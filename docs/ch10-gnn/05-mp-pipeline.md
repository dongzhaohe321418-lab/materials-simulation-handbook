# 10.5 A Materials Project Pipeline

Sections 10.3 and 10.4 covered architectures; this section covers
*data*. The CGCNN implementation we wrote is generic — it consumes any
list of pymatgen `Structure` objects with labels — but every choice
made between the database and the dataloader can shift the apparent
accuracy by a factor of two. Random splits versus structurally disjoint
splits, deduplicated entries versus polymorph duplicates, the inclusion
of high-energy metastable phases versus only ground states: each is a
seemingly mundane decision that meaningfully changes what the trained
model is actually doing.

This section walks through a complete pipeline from scratch — register
for an API key, pull 5000 oxides, deduplicate, split honestly, train,
evaluate, plot — with commentary at each step on the easily-missed
pitfalls.

!!! note "Why this section, not a Jupyter notebook?"
    The pipeline below could be — and should be — a notebook in your own
    workflow. Read this section as the *commentary* that an experienced
    practitioner would supply at each cell of that notebook, pointing
    out the silent ways things go wrong.

## 10.5.1 Setting up access to the Materials Project

The Materials Project distributes its data through the `mp-api` Python
client. To get an API key, navigate to
[https://next-gen.materialsproject.org/api](https://next-gen.materialsproject.org/api),
log in (Google or institutional SSO), and copy the generated key. The
key is a thirty-two-character alphanumeric string. Store it as an
environment variable rather than hard-coding it:

```bash
export MP_API_KEY="paste-your-key-here"
```

In Python:

```python
import os
from mp_api.client import MPRester

API_KEY = os.environ["MP_API_KEY"]
```

!!! note "Why an environment variable?"
    Hard-coding a secret in a script that you might check into git is a
    standard route to a leaked credential. The Materials Project key
    grants access to your usage logs and rate limits; revoking and
    reissuing if leaked is a chore best avoided. Environment variables
    keep the key out of source control while still making it available
    to any process you launch from the same shell. Production
    deployments use a secrets manager (HashiCorp Vault, AWS Secrets
    Manager) for the same reason. For our purposes,
    `export MP_API_KEY=...` in your `~/.bashrc` or `~/.zshrc` is
    enough.

The client is rate-limited (around 1000 requests per minute as of
2026) and ships responses in compressed batches; you should not need
to think about either limit for the dataset sizes considered here.

## 10.5.2 Querying for oxides

We want all binary, ternary and quaternary oxides — anything containing
oxygen — with a formation energy and DFT-computed structure. The
relevant endpoint is `mpr.materials.summary.search`:

```python
from mp_api.client import MPRester

with MPRester(API_KEY) as mpr:
    docs = mpr.materials.summary.search(
        elements=["O"],
        num_elements=(2, 4),
        fields=[
            "material_id",
            "structure",
            "formula_pretty",
            "formation_energy_per_atom",
            "energy_above_hull",
            "nsites",
        ],
    )
```

This returns roughly 60 000 documents. We will cap the dataset at 5000
randomly sampled entries; for a published study you would take them
all, but the smaller subset finishes training in a reasonable time on a
laptop GPU.

!!! example "Anatomy of a query"
    The `summary.search` endpoint is the workhorse: a single call
    returns every cached scalar property and the relaxed structure for
    each matching entry. Two arguments deserve attention.

    `fields=[...]`: explicitly list only the fields you need. Each
    additional field increases response size and time. For training, we
    need the structure and the target (formation energy); for splitting
    and analysis, the pretty formula and atom count are useful too.

    `elements=["O"]` and `num_elements=(2, 4)`: combined, these say
    "contains oxygen, has between two and four distinct elements".
    Other useful filters include `is_stable=True` (only entries on the
    convex hull), `band_gap=(0.5, 3.0)` (semiconductor candidates),
    `theoretical=False` (exclude purely theoretical predictions). The
    full filter set is in the `mp-api` docs; combining filters at the
    server side is far cheaper than downloading then filtering
    locally.

    For our example we want **all** oxides for the GNN to see broadly;
    we then filter and sample after the download.

    A more advanced workflow paginates the query in chunks to avoid
    server timeouts on large result sets. `summary.search` handles this
    automatically for result sets up to $\sim 10^5$; beyond that, use
    `num_chunks` and `chunk_size` to manage memory.

```python
import numpy as np

rng = np.random.default_rng(seed=42)
keep_indices = rng.choice(len(docs), size=5000, replace=False)
docs = [docs[i] for i in keep_indices]
```

A first sanity check: how is the formation energy distributed?

```python
import matplotlib.pyplot as plt

energies = np.array([d.formation_energy_per_atom for d in docs])
plt.hist(energies, bins=50)
plt.xlabel("formation energy per atom (eV)")
plt.ylabel("count")
plt.show()
```

The histogram peaks near $-2$ eV/atom for typical metal oxides, with a
tail extending to $-4$ eV/atom for the most stable rare-earth oxides
and another tail to near zero or slightly positive for the highest-
energy metastable polymorphs. The full range is roughly $-4$ to $+1$
eV/atom; a model with MAE $0.1$ eV/atom is therefore capturing roughly
98% of the variance.

## 10.5.3 The polymorph problem

The single most common mistake in materials ML is to treat a database
as a collection of independent samples. It is not. The Materials
Project contains many *polymorphs* of the same composition — structures
with identical chemical formula but different crystal symmetry — and
many *near-duplicates* generated by small distortions of a parent
phase. Anatase and rutile TiO$_2$, $\alpha$- and $\beta$-quartz,
zincblende and wurtzite ZnS: all are listed as separate entries with
near-identical formation energies and very similar atomic
neighbourhoods.

If a random 80/10/10 split places anatase in the training set and
rutile in the test set, the model sees a structure in the test set
that is in some sense already familiar — a different polymorph of an
already-seen composition. The reported test MAE then overstates the
true generalisation, sometimes by a factor of two or three.

!!! note "Derivation: why random splits inflate accuracy"
    Model the formation energy as $y = f(\mathbf{x}) + \eta$ where
    $\mathbf{x}$ is the structure descriptor, $f$ the true mapping, and
    $\eta$ measurement-plus-DFT noise. For polymorphs $\mathbf{x}_A,
    \mathbf{x}_B$ of the same composition, the descriptors are *close*
    in feature space (similar atomic environments) and $f(\mathbf{x}_A)
    \approx f(\mathbf{x}_B)$. With a random split, $\mathbf{x}_A$ may
    be in train and $\mathbf{x}_B$ in test. The model sees a training
    example near every test point and effectively does 1-nearest-
    neighbour lookup, achieving low test error not by generalising but
    by memorising. The empirical inflation factor depends on the
    polymorph density of the database; on Materials Project oxides we
    measure roughly $2\times$. The composition-disjoint split removes
    this leakage by ensuring all polymorphs of a formula share a single
    split.

There are two principled fixes.

**Composition-disjoint split.** Group structures by chemical formula
and split *groups* rather than entries. All polymorphs of TiO$_2$ go
into the same split (train, val, or test). This avoids the polymorph
contamination but allows the model to see TiO$_2$ alongside SiO$_2$,
which is still chemically related.

**Structure-disjoint split.** Cluster all structures by similarity in
some descriptor space (a featuriser like Magpie, or the
last-layer embeddings of a pre-trained GNN), then split by cluster.
This is the most honest but also the most laborious; the field has not
agreed on a single recipe.

For this section we use the composition-disjoint split.

```python
from collections import defaultdict

formula_to_docs: dict[str, list] = defaultdict(list)
for d in docs:
    formula_to_docs[d.formula_pretty].append(d)

formulae = list(formula_to_docs.keys())
rng.shuffle(formulae)

n_train_groups = int(0.8 * len(formulae))
n_val_groups = int(0.1 * len(formulae))

train_formulae = formulae[:n_train_groups]
val_formulae = formulae[n_train_groups:n_train_groups + n_val_groups]
test_formulae = formulae[n_train_groups + n_val_groups:]

def gather(group: list[str]) -> list[StructureRecord]:
    out: list[StructureRecord] = []
    for f in group:
        for d in formula_to_docs[f]:
            out.append(StructureRecord(d.structure, d.formation_energy_per_atom))
    return out

train_records = gather(train_formulae)
val_records = gather(val_formulae)
test_records = gather(test_formulae)
```

To make the cost of doing it the wrong way concrete, run the same code
twice — once with a random split and once with this composition-
disjoint split — and compare the test MAEs. The random split typically
reports $\sim 0.03$ eV/atom, while the composition-disjoint split
reports $\sim 0.06$ eV/atom on this dataset. The model has not
become worse; the evaluation has become honest.

## 10.5.4 Deduplication

Even within a single chemical formula, the database occasionally
contains exact-duplicate or near-duplicate structures generated by
different calculations of the same phase. Deduplication is a separate
operation from splitting.

The pymatgen `StructureMatcher` is the standard tool. It checks for
structural equivalence under arbitrary lattice transformations and
small atomic displacements:

```python
from pymatgen.analysis.structure_matcher import StructureMatcher

matcher = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5.0)

def deduplicate(records: list[StructureRecord]) -> list[StructureRecord]:
    kept: list[StructureRecord] = []
    for rec in records:
        if any(matcher.fit(rec.structure, k.structure) for k in kept):
            continue
        kept.append(rec)
    return kept
```

!!! note "What 'the same structure' means to `StructureMatcher`"
    `StructureMatcher` considers two structures equivalent if there
    exists a transformation of one onto the other that:

    - Permutes atom labels (so reordering does not matter).
    - Applies a rigid rotation, reflection, or lattice transformation
      (so different lattice settings of the same crystal are
      equivalent).
    - Displaces each atom by no more than `stol` fractional units (so
      small DFT-relaxation jitter does not produce false negatives).
    - Distorts the lattice parameters by no more than `ltol` relative
      and angles by no more than `angle_tol` degrees.

    The default tolerances (`ltol=0.2, stol=0.3, angle_tol=5`) are
    reasonable for DFT-relaxed structures; tighter tolerances reject
    legitimate duplicates, looser ones merge distinct polymorphs.
    `StructureMatcher.fit` returns boolean equivalence; `fit_with_electrostatics`
    is a faster shortcut that pre-screens by composition and density
    before running the full match.

`StructureMatcher` is $O(n^2)$ in the number of structures and slow,
so deduplicate *within* each composition group rather than globally:

```python
deduplicated = []
for f in formula_to_docs:
    deduplicated.extend(deduplicate(
        [StructureRecord(d.structure, d.formation_energy_per_atom)
         for d in formula_to_docs[f]]
    ))
```

On our 5000-entry sample this removes roughly 200 near-duplicates.

## 10.5.5 Training

With train, val and test sets in hand, training is identical to the
fifty-structure demonstration in §10.3.4 — the same `train_cgcnn`
function call, more data, more epochs:

```python
model = train_cgcnn(
    train_records=train_records,
    val_records=val_records,
    n_epochs=200,
    batch_size=64,
    lr=5e-4,
)
```

On a single mid-range GPU (RTX 3060 / A4000 class) this finishes in
roughly two hours. On CPU, plan on a day. Training MAE drops to about
0.03 eV/atom; validation MAE plateaus near 0.06 eV/atom under the
composition-disjoint split.

A more careful protocol adds:

- **Cosine learning-rate schedule.** `torch.optim.lr_scheduler.CosineAnnealingLR`
  over the full training run, decaying from $5 \times 10^{-4}$ to
  $10^{-6}$.
- **Validation-based checkpointing.** Save the model state with the
  lowest validation loss; restore it for test-set evaluation.
- **Early stopping.** If validation loss has not improved for 30
  epochs, terminate.
- **Gradient clipping.** Clip gradient norms to 1.0 to avoid instabilities
  with sparse high-charge species.

These are standard ML hygiene; we omitted them in §10.3 for brevity but
the production training script in the supplementary code includes them.

## 10.5.6 Evaluation: parity and MAE breakdown

After training, the test-set evaluation is the moment of truth.

```python
test_ds = CrystalGraphDataset(test_records)
test_loader = DataLoader(test_ds, batch_size=64)

model.eval()
preds: list[float] = []
targets: list[float] = []
with torch.no_grad():
    for batch in test_loader:
        batch = batch.to(next(model.parameters()).device)
        preds.extend(model(batch).cpu().numpy().tolist())
        targets.extend(batch.y.view(-1).cpu().numpy().tolist())

preds_arr = np.array(preds)
targets_arr = np.array(targets)

mae = float(np.mean(np.abs(preds_arr - targets_arr)))
rmse = float(np.sqrt(np.mean((preds_arr - targets_arr) ** 2)))
print(f"test MAE  = {mae:.4f} eV/atom")
print(f"test RMSE = {rmse:.4f} eV/atom")
```

The parity plot — predicted versus true, with the $y = x$ line —
visualises systematic biases that the scalar MAE hides.

```python
fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(targets_arr, preds_arr, s=8, alpha=0.5)
lo, hi = targets_arr.min(), targets_arr.max()
ax.plot([lo, hi], [lo, hi], "k--", lw=1)
ax.set_xlabel("DFT formation energy (eV/atom)")
ax.set_ylabel("CGCNN prediction (eV/atom)")
ax.set_aspect("equal")
plt.tight_layout()
plt.show()
```

Look for two failure modes:

- **Curvature.** A predicted-versus-true relationship that bows away
  from the diagonal indicates a systematic bias — typically over-
  prediction at the extremes. Caused by mean-regression: the model
  defaults to the training mean when uncertain.
- **Outlier classes.** A cluster of points far from the diagonal
  often corresponds to a chemistry that the model has not seen
  enough of — rare lanthanide oxides, for instance. Inspect them by
  formula and decide whether the model needs more such examples or
  whether they should be excluded from the application domain.

## 10.5.7 Caveats and where to go from here

The pipeline above is honest enough to be useful but still simplified.
For a serious campaign:

**Out-of-domain test sets.** The composition-disjoint split is a
reasonable approximation to "what would the model predict on a new
chemistry?", but it leaves alkali metal oxides closer to each other
than to transition-metal oxides. A *true* out-of-domain test might hold
out, say, all lanthanide-containing oxides as the test set, and
measure how badly the model degrades. Such tests routinely reveal
errors three to ten times worse than the in-distribution MAE.

**Energy-above-hull versus formation energy.** Formation energy is the
default target; energy-above-hull (distance from the convex stability
hull) is the more useful prediction for screening but is much harder
to learn because it can be a *small difference of large numbers*. A
model that achieves 0.05 eV/atom on formation energy might still rank
polymorphs incorrectly because the relevant differences are smaller
than the model error.

**Uncertainty estimation.** A point prediction is rarely enough for
decision-making. Standard methods (deep ensembles, MC dropout, last-
layer Bayesian heads, or Gaussian-process surrogates over CGCNN
embeddings) attach a per-prediction uncertainty. Chapter 11 will
develop these for Bayesian optimisation; for property regression they
are equally vital.

**Foundation model as backbone.** As foreshadowed in §10.4 and detailed
in Chapter 12, the modern alternative to training CGCNN from scratch is
to load pre-trained M3GNet or CHGNet weights, freeze the backbone,
attach a small property-regression head, and fine-tune on the labelled
oxide dataset. With 5000 examples this often reaches lower MAE than
training a CGCNN from scratch, and the inductive bias of the
foundation model regularises against the small-data overfitting.

!!! tip "Forward reference"
    Chapter 11 uses the CGCNN pipeline of this section as the *cheap
    oracle* in a Bayesian optimisation campaign: predict formation
    energies for a million candidates with the trained CGCNN; run DFT
    on the top hundred selected by an acquisition function. The
    composition-disjoint split is what makes the CGCNN's uncertainty
    estimates trustworthy in that downstream loop.

!!! note "A common follow-up: predicting band gaps and beyond"
    The same pipeline retargets to band gaps by replacing
    `formation_energy_per_atom` with `band_gap` in the query and the
    dataset constructor. Two practical differences emerge.

    First, band gaps are zero for all metallic compositions, producing a
    long zero-spike in the histogram. The model can easily learn to
    output zero — high accuracy on metals but useless on semiconductors.
    The fix is to split metals and non-metals into two heads (a
    classifier for metal/non-metal and a regressor for the band gap of
    non-metals) or to train two separate models.

    Second, the band gap from DFT systematically underestimates the
    experimental value by roughly 30–50%. A model trained on DFT band
    gaps therefore predicts DFT band gaps, not measured ones. Match
    the training labels to the downstream use; if you want experimental
    band gaps, train on the (much smaller) experimental database such
    as MatBench's `expt_band_gap` task.

The pipeline you have just built, however, is the right ground truth.
Before reaching for any foundation model, you should have run *some*
version of this pipeline from scratch — to understand exactly what data
flows where and to develop an intuition for which numbers are
trustworthy. With that experience in hand, Chapter 11's question
becomes much sharper: given a fast property predictor and an
uncertainty estimate, how should we *use* them to guide the next
experiment?

## 10.5.7a Uncertainty and ensembles

A single CGCNN gives a *point* prediction. For downstream applications
(Bayesian optimisation, decision-making under uncertainty) we need a
*distribution* over predictions. The cheapest reliable route is a *deep
ensemble*: train $M$ CGCNNs from different random initialisations on
the same data, and report the ensemble mean as the prediction and the
ensemble standard deviation as the uncertainty.

```python
# Train M independently initialised CGCNNs.
ensemble = []
for seed in range(5):
    torch.manual_seed(seed)
    model = train_cgcnn(train_records, val_records, n_epochs=200)
    ensemble.append(model)

def ensemble_predict(model_list, data):
    preds = []
    for m in model_list:
        m.eval()
        with torch.no_grad():
            preds.append(m(data).cpu().numpy())
    preds = np.stack(preds)         # (M, batch)
    return preds.mean(0), preds.std(0)
```

The ensemble standard deviation tracks the *epistemic* uncertainty —
the variation in predictions arising from finite training data and
random initialisation. It is well-calibrated in practice: structures
with high ensemble disagreement turn out to have high *actual* errors.
The cost is $M\times$ training compute; with $M = 5$ on a laptop GPU,
roughly ten hours.

For BO in Chapter 11 we will use exactly this construction: a CGCNN
ensemble providing $(\mu, \sigma)$ at every candidate, fed into an
acquisition function.

## 10.5.7b Multi-task learning

A single CGCNN trained on formation energy *only* discards information
that other property labels carry. Multi-task learning trains one
network on several properties simultaneously, with a shared backbone
and per-property heads. The intuition is that the chemistry that
governs formation energy also governs band gap, bulk modulus and many
others; learning all of them jointly regularises the backbone and
improves each individual prediction.

For Materials Project, a multi-task CGCNN with heads for
$\{E_\text{form}, E_\text{hull}, \text{band gap}, \text{bulk modulus}, V_\text{atom}\}$
reaches lower MAE on each individual target than five single-task
models — typically by 5–15%. The loss is a weighted sum of per-task
losses; the weights matter, and standard schemes (equal weighting,
inverse-variance weighting, learned weighting via Kendall et al. 2018's
homoscedastic uncertainty trick) all work.

The relevant change to the §10.3 implementation is to replace the
single output head with $K$ heads:

```python
self.heads = nn.ModuleList([
    nn.Sequential(nn.Linear(atom_dim, hidden_dim), nn.Softplus(),
                  nn.Linear(hidden_dim, 1))
    for _ in range(K)
])
```

and to compute the multi-task loss as $\sum_k w_k \cdot \mathrm{MAE}_k$.

## 10.5.8 A reproducible script

For completeness, the full pipeline from this section — query,
deduplicate, split, train, evaluate, plot — fits in a single
two-hundred-line script. The exercises distribute it in pieces; the
reader who wants a single-file working example will find it in the
chapter's code repository as `mp_pipeline.py`. Run it once with the
random split, once with the composition-disjoint split, and the
roughly twofold inflation of accuracy on the random split will make
the polymorph-leakage point indelible.

!!! note "Sanity checklist before reporting a number"
    Before quoting a CGCNN MAE to a colleague or paper, verify:

    1. The train/val/test split is composition-disjoint (or
       structure-disjoint), not random.
    2. The test set was held out from *every* preprocessing step
       (normalisation statistics computed on train only, not on
       train + test).
    3. The reported MAE is the *test*-set value at the *best
       validation* checkpoint, not the final-epoch value, and not the
       validation MAE itself.
    4. The Materials Project database version is recorded (or commit
       hash of your local CIF cache).
    5. The numbers are averaged over $\geq 3$ random seeds, with a
       reported standard deviation.

    Skipping any of these is the difference between a number you can
    defend at a conference and a number that quietly drifts up by 30%
    when somebody re-runs the script. Chapter 12 will return to this
    "reproducibility ethic" in the context of foundation models, where
    the surface area for hidden contamination is even larger.

### Section summary

- The Materials Project pipeline is a query → deduplicate → split →
  train → evaluate workflow; each step has a silent failure mode
  worth understanding.
- Random splits *inflate* test-set accuracy by a factor of $\sim 2$
  on polymorph-rich databases; composition-disjoint splits are the
  honest baseline.
- `StructureMatcher` is the standard tool for deduplication; tune
  `ltol, stol, angle_tol` to match your tolerance for "the same
  structure".
- Deep ensembles ($M = 5$) give calibrated uncertainty estimates for
  use in BO downstream (Chapter 11).
- Multi-task learning (joint heads for several MP properties)
  regularises the backbone and lowers per-task MAE by 5–15%.
