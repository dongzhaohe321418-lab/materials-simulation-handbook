# Methods — High-throughput band-gap screening

The protocol has five stages. Read all of them before starting; the
data-cleaning decisions in stage 1 affect what your model can do in
stage 3.

---

## Stage 0 — Setup

Environment:

```bash
conda create -n bandgap python=3.11
conda activate bandgap
pip install mp-api pymatgen ase pandas pyarrow scikit-learn matplotlib
pip install torch torch-geometric  # match versions for your CUDA
pip install cgcnn  # there are several forks; pin one
# QE for verification:
conda install -c conda-forge qe
```

Register for a Materials Project API key at
[materialsproject.org/api](https://next-gen.materialsproject.org/api).
Export it:

```bash
export MP_API_KEY=...
```

---

## Stage 1 — Query and clean the Materials Project oxide subset

### 1.1 Query

Use `mp-api` (the modern client):

```python
from mp_api.client import MPRester

with MPRester() as mpr:
    docs = mpr.materials.summary.search(
        elements=["O"],
        num_elements=(2, 4),
        nsites_max=30,
        fields=["material_id", "formula_pretty", "structure",
                "band_gap", "energy_above_hull", "is_stable",
                "is_metal", "symmetry"],
    )
```

This returns roughly 30 000–50 000 oxides depending on the MP snapshot
date. Cache to a Parquet file *immediately*:

```python
import pandas as pd
rows = [{
    "material_id": d.material_id,
    "formula": d.formula_pretty,
    "band_gap": d.band_gap,
    "energy_above_hull": d.energy_above_hull,
    "is_stable": d.is_stable,
    "spacegroup": d.symmetry.symbol,
    "nsites": len(d.structure),
    "structure_json": d.structure.to_json(),
} for d in docs]
df = pd.DataFrame(rows)
df.to_parquet("data/mp_oxides_raw.parquet")
```

The structure-as-JSON is a slight inefficiency but lets you re-load
without re-querying.

### 1.2 Clean and filter

Apply these filters in this order:

1. Remove entries with `band_gap` missing or `< 0` (data errors).
2. Remove entries with `nsites > 30` (training cost; you can add more
   later).
3. Remove entries with `energy_above_hull > 0.2 eV/atom` (unlikely
   synthesisable; reduces noise).
4. Deduplicate by (reduced formula, space group) — keeping the entry
   with lowest `energy_above_hull`.
5. Remove obvious magnetic-noisy systems: anything containing Fe, Mn,
   Co, Ni, Cr with `band_gap == 0` *and* `is_metal == False` (these
   are likely flagged-as-metal-by-DFT but should be insulators; either
   keep with caution or drop).

After cleaning, expect ≈ 5000–8000 entries. Save as
`data/mp_oxides_clean.parquet`.

### 1.3 Train/val/test split

Use a *structure-aware* split:

- Group by reduced formula.
- Random-assign each formula to one of `train`, `val`, or `test` in
  proportion 80:10:10.
- This means that all polymorphs of, e.g., TiO$_2$ end up in the same
  partition.

Save the splits as `data/{train,val,test}.json` (lists of MP IDs).

---

## Stage 2 — Train the CGCNN

### 2.1 Choose your CGCNN implementation

Recommended:

- The original `cgcnn` (Xie 2018) — simple, readable, runs without
  PyTorch Geometric. https://github.com/txie-93/cgcnn
- The `matminer` or `aviary` re-implementations if you want a more
  modern setup.

For an undergraduate project, the original CGCNN is preferred for
clarity. The "fork to pin" is up to you; document the commit hash.

### 2.2 Convert your data

CGCNN's data format is:

- A directory of CIF files (one per material).
- A `id_prop.csv` mapping `material_id, band_gap_in_eV`.
- An `atom_init.json` with one-hot atom features.

Use the `pymatgen` `Structure.to(filename=...)` writer to dump CIF
files, then write the csv.

### 2.3 Hyperparameters

For an initial model:

| Parameter | Value |
| --- | --- |
| atom_fea_len | 64 |
| h_fea_len | 128 |
| n_conv | 3 |
| n_h | 1 |
| radial cutoff | 8 Å |
| max_num_nbr | 12 |
| batch_size | 256 |
| epochs | 60 |
| lr | 0.01 (Adam) |
| step_size | 20 (LR decay) |

### 2.4 Train

```bash
python cgcnn/main.py data/train --val-ratio 0.1 --test-ratio 0.0 \
  --epochs 60 --batch-size 256 --lr 0.01 \
  --workers 4 --train-size 4000
```

Or use the supplied wrapper in `starter-code.md`. Monitor the
val-loss curve; expect 0.6–0.8 eV MAE after 10 epochs and 0.4–0.5 eV
MAE at convergence.

### 2.5 Ensemble

Train 5 models with seeds `{0, 1, 2, 3, 4}`. Save each checkpoint.
The ensemble mean is your prediction; the ensemble standard deviation
is your uncertainty.

### 2.6 Validation

On the held-out test set, compute:

- Overall MAE.
- MAE stratified by gap range: [0, 0.5), [0.5, 1.0), [1.0, 2.0),
  [2.0, 3.0), [3.0, ∞) eV.
- Parity plot of predicted vs MP gap.
- For the ensemble, the calibration plot (predicted standard
  deviation vs absolute error). A well-calibrated ensemble has these
  scaled proportionally.

Target: overall MAE < 0.5 eV. If stratified MAE in [1.5, 2.5] eV is
> 0.7 eV, your model is not useful for screening in this window;
consider augmenting the training set with more samples in that range
or moving to MEGNet/ALIGNN.

---

## Stage 3 — Build the candidate screen set

### 3.1 Where do candidates come from?

Three options, in increasing exoticism:

1. **MP itself, but a different subset.** Query for all oxides with
   `nsites ≤ 50` (a slightly larger upper bound), remove your training
   set, predict. This identifies "known but not yet flagged"
   candidates.
2. **Hypothetical structures from prototype-substitution.** Take a
   small set of known oxide prototype structures (e.g., perovskite
   ABO$_3$, spinel AB$_2$O$_4$, wurtzite AO) and substitute the cation
   sites with a fresh combinatorial set of elements. Use `pymatgen`'s
   `Substitutor` or hand-coded substitution.
3. **OQMD or AFLOW.** Other DFT databases with oxides not in MP.
   Requires more glue code.

For this project, recommend (1) + a *small* (1000-candidate)
substitution-based extension (2) for diversity. Total target: 20 000
candidates.

### 3.2 Predict

Run the ensemble on the screen set. Record for each candidate:

- Ensemble mean predicted gap.
- Ensemble standard deviation.
- Formula, spacegroup, source (MP-not-in-training vs prototype).
- If available, `energy_above_hull` (for MP candidates) or a quick
  surrogate-model stability estimate.

Save as `screen/candidates.parquet`.

---

## Stage 4 — Rank and shortlist

### 4.1 Filter

Apply, in order:

1. Predicted gap in $[1.3, 2.7]$ eV (slightly wider than the target
   window to absorb model error).
2. Predicted uncertainty $\sigma < 0.3$ eV.
3. Stability: `energy_above_hull < 0.1 eV/atom` where available; for
   prototype substitutions where it is not, use a Magpie-feature
   stability classifier or skip and flag.
4. Element-availability filter: remove rare-earth-only compositions
   if your study targets earth-abundant photocatalysts.

After filtering, expect to have 50–500 candidates.

### 4.2 Rank

Rank by *predicted distance to the target gap centre* (2.0 eV)
weighted by predicted uncertainty:

$$
\text{score} = (\mu - 2.0)^2 + \lambda \sigma^2.
$$

with $\lambda \approx 1$. The intuition: high-uncertainty candidates
are penalised even if their mean prediction looks good. Other
acquisition-style functions (UCB-style $\mu - \kappa \sigma$) are
equally defensible — pick one and explain your choice.

### 4.3 Shortlist

Take the top 10. Inspect each one manually: is the formula
plausible? Is the space group reasonable? Has anyone studied it
already? A 30-minute Google Scholar check per candidate is required
diligence.

Save as `shortlist/top10.csv` with columns:

```
mp_id, formula, spacegroup, mu_eV, sigma_eV, e_above_hull, source, notes
```

---

## Stage 5 — DFT verification

For at least 3 of the top 10:

### 5.1 Relax

Take the structure (from MP, or build from the prototype). Relax with
PBE at modest **k**-mesh (e.g., $3 \times 3 \times 3$ for a 20-atom
cell) and `ecutwfc = 60` Ry, `ecutrho = 480` Ry, `vc-relax`
calculation (`cell_dofree = 'all'`).

### 5.2 Self-consistent SCF

After relaxation, run a converged SCF with denser **k**-mesh and a
finer convergence threshold.

### 5.3 Band-structure or NSCF on dense mesh

Either:

- A non-self-consistent calculation on a dense uniform mesh
  (e.g., $8 \times 8 \times 8$), then extract the band gap from
  `bands.x` or directly from the eigenvalues, or
- A high-symmetry **k**-path band-structure calculation.

The first is faster and is what MP itself does for the `band_gap`
field. Use it for direct comparison.

### 5.4 (Optional) HSE06

If you have a VASP licence and a GPU partition, run an HSE06
band-gap calculation at the relaxed PBE structure. Expect the HSE
gap to be 0.5–1.5 eV larger than the PBE gap.

### 5.5 Report

For each verified candidate, report:

- CGCNN predicted gap ± uncertainty.
- PBE relaxed gap (your own DFT).
- (Optional) HSE06 gap.
- Discrepancy CGCNN vs PBE. Discuss whether it exceeds the
  reported CGCNN uncertainty.

---

## Pitfalls

1. **Including unrelaxed MP structures in training.** All MP
   structures *are* relaxed; do not "re-relax" before training. But
   if you build prototype-substituted candidates, you must relax
   them before predicting (different relaxation level = different
   feature distribution).
2. **Reusing the MP band gap as ground truth for prototype
   candidates.** Prototype candidates do not have an MP gap; treating
   them as "ground truth = 0" is a leakage error in stage 4. Filter
   them in stage 3.
3. **Forgetting magnetic ordering.** MP's band gap is computed at
   the *ground-state* spin configuration MP found, which for some
   oxides is non-trivial. If your CGCNN predicts a gap for a candidate
   for which the ground-state spin is ambiguous, the prediction is
   noisier than the ensemble standard deviation suggests.
4. **Top-K selection bias.** Picking the top 10 by predicted gap
   amplifies model bias: your shortlist will be enriched in
   compositions the model was overconfident on. The uncertainty
   filter mitigates this; the manual literature check completes it.
5. **Counting "novel" compounds without checking ICSD.** A compound
   not in MP is not necessarily novel; ICSD may already have it. A
   thorough study cross-checks ICSD before claiming novelty.
6. **HSE06 verification without lattice re-relaxation.** HSE06
   energies are typically computed at the PBE-relaxed structure; you
   do not relax at HSE06. State this clearly in your report.

---

## What "done" looks like

You have:

- A trained ensemble of CGCNNs with stratified MAE reported.
- A ranked, filtered shortlist of 10 candidates.
- DFT verification of at least 3, with the discrepancy honestly
  reported.
- A discussion of how many shortlisted candidates are novel and how
  many are already photocatalyst-literature staples.
