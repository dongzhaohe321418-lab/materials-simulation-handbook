# Methods — Training an MLIP from scratch

The protocol has six phases:

1. Baseline (does a foundation model already work? if so, fine-tune
   instead).
2. Sampling strategy.
3. Initial training (MACE v1).
4. Validation suite.
5. Active-learning round(s).
6. Production MD.

Read all six before starting.

---

## Phase 0 — Setup

```bash
conda create -n mlip-scratch python=3.11
conda activate mlip-scratch
conda install -c conda-forge ase pymatgen qe matplotlib numpy phonopy
pip install mace-torch lammps
```

Document your environment: `pip freeze > requirements.lock` and commit
it. MLIP libraries change quickly; pinning matters.

---

## Phase 1 — Baseline foundation-model check

Before training from scratch, confirm that an off-the-shelf MACE-MP-0
does not already solve your problem.

```python
from mace.calculators import mace_mp
calc = mace_mp(model="medium", device="cuda")
# Apply to a small set of your target structures and check forces.
```

Compute forces on:

- The equilibrium geometry of your system.
- A 10 % strained version.
- A 1500-K-rattled version.

If MACE-MP-0 gives MAE < 100 meV/Å on all three, fine-tune it for
your application rather than training from scratch — fine-tuning is
much cheaper. Document the baseline metric in your report.

If MACE-MP-0 fails on any of these tests (large MAE, unphysical
energies, or unstable forces), proceed to train from scratch.

---

## Phase 2 — Sampling strategy

Your training set must cover:

| Component | Why | How |
| --- | --- | --- |
| Equilibrium geometry | Anchors the PES minimum | DFT relaxation; include with high weight |
| Thermal samples (low T) | Solid-state harmonic region | AIMD at 100–300 K, 2 ps |
| Thermal samples (high T) | Anharmonic + diffusive | AIMD at 600–1500 K (or above $T_m$ if liquid is target), 2–4 ps |
| Rattled cells | Random off-equilibrium | `ase.atoms.rattle(stdev=0.05–0.15)` on relaxed structures |
| Strained cells | Stress tensor accuracy | ±2 %, ±5 % isotropic and uniaxial strain on relaxed; rattle each |
| Defects / boundary configs | If your application has these | Vacancy, swap, surface depending on application |

For a 30-atom system, this typically produces 400–600 frames. Save
all of them as an extended-XYZ file with `info["energy"]`,
`arrays["forces"]`, `info["stress"]`, and `info["config_type"]`
populated.

### Sampling density rule

For each `config_type`, ensure at least 50 frames are present. If a
single AIMD trajectory dominates with 200 correlated frames, subsample
to every 5th or 10th step to reduce correlation.

### Reserved held-out frames

Carve out 10 % of each `config_type` as a held-out test set
*before training begins*. This is not the validation set used during
training; it is the truth set for the final report.

---

## Phase 3 — Initial training (MACE v1)

Use mostly defaults; resist temptation to tune.

| Parameter | Value | Note |
| --- | --- | --- |
| r_max | 5.0 Å | Adjust to system; cover at least the second NN shell |
| max_L | 2 | |
| correlation | 3 | |
| hidden_irreps | `128x0e+128x1o` | |
| num_interactions | 2 | |
| batch_size | 4–8 | smaller cells allow larger batches |
| epochs | 250 | |
| forces_weight | 100 | |
| energy_weight | 1 | |
| stress_weight | 1 | leave on, even for small molecules |

Launch (see `starter-code.md` for the wrapper).

Train two models with different random seeds. Both will become the
committee for active learning later.

---

## Phase 4 — Validation suite

This is the most important phase. A model that passes parity is not
necessarily useful. Run each test below.

### 4.1 Parity plot

Per-atom energy and per-component force, on the held-out test set.
Target: F MAE < 5 % of typical force magnitude. Report per-element
breakdown if more than one element.

### 4.2 Spot-check on hand-picked configurations

Generate three or four configurations *manually*: a strain that you
didn't train on (e.g., 7 % shear), an unusual coordination, an
extreme bond length. Predict; compute DFT on the same; compare.
This catches the "out-of-distribution" failure mode.

### 4.3 MD stability test

Run NVE (no thermostat) MD at 300 K for 50 ps. Plot total energy
vs time. A well-trained MLIP shows drift below ≈ 0.5 meV/atom over
50 ps; drift above 5 meV/atom indicates a problem. Also plot the
minimum interatomic distance vs time — sudden drops below 1 Å
mean atom overlaps and a broken potential.

### 4.4 RDF check

Run NVT at the same temperatures as your training set. Compute the
RDF (`ase.geometry.analysis.Analysis.get_rdf` or `pyscal`). Compare
against the AIMD-derived RDF from your training trajectories. The
first-peak position should agree to ≈ 0.02 Å; peak heights to ≈ 10 %.

### 4.5 Phonon dispersion

Use `phonopy` to compute the harmonic phonon spectrum of the
equilibrium structure with both the MLIP and DFT (the DFT one is
expensive — a small cell is acceptable for this test, e.g., a 6-atom
primitive cell of MoSSe).

Compare the two spectra over the high-symmetry **k**-path. Acceptable
agreement: ≈ 5–10 % deviation in frequencies, no imaginary modes in
the MLIP that are not in the DFT.

### 4.6 Stress / elastic constants (if periodic)

Compute the bulk modulus by fitting a Birch–Murnaghan EoS with MLIP
and DFT. Compare.

---

## Phase 5 — Active learning

Almost certainly, one or more validation tests will fail. The standard
remedy is active learning.

### 5.1 Identify uncertainty during MD

Run a longer MD trajectory (say 50 ps NVT at the high training
temperature). At every 50 fs, evaluate forces using *both* MACE v1a
and v1b (the two committee members). Compute the per-atom force
standard deviation:

$$
\sigma_i = \frac{1}{\sqrt{3}} \|\mathbf{F}_i^{(a)} - \mathbf{F}_i^{(b)}\|.
$$

(A more robust formulation uses ≥ 3 committee members; 2 is the
minimum.)

Configurations with $\max_i \sigma_i$ exceeding a threshold (commonly
≈ 100 meV/Å) are flagged.

### 5.2 Pick AL candidates

From the flagged configurations, pick a *diverse* subset:

- Greedy selection by maximum $\sigma$ tends to give correlated
  picks (e.g., 20 consecutive MD frames around a single transition).
- A better strategy: cluster the flagged frames by descriptor (a
  simple Magpie or SOAP fingerprint) and pick one per cluster.

Target: 30–50 new frames per AL round.

### 5.3 Relabel with DFT

Run single-point SCF on each selected frame with the same DFT
settings as the training set. Extract energies, forces, and
stresses.

### 5.4 Retrain

Concatenate the new frames with the original training set. Retrain
MACE v2 from scratch (cheaper) or by warm-starting from MACE v1
weights (slightly faster).

### 5.5 Re-validate

Re-run the full validation suite. Document the changes — typically
parity MAE is largely unchanged (training set is dominated by the
old frames) but MD stability improves dramatically.

### 5.6 Repeat or stop

Run a second AL round if MD still fails. After 2–3 rounds, returns
diminish. Stop and report.

---

## Phase 6 — Production MD

With the final validated model:

- Run NVT at two temperatures (e.g., 300 K and 600 K).
- ≥ 100 ps per temperature.
- Compute one structural quantity (RDF, MSD, or dihedral
  distribution depending on system).
- Plot the per-step total energy. Drift should be undetectable.

If the production MD fails, return to phase 5.

---

## Pitfalls

1. **AIMD trajectory length confused for "training set size".** Each
   AIMD frame at 1-fs spacing is highly correlated with the next.
   Subsample by 10–20 fs minimum.
2. **Skipping the foundation-model baseline.** If MACE-MP-0 already
   does the job, you are wasting weeks. Always test the baseline.
3. **Tuning hyperparameters before sampling.** Spending a week tuning
   `max_L = 1, 2, 3` while your training set is missing the relevant
   configurations is pure ritual. Sample first.
4. **Validating only on parity.** A model with E MAE 2 meV/atom and
   F MAE 30 meV/Å can still produce nonsense MD if the training set
   misses the transition-state region.
5. **AL with two committee members.** Two is the minimum and gives
   noisy uncertainty estimates. If you can afford to train 3–5
   committee members, do so.
6. **Forgetting to fix random seeds.** Without fixed seeds, your AL
   "improvement" cannot be attributed to the new data — it might just
   be initialisation luck.
7. **AL on configurations that DFT can't converge.** Heavily distorted
   AL candidates sometimes fail to converge in DFT. Skip them and pick
   the next-lowest-uncertainty replacement.

---

## What "done" looks like

You have:

- A documented justification for training from scratch (vs
  fine-tuning).
- A final MACE potential that passes parity, MD stability, RDF, and
  phonon tests.
- At least one AL round with quantitative before/after numbers.
- A 100-ps production trajectory.
- An honest discussion of the remaining limitations.
