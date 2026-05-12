# Project 4 — Training an MLIP from scratch on a novel system

## Research question

Pick a small, well-defined materials system that is *not* covered well
by existing pre-trained MLIPs (a 2D material, a Janus dichalcogenide,
a small peptide, a metal–organic-framework fragment, or similar).
Generate your own DFT training data, train a MACE potential, run a
validation suite, and use the trained MLIP for a short (≥ 100 ps) MD
trajectory. Document the failure modes you encounter and the active-
learning rounds you needed to fix them.

This is the "deep dive" project: where Project 2 was a clean recipe
for a single element, this one asks you to act as a real practitioner
making — and recording — every decision.

## Why this project

Most MLIP papers report a clean training-set construction, a clean
parity plot, and a clean MD trajectory. The reality of training an
MLIP from scratch on a system you have not seen before is messier:

- Your first random sampling will miss configurations the MD visits.
- Your first model will explode 5 ps into the production MD.
- Your phonon dispersion will have imaginary modes near the zone
  boundary.

These are the experiences that teach you what an MLIP *is*, not what
the marketing literature says it is. The project is designed to give
you those experiences in a contained, eight-to-ten-week window.

## Suggested target systems

Pick **one**:

1. **A Janus dichalcogenide MoSSe monolayer.** A 2D semiconductor
   with a mirror-symmetry-broken sandwich (S–Mo–Se rather than
   S–Mo–S or Se–Mo–Se). Cell of 6 atoms (1×1×1) or 24 atoms (2×2×1
   for sampling). Few existing MLIPs cover this exactly.
2. **A 2D borophene polymorph (e.g., β₁₂-borophene).** A flat boron
   sheet with a non-trivial honeycomb topology. Sub-stochiometric
   boron sheets are not in MP. 24-atom cell.
3. **A small peptide (Ala–Pro–Ala) in vacuum.** ~30 atoms. Periodic
   PBC with a large vacuum buffer.
4. **A Ti–Zr–Hf MAX-phase precursor (e.g., Ti$_2$AlC alloyed with
   Zr).** Substitutional alloy chemistry; MP coverage is sparse for
   substitutional disorder.

If you have a research lab connection, pick a system that aligns
with their work. The methodology is identical.

## Expected outcomes

By the end of the project you will deliver:

- A trained MACE potential with documented force MAE on a held-out
  test set. Target: force MAE within 5 % of the typical force
  magnitude in the dataset (typically < 80 meV/Å for a covalently
  bonded system, < 50 meV/Å for a metallic one).
- A validation suite consisting of, at minimum:
  - Energy parity plot.
  - Force parity plot (separately, per element if more than one).
  - MD stability test: a ≥ 100 ps NVE/NVT trajectory at two
    temperatures, verifying no energy drift and no atom-overlap
    explosions.
  - Radial distribution function (RDF) comparison against the
    AIMD reference.
  - Phonon spectrum at $\Gamma$ (or, for periodic systems, the full
    dispersion along a high-symmetry path). The MLIP phonons should
    agree with DFT phonons to within ≈ 10 % over the spectrum, and
    no imaginary modes should appear that are not also in the DFT.
- An active-learning log: at least one round in which you identified
  high-uncertainty configurations during MD, re-labelled with DFT,
  re-trained, and demonstrated a reduction in the MD failure rate.
- A 4–6 page report describing the system, the choices made, the
  failures observed, and the final validated potential.

## Time estimate

Ten weeks, part-time.

| Week | Activity |
| --- | --- |
| 1 | Background reading. Pick the system. Compute DFT lattice parameter and basic properties. |
| 2 | Initial data sampling: short AIMD at two temperatures + rattled cells + strained cells. ≈ 300 frames. |
| 3 | Train MACE v1. Compute parity. |
| 4 | Run first production-like MD with MACE v1. *Expect failures.* Identify what went wrong. |
| 5 | Active-learning round 1: identify high-uncertainty frames (committee of 2 MACEs), re-label with DFT, augment training set. |
| 6 | Train MACE v2. Compare validation metrics. |
| 7 | Repeat AL if needed. Phonon spectrum check. |
| 8 | Final production MD: ≥ 100 ps at two temperatures. RDF check. |
| 9 | Write-up. |
| 10 | Buffer for cluster outages, model debugging, and revisions. |

## Compute budget

Variable, depending on system size. Plan for:

- ≈ 1500–3000 CPU-hours for DFT data generation, distributed across
  AIMD runs + AL relabelling. A 30-atom cell at 100 fs/step on 4
  cores is roughly 1 hour per ps.
- ≈ 30 GPU-hours for training MACE (multiple rounds).
- ≈ 200 GPU-hours for production MD if using MACE+LAMMPS at 30 fs
  per simulated ps on a single GPU.

The 10-week timeline is genuinely tight; if your cluster queues are
slow, plan to push the buffer week.

## Prerequisites

- [Chapter 6 — Running DFT in practice](../../ch06-running-dft/index.md)
  — you must be comfortable running AIMD in QE and extracting
  forces with full accuracy.
- [Chapter 9 — Machine-learning interatomic potentials](../../ch09-mlip/index.md)
  — the architectural rationale for MACE; the role of body order and
  cutoff radius; the energy/force-weight trade-off in training.
- [Chapter 8 — Statistical mechanics for simulation](../../ch08-statmech/index.md)
  — the framework for an MD-derived property (RDF, MSD, phonons) as
  a validation target.

You should have completed Project 2 *or* have equivalent experience
with at least one MLIP training pipeline. This is not a first-MLIP
project.

## What "covered well" means

A pre-trained MLIP "covers" your system if:

- Its training set includes elements present in your system at
  similar coordination environments.
- The reported test-set MAE on your system's chemistry is reasonable
  (typically < 100 meV/Å for a covalently bonded material).
- The MLIP gives stable MD at relevant temperatures.

A foundation model like MACE-MP-0 covers a substantial fraction of
crystalline materials. The point of this project is to pick a system
where it *does not*. Verify this by attempting to run MACE-MP-0 on
your target — if it gives unreasonable forces or unstable MD, you
have the right system. Document this baseline failure as the
justification for training from scratch.

## Pitfalls flagged up front

1. **Picking too big a system.** Forty atoms is generous; sixty atoms
   is challenging; a hundred atoms is a research project, not an
   undergraduate one. A six-atom cell is fine if the chemistry is
   the focus.
2. **AIMD trajectory length.** Two picoseconds at one temperature is
   not "training data"; it is a single statistically-dependent block.
   You need multiple short trajectories at different temperatures,
   plus rattled cells, plus strained cells.
3. **Force-MAE vs MD stability.** A force MAE of 100 meV/Å on a
   parity test set does not guarantee MD stability. Some of the
   worst-error configurations may not be in your test set.
4. **Phonon imaginary modes.** Small imaginary modes ($< 5$ cm$^{-1}$)
   near $\Gamma$ are usually numerical artefacts of the finite-
   difference Hessian. Larger imaginary modes elsewhere point to a
   genuinely deficient potential — usually because the relevant
   strain or twist configuration is absent from training.
5. **Cherry-picking the validation set.** It is tempting to drop the
   worst few test-set predictions as "outliers". Don't. Outliers are
   informative; they are the configurations the model is wrong on.
6. **Active-learning hubris.** A single AL round usually helps. A
   tenth AL round on the same system usually doesn't. Know when to
   stop and report the limitation.

## Deliverables checklist

- [ ] `data/dft/` — extended-XYZ training, val, test datasets, all
      versions.
- [ ] `models/v{1,2,...}/` — MACE checkpoints from each iteration.
- [ ] `validation/` — parity plots, RDF, phonon, MD-stability logs.
- [ ] `active-learning/round{1,...}/` — uncertainty-flagged frames,
      DFT relabelling, retraining diff.
- [ ] `production/md.traj` — final long-MD trajectory and analysis.
- [ ] `report.pdf`.
