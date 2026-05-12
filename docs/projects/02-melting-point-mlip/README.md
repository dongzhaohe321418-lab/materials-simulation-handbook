# Project 2 — Melting point of copper via MLIP-driven MD

## Research question

Can you train a machine-learning interatomic potential (MLIP) on
density-functional-theory data for copper and use it to compute the
melting temperature $T_m$ to within ≈ 100 K of the experimental value
(1357.77 K, sometimes rounded to 1358 K)? You will compare the MLIP's
prediction against the experimental value and, optionally, against a
classical embedded-atom-method (EAM) potential as a baseline.

This is the prototypical "train then deploy" MLIP project. It exposes
you to the entire pipeline: DFT data generation, force/energy
parity validation, MD equilibration, and a thermodynamic property
calculation (a phase transition) that no fitting target ever sees
directly.

## Why this project

Copper is the friendliest metal in the periodic table for this
exercise:

- One element. No alloy chemistry to worry about.
- Closed-shell d-electrons; spin-unpolarised DFT works. No magnetism.
- A well-tested classical EAM potential (Foiles 1986) exists, so you
  have a sanity-check baseline.
- Cubic-close-packed solid (fcc). Only one ordered phase up to melting.
- Melting temperature is far from any structural phase transition, so
  you do not need to worry about polymorphism.

The melting transition is not in your training set — DFT-MD trajectories
at 300 K, 800 K, and 1500 K sample the solid, the warm solid, and the
hot liquid respectively, but no two-phase coexistence is ever fitted.
Predicting $T_m$ is a genuine extrapolation test.

## Expected outcomes

By the end of the project you will deliver:

- A trained MACE potential for Cu with reported parity errors
  for energies, forces, and stresses on a held-out test set. Target
  metrics: energy MAE < 5 meV/atom, force MAE < 50 meV/Å, stress MAE
  < 0.5 kbar.
- A melting-point estimate from the two-phase coexistence method,
  with statistical uncertainty estimated by independent runs at the
  same temperature.
- A "hysteresis bracket" estimate from heating and cooling at constant
  rate, used as a complementary (and looser) bound.
- A comparison against the EAM-Foiles potential and a brief discussion
  of *why* MLIPs typically beat EAMs for this task — namely that EAMs
  underestimate the liquid-state cohesive energy.
- A short report (4–6 pages) presenting the pipeline, training curves,
  parity plots, and the $T_m$ estimate with its error bar.

## Time estimate

Eight weeks, part-time (≈ 15 hours per week).

| Week | Activity |
| --- | --- |
| 1 | Background reading. Set up QE/VASP and MACE environments. Verify EAM-Foiles via LAMMPS on a small Cu cell. |
| 2 | Phase 1: generate DFT training data via short AIMD at 300, 800, 1500 K and on strained cells. |
| 3 | Phase 2: train MACE on the generated set. Compute parity errors. Tune hyperparameters once. |
| 4 | Phase 3: equilibrate solid Cu at 1100 K, 1300 K, 1500 K and liquid Cu at 1700 K under NPT with the MLIP. |
| 5 | Phase 4: two-phase coexistence cells; run NPT at 1300, 1350, 1400, 1450 K. |
| 6 | Identify $T_m$ from the temperature at which the solid–liquid interface neither advances nor retreats. |
| 7 | Repeat the coexistence run at the bracketed $T_m$ with two more random seeds for an uncertainty estimate. Hysteresis comparison. |
| 8 | Write-up. |

## Compute budget

You will need:

- Approximately 800 CPU-hours for the DFT-MD data generation (10 ps
  total at 32-atom Cu cell, or equivalent). Plan for a small partition
  on a university cluster.
- A single modern GPU (e.g., A100, RTX 4090, or even an L4) for MACE
  training. Training a small MACE on ≈ 500 frames takes 6–12 hours.
- Approximately 200 GPU-hours for the production MD. Each
  coexistence run is a 6000-atom NPT trajectory of ≈ 200 ps; on a
  single GPU with MACE this is ~ 12 hours. You will run six of these.

If you do not have GPU access, you can run MACE on CPU; expect MD to be
≈ 30× slower. Plan accordingly.

## Prerequisites

- [Chapter 6 — Running DFT in practice](../../ch06-running-dft/index.md)
  for ab initio molecular dynamics in QE. You should know what
  `nstep`, `dt`, and `tempw` do.
- [Chapter 9 — Machine-learning interatomic potentials](../../ch09-mlip/index.md)
  for the message-passing formalism, the role of body order and the
  receptive field, and the meaning of `max_L` in MACE.
- [Chapter 8 — Statistical mechanics for simulation](../../ch08-statmech/index.md)
  for the two-phase coexistence method, the limitations of the
  hysteresis bracket, and Gibbs–Duhem corrections.

You should already be comfortable with LAMMPS or ASE-driven MD, with
a small (32-atom) AIMD calculation, and with submitting GPU jobs to
SLURM.

## Software you will install

- Quantum ESPRESSO 7.2 (or VASP if you have a licence) for the DFT
  training data.
- `mace-torch` (the PyTorch implementation of MACE; pip-installable).
- LAMMPS with the MACE–LAMMPS interface for production MD.
- `ase`, `pymatgen`, `matplotlib`, `numpy` for plumbing and analysis.

A working `conda` environment recipe is included in
[`methods.md`](methods.md).

## Two competing $T_m$ methods — pick coexistence as primary

You will use the **two-phase coexistence method** (sometimes called
"interface velocity" or "voronoi" method) as your primary $T_m$
estimator. The hysteresis method (heating an ordered solid until it
melts, cooling a liquid until it crystallises) is used as a
complementary, lower-quality bound.

Why coexistence is primary: the hysteresis method routinely
overestimates $T_m$ by 200–400 K because superheating of a defect-free
solid is enormous. Coexistence places a solid–liquid interface in the
cell at the start, removing the nucleation barrier and giving the true
thermodynamic $T_m$ once the interface stops moving.

Why hysteresis is still useful: it requires no special cell setup and
provides a quick reality check. If the hysteresis bracket and the
coexistence estimate disagree by less than 200 K, your pipeline is
behaving sensibly.

## Pitfalls flagged up front

1. **Training-data temperature gap.** If you only train on 300 K data
   you will get good elastic constants and a useless liquid; if you
   only train on 1500 K data the solid will be soft. You must include
   all relevant thermodynamic states in the training set.
2. **DFT-MD timestep.** A 0.5–1.0 fs timestep is appropriate for Cu;
   2 fs causes detectable energy drift in the NVE check.
3. **Strain sampling.** Without strained cells, your MLIP will get the
   ambient elastic constants wrong, which directly affects the solid's
   thermal expansion and therefore the coexistence pressure.
4. **Cell size at coexistence.** A 6000-atom slab cell is the
   sensible minimum; smaller cells produce interface-velocity noise
   that dwarfs the signal.
5. **Pressure choice.** Run coexistence at 1 bar (NPT with $p = 1$
   bar). Anyone who reports $T_m$ at zero pressure for a real liquid
   metal is reporting a different number.
6. **Force MAE in the liquid.** Liquid Cu has typical force magnitudes
   of ≈ 1 eV/Å; a 50 meV/Å MAE is roughly 5 % relative error, which is
   the threshold at which dynamical properties become quantitatively
   correct.

## Deliverables checklist

- [ ] `data/train.xyz`, `data/val.xyz`, `data/test.xyz` — extended-XYZ
      files with DFT energies, forces, and (where available) stresses.
- [ ] `mace/model.pt` — the trained MACE checkpoint.
- [ ] `parity/parity_energy.png`, `parity/parity_force.png` — held-out
      parity plots.
- [ ] `md/coexistence/T*.lammps` — coexistence-cell LAMMPS inputs.
- [ ] `analysis/interface_velocity.png` — interface velocity vs $T$;
      $T_m$ is the zero-crossing.
- [ ] `analysis/hysteresis.png` — heating and cooling curves with the
      bracket.
- [ ] `report.pdf`.
