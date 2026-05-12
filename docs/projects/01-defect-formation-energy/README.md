# Project 1 — Defect Formation Energy in Silicon

## Research question

What is the formation energy of a neutral monovacancy in crystalline silicon,
$E_f^{V_\mathrm{Si}^{0}}$, and how does the computed value converge as the
periodic supercell is enlarged? The point defect we target is among the
oldest problems in computational solid-state physics, yet it remains a
non-trivial test of every part of the modern DFT workflow: pseudopotentials,
plane-wave cut-offs, **k**-point sampling, ionic relaxation, and — most
importantly — the finite-size scaling needed to extrapolate to the dilute
limit. This project is the canonical "first real defect calculation" that
every electronic-structure practitioner should perform at least once.

## Why this project

A handful of reasons motivate spending six weeks on a single number:

1. The neutral Si vacancy has been studied for sixty years; comparison
   with experiment (positron annihilation, EPR on the negative state) and
   with a long literature of calculations gives you an unusually rich set
   of benchmarks.
2. The formation energy depends sensitively on supercell size because the
   strain field of the vacancy decays slowly. You will see, with your own
   numbers, why a single 64-atom cell is *not* converged and why
   extrapolation is required.
3. The neutral vacancy avoids — for now — the additional headache of
   image-charge corrections, but once your workflow is in place you can
   readily extend it to $V^{-1}_\mathrm{Si}$ and $V^{+1}_\mathrm{Si}$ and
   compare with the Freysoldt–Neugebauer–Van de Walle (FNV) scheme.
4. Almost every method you touch — `pw.x` from Quantum ESPRESSO, the
   Atomic Simulation Environment (ASE), the SSSP pseudopotential library,
   and a basic finite-size extrapolation — reappears in later chapters.

## Expected outcomes

By the end of the project you will deliver:

- A converged value of $E_f^{V_\mathrm{Si}^{0}}$ extrapolated to the
  infinite-supercell limit, with an uncertainty estimate of at most
  0.05 eV. The accepted PBE value is in the range 3.5–3.7 eV; LDA gives
  3.2–3.4 eV. Aim for agreement with the literature to within your
  reported uncertainty.
- A finite-size scaling plot of $E_f$ versus inverse supercell volume
  $N^{-1}$ (and, separately, versus $L^{-3}$ where $L$ is the cubic
  lattice parameter of the supercell), showing the linear behaviour
  expected from elastic-image arguments.
- A convergence study for plane-wave cut-off and **k**-point mesh on the
  smallest cell, used to fix parameters for production runs.
- A symmetry-broken relaxed geometry of the neutral vacancy showing the
  Jahn–Teller distortion (the four neighbouring atoms break $T_d$
  symmetry to $D_{2d}$ — verify this).
- A short (≈ 4 page) write-up presenting the workflow, the convergence
  plot, the relaxed geometry, and a comparison with at least two
  literature values.

## Time estimate

Six weeks, part-time (roughly 15 hours per week). The expected
distribution is:

| Week | Activity |
| --- | --- |
| 1 | Read background. Install QE locally. Run bulk Si SCF and verify lattice constant. Pick pseudopotential. |
| 2 | Plane-wave cut-off and **k**-mesh convergence on 8-atom cubic cell. |
| 3 | Build 64-atom supercell, run perfect and vacancy SCF, first crude $E_f$. |
| 4 | Geometry relaxation of the vacancy; check Jahn–Teller distortion. |
| 5 | Repeat for 216-atom and (if feasible) 512-atom supercells. |
| 6 | Finite-size extrapolation, write-up, comparison with literature. |

Be realistic: the 512-atom cell is the bottleneck. If your cluster only
gives you single-node throughput, plan to substitute a 343-atom cell or
to drop the 512-atom point altogether. The methodological lesson — that
you must scale — is preserved either way.

## Compute budget

Approximately 500 CPU-hours total, broken down as follows:

| Calculation | Cell | Cost (CPU-h) |
| --- | --- | --- |
| Convergence (cut-off + k-mesh) | 8-atom Si | 5 |
| 64-atom perfect Si SCF | 64 | 5 |
| 64-atom vacancy relaxation | 63 | 30 |
| 216-atom perfect Si SCF | 216 | 40 |
| 216-atom vacancy relaxation | 215 | 200 |
| 512-atom perfect + vacancy (optional) | 512/511 | 200 |

This is achievable on a small university partition (e.g. 4 nodes of 32
cores each for ≈ 4 hours per relaxation). If you only have access to a
laptop, the 64-atom point is doable in a weekend, and the 216-atom
point is realistic over a week, but the 512-atom calculation is not
practical without a cluster.

## Prerequisites

You must have worked through the following book chapters before starting:

- [Chapter 5 — Density functional theory](../../ch05-dft/index.md)
  for the meaning of exchange-correlation functionals and Kohn–Sham
  theory.
- [Chapter 6 — Running DFT in practice](../../ch06-running-dft/index.md)
  for the structure of a `pw.x` input file, pseudopotentials, plane-wave
  cut-offs, k-point convergence, and defect calculations.

You should be comfortable with bash, with submitting jobs to a SLURM
queue, and with at least one of NumPy/Matplotlib. A reading-level
familiarity with point-group symmetry ($T_d$, $D_{2d}$) helps when
interpreting the relaxed geometry but is not required.

## Software you will install

- Quantum ESPRESSO 7.2 or later, compiled with MPI.
- The Atomic Simulation Environment (`ase>=3.22`).
- The SSSP "efficiency" pseudopotential library (PBE) — download the
  tarball from [materialscloud.org/discover/sssp](https://www.materialscloud.org/discover/sssp).
- `pymatgen` for symmetry analysis of the relaxed structure.

Installation instructions and a verification script are given in
[`methods.md`](methods.md).

## Pitfalls flagged up front

1. **Pseudopotential transferability.** Do not mix pseudopotentials
   between the bulk and defect calculations. The same `Si.UPF` must be
   used for both, otherwise the energy difference is meaningless.
2. **Charge neutrality at fixed cell.** For the neutral defect, the cell
   stays charge-neutral by construction. The relaxation is performed at
   fixed volume — never relax the lattice parameter of a defect cell
   against vacuum, or you will accidentally see a "negative formation
   volume" artefact.
3. **k-mesh consistency.** The bulk and defect cells must sample the
   Brillouin zone equivalently. A $4 \times 4 \times 4$ mesh on the
   conventional 8-atom cell maps to a $2 \times 2 \times 2$ mesh on the
   $2 \times 2 \times 2 = 64$-atom supercell and to a $1 \times 1 \times 1$
   (Γ-only) mesh on the $3 \times 3 \times 3 = 216$-atom cell.
4. **Symmetry locking.** QE's default symmetry-detection will refuse to
   break the $T_d$ point group of the perfect crystal. You must set
   `nosym = .true.` or perturb the initial atomic positions by a small
   random displacement to allow the Jahn–Teller distortion to emerge.
5. **Smearing.** Si has a gap in the bulk but the gapped vacancy
   spectrum may close at finite size. Use a small Gaussian smearing
   (`degauss = 0.005 Ry`) for both calculations — never compare a smeared
   calculation against an unsmeared one.

## Deliverables checklist

- [ ] `convergence/` — cut-off and k-mesh convergence data, plots, README.
- [ ] `bulk/` — converged bulk Si SCF input/output and the bulk energy
      per atom $\mu_\mathrm{Si}$.
- [ ] `defect-N64/`, `defect-N216/`, optionally `defect-N512/` — each
      with `perfect.in`, `perfect.out`, `vacancy.in`, `vacancy.out`,
      and a `summary.json` that records the supercell size, the SCF and
      relaxation energies, and the computed $E_f$.
- [ ] `analysis/scaling.png` — $E_f$ vs $N^{-1}$ with a linear fit.
- [ ] `report.pdf` — at most six pages including figures.

When you are done, your supervisor should be able to reproduce your
plot in under an hour by re-running `analysis/scaling.py` on your
`summary.json` files.
