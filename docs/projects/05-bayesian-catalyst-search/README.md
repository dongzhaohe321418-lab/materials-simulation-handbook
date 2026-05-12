# Project 5 — Bayesian optimisation for catalyst composition

## Research question

Within the ternary alloy composition space Pt–Pd–Ag, can you identify
compositions whose computed hydrogen-adsorption Gibbs free energy
$\Delta G_\mathrm{H}$ is closest to zero (the Nørskov "ideal HER
catalyst" descriptor), using ≤ 50 oracle (DFT) calls? Compare a
Bayesian-optimisation (BO) strategy against a random-search baseline
of the same budget.

This project teaches you the closed loop that drives modern
"autonomous experimentation" for catalyst design: a surrogate model
(Gaussian process), an acquisition function (qEI, UCB, ...), an
oracle (DFT or an MLIP), and an iteration policy.

## Why this project

Hydrogen evolution catalysis is the cleanest test case for the
descriptor approach. Nørskov's "volcano" plot for HER plots
exchange-current density against $\Delta G_\mathrm{H}$, and the
optimum is at $\Delta G_\mathrm{H} \approx 0$. Pt sits near the top
of the volcano; many noble-metal alloys cluster nearby. The
*combinatorial* question — which Pt-rich alloy beats pure Pt? — is
the discovery question, and it is small enough to fit BO comfortably.

What makes the project pedagogically rich:

- The objective is scalar and reasonably smooth in composition space,
  but not trivial — different sites on the slab give different
  $\Delta G_\mathrm{H}$, and the choice of representative site
  matters.
- Compositions are *discrete* (a slab has a finite number of atoms),
  which forces you to deal with categorical / mixed-integer BO.
- The oracle is genuinely expensive (a DFT slab calculation),
  motivating the BO framing — random search of 50 candidates would
  cost more than 50 GPU-hours, and you want to do better.

## Expected outcomes

By the end of the project you will deliver:

- A posterior heatmap (or ternary-diagram visualisation) of predicted
  $\Delta G_\mathrm{H}$ over the Pt–Pd–Ag composition space at the
  end of the BO loop.
- A *regret curve*: the best $|\Delta G_\mathrm{H}|$ found so far vs
  iteration number, for both BO and random-search baselines, averaged
  over multiple BO seeds.
- A Pareto-style ranking of the top 5 compositions with their DFT
  $\Delta G_\mathrm{H}$ values and ensemble uncertainty.
- A short discussion of (a) whether BO beat random search, (b) the
  failure modes you encountered, and (c) the role of the *oracle*
  choice (DFT vs MLIP surrogate).
- A 4–6 page report.

## Time estimate

Six weeks, part-time.

| Week | Activity |
| --- | --- |
| 1 | Background reading. Install BoTorch. Set up the oracle (DFT slab template). |
| 2 | Build the featurisation (composition vector → Magpie or simple Pt/Pd/Ag fractions). Validate the oracle on pure Pt, Pd, Ag (literature comparison). |
| 3 | Implement the BO loop with a GP surrogate, Matérn 5/2 kernel, qEI acquisition. Test on a synthetic 2D function first. |
| 4 | Run BO on the real oracle (DFT or MLIP). 50 iterations. |
| 5 | Run random-search baselines (3 seeds × 50 iterations). Compare. |
| 6 | Write-up. |

## Compute budget

Two oracle choices, each with different cost:

**Oracle A: DFT slab calculation.** Per-call cost: ≈ 4–10 CPU-hours
on 16 cores. Total: 50 calls × 8 CPU-hours = 400 CPU-hours, plus
random-search baselines.

**Oracle B: MLIP-based surrogate.** Run a pre-trained MLIP (e.g.,
MACE-MP-0) on the same slab geometry. Per-call cost: seconds. Total:
< 1 GPU-hour. This is faster but introduces MLIP error as an extra
layer; you must validate the MLIP on a subset of compositions with
DFT.

For an undergraduate project, oracle B with DFT validation of the
final top 5 is the recommended balance. If you have generous cluster
access, oracle A is preferred for purity.

## Prerequisites

- [Chapter 6 — Running DFT in practice](../../ch06-running-dft/index.md)
  — slab geometry, dipole correction, **k**-mesh on a slab.
- [Chapter 8 — Statistical mechanics for simulation](../../ch08-statmech/index.md)
  — the meaning of a Gibbs free energy of adsorption and the
  computational hydrogen electrode reference.
- [Chapter 11 — Active learning and Bayesian optimisation for
  materials](../../ch11-active/index.md) — GPs, kernels, acquisition
  functions, the BO loop, and the regret framework.

You should be comfortable with `torch`, `gpytorch`, and `botorch`, or
prepared to learn them. You should already understand what a Gaussian
process *is*; this project will not teach the theory from scratch.

## What is $\Delta G_\mathrm{H}$ exactly

For HER, the Nørskov descriptor is

$$
\Delta G_\mathrm{H} = E_\mathrm{slab+H} - E_\mathrm{slab} - \tfrac{1}{2} E_\mathrm{H_2} + \Delta E_\mathrm{ZPE} + T \Delta S,
$$

where the ZPE + entropy correction is typically taken as
$+ 0.24$ eV at $T = 298$ K. Your oracle therefore computes:

1. The slab energy without H ($E_\mathrm{slab}$).
2. The slab energy with one H atom adsorbed at the most-stable hollow
   or top site ($E_\mathrm{slab+H}$).
3. The H$_2$ molecule energy in vacuum ($E_\mathrm{H_2}$), computed
   once.

You then combine with the standard ZPE + entropy correction. The
output of one oracle call is a single $\Delta G_\mathrm{H}$ value.

## Pitfalls flagged up front

1. **Random-alloy slabs without averaging.** A Pt$_x$Pd$_y$Ag$_z$
   slab with a *specific* atomic arrangement is not the same as the
   alloy "Pt$_x$Pd$_y$Ag$_z$". You must either (a) average over
   multiple random arrangements, (b) use the special quasi-random
   structure (SQS) approach, or (c) use a fixed *coverage-only*
   parametrisation. Pick one and stick with it.
2. **Site choice.** The hydrogen adsorption energy depends on the
   site (top, bridge, hollow, sub-surface). The *most-favourable*
   site for a given composition may not be the site you chose by hand.
   For pure Pt(111), the fcc hollow is the standard choice;
   document the convention and apply it consistently.
3. **Discrete composition space.** A slab of 32 surface atoms can
   have compositions 0/32, 1/32, ..., 32/32 of each species. The BO
   variables are *integer-valued* counts subject to a sum constraint.
   BoTorch supports mixed-integer search via SAASBO or via rounding;
   the simpler route is to relax to continuous compositions, then
   round to the nearest discrete composition before each oracle call.
4. **Oracle reliability.** A single failed DFT calculation in the BO
   loop poisons the GP posterior. Always check that each oracle call
   returned a sensible number before incorporating it.
5. **Comparing BO vs random search.** A single BO seed vs a single
   random seed is not a comparison. Run ≥ 3 seeds for each and report
   means with error bars.
6. **Kernel choice.** A Matérn 5/2 kernel on the composition simplex
   is sensible; an RBF kernel often gives misleadingly smooth
   posteriors. Justify your choice.
7. **Confusing $\Delta G$ with $\Delta E$.** Reporting an adsorption
   *energy* without the ZPE + entropy correction is not the same
   quantity as $\Delta G_\mathrm{H}$ on the Nørskov volcano. State
   the convention explicitly.

## Deliverables checklist

- [ ] `oracle/` — the DFT or MLIP oracle wrapper, with a `single_call`
      function and a test on Pt(111) reproducing $\Delta G_\mathrm{H}
      \approx -0.1$ eV.
- [ ] `bo/loop.py` — the BO loop.
- [ ] `bo/results/seed{0..2}.json` — the BO histories for three
      different random seeds.
- [ ] `random/results/seed{0..2}.json` — the random-search baselines.
- [ ] `analysis/regret.png` — the regret curve with both methods and
      error bars.
- [ ] `analysis/posterior.png` — the final GP posterior heatmap.
- [ ] `analysis/top5.csv` — the top 5 compositions with their oracle
      values.
- [ ] `report.pdf`.
