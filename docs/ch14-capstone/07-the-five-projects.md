# The Five Projects

This section is a guide to the five undergraduate-thesis-scale projects
provided in this handbook's repository. Each project lives in
`/projects/0X-name/` and has its own detailed `README.md` covering
background, step-by-step instructions, expected outcomes, and example
analysis scripts.

The five projects span the major techniques in this handbook. They are
designed to be:

- **Tractable**: each fits in 6-10 weeks of part-time work, or 3-5 weeks
  full-time.
- **Self-contained**: the project README has everything you need to
  start.
- **Open-ended at the top**: the basic version is achievable; extensions
  for ambitious students are sketched.
- **Properly bounded**: each has a clear "minimum viable result" so that
  even a struggling project produces a defensible thesis.

Pick *one* project. Do it well. Do not try two.

## Project 1: Defect Formation Energy in Silicon

**Project folder:** [`projects/01-defect-formation-energy/`](../projects/01-defect-formation-energy/README.md)
**Duration:** ~6 weeks
**Difficulty:** Foundational (you have all the prerequisites)
**Methods used:** DFT only

### What you will do

Compute the formation energy of vacancy and substitutional defects in
crystalline silicon using DFT. Carry out convergence studies for cutoff
energy, k-point grid, and supercell size. Compare your results to
published values. Extend (if time allows) to charged defects with
image-charge corrections.

### Why it is a good thesis project

Silicon is the most-studied semiconductor on Earth. Published values
exist for every defect you might compute, so you have many benchmarks.
The system is simple enough to converge cleanly on modest hardware.
You will learn the full DFT workflow — geometry, k-points, defects,
formation energies — without getting bogged down in pathological
materials physics.

### Expected outcomes

- Converged formation energies for the silicon monovacancy and one
  substitutional dopant.
- Convergence study plots for $E_\mathrm{cut}$, k-grid, and supercell.
- Comparison table with two or three published DFT values.
- (Stretch) Charged-defect formation energies with Freysoldt correction.

### Prerequisites

- [Chapter 5](../ch05-dft/index.md) (DFT)
- [Chapter 6](../ch06-running-dft/index.md) (running DFT in practice)
- [Section 4 of this chapter](04-convergence-and-validation.md)

### Computational requirements

A modest workstation (16 cores) or a small cluster allocation is
sufficient. Total compute: roughly 200-400 CPU-hours including
convergence studies.

### What can go wrong

Charged defects are surprisingly delicate; the image-charge correction
is a subtle topic. If you find yourself stuck, fall back to neutral
defects only — they are perfectly thesis-worthy.

---

## Project 2: Melting Point of Copper via MLIP-Driven MD

**Project folder:** [`projects/02-melting-point-mlip/`](../projects/02-melting-point-mlip/README.md)
**Duration:** ~8 weeks
**Difficulty:** Intermediate
**Methods used:** MD + pre-trained MLIP (or classical EAM as fallback)

### What you will do

Compute the melting point of copper by molecular dynamics, using either
a pre-trained MLIP (CHGNet, MACE-MP-0, or similar foundation model) or
a classical EAM potential as a fallback. Use the two-phase
(solid-liquid coexistence) method to bracket the melting temperature.
Compare to the experimental value (1357 K).

### Why it is a good thesis project

Melting points are notoriously hard to compute reliably. Doing this
project teaches you finite-temperature MD, the two-phase method
(superheating, supercooling, free-energy considerations), proper
thermostatting, and how to assess the quality of a force field at a
phase transition. The experimental value is known, so you have a
hard validation target.

### Expected outcomes

- A clear estimate of the melting temperature with error bars.
- Equilibration plots, energy histograms, and a phase-coexistence
  trajectory.
- Comparison of the MLIP result and the EAM result with experiment.
- Discussion of why each method succeeds or fails.

### Prerequisites

- [Chapter 7](../ch07-md/index.md) (MD)
- [Chapter 8](../ch08-statmech/index.md) (stat mech)
- [Chapter 12](../ch12-foundation/index.md) (foundation models, for the
  pre-trained MLIP path)

### Computational requirements

GPU or moderate CPU resources. With a pre-trained MLIP, total compute is
~100-300 GPU-hours. With EAM, much less; can run on a workstation.

### What can go wrong

The two-phase method requires careful setup — you need to construct a
solid-liquid interface and watch which side grows. Equilibration is
the dominant time cost. Be patient and run multiple long trajectories.

---

## Project 3: High-Throughput Band Gap Screening

**Project folder:** [`projects/03-band-gap-screening/`](../projects/03-band-gap-screening/README.md)
**Duration:** ~8 weeks
**Difficulty:** Intermediate
**Methods used:** Database query + GNN inference + DFT validation

### What you will do

Use the Materials Project database to identify candidate materials with
band gaps in a target range (e.g. 1.5-2.0 eV, useful for photovoltaics).
Use a pre-trained GNN (e.g. MEGNet, ALIGNN) to predict band gaps for a
much larger set of compositions. Validate the top candidates with
DFT calculations.

### Why it is a good thesis project

This is a complete data-science / high-throughput materials discovery
workflow in miniature: a database query, an ML prediction step, an
expensive validation step, and an analysis of how well the cheap ML
filter performed. You will learn database APIs, GNN inference, and how
to handle a high-throughput pipeline.

### Expected outcomes

- A list of 10-20 candidate materials with predicted gaps in your
  target range.
- DFT-validated gaps for 5 of these candidates.
- A precision-recall analysis of the GNN filter.
- Discussion of the failure modes (where the GNN predicts well, where
  it fails).

### Prerequisites

- [Chapter 10](../ch10-gnn/index.md) (GNNs)
- [Chapter 5](../ch05-dft/index.md) (DFT)
- Some Python and command-line comfort (Materials Project API,
  Pymatgen)

### Computational requirements

GNN inference is cheap (can run on a laptop). DFT validation is the
cost: 5 calculations × ~20 CPU-hours each, plus convergence studies.
Total ~200 CPU-hours.

### What can go wrong

PBE band gaps systematically underestimate experimental gaps by 30-50%.
Be clear in your thesis whether you are reporting PBE gaps, gaps with a
known scissor correction, or experimentally-comparable gaps from a
hybrid functional (which adds cost). Choose your story before launching.

---

## Project 4: Training an MLIP From Scratch on a Novel System

**Project folder:** [`projects/04-mlip-from-scratch/`](../projects/04-mlip-from-scratch/README.md)
**Duration:** ~10 weeks
**Difficulty:** Advanced
**Methods used:** DFT (data generation) + MLIP training + MD validation

### What you will do

Pick a system that does not yet have a publicly available accurate MLIP
(your supervisor will help you identify one — many surfaces, defects,
or alloy compositions qualify). Generate DFT training data covering the
configurations of interest. Train an MLIP (MACE, NequIP, or simpler
choices like SOAP-GAP) on the data. Validate against held-out DFT
configurations and against a target property (lattice parameter,
phonon spectrum, defect formation energy) where you have a reference.

### Why it is a good thesis project

Training an MLIP from scratch is, in a real sense, *the* defining skill
of modern computational materials science. Doing this project gives
you experience with data generation, hyper-parameter selection,
training-set curation, train/test splits, and out-of-distribution
validation. You will leave the project knowing how to build an MLIP for
*any* system you encounter later.

### Expected outcomes

- A trained MLIP with validation MAE quoted and contextualised.
- A training/validation curve showing convergence with dataset size.
- Validation of one physical property against DFT.
- Discussion of where the MLIP can and cannot be trusted.

### Prerequisites

- [Chapter 9](../ch09-mlip/index.md) (MLIPs)
- [Chapter 11](../ch11-active/index.md) (active learning, helpful for
  data selection)
- Comfort with PyTorch or JAX.

### Computational requirements

Substantial. DFT data generation is the dominant cost: maybe 500-2000
CPU-hours. Training: 50-200 GPU-hours. Best done with access to a
small GPU cluster.

### What can go wrong

The most common failure: the MLIP fits the training data well but
fails outside it. Mitigate with proper held-out test sets and, if
time permits, an active-learning loop. Be honest in your thesis about
the regime of validity.

---

## Project 5: Bayesian Optimisation for Catalyst Composition

**Project folder:** [`projects/05-bayesian-catalyst-search/`](../projects/05-bayesian-catalyst-search/README.md)
**Duration:** ~6 weeks
**Difficulty:** Intermediate
**Methods used:** Bayesian optimisation + GNN oracle (or DFT)

### What you will do

Set up a Bayesian-optimisation loop to identify the binary or ternary
alloy composition that maximises (or minimises) a catalytic property
— for example, hydrogen adsorption energy near an ideal value for the
hydrogen evolution reaction. Use a pre-trained GNN as the surrogate
"oracle" for the property, or DFT for a smaller search. Compare
random-search baselines.

### Why it is a good thesis project

Active learning and Bayesian optimisation are central to modern
materials discovery. This project gives you hands-on experience with
the BO loop, acquisition functions (expected improvement, upper
confidence bound), Gaussian processes, and exploration-exploitation
trade-offs. The scope is naturally bounded; the loop either converges
or it doesn't, and either outcome is writable.

### Expected outcomes

- A BO loop that identifies an optimal composition in fewer steps than
  random search.
- Visualisation of the acquisition function and the model uncertainty
  over the composition space.
- A comparison of two or three acquisition strategies.
- Discussion of the limits: where BO works, where it gets stuck.

### Prerequisites

- [Chapter 11](../ch11-active/index.md) (active learning and BO)
- Some Python comfort (you will use `scikit-learn`,
  `gpytorch`, or `BoTorch`).

### Computational requirements

Modest. The BO loop is light; the oracle dominates the cost. With a
pre-trained GNN, the entire project can run on a workstation. With
DFT, ~100 CPU-hours.

### What can go wrong

The BO loop can converge to a local optimum and miss the global
one — this is a feature of the algorithm, not a bug. Discuss it
honestly. The other common failure: the oracle is biased, so the
"optimum" found is an artefact of the surrogate model.

---

## A comparative table

| # | Project | Methods | Duration | Compute | Difficulty |
|---|---|---|---|---|---|
| 1 | Si defects | DFT | 6 weeks | Low | Foundational |
| 2 | Cu melting | MD + MLIP | 8 weeks | Moderate (GPU) | Intermediate |
| 3 | Band gap screening | GNN + DFT | 8 weeks | Moderate | Intermediate |
| 4 | MLIP from scratch | DFT + MLIP training | 10 weeks | High | Advanced |
| 5 | BO for catalyst | BO + GNN oracle | 6 weeks | Low | Intermediate |

## How to choose

A few practical considerations.

**If you have only used DFT before**: Project 1 is the natural
extension and gives you the cleanest learning curve.

**If you enjoyed MD or are interested in finite-temperature
properties**: Project 2.

**If you are excited about data science and want to learn the
high-throughput workflow**: Project 3.

**If you want the deepest technical project that prepares you for a
PhD in computational materials**: Project 4.

**If you are interested in algorithms and optimisation**: Project 5.

**If your supervisor has a strong preference**: take their preference,
unless it conflicts violently with your own. They know more than you do
about what will succeed in their group.

## How to use the project READMEs

Each project README is structured as follows:

1. **Background** — short scientific context.
2. **Question** — the specific question you will answer.
3. **Plan** — a week-by-week (or stage-by-stage) plan.
4. **Setup** — software to install, data to download.
5. **Tasks** — concrete tasks, each with expected output.
6. **Validation** — how to know you got the right answer.
7. **Extensions** — for if you have time spare.
8. **References** — starter reading.

Treat the README as a *starting* document, not a complete recipe. Your
thesis will deviate from it as you find unexpected results, hit
roadblocks, or discover that a particular step needs more attention
than the README anticipated.

That deviation is the project. The README gives you a map; you walk the
terrain.

## A closing word

You have now read fourteen chapters of methodology and finished with
five concrete project options. The remaining ingredient — the one that
no handbook can supply — is the *patient discipline* of seeing a
project through from start to finish.

The end of the road is a thesis, written by you, defending claims that
you have personally established by computation. That is research.

[Exercises](exercises.md) — final exercises for the chapter.
