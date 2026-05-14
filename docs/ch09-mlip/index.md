# Chapter 9 — Machine Learning Interatomic Potentials

[Open in Colab](https://colab.research.google.com/github/dongzhaohe321418-lab/materials-simulation-handbook/blob/main/notebooks/ch09-mlip-colab.ipynb){ .md-button .md-button--primary }

For seventy years the practitioner of atomistic simulation has faced an
uncomfortable choice. Quantum-mechanical methods such as density functional
theory (Chapter 5) deliver chemically accurate energies and forces, but the
$O(N^3)$ scaling and large prefactor limit them to systems of a few hundred
atoms and timescales of picoseconds. Classical force fields (Chapter 7),
written as hand-tuned functional forms in interatomic distances and angles,
scale linearly and run on a laptop, but they capture only a narrow slice of
chemistry: change the bonding environment and the parameters break.

The gap between these two regimes — call it the *accuracy–cost gap* — has
shaped what materials scientists could and could not study. Reactive
chemistry in a million-atom electrolyte, the nucleation of a solid from its
melt, the diffusion of a defect over microseconds: all sat in the no-man's
land between DFT and Lennard-Jones.

Machine learning interatomic potentials (MLIPs) close that gap. The idea is
brutally simple: generate a few thousand DFT calculations on small reference
configurations, fit a flexible regression model to the resulting
energy-force-stress data, and use the fitted model in place of DFT inside
your molecular dynamics or geometry-optimisation engine. If the model
generalises, you obtain DFT-quality forces at the cost of evaluating a few
million floating-point operations per atom — roughly four orders of
magnitude cheaper than the underlying DFT.

What makes this work is not the universal-approximator theorem (regression
on a sufficiently rich basis can fit anything; the question is whether it
generalises) but a careful encoding of physics. An interatomic potential is
not an arbitrary function of $3N$ coordinates: it is invariant under
translation, rotation, and permutation of identical atoms; it must be
smooth so that forces are continuous; it must be local so that energy
decomposes into atomic contributions that depend only on a neighbourhood.
Every successful MLIP architecture — Behler–Parrinello neural networks,
Gaussian Approximation Potentials, SchNet, NequIP, MACE — bakes these
constraints into its representation. Models that ignore them require
orders of magnitude more training data, and still extrapolate badly.

This chapter develops the field from first principles. We begin with the
motivation: why MLIPs exist and what they can do (Section 9.1). We then
catalogue the symmetries any potential must respect, and distinguish
*invariant* features (scalars under rotation) from *equivariant* features
(vectors and tensors that transform predictably) — a distinction that
separates first-generation MLIPs from the modern equivariant networks
(Section 9.2). With those constraints in hand we survey the three
descriptor families that dominate the literature: Behler–Parrinello
symmetry functions, the Smooth Overlap of Atomic Positions (SOAP), and the
Atomic Cluster Expansion (ACE) (Section 9.3).

Section 9.4 turns to two complete architectures built on these descriptors:
Behler–Parrinello neural networks and Gaussian Approximation Potentials
(GAP). These are the workhorses of the 2010s and remain competitive on many
problems. Section 9.5 introduces the equivariant revolution — NequIP and
MACE — and explains why operating on irreducible representations of
$\mathrm{O}(3)$ yields models that learn roughly twenty times faster than
their invariant predecessors. Section 9.6 is hands-on: we install
`mace-torch`, train a MACE potential on a water dataset, and wire the
trained model into ASE for molecular dynamics. The chapter closes with
eight exercises that take the reader from deriving forces by autograd to
fitting and validating their own potential.

**Prerequisites.** You should be comfortable with the DFT workflow of
Chapter 5 (we use it to generate labels), with the MD integrators of
Chapter 7 (we will feed the MLIP into them), and with the statistical
mechanics of Chapter 8 (we use ensembles and free energies as validation
tools). On the machine-learning side we assume only that you have
encountered a feed-forward neural network: every architecture in this
chapter is built up explicitly from linear layers, nonlinearities, and a
small amount of group theory which we develop as needed.

**What to expect.** By the end of the chapter you will have:

- a precise statement of the symmetries any interatomic potential must
  satisfy, and the tools to check them;
- a working implementation of a Behler–Parrinello $G^2$ descriptor;
- the ability to read and critique papers in the SOAP/ACE/MACE family;
- a trained MACE potential for liquid water with a force MAE of roughly
  $30\,\mathrm{meV}/\text{\AA}$;
- and the connections to Chapter 10 (graph neural networks generalise the
  message-passing structure introduced here), Chapter 11 (active learning
  closes the loop by querying DFT only where the MLIP is uncertain), and
  Chapter 12 (foundation models such as MACE-MP-0 amortise the training
  cost across all of materials chemistry).

The chapter is opinionated. We have chosen MACE as the working example
because, at the time of writing (2026), it represents the best trade-off
between accuracy, data-efficiency, and inference cost for a single
researcher with one GPU. The principles transfer directly to NequIP,
Allegro, SevenNet, and the other equivariant architectures you will
encounter in the literature.
