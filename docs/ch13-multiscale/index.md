# Chapter 13 — Multiscale and Coupled Methods

!!! abstract "Where this chapter sits"
    You have now met density functional theory ([Chapter 5](../ch05-dft/index.md)),
    molecular dynamics ([Chapter 7](../ch07-md/index.md)), machine-learning
    interatomic potentials ([Chapter 9](../ch09-mlip/index.md)), graph neural
    networks ([Chapter 10](../ch10-gnn/index.md)), active learning
    ([Chapter 11](../ch11-active/index.md)) and foundation models
    ([Chapter 12](../ch12-foundation/index.md)). Each method has a window of
    competence. This chapter is about what to do when your problem refuses to
    sit inside a single window.

## The single-method delusion

A common temptation, especially when one has just learnt a powerful technique,
is to believe that the technique is enough. DFT can do anything if you have
enough cores. MD can do anything if you have a good potential. An MLIP can do
anything if you have enough training data. None of this is true.

Real materials problems have a habit of spanning many decades of length and
time. Consider a few concrete examples:

- **Stress-corrosion cracking of a steel turbine blade.** The chemistry of
  hydrogen embrittlement happens at the electronic-structure scale (an angstrom,
  a femtosecond). The crack tip is a defect of nanometres. The blade is tens of
  centimetres. The engineering question — how long until it fails — lives at the
  scale of years.
- **Solid-electrolyte interphase formation in a battery.** Bond-breaking is
  quantum; ion transport in the bulk is classical MD or beyond; the SEI grows
  over thousands of cycles, which is hours of laboratory time.
- **Sintering of a structural ceramic.** Atomic-scale grain boundary motion
  feeds grain-growth kinetics, which feed microstructure evolution, which feed
  the macroscopic strength reported in a datasheet.

No single computational method spans even two of these scales, never mind four.
The job of the multiscale practitioner is to choose which methods sit at which
scales, and to make sure that information flows correctly between them.

## What "multiscale" actually means

The word is overused. In this handbook we will be careful to distinguish two
quite different ideas.

**Multiscale** modelling uses *different physics* at different scales. The
quantum mechanical region near a defect uses DFT; the elastic far-field uses a
classical force field or a continuum model. The methods are not interchangeable.
Each is doing a job that the other cannot do.

**Multifidelity** modelling uses *different cost versions of the same physics*.
A coarse mesh and a fine mesh. A small basis set and a large basis set. A
smaller MLIP and a larger one. The goal is statistical efficiency: use the
cheap model often, the expensive one rarely, and combine them carefully.

The two ideas are sometimes blended, especially in machine learning workflows,
and we will see this in [Section 1](01-coupling-schemes.md). But the
distinction matters when you choose your method.

## A roadmap for the chapter

[Section 1](01-coupling-schemes.md) sets out the two fundamental ways to couple
methods: *sequentially*, where one feeds into the next, and *concurrently*,
where two run side by side on different parts of the same system. We discuss
the hand-shake region — the place where errors live.

[Section 2](02-qm-mm.md) is the canonical concurrent multiscale technique:
QM/MM. We derive both the subtractive (ONIOM) and additive formulations, work
through how electrostatic embedding actually goes through the equations, and
write a working ASE example you can run on a laptop.

[Section 3](03-coarse-graining.md) goes the other direction: throwing away
degrees of freedom to reach longer time and length scales. We cover the MARTINI
force field for biomolecules, bottom-up CG (iterative Boltzmann inversion,
force matching, relative entropy), and the recent wave of machine-learned CG
models.

[Section 4](04-kmc-phase-field.md) covers two methods that bridge atomistic
and continuum without explicitly representing every atom: kinetic Monte Carlo
for rare-event dynamics on a lattice, and phase-field for microstructure
evolution. Both consume atomistic inputs (barriers, free energies) and produce
predictions at scales that MD will never reach.

[Section 5](05-finite-element-bridge.md) connects atomistic simulation to the
engineering world via finite-element analysis. We walk through how DFT elastic
constants become FEM stiffness matrices, and we are clear about what FEM
fundamentally cannot see.

The [exercises](exercises.md) give you a chance to practise picking the right
method for a problem, and to spot what is wrong when methods are stitched
together badly.

## Two things to keep in your head

Before you start, internalise two habits of mind. They will save you a great
deal of grief.

First, **information flow has a direction**. When you couple methods, you must
be explicit about which method is the boss at each scale, and which is the
client. If your CG model and your atomistic model disagree, which one is right?
The honest answer is that they are both right within their domain of validity,
and the multiscale worker's job is to know where those domains end.

Second, **errors compound**. If your DFT energies are good to 1 meV/atom, your
fitted force field reproduces them to 5 meV/atom, your MD samples free energies
to 20 meV/atom, your phase-field free energy interpolation introduces another
50 meV/atom, and your FEM uses macroscopic averaging — the error budget for
the predicted strength may be 100% before you have done any wet experiment. A
good multiscale workflow tracks this budget. A bad one ignores it and is
surprised when experiments disagree.

We will return to both of these themes throughout the chapter.
