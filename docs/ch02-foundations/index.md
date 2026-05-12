# Chapter 2 — What is Materials Simulation?

Materials simulation is the art of replacing an experiment with a computation. A century ago, if you wanted to know the melting point of a new alloy, you melted it. Today, in many cases, you can predict it from the periodic table, a cluster of cores, and a few hours of wall time. The same shift has happened, with varying degrees of success, for elastic constants, band gaps, catalytic activities, phonon spectra, defect formation energies, and a dozen other quantities that used to belong exclusively to the laboratory. This chapter asks what kind of activity that is, what it can and cannot deliver, and what the working tools of the trade look like in 2026.

## Why a separate chapter on foundations?

Most textbooks plunge directly into a particular method — molecular dynamics, density functional theory, phase-field modelling — and the reader emerges, several hundred pages later, with a working knowledge of that method and almost no sense of the wider landscape. This is fine if you already know which method you need. It is fatal if you do not. A graduate student arriving in a computational group typically faces a triage problem: a supervisor mentions a material, a phenomenon, and a deadline, and the student must decide whether the right tool is DFT, classical molecular dynamics, a kinetic Monte Carlo simulation, or a graph neural network. Picking the wrong one wastes months.

The aim of this chapter is therefore orientation. We map the scale ladder from electrons to engineering components, we describe what each rung produces and at what cost, and we offer an honest appraisal of the field's current capabilities. We finish with a survey of the software ecosystem you are likely to encounter, and a recommended starter stack for someone who wants to run their first calculation by the end of the week.

## What you will know by the end

After working through Chapter 2 you should be able to:

- name the principal length and time scales of materials simulation and the dominant method at each;
- match a research question (*what is the band gap of MoS$_2$?*; *how does this polymer melt?*) to an appropriate computational approach;
- describe the trade-off between accuracy and system size that drives almost every methodological choice;
- read the three diagrams that pervade the materials literature — phase diagrams, band structures, and radial distribution functions — and explain what each tells you;
- recognise the major open-source and commercial codes by name, and place them on the scale ladder;
- articulate a realistic view of what materials simulation can deliver in 2026, including the role of machine-learning interatomic potentials and foundation models.

## Roadmap

The chapter contains five thematic sections and an exercise set:

1. **[The scale ladder](01-scales.md)** — From sub-ångström electronic structure to metre-scale finite-element analysis. We introduce DFT, molecular dynamics, mesoscale methods, and continuum mechanics on a single diagram and discuss multiscale coupling.

2. **[What simulation can and cannot do](02-capabilities.md)** — A 2026 audit. We list the predictions you can trust, the predictions you should not trust, and the open frontiers where the answer is *we are working on it*.

3. **[Three diagrams every materials scientist reads](03-three-diagrams.md)** — Phase diagrams, band structures, radial distribution functions. The visual vocabulary of the field, with worked examples and Python snippets.

4. **[The software ecosystem](04-ecosystem.md)** — DFT codes, MD codes, workflow managers, ML stacks, databases, visualisation tools. Names, niches, and a recommended starter stack.

5. **[Exercises](exercises.md)** — Six mostly conceptual problems with full worked solutions.

!!! note "Prerequisites"
    We assume you have worked through Chapters 0 and 1 — that is, you are comfortable with linear algebra at the level of eigenvalue problems, you have seen Fourier transforms once, and you can write small NumPy scripts. No prior exposure to materials science is required; we introduce the relevant concepts as they arise.

!!! tip "Reading strategy"
    If you are short of time, read `01-scales.md` and `02-capabilities.md` carefully and skim the rest. These two sections frame every subsequent chapter. The other three sections are reference material you will return to as your projects demand.

## A word on perspective

It is easy, when reading any textbook in this field, to come away with the impression that materials simulation is a solved problem — that you point a code at a structure, wait, and receive the truth. This is not the case. Every method described in this handbook makes approximations, and every approximation fails somewhere. The most valuable skill a computational materials scientist can develop is not the ability to run a code, but the ability to ask, of any result, *do I believe this?* Chapter 2 is the first step in cultivating that scepticism. The capabilities section in particular is written in a deliberately critical voice. Take it seriously.

With the landscape in view, we will spend the rest of the book picking individual peaks and climbing them properly. By Chapter 12 you will have a foundation-model perspective on the whole range. For now, we walk to the base camp.
