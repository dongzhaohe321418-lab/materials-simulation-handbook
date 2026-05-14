# Chapter 12 — Foundation Models for Materials

[Open in Colab](https://colab.research.google.com/github/dongzhaohe321418-lab/materials-simulation-handbook/blob/main/notebooks/ch12-foundation-colab.ipynb){ .md-button .md-button--primary }

Chapter 11 closed with a workflow that, by 2024 standards, was already
the state of the art: a Bayesian-optimisation loop driving a sequence of
DFT calculations, guided by a graph neural network trained on a few
thousand structures. The bottleneck, as ever, was the surrogate. Each
new chemical system demanded a new training set, a new model, a new
round of hyperparameter sweeps. The pipeline worked, but it did not
*scale* across the periodic table.

This chapter is about the architectural and conceptual shift that broke
that ceiling: the arrival of *foundation models* in materials science.
The term is borrowed deliberately from the language-model literature,
where a single network — pre-trained once, at great expense, on a vast
and heterogeneous corpus — turns out to be useful, with little or no
further training, for an astonishing range of downstream tasks. GPT,
CLIP and their descendants displaced an entire ecosystem of task-specific
NLP models. The question this chapter asks is whether the same thing is
happening, right now, in computational materials science.

The answer, as of 2026, is a qualified *yes*. Three lines of work make
the case.

First, **universal machine-learning interatomic potentials**. A single
MACE-MP-0, CHGNet or SevenNet model, trained on relaxation trajectories
spanning most of the periodic table, can run molecular dynamics on
chemistries it has never seen during training — and frequently match
the accuracy of bespoke potentials fitted to that exact system. The
"zero-shot universal MLIP" is now a routine starting point for any
new simulation, in much the way that one starts a chemistry calculation
by reaching for PBE rather than deriving a new functional.

Second, **generative models for crystal structures**. MatterGen, CDVAE,
DiffCSP and related diffusion-based architectures invert the standard
forward pipeline: instead of predicting a property from a structure,
they generate candidate structures conditioned on a target property.
The classical inverse-design problem, which dominated a generation of
PhD theses, is being repackaged as a sampling problem in a learned
manifold of stable crystals.

Third, and most speculatively, **multimodal models** that tie together
structure, computed spectra, experimental measurements and the text of
the scientific literature. The picture here is less mature, but the
direction is clear: a single embedding space in which a query like
"a transparent conductor with band gap above 3 eV, made of
earth-abundant elements" returns a ranked list of candidate compositions
together with predicted XRD patterns and citations to the closest
prior work.

The chapter is organised as follows. Section 12.1 lays out the
foundation-model paradigm in general terms — what it borrows from NLP
and computer vision, why materials science has the structural
properties (a small alphabet, a shared underlying physics) that make
the paradigm viable, and how to think about zero-shot, few-shot and
fine-tuning regimes in this context. We will look at the practical
question that every reader will eventually face: given a small,
domain-specific dataset, when does one fine-tune a universal model and
when does one train from scratch?

Section 12.2 turns to MACE-MP-0 in detail. We will install the model,
run a hundred-step molecular dynamics trajectory on a binary oxide with
zero further training, and inspect its energy and force errors against
single-point DFT. Then we will fine-tune it on a hundred structures of
a domain-specific system — a perovskite catalyst, say — and observe
the cost–accuracy curve. The section closes with a comparison table
of the universal-MLIP zoo as it stood at the end of 2025: MACE-MP-0,
CHGNet, M3GNet, SevenNet, Orb, and several smaller specialised
descendants.

Section 12.3 covers generative models. We will sketch the architecture
of MatterGen — a graph-based diffusion model with equivariance built in
— and walk through a complete validation pipeline: generate candidate
structures from the model, relax each with a universal MLIP, screen the
relaxed structures by DFT, and finally hand-check the most promising
against the synthesizability heuristics of Chapter 11. Limitations are
discussed candidly. The model knows nothing about kinetic accessibility
and cheerfully proposes phases that are thermodynamically reasonable
but practically unmakeable.

Section 12.4 looks forward. We discuss autonomous laboratories (A-Lab
at LBL, the Self-Driving Lab consortium), the closing of the
simulation–experiment loop, and the open problems that the current
generation of foundation models does not solve: long-range electrostatics,
charge transfer, excited states, magnetic ordering in MLIPs, and the
question of whether the universal models truly generalise across the
periodic table or merely interpolate within the training distribution.
A short reading list of papers from 2024–2026 points to the work most
worth following.

By the end of the chapter you will have run, fine-tuned and
critically assessed a universal MLIP, generated candidate structures
with a diffusion model, and understood enough of the landscape to
read the literature as it appears. The appendix that follows the
chapter consolidates the mathematical, computational and
bibliographic resources used throughout the book.

A final remark. Foundation models do not abolish the methods of
Chapters 5 through 9. DFT remains the reference. Classical MD remains
the workhorse for million-atom simulations. What changes is the entry
point: most new projects, in 2026, begin with a pre-trained model and
add specificity only where the data demand it. The art lies in knowing
when to trust the foundation and when to bypass it.
