# Background reading — Training an MLIP from scratch

Read in order. Each annotation explains how the paper informs the
project.

---

## 1. Behler (2017) — *First Principles Neural Network Potentials for Reactive Simulations of Large Molecular and Condensed Systems*

> Behler, J. *Angew. Chem. Int. Ed.* **56**, 12828 (2017).

A practitioner's review of how to build a high-quality neural network
potential. Predates MACE but the practical advice transfers verbatim.

**Why read this.** The discussion of training-set construction — what
configurations to include, how many, at what temperatures — is the
single most important piece of advice for this project. The MLIP
architecture is the cheap part; sampling is the expensive part.

**Extract.** The recommended training-set composition (equilibrium
geometries, thermal samples, transition-state-like geometries, and
strained geometries); the rule of thumb that ≈ 10 frames per
chemical environment is the minimum; the warning that "the model
will not extrapolate".

---

## 2. Batatia, Kovács, Simm, Ortner, and Csányi (2022) — *MACE: Higher Order Equivariant Message Passing Neural Networks*

> Batatia, I. et al. *NeurIPS* 2022.

The MACE architecture paper.

**Why read this.** You are using this exact model. As in Project 2,
you must know the meaning of `max_L`, the body order, and the
relationship between the receptive field and the cutoff radius.

**Extract.** The architectural details and recommended defaults; the
typical sensitivity to hyperparameters (forces are most sensitive to
`max_L` and `hidden_irreps`; energies are most sensitive to the
training-set size).

---

## 3. Vandermause, Torrisi, Batzner et al. (2020) — *On-the-fly active learning of interatomic potentials for large-scale atomistic simulations*

> Vandermause, J. et al. *npj Comput. Mater.* **6**, 20 (2020).

The FLARE paper. Introduces a practical on-the-fly active learning
loop: a GP-based MLIP whose Bayesian variance is monitored during
MD; when variance exceeds a threshold, the configuration is
re-labelled with DFT and the model is updated.

**Why read this.** FLARE itself uses a GP; you will use a MACE
committee. The methodological lesson — *uncertainty-triggered
relabelling during MD* — is the same.

**Extract.** The trigger metric (per-atom force standard deviation
or per-cell predicted variance); the cool-down logic (do not retrigger
on every step); the typical AL rounds needed (3–10 for a new system).

---

## 4. Smith, Nebgen, Lubbers et al. (2018) — *Less is more: sampling chemical space with active learning*

> Smith, J. S. et al. *J. Chem. Phys.* **148**, 241733 (2018).

Active learning for ANI (organic molecules), with extensive
diagnostics on training-set size vs accuracy.

**Why read this.** Particularly useful if your target is a molecular
system. The curve of "accuracy vs training-set size after AL" is
plotted in their figure 5 and is approximately the same shape for
all MLIPs.

**Extract.** The plateau of accuracy past ≈ 1000 well-chosen frames;
the importance of *diversity* over *quantity* in the training set;
the cost-benefit analysis of AL vs brute-force sampling.

---

## 5. Cheng, Engel, Behler, Dellago, and Ceriotti (2019) — *Ab initio thermodynamics of liquid and solid water*

> Cheng, B. et al. *PNAS* **116**, 1110 (2019).

Trains a Behler–Parrinello NNP on water with thorough validation
(structure, dynamics, phase behaviour). Demonstrates the *whole
pipeline* you are attempting.

**Why read this.** The most thorough single-paper example of MLIP
training + validation + production. Especially useful for the
validation section: which thermodynamic observables to compute and
how to compare them to AIMD.

**Extract.** The validation hierarchy (parity → RDF → dynamical
quantities → phase behaviour); the realisation that energy errors
of < 1 meV/atom are usually sufficient if forces are good; the
discussion of what "agreement with AIMD" should mean quantitatively.

---

## 6. Kovács, Batatia, Arany, and Csányi (2023) — *Evaluation of the MACE Force Field Architecture*

> Kovács, D. P. et al. *J. Chem. Phys.* **159**, 044118 (2023).

Practical benchmarks of MACE across many systems with default
hyperparameters.

**Why read this.** Tells you what hyperparameters *not* to tune. The
defaults are good. Resist the urge to grid-search.

**Extract.** The recommended defaults; the sensitivity surface; the
training-set-size scaling.

---

## 7. Batatia, Benner, Chiang, Elena, Kovács, Riebesell et al. (2023) — *A foundation model for atomistic materials chemistry* (MACE-MP-0)

> Batatia, I. et al. arXiv:2401.00096 (2023).

The MACE-MP-0 foundation model trained on Materials Project. Useful
as a *baseline* (does it cover your system?) and as a *fine-tuning
starting point* if you choose to fine-tune rather than train from
scratch.

**Why read this.** Establishes when fine-tuning a foundation model
beats training from scratch (usually when training-set is small,
< 200 frames). For an undergraduate project, fine-tuning the
foundation model is a defensible alternative; document the choice.

**Extract.** The reported per-system MAE table; the chemistry
coverage; the fine-tuning recipe.

---

## 8. Csányi, Albaret, Payne, and De Vita (2004) — *"Learn on the fly": A hybrid classical and quantum-mechanical molecular dynamics simulation*

> Csányi, G. et al. *Phys. Rev. Lett.* **93**, 175503 (2004).

The original "learn on the fly" paper, conceptually the ancestor of
modern active learning.

**Why read this.** Helps you place AL in historical context and
emphasises that AL is not new: it is the principled solution to the
sampling problem.

**Extract.** The decomposition of the simulation cost into "cheap
when uncertain, expensive when needed"; the conceptual basis for
modern uncertainty-aware training.

---

## 9. Stocker, Gasteiger, Becker, Günnemann, and Margraf (2022) — *How robust are modern graph neural network potentials in long and hot molecular dynamics simulations?*

> Stocker, S. et al. *Mach. Learn.: Sci. Technol.* **3**, 045010 (2022).

A systematic stress-test of GNN potentials under long, hot MD. Most
GNN-MLIPs fail by ≈ 100 ps at high temperature. Quantifies the
failure modes and the relationship to training-set distribution.

**Why read this.** This paper validates the failure modes you will
*experience*: the model gives good parity statistics but explodes
in MD. Read it before your first failure so you are prepared.

**Extract.** The failure-mode taxonomy (energy drift, atom overlaps,
unphysical reconstruction); the recommended remedies (more diverse
training, lower temperature for the first MD, AL).

---

## 10. Bartók, De, Poelking et al. (2017) — *Machine learning unifies the modeling of materials and molecules*

> Bartók, A. P. et al. *Sci. Adv.* **3**, e1701816 (2017).

A panoramic review of MLIPs across systems. Gives you perspective on
what a "good" MLIP looks like across very different chemistries.

**Why read this.** Calibrates your expectations. A well-trained MLIP
on a 2D dichalcogenide will not have the same MAE as on bulk Si; the
two are not directly comparable.

**Extract.** The cross-system comparison; the relationship between
chemical complexity (number of elements, coordination diversity) and
required training-set size.

---

## Optional eleventh — if you are doing a molecular target

For molecular systems, see Smith, Isayev, and Roitberg (2017) on
ANI-1, *Chem. Sci.* **8**, 3192. Useful for understanding the
sampling philosophy on molecular conformations.

---

## Synthesis exercise

After reading, write a one-page memo:

1. Why is your chosen system not adequately covered by MACE-MP-0
   (or whatever foundation model you tested)? Show evidence.
2. What sampling protocol will you use, and what does each component
   contribute?
3. What MD-failure-mode is most likely for your system, and what
   active-learning trigger will catch it?
4. What is your validation hierarchy, in order of increasing
   stringency?

If you cannot answer these, re-read papers 1, 3, and 9 before
generating any data.
