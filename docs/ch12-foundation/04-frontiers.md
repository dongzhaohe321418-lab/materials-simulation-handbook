# 12.4 Frontiers — What Comes Next

The two preceding sections described tools that, by the end of 2025,
had matured to the point of routine use. This section is more
speculative. It surveys the directions in which the foundation-model
programme is actively expanding, the open problems that remain
unsolved, and the recent literature that a reader who wants to follow
this field should know about.

## Multimodal models

The foundation models of Sections 12.2 and 12.3 operate on a single
modality: atomic structure. A growing body of work argues that the
real leverage will come from models that simultaneously ingest
multiple representations of the same material — structure, computed
spectra, experimental measurements, text from the literature — and
build a joint embedding space across them.

The motivating analogy is CLIP. By training on $400$ million
image–caption pairs with a contrastive loss, CLIP produces a joint
embedding space in which an image of a cat and the text "a cat" land
in nearby points. Once such a space exists, an enormous range of
downstream tasks becomes trivial: zero-shot classification (compare
to text-template embeddings), retrieval (find images whose embedding
is closest to a query), even generation (used as a critic in
DALL-E / Stable Diffusion).

In materials science the analogous data would be:

- **Structure** — the same graph representation used by MACE or
  MatterGen.
- **Computed properties** — DFT-derived band structure, density of
  states, phonon spectrum, optical response.
- **Computed spectra** — simulated XRD, infrared, Raman, NMR, XANES.
- **Experimental measurements** — measured XRD patterns from ICSD or
  in-house databases; reported band gaps; magnetic susceptibilities.
- **Literature text** — abstracts and full text of materials papers,
  including the human descriptions of the systems studied.

A foundation model trained to align these modalities — in the manner
of CLIP, or of more recent multimodal transformers — would in
principle let a user issue queries of the form "find structures whose
XRD pattern matches this measurement, predicted to have a band gap
above $2$ eV, that have not been previously reported". Several
research groups (notably at MIT, Microsoft and DeepMind) have
prototypes in this direction, and the public release of such a model
seems likely within the next two years.

What makes the problem genuinely hard is the *misalignment* of the
modalities. The same crystal has many possible structural
representations (primitive cell, conventional cell, supercell), each
giving the same physical predictions but different graph inputs. Its
computed properties depend on the level of theory used. Its
experimental measurements depend on the sample and on the
measurement protocol. Its literature description is in human-written
prose, with terminology that varies between sub-fields. A useful
multimodal model must be invariant to all of these incidental
variations while remaining sensitive to the underlying physical
identity. This is harder than the corresponding problem for natural
images, where the modes of variation are more constrained.

## Autonomous laboratories and the closing loop

A parallel line of work, less visible to the modelling community but
arguably more consequential, is the construction of *autonomous
laboratories* — robotic platforms that synthesise, characterise and
return data on materials proposed by a computational pipeline, with
no human in the inner loop.

The exemplar is **A-Lab** at Lawrence Berkeley National Laboratory,
which in late 2023 reported the autonomous synthesis of $41$ novel
inorganic compounds in $17$ days, with structures proposed by a
combination of generative models and active-learning surrogates. The
A-Lab pipeline integrates powder XRD, mass spectrometry, and a
robotic precursor-mixing platform; failures (most attempts do fail)
are fed back into the model as negative examples.

Several similar platforms exist or are under construction:

- **The Self-Driving Lab** consortium (University of Toronto, NRC
  Canada) — focused on optoelectronic and catalytic materials, with
  several active loops on organic photovoltaics and reductive
  electrocatalysts.
- **MAP** at the University of Liverpool — focused on porous
  framework materials, with a particularly strong integration of
  computational screening and synthesis.
- **Aroyo** at Argonne — focused on battery cathodes, with in situ
  electrochemical characterisation.

The conceptual significance of these platforms is that they close
the simulation-experiment loop in a way that has been promised since
the inception of the Materials Genome Initiative in 2011. A
generative model proposes; an MLIP screens; DFT verifies; a robot
synthesises; characterisation returns measured properties; the loop
re-iterates. The cycle time, on these platforms, is days to weeks,
not the years that a graduate-student-mediated loop typically takes.

The bottleneck, as of 2026, is no longer the modelling. It is the
*synthesis* — the long tail of compositions for which no robust
robotic protocol exists. Most autonomous labs have a relatively
narrow window of chemistry within which they can operate (specific
solvent classes, specific temperature ranges, specific characterisation
tools), and the foundation models must learn to propose candidates
inside that window. The mismatch between what the models think is
"interesting" and what the lab is "capable of" is now an active area
of methodological work.

## What the foundation models still cannot do

Notwithstanding the impressive recent progress, several physical
phenomena remain stubbornly outside the reach of universal MLIPs and
generative models. Honest accounting of these gaps is essential.

### Long-range interactions

All current universal MLIPs use a finite cutoff (typically $5$–$6$ Å).
This is adequate for short-range covalent and metallic bonding but
fails to capture phenomena that depend on long-range Coulomb forces:
the dielectric response of ionic crystals, ferroelectric phase
transitions, the energetics of charged defects, polaronic states. The
LO–TO splitting in polar phonon spectra, for instance, requires the
non-analytic correction associated with the macroscopic dipole — and
no purely short-range MLIP can produce it.

Active research is exploring three remedies:

- **Hybrid models** that combine a short-range MLIP with an explicit
  Ewald-summed Coulomb term, with charges either fixed (using
  oxidation-state assignments) or predicted by the model itself.
- **Long-range MLIPs** that pass messages over many more layers,
  enlarging the effective receptive field. This has the cost of
  increased inference time and reduced numerical stability.
- **Foundation models with self-consistent charge equilibration**,
  in which atomic charges are themselves variables that respond to
  the local environment. The 4G-HDNNP architecture (Ko et al.,
  2021) and its successors are the canonical references.

None of these has been fully integrated into a public universal MLIP
as of mid-2026, but the field is moving quickly.

### Charge transfer

A closely related problem. Reactions involving electron transfer
between atoms — redox processes in batteries, electrochemical catalysis,
photoexcitation — depend on the redistribution of electronic charge,
which a purely local MLIP cannot represent. CHGNet's magnetic-moment
head is a partial concession to this problem; explicit charge models
are the natural next step.

### Excited states

Almost every MLIP discussed in this book is fitted to ground-state
DFT data, and the resulting potential is appropriate for ground-state
dynamics only. Excited-state phenomena — photochemistry, exciton
dynamics, non-adiabatic relaxation — are outside its remit. A
foundation model for excited-state dynamics would require training on
TDDFT, MRCI or $GW$+BSE data; such datasets exist but are perhaps
three orders of magnitude smaller than MPtrj. Progress here is gated
by data, not by architecture.

### Magnetic ordering

Universal MLIPs typically treat the magnetic state as a property of
the configuration (CHGNet) or simply ignore it (MACE-MP-0). Neither
approach captures the full complexity of frustrated magnetism, where
the energy landscape depends sensitively on the spin configuration as
much as on the atomic positions. A magnet-aware foundation MLIP would
need either to enumerate spin states explicitly or to predict them as
part of its output. The Heisenberg-MACE work of late 2024 is a step
in this direction; the field is unsettled.

### True generalisation across the periodic table

The implicit hope of the universal-MLIP programme is that a model
trained on a comprehensive corpus will generalise to *any* chemistry
within the periodic table. Empirical evidence is mixed. On elements
well-represented in the training data, generalisation is strong; on
under-represented elements (actinides, several lanthanide oxidation
states), the model is closer to an interpolator confined to the
training distribution.

A more honest statement is that the current foundation MLIPs
generalise well within the *chemistry* of their training data and
moderately well across the *configurations* (different temperatures,
pressures, defect states) of those chemistries. They do not yet
generalise across genuinely new chemistry, and there is no a priori
reason to expect that they will without further data.

## A reading list, 2024–2026

The literature in this area is moving fast enough that any printed
list is partially obsolete by the time it appears. The following are
the papers that, as of mid-2026, the editors believe are most worth
reading. Each is annotated with one or two sentences indicating what
the reader should expect to learn.

### Universal MLIPs

- **Batatia, I. et al.** *A foundation model for atomistic materials
  chemistry.* arXiv:2401.00096 (2024). The MACE-MP-0 paper.
  Definitive reference for the model, the training data, and the
  out-of-the-box performance benchmarks.

- **Deng, B. et al.** *CHGNet as a pretrained universal neural network
  potential for charge-informed atomistic modelling.* Nature Machine
  Intelligence (2023). The CHGNet paper. Demonstrates the value of
  the charge head for redox-sensitive applications.

- **Park, S. et al.** *Scalable parallel algorithm for graph neural
  network interatomic potentials in molecular dynamics simulations.*
  Journal of Chemical Theory and Computation (2024). The SevenNet
  paper. Strong on the engineering side.

- **Neumann, M. et al.** *Orb: A fast, scalable neural network
  potential.* arXiv:2410.22570 (2024). The Orb paper, with strong
  out-of-distribution results.

- **Barroso-Luque, L. et al.** *Open Materials 2024 (OMat24):
  Inorganic materials dataset and models.* arXiv:2410.12771 (2024).
  Describes the OMat24 dataset and the EquiformerV2-OMat model.

### Generative models

- **Zeni, C. et al.** *MatterGen: A generative model for inorganic
  materials design.* Nature 639, 624–632 (2025). The MatterGen paper.

- **Xie, T. et al.** *Crystal diffusion variational autoencoder for
  periodic material generation.* ICLR (2022). The CDVAE paper,
  predecessor to MatterGen.

- **Jiao, R. et al.** *Crystal structure prediction by joint
  equivariant diffusion.* NeurIPS (2023). The DiffCSP paper.

- **Antunes, L. M. et al.** *Crystal structure generation with
  autoregressive large language modeling.* Nature Communications
  (2024). A surprising demonstration that simple text-based models
  on CIF strings work for many problems.

### Autonomous labs and the experimental loop

- **Szymanski, N. J. et al.** *An autonomous laboratory for the
  accelerated synthesis of novel materials.* Nature 624, 86–91
  (2023). The A-Lab paper.

- **Burger, B. et al.** *A mobile robotic chemist.* Nature 583,
  237–241 (2020). An early but conceptually important demonstration
  in chemistry.

- **MacLeod, B. P. et al.** *Self-driving laboratory for accelerated
  discovery of thin-film materials.* Science Advances 6, eaaz8867
  (2020).

### Long-range and charge-transfer extensions

- **Ko, T. W. et al.** *A fourth-generation high-dimensional neural
  network potential with accurate electrostatics including non-local
  charge transfer.* Nature Communications 12, 398 (2021). 4G-HDNNP.

- **Cheng, B.** *Cartesian atomic cluster expansion for machine
  learning interatomic potentials.* npj Computational Materials
  (2024). Discussion of long-range corrections in the ACE
  framework.

### Reviews

- **Schmidt, J. et al.** *Recent advances and applications of machine
  learning in solid-state materials science.* npj Computational
  Materials 5, 83 (2019). The standard mid-decade review.

- **Choudhary, K. et al.** *Recent advances and applications of deep
  learning methods in materials science.* npj Computational Materials
  8, 59 (2022). Updated and broader.

- **Friederich, P. et al.** *Machine-learned potentials for next-
  generation matter simulations.* Nature Materials 20, 750–761 (2021).
  A balanced review of where the field stood at the dawn of the
  foundation-model era.

### Critiques and cautionary notes

- **Riebesell, J. et al.** *Matbench Discovery — A framework to
  evaluate machine learning crystal stability predictions.* arXiv:
  2308.14920 (2023, updated 2024). The leaderboard, but more
  importantly the careful methodology that exposed several previously
  popular models as having serious distributional weaknesses.

- **Stocker, S. et al.** *How robust are modern graph neural network
  potentials in long and hot molecular dynamics simulations?* Machine
  Learning: Science and Technology (2022). An early but still
  pertinent critique of stability under aggressive sampling.

## A summary of where we are

The honest summary, as the chapter closes, is something like this.

Foundation models for materials science exist, they work, and they
have already changed the practical workflow of computational
materials research. A typical 2026 study begins with a pre-trained
MLIP, fine-tunes it on a small system-specific dataset, uses it to
explore a phase space or candidate set, and verifies the most
interesting findings with DFT. This pipeline is faster, more general,
and more accessible than what was possible even three years ago.

At the same time, the field is not done. Long-range interactions,
charge transfer, excited states, magnetic ordering, and genuinely
out-of-distribution chemistry remain open. Generative models produce
candidates, but synthesisability is not learned. Autonomous
laboratories close one part of the loop, but only within the chemistry
they can robotically handle. The next decade of materials simulation
will be defined less by individual algorithms than by the integration
of these pieces into reliable, end-to-end pipelines — and by the
careful identification of the cases where the foundation breaks down.

The book closes here, but the work does not. The appendices that
follow consolidate the mathematical, computational and bibliographic
resources used throughout. The reader who has followed Chapters 1
through 12 is now equipped, we hope, both to use the tools described
and to read the literature critically as it appears.
