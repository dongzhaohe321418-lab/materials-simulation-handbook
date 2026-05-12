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

## Three open problems in detail

The brief survey above identified several gaps in the current
foundation-model toolkit. Three of them deserve closer scrutiny
because they are *physical* limitations of the dominant
architectures, not merely engineering debt that will be cleared
by larger datasets. Each is the subject of substantial recent
work, none has a clean solution, and each is likely to define a
sub-field of materials machine learning for the next several years.

### Long-range Coulomb in MLIPs

**The problem.** Every universal MLIP described in this book uses a
finite cutoff radius — typically $5$–$6$ Å — beyond which atomic
interactions are simply truncated. For short-range covalent and
metallic chemistry this is excellent: the energy landscape is set
by the local bonding environment, and the truncation error decays
exponentially with the cutoff. For systems where long-range
electrostatics matters, the truncation is qualitatively wrong.

The clearest diagnostic is the dielectric screening of a polar
crystal. The macroscopic Born effective charges and the LO–TO
splitting in the phonon spectrum both arise from the non-analytic
behaviour of the dynamical matrix at small wavevector,
$\mathbf{q} \to 0$, where Coulomb forces are unscreened. A purely
local MLIP cannot reproduce this: the small-$\mathbf{q}$ dynamical
matrix is built from far-apart-atom-pair contributions that lie
entirely outside the model's receptive field. Polar phonon spectra,
ferroelectric instabilities, dielectric constants, and the energetics
of charged defects in ionic solids all carry the same signature
error.

**Why current methods fail.** The dominant message-passing
architectures (MACE, NequIP, SevenNet) gather information through
local convolutions over a graph defined by the cutoff. Stacking
more layers extends the *effective* receptive field but at quadratic
cost in number of operations and worsening numerical stability of
the gradient through many message-passing steps. To capture
Coulomb forces with the correct $1/r$ tail requires either an
analytical long-range term — an Ewald-summed Coulomb interaction
layered on top of the MLIP, with charges as additional outputs —
or a fundamentally different architecture (e.g., one that operates
in reciprocal space).

The 4G-HDNNP architecture of Ko et al. (2021) handles long-range
electrostatics by *learning* the atomic charges and adding an
explicit Ewald sum. Recent universal-MLIP descendants (MACE-LR,
Cheng's Cartesian-ACE long-range extension, the eqV2 long-range
variant in development) follow the same pattern. None of these is
yet in widespread use as a drop-in foundation model; the engineering
of charge equilibration in a parameter-sharing setting across $96$
elements is genuinely difficult.

**What success would look like.** A universal MLIP that reproduces,
within DFT accuracy, the Born effective charges and LO–TO splittings
of a few canonical polar oxides (MgO, SrTiO$_3$, BaTiO$_3$, PbTiO$_3$),
the dielectric constants of layered materials, and the formation
energies of charged defects relative to neutral references. As of
mid-2026 no public model has cleared this bar; the published
benchmarks compare neutral, near-equilibrium configurations almost
exclusively.

**Recent attempts.** Anstine and Isayev (2023) survey the field;
Loche et al. (2024) report a long-range MACE variant that reproduces
the dielectric response of BaTiO$_3$ within $\sim 15$%; the
Cheng-group Cartesian-ACE long-range work (npj Computational
Materials, 2024) shows promising results on simple ionic solids but
has not yet been scaled to a universal model.

### Charged systems and electrochemistry

**The problem.** Every universal MLIP discussed here treats the
total charge of the simulation cell as fixed at zero. This is
implicit in the architecture — atoms have positions and species, no
spin or charge — and in the training data (DFT calculations almost
always at charge neutrality). For neutral, isolated systems, this is
the right choice. For electrochemistry, where an electrode is in
contact with an electrolyte at a controlled electrochemical
potential, it is fundamentally wrong.

Electrochemical interfaces have three properties no current universal
MLIP can represent. First, the electrode can carry a net charge,
balanced by the diffuse double-layer charge in the electrolyte;
energetics depend on this charge. Second, the electrochemical
potential, not the chemical composition, is the natural control
variable, requiring grand-canonical sampling in the electronic
degrees of freedom. Third, the system can transfer electrons across
the interface during reaction events, which neither a fixed-charge
nor a per-atom-charge model captures.

**Why current methods fail.** The natural way to extend an MLIP to
charged systems is to add an extra input channel (the total charge,
or the chemical potential) and retrain. The problem is that DFT
calculations of charged supercells require careful handling of the
compensating uniform-jellium background and finite-size corrections
(Makov-Payne, Freysoldt-Neugebauer); a foundation MLIP trained on
charged-system DFT data must learn to undo these corrections rather
than fold them into a meaningful physical signal. No public dataset
of charged-supercell DFT energies at the scale of MPtrj exists.

A complementary issue: the per-atom-charge-prediction approach
(CHGNet's magmom head extended to electronic charge) cannot represent
electron transfer because it does not change the *total* number of
electrons. Some configurations have only one charge-state assignment
consistent with chemistry; others (mixed-valent oxides, partially
reduced oxides) admit several, and the model has no principled way
to choose.

**What success would look like.** A universal MLIP that can simulate
an electrode in contact with an electrolyte at a specified potential
versus a reference (SHE or Ag/AgCl), correctly predicting the
double-layer capacitance, the surface-charge density as a function
of potential, and the redox potentials of small molecules at the
interface within $0.2$–$0.3$ V of experiment. None exists as of
mid-2026.

**Recent attempts.** The CENT model (Ghasemi et al. 2015, with
recent extensions) uses charge equilibration with electronegativity
parameters and is in principle applicable to charged systems, though
not at universal-MLIP scale. The constant-potential MD method of
Bonnet, Otani and Sugino (and recent extensions) is an explicit
DFT-based approach; a learned surrogate is an obvious target. The
group around Behler has reported preliminary work on
charge-conditioned 4G-HDNNP at mid-scale; the GitHub repository for
universal charge-aware models is still small. This is, in our view,
the most consequential open problem for foundation MLIPs and is
likely to attract the next generation of architectural innovation.

### Excited states and TDDFT-ML

**The problem.** Every MLIP discussed in this book is trained on
ground-state DFT, which is exact in principle for the ground-state
energy and density of a many-electron system. Many of the most
interesting physical phenomena — photoexcitation, exciton transport,
non-radiative relaxation, photochemistry — involve excited electronic
states. The ground-state PES is irrelevant for these processes; what
is needed is at least the lowest few excited-state PESs and the
non-adiabatic couplings between them.

**Why current methods fail.** The training data for an excited-state
MLIP would need to come from a method capable of computing excited
states accurately. TDDFT (the most accessible option) is roughly an
order of magnitude more expensive than ground-state DFT for the
same system. CASSCF and CASPT2 (for stronger correlation) are two to
three orders more expensive. $GW$+BSE (the closest to experiment for
solid-state spectra) is two to four orders more expensive and has
strong system-size limitations. Generating a dataset comparable to
MPtrj — $1.6 \times 10^6$ configurations — at any of these levels is
prohibitively expensive, and would have to be redone if the target
spectroscopy changes.

The architectural problem is also non-trivial. An MLIP predicts a
*scalar* energy per configuration. An excited-state model would need
to predict a *set* of energies (one per state), the gradients with
respect to atomic positions (one per state), and ideally the
non-adiabatic couplings (matrices between states). The output
dimension grows linearly with the number of states tracked, and the
ordering of states can change with configuration (state crossings,
conical intersections), making the prediction problem
combinatorially harder than the ground-state case.

**What success would look like.** A foundation model that, given a
molecular or crystalline configuration, returns the lowest $N$
adiabatic state energies and gradients with chemical accuracy
($\sim 0.1$ eV) for at least the chemistries of small organic
chromophores, transition-metal complexes, and a few canonical
inorganic semiconductors. With this in hand, photochemistry
simulations would become tractable at MLIP cost, and the screening
of light-harvesting materials, photovoltaic absorbers and
non-linear-optical compounds would accelerate by orders of
magnitude.

**Recent attempts.** SchNarc (Westermayr et al. 2020) and its
descendants have demonstrated excited-state neural-network
potentials for small molecules, with two to ten states tracked and
non-adiabatic dynamics simulated end-to-end. SPaiNN (Kornbluth et
al. 2024) extends this to equivariant architectures. The 2025
review by Kasper et al. (J. Phys. Chem. Lett.) surveys the
landscape; the universal-foundation-model regime is years away, and
even per-system models are far from the maturity of ground-state
MLIPs. This is a problem where the architecture is interesting but
the data is the limiting factor.

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
