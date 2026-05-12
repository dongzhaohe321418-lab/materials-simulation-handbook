# Appendix D — Further Reading

An annotated bibliography of textbooks, major review articles and
open courseware that the editors consider most useful for the reader
who has finished this book and wants to go further. The list is
organised by topic, following the structure of the book itself.

For each entry we have included a one- to three-sentence note about
what the work is best for and at what level. The classification
"introductory" means a student who has just finished an
undergraduate quantum-mechanics or solid-state course can usefully
read it; "intermediate" means a finishing PhD student or postdoc
will get the most from it; "advanced" means the work is genuinely a
specialist reference.

The full bibliographic details — authors, year, publisher, ISBN where
applicable — are in Appendix E.

---

## D.1 Quantum mechanics

The foundation of everything else in this book. Three textbooks span
the difficulty spectrum from undergraduate to graduate.

**Griffiths, *Introduction to Quantum Mechanics* (3rd edition, 2018).**
Introductory. The standard first textbook on the subject. Clear,
discursive, deliberately physical in its approach. Read it before
attempting any of the others if your quantum mechanics is shaky.

**Shankar, *Principles of Quantum Mechanics* (2nd edition, 1994).**
Introductory to intermediate. More formal than Griffiths, with
careful treatment of the postulates and the linear-algebra structure.
The chapter on path integrals is the best introductory treatment we
know.

**Cohen-Tannoudji, Diu and Laloë, *Quantum Mechanics* (Vols. I-III,
2nd edition, 2019).** Intermediate to advanced. The encyclopaedic
French school text. Every topic is treated thoroughly, every formula
is derived. Useful as a reference even if one does not read it
cover to cover. Volume III, on relativistic and many-body methods,
is particularly relevant to electronic-structure theory.

**Sakurai and Napolitano, *Modern Quantum Mechanics* (3rd edition,
2020).** Intermediate. The graduate-level companion to Griffiths.
The chapters on angular momentum and rotation operators are
essential preparation for any work on equivariance.

---

## D.2 Solid-state physics

For the language of band structures, Brillouin zones, phonons and
all the surrounding crystallography.

**Ashcroft and Mermin, *Solid State Physics* (1976).** Intermediate.
Despite its age, still the standard reference. The treatments of
the free-electron and tight-binding models, of phonons, and of
the de Haas–van Alphen effect remain unsurpassed. The notation is
slightly old-fashioned but mostly faithful to current usage.

**Marder, *Condensed Matter Physics* (2nd edition, 2010).**
Intermediate to advanced. A modernised Ashcroft-and-Mermin with
better treatment of magnetism, superconductivity and soft-matter
topics. The chapters on electronic structure and on phonons are
the natural complement to Chapters 5 and 8 of this book.

**Kittel, *Introduction to Solid State Physics* (8th edition,
2005).** Introductory. Wider in scope, lighter in mathematical
detail than Ashcroft-and-Mermin. Useful for orientation when one
encounters an unfamiliar phenomenon and needs the standard
elementary explanation.

**Grosso and Pastori Parravicini, *Solid State Physics* (2nd
edition, 2014).** Intermediate. Italian school text with strong
emphasis on electronic structure, semiconductors and group-
theoretic methods. The treatment of the **k·p** method and of
spin–orbit coupling is the cleanest we know.

---

## D.3 Density functional theory

The methodological core of Chapters 5 and 6.

**Martin, *Electronic Structure: Basic Theory and Practical
Methods* (2nd edition, 2020).** Advanced. The standard graduate-
level text. Comprehensive, careful, and unafraid to discuss the
practical difficulties and limitations of DFT calculations. The
companion volume *Interacting Electrons* (2016, with Reining and
Ceperley) covers post-DFT methods (GW, BSE, DMFT) at the same
level.

**Parr and Yang, *Density-Functional Theory of Atoms and
Molecules* (1989).** Intermediate to advanced. The quantum-chemistry
view of DFT: rigorous, mathematically dense, and surprisingly
relevant to materials practitioners. The chapter on the chemical
potential is the most insightful discussion of the conceptual content
of DFT we have read.

**Sholl and Steckel, *Density Functional Theory: A Practical
Introduction* (2009).** Introductory. The right book for a graduate
student running their first DFT calculation. Covers VASP-style
plane-wave pseudopotential calculations, convergence testing, and
the standard pitfalls. Light on theory, heavy on practice.

**Burke and Wagner, *DFT in a Nutshell* (International Journal of
Quantum Chemistry, 2013).** Introductory review article. Twelve
pages, free online, and the best short introduction to the conceptual
content of DFT we know. A surprising number of senior researchers
have learned the essentials from this article.

---

## D.4 Molecular dynamics and statistical mechanics

The methodological core of Chapters 7 and 8.

**Frenkel and Smit, *Understanding Molecular Simulation* (3rd
edition, 2023).** Intermediate. The standard text on MD and Monte
Carlo for materials and chemistry. Every common algorithm — Verlet,
Nosé–Hoover, replica exchange, free-energy methods — is presented
with both physical motivation and ready-to-implement detail. The
chapter on free energy is alone worth the price.

**Allen and Tildesley, *Computer Simulation of Liquids* (2nd
edition, 2017).** Intermediate. The classical reference, focused on
liquids and soft matter. Particularly strong on transport properties
and on the practical details of pair correlation functions and
structure factors.

**Tuckerman, *Statistical Mechanics: Theory and Molecular
Simulation* (2nd edition, 2023).** Advanced. Combines a careful
statistical-mechanics development with the corresponding simulation
algorithms. The reader who wants to understand the link between an
ensemble in the textbook sense and a thermostatted MD trajectory
will find it here.

**Chandler, *Introduction to Modern Statistical Mechanics* (1987).**
Intermediate. Slim and elegant; the right text for a working
researcher who wants to refresh the statistical-mechanical
foundations underlying Chapters 8 and 11.

---

## D.5 Machine learning for materials

The field that Chapters 9 through 12 inhabit. The literature is too
new to have settled textbooks; we recommend a mixture of review
articles and pivotal primary papers.

### Reviews

**Schmidt, Marques, Botti and Marques, "Recent advances and
applications of machine learning in solid-state materials science"
(*npj Computational Materials*, 2019).** Intermediate. The standard
mid-decade review of the field. Now somewhat dated on the MLIP
side but still the best overview of property-prediction work.

**Choudhary et al., "Recent advances and applications of deep
learning methods in materials science" (*npj Computational
Materials*, 2022).** Intermediate. Updated and broader, covering
graph neural networks, generative models and the use of foundation
models in their early form.

**Friederich et al., "Machine-learned potentials for next-
generation matter simulations" (*Nature Materials*, 2021).**
Intermediate. A balanced review of MLIPs as they stood just before
the foundation-model era.

**Behler, "Four generations of high-dimensional neural network
potentials" (*Chemical Reviews*, 2021).** Intermediate to advanced.
A pedagogical history of MLIP architectures from the perspective of
one of the founders of the field. Essential reading for anyone who
wants to understand why the modern architectures look the way they
do.

### Primary papers

The following are not pedagogical but are too important to omit. All
are technical research papers; we expect the reader to skim and
return as needed.

**Behler and Parrinello, "Generalized neural-network representation
of high-dimensional potential-energy surfaces" (*Physical Review
Letters*, 2007).** The paper that launched modern MLIPs.

**Bartók, Payne, Kondor and Csányi, "Gaussian approximation
potentials" (*Physical Review Letters*, 2010).** The first
demonstration of GP regression for atomistic energy surfaces.

**Schütt et al., "SchNet: A continuous-filter convolutional neural
network for modeling quantum interactions" (*NeurIPS*, 2017).** The
first deep MLIP, and the architectural ancestor of much of what
followed.

**Batzner et al., "E(3)-equivariant graph neural networks for
data-efficient and accurate interatomic potentials" (*Nature
Communications*, 2022).** The NequIP paper. The first really
convincing demonstration that strict equivariance dramatically
improves data efficiency.

**Batatia et al., "MACE: Higher order equivariant message passing
neural networks for fast and accurate force fields" (*NeurIPS*,
2022).** The MACE paper.

**Batatia et al., "A foundation model for atomistic materials
chemistry" (*arXiv:2401.00096*, 2024).** The MACE-MP-0 paper.

**Zeni et al., "MatterGen: a generative model for inorganic
materials design" (*Nature*, 2025).** The MatterGen paper.

---

## D.6 Foundational data and infrastructure

Three resources that should be in every materials-ML practitioner's
toolkit.

**The Materials Project (materialsproject.org).** Free, open, well
documented. The starting point for almost any high-throughput study
of inorganic materials. The accompanying paper (Jain et al., *APL
Materials*, 2013) is short and well worth reading.

**pymatgen (pymatgen.org).** The Python interface to the Materials
Project and to a large fraction of the related computational
ecosystem. The online documentation is the right entry point.

**ASE — Atomic Simulation Environment (wiki.fysik.dtu.dk/ase).** The
Python framework for atomistic calculations and for integrating
different electronic-structure codes. The standard glue between
DFT codes and MLIPs.

---

## D.7 Open courseware

For readers who learn best through structured lectures.

**Marzari, "Computational Methods in Physics, Chemistry and
Materials Science" (EPFL).** Recorded lecture course, freely
available online. The most pedagogical complete introduction to
periodic DFT we know.

**Cohen-Tannoudji's lectures at the Collège de France.** In French,
but transcribed lecture notes (in French and partly translated)
remain a treasure for advanced topics in quantum mechanics. Search
the Collège de France video archive.

**MIT OpenCourseWare 3.320, "Atomistic Computer Modeling of
Materials".** Older but still excellent introduction to DFT and MD
oriented to materials applications. Lecture notes and problem sets
are available.

**The Coursera course "Materials Informatics" (Schmidt et al.).**
Recent and pitched at the right level for the modern reader. Covers
much of the same ground as Chapters 9–12 of this book.

**Berkeley's "Hands-on tutorials" at the Molecular Foundry.** Short
video tutorials on the practical use of VASP, Quantum ESPRESSO and
the Materials Project. Hosted on YouTube; search for "Foundry
materials tutorial".

---

## D.8 Mathematical and computational background

For the reader who wants to fill in the mathematical or
computational prerequisites.

**Strang, *Linear Algebra and Its Applications* (4th edition, 2006),
and the corresponding MIT OpenCourseWare 18.06.** The standard.
The MIT video lectures are exceptional.

**Trefethen and Bau, *Numerical Linear Algebra* (1997).** For the
reader who wants to understand the numerical methods (Lanczos,
Davidson, conjugate gradients) that lie under the hood of every
DFT code.

**Press, Teukolsky, Vetterling, Flannery, *Numerical Recipes* (3rd
edition, 2007).** A general-purpose reference for the working
scientist who needs an algorithm now. Idiosyncratic and sometimes
opinionated, but vast in scope.

**Bishop, *Pattern Recognition and Machine Learning* (2006).**
Intermediate. Still the best introduction to the probabilistic view
of machine learning, including Gaussian processes (relevant to GAP)
and to the kernel methods underlying much of pre-2020 materials ML.

**Goodfellow, Bengio and Courville, *Deep Learning* (2016).**
Intermediate. The standard introduction to deep learning. Some
chapters (especially on regularisation and on optimisation) are
directly useful for the MLIP practitioner.

**Murphy, *Probabilistic Machine Learning: An Introduction* (2022)
and *Probabilistic Machine Learning: Advanced Topics* (2023).**
Intermediate to advanced. More current than Bishop, and covers
modern topics — diffusion models, transformers, normalising flows —
absent from older texts. The two volumes together are the closest
thing to a modern reference for the field.

---

## D.9 Personal recommendations from the editors

If asked to nominate the *single best companion* to this book for
each major topic, our votes are:

- **Quantum mechanics:** Shankar.
- **Solid-state physics:** Ashcroft and Mermin, supplemented by
  Marder for modern topics.
- **DFT:** Martin's *Electronic Structure* for theory; Sholl and
  Steckel for first practice.
- **MD and statistical mechanics:** Frenkel and Smit, with
  Tuckerman alongside for the foundations.
- **Machine learning for materials:** Schmidt et al. (2019) plus
  Behler (2021) as a starting bibliography; supplement with the
  primary papers in Appendix E.
- **Mathematical background:** Strang for linear algebra; Murphy
  for machine learning.

A reader who has worked through this book and seriously engaged with
these six recommendations should be able to read essentially any
paper in the modern computational materials literature.
