# 2.4 The Software Ecosystem

A working computational materials scientist in 2026 uses roughly a dozen pieces of software in a typical workflow. None of them does everything. None of them is universally adored. Choosing the right stack — and learning enough of each tool to use it well — is one of the unglamorous skills that separates a productive researcher from someone who spends six months wrestling with file formats. This section is a tour of the landscape: who the major players are, what they do well, where they fall short, and what we suggest you install first.

We will move from low scales to high, mirroring Section 2.1, and finish with a recommended starter stack.

## DFT codes

The DFT ecosystem is the most fragmented part of the field. There are roughly twenty actively maintained codes, each with a different combination of basis set, pseudopotential approach, parallelisation strategy, and licence. The choice matters: a calculation that takes a day in one code may take a week in another, and certain physical questions favour particular basis sets.

### Plane-wave codes

Plane-wave codes expand the wavefunction in a basis of plane waves $e^{i \mathbf{G} \cdot \mathbf{r}}$ defined by the reciprocal lattice. They are the natural choice for periodic solids and dominate solid-state DFT.

**VASP** (Vienna Ab initio Simulation Package). Commercial, licensed per group. The most widely cited DFT code in solid-state physics. Mature, well-tested pseudopotential library (PAW), excellent performance on transition metals, extensive feature set (hybrid functionals, GW, BSE, AIMD, NEB, DFPT). Cons: licence cost, closed source. If your supervisor uses VASP, you will probably use VASP.

**Quantum ESPRESSO** (QE). Open-source, GPL-licensed, GPU-accelerated. The standard open-source plane-wave code. Modular: `pw.x` does self-consistent DFT, `ph.x` does phonons, `dos.x` and `bands.x` do post-processing. Slightly less polished than VASP but free, and the community is large.

**ABINIT**. Open-source, well-developed for GW and BSE. Strong for response functions and excited states. The community is smaller but the code is robust.

**CASTEP**. Commercial but free to UK academics. Closely tied to the Cambridge solid-state community.

### Local-orbital and Gaussian-basis codes

Local-orbital codes use atom-centred basis functions, which are economical for systems with vacuum (clusters, molecules, surfaces with thick vacuum layers) but require careful basis-set convergence.

**FHI-aims**. Numeric atom-centred orbitals; all-electron; full relativistic options. Excellent for surfaces, molecules, and systems containing heavy elements. Commercial licence with academic version.

**SIESTA**. Numerical atomic orbitals, norm-conserving pseudopotentials, linear-scaling option. Open source. Strong in the Spanish-speaking community.

**CP2K**. Gaussian basis with a plane-wave auxiliary basis (the GAPW or GPW method). Particularly fast for AIMD of liquids and biomolecules; widely used for water, electrolytes, and large organic systems. Open source.

**Gaussian**. The standard for molecular quantum chemistry. Commercial. Periodic capabilities are limited.

### Real-space and other codes

**GPAW**. Real-space grids combined with the projector-augmented-wave method. Open source, Python-friendly, integrates tightly with ASE (which it was originally written alongside). Great for prototyping and small systems. Less performant than VASP or QE for large simulations.

**OCTOPUS**. Real-space TDDFT. The standard for excited-state dynamics in finite systems.

### Which to choose

For most solid-state problems, a plane-wave code is the default. If you have a VASP licence, use VASP. If not, use Quantum ESPRESSO. If you are doing GW or BSE seriously, consider ABINIT. If your system has substantial vacuum or you care about all-electron accuracy, consider FHI-aims. If you are doing AIMD of a liquid, consider CP2K. If you want to script everything from Python and don't mind a performance penalty, use GPAW.

## Molecular dynamics codes

The MD ecosystem is less fragmented but more domain-specialised than DFT. Different communities — soft matter, biomolecular, hard-condensed matter — favour different codes.

**LAMMPS** (Large-scale Atomic/Molecular Massively Parallel Simulator). Open-source, BSD-licensed. The dominant code for hard-condensed-matter MD. Supports an enormous range of force-field types (Lennard-Jones, EAM, Tersoff, ReaxFF, COMB, and crucially the modern MLIPs via the ML-IAP package, ML-PACE, ML-MACE plugins). Scales to billions of atoms on clusters. The configuration language is idiosyncratic but powerful.

**GROMACS**. Open-source. The dominant code for biomolecular MD. Highly optimised for the standard biomolecular force fields (AMBER, CHARMM, OPLS, GROMOS). Used routinely for protein folding and membrane simulations. Less natural for crystalline solids.

**AMBER**. Both a force field and a code. Standard in computational chemistry for nucleic acids. Commercial.

**NAMD**. Open-source, scalable, biomolecular focus. Strong for very large biomolecular systems (whole viruses).

**OpenMM**. Python-friendly biomolecular MD with GPU acceleration. Often used as a backend for higher-level Python tools.

**JAX-MD** and **TorchSim**. Differentiable MD frameworks. Useful for end-to-end differentiable workflows (training potentials, differentiable design). Performance has improved dramatically but is still below LAMMPS for routine work.

**ASE**. Not a code in the same sense — see below — but ASE provides a Python interface to many of the above and can run MD natively for small systems.

### Which to choose

For solids and inorganic materials, LAMMPS is the standard. For proteins and aqueous systems, GROMACS. For Python-native workflows on small systems, ASE or OpenMM. For differentiable MD, JAX-MD.

## Workflow tools

Workflow tools sit above individual codes and provide a unified interface for setting up, running, and analysing calculations. They are not strictly necessary for one-off calculations, but they are essential for any serious project that runs more than a handful of jobs.

**ASE** (Atomic Simulation Environment). A Python library that provides an `Atoms` object representing a molecule or crystal, and *calculators* that connect to dozens of DFT and MD codes. Builds structures (bulk, surface, molecule), writes input files, parses output, and exposes a consistent API for energies, forces, and stresses. ASE is the lingua franca of the Python materials world. Every modern handbook example, including most of this one, uses ASE.

**pymatgen** (Python Materials Genomics). A Python library developed at LBNL alongside the Materials Project database. Strong on structure manipulation, symmetry analysis, phase diagrams, electronic-structure analysis (especially for VASP). Overlaps with ASE in functionality but has different idioms. Many practitioners use both.

**AiiDA**. A workflow management system that records provenance — every calculation's inputs, outputs, and the chain of dependencies — in a database. Built for reproducibility at scale. Used by the European MaX centre and the Materials Cloud. Heavier learning curve than ASE; pays off when you run thousands of calculations.

**FireWorks**. A workflow engine from LBNL, used internally by the Materials Project. Less common outside that orbit.

**atomate** and **atomate2**. Pre-built workflows on top of pymatgen and FireWorks/jobflow for common DFT tasks.

### Recommendation

Learn ASE first; it is the gentlest introduction to programmatic materials simulation. Add pymatgen when you start working with Materials Project data. Adopt AiiDA only if you are running large, long-lived workflows and need provenance.

## Machine-learning stacks

The ML side of materials simulation has changed dramatically in five years. The current landscape splits into general-purpose ML libraries and materials-specific frameworks built on top of them.

### General-purpose

**PyTorch**. The dominant deep-learning framework in research. Most modern MLIPs (NequIP, Allegro, MACE, ORB) are written in PyTorch.

**JAX**. Functional, differentiable numerical computation with strong autodiff and parallelisation primitives. The basis of several recent MLIPs (Allegro-JAX, SchNet-JAX) and of differentiable MD frameworks. JAX is gaining ground especially in physics-aware ML.

**PyTorch Geometric** and **DGL**. Graph neural network libraries used as building blocks for graph-based MLIPs and property predictors.

**e3nn** (PyTorch) and **e3nn-jax**. Libraries for $E(3)$-equivariant neural networks. The foundation of modern equivariant MLIPs.

### MLIP frameworks

**MACE**. State-of-the-art equivariant MLIP based on the Atomic Cluster Expansion (ACE). The foundation model **MACE-MP-0** covers the periodic table for general inorganic chemistry; MACE-OFF covers organic chemistry. Open source, used in research and increasingly in industry.

**NequIP** and **Allegro**. Earlier-generation equivariant MLIPs from the same family. NequIP is fully equivariant; Allegro is strictly local and scales better to large systems.

**SchNet**, **DimeNet**, **GemNet**. Older invariant or message-passing GNNs; still useful and well-documented.

**SevenNet**, **ORB**. Newer entrants competing with MACE on benchmarks; ORB in particular is positioned as a fast foundation MLIP.

**FAIR-Chem** (formerly Open Catalyst Project / ocp). Meta's stack and the OC20/OC22 benchmark and models. The largest catalysis dataset and the source of several leading models (EquiformerV2, etc.).

### Materials databases

**Materials Project** (MP). Roughly 200000 inorganic crystals with DFT-computed properties, an API, and a web interface. The single most important materials database. Open access.

**OQMD** (Open Quantum Materials Database). Northwestern; a similar effort with different functional choices.

**AFLOW**. Duke; high-throughput DFT with strong electronic-property focus.

**NOMAD**. European repository accepting raw outputs from many DFT codes; focused on FAIR data principles.

**Alexandria**. Newer database with $\sim 4$ million DFT-relaxed structures using PBE and PBEsol; valuable for training ML models.

**OC20** and **OC22**. Open Catalyst datasets containing DFT-computed adsorption energies and trajectories for catalysis.

For most users, Materials Project is the entry point, accessed via its Python API (`mp-api`) or pymatgen.

## Visualisation

**VESTA**. The standard for visualising crystal structures and volumetric data (charge densities, ELF). GUI-driven, lightweight, free for non-commercial use. Most figures of crystal structures in the literature were made with VESTA.

**OVITO**. The standard for visualising MD trajectories. Powerful analysis tools (common neighbour analysis, dislocation extraction, defect identification). Free basic version; pro version for commercial use.

**ASE GUI**. The simplest viewer; bundled with ASE, runs from `ase gui structure.xyz`. Useful for quick inspection.

**Jmol**, **Avogadro**, **PyMOL**. Other viewers; PyMOL dominates biomolecular visualisation.

**py3Dmol** and **NGLView**. Browser-based viewers for embedding structures in Jupyter notebooks. Essential for any tutorial that lives in a notebook.

## A recommended starter stack

If you are new to the field and want a starter stack that will let you do useful work within a week, here is our recommendation.

**Operating system.** Linux. Windows Subsystem for Linux is acceptable; macOS works for development and small jobs but most clusters run Linux.

**Python distribution.** Miniforge or pixi for environment management. Avoid the system Python.

**Core Python libraries.**

```bash
pip install ase pymatgen numpy scipy matplotlib jupyter
```

**DFT.** Install Quantum ESPRESSO from your package manager (`apt`, `dnf`, or conda) for local prototyping. Use the production cluster's installation for real calculations.

**MD.** Install LAMMPS, ideally with the ML-IAP plugin compiled in. Conda packages exist.

**MLIPs.** Install MACE for the foundation model:

```bash
pip install mace-torch
```

You can then load MACE-MP-0 and run DFT-quality MD from Python:

```python
from ase.build import bulk
from ase.md.langevin import Langevin
from ase import units
from mace.calculators import mace_mp

atoms = bulk("Cu", cubic=True).repeat((4, 4, 4))
atoms.calc = mace_mp(model="medium", default_dtype="float32")
dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=300, friction=0.01)
dyn.run(1000)
```

**Visualisation.** VESTA for crystals, OVITO for MD trajectories, ASE GUI for quick checks, py3Dmol for notebooks.

**Database access.** Register at materialsproject.org and install `mp-api`.

**Version control.** Git. Every script you write should live in a git repository from day one.

This stack will run on a laptop and scale, with appropriate cluster access, to most research workflows. The total install time is an hour or two. You will outgrow individual pieces — perhaps you will want VASP for hybrid functionals, or AiiDA for provenance, or JAX for differentiable workflows — but the starter stack will see you through the first year.

!!! tip "Don't pre-emptively learn everything"
    A common failure mode for new students is trying to learn five codes simultaneously before doing any science. Resist this. Install the starter stack, run the tutorials at the end of each subsequent chapter, and add new tools only when a project demands them. Depth in one or two tools beats shallow familiarity with ten.

## Beyond the standard stack

Several emerging tools and platforms are worth knowing about even if you do not use them daily.

**Foundation MLIP services.** Several groups now host APIs serving foundation MLIPs: you POST a structure, get back energies and forces. This is the early form of *materials-as-a-service* and may become the default interface for many users in a few years. Hugging Face is hosting foundation model checkpoints; Orbital Materials, FAIR-Chem, and others host inference endpoints.

**Cloud and credit-based compute.** Several commercial platforms (Modal, RunPod, Lambda Labs, the major cloud providers) make GPU compute accessible without an HPC allocation. For prototyping with MLIPs this is often the fastest path, though it remains expensive for large production runs.

**Reproducibility infrastructure.** AiiDA's archive format, the Materials Cloud, NOMAD's open-data repository, and HuggingFace's datasets API are all converging toward sharing both raw inputs and outputs in standard forms. Increasingly, journals expect this.

**Autonomous laboratories.** A-Lab at LBNL, Berkeley's autonomous catalysis platforms, and similar facilities elsewhere couple simulation to robotic experiment with active learning. These are not yet standard, but they represent a clear direction.

## A historical note

The DFT community is older than the MD community in materials, which is older than the modern ML community. Their software cultures reflect this. DFT codes are mostly Fortran 90 with grumpy mailing lists. MD codes are C++ with friendlier ones. ML codes are Python on Discord and Slack. Bridging these cultures — running a workflow that touches a DFT calculation, an MD simulation, and a neural network — is what tools like ASE and AiiDA exist to do, and what we expect every practising computational materials scientist to be doing routinely within a few years.

With the ecosystem mapped, we close Chapter 2 with [exercises](exercises.md) that test your grasp of the landscape, the scale ladder, and the three diagrams.
