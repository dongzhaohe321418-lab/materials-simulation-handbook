# Chapter 6 — Running DFT in Practice

[Open in Colab](https://colab.research.google.com/github/dongzhaohe321418-lab/materials-simulation-handbook/blob/main/notebooks/ch06-running-dft-colab.ipynb){ .md-button .md-button--primary }

In [Chapter 5](../ch05-dft/index.md) we set up the machinery of Kohn-Sham density functional theory: the energy functional, the self-consistency loop, the exchange-correlation approximations. We treated the calculation abstractly — give it a set of atomic positions, get back an energy.

That abstraction breaks the moment you sit at a terminal. A real DFT calculation requires you to decide on a pseudopotential, a basis set, a plane-wave cutoff, a k-point mesh, a smearing scheme, a starting magnetisation, an SCF convergence threshold, and a dozen other knobs whose default values are wrong about half the time. None of these are theory. All of them affect the answer.

This chapter is a working tutorial. By the end you will have used [Quantum ESPRESSO](https://www.quantum-espresso.org/) to compute the band structure of silicon from scratch, and you will know — quantitatively — how converged that band structure is. The same workflow extends to almost every solid you might want to study.

## What you will build

We use silicon as our running example because it is small (2-atom primitive cell), uncontroversial (PBE works fine), and well-documented (the experimental band structure has been measured to high precision). The deliverable at the end of this chapter is the silicon band diagram, computed by you on your own machine, with convergence parameters that you have justified.

The path there is:

1. **Pseudopotentials and basis sets** ([§6.1](01-pseudo-basis.md)). Core electrons are inert; we replace them by smooth effective potentials. Valence electrons we expand in plane waves. Two numbers — the cutoff $E_\mathrm{cut}$ and the k-point grid — control the entire calculation.
2. **Your first calculation** ([§6.2](02-first-qe.md)). Install Quantum ESPRESSO, write a `pw.x` input file for silicon, run it, read the output. Then drive the same calculation from Python using ASE.
3. **Convergence testing** ([§6.3](03-convergence.md)). Sweep $E_\mathrm{cut}$ and the k-grid; plot total energy versus each; pick the smallest values that give 1 meV/atom convergence. This is the step that gets papers rejected when skipped.
4. **Band structures and DOS** ([§6.4](04-bands-dos.md)). The scf $\to$ nscf $\to$ bands pipeline. High-symmetry paths. Plotting in Python. Interpreting Si's indirect gap.
5. **Defects and formation energies** ([§6.5](05-defects.md)). A 64-atom supercell, one Si removed, the vacancy formation energy. The bridge between DFT and materials science.
6. **Exercises** ([§6.6](exercises.md)). Seven problems, solutions inline.

## What this chapter is not

It is not a comprehensive QE manual — that document is over 200 pages and you should keep it bookmarked. It is not a survey of DFT codes; VASP, ABINIT, CP2K, GPAW, SIESTA, FHI-aims and Castep would each deserve a chapter. It is not the place to debate functionals; we use PBE throughout because it is the lowest-cost defensible choice for solids, and we will revisit hybrids, meta-GGAs, and DFT+U in later chapters.

If you came here for ab-initio molecular dynamics — using DFT forces to propagate atoms through time — that is [Chapter 7](../ch07-md/index.md). Everything in this chapter is static: a fixed nuclear configuration, an electronic ground state.

## Practical prerequisites

You will need:

- A Linux or macOS machine (Windows users: WSL2 is fine). 8 GB RAM is enough for the silicon examples; 16 GB is comfortable.
- Quantum ESPRESSO 7.2 or later. Installation is in [§6.2](02-first-qe.md).
- Python 3.11 with `numpy`, `matplotlib`, and `ase>=3.23`.
- About 2 GB of disk space for pseudopotentials and scratch directories.
- Patience for SCF cycles. A converged silicon calculation takes 30 seconds on a laptop; a 64-atom defect supercell takes 10-30 minutes.

## A note on units

QE uses Rydberg atomic units internally: energies in Ry, lengths in Bohr radii ($a_0 = 0.529$ Å). Outputs include both Ry and eV. We will write energies in eV throughout the prose (because the literature does) but copy QE inputs verbatim with their native units. The conversion is $1\,\mathrm{Ry} = 13.6057\,\mathrm{eV}$.

## A note on reproducibility

Every command in this chapter has been run on a 2023 MacBook Pro (M2, 16 GB) under macOS, and verified on Ubuntu 22.04. Pseudopotential versions are pinned to the SSSP 1.3.0 PBE-efficiency release. Quantum ESPRESSO version is 7.3.1. Numerical values quoted to 1 meV will be reproducible; values quoted to 0.01 meV may not be, and depend on the BLAS library your QE was linked against. If you see deviations at the third decimal, do not panic.

Let us begin with the question that has to be answered before anything else: what, exactly, is the basis in which we represent the wavefunctions, and what potential do the valence electrons feel?
