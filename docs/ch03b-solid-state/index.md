# Chapter 3b — Solid State Physics Prerequisites

[Open in Jupyter (browser)](/materials-simulation-handbook/lite/lab/index.html?path=ch03b-solid-state.ipynb){ .md-button .md-button--primary }

> *"A solid is the densest, most boring, and most useful state of matter."* — anonymous condensed matter lecturer

Chapter 3 took you from isolated atoms to bonded crystals: you learnt what a unit cell is, drew a few Bravais lattices, met reciprocal space, and computed a structure factor. Chapter 4 then taught you quantum mechanics in finite systems — a particle in a box, a harmonic oscillator, the many-electron problem and the exponential wall. What is glaringly absent is the bridge between the two: *how does quantum mechanics behave in an infinite periodic potential?* That bridge is solid state physics, and without it nothing in the rest of the book makes sense. You cannot read a DFT output, configure a `k`-point grid, interpret a band gap, fit a machine-learning potential to a phonon spectrum, or design a graph neural network that respects translational symmetry, unless you know the four or five ideas that this chapter installs.

A note on reading order. This chapter is labelled `3b` because conceptually it expands Chapter 3 (descriptive crystallography) into its quantum counterpart, but it does use a small amount of Chapter 4 material — chiefly the time-independent Schrödinger equation, eigenvalues of Hermitian operators, and bra–ket notation. If you are following the linear path (Path A of the [learning path](../learning-path.md)), the recommended order is **Ch 3 → Ch 4 → Ch 3.5 → Ch 5**: this is the order that uses the least forward-referenced material. If you prefer the original "geometry first" ordering Ch 3 → Ch 3.5 → Ch 4, the bits of Ch 4 we lean on are brief enough that you can skim the relevant single-particle Schrödinger sections (4.2 and 4.3) and come back; nothing in this chapter requires the many-electron material from Chapter 4.5 onwards.

## Why a working materials simulator needs this chapter

Every method in Tier 1 of this book — DFT (Ch 5, Ch 6), molecular dynamics (Ch 7), MLIPs (Ch 9), graph neural networks (Ch 10) — is shaped by *one* mathematical fact: the potential a particle sees inside a crystal is *periodic*. From that single observation flow:

- **Bloch's theorem**, which says you only ever need to solve the Schrödinger equation inside one unit cell, with a `k`-dependent boundary condition. This is *exactly* the structure every plane-wave DFT code (VASP, Quantum ESPRESSO, ABINIT, CP2K) exploits — and the reason `k`-point convergence is the single most common parameter you will tune in Chapter 6.

- **Band structure**, which tells you whether a material is a metal, a semiconductor, or an insulator. We will derive band structures in two complementary limits — nearly free electrons (weak potential) and tight binding (strong, localised potential). Both limits appear in real codes: plane-wave DFT is the nearly-free-electron extreme; localised-orbital DFT (SIESTA, FHI-aims, NWChem-LCAO) is the tight-binding extreme. MLIPs that learn from local environments — SchNet, NequIP, MACE — are effectively learning an implicit tight-binding model, as we shall see in Chapter 9.

- **Phonons**, the quantised lattice vibrations. In Chapter 7 you will run molecular dynamics: atoms wobble around their equilibrium positions in a thermostatted simulation. The Fourier transform of that wobble *is* a phonon spectrum. In Chapter 8 the vibrational free energy — built from phonon densities of states — is what makes a phase stable at finite temperature. And in Chapter 9 you will benchmark MLIPs against DFT phonon dispersions, because matching phonons at the $\Gamma$ point is the strictest test of an interatomic potential.

- **The free-electron and Sommerfeld models**, which are conceptually the simplest band structure of all: $E(\mathbf k) = \hbar^2 k^2/2m$, no lattice at all. This is the "jellium" that the local density approximation in DFT (Ch 5) takes as its reference. Every time you write `xc=LDA` in an input file you are betting that the local electron density looks, to leading order, like a uniform gas. Knowing why that bet works — and when it fails — is impossible without the Sommerfeld machinery.

- **Defects and doping**, which are how real materials acquire their useful properties. In Chapter 6 you will compute defect formation energies; in the capstone project you will screen dopants for a target band gap. Both rest on the hydrogenic donor/acceptor picture developed here.

## Roadmap

This chapter has seven content sections and an exercise set.

1. **Bloch's theorem** (§3b.1). The cornerstone. We prove it from scratch starting from translational symmetry, derive the crystal momentum, define the first Brillouin zone, and introduce the band index. About 4000 words.

2. **Nearly free electrons** (§3b.2). Weak periodic potential, empty lattice, degenerate perturbation theory at zone boundaries, gap opening, physical picture in terms of standing waves.

3. **Tight binding** (§3b.3). The opposite limit. LCAO ansatz, 1D monatomic chain ($E = -2t\cos ka$), then graphene with its Dirac cones. Working Python code that diagonalises the graphene Hamiltonian along K–$\Gamma$–M–K.

4. **The free electron gas and Sommerfeld** (§3b.4). 3D density of states, Fermi energy, finite-temperature corrections, electronic specific heat. Numerical evaluation for copper.

5. **Phonons** (§3b.5). Classical lattice dynamics: 1D monatomic chain ($\omega = 2\sqrt{K/m}|\sin(ka/2)|$), then 1D diatomic chain with its acoustic and optical branches, then the 3D dynamical matrix.

6. **Specific heat: Einstein and Debye** (§3b.6). Why Dulong–Petit fails, why Einstein's tweak is not enough, why Debye's $T^3$ law works at low temperature. The Debye temperature is the single most useful number you will compute for a new material.

7. **Defects and band engineering** (§3b.7). Point defects, shallow donors, the effective mass approximation, alloying and strain. The shortest section, but a forward reference to almost every applied chapter.

8. **Exercises** (§3b.8). Eight problems, with solutions in admonitions, half analytical and half numerical.

## What you need

From Chapter 0: linear algebra (eigenvalues of $2\times2$ Hermitian matrices), Fourier series and Fourier transforms, the gradient and Laplacian, partial derivatives, complex exponentials. From Chapter 3: Bravais lattices, the reciprocal lattice, the first Brillouin zone, the structure factor. From Chapter 4: the time-independent Schrödinger equation, eigenvalues of Hermitian operators, bra-ket notation, the harmonic oscillator, and the fact that commuting Hermitian operators share an eigenbasis.

## What you do *not* need

We will not quantise the phonon field in second-quantised language. We will not develop the full theory of the magnetic Brillouin zone, time-reversal symmetry, or topological invariants. We will not do anything with phonon–phonon scattering, polarons, superconductivity, or strong correlations. Those belong in a serious solid-state course; we want only the minimum that makes the rest of the handbook readable. Five pages of Ashcroft and Mermin, faithfully.

When you finish this chapter, you will be able to read a band-structure plot, recognise an acoustic branch, write a tight-binding Hamiltonian, and explain — at three in the morning, to a sceptical reviewer — why your DFT calculation needed an $8\times 8\times 8$ Monkhorst–Pack mesh. Onwards.
