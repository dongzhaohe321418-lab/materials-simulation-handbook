# Chapter 5 — Density Functional Theory

[Open in Jupyter (browser)](/materials-simulation-handbook/lite/lab/index.html?path=ch05-dft.ipynb){ .md-button .md-button--primary }

In Chapter 4 we faced an uncomfortable truth. The non-relativistic, Born–Oppenheimer many-electron Hamiltonian

$$
\hat{H} = -\frac{\hbar^2}{2m}\sum_{i=1}^{N}\nabla_i^{2}
\;+\; \sum_{i=1}^{N} v_\mathrm{ext}(\mathbf r_i)
\;+\; \frac{1}{2}\sum_{i\neq j}\frac{e^{2}}{4\pi\varepsilon_0\,|\mathbf r_i-\mathbf r_j|}
$$

is, in principle, the answer to almost every question we want to ask about a material at zero temperature. In practice, the wavefunction $\Psi(\mathbf r_1,\sigma_1,\dots,\mathbf r_N,\sigma_N)$ lives in a $3N$-dimensional configuration space. Store it on a coarse $10^3$ grid for one electron and the cost is $10^3$ floats; for $N$ electrons, $10^{3N}$. For a single atom of iron — $N=26$ — that is $10^{78}$ numbers, more than the atoms in the visible universe. Hartree–Fock cuts this catastrophic scaling by restricting $\Psi$ to a single Slater determinant, but it pays for that simplicity in lost correlation energy that is qualitatively important for almost every chemical bond.

!!! note "In plain language"
    The exact answer to "how does this material behave?" is a single huge object — the many-electron wavefunction — and it is far too big to ever write down or store, even for one iron atom. Density functional theory makes a bold bet: that a much smaller and more familiar object, the electron density $n(\mathbf r)$ (just "how much electron charge sits at each point in space"), already contains everything we actually need to know. This chapter builds that idea up carefully, from why such a bet could possibly be allowed to how we cash it out in a calculation we can run.

This chapter develops the framework that, for the past forty years, has dominated electronic structure calculations in materials science and chemistry: **density functional theory** (DFT). The central idea is breathtakingly economical. Rather than the $3N$-dimensional wavefunction, take the **electron density**

$$
n(\mathbf r) = N\int|\Psi(\mathbf r,\mathbf r_2,\dots,\mathbf r_N)|^{2}\,\mathrm d\mathbf r_2\cdots\mathrm d\mathbf r_N,
$$

a single three-dimensional scalar field, as the fundamental variable. Hohenberg and Kohn proved in 1964 that this is enough: the ground-state density determines, in principle, every property of the system. Kohn and Sham then provided, in 1965, a practical scheme that turns the interacting many-body problem into a set of self-consistent single-particle equations of the same formal cost as Hartree–Fock — but with correlation included, at least in principle.

What follows is a careful, derivation-first development.

!!! tip "How to read this chapter (undergraduate)"
    This chapter is dense, but it rewards a layered reading — match it to the three layers from the [undergraduate guide](../undergraduate/learning-paths.md):

    - **Layer 1 — get the intuition.** Aim only to come away with three pictures: *why* the density is enough to work with, *what* the Kohn–Sham trick does (replace the hard interacting problem with an easier non-interacting one that has the same density), and *what* the self-consistent field (SCF) loop is (guess a density, build a potential, solve, get a new density, repeat until it stops changing). Sections 5.1, 5.3, and 5.5 carry these ideas.
    - **Layer 2 — the derivations.** The full Hohenberg–Kohn proofs (Section 5.2) and the uniform-electron-gas derivation of LDA exchange (Section 5.4) are where the rigour lives. It is completely fine to skim these the first time and return once the intuition is solid.
    - **Layer 3 — the judgement.** Choosing an exchange–correlation functional (Section 5.4) and knowing where DFT fails (Section 5.6) is the practitioner's craft. This matters most when you start running your own calculations in Chapter 6.

    When the symbols pile up, lean on the [beginner glossary](../undergraduate/glossary-for-beginners.md) for unfamiliar terms and the [formula-reading guide](../undergraduate/formula-reading-guide.md) for decoding equations piece by piece.

## What you will find in this chapter

- **§5.1 The Thomas–Fermi idea.** The first attempt, in 1927, to do electronic structure with the density alone. We derive the Thomas–Fermi kinetic energy functional, see why the resulting energy is variational, and understand precisely why it fails: it misses shell structure, and Teller proved that it predicts no molecule to be stable. The failure is instructive — it sharpens the question of what a *correct* density functional would need to capture.
- **§5.2 The Hohenberg–Kohn theorems.** Full proofs, every step shown, of the two theorems that put DFT on a rigorous footing. Theorem I: the ground-state density determines the external potential up to an additive constant. Theorem II: an energy functional $E[n]$ exists and is minimised by the true ground-state density. We define the universal functional $F[n]$ and confess that we cannot write it down.
- **§5.3 The Kohn–Sham construction.** The crucial reformulation that made DFT computable. We map the interacting system onto an auxiliary non-interacting system with the same density, derive the Kohn–Sham equations, and are careful about what the KS orbitals and eigenvalues do — and do not — mean physically.
- **§5.4 Exchange–correlation functionals.** Everything we do not know is hidden inside $E_{xc}[n]$. Perdew's Jacob's ladder organises the approximations: LDA, GGA, meta-GGA, hybrids, and beyond. We derive the LDA exchange energy from the uniform electron gas, explain when GGAs help, where hybrids buy accuracy with cost, and how van der Waals corrections plug a real hole. The section closes with a decision table for choosing a functional.
- **§5.5 The self-consistent field.** The KS equations are nonlinear because the potential depends on the density that the equations themselves produce. We work through the SCF loop, see why naive iteration oscillates, and implement Pulay/DIIS mixing in a complete, runnable Python program that solves a 1D model hydrogen chain end-to-end.
- **§5.6 Where DFT fails.** An honest tour of band gap underestimation, the dispersion problem, strong correlation, self-interaction error, and the inaccessibility of excited states from a ground-state theory. We point to GW, BSE, DMFT, and coupled cluster as the methods to reach for when DFT will not do.
- **§5.7 Exercises.** Eight problems, with worked solutions, ranging from re-deriving Hohenberg–Kohn in your own words to extending the SCF code.

## Why DFT became the workhorse

A typical 100-atom DFT calculation today runs on a laptop in minutes. The same calculation at coupled cluster CCSD(T) accuracy would, if it were even tractable, take weeks on a cluster. DFT sits at a sweet spot: roughly the accuracy of post-Hartree–Fock for a fraction of a percent of the cost. That is why almost every database of computed material properties — the Materials Project, OQMD, AFLOW, NOMAD — is built on DFT, and why almost every machine-learning interatomic potential (Chapter 9) is trained on DFT data. The remainder of this book takes that for granted; this chapter explains why we are entitled to.

In Chapter 6 we move from the theory to the practice of *running* DFT calculations — plane waves and pseudopotentials, $k$-point sampling, convergence testing, and the choice of code. Read this chapter first.
