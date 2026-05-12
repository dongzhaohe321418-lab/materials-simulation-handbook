# 4.7 Hartree–Fock, briefly

Within the Born–Oppenheimer approximation the central remaining task is to solve the electronic Schrödinger equation

$$\hat{H}_{\mathrm e}(\mathbf r; \mathbf R)\, \psi(\mathbf r; \mathbf R) = E(\mathbf R)\, \psi(\mathbf r; \mathbf R), \tag{4.7.1}$$

for a system of $N$ interacting electrons in the external potential of fixed nuclei. As we saw in §4.5, this problem is exponentially hard. Hartree–Fock (HF) is the simplest serious attempt to make it polynomial. The idea, conceptually, is breathtaking: assume the many-electron wavefunction is a *single* Slater determinant built from $N$ one-electron orbitals, then variationally choose those orbitals to minimise the energy. The result is a set of self-consistent one-electron equations of remarkable structure — they capture exchange exactly but neglect correlation entirely. HF is rarely used as a final method in modern materials science, but it is the conceptual scaffold on which density functional theory (Chapter 5) is built, and every electronic-structure code in the world traces some lineage to it.

This section sketches the HF construction. We will not derive the equations in full painful detail — that is a long calculation done correctly in any quantum chemistry textbook — but we will state every essential ingredient, identify the structure of the resulting equations, and pinpoint exactly what HF gets wrong.

## 4.7.1 The variational principle

The mathematical engine of HF (and DFT, and many other electronic-structure methods) is the **variational principle**: for any normalised trial wavefunction $\Psi$,

$$\langle \Psi | \hat{H} | \Psi \rangle \geq E_0, \tag{4.7.2}$$

with equality if and only if $\Psi$ is the exact ground state. We met this idea in Chapter 0.3 in the context of finding minimum-energy configurations; here it becomes the cornerstone of approximate quantum mechanics.

**Proof sketch.** Expand $|\Psi\rangle$ in the orthonormal eigenbasis $\{|\Phi_n\rangle\}$ of $\hat{H}$, $|\Psi\rangle = \sum_n c_n |\Phi_n\rangle$, with eigenvalues $E_0 \leq E_1 \leq E_2 \leq \ldots$. Normalisation gives $\sum_n |c_n|^2 = 1$. Then

$$\langle \Psi|\hat{H}|\Psi\rangle = \sum_n |c_n|^2 E_n \geq E_0 \sum_n |c_n|^2 = E_0. \quad\blacksquare$$

The strategy: choose a parameterised family $\Psi_\lambda$ of trial wavefunctions, compute $E(\lambda) = \langle\Psi_\lambda|\hat{H}|\Psi_\lambda\rangle$, and minimise over $\lambda$. The minimum is an upper bound on the true ground-state energy, and a *good* family produces a tight bound. The art is choosing a family that is rich enough to capture the relevant physics but simple enough to be computationally tractable.

## 4.7.2 The Hartree approximation

The first attempt — chronologically and pedagogically — is to assume the electrons are independent: write the wavefunction as a product of one-electron orbitals,

$$\Psi_{\mathrm H}(\mathbf r_1, \ldots, \mathbf r_N) = \phi_1(\mathbf r_1)\, \phi_2(\mathbf r_2)\,\cdots\,\phi_N(\mathbf r_N). \tag{4.7.3}$$

This is the **Hartree ansatz**. Each electron lives in its own orbital, and the orbitals are determined self-consistently by demanding that each $\phi_i$ feels the average ("mean-field") Coulomb repulsion from all the others.

Substituting (4.7.3) into the energy expectation value $\langle\hat{H}_{\mathrm e}\rangle$ and minimising with respect to each $\phi_i$ (subject to orthonormality) gives the **Hartree equations**:

$$\left[-\frac{\hbar^2}{2m_{\mathrm e}}\nabla^2 + v_{\mathrm{ext}}(\mathbf r) + v_{\mathrm H}^{(i)}(\mathbf r)\right] \phi_i(\mathbf r) = \varepsilon_i\, \phi_i(\mathbf r), \tag{4.7.4}$$

where $v_{\mathrm{ext}}$ is the electron–nuclear attraction and the **Hartree potential** for electron $i$ is

$$v_{\mathrm H}^{(i)}(\mathbf r) = \sum_{j\neq i} \int \frac{e^2 |\phi_j(\mathbf r')|^2}{|\mathbf r - \mathbf r'|}\, d^3 r'. \tag{4.7.5}$$

Each electron sees the classical electrostatic potential from the charge density $|\phi_j|^2$ of every other electron. The equations are coupled (the potential felt by $\phi_i$ depends on all the other $\phi_j$) and must be solved iteratively: guess the orbitals, build the potential, re-solve, repeat to convergence. This is the **self-consistent field** (SCF) procedure.

The Hartree picture has the right qualitative idea — replace the intractable two-body interaction with an averaged one-body field — but it has a fatal flaw: **it ignores Pauli antisymmetry**. The product (4.7.3) is symmetric under exchange of $\mathbf r_i$ and $\mathbf r_j$, not antisymmetric. The Hartree ansatz violates the Pauli principle, and as a consequence it does not enforce the exclusion principle: in (4.7.5), if electron $i$ happens to coincide with electron $j$ in the same spatial orbital, nothing stops them.

This is more than a formal complaint. The Hartree approximation systematically gets bond energies wrong by tens of eV, predicts the wrong atomic shell structure, and assigns molecules to the wrong ground-state spin. We need a properly antisymmetric ansatz.

## 4.7.3 The Hartree–Fock ansatz

The simplest properly antisymmetric trial wavefunction is a single **Slater determinant** of one-electron *spin-orbitals*. A spin-orbital $\chi_i(\mathbf x) = \phi_i(\mathbf r)\, \sigma_i(s)$ is a product of a spatial orbital and a spin function (up or down), and $\mathbf x = (\mathbf r, s)$ collects spatial and spin coordinates. The HF ansatz is

$$\Psi_{\mathrm{HF}}(\mathbf x_1, \ldots, \mathbf x_N) = \frac{1}{\sqrt{N!}}\det[\chi_i(\mathbf x_j)]. \tag{4.7.6}$$

The determinant changes sign under exchange of any two electrons (rows), satisfying (4.5.6); it vanishes if any two spin-orbitals are equal (columns), enforcing exclusion.

Now compute the energy expectation value $E_{\mathrm{HF}} = \langle\Psi_{\mathrm{HF}}|\hat{H}_{\mathrm e}|\Psi_{\mathrm{HF}}\rangle$. The calculation is tedious but elementary; the result is

$$E_{\mathrm{HF}} = \sum_i h_{ii} + \frac{1}{2}\sum_{ij}\bigl(J_{ij} - K_{ij}\bigr), \tag{4.7.7}$$

with the one-electron integrals

$$h_{ii} = \int \chi_i^*(\mathbf x)\!\left[-\tfrac{\hbar^2}{2m_{\mathrm e}}\nabla^2 + v_{\mathrm{ext}}\right]\!\chi_i(\mathbf x)\, d\mathbf x, \tag{4.7.8}$$

the **Coulomb integrals**

$$J_{ij} = \int\!\!\int \chi_i^*(\mathbf x_1)\chi_j^*(\mathbf x_2)\,\frac{e^2}{r_{12}}\,\chi_i(\mathbf x_1)\chi_j(\mathbf x_2)\,d\mathbf x_1 d\mathbf x_2, \tag{4.7.9}$$

and the **exchange integrals**

$$K_{ij} = \int\!\!\int \chi_i^*(\mathbf x_1)\chi_j^*(\mathbf x_2)\,\frac{e^2}{r_{12}}\,\chi_j(\mathbf x_1)\chi_i(\mathbf x_2)\,d\mathbf x_1 d\mathbf x_2. \tag{4.7.10}$$

The Coulomb term $J_{ij}$ is the classical electrostatic repulsion between the charge densities $|\chi_i|^2$ and $|\chi_j|^2$ — the Hartree potential of (4.7.5) is hidden inside it. The exchange term $K_{ij}$ is a *purely quantum* contribution with no classical analogue: it arises because the determinantal wavefunction correlates the positions of same-spin electrons, reducing their mutual repulsion. Exchange is zero between opposite-spin orbitals (the spin integration kills it), and it always *lowers* the energy.

!!! note "Self-interaction cancellation"
    Notice that $K_{ii} = J_{ii}$ exactly, so the sum $J_{ii} - K_{ii}$ vanishes. This means an electron does not interact with itself — the spurious self-interaction that lurks in (4.7.5) is cancelled exactly by the corresponding exchange term. This is one of the most beautiful features of Hartree–Fock, and one that approximate DFT functionals struggle to reproduce (the "self-interaction error", to which we will return in Chapter 5).

## 4.7.4 The HF equations

Minimise $E_{\mathrm{HF}}$ with respect to the spin-orbitals $\chi_i$, subject to orthonormality $\langle\chi_i|\chi_j\rangle = \delta_{ij}$. Use Lagrange multipliers $\varepsilon_{ij}$ to enforce the constraints, take the variation, and diagonalise the multiplier matrix. The result is the canonical form of the **Hartree–Fock equations**:

$$\boxed{\; \hat F\, \chi_i(\mathbf x) = \varepsilon_i\, \chi_i(\mathbf x), \;} \tag{4.7.11}$$

where the **Fock operator** is

$$\hat F = -\frac{\hbar^2}{2m_{\mathrm e}}\nabla^2 + v_{\mathrm{ext}}(\mathbf r) + \hat J(\mathbf r) - \hat K(\mathbf r). \tag{4.7.12}$$

The **Coulomb operator** $\hat J$ acts as a multiplicative potential,

$$\hat J(\mathbf r_1)\, \chi_i(\mathbf x_1) = \left[\sum_j \int \frac{e^2 |\chi_j(\mathbf x_2)|^2}{r_{12}}\, d\mathbf x_2\right] \chi_i(\mathbf x_1), \tag{4.7.13}$$

and the **exchange operator** $\hat K$ is non-local (depends on $\chi_i$ at $\mathbf r_2$, not $\mathbf r_1$):

$$\hat K(\mathbf r_1)\, \chi_i(\mathbf x_1) = \sum_j \chi_j(\mathbf x_1) \int \frac{e^2 \chi_j^*(\mathbf x_2) \chi_i(\mathbf x_2)}{r_{12}}\, d\mathbf x_2. \tag{4.7.14}$$

Equation (4.7.11) looks like a single-particle Schrödinger equation, but it is *nonlinear*: the operators $\hat J$ and $\hat K$ depend on the very orbitals $\chi_i$ we are trying to solve for. Like the Hartree equations, it must be solved iteratively — the **SCF cycle**.

**The SCF cycle.** A modern HF (or DFT) calculation looks like this:

1. Guess an initial set of orbitals $\{\chi_i^{(0)}\}$ (often from atomic orbitals on each nucleus or from a simpler theory).
2. Build the Fock operator $\hat F^{(n)}$ from the current orbitals.
3. Diagonalise $\hat F^{(n)}$ to obtain new orbitals $\{\chi_i^{(n+1)}\}$ and energies $\{\varepsilon_i^{(n+1)}\}$.
4. Test for convergence: compare $\chi^{(n+1)}$ with $\chi^{(n)}$ (or the total energy, or the density). If converged, stop. Otherwise go to step 2.

A typical molecular HF calculation reaches convergence in 10–50 iterations. The same loop, with $\hat F$ replaced by the Kohn–Sham Hamiltonian, drives every DFT code in Chapter 6.

## 4.7.5 What HF means and where it fails

Hartree–Fock has a clear physical interpretation:

- Each electron moves in the *average* electrostatic field created by all the others (the $\hat J$ term — same as Hartree).
- Additional reduction in energy comes from *Pauli exchange* between same-spin electrons (the $-\hat K$ term — unique to HF). Same-spin electrons avoid each other automatically because the Slater determinant correlates their positions; this is sometimes pictured as an "exchange hole" surrounding each electron.

What HF crucially *fails* to capture is **correlation between opposite-spin electrons**. Two electrons of opposite spin can, according to HF, sit on top of each other without paying any energy beyond the average Coulomb repulsion. In reality they avoid each other dynamically, and the energy cost of pretending otherwise — the **correlation energy** — is the gap between the HF energy and the exact non-relativistic ground-state energy:

$$E_{\mathrm{corr}} \equiv E_{\mathrm{exact}} - E_{\mathrm{HF}}. \tag{4.7.15}$$

Correlation is small in *absolute* terms — typically 1% of the total energy of a heavy atom, less for a small molecule — but it is enormous in *chemical* terms. Atomisation energies, reaction barriers, conformational preferences, hydrogen bonds, dispersion forces, magnetism: all are correlation-dominated phenomena. HF, by missing correlation, systematically underbinds molecules by tens of kcal/mol, overestimates bond lengths, and predicts entirely wrong reaction energetics.

The post-HF hierarchy of quantum chemistry — Møller–Plesset perturbation theory (MP2, MP4), configuration interaction (CISD, CCSD, CCSD(T)), multi-reference methods (CASSCF, CASPT2) — exists precisely to recover correlation systematically on top of an HF reference. The cost scales steeply, however: $\mathcal O(N^5)$ for MP2, $\mathcal O(N^7)$ for CCSD(T). HF itself nominally scales as $\mathcal O(N^4)$ from the four-index integrals (4.7.9)–(4.7.10), though linear-scaling algorithms exist for large systems.

## 4.7.6 Why this matters for the rest of the book

We will not actually *run* Hartree–Fock calculations in this book — DFT (Chapter 5) has displaced HF as the workhorse of materials simulation, for good reason. But the structural lessons of HF are essential.

1. **The single Slater determinant is a powerful organising idea.** DFT, in the Kohn–Sham formulation, also represents the ground state as a single Slater determinant of *fictitious* one-electron orbitals — though the orbitals are interpreted differently.

2. **The SCF cycle is universal.** Every electronic-structure code you will encounter — from a hand-written `pyscf` script to a national-laboratory plane-wave code — runs an SCF loop with essentially the same logic as §4.7.4.

3. **The Fock operator splits naturally into kinetic, external, Hartree (classical Coulomb) and exchange parts.** This same partitioning structures the Kohn–Sham Hamiltonian, with the addition of a *correlation* term:
$$\hat{H}_{\mathrm{KS}} = -\tfrac12 \nabla^2 + v_{\mathrm{ext}} + v_{\mathrm H} + v_{\mathrm{xc}}.$$
The first three pieces are the same as in HF (the Hartree part now including self-interaction, since we represent the system by a density, not orbitals); the last is the **exchange–correlation potential**, where the magic of DFT lives.

4. **Exchange is essentially free; correlation is hard.** HF treats exchange exactly. DFT, in approximate functionals, treats both exchange and correlation approximately — but does so in a way that captures most of the correlation as well, at HF cost.

5. **Self-interaction is a hidden danger.** HF cancels it exactly by construction. Approximate DFT functionals do not, and this is one of the major systematic errors of practical DFT — visible as too-narrow band gaps in semiconductors, over-delocalisation of charged defects, and spurious fractional charges in stretched bonds.

## 4.7.7 A bridge to DFT

Here is the conceptual leap that makes DFT — and the rest of this book — possible. Hohenberg and Kohn proved in 1964 that the ground-state electron density $n(\mathbf r) = \sum_i |\chi_i(\mathbf r)|^2$ contains all the information of the wavefunction. The exact ground-state energy is a functional of $n$ alone, $E[n]$, even though we do not know its form. Kohn and Sham proposed in 1965 to write

$$E[n] = T_{\mathrm s}[n] + \int v_{\mathrm{ext}}(\mathbf r) n(\mathbf r) d^3 r + E_{\mathrm H}[n] + E_{\mathrm{xc}}[n], \tag{4.7.16}$$

where $T_{\mathrm s}$ is the kinetic energy of a fictitious non-interacting system with the same density, $E_{\mathrm H}$ is the classical Hartree energy, and $E_{\mathrm{xc}}$ — the **exchange–correlation functional** — absorbs everything else. Minimising (4.7.16) leads to a set of one-electron equations indistinguishable in structure from HF (4.7.11), but with the non-local exchange operator $\hat K$ replaced by the local exchange–correlation potential $v_{\mathrm{xc}}(\mathbf r) = \delta E_{\mathrm{xc}}/\delta n(\mathbf r)$.

That single replacement — *non-local exchange becomes a local exchange–correlation potential of the density* — converts HF into Kohn–Sham DFT. The accuracy of the resulting theory is entirely controlled by the quality of $E_{\mathrm{xc}}[n]$, which must be modelled. Chapter 5 is the story of those models.

For now: you have all the conceptual scaffolding you need to read the rest of the book. The many-electron Schrödinger equation is exact but intractable; the Born–Oppenheimer approximation cleaves nuclei from electrons; Hartree–Fock gives a single-determinant variational ansatz with exact exchange but no correlation; and DFT, the workhorse of Chapter 5 onwards, repackages the same one-electron structure with all the missing physics bundled into a universal density functional. Turn the page.
