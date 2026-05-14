# 4.7 Hartree–Fock, briefly

Within the Born–Oppenheimer approximation the central remaining task is to solve the electronic Schrödinger equation

$$\hat{H}_{\mathrm e}(\mathbf r; \mathbf R)\, \psi(\mathbf r; \mathbf R) = E(\mathbf R)\, \psi(\mathbf r; \mathbf R), \tag{4.7.1}$$

for a system of $N$ interacting electrons in the external potential of fixed nuclei. As we saw in §4.5, this problem is exponentially hard. Hartree–Fock (HF) is the simplest serious attempt to make it polynomial. The idea, conceptually, is breathtaking: assume the many-electron wavefunction is a *single* Slater determinant built from $N$ one-electron orbitals, then variationally choose those orbitals to minimise the energy. The result is a set of self-consistent one-electron equations of remarkable structure — they capture exchange exactly but neglect correlation entirely. HF is rarely used as a final method in modern materials science, but it is the conceptual scaffold on which density functional theory (Chapter 5) is built, and every electronic-structure code in the world traces some lineage to it.

This section sketches the HF construction. We will not derive the equations in full painful detail — that is a long calculation done correctly in any quantum chemistry textbook — but we will state every essential ingredient, identify the structure of the resulting equations, and pinpoint exactly what HF gets wrong.

## 4.7.0 The plan

The strategy of Hartree–Fock has three logical steps. Knowing them in advance makes the algebra below much easier to navigate.

1. **Choose a trial form.** We restrict the $N$-electron wavefunction to a single Slater determinant of $N$ unknown one-electron spin-orbitals $\chi_i$. This builds in antisymmetry but throws away all interelectronic correlation beyond what one determinant can express.
2. **Compute the energy.** Plug the determinant into $\langle\hat H_{\mathrm e}\rangle$; the result is a sum of one-electron integrals and two-electron Coulomb–exchange integrals over the orbitals.
3. **Minimise.** Vary the orbitals $\chi_i$ to make the energy stationary, subject to orthonormality. The Euler–Lagrange equations are the Hartree–Fock equations, a nonlinear eigenvalue problem for the orbitals.

The pay-off, identified by Hartree and Fock in the late 1920s and made systematic by Roothaan and Hall in the 1950s, is that the resulting one-electron equations look almost exactly like the single-particle Schrödinger equation we have already learned to solve numerically. The cost has shrunk from exponential in $N$ to polynomial — at the price of throwing away the correlation between opposite-spin electrons.

## 4.7.1 The variational principle

The mathematical engine of HF (and DFT, and many other electronic-structure methods) is the **variational principle**: for any normalised trial wavefunction $\Psi$,

$$\langle \Psi | \hat{H} | \Psi \rangle \geq E_0, \tag{4.7.2}$$

with equality if and only if $\Psi$ is the exact ground state. We met this idea in Chapter 0.3 in the context of finding minimum-energy configurations; here it becomes the cornerstone of approximate quantum mechanics.

**Proof sketch.** Expand $|\Psi\rangle$ in the orthonormal eigenbasis $\{|\Phi_n\rangle\}$ of $\hat{H}$, $|\Psi\rangle = \sum_n c_n |\Phi_n\rangle$, with eigenvalues $E_0 \leq E_1 \leq E_2 \leq \ldots$. Normalisation gives $\sum_n |c_n|^2 = 1$. Then

$$\langle \Psi|\hat{H}|\Psi\rangle = \sum_n |c_n|^2 E_n \geq E_0 \sum_n |c_n|^2 = E_0. \quad\blacksquare$$

The strategy: choose a parameterised family $\Psi_\lambda$ of trial wavefunctions, compute $E(\lambda) = \langle\Psi_\lambda|\hat{H}|\Psi_\lambda\rangle$, and minimise over $\lambda$. The minimum is an upper bound on the true ground-state energy, and a *good* family produces a tight bound. The art is choosing a family that is rich enough to capture the relevant physics but simple enough to be computationally tractable.

!!! tip "The variational principle is the engine of approximate quantum mechanics"
    Almost every method in this book — HF, DFT, CI, CC, variational Monte Carlo, the neural network ansätze of FermiNet and PauliNet — is variational at heart. You pick an ansatz with parameters, compute the energy as a function of those parameters, and minimise. The difference between methods is the choice of ansatz. The full development of the variational principle, including the Hylleraas–Undheim theorem (that excited-state energies are also bounded from above by the corresponding eigenvalues of the projected Hamiltonian), is treated in Chapter 0.3.

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

### Deriving the Hartree equations: a worked walkthrough

To see the mean-field idea at work without yet worrying about antisymmetry, plug (4.7.3) into the energy expectation value:

$$E_{\mathrm H} = \langle\Psi_{\mathrm H}|\hat H_{\mathrm e}|\Psi_{\mathrm H}\rangle = \sum_i \int \phi_i^*(\mathbf r)\,\hat h(\mathbf r)\,\phi_i(\mathbf r)\,d\mathbf r + \frac{1}{2}\sum_{i \neq j}\int\!\!\int \frac{e^2|\phi_i(\mathbf r_1)|^2 |\phi_j(\mathbf r_2)|^2}{|\mathbf r_1 - \mathbf r_2|}\,d\mathbf r_1 d\mathbf r_2,$$

where $\hat h = -\hbar^2\nabla^2/(2m_{\mathrm e}) + v_{\mathrm{ext}}$. The first term is one-electron, summing kinetic and external-potential expectation values across the $N$ orbitals; the second is two-electron, the Coulomb energy between the *charge densities* $|\phi_i|^2$ and $|\phi_j|^2$.

!!! note "Why this step? — the mean field appears here"
    The pairwise Coulomb interaction $e^2/|\mathbf r_1 - \mathbf r_2|$ in the original Hamiltonian (4.5.3) couples the coordinates of two specific electrons. In a product wavefunction the integration *factorises*: the two-electron integral becomes a product of one-electron integrals, with each electron seeing only the *average* charge density of the others. This is the mean-field reduction — the two-body operator collapses into a sum of one-body operators, each depending on the orbitals through the density.

Minimising $E_{\mathrm H}$ over the orbitals subject to $\langle\phi_i|\phi_j\rangle = \delta_{ij}$ (via Lagrange multipliers) gives the Hartree equations (4.7.4). Each $\phi_i$ obeys a one-electron Schrödinger equation, but the potential it feels depends on all the other orbitals — hence the SCF iteration.

The Hartree approximation has historical importance (Hartree introduced it in 1928 to compute atomic structure) and is still occasionally useful as a starting point. But it is fundamentally a *bosonic* ansatz applied to fermions, and the resulting predictions are not even qualitatively right beyond the simplest atoms.

## 4.7.3 The Hartree–Fock ansatz

The simplest properly antisymmetric trial wavefunction is a single **Slater determinant** of one-electron *spin-orbitals*. A spin-orbital $\chi_i(\mathbf x) = \phi_i(\mathbf r)\, \sigma_i(s)$ is a product of a spatial orbital and a spin function (up or down), and $\mathbf x = (\mathbf r, s)$ collects spatial and spin coordinates. The HF ansatz is

$$\Psi_{\mathrm{HF}}(\mathbf x_1, \ldots, \mathbf x_N) = \frac{1}{\sqrt{N!}}\det[\chi_i(\mathbf x_j)]. \tag{4.7.6}$$

The determinant changes sign under exchange of any two electrons (rows), satisfying (4.5.6); it vanishes if any two spin-orbitals are equal (columns), enforcing exclusion.

### Slater determinants written out

For two electrons in spin-orbitals $\chi_a, \chi_b$, the Slater determinant is

$$\Psi_{\mathrm{HF}}(\mathbf x_1, \mathbf x_2) = \frac{1}{\sqrt 2}\begin{vmatrix}\chi_a(\mathbf x_1) & \chi_b(\mathbf x_1)\\ \chi_a(\mathbf x_2) & \chi_b(\mathbf x_2)\end{vmatrix} = \frac{1}{\sqrt 2}\bigl[\chi_a(\mathbf x_1)\chi_b(\mathbf x_2) - \chi_b(\mathbf x_1)\chi_a(\mathbf x_2)\bigr],$$

which is exactly the antisymmetric two-electron combination we met in §4.5.3. For three electrons,

$$\Psi_{\mathrm{HF}}(\mathbf x_1, \mathbf x_2, \mathbf x_3) = \frac{1}{\sqrt 6}\begin{vmatrix}\chi_a(\mathbf x_1) & \chi_b(\mathbf x_1) & \chi_c(\mathbf x_1)\\ \chi_a(\mathbf x_2) & \chi_b(\mathbf x_2) & \chi_c(\mathbf x_2)\\ \chi_a(\mathbf x_3) & \chi_b(\mathbf x_3) & \chi_c(\mathbf x_3)\end{vmatrix},$$

which when expanded by cofactors becomes a sum of six terms (one for each permutation of the three labels), with alternating signs to ensure antisymmetry under every pairwise exchange.

In general, the determinant of an $N\times N$ matrix is a sum of $N!$ signed products: each is a particular permutation of the orbital labels among the electrons, with sign $(-1)^P$ where $P$ is the parity of the permutation. This is what makes Slater determinants intrinsically antisymmetric — every term is generated from every other by an exchange, with a controlled sign.

!!! note "Computational implication"
    For $N$ electrons, computing the value of a Slater determinant at a single point in configuration space naively requires evaluating $N!$ terms. Numerical methods (LU decomposition) compute the determinant in $\mathcal O(N^3)$ time, which is one reason HF and Kohn–Sham DFT scale polynomially with system size rather than factorially. The $N!$ count is what kills FCI; the $\mathcal O(N^3)$ count is what makes mean-field theory practical.

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

### Koopmans' theorem

The eigenvalues $\varepsilon_i$ of the Fock operator are not arbitrary Lagrange multipliers — they have a direct physical interpretation. **Koopmans' theorem** (1934) states that, within HF and assuming the orbitals do not relax upon ionisation,

$$\boxed{\;-\varepsilon_i = I_i \quad \text{(ionisation energy from orbital $i$)},\;} \tag{4.7.K}$$

and conversely $-\varepsilon_a = A_a$ for an unoccupied (virtual) orbital, where $A$ is the electron affinity.

**Derivation sketch.** Compute the energy difference $E(\Psi_N) - E(\Psi_{N-1})$ where $\Psi_{N-1}$ is the Slater determinant obtained by removing the electron from orbital $i$. Using the energy expression (4.7.7) and the fact that the $N-1$ orbitals in $\Psi_{N-1}$ are identical to the corresponding ones in $\Psi_N$ (the "frozen-orbital" assumption),

$$E(\Psi_N) - E(\Psi_{N-1}) = h_{ii} + \sum_{j \neq i}(J_{ij} - K_{ij}) = \varepsilon_i,$$

where the last equality follows by inspecting the HF equation $\hat F \chi_i = \varepsilon_i \chi_i$ — the eigenvalue $\varepsilon_i$ is precisely the diagonal matrix element of the Fock operator, which equals $h_{ii}$ plus the Coulomb and exchange contributions from all *other* electrons. Hence the energy lost in removing the electron is $\varepsilon_i$, and the ionisation energy is $-\varepsilon_i > 0$ (since $\varepsilon_i < 0$ for bound electrons).

!!! note "Why this step?"
    The crucial input is the *frozen orbital* assumption: the remaining $N-1$ electrons are not allowed to relax in response to the missing electron. In reality they *do* relax — the remaining electrons collapse inward toward the nucleus once the screening from electron $i$ is removed — and the true ionisation energy is slightly less than $-\varepsilon_i$ (by the "orbital relaxation energy"). Koopmans is therefore an *approximation*, but a remarkably good one for valence ionisations (accurate to within $\sim 0.5$ eV for many molecules). The errors partly cancel for HF: orbital relaxation lowers the ionisation energy, but correlation typically raises it, and the two cancel by symmetry. Koopmans fails badly for core ionisation, where relaxation is enormous.

Koopmans' theorem gives HF a direct interpretation in photoelectron spectroscopy: peaks in the UPS/XPS spectrum correspond to $-\varepsilon_i$ values, suitably labelled by the orbital character. It is one of the cleanest links between calculation and experiment in the entire theory.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. Why does the Hartree product wavefunction violate the Pauli principle, and how does a Slater determinant fix this?
    2. What is the physical difference between the Coulomb integral $J_{ij}$ and the exchange integral $K_{ij}$, and why does $J_{ii} - K_{ii} = 0$ matter?
    3. The Hartree–Fock equations look like single-particle Schrödinger equations but must be solved iteratively — why, and what is this iteration called?

    If any of these is shaky, re-read the preceding section before continuing.

## 4.7.5 What HF means and where it fails

Hartree–Fock has a clear physical interpretation:

- Each electron moves in the *average* electrostatic field created by all the others (the $\hat J$ term — same as Hartree).
- Additional reduction in energy comes from *Pauli exchange* between same-spin electrons (the $-\hat K$ term — unique to HF). Same-spin electrons avoid each other automatically because the Slater determinant correlates their positions; this is sometimes pictured as an "exchange hole" surrounding each electron.

What HF crucially *fails* to capture is **correlation between opposite-spin electrons**. Two electrons of opposite spin can, according to HF, sit on top of each other without paying any energy beyond the average Coulomb repulsion. In reality they avoid each other dynamically, and the energy cost of pretending otherwise — the **correlation energy** — is the gap between the HF energy and the exact non-relativistic ground-state energy:

$$E_{\mathrm{corr}} \equiv E_{\mathrm{exact}} - E_{\mathrm{HF}}. \tag{4.7.15}$$

Correlation is small in *absolute* terms — typically 1% of the total energy of a heavy atom, less for a small molecule — but it is enormous in *chemical* terms. Atomisation energies, reaction barriers, conformational preferences, hydrogen bonds, dispersion forces, magnetism: all are correlation-dominated phenomena. HF, by missing correlation, systematically underbinds molecules by tens of kcal/mol, overestimates bond lengths, and predicts entirely wrong reaction energetics.

### The size of correlation energy: a concrete table

For first-row atoms the correlation energy is roughly $-1$ eV per pair of opposite-spin electrons. Some representative numbers (from coupled-cluster calculations):

| Species | $E_{\mathrm{HF}}$ (Ha) | $E_{\mathrm{exact}}$ (Ha) | $E_{\mathrm{corr}}$ (eV) |
|---|---:|---:|---:|
| He | $-2.8617$ | $-2.9037$ | $-1.144$ |
| Li | $-7.4327$ | $-7.4781$ | $-1.234$ |
| Be | $-14.5730$ | $-14.6674$ | $-2.568$ |
| C | $-37.6886$ | $-37.8450$ | $-4.254$ |
| Ne | $-128.5471$ | $-128.9376$ | $-10.626$ |
| H$_2$ at $R_e$ | $-1.1336$ | $-1.1745$ | $-1.114$ |
| H$_2$O | $-76.0671$ | $-76.4380$ | $-10.095$ |

In atomic units (Hartree), a correlation energy of $-0.4$ Ha is about $-10.9$ eV. Compare to typical chemical energies: atomisation energy of H$_2$ is 4.75 eV; reaction barriers are 0.1–2 eV; hydrogen-bond energies are 0.1–0.3 eV. The correlation energy is *not small* by any chemical standard, and recovering it accurately is the central business of post-HF and DFT methods.

!!! warning "What HF gets very wrong: dissociation"
    A striking failure of HF is bond dissociation. For H$_2$ stretched to large $R$, the single-determinant restricted HF wavefunction does not dissociate to two neutral H atoms; instead, it dissociates to an unphysical 50/50 mixture of H$^+$+H$^-$ and 2H. The error in the energy at $R = \infty$ is several eV — entirely due to correlation. Unrestricted HF (different orbitals for different spins) cures this in part but breaks spin symmetry. The correct treatment requires a *multi-reference* wavefunction with two determinants, which is the entry point to CASSCF and beyond.

### What correlation energy *is*, physically

It is worth being precise about what correlation means in this context. Two complementary pictures:

- **Dynamic correlation** is the short-range avoidance of electrons due to Coulomb repulsion. Even in the helium ground state, where both electrons are in $1s$, they are not independently distributed — finding electron 1 at $\mathbf r$ depresses the probability of finding electron 2 nearby. The Hartree–Fock determinant misses this entirely (between opposite-spin electrons), since the spatial part of the singlet wavefunction is $\phi(\mathbf r_1)\phi(\mathbf r_2)$ — uncorrelated. Capturing dynamic correlation requires admixing excited determinants (MP2, CCSD, etc.).
- **Static (non-dynamic) correlation** arises when more than one Slater determinant has comparable weight in the true wavefunction — for example, stretched bonds where bonding and antibonding configurations are nearly degenerate. HF assumes a single dominant determinant and fails when this is not true. Multi-reference methods (CASSCF) address this directly.

For ground-state equilibrium geometries of well-behaved molecules, dynamic correlation dominates. For transition states, photochemistry, transition-metal complexes and stretched bonds, static correlation can be larger.

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

## 4.7.6a Restricted vs unrestricted HF

A practical subtlety worth mentioning: there are two flavours of HF in common use.

**Restricted HF (RHF).** Each spatial orbital is occupied by *two* electrons, one spin-up and one spin-down. This is the natural choice for closed-shell molecules (even number of electrons, total spin $S = 0$) and produces a pure spin singlet by construction. The orbital count is halved (we solve for $N/2$ spatial orbitals rather than $N$ spin-orbitals), and computational cost roughly halves.

**Unrestricted HF (UHF).** Spin-up and spin-down electrons are described by *different* spatial orbitals. This is necessary for open-shell systems (radicals, transition metals, dissociating bonds) where the two spin channels feel different effective potentials. The trade-off is that the resulting determinant is not a pure spin state — it is a mixture of $S, S+1, S+2, \ldots$ states, and the "spin contamination" $\langle\hat S^2\rangle - S(S+1)$ is a diagnostic that should be reported in any UHF calculation.

A third option, **restricted open-shell HF (ROHF)**, forces a pure spin state at the cost of giving up the full variational freedom of UHF. It is conceptually cleaner but less commonly used in practice.

For most materials science applications (DFT calculations of crystals, molecules, surfaces) the analogous choice — spin-restricted vs spin-unrestricted Kohn–Sham — is the same. Magnetic systems require unrestricted treatments; non-magnetic insulators do not.

## 4.7.6b The Roothaan–Hall equations: HF on a basis

A purely formal point of view treats the HF equations as a non-linear differential equation for the orbitals. In practice, every HF (and DFT) calculation expands the orbitals in a finite basis:

$$\chi_i(\mathbf r) = \sum_\mu C_{\mu i}\, \varphi_\mu(\mathbf r),$$

where $\{\varphi_\mu\}$ are a fixed set of basis functions (Gaussian-type orbitals, Slater-type orbitals, plane waves, …) and $C_{\mu i}$ are coefficients to be determined. Substituting into the HF equation $\hat F \chi_i = \varepsilon_i \chi_i$ and projecting onto the basis $\{\varphi_\nu\}$ gives the **Roothaan–Hall equations**:

$$\mathbf F\,\mathbf C = \mathbf S\,\mathbf C\,\boldsymbol\varepsilon, \tag{4.7.RH}$$

where $\mathbf F_{\mu\nu} = \langle\varphi_\mu|\hat F|\varphi_\nu\rangle$ is the Fock matrix, $\mathbf S_{\mu\nu} = \langle\varphi_\mu|\varphi_\nu\rangle$ is the overlap matrix (identity for orthonormal bases), and $\boldsymbol\varepsilon$ is the diagonal matrix of eigenvalues. This is a **generalised eigenvalue problem** for the coefficient matrix $\mathbf C$ — exactly the kind that `scipy.linalg.eigh(F, S)` solves with one line of code.

The Roothaan–Hall recasting turns HF into a problem in linear algebra on a finite-dimensional vector space. The continuum spatial problem is gone; what remains is a $K\times K$ matrix problem where $K$ is the basis-set size. Standard chemistry basis sets give $K \sim 10$–$100$ per atom; the diagonalisation costs $\mathcal O(K^3)$. The four-index two-electron integrals (for $\mathbf F$) are the bottleneck and scale formally as $\mathcal O(K^4)$, although locality and density fitting reduce this in practice.

This is the form in which every quantum-chemistry code from the 1950s onwards has actually implemented HF. The conceptual structure — Fock operator, SCF iteration — is unchanged; the differential equation is replaced by a matrix eigenvalue problem.

## 4.7.7 A bridge to DFT

Here is the conceptual leap that makes DFT — and the rest of this book — possible. Hohenberg and Kohn proved in 1964 that the ground-state electron density $n(\mathbf r) = \sum_i |\chi_i(\mathbf r)|^2$ contains all the information of the wavefunction. The exact ground-state energy is a functional of $n$ alone, $E[n]$, even though we do not know its form. Kohn and Sham proposed in 1965 to write

$$E[n] = T_{\mathrm s}[n] + \int v_{\mathrm{ext}}(\mathbf r) n(\mathbf r) d^3 r + E_{\mathrm H}[n] + E_{\mathrm{xc}}[n], \tag{4.7.16}$$

where $T_{\mathrm s}$ is the kinetic energy of a fictitious non-interacting system with the same density, $E_{\mathrm H}$ is the classical Hartree energy, and $E_{\mathrm{xc}}$ — the **exchange–correlation functional** — absorbs everything else. Minimising (4.7.16) leads to a set of one-electron equations indistinguishable in structure from HF (4.7.11), but with the non-local exchange operator $\hat K$ replaced by the local exchange–correlation potential $v_{\mathrm{xc}}(\mathbf r) = \delta E_{\mathrm{xc}}/\delta n(\mathbf r)$.

That single replacement — *non-local exchange becomes a local exchange–correlation potential of the density* — converts HF into Kohn–Sham DFT. The accuracy of the resulting theory is entirely controlled by the quality of $E_{\mathrm{xc}}[n]$, which must be modelled. Chapter 5 is the story of those models.

For now: you have all the conceptual scaffolding you need to read the rest of the book. The many-electron Schrödinger equation is exact but intractable; the Born–Oppenheimer approximation cleaves nuclei from electrons; Hartree–Fock gives a single-determinant variational ansatz with exact exchange but no correlation; and DFT, the workhorse of Chapter 5 onwards, repackages the same one-electron structure with all the missing physics bundled into a universal density functional. Turn the page.

!!! tip "Forward reference: hybrid functionals"
    A useful preview. Modern hybrid DFT functionals (B3LYP, PBE0, HSE) mix a *fraction* of exact HF exchange into the Kohn–Sham potential:
    $$v_{\mathrm{xc}}^{\mathrm{hybrid}} = (1 - a)\,v_{\mathrm x}^{\mathrm{DFT}} + a\,v_{\mathrm x}^{\mathrm{HF}} + v_{\mathrm c}^{\mathrm{DFT}}.$$
    The HF-exchange admixture corrects for the self-interaction error of pure DFT functionals and improves band gaps in semiconductors, the description of charge-transfer states, and reaction barrier heights. The price is computational: the non-local HF exchange operator must be evaluated, which makes hybrid DFT roughly 10× more expensive than pure GGA DFT. This is the practical sense in which Hartree–Fock lives on — not as a stand-alone method, but as a component of every accurate modern density functional. Chapter 5 will develop this lineage in detail.
