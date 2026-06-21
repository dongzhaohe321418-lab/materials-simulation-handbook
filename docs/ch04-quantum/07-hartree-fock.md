# 4.7 Hartree–Fock, briefly

Within the Born–Oppenheimer approximation the central remaining task is to solve the electronic Schrödinger equation

$$\hat{H}_{\mathrm e}(\mathbf r; \mathbf R)\, \psi(\mathbf r; \mathbf R) = E(\mathbf R)\, \psi(\mathbf r; \mathbf R), \tag{4.7.1}$$

for a system of $N$ interacting electrons in the external potential of fixed nuclei. As we saw in §4.5, this problem is exponentially hard. Hartree–Fock (HF) is the simplest serious attempt to make it polynomial. The idea, conceptually, is breathtaking: assume the many-electron wavefunction is a *single* Slater determinant built from $N$ one-electron orbitals, then variationally choose those orbitals to minimise the energy. The result is a set of self-consistent one-electron equations of remarkable structure — they capture exchange exactly but neglect correlation entirely. HF is rarely used as a final method in modern materials science, but it is the conceptual scaffold on which density functional theory (Chapter 5) is built, and every electronic-structure code in the world traces some lineage to it.

This section sketches the HF construction. We will not derive the equations in full painful detail — that is a long calculation done correctly in any quantum chemistry textbook — but we will state every essential ingredient, identify the structure of the resulting equations, and pinpoint exactly what HF gets wrong.

!!! info "What problem are we solving?"
    The exact electronic Schrödinger equation (4.7.1) asks for one function $\psi(\mathbf r_1, \ldots, \mathbf r_N)$ of *all* the electron coordinates at once. Because the electrons repel one another, this function cannot be split into independent one-electron pieces, and storing or solving for it costs an effort that grows exponentially with the number of electrons $N$ (as §4.5 showed). We want a way to turn this one impossible $N$-electron problem into $N$ manageable *one*-electron problems — accepting some loss of accuracy in exchange for a calculation we can actually run. Hartree–Fock is the simplest principled way to do that.

!!! note "Plain-language version: a mean field"
    Imagine you are one electron in a crowd of $N$. Tracking the exact instantaneous position of every other electron, and dodging each one individually, is hopeless. So instead you pretend the other electrons are a *smooth, frozen cloud* of negative charge — a static "field" — and you just solve for your own orbital in that averaged cloud. Every electron does the same. Of course your cloud depends on everyone else's orbitals, and theirs on yours, so you have to iterate: guess everyone's orbitals, build the averaged cloud, re-solve, repeat until nothing changes. This "each particle in the averaged field of all the others" idea is called a **mean-field** approximation, and it is the heart of both Hartree and Hartree–Fock. The single extra ingredient that turns Hartree into Hartree–Fock is making the wavefunction properly *antisymmetric* — obeying the Pauli principle — which adds one genuinely quantum term (exchange) on top of the classical averaged repulsion.

!!! note "Physical picture: the averaged field plus the exchange hole"
    Two things are happening to each electron in HF.

    *The averaged field.* Each electron feels the pull of the nuclei plus the smeared-out electrostatic repulsion of the averaged charge cloud of all the other electrons. This is purely classical electrostatics applied to a cloud — it is the **Hartree** part, and it is what the $\hat J$ operator below encodes.

    *The exchange hole.* On top of this, the Pauli principle forbids two same-spin electrons from being at the same place. The antisymmetric (determinant) wavefunction builds in an automatic avoidance between same-spin electrons: around each electron there is a small region — the **exchange hole** — that other same-spin electrons keep out of. Because they keep their distance, their mutual repulsion is *reduced*, which *lowers* the energy. This reduction is the **exchange** energy, the $-\hat K$ term, and it has no classical counterpart at all: it comes purely from the antisymmetry of the wavefunction, not from any force in the Hamiltonian. What HF still misses is that *opposite*-spin electrons also avoid each other dynamically (they too repel), and a single determinant cannot represent that avoidance — that missing piece is called correlation.

| Symbol | Meaning | Units (SI) |
|---|---|---|
| $\mathbf x_i = (\mathbf r_i, s_i)$ | combined space–spin coordinate of electron $i$ | — |
| $\chi_i(\mathbf x)$ | spin-orbital: one-electron state (space $\times$ spin) | $\mathrm{m^{-3/2}}$ |
| $\phi_i(\mathbf r)$ | spatial part of a spin-orbital | $\mathrm{m^{-3/2}}$ |
| $\sigma_i(s)$ | spin part, $\alpha$ (up) or $\beta$ (down) | — |
| $\Psi$ | the full $N$-electron wavefunction (a determinant in HF) | $\mathrm{m^{-3N/2}}$ |
| $\hat h(\mathbf r)$ | one-electron operator: kinetic $+$ external (nuclear) potential | J |
| $\hat F$ | Fock operator (effective one-electron Hamiltonian) | J |
| $\varepsilon_i$ | orbital energy (eigenvalue of $\hat F$) | J |
| $J_{ij},\ \hat J$ | Coulomb integral / operator (classical averaged repulsion) | J |
| $K_{ij},\ \hat K$ | exchange integral / operator (quantum, same-spin only) | J |
| $\varepsilon_{ij}$ | Lagrange multipliers enforcing orbital orthonormality | J |
| $r_{12} = |\mathbf r_1 - \mathbf r_2|$ | distance between two electrons | m |
| $E_{\mathrm{corr}}$ | correlation energy, $E_{\mathrm{exact}} - E_{\mathrm{HF}}$ | J (or Ha) |

For any unfamiliar word — *operator*, *eigenvalue*, *wavefunction*, *self-consistent field* — see the [beginner glossary](../undergraduate/glossary-for-beginners.md), which defines each slowly before the formal version.

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

!!! note "Plain-language version: never undershoot the ground state"
    The variational principle says something reassuringly simple: *whatever trial wavefunction you guess, the energy you compute from it can never be lower than the true ground-state energy.* The true ground state is the bottom of the energy well; any guess sits at or above it. So if you have two guesses, the one giving the *lower* energy is the better one — and "minimise the energy over my adjustable knobs" becomes a rigorous recipe for getting as close to the truth as your trial form allows. It can never cheat by going below the real answer.

??? note "Full derivation: why the variational bound holds"
    We want to show $\langle\Psi|\hat H|\Psi\rangle \ge E_0$ for any normalised $\Psi$.

    Step 1 — expand in the exact eigenbasis. The Hamiltonian $\hat H$ has a complete orthonormal set of eigenstates $\{|\Phi_n\rangle\}$ with $\hat H|\Phi_n\rangle = E_n|\Phi_n\rangle$ and energies ordered $E_0 \le E_1 \le E_2 \le \cdots$. Any state can be written as a superposition of them:
    $$|\Psi\rangle = \sum_n c_n |\Phi_n\rangle, \qquad c_n = \langle\Phi_n|\Psi\rangle.$$

    Step 2 — use normalisation. Because $\langle\Phi_m|\Phi_n\rangle = \delta_{mn}$,
    $$\langle\Psi|\Psi\rangle = \sum_{m,n} c_m^* c_n \langle\Phi_m|\Phi_n\rangle = \sum_n |c_n|^2 = 1.$$

    Step 3 — evaluate the energy. Apply $\hat H$ to each $|\Phi_n\rangle$ inside the bracket:
    $$\langle\Psi|\hat H|\Psi\rangle = \sum_{m,n} c_m^* c_n \langle\Phi_m|\hat H|\Phi_n\rangle = \sum_{m,n} c_m^* c_n E_n \langle\Phi_m|\Phi_n\rangle = \sum_n |c_n|^2 E_n.$$

    Step 4 — bound it. Every eigenvalue obeys $E_n \ge E_0$, and every weight $|c_n|^2 \ge 0$, so replacing each $E_n$ by the smallest one $E_0$ can only lower the sum:
    $$\sum_n |c_n|^2 E_n \ \ge\ \sum_n |c_n|^2 E_0 = E_0 \sum_n |c_n|^2 = E_0.$$
    Equality holds only when all the weight sits on the ground state ($c_0 = 1$, all other $c_n = 0$), i.e. when $\Psi = \Phi_0$ exactly. $\blacksquare$

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

??? note "Full derivation: the Hartree energy, term by term"
    Start from the electronic Hamiltonian (4.5.3), written as a sum of one-electron operators plus pairwise repulsion:
    $$\hat H_{\mathrm e} = \sum_{k} \hat h(\mathbf r_k) + \tfrac12 \sum_{k \ne l} \frac{e^2}{|\mathbf r_k - \mathbf r_l|}, \qquad \hat h(\mathbf r) = -\frac{\hbar^2}{2m_{\mathrm e}}\nabla^2 + v_{\mathrm{ext}}(\mathbf r).$$
    The factor $\tfrac12$ stops us counting each pair $\{k,l\}$ twice. We sandwich this between the product state $\Psi_{\mathrm H} = \phi_1(\mathbf r_1)\cdots\phi_N(\mathbf r_N)$, assuming the orbitals are orthonormal, $\int \phi_i^*\phi_j\, d\mathbf r = \delta_{ij}$.

    *One-electron terms.* Take the term $\hat h(\mathbf r_k)$. It acts only on coordinate $\mathbf r_k$, so in
    $$\Big\langle \prod_m \phi_m(\mathbf r_m) \Big| \hat h(\mathbf r_k) \Big| \prod_n \phi_n(\mathbf r_n) \Big\rangle$$
    every integral except the one over $\mathbf r_k$ is just $\int |\phi_m|^2\, d\mathbf r_m = 1$ (normalisation). What survives is $\int \phi_k^*(\mathbf r_k)\, \hat h\, \phi_k(\mathbf r_k)\, d\mathbf r_k$. Summing over the $N$ choices of $k$:
    $$\sum_k \int \phi_k^*(\mathbf r)\,\hat h(\mathbf r)\,\phi_k(\mathbf r)\, d\mathbf r.$$

    *Two-electron terms.* Take the pair term for electrons $k$ and $l$. It acts only on $\mathbf r_k$ and $\mathbf r_l$; all other coordinates integrate to 1. What survives is
    $$\int\!\!\int \phi_k^*(\mathbf r_k)\phi_l^*(\mathbf r_l)\,\frac{e^2}{|\mathbf r_k - \mathbf r_l|}\,\phi_k(\mathbf r_k)\phi_l(\mathbf r_l)\, d\mathbf r_k\, d\mathbf r_l = \int\!\!\int \frac{e^2|\phi_k(\mathbf r_1)|^2 |\phi_l(\mathbf r_2)|^2}{|\mathbf r_1 - \mathbf r_2|}\, d\mathbf r_1 d\mathbf r_2,$$
    where the last step just renames the dummy integration variables $\mathbf r_k \to \mathbf r_1$, $\mathbf r_l \to \mathbf r_2$. Each electron's wavefunction modulus-squared $|\phi|^2$ is a charge density, so this is the classical Coulomb repulsion between two charge clouds. Reinstating the $\tfrac12 \sum_{k\ne l}$ prefactor gives the second term of $E_{\mathrm H}$ above. (Note the product ansatz has no $k=l$ self-term in the Hamiltonian, so the spurious self-repulsion never appears here — but, crucially, neither is there any antisymmetry.)

!!! tip "New vocabulary"
    - **Spin-orbital** — a one-electron state that carries *both* a spatial part $\phi(\mathbf r)$ and a spin label (up/down), written $\chi(\mathbf x) = \phi(\mathbf r)\sigma(s)$. An ordinary spatial orbital can hold two electrons (one up, one down); a spin-orbital holds exactly one.
    - **Slater determinant** — the antisymmetric many-electron wavefunction you build by arranging $N$ spin-orbitals into an $N\times N$ determinant. Swapping two electrons swaps two rows and flips the sign, which is exactly the antisymmetry the Pauli principle demands.
    - **Self-consistent field (SCF)** — the guess-build-resolve-repeat iteration used to solve the (nonlinear) mean-field equations. Defined slowly in the [beginner glossary](../undergraduate/glossary-for-beginners.md).

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

!!! info "What problem are we solving here?"
    We have a properly antisymmetric trial wavefunction (the Slater determinant) and the full electronic Hamiltonian. We now want a *formula for its energy* in terms of the unknown orbitals, so that — in the next step — we can minimise that energy over the orbitals. The pay-off of doing the algebra carefully is that the two-electron repulsion splits into exactly two pieces: a classical Coulomb piece $J$ that the Hartree picture already had, and an extra exchange piece $K$ that is *new*, that only couples same-spin electrons, and that has no classical interpretation. Seeing exactly where the minus sign on $K$ comes from is the whole point of the derivation below.

??? note "Full derivation: the determinant energy gives $J$ and $K$"
    We compute $E_{\mathrm{HF}} = \langle\Psi_{\mathrm{HF}}|\hat H_{\mathrm e}|\Psi_{\mathrm{HF}}\rangle$ for the normalised $N$-electron determinant. We use the two-electron version explicitly and then state the general pattern; the full $N$-electron result follows by the same bookkeeping (the *Slater–Condon rules*).

    **The two-electron case.** Take
    $$\Psi(\mathbf x_1,\mathbf x_2) = \frac{1}{\sqrt2}\big[\chi_a(\mathbf x_1)\chi_b(\mathbf x_2) - \chi_b(\mathbf x_1)\chi_a(\mathbf x_2)\big],$$
    and the Hamiltonian $\hat H_{\mathrm e} = \hat h(1) + \hat h(2) + g(1,2)$, with $g(1,2) = e^2/r_{12}$.

    *One-electron part.* Consider $\langle\Psi|\hat h(1)|\Psi\rangle$. Expanding the bra and ket each into two terms gives four products. Using orthonormality $\langle\chi_a|\chi_b\rangle = \delta_{ab}$ on the coordinate-2 integral, the two "cross" terms vanish (they leave a factor $\langle\chi_a|\chi_b\rangle = 0$), and the two "direct" terms each give $\tfrac12 \langle\chi_a|\hat h|\chi_a\rangle$ or $\tfrac12 \langle\chi_b|\hat h|\chi_b\rangle$. Adding $\hat h(1)$ and $\hat h(2)$ and collecting:
    $$\langle\Psi|\hat h(1)+\hat h(2)|\Psi\rangle = h_{aa} + h_{bb}, \qquad h_{ii} \equiv \langle\chi_i|\hat h|\chi_i\rangle.$$
    So the one-electron energy is simply the sum of one-electron expectation values over the occupied spin-orbitals — equation (4.7.8) summed over $i$.

    *Two-electron part.* Now $\langle\Psi|g(1,2)|\Psi\rangle$. Write $\Psi = \tfrac{1}{\sqrt2}(ab - ba)$ as shorthand, where $ab$ means $\chi_a(\mathbf x_1)\chi_b(\mathbf x_2)$. Then
    $$\langle\Psi|g|\Psi\rangle = \tfrac12\big[\langle ab|g|ab\rangle - \langle ab|g|ba\rangle - \langle ba|g|ab\rangle + \langle ba|g|ba\rangle\big].$$
    Because $g$ is symmetric under $1\leftrightarrow2$, the first and last terms are equal, and the two middle terms are equal, so
    $$\langle\Psi|g|\Psi\rangle = \langle ab|g|ab\rangle - \langle ab|g|ba\rangle.$$
    Writing these out as integrals:
    $$\langle ab|g|ab\rangle = \int\!\!\int \chi_a^*(\mathbf x_1)\chi_b^*(\mathbf x_2)\,\frac{e^2}{r_{12}}\,\chi_a(\mathbf x_1)\chi_b(\mathbf x_2)\, d\mathbf x_1 d\mathbf x_2 = J_{ab},$$
    $$\langle ab|g|ba\rangle = \int\!\!\int \chi_a^*(\mathbf x_1)\chi_b^*(\mathbf x_2)\,\frac{e^2}{r_{12}}\,\chi_b(\mathbf x_1)\chi_a(\mathbf x_2)\, d\mathbf x_1 d\mathbf x_2 = K_{ab}.$$
    So the two-electron energy is $J_{ab} - K_{ab}$. The **first** term is the ordinary Coulomb repulsion between charge clouds $|\chi_a|^2$ and $|\chi_b|^2$; the **second** — the *exchange* term — arose entirely from the cross-term in the antisymmetrised product, i.e. from the minus sign in the determinant. It has no analogue in classical electrostatics because there is no "charge density" being integrated against another: in $K$ the two orbitals are *swapped* between the bra and ket, so the integrand is $\chi_a^*(\mathbf x_1)\chi_b(\mathbf x_1)$ — a product of two *different* orbitals at the same point, which is not a density.

    **Why exchange is same-spin only.** Each spin-orbital factorises as $\chi_i(\mathbf x) = \phi_i(\mathbf r)\sigma_i(s)$ with $\sigma$ either $\alpha$ (up) or $\beta$ (down), and the spins are orthonormal: $\sum_s \alpha^*(s)\alpha(s) = 1$, $\sum_s \alpha^*(s)\beta(s) = 0$. In $K_{ab}$ the spin sum over coordinate 1 is $\sum_{s_1}\sigma_a^*(s_1)\sigma_b(s_1)$ and over coordinate 2 is $\sum_{s_2}\sigma_b^*(s_2)\sigma_a(s_2)$. If $a$ and $b$ have *opposite* spin, each of these factors is zero, so $K_{ab} = 0$. In $J_{ab}$, by contrast, the spin sums are $\sum_{s_1}|\sigma_a(s_1)|^2 = 1$ and $\sum_{s_2}|\sigma_b(s_2)|^2 = 1$ regardless of spin, so $J$ survives for any spin pairing. Hence: Coulomb acts between all pairs; exchange acts only between same-spin pairs.

    **General $N$.** For an $N\times N$ determinant the same logic — keep only terms that survive orthonormality — gives, for the two-electron part, a sum over all *ordered* pairs of every $J_{ij}$ minus every $K_{ij}$, with a $\tfrac12$ to avoid double-counting:
    $$E_{\mathrm{HF}} = \sum_i h_{ii} + \frac12 \sum_{i,j} \big(J_{ij} - K_{ij}\big),$$
    which is (4.7.7). The double sum includes $i = j$: there $J_{ii} = K_{ii}$ (the exchange-of-identical-orbitals integral equals the Coulomb-of-identical-orbitals integral), so $J_{ii} - K_{ii} = 0$ and the electron contributes no spurious self-repulsion — the self-interaction cancels exactly, as noted below.

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

!!! info "What problem are we solving?"
    We have the energy (4.7.7) as a functional of the orbitals. The variational principle tells us the *best* single determinant is the one whose orbitals make this energy as small as possible. So we now minimise $E_{\mathrm{HF}}$ over the orbitals. The one complication is that we are not free to vary the orbitals arbitrarily — they must stay orthonormal (a determinant of non-orthonormal orbitals is not normalised, and equal orbitals make it vanish). Minimising a quantity subject to constraints is exactly what **Lagrange multipliers** are for, and the output of that minimisation is the set of Hartree–Fock equations.

!!! note "Plain-language version: orbitals that don't lower the energy any further"
    The Hartree–Fock equations are just the statement "I have wiggled every orbital in every possible way, keeping them orthonormal, and the energy no longer goes down." At that stationary point each orbital satisfies a Schrödinger-like equation $\hat F\chi_i = \varepsilon_i\chi_i$, where the effective Hamiltonian $\hat F$ (the Fock operator) already contains the averaged repulsion from all the orbitals. Because $\hat F$ is built from the very orbitals we are solving for, the equation has to be solved by iteration — guess, build $\hat F$, solve, repeat.

Minimise $E_{\mathrm{HF}}$ with respect to the spin-orbitals $\chi_i$, subject to orthonormality $\langle\chi_i|\chi_j\rangle = \delta_{ij}$. Use Lagrange multipliers $\varepsilon_{ij}$ to enforce the constraints, take the variation, and diagonalise the multiplier matrix. The result is the canonical form of the **Hartree–Fock equations**:

??? note "Full derivation: from energy minimisation to $\hat F\chi_i = \varepsilon_i\chi_i$"
    **Set up the constrained problem.** We minimise $E_{\mathrm{HF}}[\{\chi_i\}]$ of (4.7.7) subject to $\langle\chi_i|\chi_j\rangle = \delta_{ij}$ for all $i,j$. Introduce one Lagrange multiplier $\varepsilon_{ji}$ for each constraint and form the Lagrangian
    $$\mathcal L = E_{\mathrm{HF}} - \sum_{i,j} \varepsilon_{ji}\big(\langle\chi_i|\chi_j\rangle - \delta_{ij}\big).$$
    At the minimum, $\mathcal L$ is stationary against any small variation $\chi_i \to \chi_i + \delta\chi_i$.

    **Vary the energy.** Treat $\chi_i^*$ and $\chi_i$ as independent (the standard trick for complex functions: making $\mathcal L$ stationary in $\chi_i^*$ also makes it stationary in $\chi_i$, by taking the complex conjugate). Vary $\chi_i^* \to \chi_i^* + \delta\chi_i^*$. From the one-electron part $\sum_i h_{ii}$:
    $$\delta\Big(\sum_i h_{ii}\Big) = \langle\delta\chi_i|\hat h|\chi_i\rangle.$$
    From the two-electron part $\tfrac12\sum_{j,k}(J_{jk} - K_{jk})$, only the $j=i$ and $k=i$ terms involve $\chi_i^*$. Because the expression is symmetric in its two indices, the factor of $\tfrac12$ and the two equal contributions ($j=i$ and $k=i$) cancel, leaving
    $$\delta\Big(\tfrac12\sum_{j,k}(J_{jk}-K_{jk})\Big) = \sum_j \Big\langle\delta\chi_i\Big|\,\hat J_j - \hat K_j\,\Big|\chi_i\Big\rangle,$$
    where $\hat J_j$ and $\hat K_j$ are the one-electron Coulomb and exchange operators built from orbital $\chi_j$ (defined explicitly in (4.7.13)–(4.7.14), summed over $j$). The variation of the constraint term is $\sum_j \varepsilon_{ji}\langle\delta\chi_i|\chi_j\rangle$.

    **Collect terms.** Setting $\delta\mathcal L = 0$ for arbitrary $\delta\chi_i$ means the bracketed coefficient must vanish:
    $$\Big[\hat h + \sum_j(\hat J_j - \hat K_j)\Big]\chi_i = \sum_j \varepsilon_{ji}\,\chi_j.$$
    Define the **Fock operator** $\hat F = \hat h + \sum_j(\hat J_j - \hat K_j)$. Then
    $$\hat F\,\chi_i = \sum_j \varepsilon_{ji}\,\chi_j. \tag{4.7.4a}$$
    This is the Hartree–Fock equation in its raw form: the Fock operator acting on one orbital gives a *mixture* of all the orbitals, weighted by the multiplier matrix $\varepsilon_{ji}$.

    **Diagonalise.** The multiplier matrix $\boldsymbol\varepsilon = (\varepsilon_{ji})$ is Hermitian (it can be shown that $\varepsilon_{ji} = \varepsilon_{ij}^*$ because $\hat F$ is Hermitian and the $\chi_i$ orthonormal). Any Hermitian matrix can be diagonalised by a unitary transformation $\mathbf U$ of the orbitals, $\chi_i' = \sum_k U_{ki}\chi_k$. Crucially, a *unitary mixing of the occupied orbitals leaves the determinant unchanged* (up to an overall phase) and therefore leaves both $\hat F$ and the total energy unchanged — so we are free to choose the rotation that makes $\boldsymbol\varepsilon$ diagonal, $\varepsilon_{ji} \to \varepsilon_i\delta_{ji}$. In this **canonical** basis (4.7.4a) collapses to
    $$\hat F\,\chi_i = \varepsilon_i\,\chi_i,$$
    a genuine eigenvalue equation — equation (4.7.11). The diagonal multipliers $\varepsilon_i$ are now the orbital energies, and Koopmans' theorem (below) gives them physical meaning.

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

??? note "Full derivation: $\varepsilon_i$ equals the energy of removing electron $i$"
    Two facts combine to give Koopmans' theorem. First we show $\varepsilon_i = h_{ii} + \sum_{j\ne i}(J_{ij}-K_{ij})$; then we show the energy difference equals the same thing.

    **Step 1 — the orbital energy is the diagonal Fock element.** Left-multiply the canonical HF equation $\hat F\chi_i = \varepsilon_i\chi_i$ by $\chi_i^*$ and integrate. Using $\langle\chi_i|\chi_i\rangle = 1$,
    $$\varepsilon_i = \langle\chi_i|\hat F|\chi_i\rangle = \langle\chi_i|\hat h|\chi_i\rangle + \sum_j \langle\chi_i|\hat J_j - \hat K_j|\chi_i\rangle.$$
    The first term is $h_{ii}$. From the definitions (4.7.13)–(4.7.14), $\langle\chi_i|\hat J_j|\chi_i\rangle = J_{ij}$ and $\langle\chi_i|\hat K_j|\chi_i\rangle = K_{ij}$. The $j=i$ term contributes $J_{ii}-K_{ii} = 0$ (self-interaction cancels). Hence
    $$\varepsilon_i = h_{ii} + \sum_{j\ne i}(J_{ij}-K_{ij}). \tag{4.7.Ka}$$
    In words: an orbital energy is its own one-electron energy *plus* its Coulomb–exchange interaction with every *other* occupied orbital.

    **Step 2 — the energy of removing electron $i$.** Write the total HF energy (4.7.7) by separating the orbital $i$ from the rest. Splitting the double sum into the $i$-row, the $i$-column, and the remainder:
    $$E(\Psi_N) = \underbrace{\sum_{k} h_{kk}}_{\text{includes }h_{ii}} + \tfrac12\sum_{k,l}(J_{kl}-K_{kl}).$$
    Now remove orbital $i$ (the frozen-orbital assumption keeps all other $\chi_j$ unchanged). The new energy $E(\Psi_{N-1})$ is the same sums but with $i$ deleted from the orbital list. Subtracting,
    $$E(\Psi_N) - E(\Psi_{N-1}) = h_{ii} + \tfrac12\Big[\underbrace{\sum_{j}(J_{ij}-K_{ij})}_{k=i\text{ row}} + \underbrace{\sum_{j}(J_{ji}-K_{ji})}_{l=i\text{ column}}\Big].$$
    The $i$–$i$ self-term ($J_{ii}-K_{ii}=0$) drops out, and since $J_{ij}=J_{ji}$ and $K_{ij}=K_{ji}$ (the integrals are symmetric in their two indices) the row sum and column sum are equal. The two halves combine, the $\tfrac12$ cancels, and
    $$E(\Psi_N) - E(\Psi_{N-1}) = h_{ii} + \sum_{j\ne i}(J_{ij}-K_{ij}) = \varepsilon_i,$$
    using (4.7.Ka) for the last equality. The ionisation energy is the energy you must *put in* to remove the electron, $I_i = E(\Psi_{N-1}) - E(\Psi_N) = -\varepsilon_i$. Because bound orbitals have $\varepsilon_i < 0$, this is positive, as it must be. $\blacksquare$

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

!!! warning "Common misunderstandings"
    - **"Hartree–Fock includes electron correlation."** It does not — and this is true *by definition*. The standard definition of the correlation energy is $E_{\mathrm{corr}} = E_{\mathrm{exact}} - E_{\mathrm{HF}}$ (equation (4.7.15)): correlation is *precisely the part of the energy that HF leaves out*. A single Slater determinant treats the averaged repulsion (Coulomb) and the same-spin avoidance (exchange) exactly, but it cannot represent opposite-spin electrons dodging each other instant by instant. So saying "HF includes correlation" is a contradiction in terms. HF does include *exchange* exactly — do not confuse the two; exchange is the same-spin Pauli effect, correlation is everything beyond the single-determinant mean field.
    - **"Hartree–Fock is just DFT (or DFT is just HF)."** They look structurally similar — both produce one-electron equations solved self-consistently — but they are different theories. HF is an *approximation*: a single determinant that is provably missing correlation. Kohn–Sham DFT is, *in principle, exact*: there exists an exact exchange–correlation functional $E_{\mathrm{xc}}[n]$ that would give the true ground-state energy and density. The catch is that we do not know that functional and must approximate it. So HF is "exact ansatz solved exactly but missing physics"; DFT is "exact theory solved with an approximate ingredient". They also differ in the exchange term: HF uses the exact *non-local* exchange operator $\hat K$, whereas DFT folds exchange and correlation together into a (usually) *local* potential $v_{\mathrm{xc}}(\mathbf r)$. See the forward-reference box at the end of this section.
    - **"The orbital energies $\varepsilon_i$ add up to the total energy."** They do not. Summing (4.7.Ka) over all occupied $i$ double-counts every electron–electron interaction, so $\sum_i\varepsilon_i = \sum_i h_{ii} + \sum_{ij}(J_{ij}-K_{ij})$ — which exceeds the true $E_{\mathrm{HF}}$ of (4.7.7) by $\tfrac12\sum_{ij}(J_{ij}-K_{ij})$. The total energy must be computed from (4.7.7), not from the eigenvalue sum.

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

!!! tip "Where this appears later"
    Everything in this section reappears, repackaged, in [Chapter 5 (DFT)](../ch05-dft/index.md):

    - The **single Slater determinant** returns as the Kohn–Sham determinant of fictitious one-electron orbitals (Section 5.2 there).
    - The **Fock operator's structure** (kinetic + external + Hartree + exchange) becomes the Kohn–Sham Hamiltonian, with exchange and correlation bundled into $v_{\mathrm{xc}}$.
    - The **SCF cycle** of Section 4.7.4 is run in practice for real materials in [Chapter 6 (running DFT)](../ch06-running-dft/index.md).
    - The **self-interaction cancellation** that HF achieves exactly becomes the "self-interaction error" that approximate functionals struggle with — a recurring theme in Chapter 5.

    The key contrast to carry forward: Kohn–Sham DFT has the *same one-electron structure* as HF but is, in principle, *exact* — provided one knew the exact exchange–correlation functional. HF is exact-by-construction for exchange but provably missing correlation; DFT is exact-in-principle but limited in practice by the functional we choose.

Here is the conceptual leap that makes DFT — and the rest of this book — possible. Hohenberg and Kohn proved in 1964 that the ground-state electron density $n(\mathbf r) = \sum_i |\chi_i(\mathbf r)|^2$ contains all the information of the wavefunction. The exact ground-state energy is a functional of $n$ alone, $E[n]$, even though we do not know its form. Kohn and Sham proposed in 1965 to write

$$E[n] = T_{\mathrm s}[n] + \int v_{\mathrm{ext}}(\mathbf r) n(\mathbf r) d^3 r + E_{\mathrm H}[n] + E_{\mathrm{xc}}[n], \tag{4.7.16}$$

where $T_{\mathrm s}$ is the kinetic energy of a fictitious non-interacting system with the same density, $E_{\mathrm H}$ is the classical Hartree energy, and $E_{\mathrm{xc}}$ — the **exchange–correlation functional** — absorbs everything else. Minimising (4.7.16) leads to a set of one-electron equations indistinguishable in structure from HF (4.7.11), but with the non-local exchange operator $\hat K$ replaced by the local exchange–correlation potential $v_{\mathrm{xc}}(\mathbf r) = \delta E_{\mathrm{xc}}/\delta n(\mathbf r)$.

That single replacement — *non-local exchange becomes a local exchange–correlation potential of the density* — converts HF into Kohn–Sham DFT. The accuracy of the resulting theory is entirely controlled by the quality of $E_{\mathrm{xc}}[n]$, which must be modelled. Chapter 5 is the story of those models.

For now: you have all the conceptual scaffolding you need to read the rest of the book. The many-electron Schrödinger equation is exact but intractable; the Born–Oppenheimer approximation cleaves nuclei from electrons; Hartree–Fock gives a single-determinant variational ansatz with exact exchange but no correlation; and DFT, the workhorse of Chapter 5 onwards, repackages the same one-electron structure with all the missing physics bundled into a universal density functional. Turn the page.

!!! tip "Forward reference: hybrid functionals"
    A useful preview. Modern hybrid DFT functionals (B3LYP, PBE0, HSE) mix a *fraction* of exact HF exchange into the Kohn–Sham potential:
    $$v_{\mathrm{xc}}^{\mathrm{hybrid}} = (1 - a)\,v_{\mathrm x}^{\mathrm{DFT}} + a\,v_{\mathrm x}^{\mathrm{HF}} + v_{\mathrm c}^{\mathrm{DFT}}.$$
    The HF-exchange admixture corrects for the self-interaction error of pure DFT functionals and improves band gaps in semiconductors, the description of charge-transfer states, and reaction barrier heights. The price is computational: the non-local HF exchange operator must be evaluated, which makes hybrid DFT roughly 10× more expensive than pure GGA DFT. This is the practical sense in which Hartree–Fock lives on — not as a stand-alone method, but as a component of every accurate modern density functional. Chapter 5 will develop this lineage in detail.

## 4.7.8 Check yourself

!!! question "Check yourself"
    1. The two-electron part of the HF energy is $J_{ij} - K_{ij}$. Which of these two terms has a classical electrostatic interpretation, and which has none? In one sentence, where did the term *without* a classical interpretation come from in the derivation?
    2. Show in one line that an electron does not spuriously repel itself in HF. (What is $J_{ii} - K_{ii}$, and why?)
    3. The exchange integral $K_{ij}$ is zero when orbitals $i$ and $j$ have opposite spin. Which integral in the derivation forced this, and what mathematical fact about the spin functions made it vanish?
    4. The Hartree–Fock equation $\hat F\chi_i = \varepsilon_i\chi_i$ looks like a one-electron Schrödinger equation. Why can you *not* simply diagonalise $\hat F$ once and read off the answer? What is the name of the iteration you must run instead?
    5. A colleague says "Hartree–Fock already includes electron correlation, since the determinant correlates the electrons' positions." Correct them precisely: what *does* the determinant correlate, what does it *not*, and how is the correlation energy defined?

    ??? success "Answer"
        1. $J_{ij}$ is classical — it is the Coulomb repulsion between the charge clouds $|\chi_i|^2$ and $|\chi_j|^2$. $K_{ij}$ has *no* classical interpretation; in its integrand the two orbitals are swapped between bra and ket ($\chi_i^*(\mathbf x_1)\chi_j(\mathbf x_1)$), so it is not the interaction of two densities. It arose from the *cross-term* in the antisymmetrised (determinant) product — i.e. directly from the minus sign that enforces the Pauli principle.
        2. $J_{ii} - K_{ii} = 0$, because setting $i=j$ in the definitions (4.7.9) and (4.7.10) makes the two integrals identical (the swap $j\to i$ does nothing when $j$ already equals $i$). Their difference vanishes, so an electron has no self-interaction in HF.
        3. The exchange integral $K_{ij}$. Writing $\chi = \phi\,\sigma$, the spin sums in $K_{ij}$ are $\sum_{s}\sigma_i^*(s)\sigma_j(s)$, which equals zero when the spins differ because the up and down spin functions are orthonormal ($\langle\alpha|\beta\rangle = 0$). In $J_{ij}$ the spin sums are instead $\sum_s|\sigma(s)|^2 = 1$, so $J$ survives for any spin pairing.
        4. Because $\hat F$ is built from the Coulomb and exchange operators $\hat J$ and $\hat K$, which themselves depend on *all* the occupied orbitals $\chi_j$ — the equation is nonlinear. You do not know $\hat F$ until you know its own eigenvectors. So you guess the orbitals, build $\hat F$, diagonalise to get new orbitals, and repeat to convergence. This is the **self-consistent field (SCF)** iteration.
        5. The determinant correlates the positions of *same-spin* electrons only — that is the exchange (Pauli) effect, and HF captures it exactly. It does *not* correlate opposite-spin electrons: they may, within HF, sit on top of one another paying only the average Coulomb cost. The missing dynamical avoidance of opposite-spin (and the residual same-spin) electrons is the **correlation energy**, defined as $E_{\mathrm{corr}} = E_{\mathrm{exact}} - E_{\mathrm{HF}}$ — i.e. exactly what HF leaves out, by definition.

    ??? note "Hint"
        For 1–3, look at the "Full derivation: the determinant energy gives $J$ and $K$" box: the exchange term comes from the cross-term in the antisymmetrised product, and the spin sums decide whether $K$ survives. For 4, recall that $\hat F$ in (4.7.12) contains $\hat J$ and $\hat K$, which are themselves built from the orbitals. For 5, the definition of $E_{\mathrm{corr}}$ in (4.7.15) is the whole answer.
