# 4.5 Many electrons and the exponential wall

So far we have considered a single quantum particle in a one-dimensional potential. Diagonalising a 400×400 tridiagonal matrix on a laptop returns the lowest several eigenstates of the particle-in-a-box or the harmonic oscillator to seven significant figures in a fraction of a second. It is tempting to assume that the same approach, scaled up, will work for a real material — a few hundred electrons in three dimensions on a finer grid. This section explains why that assumption is wrong by an absurd margin.

The "many-body problem" is the central computational difficulty of quantum chemistry and condensed-matter physics. Two ingredients together produce it. First, the wavefunction of $N$ particles is a function not of three coordinates but of $3N$; the configuration space scales linearly with $N$, but the number of grid points needed to *sample* it scales exponentially. Second, the Pauli principle forces the wavefunction to be antisymmetric under particle exchange, which both rules out the simplest product ansatz and dictates a particular algebraic structure (Slater determinants). The interplay is what makes interacting fermion systems hard.

## 4.5.0 Warm-up: non-interacting electrons

Before writing the full Hamiltonian, consider the *limit* in which electron–electron interaction is switched off. The Hamiltonian for $N$ non-interacting electrons in an external potential $v_{\mathrm{ext}}(\mathbf r)$ is

$$\hat H_0 = \sum_{i=1}^N \hat h(\mathbf r_i), \qquad \hat h(\mathbf r) = -\frac{\hbar^2}{2m_{\mathrm e}}\nabla^2 + v_{\mathrm{ext}}(\mathbf r).$$

Because $\hat H_0$ is a *sum* of one-electron operators, each acting on a different coordinate, the eigenvalue equation $\hat H_0 \Psi = E \Psi$ separates. Try the product ansatz

$$\Psi(\mathbf r_1, \ldots, \mathbf r_N) = \phi_{a_1}(\mathbf r_1)\,\phi_{a_2}(\mathbf r_2)\cdots\phi_{a_N}(\mathbf r_N),$$

where each $\phi_{a_i}$ is a single-particle eigenstate of $\hat h$ with eigenvalue $\varepsilon_{a_i}$. Substituting,

$$\hat H_0\,\Psi = \sum_i \hat h(\mathbf r_i) \prod_j \phi_{a_j}(\mathbf r_j) = \sum_i \varepsilon_{a_i}\,\Psi = E\,\Psi,$$

so $E = \sum_i \varepsilon_{a_i}$ is the sum of one-electron eigenvalues. **The non-interacting $N$-electron problem reduces to $N$ independent single-electron problems** — exactly the kind of problem we solved in §4.3 and §4.4.

This is a wonderful state of affairs and it is the structure that underpins band theory (Chapter 3), tight-binding models, and ultimately the *Kohn–Sham* formulation of DFT, which represents an interacting system *as if it were* a non-interacting one in a cleverly chosen effective potential. The trouble, as we are about to see, is that real electrons interact through Coulomb repulsion, and the moment we turn that interaction back on the separability is destroyed.

!!! note "Why this step?"
    The factorisation $\hat H_0 = \sum_i \hat h_i$ works precisely because each $\hat h_i$ acts on coordinate $\mathbf r_i$ only, leaving the others untouched. The product wavefunction is then an eigenfunction by elementary calculus (a partial derivative with respect to $\mathbf r_i$ touches only $\phi_{a_i}$). Adding any term that couples two coordinates — even a single $V(\mathbf r_i, \mathbf r_j)$ for one pair — destroys separability and forces us to deal with the joint $3N$-dimensional problem.

## 4.5.1 The full Hamiltonian for a real material

Consider an arbitrary molecule or solid: $N_{\mathrm e}$ electrons (mass $m_{\mathrm e}$, charge $-e$, positions $\mathbf r_i$) and $N_{\mathrm n}$ nuclei (mass $M_I$, charge $+Z_I e$, positions $\mathbf R_I$). The non-relativistic Hamiltonian, in Gaussian units for compactness, is

$$\hat{H} = \hat T_{\mathrm e} + \hat T_{\mathrm n} + \hat V_{\mathrm{ee}} + \hat V_{\mathrm{en}} + \hat V_{\mathrm{nn}}, \tag{4.5.1}$$

with the five terms

$$\hat T_{\mathrm e} = -\frac{\hbar^2}{2m_{\mathrm e}}\sum_{i=1}^{N_{\mathrm e}} \nabla_i^2, \qquad \hat T_{\mathrm n} = -\sum_{I=1}^{N_{\mathrm n}} \frac{\hbar^2}{2M_I} \nabla_I^2, \tag{4.5.2}$$

$$\hat V_{\mathrm{ee}} = \frac{1}{2}\sum_{i\neq j} \frac{e^2}{|\mathbf r_i - \mathbf r_j|}, \qquad \hat V_{\mathrm{en}} = -\sum_{i, I}\frac{Z_I e^2}{|\mathbf r_i - \mathbf R_I|}, \qquad \hat V_{\mathrm{nn}} = \frac{1}{2}\sum_{I\neq J}\frac{Z_I Z_J e^2}{|\mathbf R_I - \mathbf R_J|}. \tag{4.5.3}$$

!!! note "Reading the indices carefully"
    The sums are over different kinds of indices and it is worth pausing to make sure you parse them correctly.
    - $\hat T_{\mathrm e}$ sums over electron labels $i = 1, \ldots, N_{\mathrm e}$. Each term involves the Laplacian $\nabla_i^2$ with respect to the position $\mathbf r_i$ of electron $i$.
    - $\hat T_{\mathrm n}$ sums over nuclear labels $I = 1, \ldots, N_{\mathrm n}$, with $\nabla_I^2$ acting on $\mathbf R_I$. Capital $I$ versus lower-case $i$ is the convention that distinguishes nuclei from electrons throughout this book.
    - $\hat V_{\mathrm{ee}}$ is a double sum over *pairs* of electrons. The factor of $1/2$ corrects double-counting since $(i,j)$ and $(j,i)$ refer to the same pair; $i \neq j$ excludes self-interaction.
    - $\hat V_{\mathrm{en}}$ couples every electron $i$ to every nucleus $I$. There is no factor of $1/2$ because $i$ and $I$ run over different sets.
    - $\hat V_{\mathrm{nn}}$ is a double sum over nuclei pairs; the factor of $1/2$ and $I \neq J$ play the same role as in $\hat V_{\mathrm{ee}}$.
    Memorise this index convention now. The same letters will appear in HF (§4.7) and in every DFT discussion thereafter.

That is it. *Every* property of every material — bond lengths, lattice constants, elastic moduli, band gaps, magnetisation, superconductivity, ferroelectricity, thermal conductivity — is encoded in solving the eigenvalue equation $\hat{H} \Psi = E \Psi$ for this operator. Dirac, having written down a similar Hamiltonian in 1929, declared:

> "The underlying physical laws necessary for the mathematical theory of a large part of physics and the whole of chemistry are thus completely known, and the difficulty is only that the exact application of these laws leads to equations much too complicated to be soluble."

The history of computational materials science is the history of grappling with that "only". In this section we measure the size of the problem and convince ourselves that it cannot be solved by brute force.

## 4.5.2 The wavefunction lives in $3N$-dimensional space

In §4.3 we discretised the wavefunction of a single particle on a 1D grid of, say, $N_g = 400$ points. The wavefunction was represented by 400 complex numbers. The Hamiltonian was a $400\times 400$ matrix — 160 000 entries, of which only $\sim 1200$ were non-zero in the tridiagonal structure. Trivial.

Now consider $N$ electrons in three dimensions. The many-electron wavefunction is a function

$$\Psi(\mathbf r_1, \mathbf r_2, \ldots, \mathbf r_N) \tag{4.5.4}$$

of $3N$ continuous coordinates (we are ignoring spin for the moment). To sample $\Psi$ on a grid with $N_g$ points per spatial direction, we need

$$N_g^{3N} \tag{4.5.5}$$

complex numbers. The Hamiltonian becomes a matrix of size $N_g^{3N} \times N_g^{3N}$. This is the **exponential wall**.

Let us put numbers on it. Take a coarse spatial grid of $N_g = 10$ points per dimension — laughably under-resolved for any real chemistry, but a good lower bound. Then:

| $N$ electrons | grid coords $3N$ | basis states $10^{3N}$ |
|---:|---:|---:|
| 1 | 3 | $10^{3}$ |
| 2 | 6 | $10^{6}$ |
| 3 | 9 | $10^{9}$ |
| 5 | 15 | $10^{15}$ |
| 10 | 30 | $10^{30}$ |
| 20 | 60 | $10^{60}$ |
| 100 | 300 | $10^{300}$ |

The third row, three electrons, already gives a billion-dimensional Hilbert space. The fifth row, ten electrons — roughly a water molecule — gives $10^{30}$ basis functions. A single double-precision complex number occupies 16 bytes, so storing $\Psi$ for ten electrons on a 10×10×10 grid would require $1.6 \times 10^{31}$ bytes. The total digital data created by humanity to date is approximately $10^{23}$ bytes. Just *writing down* the wavefunction of ten electrons on a coarse 1000-point grid would require **a hundred million times the entire world's data storage**.

!!! example "Putting the absurdity in context"
    Take 10 electrons (a water molecule) on a $10\times 10\times 10$ grid: $10^{30}$ basis states. The number of atoms in the Earth is approximately $1.3\times 10^{50}$. So the wavefunction of *water* on a 1000-point coarse grid has fewer entries than there are atoms in the Earth — but only by twenty orders of magnitude, which on a logarithmic scale leaves us essentially nowhere. Refining the grid to chemical resolution (100 points per direction) brings the basis up to $10^{60}$, which is more than the atoms in the *solar system* ($\sim 10^{57}$). The brute-force approach is not just impractical, it is *cosmologically* impossible.

!!! warning "It gets worse"
    Even if we could store $\Psi$, we would still need to diagonalise the Hamiltonian. Dense diagonalisation costs $\mathcal O(D^3)$ operations for a $D\times D$ matrix; sparse Lanczos methods can reach $\mathcal O(D \cdot \text{iterations})$ but still need to apply $\hat{H}$ to a vector of length $D$. For our ten-electron problem, $D = 10^{30}$ — and the fastest supercomputers in 2026 perform of order $10^{18}$ floating-point operations per second. A single matrix–vector product would take $10^{12}$ seconds, roughly thirty thousand years.

This is **not** a problem that will be solved by Moore's law. Doubling our compute power every two years lets us add one electron to the calculation every six years or so. At that rate, going from 10 electrons to 30 (a small organic molecule) would take 120 years of hardware improvement. The Schrödinger equation must be tackled differently.

## 4.5.3 Pauli antisymmetry

Worse: a generic function $\Psi(\mathbf r_1, \ldots, \mathbf r_N)$ is not even an admissible wavefunction for electrons. The Pauli principle, stated formally, is the postulate that the wavefunction of identical *fermions* (electrons, protons, neutrons, ...) must change sign under exchange of any two particles:

$$\Psi(\ldots, \mathbf r_i, \ldots, \mathbf r_j, \ldots) = -\Psi(\ldots, \mathbf r_j, \ldots, \mathbf r_i, \ldots). \tag{4.5.6}$$

(Strictly the exchange acts on combined space-and-spin coordinates $\mathbf x_i = (\mathbf r_i, \sigma_i)$, but in this section we will keep spin implicit.) Identical *bosons* (photons, $^4$He nuclei, …) take a plus sign instead. This is the **spin-statistics theorem**, derived in relativistic quantum field theory and accepted as a postulate in non-relativistic quantum mechanics.

An immediate corollary: if two electrons are in the same single-particle state ($\mathbf r_i = \mathbf r_j$ and same spin), then (4.5.6) demands $\Psi = -\Psi$, so $\Psi = 0$. This is the exclusion principle: *no two electrons can occupy the same one-electron state*. Stuffing more than two electrons into the lowest level of a particle-in-a-box (one spin-up, one spin-down) is forbidden; the third electron must go in $n = 2$. This is what makes atoms beyond hydrogen possess shell structure, and what makes the periodic table look the way it does.

### Why antisymmetry — a two-electron example

It is worth seeing antisymmetry in action on the smallest possible system. Take two electrons in single-particle orbitals $\phi_a, \phi_b$ (orthonormal, both spatial only — we suppress spin for now). The product ansatz $\Psi_{\mathrm{prod}}(\mathbf r_1, \mathbf r_2) = \phi_a(\mathbf r_1) \phi_b(\mathbf r_2)$ treats the electrons as distinguishable: electron 1 is in $\phi_a$, electron 2 is in $\phi_b$. But electrons are *identical* — there is no fact of the matter about "which one is which" — and the wavefunction must reflect this.

The two possibilities consistent with identicalness are the *symmetric* and *antisymmetric* combinations:

$$\Psi_\pm(\mathbf r_1, \mathbf r_2) = \frac{1}{\sqrt 2}\bigl[\phi_a(\mathbf r_1)\phi_b(\mathbf r_2) \pm \phi_b(\mathbf r_1)\phi_a(\mathbf r_2)\bigr].$$

Under exchange $\mathbf r_1 \leftrightarrow \mathbf r_2$, $\Psi_+$ is unchanged (bosons) and $\Psi_-$ acquires a minus sign (fermions, including electrons).

Now evaluate the joint probability density at $\mathbf r_1 = \mathbf r_2 \equiv \mathbf r$. For the antisymmetric combination,

$$\Psi_-(\mathbf r, \mathbf r) = \frac{1}{\sqrt 2}\bigl[\phi_a(\mathbf r)\phi_b(\mathbf r) - \phi_b(\mathbf r)\phi_a(\mathbf r)\bigr] = 0.$$

The probability of finding both electrons at the same point is exactly zero, regardless of the orbitals $\phi_a, \phi_b$. This is the "exchange hole" — a quantum-mechanical avoidance, present *even without Coulomb repulsion*, between electrons of the same spin. It is the origin of why electrons in different shells of an atom don't collapse together: Pauli, not Coulomb, holds them apart.

For the symmetric (bosonic) combination, by contrast, $\Psi_+(\mathbf r, \mathbf r) = \sqrt 2 \phi_a(\mathbf r)\phi_b(\mathbf r) \neq 0$ — bosons *like* to be together (bunching).

This single calculation, no more elaborate than the algebra above, is the seed from which exchange interactions in magnetism, the Slater determinant in §4.5.4, and the Fermi sea in solid-state physics all grow.

!!! tip "Spin enters the picture"
    With spin included, the full wavefunction is antisymmetric under exchange of the combined space-spin coordinates. This permits two electrons of *opposite* spin to occupy the same spatial orbital (the spatial part is symmetric, the spin part antisymmetric — the *singlet*); but two electrons of *parallel* spin must have an antisymmetric spatial part (and a symmetric *triplet* spin), forcing them into different spatial orbitals. The pattern "two opposite-spin electrons per spatial orbital" is what makes the periodic table look the way it does.

For numerical purposes, antisymmetry restricts the wavefunction to a subspace — the antisymmetric subspace of $L^2(\mathbb R^{3N})$. Its dimension is much smaller than $N_g^{3N}$ (only $\binom{N_g^3}{N}$ ways to occupy $N_g^3$ orbitals with $N$ electrons), but for any non-trivial $N$ this is still astronomical.

### Worked example: helium singlet–triplet splitting

The simplest place in nature where antisymmetry has a *measurable* consequence — beyond the structure of the periodic table itself — is the excited-state spectrum of helium. He has two electrons; in its ground state both occupy the $1s$ orbital, one with spin up and one with spin down. Excite one electron to the $2s$ orbital and you have an electron in $1s$ and an electron in $2s$, with two possible spin configurations:

- **Para-helium (singlet, $S = 0$).** Spatial part symmetric, $\Psi_{\mathrm{sym}}^{\mathrm{space}}(\mathbf r_1, \mathbf r_2) = \tfrac{1}{\sqrt 2}[\phi_{1s}(\mathbf r_1)\phi_{2s}(\mathbf r_2) + \phi_{2s}(\mathbf r_1)\phi_{1s}(\mathbf r_2)]$; spin part antisymmetric, $\tfrac{1}{\sqrt 2}(|\!\uparrow\downarrow\rangle - |\!\downarrow\uparrow\rangle)$.
- **Ortho-helium (triplet, $S = 1$).** Spatial part antisymmetric, $\Psi_{\mathrm{anti}}^{\mathrm{space}} = \tfrac{1}{\sqrt 2}[\phi_{1s}(\mathbf r_1)\phi_{2s}(\mathbf r_2) - \phi_{2s}(\mathbf r_1)\phi_{1s}(\mathbf r_2)]$; spin part one of the three symmetric combinations $|\!\uparrow\uparrow\rangle$, $|\!\downarrow\downarrow\rangle$, $\tfrac{1}{\sqrt 2}(|\!\uparrow\downarrow\rangle + |\!\downarrow\uparrow\rangle)$.

Now compute the Coulomb repulsion expectation value $\langle \Psi | e^2/r_{12} | \Psi\rangle$ for the two cases. After a few lines of algebra,

$$E_\pm = J \pm K, \qquad J = \int |\phi_{1s}(\mathbf r_1)|^2 |\phi_{2s}(\mathbf r_2)|^2\,\frac{e^2}{r_{12}}\,d\mathbf r_1 d\mathbf r_2,$$
$$K = \int \phi_{1s}^*(\mathbf r_1)\phi_{2s}(\mathbf r_1)\,\frac{e^2}{r_{12}}\,\phi_{2s}^*(\mathbf r_2)\phi_{1s}(\mathbf r_2)\,d\mathbf r_1 d\mathbf r_2.$$

The "+" applies to the singlet (symmetric spatial), the "−" to the triplet (antisymmetric spatial). The first integral $J$ is the classical Coulomb energy; the second integral $K$ is the **exchange integral** and is purely quantum-mechanical. Both $J$ and $K$ are positive.

The triplet, $E_- = J - K$, lies *below* the singlet, $E_+ = J + K$, by $2K \approx 0.8$ eV in helium — an enormous splitting, measurable spectroscopically (the 1s2s $^1S$ and $^3S$ levels of He) and entirely due to exchange. The triplet wins because the antisymmetric spatial wavefunction has a node at $\mathbf r_1 = \mathbf r_2$ (we saw this in §4.5.3), so the two electrons avoid each other and pay less Coulomb energy.

!!! tip "Hund's rule"
    This is the simplest example of **Hund's rule**: for a given electronic configuration, the term with the *highest spin* lies lowest in energy, because parallel-spin electrons are kept apart by antisymmetry and so pay less Coulomb repulsion. The same mechanism makes magnetism possible in transition metals and rare earths.

## 4.5.4 Slater determinants

The simplest antisymmetric wavefunctions are **Slater determinants**. Given $N$ orthonormal one-electron orbitals $\phi_1, \phi_2, \ldots, \phi_N$ (each a function of one spatial coordinate, and for now a single spin), build

$$\Phi(\mathbf r_1, \ldots, \mathbf r_N) = \frac{1}{\sqrt{N!}}\,\det\!\begin{pmatrix} \phi_1(\mathbf r_1) & \phi_2(\mathbf r_1) & \cdots & \phi_N(\mathbf r_1)\\ \phi_1(\mathbf r_2) & \phi_2(\mathbf r_2) & \cdots & \phi_N(\mathbf r_2)\\ \vdots & \vdots & \ddots & \vdots\\ \phi_1(\mathbf r_N) & \phi_2(\mathbf r_N) & \cdots & \phi_N(\mathbf r_N)\end{pmatrix}. \tag{4.5.7}$$

Two structural facts make Slater determinants the workhorse of quantum chemistry.

1. **They are automatically antisymmetric.** Swapping any two electrons $i, j$ swaps two rows of the determinant, which flips its sign — exactly Pauli's rule (4.5.6).

2. **The exclusion principle is built in.** If two orbitals are equal ($\phi_a = \phi_b$), two columns of the determinant are equal and $\det = 0$. You cannot put two electrons (with the same spin) into the same orbital.

Slater determinants form a complete basis for the antisymmetric Hilbert space: the most general $N$-electron wavefunction can be written as a linear combination of all possible Slater determinants built from any complete one-electron basis $\{\chi_i\}$. The number of such determinants for $N$ electrons in $N_b$ basis orbitals is $\binom{N_b}{N}$ (per spin channel), which is still enormous — but it is the *correct* counting of the antisymmetric subspace.

This expansion is called **full configuration interaction** (FCI). With $N = 20$ electrons and $N_b = 100$ basis orbitals, the number of determinants is $\binom{100}{20} \approx 5\times 10^{20}$. FCI scales factorially with system size and is feasible only for molecules with under ~20 electrons in modest basis sets. For materials, FCI is hopeless.

!!! note "What a Slater determinant *cannot* do"
    A single Slater determinant treats electrons as moving in some average potential created by the others (Hartree–Fock, §4.7). It does not capture *correlation*: the fact that an electron at $\mathbf r$ actively *repels* nearby electrons because of $1/|\mathbf r - \mathbf r'|$, regardless of spin. The energy gap between the best single-determinant solution and the true ground state is called the **correlation energy** and is typically 1% of the total energy — but 100% of the chemistry. Bond-making is correlation. Reaction barriers are correlation. Magnetism is correlation. Reproducing correlation is the central challenge of every electronic-structure method.

## 4.5.5 An order-of-magnitude budget

Let us be concrete. Imagine we want to compute the ground-state wavefunction of a single benzene molecule, $\mathrm C_6 \mathrm H_6$. It has 42 electrons (each carbon contributes 6, each hydrogen 1).

- **Real-space grid.** A modest 3D grid with $N_g = 50$ points per direction (1 Å resolution over a 5 Å box per electron, far too coarse for chemistry) gives configuration space of dimension $50^{3\cdot 42} = 50^{126} \approx 10^{214}$.

- **Atomic-orbital basis.** Standard chemistry basis sets (cc-pVDZ) use about 20 contracted basis functions per heavy atom; benzene has roughly $6\cdot 20 + 6\cdot 5 = 150$ basis functions. The number of Slater determinants for 42 electrons in 150 orbitals (per spin) is $\binom{150}{21}^2 \approx 10^{49}$. FCI in this basis is impossible.

- **CCSD(T)**, the "gold standard" of molecular quantum chemistry, scales as $\mathcal O(N_b^7)$. For benzene with $N_b = 150$, this is $150^7 \approx 1.7 \times 10^{15}$ operations — feasible on a workstation. Result: chemical accuracy (1 kcal/mol) for the benzene total energy.

- **DFT** (Chapter 5) scales as $\mathcal O(N_b^3)$. For benzene with $N_b = 150$ this is $\sim 3 \times 10^6$ operations per SCF step — runs in seconds on a laptop. Result: bond lengths correct to about 0.01 Å.

The contrast is the entire point. Solving the Schrödinger equation directly is impossible. Reducing the problem to a well-chosen mean-field theory (Hartree–Fock) or to a functional of the electron density (DFT) makes routine chemistry tractable. The trade-off is that you no longer compute the full wavefunction; you compute an effective single-particle theory and hope that the part you have thrown away — correlation — is small or systematically correctable.

## 4.5.6 Two strategies for escape

The remaining sections of the chapter introduce two complementary approximations that begin to chip away at the exponential wall.

**The Born–Oppenheimer separation** (§4.6) recognises that nuclei are 1836 times heavier than electrons. We can freeze the nuclei, solve the *electronic* Schrödinger equation at fixed $\{\mathbf R_I\}$, and treat the resulting energy as a potential for the nuclei. The wavefunction factorises:

$$\Psi(\mathbf r, \mathbf R) \approx \chi(\mathbf R)\, \psi(\mathbf r; \mathbf R).$$

This does not reduce the electronic problem — that is still $N_{\mathrm e}$-dimensional — but it removes the nuclear coordinates from the worst-scaling part of the calculation. The PES $E_{\mathrm{BO}}(\mathbf R)$ then becomes the central object of all atomistic simulation, including the classical and machine-learning approaches of Chapters 7–9.

**Hartree–Fock** (§4.7) and its descendants tackle the electronic problem itself by restricting the wavefunction to the simplest antisymmetric form — a single Slater determinant — and minimising the energy variationally. This reduces the many-body problem to a nonlinear set of self-consistent one-electron equations, scaling polynomially rather than exponentially. The missing physics (correlation) must then be reintroduced by post-Hartree–Fock methods or absorbed into an exchange-correlation functional in DFT.

The two strategies are independent and composable: every DFT calculation you will run in Chapter 6 uses Born–Oppenheimer to freeze the nuclei *and* an effective single-determinant ansatz for the electronic ground state. Together they convert an impossible $10^{300}$-dimensional eigenvalue problem into a tractable $10^3$-dimensional one, with errors that are typically a fraction of an electron-volt per atom.

## 4.5.6a The hierarchy of electronic-structure methods

It is useful, before closing this section, to sketch the landscape of methods that have been built to evade the exponential wall. Each method is a different compromise between cost and accuracy, and each is appropriate for a different class of problem.

| Method | Scaling | Captures correlation? | Where it lives |
|---|---|---|---|
| Hartree (mean-field, no exchange) | $\mathcal O(N^3)$ | No | Historical curiosity |
| Hartree–Fock | $\mathcal O(N^4)$ | Exchange only | §4.7; basis for post-HF |
| MP2 | $\mathcal O(N^5)$ | Perturbative | Small molecules |
| CCSD | $\mathcal O(N^6)$ | High accuracy | Benchmark chemistry |
| CCSD(T) | $\mathcal O(N^7)$ | "Gold standard" | $\lesssim 30$ atoms |
| FCI | $\mathcal O(N!)$ | Exact in basis | $\lesssim 15$ electrons |
| DFT (KS) | $\mathcal O(N^3)$ | Approximate (functional) | Chapter 5; workhorse |
| QMC (DMC/VMC) | $\mathcal O(N^{3-4})$ | Exact in principle, stochastic | High accuracy on demand |
| ML potentials | $\mathcal O(N)$ | Inherited from training data | Chapter 9; large MD |

Two observations: (1) every accurate method scales at least as the cube of system size, which is why a method that scales as $N^4$ or worse is described as "expensive"; (2) the only method on the list that scales *exactly* as the exponential wall is FCI, the brute-force diagonalisation we have just argued is impossible. Every other method is some clever rearrangement that captures the dominant correlations in polynomial time.

The trade-off is what we have stressed throughout: solving the Schrödinger equation exactly is forbidden, but solving an *approximate* version of it cheaply is the entire science of electronic structure. The exponential wall is the boundary between physics and computational physics.

## 4.5.7 Take-aways

Three messages from this section should stick.

1. **The Schrödinger equation for any non-trivial material is too large to solve.** The wavefunction lives in an exponentially large space, and storing it is physically impossible even in principle.

2. **Pauli antisymmetry both helps and hurts.** It restricts the wavefunction to a smaller subspace (helpful) but forces a particular algebraic structure — Slater determinants — that complicates calculations (hurtful). Antisymmetry is also the origin of exchange, the most important purely quantum interaction in the many-electron problem.

3. **All electronic-structure methods are approximations.** There is no "exact" method one could in principle run on a sufficiently large computer; the cost of an exact method scales exponentially with system size, full stop. Density functional theory, coupled-cluster theory, quantum Monte Carlo, and the modern neural-network ansätze of variational Monte Carlo are all attempts to capture the most important physics in a polynomially-scaling representation.

The rest of this chapter, and indeed the rest of the book, is the long story of those approximations.

!!! tip "How to think about approximations"
    A useful taxonomy: an *uncontrolled* approximation is one for which we have no a priori estimate of the error and no systematic way to improve it (most density functionals fall here); a *controlled* approximation has a known error bound that shrinks as a parameter is tuned (e.g.\ FCI in a finite basis, with basis-set convergence). The most powerful methods combine the two — for instance, **embedding theories** that solve a small, strongly correlated region with a controlled method (CCSD(T) or FCI) and treat the surrounding environment with a cheaper, uncontrolled one (DFT). The art of modern electronic structure is the *composition* of approximations, each correcting the errors of another, with the system size scaling tamed by careful exploitation of locality, sparsity, and renormalisation.

The remainder of Chapter 4 lays the groundwork. In §4.6 we make the first great simplification — the Born–Oppenheimer separation — which decouples electrons from nuclei. In §4.7 we take the first serious crack at the resulting electronic problem with Hartree–Fock. Then in Chapter 5 we replace the explicit wavefunction with an electron density, the conceptual leap that opened the door to *practical* materials simulation.
