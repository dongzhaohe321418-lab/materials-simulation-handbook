# 4.8 Exercises

Difficulty markers: $\star$ = warm-up, $\star\star$ = standard, $\star\star\star$ = harder. Full solutions are given in collapsed admonitions. Try each exercise *before* opening the solution.

!!! info "What problem are we solving?"
    These exercises check that you can *reproduce*, not just recognise, the
    central results of Chapter 4: that the textbook wavefunctions really are
    normalised and orthogonal, that the harmonic-oscillator and
    particle-in-a-box energies come out where the chapter said, and that the
    "exponential wall" of many-electron quantum mechanics is a countable,
    concrete number rather than a slogan. Every formula you need is restated
    inside the problem, so you can work each one without flipping back through
    the chapter.

!!! note "How these exercises are graded"
    The original difficulty stars ($\star$ to $\star\star\star$) are kept. On
    top of them, each problem is tagged with a **level** that says what kind of
    thinking it asks for, so you can build up gradually:

    - **Level A — recall.** Reproduce a definition or a one-line result.
    - **Level B — single calculation.** Turn one formula into one number or one
      function.
    - **Level C — multi-step.** Chain two or three steps (an integral, then a
      substitution, then a limit).
    - **Level D — apply / code.** Run or modify a calculation and read the result.
    - **Level E — critique / estimate.** Judge an approximation, estimate an
      order of magnitude, or explain *why* a method scales the way it does.

    A rough order of attack for a first pass: do all the Level A and B problems
    first (4.1, 4.2, 4.3, the recall parts of 4.4), then the Level C ones (4.4,
    4.5), then the Level D/E ones (4.6, 4.7, 4.8).

!!! tip "New vocabulary"
    - **Normalisation** — scaling a wavefunction so that $\int|\psi|^2\,dx = 1$,
      i.e. the particle is somewhere with total probability one. See the
      [beginner glossary](../undergraduate/glossary-for-beginners.md) entry on
      *Wavefunction*.
    - **Orthogonality** — two states with $\langle\phi_m|\phi_n\rangle = 0$; they
      have "no overlap" and represent independent measurement outcomes.
    - **Eigenstate / eigenvalue** — a state the operator merely rescales,
      $\hat A\psi = a\psi$; the number $a$ is the measurable value. See the
      glossary entries on *Eigenvalue* and *Eigenvector*.
    - **Expectation value** — the probability-weighted average of a quantity in
      a given state, $\langle x\rangle = \int \psi^* x\,\psi\,dx$.

The symbols below recur across several exercises; this table collects them once
so none appears undefined.

| Symbol | Meaning | Units (SI) |
|---|---|---|
| $L$ | width of the 1-D box | m |
| $n$ | quantum number (state label), $n = 1, 2, \dots$ (box) or $0, 1, 2, \dots$ (oscillator) | — |
| $m$, $m_e$ | particle mass; electron mass $9.109\times10^{-31}$ kg | kg |
| $\hbar$ | reduced Planck constant, $1.055\times10^{-34}$ | J s |
| $e$ | elementary charge, $1.602\times10^{-19}$ | C (1 eV $= e$ joule) |
| $\omega$ | angular frequency of the oscillator | rad s$^{-1}$ |
| $\ell$ | oscillator length $\sqrt{\hbar/(m\omega)}$ | m |
| $\xi$ | dimensionless oscillator coordinate $x/\ell$ | — |
| $H_n(\xi)$ | $n$th Hermite polynomial | — |
| $\hat A$, $\hat H$ | a Hermitian operator; the Hamiltonian (energy operator) | varies; J |
| $N_g$, $N$ | number of grid points; number of electrons | — |

---

## Exercise 4.1 — Normalisation of particle-in-a-box states ($\star$, Level B)

Verify directly that the analytical eigenfunctions

$$\psi_n(x) = \sqrt{\frac{2}{L}}\,\sin\!\left(\frac{n\pi x}{L}\right), \qquad n = 1, 2, 3, \ldots,$$

are normalised on $[0, L]$. *Self-contained statement:* "normalised on $[0,L]$"
means $\int_0^L |\psi_n(x)|^2\,dx = 1$, i.e. the probability of finding the
particle *somewhere* in the box is exactly one. You will need the double-angle
identity $\sin^2\theta = \tfrac12 - \tfrac12\cos 2\theta$.

??? note "Hint"
    Pull the constant $2/L$ outside the integral, replace $\sin^2$ by the
    double-angle identity, and integrate term by term. The cosine term
    integrates to a sine that vanishes at both ends because its argument is a
    whole number of half-periods. You are left with $\tfrac12 \times L$.

??? success "Solution"
    Compute
    $$\int_0^L |\psi_n|^2\, dx = \frac{2}{L}\int_0^L \sin^2\!\left(\frac{n\pi x}{L}\right) dx.$$
    Using the identity $\sin^2\theta = \tfrac12 - \tfrac12 \cos 2\theta$,
    $$\int_0^L \sin^2\!\left(\frac{n\pi x}{L}\right) dx = \int_0^L \!\left[\tfrac12 - \tfrac12 \cos\!\left(\tfrac{2n\pi x}{L}\right)\right] dx = \tfrac{L}{2} - \tfrac12 \cdot \frac{L}{2n\pi}\sin\!\left(\tfrac{2n\pi x}{L}\right)\Big|_0^L = \tfrac{L}{2},$$
    since $\sin(2n\pi) = \sin(0) = 0$ for any integer $n$. Therefore $\int_0^L |\psi_n|^2 dx = (2/L)(L/2) = 1$. $\blacksquare$

---

## Exercise 4.2 — Orthogonality of non-degenerate eigenstates ($\star$, Level C)

Prove from first principles that two eigenfunctions of a Hermitian operator $\hat A$ belonging to *different* eigenvalues are orthogonal. (This was proved as a theorem in §4.2.6; rewrite the proof in your own words, without consulting the text.)

*Self-contained statement:* you may use two facts established earlier in the
chapter. (i) **Hermiticity** means $\langle \phi_m | \hat A \phi_n\rangle =
\langle \hat A \phi_m | \phi_n\rangle$ for all states — the operator can act
"to the left" or "to the right" inside an inner product with the same result.
(ii) The eigenvalues of a Hermitian operator are **real**, so $a_m^* = a_m$.
"Orthogonal" means $\langle \phi_m | \phi_n\rangle = 0$.

??? note "Hint"
    Evaluate the *same* number $\langle \phi_m | \hat A \phi_n\rangle$ two ways:
    once letting $\hat A$ act on $\phi_n$ (giving $a_n$), once using Hermiticity
    to let it act on $\phi_m$ (giving $a_m$). Subtract. You get
    $(a_n - a_m)\langle \phi_m | \phi_n\rangle = 0$; now use $a_m \ne a_n$.

??? success "Solution"
    Let $\hat A\phi_m = a_m\phi_m$ and $\hat A\phi_n = a_n\phi_n$ with $a_m \neq a_n$. Compute
    $$\langle \phi_m | \hat A \phi_n\rangle = a_n \langle \phi_m | \phi_n\rangle.$$
    By Hermiticity,
    $$\langle \phi_m | \hat A \phi_n\rangle = \langle \hat A \phi_m | \phi_n\rangle = a_m^* \langle \phi_m | \phi_n\rangle = a_m \langle \phi_m | \phi_n\rangle,$$
    where we used the fact (also proved in §4.2.5) that Hermitian eigenvalues are real, so $a_m^* = a_m$. Equating the two expressions,
    $$(a_n - a_m)\langle \phi_m | \phi_n\rangle = 0,$$
    and since $a_m \neq a_n$, the inner product must vanish: $\langle \phi_m | \phi_n\rangle = 0$. $\blacksquare$

---

## Exercise 4.3 — Ground-state energy of an electron in a 1 Å box ($\star$, Level B)

Estimate the ground-state energy of an electron confined to a 1D box of width $L = 1\ \mathrm{\AA} = 10^{-10}$ m. Give the answer in eV.

*Self-contained statement:* the particle-in-a-box energy levels (derived in the
chapter as equation (4.3.7)) are

$$E_n = \frac{n^2\pi^2\hbar^2}{2 m L^2}, \qquad n = 1, 2, 3, \ldots \tag{4.8.1}$$

so the ground state ($n=1$) is $E_1 = \pi^2\hbar^2/(2 m_e L^2)$. Use
$\hbar = 1.055\times10^{-34}\ \mathrm{J\,s}$,
$m_e = 9.109\times10^{-31}\ \mathrm{kg}$, and convert joules to electronvolts by
dividing by $e = 1.602\times10^{-19}$ (since $1\ \mathrm{eV} = e$ joule).

??? note "Hint"
    Square $\hbar$ and $L$ first, keep track of the powers of ten, and only
    convert to eV at the very end. A useful sanity check: the answer should be
    *tens* of eV, comparable to atomic binding energies, not the milli-eV of
    thermal motion.

??? success "Solution"
    From (4.3.7), $E_1 = \pi^2 \hbar^2/(2 m_e L^2)$. Plug in:
    $$E_1 = \frac{\pi^2 \cdot (1.055\times 10^{-34})^2}{2 \cdot 9.109\times 10^{-31} \cdot (10^{-10})^2}\ \mathrm{J} = 6.025\times 10^{-18}\ \mathrm{J}.$$
    Divide by $e = 1.602\times 10^{-19}$ C to convert to eV:
    $$E_1 \approx 37.6\ \mathrm{eV}.$$
    For comparison, the ionisation energy of hydrogen is 13.6 eV; the particle-in-a-box ground state at $L = 1$ Å is in the same ballpark, which is why the box is a halfway-reasonable cartoon of an atom.

---

## Exercise 4.4 — First three Hermite polynomials from the recursion ($\star\star$, Level C)

Starting from $H_0(\xi) = 1$ and the recursion

$$H_{n+1}(\xi) = 2\xi H_n(\xi) - 2n H_{n-1}(\xi), \tag{4.8.2}$$

derive $H_1(\xi)$, $H_2(\xi)$, and $H_3(\xi)$ explicitly. Verify that $H_2(\xi)$ has two real roots and identify their locations.

*Self-contained statement:* the Hermite polynomials $H_n(\xi)$ are the
polynomial factors in the harmonic-oscillator eigenfunctions
$\psi_n(x) \propto H_n(\xi)\,e^{-\xi^2/2}$ with $\xi = x/\ell$ (Section 4.4).
The recursion (4.8.2) generates each one from the two below it. The roots of
$H_n$ are the **nodes** (zero crossings) of the $n$th excited state, so $H_2$
having two roots says $\psi_2$ crosses zero twice.

??? note "Hint 1 — getting started"
    The recursion needs *two* lower polynomials. You are given $H_0 = 1$. To use
    (4.8.2) for $n = 0$ you also need $H_{-1}$: take the standard convention
    $H_{-1} \equiv 0$ (its coefficient $2n = 0$ at $n=0$ kills it anyway). Then
    the $n=0$ case gives $H_1$ directly.

??? note "Hint 2 — climbing the ladder"
    Substitute $n = 0$ into (4.8.2) to get $H_1$, then $n = 1$ to get $H_2$
    (it needs $H_1$ and $H_0$), then $n = 2$ to get $H_3$ (it needs $H_2$ and
    $H_1$). Expand the brackets fully and collect like powers of $\xi$ at each
    step — do not leave $2\xi(4\xi^2 - 2)$ unmultiplied.

??? note "Hint 3 — the roots of $H_2$"
    Set your expression for $H_2$ equal to zero. It is a quadratic in $\xi$ with
    no linear term, so solve $4\xi^2 - 2 = 0$ for $\xi^2$, then take both square
    roots. Two real roots, symmetric about the origin, confirm the two nodes.

??? success "Solution"
    Apply the recursion with $H_{-1} \equiv 0$ and $H_0 = 1$.
    - $H_1 = 2\xi \cdot H_0 - 2\cdot 0 \cdot H_{-1} = 2\xi$.
    - $H_2 = 2\xi \cdot H_1 - 2\cdot 1 \cdot H_0 = 2\xi(2\xi) - 2 = 4\xi^2 - 2$.
    - $H_3 = 2\xi \cdot H_2 - 2\cdot 2 \cdot H_1 = 2\xi(4\xi^2 - 2) - 4(2\xi) = 8\xi^3 - 4\xi - 8\xi = 8\xi^3 - 12\xi$.
    Solving $H_2(\xi) = 4\xi^2 - 2 = 0$ gives $\xi^2 = 1/2$, hence $\xi = \pm 1/\sqrt{2}$. These are the two nodes of $\psi_2(x)$ — the wavefunction of the second excited state of the harmonic oscillator, which has two zero crossings, as expected.

---

## Exercise 4.5 — Expectation value of position for an SHO ground state ($\star\star$, Level C)

Compute $\langle x\rangle$ and $\langle x^2\rangle$ in the ground state $\psi_0$ of the 1D harmonic oscillator, and use them to evaluate $\Delta x = \sqrt{\langle x^2\rangle - \langle x\rangle^2}$. Express your answer in units of the oscillator length $\ell = \sqrt{\hbar/(m\omega)}$.

*Self-contained statement:* the harmonic-oscillator ground state (equation
(4.4.9), written here in terms of $\ell$) is the normalised Gaussian

$$\psi_0(x) = \frac{1}{(\pi\ell^2)^{1/4}}\,\exp\!\left(-\frac{x^2}{2\ell^2}\right). \tag{4.8.3}$$

An expectation value is the probability-weighted average
$\langle f\rangle = \int_{-\infty}^{\infty}\psi_0^*\,f\,\psi_0\,dx
= \int_{-\infty}^{\infty} f(x)\,|\psi_0(x)|^2\,dx$, and $\Delta x$ is the
root-mean-square spread about that average. You may use the standard Gaussian
integral

$$\int_{-\infty}^{\infty} x^2\,e^{-x^2/\ell^2}\,dx = \frac{\ell^3\sqrt{\pi}}{2}. \tag{4.8.4}$$

??? note "Hint"
    Do not integrate $\langle x\rangle$ by brute force: $|\psi_0|^2$ is an
    *even* function and $x$ is *odd*, so the integrand is odd and the integral
    over the symmetric range $(-\infty,\infty)$ is zero by symmetry — hence
    $\langle x\rangle = 0$. For $\langle x^2\rangle$, note
    $|\psi_0|^2 = (\pi\ell^2)^{-1/2} e^{-x^2/\ell^2}$ and plug straight into
    (4.8.4). With $\langle x\rangle = 0$, $\Delta x = \sqrt{\langle x^2\rangle}$.

??? success "Solution"
    The ground state from (4.4.9) is
    $$\psi_0(x) = \left(\frac{m\omega}{\pi\hbar}\right)^{1/4}\exp\!\left(-\frac{m\omega x^2}{2\hbar}\right) = \frac{1}{(\pi\ell^2)^{1/4}}\exp\!\left(-\frac{x^2}{2\ell^2}\right).$$
    Since $|\psi_0|^2$ is even, $\langle x\rangle = 0$ by symmetry. For $\langle x^2\rangle$ use the standard Gaussian integral
    $$\int_{-\infty}^\infty x^2 e^{-x^2/\ell^2}\, dx = \frac{\ell^3 \sqrt\pi}{2}.$$
    Then
    $$\langle x^2\rangle = \frac{1}{\ell \sqrt\pi}\int_{-\infty}^\infty x^2 e^{-x^2/\ell^2}\, dx = \frac{1}{\ell\sqrt\pi}\cdot \frac{\ell^3 \sqrt\pi}{2} = \frac{\ell^2}{2}.$$
    Hence $\Delta x = \ell/\sqrt{2}$. By an analogous calculation $\Delta p = \hbar/(\sqrt{2}\ell)$, and the product saturates the uncertainty bound: $\Delta x\, \Delta p = \hbar/2$. The harmonic-oscillator ground state is a *minimum-uncertainty* state.

---

## Exercise 4.6 — Modifying the FD code for a double well ($\star\star$, Level D)

Take the finite-difference code of §4.3 and modify it to solve the double-well potential

$$V(x) = V_0 \left[\left(\frac{x}{a}\right)^2 - 1\right]^2, \tag{4.8.5}$$

with $V_0 = 1.0$ eV and $a = 0.5$ nm. Plot the lowest four eigenstates and their probability densities, and explain the pattern of (near-)degeneracies you see in terms of tunnelling between the two wells.

*Self-contained statement:* the finite-difference method (Section 4.3)
discretises $-\tfrac{\hbar^2}{2m}\psi'' + V(x)\psi = E\psi$ on a grid of spacing
$h$, turning the Hamiltonian into a tridiagonal matrix: $2\,\mathrm{pre} + V_i$
on the diagonal and $-\mathrm{pre}$ on the two neighbours, with
$\mathrm{pre} = \hbar^2/(2 m h^2)$. Diagonalising it (`np.linalg.eigh`) returns
the eigenvalues (energies) and eigenvectors (wavefunctions). The potential
(4.8.5) is a "W" shape: two minima at $x = \pm a$ separated by a hump of height
$V_0$ at $x = 0$. The runtime is plain NumPy/Matplotlib, so this one *does* run
live in the browser.

??? note "Hint"
    You only need to change the *potential*, not the solver. Replace the box or
    oscillator $V$ with `V = V0 * ((x/a)**2 - 1)**2` (with `V0` and `a` in SI
    units), keep a grid wide enough to contain both wells with room to spare
    (a few nm), and diagonalise as before. After solving, look at the spacing of
    the first four eigenvalues: you should find them in two close pairs, not
    four evenly spaced levels.

??? success "Solution"
    A complete script:
    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    HBAR = 1.054_571_817e-34
    M_E  = 9.109_383_7e-31
    EV   = 1.602_176_634e-19

    def build_H(x, V, mass=M_E):
        h = x[1] - x[0]
        pre = HBAR**2 / (2 * mass * h**2)
        n = x.size
        H = (np.diag(2*pre*np.ones(n) + V)
             + np.diag(-pre*np.ones(n-1), k=1)
             + np.diag(-pre*np.ones(n-1), k=-1))
        return H

    L = 3.0e-9
    n_grid = 1200
    x = np.linspace(-L/2, L/2, n_grid)

    V0 = 1.0 * EV
    a  = 0.5e-9
    V  = V0 * ((x/a)**2 - 1)**2

    H = build_H(x, V)
    eigvals, eigvecs = np.linalg.eigh(H)
    eigvecs /= np.sqrt(x[1] - x[0])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x*1e9, V/EV, "k", lw=1.5)
    scale = 0.15
    for n in range(4):
        ax.plot(x*1e9, eigvals[n]/EV + scale*eigvecs[:, n]/np.max(np.abs(eigvecs[:, n])),
                label=f"n={n}, E={eigvals[n]/EV:.4f} eV")
        ax.axhline(eigvals[n]/EV, color="gray", ls=":", lw=0.5)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("E (eV)")
    ax.set_ylim(-0.05, 1.5)
    ax.legend()
    plt.tight_layout()
    plt.show()
    ```
    **What you should see.** The lowest two eigenvalues come in a closely-spaced *pair* (symmetric and antisymmetric combinations of "left-well" and "right-well" states), and similarly the third and fourth. The splitting within each pair is the **tunnel splitting**: even though the two wells are separated by a barrier of height $V_0$, the wavefunction has a small but non-zero amplitude in the classically forbidden region, allowing the particle to be in a coherent superposition of "in the left well" and "in the right well". The smaller the barrier or the lighter the particle, the larger the tunnel splitting. This is the prototype of every quantum-mechanical double well — ammonia inversion, hydrogen bonds, ferroelectric phase transitions.

---

## Exercise 4.7 — Energies via Hellmann–Feynman ($\star\star$, Level C)

For a normalised energy eigenstate $|\psi_\lambda\rangle$ of a Hamiltonian $\hat{H}(\lambda)$ depending on a parameter $\lambda$, prove the Hellmann–Feynman theorem:

$$\frac{dE(\lambda)}{d\lambda} = \langle\psi_\lambda | \frac{\partial \hat{H}}{\partial \lambda} | \psi_\lambda\rangle.$$

*Self-contained statement:* "normalised" means $\langle\psi_\lambda|\psi_\lambda\rangle = 1$
for *every* $\lambda$, and the eigenvalue equation is
$\hat H(\lambda)|\psi_\lambda\rangle = E(\lambda)|\psi_\lambda\rangle$. The
notation $\partial_\lambda \equiv \partial/\partial\lambda$. The point of the
theorem is that to get $dE/d\lambda$ you only differentiate the *operator*, never
the (hard-to-find) wavefunction.

??? note "Hint"
    Start from $E(\lambda) = \langle\psi_\lambda|\hat H|\psi_\lambda\rangle$ and
    differentiate using the product rule — three terms appear, because
    $\langle\psi_\lambda|$, $\hat H$, and $|\psi_\lambda\rangle$ each depend on
    $\lambda$. In the term where $\hat H$ hits $|\partial_\lambda\psi_\lambda\rangle$,
    use the eigenvalue equation; in the term where it hits
    $\langle\partial_\lambda\psi_\lambda|$, use Hermiticity to move it onto
    $|\psi_\lambda\rangle$. The two wavefunction-derivative terms then combine
    into $E\,\partial_\lambda\langle\psi_\lambda|\psi_\lambda\rangle =
    E\,\partial_\lambda(1) = 0$.

??? success "Solution"
    By assumption, $\hat{H}(\lambda)|\psi_\lambda\rangle = E(\lambda)|\psi_\lambda\rangle$ with $\langle\psi_\lambda|\psi_\lambda\rangle = 1$. Take the inner product with $\langle\psi_\lambda|$:
    $$E(\lambda) = \langle\psi_\lambda | \hat{H} | \psi_\lambda\rangle.$$
    Differentiate with respect to $\lambda$:
    $$\frac{dE}{d\lambda} = \langle \partial_\lambda \psi_\lambda | \hat{H} | \psi_\lambda\rangle + \langle\psi_\lambda | \partial_\lambda \hat{H} | \psi_\lambda\rangle + \langle\psi_\lambda | \hat{H} | \partial_\lambda \psi_\lambda\rangle.$$
    Use the eigenvalue equation in the first and third terms, and the Hermiticity of $\hat{H}$:
    $$\langle \partial_\lambda \psi_\lambda | \hat{H} | \psi_\lambda\rangle + \langle\psi_\lambda | \hat{H} | \partial_\lambda \psi_\lambda\rangle = E(\lambda)\bigl[\langle \partial_\lambda \psi_\lambda | \psi_\lambda\rangle + \langle\psi_\lambda | \partial_\lambda \psi_\lambda\rangle\bigr] = E(\lambda)\, \partial_\lambda \langle\psi_\lambda|\psi_\lambda\rangle = E(\lambda)\cdot 0 = 0,$$
    using normalisation. Hence
    $$\frac{dE}{d\lambda} = \langle\psi_\lambda | \partial_\lambda \hat{H} | \psi_\lambda\rangle. \quad\blacksquare$$
    **Application.** Taking $\lambda \to \mathbf R_I$ (a nuclear position) gives the force formula (4.6.14) that we used in the BO section. The force on a nucleus is the expectation value of the gradient of the Hamiltonian — no need to differentiate the wavefunction, which is the gift that makes ab-initio MD feasible.

---

## Exercise 4.8 — Estimating the size of a many-electron Hilbert space ($\star\star\star$, Level E)

Consider $N$ spinless electrons confined to a 1D box of width $L$, discretised on a grid of $N_g$ points. How many *antisymmetric* basis states are there? Evaluate this number for $N = 20$ electrons and $N_g = 100$ grid points, and compare with $N_g^N$ (the dimension if antisymmetry were not enforced).

*Self-contained statement:* a many-electron wavefunction must be **antisymmetric**
under swapping any two electrons (the Pauli principle, Section 4.5); for spinless
electrons this forbids two electrons sitting on the same grid point. The number
of ways to choose $N$ distinct points out of $N_g$, order not mattering, is the
**binomial coefficient**

$$\binom{N_g}{N} = \frac{N_g!}{N!\,(N_g - N)!}. \tag{4.8.6}$$

Without antisymmetry, each of the $N$ electrons could independently sit on any of
the $N_g$ points, giving $N_g^N$ arrangements.

??? note "Hint 1 — why a binomial coefficient"
    Antisymmetry (a Slater determinant) is built from $N$ *distinct*
    single-particle states, and swapping two of them only flips the overall
    sign — it is not a new physical state. So you are *choosing* an unordered
    set of $N$ occupied grid points out of $N_g$. "Choose $N$ from $N_g$,
    order irrelevant" is exactly $\binom{N_g}{N}$.

??? note "Hint 2 — evaluating $\binom{100}{20}$ without overflow"
    Do not compute $100!$ directly — it overflows ordinary floating point.
    Instead build the product incrementally,
    $\binom{N_g}{N} = \prod_{i=0}^{N-1}\frac{N_g - i}{i + 1}$, multiplying and
    dividing as you go so the running value stays moderate. In Python,
    `math.comb(100, 20)` gives the exact integer
    ($\approx 5.36\times10^{20}$).

??? note "Hint 3 — the comparison"
    Compare $\binom{100}{20}\approx 5.4\times10^{20}$ with
    $N_g^N = 100^{20} = 10^{40}$. The ratio is about $1.9\times10^{19}$, which
    is of order $N! = 20! \approx 2.4\times10^{18}$ (a few times larger). So
    enforcing antisymmetry removes roughly a factor of $N!$ — a huge saving in
    relative terms, yet the surviving count is *still* astronomically large.
    That is the quantitative face of the "exponential wall".

??? success "Solution"
    Place a single electron on the grid in one of $N_g$ ways. With $N$ electrons and antisymmetry, two electrons cannot share the same grid point (Pauli, since we are taking them to be spinless or all of one spin). The number of antisymmetric basis states is therefore the number of ways to *choose* $N$ distinct grid points out of $N_g$ — independent of order, because the determinantal structure fixes the sign:
    $$\dim \mathcal H_{\mathrm{antisym}} = \binom{N_g}{N}.$$
    For $N_g = 100$, $N = 20$:
    $$\binom{100}{20} = \frac{100!}{20!\,80!} \approx 5.36\times 10^{20}.$$
    Without antisymmetry, the unrestricted Hilbert space would have $N_g^N = 100^{20} = 10^{40}$ basis states.
    Antisymmetry reduces the count by a factor of $N! \cdot (\text{double-occupancy combinatorics}) \sim N!$, but the resulting number is *still* astronomical — 20 electrons on a 100-point grid require half a sextillion basis states. This is a quantitative version of the exponential wall of §4.5, and explains why FCI is impossible for any system with more than ~20 electrons in modest basis sets.
    **Take-away.** Pauli antisymmetry helps — it removes a factor of $\sim N!$ from the basis count — but it does not save us from exponential scaling. Approximate methods (HF, DFT, MP2, CCSD, …) that scale polynomially in $N$ are the only way forward.

---

## Check yourself

Before moving on, make sure you can answer these without looking back at the
solutions. They test the threads running through the whole exercise set.

!!! warning "Common misunderstandings"
    - **Normalisation is about $|\psi|^2$, not $\psi$.** In Exercise 4.1 the
      condition is $\int|\psi_n|^2\,dx = 1$. The wavefunction $\psi$ itself can
      be negative; only its squared modulus is a probability density.
    - **"Orthogonal" needs *different* eigenvalues.** The proof in Exercise 4.2
      relies on $a_m \ne a_n$. For *degenerate* states (equal eigenvalues) the
      theorem says nothing — you have to orthogonalise them by hand.
    - **A node count is not an energy.** In Exercise 4.4, $H_2$ having two roots
      tells you $\psi_2$ has two nodes; it does not directly give the energy.
      More nodes does correlate with higher energy, but the energy comes from
      the eigenvalue, not from counting zeros.
    - **Antisymmetry helps but does not rescue you.** In Exercise 4.8 the
      Pauli principle shrinks the basis by a factor of order $N!$, yet
      $\binom{100}{20}$ is still $\sim10^{20}$. The wall is *exponential*; a
      factorial prefactor does not defeat it.

!!! question "Check yourself"
    1. In Exercise 4.1, where does the factor $\sqrt{2/L}$ in
       $\psi_n = \sqrt{2/L}\,\sin(n\pi x/L)$ come from, and what would the
       integral $\int_0^L|\psi_n|^2\,dx$ equal if you forgot it?
    2. State the two facts about a Hermitian operator that the orthogonality
       proof of Exercise 4.2 relies on.
    3. The ground-state energy of an electron in a $1\,\mathrm{\AA}$ box
       (Exercise 4.3) came out near $38\,\mathrm{eV}$. If you doubled the box
       width to $2\,\mathrm{\AA}$, by what factor would $E_1$ change?
    4. Using the recursion (4.8.2), how many real roots does $H_3(\xi)$ have,
       and what does that tell you about $\psi_3$?
    5. For the oscillator ground state (Exercise 4.5), why is $\langle x\rangle$
       zero without doing any integral?
    6. In Exercise 4.8, why is the antisymmetric count a *binomial* coefficient
       $\binom{N_g}{N}$ rather than $N_g^{\,N}$?

    ??? success "Answer"
        1. It is the **normalisation constant**: $\int_0^L \sin^2(n\pi x/L)\,dx
           = L/2$, so to make $\int_0^L|\psi_n|^2\,dx = 1$ we need
           $(2/L)(L/2) = 1$. Without the factor, the integral would equal
           $L/2$, not $1$ — the state would not be normalised.
        2. (i) Hermiticity, $\langle\phi_m|\hat A\phi_n\rangle =
           \langle\hat A\phi_m|\phi_n\rangle$; and (ii) the eigenvalues are
           real, $a_m^* = a_m$.
        3. From $E_1 = \pi^2\hbar^2/(2 m_e L^2)\propto 1/L^2$, doubling $L$
           divides $E_1$ by $2^2 = 4$, giving roughly $9.4\,\mathrm{eV}$. Wider
           box, lower confinement energy.
        4. $H_3(\xi) = 8\xi^3 - 12\xi = 4\xi(2\xi^2 - 3)$ has three real roots,
           at $\xi = 0$ and $\xi = \pm\sqrt{3/2}$. So $\psi_3$ has three nodes —
           the standard "$n$ nodes in the $n$th excited state" pattern.
        5. Because $|\psi_0|^2$ is **even** in $x$ while $x$ is **odd**, so the
           integrand $x\,|\psi_0|^2$ is odd and integrates to zero over the
           symmetric range $(-\infty,\infty)$. The state is symmetric about the
           origin, so its average position is the centre.
        6. Antisymmetry (a Slater determinant) is built from $N$ *distinct*
           occupied points, and reordering them only flips a sign rather than
           giving a new state. So you *choose* an unordered set of $N$ points
           from $N_g$ — exactly $\binom{N_g}{N}$. The $N_g^{\,N}$ count would
           treat ordered, possibly repeated occupations as distinct, which
           antisymmetry forbids.

## Going further

If you want more practice, try:

- Solve numerically for an electron in a finite square well of depth $V_0 = 5$ eV and width $L = 1$ nm. Count the number of bound states and compare with the analytical transcendental-equation result.
- Repeat the SHO calculation in *atomic units* ($\hbar = m_e = 1$). Confirm that $E_n = \omega(n + 1/2)$ with $\omega = 1$ gives eigenvalues $0.5, 1.5, 2.5, \ldots$ hartree.
- For the double-well potential of Exercise 4.6, plot the tunnel splitting $\Delta E_{01} = E_1 - E_0$ as a function of barrier height $V_0$. Verify that it decays exponentially with $\sqrt{V_0}$ — the WKB prediction.
- Compute the Hartree–Fock ground-state energy of the helium atom analytically using a single-zeta Slater orbital $\phi(\mathbf r) = (\zeta^3/\pi)^{1/2} e^{-\zeta r}$, minimising over $\zeta$. Compare with the experimental value (-79.0 eV) and identify the correlation energy.

In Chapter 5 we will leave wavefunctions behind and reformulate everything in terms of the electron density. The mathematical machinery developed in this chapter — Hermitian operators, eigenvalue problems, finite-difference discretisations, SCF cycles — will all reappear, but in a form that can finally be applied to real materials.
