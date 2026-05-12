# 4.8 Exercises

Difficulty markers: $\star$ = warm-up, $\star\star$ = standard, $\star\star\star$ = harder. Full solutions are given in collapsed admonitions. Try each exercise *before* opening the solution.

---

## Exercise 4.1 — Normalisation of particle-in-a-box states ($\star$)

Verify directly that the analytical eigenfunctions

$$\psi_n(x) = \sqrt{\frac{2}{L}}\,\sin\!\left(\frac{n\pi x}{L}\right), \qquad n = 1, 2, 3, \ldots,$$

are normalised on $[0, L]$.

??? success "Solution"
    Compute
    $$\int_0^L |\psi_n|^2\, dx = \frac{2}{L}\int_0^L \sin^2\!\left(\frac{n\pi x}{L}\right) dx.$$
    Using the identity $\sin^2\theta = \tfrac12 - \tfrac12 \cos 2\theta$,
    $$\int_0^L \sin^2\!\left(\frac{n\pi x}{L}\right) dx = \int_0^L \!\left[\tfrac12 - \tfrac12 \cos\!\left(\tfrac{2n\pi x}{L}\right)\right] dx = \tfrac{L}{2} - \tfrac12 \cdot \frac{L}{2n\pi}\sin\!\left(\tfrac{2n\pi x}{L}\right)\Big|_0^L = \tfrac{L}{2},$$
    since $\sin(2n\pi) = \sin(0) = 0$ for any integer $n$. Therefore $\int_0^L |\psi_n|^2 dx = (2/L)(L/2) = 1$. $\blacksquare$

---

## Exercise 4.2 — Orthogonality of non-degenerate eigenstates ($\star$)

Prove from first principles that two eigenfunctions of a Hermitian operator $\hat A$ belonging to *different* eigenvalues are orthogonal. (This was proved as a theorem in §4.2.6; rewrite the proof in your own words, without consulting the text.)

??? success "Solution"
    Let $\hat A\phi_m = a_m\phi_m$ and $\hat A\phi_n = a_n\phi_n$ with $a_m \neq a_n$. Compute
    $$\langle \phi_m | \hat A \phi_n\rangle = a_n \langle \phi_m | \phi_n\rangle.$$
    By Hermiticity,
    $$\langle \phi_m | \hat A \phi_n\rangle = \langle \hat A \phi_m | \phi_n\rangle = a_m^* \langle \phi_m | \phi_n\rangle = a_m \langle \phi_m | \phi_n\rangle,$$
    where we used the fact (also proved in §4.2.5) that Hermitian eigenvalues are real, so $a_m^* = a_m$. Equating the two expressions,
    $$(a_n - a_m)\langle \phi_m | \phi_n\rangle = 0,$$
    and since $a_m \neq a_n$, the inner product must vanish: $\langle \phi_m | \phi_n\rangle = 0$. $\blacksquare$

---

## Exercise 4.3 — Ground-state energy of an electron in a 1 Å box ($\star$)

Estimate the ground-state energy of an electron confined to a 1D box of width $L = 1\ \mathrm{\AA} = 10^{-10}$ m. Give the answer in eV.

??? success "Solution"
    From (4.3.7), $E_1 = \pi^2 \hbar^2/(2 m_e L^2)$. Plug in:
    $$E_1 = \frac{\pi^2 \cdot (1.055\times 10^{-34})^2}{2 \cdot 9.109\times 10^{-31} \cdot (10^{-10})^2}\ \mathrm{J} = 6.025\times 10^{-18}\ \mathrm{J}.$$
    Divide by $e = 1.602\times 10^{-19}$ C to convert to eV:
    $$E_1 \approx 37.6\ \mathrm{eV}.$$
    For comparison, the ionisation energy of hydrogen is 13.6 eV; the particle-in-a-box ground state at $L = 1$ Å is in the same ballpark, which is why the box is a halfway-reasonable cartoon of an atom.

---

## Exercise 4.4 — First three Hermite polynomials from the recursion ($\star\star$)

Starting from $H_0(\xi) = 1$ and the recursion

$$H_{n+1}(\xi) = 2\xi H_n(\xi) - 2n H_{n-1}(\xi),$$

derive $H_1(\xi)$, $H_2(\xi)$, and $H_3(\xi)$ explicitly. Verify that $H_2(\xi)$ has two real roots and identify their locations.

??? success "Solution"
    Apply the recursion with $H_{-1} \equiv 0$ and $H_0 = 1$.
    - $H_1 = 2\xi \cdot H_0 - 2\cdot 0 \cdot H_{-1} = 2\xi$.
    - $H_2 = 2\xi \cdot H_1 - 2\cdot 1 \cdot H_0 = 2\xi(2\xi) - 2 = 4\xi^2 - 2$.
    - $H_3 = 2\xi \cdot H_2 - 2\cdot 2 \cdot H_1 = 2\xi(4\xi^2 - 2) - 4(2\xi) = 8\xi^3 - 4\xi - 8\xi = 8\xi^3 - 12\xi$.
    Solving $H_2(\xi) = 4\xi^2 - 2 = 0$ gives $\xi^2 = 1/2$, hence $\xi = \pm 1/\sqrt{2}$. These are the two nodes of $\psi_2(x)$ — the wavefunction of the second excited state of the harmonic oscillator, which has two zero crossings, as expected.

---

## Exercise 4.5 — Expectation value of position for an SHO ground state ($\star\star$)

Compute $\langle x\rangle$ and $\langle x^2\rangle$ in the ground state $\psi_0$ of the 1D harmonic oscillator, and use them to evaluate $\Delta x = \sqrt{\langle x^2\rangle - \langle x\rangle^2}$. Express your answer in units of the oscillator length $\ell = \sqrt{\hbar/(m\omega)}$.

??? success "Solution"
    The ground state from (4.4.9) is
    $$\psi_0(x) = \left(\frac{m\omega}{\pi\hbar}\right)^{1/4}\exp\!\left(-\frac{m\omega x^2}{2\hbar}\right) = \frac{1}{(\pi\ell^2)^{1/4}}\exp\!\left(-\frac{x^2}{2\ell^2}\right).$$
    Since $|\psi_0|^2$ is even, $\langle x\rangle = 0$ by symmetry. For $\langle x^2\rangle$ use the standard Gaussian integral
    $$\int_{-\infty}^\infty x^2 e^{-x^2/\ell^2}\, dx = \frac{\ell^3 \sqrt\pi}{2}.$$
    Then
    $$\langle x^2\rangle = \frac{1}{\ell \sqrt\pi}\int_{-\infty}^\infty x^2 e^{-x^2/\ell^2}\, dx = \frac{1}{\ell\sqrt\pi}\cdot \frac{\ell^3 \sqrt\pi}{2} = \frac{\ell^2}{2}.$$
    Hence $\Delta x = \ell/\sqrt{2}$. By an analogous calculation $\Delta p = \hbar/(\sqrt{2}\ell)$, and the product saturates the uncertainty bound: $\Delta x\, \Delta p = \hbar/2$. The harmonic-oscillator ground state is a *minimum-uncertainty* state.

---

## Exercise 4.6 — Modifying the FD code for a double well ($\star\star$)

Take the finite-difference code of §4.3 and modify it to solve the double-well potential

$$V(x) = V_0 \left[\left(\frac{x}{a}\right)^2 - 1\right]^2,$$

with $V_0 = 1.0$ eV and $a = 0.5$ nm. Plot the lowest four eigenstates and their probability densities, and explain the pattern of (near-)degeneracies you see in terms of tunnelling between the two wells.

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

## Exercise 4.7 — Energies via Hellmann–Feynman ($\star\star$)

For a normalised energy eigenstate $|\psi_\lambda\rangle$ of a Hamiltonian $\hat{H}(\lambda)$ depending on a parameter $\lambda$, prove the Hellmann–Feynman theorem:

$$\frac{dE(\lambda)}{d\lambda} = \langle\psi_\lambda | \frac{\partial \hat{H}}{\partial \lambda} | \psi_\lambda\rangle.$$

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

## Exercise 4.8 — Estimating the size of a many-electron Hilbert space ($\star\star\star$)

Consider $N$ spinless electrons confined to a 1D box of width $L$, discretised on a grid of $N_g$ points. How many *antisymmetric* basis states are there? Evaluate this number for $N = 20$ electrons and $N_g = 100$ grid points, and compare with $N_g^N$ (the dimension if antisymmetry were not enforced).

??? success "Solution"
    Place a single electron on the grid in one of $N_g$ ways. With $N$ electrons and antisymmetry, two electrons cannot share the same grid point (Pauli, since we are taking them to be spinless or all of one spin). The number of antisymmetric basis states is therefore the number of ways to *choose* $N$ distinct grid points out of $N_g$ — independent of order, because the determinantal structure fixes the sign:
    $$\dim \mathcal H_{\mathrm{antisym}} = \binom{N_g}{N}.$$
    For $N_g = 100$, $N = 20$:
    $$\binom{100}{20} = \frac{100!}{20!\,80!} \approx 5.36\times 10^{20}.$$
    Without antisymmetry, the unrestricted Hilbert space would have $N_g^N = 100^{20} = 10^{40}$ basis states.
    Antisymmetry reduces the count by a factor of $N! \cdot (\text{double-occupancy combinatorics}) \sim N!$, but the resulting number is *still* astronomical — 20 electrons on a 100-point grid require half a sextillion basis states. This is a quantitative version of the exponential wall of §4.5, and explains why FCI is impossible for any system with more than ~20 electrons in modest basis sets.
    **Take-away.** Pauli antisymmetry helps — it removes a factor of $\sim N!$ from the basis count — but it does not save us from exponential scaling. Approximate methods (HF, DFT, MP2, CCSD, …) that scale polynomially in $N$ are the only way forward.

---

## Going further

If you want more practice, try:

- Solve numerically for an electron in a finite square well of depth $V_0 = 5$ eV and width $L = 1$ nm. Count the number of bound states and compare with the analytical transcendental-equation result.
- Repeat the SHO calculation in *atomic units* ($\hbar = m_e = 1$). Confirm that $E_n = \omega(n + 1/2)$ with $\omega = 1$ gives eigenvalues $0.5, 1.5, 2.5, \ldots$ hartree.
- For the double-well potential of Exercise 4.6, plot the tunnel splitting $\Delta E_{01} = E_1 - E_0$ as a function of barrier height $V_0$. Verify that it decays exponentially with $\sqrt{V_0}$ — the WKB prediction.
- Compute the Hartree–Fock ground-state energy of the helium atom analytically using a single-zeta Slater orbital $\phi(\mathbf r) = (\zeta^3/\pi)^{1/2} e^{-\zeta r}$, minimising over $\zeta$. Compare with the experimental value (-79.0 eV) and identify the correlation energy.

In Chapter 5 we will leave wavefunctions behind and reformulate everything in terms of the electron density. The mathematical machinery developed in this chapter — Hermitian operators, eigenvalue problems, finite-difference discretisations, SCF cycles — will all reappear, but in a form that can finally be applied to real materials.
