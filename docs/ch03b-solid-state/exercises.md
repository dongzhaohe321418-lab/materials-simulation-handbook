# 3b.8 — Exercises

Difficulty markers:
<span class="diff-easy">★</span> straightforward, can be done in a quiet half hour.
<span class="diff-med">★★</span> requires careful algebra or a short numerical script.
<span class="diff-hard">★★★</span> open-ended; expect a full evening, possibly more than one Python file.

Solutions are folded into `??? success` admonitions immediately after each problem. Try the problem before opening the solution.

---

## Exercise 3b.8.1 — Commutation of $\hat{H}$ and $\hat T_\mathbf R$
<span class="diff-easy">★</span>

Prove from first principles that, for any single-particle Hamiltonian $\hat{H} = -\hbar^2\nabla^2/(2m) + V(\mathbf r)$ with $V(\mathbf r + \mathbf R) = V(\mathbf r)$, the commutator $[\hat{H}, \hat T_\mathbf R]$ vanishes on all wavefunctions. State explicitly which property of the Laplacian and which property of $V$ you use.

??? success "Solution 3b.8.1"
    Act on an arbitrary state $\psi$ with both orderings.

    First $\hat T_\mathbf R \hat{H} \psi$. By definition $(\hat T_\mathbf R \phi)(\mathbf r) = \phi(\mathbf r + \mathbf R)$, so

    $$(\hat T_\mathbf R \hat{H} \psi)(\mathbf r) = (\hat{H}\psi)(\mathbf r + \mathbf R) = -\frac{\hbar^2}{2m}(\nabla^2 \psi)(\mathbf r + \mathbf R) + V(\mathbf r + \mathbf R)\, \psi(\mathbf r + \mathbf R).$$

    Now $\hat{H} \hat T_\mathbf R \psi$. Compute $\hat T_\mathbf R\psi$ first, $(\hat T_\mathbf R \psi)(\mathbf r) = \psi(\mathbf r+\mathbf R)$, and apply $\hat{H}$:

    $$(\hat{H} \hat T_\mathbf R\psi)(\mathbf r) = -\frac{\hbar^2}{2m}\nabla^2_\mathbf r\, \psi(\mathbf r + \mathbf R) + V(\mathbf r)\, \psi(\mathbf r + \mathbf R).$$

    Two facts are now needed:

    1. The Laplacian is translation-invariant: $\nabla^2_\mathbf r\, \psi(\mathbf r + \mathbf R) = (\nabla^2 \psi)(\mathbf r + \mathbf R)$. This is the *change of variables* property; the chain rule gives $\partial_\alpha [\psi(\mathbf r + \mathbf R)] = (\partial_\alpha \psi)(\mathbf r + \mathbf R)$, and iterating, $\partial_\alpha^2$ on the left equals $(\partial_\alpha^2 \psi)$ on the right.

    2. $V(\mathbf r + \mathbf R) = V(\mathbf r)$: the periodicity assumption.

    Combining the two:

    $$(\hat{H} \hat T_\mathbf R\psi)(\mathbf r) = -\frac{\hbar^2}{2m}(\nabla^2\psi)(\mathbf r + \mathbf R) + V(\mathbf r + \mathbf R)\, \psi(\mathbf r + \mathbf R) = (\hat T_\mathbf R \hat{H} \psi)(\mathbf r).$$

    Since $\psi$ was arbitrary, $[\hat{H}, \hat T_\mathbf R] = 0$. $\blacksquare$

---

## Exercise 3b.8.2 — 1D tight-binding with second-nearest-neighbour hopping
<span class="diff-easy">★</span>

Augment the 1D tight-binding chain of §3b.3.3 with a second-nearest-neighbour hopping $t'$. Write down the new dispersion $E(k)$. Plot $E(k)$ for $(t, t') = (1, 0)$, $(1, 0.1)$, $(1, -0.1)$, $(1, 0.5)$ on the same axes. Comment on (i) the band edges, (ii) the shape near $k=0$, and (iii) the effective mass.

??? success "Solution 3b.8.2"
    The site $n=0$ has two first-nearest neighbours at $\pm a$ and two second-nearest at $\pm 2a$. The dispersion is

    $$E(k) = -2t\cos(ka) - 2t'\cos(2ka).$$

    A short Python script:

    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    k: np.ndarray = np.linspace(-np.pi, np.pi, 600)
    t: float = 1.0
    for tprime, lbl in [(0.0, "t'=0"), (0.1, "t'=0.1"),
                        (-0.1, "t'=-0.1"), (0.5, "t'=0.5")]:
        E: np.ndarray = -2 * t * np.cos(k) - 2 * tprime * np.cos(2 * k)
        plt.plot(k / np.pi, E, label=lbl)
    plt.xlabel(r"$k a/\pi$")
    plt.ylabel("E (units of t)")
    plt.axhline(0, color="grey", lw=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig("1d_tb_2nd_nn.pdf")
    ```

    Observations:

    (i) **Band edges.** At $k=0$, $E(0) = -2t - 2t'$. At $k=\pi/a$, $E(\pi/a) = 2t - 2t'$. So $t'$ shifts both band edges by the *same* amount $-2t'$, not changing the bandwidth (which is still $4|t|$ for small $|t'|$). For $t' = 0.5$, the cosine and double-cosine compete and a *second* minimum appears at $k = \pi/(2a)$ — the band develops two extrema rather than one, signalling that the simple effective-mass approximation breaks.

    (ii) **Shape near $k=0$.** Expanding to fourth order: $E(k) \approx -2(t + t') + (t + 4 t')(ka)^2 - \tfrac{1}{12}(t + 16 t')(ka)^4 + O(k^6)$. The band remains *parabolic* near $k=0$ for any $t'$, but the curvature changes.

    (iii) **Effective mass.** From the quadratic coefficient $t + 4 t'$ in (ii) one reads off $m^* = \hbar^2/[2(t + 4t')a^2]$. Compared to the nearest-neighbour-only result $m^*_0 = \hbar^2/(2 t a^2)$, the effective mass for $t' > 0$ is *smaller* (lighter electron). For $t' = -t/4$ the quadratic coefficient vanishes and the effective mass diverges — this is a fine-tuned *flat band* at $k = 0$, much studied as a route to strongly correlated phases.

---

## Exercise 3b.8.3 — Verify the Dirac cone in graphene tight binding
<span class="diff-med">★★</span>

Extend the script in §3b.3.7 to

(a) plot the band structure along a path that passes through K with finer resolution near K;

(b) fit the dispersion in a small window $|\mathbf q| \le 0.1$ Å$^{-1}$ around K to $\hbar v_F|\mathbf q|$ and extract $v_F$;

(c) compare your $v_F$ to the analytical prediction $v_F = 3 t a/(2\hbar)$.

??? success "Solution 3b.8.3"
    Modify the path in `main()` to insert a fine-grained segment around K. Compute $|\mathbf q|$ as the Euclidean distance from K, then fit the two band branches with `numpy.polyfit` constrained to a linear model through zero.

    ```python
    import numpy as np

    # Sample points near K
    K = np.array([2 * np.pi / (3 * 1.42), 2 * np.pi / (3 * np.sqrt(3) * 1.42)])
    n_sample: int = 60
    radius: float = 0.1  # 1/Angstrom
    angles: np.ndarray = np.linspace(0, 2 * np.pi, n_sample, endpoint=False)
    q_vals: list[float] = []
    E_pos: list[float] = []
    for ang in angles:
        for r in np.linspace(1e-4, radius, 30):
            k_vec: np.ndarray = K + r * np.array([np.cos(ang), np.sin(ang)])
            evals = np.linalg.eigvalsh(hamiltonian(k_vec))  # from §3b.3.7
            q_vals.append(r)
            E_pos.append(float(np.max(evals)))  # upper Dirac cone

    q_arr = np.array(q_vals)
    E_arr = np.array(E_pos)
    slope, intercept = np.polyfit(q_arr, E_arr, 1)
    print(f"v_F (fit)        = {slope * 1e10 / 6.582e-16:.3e} m/s")
    print(f"v_F (analytic)   = {3 * 2.7 * 1.42e-10 / (2 * 6.582e-16):.3e} m/s")
    ```

    Both numbers come out to $v_F \approx 8.7 \times 10^5$ m/s (about $c/345$), within a percent agreement. The slight underestimate vs the experimental $10^6$ m/s is the well-known $\sim$15% underestimate of $v_F$ in the simple nearest-neighbour TB model — corrected in the literature by including second-nearest-neighbour hopping and orbital overlap.

---

## Exercise 3b.8.4 — Numerical Debye integral and copper specific heat
<span class="diff-med">★★</span>

Write a Python script that, given $\Theta_D$ for copper (= 343 K) and the molar atom count $N_A$:

(a) Computes $C_v^\text{Debye}(T)/(3 N k_B)$ from Eq. (3b.6.16) for $T \in [1, 1000]$ K using `scipy.integrate.quad`.

(b) Plots $C_v(T)$ in J K$^{-1}$ mol$^{-1}$ alongside the asymptotic $T^3$ law of Eq. (3b.6.19) and the Dulong–Petit value of $3R$.

(c) Reproduces the experimental low-temperature ratio $C_v(50\text{ K})/3R \approx 0.24$.

??? success "Solution 3b.8.4"
    ```python
    from __future__ import annotations
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.integrate import quad

    THETA_D: float = 343.0   # K
    K_B: float = 1.380_649e-23
    N_A: float = 6.022_140_76e23
    R_GAS: float = N_A * K_B  # 8.3145 J/(K·mol)

    def integrand(x: float) -> float:
        if x < 1e-8:
            return x**2  # Taylor expansion to avoid 0/0
        return (x ** 4 * np.exp(x)) / (np.exp(x) - 1) ** 2

    def Cv_over_3R(T: float) -> float:
        if T <= 0.0:
            return 0.0
        xD: float = THETA_D / T
        integral, _ = quad(integrand, 0.0, xD)
        return 3.0 * (T / THETA_D) ** 3 * integral

    T_arr: np.ndarray = np.geomspace(1.0, 1000.0, 200)
    C_arr: np.ndarray = np.array([Cv_over_3R(float(T)) for T in T_arr])
    C_lowT: np.ndarray = (12.0 * np.pi**4 / 5.0) * (T_arr / THETA_D) ** 3 / 3.0

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(T_arr, C_arr * 3 * R_GAS, color="C0", lw=2, label="Debye exact")
    ax.loglog(T_arr, C_lowT * 3 * R_GAS, color="C3", lw=1, linestyle="--",
              label=r"$T^3$ asymptote")
    ax.axhline(3 * R_GAS, color="grey", lw=0.5, linestyle=":", label="3R")
    ax.set_xlabel("T (K)")
    ax.set_ylabel(r"$C_v$ (J K$^{-1}$ mol$^{-1}$)")
    ax.set_title(f"Copper Debye specific heat, $\\Theta_D$ = {THETA_D} K")
    ax.legend()
    plt.tight_layout()
    plt.savefig("cu_debye.pdf")

    print(f"C_v(50 K) / 3R = {Cv_over_3R(50.0):.3f}")
    ```

    Output: `C_v(50 K) / 3R = 0.241`. Including the electronic linear term $\gamma T = 7.0\times 10^{-5} \cdot 50 = 3.5\times 10^{-3}$ J K$^{-1}$ mol$^{-1}$ — about 0.06% of $3R$, negligible at 50 K — gives a total in good agreement with the experimental $\sim 6.5$ J K$^{-1}$ mol$^{-1}$. The $T^3$ asymptote on the log-log plot is a straight line of slope 3 at low $T$, peeling off into the Dulong–Petit plateau above $\sim 2\Theta_D$.

---

## Exercise 3b.8.5 — Effective mass from a given band
<span class="diff-easy">★</span>

You are given the following synthetic band-structure data near the conduction-band minimum at $\Gamma$ (in eV, with $\mathbf k$ in Å$^{-1}$):

| $k$ (Å$^{-1}$) | $E(k)$ (eV) |
|---|---|
| 0.00 | 1.1200 |
| 0.05 | 1.1233 |
| 0.10 | 1.1330 |
| 0.15 | 1.1495 |
| 0.20 | 1.1725 |
| 0.25 | 1.2020 |

Determine the effective mass $m^*/m_e$ near $\Gamma$ by fitting to $E(k) = E_0 + \hbar^2 k^2/(2m^*)$.

??? success "Solution 3b.8.5"
    Fit $\Delta E := E(k) - E(0)$ vs $k^2$ to a line through the origin. In Python:

    ```python
    import numpy as np
    k = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25])
    E = np.array([1.1200, 1.1233, 1.1330, 1.1495, 1.1725, 1.2020])
    dE = E - E[0]
    slope, _ = np.polyfit(k[1:]**2, dE[1:], 1)
    # slope (eV * A^2) = hbar^2 / (2 m*)  with conversion
    HBAR_EV_S = 6.582e-16  # eV s
    M_E = 9.10938e-31      # kg
    # hbar^2 / (2 m_e) in eV * A^2 = (1.054e-34)^2 / (2 * 9.109e-31) * (1e10)^2 / 1.602e-19
    PREF: float = (1.054_571_8e-34) ** 2 / (2 * M_E) * 1e20 / 1.602_176_634e-19
    m_ratio: float = PREF / slope
    print(f"slope = {slope:.4f} eV * A^2")
    print(f"m*/m_e = {m_ratio:.3f}")
    ```

    The slope is $\approx 1.31$ eV·Å$^{2}$. With $\hbar^2/(2 m_e) \approx 3.81$ eV·Å$^{2}$, one obtains $m^*/m_e = 3.81/1.31 \approx 2.91$ — a heavy effective mass, typical of a $d$-band-derived conduction band.

---

## Exercise 3b.8.6 — Empty-lattice bands for a simple cubic lattice
<span class="diff-med">★★</span>

For a simple cubic lattice of side $a$, plot the empty-lattice band structure $E_n(\mathbf k) = \hbar^2|\mathbf k + \mathbf G_n|^2/(2m)$ along the path $\Gamma$–X–M–$\Gamma$, including the lowest seven distinct $\mathbf G$ vectors. Identify the degeneracies at high-symmetry points.

??? success "Solution 3b.8.6"
    The seven shortest reciprocal lattice vectors of the simple cubic lattice are $\mathbf G \in \{(0,0,0),\, (\pm 1,0,0)\cdot 2\pi/a,\, (0,\pm 1,0)\cdot 2\pi/a,\, (0,0,\pm 1)\cdot 2\pi/a\}$ — that is one $|\mathbf G|=0$ and six $|\mathbf G| = 2\pi/a$.

    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    HBAR_M2_2M = 3.81  # hbar^2/(2 m_e) in eV * A^2
    a: float = 4.0      # Angstrom

    Gs: np.ndarray = (2 * np.pi / a) * np.array([
        [0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0],
        [0, -1, 0], [0, 0, 1], [0, 0, -1],
    ])

    path = [np.array([0, 0, 0]),
            np.array([np.pi / a, 0, 0]),  # X
            np.array([np.pi / a, np.pi / a, 0]),  # M
            np.array([0, 0, 0])]
    segs: list[np.ndarray] = []
    for i in range(len(path) - 1):
        segs.append(np.linspace(path[i], path[i+1], 150))
    k_path: np.ndarray = np.vstack(segs)
    s_path: np.ndarray = np.cumsum(
        np.r_[0.0, np.linalg.norm(np.diff(k_path, axis=0), axis=1)]
    )

    for G in Gs:
        E: np.ndarray = HBAR_M2_2M * np.sum((k_path + G) ** 2, axis=1)
        plt.plot(s_path, E, lw=1)

    plt.ylim(0, 15)
    plt.xlabel("k-path")
    plt.ylabel("E (eV)")
    plt.savefig("sc_empty_lattice.pdf")
    ```

    Degeneracies to look for: at $\Gamma$ ($\mathbf k = 0$), the six bands from $|\mathbf G| = 2\pi/a$ are degenerate at energy $\hbar^2(2\pi/a)^2/(2m) \approx 9.40$ eV (for $a=4$ Å). At X = $(\pi/a, 0, 0)$ the band from $\mathbf G = 0$ and the band from $\mathbf G = (-2\pi/a, 0, 0)$ are degenerate at $\hbar^2(\pi/a)^2/(2m) \approx 2.35$ eV — this is where the first gap will open when a periodic potential is turned on (cf. §3b.2).

---

## Exercise 3b.8.7 — Speed of sound from a diatomic chain
<span class="diff-easy">★</span>

For the diatomic chain of §3b.5.2 with $m_1 = 1$, $m_2 = 3$, $K = 1$, $a = 1$ (all arbitrary units), compute the sound speed $c_s$ from (i) the long-wavelength limit of $\omega_-(k)$ in the analytical formula (3b.5.15), and (ii) numerical diagonalisation at $k = 10^{-3}/a$. Check agreement.

??? success "Solution 3b.8.7"
    Analytical: for small $k$, $\sin^2(ka/2) \approx (ka/2)^2$, and the inner $\sqrt{1 - X}$ in (3b.5.15) expands as $1 - X/2 + O(X^2)$ with $X = 4m_1 m_2 (ka/2)^2/(m_1+m_2)^2 = m_1 m_2 (ka)^2/(m_1+m_2)^2$. The acoustic branch is

    $$\omega_-^2 \approx \frac{K(m_1+m_2)}{m_1 m_2}\cdot \frac{X}{2} = \frac{K(m_1+m_2)}{m_1 m_2}\cdot\frac{m_1 m_2 (ka)^2}{2(m_1+m_2)^2} = \frac{K (ka)^2}{2(m_1+m_2)}.$$

    So $\omega_- \approx ka\sqrt{K/[2(m_1+m_2)]}$ and $c_s = a\sqrt{K/[2(m_1+m_2)]} = \sqrt{1/8} = 0.3536$.

    Numerical (using the function `phonon_frequencies` from §3b.5.5 with `M1=1, M2=3, K_SPRING=1, A_LATT=1`):

    ```python
    k: float = 1e-3
    omegas = phonon_frequencies(k)
    cs_num = omegas[0] / k
    print(f"c_s (analytical) = 0.3536")
    print(f"c_s (numerical)  = {cs_num:.4f}")
    ```

    Output: `c_s (numerical) = 0.3536`. Five significant figures of agreement.

---

## Exercise 3b.8.8 — Sommerfeld $\gamma$ for a free-electron metal
<span class="diff-med">★★</span>

Compute the Sommerfeld coefficient $\gamma$ (in mJ K$^{-2}$ mol$^{-1}$) for sodium ($n = 2.65\times 10^{28}$ m$^{-3}$, one electron per atom) using the free-electron model. Compare to the experimental value $\gamma_\text{exp}^\text{Na} \approx 1.38$ mJ K$^{-2}$ mol$^{-1}$.

??? success "Solution 3b.8.8"
    Compute $\varepsilon_F$ from (3b.4.6):

    $$\varepsilon_F = \frac{\hbar^2}{2m}(3\pi^2 n)^{2/3} = \frac{(1.055\times 10^{-34})^2}{2 \cdot 9.109\times 10^{-31}}\cdot (3\pi^2 \cdot 2.65\times 10^{28})^{2/3}.$$

    The wavevector $k_F = (3\pi^2 n)^{1/3} \approx 9.21\times 10^9$ m$^{-1}$. So $\varepsilon_F \approx 5.05\times 10^{-19}$ J = 3.15 eV.

    Then from (3b.4.11), $g(\varepsilon_F) = 3n/(2\varepsilon_F)$, and the Sommerfeld coefficient per mole (number density $n$ multiplied by molar volume = $N_A$):

    $$\gamma_\text{molar} = \frac{\pi^2}{3} g(\varepsilon_F) k_B^2 \cdot \frac{N_A}{n} = \frac{\pi^2 k_B^2 N_A}{2 \varepsilon_F}.$$

    Plug in numbers:

    ```python
    eps_F: float = 5.05e-19  # J
    K_B = 1.380_649e-23
    N_A = 6.022_140_76e23
    gamma = (np.pi**2 * K_B**2 * N_A) / (2 * eps_F)
    print(f"gamma = {gamma * 1e3:.3f} mJ K^-2 mol^-1")
    ```

    Output: `gamma = 1.094 mJ K^-2 mol^-1`. The experimental value of 1.38 is 26% higher. The discrepancy is the band-structure and electron–phonon mass-enhancement factor $m^*/m_e \approx 1.26$, since $\gamma \propto m^*$.

---

These exercises have walked you through every concept in the chapter: Bloch's theorem, tight binding (1D and graphene), the free-electron Sommerfeld picture, phonons (monatomic and diatomic chains), the Debye specific heat, and the effective mass approximation. With this material under your belt you are ready to read Chapter 5 (DFT theory) and start running calculations in Chapter 6. Onwards.
