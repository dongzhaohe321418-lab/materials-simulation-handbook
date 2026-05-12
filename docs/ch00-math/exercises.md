# Chapter 0 — Exercises

Eight problems covering the material of Sections 0.1–0.5. Difficulty is marked with stars: <span class="diff-easy">★</span> easy, <span class="diff-med">★★</span> medium, <span class="diff-hard">★★★</span> hard. Worked solutions are collapsed below each problem; attempt the question first.

---

## Exercise 1 <span class="diff-easy">★</span> — Arrhenius slopes

A diffusion experiment on aluminium-doped silicon yields the following rate constants at two temperatures:

| $T$ / K | $k$ / s$^{-1}$ |
|---|---|
| 800 | $1.2 \times 10^{-4}$ |
| 1000 | $4.5 \times 10^{-3}$ |

Estimate the activation energy $E_\mathrm{a}$ in eV, assuming an Arrhenius rate law. Use $k_\mathrm{B} = 8.617 \times 10^{-5}\,\mathrm{eV/K}$.

??? success "Solution"
    From (0.1.4), $k(T) = A \exp(-E_\mathrm{a}/(k_\mathrm{B} T))$. Taking the ratio of the two measurements eliminates the prefactor:
    $$
    \frac{k_2}{k_1} = \exp\!\left[ \frac{E_\mathrm{a}}{k_\mathrm{B}} \left( \frac{1}{T_1} - \frac{1}{T_2} \right) \right].
    $$
    Taking logarithms,
    $$
    E_\mathrm{a} = k_\mathrm{B}\, \frac{\ln(k_2 / k_1)}{1/T_1 - 1/T_2}.
    $$
    With $k_2 / k_1 = 37.5$, $\ln 37.5 \approx 3.624$, and $1/800 - 1/1000 = 2.5 \times 10^{-4}\,\mathrm{K^{-1}}$,
    $$
    E_\mathrm{a} \approx (8.617 \times 10^{-5})\, \frac{3.624}{2.5 \times 10^{-4}} \approx 1.25\,\mathrm{eV}.
    $$
    A value in the right ballpark for substitutional diffusion in silicon.

---

## Exercise 2 <span class="diff-easy">★</span> — Domain and range

State the natural domain and range of each function, treating it as a function $\mathbb{R} \to \mathbb{R}$.

(a) $f(x) = \ln(1 - x^2)$.
(b) $g(x) = \dfrac{1}{x^2 + 1}$.
(c) $h(x) = \sqrt{x - 3}$.

??? success "Solution"
    (a) The argument of $\ln$ must be positive: $1 - x^2 > 0$, i.e. $x \in (-1, 1)$. As $x \to 0$, $\ln 1 = 0$; as $x \to \pm 1$, $\ln(1 - x^2) \to -\infty$. The range is $(-\infty, 0]$.

    (b) Domain is all of $\mathbb{R}$ since $x^2 + 1 \ge 1 > 0$. The denominator ranges over $[1, \infty)$, so $g$ ranges over $(0, 1]$, with $g(0) = 1$ the maximum.

    (c) The radicand must be non-negative: $x \ge 3$. Domain is $[3, \infty)$, range $[0, \infty)$.

---

## Exercise 3 <span class="diff-easy">★</span> — Dot and cross products

Let $\mathbf{a} = (1, 0, 0)$ and $\mathbf{b} = (1, 1, 0)$.

(a) Compute $\mathbf{a} \cdot \mathbf{b}$ and the angle between them.
(b) Compute $\mathbf{a} \times \mathbf{b}$ and confirm its magnitude equals the area of the parallelogram they span.
(c) Verify in NumPy.

??? success "Solution"
    (a) $\mathbf{a} \cdot \mathbf{b} = 1 \cdot 1 + 0 + 0 = 1$. Norms: $\lVert \mathbf{a} \rVert = 1$, $\lVert \mathbf{b} \rVert = \sqrt 2$. From (0.2.6), $\cos\theta = 1/\sqrt 2$, so $\theta = 45^\circ$.

    (b) From (0.2.7),
    $$
    \mathbf{a} \times \mathbf{b} = (0 \cdot 0 - 0 \cdot 1,\; 0 \cdot 1 - 1 \cdot 0,\; 1 \cdot 1 - 0 \cdot 1) = (0, 0, 1).
    $$
    Magnitude $1$. The parallelogram has base $\lVert\mathbf{a}\rVert = 1$ and height $\lVert\mathbf{b}\rVert \sin\theta = \sqrt 2 \cdot (1/\sqrt 2) = 1$, area $1$. Agreement.

    (c) Code:
    ```python
    import numpy as np
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([1.0, 1.0, 0.0])
    print(np.dot(a, b))             # 1.0
    print(np.degrees(np.arccos(np.dot(a, b) /
                               (np.linalg.norm(a) * np.linalg.norm(b)))))  # 45.0
    print(np.cross(a, b))           # [0. 0. 1.]
    print(np.linalg.norm(np.cross(a, b)))  # 1.0
    ```

---

## Exercise 4 <span class="diff-med">★★</span> — A 3x3 eigenvalue problem

Find the eigenvalues and eigenvectors of

$$
A = \begin{pmatrix} 2 & 0 & 0 \\ 0 & 3 & 1 \\ 0 & 1 & 3 \end{pmatrix}.
$$

Identify whether $A$ is positive definite, and verify numerically.

??? success "Solution"
    The block structure of $A$ — a $1 \times 1$ block of value $2$ and a $2 \times 2$ block in the lower-right — means the eigenvalues are the eigenvalues of the blocks taken separately.

    The first block gives eigenvalue $\lambda_1 = 2$ with eigenvector $\mathbf{v}_1 = (1, 0, 0)$.

    The second block is $\begin{pmatrix} 3 & 1 \\ 1 & 3 \end{pmatrix}$. By the recipe of Section 0.2, the characteristic polynomial is $(3 - \lambda)^2 - 1 = \lambda^2 - 6\lambda + 8 = 0$, giving $\lambda = 2$ or $\lambda = 4$. So $\lambda_2 = 2$ with $2$-block eigenvector $(1, -1)/\sqrt 2$, hence full eigenvector $\mathbf{v}_2 = (0, 1, -1)/\sqrt 2$; and $\lambda_3 = 4$ with $\mathbf{v}_3 = (0, 1, 1)/\sqrt 2$.

    All eigenvalues are strictly positive, so $A$ is positive definite. The eigenvalue $2$ has multiplicity $2$; any vector in the plane spanned by $\mathbf{v}_1$ and $\mathbf{v}_2$ is an eigenvector with eigenvalue $2$.

    ```python
    import numpy as np
    A = np.array([[2.0, 0.0, 0.0],
                  [0.0, 3.0, 1.0],
                  [0.0, 1.0, 3.0]])
    eigvals, eigvecs = np.linalg.eigh(A)
    print(eigvals)        # [2. 2. 4.]
    print(eigvecs)
    ```

---

## Exercise 5 <span class="diff-med">★★</span> — Taylor expansion and a harmonic potential

A diatomic molecule has potential energy

$$
U(r) = D\left[ 1 - e^{-\alpha (r - r_0)} \right]^2
$$

(the Morse potential), with $D, \alpha, r_0 > 0$.

(a) Show that $r = r_0$ is the unique minimum.
(b) Expand $U$ to second order about $r_0$ and identify the effective spring constant $k$.
(c) Give the harmonic angular frequency $\omega = \sqrt{k/\mu}$ in terms of the reduced mass $\mu$.

??? success "Solution"
    (a) Let $u = e^{-\alpha(r - r_0)}$. Then $U = D(1 - u)^2$ and
    $$
    \frac{\mathrm{d} U}{\mathrm{d} r} = 2 D (1 - u) \cdot (-1) \cdot \frac{\mathrm{d} u}{\mathrm{d} r} = 2 D (1 - u) \cdot \alpha u.
    $$
    The factor $\alpha u > 0$ always, so the derivative vanishes only when $u = 1$, i.e. $r = r_0$. Substituting back, $U(r_0) = 0$. To confirm it is a minimum, note $U \ge 0$ everywhere and equals $0$ only at $r_0$.

    (b) The Taylor expansion to second order needs $U(r_0)$, $U'(r_0)$ and $U''(r_0)$. We have $U(r_0) = 0$ and $U'(r_0) = 0$. For $U''$, differentiate the expression for $U'$:
    $$
    U' = 2 D \alpha\, u\,(1 - u), \qquad u = e^{-\alpha(r - r_0)},\; u' = -\alpha u.
    $$
    Product rule:
    $$
    U'' = 2 D \alpha \big[ u'(1 - u) + u \cdot (-u') \big] = 2 D \alpha\, u'\,(1 - 2u).
    $$
    At $r = r_0$, $u = 1$ and $u' = -\alpha$:
    $$
    U''(r_0) = 2 D \alpha\, (-\alpha)\,(1 - 2) = 2 D \alpha^2.
    $$
    Hence
    $$
    U(r) \approx \tfrac{1}{2} \cdot 2 D \alpha^2 \cdot (r - r_0)^2 = D \alpha^2 (r - r_0)^2,
    $$
    so the effective spring constant is $k = 2 D \alpha^2$.

    (c) $\omega = \sqrt{k/\mu} = \alpha \sqrt{2 D / \mu}$. This is the harmonic frequency that you would extract from the eigenvalues of a vibrational Hessian in Chapter 8.

---

## Exercise 6 <span class="diff-med">★★</span> — Central differences and convergence

Implement central differences for $f(x) = \sin x$ at $x = 1$, sweep the step size $h$ from $10^{-1}$ down to $10^{-15}$, and plot the absolute error against $h$ on a log-log scale.

(a) For what range of $h$ does the truncation error of $O(h^2)$ dominate?
(b) At what $h$ does round-off start to dominate, and why?
(c) What is roughly the minimum achievable error?

??? success "Solution"
    The reference value is $\cos 1 = 0.5403023\ldots$.

    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    x0: float = 1.0
    hs: np.ndarray = np.logspace(-15, -1, 200)
    errs: np.ndarray = np.abs(
        (np.sin(x0 + hs) - np.sin(x0 - hs)) / (2 * hs) - np.cos(x0)
    )

    plt.loglog(hs, errs)
    plt.loglog(hs, hs**2, '--', label=r'$h^2$ reference')
    plt.xlabel('h'); plt.ylabel('absolute error'); plt.legend(); plt.show()
    ```

    (a) For $h \gtrsim 10^{-5}$ or so, the error scales like $h^2$ — the slope on the log-log plot is $2$.

    (b) Below roughly $h \sim 10^{-5}$ the error starts to rise again because the numerator $\sin(x+h) - \sin(x-h)$ is the difference of two nearly equal floating-point numbers, suffering relative round-off of order $\varepsilon_\mathrm{mach}/h \sim 10^{-16}/h$.

    (c) The two error sources balance when $h^2 \sim \varepsilon_\mathrm{mach}/h$, i.e. $h \sim \varepsilon_\mathrm{mach}^{1/3} \sim 10^{-5}$. The minimum error is then $\sim \varepsilon_\mathrm{mach}^{2/3} \sim 10^{-11}$. This is the floor for naive central differences in double precision; achieving better accuracy needs Richardson extrapolation or complex-step differentiation.

---

## Exercise 7 <span class="diff-med">★★</span> — FFT of a Gaussian

Take $f(x) = e^{-x^2/(2 \sigma^2)}$ with $\sigma = 1$, sampled on $x \in [-10, 10]$ at $N = 2048$ equally spaced points. Compute the discrete Fourier transform with `np.fft.fft` and compare its magnitude with the analytical result $\tilde f(k) = \sigma\sqrt{2\pi}\, e^{-\sigma^2 k^2 / 2}$.

??? success "Solution"
    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    sigma: float = 1.0
    L: float = 20.0
    N: int = 2048
    x: np.ndarray = np.linspace(-L/2, L/2, N, endpoint=False)
    dx: float = x[1] - x[0]
    f: np.ndarray = np.exp(-x**2 / (2 * sigma**2))

    # NumPy's FFT returns sum_n f_n exp(-i 2 pi m n / N); to approximate
    # the continuous Fourier transform we multiply by dx and account for
    # the origin shift with fftshift.
    F: np.ndarray = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(f))) * dx
    k: np.ndarray = np.fft.fftshift(np.fft.fftfreq(N, d=dx)) * 2 * np.pi

    F_analytic: np.ndarray = sigma * np.sqrt(2 * np.pi) * np.exp(-sigma**2 * k**2 / 2)

    plt.plot(k, np.real(F), label='FFT (real part)')
    plt.plot(k, F_analytic, '--', label='analytic')
    plt.xlim(-6, 6); plt.xlabel('k'); plt.legend(); plt.show()

    print(np.max(np.abs(np.real(F) - F_analytic)))  # ~1e-10
    ```

    The Gaussian is its own Fourier-transform pair (up to scaling): a sharp Gaussian in $x$ gives a broad one in $k$, and vice versa — the position–momentum uncertainty relation made manifest. The two `*shift` calls are essential because NumPy's FFT places the zero-frequency component at index $0$ and treats the signal as $L$-periodic; `fftshift` rearranges to physicists' centred convention.

---

## Exercise 8 <span class="diff-hard">★★★</span> — From variational principle to eigenvalue problem

A trial wavefunction $\psi$ in a two-dimensional Hilbert space is expanded in an orthonormal basis $\{\phi_1, \phi_2\}$ as $\psi = c_1 \phi_1 + c_2 \phi_2$ with $c_1, c_2 \in \mathbb{C}$. Define the energy functional

$$
E[\psi] = \frac{\langle \psi | \hat H | \psi \rangle}{\langle \psi | \psi \rangle}.
$$

(a) Let $H_{ij} = \langle \phi_i | \hat H | \phi_j \rangle$ be the matrix elements of a Hermitian Hamiltonian ($H_{ji} = \overline{H_{ij}}$). Write $E$ as a function of $c_1, c_2$ in terms of the $H_{ij}$.

(b) Minimise $E$ over $c_1, c_2$ subject to the constraint $|c_1|^2 + |c_2|^2 = 1$. Show that the stationary points are eigenvectors of the matrix $H = (H_{ij})$ and that the corresponding stationary values are its eigenvalues.

(c) Specialise to $H = \begin{pmatrix} 0 & t \\ t & \Delta \end{pmatrix}$ with $t, \Delta > 0$ real. Compute the two eigenvalues and discuss the limits $t \ll \Delta$ and $t \gg \Delta$. This is the textbook two-level problem (bonding/antibonding orbitals, level repulsion).

??? success "Solution"
    (a) The numerator is
    $$
    \langle \psi | \hat H | \psi \rangle = \sum_{i,j} \overline{c_i} c_j H_{ij},
    $$
    which in matrix form is $\mathbf{c}^\dagger H \mathbf{c}$ with $\mathbf{c} = (c_1, c_2)^\top$. The denominator is $\mathbf{c}^\dagger \mathbf{c} = |c_1|^2 + |c_2|^2$. So
    $$
    E(\mathbf{c}) = \frac{\mathbf{c}^\dagger H \mathbf{c}}{\mathbf{c}^\dagger \mathbf{c}}.
    $$

    (b) Use the Lagrangian
    $$
    L(\mathbf{c}, \lambda) = \mathbf{c}^\dagger H \mathbf{c} - \lambda \big( \mathbf{c}^\dagger \mathbf{c} - 1 \big).
    $$
    Treat $\mathbf{c}$ and $\mathbf{c}^\dagger$ as independent (standard Wirtinger-calculus trick — see the discussion at the end of Section 0.3). The stationarity condition $\partial L / \partial \mathbf{c}^\dagger = 0$ gives
    $$
    H \mathbf{c} = \lambda \mathbf{c},
    $$
    exactly an eigenvalue equation. Multiplying on the left by $\mathbf{c}^\dagger$ and using the constraint $\mathbf{c}^\dagger \mathbf{c} = 1$,
    $$
    E = \mathbf{c}^\dagger H \mathbf{c} = \mathbf{c}^\dagger (\lambda \mathbf{c}) = \lambda,
    $$
    so the stationary value of $E$ is the eigenvalue. The lowest eigenvalue is the ground-state energy; this is the discrete two-state version of the variational theorem (0.3.19).

    (c) The characteristic polynomial of $H$ is
    $$
    \det(H - \lambda I) = -\lambda(\Delta - \lambda) - t^2 = \lambda^2 - \Delta \lambda - t^2 = 0,
    $$
    giving
    $$
    \lambda_{\pm} = \frac{\Delta}{2} \pm \sqrt{\left(\frac{\Delta}{2}\right)^2 + t^2}.
    $$
    The level splitting is $\lambda_+ - \lambda_- = 2 \sqrt{(\Delta/2)^2 + t^2}$, always at least $\Delta$ — the hallmark **level repulsion** of quantum mechanics. In the weak-coupling limit $t \ll \Delta$:
    $$
    \lambda_+ \approx \Delta + \frac{t^2}{\Delta}, \qquad \lambda_- \approx -\frac{t^2}{\Delta},
    $$
    a tiny perturbative shift. In the strong-coupling limit $t \gg \Delta$:
    $$
    \lambda_{\pm} \approx \frac{\Delta}{2} \pm t,
    $$
    the symmetric bonding/antibonding splitting familiar from molecular orbital theory. The latter case is exactly the LCAO picture of the $\mathrm{H}_2^+$ ion you will meet in Chapter 4.

    Numerical sanity check:
    ```python
    import numpy as np
    for t, Delta in [(0.1, 1.0), (1.0, 0.1)]:
        H = np.array([[0.0, t], [t, Delta]])
        print(t, Delta, np.linalg.eigvalsh(H))
    ```

---

That is the toolkit. Solving Exercise 8 by hand is genuine progress: you have just rederived the variational origin of the configuration-interaction matrix, the same equation Chapter 5 will solve with thousands of basis functions rather than two. Onward to the computational side of the book.
