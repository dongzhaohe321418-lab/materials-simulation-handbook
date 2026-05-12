# 4.3 Solving it numerically — particle in a box

<figure markdown>
![Particle in a box eigenstates](../assets/figures/ch04/fig_particle_in_box.png){ width="700" }
<figcaption>Figure 4.3.1. The first four eigenstates of the 1D infinite square well. Energies scale as \(n^2\) (left panel); wavefunctions are sinusoids that vanish at the walls (right panel). The number of nodes increases by one with each level.</figcaption>
</figure>

We have a postulate (the Schrödinger equation) and a framework (Hermitian operators, eigenvalue problems, bra-ket notation). It is time to *solve* something. We pick the simplest non-trivial problem in quantum mechanics: a single particle confined to a one-dimensional region, with infinite potential walls. This is the "particle in a box", and it is the quantum analogue of a guitar string fixed at both ends.

The motivation is partly pedagogical and partly practical. Pedagogically, the box is the first place a reader meets quantised energies emerging *from a calculation* rather than as a postulate — the discreteness is forced on us by the boundary conditions, not put in by hand. Practically, the box is a surprisingly good caricature of an electron in a quantum well (a thin semiconductor layer) or in a long conjugated molecule like a polyene, and one can already make order-of-magnitude predictions about light absorption from this model.

We will solve the problem twice: once analytically with paper and pencil, and once numerically by turning the Hamiltonian into a matrix and diagonalising it. The numerical method we develop here — finite differences plus `scipy.linalg.eigh` — is exactly the method we will reuse in §4.4 for the harmonic oscillator and which, in spirit, underlies modern plane-wave electronic-structure codes.

## 4.3.1 The model

Consider a single particle of mass $m$ in one dimension, with potential

$$V(x) = \begin{cases} 0, & 0 < x < L,\\ \infty, & \text{otherwise.}\end{cases} \tag{4.3.1}$$

Inside the box the particle is free; outside, the infinite potential forbids any wavefunction amplitude. Continuity of $\psi$ therefore demands

$$\psi(0) = \psi(L) = 0. \tag{4.3.2}$$

Inside the box the time-independent Schrödinger equation (4.2.6) reads

$$-\frac{\hbar^2}{2m}\frac{d^2 \psi}{dx^2} = E\, \psi. \tag{4.3.3}$$

## 4.3.2 Analytical solution

Equation (4.3.3) is a linear second-order ODE with constant coefficients — the same equation that governs a simple harmonic oscillator in classical mechanics, with the spatial coordinate playing the role of time. Define

$$k^2 \equiv \frac{2mE}{\hbar^2}, \tag{4.3.4}$$

so that (4.3.3) becomes $\psi'' + k^2 \psi = 0$. The general real solution is

$$\psi(x) = A\sin(kx) + B\cos(kx). \tag{4.3.5}$$

Apply the boundary conditions. At $x = 0$,

$$\psi(0) = B = 0,$$

so $B = 0$ and the wavefunction is a pure sine. At $x = L$,

$$\psi(L) = A\sin(kL) = 0.$$

Either $A = 0$ (the trivial, unnormalisable solution we discard) or $\sin(kL) = 0$, i.e.\ $kL = n\pi$ for some integer $n$. Thus the allowed wavenumbers are

$$k_n = \frac{n\pi}{L}, \quad n = 1, 2, 3, \ldots \tag{4.3.6}$$

(Negative $n$ give the same wavefunction up to an overall sign and are discarded; $n = 0$ gives $\psi \equiv 0$.) The corresponding energies follow from (4.3.4):

$$\boxed{\; E_n = \frac{\hbar^2 k_n^2}{2m} = \frac{n^2 \pi^2 \hbar^2}{2 m L^2}, \quad n = 1, 2, 3, \ldots \;} \tag{4.3.7}$$

This is our first quantised spectrum. Three features deserve note.

- The energies scale as $n^2$: the levels are not equally spaced. Adjacent gaps grow with $n$.
- The lowest allowed energy is $E_1 = \pi^2 \hbar^2/(2mL^2) > 0$. Even in the ground state the particle cannot be at rest, as a classical particle could. This *zero-point energy* is a direct consequence of the uncertainty principle: confining the particle to a region of width $L$ forces $\Delta p \gtrsim \hbar/L$ and hence $E \gtrsim \hbar^2/(2mL^2)$. We will meet a closely-related zero-point energy in §4.4.
- The energies scale as $1/L^2$: a smaller box gives more widely spaced levels. This is why nanoscale confinement (quantum dots, quantum wells) produces tunable optical properties.

The wavefunctions are

$$\psi_n(x) = A_n \sin\!\left(\frac{n\pi x}{L}\right). \tag{4.3.8}$$

The amplitude $A_n$ is fixed by normalisation, equation (4.2.4) in one dimension:

$$1 = \int_0^L |\psi_n(x)|^2\, dx = A_n^2 \int_0^L \sin^2\!\left(\frac{n\pi x}{L}\right) dx.$$

Using $\sin^2\theta = \tfrac12(1 - \cos 2\theta)$,

$$\int_0^L \sin^2\!\left(\frac{n\pi x}{L}\right) dx = \frac{L}{2},$$

so $A_n^2 \cdot L/2 = 1$, hence $A_n = \sqrt{2/L}$. The normalised eigenfunctions are

$$\boxed{\; \psi_n(x) = \sqrt{\frac{2}{L}}\, \sin\!\left(\frac{n\pi x}{L}\right). \;} \tag{4.3.9}$$

One quick sanity check: the eigenfunctions are orthogonal. For $m \neq n$,

$$\int_0^L \psi_m^* \psi_n\, dx = \frac{2}{L}\int_0^L \sin\!\left(\frac{m\pi x}{L}\right)\sin\!\left(\frac{n\pi x}{L}\right) dx = 0,$$

using the standard sine-sine integral. This is the orthogonality theorem of §4.2.6 made explicit.

!!! example "Numerical scale"
    For an electron ($m_e = 9.109 \times 10^{-31}$ kg) in a box of $L = 1$ nm, the ground-state energy is
    $$E_1 = \frac{\pi^2 (1.055 \times 10^{-34})^2}{2 \cdot 9.109 \times 10^{-31} \cdot (10^{-9})^2} \approx 6.0 \times 10^{-20}\ \mathrm{J} \approx 0.376\ \mathrm{eV}.$$
    The first excited state is at $4 E_1 \approx 1.5$ eV, and the $1 \to 2$ transition occurs at a wavelength of $\sim 1100$ nm — the near-infrared. Make the box 0.5 nm and the transition shifts into the visible. This is the physics of quantum-confined optical materials.

## 4.3.3 Discretising the Hamiltonian

We now solve exactly the same problem numerically, with the explicit aim that the method should generalise to any 1D potential $V(x)$. The strategy is:

1. Replace the continuous coordinate $x \in [0, L]$ by a discrete grid of $N$ points.
2. Replace the second-derivative operator by a finite-difference approximation, turning $\hat H$ into a finite-size matrix.
3. Diagonalise the matrix to obtain approximate eigenvalues and eigenvectors of $\hat H$.

**The grid.** Place $N$ equally spaced points $x_1, x_2, \ldots, x_N$ inside the box, with spacing $h = L/(N+1)$ and positions $x_i = i\, h$ for $i = 1, \ldots, N$. The endpoints $x_0 = 0$ and $x_{N+1} = L$ are *not* part of the grid; the boundary conditions $\psi(0) = \psi(L) = 0$ are imposed by simply not including those points.

**The second derivative.** A Taylor expansion gives

$$\psi(x + h) = \psi(x) + h\psi'(x) + \frac{h^2}{2}\psi''(x) + \frac{h^3}{6}\psi'''(x) + \mathcal O(h^4),$$
$$\psi(x - h) = \psi(x) - h\psi'(x) + \frac{h^2}{2}\psi''(x) - \frac{h^3}{6}\psi'''(x) + \mathcal O(h^4).$$

Adding and rearranging,

$$\psi''(x) = \frac{\psi(x+h) - 2\psi(x) + \psi(x-h)}{h^2} + \mathcal O(h^2). \tag{4.3.10}$$

This is the **central second-difference** formula. On the grid, with $\psi_i \equiv \psi(x_i)$, it reads

$$\psi''(x_i) \approx \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{h^2}.$$

**The Hamiltonian matrix.** Inside the box $V = 0$, so $\hat H = -\frac{\hbar^2}{2m}\partial_x^2$, and the discrete Hamiltonian is the $N\times N$ matrix

$$H_{ij} = -\frac{\hbar^2}{2m h^2}\, \begin{cases} -2, & i = j,\\ 1, & |i - j| = 1,\\ 0, & \text{otherwise.}\end{cases} \tag{4.3.11}$$

In matrix form,

$$\hat H = \frac{\hbar^2}{2m h^2}\, \begin{pmatrix} 2 & -1 & & & \\ -1 & 2 & -1 & & \\ & -1 & 2 & -1 & \\ & & \ddots & \ddots & \ddots \\ & & & -1 & 2\end{pmatrix}. \tag{4.3.12}$$

This is a real symmetric tridiagonal matrix. Real symmetric matrices have real eigenvalues and orthogonal eigenvectors — the discrete analogue of our continuum theorem in §4.2.5–6. (The continuum operator is Hermitian; its finite-difference approximation is *symmetric*, which is the real version of the same condition.)

To include an arbitrary potential $V(x)$ we simply add a diagonal matrix:

$$H_{ii} \to H_{ii} + V(x_i). \tag{4.3.13}$$

The off-diagonal kinetic-energy part is the same for every problem. This is what makes the method so general.

**Boundary conditions.** Note that at the first grid point $i = 1$, the second-difference formula involves $\psi_0 \equiv \psi(0)$, which the Dirichlet boundary condition sets to zero — and so the term $-\psi_0/h^2$ simply does not contribute. Similarly at $i = N$. The matrix (4.3.11) implicitly enforces $\psi(0) = \psi(L) = 0$. Other boundary conditions (periodic, von Neumann, …) would modify the corners of the matrix.

## 4.3.4 A complete Python implementation

The script below solves the particle-in-a-box numerically and compares with the analytical answer. It uses SI units and is parameterised on $m$ and $L$, so you can change the mass or the box width with a single line.

```python
"""particle_in_a_box.py — Solve the 1D infinite square well by finite differences.

Reference: §4.3 of the Materials Simulation Handbook.
Requires: numpy, scipy, matplotlib.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
HBAR: float = 1.054_571_817e-34   # J s
M_E: float = 9.109_383_7e-31      # kg
EV: float = 1.602_176_634e-19     # J per eV


def build_hamiltonian(
    n_grid: int,
    box_length: float,
    mass: float = M_E,
    potential: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Build the finite-difference Hamiltonian on an interior grid.

    Parameters
    ----------
    n_grid : int
        Number of interior grid points (excluding the two walls).
    box_length : float
        Width L of the box, in metres.
    mass : float, optional
        Particle mass in kg (default: electron mass).
    potential : array_like, optional
        Optional V(x) sampled at the grid points, in joules. If omitted,
        the potential is zero inside the box.

    Returns
    -------
    x : np.ndarray, shape (n_grid,)
        Interior grid positions in metres.
    H : np.ndarray, shape (n_grid, n_grid)
        The Hamiltonian matrix, ready for diagonalisation.
    """
    h = box_length / (n_grid + 1)
    x = np.linspace(h, box_length - h, n_grid)

    # Kinetic part: tridiagonal (-1, 2, -1) scaled by hbar^2 / (2 m h^2).
    prefactor = HBAR**2 / (2.0 * mass * h**2)
    main = 2.0 * prefactor * np.ones(n_grid)
    off = -prefactor * np.ones(n_grid - 1)
    H = np.diag(main) + np.diag(off, k=1) + np.diag(off, k=-1)

    if potential is not None:
        if potential.shape != (n_grid,):
            raise ValueError("potential must have shape (n_grid,)")
        H = H + np.diag(potential)

    return x, H


def analytic_box(
    n: int,
    x: np.ndarray,
    box_length: float,
    mass: float = M_E,
) -> tuple[float, np.ndarray]:
    """Analytical eigenstate of the infinite square well.

    Returns the energy in joules and the (real, normalised) wavefunction
    sampled at x.
    """
    energy = (n**2 * np.pi**2 * HBAR**2) / (2.0 * mass * box_length**2)
    psi = np.sqrt(2.0 / box_length) * np.sin(n * np.pi * x / box_length)
    return energy, psi


def solve_and_plot(n_grid: int = 400, box_length: float = 1.0e-9) -> None:
    """Diagonalise H and compare the first four eigenstates with theory."""
    x, H = build_hamiltonian(n_grid, box_length)

    # Use a dense eigensolver here: 400 x 400 is trivial. For larger
    # problems use scipy.sparse.linalg.eigsh on a sparse matrix.
    eigvals, eigvecs = np.linalg.eigh(H)

    # Discrete eigenvectors must be rescaled so that sum |psi_i|^2 dx = 1.
    h = box_length / (n_grid + 1)
    eigvecs = eigvecs / np.sqrt(h)

    print(f"{'n':>3} {'E_num (eV)':>14} {'E_ana (eV)':>14} {'rel err':>10}")
    for n in range(1, 5):
        e_ana, _ = analytic_box(n, x, box_length)
        e_num = eigvals[n - 1]
        rel = abs(e_num - e_ana) / e_ana
        print(f"{n:>3d} {e_num/EV:>14.6f} {e_ana/EV:>14.6f} {rel:>10.2e}")

    fig, axes = plt.subplots(2, 2, figsize=(9, 6), sharex=True)
    for n, ax in zip(range(1, 5), axes.ravel()):
        e_ana, psi_ana = analytic_box(n, x, box_length)
        psi_num = eigvecs[:, n - 1]
        # Eigenvectors are defined up to a sign: align with the analytic.
        if np.dot(psi_num, psi_ana) < 0:
            psi_num = -psi_num
        ax.plot(x * 1e9, psi_ana, "k-", lw=2, label="analytic")
        ax.plot(x * 1e9, psi_num, "r--", lw=1.2, label="FD")
        ax.set_title(f"n = {n},  E = {eigvals[n-1]/EV:.3f} eV")
        ax.set_xlabel("x (nm)")
        ax.set_ylabel(r"$\psi_n(x)$")
        ax.legend()
    fig.tight_layout()
    plt.savefig("particle_in_a_box.png", dpi=140)


if __name__ == "__main__":
    solve_and_plot()
```

Run the script. Typical output (for $N = 400$ grid points and $L = 1$ nm) is:

```
  n      E_num (eV)      E_ana (eV)    rel err
  1        0.376024        0.376033   2.55e-05
  2        1.504099        1.504133   2.27e-05
  3        3.384067        3.384300   6.89e-05
  4        6.015977        6.016531   9.21e-05
```

Four-decimal agreement with theory on a 400-point grid — and the relative error scales as $h^2$, the order of the finite-difference truncation in (4.3.10). The plot shows the numerical eigenfunctions overlaid on the analytical sines: indistinguishable to the eye.

!!! note "What you have just done"
    You have solved a quantum mechanical eigenvalue problem with general-purpose linear algebra. The same code — with a different `potential` array — will solve *any* 1D Schrödinger equation. In §4.4 we will reuse it verbatim for the harmonic oscillator. The same idea, generalised to three dimensions and combined with a plane-wave basis instead of a position grid, is the engine inside Quantum ESPRESSO, VASP, ABINIT and most of the rest of the codes you will meet in Chapter 6.

## 4.3.5 Convergence and pitfalls

A few practical remarks.

**Grid spacing.** The error in the second-difference formula (4.3.10) scales as $h^2$, so halving $h$ should reduce the energy error by a factor of four. You can verify this empirically by running the script with `n_grid = 100, 200, 400, 800` and tabulating $E_1$.

**High-$n$ states.** Finite-difference methods are accurate for the *low-energy* eigenstates whose wavelengths span many grid points, but error grows rapidly when the wavelength approaches the grid spacing. With $N$ grid points you can trust roughly the first $N/10$ eigenstates. For the box, $\psi_n$ has $n$ half-wavelengths fitting into $L$, so $\lambda_n = 2L/n$. For the formula to be accurate we need $\lambda_n \gg h$, i.e.\ $n \ll 2L/h = 2(N+1)$.

**Sparse storage.** Our dense `np.diag` construction is wasteful: the Hamiltonian is tridiagonal and has only $3N$ non-zeros, not $N^2$. For large $N$ replace `np.diag` constructions with `scipy.sparse.diags` and use `scipy.sparse.linalg.eigsh` (Lanczos) to compute the lowest few eigenpairs. This is essential in 2D and 3D, where naive dense storage of an $N^3 \times N^3$ matrix would require terabytes.

**Units.** We have worked in SI throughout. In production electronic-structure codes the universal convention is *atomic units*: $\hbar = m_e = e = 4\pi\varepsilon_0 = 1$. Energies are then in *hartrees* ($1\ \mathrm{Ha} = 27.211$ eV) and lengths in *bohrs* ($1\ \mathrm{a_0} = 0.529$ Å). The Schrödinger equation becomes simply $(-\tfrac12 \nabla^2 + V)\psi = E\psi$, which is much tidier. We will switch to atomic units in Chapter 5.

**Boundary conditions matter.** Different physics calls for different boundary conditions. Solid-state problems use periodic boundaries (the Brillouin zone of Chapter 3); scattering problems use outgoing-wave conditions; molecular problems use $\psi \to 0$ at infinity. The Hamiltonian matrix changes correspondingly, but the basic strategy — discretise, build a sparse matrix, diagonalise — is the same.

## 4.3.6 Looking ahead

We have solved the simplest quantum mechanical problem twice over and met every ingredient that will reappear in more elaborate settings:

- a Hamiltonian operator,
- boundary conditions,
- a discrete spectrum,
- orthonormal eigenfunctions,
- a numerical scheme (finite differences) that turns the spectral problem into matrix diagonalisation.

In §4.4 we keep the numerical machinery exactly as it is and substitute a different potential — the harmonic oscillator. The analytical solution is more elaborate (Hermite polynomials), but the *code* is the same, with two lines changed. That is the point of working numerically: once the infrastructure is in place, every new physical problem reduces to specifying $V(x)$.
