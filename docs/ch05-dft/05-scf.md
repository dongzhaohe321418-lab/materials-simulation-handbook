# 5.5 The Self-Consistent Field Loop

```mermaid
stateDiagram-v2
    [*] --> Guess : initial guess n⁰(r)
    Guess --> Veff : build v_KS[n] = v_ext + v_H[n] + v_xc[n]
    Veff --> Solve : solve (−½∇² + v_KS) φᵢ = εᵢ φᵢ
    Solve --> NewDensity : n_out(r) = Σ |φᵢ|² (occupied)
    NewDensity --> Mix : mix n_in and n_out (Pulay/Broyden)
    Mix --> Check : converged?
    Check --> Veff : no — iterate
    Check --> [*] : yes — output E, forces, ρ
```
*State diagram of the Kohn–Sham SCF loop: guess the density, build the effective potential, solve the one-electron equations, mix, repeat until self-consistent.*

<figure markdown>
![SCF total energy convergence](../assets/figures/ch05/fig_scf_convergence.png){ width="750" }
<figcaption>Figure 5.5.1. Typical SCF convergence behaviour (synthetic example). The total energy approaches the converged value approximately exponentially (left), and the energy change per step \(|\Delta E|\) drops below the user-specified threshold (here \(10^{-6}\) Ry) after a few tens of iterations (right). Real calculations may oscillate before locking in, especially for metals and magnetic systems.</figcaption>
</figure>

The Kohn–Sham equations,

$$
\Big[-\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}[n](\mathbf r)\Big]\phi_i(\mathbf r) = \varepsilon_i\,\phi_i(\mathbf r),
\qquad
n = \sum_i^\mathrm{occ}|\phi_i|^{2},
\qquad
v_\mathrm{KS}[n] = v_\mathrm{ext} + v_H[n] + v_{xc}[n],
$$

are nonlinear: the operator on the left depends, through $v_\mathrm{KS}$, on the density that the eigenfunctions themselves produce. The standard way to solve a nonlinear equation $n = \mathcal F[n]$ is fixed-point iteration: take an initial guess $n^{(0)}$, evaluate $n^{(1)} = \mathcal F[n^{(0)}]$, and repeat until $\|n^{(k+1)} - n^{(k)}\|$ is small.

In DFT this is the **self-consistent field** (SCF) loop. Naive iteration almost always fails to converge for systems beyond hydrogen-like atoms; the loop oscillates between charge-rich and charge-poor solutions, sometimes diverging outright. Practical SCF codes use *mixing* schemes — careful linear combinations of densities (or potentials) from successive iterations — to suppress these oscillations.

This section walks through the algorithm, explains why naive iteration fails, develops the mixing schemes (linear, Pulay/DIIS, Anderson), and ends with a complete Python implementation that solves a 1D model "hydrogen chain" using LDA exchange and finite differences.

## 5.5.1 The basic SCF algorithm

The textbook KS-SCF loop is:

1. **Initial guess.** Construct an initial density $n^{(0)}(\mathbf r)$. Common choices: superposition of free-atom densities, the previous SCF solution at a nearby geometry, or — for very simple systems — the uniform density.
2. **Build the potential.** Compute $v_\mathrm{KS}[n^{(k)}] = v_\mathrm{ext} + v_H[n^{(k)}] + v_{xc}[n^{(k)}]$.
3. **Diagonalise.** Solve the eigenvalue problem $\hat{H}_\mathrm{KS}\phi_i = \varepsilon_i\phi_i$ for the lowest $N_\mathrm{occ}$ eigenpairs.
4. **Form the new density.** $n_\mathrm{out}^{(k)}(\mathbf r) = \sum_i^\mathrm{occ}|\phi_i(\mathbf r)|^{2}$.
5. **Check convergence.** Compute residuals — change in density, change in energy, maximum force. If below tolerance, stop.
6. **Mix.** $n^{(k+1)} = \mathcal M(n^{(k)}, n_\mathrm{out}^{(k)}; \text{history})$. Go to step 2.

The interesting step is 6.

## 5.5.2 Why naive iteration fails

The naive mixing scheme is $n^{(k+1)} = n_\mathrm{out}^{(k)}$ — accept the output density unchanged. This is fixed-point iteration on the map $\mathcal F: n \mapsto n_\mathrm{out}[n]$.

Convergence of fixed-point iteration requires that the Jacobian (the *dielectric response* of the system) have spectral radius below unity in some norm. For metallic or polarisable systems it typically does not. Physically: if the density at one iteration has slightly too much charge in region $A$, the new Hartree potential pushes electrons out of $A$. The screening response sends *more* charge out of $A$ than the original excess — overshoot — and the next iteration has too little charge in $A$. The system oscillates with growing amplitude. This is **charge sloshing**.

The cure is to dampen the iteration: take only a fraction of the new density and combine it with the previous one. This is *linear mixing*:

$$
n^{(k+1)} = (1-\alpha)\,n^{(k)} + \alpha\,n_\mathrm{out}^{(k)},
\qquad \alpha \in (0,1].
\tag{5.42}
$$

Small $\alpha$ (e.g., $\alpha = 0.1$) almost always converges but does so slowly — convergence rate scales as $1 - \alpha$ per iteration. Large $\alpha$ converges fast when it converges and oscillates when it does not. For typical insulators $\alpha = 0.3$ is reasonable; for metals one often needs $\alpha = 0.05$.

Linear mixing is robust but slow. Modern codes use **acceleration schemes** based on the history of recent densities.

## 5.5.3 Pulay / DIIS mixing

The **Direct Inversion in the Iterative Subspace** (DIIS) method, due to Péter Pulay (1980), is a powerful workhorse. The idea: at iteration $k$, we have a history $\{n^{(j)}, n_\mathrm{out}^{(j)}\}_{j=k-m+1}^{k}$ of $m$ recent inputs and outputs. Define the **residual** of each:

$$
r^{(j)} \equiv n_\mathrm{out}^{(j)} - n^{(j)}.
$$

At self-consistency, $r = 0$. Search for the linear combination

$$
\bar n = \sum_j c_j\, n^{(j)},
\qquad \sum_j c_j = 1,
$$

that minimises the residual norm $\|\sum_j c_j r^{(j)}\|^{2}$. The solution is obtained by setting up the matrix

$$
B_{jk} = \langle r^{(j)}|r^{(k)}\rangle
$$

(here $\langle\cdot|\cdot\rangle$ is the $L^{2}$ inner product on the real-space grid) and solving the constrained linear system

$$
\begin{pmatrix} B & \mathbf 1 \\ \mathbf 1^{T} & 0 \end{pmatrix}
\begin{pmatrix} \mathbf c \\ -\lambda \end{pmatrix}
= \begin{pmatrix} \mathbf 0 \\ 1 \end{pmatrix},
\tag{5.43}
$$

where $\lambda$ is the Lagrange multiplier for the constraint $\sum c_j = 1$. The next-iteration density is

$$
n^{(k+1)} = \sum_j c_j\,n_\mathrm{out}^{(j)}
\qquad(\text{or}\;\;\sum_j c_j\,n^{(j)} + \alpha\sum_j c_j r^{(j)},\text{ DIIS with damping}).
\tag{5.44}
$$

In practice DIIS converges much faster than linear mixing — often quadratically near the solution. It needs a small history (typically $m = 6$–$10$). Far from convergence DIIS can be unstable; codes typically start with several linear-mixing steps before switching to DIIS, or fall back to linear mixing if DIIS diverges.

### Anderson mixing

Anderson mixing (1965) is an older, closely related method. At iteration $k$ with output residual $r^{(k)} = n_\mathrm{out}^{(k)} - n^{(k)}$, set

$$
n^{(k+1)} = n^{(k)} + \alpha\,r^{(k)} - \beta\big(r^{(k)} - r^{(k-1)}\big),
$$

with $\alpha,\beta$ chosen to minimise $\|r^{(k+1)}\|$ over the affine subspace. The general $m$-step Anderson method is essentially equivalent to DIIS; many modern codes use this formulation.

### Broyden's second method

Broyden mixing is a quasi-Newton scheme that approximates the inverse Jacobian of the SCF map from the iteration history. It generalises both Anderson and DIIS and is the default in some plane-wave codes (e.g., VASP). The implementation is more involved but the convergence behaviour is similar to DIIS for most problems.

## 5.5.4 Convergence criteria

Several quantities can be monitored:

- **Density residual** $\|n_\mathrm{out} - n_\mathrm{in}\| = \int|n_\mathrm{out} - n_\mathrm{in}|^{2}\,\mathrm d\mathbf r$ or $\max|n_\mathrm{out} - n_\mathrm{in}|$. The most rigorous criterion.
- **Energy difference** $|E^{(k+1)} - E^{(k)}|$. Easy to compute; typical tolerance $10^{-5}$ to $10^{-8}$ Ha. Beware: small energy changes do not always mean converged densities.
- **Force / stress changes**: critical for geometry optimisations. Want forces converged to $\sim 10^{-3}$ Ha/Bohr or better before trusting them.

A common protocol: require *both* an energy tolerance and a density tolerance to be satisfied for two consecutive iterations.

## 5.5.5 A complete Python implementation

We now solve the Kohn–Sham equations for a one-dimensional "hydrogen chain" model: $N_\mathrm{atom}$ protons placed on a line, with the electron–nucleus interaction softened to avoid the 1D Coulomb singularity, periodic boundary conditions on a finite box, and LDA exchange. We use a real-space finite-difference discretisation and direct diagonalisation. The code is ~150 lines, fully type-hinted, runnable on a laptop in under a second per iteration.

```python
"""scf_1d_hchain.py
A minimal Kohn--Sham SCF solver for a 1D hydrogen chain using LDA exchange.

Real-space finite-difference Hamiltonian on a uniform grid with periodic
boundary conditions. Linear and Pulay/DIIS mixing.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigh, solve
import matplotlib.pyplot as plt


@dataclass
class Grid:
    """Uniform 1D real-space grid with periodic boundary conditions."""
    n: int                # number of grid points
    L: float              # box length (Bohr)

    @property
    def dx(self) -> float:
        return self.L / self.n

    @property
    def x(self) -> NDArray[np.float64]:
        return np.arange(self.n) * self.dx


def kinetic_matrix(g: Grid) -> NDArray[np.float64]:
    """Second-order central finite-difference kinetic operator -1/2 d^2/dx^2."""
    n, dx = g.n, g.dx
    main = np.full(n, 1.0 / dx ** 2)
    off = np.full(n - 1, -0.5 / dx ** 2)
    T = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    # periodic boundary conditions
    T[0, -1] = T[-1, 0] = -0.5 / dx ** 2
    return T


def soft_coulomb(x: NDArray[np.float64], x0: float, L: float,
                 a: float = 1.0) -> NDArray[np.float64]:
    """Soft-Coulomb e-N attraction, periodic in box of length L.
    v(x) = -1/sqrt((x-x0)^2 + a^2), summed over periodic images."""
    v = np.zeros_like(x)
    for m in range(-2, 3):  # nearest images suffice
        dx_arr = x - x0 - m * L
        v += -1.0 / np.sqrt(dx_arr ** 2 + a ** 2)
    return v


def external_potential(g: Grid, positions: list[float]) -> NDArray[np.float64]:
    """Total e-N potential for protons at given positions."""
    v = np.zeros(g.n)
    for x0 in positions:
        v += soft_coulomb(g.x, x0, g.L)
    return v


def hartree_potential(g: Grid, n: NDArray[np.float64],
                      a: float = 1.0) -> NDArray[np.float64]:
    """Hartree potential from density n(x), using soft Coulomb kernel."""
    vH = np.zeros(g.n)
    dx = g.dx
    for i, xi in enumerate(g.x):
        # sum over periodic images
        contrib = 0.0
        for m in range(-2, 3):
            d = g.x - xi - m * g.L
            contrib += np.sum(n / np.sqrt(d ** 2 + a ** 2)) * dx
        vH[i] = contrib
    return vH


def lda_exchange_potential(n: NDArray[np.float64]) -> NDArray[np.float64]:
    """LDA exchange potential, 3D form applied to 1D effective density.
    v_x = -(3/pi)^{1/3} n^{1/3}. Pedagogical, not 1D-exact."""
    return -((3.0 / np.pi) ** (1.0 / 3.0)) * np.cbrt(np.maximum(n, 1e-12))


def lda_exchange_energy(n: NDArray[np.float64], dx: float) -> float:
    """LDA exchange energy per (5.36), pedagogical 3D form."""
    cx = -0.75 * (3.0 / np.pi) ** (1.0 / 3.0)
    return cx * float(np.sum(n ** (4.0 / 3.0)) * dx)


def build_density(orbitals: NDArray[np.float64], n_occ: int,
                  dx: float) -> NDArray[np.float64]:
    """Density from doubly occupied lowest orbitals; orbitals are columns."""
    psi = orbitals[:, :n_occ]
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2, axis=0) * dx)  # normalise
    return 2.0 * np.sum(np.abs(psi) ** 2, axis=1)


def total_energy(eigvals: NDArray[np.float64], n: NDArray[np.float64],
                 vH: NDArray[np.float64], v_ext: NDArray[np.float64],
                 dx: float, n_occ: int) -> float:
    """Kohn--Sham total energy from band-sum (5.28) form, LDA-X only."""
    band = 2.0 * float(np.sum(eigvals[:n_occ]))
    # double-counting corrections
    eH = 0.5 * float(np.sum(n * vH) * dx)
    vx = lda_exchange_potential(n)
    ex = lda_exchange_energy(n, dx)
    e_dc = -eH + ex - float(np.sum(n * vx) * dx)
    return band + e_dc


def pulay_mix(n_in_hist: list[NDArray[np.float64]],
              n_out_hist: list[NDArray[np.float64]],
              m: int = 6) -> NDArray[np.float64]:
    """Pulay/DIIS mixing using last m iterates."""
    k = len(n_in_hist)
    use = min(k, m)
    inputs = np.array(n_in_hist[-use:])
    outputs = np.array(n_out_hist[-use:])
    res = outputs - inputs  # residuals
    B = res @ res.T
    A = np.zeros((use + 1, use + 1))
    A[:use, :use] = B
    A[use, :use] = 1.0
    A[:use, use] = 1.0
    rhs = np.zeros(use + 1)
    rhs[-1] = 1.0
    try:
        sol = solve(A, rhs)
    except np.linalg.LinAlgError:
        return 0.7 * n_in_hist[-1] + 0.3 * n_out_hist[-1]
    coeffs = sol[:use]
    return coeffs @ outputs


def scf(positions: list[float], n_electrons: int, grid: Grid,
        alpha: float = 0.3, tol: float = 1e-6, max_iter: int = 100,
        use_pulay_after: int = 3) -> dict:
    """Run the SCF loop. Returns dict with density, eigvals, energies, etc."""
    n_occ = n_electrons // 2
    T = kinetic_matrix(grid)
    v_ext = external_potential(grid, positions)
    # initial guess: uniform density
    n = np.full(grid.n, n_electrons / grid.L)
    n_in_hist: list[NDArray[np.float64]] = []
    n_out_hist: list[NDArray[np.float64]] = []
    energies: list[float] = []
    for it in range(max_iter):
        vH = hartree_potential(grid, n)
        vx = lda_exchange_potential(n)
        v_ks = v_ext + vH + vx
        H = T + np.diag(v_ks)
        eigvals, eigvecs = eigh(H)
        n_out = build_density(eigvecs, n_occ, grid.dx)
        # renormalise to N electrons (FD discretisation drifts slightly)
        n_out *= n_electrons / (np.sum(n_out) * grid.dx)
        E = total_energy(eigvals, n_out, vH, v_ext, grid.dx, n_occ)
        energies.append(E)
        residual = float(np.max(np.abs(n_out - n)))
        print(f"iter {it:3d}  E = {E: .6f} Ha   |dn|_inf = {residual:.2e}")
        if residual < tol:
            print(f"Converged in {it+1} iterations.")
            break
        n_in_hist.append(n.copy())
        n_out_hist.append(n_out.copy())
        if it < use_pulay_after:
            n = (1 - alpha) * n + alpha * n_out
        else:
            n = pulay_mix(n_in_hist, n_out_hist, m=6)
    return {"density": n, "x": grid.x, "eigvals": eigvals,
            "energies": energies, "v_ks": v_ks}


def main() -> None:
    grid = Grid(n=256, L=20.0)
    positions = [4.0, 8.0, 12.0, 16.0]   # H4 chain, ~4 Bohr spacing
    n_electrons = 4
    result = scf(positions, n_electrons, grid, alpha=0.3, tol=1e-6,
                 max_iter=80)
    fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    axes[0].plot(result["x"], result["density"], lw=2)
    for x0 in positions:
        axes[0].axvline(x0, color="grey", ls=":", alpha=0.5)
    axes[0].set_ylabel(r"$n(x)$ (Bohr$^{-1}$)")
    axes[0].set_title("Converged Kohn-Sham density, 1D H$_4$ chain (LDA-X)")
    axes[1].plot(result["x"], result["v_ks"], lw=2)
    for x0 in positions:
        axes[1].axvline(x0, color="grey", ls=":", alpha=0.5)
    axes[1].set_xlabel("x (Bohr)")
    axes[1].set_ylabel(r"$v_\mathrm{KS}(x)$ (Ha)")
    fig.tight_layout()
    fig.savefig("scf_1d_hchain.png", dpi=150)
    print("Eigenvalues (Ha):", result["eigvals"][:6])


if __name__ == "__main__":
    main()
```

### How to read this code

- **Grid and operators.** A uniform real-space grid of $n=256$ points on a 20-Bohr box, with periodic boundary conditions. The kinetic operator is the second-order central difference $T = -\tfrac{1}{2}D^{2}$ assembled as a dense matrix; for larger systems one would use a sparse representation.
- **Soft Coulomb.** A 1D Coulomb $-1/|x-x_0|$ is singular; replacing $1/|x|$ with $1/\sqrt{x^{2}+a^{2}}$ regularises it and makes the model physically reasonable. The Hartree kernel uses the same softening.
- **LDA exchange.** We use the 3D LDA exchange potential $v_x \propto -n^{1/3}$ applied to our 1D density. This is *not* the rigorous 1D exchange (which is different functionally), but it is pedagogically standard and produces qualitatively correct behaviour.
- **Mixing.** The first few iterations use linear mixing with $\alpha = 0.3$ to stabilise the history; subsequent iterations switch to Pulay/DIIS using the most recent six densities.
- **Total energy.** Equation (5.28) is computed via the band-sum form: $E = 2\sum_i^\mathrm{occ}\varepsilon_i - U_H + E_x - \int n v_x$, accounting for double-counting between the band-sum $2\sum\varepsilon_i$ and the explicit Hartree/exchange energies.

### What you should see

Running `python scf_1d_hchain.py` produces:

- Convergence in roughly 10–20 iterations to $|\Delta n|_\infty < 10^{-6}$.
- A density peaked at each proton site, with smooth bonding charge in between — a one-dimensional analogue of the bond-charge build-up in a real H chain.
- Eigenvalues that group into a "band" of four nearly-degenerate states for the four-atom chain — the discrete analogue of the bonding band of an infinite 1D chain.
- A figure `scf_1d_hchain.png` showing the converged density and Kohn–Sham potential.

!!! note "Pedagogical, not production"
    This code is a teaching tool, not a research code. The Hartree integral is $\mathcal O(N^{2})$ in grid size (no FFT), the LDA exchange is the 3D form applied to a 1D problem, and there is no provision for spin polarisation, correlation, or geometry optimisation. Chapter 6 introduces production codes (Quantum ESPRESSO, VASP, ABINIT) that handle all of this rigorously.

## 5.5.6 Practical tips for SCF convergence

When your SCF fails to converge, here are the levers to pull, in approximate order of how often they help.

1. **Smear the occupations.** For metallic systems, the Fermi level sits in a band of states and small changes in the potential cause discontinuous re-occupation between iterations. Replace the integer occupations $f_i \in \{0,2\}$ with smooth Fermi–Dirac (or Gaussian, or Methfessel–Paxton) occupations $f_i = f(\varepsilon_i, \mu, T)$ at an artificial electronic temperature $T$ (typically $kT = 0.05$–$0.2$ eV). This is *always* needed for metals and often helpful for small-gap systems.
2. **Reduce the mixing parameter $\alpha$.** Try $\alpha = 0.1$ or even $0.05$ if oscillating.
3. **Increase the history depth $m$ for DIIS/Anderson.** From 6 to 10 or more.
4. **Mix the potential, not the density** (or vice versa). Some codes mix $v_\mathrm{KS}$; others mix $n$; for hard problems one can flip the choice.
5. **Pre-conditioning.** Long-wavelength density oscillations (charge sloshing in metals) decouple from short-wavelength ones. The **Kerker preconditioner** multiplies the residual in reciprocal space by $|\mathbf q|^{2}/(|\mathbf q|^{2}+q_0^{2})$ to dampen the long-wavelength components that drive oscillation.
6. **Better initial guess.** Restart from the previous geometry's wavefunctions in a relaxation; superposition of atomic densities is good; uniform is the worst.
7. **Tighten the grid / basis.** Sometimes SCF instabilities are artefacts of an under-converged basis or **k**-point sampling.

## 5.5.7 Closing the loop

We now have a complete picture: choose a functional (§5.4), discretise the wavefunctions on a grid (or basis), iterate the SCF loop with sensible mixing, and read out the converged density, total energy, eigenvalues, and (with the Hellmann–Feynman theorem) forces. From these, all the standard observables follow: cohesive energies, lattice constants, elastic moduli, vibrational spectra, band structures.

The next section, §5.6, gives an honest account of where this machine fails and what to do about it. Chapter 6 turns the cogs of a production calculation: plane-wave basis sets, pseudopotentials, **k**-point sampling, and practical convergence testing.
