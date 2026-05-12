"""Minimal Kohn-Sham SCF solver for a 1D hydrogen chain (LDA exchange).

Extracted from docs/ch05-dft/05-scf.md (S 5.5.5).
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigh, solve


@dataclass
class Grid:
    n: int
    L: float

    @property
    def dx(self) -> float:
        return self.L / self.n

    @property
    def x(self) -> NDArray[np.float64]:
        return np.arange(self.n) * self.dx


def kinetic_matrix(g: Grid) -> NDArray[np.float64]:
    """Periodic-BC second-order finite-difference -1/2 d^2/dx^2."""
    n, dx = g.n, g.dx
    main = np.full(n, 1.0 / dx ** 2)
    off = np.full(n - 1, -0.5 / dx ** 2)
    T = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    T[0, -1] = T[-1, 0] = -0.5 / dx ** 2
    return T


def soft_coulomb(
    x: NDArray[np.float64], x0: float, L: float, a: float = 1.0
) -> NDArray[np.float64]:
    v = np.zeros_like(x)
    for m in range(-2, 3):
        dx_arr = x - x0 - m * L
        v += -1.0 / np.sqrt(dx_arr ** 2 + a ** 2)
    return v


def external_potential(g: Grid, positions: list[float]) -> NDArray[np.float64]:
    v = np.zeros(g.n)
    for x0 in positions:
        v += soft_coulomb(g.x, x0, g.L)
    return v


def hartree_potential(
    g: Grid, n: NDArray[np.float64], a: float = 1.0
) -> NDArray[np.float64]:
    vH = np.zeros(g.n)
    dx = g.dx
    for i, xi in enumerate(g.x):
        contrib = 0.0
        for m in range(-2, 3):
            d = g.x - xi - m * g.L
            contrib += np.sum(n / np.sqrt(d ** 2 + a ** 2)) * dx
        vH[i] = contrib
    return vH


def lda_exchange_potential(n: NDArray[np.float64]) -> NDArray[np.float64]:
    return -((3.0 / np.pi) ** (1.0 / 3.0)) * np.cbrt(np.maximum(n, 1e-12))


def lda_exchange_energy(n: NDArray[np.float64], dx: float) -> float:
    cx = -0.75 * (3.0 / np.pi) ** (1.0 / 3.0)
    return cx * float(np.sum(n ** (4.0 / 3.0)) * dx)


def build_density(
    orbitals: NDArray[np.float64], n_occ: int, dx: float
) -> NDArray[np.float64]:
    psi = orbitals[:, :n_occ].copy()
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2, axis=0) * dx)
    return 2.0 * np.sum(np.abs(psi) ** 2, axis=1)


def total_energy(
    eigvals: NDArray[np.float64],
    n: NDArray[np.float64],
    vH: NDArray[np.float64],
    v_ext: NDArray[np.float64],
    dx: float,
    n_occ: int,
) -> float:
    band = 2.0 * float(np.sum(eigvals[:n_occ]))
    eH = 0.5 * float(np.sum(n * vH) * dx)
    vx = lda_exchange_potential(n)
    ex = lda_exchange_energy(n, dx)
    e_dc = -eH + ex - float(np.sum(n * vx) * dx)
    return band + e_dc


def pulay_mix(
    n_in_hist: list[NDArray[np.float64]],
    n_out_hist: list[NDArray[np.float64]],
    m: int = 6,
) -> NDArray[np.float64]:
    k = len(n_in_hist)
    use = min(k, m)
    inputs = np.array(n_in_hist[-use:])
    outputs = np.array(n_out_hist[-use:])
    res = outputs - inputs
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


def scf(
    positions: list[float],
    n_electrons: int,
    grid: Grid,
    alpha: float = 0.3,
    tol: float = 1e-6,
    max_iter: int = 100,
    use_pulay_after: int = 3,
    verbose: bool = False,
) -> dict:
    n_occ = n_electrons // 2
    T = kinetic_matrix(grid)
    v_ext = external_potential(grid, positions)
    n = np.full(grid.n, n_electrons / grid.L)
    n_in_hist: list[NDArray[np.float64]] = []
    n_out_hist: list[NDArray[np.float64]] = []
    energies: list[float] = []
    eigvals = np.zeros(grid.n)
    eigvecs = np.zeros((grid.n, grid.n))
    v_ks = np.zeros(grid.n)
    converged = False
    for it in range(max_iter):
        vH = hartree_potential(grid, n)
        vx = lda_exchange_potential(n)
        v_ks = v_ext + vH + vx
        H = T + np.diag(v_ks)
        eigvals, eigvecs = eigh(H)
        n_out = build_density(eigvecs, n_occ, grid.dx)
        n_out *= n_electrons / (np.sum(n_out) * grid.dx)
        E = total_energy(eigvals, n_out, vH, v_ext, grid.dx, n_occ)
        energies.append(E)
        residual = float(np.max(np.abs(n_out - n)))
        if verbose:
            print(f"iter {it:3d}  E = {E: .6f} Ha   |dn|_inf = {residual:.2e}")
        if residual < tol:
            converged = True
            break
        n_in_hist.append(n.copy())
        n_out_hist.append(n_out.copy())
        if it < use_pulay_after:
            n = (1 - alpha) * n + alpha * n_out
        else:
            n = pulay_mix(n_in_hist, n_out_hist, m=6)
    return {
        "density": n,
        "x": grid.x,
        "eigvals": eigvals,
        "energies": energies,
        "v_ks": v_ks,
        "converged": converged,
        "n_iter": len(energies),
    }


if __name__ == "__main__":
    grid = Grid(n=128, L=20.0)
    positions = [4.0, 8.0, 12.0, 16.0]
    result = scf(positions, 4, grid, tol=1e-6, max_iter=80, verbose=True)
    print("converged:", result["converged"])
