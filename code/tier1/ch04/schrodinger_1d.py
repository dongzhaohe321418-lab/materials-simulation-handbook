"""Solve the 1D infinite square well by finite differences.

Extracted from docs/ch04-quantum/03-particle-in-box.md (S 4.3.4).
"""
from __future__ import annotations

import numpy as np

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

    Returns the interior grid positions and the dense Hamiltonian matrix.
    """
    h = box_length / (n_grid + 1)
    x = np.linspace(h, box_length - h, n_grid)

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

    Returns (energy in joules, normalised wavefunction on x).
    """
    energy = (n**2 * np.pi**2 * HBAR**2) / (2.0 * mass * box_length**2)
    psi = np.sqrt(2.0 / box_length) * np.sin(n * np.pi * x / box_length)
    return energy, psi


def solve_box(
    n_grid: int = 400,
    box_length: float = 1.0e-9,
    mass: float = M_E,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Diagonalise H and return grid, eigenvalues (J), normalised eigenvectors."""
    x, H = build_hamiltonian(n_grid, box_length, mass=mass)
    eigvals, eigvecs = np.linalg.eigh(H)
    h = box_length / (n_grid + 1)
    eigvecs = eigvecs / np.sqrt(h)
    return x, eigvals, eigvecs


if __name__ == "__main__":
    x, eigvals, _ = solve_box()
    print(f"{'n':>3} {'E_num (eV)':>14} {'E_ana (eV)':>14} {'rel err':>10}")
    for n in range(1, 5):
        e_ana, _ = analytic_box(n, x, 1.0e-9)
        e_num = eigvals[n - 1]
        rel = abs(e_num - e_ana) / e_ana
        print(f"{n:>3d} {e_num/EV:>14.6f} {e_ana/EV:>14.6f} {rel:>10.2e}")
