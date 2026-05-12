"""Solve the 1D quantum SHO by finite differences.

Extracted from docs/ch04-quantum/04-harmonic-oscillator.md (S 4.4.4).
"""
from __future__ import annotations

import numpy as np

HBAR: float = 1.054_571_817e-34
M_E: float = 9.109_383_7e-31
EV: float = 1.602_176_634e-19


def build_hamiltonian(
    x: np.ndarray,
    mass: float,
    potential: np.ndarray,
) -> np.ndarray:
    """1D finite-difference Hamiltonian on a regular grid x with V(x)."""
    h = x[1] - x[0]
    prefactor = HBAR**2 / (2.0 * mass * h**2)
    n = x.size
    main = 2.0 * prefactor * np.ones(n) + potential
    off = -prefactor * np.ones(n - 1)
    return np.diag(main) + np.diag(off, k=1) + np.diag(off, k=-1)


def solve_harmonic(
    omega: float = 1.0e15,
    mass: float = M_E,
    box_half_width: float = 4.0e-9,
    n_grid: int = 800,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve the SHO numerically.

    Returns (x grid, eigenvalues in joules, normalised eigenvectors).
    """
    x = np.linspace(-box_half_width, box_half_width, n_grid)
    h = x[1] - x[0]
    V = 0.5 * mass * omega**2 * x**2
    H = build_hamiltonian(x, mass, V)
    eigvals, eigvecs = np.linalg.eigh(H)
    eigvecs = eigvecs / np.sqrt(h)
    return x, eigvals, eigvecs


if __name__ == "__main__":
    omega = 1.0e15
    x, eigvals, _ = solve_harmonic(omega=omega)
    quantum = HBAR * omega
    print(f"hbar*omega = {quantum/EV:.6f} eV")
    for n in range(4):
        e_ana = quantum * (n + 0.5)
        e_num = eigvals[n]
        rel = abs(e_num - e_ana) / e_ana
        print(f"n={n} E_num={e_num/EV:.6f} eV  E_ana={e_ana/EV:.6f} eV  rel={rel:.2e}")
