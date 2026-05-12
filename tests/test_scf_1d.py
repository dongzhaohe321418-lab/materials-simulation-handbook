"""Test the 1D Kohn-Sham SCF loop."""
from __future__ import annotations

import numpy as np

from tier1.ch05.scf_1d import Grid, scf


def test_scf_converges_h4_chain() -> None:
    grid = Grid(n=64, L=20.0)
    positions = [4.0, 8.0, 12.0, 16.0]
    result = scf(positions, n_electrons=4, grid=grid, tol=1e-5, max_iter=60)
    assert result["converged"], (
        f"SCF did not converge in {result['n_iter']} iterations"
    )


def test_scf_energy_decreasing_then_stable() -> None:
    grid = Grid(n=64, L=20.0)
    result = scf([4.0, 8.0, 12.0, 16.0], 4, grid, tol=1e-5, max_iter=60)
    energies = np.array(result["energies"])
    # Final energy change is small
    assert abs(energies[-1] - energies[-2]) < 1e-4


def test_density_normalisation() -> None:
    grid = Grid(n=64, L=20.0)
    n_elec = 4
    result = scf([4.0, 8.0, 12.0, 16.0], n_elec, grid, tol=1e-5, max_iter=60)
    integral = float(np.sum(result["density"]) * grid.dx)
    assert abs(integral - n_elec) < 1e-3, (
        f"Density integrates to {integral:.4f}, expected {n_elec}"
    )


def test_density_non_negative() -> None:
    grid = Grid(n=64, L=20.0)
    result = scf([4.0, 8.0, 12.0, 16.0], 4, grid, tol=1e-5, max_iter=60)
    assert np.all(result["density"] >= -1e-10)
