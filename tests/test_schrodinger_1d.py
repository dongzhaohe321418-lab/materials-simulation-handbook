"""Test the 1D finite-difference particle-in-a-box solver."""
from __future__ import annotations

import numpy as np
import pytest

from tier1.ch04.schrodinger_1d import (
    HBAR,
    M_E,
    analytic_box,
    build_hamiltonian,
    solve_box,
)


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_eigenvalues_match_analytical(n: int) -> None:
    """Lowest four FD eigenvalues match n^2 pi^2 hbar^2 / (2 m L^2) within 1%."""
    box_length = 1.0e-9
    x, eigvals, _ = solve_box(n_grid=400, box_length=box_length)
    e_ana, _ = analytic_box(n, x, box_length)
    rel_err = abs(eigvals[n - 1] - e_ana) / e_ana
    assert rel_err < 1e-2, (
        f"State n={n}: numerical {eigvals[n-1]:.3e} vs analytical "
        f"{e_ana:.3e}, rel_err={rel_err:.2e}"
    )


def test_hamiltonian_symmetric() -> None:
    _, H = build_hamiltonian(50, 1.0e-9)
    assert np.allclose(H, H.T)


def test_potential_shape_validation() -> None:
    with pytest.raises(ValueError):
        build_hamiltonian(50, 1.0e-9, potential=np.zeros(49))


def test_eigenvalues_ordered_ascending() -> None:
    _, eigvals, _ = solve_box(n_grid=200, box_length=1.0e-9)
    assert np.all(np.diff(eigvals[:10]) > 0)


def test_eigenvalues_scale_as_n_squared() -> None:
    """E_n / E_1 should be ~n^2 for low states."""
    _, eigvals, _ = solve_box(n_grid=400, box_length=1.0e-9)
    for n in (2, 3, 4):
        ratio = eigvals[n - 1] / eigvals[0]
        assert abs(ratio - n ** 2) / (n ** 2) < 1e-2
