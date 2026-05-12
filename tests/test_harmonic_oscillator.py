"""Test the 1D quantum harmonic-oscillator FD solver."""
from __future__ import annotations

import numpy as np

from tier1.ch04.harmonic_oscillator import HBAR, M_E, solve_harmonic


def test_ground_state_is_half_hbar_omega() -> None:
    omega = 1.0e15
    _, eigvals, _ = solve_harmonic(omega=omega, n_grid=600)
    e0_analytic = 0.5 * HBAR * omega
    rel_err = abs(eigvals[0] - e0_analytic) / e0_analytic
    assert rel_err < 1e-3, f"E0 rel_err = {rel_err:.2e}"


def test_levels_equally_spaced() -> None:
    """Adjacent quantum HO levels differ by hbar*omega."""
    omega = 1.0e15
    _, eigvals, _ = solve_harmonic(omega=omega, n_grid=800)
    spacings = np.diff(eigvals[:5])
    expected = HBAR * omega
    for s in spacings:
        rel = abs(s - expected) / expected
        assert rel < 1e-3, f"Spacing {s:.3e} vs {expected:.3e} rel_err {rel:.2e}"


def test_ground_state_node_count() -> None:
    """The ground state is nodeless wherever the amplitude is appreciable."""
    _, _, eigvecs = solve_harmonic(omega=1.0e15, n_grid=400)
    psi0 = eigvecs[:, 0]
    if psi0[len(psi0) // 2] < 0:
        psi0 = -psi0
    # Look only at the region where |psi0| is above 1% of its peak —
    # outside that, values are at numerical noise level and can flip sign.
    peak = float(np.max(np.abs(psi0)))
    significant = np.abs(psi0) > 0.01 * peak
    assert np.all(psi0[significant] > 0), (
        "Ground state should be nodeless where amplitude is appreciable"
    )
