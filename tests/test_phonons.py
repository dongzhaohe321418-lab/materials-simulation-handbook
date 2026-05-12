"""Test the 1D diatomic-chain phonon dispersion."""
from __future__ import annotations

import numpy as np

from tier1.ch03b.phonons_diatomic import (
    M1,
    M2,
    K_SPRING,
    analytical_branches,
    compute_bands,
    phonon_frequencies,
)


def test_acoustic_branch_starts_at_zero() -> None:
    """omega_- (k=0) should vanish for any mass ratio."""
    om_acoustic = phonon_frequencies(0.0)[0]
    assert om_acoustic < 1e-10, f"omega_- (k=0) = {om_acoustic:.3e}"


def test_optical_branch_finite_at_gamma() -> None:
    """omega_+ (k=0) = sqrt(2K/mu) where mu is the reduced mass."""
    mu = M1 * M2 / (M1 + M2)
    expected = np.sqrt(2 * K_SPRING / mu)
    om_optical = phonon_frequencies(0.0)[1]
    assert abs(om_optical - expected) / expected < 1e-8


def test_numerical_matches_analytical() -> None:
    k_vals, omega = compute_bands(200)
    om_m, om_p = analytical_branches(k_vals)
    assert np.max(np.abs(omega[:, 0] - om_m)) < 1e-8
    assert np.max(np.abs(omega[:, 1] - om_p)) < 1e-8


def test_optical_above_acoustic() -> None:
    _, omega = compute_bands(64)
    assert np.all(omega[:, 1] >= omega[:, 0] - 1e-12)


def test_zone_boundary_gap() -> None:
    """At k=pi/a: omega_+^2 = 2K/min(m1,m2), omega_-^2 = 2K/max(m1,m2)."""
    om = phonon_frequencies(np.pi)
    assert abs(om[0] ** 2 - 2 * K_SPRING / max(M1, M2)) < 1e-8
    assert abs(om[1] ** 2 - 2 * K_SPRING / min(M1, M2)) < 1e-8
