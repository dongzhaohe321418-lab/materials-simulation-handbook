"""Test the graphene tight-binding band structure."""
from __future__ import annotations

import numpy as np

from tier1.ch03b.graphene_tb import (
    bands_along_path,
    hamiltonian,
    high_symmetry_points,
)


def test_bands_touch_zero_at_K() -> None:
    """Both bands vanish at the K point (Dirac point)."""
    K = high_symmetry_points()["K"]
    eigs = np.linalg.eigvalsh(hamiltonian(K))
    assert np.max(np.abs(eigs)) < 1e-10, (
        f"|E| at K = {np.max(np.abs(eigs)):.3e}, expected ~0"
    )


def test_hamiltonian_hermitian() -> None:
    rng = np.random.default_rng(0)
    for _ in range(5):
        k = rng.normal(size=2)
        H = hamiltonian(k)
        assert np.allclose(H, H.conj().T)


def test_bands_symmetric_about_zero() -> None:
    """The two graphene bands E_+ and E_- are negatives of each other."""
    pts = high_symmetry_points()
    path = np.array([pts["K"], pts["Gamma"], pts["M"], pts["K"]])
    _, E = bands_along_path(path, n_per_seg=50)
    assert np.allclose(E[:, 0], -E[:, 1], atol=1e-10)


def test_bandwidth_at_gamma() -> None:
    """At Gamma, |f| = 3t so the two eigenvalues are +-3t."""
    Gamma = np.array([0.0, 0.0])
    eigs = np.linalg.eigvalsh(hamiltonian(Gamma))
    assert abs(abs(eigs[1]) - 3 * 2.7) < 1e-8
