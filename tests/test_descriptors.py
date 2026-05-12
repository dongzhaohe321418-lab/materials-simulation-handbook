"""Test the Behler G^2 descriptor implementation."""
from __future__ import annotations

import numpy as np

from tier1.ch09.descriptors_g2 import cosine_cutoff, g2_descriptor


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    """Random 3D rotation via QR decomposition."""
    A = rng.normal(size=(3, 3))
    Q, R = np.linalg.qr(A)
    # Ensure proper rotation (determinant +1).
    Q *= np.sign(np.diag(R))
    if np.linalg.det(Q) < 0:
        Q[:, 0] = -Q[:, 0]
    return Q


def test_cosine_cutoff_zero_at_rc() -> None:
    r = np.array([0.0, 3.0, 6.0, 7.0])
    fc = cosine_cutoff(r, r_c=6.0)
    assert fc[0] == 1.0
    assert fc[-1] == 0.0  # beyond cutoff
    # at r=6 the cosine becomes -1 so 0.5*((-1)+1) = 0, but `r < r_c` is False -> 0.
    assert fc[2] == 0.0


def test_g2_invariant_under_rotation() -> None:
    rng = np.random.default_rng(123)
    # Cubic supercell large enough that rotating positions stays inside box image
    cell = 20.0 * np.eye(3)
    positions = rng.uniform(7.0, 13.0, size=(6, 3))  # well inside
    eta = np.array([2.0, 4.0])
    r_s = np.array([2.0, 3.0])

    G1 = g2_descriptor(positions, cell, eta, r_s, r_c=5.0)
    R = _random_rotation(rng)
    centre = positions.mean(axis=0)
    rotated = (positions - centre) @ R.T + centre
    G2 = g2_descriptor(rotated, cell, eta, r_s, r_c=5.0)
    # Per-atom descriptors should match up to a permutation; since we did not
    # permute, they should match directly.
    assert np.allclose(np.sort(G1, axis=0), np.sort(G2, axis=0), atol=1e-8)


def test_g2_invariant_under_translation() -> None:
    rng = np.random.default_rng(0)
    cell = 20.0 * np.eye(3)
    positions = rng.uniform(5.0, 15.0, size=(4, 3))
    eta = np.array([2.0])
    r_s = np.array([2.5])
    G1 = g2_descriptor(positions, cell, eta, r_s, r_c=5.0)
    shift = np.array([1.3, -0.7, 0.5])
    G2 = g2_descriptor(positions + shift, cell, eta, r_s, r_c=5.0)
    assert np.allclose(G1, G2, atol=1e-10)


def test_g2_shape() -> None:
    cell = 10.0 * np.eye(3)
    positions = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]])
    eta = np.array([2.0, 4.0, 8.0])
    r_s = np.array([1.5, 2.0, 2.5])
    G = g2_descriptor(positions, cell, eta, r_s, r_c=5.0)
    assert G.shape == (2, 3)


def test_g2_nonnegative() -> None:
    rng = np.random.default_rng(2)
    cell = 10.0 * np.eye(3)
    positions = rng.uniform(0, 10, size=(10, 3))
    eta = np.array([2.0, 4.0])
    r_s = np.array([1.5, 3.0])
    G = g2_descriptor(positions, cell, eta, r_s, r_c=5.0)
    assert np.all(G >= 0.0)
