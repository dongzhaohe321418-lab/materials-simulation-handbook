"""Behler G^2 radial symmetry function (descriptors).

Extracted from docs/ch09-mlip/03-descriptors.md (S 9.3.1).
"""
from __future__ import annotations
import numpy as np
from numpy.typing import NDArray


def cosine_cutoff(r: NDArray[np.float64], r_c: float) -> NDArray[np.float64]:
    # Clip r before cos to avoid invalid values for r = inf (the diagonal).
    r_safe = np.where(r < r_c, r, r_c)
    fc = 0.5 * (np.cos(np.pi * r_safe / r_c) + 1.0)
    return np.where(r < r_c, fc, 0.0)


def g2_descriptor(
    positions: NDArray[np.float64],
    cell: NDArray[np.float64],
    eta: NDArray[np.float64],
    r_s: NDArray[np.float64],
    r_c: float = 6.0,
) -> NDArray[np.float64]:
    """Compute Behler G^2 radial symmetry functions for every atom.

    positions: (N, 3) Cartesian coordinates (A).
    cell:      (3, 3) lattice vectors (rows).
    eta, r_s:  (K,) Gaussian widths and centres.
    """
    n_atoms = positions.shape[0]
    n_param = eta.shape[0]
    G = np.zeros((n_atoms, n_param))

    inv_cell = np.linalg.inv(cell)
    frac = positions @ inv_cell
    dfrac = frac[:, None, :] - frac[None, :, :]
    dfrac -= np.round(dfrac)
    dxyz = dfrac @ cell
    r = np.linalg.norm(dxyz, axis=-1)
    np.fill_diagonal(r, np.inf)

    mask = r < r_c
    fc = cosine_cutoff(r, r_c)

    for k in range(n_param):
        contrib = np.exp(-eta[k] * (r - r_s[k]) ** 2) * fc * mask
        G[:, k] = contrib.sum(axis=1)
    return G


if __name__ == "__main__":
    cell = 10.0 * np.eye(3)
    positions = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [0.0, 2.5, 0.0]])
    eta = np.array([4.0, 4.0])
    r_s = np.array([2.0, 3.0])
    G = g2_descriptor(positions, cell, eta, r_s, r_c=6.0)
    print(G)
