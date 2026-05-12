"""Tight-binding band structure of graphene along K-Gamma-M-K.

Extracted from docs/ch03b-solid-state/03-tight-binding.md (S 3b.3.7).
"""
from __future__ import annotations
import numpy as np
import numpy.typing as npt

T_HOP: float = 2.7       # nearest-neighbour hopping amplitude (eV)
A_CC: float = 1.42       # carbon-carbon bond length (Angstrom)

DELTAS: npt.NDArray[np.float64] = A_CC * np.array([
    [1.0, 0.0],
    [-0.5,  0.5 * np.sqrt(3.0)],
    [-0.5, -0.5 * np.sqrt(3.0)],
])


def hamiltonian(k: npt.NDArray[np.float64]) -> npt.NDArray[np.complex128]:
    """Return the 2x2 tight-binding Hamiltonian at wavevector k (1/A)."""
    f: complex = -T_HOP * np.sum(np.exp(1j * DELTAS @ k))
    H = np.array(
        [[0.0 + 0j, f],
         [np.conj(f), 0.0 + 0j]],
        dtype=np.complex128,
    )
    return H


def bands_along_path(
    path: npt.NDArray[np.float64],
    n_per_seg: int = 200,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Diagonalise H(k) along a piecewise-linear k-path."""
    s_list: list[float] = []
    E_list: list[npt.NDArray[np.float64]] = []
    s_cum: float = 0.0
    for i in range(len(path) - 1):
        k_a, k_b = path[i], path[i + 1]
        seg_len: float = float(np.linalg.norm(k_b - k_a))
        for j in range(n_per_seg):
            frac = j / n_per_seg
            k = (1 - frac) * k_a + frac * k_b
            evals = np.linalg.eigvalsh(hamiltonian(k))
            E_list.append(evals)
            s_list.append(s_cum + frac * seg_len)
        s_cum += seg_len
    evals_end = np.linalg.eigvalsh(hamiltonian(path[-1]))
    E_list.append(evals_end)
    s_list.append(s_cum)
    return np.array(s_list), np.array(E_list)


def high_symmetry_points() -> dict[str, npt.NDArray[np.float64]]:
    return {
        "Gamma": np.array([0.0, 0.0]),
        "K": np.array(
            [2 * np.pi / (3 * A_CC), 2 * np.pi / (3 * np.sqrt(3.0) * A_CC)]
        ),
        "M": np.array([2 * np.pi / (3 * A_CC), 0.0]),
    }


if __name__ == "__main__":
    pts = high_symmetry_points()
    path = np.array([pts["K"], pts["Gamma"], pts["M"], pts["K"]])
    s, E = bands_along_path(path, n_per_seg=300)
    print(f"E range: [{E.min():.3f}, {E.max():.3f}] eV")
    print(f"|E| at K point: {np.abs(np.linalg.eigvalsh(hamiltonian(pts['K']))).max():.2e}")
