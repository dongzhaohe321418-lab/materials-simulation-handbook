"""Phonon dispersion of a 1D diatomic chain.

Extracted from docs/ch03b-solid-state/05-phonons.md (S 3b.5.5).
"""
from __future__ import annotations
import numpy as np
import numpy.typing as npt

M1: float = 1.0      # mass of atom A
M2: float = 3.0      # mass of atom B (heavier)
K_SPRING: float = 1.0
A_LATT: float = 1.0


def dynamical_matrix(
    k: float,
    m1: float = M1,
    m2: float = M2,
    K: float = K_SPRING,
    a: float = A_LATT,
) -> npt.NDArray[np.complex128]:
    """Return the 2x2 dynamical matrix at wavevector k."""
    D = np.zeros((2, 2), dtype=np.complex128)
    D[0, 0] = (2 * K) / m1
    D[1, 1] = (2 * K) / m2
    coupling = -K * (1.0 + np.exp(-1j * k * a))
    D[0, 1] = coupling / np.sqrt(m1 * m2)
    D[1, 0] = np.conj(D[0, 1])
    return D


def phonon_frequencies(
    k: float,
    m1: float = M1,
    m2: float = M2,
    K: float = K_SPRING,
    a: float = A_LATT,
) -> npt.NDArray[np.float64]:
    """Sqrt of eigenvalues (clipped at 0) of the dynamical matrix."""
    eigs = np.linalg.eigvalsh(dynamical_matrix(k, m1, m2, K, a))
    return np.sqrt(np.clip(eigs, a_min=0.0, a_max=None))


def analytical_branches(
    k: npt.NDArray[np.float64],
    m1: float = M1,
    m2: float = M2,
    K: float = K_SPRING,
    a: float = A_LATT,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    pref = K * (m1 + m2) / (m1 * m2)
    discr = 1.0 - (4 * m1 * m2 / (m1 + m2) ** 2) * np.sin(k * a / 2) ** 2
    omega_minus = np.sqrt(pref * (1.0 - np.sqrt(discr)))
    omega_plus = np.sqrt(pref * (1.0 + np.sqrt(discr)))
    return omega_minus, omega_plus


def compute_bands(
    n_k: int = 200,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Return (k_values, omega) where omega has shape (n_k, 2)."""
    k_vals = np.linspace(-np.pi / A_LATT, np.pi / A_LATT, n_k)
    omega = np.array([phonon_frequencies(float(k)) for k in k_vals])
    return k_vals, omega


if __name__ == "__main__":
    k_vals, omega = compute_bands(200)
    om_m, om_p = analytical_branches(k_vals)
    print(f"max |num - ana| acoustic: {np.max(np.abs(omega[:, 0] - om_m)):.2e}")
    print(f"max |num - ana| optical:  {np.max(np.abs(omega[:, 1] - om_p)):.2e}")
    print(f"omega- at k=0: {omega[len(k_vals)//2, 0]:.3e}")
