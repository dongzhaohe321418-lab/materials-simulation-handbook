"""Trajectory analysis helpers: MSD, RDF, VACF/VDOS.

Extracted from docs/ch07-md/06-analysis.md.
"""
from __future__ import annotations
import numpy as np


def msd(positions: np.ndarray) -> np.ndarray:
    """Time-origin-averaged MSD.

    positions: (n_frames, n_atoms, 3) of unwrapped positions.
    Returns array of length n_frames.
    """
    n_frames = positions.shape[0]
    msd_out = np.zeros(n_frames, dtype=np.float64)
    for tau in range(n_frames):
        diffs = positions[tau:] - positions[: n_frames - tau]
        msd_out[tau] = np.mean(np.sum(diffs ** 2, axis=-1))
    return msd_out


def diffusion_coefficient(
    t: np.ndarray, msd_t: np.ndarray, fit_range: tuple[float, float]
) -> float:
    """Linear-regression D = slope/6 over a chosen time window."""
    mask = (t >= fit_range[0]) & (t <= fit_range[1])
    slope, _ = np.polyfit(t[mask], msd_t[mask], 1)
    return slope / 6.0


def rdf(
    positions: np.ndarray,
    box: np.ndarray,
    r_max: float,
    n_bins: int = 200,
) -> tuple[np.ndarray, np.ndarray]:
    """Radial distribution function under minimum-image convention.

    positions: (n_frames, n_atoms, 3) wrapped or unwrapped.
    box: (3,) orthorhombic box dimensions.
    """
    n_frames, n_atoms, _ = positions.shape
    rho = n_atoms / np.prod(box)
    dr = r_max / n_bins
    hist = np.zeros(n_bins, dtype=np.int64)

    for frame in positions:
        for i in range(n_atoms - 1):
            dvec = frame[i + 1:] - frame[i]
            dvec -= box * np.round(dvec / box)
            d = np.linalg.norm(dvec, axis=-1)
            d = d[d < r_max]
            idx = (d / dr).astype(np.int64)
            np.add.at(hist, idx, 1)

    r = (np.arange(n_bins) + 0.5) * dr
    shell_vol = 4 * np.pi * r ** 2 * dr
    norm = shell_vol * rho * n_atoms * n_frames / 2.0
    g = hist / norm
    return r, g


def vacf_and_vdos(
    velocities: np.ndarray, dt: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Velocity autocorrelation and vibrational DOS via FFT.

    velocities: (n_frames, n_atoms, 3).
    """
    n_frames = velocities.shape[0]
    v_flat = velocities.reshape(n_frames, -1)
    F = np.fft.rfft(v_flat, axis=0, n=2 * n_frames)
    acf = np.fft.irfft(F * np.conj(F), axis=0)[:n_frames]
    acf = acf.sum(axis=1)
    norm = (n_frames - np.arange(n_frames)) * v_flat.shape[1]
    acf = acf / norm
    acf /= acf[0]
    t = np.arange(n_frames) * dt
    vdos = np.abs(np.fft.rfft(acf))
    return t, acf, vdos


if __name__ == "__main__":
    rng = np.random.default_rng(0)
    n_frames, n_atoms = 200, 50
    # Synthetic Brownian-like walk
    steps = rng.normal(0.0, 0.1, size=(n_frames, n_atoms, 3))
    positions = np.cumsum(steps, axis=0)
    m = msd(positions)
    t = np.arange(n_frames) * 1.0
    D = diffusion_coefficient(t, m, (10.0, 150.0))
    print(f"diffusion D = {D:.4f}")
