"""Test MSD, RDF, and VACF helpers."""
from __future__ import annotations

import numpy as np

from tier1.ch07.trajectory_analysis import (
    diffusion_coefficient,
    msd,
    rdf,
    vacf_and_vdos,
)


def test_msd_zero_lag_is_zero() -> None:
    rng = np.random.default_rng(0)
    pos = rng.normal(size=(50, 10, 3))
    out = msd(pos)
    assert out[0] == 0.0


def test_msd_brownian_slope_gives_diffusion() -> None:
    """For Brownian motion with step variance s^2 per frame and dt=1,
    MSD(t) = 3 * s^2 * t in 3D so D = slope/6 = s^2/2."""
    rng = np.random.default_rng(42)
    n_frames, n_atoms = 400, 200
    step_var = 0.04
    steps = rng.normal(0.0, np.sqrt(step_var), size=(n_frames, n_atoms, 3))
    positions = np.cumsum(steps, axis=0)
    t = np.arange(n_frames, dtype=np.float64)
    m = msd(positions)
    D = diffusion_coefficient(t, m, (10.0, 200.0))
    D_expected = step_var / 2.0
    rel_err = abs(D - D_expected) / D_expected
    assert rel_err < 0.10, f"D = {D:.4f} vs {D_expected:.4f} rel_err {rel_err:.2e}"


def test_msd_monotonic_for_random_walk() -> None:
    rng = np.random.default_rng(1)
    steps = rng.normal(size=(100, 30, 3)) * 0.1
    pos = np.cumsum(steps, axis=0)
    m = msd(pos)
    # Allow tiny non-monotonic blips from finite statistics in early lag
    # but the overall trend must be increasing.
    assert m[-1] > m[10]


def test_rdf_uniform_random_approaches_one() -> None:
    """An ideal-gas configuration should give g(r) ~ 1 at intermediate r."""
    rng = np.random.default_rng(7)
    box = np.array([10.0, 10.0, 10.0])
    n_atoms = 800
    positions = rng.uniform(0.0, 10.0, size=(5, n_atoms, 3))
    r, g = rdf(positions, box, r_max=4.0, n_bins=40)
    # Look at the middle of the r range away from r=0 noise
    mid = g[5:35]
    assert abs(mid.mean() - 1.0) < 0.1, f"mean g(r) = {mid.mean():.3f}"


def test_vacf_first_value_is_one() -> None:
    rng = np.random.default_rng(0)
    vel = rng.normal(size=(64, 20, 3))
    _, acf, _ = vacf_and_vdos(vel, dt=1.0)
    assert abs(acf[0] - 1.0) < 1e-10
