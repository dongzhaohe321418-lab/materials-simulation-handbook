"""Test the velocity-Verlet integrator on a 1D harmonic oscillator."""
from __future__ import annotations

import numpy as np

from tier1.ch07.velocity_verlet import Oscillator, velocity_verlet


def test_energy_conservation_over_1000_steps() -> None:
    osc = Oscillator(mass=1.0, k=1.0)
    t, x, v = velocity_verlet(
        x0=1.0, v0=0.0, force_fn=osc.force, mass=osc.mass, dt=0.05, n_steps=1000
    )
    E = 0.5 * osc.mass * v ** 2 + osc.potential(x)
    drift = (E.max() - E.min()) / E[0]
    assert drift < 1e-3, f"energy spread {drift:.2e} exceeds 1e-3"


def test_no_secular_drift() -> None:
    """A symplectic integrator should not drift monotonically."""
    osc = Oscillator()
    _, x, v = velocity_verlet(1.0, 0.0, osc.force, osc.mass, 0.05, 2000)
    E = 0.5 * osc.mass * v ** 2 + osc.potential(x)
    # The mean of the second half should be close to the mean of the first half.
    half = len(E) // 2
    mean1 = E[:half].mean()
    mean2 = E[half:].mean()
    assert abs(mean2 - mean1) / E[0] < 1e-4


def test_oscillation_period() -> None:
    """The recovered period should be 2 pi / omega for omega=1."""
    osc = Oscillator(mass=1.0, k=1.0)
    t, x, _ = velocity_verlet(1.0, 0.0, osc.force, osc.mass, 0.01, 2000)
    # Find first negative-going zero crossing after the start.
    # x starts at +1, so first zero crossing is at t = pi/2.
    sign_change = np.where(np.diff(np.signbit(x)))[0]
    assert len(sign_change) > 0
    t_zero = t[sign_change[0]]
    expected = np.pi / 2
    assert abs(t_zero - expected) / expected < 0.05


def test_returns_correct_shapes() -> None:
    osc = Oscillator()
    n = 100
    t, x, v = velocity_verlet(0.5, 0.0, osc.force, osc.mass, 0.1, n)
    assert t.shape == (n + 1,)
    assert x.shape == (n + 1,)
    assert v.shape == (n + 1,)
