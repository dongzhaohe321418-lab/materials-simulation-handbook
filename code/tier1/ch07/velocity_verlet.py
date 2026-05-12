"""Velocity-Verlet integration of a 1D harmonic oscillator.

Extracted from docs/ch07-md/01-integration.md.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import numpy as np


@dataclass
class Oscillator:
    """1D harmonic oscillator U(x) = (1/2) k x^2."""
    mass: float = 1.0
    k: float = 1.0

    def force(self, x: float) -> float:
        return -self.k * x

    def potential(self, x: float) -> float:
        return 0.5 * self.k * x * x

    @property
    def omega(self) -> float:
        return float(np.sqrt(self.k / self.mass))


def velocity_verlet(
    x0: float,
    v0: float,
    force_fn: Callable[[float], float],
    mass: float,
    dt: float,
    n_steps: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Integrate Newton's equations with velocity-Verlet."""
    t = np.zeros(n_steps + 1, dtype=np.float64)
    x = np.zeros(n_steps + 1, dtype=np.float64)
    v = np.zeros(n_steps + 1, dtype=np.float64)
    x[0], v[0] = x0, v0
    a = force_fn(x0) / mass
    for n in range(n_steps):
        v_half = v[n] + 0.5 * dt * a
        x[n + 1] = x[n] + dt * v_half
        a_new = force_fn(x[n + 1]) / mass
        v[n + 1] = v_half + 0.5 * dt * a_new
        a = a_new
        t[n + 1] = t[n] + dt
    return t, x, v


if __name__ == "__main__":
    osc = Oscillator()
    t, x, v = velocity_verlet(1.0, 0.0, osc.force, osc.mass, 0.05, 1000)
    E = 0.5 * osc.mass * v ** 2 + osc.potential(x)
    print(f"E[0] = {E[0]:.10f}  E[-1] = {E[-1]:.10f}")
    print(f"spread = {(E.max() - E.min()) / E[0]:.2e}")
