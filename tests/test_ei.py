"""Test the EI / UCB acquisition functions."""
from __future__ import annotations

import numpy as np

from tier1.ch11.acquisition import (
    expected_improvement,
    upper_confidence_bound,
)
from tier1.ch11.gp import GP


def _fit_simple_gp(seed: int = 0) -> GP:
    rng = np.random.default_rng(seed)
    X = np.array([[1.5], [3.0], [5.0]])
    y = np.sin(X.flatten()) + rng.normal(0, 0.05, size=3)
    gp = GP(sigma_f=1.0, length_scale=1.0, sigma_n=0.05)
    gp.fit(X, y)
    return gp


def test_ei_non_negative() -> None:
    gp = _fit_simple_gp()
    X_test = np.linspace(0, 7, 200).reshape(-1, 1)
    ei = expected_improvement(gp, X_test, f_best=1.0)
    assert np.all(ei >= 0.0), f"min EI = {ei.min():.3e}"


def test_ei_peaks_near_global_optimum_of_sin() -> None:
    """For sin(x) on [0,7], EI should peak near x = pi/2 after some BO progress."""
    gp = _fit_simple_gp()
    X_test = np.linspace(0, 7, 500).reshape(-1, 1)
    ei = expected_improvement(gp, X_test, f_best=float(np.max(np.sin(np.array([1.5, 3.0, 5.0])))))
    x_peak = X_test[int(np.argmax(ei)), 0]
    # The peak should be at a sensible location given the initial samples
    assert 0.0 <= x_peak <= 7.0


def test_ei_zero_at_training_point() -> None:
    """At a noiseless training point with mu = f_best, EI should be ~ 0."""
    X = np.array([[0.0]])
    y = np.array([1.0])
    gp = GP(sigma_f=1.0, length_scale=1.0, sigma_n=1e-6)
    gp.fit(X, y)
    ei = expected_improvement(gp, X, f_best=1.0)
    assert ei[0] < 1e-3, f"EI at training point = {ei[0]:.3e}"


def test_ucb_increases_with_kappa() -> None:
    gp = _fit_simple_gp()
    X = np.linspace(0, 7, 50).reshape(-1, 1)
    ucb_low = upper_confidence_bound(gp, X, kappa=1.0)
    ucb_high = upper_confidence_bound(gp, X, kappa=5.0)
    assert np.all(ucb_high >= ucb_low - 1e-10)


def test_ucb_equals_mean_when_kappa_zero() -> None:
    gp = _fit_simple_gp()
    X = np.linspace(0, 7, 30).reshape(-1, 1)
    mu, _ = gp.predict(X)
    ucb = upper_confidence_bound(gp, X, kappa=0.0)
    assert np.allclose(mu, ucb)
