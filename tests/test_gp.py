"""Test the NumPy Gaussian process regression."""
from __future__ import annotations

import numpy as np

from tier1.ch11.gp import GP, rbf_kernel


def test_rbf_kernel_self_is_sigma_f_squared() -> None:
    x = np.array([[0.0], [1.0]])
    K = rbf_kernel(x, x, sigma_f=2.0, length_scale=1.0)
    assert abs(K[0, 0] - 4.0) < 1e-10
    assert abs(K[1, 1] - 4.0) < 1e-10


def test_rbf_kernel_symmetric() -> None:
    rng = np.random.default_rng(0)
    x = rng.normal(size=(5, 2))
    K = rbf_kernel(x, x, 1.0, 1.0)
    assert np.allclose(K, K.T)


def test_gp_interpolates_noiseless_training_points() -> None:
    X = np.array([[0.0], [1.0], [2.0], [3.0]])
    y = np.sin(X.flatten())
    gp = GP(sigma_f=1.0, length_scale=1.0, sigma_n=1e-6)
    gp.fit(X, y)
    mu, var = gp.predict(X)
    assert np.allclose(mu, y, atol=1e-3), f"|mu - y| max = {np.max(np.abs(mu - y)):.2e}"
    # Variance at training points should be near zero (essentially the noise).
    assert np.all(var < 1e-3)


def test_gp_uncertainty_grows_away_from_data() -> None:
    X = np.array([[0.0], [1.0]])
    y = np.array([0.0, 1.0])
    gp = GP(sigma_f=1.0, length_scale=0.3, sigma_n=0.01)
    gp.fit(X, y)
    _, var_near = gp.predict(np.array([[0.5]]))
    _, var_far = gp.predict(np.array([[10.0]]))
    assert var_far[0] > var_near[0]


def test_gp_hyperparameter_optimisation_runs() -> None:
    rng = np.random.default_rng(42)
    X = np.linspace(0, 5, 10).reshape(-1, 1)
    y = np.sin(X.flatten()) + rng.normal(0, 0.05, size=10)
    gp = GP()
    gp.optimise_hyperparameters(X, y, n_restarts=3)
    assert gp.length_scale > 0
    assert gp.sigma_f > 0
    assert gp.sigma_n > 0
