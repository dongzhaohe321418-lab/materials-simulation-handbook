"""Acquisition functions (EI, UCB) and a short BO loop.

Extracted from docs/ch11-active/03-acquisition.md (S 11.3.7).
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.stats import norm

from .gp import GP


def expected_improvement(
    gp: GP,
    X: NDArray[np.float64],
    f_best: float,
    xi: float = 0.01,
) -> NDArray[np.float64]:
    """Expected improvement for a maximisation objective."""
    mu, var = gp.predict(X)
    sigma = np.sqrt(np.maximum(var, 1e-12))
    z = (mu - f_best - xi) / sigma
    ei = sigma * (z * norm.cdf(z) + norm.pdf(z))
    ei[sigma < 1e-9] = 0.0
    return ei


def upper_confidence_bound(
    gp: GP,
    X: NDArray[np.float64],
    kappa: float = 2.0,
) -> NDArray[np.float64]:
    """UCB acquisition."""
    mu, var = gp.predict(X)
    return mu + kappa * np.sqrt(np.maximum(var, 1e-12))


def bo_loop_sin(
    n_iter: int = 10,
    noise_level: float = 0.05,
    seed: int = 0,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Run a BO loop maximising sin(x) on [0, 7]."""
    rng = np.random.default_rng(seed)
    true_f = np.sin

    X_train = np.array([[1.5], [3.0], [5.0]])
    y_train = true_f(X_train.flatten()) + rng.normal(0, noise_level, size=3)

    X_cand = np.linspace(0, 7, 500).reshape(-1, 1)

    for _ in range(n_iter):
        gp = GP()
        gp.optimise_hyperparameters(X_train, y_train)
        f_best = float(np.max(y_train))
        ei = expected_improvement(gp, X_cand, f_best, xi=0.01)
        x_next = X_cand[int(np.argmax(ei))].reshape(1, -1)
        y_next = true_f(float(x_next)) + rng.normal(0, noise_level)
        X_train = np.vstack([X_train, x_next])
        y_train = np.append(y_train, y_next)
    return X_train, y_train


if __name__ == "__main__":
    X, y = bo_loop_sin(n_iter=10)
    print(f"best y = {y.max():.3f} at x = {float(X[int(np.argmax(y))]):.3f}")
