"""Minimal Gaussian process regression in pure NumPy.

Extracted from docs/ch11-active/02-gp.md (S 11.2.6).
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize


def rbf_kernel(
    x1: NDArray[np.float64],
    x2: NDArray[np.float64],
    sigma_f: float,
    length_scale: float,
) -> NDArray[np.float64]:
    """RBF kernel between two sets of inputs of shape (n, d)."""
    sq = (
        np.sum(x1 ** 2, axis=1, keepdims=True)
        + np.sum(x2 ** 2, axis=1)
        - 2.0 * x1 @ x2.T
    )
    sq = np.maximum(sq, 0.0)
    return sigma_f ** 2 * np.exp(-0.5 * sq / length_scale ** 2)


class GP:
    """Gaussian process regression with RBF kernel and Gaussian noise."""

    def __init__(
        self,
        sigma_f: float = 1.0,
        length_scale: float = 1.0,
        sigma_n: float = 0.1,
    ) -> None:
        self.sigma_f = sigma_f
        self.length_scale = length_scale
        self.sigma_n = sigma_n
        self.X: NDArray[np.float64] | None = None
        self.y: NDArray[np.float64] | None = None
        self.alpha: NDArray[np.float64] | None = None
        self.L: NDArray[np.float64] | None = None

    def _prepare_X(self, X: NDArray[np.float64]) -> NDArray[np.float64]:
        X = np.atleast_2d(X)
        return X

    def fit(
        self, X: NDArray[np.float64], y: NDArray[np.float64]
    ) -> None:
        X = self._prepare_X(X)
        if X.shape[0] != y.shape[0]:
            X = X.T
        self.X = X
        self.y = y.astype(np.float64)
        K = rbf_kernel(self.X, self.X, self.sigma_f, self.length_scale)
        K += self.sigma_n ** 2 * np.eye(len(self.X))
        self.L, lower = cho_factor(K, lower=True)
        self._lower = lower
        self.alpha = cho_solve((self.L, lower), self.y)

    def predict(
        self, X_star: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        assert self.X is not None and self.alpha is not None
        X_star = np.atleast_2d(X_star)
        if X_star.shape[1] != self.X.shape[1]:
            X_star = X_star.T
        K_star = rbf_kernel(self.X, X_star, self.sigma_f, self.length_scale)
        mu = K_star.T @ self.alpha
        v = cho_solve((self.L, self._lower), K_star)
        K_starstar = rbf_kernel(
            X_star, X_star, self.sigma_f, self.length_scale
        )
        sigma2 = np.diag(K_starstar) - np.sum(K_star * v, axis=0)
        sigma2 = np.maximum(sigma2, 1e-12)
        return mu, sigma2

    def log_marginal_likelihood(
        self,
        params: NDArray[np.float64],
        X: NDArray[np.float64],
        y: NDArray[np.float64],
    ) -> float:
        sigma_f, length_scale, sigma_n = np.exp(params)
        K = rbf_kernel(X, X, sigma_f, length_scale)
        K += sigma_n ** 2 * np.eye(len(X))
        try:
            L, lower = cho_factor(K, lower=True)
        except np.linalg.LinAlgError:
            return 1e10
        alpha = cho_solve((L, lower), y)
        nll = 0.5 * y @ alpha
        nll += np.sum(np.log(np.diag(L)))
        nll += 0.5 * len(X) * np.log(2 * np.pi)
        return float(nll)

    def optimise_hyperparameters(
        self,
        X: NDArray[np.float64],
        y: NDArray[np.float64],
        n_restarts: int = 5,
        seed: int = 0,
    ) -> None:
        X = np.atleast_2d(X)
        if X.shape[0] != y.shape[0]:
            X = X.T
        best_nll = np.inf
        best_params = None
        rng = np.random.default_rng(seed)
        for _ in range(n_restarts):
            init = rng.normal(0.0, 1.0, size=3)
            res = minimize(
                self.log_marginal_likelihood,
                init,
                args=(X, y),
                method="L-BFGS-B",
            )
            if res.fun < best_nll:
                best_nll = res.fun
                best_params = res.x
        assert best_params is not None
        self.sigma_f, self.length_scale, self.sigma_n = np.exp(best_params)
        self.fit(X, y)


if __name__ == "__main__":
    rng = np.random.default_rng(42)
    X_train = np.array([[0.5], [1.0], [2.0], [3.5], [4.0], [5.5]])
    y_train = np.sin(X_train.flatten()) + rng.normal(0, 0.1, size=6)
    gp = GP()
    gp.optimise_hyperparameters(X_train, y_train)
    print(f"sigma_f={gp.sigma_f:.3f}  L={gp.length_scale:.3f}  noise={gp.sigma_n:.3f}")
    mu, var = gp.predict(np.array([[1.0]]))
    print(f"mu(1.0) = {mu[0]:.3f}  var = {var[0]:.3e}")
