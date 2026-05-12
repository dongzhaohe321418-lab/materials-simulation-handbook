#!/usr/bin/env python3
"""GP prior samples and posterior given 3 noiseless data points."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS, PALETTE
import numpy as np
import matplotlib.pyplot as plt

setup()
rng = np.random.default_rng(0)

def rbf(x1, x2, l=1.0, sf=1.0):
    d = x1[:, None] - x2[None, :]
    return sf**2 * np.exp(-0.5*(d/l)**2)

x = np.linspace(-5, 5, 250)
K = rbf(x, x, l=1.0) + 1e-8*np.eye(len(x))
L = np.linalg.cholesky(K)

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
for i in range(5):
    sample = L @ rng.standard_normal(len(x))
    axes[0].plot(x, sample, lw=1.5, alpha=0.85, color=PALETTE[i])
axes[0].fill_between(x, -2, 2, color=COLOURS["grey"], alpha=0.15,
                     label=r"$\pm 2\sigma$ prior")
axes[0].set_title(r"GP prior: samples from $\mathcal{N}(0, K)$ (RBF kernel)")
axes[0].set_xlabel(r"$x$"); axes[0].set_ylabel(r"$f(x)$")
axes[0].set_xlim(-5, 5); axes[0].set_ylim(-3, 3)
axes[0].legend(frameon=False, loc="upper right")
axes[0].grid(alpha=0.3, linestyle="--")

xd = np.array([-3.0, 0.5, 2.5])
yd = np.array([-1.0, 1.3, -0.4])
Kxx = rbf(xd, xd) + 1e-6*np.eye(3)
Kxs = rbf(xd, x)
alpha = np.linalg.solve(Kxx, yd)
mu = Kxs.T @ alpha
v = np.linalg.solve(Kxx, Kxs)
cov = rbf(x, x) - Kxs.T @ v
sd = np.sqrt(np.maximum(np.diag(cov), 0.0))

axes[1].fill_between(x, mu-2*sd, mu+2*sd, color=COLOURS["blue"], alpha=0.2,
                     label=r"$\pm 2\sigma$ posterior")
axes[1].plot(x, mu, color=COLOURS["blue"], lw=2.2, label=r"posterior mean")
Lpost = np.linalg.cholesky(cov + 1e-8*np.eye(len(x)))
for i in range(3):
    axes[1].plot(x, mu + Lpost @ rng.standard_normal(len(x)),
                 lw=1.0, alpha=0.6, color=COLOURS["red"])
axes[1].scatter(xd, yd, color="k", s=60, zorder=5, label="observations")

axes[1].annotate("uncertainty\nshrinks at data", xy=(xd[1], 1.3),
                 xytext=(-2.4, 2.3), fontsize=11,
                 arrowprops=dict(arrowstyle="->", color="k", lw=0.9))

axes[1].set_title(r"GP posterior given 3 noiseless observations")
axes[1].set_xlabel(r"$x$"); axes[1].set_ylabel(r"$f(x)$")
axes[1].set_xlim(-5, 5); axes[1].set_ylim(-3, 3)
axes[1].legend(frameon=False, loc="lower right")
axes[1].grid(alpha=0.3, linestyle="--")

fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch11/fig_gaussian_kernel.png")
save(fig, out)
