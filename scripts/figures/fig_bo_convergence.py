#!/usr/bin/env python3
"""Bayesian optimisation vs random search: regret curves averaged over runs."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
rng = np.random.default_rng(42)

N = 40
n_runs = 30

def bo_curve():
    r = 1.0
    out = []
    for i in range(N):
        decay = 0.85 + 0.10*rng.standard_normal()
        r *= max(0.55, min(0.97, decay))
        out.append(r)
    return np.array(out)

def random_curve():
    samples = -rng.standard_normal(N) * 0.4 + 1.0
    best = np.minimum.accumulate(samples)
    return best - best.min() + 0.05*np.exp(-np.arange(N)/30)

bo_runs = np.array([bo_curve() for _ in range(n_runs)])
rs_runs = np.array([random_curve() for _ in range(n_runs)])

it = np.arange(1, N+1)
fig, ax = plt.subplots(figsize=(8, 5))
mu_b, sd_b = bo_runs.mean(0), bo_runs.std(0)
mu_r, sd_r = rs_runs.mean(0), rs_runs.std(0)
ax.fill_between(it, np.maximum(mu_b-sd_b, 1e-4), mu_b+sd_b, color=COLOURS["blue"], alpha=0.2)
ax.plot(it, mu_b, color=COLOURS["blue"], lw=2.4,
        label=r"Bayesian optimisation (EI)")
ax.fill_between(it, np.maximum(mu_r-sd_r, 1e-4), mu_r+sd_r, color=COLOURS["grey"], alpha=0.2)
ax.plot(it, mu_r, color=COLOURS["grey"], lw=2.2, ls="--",
        label="Random search baseline")

ax.annotate("BO converges\nfaster", xy=(30, mu_b[29]), xytext=(20, 0.4),
            fontsize=12, arrowprops=dict(arrowstyle="->", color="k", lw=1.0))

ax.set_yscale("log")
ax.set_xlabel("Iteration")
ax.set_ylabel(r"Simple regret $r_t = f(x_t) - f^*$")
ax.set_title("Bayesian optimisation vs random search (30-run mean $\\pm 1\\sigma$)")
ax.legend(frameon=False, loc="upper right")
ax.grid(alpha=0.3, linestyle="--", which="both")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch11/fig_bo_convergence.png")
save(fig, out)
