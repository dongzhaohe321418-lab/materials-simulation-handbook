#!/usr/bin/env python3
"""GP posterior + Expected Improvement (EI) acquisition function in 1D."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

setup()

def rbf(x1, x2, l=1.0, sf=1.0):
    d = x1[:, None] - x2[None, :]
    return sf**2 * np.exp(-0.5*(d/l)**2)

def f(x):
    return np.sin(1.4*x) + 0.3*np.cos(3*x) + 0.05*x

x = np.linspace(-4, 4, 500)
xd = np.array([-3.5, -1.5, 0.7, 2.5])
yd = f(xd)
Kxx = rbf(xd, xd) + 1e-6*np.eye(len(xd))
Kxs = rbf(xd, x)
alpha = np.linalg.solve(Kxx, yd)
mu = Kxs.T @ alpha
v = np.linalg.solve(Kxx, Kxs)
cov = rbf(x, x) - Kxs.T @ v
sd = np.sqrt(np.maximum(np.diag(cov), 1e-12))

f_best = yd.min()
imp = f_best - mu - 0.01
Z = imp / sd
ei = imp*norm.cdf(Z) + sd*norm.pdf(Z)
ei = np.maximum(ei, 0)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                gridspec_kw=dict(height_ratios=[2, 1]))
ax1.plot(x, f(x), color="k", lw=1.2, ls=":", label=r"true $f$ (hidden)")
ax1.fill_between(x, mu-2*sd, mu+2*sd, color=COLOURS["blue"], alpha=0.2,
                 label=r"$\pm 2\sigma$")
ax1.plot(x, mu, color=COLOURS["blue"], lw=2.2, label=r"posterior mean")
ax1.scatter(xd, yd, color="k", s=60, zorder=5, label="observations")
ax1.axhline(f_best, color=COLOURS["red"], lw=1.2, ls="--",
            label=r"$f_{\rm best}$ so far")
ax1.set_ylabel(r"$f(x)$")
ax1.set_title("GP posterior")
ax1.legend(frameon=False, loc="upper right", fontsize=10)
ax1.grid(alpha=0.3, linestyle="--")

ax2.fill_between(x, 0, ei, color=COLOURS["green"], alpha=0.4)
ax2.plot(x, ei, color=COLOURS["green"], lw=2.2, label=r"$\mathrm{EI}(x)$")
x_next = x[np.argmax(ei)]
ax2.axvline(x_next, color=COLOURS["red"], lw=1.6, ls="--",
            label=fr"next query $x^* = {x_next:.2f}$")
ax2.annotate("argmax of EI", xy=(x_next, ei.max()), xytext=(x_next+0.5, ei.max()*0.7),
             fontsize=11, arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
ax2.set_ylabel(r"$\mathrm{EI}(x)$")
ax2.set_xlabel(r"$x$")
ax2.set_title("Expected Improvement acquisition")
ax2.legend(frameon=False, loc="upper right", fontsize=11)
ax2.grid(alpha=0.3, linestyle="--")

fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch11/fig_ei_acquisition.png")
save(fig, out)
