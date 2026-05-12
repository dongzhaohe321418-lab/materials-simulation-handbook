#!/usr/bin/env python3
"""Lennard-Jones 12-6 potential vs r, showing minimum at r=r_min."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
sigma, eps = 1.0, 1.0
r = np.linspace(0.92, 3.0, 600)
V = 4*eps*((sigma/r)**12 - (sigma/r)**6)
r_min = 2**(1/6) * sigma

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(r, V, color=COLOURS["blue"], lw=2.2, label=r"$V_{LJ}(r) = 4\varepsilon\,[(\sigma/r)^{12} - (\sigma/r)^{6}]$")
ax.axhline(0, color="k", lw=0.6)
ax.axhline(-eps, color=COLOURS["grey"], lw=0.8, ls="--")
ax.axvline(sigma, color=COLOURS["orange"], lw=0.8, ls=":")
ax.axvline(r_min, color=COLOURS["green"], lw=0.8, ls=":")
ax.plot([r_min], [-eps], 'o', color=COLOURS["green"], ms=8, zorder=5)

ax.annotate(r"minimum: $r_{\min}=2^{1/6}\sigma$", xy=(r_min, -eps),
            xytext=(1.55, -0.55), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="k", lw=1.0))
ax.annotate(r"zero crossing: $V(\sigma)=0$", xy=(sigma, 0),
            xytext=(1.25, 0.55), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="k", lw=1.0))
ax.text(2.55, -eps-0.08, r"well depth $-\varepsilon$",
        fontsize=12, color=COLOURS["grey"])

ax.set_xlabel(r"$r / \sigma$")
ax.set_ylabel(r"$V(r) / \varepsilon$")
ax.set_title("Lennard-Jones 12-6 potential")
ax.set_ylim(-1.25, 1.5)
ax.set_xlim(0.9, 3.0)
ax.grid(alpha=0.3, linestyle="--")
ax.legend(frameon=False, loc="upper right")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch07/fig_lj_potential.png")
save(fig, out)
