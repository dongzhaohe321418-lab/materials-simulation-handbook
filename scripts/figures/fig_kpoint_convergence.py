#!/usr/bin/env python3
"""Synthetic k-grid convergence: E vs N for NxNxN Monkhorst-Pack."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
Ns = np.arange(2, 17)
E_inf = -157.4324
E = E_inf + 0.02/Ns**2 * (1 + 0.6*np.cos(Ns*1.3))
dE_mev = np.abs((E - E_inf) * 1000 * 13.6057)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(Ns, dE_mev, "-o", color=COLOURS["blue"], lw=2.2, ms=8,
        label=r"$|E_{\rm tot}(N) - E_{\infty}|$")
ax.axhline(1, color=COLOURS["green"], ls="--", lw=1.6,
           label="1 meV/atom tolerance")
ax.axvline(10, color=COLOURS["grey"], ls=":", lw=1.2)
ax.annotate("converged grid\n($\\approx$ 1 meV/atom)", xy=(10, 4),
            xytext=(3.5, 0.9),
            fontsize=11, arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
ax.set_xlabel(r"$N$ for $N \times N \times N$ Monkhorst-Pack grid")
ax.set_ylabel(r"$|E - E_{\infty}|$ (meV/atom)")
ax.set_title("$\\mathbf{k}$-point convergence (synthetic example)")
ax.set_yscale("log")
ax.legend(frameon=False, loc="upper right")
ax.grid(alpha=0.3, linestyle="--", which="both")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch06/fig_kpoint_convergence.png")
save(fig, out)
