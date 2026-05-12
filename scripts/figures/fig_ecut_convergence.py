#!/usr/bin/env python3
"""Synthetic plane-wave cutoff (E_cut) convergence curve."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
ecuts = np.array([20, 30, 40, 50, 60, 70, 80, 100, 120])
E_inf = -157.4324
# Monotonic exponential convergence reaching < 1 meV/atom near 80 Ry
E = E_inf + 0.05 * np.exp(-(ecuts-20)/12.0)
dE_mev = np.abs((E - E_inf) * 1000 * 13.6057)

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(ecuts, dE_mev, "-o", color=COLOURS["blue"], lw=2.2, ms=8,
        label=r"$|E_{\rm tot}(E_{\rm cut}) - E_{\infty}|$")
ax.axhline(1, color=COLOURS["green"], ls="--", lw=1.6,
           label="1 meV/atom tolerance")
ax.axvline(100, color=COLOURS["grey"], ls=":", lw=1.2)
ax.annotate("converged cutoff", xy=(100, 1.0), xytext=(60, 8), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
ax.set_xlabel(r"$E_{\rm cut}$ (Ry)")
ax.set_ylabel(r"$|E - E_{\infty}|$ (meV/atom)")
ax.set_title("Plane-wave cutoff convergence (synthetic example)")
ax.set_yscale("log")
ax.set_ylim(0.1, 1000)
ax.legend(frameon=False, loc="upper right")
ax.grid(alpha=0.3, linestyle="--", which="both")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch06/fig_ecut_convergence.png")
save(fig, out)
