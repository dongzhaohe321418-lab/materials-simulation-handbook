#!/usr/bin/env python3
"""Synthetic total energy vs SCF iteration showing exponential convergence."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
n = 24
it = np.arange(1, n+1)
E_final = -157.4324
# Faster exponential decay so |dE| reaches < 1e-6 by iteration ~20
E = E_final + 0.5 * np.exp(-it/1.8) * (1 + 0.05*np.cos(it/1.5))
dE = np.abs(np.diff(E, prepend=E[0]+1))
dE[0] = 0.5

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].plot(it, E, "-o", color=COLOURS["blue"], lw=2, ms=6)
axes[0].axhline(E_final, color=COLOURS["red"], ls="--", lw=1.4,
                label=fr"converged: $E_{{\rm tot}} = {E_final:.4f}$ Ry")
axes[0].set_xlabel("SCF iteration")
axes[0].set_ylabel(r"$E_{\rm tot}$ (Ry)")
axes[0].set_title("Total energy vs SCF iteration")
axes[0].legend(frameon=False, loc="upper right")
axes[0].grid(alpha=0.3, linestyle="--")

axes[1].semilogy(it, dE, "-o", color=COLOURS["red"], lw=2, ms=6)
axes[1].axhline(1e-6, color=COLOURS["green"], ls="--", lw=1.4,
                label=r"$\mathrm{conv\_thr} = 10^{-6}$ Ry")
axes[1].annotate("convergence reached", xy=(20, 1e-6), xytext=(8, 1e-4),
                 fontsize=11,
                 arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
axes[1].set_xlabel("SCF iteration")
axes[1].set_ylabel(r"$|\Delta E|$ (Ry, log scale)")
axes[1].set_title("Energy change per step")
axes[1].legend(frameon=False, loc="upper right")
axes[1].grid(alpha=0.3, linestyle="--", which="both")

fig.suptitle("SCF convergence (synthetic example)", fontsize=15, y=1.02)
fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch05/fig_scf_convergence.png")
save(fig, out)
