#!/usr/bin/env python3
"""Comparison of Dulong-Petit, Einstein, and Debye heat capacities."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
T_TD = np.linspace(0.02, 2.0, 300)  # T / Theta_D

def debye_int(x_max, n=400):
    x = np.linspace(1e-6, x_max, n)
    y = x**4 * np.exp(x) / (np.exp(x)-1)**2
    return np.trapz(y, x)

Cv_debye = np.array([9 * T**3 * debye_int(1.0/T) for T in T_TD])
x = 1.0 / T_TD
Cv_einstein = (x**2 * np.exp(x)) / (np.exp(x)-1)**2
Cv_DP = np.ones_like(T_TD)

fig, axes = plt.subplots(1, 2, figsize=(11, 5))
for ax in axes:
    ax.plot(T_TD, Cv_DP, color=COLOURS["grey"], lw=2.0, ls="--",
            label=r"Dulong-Petit ($3Nk_B$)")
    ax.plot(T_TD, Cv_einstein, color=COLOURS["red"], lw=2.2, label="Einstein")
    ax.plot(T_TD, Cv_debye, color=COLOURS["blue"], lw=2.2, label="Debye")
    ax.set_xlabel(r"$T / \Theta_D$")
    ax.set_ylabel(r"$C_V / 3Nk_B$")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(alpha=0.3, linestyle="--")

axes[0].set_xlim(0, 2.0); axes[0].set_ylim(0, 1.15)
axes[0].set_title("Linear scale: classical limit")
axes[0].annotate(r"Dulong-Petit limit", xy=(1.8, 1.0), xytext=(0.9, 1.08),
                 fontsize=11, arrowprops=dict(arrowstyle="->", color="k", lw=0.9))

# Inset low-T behaviour
axes[1].set_xscale("log"); axes[1].set_yscale("log")
mask = T_TD < 0.6
axes[1].plot(T_TD[mask], 234.0/(1/T_TD[mask])**3 * (T_TD[mask] >= 0)*0 + 234.0*T_TD[mask]**3 * (1/(1)),
             color="k", lw=1, ls=":", label=r"Debye low-$T$: $\propto T^3$")
# overlay actual debye
axes[1].plot(T_TD, Cv_debye, color=COLOURS["blue"], lw=2.2, label="Debye")
axes[1].plot(T_TD, Cv_einstein, color=COLOURS["red"], lw=2.2, label="Einstein")
axes[1].set_xlim(0.05, 2.0); axes[1].set_ylim(1e-4, 1.5)
axes[1].set_title(r"Log scale: low-$T$ scaling (Debye $\propto T^3$)")
axes[1].legend(frameon=False, loc="lower right")

fig.suptitle("Lattice heat capacity: three models compared", fontsize=15, y=1.02)
fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch03b/fig_debye_einstein_cv.png")
save(fig, out)
