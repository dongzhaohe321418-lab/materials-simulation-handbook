#!/usr/bin/env python3
"""1D tight-binding band E(k) = -2t cos(ka) with effective-mass parabola near k=0."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
a, t = 1.0, 1.0
k = np.linspace(-np.pi/a, np.pi/a, 600)
E = -2*t*np.cos(k*a)
E_eff = -2*t + t*a**2*k**2

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(k*a/np.pi, E, color=COLOURS["blue"], lw=2.2,
        label=r"$E(k) = -2t\cos(ka)$")
mask = np.abs(k*a) < 1.0
ax.plot(k[mask]*a/np.pi, E_eff[mask], color=COLOURS["red"], lw=2.0, ls="--",
        label=r"Effective mass: $E \approx -2t + \dfrac{\hbar^2 k^2}{2 m^*}$")
ax.axhline(-2*t, color=COLOURS["grey"], lw=0.6, ls=":")
ax.axhline(2*t,  color=COLOURS["grey"], lw=0.6, ls=":")

ax.annotate(r"band minimum ($k=0$)", xy=(0, -2*t),
            xytext=(-0.5, -1.4), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
ax.annotate(r"band maximum ($k=\pm\pi/a$)", xy=(1, 2*t),
            xytext=(-0.2, 1.5), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
ax.text(0.05, 1.7, r"bandwidth $W = 4t$", fontsize=11, color=COLOURS["grey"])
ax.annotate("", xy=(0.9, -2*t), xytext=(0.9, 2*t),
            arrowprops=dict(arrowstyle="<->", color=COLOURS["grey"], lw=0.9))

ax.set_xlabel(r"$ka/\pi$")
ax.set_ylabel(r"$E / t$")
ax.set_title(r"1D tight-binding band and effective-mass approximation")
ax.legend(frameon=False, loc="lower center", fontsize=11)
ax.set_xlim(-1, 1)
ax.grid(alpha=0.3, linestyle="--")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch03b/fig_1d_tb_bands.png")
save(fig, out)
