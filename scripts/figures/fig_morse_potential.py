#!/usr/bin/env python3
"""Morse potential with harmonic approximation overlay near the minimum."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
De, a, re = 1.0, 1.5, 1.0
r = np.linspace(0.5, 4.0, 600)
V_morse = De * (1 - np.exp(-a*(r-re)))**2 - De
k = 2*De*a**2  # harmonic curvature at minimum
V_harm = 0.5*k*(r-re)**2 - De

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(r, V_morse, color=COLOURS["blue"], lw=2.2,
        label=r"Morse: $V(r) = D_e\,[1 - e^{-a(r-r_e)}]^2 - D_e$")
ax.plot(r, V_harm, color=COLOURS["red"], lw=2.0, ls="--",
        label=r"Harmonic: $V(r) \approx \frac{1}{2} k(r-r_e)^2 - D_e$")
ax.axhline(0, color="k", lw=0.5)
ax.axhline(-De, color=COLOURS["grey"], lw=0.6, ls=":")

ax.annotate(r"dissociation energy $D_e$", xy=(3.5, 0), xytext=(2.6, -0.4),
            arrowprops=dict(arrowstyle="<->", color="k", lw=1.0), fontsize=12)
ax.text(2.0, -0.95, r"$-D_e$", color=COLOURS["grey"], fontsize=12)
ax.annotate(r"harmonic fails for large $|r-r_e|$", xy=(2.2, 0.4),
            xytext=(2.4, 1.0),
            arrowprops=dict(arrowstyle="->", color="k", lw=1.0), fontsize=11)

ax.set_xlabel(r"$r$ (arbitrary units)")
ax.set_ylabel(r"$V(r)$ (arbitrary units)")
ax.set_title("Morse potential vs harmonic approximation")
ax.set_ylim(-1.2, 1.5)
ax.set_xlim(0.5, 4.0)
ax.grid(alpha=0.3, linestyle="--")
ax.legend(frameon=False, loc="upper right")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch04/fig_morse_potential.png")
save(fig, out)
