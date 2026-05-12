#!/usr/bin/env python3
"""Synthetic g(r) for crystal, liquid, and gas to illustrate phase signatures."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
r = np.linspace(0.0, 6.0, 800)

def gauss(x, mu, sig):
    return np.exp(-0.5*((x-mu)/sig)**2)

crystal = np.zeros_like(r)
shells = [1.0, 1.41, 1.73, 2.0, 2.24, 2.45, 2.83, 3.16, 3.46, 3.74, 4.0, 4.24]
heights = [3.5, 2.6, 2.2, 2.0, 1.9, 1.7, 1.6, 1.55, 1.5, 1.45, 1.4, 1.35]
for s, h in zip(shells, heights):
    crystal += (h-1) * gauss(r, s, 0.03)
crystal += 1.0
crystal[r<0.85] = 0.0

liquid = np.ones_like(r)
liquid += 2.5*gauss(r, 1.05, 0.10)
liquid += 0.6*gauss(r, 2.0,  0.18)
liquid += 0.25*gauss(r, 3.0, 0.25)
liquid[r<0.85] = 0.0

gas = 1.0 - np.exp(-(r/1.5)**6)
gas[r<0.85] = 0.0

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(r, crystal, color=COLOURS["blue"], lw=2.2, label="Crystal")
ax.plot(r, liquid,  color=COLOURS["red"], lw=2.2, label="Liquid")
ax.plot(r, gas,     color=COLOURS["green"], lw=2.2, label="Gas")
ax.axhline(1, color=COLOURS["grey"], lw=0.7, ls=":")
ax.annotate(r"sharp peaks $\Rightarrow$ long-range order",
            xy=(1.0, 3.5), xytext=(1.6, 4.0), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
ax.annotate(r"first coordination shell",
            xy=(1.05, 3.5), xytext=(0.05, 3.0), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
ax.annotate(r"$g(r) \to 1$ (no correlations)", xy=(5, 1.05),
            xytext=(3.8, 1.7), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.9))

ax.set_xlabel(r"$r$ (reduced units)")
ax.set_ylabel(r"$g(r)$")
ax.set_title("Radial distribution function: phase signatures (schematic)")
ax.set_xlim(0, 6)
ax.set_ylim(0, 4.7)
ax.legend(frameon=False, loc="upper right")
ax.grid(alpha=0.3, linestyle="--")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch07/fig_rdf_phases.png")
save(fig, out)
