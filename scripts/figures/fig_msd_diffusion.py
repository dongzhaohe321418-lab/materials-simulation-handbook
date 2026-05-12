#!/usr/bin/env python3
"""MSD vs time for solid, liquid, gas."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
t = np.linspace(0, 10, 500)
msd_solid = 0.12*(1 - np.exp(-3*t))
D_liq = 0.25
msd_liq_ballistic = (1.5*t)**2
msd_liq_diff = 6*D_liq*t
with np.errstate(divide="ignore", invalid="ignore"):
    msd_liq = 1.0/(1.0/msd_liq_ballistic[1:] + 1.0/msd_liq_diff[1:])
msd_liq = np.concatenate(([0], msd_liq))
D_gas = 1.0
msd_gas = 6*D_gas*t

fig, ax = plt.subplots(figsize=(8.5, 5.5))
ax.plot(t, msd_solid, color=COLOURS["blue"], lw=2.2,
        label=r"Solid: $\langle r^2 \rangle \to$ constant (bounded vibration)")
ax.plot(t, msd_liq,   color=COLOURS["red"], lw=2.2,
        label=fr"Liquid: linear, slope $= 6D$, $D \approx {D_liq}$")
ax.plot(t, msd_gas,   color=COLOURS["green"], lw=2.2,
        label=fr"Gas: large $D \approx {D_gas}$")

ax.annotate(r"slope $= 6D$  ($d=3$)", xy=(7.5, 6*D_liq*7.5),
            xytext=(5.0, 9.5), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="k", lw=1.0))
ax.annotate("ballistic regime\n$\\langle r^2\\rangle \\propto t^2$",
            xy=(0.3, msd_liq[15]), xytext=(0.6, 4.0), fontsize=11,
            arrowprops=dict(arrowstyle="->", color="k", lw=0.9))

ax.set_xlabel(r"$t$ (ps)")
ax.set_ylabel(r"$\langle |\mathbf{r}(t) - \mathbf{r}(0)|^2 \rangle$ (Å$^2$)")
ax.set_title("Mean-squared displacement: phase signatures (schematic)")
ax.set_xlim(0, 10); ax.set_ylim(0, 18)
ax.legend(frameon=False, loc="upper left")
ax.grid(alpha=0.3, linestyle="--")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch07/fig_msd_diffusion.png")
save(fig, out)
