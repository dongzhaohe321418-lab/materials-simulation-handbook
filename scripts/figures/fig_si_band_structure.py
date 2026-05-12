#!/usr/bin/env python3
"""Schematic Si band structure: VBM at Gamma, CBM near X => indirect gap."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()

nk = 250
seg1 = np.linspace(0, 1, nk)
seg2 = np.linspace(1, 2.2, nk)
seg3 = np.linspace(2.2, 3.0, nk)
seg4 = np.linspace(3.0, 4.2, nk)
k_all = np.concatenate([seg1, seg2, seg3, seg4])

def parabola(k, k0, c, top):
    return top - c*(k-k0)**2

VB_top = parabola(k_all, 1.0, 1.5, 0.0)
VB_mid = parabola(k_all, 1.0, 1.6, -0.8)
VB_low = parabola(k_all, 1.0, 1.2, -3.0)
CB1 = parabola(k_all, 2.0, -1.0, 1.15)  # CBM near X
CB2 = parabola(k_all, 1.0, -2.5, 3.4)   # direct gap at Gamma higher

fig, ax = plt.subplots(figsize=(9, 6))
for band, col, lab in [(VB_top, COLOURS["blue"], "valence"),
                       (VB_mid, COLOURS["blue"], None),
                       (VB_low, COLOURS["blue"], None),
                       (CB1, COLOURS["red"], "conduction"),
                       (CB2, COLOURS["red"], None)]:
    ax.plot(k_all, band, color=col, lw=2.2, label=lab)

# Indirect gap arrow
ax.annotate("", xy=(2.0, 1.15), xytext=(1.0, 0.0),
            arrowprops=dict(arrowstyle="->", color=COLOURS["green"], lw=2.6))
ax.text(1.45, 0.45, r"indirect gap $\approx 1.17$ eV",
        color=COLOURS["green"], fontsize=13)

# Direct gap arrow
ax.annotate("", xy=(1.0, 3.4), xytext=(1.0, 0.0),
            arrowprops=dict(arrowstyle="->", color=COLOURS["orange"], lw=1.8, ls="--"))
ax.text(0.7, 1.8, r"direct gap $\approx 3.4$ eV",
        color=COLOURS["orange"], fontsize=12, rotation=90)

ax.text(1.0, -3.5, "VBM", fontsize=11, ha="center")
ax.text(2.0, 1.4,  "CBM", fontsize=11, ha="center", color=COLOURS["red"])

ticks = [0, 1, 2.2, 3.0, 4.2]
labels = ["L", r"$\Gamma$", "X", "U,K", r"$\Gamma$"]
ax.set_xticks(ticks); ax.set_xticklabels(labels, fontsize=13)
for t in ticks:
    ax.axvline(t, color=COLOURS["grey"], lw=0.6, ls=":")
ax.axhline(0, color="k", lw=0.8, ls="--")

ax.set_ylabel(r"$E - E_{\rm VBM}$ (eV)")
ax.set_title("Silicon band structure (schematic): indirect gap")
ax.set_xlim(0, 4.2); ax.set_ylim(-4, 5)
ax.legend(frameon=False, loc="upper right")
ax.grid(alpha=0.3, linestyle="--", axis="y")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch06/fig_si_band_structure.png")
save(fig, out)
