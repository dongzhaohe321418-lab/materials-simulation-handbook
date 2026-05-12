#!/usr/bin/env python3
"""Empty lattice bands and gap-opening at zone boundary for nearly-free electrons."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
a = 1.0
k = np.linspace(-np.pi/a, np.pi/a, 700)
E1 = (k)**2
E2 = (k + 2*np.pi/a)**2
E3 = (k - 2*np.pi/a)**2

V = 1.5
def mix(Ea, Eb, V):
    avg = 0.5*(Ea+Eb); diff = 0.5*(Ea-Eb)
    s = np.sqrt(diff**2 + V**2)
    return avg - s, avg + s

band_low_pos = mix(E1, E3, V)[0]
band_low_neg = mix(E1, E2, V)[0]
band_high_pos = mix(E1, E3, V)[1]
band_high_neg = mix(E1, E2, V)[1]
band_low  = np.where(k>=0, band_low_pos,  band_low_neg)
band_high = np.where(k>=0, band_high_pos, band_high_neg)

fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
axes[0].plot(k*a/np.pi, E1, color=COLOURS["blue"], lw=2)
axes[0].plot(k*a/np.pi, E2, color=COLOURS["blue"], lw=2)
axes[0].plot(k*a/np.pi, E3, color=COLOURS["blue"], lw=2)
axes[0].axvline(1, color=COLOURS["grey"], lw=0.8, ls=":")
axes[0].axvline(-1, color=COLOURS["grey"], lw=0.8, ls=":")
axes[0].set_title(r"Empty lattice: $V = 0$")
axes[0].set_xlabel(r"$ka/\pi$"); axes[0].set_ylabel(r"$E$ (arb. units)")
axes[0].set_ylim(0, 35); axes[0].set_xlim(-1, 1)
axes[0].text(0, 32, r"bands cross at zone boundary", ha="center", fontsize=11)
axes[0].grid(alpha=0.3, linestyle="--")

axes[1].plot(k*a/np.pi, band_low,  color=COLOURS["red"], lw=2.2)
axes[1].plot(k*a/np.pi, band_high, color=COLOURS["red"], lw=2.2)
axes[1].axvline(1, color=COLOURS["grey"], lw=0.8, ls=":")
axes[1].axvline(-1, color=COLOURS["grey"], lw=0.8, ls=":")
gap_top = band_high[np.argmin(np.abs(k-np.pi/a))]
gap_bot = band_low[np.argmin(np.abs(k-np.pi/a))]
axes[1].annotate("", xy=(1, gap_top), xytext=(1, gap_bot),
                 arrowprops=dict(arrowstyle="<->", color="k", lw=1.4))
axes[1].text(0.55, 0.5*(gap_top+gap_bot), r"gap $= 2|V_G|$", fontsize=12)
axes[1].set_title(r"NFE: $V \ne 0$ opens gap at $k = \pi/a$")
axes[1].set_xlabel(r"$ka/\pi$")
axes[1].set_xlim(-1, 1)
axes[1].grid(alpha=0.3, linestyle="--")

fig.suptitle("Nearly-free-electron model: opening band gaps at zone boundaries",
             fontsize=15, y=1.02)
fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch03b/fig_nfe_band_gap.png")
save(fig, out)
