#!/usr/bin/env python3
"""First four eigenfunctions and energy levels of the 1D infinite square well.

Two-panel: (left) energy levels diagram, (right) wavefunctions overlayed on levels.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS, PALETTE
import numpy as np
import matplotlib.pyplot as plt

setup()
L = 1.0
x = np.linspace(0, L, 500)
n_max = 4

fig, axes = plt.subplots(1, 2, figsize=(11, 5.5))

# LEFT: energy ladder
axL = axes[0]
for n in range(1, n_max+1):
    E = n**2
    axL.hlines(E, 0, 1, colors=PALETTE[n-1], linewidth=2.6)
    axL.text(1.04, E, fr"$E_{{{n}}} = {n**2}\,E_1$", va="center", fontsize=13,
             color=PALETTE[n-1])
axL.set_xlim(-0.05, 1.6)
axL.set_ylim(0, 18)
axL.set_xticks([])
axL.set_ylabel(r"$E / E_1$  where  $E_1 = \dfrac{\pi^2 \hbar^2}{2 m L^2}$")
axL.set_title("Energy levels: $E_n \\propto n^2$")
axL.grid(alpha=0.3, linestyle="--", axis="y")
axL.spines["bottom"].set_visible(False)

# RIGHT: wavefunctions on potential well
axR = axes[1]
axR.axvline(0, color="k", lw=2.5)
axR.axvline(L, color="k", lw=2.5)
for n in range(1, n_max+1):
    psi = np.sqrt(2/L) * np.sin(n*np.pi*x/L)
    E = n**2
    c = PALETTE[n-1]
    axR.fill_between(x, E, E + 0.5*psi, color=c, alpha=0.35)
    axR.plot(x, E + 0.5*psi, color=c, lw=2.0, label=fr"$n={n}$")
    axR.hlines(E, 0, L, colors=c, linewidth=0.8, linestyles=":")
axR.text(0.02, 16.5, r"$\psi_n(x) = \sqrt{2/L}\,\sin(n\pi x/L)$", fontsize=12)
axR.set_xlabel(r"$x / L$")
axR.set_ylabel(r"$E / E_1$  (wavefunctions $\psi_n$ overlayed)")
axR.set_title("Eigenfunctions on the infinite square well")
axR.set_xlim(-0.05, 1.05)
axR.set_ylim(0, 18)
axR.legend(frameon=False, loc="upper right")
axR.grid(alpha=0.3, linestyle="--", axis="y")

fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch04/fig_particle_in_box.png")
save(fig, out)
