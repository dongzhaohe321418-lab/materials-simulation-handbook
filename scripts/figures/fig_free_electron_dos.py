#!/usr/bin/env python3
"""Free-electron DOS g(E) ~ sqrt(E) with E_F for Cu marked and filled states shaded."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
E = np.linspace(0, 12, 600)  # eV
g = np.sqrt(E)
E_F = 7.0  # Cu approximate

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(E, g, color=COLOURS["blue"], lw=2.2,
        label=r"$g(\varepsilon) = \dfrac{V}{2\pi^2}\left(\dfrac{2m}{\hbar^2}\right)^{3/2}\sqrt{\varepsilon}$")
ax.fill_between(E[E<=E_F], g[E<=E_F], color=COLOURS["blue"], alpha=0.3,
                label=r"Filled states at $T=0$")
ax.axvline(E_F, color=COLOURS["red"], ls="--", lw=2.0,
           label=fr"$\varepsilon_F \approx {E_F:.0f}$ eV (Cu)")
ax.annotate(r"Fermi energy $\varepsilon_F$", xy=(E_F, np.sqrt(E_F)),
            xytext=(8.5, 2.4), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="k", lw=1.0))
ax.set_xlabel(r"$\varepsilon$ (eV)")
ax.set_ylabel(r"$g(\varepsilon)$ (arb. units)")
ax.set_title(r"Free-electron density of states: $g(\varepsilon) \propto \sqrt{\varepsilon}$")
ax.set_xlim(0, 12)
ax.legend(frameon=False, loc="upper left", fontsize=11)
ax.grid(alpha=0.3, linestyle="--")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch03b/fig_free_electron_dos.png")
save(fig, out)
