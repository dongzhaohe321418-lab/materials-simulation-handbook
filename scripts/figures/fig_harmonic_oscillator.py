#!/usr/bin/env python3
"""First four eigenfunctions of the 1D quantum harmonic oscillator on parabolic potential."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS, PALETTE
import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.hermite import hermval
from math import factorial

setup()
x = np.linspace(-4.5, 4.5, 600)
V = 0.5 * x**2  # natural units: hbar=omega=m=1

def psi_n(n, x):
    coef = np.zeros(n+1); coef[n] = 1.0
    H = hermval(x, coef)
    norm = 1.0/np.sqrt(2**n * factorial(n) * np.sqrt(np.pi))
    return norm * H * np.exp(-x**2/2)

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(x, V, color="k", lw=2.0, label=r"$V(x) = \frac{1}{2} m\omega^2 x^2$")

for n in range(4):
    E = n + 0.5
    psi = psi_n(n, x)
    c = PALETTE[n]
    ax.hlines(E, x[0], x[-1], colors=c, linewidth=0.8, linestyles=":")
    ax.plot(x, E + 0.55*psi, color=c, lw=2.0, label=fr"$n={n}$, $E_{{{n}}}={E:.1f}\,\hbar\omega$")
    ax.fill_between(x, E, E + 0.55*psi, color=c, alpha=0.3)

# Annotate zero-point energy
ax.annotate(r"zero-point energy $\frac{1}{2}\hbar\omega$",
            xy=(-1.0, 0.5), xytext=(-4.2, 1.6), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="k", lw=1.0))
# Annotate even spacing
ax.annotate(r"$\Delta E = \hbar\omega$", xy=(3.6, 2.5), xytext=(2.8, 4.4),
            fontsize=12, arrowprops=dict(arrowstyle="->", color="k", lw=1.0))

ax.set_xlabel(r"$x$ (natural units)")
ax.set_ylabel(r"$E / \hbar\omega$  (with $\psi_n$ overlayed)")
ax.set_title("Quantum harmonic oscillator: first four eigenstates")
ax.set_xlim(-4.5, 4.5)
ax.set_ylim(0, 6)
ax.legend(frameon=False, loc="upper right", fontsize=11)
ax.grid(alpha=0.3, linestyle="--", axis="y")

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch04/fig_harmonic_oscillator.png")
save(fig, out)
