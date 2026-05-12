#!/usr/bin/env python3
"""Graphene nearest-neighbour tight-binding bands along K-Gamma-M-K with Dirac point."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
a, t = 1.0, 1.0

def E_pm(k):
    kx, ky = k[..., 0], k[..., 1]
    f = 1 + 4*np.cos(np.sqrt(3)*ky*a/2)*np.cos(3*kx*a/2) + 4*np.cos(np.sqrt(3)*ky*a/2)**2
    f = np.maximum(f, 0.0)
    return -t*np.sqrt(f), t*np.sqrt(f)

G = np.array([0.0, 0.0])
M = np.array([0.0, 2*np.pi/(np.sqrt(3)*a)])
K = np.array([2*np.pi/(3*a), 2*np.pi/(3*np.sqrt(3)*a)])

def path(p1, p2, n=160):
    return np.array([p1 + (p2-p1)*t for t in np.linspace(0, 1, n)])
def dists(p1, p2, n):
    return np.linspace(0, np.linalg.norm(p2-p1), n)

seg1 = path(K, G); seg2 = path(G, M); seg3 = path(M, K)
ks = np.concatenate([seg1, seg2[1:], seg3[1:]])
d1 = dists(K, G, 160)
d2 = dists(G, M, 160) + d1[-1]
d3 = dists(M, K, 160) + d2[-1]
xax = np.concatenate([d1, d2[1:], d3[1:]])

Em, Ep = E_pm(ks)

fig, ax = plt.subplots(figsize=(8.5, 5.5))
ax.plot(xax, Ep, color=COLOURS["red"], lw=2.2, label=r"$\pi^*$ (anti-bonding)")
ax.plot(xax, Em, color=COLOURS["blue"], lw=2.2, label=r"$\pi$ (bonding)")
ax.axhline(0, color="k", lw=0.8, ls="--", label=r"Fermi level $\varepsilon_F$")

ticks = [0, d1[-1], d2[-1], d3[-1]]
ax.set_xticks(ticks)
ax.set_xticklabels(["K", r"$\Gamma$", "M", "K"])
for x in ticks:
    ax.axvline(x, color=COLOURS["grey"], lw=0.6, ls=":")
ax.set_ylabel(r"$E / t$")
ax.set_title("Graphene tight-binding bands (touching at Dirac point K)")
ax.legend(frameon=False, loc="upper right")
ax.grid(alpha=0.3, linestyle="--", axis="y")

ax.annotate("Dirac point\n(bands touch at $\\varepsilon_F$)",
            xy=(0, 0), xytext=(0.5, 1.8), fontsize=12,
            arrowprops=dict(arrowstyle="->", color="k", lw=1.0))

out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch03b/fig_graphene_bands.png")
save(fig, out)
