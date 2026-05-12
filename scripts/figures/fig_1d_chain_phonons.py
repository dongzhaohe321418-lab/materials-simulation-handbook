#!/usr/bin/env python3
"""Monatomic and diatomic 1D chain phonon dispersions side-by-side."""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS
import numpy as np
import matplotlib.pyplot as plt

setup()
fig, axes = plt.subplots(1, 2, figsize=(11, 5))

# Monatomic
a, K, m = 1.0, 1.0, 1.0
k = np.linspace(-np.pi/a, np.pi/a, 600)
omega = 2*np.sqrt(K/m)*np.abs(np.sin(k*a/2))
axes[0].plot(k*a/np.pi, omega, color=COLOURS["blue"], lw=2.2,
             label=r"$\omega(k) = 2\sqrt{K/m}\,|\sin(ka/2)|$")
axes[0].axhline(2*np.sqrt(K/m), color=COLOURS["grey"], lw=0.8, ls="--")
axes[0].annotate(r"$\omega_{\max} = 2\sqrt{K/m}$", xy=(0.0, 2.0),
                 xytext=(-0.5, 2.15), fontsize=11,
                 arrowprops=dict(arrowstyle="->", color="k", lw=0.9))
# Linear (sound) regime near k=0
v_s = np.sqrt(K/m)*a
k_lin = np.linspace(-0.4, 0.4, 50)/a
axes[0].plot(k_lin*a/np.pi, v_s*np.abs(k_lin), color=COLOURS["red"], lw=1.4, ls=":",
             label=r"$\omega \approx v_s |k|$ (acoustic limit)")
axes[0].set_xlabel(r"$k a / \pi$")
axes[0].set_ylabel(r"$\omega / \sqrt{K/m}$")
axes[0].set_title("Monatomic chain: single acoustic branch")
axes[0].set_xlim(-1, 1); axes[0].set_ylim(0, 2.4)
axes[0].legend(frameon=False, loc="lower center", fontsize=10)
axes[0].grid(alpha=0.3, linestyle="--")

# Diatomic
m1, m2 = 1.0, 3.0
A = K*(1/m1 + 1/m2)
B = 4*K**2/(m1*m2)
omega_minus = np.sqrt(A - np.sqrt(A**2 - B*np.sin(k*a/2)**2))
omega_plus  = np.sqrt(A + np.sqrt(A**2 - B*np.sin(k*a/2)**2))
axes[1].plot(k*a/np.pi, omega_minus, color=COLOURS["blue"], lw=2.2, label="Acoustic")
axes[1].plot(k*a/np.pi, omega_plus,  color=COLOURS["red"],  lw=2.2, label="Optical")
gap_top = omega_plus[np.argmin(np.abs(k - np.pi/a))]
gap_bot = omega_minus[np.argmin(np.abs(k - np.pi/a))]
axes[1].annotate("", xy=(1.0, gap_top), xytext=(1.0, gap_bot),
                 arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
axes[1].text(0.5, 0.5*(gap_top+gap_bot), "phonon\nband gap", fontsize=11)
axes[1].set_xlabel(r"$k a / \pi$")
axes[1].set_ylabel(r"$\omega$ (arb. units)")
axes[1].set_title(fr"Diatomic chain ($m_2/m_1 = {m2/m1:.0f}$)")
axes[1].set_xlim(-1, 1)
axes[1].legend(frameon=False, loc="center right")
axes[1].grid(alpha=0.3, linestyle="--")

fig.tight_layout()
out = os.path.join(os.path.dirname(__file__), "../../docs/assets/figures/ch03b/fig_1d_chain_phonons.png")
save(fig, out)
