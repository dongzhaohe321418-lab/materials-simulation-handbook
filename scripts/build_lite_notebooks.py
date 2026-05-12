"""Generate JupyterLite-compatible notebooks for the handbook.

Each notebook mirrors the meatiest, Pyodide-compatible code from a chapter and
includes a markdown header linking back to the corresponding page on the main
site. Run from the repository root::

    python scripts/build_lite_notebooks.py

This script is deterministic — re-running it overwrites the generated
notebooks/*.ipynb files in place. Hand-edits to those files will be lost.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO_ROOT / "notebooks"
SITE_BASE = "https://dongzhaohe321418-lab.github.io/materials-simulation-handbook"


def md(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.splitlines(keepends=True) or [""],
    }


def code(source: str) -> dict:
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": source.splitlines(keepends=True) or [""],
    }


def notebook(cells: Iterable[dict]) -> dict:
    return {
        "cells": list(cells),
        "metadata": {
            "kernelspec": {
                "name": "python",
                "display_name": "Python (Pyodide)",
                "language": "python",
            },
            "language_info": {"name": "python"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def write_notebook(name: str, cells: list[dict]) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / name
    path.write_text(json.dumps(notebook(cells), indent=1) + "\n")
    return path


def header(chapter_title: str, slug: str, summary: str) -> dict:
    url = f"{SITE_BASE}/{slug}/"
    text = (
        f"# {chapter_title} — live notebook\n\n"
        f"This notebook runs entirely in your browser via the Pyodide kernel.\n"
        f"All state is saved to your browser's local storage; nothing is sent to a server.\n\n"
        f"Back to the chapter: <{url}>\n\n"
        f"{summary}\n\n"
        f"Only Pyodide-compatible packages are used (numpy, scipy, matplotlib, ipywidgets).\n"
    )
    return md(text)


# ---------------------------------------------------------------------------
# Chapter 0 — Mathematics from Scratch
# ---------------------------------------------------------------------------
CH00 = [
    header(
        "Chapter 0 · Mathematics from Scratch",
        "ch00-math",
        "We touch the four threads of the chapter: linear algebra, calculus, "
        "Fourier analysis and probability. Each section is a few self-contained "
        "cells you can edit and re-run.",
    ),
    md("## 0.1 Linear algebra — eigenvalues of a symmetric matrix\n\n"
       "The single most important operation in computational materials science. "
       "Every Hamiltonian diagonalisation, every principal-component analysis, "
       "every spectral graph technique reduces to this."),
    code("import numpy as np\n\n"
         "rng = np.random.default_rng(0)\n"
         "A = rng.standard_normal((5, 5))\n"
         "A = (A + A.T) / 2  # symmetrise\n"
         "eigvals, eigvecs = np.linalg.eigh(A)\n"
         "print('eigenvalues:', np.round(eigvals, 4))\n"
         "print('orthonormality check:', np.allclose(eigvecs.T @ eigvecs, np.eye(5)))\n"),
    md("## 0.2 Calculus — gradient descent on a paraboloid\n\n"
       "The workhorse of geometry optimisation and machine-learning training. "
       "Tweak the learning rate and watch convergence change."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "def f(x):\n"
         "    return 0.5 * (x[0] ** 2 + 10.0 * x[1] ** 2)\n\n"
         "def grad(x):\n"
         "    return np.array([x[0], 10.0 * x[1]])\n\n"
         "x = np.array([3.0, 1.5])\n"
         "lr = 0.1\n"
         "trail = [x.copy()]\n"
         "for _ in range(40):\n"
         "    x = x - lr * grad(x)\n"
         "    trail.append(x.copy())\n"
         "trail = np.array(trail)\n\n"
         "xs, ys = np.meshgrid(np.linspace(-3.5, 3.5, 80), np.linspace(-2, 2, 80))\n"
         "zs = 0.5 * (xs ** 2 + 10.0 * ys ** 2)\n"
         "fig, ax = plt.subplots(figsize=(6, 4))\n"
         "ax.contour(xs, ys, zs, levels=20, colors='lightgrey')\n"
         "ax.plot(trail[:, 0], trail[:, 1], 'o-', color='C3', ms=3)\n"
         "ax.set_xlabel('x')\n"
         "ax.set_ylabel('y')\n"
         "ax.set_title('Gradient descent on f(x,y) = (x^2 + 10y^2)/2')\n"
         "plt.show()\n"),
    md("## 0.3 Fourier — a square wave from sines\n\n"
       "Reciprocal-space thinking is built on this picture: any well-behaved "
       "periodic function is a sum of sinusoids."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "x = np.linspace(0, 2 * np.pi, 1000)\n"
         "target = np.sign(np.sin(x))\n\n"
         "fig, ax = plt.subplots(figsize=(7, 4))\n"
         "ax.plot(x, target, 'k', lw=1.5, label='target')\n"
         "approx = np.zeros_like(x)\n"
         "for k in range(1, 12, 2):\n"
         "    approx += (4 / (np.pi * k)) * np.sin(k * x)\n"
         "    if k in (1, 3, 7, 11):\n"
         "        ax.plot(x, approx, label=f'up to k={k}', lw=1)\n"
         "ax.legend()\n"
         "ax.set_xlabel('x')\n"
         "ax.set_ylabel('partial Fourier sum')\n"
         "plt.show()\n"),
    md("## 0.4 Probability — the central limit theorem in action\n\n"
       "Sums of independent random variables converge to a Gaussian. This is "
       "why noise in everything from thermostats to GP priors is normally "
       "distributed."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "rng = np.random.default_rng(42)\n"
         "fig, axes = plt.subplots(1, 3, figsize=(10, 3), sharey=True)\n"
         "for ax, n in zip(axes, [1, 5, 50]):\n"
         "    samples = rng.uniform(-1, 1, size=(20000, n)).mean(axis=1)\n"
         "    ax.hist(samples, bins=60, density=True, alpha=0.7)\n"
         "    ax.set_title(f'mean of {n} uniforms')\n"
         "fig.tight_layout()\n"
         "plt.show()\n"),
]


# ---------------------------------------------------------------------------
# Chapter 1 — Python and Scientific Computing
# ---------------------------------------------------------------------------
CH01 = [
    header(
        "Chapter 1 · Python and Scientific Computing",
        "ch01-python",
        "Hands-on NumPy and Matplotlib drills that mirror the chapter. "
        "Every cell is editable; experiment freely.",
    ),
    md("## 1.1 NumPy basics — vectorisation beats loops\n"),
    code("import numpy as np\n"
         "import time\n\n"
         "N = 1_000_000\n"
         "x = np.arange(N, dtype=float)\n\n"
         "t0 = time.perf_counter()\n"
         "total = 0.0\n"
         "for v in x[:50_000]:\n"
         "    total += v * v\n"
         "loop_t = time.perf_counter() - t0\n\n"
         "t0 = time.perf_counter()\n"
         "_ = (x * x).sum()\n"
         "vec_t = time.perf_counter() - t0\n\n"
         "print(f'loop (50k items):  {loop_t * 1000:.1f} ms')\n"
         "print(f'vectorised ({N}): {vec_t * 1000:.2f} ms')\n"),
    md("## 1.2 Broadcasting — pairwise distance matrix without a loop\n"),
    code("import numpy as np\n\n"
         "rng = np.random.default_rng(0)\n"
         "pts = rng.uniform(0, 1, size=(6, 3))\n"
         "diff = pts[:, None, :] - pts[None, :, :]\n"
         "dist = np.linalg.norm(diff, axis=-1)\n"
         "print(np.round(dist, 3))\n"),
    md("## 1.3 Matplotlib — a publication-grade figure\n"),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "x = np.linspace(0, 2 * np.pi, 400)\n"
         "fig, ax = plt.subplots(figsize=(6, 4))\n"
         "for k in (1, 2, 3):\n"
         "    ax.plot(x, np.sin(k * x) / k, label=f'sin({k}x)/{k}')\n"
         "ax.set_xlabel('x')\n"
         "ax.set_ylabel('amplitude')\n"
         "ax.set_title('A small family of sinusoids')\n"
         "ax.legend(frameon=False)\n"
         "ax.grid(alpha=0.3)\n"
         "fig.tight_layout()\n"
         "plt.show()\n"),
    md("## 1.4 An ipywidgets slider — interactive learning\n\n"
       "Move the slider and re-run the cell to change the frequency. "
       "If the slider does not appear, run the cell once more after the "
       "kernel finishes loading."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n"
         "from ipywidgets import interact, FloatSlider\n\n"
         "x = np.linspace(0, 2 * np.pi, 400)\n\n"
         "@interact(k=FloatSlider(min=0.5, max=5.0, step=0.1, value=1.0))\n"
         "def plot_sine(k):\n"
         "    fig, ax = plt.subplots(figsize=(6, 3))\n"
         "    ax.plot(x, np.sin(k * x))\n"
         "    ax.set_ylim(-1.1, 1.1)\n"
         "    ax.set_title(f'sin({k:.1f} x)')\n"
         "    plt.show()\n"),
]


# ---------------------------------------------------------------------------
# Chapter 3.5 — Solid State Physics Prerequisites
# ---------------------------------------------------------------------------
CH035 = [
    header(
        "Chapter 3.5 · Solid State Physics Prerequisites",
        "ch03b-solid-state",
        "Two cornerstones of solid-state physics implemented from scratch: "
        "the tight-binding band of a 1D chain, and the dispersion of a "
        "monatomic phonon chain.",
    ),
    md("## 3.5.3 Tight-binding 1D chain\n\n"
       "With on-site energy $\\varepsilon$ and nearest-neighbour hopping $-t$, "
       "the band structure is $E(k) = \\varepsilon - 2t \\cos(ka)$. "
       "Diagonalise a finite chain to recover the same shape."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "N = 200          # number of sites (periodic chain)\n"
         "eps = 0.0\n"
         "t = 1.0\n"
         "a = 1.0\n\n"
         "H = eps * np.eye(N) - t * (np.eye(N, k=1) + np.eye(N, k=-1))\n"
         "H[0, -1] = H[-1, 0] = -t  # periodic boundary\n"
         "energies = np.sort(np.linalg.eigvalsh(H))\n\n"
         "k = np.linspace(-np.pi / a, np.pi / a, 400)\n"
         "E_analytic = eps - 2 * t * np.cos(k * a)\n\n"
         "fig, ax = plt.subplots(figsize=(6, 4))\n"
         "ax.plot(k, E_analytic, 'k-', lw=1.5, label='analytic')\n"
         "ax.plot(np.linspace(-np.pi, np.pi, N), energies, 'o', ms=3, alpha=0.6, label='finite chain')\n"
         "ax.set_xlabel('k a')\n"
         "ax.set_ylabel('E / t')\n"
         "ax.legend()\n"
         "ax.set_title(f'Tight-binding band, N = {N}')\n"
         "plt.show()\n"),
    md("## 3.5.5 Phonons of a 1D monatomic chain\n\n"
       "Spring constant $K$, mass $m$. Dispersion: "
       "$\\omega(k) = 2\\sqrt{K/m}\\,|\\sin(ka/2)|$."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "K = 1.0\n"
         "m = 1.0\n"
         "a = 1.0\n"
         "k = np.linspace(-np.pi / a, np.pi / a, 400)\n"
         "omega = 2 * np.sqrt(K / m) * np.abs(np.sin(k * a / 2))\n\n"
         "fig, ax = plt.subplots(figsize=(6, 4))\n"
         "ax.plot(k * a, omega)\n"
         "ax.set_xlabel('k a')\n"
         "ax.set_ylabel('omega (sqrt(K/m))')\n"
         "ax.set_title('Monatomic phonon chain')\n"
         "ax.grid(alpha=0.3)\n"
         "plt.show()\n"),
    md("Compare with the long-wavelength sound speed $c_s = a\\sqrt{K/m}$ "
       "(slope of the curve at $k = 0$)."),
]


# ---------------------------------------------------------------------------
# Chapter 4 — Quantum Mechanics: particle in a box & harmonic oscillator
# ---------------------------------------------------------------------------
CH04 = [
    header(
        "Chapter 4 · Quantum Mechanics for Materials",
        "ch04-quantum",
        "Finite-difference solution of two textbook problems: the particle in "
        "a box (analytical comparison) and the quantum harmonic oscillator "
        "(spotting the equal level spacing).",
    ),
    md("## 4.3 Particle in a 1D box — finite-difference diagonalisation\n\n"
       "Build the Hamiltonian as a tridiagonal matrix and diagonalise. "
       "Compare against the closed form $E_n = n^2 \\pi^2 \\hbar^2 / (2 m L^2)$."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "HBAR = 1.054_571_817e-34\n"
         "M_E = 9.109_383_7e-31\n"
         "EV = 1.602_176_634e-19\n\n"
         "def build_hamiltonian(n_grid, box_length, potential=None, mass=M_E):\n"
         "    h = box_length / (n_grid + 1)\n"
         "    x = np.linspace(h, box_length - h, n_grid)\n"
         "    pref = HBAR ** 2 / (2.0 * mass * h ** 2)\n"
         "    main = 2.0 * pref * np.ones(n_grid)\n"
         "    off = -pref * np.ones(n_grid - 1)\n"
         "    H = np.diag(main) + np.diag(off, k=1) + np.diag(off, k=-1)\n"
         "    if potential is not None:\n"
         "        H = H + np.diag(potential)\n"
         "    return x, H\n\n"
         "L = 1.0e-9\n"
         "N = 400\n"
         "x, H = build_hamiltonian(N, L)\n"
         "eigvals, eigvecs = np.linalg.eigh(H)\n"
         "h = L / (N + 1)\n"
         "eigvecs = eigvecs / np.sqrt(h)\n\n"
         "print(f'{\"n\":>3} {\"E_num (eV)\":>12} {\"E_ana (eV)\":>12} {\"rel err\":>10}')\n"
         "for n in range(1, 5):\n"
         "    e_ana = (n ** 2 * np.pi ** 2 * HBAR ** 2) / (2.0 * M_E * L ** 2)\n"
         "    e_num = eigvals[n - 1]\n"
         "    rel = abs(e_num - e_ana) / e_ana\n"
         "    print(f'{n:>3d} {e_num / EV:>12.6f} {e_ana / EV:>12.6f} {rel:>10.2e}')\n"),
    code("fig, axes = plt.subplots(2, 2, figsize=(8, 5), sharex=True)\n"
         "for n, ax in zip(range(1, 5), axes.ravel()):\n"
         "    psi_ana = np.sqrt(2.0 / L) * np.sin(n * np.pi * x / L)\n"
         "    psi_num = eigvecs[:, n - 1]\n"
         "    if np.dot(psi_num, psi_ana) < 0:\n"
         "        psi_num = -psi_num\n"
         "    ax.plot(x * 1e9, psi_ana, 'k-', lw=2, label='analytic')\n"
         "    ax.plot(x * 1e9, psi_num, 'r--', lw=1.2, label='finite-diff')\n"
         "    ax.set_title(f'n = {n}, E = {eigvals[n - 1] / EV:.3f} eV')\n"
         "    ax.set_xlabel('x (nm)')\n"
         "    ax.legend(fontsize=8)\n"
         "fig.tight_layout()\n"
         "plt.show()\n"),
    md("## 4.4 Harmonic oscillator — same code, different potential\n\n"
       "Replace V(x) with $\\tfrac{1}{2} m \\omega^2 x^2$; expect equally spaced "
       "levels $E_n = \\hbar \\omega (n + 1/2)$."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "# Atomic-ish units: hbar = m = omega = 1\n"
         "N = 600\n"
         "L = 20.0       # box from -L/2 to L/2\n"
         "h = L / (N + 1)\n"
         "x = np.linspace(-L / 2 + h, L / 2 - h, N)\n"
         "V = 0.5 * x ** 2\n"
         "pref = 1.0 / (2.0 * h ** 2)\n"
         "main = 2.0 * pref + V\n"
         "off = -pref * np.ones(N - 1)\n"
         "H = np.diag(main) + np.diag(off, k=1) + np.diag(off, k=-1)\n"
         "vals, vecs = np.linalg.eigh(H)\n"
         "print('first six levels:', np.round(vals[:6], 4))\n"
         "print('expected n+1/2  :', np.arange(6) + 0.5)\n\n"
         "fig, ax = plt.subplots(figsize=(6, 4))\n"
         "for n in range(4):\n"
         "    psi = vecs[:, n] / np.sqrt(h)\n"
         "    ax.plot(x, psi + vals[n], label=f'n={n}')\n"
         "ax.plot(x, 0.5 * x ** 2, 'k', lw=0.7)\n"
         "ax.set_ylim(-0.2, 4.5)\n"
         "ax.set_xlim(-5, 5)\n"
         "ax.set_xlabel('x')\n"
         "ax.set_ylabel('E (and shifted psi)')\n"
         "ax.legend(fontsize=8)\n"
         "ax.set_title('Quantum harmonic oscillator')\n"
         "plt.show()\n"),
]


# ---------------------------------------------------------------------------
# Chapter 5 — DFT (the SCF iteration as a toy 1D Kohn-Sham loop)
# ---------------------------------------------------------------------------
CH05 = [
    header(
        "Chapter 5 · Density Functional Theory",
        "ch05-dft",
        "A toy self-consistent-field loop in 1D. We pretend the "
        "exchange-correlation functional is a simple LDA-like local term and "
        "watch the density converge under linear mixing.",
    ),
    md("## A miniature 1D Kohn–Sham SCF\n\n"
       "Two electrons in a soft Coulombic external well. The Hartree term is "
       "computed from the running density by direct integration; XC is "
       "approximated by $v_{xc}(r) = -c \\rho^{1/3}$ — a 1D caricature of LDA. "
       "Watch the convergence of the total density."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "N = 300\n"
         "L = 14.0\n"
         "x = np.linspace(-L / 2, L / 2, N)\n"
         "dx = x[1] - x[0]\n\n"
         "def kinetic_matrix(n, dx):\n"
         "    pref = 0.5 / dx ** 2\n"
         "    return np.diag(2 * pref * np.ones(n)) - np.diag(pref * np.ones(n - 1), 1) \\\n"
         "        - np.diag(pref * np.ones(n - 1), -1)\n\n"
         "def soft_coulomb(x, alpha=1.0, x0=0.0):\n"
         "    return -1.0 / np.sqrt((x - x0) ** 2 + alpha ** 2)\n\n"
         "v_ext = soft_coulomb(x)\n"
         "T = kinetic_matrix(N, dx)\n\n"
         "def trapz(y, x):\n"
         "    return float(0.5 * np.sum((y[1:] + y[:-1]) * np.diff(x)))\n\n"
         "def hartree(rho, x):\n"
         "    # 1D soft Coulomb interaction\n"
         "    return np.array([trapz(rho / np.sqrt((x - xi) ** 2 + 1.0), x) for xi in x])\n\n"
         "def vxc(rho, c=0.5):\n"
         "    return -c * np.cbrt(np.clip(rho, 1e-12, None))\n\n"
         "# Initial guess: ground state of the bare external potential\n"
         "H0 = T + np.diag(v_ext)\n"
         "vals, vecs = np.linalg.eigh(H0)\n"
         "psi = vecs[:, 0] / np.sqrt(dx)\n"
         "rho = 2.0 * psi ** 2  # two electrons, spin-paired\n\n"
         "history = [rho.copy()]\n"
         "mix = 0.3\n"
         "for it in range(40):\n"
         "    v_eff = v_ext + hartree(rho, x) + vxc(rho)\n"
         "    H = T + np.diag(v_eff)\n"
         "    vals, vecs = np.linalg.eigh(H)\n"
         "    psi_new = vecs[:, 0] / np.sqrt(dx)\n"
         "    rho_new = 2.0 * psi_new ** 2\n"
         "    delta = np.max(np.abs(rho_new - rho))\n"
         "    rho = mix * rho_new + (1 - mix) * rho\n"
         "    history.append(rho.copy())\n"
         "    if delta < 1e-6:\n"
         "        break\n"
         "print(f'converged in {it + 1} iterations, eigenvalue = {vals[0]:.4f} Ha')\n"),
    code("fig, ax = plt.subplots(figsize=(7, 4))\n"
         "for k, r in enumerate(history[:: max(1, len(history) // 6)]):\n"
         "    ax.plot(x, r, label=f'iter {k * max(1, len(history) // 6)}', alpha=0.7)\n"
         "ax.set_xlabel('x (bohr)')\n"
         "ax.set_ylabel('density rho(x)')\n"
         "ax.set_title('SCF convergence under linear mixing')\n"
         "ax.legend(fontsize=8)\n"
         "plt.show()\n"),
    md("Experiment: shrink `mix` to 0.1 — convergence slows. "
       "Push it past 1.0 and you see oscillations or divergence, exactly the "
       "pathology DIIS and Pulay mixing are designed to avoid."),
]


# ---------------------------------------------------------------------------
# Chapter 7 — Molecular Dynamics: Euler vs velocity-Verlet
# ---------------------------------------------------------------------------
CH07 = [
    header(
        "Chapter 7 · Molecular Dynamics",
        "ch07-md",
        "The defining demonstration of why symplectic integrators win: a "
        "harmonic oscillator integrated with Euler vs velocity-Verlet, looking "
        "at energy drift over many periods.",
    ),
    md("## 7.1 Euler vs velocity-Verlet on a 1D oscillator\n\n"
       "Equation of motion: $\\ddot{x} = -x$. Total energy "
       "$E = \\tfrac{1}{2}(\\dot{x}^2 + x^2)$ should be conserved. "
       "Euler grows it exponentially; Verlet keeps it bounded forever."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "def euler(x0, v0, dt, n_steps):\n"
         "    x = np.empty(n_steps + 1)\n"
         "    v = np.empty(n_steps + 1)\n"
         "    x[0], v[0] = x0, v0\n"
         "    for n in range(n_steps):\n"
         "        a = -x[n]\n"
         "        x[n + 1] = x[n] + dt * v[n]\n"
         "        v[n + 1] = v[n] + dt * a\n"
         "    return x, v\n\n"
         "def velocity_verlet(x0, v0, dt, n_steps):\n"
         "    x = np.empty(n_steps + 1)\n"
         "    v = np.empty(n_steps + 1)\n"
         "    x[0], v[0] = x0, v0\n"
         "    a = -x[0]\n"
         "    for n in range(n_steps):\n"
         "        x[n + 1] = x[n] + dt * v[n] + 0.5 * dt ** 2 * a\n"
         "        a_new = -x[n + 1]\n"
         "        v[n + 1] = v[n] + 0.5 * dt * (a + a_new)\n"
         "        a = a_new\n"
         "    return x, v\n\n"
         "dt = 0.1\n"
         "N = 5000\n"
         "xe, ve = euler(1.0, 0.0, dt, N)\n"
         "xv, vv = velocity_verlet(1.0, 0.0, dt, N)\n"
         "E_e = 0.5 * (ve ** 2 + xe ** 2)\n"
         "E_v = 0.5 * (vv ** 2 + xv ** 2)\n"
         "t = np.arange(N + 1) * dt\n\n"
         "fig, ax = plt.subplots(figsize=(7, 4))\n"
         "ax.semilogy(t, E_e, label='forward Euler')\n"
         "ax.semilogy(t, E_v, label='velocity-Verlet')\n"
         "ax.set_xlabel('time')\n"
         "ax.set_ylabel('total energy (log scale)')\n"
         "ax.set_title('Energy conservation: Euler vs Verlet')\n"
         "ax.legend()\n"
         "plt.show()\n"),
    md("## 7.2 Convergence of Verlet energy fluctuation with step size\n\n"
       "Halving $\\Delta t$ should quarter the amplitude of the bounded "
       "fluctuation — the hallmark of a second-order symplectic scheme."),
    code("import numpy as np\n\n"
         "for dt in [0.4, 0.2, 0.1, 0.05]:\n"
         "    N = int(200 / dt)\n"
         "    xv, vv = velocity_verlet(1.0, 0.0, dt, N)\n"
         "    E = 0.5 * (vv ** 2 + xv ** 2)\n"
         "    print(f'dt = {dt:5.2f}   peak-to-peak dE = {E.max() - E.min():.5f}')\n"),
]


# ---------------------------------------------------------------------------
# Chapter 8 — Statistical Mechanics from Simulations
# ---------------------------------------------------------------------------
CH08 = [
    header(
        "Chapter 8 · Statistical Mechanics from Simulations",
        "ch08-statmech",
        "Two short experiments: sampling a Boltzmann distribution with "
        "Metropolis Monte Carlo, and extracting a diffusion coefficient from a "
        "1D random walk via the Einstein relation.",
    ),
    md("## 8.1 Metropolis sampling of a 1D Boltzmann distribution\n\n"
       "Target: $\\rho(x) \\propto \\exp(-\\beta U(x))$ with a double-well "
       "potential. We plot the histogram against the analytic distribution."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "rng = np.random.default_rng(7)\n\n"
         "def U(x):\n"
         "    return (x ** 2 - 1) ** 2\n\n"
         "beta = 5.0\n"
         "x = 0.0\n"
         "samples = []\n"
         "for _ in range(200_000):\n"
         "    xp = x + rng.normal(0, 0.4)\n"
         "    if rng.random() < np.exp(-beta * (U(xp) - U(x))):\n"
         "        x = xp\n"
         "    samples.append(x)\n"
         "samples = np.array(samples[5000:])  # discard burn-in\n\n"
         "grid = np.linspace(-2, 2, 400)\n"
         "p = np.exp(-beta * U(grid))\n"
         "p /= 0.5 * np.sum((p[1:] + p[:-1]) * np.diff(grid))\n\n"
         "fig, ax = plt.subplots(figsize=(7, 4))\n"
         "ax.hist(samples, bins=80, density=True, alpha=0.6, label='Metropolis')\n"
         "ax.plot(grid, p, 'k', lw=2, label='analytic Boltzmann')\n"
         "ax.set_xlabel('x')\n"
         "ax.set_ylabel('rho(x)')\n"
         "ax.set_title('Double-well at beta = 5')\n"
         "ax.legend()\n"
         "plt.show()\n"),
    md("## 8.2 Diffusion from a random walk — the Einstein relation\n\n"
       "Mean-squared displacement of an unbiased walk grows linearly in time, "
       "with slope $2 D$ in 1D."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "rng = np.random.default_rng(11)\n"
         "n_walkers = 2000\n"
         "n_steps = 1500\n"
         "step_size = 1.0\n"
         "steps = rng.normal(0.0, step_size, size=(n_walkers, n_steps))\n"
         "x = np.cumsum(steps, axis=1)\n"
         "msd = (x ** 2).mean(axis=0)\n"
         "t = np.arange(1, n_steps + 1)\n"
         "D_fit = np.polyfit(t, msd, 1)[0] / 2\n\n"
         "fig, ax = plt.subplots(figsize=(6, 4))\n"
         "ax.plot(t, msd, label='measured MSD')\n"
         "ax.plot(t, 2 * D_fit * t, 'r--', label=f'fit: D = {D_fit:.3f}')\n"
         "ax.set_xlabel('step')\n"
         "ax.set_ylabel('<x^2>')\n"
         "ax.legend()\n"
         "ax.set_title('Einstein relation in 1D')\n"
         "plt.show()\n"),
]


# ---------------------------------------------------------------------------
# Chapter 11 — Active Learning and Bayesian Optimisation
# ---------------------------------------------------------------------------
CH11 = [
    header(
        "Chapter 11 · Active Learning and Bayesian Optimisation",
        "ch11-active",
        "A Gaussian process regressor with an RBF kernel, and a Bayesian "
        "optimisation loop on a toy 1D objective using Expected Improvement.",
    ),
    md("## 11.1 Gaussian process regression on noisy samples\n"),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n\n"
         "rng = np.random.default_rng(3)\n\n"
         "def rbf(x1, x2, length=0.4, sigma=1.0):\n"
         "    d2 = (x1[:, None] - x2[None, :]) ** 2\n"
         "    return sigma ** 2 * np.exp(-0.5 * d2 / length ** 2)\n\n"
         "def gp_predict(x_train, y_train, x_test, noise=0.05, length=0.4):\n"
         "    K = rbf(x_train, x_train, length) + noise ** 2 * np.eye(len(x_train))\n"
         "    Ks = rbf(x_train, x_test, length)\n"
         "    Kss = rbf(x_test, x_test, length)\n"
         "    L = np.linalg.cholesky(K)\n"
         "    alpha = np.linalg.solve(L.T, np.linalg.solve(L, y_train))\n"
         "    mean = Ks.T @ alpha\n"
         "    v = np.linalg.solve(L, Ks)\n"
         "    var = np.diag(Kss) - np.sum(v ** 2, axis=0)\n"
         "    return mean, np.sqrt(np.clip(var, 0, None))\n\n"
         "def true_f(x):\n"
         "    return np.sin(3 * x) + 0.3 * x\n\n"
         "x_train = rng.uniform(-2, 2, 8)\n"
         "y_train = true_f(x_train) + rng.normal(0, 0.1, len(x_train))\n"
         "x_test = np.linspace(-2.5, 2.5, 300)\n"
         "mu, std = gp_predict(x_train, y_train, x_test)\n\n"
         "fig, ax = plt.subplots(figsize=(7, 4))\n"
         "ax.plot(x_test, true_f(x_test), 'k--', lw=1, label='truth')\n"
         "ax.plot(x_train, y_train, 'ko', label='observations')\n"
         "ax.plot(x_test, mu, 'C0', label='GP mean')\n"
         "ax.fill_between(x_test, mu - 2 * std, mu + 2 * std, alpha=0.2, color='C0', label='+/- 2 sigma')\n"
         "ax.legend(fontsize=8)\n"
         "ax.set_title('Gaussian process regression')\n"
         "plt.show()\n"),
    md("## 11.2 Bayesian optimisation with Expected Improvement\n\n"
       "Maximise a 1D objective by alternating GP fits and EI-driven queries. "
       "Watch the algorithm zero in on the optimum within a handful of "
       "evaluations."),
    code("import numpy as np\n"
         "import matplotlib.pyplot as plt\n"
         "from math import erf, sqrt\n\n"
         "def phi(z):\n"
         "    return np.exp(-0.5 * z ** 2) / np.sqrt(2 * np.pi)\n\n"
         "def Phi(z):\n"
         "    return 0.5 * (1 + np.vectorize(lambda x: erf(x / sqrt(2)))(z))\n\n"
         "def ei(mu, std, best, xi=0.01):\n"
         "    std = np.clip(std, 1e-9, None)\n"
         "    z = (mu - best - xi) / std\n"
         "    return (mu - best - xi) * Phi(z) + std * phi(z)\n\n"
         "def objective(x):\n"
         "    return -((x - 0.4) ** 2) + 0.3 * np.sin(8 * x)\n\n"
         "rng = np.random.default_rng(0)\n"
         "x_grid = np.linspace(-1, 1, 400)\n"
         "x_obs = rng.uniform(-1, 1, 3)\n"
         "y_obs = objective(x_obs)\n\n"
         "fig, axes = plt.subplots(2, 3, figsize=(11, 5), sharex=True)\n"
         "for step, ax in enumerate(axes.ravel()):\n"
         "    mu, std = gp_predict(x_obs, y_obs, x_grid, noise=0.01, length=0.15)\n"
         "    acq = ei(mu, std, y_obs.max())\n"
         "    x_next = x_grid[int(np.argmax(acq))]\n"
         "    ax.plot(x_grid, objective(x_grid), 'k--', lw=0.8, label='truth')\n"
         "    ax.plot(x_grid, mu, 'C0')\n"
         "    ax.fill_between(x_grid, mu - 2 * std, mu + 2 * std, alpha=0.15, color='C0')\n"
         "    ax.plot(x_obs, y_obs, 'ko')\n"
         "    ax.axvline(x_next, color='C3', lw=0.8)\n"
         "    ax.set_title(f'after {step} suggestions')\n"
         "    x_obs = np.append(x_obs, x_next)\n"
         "    y_obs = np.append(y_obs, objective(x_next))\n"
         "fig.suptitle('Bayesian optimisation under Expected Improvement')\n"
         "fig.tight_layout()\n"
         "plt.show()\n"
         "print(f'best objective found: {y_obs.max():.4f} at x = {x_obs[int(np.argmax(y_obs))]:.4f}')\n"),
]


NOTEBOOK_SPECS = {
    "ch00-math.ipynb": CH00,
    "ch01-python.ipynb": CH01,
    "ch03b-solid-state.ipynb": CH035,
    "ch04-quantum.ipynb": CH04,
    "ch05-dft.ipynb": CH05,
    "ch07-md.ipynb": CH07,
    "ch08-statmech.ipynb": CH08,
    "ch11-active.ipynb": CH11,
}


def main() -> None:
    for name, cells in NOTEBOOK_SPECS.items():
        path = write_notebook(name, cells)
        print(f"wrote {path.relative_to(REPO_ROOT)}  ({len(cells)} cells)")


if __name__ == "__main__":
    main()
