#!/usr/bin/env python3
"""Real silicon band structure from the sp3d5s* tight-binding model.

Slater-Koster parameters for Si from
    Jancu, Scholz, Beltram, Bassani, PRB 57, 6493 (1998).

Bands are diagonalised along the standard FCC path L -> Gamma -> X -> U,K
-> Gamma at ~50 points per segment.  The output is the canonical Si band
structure with an indirect gap of ~1.17 eV, valence-band maximum at Gamma
and conduction-band minimum about 85% of the way from Gamma to X.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS

# ----------------------------------------------------------------------
# Jancu et al. sp3d5s* parameters for Si (eV)
# ----------------------------------------------------------------------
# On-site energies
E_s = -2.15168
E_p = 4.22925
E_sstar = 19.11650
E_d = 13.78950

# Two-centre Slater-Koster integrals (eV)
V_sss = -1.95933
V_sps = 3.02562
V_sds = -2.28485
V_pps = 4.10364
V_ppp = -1.51801
V_pds = -1.35554
V_pdp = 2.38479
V_dds = -1.68136
V_ddp = 2.58880
V_ddd = -1.81400
V_sstar_sstar_s = -4.24135
V_sstar_p_s = 3.15565
V_sstar_d_s = -0.80993
V_s_sstar_s = -1.52230   # between s on one atom and s* on the other

# Conventional Si lattice constant (Angstrom).  Bands are insensitive to
# this within the tight-binding model (it only sets the bond direction
# pattern), so the exact value matters mostly for the k-point coordinates.
A_LATT: float = 5.431

# Orbital ordering within an atom
ORB = ["s", "px", "py", "pz", "xy", "yz", "zx", "x2-y2", "3z2-r2", "s*"]
N_ORB = len(ORB)        # 10
N_ATOMS_PER_CELL = 2
N_BASIS = N_ORB * N_ATOMS_PER_CELL  # 20


# ----------------------------------------------------------------------
# Slater-Koster table: hopping between orbitals on neighbouring atoms
# along a bond with direction cosines (l, m, n).
# Reference: Slater & Koster, Phys. Rev. 94, 1498 (1954), Table I.
# All entries are <orbital_alpha (atom A) | H | orbital_beta (atom B)>.
# ----------------------------------------------------------------------
def sk_hop(alpha: str, beta: str, l: float, m: float, n: float) -> float:
    """Slater-Koster two-centre integral for one bond direction.

    alpha, beta : orbital labels from ORB.
    l, m, n     : direction cosines pointing from atom A to atom B.
    """
    sqrt3 = np.sqrt(3.0)

    # ----- s-s -----
    if alpha == "s" and beta == "s":
        return V_sss
    # ----- s-s* and s*-s -----
    if (alpha, beta) in (("s", "s*"), ("s*", "s")):
        return V_s_sstar_s
    # ----- s*-s* -----
    if alpha == "s*" and beta == "s*":
        return V_sstar_sstar_s

    # ----- s - p   (sign convention: pos l for px when alpha=s, beta=px) -----
    if alpha == "s" and beta in ("px", "py", "pz"):
        dc = {"px": l, "py": m, "pz": n}[beta]
        return dc * V_sps
    if alpha in ("px", "py", "pz") and beta == "s":
        dc = {"px": l, "py": m, "pz": n}[alpha]
        return -dc * V_sps     # <p|H|s> = -<s|H|p> for the same bond

    # ----- s* - p -----
    if alpha == "s*" and beta in ("px", "py", "pz"):
        dc = {"px": l, "py": m, "pz": n}[beta]
        return dc * V_sstar_p_s
    if alpha in ("px", "py", "pz") and beta == "s*":
        dc = {"px": l, "py": m, "pz": n}[alpha]
        return -dc * V_sstar_p_s

    # ----- s - d -----
    if alpha == "s" and beta.startswith(("xy", "yz", "zx", "x2", "3z2")):
        return _s_d(beta, l, m, n, V_sds)
    if alpha.startswith(("xy", "yz", "zx", "x2", "3z2")) and beta == "s":
        return _s_d(alpha, l, m, n, V_sds)

    # ----- s* - d -----
    if alpha == "s*" and beta.startswith(("xy", "yz", "zx", "x2", "3z2")):
        return _s_d(beta, l, m, n, V_sstar_d_s)
    if alpha.startswith(("xy", "yz", "zx", "x2", "3z2")) and beta == "s*":
        return _s_d(alpha, l, m, n, V_sstar_d_s)

    # ----- p - p -----
    p_set = {"px", "py", "pz"}
    if alpha in p_set and beta in p_set:
        dc = {"px": l, "py": m, "pz": n}
        a = dc[alpha]; b = dc[beta]
        if alpha == beta:
            return a * a * V_pps + (1.0 - a * a) * V_ppp
        return a * b * (V_pps - V_ppp)

    # ----- p - d -----
    if alpha in p_set and beta.startswith(("xy", "yz", "zx", "x2", "3z2")):
        return _p_d(alpha, beta, l, m, n)
    if alpha.startswith(("xy", "yz", "zx", "x2", "3z2")) and beta in p_set:
        # <d|H|p>(R) = <p|H|d>(-R) under sign flip of direction cosines.
        return _p_d(beta, alpha, -l, -m, -n)

    # ----- d - d -----
    if (alpha.startswith(("xy", "yz", "zx", "x2", "3z2")) and
            beta.startswith(("xy", "yz", "zx", "x2", "3z2"))):
        return _d_d(alpha, beta, l, m, n)

    raise ValueError(f"unhandled SK pair: {alpha}, {beta}")


def _s_d(d_orb: str, l: float, m: float, n: float, V: float) -> float:
    """Slater-Koster <s | H | d_orb> along (l, m, n)."""
    sqrt3 = np.sqrt(3.0)
    if d_orb == "xy":
        return sqrt3 * l * m * V
    if d_orb == "yz":
        return sqrt3 * m * n * V
    if d_orb == "zx":
        return sqrt3 * n * l * V
    if d_orb == "x2-y2":
        return 0.5 * sqrt3 * (l * l - m * m) * V
    if d_orb == "3z2-r2":
        return (n * n - 0.5 * (l * l + m * m)) * V
    raise ValueError(d_orb)


def _p_d(p_orb: str, d_orb: str, l: float, m: float, n: float) -> float:
    """Slater-Koster <p_orb (A) | H | d_orb (B)> along (l, m, n) from A to B.

    Uses the standard table from Slater & Koster (1954) with sigma and pi
    bond strengths V_pds, V_pdp.
    """
    sqrt3 = np.sqrt(3.0)
    pidx = {"px": l, "py": m, "pz": n}[p_orb]
    if d_orb == "xy":
        # E_{x,xy} = sqrt3 l^2 m V_pds + m (1 - 2 l^2) V_pdp
        if p_orb == "px":
            return sqrt3 * l * l * m * V_pds + m * (1.0 - 2.0 * l * l) * V_pdp
        if p_orb == "py":
            return sqrt3 * m * m * l * V_pds + l * (1.0 - 2.0 * m * m) * V_pdp
        if p_orb == "pz":
            return sqrt3 * l * m * n * V_pds - 2.0 * l * m * n * V_pdp
    if d_orb == "yz":
        if p_orb == "px":
            return sqrt3 * l * m * n * V_pds - 2.0 * l * m * n * V_pdp
        if p_orb == "py":
            return sqrt3 * m * m * n * V_pds + n * (1.0 - 2.0 * m * m) * V_pdp
        if p_orb == "pz":
            return sqrt3 * n * n * m * V_pds + m * (1.0 - 2.0 * n * n) * V_pdp
    if d_orb == "zx":
        if p_orb == "px":
            return sqrt3 * l * l * n * V_pds + n * (1.0 - 2.0 * l * l) * V_pdp
        if p_orb == "py":
            return sqrt3 * l * m * n * V_pds - 2.0 * l * m * n * V_pdp
        if p_orb == "pz":
            return sqrt3 * n * n * l * V_pds + l * (1.0 - 2.0 * n * n) * V_pdp
    if d_orb == "x2-y2":
        # standard SK entries
        if p_orb == "px":
            return (0.5 * sqrt3 * l * (l * l - m * m) * V_pds
                    + l * (1.0 - l * l + m * m) * V_pdp)
        if p_orb == "py":
            return (0.5 * sqrt3 * m * (l * l - m * m) * V_pds
                    - m * (1.0 + l * l - m * m) * V_pdp)
        if p_orb == "pz":
            return (0.5 * sqrt3 * n * (l * l - m * m) * V_pds
                    - n * (l * l - m * m) * V_pdp)
    if d_orb == "3z2-r2":
        if p_orb == "px":
            return (l * (n * n - 0.5 * (l * l + m * m)) * V_pds
                    - sqrt3 * l * n * n * V_pdp)
        if p_orb == "py":
            return (m * (n * n - 0.5 * (l * l + m * m)) * V_pds
                    - sqrt3 * m * n * n * V_pdp)
        if p_orb == "pz":
            return (n * (n * n - 0.5 * (l * l + m * m)) * V_pds
                    + sqrt3 * n * (l * l + m * m) * V_pdp)
    raise ValueError(f"_p_d unhandled: {p_orb}, {d_orb}")


def _d_d(d1: str, d2: str, l: float, m: float, n: float) -> float:
    """Slater-Koster <d1 | H | d2> along (l, m, n)."""
    sqrt3 = np.sqrt(3.0)
    # Use the symmetric table from Slater-Koster (1954).
    # Diagonal d-d:
    if d1 == d2 == "xy":
        return (3.0 * l * l * m * m * V_dds
                + (l * l + m * m - 4.0 * l * l * m * m) * V_ddp
                + (n * n + l * l * m * m) * V_ddd)
    if d1 == d2 == "yz":
        return (3.0 * m * m * n * n * V_dds
                + (m * m + n * n - 4.0 * m * m * n * n) * V_ddp
                + (l * l + m * m * n * n) * V_ddd)
    if d1 == d2 == "zx":
        return (3.0 * n * n * l * l * V_dds
                + (n * n + l * l - 4.0 * n * n * l * l) * V_ddp
                + (m * m + n * n * l * l) * V_ddd)
    if d1 == d2 == "x2-y2":
        return (0.75 * (l * l - m * m) ** 2 * V_dds
                + (l * l + m * m - (l * l - m * m) ** 2) * V_ddp
                + (n * n + 0.25 * (l * l - m * m) ** 2) * V_ddd)
    if d1 == d2 == "3z2-r2":
        f = n * n - 0.5 * (l * l + m * m)
        return (f * f * V_dds
                + 3.0 * n * n * (l * l + m * m) * V_ddp
                + 0.75 * (l * l + m * m) ** 2 * V_ddd)
    # Off-diagonal: symmetric, swap if needed
    if d1 == "xy" and d2 == "yz":
        return (3.0 * l * m * m * n * V_dds
                + l * n * (1.0 - 4.0 * m * m) * V_ddp
                + l * n * (m * m - 1.0) * V_ddd)
    if d1 == "xy" and d2 == "zx":
        return (3.0 * l * l * m * n * V_dds
                + m * n * (1.0 - 4.0 * l * l) * V_ddp
                + m * n * (l * l - 1.0) * V_ddd)
    if d1 == "yz" and d2 == "zx":
        return (3.0 * l * m * n * n * V_dds
                + l * m * (1.0 - 4.0 * n * n) * V_ddp
                + l * m * (n * n - 1.0) * V_ddd)
    # symmetric counterparts:
    if (d1, d2) in (("yz", "xy"), ("zx", "xy"), ("zx", "yz")):
        return _d_d(d2, d1, l, m, n)

    if d1 == "xy" and d2 == "x2-y2":
        return (1.5 * l * m * (l * l - m * m) * V_dds
                + 2.0 * l * m * (m * m - l * l) * V_ddp
                + 0.5 * l * m * (l * l - m * m) * V_ddd)
    if d1 == "yz" and d2 == "x2-y2":
        return (1.5 * m * n * (l * l - m * m) * V_dds
                - m * n * (1.0 + 2.0 * (l * l - m * m)) * V_ddp
                + m * n * (1.0 + 0.5 * (l * l - m * m)) * V_ddd)
    if d1 == "zx" and d2 == "x2-y2":
        return (1.5 * n * l * (l * l - m * m) * V_dds
                + n * l * (1.0 - 2.0 * (l * l - m * m)) * V_ddp
                - n * l * (1.0 - 0.5 * (l * l - m * m)) * V_ddd)
    if (d1, d2) in (("x2-y2", "xy"), ("x2-y2", "yz"), ("x2-y2", "zx")):
        return _d_d(d2, d1, l, m, n)

    if d1 == "xy" and d2 == "3z2-r2":
        return (sqrt3 * l * m * (n * n - 0.5 * (l * l + m * m)) * V_dds
                - 2.0 * sqrt3 * l * m * n * n * V_ddp
                + 0.5 * sqrt3 * l * m * (1.0 + n * n) * V_ddd)
    if d1 == "yz" and d2 == "3z2-r2":
        return (sqrt3 * m * n * (n * n - 0.5 * (l * l + m * m)) * V_dds
                + sqrt3 * m * n * (l * l + m * m - n * n) * V_ddp
                - 0.5 * sqrt3 * m * n * (l * l + m * m) * V_ddd)
    if d1 == "zx" and d2 == "3z2-r2":
        return (sqrt3 * l * n * (n * n - 0.5 * (l * l + m * m)) * V_dds
                + sqrt3 * l * n * (l * l + m * m - n * n) * V_ddp
                - 0.5 * sqrt3 * l * n * (l * l + m * m) * V_ddd)
    if (d1, d2) in (("3z2-r2", "xy"), ("3z2-r2", "yz"), ("3z2-r2", "zx")):
        return _d_d(d2, d1, l, m, n)

    if d1 == "x2-y2" and d2 == "3z2-r2":
        return (0.5 * sqrt3 * (l * l - m * m) *
                (n * n - 0.5 * (l * l + m * m)) * V_dds
                + sqrt3 * n * n * (m * m - l * l) * V_ddp
                + 0.25 * sqrt3 * (1.0 + n * n) * (l * l - m * m) * V_ddd)
    if d1 == "3z2-r2" and d2 == "x2-y2":
        return _d_d(d2, d1, l, m, n)

    raise ValueError(f"_d_d unhandled: {d1}, {d2}")


# ----------------------------------------------------------------------
# Bloch Hamiltonian construction
# ----------------------------------------------------------------------
def hop_block(l: float, m: float, n: float) -> np.ndarray:
    """10x10 hopping block from atom A to atom B along direction (l,m,n)."""
    H = np.zeros((N_ORB, N_ORB))
    for i, a in enumerate(ORB):
        for j, b in enumerate(ORB):
            H[i, j] = sk_hop(a, b, l, m, n)
    return H


@dataclass
class TBModel:
    """Pre-computed onsite blocks and the four NN hopping blocks."""

    onsite: np.ndarray              # (10,10) diagonal of orbital energies
    nn_blocks: list[np.ndarray]      # four (10,10) NN hopping matrices
    nn_vecs: np.ndarray              # (4, 3) NN displacement vectors (A)


def build_model() -> TBModel:
    """Assemble onsite block and four NN hopping blocks for Si diamond."""
    onsite = np.diag([
        E_s, E_p, E_p, E_p,
        E_d, E_d, E_d, E_d, E_d,
        E_sstar,
    ])
    # Nearest-neighbour vectors from atom A at (0,0,0) to atom B
    # at (a/4)(1,1,1) plus the three equivalents.  In conventional units:
    a = A_LATT
    d = a / 4.0
    nn_vecs = np.array([
        [ d,  d,  d],
        [ d, -d, -d],
        [-d,  d, -d],
        [-d, -d,  d],
    ])
    # direction cosines (unit vectors along each bond)
    nn_blocks = []
    for v in nn_vecs:
        norm = np.linalg.norm(v)
        l, m, n = v / norm
        nn_blocks.append(hop_block(l, m, n))
    return TBModel(onsite=onsite, nn_blocks=nn_blocks, nn_vecs=nn_vecs)


def bloch_hamiltonian(k: np.ndarray, model: TBModel) -> np.ndarray:
    """Construct the 20x20 Bloch Hamiltonian at crystal momentum k.

    k is in Cartesian reciprocal-space coordinates (units 1/A).
    """
    H = np.zeros((N_BASIS, N_BASIS), dtype=complex)
    # On-site
    H[:N_ORB, :N_ORB] = model.onsite
    H[N_ORB:, N_ORB:] = model.onsite
    # Inter-sublattice hopping with the Bloch phase
    H_AB = np.zeros((N_ORB, N_ORB), dtype=complex)
    for v, block in zip(model.nn_vecs, model.nn_blocks):
        phase = np.exp(1j * np.dot(k, v))
        H_AB += block * phase
    H[:N_ORB, N_ORB:] = H_AB
    H[N_ORB:, :N_ORB] = H_AB.conj().T
    return H


# ----------------------------------------------------------------------
# k-path along L -> Gamma -> X -> U,K -> Gamma in 2*pi/a units
# ----------------------------------------------------------------------
def build_kpath(n_per_seg: int = 50) -> tuple[np.ndarray, np.ndarray,
                                              list[tuple[float, str]]]:
    """Return Cartesian k-points (1/A), cumulative k-distance, and labels."""
    two_pi_a = 2.0 * np.pi / A_LATT
    # high-symmetry points in Cartesian (units of 2pi/a)
    pts = {
        "L":     np.array([0.5, 0.5, 0.5]),
        "G":     np.array([0.0, 0.0, 0.0]),
        "X":     np.array([0.0, 1.0, 0.0]),
        "U":     np.array([0.25, 1.0, 0.25]),
        "K":     np.array([0.75, 0.75, 0.0]),
    }
    # path: L - G - X - U / K - G (with jump from U to K)
    segs = [("L", "G"), ("G", "X"), ("X", "U"), ("K", "G")]
    labels: list[tuple[float, str]] = []
    k_pts: list[np.ndarray] = []
    k_cum: list[float] = []
    cum = 0.0
    for i, (a, b) in enumerate(segs):
        start = pts[a] * two_pi_a
        end = pts[b] * two_pi_a
        labels.append((cum, a if a != "K" else "U,K"))
        seg_pts = np.linspace(start, end, n_per_seg, endpoint=False)
        for p in seg_pts:
            k_pts.append(p)
            k_cum.append(cum)
            cum_diff = np.linalg.norm(end - start) / n_per_seg
            cum += cum_diff
    # append the final endpoint of the last segment
    k_pts.append(pts["G"] * two_pi_a)
    k_cum.append(cum)
    labels.append((cum, r"$\Gamma$"))
    # Replace G with Gamma symbol in labels
    labels = [(x, (r"$\Gamma$" if name == "G" else name)) for (x, name) in labels]
    return np.array(k_pts), np.array(k_cum), labels


def main() -> None:
    setup()
    model = build_model()
    k_pts, k_cum, labels = build_kpath(n_per_seg=50)
    n_k = k_pts.shape[0]
    bands = np.zeros((n_k, N_BASIS))
    for i, k in enumerate(k_pts):
        H = bloch_hamiltonian(k, model)
        w = np.linalg.eigvalsh(H)
        bands[i, :] = w

    # Determine VBM and CBM.  Si has 4 valence electrons per atom * 2 atoms = 8
    # electrons per primitive cell, so 4 bands are doubly-occupied.  The
    # valence-band top is the maximum of band index 3 (zero-indexed) over k,
    # the conduction-band minimum is the minimum of band index 4 over k.
    VB_top = bands[:, 3].max()
    CB_bot = bands[:, 4].min()
    indirect_gap = CB_bot - VB_top
    # Locate CBM
    k_cbm = k_cum[np.argmin(bands[:, 4])]
    # Direct gap at Gamma: indices of Gamma points are at labels positions [1] and [-1]
    gamma_idx = np.argmin(np.linalg.norm(k_pts - 0.0, axis=1))
    direct_gap = bands[gamma_idx, 4] - bands[gamma_idx, 3]
    print(f"Indirect gap = {indirect_gap:.3f} eV "
          f"(CBM at k = {k_cbm:.3f}, VBM at Gamma)")
    print(f"Direct gap at Gamma = {direct_gap:.3f} eV")

    # ------- Plot -------
    fig, ax = plt.subplots(figsize=(9, 6))
    # Plot the lowest 16 bands relative to the VBM
    e_shift = VB_top
    n_plot = 16
    for b in range(n_plot):
        col = COLOURS["blue"] if b <= 3 else COLOURS["red"]
        ax.plot(k_cum, bands[:, b] - e_shift, color=col, lw=1.5,
                label=("valence" if b == 0 else
                       ("conduction" if b == 4 else None)))

    # Fermi level: top of valence band
    ax.axhline(0.0, color="k", lw=0.8, ls="--", label=r"$E_F$ (VBM)")

    # High-symmetry tick marks
    tick_x = [labels[0][0], labels[1][0], labels[2][0], labels[3][0], labels[4][0]]
    tick_l = [labels[0][1], labels[1][1], labels[2][1], labels[3][1], labels[4][1]]
    ax.set_xticks(tick_x)
    ax.set_xticklabels(tick_l, fontsize=13)
    for tx in tick_x:
        ax.axvline(tx, color=COLOURS["grey"], lw=0.6, ls=":")

    ax.set_xlim(k_cum.min(), k_cum.max())
    ax.set_ylim(-13, 8)
    ax.set_ylabel(r"$E - E_{\rm VBM}$ (eV)")
    ax.set_title(r"Silicon band structure (sp$^3$d$^5$s$^*$ tight binding, Jancu 1998)")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(alpha=0.3, linestyle="--", axis="y")

    out = os.path.join(
        os.path.dirname(__file__),
        "../../docs/assets/figures/ch06/fig_si_bands_real.png",
    )
    save(fig, out)


if __name__ == "__main__":
    main()
