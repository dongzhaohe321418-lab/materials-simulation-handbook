#!/usr/bin/env python3
"""Real Lennard-Jones MD simulation of argon producing g(r) for three phases.

Pure NumPy molecular dynamics on N=500 atoms in a cubic periodic box.
The Lennard-Jones potential with the standard Ar parameters
(sigma = 3.405 A, eps/kB = 119.8 K) is evolved with velocity Verlet at
dt = 5 fs.  Velocity rescaling (Berendsen-style) is used for 2000
equilibration steps at each target temperature, followed by 5000
production steps during which g(r) is sampled every 10 steps.

The resulting g(r) curves are characteristic of a solid (T = 50 K),
a liquid (T = 150 K), and a gas (T = 600 K).
"""
from __future__ import annotations

import os
import sys
import time
from dataclasses import dataclass

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from _style import setup, save, COLOURS

# ----------------------------------------------------------------------
# Physical constants and Ar parameters
# ----------------------------------------------------------------------
KB_EV: float = 8.617333262e-5          # eV / K
AMU_TO_EV_PS2_PER_A2: float = 1.0364e-4  # convert amu * A^2 / ps^2 -> eV
# i.e. KE [eV] = 0.5 * m[amu] * v^2[A^2/ps^2] * AMU_TO_EV_PS2_PER_A2

SIGMA_A: float = 3.405                  # Angstrom
EPS_K: float = 119.8                    # K  (eps / kB)
EPS_EV: float = EPS_K * KB_EV           # eV
MASS_AMU: float = 39.948                # Ar atomic mass


@dataclass
class MDState:
    """Container for current MD state in Angstroms, ps, eV units."""

    pos: np.ndarray                     # (N, 3) positions [A]
    vel: np.ndarray                     # (N, 3) velocities [A/ps]
    forces: np.ndarray                  # (N, 3) forces [eV/A]
    box: float                          # cubic box length [A]


# ----------------------------------------------------------------------
# Initial configuration: FCC lattice in a cubic box of n^3 conventional cells
# ----------------------------------------------------------------------
def fcc_positions(n_cells: int, a: float) -> np.ndarray:
    """Generate FCC atomic positions in a cubic box of side n_cells*a."""
    basis = np.array([[0.0, 0.0, 0.0],
                      [0.5, 0.5, 0.0],
                      [0.5, 0.0, 0.5],
                      [0.0, 0.5, 0.5]])
    coords: list[np.ndarray] = []
    for ix in range(n_cells):
        for iy in range(n_cells):
            for iz in range(n_cells):
                offset = np.array([ix, iy, iz], dtype=float)
                coords.append((basis + offset) * a)
    return np.vstack(coords)


def maxwell_boltzmann(n_atoms: int, mass_amu: float, T: float,
                      rng: np.random.Generator) -> np.ndarray:
    """Draw initial velocities from a Maxwell-Boltzmann distribution.

    Units: A / ps.  v ~ N(0, sqrt(kT/m)).
    """
    # kT in eV; mass in amu -> velocity in A/ps via the conversion constant
    sigma_v = np.sqrt(KB_EV * T / (mass_amu * AMU_TO_EV_PS2_PER_A2))
    v = rng.normal(0.0, sigma_v, size=(n_atoms, 3))
    v -= v.mean(axis=0)        # remove COM drift
    return v


# ----------------------------------------------------------------------
# Lennard-Jones forces with PBC, minimum-image, cutoff at 2.5 sigma with
# energy shift.  Vectorised with NumPy via pairwise broadcast.
# ----------------------------------------------------------------------
def lj_forces(pos: np.ndarray, box: float,
              rcut: float, sigma: float, eps: float) -> tuple[np.ndarray, float]:
    """Compute LJ forces and (shifted) potential energy under PBC.

    Returns
    -------
    forces : (N, 3) eV / A
    pe     : total potential energy in eV
    """
    n = pos.shape[0]
    # pairwise displacements (N, N, 3); disp[i, j] = r_i - r_j
    disp = pos[:, None, :] - pos[None, :, :]
    disp -= box * np.round(disp / box)
    r2 = np.sum(disp * disp, axis=-1)
    # avoid self-interaction and any out-of-cutoff division: clamp r2
    mask = (r2 < rcut * rcut) & (r2 > 0.0)
    safe_r2 = np.where(mask, r2, 1.0)         # placeholder for masked entries

    sig2 = sigma * sigma
    inv_r2 = sig2 / safe_r2
    inv_r6 = inv_r2 ** 3
    inv_r12 = inv_r6 ** 2

    # potential (with energy shift to zero at rcut)
    sr6_cut = (sigma / rcut) ** 6
    u_cut = 4.0 * eps * (sr6_cut ** 2 - sr6_cut)
    u_pair = np.where(mask, 4.0 * eps * (inv_r12 - inv_r6) - u_cut, 0.0)
    pe = 0.5 * u_pair.sum()

    # force magnitude factor:  F_i_from_j = -dU/dr * r_hat
    #   dU/dr = (24 eps / r) * ( (sigma/r)^6 - 2 (sigma/r)^12 )
    #   so the vector force on i from j is
    #     F = (24 eps / r^2) * (2 (sigma/r)^12 - (sigma/r)^6) * (r_i - r_j)
    fmag = np.where(mask,
                    24.0 * eps * (2.0 * inv_r12 - inv_r6) / safe_r2,
                    0.0)
    forces = (fmag[:, :, None] * disp).sum(axis=1)
    # convert: factor 24 eps / r^2 has units eV/A^2;
    # multiplying by disp [A] gives eV/A. Good.
    return forces, pe


# ----------------------------------------------------------------------
# Velocity-Verlet step
# ----------------------------------------------------------------------
def velocity_verlet(state: MDState, dt: float, mass_amu: float,
                    rcut: float, sigma: float, eps: float) -> float:
    """Advance state by one velocity-Verlet step. Returns potential energy."""
    # convert force to acceleration in A/ps^2:
    # a [A/ps^2] = F [eV/A] / m [amu] / AMU_TO_EV_PS2_PER_A2
    inv_m_scale = 1.0 / (mass_amu * AMU_TO_EV_PS2_PER_A2)

    state.vel += 0.5 * dt * state.forces * inv_m_scale
    state.pos += dt * state.vel
    # wrap into box
    state.pos -= state.box * np.floor(state.pos / state.box)

    state.forces, pe = lj_forces(state.pos, state.box, rcut, sigma, eps)
    state.vel += 0.5 * dt * state.forces * inv_m_scale
    return pe


def current_temperature(vel: np.ndarray, mass_amu: float) -> float:
    """Instantaneous temperature in K from the kinetic energy."""
    n = vel.shape[0]
    ke_ev = 0.5 * mass_amu * AMU_TO_EV_PS2_PER_A2 * np.sum(vel * vel)
    # equipartition: KE = (3/2) N kT  => T = 2 KE / (3 N kB)
    return 2.0 * ke_ev / (3.0 * n * KB_EV)


def berendsen_rescale(state: MDState, T_target: float, mass_amu: float,
                      tau_ratio: float = 0.1) -> None:
    """Scale velocities towards T_target (mild Berendsen-style)."""
    T = current_temperature(state.vel, mass_amu)
    if T <= 0.0:
        return
    lam = np.sqrt(1.0 + tau_ratio * (T_target / T - 1.0))
    state.vel *= lam


# ----------------------------------------------------------------------
# Radial distribution function accumulation
# ----------------------------------------------------------------------
def accumulate_rdf(pos: np.ndarray, box: float, n_bins: int,
                   r_max: float, hist: np.ndarray) -> None:
    """Accumulate pair distances into the histogram (in-place)."""
    n = pos.shape[0]
    dr = r_max / n_bins
    for i in range(n - 1):
        disp = pos[i + 1:] - pos[i]
        disp -= box * np.round(disp / box)
        d = np.sqrt(np.sum(disp * disp, axis=-1))
        d = d[d < r_max]
        idx = (d / dr).astype(np.int64)
        np.add.at(hist, idx, 1)


def finalise_rdf(hist: np.ndarray, n_frames: int, n_atoms: int,
                 box: float, r_max: float) -> tuple[np.ndarray, np.ndarray]:
    """Normalise the pair histogram into g(r)."""
    n_bins = hist.size
    dr = r_max / n_bins
    r = (np.arange(n_bins) + 0.5) * dr
    shell_vol = 4.0 * np.pi * r * r * dr
    rho = n_atoms / box ** 3
    norm = shell_vol * rho * n_atoms * n_frames / 2.0
    g = hist / norm
    return r, g


# ----------------------------------------------------------------------
# Top-level driver
# ----------------------------------------------------------------------
def run_phase(T_target: float, box: float, n_cells: int, a_fcc: float,
              n_equil: int, n_prod: int, dt: float,
              rcut: float, sigma: float, eps: float,
              rng: np.random.Generator,
              label: str) -> tuple[np.ndarray, np.ndarray]:
    """Run one MD trajectory at T_target and return its g(r)."""
    pos = fcc_positions(n_cells, a_fcc)
    pos -= box * np.floor(pos / box)
    vel = maxwell_boltzmann(pos.shape[0], MASS_AMU, T_target, rng)
    forces, _ = lj_forces(pos, box, rcut, sigma, eps)
    state = MDState(pos=pos, vel=vel, forces=forces, box=box)

    # equilibration with periodic rescaling
    t0 = time.time()
    for step in range(n_equil):
        velocity_verlet(state, dt, MASS_AMU, rcut, sigma, eps)
        if step % 10 == 0:
            berendsen_rescale(state, T_target, MASS_AMU, tau_ratio=0.1)

    # production with continued (much milder) rescaling and rdf sampling
    n_bins = 200
    r_max = box / 2.0 - 1e-6
    hist = np.zeros(n_bins, dtype=np.int64)
    n_frames = 0
    for step in range(n_prod):
        velocity_verlet(state, dt, MASS_AMU, rcut, sigma, eps)
        if step % 10 == 0:
            berendsen_rescale(state, T_target, MASS_AMU, tau_ratio=0.02)
            accumulate_rdf(state.pos, state.box, n_bins, r_max, hist)
            n_frames += 1

    r, g = finalise_rdf(hist, n_frames, state.pos.shape[0], state.box, r_max)
    T_final = current_temperature(state.vel, MASS_AMU)
    elapsed = time.time() - t0
    print(f"  {label}: target {T_target:.0f} K, final {T_final:.1f} K, "
          f"{n_frames} samples, {elapsed:.1f} s")
    return r, g


def main() -> None:
    setup()
    rng = np.random.default_rng(20260512)

    # FCC initial geometry: 5x5x5 cells = 500 atoms.
    n_cells = 5
    n_atoms = 4 * n_cells ** 3
    assert n_atoms == 500

    # MD parameters
    dt = 5e-3        # 5 fs in ps
    n_equil = 2000
    n_prod = 5000
    rcut = 2.5 * SIGMA_A

    # Each phase is run at BOTH its characteristic temperature AND its
    # characteristic density. Reusing one density for all three would
    # just give the same fluid at three temperatures — a 600 K run at
    # liquid density is a supercritical fluid, not a vapour. The reduced
    # densities rho* = rho * sigma^3 of ~1.0 / 0.84 / 0.05 are the
    # textbook solid / triple-point-liquid / dilute-gas values for the
    # Lennard-Jones system.
    sigma3 = SIGMA_A ** 3
    runs = [
        # label,                T (K),  rho* (reduced), colour
        ("Solid (T = 50 K)",     50.0,  1.00, COLOURS["blue"]),
        ("Liquid (T = 150 K)",  150.0,  0.84, COLOURS["red"]),
        ("Gas (T = 600 K)",     600.0,  0.05, COLOURS["green"]),
    ]

    print(f"N = {n_atoms} atoms, sigma = {SIGMA_A} A, eps/kB = {EPS_K} K")

    fig, ax = plt.subplots(figsize=(9, 5))
    max_half_box = 0.0
    for label, T, rho_star, colour in runs:
        rho = rho_star / sigma3                  # number density in A^-3
        box = (n_atoms / rho) ** (1.0 / 3.0)
        a_fcc = box / n_cells
        max_half_box = max(max_half_box, box / 2.0)
        print(f"  {label}: rho* = {rho_star:.2f}, rho = {rho:.4f} A^-3, "
              f"box = {box:.2f} A")
        r, g = run_phase(
            T_target=T, box=box, n_cells=n_cells, a_fcc=a_fcc,
            n_equil=n_equil, n_prod=n_prod, dt=dt,
            rcut=rcut, sigma=SIGMA_A, eps=EPS_EV,
            rng=rng, label=label,
        )
        ax.plot(r / SIGMA_A, g, color=colour, lw=2.0, label=label)
        # quick numerical sanity print: first peak
        if T == 150.0:
            idx_peak = np.argmax(g)
            print(f"  liquid first-peak at r/sigma = {r[idx_peak]/SIGMA_A:.3f}, "
                  f"g_max = {g[idx_peak]:.2f}")

    ax.axhline(1.0, color=COLOURS["grey"], lw=0.7, ls=":")
    ax.set_xlabel(r"$r / \sigma$")
    ax.set_ylabel(r"$g(r)$")
    ax.set_title(r"Radial distribution function from LJ-MD of 500 Ar atoms")
    ax.set_xlim(0.0, 6.0)
    ax.set_ylim(0.0, None)
    ax.legend(frameon=False, loc="upper right")
    ax.grid(alpha=0.3, linestyle="--")

    out = os.path.join(
        os.path.dirname(__file__),
        "../../docs/assets/figures/ch07/fig_rdf_real.png",
    )
    save(fig, out)


if __name__ == "__main__":
    main()
