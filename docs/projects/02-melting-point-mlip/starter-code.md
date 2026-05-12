# Starter code — Melting point of copper via MLIP-driven MD

All scripts are Python 3.11 with type hints, organised as a package
`project02/`. Together they implement: AIMD data generation, training
config, parity validation, coexistence-cell construction, interface
analysis, and a small CLI to orchestrate the four phases.

---

## `project02/config.py` — central configuration

```python
"""Shared configuration for the Cu melting-point project."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DFTConfig:
    """DFT settings for QE."""
    pseudo_cu: str = "Cu.pbe-spn-rrkjus_psl.1.0.0.UPF"
    pseudo_dir: Path = Path.home() / "pseudo" / "sssp_efficiency_PBE"
    ecutwfc: float = 60.0       # Ry
    ecutrho: float = 720.0      # Ry
    smearing: str = "mv"        # Marzari-Vanderbilt
    degauss: float = 0.01       # Ry
    mixing_beta: float = 0.3
    conv_thr: float = 1.0e-7


@dataclass(frozen=True)
class MACEConfig:
    """MACE training hyperparameters."""
    r_max: float = 5.0
    max_L: int = 2
    correlation: int = 3
    hidden_irreps: str = "128x0e+128x1o"
    num_interactions: int = 2
    batch_size: int = 8
    epochs: int = 200
    lr: float = 0.01
    weight_decay: float = 5.0e-7
    energy_weight: float = 1.0
    forces_weight: float = 100.0
    stress_weight: float = 1.0


@dataclass(frozen=True)
class MDConfig:
    """MD parameters used in equilibration and coexistence."""
    timestep_fs: float = 1.0
    thermostat_tau_ps: float = 0.1
    barostat_tau_ps: float = 1.0
    pressure_bar: float = 1.0
    equil_ps: float = 30.0
    sample_ps: float = 170.0    # so total = 200 ps


@dataclass(frozen=True)
class ProjectConfig:
    a0_cu: float = 3.615        # Å, PBE eq. lattice parameter
    work_dir: Path = Path("./runs").resolve()
    dft: DFTConfig = DFTConfig()
    mace: MACEConfig = MACEConfig()
    md: MDConfig = MDConfig()


CFG = ProjectConfig()
```

---

## `project02/build_cells.py` — bulk and coexistence cells

```python
"""Construct the bulk fcc Cu cell and the two-phase coexistence cell."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import write

from .config import CFG


def fcc_cu(n: tuple[int, int, int]) -> Atoms:
    """n[i] copies along axis i of the 4-atom conventional fcc cell."""
    conv = bulk("Cu", crystalstructure="fcc", a=CFG.a0_cu, cubic=True)
    return conv.repeat(n)


def coexistence_initial(n_xy: int = 8, n_z: int = 32) -> Atoms:
    """Return the elongated coexistence cell (all solid; melting handled later)."""
    return fcc_cu((n_xy, n_xy, n_z))


def split_solid_liquid_mask(atoms: Atoms) -> np.ndarray:
    """Boolean mask: True for atoms in the half to be melted (z > L_z/2)."""
    z = atoms.positions[:, 2]
    return z > 0.5 * atoms.cell[2, 2]


if __name__ == "__main__":
    out = Path("structures")
    out.mkdir(exist_ok=True)
    write(out / "Cu_864.xyz", fcc_cu((6, 6, 6)))
    write(out / "Cu_coex_init.xyz", coexistence_initial())
    print(f"864-atom cell and coexistence cell ({8 * 8 * 32 * 4} atoms) written.")
```

---

## `project02/aimd_qe.py` — generate DFT-MD training data

```python
"""Build QE input files for the AIMD training trajectories."""
from __future__ import annotations

from pathlib import Path
from typing import Mapping

from ase import Atoms
from ase.calculators.espresso import Espresso

from .config import CFG


def aimd_calc(prefix: str,
              temperature_K: float,
              nstep: int,
              workdir: Path) -> Espresso:
    """Return an Espresso calculator configured for BO-MD on Cu."""
    workdir.mkdir(parents=True, exist_ok=True)
    input_data: Mapping[str, Mapping[str, object]] = {
        "control": {
            "calculation": "md",
            "prefix": prefix,
            "pseudo_dir": str(CFG.dft.pseudo_dir),
            "outdir": str(workdir),
            "nstep": nstep,
            "dt": 41.34,  # 1 fs in Ry atomic units
            "verbosity": "low",
        },
        "system": {
            "ecutwfc": CFG.dft.ecutwfc,
            "ecutrho": CFG.dft.ecutrho,
            "occupations": "smearing",
            "smearing": CFG.dft.smearing,
            "degauss": CFG.dft.degauss,
        },
        "electrons": {
            "conv_thr": CFG.dft.conv_thr,
            "mixing_beta": CFG.dft.mixing_beta,
        },
        "ions": {
            "ion_dynamics": "verlet",
            "ion_temperature": "rescale-v",
            "tempw": temperature_K,
            "nraise": 50,
        },
    }
    return Espresso(
        input_data=input_data,
        pseudopotentials={"Cu": CFG.dft.pseudo_cu},
        kpts=(2, 2, 2),
        koffset=(0, 0, 0),
        directory=str(workdir),
    )


def run_one(atoms: Atoms, prefix: str, temperature_K: float, ps: float) -> None:
    """Drive a single AIMD trajectory and let QE write the .pwo log."""
    nstep = int(ps * 1000)  # 1 fs steps
    atoms = atoms.copy()
    atoms.calc = aimd_calc(prefix, temperature_K, nstep,
                           CFG.work_dir / "aimd" / prefix)
    atoms.get_potential_energy()  # triggers the run
```

---

## `project02/build_dataset.py` — convert QE outputs to extended-XYZ

```python
"""Read AIMD trajectories and write a single extxyz training file."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase.io import read, write


def collect_frames(label_to_path: dict[str, Path],
                   stride: int = 20) -> list:
    frames = []
    for label, path in label_to_path.items():
        traj = read(path, index=f"::{stride}")
        for f in traj:
            f.info["config_type"] = label
        frames.extend(traj)
    return frames


def write_split(frames: list,
                out_dir: Path,
                fractions: tuple[float, float, float] = (0.8, 0.1, 0.1),
                seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    idx = np.arange(len(frames))
    rng.shuffle(idx)
    n_train = int(fractions[0] * len(frames))
    n_val = int(fractions[1] * len(frames))
    splits = {
        "train": [frames[i] for i in idx[:n_train]],
        "val": [frames[i] for i in idx[n_train:n_train + n_val]],
        "test": [frames[i] for i in idx[n_train + n_val:]],
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, fs in splits.items():
        write(out_dir / f"cu_{name}.xyz", fs)
        print(f"{name}: {len(fs)} frames")


if __name__ == "__main__":
    base = Path("runs/aimd")
    label_paths = {
        "solid_300K": base / "cu_solid_300K" / "espresso.pwo",
        "solid_800K": base / "cu_solid_800K" / "espresso.pwo",
        "liquid_1500K": base / "cu_liquid_1500K" / "espresso.pwo",
    }
    frames = collect_frames(label_paths)
    write_split(frames, Path("data"))
```

---

## `project02/train_mace.py` — train the MLIP

```python
"""Thin wrapper around mace-torch's CLI for reproducible training."""
from __future__ import annotations

import subprocess
from pathlib import Path

from .config import CFG


def train(name: str = "cu_mace_v1",
          train_xyz: Path = Path("data/cu_train.xyz"),
          val_xyz: Path = Path("data/cu_val.xyz"),
          e0_cu_ev: float = -1234.567,  # placeholder; compute for your pseudopotential
          device: str = "cuda") -> None:
    cmd = [
        "mace_run_train",
        f"--name={name}",
        f"--train_file={train_xyz}",
        f"--valid_file={val_xyz}",
        f"--E0s={{\"Cu\": {e0_cu_ev}}}",
        f"--r_max={CFG.mace.r_max}",
        f"--max_L={CFG.mace.max_L}",
        f"--correlation={CFG.mace.correlation}",
        f"--hidden_irreps={CFG.mace.hidden_irreps}",
        f"--num_interactions={CFG.mace.num_interactions}",
        f"--batch_size={CFG.mace.batch_size}",
        f"--max_num_epochs={CFG.mace.epochs}",
        f"--lr={CFG.mace.lr}",
        f"--weight_decay={CFG.mace.weight_decay}",
        f"--energy_weight={CFG.mace.energy_weight}",
        f"--forces_weight={CFG.mace.forces_weight}",
        f"--stress_weight={CFG.mace.stress_weight}",
        f"--device={device}",
    ]
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    train()
```

---

## `project02/parity.py` — validate the trained potential

```python
"""Compute parity-plot statistics on the held-out test set."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from mace.calculators import MACECalculator


def evaluate(model_path: Path, test_xyz: Path
             ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    calc = MACECalculator(model_paths=[str(model_path)], device="cuda")
    frames = read(str(test_xyz), index=":")
    e_dft, e_mlip, f_dft, f_mlip = [], [], [], []
    for f in frames:
        n = len(f)
        e_dft.append(float(f.info["energy"]) / n)
        f_dft.extend(np.asarray(f.arrays["forces"]).ravel().tolist())
        f.calc = calc
        e_mlip.append(float(f.get_potential_energy()) / n)
        f_mlip.extend(np.asarray(f.get_forces()).ravel().tolist())
    return (np.array(e_dft), np.array(e_mlip),
            np.array(f_dft), np.array(f_mlip))


def parity_plots(e_dft: np.ndarray, e_mlip: np.ndarray,
                 f_dft: np.ndarray, f_mlip: np.ndarray,
                 out_dir: Path = Path("parity")) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    e_mae = 1000.0 * float(np.mean(np.abs(e_dft - e_mlip)))   # meV/atom
    f_mae = 1000.0 * float(np.mean(np.abs(f_dft - f_mlip)))   # meV/Å
    print(f"E MAE: {e_mae:.2f} meV/atom")
    print(f"F MAE: {f_mae:.2f} meV/Å")

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].scatter(e_dft, e_mlip, s=8, alpha=0.4)
    lo, hi = min(e_dft.min(), e_mlip.min()), max(e_dft.max(), e_mlip.max())
    ax[0].plot([lo, hi], [lo, hi], "k--")
    ax[0].set_xlabel("DFT energy (eV/atom)")
    ax[0].set_ylabel("MLIP energy (eV/atom)")
    ax[0].set_title(f"E MAE = {e_mae:.1f} meV/atom")

    ax[1].scatter(f_dft, f_mlip, s=4, alpha=0.2)
    lo, hi = min(f_dft.min(), f_mlip.min()), max(f_dft.max(), f_mlip.max())
    ax[1].plot([lo, hi], [lo, hi], "k--")
    ax[1].set_xlabel("DFT force (eV/Å)")
    ax[1].set_ylabel("MLIP force (eV/Å)")
    ax[1].set_title(f"F MAE = {f_mae:.1f} meV/Å")
    fig.tight_layout()
    fig.savefig(out_dir / "parity.png", dpi=200)


if __name__ == "__main__":
    e_d, e_m, f_d, f_m = evaluate(Path("cu_mace_v1.model"),
                                  Path("data/cu_test.xyz"))
    parity_plots(e_d, e_m, f_d, f_m)
```

---

## `project02/coex_md.py` — build and run the coexistence cell

```python
"""Drive coexistence MD with MACE via ASE.

Uses ASE's MD machinery rather than LAMMPS for ease of use; for
production-scale runs, switch to LAMMPS with pair_style mace."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.npt import NPT
from ase import units
from mace.calculators import MACECalculator

from .build_cells import coexistence_initial, split_solid_liquid_mask
from .config import CFG


def melt_half(atoms: Atoms, t_hot_K: float = 2500.0,
              ps: float = 5.0, dt_fs: float = 1.0) -> Atoms:
    """Heat the upper half to t_hot_K while freezing the lower half."""
    mask = split_solid_liquid_mask(atoms)
    from ase.constraints import FixAtoms
    atoms_local = atoms.copy()
    atoms_local.set_constraint(FixAtoms(indices=np.where(~mask)[0]))
    MaxwellBoltzmannDistribution(atoms_local, temperature_K=t_hot_K)
    dyn = Langevin(atoms_local, dt_fs * units.fs,
                   temperature_K=t_hot_K, friction=0.01)
    dyn.run(int(ps * 1000))
    atoms_local.set_constraint()
    return atoms_local


def run_coex(model_path: Path,
             temperature_K: float,
             ps_total: float = 200.0,
             out_traj: Path = Path("coex.traj")) -> None:
    init = coexistence_initial()
    init.calc = MACECalculator(model_paths=[str(model_path)], device="cuda")
    melted = melt_half(init)
    # Switch to NPT for production.
    MaxwellBoltzmannDistribution(melted, temperature_K=temperature_K)
    npt = NPT(melted, timestep=CFG.md.timestep_fs * units.fs,
              temperature_K=temperature_K,
              externalstress=CFG.md.pressure_bar * units.bar,
              ttime=CFG.md.thermostat_tau_ps * 1000 * units.fs,
              pfactor=(CFG.md.barostat_tau_ps * 1000 * units.fs) ** 2
                      * units.GPa)
    from ase.io.trajectory import Trajectory
    traj = Trajectory(str(out_traj), "w", melted)
    npt.attach(traj.write, interval=100)  # every 100 fs
    npt.run(int(ps_total * 1000))


if __name__ == "__main__":
    for T in (1300.0, 1350.0, 1400.0, 1450.0):
        run_coex(Path("cu_mace_v1.model"), T,
                 out_traj=Path(f"coex_T{int(T)}.traj"))
```

---

## `project02/interface.py` — locate the solid–liquid interface

```python
"""Detect the interface position via a z-resolved Steinhardt q6 profile."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read


def q6_profile(atoms, n_bins: int = 40) -> tuple[np.ndarray, np.ndarray]:
    """Approximate per-atom order parameter via coordination-based proxy.

    A full Steinhardt analysis is more accurate; for an undergraduate
    project this proxy (count of neighbours within 3.5 Å) suffices to
    distinguish solid (12 NN) from liquid (~ 10 NN with noise)."""
    pos = atoms.positions
    cell = np.asarray(atoms.cell)
    inv_cell = np.linalg.inv(cell)
    n_atoms = len(atoms)
    coord = np.zeros(n_atoms)
    cutoff = 3.5
    for i in range(n_atoms):
        diffs = pos - pos[i]
        f = diffs @ inv_cell
        f -= np.round(f)
        d = np.linalg.norm(f @ cell, axis=1)
        coord[i] = np.sum((d > 0.1) & (d < cutoff))
    z_edges = np.linspace(0, cell[2, 2], n_bins + 1)
    z_centres = 0.5 * (z_edges[:-1] + z_edges[1:])
    profile = np.zeros(n_bins)
    z = pos[:, 2] % cell[2, 2]
    for b in range(n_bins):
        in_bin = (z >= z_edges[b]) & (z < z_edges[b + 1])
        profile[b] = coord[in_bin].mean() if in_bin.any() else np.nan
    return z_centres, profile


def interface_position(z: np.ndarray, profile: np.ndarray) -> float:
    """Return the z at which the profile crosses the midpoint."""
    valid = ~np.isnan(profile)
    z_v, p_v = z[valid], profile[valid]
    midpoint = 0.5 * (np.nanmin(profile) + np.nanmax(profile))
    above = p_v > midpoint
    # First downward crossing index.
    for i in range(len(p_v) - 1):
        if above[i] and not above[i + 1]:
            frac = (p_v[i] - midpoint) / (p_v[i] - p_v[i + 1] + 1e-9)
            return float(z_v[i] + frac * (z_v[i + 1] - z_v[i]))
    return float(np.nan)


def interface_velocity(traj_path: Path, equil_ps: float = 30.0,
                       sample_every_ps: float = 1.0) -> tuple[float, float]:
    frames = read(str(traj_path), index=":")
    # We assume one frame per 100 fs from coex_md.py.
    times_ps = np.arange(len(frames)) * 0.1
    positions = []
    for f in frames:
        z, p = q6_profile(f)
        positions.append(interface_position(z, p))
    positions = np.array(positions)
    keep = times_ps >= equil_ps
    t = times_ps[keep]
    z = positions[keep]
    valid = ~np.isnan(z)
    slope, intercept = np.polyfit(t[valid], z[valid], 1)
    return float(slope), float(intercept)


def melting_point_from_velocities(temps: list[float],
                                  velocities: list[float]) -> float:
    """Linear interpolation through the v(T) = 0 crossing."""
    slope, intercept = np.polyfit(temps, velocities, 1)
    return float(-intercept / slope)


if __name__ == "__main__":
    temps = [1300.0, 1350.0, 1400.0, 1450.0]
    vs = []
    for T in temps:
        v, _ = interface_velocity(Path(f"coex_T{int(T)}.traj"))
        vs.append(v)
        print(f"T = {T:.0f} K  v = {v:.3f} Å/ps")
    tm = melting_point_from_velocities(temps, vs)
    print(f"T_m = {tm:.0f} K")
    fig, ax = plt.subplots()
    ax.plot(temps, vs, "o-")
    ax.axhline(0, color="k", ls=":")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("Interface velocity (Å/ps)")
    fig.savefig("analysis/interface_velocity.png", dpi=200)
```

---

## Putting it together

```bash
# Phase 1
python -m project02.aimd_qe          # custom driver around `run_one`
python -m project02.build_dataset
# Phase 2
python -m project02.train_mace
python -m project02.parity
# Phase 3 + 4
python -m project02.coex_md
python -m project02.interface
```

The whole pipeline is several days of wall time, dominated by AIMD
and by the four coexistence runs. Plan accordingly.

---

## Notes

- The `q6_profile` here is a coordination-number proxy, not a true
  Steinhardt $q_6$. For publication-grade analysis switch to
  `pyscal` or to LAMMPS's `compute orientorder/atom`.
- `MACECalculator`'s API changes between mace-torch releases; check
  the keyword names against your installed version.
- All MD here uses ASE for portability. For large coexistence cells
  (≥ 8000 atoms), switch to LAMMPS with `pair_style mace` for a
  ≈ 10× speedup.
- `dt_fs = 1.0` is appropriate for Cu; do not exceed 2 fs.
