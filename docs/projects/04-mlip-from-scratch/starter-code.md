# Starter code — Training an MLIP from scratch

Python 3.11 with type hints throughout. The package `project04/`
contains the sampling driver, the training launcher, the validation
suite, and the active-learning loop.

The starter code is deliberately *less abstract* than Project 2:
because you may need to inspect and edit the loops yourself during
debugging, the code keeps the abstraction shallow.

---

## `project04/config.py`

```python
"""Shared configuration for the from-scratch MLIP project."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SamplingConfig:
    aimd_temperatures_K: tuple[float, ...] = (100.0, 300.0, 800.0, 1500.0)
    aimd_ps_per_T: float = 2.0
    aimd_sample_every_fs: float = 20.0
    rattle_stdev_ang: tuple[float, ...] = (0.05, 0.10, 0.15)
    strain_fractions: tuple[float, ...] = (-0.05, -0.02, 0.02, 0.05)
    n_rattled_per_strain: int = 3


@dataclass(frozen=True)
class MACEConfig:
    r_max: float = 5.0
    max_L: int = 2
    correlation: int = 3
    hidden_irreps: str = "128x0e+128x1o"
    num_interactions: int = 2
    batch_size: int = 8
    epochs: int = 250
    lr: float = 0.01
    weight_decay: float = 5.0e-7
    energy_weight: float = 1.0
    forces_weight: float = 100.0
    stress_weight: float = 1.0


@dataclass(frozen=True)
class ALConfig:
    md_ps: float = 50.0
    check_every_fs: float = 50.0
    sigma_threshold_eV_per_A: float = 0.1
    n_select_per_round: int = 40


@dataclass(frozen=True)
class ProjectConfig:
    work_dir: Path = Path("./runs").resolve()
    data_dir: Path = Path("./data").resolve()
    models_dir: Path = Path("./models").resolve()
    sampling: SamplingConfig = SamplingConfig()
    mace: MACEConfig = MACEConfig()
    al: ALConfig = ALConfig()


CFG = ProjectConfig()
```

---

## `project04/baseline_mp0.py` — does the foundation model already work?

```python
"""Compute baseline MAE of MACE-MP-0 on a small probe set."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read


def baseline_metrics(probe_xyz: Path,
                     model: str = "medium") -> tuple[float, float]:
    """Return (E MAE meV/atom, F MAE meV/Å) of MACE-MP-0 on probe set."""
    from mace.calculators import mace_mp
    calc = mace_mp(model=model, device="cuda")
    frames = read(str(probe_xyz), index=":")
    e_truth, e_pred, f_truth, f_pred = [], [], [], []
    for f in frames:
        n = len(f)
        e_truth.append(float(f.info["energy"]) / n)
        f_truth.extend(np.asarray(f.arrays["forces"]).ravel().tolist())
        f.calc = calc
        e_pred.append(float(f.get_potential_energy()) / n)
        f_pred.extend(np.asarray(f.get_forces()).ravel().tolist())
    e_mae_meV = 1000.0 * float(np.mean(np.abs(np.asarray(e_truth) - np.asarray(e_pred))))
    f_mae_meV = 1000.0 * float(np.mean(np.abs(np.asarray(f_truth) - np.asarray(f_pred))))
    return e_mae_meV, f_mae_meV


if __name__ == "__main__":
    e, f = baseline_metrics(Path("data/probe.xyz"))
    print(f"MACE-MP-0 baseline: E MAE = {e:.1f} meV/atom, F MAE = {f:.1f} meV/Å")
    if f < 100:
        print("Foundation model appears adequate — consider fine-tuning rather than training from scratch.")
```

---

## `project04/sample.py` — generate the training configurations

```python
"""Build the sampling set: equilibrium + rattled + strained + AIMD frames."""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms
from ase.io import read, write

from .config import CFG


def rattle_set(equilibrium: Atoms,
               stdevs: Iterable[float] = CFG.sampling.rattle_stdev_ang,
               n_per_stdev: int = 5,
               seed: int = 0) -> list[Atoms]:
    rng = np.random.default_rng(seed)
    out: list[Atoms] = []
    for stdev in stdevs:
        for _ in range(n_per_stdev):
            a = equilibrium.copy()
            a.rattle(stdev=stdev, seed=int(rng.integers(0, 2**31)))
            a.info["config_type"] = f"rattle_{stdev:.2f}"
            out.append(a)
    return out


def strain_set(equilibrium: Atoms,
               fractions: Iterable[float] = CFG.sampling.strain_fractions,
               rattle_per: int = CFG.sampling.n_rattled_per_strain,
               seed: int = 1) -> list[Atoms]:
    rng = np.random.default_rng(seed)
    out: list[Atoms] = []
    for frac in fractions:
        scaled = equilibrium.copy()
        scaled.set_cell(equilibrium.cell * (1.0 + frac), scale_atoms=True)
        scaled.info["config_type"] = f"strain_{frac:+.2f}"
        out.append(scaled.copy())
        for _ in range(rattle_per):
            a = scaled.copy()
            a.rattle(stdev=0.05, seed=int(rng.integers(0, 2**31)))
            a.info["config_type"] = f"strain_{frac:+.2f}_rattle"
            out.append(a)
    return out


def aimd_frames(aimd_xyz: Path) -> list[Atoms]:
    """Read AIMD frames produced externally by run_aimd.py."""
    return read(str(aimd_xyz), index=":")


def build_master(equilibrium: Atoms,
                 aimd_paths: list[Path],
                 out_xyz: Path = Path("data/master.xyz")) -> None:
    frames: list[Atoms] = [equilibrium.copy()]
    frames[-1].info["config_type"] = "equilibrium"
    frames += rattle_set(equilibrium)
    frames += strain_set(equilibrium)
    for p in aimd_paths:
        frames += aimd_frames(p)
    out_xyz.parent.mkdir(parents=True, exist_ok=True)
    write(str(out_xyz), frames)
    print(f"Wrote {len(frames)} configurations to {out_xyz}")


if __name__ == "__main__":
    eq = read("structures/equilibrium.xyz")
    build_master(eq, list((CFG.data_dir / "aimd").glob("*.xyz")))
```

---

## `project04/run_aimd.py` — drive an AIMD trajectory with QE

```python
"""Run a single AIMD trajectory; the build script then converts .pwo -> .xyz."""
from __future__ import annotations

from pathlib import Path

from ase import Atoms
from ase.calculators.espresso import Espresso

from .config import CFG


def aimd_calc(atoms: Atoms, prefix: str, temperature_K: float,
              ps: float, workdir: Path) -> Espresso:
    workdir.mkdir(parents=True, exist_ok=True)
    nstep = int(ps * 1000)  # 1 fs steps
    return Espresso(
        input_data={
            "control": {
                "calculation": "md", "prefix": prefix,
                "outdir": str(workdir), "verbosity": "low",
                "nstep": nstep, "dt": 41.34,
            },
            "system": {
                "ecutwfc": 60.0, "ecutrho": 480.0,
                "occupations": "smearing", "smearing": "mv", "degauss": 0.01,
            },
            "electrons": {"conv_thr": 1.0e-7, "mixing_beta": 0.3},
            "ions": {"ion_dynamics": "verlet",
                     "ion_temperature": "rescale-v",
                     "tempw": temperature_K, "nraise": 50},
        },
        pseudopotentials={},   # populate per-element in caller
        kpts=(2, 2, 1),         # adjust to cell
        directory=str(workdir),
    )


def run(atoms: Atoms, prefix: str, T_K: float, ps: float,
        pseudos: dict[str, str], kpts: tuple[int, int, int]) -> None:
    atoms = atoms.copy()
    workdir = CFG.work_dir / "aimd" / prefix
    calc = aimd_calc(atoms, prefix, T_K, ps, workdir)
    calc.parameters["pseudopotentials"] = pseudos
    calc.parameters["kpts"] = kpts
    atoms.calc = calc
    atoms.get_potential_energy()


if __name__ == "__main__":
    from ase.io import read
    eq = read("structures/equilibrium.xyz")
    for T in CFG.sampling.aimd_temperatures_K:
        run(eq, prefix=f"aimd_T{int(T)}", T_K=T,
            ps=CFG.sampling.aimd_ps_per_T,
            pseudos={"Mo": "Mo.UPF", "S": "S.UPF", "Se": "Se.UPF"},
            kpts=(2, 2, 1))
```

---

## `project04/train.py` — wrapper around mace-torch

```python
"""Train MACE on the master training set."""
from __future__ import annotations

import subprocess
from pathlib import Path

from .config import CFG


def train(name: str,
          train_xyz: Path = Path("data/train.xyz"),
          val_xyz: Path = Path("data/val.xyz"),
          E0s_json: str = '{"Mo": -1234.0, "S": -456.0, "Se": -789.0}',
          device: str = "cuda") -> None:
    CFG.models_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mace_run_train",
        f"--name={name}",
        f"--train_file={train_xyz}",
        f"--valid_file={val_xyz}",
        f"--E0s={E0s_json}",
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


def train_committee(version: str = "v1", n: int = 2) -> None:
    for seed in range(n):
        # MACE CLI accepts --seed; if not, set the env var.
        import os
        os.environ["PYTHONHASHSEED"] = str(seed)
        train(name=f"mace_{version}_seed{seed}")


if __name__ == "__main__":
    train_committee("v1", n=2)
```

---

## `project04/validation.py` — the full validation suite

```python
"""Validation suite: parity, MD stability, RDF, phonon spectrum."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms, units
from ase.io import read
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


def load_mace(model_path: Path):
    from mace.calculators import MACECalculator
    return MACECalculator(model_paths=[str(model_path)], device="cuda")


def parity(model_path: Path, test_xyz: Path,
           out_dir: Path = Path("validation/parity")) -> tuple[float, float]:
    out_dir.mkdir(parents=True, exist_ok=True)
    calc = load_mace(model_path)
    frames = read(str(test_xyz), index=":")
    e_d, e_p, f_d, f_p = [], [], [], []
    for f in frames:
        n = len(f)
        e_d.append(float(f.info["energy"]) / n)
        f_d.extend(np.asarray(f.arrays["forces"]).ravel().tolist())
        f.calc = calc
        e_p.append(float(f.get_potential_energy()) / n)
        f_p.extend(np.asarray(f.get_forces()).ravel().tolist())
    e_d, e_p = np.array(e_d), np.array(e_p)
    f_d, f_p = np.array(f_d), np.array(f_p)
    e_mae = 1000.0 * float(np.mean(np.abs(e_d - e_p)))
    f_mae = 1000.0 * float(np.mean(np.abs(f_d - f_p)))

    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].scatter(e_d, e_p, s=8, alpha=0.5)
    lo, hi = min(e_d.min(), e_p.min()), max(e_d.max(), e_p.max())
    ax[0].plot([lo, hi], [lo, hi], "k--")
    ax[0].set_xlabel("DFT (eV/atom)"); ax[0].set_ylabel("MACE (eV/atom)")
    ax[0].set_title(f"E MAE = {e_mae:.2f} meV/atom")
    ax[1].scatter(f_d, f_p, s=3, alpha=0.2)
    lo, hi = min(f_d.min(), f_p.min()), max(f_d.max(), f_p.max())
    ax[1].plot([lo, hi], [lo, hi], "k--")
    ax[1].set_xlabel("DFT (eV/Å)"); ax[1].set_ylabel("MACE (eV/Å)")
    ax[1].set_title(f"F MAE = {f_mae:.1f} meV/Å")
    fig.tight_layout()
    fig.savefig(out_dir / "parity.png", dpi=200)
    return e_mae, f_mae


def md_stability(model_path: Path, equilibrium: Atoms,
                 T_K: float = 300.0, ps: float = 50.0,
                 out_dir: Path = Path("validation/md")) -> dict[str, float]:
    """Run NVT MD; record energy drift and minimum interatomic distance."""
    out_dir.mkdir(parents=True, exist_ok=True)
    atoms = equilibrium.copy()
    atoms.calc = load_mace(model_path)
    MaxwellBoltzmannDistribution(atoms, temperature_K=T_K)
    dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=T_K, friction=0.01)
    energies: list[float] = []
    min_dist: list[float] = []

    def record() -> None:
        energies.append(float(atoms.get_total_energy()))
        pos = atoms.get_positions()
        d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)
        np.fill_diagonal(d, np.inf)
        min_dist.append(float(d.min()))

    dyn.attach(record, interval=50)
    dyn.run(int(ps * 1000))
    e = np.array(energies)
    drift_meV_per_atom = 1000.0 * float((e[-1] - e[0]) / len(atoms))
    fig, ax = plt.subplots(2, 1, figsize=(6, 4), sharex=True)
    t = np.arange(len(energies)) * 0.05  # ps
    ax[0].plot(t, e - e[0])
    ax[0].set_ylabel("E - E(0) (eV)")
    ax[1].plot(t, min_dist)
    ax[1].set_ylabel("min d (Å)")
    ax[1].set_xlabel("t (ps)")
    fig.tight_layout()
    fig.savefig(out_dir / f"md_T{int(T_K)}.png", dpi=200)
    return {"drift_meV_per_atom": drift_meV_per_atom,
            "min_distance_angstrom": float(min(min_dist))}


if __name__ == "__main__":
    eq = read("structures/equilibrium.xyz")
    e_mae, f_mae = parity(Path("models/mace_v1_seed0.model"),
                          Path("data/test.xyz"))
    print(f"Parity: E {e_mae:.2f} meV/atom, F {f_mae:.1f} meV/Å")
    md = md_stability(Path("models/mace_v1_seed0.model"), eq, T_K=300.0)
    print(f"MD: drift {md['drift_meV_per_atom']:.3f} meV/atom, "
          f"min d {md['min_distance_angstrom']:.2f} Å")
```

---

## `project04/active_learning.py` — uncertainty-triggered relabelling

```python
"""Active-learning loop driven by a 2-member MACE committee."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms, units
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from .config import CFG
from .validation import load_mace


def committee_uncertainty(atoms: Atoms,
                          calcs: list) -> float:
    """Max per-atom force standard deviation across committee members."""
    fs = []
    for c in calcs:
        atoms.calc = c
        fs.append(np.asarray(atoms.get_forces()))
    arr = np.stack(fs, axis=0)
    sigma = arr.std(axis=0)            # shape (n_atoms, 3)
    per_atom = np.linalg.norm(sigma, axis=-1) / np.sqrt(3)
    return float(per_atom.max())


def collect_candidates(equilibrium: Atoms,
                       committee_paths: list[Path],
                       T_K: float = 600.0,
                       ps: float = CFG.al.md_ps,
                       sigma_threshold: float = CFG.al.sigma_threshold_eV_per_A,
                       out_xyz: Path = Path("active-learning/round1/candidates.xyz")) -> int:
    out_xyz.parent.mkdir(parents=True, exist_ok=True)
    calcs = [load_mace(p) for p in committee_paths]
    atoms = equilibrium.copy()
    atoms.calc = calcs[0]   # propagate with the first
    MaxwellBoltzmannDistribution(atoms, temperature_K=T_K)
    dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=T_K, friction=0.01)

    flagged: list[Atoms] = []
    cooldown = 0
    check_interval = int(CFG.al.check_every_fs)

    def check() -> None:
        nonlocal cooldown
        if cooldown > 0:
            cooldown -= 1
            return
        sigma = committee_uncertainty(atoms, calcs)
        if sigma > sigma_threshold:
            snapshot = atoms.copy()
            snapshot.info["sigma_eV_per_A"] = sigma
            snapshot.info["config_type"] = f"al_T{int(T_K)}"
            flagged.append(snapshot)
            cooldown = 50   # 50 * 50 fs = 2.5 ps cooldown

    dyn.attach(check, interval=check_interval)
    dyn.run(int(ps * 1000))

    # Diversity-pick by simple sort+stride; replace with clustering for a better job.
    flagged.sort(key=lambda a: -a.info["sigma_eV_per_A"])
    n = min(CFG.al.n_select_per_round, len(flagged))
    selected = flagged[:n]
    write(str(out_xyz), selected)
    return len(selected)


if __name__ == "__main__":
    n = collect_candidates(
        equilibrium=read("structures/equilibrium.xyz"),
        committee_paths=[Path("models/mace_v1_seed0.model"),
                         Path("models/mace_v1_seed1.model")],
    )
    print(f"Flagged {n} candidates for DFT relabelling.")
```

---

## `project04/phonons.py` — phonon-spectrum sanity check

```python
"""Compute MLIP phonon dispersion via phonopy and compare to a DFT reference."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.io import read

from .validation import load_mace


def compute_phonons(equilibrium: Atoms,
                    model_path: Path,
                    supercell: tuple[int, int, int] = (2, 2, 2),
                    displacement_ang: float = 0.01) -> "Phonopy":
    from phonopy import Phonopy
    from phonopy.structure.atoms import PhonopyAtoms

    cell = PhonopyAtoms(symbols=equilibrium.get_chemical_symbols(),
                        positions=equilibrium.get_positions(),
                        cell=equilibrium.get_cell())
    phonon = Phonopy(cell, np.diag(supercell))
    phonon.generate_displacements(distance=displacement_ang)
    supercells = phonon.get_supercells_with_displacements()

    calc = load_mace(model_path)
    force_sets = []
    for sc in supercells:
        ase_atoms = Atoms(symbols=sc.get_chemical_symbols(),
                          positions=sc.get_positions(),
                          cell=sc.get_cell(), pbc=True)
        ase_atoms.calc = calc
        force_sets.append(ase_atoms.get_forces())
    phonon.set_forces(force_sets)
    phonon.produce_force_constants()
    return phonon


def plot_dispersion(phonon, path_qpoints: list[list[float]],
                    out_png: Path) -> None:
    phonon.set_band_structure(path_qpoints)
    bs = phonon.get_band_structure_dict()
    fig, ax = plt.subplots(figsize=(6, 4))
    for freqs in bs["frequencies"]:
        ax.plot(freqs, color="C0", alpha=0.7)
    ax.set_xlabel("k-path index")
    ax.set_ylabel("Frequency (THz)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)


if __name__ == "__main__":
    eq = read("structures/equilibrium.xyz")
    phonon = compute_phonons(eq, Path("models/mace_v1_seed0.model"))
    plot_dispersion(phonon, [[[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0]]],
                    Path("validation/phonon.png"))
```

---

## Running the loop

```bash
# Phase 1
python -m project04.baseline_mp0
# Phase 2
python -m project04.run_aimd      # multiple temperatures
python -m project04.sample
# Phase 3
python -m project04.train         # committee of 2
# Phase 4
python -m project04.validation
# Phase 5
python -m project04.active_learning
# ... DFT relabel selected frames externally, then:
python -m project04.train         # MACE v2
# Repeat as needed.
```

---

## Notes

- The DFT relabelling step is left as a manual prompt: take the
  selected candidates, run single-point SCF with QE, parse the output
  into the same extended-XYZ format as the original training set,
  and concatenate.
- The phonon code assumes you have `phonopy` installed; install via
  `pip install phonopy`.
- All `subprocess` calls in `train.py` are tied to the `mace_run_train`
  CLI as of mace-torch 0.3.x; flag names occasionally change between
  releases — adjust as needed.
- The committee here is two members; if you have compute to spare,
  use four or five for better uncertainty estimates.
