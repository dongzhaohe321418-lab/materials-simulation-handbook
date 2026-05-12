# Starter code — Defect formation energy in silicon

The scripts below are working Python 3.11 with type hints. They are
deliberately written for clarity, not for maximum compactness. Treat
them as a scaffold: read them, run them, then modify them. The
"runnable" status assumes you have `ase`, `numpy`, `matplotlib`, and a
working `pw.x` on your `PATH`, plus a valid `ESPRESSO_PSEUDO`
environment variable.

All files below live under `project01/`.

---

## `project01/config.py` — central configuration

```python
"""Single source of truth for parameters across the project."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class QEConfig:
    """Quantum ESPRESSO parameters."""
    pseudo_dir: Path = Path.home() / "pseudo" / "sssp_efficiency_PBE"
    pseudo_si: str = "Si.pbe-n-rrkjus_psl.1.0.0.UPF"
    ecutwfc: float = 40.0       # Ry
    ecutrho: float = 320.0      # Ry
    conv_thr: float = 1.0e-9
    mixing_beta: float = 0.5
    degauss: float = 0.005      # Ry, Gaussian smearing
    forc_conv_thr: float = 1.0e-4   # Ry/bohr
    etot_conv_thr: float = 1.0e-5   # Ry


@dataclass(frozen=True)
class ProjectConfig:
    """Top-level paths and physical constants."""
    a0_si: float = 5.469            # Å, PBE-equilibrium lattice parameter
    work_dir: Path = Path("./runs").resolve()
    qe: QEConfig = QEConfig()


RY_TO_EV: float = 13.605693122994
BOHR_TO_ANG: float = 0.529177210903


CFG = ProjectConfig()
```

---

## `project01/build_cells.py` — build perfect and vacancy supercells

```python
"""Construct perfect and vacancy supercells of silicon."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import bulk
from ase.io import write

from .config import CFG


def perfect_si_supercell(n: int) -> Atoms:
    """Return an n*n*n repeat of the 8-atom conventional Si cell."""
    si_conv = bulk("Si", crystalstructure="diamond", a=CFG.a0_si, cubic=True)
    return si_conv.repeat((n, n, n))


def make_vacancy(perfect: Atoms,
                 index: int = 0,
                 perturb: float = 0.05,
                 rng_seed: int = 42) -> Atoms:
    """Remove an atom and rattle its four nearest neighbours.

    The rattle breaks T_d to allow the Jahn-Teller distortion to
    relax into a D_{2d} state. Returns a new Atoms with len(perfect) - 1
    atoms."""
    rng = np.random.default_rng(rng_seed)
    centre = perfect.positions[index].copy()
    vac = perfect.copy()
    del vac[index]
    # Compute distances to the original position (centre).
    diffs = vac.positions - centre
    # Apply minimum-image convention for periodic NN search.
    cell = np.asarray(vac.cell)
    inv_cell = np.linalg.inv(cell)
    fractional = diffs @ inv_cell
    fractional -= np.round(fractional)
    diffs_mic = fractional @ cell
    d = np.linalg.norm(diffs_mic, axis=1)
    nn = np.argsort(d)[:4]
    vac.positions[nn] += perturb * rng.standard_normal(size=(4, 3))
    return vac


def build_all(sizes: tuple[int, ...] = (2, 3, 4)) -> dict[int, dict[str, Atoms]]:
    """Build perfect and vacancy cells for each n in sizes (n*n*n conv. cells).

    Returns a nested dict keyed by atom count, with 'perfect' and 'vacancy'."""
    out: dict[int, dict[str, Atoms]] = {}
    for n in sizes:
        perfect = perfect_si_supercell(n)
        vacancy = make_vacancy(perfect)
        out[len(perfect)] = {"perfect": perfect, "vacancy": vacancy}
    return out


if __name__ == "__main__":
    outdir = Path("structures")
    outdir.mkdir(exist_ok=True)
    cells = build_all()
    for natoms, pair in cells.items():
        write(outdir / f"Si_N{natoms:03d}_perfect.xyz", pair["perfect"])
        write(outdir / f"Si_N{natoms - 1:03d}_vacancy.xyz", pair["vacancy"])
        print(f"N={natoms}: perfect {len(pair['perfect'])} atoms, "
              f"vacancy {len(pair['vacancy'])} atoms")
```

---

## `project01/qe_runner.py` — drive `pw.x` via ASE

```python
"""Wrap ASE's Espresso calculator for our defect workflow."""
from __future__ import annotations

from pathlib import Path
from typing import Literal

from ase import Atoms
from ase.calculators.espresso import Espresso

from .config import CFG, RY_TO_EV


def _kpt_for(n: int) -> tuple[int, int, int]:
    """Return the equivalent k-mesh for an n*n*n conventional supercell."""
    if n == 1:
        return (8, 8, 8)
    if n == 2:
        return (4, 4, 4)
    if n == 3:
        return (3, 3, 3)
    if n == 4:
        return (2, 2, 2)
    # General fallback: ceil(8/n) rounded to nearest >=1
    k = max(1, round(8 / n))
    return (k, k, k)


def build_calc(atoms: Atoms,
               calculation: Literal["scf", "relax"] = "scf",
               n_super: int = 2,
               workdir: Path | None = None,
               prefix: str = "si") -> Espresso:
    """Return an Espresso calculator with our standard parameters.

    Parameters
    ----------
    atoms : the system to compute. Used only for species detection.
    calculation : 'scf' or 'relax'.
    n_super : the supercell index n (n*n*n conventional cells), used to
        pick the equivalent k-mesh.
    workdir : where QE writes its outdir tree.
    prefix : QE 'prefix' string; must be unique per calculation.
    """
    workdir = workdir or (CFG.work_dir / prefix)
    workdir.mkdir(parents=True, exist_ok=True)
    input_data: dict[str, dict[str, object]] = {
        "control": {
            "calculation": calculation,
            "prefix": prefix,
            "pseudo_dir": str(CFG.qe.pseudo_dir),
            "outdir": str(workdir),
            "verbosity": "low",
            "etot_conv_thr": CFG.qe.etot_conv_thr,
            "forc_conv_thr": CFG.qe.forc_conv_thr,
        },
        "system": {
            "ecutwfc": CFG.qe.ecutwfc,
            "ecutrho": CFG.qe.ecutrho,
            "occupations": "smearing",
            "smearing": "gaussian",
            "degauss": CFG.qe.degauss,
            "nosym": True,
        },
        "electrons": {
            "conv_thr": CFG.qe.conv_thr,
            "mixing_beta": CFG.qe.mixing_beta,
        },
    }
    if calculation == "relax":
        input_data["ions"] = {"ion_dynamics": "bfgs"}
    pseudos = {"Si": CFG.qe.pseudo_si}
    return Espresso(
        input_data=input_data,
        pseudopotentials=pseudos,
        kpts=_kpt_for(n_super),
        koffset=(0, 0, 0),
        directory=str(workdir),
    )


def run_and_get_energy_ev(atoms: Atoms,
                          calculation: Literal["scf", "relax"],
                          n_super: int,
                          prefix: str) -> float:
    """Run pw.x and return the total energy in eV."""
    atoms = atoms.copy()
    atoms.calc = build_calc(atoms, calculation=calculation,
                            n_super=n_super, prefix=prefix)
    e_ev: float = atoms.get_potential_energy()  # ASE returns eV already
    # ASE returns eV; we keep this docstring explicit for clarity.
    _ = RY_TO_EV  # silence the import-unused linter
    return e_ev
```

---

## `project01/run_convergence.py` — cut-off and k-mesh convergence

```python
"""Sweep ecutwfc and k-mesh on the 2-atom primitive Si cell."""
from __future__ import annotations

import csv
from pathlib import Path

from ase.build import bulk
from ase.calculators.espresso import Espresso

from .config import CFG


def _calc(ecut: float, kpts: tuple[int, int, int], prefix: str) -> Espresso:
    workdir = CFG.work_dir / "conv" / prefix
    workdir.mkdir(parents=True, exist_ok=True)
    return Espresso(
        input_data={
            "control": {"calculation": "scf", "prefix": prefix,
                        "pseudo_dir": str(CFG.qe.pseudo_dir),
                        "outdir": str(workdir), "verbosity": "low"},
            "system": {"ecutwfc": ecut, "ecutrho": 8 * ecut,
                       "occupations": "fixed"},
            "electrons": {"conv_thr": 1.0e-9, "mixing_beta": 0.5},
        },
        pseudopotentials={"Si": CFG.qe.pseudo_si},
        kpts=kpts,
        directory=str(workdir),
    )


def sweep_ecut(out_csv: Path = Path("convergence/ecut.csv")) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    si = bulk("Si", crystalstructure="diamond", a=CFG.a0_si)
    rows: list[tuple[float, float]] = []
    for ecut in [30.0, 35.0, 40.0, 45.0, 50.0, 60.0]:
        si.calc = _calc(ecut, (8, 8, 8), prefix=f"ecut_{int(ecut)}")
        e = si.get_potential_energy() / len(si)
        rows.append((ecut, e))
        print(f"ecut={ecut} Ry  E/atom={e:.6f} eV")
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ecutwfc_Ry", "E_per_atom_eV"])
        w.writerows(rows)


def sweep_kpts(out_csv: Path = Path("convergence/kpts.csv")) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    si = bulk("Si", crystalstructure="diamond", a=CFG.a0_si)
    rows: list[tuple[int, float]] = []
    for k in [4, 6, 8, 10, 12]:
        si.calc = _calc(CFG.qe.ecutwfc, (k, k, k), prefix=f"k_{k}")
        e = si.get_potential_energy() / len(si)
        rows.append((k, e))
        print(f"k={k}  E/atom={e:.6f} eV")
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["k_per_axis", "E_per_atom_eV"])
        w.writerows(rows)


if __name__ == "__main__":
    sweep_ecut()
    sweep_kpts()
```

---

## `project01/run_defect.py` — orchestrate the defect calculations

```python
"""Run perfect + vacancy for each supercell size and record E_f."""
from __future__ import annotations

import json
from pathlib import Path

from .build_cells import build_all
from .config import CFG
from .qe_runner import run_and_get_energy_ev


def formation_energy_ev(e_defect: float, e_bulk: float, n_bulk: int) -> float:
    """E_f = E_def - (N-1)/N * E_bulk for an elemental crystal."""
    mu = e_bulk / n_bulk
    return e_defect - (n_bulk - 1) * mu


def run_size(n: int) -> dict[str, float | int]:
    """Run the perfect and vacancy calculations for an n*n*n conv. supercell."""
    cells = build_all(sizes=(n,))
    natoms = next(iter(cells))
    perfect = cells[natoms]["perfect"]
    vacancy = cells[natoms]["vacancy"]

    e_bulk = run_and_get_energy_ev(perfect, "scf", n_super=n,
                                   prefix=f"perfect_N{natoms}")
    e_def = run_and_get_energy_ev(vacancy, "relax", n_super=n,
                                  prefix=f"vac_N{natoms - 1}")
    ef = formation_energy_ev(e_def, e_bulk, natoms)
    return {"n_super": n, "n_bulk_atoms": natoms,
            "E_bulk_eV": e_bulk, "E_defect_eV": e_def, "E_f_eV": ef}


def main() -> None:
    summaries: list[dict[str, float | int]] = []
    for n in (2, 3, 4):  # 64, 216, 512 atoms
        try:
            summary = run_size(n)
        except Exception as err:  # noqa: BLE001 — record and continue
            print(f"n={n} failed: {err!r}")
            continue
        summaries.append(summary)
        with (CFG.work_dir / f"summary_n{n}.json").open("w") as fh:
            json.dump(summary, fh, indent=2)
        print(summary)
    with Path("analysis/summary.json").open("w") as fh:
        json.dump(summaries, fh, indent=2)


if __name__ == "__main__":
    Path("analysis").mkdir(exist_ok=True)
    main()
```

---

## `project01/analyse.py` — finite-size extrapolation and plotting

```python
"""Fit E_f vs 1/N and plot the scaling."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def fit_linear_inverse_n(n: np.ndarray, ef: np.ndarray
                         ) -> tuple[float, float, float, float]:
    """Weighted least squares: E_f(N) = E_inf + alpha / N.

    Weights are proportional to N (larger cells trusted more).
    Returns (E_inf, alpha, se_E_inf, se_alpha).
    """
    x = 1.0 / n
    w = n.astype(float)
    W = np.diag(w)
    X = np.column_stack([np.ones_like(x), x])
    XtWX = X.T @ W @ X
    XtWy = X.T @ W @ ef
    beta = np.linalg.solve(XtWX, XtWy)
    residuals = ef - X @ beta
    dof = max(1, len(n) - 2)
    sigma2 = float(residuals.T @ W @ residuals / dof)
    cov = sigma2 * np.linalg.inv(XtWX)
    return float(beta[0]), float(beta[1]), float(np.sqrt(cov[0, 0])), float(np.sqrt(cov[1, 1]))


def main(summary_path: Path = Path("analysis/summary.json"),
         out_png: Path = Path("analysis/scaling.png")) -> None:
    data = json.loads(summary_path.read_text())
    n = np.array([row["n_bulk_atoms"] for row in data])
    ef = np.array([row["E_f_eV"] for row in data])

    e_inf, alpha, se_inf, _ = fit_linear_inverse_n(n, ef)
    print(f"E_f(N -> inf) = {e_inf:.3f} +/- {se_inf:.3f} eV")
    print(f"slope alpha   = {alpha:.3f} eV*atoms")

    x_fine = np.linspace(0.0, max(1.0 / n), 200)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(1.0 / n, ef, "o", label="DFT")
    ax.plot(x_fine, e_inf + alpha * x_fine, "-",
            label=f"fit: $E_f^\\infty = {e_inf:.2f}$ eV")
    ax.axhline(e_inf, ls=":", color="grey")
    ax.set_xlabel(r"$1 / N$")
    ax.set_ylabel(r"$E_f^{V_\mathrm{Si}^0}$ (eV)")
    ax.legend()
    ax.set_xlim(left=0)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
```

---

## `project01/check_symmetry.py` — verify Jahn–Teller

```python
"""Check that the relaxed vacancy is D_{2d} (not T_d)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read


def four_nn_distances(relaxed: Atoms,
                      vacancy_position: np.ndarray) -> np.ndarray:
    """Return the four shortest distances from the vacancy site."""
    diffs = relaxed.positions - vacancy_position
    cell = np.asarray(relaxed.cell)
    inv_cell = np.linalg.inv(cell)
    f = diffs @ inv_cell
    f -= np.round(f)
    d = np.linalg.norm(f @ cell, axis=1)
    return np.sort(d)[:4]


def report(traj_path: Path, vacancy_position: np.ndarray) -> None:
    atoms = read(traj_path)
    d = four_nn_distances(atoms, vacancy_position)
    print(f"NN distances (Å): {d}")
    spread = float(d.max() - d.min())
    if spread < 0.02:
        print("Symmetry appears to be T_d — JT distortion did NOT occur.")
    else:
        print(f"D_{2d}-like distortion: spread {spread*1000:.0f} mÅ")


if __name__ == "__main__":
    import sys
    traj = Path(sys.argv[1])
    centre = np.array(list(map(float, sys.argv[2:5])))
    report(traj, centre)
```

---

## Running it end-to-end

```bash
cd project01
python -m project01.run_convergence
python -m project01.run_defect
python -m project01.analyse
python -m project01.check_symmetry runs/vac_N063/espresso.pwo  X Y Z
```

The last command takes the original Cartesian coordinates of the
removed atom; print them from `build_cells.py` if you have not kept
track. The whole pipeline produces `convergence/*.csv`,
`runs/summary_n*.json`, `analysis/summary.json`, and
`analysis/scaling.png`.

---

## Style notes

- Every function is fully type-annotated.
- All paths are `pathlib.Path` objects, not raw strings.
- No global state beyond the immutable `CFG` dataclass.
- Each script is a `python -m project01.*` entry point — you must
  add `__init__.py` to the `project01/` directory.

If the QE runs are too slow on your machine, drop the `n=4` (512-atom)
size from `run_defect.py`'s loop; the two-point fit is degraded but
the pipeline still works.
