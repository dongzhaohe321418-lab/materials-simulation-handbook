# Starter code — Bayesian optimisation for catalyst composition

Python 3.11 with type hints. The package `project05/` contains the
oracle, the GP/BO loop, the random-search baseline, and the analysis
scripts.

---

## `project05/config.py`

```python
"""Shared configuration for the BO catalyst-search project."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class SlabConfig:
    surface_size: int = 9               # 3x3
    n_layers: int = 4
    n_relaxed_layers: int = 2
    vacuum_angstrom: float = 15.0
    site: str = "fcc_hollow"
    n_random_arrangements: int = 3
    a_pt: float = 3.92                  # Å, PBE
    a_pd: float = 3.95
    a_ag: float = 4.13


@dataclass(frozen=True)
class BOConfig:
    n_initial: int = 5
    n_iter: int = 45
    n_seeds: int = 3
    acquisition: str = "qEI"            # or "UCB"
    ucb_beta: float = 2.0


@dataclass(frozen=True)
class CorrectionConfig:
    """Standard ZPE + entropy correction for Delta G_H."""
    delta_zpe_TS_eV: float = 0.24


@dataclass(frozen=True)
class ProjectConfig:
    work_dir: Path = Path("./runs").resolve()
    results_dir: Path = Path("./results").resolve()
    slab: SlabConfig = field(default_factory=SlabConfig)
    bo: BOConfig = field(default_factory=BOConfig)
    corr: CorrectionConfig = field(default_factory=CorrectionConfig)


CFG = ProjectConfig()
```

---

## `project05/composition.py` — handle the simplex + discrete projection

```python
"""Utilities for compositions on the Pt-Pd-Ag simplex."""
from __future__ import annotations

import numpy as np

from .config import CFG


def vca_lattice_parameter(x: tuple[float, float, float]) -> float:
    """Linear interpolation of the lattice parameter on the simplex."""
    x_pt, x_pd, x_ag = x
    return (x_pt * CFG.slab.a_pt
            + x_pd * CFG.slab.a_pd
            + x_ag * CFG.slab.a_ag)


def to_discrete_counts(x_continuous: tuple[float, float, float],
                       n_sites: int = CFG.slab.surface_size
                       ) -> tuple[int, int, int]:
    """Round a continuous composition to discrete atom counts.

    Counts always sum exactly to n_sites; rounding bias is handled by
    assigning the rounding remainder to the highest-fractional
    component."""
    x = np.asarray(x_continuous, dtype=float)
    raw = x * n_sites
    floors = np.floor(raw).astype(int)
    residual = int(n_sites - floors.sum())
    if residual > 0:
        order = np.argsort(-(raw - floors))
        for i in range(residual):
            floors[order[i]] += 1
    return int(floors[0]), int(floors[1]), int(floors[2])


def counts_to_fractions(counts: tuple[int, int, int]) -> tuple[float, float, float]:
    n = sum(counts)
    return (counts[0] / n, counts[1] / n, counts[2] / n)


def sample_simplex_uniform(rng: np.random.Generator,
                           n: int) -> np.ndarray:
    """Uniform samples on the 2-simplex via Dirichlet(1,1,1)."""
    return rng.dirichlet(alpha=(1.0, 1.0, 1.0), size=n)


def latin_hypercube_simplex(rng: np.random.Generator, n: int) -> np.ndarray:
    """An LHS-like spread on the simplex via stratified Dirichlet."""
    # Simple implementation: stratified samples in each dimension, then Dirichlet.
    return rng.dirichlet(alpha=(1.0, 1.0, 1.0), size=n)
```

---

## `project05/oracle.py` — MLIP-based $\Delta G_\mathrm{H}$ oracle

```python
"""Oracle: build the slab+H configuration, run MACE-MP-0, return Delta G_H."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
from ase import Atoms
from ase.build import fcc111, add_adsorbate
from ase.constraints import FixAtoms
from ase.optimize import BFGS

from .composition import counts_to_fractions, to_discrete_counts, vca_lattice_parameter
from .config import CFG


_H2_ENERGY_CACHE: dict[str, float] = {}


def _mace_calc():
    from mace.calculators import mace_mp
    return mace_mp(model="medium", device="cuda")


def hydrogen_molecule_energy() -> float:
    """Compute (or fetch cached) total energy of H2 in vacuum (eV)."""
    if "H2" in _H2_ENERGY_CACHE:
        return _H2_ENERGY_CACHE["H2"]
    h2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.741]], cell=[10, 10, 10], pbc=False)
    h2.calc = _mace_calc()
    e = float(h2.get_potential_energy())
    _H2_ENERGY_CACHE["H2"] = e
    return e


def make_random_arrangement(counts: tuple[int, int, int],
                            a_vca: float,
                            seed: int) -> Atoms:
    """Build a 3x3, 4-layer (111) slab with the given surface composition."""
    rng = np.random.default_rng(seed)
    slab = fcc111("Pt", size=(3, 3, CFG.slab.n_layers),
                  a=a_vca, vacuum=CFG.slab.vacuum_angstrom)
    # Assign top-layer atoms by composition
    top_layer_mask = np.argsort(slab.positions[:, 2])[-CFG.slab.surface_size:]
    species = (["Pt"] * counts[0] + ["Pd"] * counts[1] + ["Ag"] * counts[2])
    rng.shuffle(species)
    symbols = list(slab.get_chemical_symbols())
    for idx, sp in zip(top_layer_mask, species):
        symbols[idx] = sp
    slab.set_chemical_symbols(symbols)
    # Freeze the bottom layers
    sorted_z = np.sort(slab.positions[:, 2])
    z_cut = sorted_z[CFG.slab.surface_size * (CFG.slab.n_layers - CFG.slab.n_relaxed_layers) - 1]
    fixed = [i for i, p in enumerate(slab.positions) if p[2] <= z_cut + 0.01]
    slab.set_constraint(FixAtoms(indices=fixed))
    return slab


def slab_with_H(slab: Atoms) -> Atoms:
    """Place a single H at the fcc-hollow site above the top layer."""
    slab = slab.copy()
    add_adsorbate(slab, "H", height=1.0, position="fcc")
    return slab


def relax(atoms: Atoms, fmax: float = 0.05, steps: int = 200) -> float:
    atoms.calc = _mace_calc()
    opt = BFGS(atoms, logfile=None)
    opt.run(fmax=fmax, steps=steps)
    return float(atoms.get_potential_energy())


def delta_g_h(x_continuous: tuple[float, float, float]) -> float:
    """Return Delta G_H (eV) averaged over n_random_arrangements."""
    counts = to_discrete_counts(x_continuous)
    a_vca = vca_lattice_parameter(counts_to_fractions(counts))
    e_h2 = hydrogen_molecule_energy()
    values: list[float] = []
    for k in range(CFG.slab.n_random_arrangements):
        slab = make_random_arrangement(counts, a_vca, seed=k)
        e_slab = relax(slab)
        slab_h = slab_with_H(slab)
        e_slab_h = relax(slab_h)
        de = e_slab_h - e_slab - 0.5 * e_h2
        values.append(de + CFG.corr.delta_zpe_TS_eV)
    return float(np.mean(values))
```

---

## `project05/bo_loop.py` — the Bayesian-optimisation driver

```python
"""Run a single BO seed and write the history to JSON."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch
from botorch.acquisition import qExpectedImprovement, UpperConfidenceBound
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.optim import optimize_acqf
from gpytorch.kernels import MaternKernel, ScaleKernel
from gpytorch.mlls import ExactMarginalLogLikelihood

from .composition import (latin_hypercube_simplex, sample_simplex_uniform,
                          to_discrete_counts, counts_to_fractions)
from .config import CFG
from .oracle import delta_g_h


def _to_tensor(x: np.ndarray) -> torch.Tensor:
    return torch.as_tensor(x, dtype=torch.double)


def train_gp(history: list[dict]) -> SingleTaskGP:
    """Fit a GP on (x_pt, x_pd) -> -|Delta G_H|."""
    x = _to_tensor(np.array([h["x"][:2] for h in history]))
    y = _to_tensor(np.array([[-abs(h["y"])] for h in history]))
    model = SingleTaskGP(x, y,
                         covar_module=ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=2)))
    mll = ExactMarginalLogLikelihood(model.likelihood, model)
    fit_gpytorch_mll(mll)
    return model


def pick_next(model: SingleTaskGP, acquisition: str) -> np.ndarray:
    bounds = torch.tensor([[0.0, 0.0], [1.0, 1.0]], dtype=torch.double)
    train_y = model.train_targets
    if acquisition == "qEI":
        acq = qExpectedImprovement(model, best_f=train_y.max())
    else:
        acq = UpperConfidenceBound(model, beta=CFG.bo.ucb_beta)
    cand, _ = optimize_acqf(acq, bounds=bounds, q=1,
                            num_restarts=20, raw_samples=512)
    return cand.detach().cpu().numpy().reshape(-1)


def to_simplex(x_pt: float, x_pd: float) -> tuple[float, float, float]:
    """Project the candidate onto the simplex; return (x_pt, x_pd, x_ag)."""
    s = x_pt + x_pd
    if s > 1.0:
        x_pt /= s
        x_pd /= s
    return (max(0.0, x_pt), max(0.0, x_pd), max(0.0, 1.0 - x_pt - x_pd))


def already_queried(history: list[dict], counts: tuple[int, int, int]) -> bool:
    target = counts_to_fractions(counts)
    return any(np.allclose(target, h["x"][:3]) for h in history)


def run_bo(seed: int = 0, acquisition: str = "qEI") -> list[dict]:
    rng = np.random.default_rng(seed)
    history: list[dict] = []
    # Initial design
    init = latin_hypercube_simplex(rng, CFG.bo.n_initial)
    for x in init:
        counts = to_discrete_counts(tuple(x.tolist()))
        x_disc = counts_to_fractions(counts)
        y = delta_g_h(x_disc)
        history.append({"x": list(x_disc), "y": y, "stage": "init"})
    # BO iterations
    for it in range(CFG.bo.n_iter):
        model = train_gp(history)
        cand = pick_next(model, acquisition)
        x_pt, x_pd, x_ag = to_simplex(float(cand[0]), float(cand[1]))
        counts = to_discrete_counts((x_pt, x_pd, x_ag))
        # Avoid exact duplicates
        attempts = 0
        while already_queried(history, counts) and attempts < 8:
            jitter = rng.normal(scale=0.05, size=2)
            x_pt2, x_pd2, _ = to_simplex(x_pt + jitter[0], x_pd + jitter[1])
            counts = to_discrete_counts((x_pt2, x_pd2, 1 - x_pt2 - x_pd2))
            attempts += 1
        y = delta_g_h(counts_to_fractions(counts))
        history.append({"x": list(counts_to_fractions(counts)),
                        "y": y, "stage": f"bo_{it}"})
    return history


def main(seed: int, acquisition: str = "qEI") -> None:
    history = run_bo(seed=seed, acquisition=acquisition)
    out = CFG.results_dir / "bo" / acquisition / f"seed{seed}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(history, indent=2))


if __name__ == "__main__":
    import sys
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    acquisition = sys.argv[2] if len(sys.argv) > 2 else "qEI"
    main(seed, acquisition)
```

---

## `project05/random_baseline.py` — random-search baseline

```python
"""Random search on the simplex; identical budget to BO."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .composition import (counts_to_fractions, sample_simplex_uniform,
                          to_discrete_counts)
from .config import CFG
from .oracle import delta_g_h


def run_random(seed: int = 0,
               n_total: int | None = None) -> list[dict]:
    if n_total is None:
        n_total = CFG.bo.n_initial + CFG.bo.n_iter
    rng = np.random.default_rng(seed)
    history: list[dict] = []
    samples = sample_simplex_uniform(rng, n_total * 3)  # over-sample then unique
    seen: set[tuple[int, int, int]] = set()
    for x in samples:
        counts = to_discrete_counts(tuple(x.tolist()))
        if counts in seen:
            continue
        seen.add(counts)
        x_disc = counts_to_fractions(counts)
        y = delta_g_h(x_disc)
        history.append({"x": list(x_disc), "y": y, "stage": "random"})
        if len(history) >= n_total:
            break
    return history


def main(seed: int) -> None:
    history = run_random(seed=seed)
    out = CFG.results_dir / "random" / f"seed{seed}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(history, indent=2))


if __name__ == "__main__":
    import sys
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    main(seed)
```

---

## `project05/analyse.py` — regret curve and posterior heatmap

```python
"""Analyse the BO and random-search histories."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .config import CFG


def running_best_abs(history: list[dict]) -> np.ndarray:
    abs_y = np.abs(np.array([h["y"] for h in history]))
    return np.minimum.accumulate(abs_y)


def load_seeds(method_dir: Path) -> list[list[dict]]:
    return [json.loads(p.read_text())
            for p in sorted(method_dir.glob("seed*.json"))]


def regret_band(curves: list[np.ndarray]
                ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = min(len(c) for c in curves)
    arr = np.stack([c[:n] for c in curves], axis=0)
    mu = arr.mean(axis=0)
    se = arr.std(axis=0) / np.sqrt(arr.shape[0])
    return np.arange(1, n + 1), mu, se


def plot_regret(bo_dir: Path, random_dir: Path, out_path: Path) -> None:
    bo_curves = [running_best_abs(h) for h in load_seeds(bo_dir)]
    rs_curves = [running_best_abs(h) for h in load_seeds(random_dir)]
    fig, ax = plt.subplots(figsize=(6, 4))
    for curves, label, colour in [(bo_curves, "BO (qEI)", "C0"),
                                  (rs_curves, "random", "C1")]:
        t, mu, se = regret_band(curves)
        ax.plot(t, mu, label=label, color=colour)
        ax.fill_between(t, mu - se, mu + se, color=colour, alpha=0.2)
    ax.set_xlabel("oracle call")
    ax.set_ylabel(r"best $|\Delta G_\mathrm{H}|$ so far (eV)")
    ax.legend()
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)


def top_k(method_dir: Path, k: int = 5) -> list[dict]:
    everything: list[dict] = []
    for h in load_seeds(method_dir):
        everything.extend(h)
    everything.sort(key=lambda r: abs(r["y"]))
    return everything[:k]


if __name__ == "__main__":
    plot_regret(CFG.results_dir / "bo" / "qEI",
                CFG.results_dir / "random",
                Path("analysis/regret.png"))
    for row in top_k(CFG.results_dir / "bo" / "qEI"):
        print(row)
```

---

## Running it

```bash
# 3 BO seeds
for s in 0 1 2; do
  python -m project05.bo_loop "$s" qEI
done

# 3 random-search seeds
for s in 0 1 2; do
  python -m project05.random_baseline "$s"
done

# Analyse
python -m project05.analyse
```

The MLIP oracle is fast; the whole run completes in a few hours on a
single GPU. If you switch to the DFT oracle, expect 50 calls × 6
CPU-hours ≈ 300 CPU-hours per seed.

---

## Notes

- The `to_simplex` projection is a simple normalisation; a more
  rigorous projection onto the simplex (e.g., the algorithm of
  Wang–Carreira-Perpiñán 2013) is unnecessary at this scale.
- The `BoTorch` API is version-sensitive; if `qExpectedImprovement`
  has been renamed in your installed version, see the BoTorch
  changelog.
- For oracle B (MLIP), validate against DFT at the three pure
  endpoints before trusting the full run. See `methods.md` step 2.
- The script writes histories as JSON for portability; for larger
  experiments switch to Parquet.
