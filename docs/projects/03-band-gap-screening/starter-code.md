# Starter code — High-throughput band-gap screening

Python 3.11 throughout, with type hints. Each script is a standalone
entry point. Together they cover the five stages of the methods
document.

---

## `project03/config.py`

```python
"""Shared configuration."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class QueryConfig:
    """Materials Project query filters."""
    elements_required: tuple[str, ...] = ("O",)
    n_elements: tuple[int, int] = (2, 4)
    nsites_max: int = 30


@dataclass(frozen=True)
class CleanConfig:
    """Post-query cleaning thresholds."""
    e_hull_max: float = 0.2          # eV/atom
    band_gap_min: float = 0.0


@dataclass(frozen=True)
class CGCNNConfig:
    """CGCNN hyperparameters."""
    atom_fea_len: int = 64
    h_fea_len: int = 128
    n_conv: int = 3
    n_h: int = 1
    radial_cutoff: float = 8.0       # Å
    max_num_nbr: int = 12
    batch_size: int = 256
    epochs: int = 60
    lr: float = 1.0e-2
    n_seeds: int = 5


@dataclass(frozen=True)
class ScreenConfig:
    """Acquisition and filtering for the screen."""
    gap_min: float = 1.3
    gap_max: float = 2.7
    sigma_max: float = 0.3
    e_hull_filter: float = 0.1
    target_centre: float = 2.0
    lambda_uncertainty: float = 1.0


@dataclass(frozen=True)
class ProjectConfig:
    data_dir: Path = Path("data")
    models_dir: Path = Path("models")
    screen_dir: Path = Path("screen")
    shortlist_dir: Path = Path("shortlist")
    query: QueryConfig = field(default_factory=QueryConfig)
    clean: CleanConfig = field(default_factory=CleanConfig)
    cgcnn: CGCNNConfig = field(default_factory=CGCNNConfig)
    screen: ScreenConfig = field(default_factory=ScreenConfig)


CFG = ProjectConfig()
```

---

## `project03/query_mp.py` — fetch and cache the MP oxide subset

```python
"""Query Materials Project for oxides and cache locally."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from mp_api.client import MPRester

from .config import CFG


def fetch_oxides() -> pd.DataFrame:
    """Query MP and return a DataFrame of summary records."""
    with MPRester() as mpr:
        docs = mpr.materials.summary.search(
            elements=list(CFG.query.elements_required),
            num_elements=CFG.query.n_elements,
            nsites_max=CFG.query.nsites_max,
            fields=["material_id", "formula_pretty", "structure",
                    "band_gap", "energy_above_hull", "is_stable",
                    "is_metal", "symmetry"],
        )
    rows: list[dict] = []
    for d in docs:
        rows.append({
            "material_id": d.material_id,
            "formula": d.formula_pretty,
            "band_gap": float(d.band_gap) if d.band_gap is not None else None,
            "e_hull": float(d.energy_above_hull) if d.energy_above_hull is not None else None,
            "is_stable": bool(d.is_stable),
            "is_metal": bool(d.is_metal),
            "spacegroup": d.symmetry.symbol if d.symmetry else None,
            "nsites": len(d.structure),
            "structure_json": d.structure.to_json(),
        })
    return pd.DataFrame(rows)


def main(out_path: Path = CFG.data_dir / "mp_oxides_raw.parquet") -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df = fetch_oxides()
    df.to_parquet(out_path)
    print(f"Wrote {len(df)} entries to {out_path}")


if __name__ == "__main__":
    main()
```

---

## `project03/clean_data.py` — filter and split

```python
"""Clean the MP query and produce train/val/test splits."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .config import CFG


def clean(df: pd.DataFrame) -> pd.DataFrame:
    """Apply the standard filters (see methods.md)."""
    n0 = len(df)
    df = df.dropna(subset=["band_gap", "e_hull"])
    df = df[df["band_gap"] >= CFG.clean.band_gap_min]
    df = df[df["e_hull"] <= CFG.clean.e_hull_max]
    df = df[df["nsites"] <= CFG.query.nsites_max]
    df = df.sort_values("e_hull").drop_duplicates(
        subset=["formula", "spacegroup"], keep="first"
    )
    print(f"Cleaning: {n0} -> {len(df)}")
    return df.reset_index(drop=True)


def split_by_formula(df: pd.DataFrame,
                     fractions: tuple[float, float, float] = (0.8, 0.1, 0.1),
                     seed: int = 0) -> dict[str, list[str]]:
    rng = np.random.default_rng(seed)
    formulas = df["formula"].unique()
    rng.shuffle(formulas)
    n_train = int(fractions[0] * len(formulas))
    n_val = int(fractions[1] * len(formulas))
    train_f = set(formulas[:n_train])
    val_f = set(formulas[n_train:n_train + n_val])
    test_f = set(formulas[n_train + n_val:])
    def pick(s: set[str]) -> list[str]:
        return df[df["formula"].isin(s)]["material_id"].tolist()
    return {"train": pick(train_f), "val": pick(val_f), "test": pick(test_f)}


def main() -> None:
    raw = pd.read_parquet(CFG.data_dir / "mp_oxides_raw.parquet")
    clean_df = clean(raw)
    clean_df.to_parquet(CFG.data_dir / "mp_oxides_clean.parquet")
    splits = split_by_formula(clean_df)
    for name, ids in splits.items():
        with (CFG.data_dir / f"{name}.json").open("w") as fh:
            json.dump(ids, fh)
        print(f"{name}: {len(ids)}")


if __name__ == "__main__":
    main()
```

---

## `project03/prepare_cgcnn.py` — write CIFs and the id_prop file

```python
"""Convert the cleaned DataFrame into CGCNN's on-disk layout."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from pymatgen.core import Structure

from .config import CFG


def write_dataset(name: str, ids: list[str], df: pd.DataFrame,
                  out_dir: Path) -> None:
    cif_dir = out_dir / name
    cif_dir.mkdir(parents=True, exist_ok=True)
    rows: list[tuple[str, float]] = []
    by_id = df.set_index("material_id")
    for mid in ids:
        if mid not in by_id.index:
            continue
        record = by_id.loc[mid]
        s = Structure.from_str(record["structure_json"], fmt="json")
        cif_path = cif_dir / f"{mid}.cif"
        s.to(fmt="cif", filename=str(cif_path))
        rows.append((mid, float(record["band_gap"])))
    pd.DataFrame(rows, columns=["material_id", "band_gap_eV"]).to_csv(
        cif_dir / "id_prop.csv", index=False, header=False
    )
    print(f"{name}: wrote {len(rows)} CIFs to {cif_dir}")


def main() -> None:
    df = pd.read_parquet(CFG.data_dir / "mp_oxides_clean.parquet")
    for name in ("train", "val", "test"):
        ids = json.loads((CFG.data_dir / f"{name}.json").read_text())
        write_dataset(name, ids, df, CFG.data_dir / "cgcnn")


if __name__ == "__main__":
    main()
```

---

## `project03/train_cgcnn.py` — train an ensemble

```python
"""Train an ensemble of CGCNN models with different seeds.

This script assumes you have the original `cgcnn` repo cloned and
installed (i.e., `cgcnn/main.py` is importable as a module). If not,
use the included subprocess-based shim."""
from __future__ import annotations

import subprocess
from pathlib import Path

from .config import CFG


def train_one(seed: int,
              cgcnn_main: Path = Path("cgcnn/main.py"),
              data_root: Path = CFG.data_dir / "cgcnn" / "train",
              out_dir: Path = CFG.models_dir) -> Path:
    """Train one CGCNN model and return the checkpoint path."""
    out_dir.mkdir(parents=True, exist_ok=True)
    ckpt = out_dir / f"cgcnn_seed{seed}.pt"
    cmd = [
        "python", str(cgcnn_main), str(data_root),
        "--epochs", str(CFG.cgcnn.epochs),
        "--batch-size", str(CFG.cgcnn.batch_size),
        "--lr", str(CFG.cgcnn.lr),
        "--atom-fea-len", str(CFG.cgcnn.atom_fea_len),
        "--h-fea-len", str(CFG.cgcnn.h_fea_len),
        "--n-conv", str(CFG.cgcnn.n_conv),
        "--n-h", str(CFG.cgcnn.n_h),
        "--workers", "4",
        "--seed", str(seed),
        "--checkpoint", str(ckpt),
    ]
    subprocess.run(cmd, check=True)
    return ckpt


def main() -> None:
    for seed in range(CFG.cgcnn.n_seeds):
        ckpt = train_one(seed)
        print(f"Seed {seed}: {ckpt}")


if __name__ == "__main__":
    main()
```

---

## `project03/predict_ensemble.py` — predict + uncertainty

```python
"""Apply the ensemble to a set of CIF files and record (mu, sigma)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import torch

from .config import CFG


def load_model(ckpt: Path) -> torch.nn.Module:
    """Load a CGCNN checkpoint."""
    state = torch.load(str(ckpt), map_location="cpu")
    from cgcnn.model import CrystalGraphConvNet
    model = CrystalGraphConvNet(
        orig_atom_fea_len=state["orig_atom_fea_len"],
        nbr_fea_len=state["nbr_fea_len"],
        atom_fea_len=CFG.cgcnn.atom_fea_len,
        h_fea_len=CFG.cgcnn.h_fea_len,
        n_conv=CFG.cgcnn.n_conv,
        n_h=CFG.cgcnn.n_h,
    )
    model.load_state_dict(state["state_dict"])
    model.eval()
    return model


def predict_dir(cif_dir: Path,
                checkpoints: list[Path]) -> pd.DataFrame:
    """Return a DataFrame with material_id, mu, sigma."""
    from cgcnn.data import CIFData, collate_pool
    from torch.utils.data import DataLoader

    data = CIFData(str(cif_dir))
    loader = DataLoader(data, batch_size=CFG.cgcnn.batch_size,
                        shuffle=False, collate_fn=collate_pool)

    all_preds: list[np.ndarray] = []
    ids: list[str] = []
    for ckpt in checkpoints:
        model = load_model(ckpt)
        preds: list[float] = []
        ids_local: list[str] = []
        with torch.no_grad():
            for inp, _, cif_ids in loader:
                out = model(*inp)
                preds.extend(out.flatten().tolist())
                ids_local.extend(cif_ids)
        if not ids:
            ids = ids_local
        all_preds.append(np.asarray(preds))
    arr = np.stack(all_preds, axis=0)
    mu = arr.mean(axis=0)
    sigma = arr.std(axis=0)
    return pd.DataFrame({"material_id": ids, "mu_eV": mu, "sigma_eV": sigma})


def main() -> None:
    checkpoints = sorted(CFG.models_dir.glob("cgcnn_seed*.pt"))
    test = predict_dir(CFG.data_dir / "cgcnn" / "test", checkpoints)
    test.to_parquet(CFG.data_dir / "test_predictions.parquet")
    screen = predict_dir(CFG.data_dir / "cgcnn" / "screen", checkpoints)
    screen.to_parquet(CFG.screen_dir / "screen_predictions.parquet")


if __name__ == "__main__":
    main()
```

---

## `project03/evaluate.py` — stratified MAE and parity

```python
"""Compute overall and stratified MAE, plus the parity plot."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import CFG


def stratified_mae(df: pd.DataFrame,
                   bins: list[float] | None = None) -> pd.DataFrame:
    """MAE broken down by true-gap bin."""
    if bins is None:
        bins = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, np.inf]
    df = df.copy()
    df["bin"] = pd.cut(df["band_gap"], bins=bins)
    summary = df.groupby("bin", observed=False).apply(
        lambda g: pd.Series({
            "n": len(g),
            "mae": float(np.mean(np.abs(g["band_gap"] - g["mu_eV"]))),
        })
    )
    return summary


def parity(df: pd.DataFrame, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(df["band_gap"], df["mu_eV"], s=6, alpha=0.4)
    lo = float(min(df["band_gap"].min(), df["mu_eV"].min()))
    hi = float(max(df["band_gap"].max(), df["mu_eV"].max()))
    ax.plot([lo, hi], [lo, hi], "k--")
    ax.set_xlabel("MP band gap (eV)")
    ax.set_ylabel("Predicted gap (eV)")
    mae = float(np.mean(np.abs(df["band_gap"] - df["mu_eV"])))
    ax.set_title(f"Test MAE = {mae:.3f} eV (n = {len(df)})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)


def main() -> None:
    df_truth = pd.read_parquet(CFG.data_dir / "mp_oxides_clean.parquet")
    df_pred = pd.read_parquet(CFG.data_dir / "test_predictions.parquet")
    merged = df_truth.merge(df_pred, on="material_id")
    print("Overall MAE (eV):",
          float(np.mean(np.abs(merged["band_gap"] - merged["mu_eV"]))))
    table = stratified_mae(merged)
    print(table)
    table.to_csv("analysis/stratified_mae.csv")
    parity(merged, Path("analysis/parity.png"))


if __name__ == "__main__":
    main()
```

---

## `project03/rank.py` — apply screen filters and produce top-10

```python
"""Filter and rank screen-set candidates."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from .config import CFG


def score(mu: float, sigma: float,
          centre: float = CFG.screen.target_centre,
          lam: float = CFG.screen.lambda_uncertainty) -> float:
    return (mu - centre) ** 2 + lam * sigma ** 2


def main() -> None:
    preds = pd.read_parquet(CFG.screen_dir / "screen_predictions.parquet")
    meta = pd.read_parquet(CFG.data_dir / "screen_meta.parquet")
    df = preds.merge(meta, on="material_id")
    df = df[(df["mu_eV"] >= CFG.screen.gap_min) &
            (df["mu_eV"] <= CFG.screen.gap_max) &
            (df["sigma_eV"] <= CFG.screen.sigma_max)]
    if "e_hull" in df.columns:
        df = df[df["e_hull"].fillna(0.0) <= CFG.screen.e_hull_filter]
    df["score"] = df.apply(lambda r: score(r["mu_eV"], r["sigma_eV"]), axis=1)
    df = df.sort_values("score").reset_index(drop=True)
    CFG.shortlist_dir.mkdir(parents=True, exist_ok=True)
    df.head(10).to_csv(CFG.shortlist_dir / "top10.csv", index=False)
    print(df.head(10))


if __name__ == "__main__":
    main()
```

---

## `project03/dft_verify.py` — make QE inputs for top candidates

```python
"""Build QE input directories for the top-K candidates."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from ase.calculators.espresso import Espresso
from ase.io import write
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor

from .config import CFG


def _calc(prefix: str, workdir: Path) -> Espresso:
    workdir.mkdir(parents=True, exist_ok=True)
    return Espresso(
        input_data={
            "control": {
                "calculation": "vc-relax",
                "prefix": prefix,
                "outdir": str(workdir),
                "verbosity": "low",
            },
            "system": {
                "ecutwfc": 60.0,
                "ecutrho": 480.0,
                "occupations": "smearing",
                "smearing": "mv",
                "degauss": 0.005,
            },
            "electrons": {"conv_thr": 1.0e-8, "mixing_beta": 0.3},
            "ions": {"ion_dynamics": "bfgs"},
            "cell": {"cell_dofree": "all"},
        },
        pseudopotentials={},   # populate per-element below
        kpts=(3, 3, 3),
        directory=str(workdir),
    )


def build_inputs(top_csv: Path = CFG.shortlist_dir / "top10.csv",
                 meta: Path = CFG.data_dir / "mp_oxides_clean.parquet",
                 out_dir: Path = Path("dft"),
                 n_verify: int = 3) -> None:
    short = pd.read_csv(top_csv).head(n_verify)
    mmeta = pd.read_parquet(meta).set_index("material_id")
    for _, row in short.iterrows():
        record = mmeta.loc[row["material_id"]]
        structure = Structure.from_str(record["structure_json"], fmt="json")
        atoms = AseAtomsAdaptor().get_atoms(structure)
        sub = out_dir / row["material_id"]
        atoms.calc = _calc(row["material_id"], sub)
        write(sub / "input.pwi", atoms, format="espresso-in",
              input_data=atoms.calc.parameters["input_data"],
              pseudopotentials=atoms.calc.parameters["pseudopotentials"],
              kpts=atoms.calc.parameters["kpts"])
        print(f"Wrote {sub / 'input.pwi'}")


if __name__ == "__main__":
    build_inputs()
```

---

## Running it

```bash
python -m project03.query_mp
python -m project03.clean_data
python -m project03.prepare_cgcnn
python -m project03.train_cgcnn       # 5 seeds, several hours each
python -m project03.predict_ensemble
python -m project03.evaluate
# Construct screen-set CIFs separately (not shown here — depends on choice of source).
python -m project03.rank
python -m project03.dft_verify
# Submit each dft/<id>/input.pwi to your cluster.
```

---

## Notes

- The `cgcnn` import paths follow the original Xie repo; if you fork
  a different one, adjust the imports in `predict_ensemble.py` and
  `train_cgcnn.py`.
- `screen_meta.parquet` is constructed by your own screen-set builder
  (not shown). The minimum columns are `material_id, formula,
  spacegroup, e_hull` (optional), `source`.
- All file paths are `pathlib.Path` objects. Do not silently use
  string concatenation.
- Pin the MP API snapshot date in your `report.pdf`; MP updates and
  IDs change.
