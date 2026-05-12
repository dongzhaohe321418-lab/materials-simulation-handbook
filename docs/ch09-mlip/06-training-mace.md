# 9.6 Training a MACE potential

We now turn from theory to practice. In this section we train a MACE
potential for liquid water using a small reference dataset of $\sim
1000$ DFT configurations, validate it on a held-out test set, and wire
the trained model into ASE to run a short MD trajectory. The
walkthrough uses `mace-torch` (Batatia et al., 2023), which is the
reference implementation. The same pipeline applies, with small
changes to the configuration file, to NequIP, Allegro, and any other
equivariant MLIP.

## 9.6.1 Choosing and preparing the data

### What a training dataset looks like

A training dataset for an MLIP consists of *configurations*, each one
a snapshot of an atomic system labelled by DFT outputs:

- **Atomic positions** $\mathbf{r}_i$ in $\text{\AA}$.
- **Atomic species** $Z_i$.
- **Simulation cell** (the three lattice vectors).
- **Total energy** $E^\mathrm{DFT}$ in eV.
- **Forces** $\mathbf{F}_i^\mathrm{DFT}$ on each atom in $\mathrm{eV}/\text{\AA}$.
- **Stress tensor** $\sigma_{\alpha\beta}^\mathrm{DFT}$ in
  $\mathrm{eV}/\text{\AA}^3$ (optional but recommended for periodic
  systems).

For our walkthrough we use a 1000-configuration water dataset drawn
from a previously published AIMD trajectory at 300 K and 1 bar, with
energies and forces computed at the revPBE0+D3 level (a hybrid
functional with dispersion correction). The dataset spans liquid
water at near-ambient conditions: it does not cover ice phases, the
gas–liquid critical point, or chemistry such as proton transfer in
hydronium. A potential trained on this data will be reliable in the
configurational neighbourhood of room-temperature liquid water and
should not be trusted outside it.

The dataset takes the form of an *extended XYZ* (extxyz) file. Each
frame is a block of the form

```
192
Lattice="12.42 0.0 0.0 0.0 12.42 0.0 0.0 0.0 12.42" Properties=species:S:1:pos:R:3:forces:R:3 energy=-2034.527 stress="..." pbc="T T T"
O 0.123 4.821 6.117  0.012  -0.054  0.119
H 0.876 4.110 6.510  0.341  -0.182  -0.075
...
```

The first line is the atom count; the second is a comment line
encoding the lattice vectors, per-atom property layout, total energy,
stress (Voigt order), and periodic boundary flags; the remaining
lines are the per-atom data. ASE reads and writes this format
natively, and `mace-torch` consumes it directly.

The alternative is HDF5, used when the dataset is too large to fit
comfortably as plain text (above $\sim 10^5$ configurations). For our
1000-frame water example, extxyz is more convenient.

### Train/validation/test split

A standard split is 80 % train, 10 % validation, 10 % test. The
validation set is used during training to monitor overfitting and to
decide when to reduce the learning rate; the test set is held out
entirely until the end and used to report final metrics. With 1000
configurations: 800 training, 100 validation, 100 test.

The split must be *random over configurations*, not stratified by
time order within the AIMD trajectory: consecutive AIMD frames are
highly correlated, and a chronologically late test set would
overestimate generalisation. Set a fixed random seed for
reproducibility.

```python
import ase.io
import numpy as np

frames = ase.io.read("water_revpbe0.xyz", index=":")
rng = np.random.default_rng(seed=0)
idx = rng.permutation(len(frames))
n_train, n_val = 800, 100
ase.io.write("train.xyz", [frames[i] for i in idx[:n_train]])
ase.io.write("val.xyz",   [frames[i] for i in idx[n_train:n_train+n_val]])
ase.io.write("test.xyz",  [frames[i] for i in idx[n_train+n_val:]])
```

## 9.6.2 Installing `mace-torch`

`mace-torch` is on PyPI:

```bash
pip install mace-torch ase
```

Confirm GPU support is operational:

```python
import torch
assert torch.cuda.is_available(), "GPU required for reasonable training time"
print(torch.cuda.get_device_name(0))
```

A modern consumer GPU (RTX 4090, A6000) trains the configuration below
in about three hours. On CPU the same training takes one to two days
and is not recommended.

## 9.6.3 Configuring the model

The minimal MACE configuration for water is:

```python
from dataclasses import dataclass

@dataclass
class MACEConfig:
    # Data
    train_file: str = "train.xyz"
    valid_file: str = "val.xyz"
    test_file:  str = "test.xyz"
    E0s: dict = None        # isolated atom energies, keyed by Z
    # Architecture
    r_max: float = 5.0      # neighbour cutoff in angstrom
    num_layers: int = 2     # message-passing layers
    hidden_irreps: str = "128x0e + 128x1o"   # 128 scalars, 128 vectors
    max_ell: int = 1
    correlation: int = 3    # body-order per layer
    num_radial_basis: int = 8
    radial_mlp: tuple = (64, 64, 64)
    # Training
    batch_size: int = 5
    max_num_epochs: int = 100
    lr: float = 0.01
    energy_weight: float = 1.0
    forces_weight: float = 100.0
    stress_weight: float = 1.0
    weight_decay: float = 5e-7
    # Output
    name: str = "water_mace"
    seed: int = 1
```

The `hidden_irreps` field deserves comment. The string
`"128x0e + 128x1o"` declares the equivariant feature space: 128
channels of even scalars (irrep $0e$) plus 128 channels of odd
vectors (irrep $1o$). For organic chemistry with hydrogen bonding,
including the $\ell = 1$ channels is helpful; pushing to $\ell = 2$
($+ \mathtt{128x2e}$) buys further accuracy at $\sim 30 \%$ extra
inference cost.

The *isolated atom energies* `E0s` are a subtraction that makes the
fitting target well-conditioned. Without them the model has to learn
a large baseline ($\sim -2000\,\mathrm{eV}$ for a 192-atom water box)
in addition to the small ($\sim 1\,\mathrm{eV}$) configurational
variation. With per-element subtractions

$$
E^\mathrm{target} = E^\mathrm{DFT} - \sum_i E_0(Z_i),
$$

the target sits within a few eV of zero. The reference $E_0$ values
come from DFT calculations on a single isolated atom in a large box;
for the revPBE0+D3 functional and our pseudopotentials,
$E_0(\mathrm{H}) \approx -13.6\,\mathrm{eV}$,
$E_0(\mathrm{O}) \approx -432.5\,\mathrm{eV}$.

## 9.6.4 The training script

A complete training script follows. The `mace-torch` package exposes
the training loop as a single function, but we expand it here to make
the flow legible.

```python
"""Train a MACE potential for liquid water."""
from __future__ import annotations
from pathlib import Path
import torch
from torch.utils.data import DataLoader
import ase.io
from mace.data import AtomicData, Collater, KeySpecification
from mace.modules import MACE, WeightedEnergyForcesStressLoss
from mace.tools import AtomicNumberTable, scatter_sum
from mace.tools.scripts_utils import get_dataset_from_xyz


def main(cfg: MACEConfig) -> None:
    torch.manual_seed(cfg.seed)
    device = torch.device("cuda")

    # ---- 1. Load and tokenise data ----
    train, valid, test = get_dataset_from_xyz(
        train_path=cfg.train_file,
        valid_path=cfg.valid_file,
        test_path=cfg.test_file,
        config_type_weights={"Default": 1.0},
        energy_key="energy",
        forces_key="forces",
        stress_key="stress",
    )
    z_table = AtomicNumberTable([1, 8])    # H, O

    collate = Collater()
    train_loader = DataLoader(
        [AtomicData.from_config(c, z_table=z_table, cutoff=cfg.r_max) for c in train],
        batch_size=cfg.batch_size, shuffle=True, collate_fn=collate,
    )
    valid_loader = DataLoader(
        [AtomicData.from_config(c, z_table=z_table, cutoff=cfg.r_max) for c in valid],
        batch_size=cfg.batch_size, shuffle=False, collate_fn=collate,
    )

    # ---- 2. Construct the model ----
    model = MACE(
        r_max=cfg.r_max,
        num_bessel=cfg.num_radial_basis,
        num_polynomial_cutoff=5,
        max_ell=cfg.max_ell,
        num_interactions=cfg.num_layers,
        num_elements=2,
        hidden_irreps=cfg.hidden_irreps,
        MLP_irreps="16x0e",
        atomic_energies=torch.tensor(
            [cfg.E0s[1], cfg.E0s[8]], dtype=torch.float64),
        avg_num_neighbors=24.0,
        atomic_numbers=[1, 8],
        correlation=cfg.correlation,
    ).to(device)

    # ---- 3. Loss and optimiser ----
    loss_fn = WeightedEnergyForcesStressLoss(
        energy_weight=cfg.energy_weight,
        forces_weight=cfg.forces_weight,
        stress_weight=cfg.stress_weight,
    )
    optimiser = torch.optim.AdamW(
        model.parameters(), lr=cfg.lr, weight_decay=cfg.weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimiser, factor=0.8, patience=5)

    # ---- 4. Training loop ----
    best_val = float("inf")
    out_dir = Path(cfg.name); out_dir.mkdir(exist_ok=True)
    for epoch in range(cfg.max_num_epochs):
        model.train()
        for batch in train_loader:
            batch = batch.to(device)
            out = model(batch.to_dict(), training=True)
            loss = loss_fn(pred=out, ref=batch)
            optimiser.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            optimiser.step()

        # Validation
        model.eval()
        val_loss = 0.0
        n_struct = 0
        with torch.enable_grad():   # forces need grad even in eval mode
            for batch in valid_loader:
                batch = batch.to(device)
                out = model(batch.to_dict(), training=False)
                val_loss += loss_fn(pred=out, ref=batch).item()
                n_struct += batch.num_graphs
        val_loss /= max(n_struct, 1)
        scheduler.step(val_loss)
        print(f"epoch {epoch:3d}  val_loss {val_loss:.4f}")

        if val_loss < best_val:
            best_val = val_loss
            torch.save(model.state_dict(), out_dir / "best.pt")

    print(f"Best validation loss: {best_val:.4f}")
```

A few non-obvious points are worth flagging.

1. **`training=True` versus `training=False`.** Setting `training=True`
   tells MACE to compute and backpropagate forces by autograd. In
   `eval`-mode validation we still need gradient flow (forces are
   produced by differentiating the energy), but we do not want
   gradients with respect to weights. The model handles this
   correctly when both `model.eval()` and `torch.enable_grad()` are
   active.
2. **Gradient clipping.** Without `clip_grad_norm_`, NaN losses
   appear roughly one time in five during the first epoch — a
   single batch with an unusual configuration can produce an
   explosive gradient that destabilises training. A clip value of
   $10$ is conservative and rarely triggers in steady state.
3. **`avg_num_neighbors`.** This is the average count of neighbours
   within $r_\mathrm{c}$, used to normalise the per-atom feature sums.
   The right value can be measured by iterating once through the
   training set; for water with $r_\mathrm{c} = 5\,\text{\AA}$ it is
   approximately $24$.

## 9.6.5 Validation

After training, evaluate on the held-out test set.

```python
import numpy as np
import matplotlib.pyplot as plt

def evaluate(model, test_path, z_table, r_max, device):
    test_frames = ase.io.read(test_path, index=":")
    e_pred, e_dft, f_pred, f_dft = [], [], [], []
    for atoms in test_frames:
        data = AtomicData.from_atoms(atoms, z_table=z_table, cutoff=r_max).to(device)
        out = model(data.to_dict(), training=False)
        n = len(atoms)
        e_pred.append(out["energy"].item() / n)
        e_dft.append(atoms.info["energy"] / n)
        f_pred.append(out["forces"].detach().cpu().numpy())
        f_dft.append(atoms.arrays["forces"])
    e_pred, e_dft = np.array(e_pred), np.array(e_dft)
    f_pred = np.concatenate(f_pred);  f_dft = np.concatenate(f_dft)
    return e_pred, e_dft, f_pred, f_dft


e_pred, e_dft, f_pred, f_dft = evaluate(model, "test.xyz", z_table,
                                        cfg.r_max, device)
energy_mae  = np.mean(np.abs(e_pred - e_dft)) * 1000      # meV/atom
energy_rmse = np.sqrt(np.mean((e_pred - e_dft)**2)) * 1000
force_mae   = np.mean(np.abs(f_pred - f_dft)) * 1000      # meV/AA
force_rmse  = np.sqrt(np.mean((f_pred - f_dft)**2)) * 1000
print(f"Energy MAE/RMSE: {energy_mae:.2f} / {energy_rmse:.2f}  meV/atom")
print(f"Force  MAE/RMSE: {force_mae:.2f} / {force_rmse:.2f}  meV/AA")
```

For the configuration above on the revPBE0+D3 water dataset, typical
results are

| Metric | Value |
|---|---|
| Energy MAE | $0.5$ meV/atom |
| Energy RMSE | $0.8$ meV/atom |
| Force MAE | $25$ meV/Å |
| Force RMSE | $40$ meV/Å |

These are competitive with the reported accuracies of bespoke
revPBE0 water potentials. The force MAE of $25\,\mathrm{meV}/\text{\AA}$
is roughly the magnitude of the residual DFT integration noise itself,
so further reduction is bounded by the data quality.

Make parity plots:

```python
fig, axes = plt.subplots(1, 2, figsize=(8, 4))
axes[0].scatter(e_dft, e_pred, s=4, alpha=0.5)
lim = [min(e_dft.min(), e_pred.min()), max(e_dft.max(), e_pred.max())]
axes[0].plot(lim, lim, "k--", lw=0.8)
axes[0].set(xlabel="DFT energy (eV/atom)", ylabel="MACE energy (eV/atom)")
axes[1].scatter(f_dft.ravel(), f_pred.ravel(), s=1, alpha=0.3)
lim = [f_dft.min(), f_dft.max()]
axes[1].plot(lim, lim, "k--", lw=0.8)
axes[1].set(xlabel="DFT force (eV/Å)", ylabel="MACE force (eV/Å)")
fig.tight_layout()
fig.savefig("parity.png", dpi=200)
```

A correct parity plot is a tight diagonal cloud with no systematic
curvature. Bowing, fanning, or off-diagonal clusters indicate a
problem: bowing usually means insufficient training data; fanning
means the model is biased low; off-diagonal clusters mean a subset of
configurations is poorly represented.

## 9.6.6 Using the trained model — ASE calculator

`mace-torch` provides an ASE `Calculator` wrapper:

```python
from mace.calculators import MACECalculator
from ase import units
from ase.io import read, write
from ase.md.langevin import Langevin

calc = MACECalculator(
    model_paths=["water_mace/best.pt"],
    device="cuda",
    default_dtype="float64",
)

atoms = read("water_revpbe0.xyz", index=0)
atoms.calc = calc

# 10 ps of Langevin MD at 300 K, dt = 0.5 fs
dyn = Langevin(
    atoms,
    timestep=0.5 * units.fs,
    temperature_K=300,
    friction=0.01 / units.fs,
)
dyn.run(20_000)
write("mace_md.xyz", atoms)
```

A few caveats. `default_dtype="float64"` is *not* the default in
`mace-torch`; the default is `float32`, which trains faster but
produces slightly worse force accuracy and can drift energy in NVE
runs. For production MD use `float64`. The Langevin friction
$\gamma = 0.01\,\mathrm{fs}^{-1}$ is on the heavy side, suitable for
fast equilibration; for production sampling reduce to
$\gamma \approx 0.001\,\mathrm{fs}^{-1}$ to perturb dynamics less.

Cross-reference Chapter 7 for a deeper treatment of MD integrators
and thermostats; the MLIP is the calculator, the rest of the
machinery is the same as for a classical force field.

## 9.6.7 Active learning — a sketch

A single training run rarely produces a production-ready potential.
The trained model is good in the configurational neighbourhood of the
training data, but configurations encountered during long MD —
particularly rare events such as defect motion, surface diffusion, or
chemical reaction — may sit outside this neighbourhood. Active
learning closes the loop by detecting these out-of-distribution
configurations and adding DFT labels for them iteratively.

The basic loop is:

1. Run MD with the current MLIP.
2. Detect uncertain configurations — either by ensemble disagreement
   (train $K=3$ MACE models with different random seeds and flag
   configurations where their force predictions disagree by more than
   a threshold) or by Bayesian uncertainty (GAP). MACE typically uses
   ensembling.
3. Run DFT on the flagged configurations.
4. Add the new labels to the training set and retrain (or fine-tune
   from the previous checkpoint).
5. Repeat until uncertainties stay below threshold for the duration
   of a target-length trajectory.

Chapter 11 develops this in detail, including how to set the
disagreement threshold, how to budget DFT calls, and how to detect
*genuine* uncertainty versus mere variance from random initialisation.

## 9.6.8 Common training failures and remedies

**NaN loss within the first few batches.**
Almost always one of: (i) a configuration with an unphysically short
bond ($r_{ij} < 0.5\,\text{\AA}$), often from an unrelaxed structure,
which produces an exploding gradient; (ii) an isolated-atom energy
$E_0$ that is wrong by tens of eV, so the regression target is huge;
or (iii) batch size $> 1$ with structures whose atom counts span an
order of magnitude, producing imbalanced loss contributions. Fix:
visually inspect a few configurations, double-check the $E_0$ values,
and confirm gradient clipping is active.

**Validation loss decreases then steadily increases.**
Classical overfitting. The training set is too small, or the model is
too large, or training has run too long. Fix: reduce model capacity
(fewer channels, fewer layers), increase weight decay, or stop
training at the validation minimum (the script above already saves
the best checkpoint).

**Force MAE is good but energy MAE is poor.**
The energy weight is too low. Try `energy_weight = 10.0` instead of
$1.0$. Conversely, if energy is great but forces are poor, lower
`energy_weight` or raise `forces_weight`.

**Test parity plot is a tight diagonal in the bulk and a tail of
outliers.**
The outliers are configurations from a region of phase space
underrepresented in training. This is the active-learning signal:
identify the outliers (by force-error), label them with DFT, and add
to training.

**Inference is slow on large systems.**
Confirm `default_dtype="float64"` is the bottleneck (try `float32`
for production MD if the slight accuracy loss is acceptable). Confirm
the system fits in GPU memory (a $1000$-atom MACE inference uses
roughly $4\,\mathrm{GB}$ for the default configuration; larger
systems need either a bigger GPU, batched domain decomposition, or a
smaller `hidden_irreps`).

**Energy drift in NVE.**
The classic symptom of insufficient cutoff smoothness or
`float32`-induced rounding. Switch to `float64`, increase
`num_polynomial_cutoff` from $5$ to $6$, and check that the cutoff
$r_\mathrm{c} = 5\,\text{\AA}$ encloses the first two coordination
shells of every species in the system.

## 9.6.9 What we have

You now have a working MACE potential for liquid water, an ASE
calculator interface, validated accuracy in line with bespoke
literature potentials, and a checklist for diagnosing training
failures. Chapter 11 closes the loop with active learning; Chapter 12
shows how the foundation-model paradigm (MACE-MP-0) makes this entire
workflow shorter by reusing a pre-trained checkpoint.

The exercises that close this chapter walk you through deriving the
core ideas yourself, from the analytical force expression of a BPNN
to a from-scratch SOAP invariance proof, ending with a tiny MACE
training run on a fifty-configuration toy dataset.
