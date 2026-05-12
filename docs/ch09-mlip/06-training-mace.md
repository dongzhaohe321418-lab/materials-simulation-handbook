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

### Every parameter explained

Each field of the configuration deserves a careful look — wrong
choices here are the leading cause of MACE training failures.

| Field | Default | What it controls |
|---|---|---|
| `r_max` | 5.0 Å | Radial cutoff. Smaller → faster but loses long-range context; larger → slower and risks correlated neighbours. Match to first two coordination shells. |
| `num_layers` | 2 | Number of message-passing iterations. Each layer extends the receptive field by `r_max`. Two is sufficient for most problems. |
| `hidden_irreps` | "128x0e + 128x1o" | Equivariant feature space; see below for syntax. |
| `max_ell` | 1 | Maximum spherical-harmonic order used in tensor products. Going to 2 buys $\sim 30\%$ accuracy gain at $\sim 30\%$ extra cost. |
| `correlation` | 3 | Per-layer body order $\nu$. Three captures four-body correlations within one layer. |
| `num_radial_basis` | 8 | Bessel functions in the radial basis. Eight is plenty; sixteen rarely helps. |
| `radial_mlp` | (64, 64, 64) | Hidden widths of the radial MLP that maps $R(r)$ to weights. Default works for almost everything. |
| `batch_size` | 5 | Structures per gradient update. Larger batches stabilise gradients but slow per-epoch wall time. |
| `lr` | 0.01 | Initial Adam learning rate. Start here; the scheduler will decay. |
| `energy_weight` | 1.0 | Loss weight for energies. Use $1.0$ as the unit. |
| `forces_weight` | 100.0 | Loss weight for forces. With energy weight $1$, this gives the standard $1:100$ ratio. |
| `stress_weight` | 1.0 | Loss weight for stresses. Use $1$–$10$ when stresses are present. |
| `weight_decay` | 5e-7 | $L^2$ regularisation on weights. Small — most regularisation comes from architecture and early stopping. |

The `hidden_irreps` field deserves further comment. The string
`"128x0e + 128x1o"` declares the equivariant feature space: 128
channels of even scalars (irrep $0e$) plus 128 channels of odd
vectors (irrep $1o$). For organic chemistry with hydrogen bonding,
including the $\ell = 1$ channels is helpful; pushing to $\ell = 2$
($+ \mathtt{128x2e}$) buys further accuracy at $\sim 30 \%$ extra
inference cost.

The e3nn naming convention: each block is `{channels}x{ell}{parity}`,
where channels is a positive integer, $\ell \in \{0, 1, 2, \dots\}$ is
the irrep order, and parity is $e$ (even) or $o$ (odd). For physical
quantities: scalars are $0e$, pseudo-scalars are $0o$, position
vectors are $1o$, axial vectors are $1e$, rank-2 symmetric tensors are
$2e$. The parity selection enforces conservation of parity through the
network (relevant if your system has chiral content).

A common scaling: doubling the channel count roughly doubles the cost
and adds $\sim 0.1\,\mathrm{meV/atom}$ accuracy at the typical
training-set sizes of $10^3$–$10^4$ structures. Beyond $256$ channels
the returns diminish rapidly, and one is better served by adding more
training data.

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

### Training curves to expect

A well-configured MACE training on the 1000-frame water dataset
should produce learning curves with the following milestones:

- **Epoch 1–5:** Energy and force losses drop by an order of magnitude
  as the model learns the rough scale of the energies and the
  approximate force directions. Validation loss tracks training loss
  closely.
- **Epoch 5–30:** Energy MAE drops from $\sim 5$ to $\sim 1\,\mathrm{meV/atom}$;
  force MAE from $\sim 200$ to $\sim 50\,\mathrm{meV}/\text{\AA}$.
  Learning-rate scheduler may step down once.
- **Epoch 30–100:** Slow refinement. Energy MAE plateaus around
  $0.5\,\mathrm{meV/atom}$; force MAE around $25$–$30\,\mathrm{meV}/\text{\AA}$.
  Train/val gap of $\sim 20\%$ is normal; larger gaps suggest
  overfitting.
- **Epoch 100+:** Diminishing returns. Force MAE continues to creep
  down at $\sim 1\,\mathrm{meV}/\text{\AA}$ per 50 epochs but
  validation loss may plateau. Save the best checkpoint and stop.

If your curves deviate substantially — e.g. losses plateau after the
first epoch, or oscillate violently — refer to the common-failure
list in §9.6.8.

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

### A full validation suite

Parity plots are necessary but not sufficient. A production-ready
MACE potential must pass four independent tests; each catches failure
modes the others miss.

**(i) Parity plots on the held-out test set.** As above. Look for
$R^2 > 0.999$ on energies, $R^2 > 0.99$ on forces. Histograms of
errors should be tight Gaussians with no heavy tails.

**(ii) MD stability at multiple temperatures.** Run NVE simulations
for $50\,\mathrm{ps}$ at $300, 600, 1500\,\mathrm{K}$. Check:

- Energy drift in NVE: should be $< 1\,\mathrm{meV/atom}$ per
  $\mathrm{ps}$, ideally $< 0.1\,\mathrm{meV/atom}$ per $\mathrm{ps}$.
- No unphysical configurations (atoms within $0.5\,\mathrm{\AA}$
  of one another, exploded structures).
- Conservation of total momentum and angular momentum (only if no
  thermostat).

The $1500\,\mathrm{K}$ test is the most informative: it pushes the
model out of the training distribution and reveals where it breaks.
Most well-trained MACE potentials for water survive $300$–$600\,\mathrm{K}$
without problem but show drift or instability at $> 1000\,\mathrm{K}$
unless the training data included hot configurations.

**(iii) Radial distribution function (RDF) agreement.** Compute the
RDF $g_\mathrm{OO}(r)$ from a $100\,\mathrm{ps}$ trajectory in the
NVT ensemble at $300\,\mathrm{K}$ and compare to the DFT-MD reference
(from the same training-data source) or to experimental neutron
scattering. The first peak position and height are the most
sensitive features. Discrepancies of $\sim 1\%$ in peak position or
$\sim 5\%$ in peak height are typical and acceptable; larger
discrepancies indicate force-fitting problems.

**(iv) Phonon spectrum sanity.** Run a phonon calculation (small
displacements, finite-difference dynamical matrix) on a reference
crystal or, for liquids, compute the vibrational density of states
from the velocity autocorrelation function. Check:

- No imaginary modes in stable phases (Γ-point phonons all positive).
- O–H stretch in the right range ($3000$–$3500\,\mathrm{cm}^{-1}$
  for water).
- Bending modes ($1600\,\mathrm{cm}^{-1}$), librational ($600\,\mathrm{cm}^{-1}$),
  and translational ($200\,\mathrm{cm}^{-1}$) bands present.

Imaginary modes are an unambiguous failure signal: the potential
predicts unstable phonons in a stable phase, meaning the energy surface
has the wrong curvature at equilibrium. Re-train with more training
data on near-equilibrium configurations.

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

## 9.6.6a Inference timing and GPU memory

For production deployment, two practical numbers matter: how fast
the trained model evaluates and how much GPU memory it needs.

For a default-configuration MACE (128x0e + 128x1o, 2 layers,
`max_ell=1`) on an NVIDIA A100:

- **192-atom water box**: $\sim 4\,\mathrm{ms}$ per energy+force
  evaluation, $\sim 800\,\mathrm{MB}$ peak GPU memory.
- **1000-atom water box**: $\sim 18\,\mathrm{ms}$ per evaluation,
  $\sim 3.5\,\mathrm{GB}$ peak GPU memory.
- **10,000-atom box**: $\sim 200\,\mathrm{ms}$ per evaluation,
  $\sim 30\,\mathrm{GB}$ peak GPU memory — close to the A100 limit.

Memory scales roughly $O(N \times C^2 \times \ell_\mathrm{max}^2)$
where $C$ is the channel count and $\ell_\mathrm{max}$ is the
maximum irrep. For very large systems, halving `hidden_irreps` from
$128$ to $64$ quarters the memory at the cost of $\sim 0.5\,\mathrm{meV/atom}$
extra energy error.

For systems beyond $\sim 30,000$ atoms, even an A100 cannot fit a
full forward pass. Domain decomposition splits the system into
spatial tiles with halo regions of width $T \times r_\mathrm{c}$
(for $T$ layers); each tile runs on a separate GPU. The `mace-torch`
package supports this via `MACECalculator(...)` integrated with the
LAMMPS pair style `mace`.

For comparison: a classical force field (TIP4P/2005 water) evaluates
a 1000-atom box in $\sim 50\,\mu\mathrm{s}$ on a single CPU core,
so $\sim 400\times$ faster than MACE. The MACE gives DFT-quality
forces; TIP4P/2005 gives parameterised approximations. The choice
between them is the topic of §9.1.4a.

## 9.6.6b Stress predictions and constant-pressure MD

For NPT simulations and equation-of-state calculations, the model
must predict not only energies and forces but the *stress tensor*

$$
\sigma_{\alpha\beta} = -\frac{1}{V} \sum_i \mathbf{r}_i^\alpha \mathbf{F}_i^\beta
                     + \text{kinetic part},
$$

where $V$ is the cell volume. MACE computes the stress automatically
via autograd on $U$ with respect to a virtual strain applied to the
cell:

$$
\sigma_{\alpha\beta} = -\frac{1}{V} \frac{\partial U}{\partial \epsilon_{\alpha\beta}}
\Big|_{\epsilon = 0}.
$$

If stress targets are provided in the training data (the `stress=...`
field in extxyz), `mace-torch` adds a stress-matching term to the
loss with weight `stress_weight`. The typical convention
$\lambda_E : \lambda_F : \lambda_\mathrm{stress} = 1 : 100 : 10$
gives roughly equal residual error on each. Without stress training,
stresses are still predicted (autograd derives them from the energy)
but their accuracy is unconstrained — adequate for monitoring but
not for high-precision EOS work.

For NPT MD with an MLIP, use the ASE `NPT` integrator (Parrinello–Rahman)
with the MACE calculator as usual. Equilibration to a target pressure
typically takes 5–10 ps; production runs of $50$–$100\,\mathrm{ps}$
suffice for converged density and bulk modulus.

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

A practical preview: a typical active-learning loop adds 50–200 new
DFT configurations per iteration and converges in 3–10 iterations
for a single chemistry. The total DFT cost is typically $\sim 10\times$
less than would be needed to "guess" the right training set upfront,
because active learning targets exactly the configurations the model
finds unfamiliar.

The two operational decisions in any active-learning loop:

- **The ensemble size $K$.** Three models is the practical minimum
  for a meaningful variance estimate; five is the typical
  production value. Each ensemble member is trained from scratch with
  a different random seed, on the same training data.
- **The disagreement threshold.** Set as a multiple of the
  training-set average force disagreement: typically $3\sigma$ to
  $5\sigma$. Too low and you label many configurations the ensemble
  is merely *slightly* uncertain about, wasting DFT calls. Too high
  and you miss genuine extrapolations.

Cross-reference Chapter 11 for code.

## 9.6.8 Common training failures and remedies

**Force MAE doesn't drop below $\sim 50\,\mathrm{meV}/\text{\AA}$.**
You have probably saturated the expressive power of an $\ell = 1$
network on a system with angular structure beyond what $\ell = 1$
can capture. Try increasing `max_ell` from $1$ to $2$ and adding
`+ 128x2e` to `hidden_irreps`. Expect a $30\%$ slowdown and a
$\sim 2\times$ reduction in force MAE. For systems with strong sp$^3$
hybridisation (water, ammonia, hydrocarbons), $\ell = 2$ is often
the sweet spot.

**MD trajectory unstable — atoms fly apart within a few ps.**
The training data does not cover the configurations encountered in
MD. Two distinct sub-failures:

- **Force noise too large.** A trained potential with force MAE
  comparable to thermal force fluctuations ($\sim k_B T / a$ where
  $a$ is a typical bond length) will produce a random-walk-like
  drift in the energy. Solution: shrink the force MAE, either by
  more training data or by raising `max_ell`.
- **Out-of-distribution configurations.** The MD visits regions of
  configuration space the training set did not cover (e.g. a
  transient compressed bond near a collision). The potential
  extrapolates wildly. Solution: train on hotter configurations
  (a $1500\,\mathrm{K}$ ab-initio MD sub-trajectory adds the
  high-energy tail to the training distribution), or run active
  learning to add the offending configurations.

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

### A summary checklist before deploying a potential

Before declaring victory and using your potential for production MD,
run through this checklist. Most published "MLIP works for system X"
papers omit one or more of these tests, and the literature suffers
for it.

- [ ] Energy MAE on held-out test set $< 1\,\mathrm{meV/atom}$.
- [ ] Force MAE on held-out test set $< 50\,\mathrm{meV}/\text{\AA}$
      (better $< 30$).
- [ ] Stress MAE $< 10^{-3}\,\mathrm{eV}/\text{\AA}^3$ if you trained
      on stresses.
- [ ] NVE energy drift $< 1\,\mathrm{meV/atom/ps}$ at the target
      temperature.
- [ ] MD trajectory stable for $\geq 100\,\mathrm{ps}$ at the target
      temperature and at $1.5\times$ the target temperature (an
      out-of-distribution stress test).
- [ ] Radial distribution function matches DFT-MD reference within
      $5\%$ in peak heights.
- [ ] No imaginary phonon modes in stable phases.
- [ ] Density (NPT equilibration) within $2\%$ of the experimental
      or DFT-MD value.
- [ ] If you use forces in any downstream analysis (e.g. infrared
      spectra from dipole-velocity correlation), verify the
      derivatives are smooth enough — increase
      `num_polynomial_cutoff` if necessary.

Only after all of these pass should the potential be considered
production-ready. The validation effort is typically 20–30% of the
total time investment in a new MLIP, and it is the right place to
spend that time.

## 9.6.10 What to read next

- Chapter 10 (enhanced sampling) — now that MD is cheap, use it
  with metadynamics, umbrella sampling, transition-path sampling.
- Chapter 11 (active learning) — automate the configuration-set
  expansion.
- Chapter 12 (foundation models) — MACE-MP-0 zero-shot evaluation
  and fine-tuning workflows.
- Behler, *Constructing high-dimensional neural network potentials*,
  Int. J. Quantum Chem. 115, 1032 (2015) — a thorough review of the
  BPNN family.
- Drautz, *Atomic cluster expansion for accurate and transferable
  interatomic potentials*, Phys. Rev. B 99, 014104 (2019) — the ACE
  paper.
- Batatia et al., *MACE: Higher Order Equivariant Message Passing
  Neural Networks for Fast and Accurate Force Fields*, NeurIPS 2022 —
  the MACE paper.
- Batatia et al., *A foundation model for atomistic materials
  chemistry*, arXiv:2401.00096 (2024) — MACE-MP-0.

The exercises that close this chapter walk you through deriving the
core ideas yourself, from the analytical force expression of a BPNN
to a from-scratch SOAP invariance proof, ending with a tiny MACE
training run on a fifty-configuration toy dataset.
