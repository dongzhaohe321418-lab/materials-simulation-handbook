# Methods — Melting point of copper via MLIP-driven MD

The protocol has four phases:

1. Generate DFT training data.
2. Train and validate the MACE potential.
3. Equilibrate solid and liquid Cu with the MLIP.
4. Run two-phase coexistence to extract $T_m$.

Read the whole document before doing any of it.

---

## Phase 0 — Environment setup

```bash
conda create -n cu-mlip python=3.11
conda activate cu-mlip
conda install -c conda-forge ase pymatgen matplotlib numpy
pip install mace-torch lammps
# QE for DFT data:
conda install -c conda-forge qe
# Optional: LAMMPS with MACE interface
# Download from https://github.com/ACEsuit/mace and follow the LAMMPS
# patch instructions; verify with `lmp -h | grep mace`.
```

Verify each component:

```bash
python -c "import mace; print(mace.__version__)"
pw.x < /dev/null   # should print QE banner and exit
lmp -h | grep mace # should list 'mace' pair_style
```

---

## Phase 1 — DFT training data

### 1.1 Choose the DFT settings

Use PBE with the SSSP "efficiency" Cu pseudopotential. Convergence
parameters:

- $E_\mathrm{cut}^\mathrm{wfc} = 60$ Ry (Cu needs more than Si because
  of the d electrons).
- $E_\mathrm{cut}^\mathrm{rho} = 720$ Ry (= 12 × wavefunction cut-off
  to be safe with ultrasoft).
- **k**-mesh: $4 \times 4 \times 4$ on the 4-atom conventional cell;
  scale to equivalent meshes on larger cells.
- Marzari–Vanderbilt cold smearing, `degauss = 0.01 Ry`.
- `mixing_beta = 0.3` for transition metals.

### 1.2 Build the training-set sampling cells

Use a 32-atom Cu cell (a $2 \times 2 \times 2$ replication of the
4-atom conventional fcc cell) for the DFT-MD trajectories.

Three thermodynamic states:

| Label | Temperature | Cell | Duration | Sampling interval |
| --- | --- | --- | --- | --- |
| `solid_300K` | 300 K | 32-atom fcc | 2 ps | every 20 fs |
| `solid_800K` | 800 K | 32-atom fcc | 2 ps | every 20 fs |
| `liquid_1500K` | 1500 K | 32-atom (pre-melted) | 4 ps | every 20 fs |

This gives 100 + 100 + 200 = 400 frames. Add 100 more from strained
cells (volumetric strain at ±1, ±2, ±3, ±5 % and isolated shear of
±2 %) on the static lattice. Total ≈ 500 frames.

### 1.3 Pre-melt for the liquid run

A 32-atom Cu cell at 1500 K will not melt within 4 ps from the
crystalline initial condition because of the superheating barrier.
Pre-melt it: take the cell to 3000 K for 1 ps (which definitely
liquefies it), then quench-equilibrate to 1500 K for 1 ps. Discard
the warm-up frames. The remaining 2-ps trajectory at 1500 K samples
the equilibrium liquid.

### 1.4 Run AIMD

Use Born–Oppenheimer MD (`calculation = 'md'`, `ion_dynamics =
'verlet'`) with a 1-fs timestep. A representative `pw.x` input:

```text
&CONTROL
  calculation = 'md'
  prefix      = 'cu_solid_300K'
  pseudo_dir  = '...'
  outdir      = './tmp_cu_300/'
  nstep       = 2000
  dt          = 41.34          ! 1 fs in atomic units
/
&SYSTEM
  ibrav = 0
  nat = 32, ntyp = 1
  ecutwfc = 60.0, ecutrho = 720.0
  occupations = 'smearing', smearing = 'mv', degauss = 0.01
/
&ELECTRONS
  conv_thr = 1.0d-7
  mixing_beta = 0.3
/
&IONS
  ion_dynamics = 'verlet'
  ion_temperature = 'rescale-v'
  tempw = 300.0
  nraise = 50
/
ATOMIC_SPECIES
  Cu  63.546  Cu.pbe-spn-rrkjus_psl.1.0.0.UPF
CELL_PARAMETERS angstrom
  ... (3.615 * 2 * <identity>)
ATOMIC_POSITIONS angstrom
  ... (32 Cu atoms)
K_POINTS automatic
  2 2 2 0 0 0
```

Submit each trajectory; check that the kinetic temperature (printed
by QE) hovers around the requested value after the first few hundred
steps.

### 1.5 Convert to extended-XYZ

Use ASE to read the QE output trajectories and write a single
`extxyz` file with energies, forces, and stresses:

```python
from ase.io import read, write
frames = []
for label in ["cu_solid_300K", "cu_solid_800K", "cu_liquid_1500K"]:
    traj = read(f"runs/{label}.pwo", index="::20")  # every 20 fs
    frames.extend(traj)
write("data/cu_all.xyz", frames)
```

Make sure each frame has `info["energy"]`, `arrays["forces"]`, and
`info["stress"]` populated. MACE's training loader expects these keys.

### 1.6 Train/val/test split

Random 80/10/10 split. Use a fixed random seed. Keep the test set
**hidden** until after hyperparameter tuning.

---

## Phase 2 — Train MACE

### 2.1 Hyperparameters

Sensible defaults for a small-system MLIP on a transition metal:

| Parameter | Value | Note |
| --- | --- | --- |
| `r_max` | 5.0 Å | Cu nearest-neighbour distance is 2.56 Å; 5 Å covers third shell. |
| `num_interactions` | 2 | Two-layer message passing → effective receptive field ≈ 10 Å. |
| `max_L` | 2 | $L = 2$ is the standard balance of accuracy and cost. |
| `hidden_irreps` | `128x0e+128x1o` | 128 invariant + 128 vector features. |
| `correlation` | 3 | Body order 4 per layer. |
| `batch_size` | 8 | |
| `lr` | 0.01 → 0.001 with cosine annealing | |
| `epochs` | 200 | |
| `loss_weights` | E:1, F:100, S:1 | Force-weight 100 is standard. |

### 2.2 Launch the training

`mace-torch` exposes a CLI `mace_run_train`. A typical config:

```bash
mace_run_train \
  --name="cu_mace_v1" \
  --train_file="data/cu_all_train.xyz" \
  --valid_file="data/cu_all_val.xyz" \
  --E0s='{"Cu": -1234.567}' \
  --r_max=5.0 \
  --max_L=2 \
  --correlation=3 \
  --hidden_irreps='128x0e+128x1o' \
  --num_interactions=2 \
  --batch_size=8 \
  --max_num_epochs=200 \
  --lr=0.01 \
  --weight_decay=5e-7 \
  --energy_weight=1.0 \
  --forces_weight=100.0 \
  --stress_weight=1.0 \
  --device=cuda
```

The `E0s` value is the isolated-atom energy of Cu in your DFT setting.
Compute it once by running a single Cu atom in a large vacuum box.

Expect 6–12 hours on an A100 for ~ 500 frames.

### 2.3 Validate

After training, generate parity plots on the test set:

```python
from ase.io import read
from mace.calculators import MACECalculator
import numpy as np

calc = MACECalculator(model_path="cu_mace_v1.model", device="cuda")
test_frames = read("data/cu_all_test.xyz", index=":")
e_dft, e_mlip, f_dft, f_mlip = [], [], [], []
for f in test_frames:
    e_dft.append(f.info["energy"] / len(f))
    f_dft.extend(f.arrays["forces"].ravel().tolist())
    f.calc = calc
    e_mlip.append(f.get_potential_energy() / len(f))
    f_mlip.extend(f.get_forces().ravel().tolist())
print("E MAE (meV/atom):", 1000 * np.mean(np.abs(np.array(e_dft) - np.array(e_mlip))))
print("F MAE (meV/Å):", 1000 * np.mean(np.abs(np.array(f_dft) - np.array(f_mlip))))
```

Target metrics:

- E MAE < 5 meV/atom on the held-out test set.
- F MAE < 50 meV/Å on the held-out test set.

If you miss the F target, the usual reason is that the liquid is
under-represented; add another 100 frames from a longer 1500 K
AIMD run and retrain.

---

## Phase 3 — Equilibrate solid and liquid

Build LAMMPS input cells using ASE's `lammpsdata` writer:

- A 6 × 6 × 6 fcc supercell of Cu: 864 atoms (this is your "small"
  test cell).
- A 8 × 8 × 32 elongated fcc supercell for coexistence (described
  below): 8192 atoms.

Run NPT for each:

- `solid_1300K`: NPT at $T = 1300$ K, $p = 1$ bar, 100 ps.
- `solid_1500K`: NPT at 1500 K, 1 bar, 100 ps.
- `liquid_1700K`: melt the 1500 K configuration at 2500 K for 20 ps,
  then equilibrate at 1700 K for 100 ps.

LAMMPS pair-style for MACE:

```text
pair_style mace no_domain_decomposition
pair_coeff * * cu_mace_v1.model-lammps.pt Cu
```

Time step: 1.0 fs. Nosé–Hoover thermostat/barostat with relaxation
times 0.1 ps and 1.0 ps respectively.

Monitor:

- Volume (should be stable to ≈ 1 % after equilibration).
- Mean-squared displacement (should be diffusive for liquid, plateau
  for solid).
- Radial distribution function (use `compute rdf` in LAMMPS).

---

## Phase 4 — Two-phase coexistence

### 4.1 Construct the coexistence cell

Build an elongated supercell, e.g., $8 \times 8 \times 32$ fcc
conventional cells along the **z** axis. Equilibrate the whole cell
at $T_\mathrm{cold} = 800$ K as solid, then take *half* the cell
along **z** (atoms with $z > L_z/2$) and re-equilibrate that half
at $T_\mathrm{hot} = 2500$ K *with positions fixed for $z < L_z/2$*.
This produces a cell that is solid in one half and liquid in the
other. Save this as `coex_init.lmp`.

### 4.2 Run NPT coexistence at trial temperatures

For each trial temperature $T \in \{1300, 1350, 1400, 1450\}$ K:

1. Read `coex_init.lmp`.
2. Run NPT at $(T, 1\,\text{bar})$ for 200 ps.
3. At every picosecond, compute the position of the solid–liquid
   interface by analysing the local Steinhardt order parameter $q_6$
   along **z** (LAMMPS `compute orientorder/atom`). The interface is
   the position where the running average of $q_6$ along **z** crosses
   the midpoint between the solid and liquid bulk values.
4. Compute the interface velocity $v(T) = dz_\mathrm{interface}/dt$
   by linear regression of position vs time after the first 30 ps
   (the equilibration period).

### 4.3 Extract $T_m$

Plot $v(T)$ versus $T$. The melting temperature is the zero-crossing,
which you obtain by linear interpolation between the two trial
temperatures that bracket $v = 0$. For Cu with a well-trained MACE,
you should find $T_m$ in the range 1300–1400 K.

Report the result with an uncertainty derived from:

- The standard error of the linear fit $v(T)$ → $T_m$.
- A repeat of the coexistence run at the bracketed $T_m$ with two
  additional random thermostat seeds, to bound the run-to-run scatter.

### 4.4 Hysteresis cross-check

In parallel, run a 6 × 6 × 6 NPT cell from 800 K to 2500 K at a
heating rate of 5 K/ps; identify the temperature at which the
volume jumps (the superheated melting temperature $T_+$). Run the
reverse cooling from 2500 K to 800 K; identify $T_-$ (the supercooled
crystallisation temperature). The true $T_m$ lies between $T_+$ and
$T_-$. For a clean Cu cell expect $T_+ \approx 1500$–1700 K and
$T_- \approx 800$–1000 K — so the hysteresis bracket is loose, but
must contain your coexistence-derived $T_m$.

---

## Pitfalls

1. **Forgetting to compute the isolated-atom energy.** Without a
   correct `E0`, the per-atom energy MAE is meaningless and the MLIP
   is offset by an uninteresting constant.
2. **Training on energies that include the PV term.** QE reports the
   total energy without the PV term by default; do not add it.
3. **Using `fix nvt` instead of `fix npt`.** A constant-volume
   coexistence will produce a meaningless pressure swing as the
   interface moves; you must use NPT.
4. **Cell too short along z.** If $L_z < 30$ Å, the periodic interface
   self-interacts. Use at least $L_z \approx 100$ Å.
5. **Coexistence run too short.** The first 30 ps are equilibration;
   any velocity estimated from a 50-ps run will be dominated by
   transients. Run at least 200 ps per temperature.
6. **Mixing single-precision MD with double-precision DFT.** MACE
   in PyTorch defaults to float32; this is fine for MD but introduces
   noise on the order of 0.1 meV/atom that you should be aware of.
7. **Reporting $T_m$ to four significant figures.** Your uncertainty
   is ≈ 50 K; report 1350 ± 50 K, not 1351 K.

---

## What "done" looks like

You have a number, $T_m = T \pm \sigma$, with $\sigma$ derived from a
proper uncertainty analysis. You can plot interface velocity versus
temperature with the zero-crossing marked. You can defend why your
training set is appropriate. You can name the dominant source of
error.
