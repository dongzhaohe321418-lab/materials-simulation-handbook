"""Generate Google-Colab-ready notebooks for the handbook's heavy chapters.

Chapters 6, 9, 10 and 12 depend on packages that cannot run under the
browser-based JupyterLite kernel (torch, MACE, ASE, pymatgen, torch-geometric)
and often want a GPU. These notebooks mirror the meatiest runnable code from
those chapters, with an install cell at the top, and are meant to be opened in
Google Colab. Run from the repository root::

    python scripts/build_colab_notebooks.py

This script is deterministic — re-running it overwrites the generated
notebooks/*-colab.ipynb files in place. Hand-edits to those files will be lost.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO_ROOT / "notebooks"
SITE_BASE = "https://dongzhaohe321418-lab.github.io/materials-simulation-handbook"


def md(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.splitlines(keepends=True) or [""],
    }


def code(source: str) -> dict:
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": source.splitlines(keepends=True) or [""],
    }


def notebook(cells: Iterable[dict]) -> dict:
    return {
        "cells": list(cells),
        "metadata": {
            "kernelspec": {
                "name": "python3",
                "display_name": "Python 3",
                "language": "python",
            },
            "language_info": {"name": "python"},
            "accelerator": "GPU",
            "colab": {"provenance": []},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def write_notebook(name: str, cells: list[dict]) -> Path:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / name
    path.write_text(json.dumps(notebook(cells), indent=1) + "\n")
    return path


def header(chapter_title: str, slug: str, summary: str) -> dict:
    url = f"{SITE_BASE}/{slug}/"
    text = (
        f"# {chapter_title} — Colab notebook\n\n"
        f"Back to the chapter: <{url}>\n\n"
        f"{summary}\n\n"
        f"This notebook is meant for **Google Colab** rather than the in-browser "
        f"JupyterLite kernel, because it needs heavy packages (and, where noted, "
        f"a GPU) that cannot run under Pyodide. Open it in Colab, and where the "
        f"install cell mentions it, switch the runtime to a GPU via "
        f"**Runtime -> Change runtime type -> GPU** before running the rest.\n"
    )
    return md(text)


# ---------------------------------------------------------------------------
# Chapter 6 — Running DFT in Practice
# ---------------------------------------------------------------------------
CH06 = [
    header(
        "Chapter 6 · Running DFT in Practice",
        "ch06-running-dft",
        "We rehearse the full DFT workflow — build a structure, attach a "
        "calculator, relax it, sweep an equation of state, plot the result. "
        "Quantum ESPRESSO is not pip-installable, so we use ASE's built-in EMT "
        "calculator as a stand-in: the *workflow* is identical to a real DFT "
        "run, only the energy model is cheaper. A clearly marked block shows "
        "how to swap EMT for the real Espresso calculator.",
    ),
    md("## Install\n\n"
       "`ase` is the only dependency here and installs in a few seconds. No GPU "
       "is needed for this notebook — EMT is a cheap analytic potential. A real "
       "DFT run would need a separate Quantum ESPRESSO or VASP installation, "
       "which is not pip-installable and is out of scope for a Colab notebook."),
    code("!pip install ase\n"),
    md("## Build a crystal structure with ASE\n\n"
       "We start exactly as a real DFT study would: construct the periodic "
       "cell. Here it is an aluminium FCC primitive cell. ASE's `bulk` builder "
       "knows the crystal structures of the elements."),
    code("from ase.build import bulk\n"
         "from ase.calculators.emt import EMT\n\n"
         "al = bulk('Al', crystalstructure='fcc', a=4.05)\n"
         "print(al)\n"
         "print('volume per atom:', al.get_volume() / len(al), 'A^3')\n"),
    md("## Attach a calculator and relax the cell\n\n"
       "We attach the EMT calculator — the stand-in for a DFT engine — and "
       "relax the atomic positions and cell with a filter plus an optimiser. "
       "This mirrors a DFT geometry optimisation step for step."),
    code("from ase.optimize import BFGS\n"
         "from ase.filters import FrechetCellFilter\n\n"
         "al.calc = EMT()\n"
         "print('initial energy:', al.get_potential_energy(), 'eV')\n\n"
         "relaxed = al.copy()\n"
         "relaxed.calc = EMT()\n"
         "opt = BFGS(FrechetCellFilter(relaxed), logfile=None)\n"
         "opt.run(fmax=0.001)\n"
         "print('relaxed energy:', relaxed.get_potential_energy(), 'eV')\n"
         "print('relaxed lattice constant:', relaxed.cell.cellpar()[0], 'A')\n"),
    md("## Sweep an equation of state\n\n"
       "Scan the lattice constant, compute the energy at each volume, and fit a "
       "Birch-Murnaghan equation of state to recover the equilibrium volume and "
       "bulk modulus. This is one of the most common DFT validation tasks."),
    code("import numpy as np\n\n"
         "def equation_of_state(a_values):\n"
         "    out = {'a': [], 'V': [], 'E': []}\n"
         "    for a in a_values:\n"
         "        atoms = bulk('Al', crystalstructure='fcc', a=a)\n"
         "        atoms.calc = EMT()\n"
         "        out['a'].append(a)\n"
         "        out['V'].append(atoms.get_volume())\n"
         "        out['E'].append(atoms.get_potential_energy())\n"
         "    return {k: np.array(v) for k, v in out.items()}\n\n"
         "a_grid = np.linspace(3.8, 4.3, 11)\n"
         "eos = equation_of_state(a_grid)\n"
         "for a, V, E in zip(eos['a'], eos['V'], eos['E']):\n"
         "    print(f'a = {a:.3f} A   V = {V:7.3f} A^3   E = {E:8.4f} eV')\n"),
    md("## Fit and plot the equation of state\n\n"
       "ASE ships a Birch-Murnaghan fitter. The fitted minimum is the "
       "predicted equilibrium volume; the curvature is the bulk modulus."),
    code("import matplotlib.pyplot as plt\n"
         "from ase.eos import EquationOfState\n\n"
         "eos_fit = EquationOfState(eos['V'], eos['E'])\n"
         "V0, E0, B = eos_fit.fit()\n"
         "print(f'equilibrium volume V0 = {V0:.3f} A^3')\n"
         "print(f'equilibrium energy E0 = {E0:.4f} eV')\n"
         "print(f'bulk modulus B = {B / 1e-3:.1f} meV/A^3')\n\n"
         "fig, ax = plt.subplots(figsize=(6, 4))\n"
         "eos_fit.plot(ax=ax)\n"
         "ax.set_title('Aluminium equation of state (EMT stand-in)')\n"
         "plt.show()\n"),
    md("## Using real DFT instead of EMT\n\n"
       "Everything above is the genuine DFT workflow; only the energy model is "
       "a stand-in. To run real density functional theory, install Quantum "
       "ESPRESSO separately (it is not pip-installable) and swap the EMT "
       "calculator for ASE's `Espresso` calculator. The rest of the notebook — "
       "`bulk`, `BFGS`, `EquationOfState` — is unchanged. The cell below is "
       "**not meant to run here**; it shows the one substitution."),
    code("# ---- TO USE REAL DFT: replace `EMT()` with the block below ----\n"
         "# Requires a working Quantum ESPRESSO install and pseudopotentials.\n"
         "#\n"
         "# from pathlib import Path\n"
         "# from ase.calculators.espresso import Espresso, EspressoProfile\n"
         "#\n"
         "# profile = EspressoProfile(\n"
         "#     command='pw.x',  # or 'mpirun -np 4 pw.x'\n"
         "#     pseudo_dir=Path.home() / 'pseudo/SSSP_1.3.0_PBE_efficiency',\n"
         "# )\n"
         "# al.calc = Espresso(\n"
         "#     profile=profile,\n"
         "#     pseudopotentials={'Al': 'Al.pbe-n-kjpaw_psl.1.0.0.UPF'},\n"
         "#     input_data={\n"
         "#         'system': {'ecutwfc': 50, 'ecutrho': 400,\n"
         "#                    'occupations': 'smearing', 'smearing': 'mv',\n"
         "#                    'degauss': 0.01},\n"
         "#         'electrons': {'conv_thr': 1e-9, 'mixing_beta': 0.4},\n"
         "#     },\n"
         "#     kpts=(8, 8, 8),\n"
         "# )\n"
         "# energy = al.get_potential_energy()  # now a genuine DFT energy\n"
         "print('See comment above: swap EMT for Espresso for production DFT.')\n"),
]


# ---------------------------------------------------------------------------
# Chapter 9 — Machine Learning Interatomic Potentials
# ---------------------------------------------------------------------------
CH09 = [
    header(
        "Chapter 9 · Machine Learning Interatomic Potentials",
        "ch09-mlip",
        "We train a tiny MACE potential end to end on a small synthetic dataset "
        "and check its parity against the reference. The dataset and model are "
        "deliberately small so the run finishes in a few minutes on a Colab "
        "GPU; the chapter's production recipe uses the same code with a larger "
        "dataset, more epochs and the official `mace_run_train` entry point.",
    ),
    md("## Install\n\n"
       "`mace-torch` pulls in the MACE library and its `e3nn` dependency. "
       "`torch` itself comes pre-installed on Colab, so it is not listed here. "
       "**Set the runtime to GPU** (Runtime -> Change runtime type -> GPU) "
       "before running the training cell — MACE training on CPU is painfully "
       "slow."),
    code("!pip install mace-torch\n"),
    md("## Check the GPU\n\n"
       "MACE training needs a GPU for a reasonable turnaround. If this cell "
       "fails, switch the Colab runtime type to GPU and re-run from the top."),
    code("import torch\n"
         "assert torch.cuda.is_available(), 'GPU required: Runtime -> Change runtime type -> GPU'\n"
         "print('device:', torch.cuda.get_device_name(0))\n"),
    md("## Build a tiny training dataset\n\n"
       "A real MACE study reads thousands of DFT-labelled structures from an "
       "extended-XYZ file. To keep this notebook self-contained and fast, we "
       "generate a small set of rattled bulk-copper cells and label them with "
       "ASE's cheap EMT calculator standing in for DFT. The workflow — write "
       "`train.xyz` / `valid.xyz` / `test.xyz` — is exactly the chapter's."),
    code("import numpy as np\n"
         "import ase.io\n"
         "from ase.build import bulk\n"
         "from ase.calculators.emt import EMT\n\n"
         "rng = np.random.default_rng(0)\n"
         "frames = []\n"
         "base = bulk('Cu', crystalstructure='fcc', a=3.61, cubic=True).repeat((2, 2, 2))\n"
         "for _ in range(120):\n"
         "    atoms = base.copy()\n"
         "    atoms.rattle(stdev=0.08, rng=rng)\n"
         "    atoms.calc = EMT()\n"
         "    atoms.info['energy'] = atoms.get_potential_energy()\n"
         "    atoms.arrays['forces'] = atoms.get_forces()\n"
         "    atoms.calc = None\n"
         "    frames.append(atoms)\n\n"
         "idx = rng.permutation(len(frames))\n"
         "ase.io.write('train.xyz', [frames[i] for i in idx[:90]])\n"
         "ase.io.write('valid.xyz', [frames[i] for i in idx[90:105]])\n"
         "ase.io.write('test.xyz',  [frames[i] for i in idx[105:]])\n"
         "print('wrote 90 train / 15 valid / 15 test structures')\n"),
    md("## Train a small MACE model\n\n"
       "We call MACE's official training entry point with a small architecture "
       "(two interaction layers, modest hidden dimension) and a short epoch "
       "budget. This is the same `mace_run_train` command the chapter uses for "
       "production, only with smaller knobs so it finishes quickly."),
    code("import sys\n"
         "from mace.cli.run_train import main as mace_run_train\n\n"
         "args = [\n"
         "    '--name', 'cu_mace',\n"
         "    '--train_file', 'train.xyz',\n"
         "    '--valid_file', 'valid.xyz',\n"
         "    '--test_file', 'test.xyz',\n"
         "    '--energy_key', 'energy',\n"
         "    '--forces_key', 'forces',\n"
         "    '--model', 'MACE',\n"
         "    '--r_max', '5.0',\n"
         "    '--num_interactions', '2',\n"
         "    '--hidden_irreps', '32x0e + 32x1o',\n"
         "    '--num_radial_basis', '8',\n"
         "    '--max_ell', '2',\n"
         "    '--correlation', '3',\n"
         "    '--batch_size', '5',\n"
         "    '--max_num_epochs', '30',\n"
         "    '--energy_weight', '1.0',\n"
         "    '--forces_weight', '100.0',\n"
         "    '--lr', '0.01',\n"
         "    '--device', 'cuda',\n"
         "    '--default_dtype', 'float64',\n"
         "    '--seed', '1',\n"
         "]\n"
         "sys.argv = ['mace_run_train'] + args\n"
         "mace_run_train()\n"),
    md("## Evaluate the trained potential\n\n"
       "We load the saved model with the `MACECalculator` — the same object "
       "you would attach to run MD — and compare its per-atom energies against "
       "the reference labels on the held-out test set."),
    code("from mace.calculators import MACECalculator\n\n"
         "calc = MACECalculator(model_paths=['cu_mace.model'], device='cuda',\n"
         "                      default_dtype='float64')\n\n"
         "test_frames = ase.io.read('test.xyz', index=':')\n"
         "e_pred, e_ref = [], []\n"
         "for atoms in test_frames:\n"
         "    ref = atoms.info['energy'] / len(atoms)\n"
         "    atoms.calc = calc\n"
         "    e_pred.append(atoms.get_potential_energy() / len(atoms))\n"
         "    e_ref.append(ref)\n"
         "e_pred, e_ref = np.array(e_pred), np.array(e_ref)\n"
         "mae = np.mean(np.abs(e_pred - e_ref)) * 1000\n"
         "print(f'energy MAE = {mae:.2f} meV/atom')\n"),
    md("## Parity plot\n\n"
       "The hallmark figure of every MLIP paper: predicted versus reference "
       "energy. Points hugging the diagonal mean the potential has learned the "
       "energy surface."),
    code("import matplotlib.pyplot as plt\n\n"
         "fig, ax = plt.subplots(figsize=(5, 5))\n"
         "ax.scatter(e_ref, e_pred, s=20, alpha=0.7)\n"
         "lim = [min(e_ref.min(), e_pred.min()), max(e_ref.max(), e_pred.max())]\n"
         "ax.plot(lim, lim, 'k--', lw=0.8)\n"
         "ax.set_xlabel('reference energy (eV/atom)')\n"
         "ax.set_ylabel('MACE energy (eV/atom)')\n"
         "ax.set_aspect('equal')\n"
         "ax.set_title('MACE parity on held-out test set')\n"
         "plt.show()\n"),
]


# ---------------------------------------------------------------------------
# Chapter 10 — Graph Neural Networks for Materials
# ---------------------------------------------------------------------------
CH10 = [
    header(
        "Chapter 10 · Graph Neural Networks for Materials",
        "ch10-gnn",
        "We turn a crystal into a graph, wrap it in a PyTorch Geometric `Data` "
        "object, and run a single forward pass of a small CGCNN-style model. "
        "This is the core of the chapter's pipeline; the full study trains the "
        "same architecture on thousands of Materials Project structures.",
    ),
    md("## Install\n\n"
       "`torch-geometric` provides the graph data structures and "
       "message-passing primitives; `pymatgen` builds and analyses the crystal "
       "structures. `torch` comes pre-installed on Colab. A GPU is optional for "
       "a single forward pass but recommended once you start training."),
    code("!pip install torch-geometric pymatgen\n"),
    md("## Build a crystal and convert it to a graph\n\n"
       "We take an SrTiO3 perovskite, find every neighbour within a cutoff "
       "radius with pymatgen, and read off the directed edge list and "
       "edge distances — the raw ingredients of a crystal graph."),
    code("import numpy as np\n"
         "import torch\n"
         "from pymatgen.core import Structure, Lattice\n\n"
         "lattice = Lattice.cubic(3.905)\n"
         "structure = Structure(\n"
         "    lattice,\n"
         "    ['Sr', 'Ti', 'O', 'O', 'O'],\n"
         "    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5],\n"
         "     [0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]],\n"
         ")\n\n"
         "cutoff = 5.0\n"
         "centres, points, _, distances = structure.get_neighbor_list(\n"
         "    r=cutoff, exclude_self=True)\n"
         "Z = np.array(structure.atomic_numbers, dtype=np.int64)\n"
         "edge_index = np.stack([centres, points], axis=0)\n"
         "print(f'{len(Z)} atoms, {edge_index.shape[1]} directed edges')\n"),
    md("## Expand distances in a Gaussian basis\n\n"
       "Raw scalar distances make poor neural-network inputs. CGCNN expands "
       "each edge distance in a fixed bank of Gaussians, turning one number "
       "into a smooth feature vector."),
    code("import torch.nn as nn\n\n"
         "class GaussianBasis(nn.Module):\n"
         "    def __init__(self, r_min=0.0, r_max=8.0, n_basis=64):\n"
         "        super().__init__()\n"
         "        centres = torch.linspace(r_min, r_max, n_basis)\n"
         "        self.register_buffer('centres', centres)\n"
         "        self.sigma = float(centres[1] - centres[0])\n\n"
         "    def forward(self, r):\n"
         "        delta = r.unsqueeze(-1) - self.centres\n"
         "        return torch.exp(-0.5 * (delta / self.sigma) ** 2)\n\n"
         "basis = GaussianBasis(0.0, 8.0, 64)\n"
         "r = torch.tensor(distances, dtype=torch.float32)\n"
         "edge_attr = basis(r)\n"
         "print('edge feature tensor:', edge_attr.shape)\n"),
    md("## Wrap it in a PyTorch Geometric `Data` object\n\n"
       "`Data` is the standard container: node features, edge index, edge "
       "features, and (in training) a target. A `batch` vector of zeros marks "
       "every atom as belonging to graph 0."),
    code("from torch_geometric.data import Data\n\n"
         "graph = Data(\n"
         "    Z=torch.tensor(Z, dtype=torch.long),\n"
         "    edge_index=torch.tensor(edge_index, dtype=torch.long),\n"
         "    edge_attr=edge_attr,\n"
         "    batch=torch.zeros(len(Z), dtype=torch.long),\n"
         ")\n"
         "print(graph)\n"),
    md("## Define a small CGCNN-style model\n\n"
       "A faithful, compact re-implementation of the Crystal Graph "
       "Convolutional Neural Network (Xie & Grossman, 2018): an element "
       "embedding, a stack of gated message-passing layers, a global pool, and "
       "a small read-out head that produces one scalar per crystal."),
    code("import torch.nn.functional as F\n"
         "from torch_geometric.nn import MessagePassing, global_mean_pool\n\n"
         "class CGCNNConv(MessagePassing):\n"
         "    def __init__(self, atom_dim, edge_dim):\n"
         "        super().__init__(aggr='add')\n"
         "        z_dim = 2 * atom_dim + edge_dim\n"
         "        self.gate_linear = nn.Linear(z_dim, atom_dim)\n"
         "        self.core_linear = nn.Linear(z_dim, atom_dim)\n"
         "        self.bn_msg = nn.BatchNorm1d(atom_dim)\n"
         "        self.bn_out = nn.BatchNorm1d(atom_dim)\n\n"
         "    def forward(self, h, edge_index, e):\n"
         "        agg = self.propagate(edge_index, h=h, e=e)\n"
         "        return self.bn_out(h + agg)\n\n"
         "    def message(self, h_i, h_j, e):\n"
         "        z = torch.cat([h_i, h_j, e], dim=-1)\n"
         "        gate = torch.sigmoid(self.gate_linear(z))\n"
         "        core = F.softplus(self.core_linear(z))\n"
         "        return self.bn_msg(gate * core)\n\n"
         "class CGCNN(nn.Module):\n"
         "    def __init__(self, n_elements=100, atom_dim=64, edge_dim=64,\n"
         "                 n_conv=3, hidden_dim=128, n_targets=1):\n"
         "        super().__init__()\n"
         "        self.embedding = nn.Embedding(n_elements, atom_dim)\n"
         "        self.convs = nn.ModuleList(\n"
         "            [CGCNNConv(atom_dim, edge_dim) for _ in range(n_conv)])\n"
         "        self.head = nn.Sequential(\n"
         "            nn.Linear(atom_dim, hidden_dim), nn.Softplus(),\n"
         "            nn.Linear(hidden_dim, n_targets))\n\n"
         "    def forward(self, data):\n"
         "        h = self.embedding(data.Z)\n"
         "        for conv in self.convs:\n"
         "            h = conv(h, data.edge_index, data.edge_attr)\n"
         "        h_G = global_mean_pool(h, data.batch)\n"
         "        return self.head(h_G).squeeze(-1)\n\n"
         "model = CGCNN(edge_dim=edge_attr.shape[1])\n"
         "print(model)\n"),
    md("## Run one forward pass\n\n"
       "With an untrained model the number is meaningless — but a clean forward "
       "pass confirms the graph, the basis expansion and the message-passing "
       "layers all fit together. Training this model on Materials Project data "
       "is the subject of the chapter's pipeline section."),
    code("model.eval()\n"
         "with torch.no_grad():\n"
         "    prediction = model(graph)\n"
         "print('predicted property (untrained):', float(prediction))\n"),
]


# ---------------------------------------------------------------------------
# Chapter 12 — Foundation Models for Materials
# ---------------------------------------------------------------------------
CH12 = [
    header(
        "Chapter 12 · Foundation Models for Materials",
        "ch12-foundation",
        "We load MACE-MP-0, a pre-trained universal interatomic potential, and "
        "use it zero-shot: no training, no fine-tuning. We evaluate energies "
        "and forces on a crystal it has never explicitly seen, then run a short "
        "molecular-dynamics trajectory — all from a downloaded checkpoint.",
    ),
    md("## Install\n\n"
       "`mace-torch` provides both the MACE library and the `mace_mp` "
       "foundation-model loader, which downloads the pre-trained MACE-MP-0 "
       "checkpoint on first use. `torch` comes pre-installed on Colab. "
       "**Set the runtime to GPU** (Runtime -> Change runtime type -> GPU) for "
       "the molecular-dynamics cell; the single-point evaluation is fine on "
       "CPU."),
    code("!pip install mace-torch\n"),
    md("## Load the MACE-MP-0 foundation model\n\n"
       "`mace_mp` is the foundation-model loader. The first call downloads the "
       "checkpoint (the `medium` model is a good default) and returns an ASE "
       "calculator ready to attach to any structure — no training step."),
    code("import torch\n"
         "from mace.calculators import mace_mp\n\n"
         "device = 'cuda' if torch.cuda.is_available() else 'cpu'\n"
         "calc = mace_mp(model='medium', dispersion=False,\n"
         "               default_dtype='float64', device=device)\n"
         "print('MACE-MP-0 loaded on', device)\n"),
    md("## Build a structure and run a zero-shot evaluation\n\n"
       "We construct an MgO rock-salt supercell and attach the foundation "
       "model. The energy and forces come straight out — zero-shot — with no "
       "system-specific training. This is the headline capability of a "
       "materials foundation model."),
    code("import numpy as np\n"
         "from ase.build import bulk\n\n"
         "mgo = bulk('MgO', crystalstructure='rocksalt', a=4.21).repeat((3, 3, 3))\n"
         "mgo.calc = calc\n\n"
         "energy = mgo.get_potential_energy()\n"
         "forces = mgo.get_forces()\n"
         "print(f'{len(mgo)} atoms')\n"
         "print(f'total energy   = {energy:.4f} eV')\n"
         "print(f'energy / atom  = {energy / len(mgo):.4f} eV')\n"
         "print(f'max |force|    = {np.abs(forces).max():.4e} eV/A  '\n"
         "      '(near zero: the ideal lattice is close to equilibrium)')\n"),
    md("## Probe a rattled structure\n\n"
       "Displace the atoms slightly and the forces become non-trivial — the "
       "foundation model returns a restoring force pushing every atom back "
       "towards its lattice site, exactly as a DFT calculation would."),
    code("rng = np.random.default_rng(0)\n"
         "rattled = mgo.copy()\n"
         "rattled.rattle(stdev=0.05, rng=rng)\n"
         "rattled.calc = calc\n"
         "print(f'rattled energy    = {rattled.get_potential_energy():.4f} eV')\n"
         "print(f'rattled max|force| = {np.abs(rattled.get_forces()).max():.4f} eV/A')\n"),
    md("## A short zero-shot molecular-dynamics run\n\n"
       "Because the foundation model returns forces, we can drive molecular "
       "dynamics with it directly. Here is a brief NVT Langevin trajectory at "
       "600 K — short enough to finish quickly, long enough to show the "
       "temperature settling around its target. Use the GPU runtime for this "
       "cell."),
    code("from ase.md.langevin import Langevin\n"
         "from ase.md.velocitydistribution import MaxwellBoltzmannDistribution\n"
         "from ase.units import fs\n\n"
         "md_atoms = mgo.copy()\n"
         "md_atoms.calc = calc\n"
         "MaxwellBoltzmannDistribution(md_atoms, temperature_K=600.0)\n"
         "dyn = Langevin(md_atoms, timestep=1.0 * fs, temperature_K=600.0,\n"
         "               friction=0.01)\n\n"
         "energies, temperatures = [], []\n"
         "def record():\n"
         "    energies.append(float(md_atoms.get_potential_energy()))\n"
         "    temperatures.append(float(md_atoms.get_temperature()))\n"
         "dyn.attach(record, interval=1)\n"
         "dyn.run(100)\n"
         "print(f'final energy = {energies[-1]:.3f} eV')\n"
         "print(f'mean T over last 50 steps = {np.mean(temperatures[-50:]):.1f} K')\n"),
    md("## Plot the trajectory\n\n"
       "Energy and temperature versus step. The thermostat draws the "
       "temperature towards 600 K while the potential energy responds — the "
       "expected signature of an equilibrating NVT run."),
    code("import matplotlib.pyplot as plt\n\n"
         "fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))\n"
         "ax1.plot(energies)\n"
         "ax1.set_xlabel('MD step')\n"
         "ax1.set_ylabel('potential energy (eV)')\n"
         "ax1.set_title('Energy')\n"
         "ax2.plot(temperatures)\n"
         "ax2.axhline(600.0, ls='--', color='grey', label='target 600 K')\n"
         "ax2.set_xlabel('MD step')\n"
         "ax2.set_ylabel('temperature (K)')\n"
         "ax2.set_title('Temperature')\n"
         "ax2.legend()\n"
         "fig.tight_layout()\n"
         "plt.show()\n"),
]


NOTEBOOK_SPECS = {
    "ch06-running-dft-colab.ipynb": CH06,
    "ch09-mlip-colab.ipynb": CH09,
    "ch10-gnn-colab.ipynb": CH10,
    "ch12-foundation-colab.ipynb": CH12,
}


def main() -> None:
    for name, cells in NOTEBOOK_SPECS.items():
        path = write_notebook(name, cells)
        print(f"wrote {path.relative_to(REPO_ROOT)}  ({len(cells)} cells)")


if __name__ == "__main__":
    main()
