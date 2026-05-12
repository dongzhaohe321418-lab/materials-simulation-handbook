# 1.1 Installing Python the right way

Almost every Python disaster in a research group can be traced to one of three root causes: using the operating system's bundled Python, installing the wrong distribution, or mixing `pip` and `conda` carelessly. This section shows you the path that avoids all three.

## Why not the Python that came with my computer?

macOS and most Linux distributions ship with a system Python. Do **not** use it for science. The system Python:

- is tied to OS updates — upgrading macOS may silently change your interpreter;
- often lacks the headers needed to build scientific packages;
- requires `sudo` to install anything globally, which corrupts the system if you ever pin a different version of a shared library;
- is missing important compiled dependencies (BLAS, LAPACK, FFTW) or links against slow defaults.

The system Python exists to run system utilities, not your DFT pipeline.

## Why not Anaconda?

Anaconda — the full distribution — bundles hundreds of packages you do not need and several you might not want (licensed, telemetry-enabled defaults, an enormous installer). Since 2020, Anaconda's commercial terms have also restricted its use in organisations with more than 200 employees.

We use **Miniconda** (or its drop-in faster sibling **Mamba**, distributed as Miniforge), which gives you:

- the `conda` package manager,
- a small base environment,
- access to the `conda-forge` channel (community-maintained, well-curated, license-clean),

and nothing else. You install only what you need.

!!! tip "Recommendation"
    Install **Miniforge3**. It is Miniconda preconfigured with `conda-forge` as the default channel and `mamba` as a faster solver. One installer, no licensing concerns, no surprises.

## Step-by-step installation

### macOS (Apple Silicon and Intel)

Open Terminal and run:

```bash
# Pick the right arch — uname -m prints arm64 on Apple Silicon, x86_64 on Intel
ARCH=$(uname -m)
curl -L -o miniforge.sh \
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-${ARCH}.sh"
bash miniforge.sh -b -p "$HOME/miniforge3"
rm miniforge.sh
```

The `-b` flag accepts the licence non-interactively; `-p` chooses the install prefix. Now initialise your shell:

```bash
"$HOME/miniforge3/bin/conda" init "$(basename "$SHELL")"
exec "$SHELL"
```

The second line restarts your shell so the new `conda` is on your `PATH`. Verify:

```bash
which conda
conda --version
python --version
```

You should see paths under `~/miniforge3` and a Python version of 3.11 or newer.

### Linux (x86_64 and aarch64)

```bash
ARCH=$(uname -m)
wget -O miniforge.sh \
  "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh"
bash miniforge.sh -b -p "$HOME/miniforge3"
rm miniforge.sh
"$HOME/miniforge3/bin/conda" init bash
exec bash
```

If you are on a cluster, install into `$HOME` rather than a shared directory. Conda environments contain symlinks that break if moved, so installing to your own home directory is the path of least resistance.

### Windows

Use the Windows Subsystem for Linux (WSL2) with Ubuntu, then follow the Linux instructions above. Native Windows Python works but is a second-class citizen for materials simulation — most DFT and MD codes assume a POSIX environment.

## What's an environment, and why have many?

An **environment** is an isolated directory containing a specific Python interpreter and a specific set of packages. Different projects can live in different environments without interfering. If you upgrade `numpy` in one and it breaks an old script in another, the old script does not care, because its environment still has the old `numpy`.

You should have:

- a *base* environment that you never touch except to update `conda` itself;
- one environment per project (or per closely related cluster of projects).

!!! warning "Never install into base"
    The moment you `pip install` something into your `base` environment, you have lost the ability to reset cleanly. Treat `base` as read-only.

## Creating the `materials-simulation` environment

Save the following as `environment.yml` in the root of a new project directory:

```yaml
name: materials-simulation
channels:
  - conda-forge
  - pytorch
dependencies:
  - python=3.11
  - numpy>=1.26
  - scipy>=1.11
  - matplotlib>=3.8
  - jupyterlab
  - ipykernel
  - ase>=3.22
  - pymatgen>=2024.1
  - pytorch>=2.2
  - cpuonly                # remove on a GPU machine, add pytorch-cuda=12.1
  - pandas
  - h5py
  - pyyaml
  - tqdm
  - pip
  - pip:
      - mace-torch>=0.3.6
      - nequip
```

Create the environment:

```bash
mamba env create -f environment.yml
conda activate materials-simulation
```

If you do not have `mamba`, swap it for `conda` — the syntax is identical, only slower. The first solve may take a few minutes; `conda-forge` is a large channel.

### Activating and deactivating

```bash
conda activate materials-simulation  # turn it on
python -c "import numpy; print(numpy.__version__)"
conda deactivate                     # turn it off
```

Your shell prompt should show `(materials-simulation)` while the environment is active. If it does not, your shell init step earlier did not take effect; rerun `conda init` and open a fresh terminal.

## `requirements.txt` vs `environment.yml`

You will see both formats in the wild. They are not interchangeable.

- **`requirements.txt`** is the pure-pip format. One package per line, optional version pins (`numpy==1.26.2`). It records only Python packages.
- **`environment.yml`** is conda's format. It records the **Python version**, **non-Python dependencies** (BLAS, MPI, CUDA), and channels. It can embed a `pip:` section for packages not on conda-forge.

For materials simulation, `environment.yml` is the right default because we depend on compiled native libraries (BLAS for NumPy, libxc for DFT, CUDA for PyTorch). A bare `requirements.txt` cannot specify these.

A pragmatic rule:

- Use `environment.yml` to *recreate* a project environment from scratch.
- Use `pip freeze > requirements.txt` to *snapshot* exact pinned versions for a finished paper.

We come back to this in [Section 1.4](04-reproducibility.md).

## Verifying the install

With the environment active, run:

```bash
python - <<'EOF'
import numpy, scipy, matplotlib, ase, pymatgen, torch
print(f"numpy      {numpy.__version__}")
print(f"scipy      {scipy.__version__}")
print(f"matplotlib {matplotlib.__version__}")
print(f"ase        {ase.__version__}")
print(f"pymatgen   {pymatgen.__version__}")
print(f"torch      {torch.__version__}  CUDA={torch.cuda.is_available()}")
EOF
```

If all lines print without a traceback, your install is healthy. A quick functional test:

```python
# verify_install.py
import numpy as np
from ase.build import bulk

si = bulk("Si", cubic=True)
print(si)                     # 8-atom Si conventional cell
print("Volume:", si.get_volume(), "Å³")
print("First-NN distance:", si.get_all_distances(mic=True)[0, 1], "Å")

# A trivial NumPy check
A = np.random.default_rng(0).normal(size=(3, 3))
A = A + A.T                    # symmetric
eigs = np.linalg.eigvalsh(A)
print("Eigenvalues:", eigs)
```

Run with `python verify_install.py`. You should see a silicon cell, a volume near 160 Å³, and three real eigenvalues.

## Common pitfalls

### "I have three Pythons. Which one is running?"

Check with:

```bash
which -a python python3
type python
```

If `which -a` lists multiple, the *first* on your `PATH` is what runs. When `conda activate` is active, this should be the one inside `~/miniforge3/envs/<your-env>/bin/`.

### `pip install` ignored my environment

If `pip` lives outside your environment, it installs to the wrong place. Always check:

```bash
which pip
pip --version
```

The path should be inside the active environment. If it is not, install `pip` *into* the environment with `conda install pip`, then close and reopen your terminal.

### "It worked yesterday and now it doesn't"

Three usual suspects:

1. You ran `pip install <something>` into `base` and it overrode a conda-managed library. Recreate the environment from `environment.yml`.
2. A system update changed your shell init file. Re-run `conda init`.
3. Your environment file was edited and re-solved; the new solve picked different versions. This is exactly why we pin versions before publishing — see [Section 1.4](04-reproducibility.md).

### Permission errors

If you see `Permission denied` or are tempted to `sudo pip install`, **stop**. You are about to corrupt your system Python. Activate your conda environment first; `pip` inside an environment never needs `sudo`.

### Mixing `pip` and `conda`

`conda` and `pip` do not know about each other's metadata. A safe rule: install everything with `conda`/`mamba` first, then *only* use `pip` for packages unavailable on conda-forge (e.g., `mace-torch`). Never use `pip` to upgrade something `conda` installed.

### Apple Silicon and Intel-only wheels

A few packages still ship only `x86_64` wheels. If you hit `no matching distribution`, try `conda-forge` first (which usually has an `osx-arm64` build); if not, run inside Rosetta with:

```bash
CONDA_SUBDIR=osx-64 mamba env create -f environment.yml
```

This forces an Intel environment on Apple Silicon; macOS will run it under Rosetta translation.

## Beyond the basics: lockfiles

For a project that absolutely must reproduce, generate a lockfile:

```bash
conda env export --no-builds > environment.lock.yml
```

This records every transitively-resolved package version. Commit it. To recreate the *exact* environment six months later:

```bash
mamba env create -f environment.lock.yml -n materials-simulation-snapshot
```

We will revisit lockfiles, hashes, and full reproducibility in [Section 1.4](04-reproducibility.md). For now, you have a working scientific Python install. Time to use it.
