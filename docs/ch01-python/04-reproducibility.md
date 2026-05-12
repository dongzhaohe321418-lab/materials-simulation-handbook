# 1.4 Reproducibility, environments, and version control

If a result cannot be reproduced, it is not a result. This sounds obvious until you try to rerun, eighteen months later, the code that produced Figure 3 of a paper you wrote, and discover that the script depends on a package version that no longer exists, in an environment you cannot reconstruct, on a laptop that has since been wiped.

This section is the cheapest insurance policy you will ever buy: a small set of habits that take five minutes a week and make the difference between work that compounds and work that evaporates.

## Why reproducibility matters

The reproducibility crisis is real and not confined to social science. In computational materials, several high-profile cases have shown how fragile results can be:

- A 2016 study attempted to rerun 176 published DFT calculations across 15 codes. Energy differences for the *same* material and the *same* exchange-correlation functional varied by more than 1 meV per atom across codes, and by much more when the codes were used with their default settings.
- A widely cited machine-learning-for-materials benchmark turned out, two years after publication, to have leaked test-set composition into training. Once corrected, several "state-of-the-art" models lost most of their advantage.
- A small but recurring fraction of *Nature* materials papers — single-digit percent in surveys — cannot be rerun even by the original authors, usually because of lost code, lost environments, or lost data.

You will not avoid all of these. But the version of your future self that has to fix a Reviewer 2 comment three months after submission will be grateful that *you* can rerun your own code.

## The four pillars

In rough order of importance:

1. **The exact code** — every script, configuration file, and notebook used to produce the figure.
2. **The exact environment** — the Python version, the package versions, the BLAS implementation, the CUDA driver.
3. **The exact inputs** — the structure files, the DFT input files, the data splits.
4. **The exact random state** — every seed, in every library that has one.

Lose any one and your work is in trouble. We address them in turn.

## Pinning versions

A version pin says: "this code was tested with exactly this version of this package." Pins live in `requirements.txt` or `environment.yml`.

Two kinds of pins:

```
# Minimum, with a permissive upper bound — for active development
numpy>=1.26,<2.0

# Exact — for archival snapshots
numpy==1.26.4
```

During development, use loose pins (`>=`) so you can upgrade. When you submit a paper, *freeze* the environment with exact pins.

### Freezing with pip

```bash
pip freeze > requirements-frozen.txt
```

The output looks like:

```
ase==3.22.1
matplotlib==3.8.2
numpy==1.26.4
pymatgen==2024.2.20
scipy==1.11.4
torch==2.2.0
```

Commit this file alongside the paper-producing scripts.

### Freezing with conda

```bash
conda env export --no-builds > environment-frozen.yml
```

`--no-builds` strips the per-platform build hashes so the file is at least *recreatable* on a different OS. For full hash-level reproducibility, omit the flag — but accept that the file is then only usable on the OS that produced it.

### Lockfiles

The next level of rigour is a *lockfile* — a pinned, hashed, fully resolved dependency graph. Two good tools:

- **`pip-tools`** — `pip-compile requirements.in` produces `requirements.txt` with all transitive dependencies pinned.
- **`conda-lock`** — `conda-lock lock -f environment.yml` produces per-platform lockfiles.

For a thesis or a long-term project, lockfiles are worth the effort. For an exploratory notebook, a frozen `requirements.txt` is enough.

## Git basics for scientists

Git is the de facto standard for source control. You will use a fraction of its features, and that fraction is enough.

### Initialise a repository

```bash
mkdir my-project
cd my-project
git init
git branch -M main
```

### The basic loop

You edit files, then *stage* and *commit* the changes:

```bash
git status                       # what has changed?
git add my_script.py             # stage one file
git add .                        # stage everything in the current dir
git commit -m "compute RDF for Si supercell"
```

Commit messages should be short and informative. "fix" is a bad message; "fix RDF normalisation: divide by N(N-1)/2 not N²" is a good one.

### Branches

A branch is a parallel line of development. When you start a risky change, create one:

```bash
git checkout -b feature/mace-potential
# ... make changes, commit ...
git checkout main
git merge feature/mace-potential
```

For small projects you may live entirely on `main`. Branches matter most when collaborating or experimenting with something you might throw away.

### Pushing to a remote

A remote is a hosted copy on GitHub, GitLab, or your group's server.

```bash
git remote add origin git@github.com:your-username/my-project.git
git push -u origin main
```

After the first push, `git push` on its own is enough.

!!! tip "Push the day you start"
    Push the empty repository, then push after every working session. A bricked laptop is annoying. A bricked laptop containing your only copy of a thesis chapter is a tragedy.

### `.gitignore` for scientific projects

You do **not** want to commit large output files, build artefacts, or environment caches. Save the following as `.gitignore`:

```gitignore
# Python
__pycache__/
*.py[cod]
*.egg-info/
.ipynb_checkpoints/
.venv/
.env

# Editors
.vscode/
.idea/
*.swp

# Operating system
.DS_Store
Thumbs.db

# Scientific outputs — keep code, drop large artefacts
*.cube
*.wfn
*.chk
*.h5
*.hdf5
*.traj
*.xyz.gz
*.npy
*.npz
*.pdf
*.png
data/raw/
data/processed/
results/
runs/
wandb/
mlruns/

# Compiled DFT outputs
OUTCAR
CHGCAR
WAVECAR
vasprun.xml
```

Adjust to taste. Anything that takes seconds to regenerate from code should not be in git; anything that takes a week of compute should be archived elsewhere (a Zenodo deposit, an institutional store), with the *recipe* committed to git.

## Structuring a scientific repository

A layout that scales from one student to a small group:

```
my-project/
├── README.md                # one paragraph: what, why, how to run
├── environment.yml          # development pins
├── environment-frozen.yml   # archival pins
├── .gitignore
├── src/
│   └── myproject/
│       ├── __init__.py
│       ├── data.py
│       ├── models.py
│       └── analysis.py
├── scripts/
│   ├── 01_prepare_structures.py
│   ├── 02_run_dft.py
│   └── 03_make_figures.py
├── notebooks/
│   └── exploratory.ipynb
├── tests/
│   └── test_data.py
├── configs/
│   └── default.yaml
├── data/                    # gitignored except for README and pointer files
└── results/                 # gitignored
```

Three rules:

- **Library in `src/`, runners in `scripts/`.** The library has no top-level side effects and is importable. The scripts call into the library and contain `if __name__ == "__main__":` blocks. Notebooks may import the library but should never define the same functions twice.
- **Scripts numbered in pipeline order.** Anyone reading the repository can see the order of operations at a glance.
- **Data and results are *not* in git.** Add a `data/README.md` that explains where the data came from and how to obtain it.

## Random seeds

Nearly every paper in machine-learning materials science contains the word "random" — initialisation, train/test splits, data shuffling, MD initial velocities. Set every seed.

```python
import os
import random
import numpy as np
import torch

def set_seed(seed: int = 0) -> None:
    """Make NumPy, Python, and PyTorch deterministic."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    # Deterministic CUDA at a small performance cost
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
```

A subtlety: NumPy now prefers a *generator object* over the global state. The modern, recommended idiom:

```python
rng = np.random.default_rng(seed=42)
x = rng.normal(size=100)
y = rng.choice(structures, size=10, replace=False)
```

Pass `rng` into every function that needs randomness, rather than reseeding globally. This makes parallel code and library code much easier to reason about.

!!! warning "Determinism is not free"
    Even with all seeds set, GPU operations are not bitwise deterministic across hardware or driver versions. You can reproduce the same *result* on the same hardware; you may not reproduce identical floating-point numbers across machines.

## Logging software versions in output files

Every output file should record what produced it. The smallest useful header:

```python
import sys
import platform
import json
from datetime import datetime, timezone
from importlib.metadata import version, PackageNotFoundError

def provenance() -> dict:
    pkgs = ["numpy", "scipy", "matplotlib", "ase", "pymatgen", "torch", "mace-torch"]
    out: dict[str, str] = {}
    for p in pkgs:
        try:
            out[p] = version(p)
        except PackageNotFoundError:
            out[p] = "not installed"
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": out,
    }

with open("results/run_001/provenance.json", "w") as f:
    json.dump(provenance(), f, indent=2)
```

A `provenance.json` next to every results directory will save you in nine out of ten "what did I run?" emergencies. For the tenth, commit the git SHA:

```python
import subprocess

def git_sha() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
```

Then your results carry both the environment *and* the exact code version that produced them.

## Jupyter notebooks versus scripts

Notebooks are excellent for exploration: tight feedback loops, rich output, prose alongside code. They are *bad* at:

- diffing in version control (the JSON format is noisy);
- enforcing execution order (cells can be run out of order);
- running unattended (no command-line interface);
- importing into other code (functions defined in a notebook are awkward to reuse).

A simple discipline:

- **Notebooks for exploration and figures.** Anything you might throw away. Always restart-and-run-all before you commit, so the saved state matches the executed state. Consider stripping outputs with `nbstripout` before committing — it makes diffs sane.
- **Scripts for everything else.** As soon as a notebook does something you will repeat, lift its functions into a module under `src/` and import them. The notebook becomes a thin orchestrator; the logic is testable and reusable.

Tools that help:

- `jupytext` — pairs a notebook with a plain `.py` file kept in sync, so diffs are readable.
- `nbstripout` — strips notebook outputs on commit. Install once, forget about it.
- `papermill` — runs a notebook with parameters from the command line, for parameter sweeps.

## FAIR principles, in one paragraph

The FAIR principles say that scientific data — and, by extension, code — should be **F**indable, **A**ccessible, **I**nteroperable, and **R**eusable. Concretely: deposit your code at a persistent identifier (a Zenodo DOI for a paper-version snapshot, a GitHub repo for the living version); deposit your structures and trained models with an open licence; use standard formats (CIF for crystals, extended XYZ for trajectories, ONNX or `torch.save` for models); and document inputs and outputs so that a reader who is not you can rerun the pipeline. None of this is hard. All of it is the difference between a paper that has impact and a paper that has only citations.

## A minimum starter checklist

Before you run a calculation worth a coffee or more:

- [ ] Project is in a git repository, pushed to a remote.
- [ ] `environment.yml` (or `requirements.txt`) is up to date and committed.
- [ ] `.gitignore` excludes large outputs.
- [ ] All random seeds are set in one place.
- [ ] A `provenance.json` is written next to every results directory.
- [ ] A `README.md` describes how to run the pipeline end to end.

Five minutes. Worth a thesis.
