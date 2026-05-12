---
title: Try in Jupyter
---

# Try in Jupyter (in your browser)

The handbook ships with a full [JupyterLab](https://jupyter.org/) environment
that runs entirely in your browser via the
[Pyodide](https://pyodide.org/) WebAssembly distribution of CPython. No
installation, no account, no server: everything executes on your machine
inside the browser tab.

[Launch JupyterLab](/materials-simulation-handbook/lite/lab/index.html){ .md-button .md-button--primary }
[Open the file browser](/materials-simulation-handbook/lite/index.html){ .md-button }

## What is preloaded

- Python 3.11 (Pyodide build)
- NumPy, SciPy, Matplotlib, ipywidgets
- One ready-to-run notebook per substantial chapter:

| Chapter | Notebook |
| --- | --- |
| Chapter 0 · Mathematics from Scratch | [ch00-math.ipynb](/materials-simulation-handbook/lite/lab/index.html?path=ch00-math.ipynb) |
| Chapter 1 · Python and Scientific Computing | [ch01-python.ipynb](/materials-simulation-handbook/lite/lab/index.html?path=ch01-python.ipynb) |
| Chapter 3.5 · Solid State Prerequisites | [ch03b-solid-state.ipynb](/materials-simulation-handbook/lite/lab/index.html?path=ch03b-solid-state.ipynb) |
| Chapter 4 · Quantum Mechanics | [ch04-quantum.ipynb](/materials-simulation-handbook/lite/lab/index.html?path=ch04-quantum.ipynb) |
| Chapter 5 · Density Functional Theory | [ch05-dft.ipynb](/materials-simulation-handbook/lite/lab/index.html?path=ch05-dft.ipynb) |
| Chapter 7 · Molecular Dynamics | [ch07-md.ipynb](/materials-simulation-handbook/lite/lab/index.html?path=ch07-md.ipynb) |
| Chapter 8 · Statistical Mechanics | [ch08-statmech.ipynb](/materials-simulation-handbook/lite/lab/index.html?path=ch08-statmech.ipynb) |
| Chapter 11 · Active Learning and Bayesian Optimisation | [ch11-active.ipynb](/materials-simulation-handbook/lite/lab/index.html?path=ch11-active.ipynb) |

## What runs and what does not

JupyterLite uses the Pyodide kernel, which is a port of CPython to
WebAssembly. The pure-Python and Fortran-free scientific stack works well
(NumPy, SciPy, Matplotlib, scikit-learn, pandas, sympy, networkx,
ipywidgets). Packages that wrap large native compiled codes do **not** run
in the browser:

- ASE, pymatgen, MACE, NequIP, PyTorch, TensorFlow, JAX
- LAMMPS, GROMACS, Quantum ESPRESSO, VASP
- anything that requires a GPU

For these chapters (6 on running DFT, 9 on MLIPs, 10 on GNNs, 12 on
foundation models) we keep the "Open in Colab" route, where a free Linux
VM with the full stack is one click away.

## Privacy and persistence

Every notebook you edit is saved to your browser's local storage. Nothing
leaves your machine. Clearing your browser data will wipe the notebooks;
download a copy from the JupyterLab file menu first if you want to keep
them.

## Reporting issues

If a notebook fails to load or a cell crashes the kernel, please open an
issue at the [GitHub repository](https://github.com/dongzhaohe321418-lab/materials-simulation-handbook/issues)
with the browser name, version and a screenshot of the JupyterLab
console.
