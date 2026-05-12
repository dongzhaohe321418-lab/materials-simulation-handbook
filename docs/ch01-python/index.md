# Chapter 1 — Python and Scientific Computing

[Open in Jupyter (browser)](/materials-simulation-handbook/lite/lab/index.html?path=ch01-python.ipynb){ .md-button .md-button--primary }

Before we can simulate a single atom, we need a working scientific Python environment and enough fluency in it to read, modify, and write the code that the rest of this book is built on. This chapter assumes **no previous programming experience**. If you have written Python before, skim the headings and skip to the exercises at the end; you will lose nothing by checking your fundamentals.

## Why Python?

There are faster languages. Fortran will out-run a careful Python loop by two orders of magnitude, and C++ underpins most production DFT codes. Yet almost every paper in machine-learning interatomic potentials, every modern materials informatics pipeline, and every active workflow engine — ASE, pymatgen, FHI-aims-tools, AiiDA, PyTorch Geometric, MACE, JAX-MD — is written in Python, or exposes a Python interface as its primary user surface.

Python won the scientific computing wars for three reasons:

1. **Readability**. A Python script reads almost like pseudocode, which lowers the barrier between an idea on a whiteboard and a runnable experiment.
2. **An ecosystem of fast libraries**. The slow parts (linear algebra, FFTs, sparse solvers, neural network kernels) are written in C, Fortran, or CUDA and merely *called* from Python. We get C-speed on the hot path and Python-comfort everywhere else.
3. **Glue**. Python is the connective tissue between data, simulation engines, and analysis. ASE talks to thirty different DFT codes through a uniform interface; pymatgen ingests CIFs and writes VASP inputs; PyTorch lets you train a potential and use it inside an MD run in the same script.

You will spend the rest of your career writing Python whether you want to or not. Better to do it well.

## What we'll cover

This chapter is organised into five short sections and an exercise set:

- **[Installation](01-install.md)** — Setting up a clean, reproducible Python environment with `miniconda` or `mamba`. We avoid the pitfalls (system Python, broken `pip`, mysterious permission errors) that consume the first week of every graduate student's life.
- **[NumPy](02-numpy.md)** — The single most important library in scientific Python. We treat arrays, broadcasting, indexing, reductions, and basic linear algebra at a level sufficient for the rest of the book.
- **[Matplotlib](03-matplotlib.md)** — Publication-quality plotting. The figure/axes mental model, subplots, log axes, twin axes, saving PDFs for a thesis.
- **[Reproducibility](04-reproducibility.md)** — Version control, environment locking, random seeds, and the difference between code that *ran once* and code that *can be run again*. This is the section that separates a working scientist from a frustrated one.
- **[Exercises](exercises.md)** — Eight problems with worked solutions, ranging from one-liners to a small modelling task.

## Tools we will use across the book

| Library | Role |
|---|---|
| **NumPy** | $n$-dimensional arrays; the foundation everything else sits on |
| **SciPy** | Optimisation, integration, sparse linear algebra, special functions |
| **Matplotlib** | Plotting |
| **ASE** | Atomic structures, calculators (DFT, MD), file I/O |
| **pymatgen** | Crystallography, symmetry, Materials Project API |
| **PyTorch** | Neural networks, autograd, GPU |
| **MACE / NequIP** | Equivariant message-passing potentials (Chapters 9–10) |

By the end of Chapter 1 you will be able to install these, write a script that loads a crystal structure and computes a property, and commit it to a git repository in a way that someone else — including your future self — can actually reproduce.

## The return on investment

Every hour spent on Python fundamentals saves a day later. The most expensive bug in a research codebase is not a wrong formula; it is a silent type coercion, a stale environment, a missing seed, a misaligned axis. The habits we install in this chapter — small, tested functions; pinned environments; version-controlled notebooks — are not optional polish. They are the difference between a paper that survives peer review and a paper that quietly retracts itself two years later.

Take it slowly. Type every example. Break them on purpose and read the traceback. Then turn the page.
