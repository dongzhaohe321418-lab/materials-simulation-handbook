# Code

Standalone, runnable reference implementations associated with the handbook.

## Where the code actually lives (for now)

The vast majority of teaching code in this handbook is **embedded inline** in the chapter markdown under `docs/`. Each chapter section that introduces an algorithm follows it with a complete, runnable Python listing in a fenced code block. To extract a listing, copy it out of the rendered page or out of the corresponding `docs/chXX-*/0Y-section.md` file.

Figure-generating scripts live under [`scripts/figures/`](../scripts/figures/) — these are the matplotlib programs that produce the PNGs in [`docs/assets/figures/`](../docs/assets/figures/) and can be re-run to regenerate or modify the figures.

The five capstone project templates each ship their own `starter-code.md` with a full scaffold; see [`docs/projects/`](../docs/projects/).

## Roadmap for this directory

As the handbook matures, the heavier worked examples — the from-scratch SCF loop, the CGCNN implementation, the MACE-on-water training pipeline — will be lifted out of the inline listings into stand-alone, pip-installable Python modules under this directory, with proper tests and CI. Pull requests adding such modules are welcome; see [CONTRIBUTING.md](../CONTRIBUTING.md).
