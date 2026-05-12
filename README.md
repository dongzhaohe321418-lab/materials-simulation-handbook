# Materials Simulation Handbook

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Text: CC BY 4.0](https://img.shields.io/badge/Text-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Deploy](https://github.com/dongzhaohe321418-lab/materials-simulation-handbook/actions/workflows/deploy.yml/badge.svg)](https://github.com/dongzhaohe321418-lab/materials-simulation-handbook/actions/workflows/deploy.yml)
[![Site](https://img.shields.io/badge/docs-online-blue)](https://dongzhaohe321418-lab.github.io/materials-simulation-handbook/)

> A three-tier open handbook on computational materials science, taking a reader from A-level / first-year-university maths and physics all the way to running modern foundation-model workflows such as MACE-MP-0.

**Read online:** <https://dongzhaohe321418-lab.github.io/materials-simulation-handbook/>

---

## What is this?

The *Materials Simulation Handbook* is a self-contained, open-source textbook with accompanying runnable code, through the core of computational materials science, and with some small projects for undergraduate/ graduate level. The book is built in three explicit tiers — *Prerequisites*,* Core*, *Research*. You can read it linearly as a textbook, skip the early tiers if you already have the background, or work backwards from one of five capstone projects.

---

## Who is this for?

The **prerequisite floor** is A-level (or international equivalent) mathematics and physics, plus a willingness to install Python and run code from the command line. Concretely: you can differentiate and integrate single-variable functions; you have seen matrices and vectors; you can read a graph of energy vs. distance; you know what a wavefunction is even if you have never solved Schrödinger's equation by yourself :)

The **comfort zone** the textbook is primarily designed for undergraduates with relevant quantum chemistry/ quantum mechanics or material sicence background. You could use external LLMs, e.g. claude or gemini, to help you better understand the ideas you are less familiar with.

If you are not sure where you sit, do the [twenty-question prerequisites checker](docs/prerequisites-checker.md). It takes about half an hour and tells you which chapter to start at and which earlier ones to skim or skip.

---

## The three tiers

**Tier 0 — Prerequisites (Chapters 0, 1, 2, 3, 3.5).** Mathematics from numbers up to gradients and Fourier; Python and scientific computing; what materials simulation is and what it isn't; atoms, bonds, crystals, reciprocal space; and a compact solid-state physics primer covering Bloch's theorem, nearly-free-electron and tight-binding models, phonons, and defects. Outcome: you can read a band structure, write a NumPy loop, and recognise every symbol in the chapters that follow.

**Tier 1 — Core (Chapters 4, 5, 6, 7, 8, 9, 10).** The textbook proper. Quantum mechanics aimed at materials; density functional theory derived twice (once historically, once cleanly); running real DFT in Quantum ESPRESSO with convergence testing and defect calculations; molecular dynamics with thermostats and LAMMPS; statistical mechanics from a simulation point of view; machine-learning interatomic potentials from descriptors through Behler–Parrinello to MACE; and graph neural networks for crystals. Outcome: you can choose between DFT, MD, and an MLIP for a given problem, run all three, and explain the trade-offs.

**Tier 2 — Research (Chapters 11, 12, 13, 14).** The contemporary frontier and the meta-skills of doing research. Active learning and Bayesian optimisation for sample-efficient exploration; foundation models for atoms (MACE-MP-0, MatterGen and friends); multiscale coupling from QM/MM through coarse-graining to phase-field and finite-element bridges; and a final chapter on how to design, scope, run and write up your own project. Outcome: you can pick a problem in the literature, identify a workable method, plan an eight- to twelve-week project around it, and write a thesis chapter on what you found.

---

## The five capstone projects

The `docs/projects/` directory holds five self-contained projects, each with its own README, prerequisite list, and grading rubric. They are designed to be doable in four to ten weeks and to produce something publishable as a workshop note or thesis chapter.

1. **Project 1 — Defect Formation Energy.** Compute the vacancy formation energy of a chosen material in DFT, with full convergence analysis. *~4 weeks.*
2. **Project 2 — Melting Point from an MLIP.** Train or fine-tune a machine-learning potential, then estimate the melting point via two-phase coexistence MD. *~6 weeks.*
3. **Project 3 — Band-Gap Screening.** Build a GNN-based screening pipeline over a Materials Project subset and validate predictions against held-out DFT. *~6 weeks.*
4. **Project 4 — MLIP from Scratch.** Implement a Behler–Parrinello-style potential end-to-end and benchmark it against a published reference. *~8 weeks.*
5. **Project 5 — Bayesian Catalyst Search.** Use Bayesian optimisation with surrogate models to search a small catalyst space with a constrained budget of DFT evaluations. *~10 weeks.*

Each project lists exactly which chapters it depends on, so [Path C — Project-driven](docs/learning-path.md#path-c--project-driven-8-weeks) is a viable route through the book.

---

## Quick start

```bash
git clone https://github.com/dongzhaohe321418-lab/materials-simulation-handbook.git
cd materials-simulation-handbook
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
mkdocs serve   # local preview at http://127.0.0.1:8000
```

---

## Repository layout

```
materials-simulation-handbook/
├── docs/                # The book itself (Markdown, rendered with MkDocs Material)
│   ├── ch00-math/       # Tier 0
│   ├── ch01-python/
│   ├── ch02-foundations/
│   ├── ch03-atoms/
│   ├── ch03b-solid-state/
│   ├── ch04-quantum/    # Tier 1
│   ├── ch05-dft/
│   ├── ch06-running-dft/
│   ├── ch07-md/
│   ├── ch08-statmech/
│   ├── ch09-mlip/
│   ├── ch10-gnn/
│   ├── ch11-active/     # Tier 2
│   ├── ch12-foundation/
│   ├── ch13-multiscale/
│   ├── ch14-capstone/
│   ├── projects/        # Five capstone projects
│   ├── appendix/
│   ├── learning-path.md
│   └── prerequisites-checker.md
├── code/                # Runnable reference implementations, organised by chapter
├── notebooks/           # Jupyter notebooks (interactive companion)
├── exercises/           # Problem sets with solutions
├── mkdocs.yml           # Site configuration
└── .github/             # CI: builds and deploys the site on every push to main
```

---

## How to cite

```bibtex
@misc{dong2026materialshandbook,
  author       = {Zhaohe Dong},
  title        = {Materials Simulation Handbook},
  year         = {2026},
  howpublished = {\url{https://github.com/dongzhaohe321418-lab/materials-simulation-handbook}},
  note         = {Open educational resource. Three-tier curriculum from prerequisites to foundation models. Text licensed CC BY 4.0; code licensed MIT.}
}
```

If you cite a specific chapter, please give the chapter title and section number alongside the repository URL, e.g. *"Materials Simulation Handbook, Ch 9 §5 (equivariant networks)."*

---

## Contributing

Issues and pull requests are welcome — typo fixes, clearer derivations, additional worked examples and corrections to the physics are all valuable. See [CONTRIBUTING.md](CONTRIBUTING.md) for details and the [Code of Conduct](CODE_OF_CONDUCT.md). For larger changes (a new section, a new project, a new chapter) please open an issue first so the structure stays coherent.

---

## Licence

- **Text** (everything under `docs/`): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). You may reuse, adapt and redistribute the text, including for commercial purposes, provided you give appropriate credit.
- **Code** (everything under `code/`, `notebooks/`, `exercises/`): [MIT](LICENSE).
