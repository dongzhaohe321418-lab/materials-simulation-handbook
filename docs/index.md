---
title: Materials Simulation Handbook
hide:
  - navigation
---

# Materials Simulation Handbook

> A three-tier, textbook-grade open handbook on computational materials science — from the maths you need, through DFT and MD, all the way to foundation models such as MACE-MP-0.

This book is for anyone who wants to understand **how we simulate matter on a computer** — the electrons inside silicon, the atoms diffusing in a battery, the alloys we have not yet discovered — and who wants to do it themselves, with real code, on a real laptop.

You do **not** need a physics or chemistry degree to start, but you do need A-level (or equivalent first-year university) mathematics and physics: you should be comfortable differentiating and integrating, you should have met matrices and vectors, and you should at least have heard of the Schrödinger equation. If any of that sounds shaky, take the [prerequisites self-test](prerequisites-checker.md) — it diagnoses in 20 questions whether you should start at Chapter 0 or skip ahead.

---

## What you will learn

By the end of this handbook you will be able to:

1.  Read a paper that says *"we performed DFT calculations with PBE and a 60 Ry plane-wave cutoff"* and know exactly what each word means.
2.  Run your own DFT calculation in Quantum ESPRESSO and reproduce the band structure of silicon.
3.  Set up a molecular dynamics simulation in LAMMPS and compute thermodynamic averages.
4.  Train a machine-learning interatomic potential (MACE) on a small dataset and use it to run MD.
5.  Build a graph neural network from scratch in PyTorch to predict materials properties.
6.  Use Bayesian optimisation to drive an active-learning loop for materials discovery.
7.  Use foundation models like MACE-MP-0 zero-shot across the periodic table.

---

## How the book is structured

The handbook is divided into five parts:

| Part | Chapters | What it covers |
|------|----------|----------------|
| **I. Prerequisites** | 0 – 1 | Mathematics and Python you need before anything else |
| **II. The big picture** | 2 – 3 | What we are simulating and why; atoms, bonds, and crystals |
| **III. Electrons** | 4 – 6 | Quantum mechanics, DFT, and running real calculations |
| **IV. Atoms in motion** | 7 – 8 | Molecular dynamics and statistical mechanics |
| **V. Modern machine learning** | 9 – 12 | MLIPs, GNNs, active learning, foundation models |

Each chapter has the same structure: a conceptual overview, the mathematics, a worked example with code you can run, exercises (with full solutions), and an annotated further-reading list.

---

## How to read it

- **If you are starting from zero**: read Chapter 0 and Chapter 1 carefully, then proceed in order.
- **If you already know first-year university maths and Python**: skim Chapter 0 and Chapter 1, then start at Chapter 2.
- **If you are a chemist who wants to learn ML**: start with Chapter 2 for orientation, then jump to Part V.
- **If you are an ML researcher who wants to learn materials**: read Chapters 2 – 9 in order; the maths chapters will go quickly.

The book is deliberately self-contained. If something appears that you do not understand, the previous chapter explained it — go back and re-read.

---

## Quick start

```bash
# 1. Clone the repository
git clone https://github.com/dongzhaohe321418-lab/materials-simulation-handbook.git
cd materials-simulation-handbook

# 2. Create the environment
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 3. Read locally (optional — the site is also online)
mkdocs serve
# Open http://127.0.0.1:8000
```

For the scientific Python stack used in the worked examples, see [Chapter 1.1](ch01-python/01-install.md).

---

## Why I wrote this

There are many excellent textbooks on density functional theory, on molecular dynamics, and increasingly on machine learning for materials. There is no single book that takes a complete beginner from "what is a vector" through to "fine-tune a foundation model for a catalyst dataset". This book tries to be that bridge.

It is also free and open. Every page can be edited; every code example can be improved. Pull requests are welcome.

---

## How to cite

If you use this handbook in a course or paper, please cite it as:

> Dong, Z. (2026). *Materials Simulation Handbook* [Open educational resource]. https://github.com/dongzhaohe321418-lab/materials-simulation-handbook

A machine-readable `CITATION.cff` is available in the repository root.

---

## Acknowledgements

This handbook stands on the shoulders of the textbooks of Martin, Frenkel & Smit, Marder, and Sholl & Steckel, and on the open-source codes — Quantum ESPRESSO, LAMMPS, ASE, pymatgen, PyTorch, MACE — without which none of this would be possible.
