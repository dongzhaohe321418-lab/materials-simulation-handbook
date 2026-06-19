# Undergraduate Learning Paths

These are slower, self-study routes through the chapters the handbook already contains. They do not add new physics; they take the existing material and lay it out as a gentle, time-budgeted plan you can follow on your own, a few evenings a week. Each path tells you which chapters to read, in which order, what to read carefully and what to skim, roughly how long it takes, and a small project to finish on.

The canonical map of the whole book is [the main Learning Path](../learning-path.md). That page shows the full dependency graph, the per-chapter hour estimates, and three routes aimed at readers who already have a strong background. The paths here are the gentler cousins of those routes: they cover fewer chapters, move more slowly, and assume you are teaching yourself with no one to ask. Where the two pages overlap they agree; where they differ, the main page is the authority on dependencies and this page is the authority on pacing for a self-taught beginner.

If you are not sure where to start, spend twenty minutes on the [prerequisites checker](../prerequisites-checker.md) first. It points you at a specific starting chapter based on what you already know, and it will tell you whether you can safely skim Chapter 0 or should read it line by line.

!!! warning "Be realistic about how long this takes"

    You cannot master density functional theory, molecular dynamics or graph neural networks in a few days, and no honest plan will tell you otherwise. Each of the paths below is **weeks of part-time study** — typically two or three focused evenings a week over two to four months. That is not slow; that is the normal pace for learning a research-level subject from scratch. If a week feels hard, that is the material being hard, not you being slow. Plan for the long haul, rest when a chapter stops going in, and come back to it. Steady part-time progress finishes the path; cramming does not.

A note on the tables below: "Read" means the sections to study that week, "Do" means the by-hand work, code, or exercises to attempt. You are not expected to finish everything in a week — attempting it and noting what you did not understand is the point.

---

## Path A — Materials Simulation Core

### Goal

Understand what density functional theory (DFT) is, why it works, and run your first electronic-structure calculation on a real material.

### Who should choose this

You want to understand how the energies and forces in modern materials science are actually computed from quantum mechanics, and you want to come out the other side having run a genuine DFT calculation yourself. You are comfortable with the idea of wavefunctions and energy levels from first-year physics but have never seen them used in anger. This is the longest of the three paths and the most conceptually demanding.

### Chapter sequence

1. [Chapter 0 — Mathematics](../ch00-math/index.md)
2. [Chapter 1 — Python & Scientific Computing](../ch01-python/index.md)
3. [Chapter 2 — What is Materials Simulation?](../ch02-foundations/index.md)
4. [Chapter 3 — Atoms, Bonds & Solids](../ch03-atoms/index.md)
5. [Chapter 3.5 — Solid-State Physics](../ch03b-solid-state/index.md)
6. [Chapter 4 — Quantum Mechanics](../ch04-quantum/index.md)
7. [Chapter 5 — Density Functional Theory](../ch05-dft/index.md)
8. [Chapter 6 — Running DFT in Practice](../ch06-running-dft/index.md)

### What to read carefully vs skim on a first pass

Read **carefully**: the linear-algebra and complex-number parts of Chapter 0 (they underpin everything quantum), the band-structure and reciprocal-space ideas in Chapter 3.5, the single-particle Schrödinger equation in Chapter 4, and the Hohenberg–Kohn and Kohn–Sham sections of Chapter 5 (Sections 5.2 and 5.3) — these are the conceptual heart of the path.

**Skim** on a first pass: the Thomas–Fermi opening of Chapter 5 (Section 5.1) and Hartree–Fock in Chapter 4 — both are historically important but not load-bearing for a first working picture, exactly as the main Learning Path notes. You can also skim the exchange–correlation zoo in Chapter 5 (Section 5.4) the first time through; learn that the choice of functional matters and come back for the details.

### Approximate study time

About **10–14 weeks** of part-time study (roughly two or three evenings a week). The quantum and DFT chapters are the slow part; do not be surprised if Chapters 4 and 5 each take two weeks on their own.

### A weekly plan

| Week | Focus | Read | Do |
|---|---|---|---|
| 1 | Maths refresh | Ch 0 (linear algebra, calculus, complex numbers) | Work the by-hand exercises; flag what is rusty |
| 2 | Python & environment | Ch 1 | Install Python, run the notebooks, change a number and rerun |
| 3 | The big picture | Ch 2, start Ch 3 | Read the index pages; sketch the simulation landscape |
| 4 | Atoms and crystals | Ch 3 | Build a unit cell on paper; do the crystal exercises |
| 5–6 | Solid-state ideas | Ch 3.5 | Reciprocal lattice, Brillouin zone, simple band picture by hand |
| 7–8 | Quantum mechanics | Ch 4 (skim Hartree–Fock) | Solve the particle-in-a-box; run the 1D solver notebook |
| 9–10 | DFT theory | Ch 5 (Sections 5.2–5.3 carefully) | Run the 1D Kohn–Sham SCF program; watch it converge |
| 11–12 | Running DFT | Ch 6 | Install a DFT code; converge a calculation on one material |
| 13–14 | Consolidate + project | Re-read Ch 5–6 index pages | Start the final project below |

### Code to run

In Chapter 5, run the self-consistent-field program that solves a 1D model hydrogen chain end to end — watching it iterate to convergence is the single best way to make Kohn–Sham theory concrete. It is NumPy/SciPy, so it runs live in the browser. In Chapter 6, run a real DFT calculation in Quantum ESPRESSO or via ASE; that code does **not** run in the browser, so follow the chapter's installation notes and run it locally or on a free cluster.

### Exercises to prioritise

From Chapter 4, the single-particle Schrödinger problems. From Chapter 5, re-deriving the Hohenberg–Kohn argument in your own words and the exercise that extends the SCF code. From Chapter 6, the convergence-testing exercises (k-points and plane-wave cutoff) — convergence testing is the practical skill that separates a usable calculation from a meaningless one.

### A final small project

Finish with a first electronic-structure project from the [undergraduate projects](undergraduate-projects.md) page — converge the total energy of a simple crystal such as silicon and compute its lattice constant, then compare with the accepted value.

---

## Path B — Molecular Dynamics Core

### Goal

Understand how atoms move in time under forces, and be able to run and analyse simple molecular-dynamics (MD) output.

### Who should choose this

You are drawn to the dynamic, finite-temperature side of simulation — melting, diffusion, atoms jostling in a box — rather than the zero-temperature quantum picture. This path is shorter and, for most beginners, gentler than Path A, because the core idea (integrate Newton's second law) is one you already half-know from mechanics.

### Chapter sequence

1. [Chapter 0 — Mathematics](../ch00-math/index.md)
2. [Chapter 1 — Python & Scientific Computing](../ch01-python/index.md)
3. [Chapter 2 — What is Materials Simulation?](../ch02-foundations/index.md)
4. [Chapter 3 — Atoms, Bonds & Solids](../ch03-atoms/index.md)
5. [Chapter 7 — Molecular Dynamics](../ch07-md/index.md)
6. [Chapter 8 — Statistical Mechanics](../ch08-statmech/index.md)

### What to read carefully vs skim on a first pass

Read **carefully**: the ordinary-differential-equation and probability parts of Chapter 0, the integration section of Chapter 7 (Section 7.1, the velocity-Verlet integrator and time-step choice), and the ensemble idea in Chapter 8 — knowing which ensemble (NVE, NVT, NPT) you are sampling is the conceptual key to the whole path.

**Skim** on a first pass: the full survey of force fields in Chapter 7 (Section 7.4) — learn that Lennard-Jones and EAM exist and why hand-tuned potentials are limited, but do not memorise each form. The Ewald-summation details for long-range Coulomb interactions (Section 7.2) can also wait; understand the minimum-image convention first.

### Approximate study time

About **7–10 weeks** of part-time study. The integrator and analysis ideas land quickly once you have written one Verlet step yourself.

### A weekly plan

| Week | Focus | Read | Do |
|---|---|---|---|
| 1 | Maths refresh | Ch 0 (ODE intuition, probability) | Work through a simple numerical ODE by hand |
| 2 | Python & environment | Ch 1 | Run the notebooks; plot a trajectory |
| 3 | Big picture + atoms | Ch 2, Ch 3 | Build a small periodic box of atoms |
| 4–5 | Integration | Ch 7 (Section 7.1 carefully) | Write a velocity-Verlet integrator in NumPy; check energy conservation |
| 6 | Boundaries & thermostats | Ch 7 (Sections 7.2–7.3) | Run a Lennard-Jones argon simulation |
| 7 | Trajectory analysis | Ch 7 (Section 7.6) | Compute a radial distribution function and mean-squared displacement |
| 8–9 | Statistical mechanics | Ch 8 | Connect temperature and ensemble to your trajectory; do the sampling exercises |
| 10 | Consolidate + project | Re-read Ch 7–8 index pages | Start the final project below |

### Code to run

The centrepiece is the velocity-Verlet integrator in Chapter 7, written from scratch in NumPy — it runs live in the browser and conserving energy over a long run is the satisfying proof that you have it right. The trajectory-analysis code (radial distribution function, mean-squared displacement) is also pure NumPy. The full LAMMPS argon and copper-melting workflows do **not** run in the browser; follow the chapter to run them locally.

### Exercises to prioritise

From Chapter 7, the time-step-stability exercise (push the step size until energy drifts) and the diffusion-coefficient calculation from a mean-squared displacement. From Chapter 8, the exercises that ask which ensemble a given setup samples — these tie the dynamics back to thermodynamics.

### A final small project

Finish with an MD project from the [undergraduate projects](undergraduate-projects.md) page — for example, run a Lennard-Jones liquid and measure its self-diffusion coefficient as a function of temperature.

---

## Path C — Machine Learning for Materials

### Goal

Understand, at an undergraduate level, what machine-learning interatomic potentials (MLIPs), graph neural networks (GNNs) and active learning are, and why they have changed the field.

### Who should choose this

You are interested in the machine-learning side of modern materials simulation and want to understand the modern pipeline rather than re-derive quantum mechanics from scratch. This is a frontier-facing path.

!!! note "Chapters 9 and 10 expect some comfort with Python and linear algebra"

    The machine-learning chapters move faster than the rest of the book and assume you can read a moderate Python script without panic and are comfortable with vectors and matrices. If that is not yet true, spend real time on Chapters 0 and 1 first — it is not wasted, and the later chapters will be far less bruising for it. There is no shame in taking an extra week here.

### Chapter sequence

1. [Chapter 0 — Mathematics](../ch00-math/index.md)
2. [Chapter 1 — Python & Scientific Computing](../ch01-python/index.md)
3. [Chapter 2 — What is Materials Simulation?](../ch02-foundations/index.md)
4. [Chapter 3 — Atoms, Bonds & Solids](../ch03-atoms/index.md)
5. [Chapter 7 — Molecular Dynamics](../ch07-md/index.md) — MD basics only (where the potentials are used)
6. [Chapter 9 — ML Interatomic Potentials](../ch09-mlip/index.md)
7. [Chapter 10 — Graph Neural Networks](../ch10-gnn/index.md)
8. [Chapter 11 — Active Learning & Bayesian Optimisation](../ch11-active/index.md)

### What to read carefully vs skim on a first pass

Read **carefully**: the linear-algebra of Chapter 0; the symmetry-and-descriptor framework in Chapter 9 (an interatomic potential must be invariant under translation, rotation and permutation — this idea recurs everywhere); the message-passing core of Chapter 10; and the active-learning loop in Chapter 11 (surrogate model → acquisition function → query → repeat).

**Skim** on a first pass: from Chapter 7, take only the basics — what a trajectory is and what a force model is for — since this path uses MD as context, not as a destination. Following the main Learning Path's guidance, skim the historical first-generation potentials in Chapter 9 (Behler–Parrinello and GAP, Section 9.4) and the comparative architecture survey in Chapter 10 (Section 10.4) on a first pass.

### Approximate study time

About **10–14 weeks** of part-time study. Chapters 9 and 10 are dense; budget two weeks each and do not rush them.

### A weekly plan

| Week | Focus | Read | Do |
|---|---|---|---|
| 1 | Maths refresh | Ch 0 (linear algebra) | Matrix and vector exercises until comfortable |
| 2 | Python & environment | Ch 1 | Run the notebooks; get comfortable reading scripts |
| 3 | Big picture + atoms | Ch 2, Ch 3 | Represent a structure as atoms with positions and types |
| 4 | MD basics (context) | Ch 7 (skim, integration idea only) | Run one short trajectory to see what data looks like |
| 5–6 | MLIPs | Ch 9 (Sections 9.2–9.3 carefully; skim 9.4) | Inspect descriptors; run the provided MLIP notebook |
| 7–8 | Graph neural networks | Ch 10 (Sections 10.1–10.3) | Run the GNN property-prediction notebook |
| 9–10 | Active learning | Ch 11 | Step through the active-learning loop; do the acquisition exercises |
| 11–12 | Consolidate + project | Re-read Ch 9–11 index pages | Start the final project below |

### Code to run

The MLIP and GNN notebooks in Chapters 9 and 10 are the core of this path, and you can launch them in Google Colab from the button at the top of each chapter. Note that these use PyTorch and related machine-learning libraries, which do **not** run live in the handbook's browser sandbox — Colab or a local GPU-free machine is the way to run them. The active-learning surrogate examples in Chapter 11 include simpler NumPy/SciPy pieces that do run in the browser.

### Exercises to prioritise

From Chapter 9, the exercises on why invariance and locality matter — they are the conceptual foundation. From Chapter 10, the exercise that runs the trained model on a held-out structure and checks the prediction. From Chapter 11, building one full pass of the active-learning loop by hand on a toy function.

### A final small project

Finish with a machine-learning project from the [undergraduate projects](undergraduate-projects.md) page — for example, train a small interatomic potential on a provided dataset and test how well it predicts forces on configurations it has never seen.

---

!!! tip "Combining or shortening the paths"

    These routes are deliberately separate, but they share a spine: Chapters 0, 1, 2 and 3 appear at the start of all three. If you intend to do more than one path, do that shared foundation **once** and keep your notes — you will not need to repeat it.

    To **shorten** a path, lean on the prerequisites checker. If it tells you your maths and Python are already solid, you can compress Chapters 0 and 1 into a quick skim (a few hours rather than two weeks) and start at the big-picture chapters. The main [Learning Path](../learning-path.md) describes faster, less hand-held routes for readers with a strong background — switch to those once a path stops feeling like a stretch.

    To **combine** paths, a natural order is A then B (DFT gives you the forces; MD puts the atoms in motion under them) and then C (machine-learned potentials replace DFT inside MD, so both earlier paths feed it). Doing all three back to back is most of the core book — expect several months of part-time study, and treat each path's final project as a checkpoint rather than racing to the next chapter.
