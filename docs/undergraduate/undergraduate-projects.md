# Undergraduate Projects

This page collects six small, self-contained projects that you can run on a normal laptop in anything from a few days to about two weeks of part-time work. They are deliberately *gentler* than the five large capstone projects in [Chapter 14 (Capstone)](../ch14-capstone/index.md). The capstones assume you can already string several methods together and run a real simulation campaign; the projects here assume only that you have read the relevant chapter overview, can run a short Python script, and are willing to be patient. If a project below excites you, treat it as a warm-up for the matching capstone rather than a competitor to it.

Each project connects back to chapters you have already met, uses small inputs, and produces something concrete you can show a supervisor. None of them require you to "master" a method — that is not realistic in two weeks, and we will not pretend otherwise. The aim is narrower and more honest: get one thing working, understand why it works, and learn what breaks it.

## How to approach a small project

The most common way a small project goes wrong is that it quietly becomes a large one. Guard against that from the start.

- **Keep it small.** Resist the urge to add features. A box that holds one particle in one dimension teaches you more, faster, than a half-finished three-dimensional solver.
- **Get *something* running first.** Before you make anything correct or pretty, make it run end to end, even if the answer is wrong. A pipeline that runs and prints a wrong number is far closer to done than elegant code that has never executed.
- **Change one thing at a time.** When you tune a parameter — a grid spacing, a time step, a cutoff — change it alone and watch what moves. If you change three things and the answer improves, you have learned nothing about which one mattered.
- **Write down what you did.** Keep a plain text log: what you ran, what number came out, what you changed next. This is not bureaucracy; it is the difference between a result you can defend and a result you cannot reproduce next week.

!!! tip "Pick your first project by what you have"

    Projects 1, 2, 3 and 6 run entirely with NumPy, SciPy, Matplotlib and scikit-learn, so they work on any laptop and (for the NumPy ones) even in the browser sandbox. Projects 4 and 5 need a real density functional theory (DFT) code installed and some patience — start with the others if you have never run scientific software before.

A note on difficulty, if it helps you choose an order:
<span class="diff-easy">★ easy</span> Projects 1, 2 ·
<span class="diff-med">★★ medium</span> Projects 3, 6 ·
<span class="diff-hard">★★★ hard</span> Projects 4, 5.

---

## Project 1 — One-dimensional particle in a box

**Goal.** Solve the time-independent Schrödinger equation for a particle confined to a one-dimensional box using finite differences, and check your numerical energy levels against the textbook formula $E_n = \dfrac{n^2 \pi^2 \hbar^2}{2 m L^2}$.

**Related chapters.** [Chapter 1 (Python)](../ch01-python/index.md) for the array tools; [Chapter 4 (Quantum)](../ch04-quantum/index.md) for the physics of bound states and wavefunctions.

**Prerequisites.** You can build a NumPy array and index it. You have seen that a second derivative can be approximated by neighbouring grid points. You recognise an eigenvalue problem as "matrix in, special vectors and numbers out" — you do not need to be able to derive it.

**Estimated time.** Two to four days, part-time.

**Tools.** NumPy only (it runs live in the browser). Matplotlib for the plots.

**Step-by-step plan.**

1. Put $N$ grid points evenly across the box from $0$ to $L$ with spacing $h = L/(N+1)$, using the interior points only (the wavefunction is zero at both walls).
2. Approximate the second derivative with the standard three-point stencil $\psi''(x_i) \approx \dfrac{\psi_{i-1} - 2\psi_i + \psi_{i+1}}{h^2}$.
3. Assemble the Hamiltonian as a matrix: the kinetic term $-\dfrac{\hbar^2}{2m}\dfrac{d^2}{dx^2}$ becomes a tridiagonal matrix with $-2$ on the diagonal and $1$ on the two off-diagonals, scaled by $-\dfrac{\hbar^2}{2 m h^2}$. Inside the box the potential is zero, so there is nothing to add.
4. Diagonalise with `numpy.linalg.eigh`, which returns eigenvalues already sorted.
5. Compare the lowest five numerical eigenvalues to $E_n = \dfrac{n^2 \pi^2 \hbar^2}{2 m L^2}$.
6. Plot the first three eigenvectors as wavefunctions and confirm they show one, two and three half-waves.

!!! tip "Use the simplest units you can"

    Set $\hbar = 1$, $m = 1$ and $L = 1$. Then the exact levels are simply $E_n = n^2 \pi^2 / 2 \approx 4.93,\ 19.7,\ 44.4,\dots$ and you can eyeball whether your solver is right before worrying about physical units.

**What to submit.** A script, a table of numerical versus analytic energies for $n = 1\ldots5$ with the percentage error, and one figure of the first three wavefunctions.

**What counts as success.** Your lowest levels match the formula to better than about 1% with a few hundred grid points, the error shrinks as you add points, and the wavefunctions have the right number of nodes.

**Stretch goal.** Add a potential: put a step or a small barrier in the middle of the box and watch the levels shift and the wavefunctions deform.

**Common failure points.** Forgetting the $1/h^2$ factor (energies off by a constant); including the wall points so the box is effectively the wrong length; comparing to $E_n$ with a misremembered factor of two.

---

## Project 2 — Harmonic-oscillator molecular dynamics

**Goal.** Integrate the motion of a single one-dimensional harmonic oscillator with the velocity-Verlet algorithm, show that the total energy stays nearly constant over many oscillations, and measure how that conservation depends on the time step.

**Related chapters.** [Chapter 1 (Python)](../ch01-python/index.md) for loops and arrays; [Chapter 7 (Molecular Dynamics)](../ch07-md/index.md) for the integrator and the idea of energy conservation as a quality check.

**Prerequisites.** You know that force is mass times acceleration, and that for a spring $F = -k x$. You can write a `for` loop that updates two numbers each step.

**Estimated time.** Two to three days, part-time.

**Tools.** NumPy only (runs live in the browser). Matplotlib for plotting.

**Step-by-step plan.**

1. Choose $m = 1$ and $k = 1$, so the exact angular frequency is $\omega = 1$ and the period is $2\pi$.
2. Write the velocity-Verlet update: advance the position by $x + v\,\Delta t + \tfrac{1}{2} a\,\Delta t^2$, compute the new acceleration from the new position, then advance the velocity by $\tfrac{1}{2}(a_{\text{old}} + a_{\text{new}})\,\Delta t$.
3. Run from rest at $x = 1$ for several periods, storing position, velocity and time each step.
4. Compute the total energy $E = \tfrac{1}{2} m v^2 + \tfrac{1}{2} k x^2$ at every step and plot it against time.
5. Repeat with a few time steps, say $\Delta t = 0.2,\ 0.1,\ 0.05$ of the period, and overlay the energy traces.
6. Plot the position against the analytic solution $x(t) = \cos(\omega t)$ for the smallest time step.

**What to submit.** A script, one figure of energy versus time for the three time steps, and one figure comparing the trajectory to the exact cosine.

**What counts as success.** The energy oscillates within a small band rather than drifting steadily up or down, the band narrows as the time step shrinks, and the trajectory tracks the cosine for the smallest step.

!!! note "Why Verlet does not drift"

    A naive Euler step leaks energy and the oscillator either spirals in or blows up. Velocity-Verlet is *symplectic*: its energy wobbles but does not march off in one direction, which is exactly why molecular dynamics codes use it. Seeing this difference yourself is the whole point of the project.

**Stretch goal.** Implement plain forward Euler alongside Verlet and plot both energy traces together — the contrast is the clearest single lesson in the project.

**Common failure points.** Updating the velocity with the *old* acceleration only (that is Euler, not Verlet); using a time step near the period so the motion is badly under-sampled; comparing to $\sin$ when you started from rest.

---

## Project 3 — Lennard-Jones potential and a tiny MD analysis

**Goal.** Plot the Lennard-Jones potential and the force derived from it, run a very small two-dimensional Lennard-Jones simulation (or analyse a trajectory you are given), and compute a radial distribution function from the particle positions.

**Related chapters.** [Chapter 7 (Molecular Dynamics)](../ch07-md/index.md) for the simulation; [Chapter 8 (Statistical Mechanics)](../ch08-statmech/index.md) for what the radial distribution function means.

**Prerequisites.** You finished Project 2 or are otherwise comfortable with a Verlet loop. You can take a derivative of a simple function by hand. You can build a histogram with NumPy.

**Estimated time.** Five days to about two weeks, part-time, depending on whether you write the simulation or analyse a supplied trajectory.

**Tools.** NumPy and Matplotlib (both run live in the browser). ASE is optional: it can run the dynamics for you and provides neighbour lists, but it does not run in the browser and is not required — a plain NumPy loop with a handful of particles is enough.

**Step-by-step plan.**

1. Write the potential $V(r) = 4\varepsilon\!\left[\left(\dfrac{\sigma}{r}\right)^{12} - \left(\dfrac{\sigma}{r}\right)^{6}\right]$ and plot it; mark the minimum at $r = 2^{1/6}\sigma$.
2. Derive the force $F(r) = -\dfrac{dV}{dr}$ analytically, code it, and plot it; check that the force is zero exactly at the minimum of the potential.
3. Place a small grid of particles, say $16$ to $36$, in a two-dimensional periodic box at a modest density.
4. Run velocity-Verlet for a few thousand steps, computing pairwise forces with the minimum-image convention for the periodic boundary.
5. Save snapshots of the positions, or load the trajectory you were given.
6. Compute the radial distribution function $g(r)$ by histogramming all pair distances and normalising by the number expected in an ideal gas at the same density; plot it.

**What to submit.** A script, the potential-and-force figure, and a $g(r)$ plot with a short paragraph explaining its shape.

**What counts as success.** The force crosses zero at the potential minimum; the simulation runs without particles flying apart; and $g(r)$ is near zero at very short range, rises to a clear first peak near the equilibrium spacing, and settles towards one at large $r$.

**Stretch goal.** Run the same system at two temperatures and show how the first peak of $g(r)$ sharpens as the system cools towards a more ordered, solid-like arrangement.

**Common failure points.** Forgetting the minimum-image convention, so particles interact across the wrong side of the box; a time step too large for the steep $r^{-12}$ wall, which makes energies explode; normalising $g(r)$ incorrectly so it never approaches one.

---

## Project 4 — Your first DFT calculation of silicon

**Goal.** Run a real density functional theory calculation on silicon: relax or scan the lattice, extract the equilibrium lattice constant, and compute a band structure — while taking convergence seriously rather than trusting the first number that appears.

**Related chapters.** [Chapter 5 (DFT)](../ch05-dft/index.md) for what the calculation is actually solving; [Chapter 6 (Running DFT)](../ch06-running-dft/index.md) for the practical mechanics of inputs, pseudopotentials and outputs.

**Prerequisites.** You have read the Chapter 5 and Chapter 6 overviews. You can edit a text input file and run a command in a terminal. You understand, at least roughly, that a plane-wave calculation has a cutoff energy and a grid of $\mathbf{k}$-points.

**Estimated time.** One to two weeks, part-time — and genuinely so. Most of that time is installation, reading error messages, and waiting, not physics.

**Tools.** A real DFT code, most simply Quantum ESPRESSO, optionally driven through ASE in Python. You will also need a silicon pseudopotential.

!!! warning "This does not run in the browser"

    Quantum ESPRESSO and ASE-driven DFT are not part of the live browser sandbox and will not run there. You need them installed locally, or an account on a cloud notebook such as Google Colab, or access to an HPC cluster. Expect installation alone to take a session, and do not be discouraged when the first run fails — that is the normal experience, not a sign you have done something wrong.

**Step-by-step plan.**

1. Install Quantum ESPRESSO locally, or open a Colab notebook that installs it, or get an account on your group's cluster.
2. Download a silicon pseudopotential and confirm it matches the exchange-correlation functional you intend to use.
3. Build the diamond-cubic silicon cell (two atoms in the primitive cell) and write a single self-consistent input with a modest cutoff and a coarse $\mathbf{k}$-point grid.
4. Get one calculation to finish and converge — this milestone alone is worth a day.
5. Scan a range of lattice constants, run a self-consistent calculation at each, and find the energy minimum to get the equilibrium lattice constant; compare to the measured value near $5.43\ \text{Å}$.
6. At the relaxed lattice constant, run a non-self-consistent calculation along a standard $\mathbf{k}$-path and plot the band structure.

**What to submit.** Your input files, a plot of energy versus lattice constant with the minimum marked, your fitted lattice constant next to the experimental value, and a band-structure figure.

**What counts as success.** A calculation that converges cleanly; a lattice constant within roughly 1–2% of experiment; and a band structure that has the expected shape for silicon, including an indirect gap (do not be alarmed that the gap is too small — standard functionals are known to underestimate it).

**Stretch goal.** Recompute the lattice constant with a second exchange-correlation functional and report how much it shifts, as a first taste of functional dependence.

**Common failure points.** A cutoff or $\mathbf{k}$-grid too coarse, so the lattice constant is wrong but stable (Project 5 exists precisely to cure this); a pseudopotential that does not match the functional; mislabelling the high-symmetry points on the band-structure path.

---

## Project 5 — Cutoff and k-point convergence study

**Goal.** Take the silicon calculation from Project 4 and make it trustworthy: systematically vary the plane-wave cutoff and the $\mathbf{k}$-point density, plot the total energy against each, and decide for yourself what "converged" means in practice.

**Related chapters.** [Chapter 6 (Running DFT)](../ch06-running-dft/index.md), which is built around exactly this kind of numerical control.

**Prerequisites.** Project 4 is running and produces a total energy you can read from the output. You can write a short loop or a shell script that runs the same input with one parameter changed.

**Estimated time.** Three days to a week, part-time, on top of Project 4.

**Tools.** The same DFT code as Project 4. NumPy and Matplotlib to collect the numbers and plot the curves.

!!! warning "Same caveat as Project 4"

    This project also depends on real DFT software and will not run in the browser. Run it locally, on Colab, or on a cluster.

**Step-by-step plan.**

1. Fix the lattice constant and the $\mathbf{k}$-grid, then run the calculation at a ladder of increasing plane-wave cutoffs.
2. Plot total energy against cutoff and watch the curve flatten.
3. Choose the cutoff at which the energy stops changing by more than your chosen tolerance — a few meV per atom is a common rule of thumb.
4. Now fix the cutoff at that converged value and run a ladder of increasingly dense $\mathbf{k}$-grids.
5. Plot total energy against $\mathbf{k}$-grid density and again pick where it flattens.
6. Re-run the lattice-constant scan from Project 4 with both parameters at their converged values and note whether your lattice constant moved.

**What to submit.** Two convergence plots (energy versus cutoff, energy versus $\mathbf{k}$-density), your chosen converged values with the tolerance you used, and a short paragraph on whether the Project 4 lattice constant changed once converged.

**What counts as success.** Both curves visibly flatten; you can state a cutoff and a $\mathbf{k}$-grid that hold the energy steady to your stated tolerance; and you can explain in one sentence why a number that looks stable can still be wrong if you never tested it.

!!! note "Convergence is the habit, not the homework"

    Every total energy, force or lattice constant a DFT code prints is only as trustworthy as the cutoff and $\mathbf{k}$-grid behind it. Doing this study once, by hand, builds the instinct to never quote a number you have not converged. That instinct matters more than the silicon result itself.

**Stretch goal.** Show that energy *differences* (such as the lattice-constant scan) converge faster than absolute total energies, because systematic errors partly cancel.

**Common failure points.** Changing the cutoff and the $\mathbf{k}$-grid at the same time, so you cannot tell which drove the change; comparing absolute total energies across runs without realising they are not directly meaningful; declaring convergence after one extra point instead of confirming the curve has truly levelled.

---

## Project 6 — A simple machine-learning model for a materials property

**Goal.** Take a small dataset of materials, turn each material into a handful of numerical features, fit a simple regression model to predict a property such as formation energy, and evaluate it with an honest train/test split.

**Related chapters.** [Chapter 9 (Machine-Learning Interatomic Potentials)](../ch09-mlip/index.md) for where machine learning meets materials; [Chapter 11 (Active Learning)](../ch11-active/index.md) for how models and data improve together once you have a baseline.

**Prerequisites.** You can load a table with NumPy or pandas. You understand that you must not test a model on the same data you trained it on. You do not need any prior machine-learning experience.

**Estimated time.** Five days to about two weeks, part-time.

**Tools.** scikit-learn, NumPy and Matplotlib. (scikit-learn does not run in the browser sandbox, but it installs easily on any laptop.) Keep the dataset and the model small on purpose.

**Step-by-step plan.**

1. Get a small dataset — a few hundred to a few thousand entries — of compositions paired with a target property such as formation energy.
2. Build simple composition-based features: for example the element-weighted averages of properties like atomic mass, electronegativity and group number. A handful of features is plenty to start.
3. Split the data into training and test sets with `train_test_split`, holding out the test set and never looking at it during fitting.
4. Fit a baseline model — ordinary linear regression, or ridge regression — on the training set only.
5. Predict on the held-out test set and report the mean absolute error and the coefficient of determination $R^2$.
6. Plot predicted against true values for the test set; a good model hugs the diagonal.

**What to submit.** A script, the test-set error and $R^2$, and a predicted-versus-true scatter plot.

**What counts as success.** Your model beats the trivial baseline of always predicting the mean; the test error is reported on data the model never saw; and you can name one feature that mattered and one limitation of your approach.

!!! warning "Be honest about the split"

    The single most common way a beginner's model looks better than it is comes from leakage: scaling features using the whole dataset, tuning until the test score improves, or letting near-duplicate entries fall on both sides of the split. Fit every transformation on the training set alone, decide your model *before* you look at the test score, and report that score once.

**Stretch goal.** Try a slightly richer model such as a small random forest and compare it fairly to your linear baseline using the same split — sometimes the simple model is competitive, which is itself a useful finding.

**Common failure points.** Reporting the training error and calling it performance; standardising features using statistics from the full dataset; chasing a higher test score by repeatedly peeking at it, which quietly turns your test set into a second training set.

---

!!! info "When you have finished one of these"

    A completed small project is the natural springboard to the matching large project. When you are ready for more, read the [Chapter 14 (Capstone)](../ch14-capstone/index.md) overview: the convergence habits from Projects 4 and 5, the dynamics from Projects 2 and 3, and the modelling discipline from Project 6 are exactly the foundations the capstones build on.
