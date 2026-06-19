# v1.4 Undergraduate-Expansion Roadmap

This roadmap proposes *where* the existing chapters of the *Materials Simulation Handbook* need more scaffolding for self-study undergraduates, and roughly *what kind* of scaffolding each one needs. It is a planning document for the **v1.4 milestone**, not an expansion in itself: it guides future work but does not add the material. When the work is actually done, it should be applied in place, chapter by chapter, using the [in-place chapter expansion template](chapter-expansion-template.md), so that the additions sit inside the existing chapters rather than in a parallel beginner's book.

The audience we are scaffolding for is the one described throughout this undergraduate layer: a materials, physics, chemistry or engineering undergraduate with first-year maths and physics, some Python, and little research experience, who tires when many symbols and acronyms arrive at once. The principle is *add, never dumb down*: intuition before formalism, definitions before derivations, toy examples before full theory, and honest signposting of what is hard. Nothing here should reduce the scientific content of the existing chapters; it should lower the activation energy for reaching it.

**Chapter 5 (DFT) has been used as the pilot** for this approach — it is the chapter where the template's pause-and-recall boxes, hidden-prerequisite call-outs and toy-before-theory ordering were first tried — so where this roadmap says "as in the Chapter 5 pilot", it means a pattern that already exists and can be copied.

One of the v1.4 milestone's explicit goals is to add **"pause-and-recall / check-yourself boxes"** — short `!!! question "Check yourself"` admonitions with a collapsible `??? success "Answer"` — at the points where a reader has just been handed a new idea and is most likely to read on without having absorbed it. This roadmap flags, per chapter, the best places for them. For honesty about study time, the per-chapter hour estimates in the [learning path](../learning-path.md) and the gentler routes in [learning-paths.md](learning-paths.md) should be treated as the baseline; scaffolding does not make a hard chapter quick.

Each subsection below uses the same labelled fields so that the work can be triaged. Priority is the practical question: *if we only get to a few chapters before the v1.4 release, which ones move the needle most?*

---

## Chapter 0 — [Mathematics from Scratch](../ch00-math/index.md)

- **Likely difficulty for undergraduates:** Low to Medium. The individual topics are familiar from school, but the chapter compresses a lot of notation that the rest of the book then assumes silently.
- **Where it probably jumps too quickly.** Eigenvalues and eigenvectors, and the variational idea in the calculus section, get less space than they carry later. Almost every hard chapter (4, 5, 3.5) leans on "diagonalise this Hermitian matrix" as if it were obvious.
- **Hidden prerequisites.** Comfort reading $\sum$ and $\int$ as *operations you could actually carry out*, not just decoration; the idea that a matrix can *act on* a vector; that a derivative is a limit, not just a rule.
- **Concepts needing slower explanation.** Eigen-decomposition geometrically (stretch directions); what "a functional" means before §5 ever needs it; the gradient $\nabla$ and Laplacian $\nabla^2$ as the same objects in 1D and 3D.
- **Suggested toy examples.** Diagonalise a single $2\times2$ symmetric matrix by hand and *see* the eigenvectors as axes; Fourier-decompose a square wave into a few sines.
- **Suggested worked examples.** One fully worked eigenvalue problem reused verbatim later in Chapter 4 (so the reader meets it twice); one worked "minimise this functional" toy that prefigures Hohenberg–Kohn.
- **Suggested code-lab improvements.** A short NumPy/Matplotlib cell that plots eigenvectors of a $2\times2$ matrix as arrows (runs live in the browser); a Fourier-reconstruction slider notebook.
- **Suggested "Common Misunderstandings" topics.** "Eigenvalues are just numbers you compute" vs. what they mean; confusing a vector's components with the vector; thinking Fourier transforms are only for signals, not for crystals.
- **Suggested exercise improvements.** Add 2–3 *graded* warm-ups before the existing problems, with `??? success` answers, plus one "you will see this again in Chapter N" tagged problem per topic.
- **Priority: High.** This is the foundation: every High-priority electronic-structure chapter inherits its notation here. A small investment in eigen/Fourier/functional intuition pays off four chapters deep, so it earns High despite its low intrinsic difficulty.

---

## Chapter 1 — [Python and Scientific Computing](../ch01-python/index.md)

- **Likely difficulty for undergraduates:** Low. The barrier is environment setup and patience, not concepts.
- **Where it probably jumps too quickly.** The jump from "I can write a loop" to "I can think in vectorised NumPy" is real and often skipped; reproducibility/version-control material can feel like chores rather than tools.
- **Hidden prerequisites.** A working terminal mindset; knowing that an error message is information, not failure; the idea of an array shape.
- **Concepts needing slower explanation.** Broadcasting; the difference between a Python list and a NumPy array; why `import` and environments matter before they have ever bitten the reader.
- **Suggested toy examples.** Re-implement one Chapter 0 maths operation (e.g. a dot product) first with a loop, then vectorised, and time both.
- **Suggested worked examples.** A single end-to-end "read numbers, compute, plot" mini-pipeline the reader can later recognise inside every later notebook.
- **Suggested code-lab improvements.** A "debugging clinic" cell with three deliberately broken snippets and `??? note "Hint"` boxes; all live-runnable since it is pure NumPy/Matplotlib.
- **Suggested "Common Misunderstandings" topics.** "NumPy is just faster lists"; shape errors; thinking a plot that renders is a plot that is correct.
- **Suggested exercise improvements.** Add small, confidence-building exercises that mirror code the reader will meet in Chapters 4–9, so the syntax is not new when the physics is.
- **Priority: Medium.** Essential but not the bottleneck; most readers self-rescue here. Worth a light pass, not a heavy one.

---

## Chapter 2 — [What is Materials Simulation?](../ch02-foundations/index.md)

- **Likely difficulty for undergraduates:** Low. This is a map chapter, mostly prose.
- **Where it probably jumps too quickly.** The "scale ladder" can introduce many method names (DFT, MD, MLIP, kMC, phase-field) before any of them mean anything.
- **Hidden prerequisites.** Almost none — which is its strength. The risk is the opposite: acronyms used here as orientation reappear as load-bearing later.
- **Concepts needing slower explanation.** What "ab initio" actually promises; the trade-off triangle of accuracy, system size and time, stated once and referred back to.
- **Suggested toy examples.** A single table mapping a real question ("why does ice float?") to the scale and method that answers it.
- **Suggested worked examples.** Not a derivation chapter; instead a worked "reading a phase diagram" walk-through tied to (Section 2.3).
- **Suggested code-lab improvements.** Keep code light here; perhaps one tiny plot of length/time scales. Do not over-engineer.
- **Suggested "Common Misunderstandings" topics.** "Simulation = exact"; "bigger model = better"; thinking each method is a rival rather than a rung on a ladder.
- **Suggested exercise improvements.** Add reflective "which method would you reach for, and why?" prompts with model answers, to rehearse judgement early.
- **Priority: Low.** Already accessible; its main job is to forward-reference, which a glossary cross-link can support cheaply.

---

## Chapter 3 — [Atoms, Bonds, and Solids](../ch03-atoms/index.md)

- **Likely difficulty for undergraduates:** Low to Medium. Familiar chemistry, with the first appearance of reciprocal space, which is genuinely hard.
- **Where it probably jumps too quickly.** Reciprocal space and the structure factor. "Reciprocal space without tears" still tends to produce a few tears.
- **Hidden prerequisites.** Fourier intuition from Chapter 0; comfort with 3D vectors and dot products; that a lattice is a set of *translations*, not just dots on paper.
- **Concepts needing slower explanation.** Why reciprocal space exists at all (it is the natural home of anything periodic); the relationship between a real lattice vector and its reciprocal partner.
- **Suggested toy examples.** A 1D "lattice" of equally spaced points and its reciprocal, worked fully before any 3D Bravais lattice appears.
- **Suggested worked examples.** One worked structure-factor calculation for a two-atom basis, every term shown, reused as the template for the exercises.
- **Suggested code-lab improvements.** A live NumPy/Matplotlib cell that draws a 2D lattice and its reciprocal side by side and lets the reader change the basis vectors.
- **Suggested "Common Misunderstandings" topics.** "Reciprocal space is imaginary/abstract nonsense"; confusing the unit cell with the primitive cell; thinking the reciprocal lattice has units of length.
- **Suggested exercise improvements.** Add a scaffolded ladder from 1D to 2D to 3D reciprocal lattices, with `??? success` checkpoints, so the dimensional jump is gradual.
- **Priority: Medium.** It feeds the High-priority Chapter 3.5; getting reciprocal space right here removes a recurring stumbling block, but the heaviest lifting is downstream.

---

## Chapter 3.5 — [Solid State Physics Prerequisites](../ch03b-solid-state/index.md)

- **Likely difficulty for undergraduates:** High. At ~16 hours it is one of the longest prerequisite chapters and packs Bloch, bands, phonons and thermal models together.
- **Where it probably jumps too quickly.** The Bloch's-theorem proof from translational symmetry; the leap from "a $2\times2$ tight-binding matrix" to graphene's Dirac cones; the diatomic-chain acoustic/optical split.
- **Hidden prerequisites.** Reciprocal space (Ch 3); eigenvalues of $2\times2$ Hermitian matrices and bra-ket notation (Ch 4); Fourier series (Ch 0); crucially, the reading-order dependency on Chapter 4 that the chapter itself flags.
- **Concepts needing slower explanation.** What crystal momentum $\mathbf{k}$ *is* (a label, not an ordinary momentum); why a band gap opens at a zone boundary; the physical meaning of an acoustic vs. optical branch.
- **Suggested toy examples.** The 1D monatomic chain ($\omega = 2\sqrt{K/m}\,|\sin(ka/2)|$) carried all the way through by hand before the diatomic and 3D cases; an "empty lattice" band folded back into the first zone.
- **Suggested worked examples.** A fully worked 1D tight-binding chain giving $E = -2t\cos ka$, with the dispersion plotted, used as the on-ramp to the graphene code.
- **Suggested code-lab improvements.** Split the graphene diagonalisation lab into a *first* pure-NumPy 1D-chain band plot (browser-runnable) before the 2D graphene case; add a phonon-dispersion plotting cell for the diatomic chain.
- **Suggested "Common Misunderstandings" topics.** "$\mathbf{k}$ is the electron's momentum"; "more $k$-points always means more accuracy regardless of cost"; confusing a band index with an energy level; thinking phonons are particles you could catch.
- **Suggested exercise improvements.** Insert pause-and-recall boxes after Bloch's theorem and after the first dispersion relation; add a guided multi-part problem that builds the diatomic chain step by step.
- **Priority: High.** This is the electronic-structure on-ramp. Misunderstandings here (especially of $\mathbf{k}$ and $k$-point sampling) directly damage comprehension of Chapters 5 and 6, where they become the most-tuned parameter. It is the single chapter most worth scaffolding.

---

## Chapter 4 — [Quantum Mechanics for Materials](../ch04-quantum/index.md)

- **Likely difficulty for undergraduates:** High. At ~18 hours it is the longest Tier-1 chapter and, for most readers, their first formal quantum mechanics.
- **Where it probably jumps too quickly.** The postulates and the Born rule arrive fast in §4.2; the "exponential wall" arithmetic in §4.5 can land as a number rather than a felt impossibility; Hartree–Fock (§4.7) is dense.
- **Hidden prerequisites.** Complex numbers and complex exponentials (Ch 0); eigenvalue problems as the *same* thing as $\hat{H}\psi=E\psi$; the gradient and Laplacian; comfort that an operator is something that *acts on* a function.
- **Concepts needing slower explanation.** What a wavefunction *is* before $|\psi|^2$ is interpreted; why Hermitian operators give real eigenvalues, stated as "that is why energies are real"; antisymmetry and the Slater determinant; the *meaning* of the $3N$-dimensional configuration space.
- **Suggested toy examples.** The particle in a 1D box solved on a tiny grid by hand-sized matrices the reader can almost diagonalise mentally, before the finite-difference code; a "two electrons on a 3-point grid" count that makes $10^{30}$ tangible.
- **Suggested worked examples.** A fully worked normalisation and a fully worked orthogonality check (these are exactly the §4.8 exercise themes — promote one of each into the body as a solved model).
- **Suggested code-lab improvements.** Keep the particle-in-a-box and harmonic-oscillator labs (pure SciPy, browser-runnable) but precede each with a "what should the answer roughly look like?" prediction box; add a double-well teaser that foreshadows bonding.
- **Suggested "Common Misunderstandings" topics.** "$\psi$ is a physical wave you could measure"; "the electron is spread out like jelly"; confusing eigenvalue with eigenfunction; thinking Born–Oppenheimer is exact; "Hartree–Fock is just DFT".
- **Suggested exercise improvements.** Add pause-and-recall boxes after the Born rule and after the exponential-wall estimate; provide a hint ladder (`??? note "Hint"`) on the harder Hermite-polynomial problem.
- **Priority: High.** Everything from Chapter 5 onward is an approximation to the equation set up here; if the reader does not *feel* why the many-electron problem is hopeless, DFT looks like an arbitrary trick rather than a necessity. Core on-ramp.

---

## Chapter 5 — [Density Functional Theory](../ch05-dft/index.md)

- **Likely difficulty for undergraduates:** High. Abstract (a functional of a density), proof-heavy (Hohenberg–Kohn), and acronym-dense (LDA/GGA/SCF/DIIS).
- **Where it probably jumps too quickly.** The Hohenberg–Kohn proofs (§5.2) for a reader new to proof by contradiction; the conceptual switch from wavefunction to density; the Jacob's-ladder zoo of functionals (§5.4).
- **Hidden prerequisites.** The many-electron problem and the exponential wall (Ch 4); Bloch's theorem and reciprocal space (Ch 3.5); the variational principle and the idea of a functional (Ch 0); self-consistency as a fixed-point idea.
- **Concepts needing slower explanation.** Why "the density determines everything" is surprising and what it does *not* claim; what the Kohn–Sham orbitals *are not* (not the real electrons); where the unknown $E_{xc}[n]$ hides; why the SCF loop is even necessary.
- **Suggested toy examples.** A one-page "energy as a functional of a 1-number density" caricature before the real $E[n]$; a hand-iterated 2-step fixed-point to motivate SCF mixing.
- **Suggested worked examples.** A guided walk-through of the existing 1D model-hydrogen-chain SCF program (§5.5), annotated line by line, as the pilot already begins to do.
- **Suggested code-lab improvements.** The 1D SCF program is NumPy/SciPy and runs live — add a "watch it oscillate without mixing, then converge with mixing" before/after cell so the reader sees *why* DIIS exists.
- **Suggested "Common Misunderstandings" topics.** "KS eigenvalues are the real electron energies"; "the KS band gap is the true band gap"; "DFT is exact" vs. "DFT is fitted"; "LDA assumes the material is a uniform gas". (Several of these are already drafted in the pilot.)
- **Suggested exercise improvements.** Keep the eight problems; add pause-and-recall boxes after each Hohenberg–Kohn theorem and after the Kohn–Sham mapping, in the pilot style.
- **Priority: High.** The conceptual heart of the electronic-structure half of the book and the chosen pilot; finishing its scaffolding both helps readers most and supplies the worked pattern every other chapter copies.

---

## Chapter 6 — [Running DFT in Practice](../ch06-running-dft/index.md)

- **Likely difficulty for undergraduates:** Medium to High. The physics is lighter than Chapter 5, but the *operational* load (installs, input files, convergence sweeps) is where self-study readers stall and quietly give up.
- **Where it probably jumps too quickly.** The number of knobs introduced at once in §6.1 (cutoff, $k$-grid, pseudopotential, smearing, magnetisation); the scf $\to$ nscf $\to$ bands pipeline as three steps that look interchangeable but are not.
- **Hidden prerequisites.** A Linux/WSL2 environment and shell comfort (Ch 1 and Appendix B); $k$-point sampling and band structure (Ch 3.5); what SCF convergence means (Ch 5); patience with non-Python tools.
- **Concepts needing slower explanation.** What a pseudopotential is actually replacing and why that is allowed; why convergence testing is non-negotiable rather than busywork; reading a QE output file as a document with structure.
- **Suggested toy examples.** A "convergence in one variable" cartoon — total energy vs. cutoff on a handful of points — plotted before the full silicon sweep, so the shape is expected.
- **Suggested worked examples.** A fully annotated silicon `pw.x` input file with every line explained in a side margin; a worked vacancy-formation-energy bookkeeping calculation (the §6.5 arithmetic shown in full).
- **Suggested code-lab improvements.** QE and ASE do **not** run in the browser — say so prominently. Add a Colab-first, install-second path; provide a pre-computed convergence dataset so a reader without QE installed can still do the *plotting and judging* in live NumPy/Matplotlib.
- **Suggested "Common Misunderstandings" topics.** "Defaults are fine"; "more $k$-points/higher cutoff is always better" (cost!); "the calculation converged, so the answer is right"; confusing SCF convergence with $k$-point/cutoff convergence.
- **Suggested exercise improvements.** Add a checklist-style "before you trust this number" self-audit; provide a hint ladder for the defect-supercell problem.
- **Priority: High.** This is where theory becomes a runnable skill and where the on-ramp either succeeds or strands the reader. The browser-can't-run-it reality makes clear signposting and a no-install fallback essential.

---

## Chapter 7 — [Molecular Dynamics](../ch07-md/index.md)

- **Likely difficulty for undergraduates:** Medium. Newtonian intuition is friendly; thermostats and PBC are where it gets subtle.
- **Where it probably jumps too quickly.** Why the Verlet integrator is preferred over naive Euler (symplecticness stated, not felt); thermostats as more than "a knob that sets temperature"; the minimum-image convention in PBC.
- **Hidden prerequisites.** ODE intuition and finite-difference stepping (Ch 0); NumPy (Ch 1); phonons help but are not required (Ch 3.5); the idea of a potential energy surface (Ch 4).
- **Concepts needing slower explanation.** Energy drift and what conserving energy means numerically; temperature as a property of an *ensemble*, not an atom; periodic images as copies, not walls.
- **Suggested toy examples.** A 1D harmonic oscillator integrated with Euler (drifts) vs. velocity-Verlet (stable), plotted, before any many-atom MD; two particles in a periodic box to show minimum image.
- **Suggested worked examples.** A worked single velocity-Verlet step by hand for one particle; a worked instantaneous-temperature calculation from a handful of velocities.
- **Suggested code-lab improvements.** The integrator comparison is pure NumPy and runs live — make it the opening lab. LAMMPS/ASE workflows do not run in-browser; flag that and offer a small browser-runnable Lennard-Jones MD as the conceptual stand-in.
- **Suggested "Common Misunderstandings" topics.** "A single atom has a temperature"; "any integrator will do if the timestep is small"; "PBC means the box has walls"; confusing thermostat-set temperature with measured temperature.
- **Suggested exercise improvements.** Add a pause-and-recall after the energy-drift demonstration; add a guided problem on choosing a stable timestep.
- **Priority: Medium.** Important and widely used, but more forgiving than the electronic-structure chapters; the toy-integrator demo is high value for low effort.

---

## Chapter 8 — [Statistical Mechanics from Simulations](../ch08-statmech/index.md)

- **Likely difficulty for undergraduates:** Medium to High. Abstract (ensembles, partition functions) and easy to learn as formulae without meaning.
- **Where it probably jumps too quickly.** The partition function appearing as a definition with no felt purpose; free-energy methods (thermodynamic integration) as machinery before the reader sees *why* free energy is hard to get directly.
- **Hidden prerequisites.** Probability and expectation values (Ch 0); MD trajectories as samples (Ch 7); the idea that a simulation samples a distribution; logarithms and exponentials with physical arguments.
- **Concepts needing slower explanation.** What an ensemble *is* (a probability distribution over microstates); why you cannot just read free energy off a trajectory; ergodicity as "averaging over time equals averaging over the ensemble".
- **Suggested toy examples.** A two-state system (just two energies) for which the partition function, average energy and free energy are computed by hand, then re-derived from a tiny Monte Carlo sample.
- **Suggested worked examples.** A worked thermodynamic-integration toy on a 1D model where the answer is known analytically, so the method can be checked.
- **Suggested code-lab improvements.** A live NumPy Monte Carlo of the two-state (or 2D Ising) system showing how a sampled average approaches the analytic one; pure NumPy, browser-runnable.
- **Suggested "Common Misunderstandings" topics.** "Free energy is just energy"; "a longer run is automatically converged"; confusing time average with ensemble average; "entropy is disorder" used as if it were a calculation.
- **Suggested exercise improvements.** Add a pause-and-recall after the partition-function definition; add a graded ladder from a two-state system up to a small Ising model.
- **Priority: Medium.** Conceptually demanding, but its hardest dependencies (Chs 0, 7) are upstream; scaffolding the two-state toy removes most of the fog cheaply.

---

## Chapter 9 — [Machine Learning Interatomic Potentials](../ch09-mlip/index.md)

- **Likely difficulty for undergraduates:** High. At ~20 hours it is the longest chapter and braids physics (symmetries, descriptors) with machine-learning machinery (NequIP, MACE).
- **Where it probably jumps too quickly.** The symmetry requirements (invariance vs. equivariance) stated abstractly; the descriptor zoo (SOAP, ACE) before the reader knows what a descriptor must achieve; the jump to equivariant networks in §9.5.
- **Hidden prerequisites.** Basic ML vocabulary (training set, loss, overfitting) that the book may assume; descriptors as fixed-length vectors (Ch 3 reciprocal-space habits help); where the training data comes from (Ch 6); what the potential is *for* (Ch 7).
- **Concepts needing slower explanation.** Why translational/rotational/permutation symmetry must be *built in* rather than learnt; the difference between an invariant energy and equivariant forces; what "local environment" means physically.
- **Suggested toy examples.** A 2-atom-in-a-box example showing that a naive coordinate list is not rotation-invariant, then fixing it with a distance; a hand-built tiny descriptor for a 3-atom cluster.
- **Suggested worked examples.** A worked symmetry check: take a small structure, rotate it, and verify the energy is unchanged but forces rotate with it.
- **Suggested code-lab improvements.** PyTorch/ASE training does **not** run in the browser — flag it clearly. Provide a tiny pure-NumPy linear-regression-on-distances "potential" the reader *can* run live, as the conceptual seed before MACE.
- **Suggested "Common Misunderstandings" topics.** "MLIPs are exact because they fit DFT"; "more parameters always means better"; confusing invariance with equivariance; "the potential will extrapolate to structures it never saw".
- **Suggested exercise improvements.** Add a pause-and-recall after the symmetry section; add a low-stakes "predict the failure mode" exercise on extrapolation.
- **Priority: Medium.** Heavily used at the frontier, but it sits past the electronic-structure on-ramp; readers who reach it are already committed. High *length* makes a milestone-by-milestone split (symmetries first) worthwhile.

---

## Chapter 10 — [Graph Neural Networks for Materials](../ch10-gnn/index.md)

- **Likely difficulty for undergraduates:** Medium to High. The graph idea is intuitive; message passing and the model-evolution survey are where readers get lost.
- **Where it probably jumps too quickly.** Message passing presented as equations before the "atoms pass notes to neighbours" picture is fixed; the CGCNN-from-scratch build assuming PyTorch fluency.
- **Hidden prerequisites.** The descriptor concept from Chapter 9; basic neural-network vocabulary (layer, weights, training); what a graph is (nodes, edges); crystals-as-periodic-graphs needing PBC intuition (Ch 7).
- **Concepts needing slower explanation.** What a "message" is and why aggregation must be permutation-invariant; how a crystal becomes a graph despite being infinite; why a GNN can respect translational symmetry "for free".
- **Suggested toy examples.** A 3-node graph with hand-computed one round of message passing; a tiny molecule drawn as a graph before any crystal.
- **Suggested worked examples.** A worked single message-passing update on a small graph, every sum shown, mirroring the CGCNN update rule.
- **Suggested code-lab improvements.** The CGCNN build needs PyTorch (not browser-runnable — say so). Add a pure-NumPy "one message-passing layer by hand" cell that *is* runnable live, as the conceptual precursor.
- **Suggested "Common Misunderstandings" topics.** "A GNN sees the whole crystal at once"; "more layers always helps"; confusing node features with edge features; "the model understands chemistry".
- **Suggested exercise improvements.** Add a pause-and-recall after the message-passing definition; add a "draw this structure as a graph" hand exercise.
- **Priority: Low to Medium.** Specialised and downstream; valuable scaffolding is the by-hand message-passing toy, but few self-study readers reach this depth, so it ranks below the on-ramp chapters.

---

## Chapter 11 — [Active Learning and Bayesian Optimisation](../ch11-active/index.md)

- **Likely difficulty for undergraduates:** Medium to High. Gaussian processes are notoriously hard to meet for the first time.
- **Where it probably jumps too quickly.** The leap from "fit a curve" to "a distribution over functions"; acquisition functions presented as formulae before the exploration–exploitation tension is felt.
- **Hidden prerequisites.** Probability, Gaussians, and covariance (Ch 0); the idea of uncertainty as a quantity you can compute; sampling intuition (Ch 8); what an MLIP costs to evaluate (Ch 9).
- **Concepts needing slower explanation.** What a Gaussian process *is* (a prior over functions, narrowing as data arrive); why uncertainty, not just the mean, drives the next experiment; the exploration–exploitation trade-off as a decision, not a formula.
- **Suggested toy examples.** A 1D GP fit to three points with the uncertainty band drawn, updated as points are added; a one-dimensional "where do I sample next?" decision made by hand from an acquisition curve.
- **Suggested worked examples.** A worked single BO step on a 1D function with a known optimum, so the reader sees the loop close.
- **Suggested code-lab improvements.** A 1D GP and an expected-improvement acquisition in pure NumPy/SciPy — browser-runnable — as the centrepiece, before any materials-scale example.
- **Suggested "Common Misunderstandings" topics.** "BO is just grid search"; "the GP mean is the answer"; "uncertainty is error"; confusing exploration with randomness.
- **Suggested exercise improvements.** Add a pause-and-recall after the GP definition; add a guided "choose the next point" exercise from a plotted acquisition function.
- **Priority: Medium.** The 1D GP visualisation is unusually high-leverage scaffolding, but the chapter is Tier 2 and optional for many readers.

---

## Chapter 12 — [Foundation Models for Materials](../ch12-foundation/index.md)

- **Likely difficulty for undergraduates:** Medium. More conceptual and current-affairs than mathematical; the risk is hype and unfamiliar named systems, not hard derivations.
- **Where it probably jumps too quickly.** The paradigm shift ("pretrain once, use everywhere") assumed rather than motivated; named systems (MACE-MP-0, MatterGen) arriving as a list.
- **Hidden prerequisites.** MLIPs and GNNs (Chs 9, 10) as the thing being scaled up; what "pretraining" and "transfer" mean; honest expectations about reliability.
- **Concepts needing slower explanation.** What "foundation model" means here vs. in language models; what a *universal* potential does and does not guarantee; what a generative materials model is actually proposing.
- **Suggested toy examples.** A before/after sketch: a potential trained per-system vs. one pretrained across many — what changes for the user.
- **Suggested worked examples.** Not a derivation chapter; instead a worked "sanity-check a universal-potential prediction" walk-through (compare to a known value).
- **Suggested code-lab improvements.** These models do not run in the browser and often need real downloads/compute — state this plainly and keep code as clearly-labelled illustrative snippets, not live cells.
- **Suggested "Common Misunderstandings" topics.** "A universal potential works for everything"; "generated structures are guaranteed stable"; "foundation model = AGI for materials"; treating a single prediction as ground truth.
- **Suggested exercise improvements.** Add reflective "when would you trust / not trust a universal model?" prompts with model answers; emphasise validation habits over computation.
- **Priority: Low.** Fast-moving and downstream; the best scaffolding is calibration of expectations and cross-links to validation, which are cheap to add.

---

## Chapter 13 — [Multiscale Methods](../ch13-multiscale/index.md)

- **Likely difficulty for undergraduates:** Medium to High. A survey of several distinct couplings (QM/MM, coarse-graining, kMC, phase-field, finite element), each with its own vocabulary.
- **Where it probably jumps too quickly.** Switching methods between subsections without a shared "what scale, what is the handoff?" frame; the finite-element bridge assuming continuum-mechanics comfort.
- **Hidden prerequisites.** DFT (Ch 5) and MD (Ch 7) as the two scales being bridged; the scale-ladder picture from Chapter 2; basic continuum/PDE intuition for the FE section.
- **Concepts needing slower explanation.** What "coupling" physically means at the boundary between two methods; why information is lost going up scales (coarse-graining); kMC as event-driven rather than time-stepped.
- **Suggested toy examples.** A 1D "fine region / coarse region" cartoon showing the handoff at the interface; a tiny kMC of one hopping atom to contrast with MD's fixed timestep.
- **Suggested worked examples.** A worked QM/MM energy-partitioning bookkeeping on a toy system, so the double-counting pitfall is visible.
- **Suggested code-lab improvements.** Most real multiscale tooling is heavyweight; keep live code to a small NumPy kMC of a single process (browser-runnable) and flag the rest as illustrative.
- **Suggested "Common Misunderstandings" topics.** "Multiscale = just a bigger simulation"; "coarse-graining loses nothing important"; confusing kMC time with MD time; QM/MM double-counting at the boundary.
- **Suggested exercise improvements.** Add a "which scale owns this question?" sorting exercise; add a pause-and-recall after the coupling-schemes overview.
- **Priority: Low.** Treated as a reference chapter even in the existing learning path; scaffolding helps but reaches few readers, so it ranks last alongside Chapter 12.

---

## Chapter 14 — [Designing Your Own Project](../ch14-capstone/index.md)

- **Likely difficulty for undergraduates:** Medium. Low on derivation, high on judgement and open-endedness — which is its own kind of hard for a beginner used to set problems.
- **Where it probably jumps too quickly.** The assumption that the reader can scope a tractable question; the convergence-and-validation material relies on every Tier-1 chapter at once.
- **Hidden prerequisites.** Working familiarity with at least one method (Chs 5–11); the convergence habits from Chapter 6; comfort reading a paper (a skill rarely taught explicitly).
- **Concepts needing slower explanation.** What makes a question *tractable* in a semester; what "validation" means in practice; how to tell a result from an artefact.
- **Suggested toy examples.** A worked "scoping" example that takes one over-ambitious idea and narrows it to something doable, shown as a before/after.
- **Suggested worked examples.** A short worked validation story — a plausible-but-wrong result caught by a convergence check — drawn from the existing common-pitfalls material.
- **Suggested code-lab improvements.** Mostly project-specific; add a reusable "project sanity-check" checklist (not code) and cross-link the five projects in the nav.
- **Suggested "Common Misunderstandings" topics.** "A bigger project is a better project"; "if it ran, it is correct"; "I need a novel idea before I can start"; underestimating validation time.
- **Suggested exercise improvements.** The existing open-ended exercises are appropriate; add scaffolded scoping prompts with model answers and a self-audit checklist, cross-linked to [undergraduate-projects.md](undergraduate-projects.md).
- **Priority: Low to Medium.** The capstone benefits most *after* the on-ramp chapters are scaffolded; its judgement-building additions are valuable but depend on the earlier work being in place.

---

## Suggested order of work

The cheapest path to the biggest improvement is to scaffold the electronic-structure on-ramp first, building outward from the existing pilot, then handle the foundation, then the rest by dependency order:

1. **Chapter 5 (DFT)** — already piloted; finish and polish it as the reference implementation of the template, so every later chapter has a worked pattern to copy.
2. **Chapter 4 (Quantum Mechanics)** — the equation Chapter 5 approximates; scaffold the exponential wall and the wavefunction concept.
3. **Chapter 3.5 (Solid State Physics)** — fix $\mathbf{k}$, bands and $k$-point intuition before they bite in Chapters 5 and 6.
4. **Chapter 6 (Running DFT)** — turn theory into a runnable skill; add the no-install fallback.
5. **Chapter 0 (Mathematics)** — back-fill the eigen/Fourier/functional intuition the four chapters above lean on; do it once the downstream chapters have shown exactly which gaps matter.
6. **Chapters 7, 8, 9** (Medium priority) — MD, statistical mechanics, MLIPs, in dependency order.
7. **Chapters 1, 2, 3** (light passes) — confidence-building and cross-linking; cheap and low-risk.
8. **Chapters 10–14** (Low priority) — the GNN/foundation/multiscale/capstone tail, scaffolded last, as few self-study readers reach them and they depend on everything above.

Throughout, apply the additions in place using the [chapter expansion template](chapter-expansion-template.md), keep the pause-and-recall / check-yourself boxes consistent with the **Chapter 5 pilot**, and respect the honest study-time estimates in the [learning path](../learning-path.md) — scaffolding lowers the activation energy, it does not shorten the chapter.
