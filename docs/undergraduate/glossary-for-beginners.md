# A Beginner's Glossary

This page collects the words and symbols that appear again and again across the handbook and
explains each one slowly, in ordinary language, before showing the formal version. It is
deliberately patient: it would rather use three short sentences than one dense one.

!!! info "How this differs from Appendix C"
    The book already has a terse **Appendix C glossary** meant for quick reference once you
    know the field. This page is the *slow, friendly* companion to it — same terms, more
    hand-holding. When a definition here feels long, that is on purpose. For a one-line
    reminder, use Appendix C; to look something up fast, use the search box at the top of the
    site (it indexes every page).

Each entry follows the same six-part shape so you always know where to look:
**Plain meaning**, **In this handbook**, **In equations**, **In code or input files**,
**Why it matters**, and one **Common misunderstanding** to watch for. The terms are grouped by
topic rather than alphabetised, so that related ideas sit together. New to the handbook? The
[learning paths](learning-paths.md) page suggests an order; the
[formula reading guide](formula-reading-guide.md) and
[code reading guide](code-reading-guide.md) go deeper on equations and code respectively.

---

## Structure and crystals

### Atom

**Plain meaning:** The basic building block of ordinary matter — a small positive nucleus
surrounded by electrons. In simulation we usually treat it as a labelled point in space with a
chemical identity (e.g. "this is a silicon atom here").

**In this handbook:** Atoms first appear properly in [Chapter 3 (atoms and structure)](../ch03-atoms/index.md)
and underlie everything afterwards, from quantum calculations to molecular dynamics.

**In equations:** Often indexed by a subscript: atom $i$ sits at position $\mathbf{r}_i$ and
carries an atomic number $Z_i$. A whole structure is then the set $\{\mathbf{r}_i\}$.

**In code or input files:** An atom is typically a row of data: an element symbol plus three
coordinates, for example `("Si", [0.0, 0.0, 0.0])`. Libraries like ASE store collections of
these in an `Atoms` object (ASE does not run live in the browser).

**Why it matters:** Almost every method in the book ultimately computes the energy and forces
for a particular arrangement of atoms, so "where are the atoms" is the central input.

**Common misunderstanding:** Beginners sometimes picture electrons orbiting like planets. In
quantum simulation the electrons are described by a smeared-out probability cloud, not tiny
particles on fixed tracks.

### Bond

**Plain meaning:** A bond is a region of strong attraction holding two atoms together because
sharing or transferring electrons lowers their combined energy. It is a useful picture, not a
physical stick.

**In this handbook:** Bonds appear when we discuss molecules and structure in
[Chapter 3 (atoms and structure)](../ch03-atoms/index.md), and reappear as graph "edges" in
[Chapter 10 (graph neural networks)](../ch10-gnn/index.md).

**In equations:** In simple force fields a bond contributes a term like
$\tfrac{1}{2}k\,(r - r_0)^2$, where $r$ is the current distance and $r_0$ the preferred bond
length — a spring, in other words.

**In code or input files:** Bonds may be listed explicitly as pairs of atom indices, e.g.
`bonds = [(0, 1), (1, 2)]`, or inferred from distances by a cutoff.

**Why it matters:** Whether two atoms are "bonded" controls how a material deforms, reacts and
vibrates, and it is how chemists reason about structure.

**Common misunderstanding:** A bond is not an all-or-nothing object. In quantum methods (DFT)
there is no explicit "bond" variable at all — bonding emerges from the electron density.

### Crystal

**Plain meaning:** A crystal is a solid whose atoms repeat in a regular, periodic pattern, like
a 3D wallpaper that tiles space forever.

**In this handbook:** Crystals are the main subject of
[Chapter 3b (solid-state)](../ch03b-solid-state/index.md) and are assumed by most DFT
calculations in [Chapter 5 (DFT)](../ch05-dft/index.md).

**In equations:** Crystallinity is expressed by translational symmetry: the structure looks
identical after a shift by any lattice vector $\mathbf{R} = n_1\mathbf{a}_1 + n_2\mathbf{a}_2 +
n_3\mathbf{a}_3$ with integer $n_i$.

**In code or input files:** A crystal is usually given as a unit cell plus a list of atoms, e.g.
a `cell` matrix and fractional coordinates in a `.cif` or POSCAR file.

**Why it matters:** Periodicity lets us simulate an effectively infinite solid using only a
small repeating piece, which makes solid-state calculations affordable.

**Common misunderstanding:** "Crystal" does not mean transparent or gem-like. A block of copper
is a crystal; glass, despite looking solid, is *not* (it lacks long-range order).

### Unit cell

**Plain meaning:** The smallest tile that, when stacked in all directions, builds the whole
crystal. Choose it once and you have described an infinite solid.

**In this handbook:** Introduced in [Chapter 3b (solid-state)](../ch03b-solid-state/index.md);
every periodic DFT input in [Chapter 6 (running DFT)](../ch06-running-dft/index.md) starts from
one.

**In equations:** Defined by three edge vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$
collected into a cell matrix; its volume is the scalar triple product
$V = \mathbf{a}_1 \cdot (\mathbf{a}_2 \times \mathbf{a}_3)$.

**In code or input files:** Often a $3\times 3$ array, e.g.
`cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]` for a cubic cell of side `a`.

**Why it matters:** The cell sets the periodic boundary and the volume, and changing it (strain)
is how you study pressure, elasticity and phase changes.

**Common misunderstanding:** The unit cell is not unique — many valid choices describe the same
crystal. The *conventional* cell (often more symmetric) and the *primitive* cell (smallest) can
differ in size.

### Lattice

**Plain meaning:** The bare grid of repeat points — an infinite set of mathematical dots marking
where the pattern repeats, before you put any atoms on them.

**In this handbook:** Discussed in [Chapter 3b (solid-state)](../ch03b-solid-state/index.md)
alongside the unit cell and basis.

**In equations:** The set of lattice vectors $\mathbf{R} = \sum_i n_i \mathbf{a}_i$ with integer
$n_i$. The lattice is the points; the atoms are added separately.

**In code or input files:** Represented by the same cell vectors as the unit cell; the lattice
itself carries no atoms, so it is usually implicit in the `cell` definition.

**Why it matters:** Separating "the grid" (lattice) from "what sits on it" (basis) is the clean
way to describe any crystal, however complicated.

**Common misunderstanding:** The lattice is not the same as the atoms. A lattice point need not
have an atom exactly on it, and one lattice point can correspond to several atoms.

### Basis

**Plain meaning:** The group of atoms attached to each lattice point — the actual "stuff" you
stamp down at every repeat position.

**In this handbook:** Appears with lattice and unit cell in
[Chapter 3b (solid-state)](../ch03b-solid-state/index.md).

**In equations:** Atom positions are $\mathbf{r} = \mathbf{R} + \boldsymbol{\tau}_j$, where
$\mathbf{R}$ is a lattice vector and $\boldsymbol{\tau}_j$ are the basis offsets within one cell.

**In code or input files:** The basis is the list of fractional coordinates inside the cell, e.g.
`positions = [[0,0,0], [0.25,0.25,0.25]]` for a two-atom basis.

**Why it matters:** Crystal = lattice + basis. Without the basis you have an empty grid; the
basis is what makes it diamond, salt or quartz.

**Common misunderstanding:** "Basis" here is unrelated to the "basis set" used in DFT (the
functions that build the wavefunction). Same word, different meaning — context tells you which.

### Reciprocal space

**Plain meaning:** A second coordinate system, paired with ordinary ("real") space, in which
periodic patterns are easier to describe. It is the natural home of waves and repeating
structures.

**In this handbook:** Introduced in [Chapter 3b (solid-state)](../ch03b-solid-state/index.md)
and used throughout DFT in [Chapter 5 (DFT)](../ch05-dft/index.md).

**In equations:** Reciprocal space is the *Fourier dual* of real space. Its basis vectors
$\mathbf{b}_i$ satisfy $\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi\,\delta_{ij}$, and points in it
are wavevectors $\mathbf{k}$ with units of inverse length.

**In code or input files:** Rarely typed directly, but it underlies fast Fourier transforms
(`numpy.fft`) and the k-point grids you specify for a calculation.

**Why it matters:** In a crystal, the physics simplifies enormously when written in reciprocal
space; band structures and plane-wave DFT live there.

**Common misunderstanding:** "Reciprocal" does not just mean "one over the distance" for a
single number. It is a whole dual space obtained by Fourier transform, with its own vectors and
geometry.

### Brillouin zone

**Plain meaning:** A specially chosen region of reciprocal space that contains every distinct
wave a crystal can support — once you know what happens inside it, you know everything
everywhere by symmetry.

**In this handbook:** Defined in [Chapter 3b (solid-state)](../ch03b-solid-state/index.md) and
sampled when computing band structures in [Chapter 5 (DFT)](../ch05-dft/index.md).

**In equations:** It is the primitive cell of the reciprocal lattice (specifically the
Wigner–Seitz cell), so any wavevector can be reduced into it: $\mathbf{k} \to \mathbf{k} +
\mathbf{G}$ for some reciprocal lattice vector $\mathbf{G}$.

**In code or input files:** You meet it through high-symmetry point labels along a band path,
e.g. `Gamma -> X -> M -> Gamma`, and through k-point density settings.

**Why it matters:** Electronic and vibrational properties of solids are catalogued by where they
sit in the Brillouin zone; the band gap, for instance, is a feature of specific points.

**Common misunderstanding:** It is not a region of real space and has nothing to do with where
the atoms physically are; it lives entirely in reciprocal (wavevector) space.

---

## Quantum mechanics

### Wavefunction

**Plain meaning:** A mathematical object that holds everything knowable about a quantum system.
For one electron it is a function whose squared size tells you how likely the electron is to be
found at each point.

**In this handbook:** Central to [Chapter 4 (quantum mechanics)](../ch04-quantum/index.md) and
the foundation for [Chapter 5 (DFT)](../ch05-dft/index.md).

**In equations:** Written $\psi(\mathbf{r})$ or as a state $|\psi\rangle$. The probability
density of finding the particle at $\mathbf{r}$ is $|\psi(\mathbf{r})|^2$, and it integrates to
one: $\int |\psi|^2\,\mathrm{d}\mathbf{r} = 1$.

**In code or input files:** Usually a complex-valued array over a spatial grid, e.g.
`psi = np.zeros(N, dtype=complex)`, with probabilities given by `np.abs(psi)**2`.

**Why it matters:** Solving for the wavefunction (or, in DFT, for the density that stands in for
it) is what "doing the quantum mechanics" means in practice.

**Common misunderstanding:** $\psi$ itself is not a probability and can be negative or complex.
Only $|\psi|^2$ is a probability density.

### Operator

**Plain meaning:** A rule that takes a function and gives back another function — a verb for
quantum mechanics. "Differentiate", "multiply by position" and "measure the energy" are all
operators.

**In this handbook:** Introduced in [Chapter 4 (quantum mechanics)](../ch04-quantum/index.md);
the most important one is the Hamiltonian.

**In equations:** Operators wear hats: $\hat{A}$, $\hat{T}$ (kinetic), $\hat{H}$ (energy). They
act on states: $\hat{A}\,\psi(\mathbf{r})$ is a new function.

**In code or input files:** On a grid, an operator often becomes a matrix, so applying it is a
matrix–vector product, e.g. `H @ psi` in NumPy.

**Why it matters:** Every physical quantity you can measure corresponds to an operator, so
operators are how questions ("what is the energy?") get encoded mathematically.

**Common misunderstanding:** Operators usually do *not* commute: $\hat{A}\hat{B}$ and
$\hat{B}\hat{A}$ can differ. Order matters, unlike ordinary number multiplication.

### Hamiltonian

**Plain meaning:** The operator that represents the total energy of a system — kinetic plus
potential. It is the master object whose solutions describe the system's allowed states.

**In this handbook:** Appears throughout [Chapter 4 (quantum mechanics)](../ch04-quantum/index.md)
and again in [Chapter 5 (DFT)](../ch05-dft/index.md), where it includes the electron–electron
interaction.

**In equations:** For one particle, $\hat{H} = \hat{T} + \hat{V} = -\dfrac{\hbar^2}{2m}\nabla^2 +
V(\mathbf{r})$. It enters the central eigenvalue equation $\hat{H}\psi = E\psi$.

**In code or input files:** Built as a matrix combining a kinetic part (a discretised Laplacian)
and a potential (a diagonal array), e.g. `H = T + np.diag(V)`.

**Why it matters:** Almost every quantum calculation reduces to "construct $\hat{H}$, then find
its eigenvalues and eigenvectors". It is the heart of the method.

**Common misunderstanding:** The Hamiltonian is the energy *operator*, not a single energy
number. The numbers (energies) come out only after you solve its eigenvalue equation.

### Eigenvalue

**Plain meaning:** When an operator acts on certain special functions, it simply rescales them
without changing their shape. The scale factor is the eigenvalue. In quantum mechanics these are
the measurable values — for $\hat{H}$, the allowed energies.

**In this handbook:** Introduced in [Chapter 0 (maths)](../ch00-math/index.md) as linear
algebra, then given physical meaning in [Chapter 4 (quantum mechanics)](../ch04-quantum/index.md).

**In equations:** The eigenvalue equation is $\hat{A}\psi = a\psi$, where $a$ is the eigenvalue
(a number) and $\psi$ the eigenvector. For energies: $\hat{H}\psi = E\psi$.

**In code or input files:** Found numerically with one call, e.g.
`evals, evecs = np.linalg.eigh(H)`; `evals` are the eigenvalues, returned in ascending order.

**Why it matters:** Quantised energy levels, band energies and vibrational frequencies are all
eigenvalues; computing them is a recurring task across the book.

**Common misunderstanding:** An eigenvalue is not "the answer" on its own — it always comes
paired with its eigenvector, and the same operator has many eigenvalue–eigenvector pairs.

### Eigenvector

**Plain meaning:** The special function (or vector) that an operator only stretches or shrinks,
never reshapes. In quantum mechanics it is the state that goes with a given measured value.

**In this handbook:** Paired with eigenvalues in [Chapter 0 (maths)](../ch00-math/index.md) and
[Chapter 4 (quantum mechanics)](../ch04-quantum/index.md); the eigenvectors of $\hat{H}$ are the
allowed wavefunctions.

**In equations:** In $\hat{A}\psi = a\psi$, the $\psi$ is the eigenvector (often called an
*eigenstate* or *eigenfunction* in this context).

**In code or input files:** Returned alongside eigenvalues, e.g. column `evecs[:, 0]` is the
eigenvector for the lowest eigenvalue `evals[0]`.

**Why it matters:** Eigenvectors are the actual quantum states — orbitals, bands and normal
modes — so they tell you *what* the system is doing, while eigenvalues tell you the energy cost.

**Common misunderstanding:** Eigenvectors are only defined up to scale (and a phase), so a
solver may return one pointing the "opposite" way; that is the same physical state.

### Boundary condition

**Plain meaning:** A rule about what happens at the edges of your simulated region — does the
wavefunction vanish there, or does the box wrap around so that leaving one side re-enters the
other?

**In this handbook:** Discussed in [Chapter 4 (quantum mechanics)](../ch04-quantum/index.md) and
crucially in [Chapter 3b (solid-state)](../ch03b-solid-state/index.md), where periodic boundary
conditions make a small cell behave like an infinite crystal.

**In equations:** For a particle in a box, $\psi(0) = \psi(L) = 0$; for a crystal, the periodic
(Bloch) condition $\psi(\mathbf{r} + \mathbf{R}) = e^{i\mathbf{k}\cdot\mathbf{R}}\,\psi(\mathbf{r})$.

**In code or input files:** Often a flag, e.g. `pbc=True` for periodic boundaries in an
atomistic simulation, or a choice of how the grid edges connect.

**Why it matters:** The boundary condition changes the allowed solutions completely; the same
Hamiltonian gives discrete levels in a box but continuous bands under periodicity.

**Common misunderstanding:** Periodic boundaries do not mean the system is truly infinite — it
is one cell repeated, so artefacts appear if the cell is too small for the physics you want.

---

## DFT

### Density

**Plain meaning:** The electron density is how much electron "charge cloud" sits at each point in
space. Add it up everywhere and you get the total number of electrons.

**In this handbook:** The starring quantity of [Chapter 5 (DFT)](../ch05-dft/index.md) — "DFT"
literally means *density* functional theory.

**In equations:** Written $n(\mathbf{r})$ or $\rho(\mathbf{r})$, built from the occupied orbitals
as $n(\mathbf{r}) = \sum_i |\psi_i(\mathbf{r})|^2$, with $\int n(\mathbf{r})\,\mathrm{d}\mathbf{r}
= N$ electrons.

**In code or input files:** Stored as a value on a 3D grid, e.g. a `density` array of shape
`(nx, ny, nz)`, and written to files such as a charge-density output.

**Why it matters:** DFT's central claim is that this single 3D function — not the full
many-electron wavefunction — is enough to determine the ground-state energy, which is what makes
it affordable.

**Common misunderstanding:** The density is a *real, positive* scalar field, not a wavefunction;
it has no phase and you cannot recover individual orbitals from it directly.

### Functional

**Plain meaning:** A function takes a number and returns a number; a *functional* takes a whole
function and returns a single number. Think of it as a machine that scores an entire curve with
one value.

**In this handbook:** The key idea behind [Chapter 5 (DFT)](../ch05-dft/index.md): the energy is
a functional of the density.

**In equations:** Written with square brackets to show the input is a function: $E[n]$ means "the
energy, given the whole density function $n(\mathbf{r})$". An example is
$E[n] = \int n(\mathbf{r})\,v(\mathbf{r})\,\mathrm{d}\mathbf{r} + \dots$.

**In code or input files:** A functional becomes a routine that takes an array (the sampled
function) and returns a scalar, e.g. `def E(n): return np.sum(n * v) * dV`.

**Why it matters:** The entire framework rests on the existence of an energy functional of the
density; "which functional" is one of the main accuracy choices you make.

**Common misunderstanding:** A functional is not just a complicated function. The defining point
is that its *input is a function* and its output is a single number.

### Exchange–correlation

**Plain meaning:** A correction term that captures the subtle quantum effects of electrons
avoiding one another — effects too complicated to write down exactly, so they are bundled into
one approximate piece of the energy.

**In this handbook:** Discussed in [Chapter 5 (DFT)](../ch05-dft/index.md) as the part of the
functional we must approximate.

**In equations:** The exchange–correlation energy $E_{\mathrm{xc}}[n]$ is one term in the total
DFT energy $E[n] = T_s[n] + E_{\mathrm{H}}[n] + E_{\mathrm{xc}}[n] + \int v_{\mathrm{ext}}\,n\,
\mathrm{d}\mathbf{r}$.

**In code or input files:** Chosen by name in an input file, e.g. an `xc` or functional keyword
set to `PBE` or `LDA`.

**Why it matters:** It is the one term DFT cannot compute exactly, so the choice of
exchange–correlation approximation is the biggest source of systematic error in a calculation.

**Common misunderstanding:** "Exchange" and "correlation" are two distinct physical effects
bundled together; and no approximation is universally best — a functional good for one property
may be poor for another.

### Self-consistent field

**Plain meaning:** A chicken-and-egg loop: the electron density depends on the potential, but the
potential depends on the density. You guess, recompute, and repeat until the two stop changing —
then you have a *self-consistent* answer.

**In this handbook:** The core algorithm of [Chapter 5 (DFT)](../ch05-dft/index.md), run in
practice in [Chapter 6 (running DFT)](../ch06-running-dft/index.md).

**In equations:** You iterate $n_{\text{in}} \to V[n_{\text{in}}] \to \hat{H} \to \{\psi_i\} \to
n_{\text{out}}$ until $\lVert n_{\text{out}} - n_{\text{in}}\rVert$ falls below a tolerance.

**In code or input files:** Appears as convergence settings, e.g. `scf_tol = 1e-6` and a maximum
iteration count; the output prints the energy at each SCF step.

**Why it matters:** Whether your DFT result is trustworthy depends on the SCF loop having
genuinely converged; unconverged runs give meaningless energies.

**Common misunderstanding:** "Self-consistent" does not mean "correct" — it means internally
consistent. A converged calculation can still be wrong if the functional, basis or cell is poor.

### Pseudopotential

**Plain meaning:** A simplified stand-in for an atom's nucleus *and* its tightly bound inner
electrons, so the simulation only has to track the chemically active outer electrons.

**In this handbook:** Introduced in [Chapter 5 (DFT)](../ch05-dft/index.md) and selected in
practice in [Chapter 6 (running DFT)](../ch06-running-dft/index.md).

**In equations:** It replaces the steep true potential near the nucleus with a smooth
$V_{\mathrm{ps}}(\mathbf{r})$ that matches the real one outside a chosen core radius
$r_{\mathrm{c}}$.

**In code or input files:** Supplied as a per-element file, referenced like
`pseudo = {"Si": "Si.UPF"}` (such DFT codes do not run live in the browser).

**Why it matters:** Core electrons barely change during chemistry but are expensive to compute;
pseudopotentials remove that cost and make plane-wave DFT feasible.

**Common misunderstanding:** Pseudopotentials are not a single fixed choice — they come in
families and must be matched to your functional and cutoff, or results degrade.

### Plane-wave cutoff

**Plain meaning:** A setting that decides how finely the electron wavefunctions and density are
resolved by limiting how wiggly the building-block waves are allowed to be. Higher cutoff means
finer detail and more cost.

**In this handbook:** A key convergence parameter in
[Chapter 6 (running DFT)](../ch06-running-dft/index.md).

**In equations:** Only plane waves with kinetic energy below the cutoff are kept:
$\dfrac{\hbar^2}{2m}|\mathbf{k} + \mathbf{G}|^2 \le E_{\mathrm{cut}}$.

**In code or input files:** A single energy value, e.g. `ecutwfc = 40` (in rydberg or eV
depending on the code).

**Why it matters:** Too low a cutoff gives inaccurate energies and forces; you must test that
results stop changing as you raise it (a convergence test).

**Common misunderstanding:** A higher cutoff is not "more correct physics" — it just reduces one
numerical error. Beyond convergence it only wastes compute time.

### k-point

**Plain meaning:** A sampling point in reciprocal space. Because a crystal's properties vary
smoothly across the Brillouin zone, you approximate sums over it by evaluating at a finite grid
of k-points.

**In this handbook:** Introduced with reciprocal space in
[Chapter 3b (solid-state)](../ch03b-solid-state/index.md) and chosen in practice in
[Chapter 6 (running DFT)](../ch06-running-dft/index.md).

**In equations:** A wavevector $\mathbf{k}$ in the Brillouin zone; integrals become weighted sums
$\int_{\mathrm{BZ}} f(\mathbf{k})\,\mathrm{d}\mathbf{k} \approx \sum_{\mathbf{k}} w_{\mathbf{k}}
f(\mathbf{k})$.

**In code or input files:** Specified as a grid, e.g. `kpoints = (8, 8, 8)` for a
Monkhorst–Pack mesh.

**Why it matters:** Too few k-points gives wrong energies and metallic properties; like the
cutoff, k-point density needs a convergence test.

**Common misunderstanding:** A denser k-grid is needed for metals than for insulators, and the
required density goes *down* as the cell gets bigger — bigger cell, sparser k-grid.

### Band structure

**Plain meaning:** A plot showing the allowed electron energies as you move through reciprocal
space. The curves ("bands") reveal whether a material is a metal, semiconductor or insulator.

**In this handbook:** A headline output of DFT, explained in
[Chapter 5 (DFT)](../ch05-dft/index.md).

**In equations:** The energy eigenvalues as a function of wavevector, $E_n(\mathbf{k})$, where
$n$ labels the band; the band gap is the energy difference between the highest filled and lowest
empty band.

**In code or input files:** Produced by solving the eigenvalue problem along a k-path and
plotting `E` versus `k`, e.g. with Matplotlib.

**Why it matters:** Conductivity, optical absorption and the band gap — the quantities that
decide whether a material is useful in electronics — are read straight off the band structure.

**Common misunderstanding:** Standard DFT systematically *underestimates* band gaps, so a
computed gap should be treated with caution rather than as an exact prediction.

### Density of states

**Plain meaning:** A count of how many electron states exist at each energy — many states packed
into a narrow energy range gives a tall peak, few states gives a low region.

**In this handbook:** Presented alongside band structure in
[Chapter 5 (DFT)](../ch05-dft/index.md).

**In equations:** $g(E) = \sum_n \int_{\mathrm{BZ}} \delta\big(E - E_n(\mathbf{k})\big)\,
\mathrm{d}\mathbf{k}$ — formally a sum of spikes, smoothed in practice.

**In code or input files:** Built as a histogram of eigenvalues, e.g.
`dos, edges = np.histogram(evals, bins=200)`, often broadened with a Gaussian.

**Why it matters:** The density of states near the highest filled level controls many
properties, from electrical conduction to chemical reactivity at surfaces.

**Common misunderstanding:** The density of states discards the k-information; two materials with
the same DOS can still have different band structures, so it is a summary, not the full picture.

---

## Molecular dynamics and statistics

### Force field

**Plain meaning:** A set of simple formulas (with fitted parameters) that give the energy and
forces of a set of atoms directly from their positions, without solving any quantum mechanics.
Also called an *interatomic potential*.

**In this handbook:** Introduced in [Chapter 7 (molecular dynamics)](../ch07-md/index.md) as the
cheap engine that drives the atoms.

**In equations:** A sum of contributions, e.g. a Lennard-Jones pair term
$V(r) = 4\varepsilon\big[(\sigma/r)^{12} - (\sigma/r)^6\big]$, plus bond, angle and other terms.

**In code or input files:** Defined by parameter blocks and a functional form, e.g. a
`pair_style lj/cut` line in a LAMMPS input (LAMMPS does not run in the browser).

**Why it matters:** Because it is cheap, a force field lets you simulate millions of atoms or
long timescales that quantum methods cannot reach.

**Common misunderstanding:** A force field is only as good as the situations it was fitted for;
used outside that range (e.g. bond breaking) it can give confidently wrong answers.

### Potential energy surface

**Plain meaning:** The landscape of energy as a function of where all the atoms are — hills,
valleys and passes. Minima are stable structures; passes between them are reaction barriers.

**In this handbook:** A unifying picture used across
[Chapter 7 (molecular dynamics)](../ch07-md/index.md) and the machine-learning chapters.

**In equations:** A single scalar function of all coordinates, $E(\mathbf{r}_1, \dots,
\mathbf{r}_N)$; the force on atom $i$ is its negative gradient,
$\mathbf{F}_i = -\nabla_{\mathbf{r}_i} E$.

**In code or input files:** Not stored explicitly — it is whatever your energy function returns;
`forces = -grad(E, positions)` samples its slope at one point.

**Why it matters:** Every structure optimisation, every MD trajectory and every reaction path is
really a journey on this surface, so the picture ties the methods together.

**Common misunderstanding:** It is enormously high-dimensional (three coordinates per atom), so
the intuitive 2D "hills and valleys" image is only a cartoon — real surfaces cannot be drawn.

### Molecular dynamics

**Plain meaning:** Simulating motion by repeatedly computing the forces on every atom and
nudging them forward a tiny time step, over and over, to watch the system evolve.

**In this handbook:** The subject of [Chapter 7 (molecular dynamics)](../ch07-md/index.md).

**In equations:** Newton's second law integrated step by step:
$\mathbf{F}_i = m_i \ddot{\mathbf{r}}_i$, advanced with a scheme like velocity Verlet over a
time step $\Delta t$.

**In code or input files:** A loop over steps, e.g.
`for step in range(nsteps): forces = compute(positions); positions += ...`; the time step
appears as `dt`.

**Why it matters:** MD links microscopic forces to measurable behaviour — diffusion, melting,
thermal expansion — by literally letting the atoms move.

**Common misunderstanding:** The time step must be small (typically around a femtosecond,
$10^{-15}$ s); too large a `dt` makes the simulation blow up because fast vibrations are missed.

### Ensemble

**Plain meaning:** The set of conditions you hold fixed during a simulation — for example
constant number of atoms, volume and temperature — together with the imagined collection of all
microstates consistent with those conditions.

**In this handbook:** A statistical-mechanics idea from
[Chapter 8 (statistical mechanics)](../ch08-statmech/index.md), applied in MD.

**In equations:** Named by their fixed quantities, e.g. NVE (constant $N$, $V$, $E$) or NVT
(constant $N$, $V$, $T$); averages are weighted by probabilities like $p \propto
e^{-E/k_{\mathrm{B}}T}$.

**In code or input files:** Chosen by keyword, e.g. an `nvt` or `npt` integrator/fix in an MD
input.

**Why it matters:** The ensemble decides which thermodynamic quantities are controlled and which
fluctuate, so it must match the experiment you are imitating.

**Common misunderstanding:** An ensemble is a collection of *many* microstates, not a single
trajectory; the link between a long time-average and the ensemble average relies on the
(usually assumed) ergodic hypothesis.

### Thermostat

**Plain meaning:** An algorithm bolted onto an MD simulation that gently adds or removes kinetic
energy so the system stays at a target temperature.

**In this handbook:** Discussed with constant-temperature MD in
[Chapter 7 (molecular dynamics)](../ch07-md/index.md) and
[Chapter 8 (statistical mechanics)](../ch08-statmech/index.md).

**In equations:** Temperature is tied to average kinetic energy via $\tfrac{3}{2}N k_{\mathrm{B}}
T = \langle \sum_i \tfrac{1}{2} m_i v_i^2 \rangle$; a thermostat steers the velocities to hit the
target $T$.

**In code or input files:** A named choice with a coupling time, e.g. `thermostat = "nose-hoover"`
and a `tau` parameter.

**Why it matters:** To simulate the NVT ensemble (constant temperature) you need a thermostat;
without one, MD conserves energy instead and the temperature drifts.

**Common misunderstanding:** Different thermostats can distort dynamics differently; some are
fine for equilibrium averages but a poor choice if you care about diffusion or transport rates.

### Free energy

**Plain meaning:** The "useful" energy that decides which state a system actually prefers at a
given temperature, balancing low energy against high disorder (entropy).

**In this handbook:** A key concept in
[Chapter 8 (statistical mechanics)](../ch08-statmech/index.md).

**In equations:** The Helmholtz free energy is $F = U - TS$ (internal energy minus temperature
times entropy); systems at fixed $T, V$ move towards minimum $F$.

**In code or input files:** Rarely a single variable — it is computed by specialised methods
(e.g. thermodynamic integration) that average over many MD snapshots.

**Why it matters:** Phase stability, solubility and reaction spontaneity are governed by free
energy, not by energy alone, because temperature and entropy matter.

**Common misunderstanding:** The lowest-energy structure is not always the stable one: at finite
temperature the system minimises *free* energy, so entropy can win and favour a higher-energy
phase.

---

## Machine learning

### Descriptor

**Plain meaning:** A fixed recipe that turns an atom and its neighbourhood into a list of numbers
a machine-learning model can read — a numerical fingerprint of "what does it look like around
here?".

**In this handbook:** Introduced in
[Chapter 9 (machine-learning interatomic potentials)](../ch09-mlip/index.md).

**In equations:** A vector-valued function of the local environment,
$\mathbf{D}_i = \mathbf{D}(\{\mathbf{r}_j - \mathbf{r}_i\})$, built to be invariant to rotation,
translation and atom relabelling.

**In code or input files:** Computed per atom into an array, e.g. `descriptors = featurise(atoms)`
of shape `(n_atoms, n_features)`.

**Why it matters:** A model is only as good as its descriptor; a well-designed descriptor with
the right symmetries is what lets the model generalise to new structures.

**Common misunderstanding:** Descriptors must respect physical symmetries — if you just feed in
raw $xyz$ coordinates, the model wrongly thinks a rotated molecule is a different one.

### Machine-learning interatomic potential

**Plain meaning:** A model trained on quantum (e.g. DFT) data to predict energies and forces
almost as accurately as DFT but far faster, acting as a learned force field.

**In this handbook:** The central topic of
[Chapter 9 (machine-learning interatomic potentials)](../ch09-mlip/index.md), often abbreviated
MLIP.

**In equations:** A parametrised energy $E_\theta(\{\mathbf{r}_i\})$ with learnable parameters
$\theta$, trained by minimising a loss over reference energies and forces,
$\mathcal{L} = \sum |E_\theta - E_{\mathrm{DFT}}|^2 + \dots$.

**In code or input files:** A trained model loaded and queried like
`energy = model.predict(atoms)` (these models use PyTorch and do not run in the browser).

**Why it matters:** MLIPs bridge the gap between DFT's accuracy and MD's speed, enabling large,
long simulations at near-quantum quality.

**Common misunderstanding:** An MLIP is reliable only near the data it was trained on;
*extrapolating* to unfamiliar structures can give large, silent errors (see active learning).

### Graph

**Plain meaning:** A way of representing a structure as dots joined by lines — atoms as dots,
their connections as lines — which is a natural fit for molecules and materials.

**In this handbook:** The data structure behind
[Chapter 10 (graph neural networks)](../ch10-gnn/index.md).

**In equations:** A graph $G = (V, E)$ is a set of nodes $V$ and edges $E$; for atoms, nodes
carry features (element, charge) and edges carry relations (distance).

**In code or input files:** Often a `Data` object with node features and an `edge_index` array
listing which atoms connect to which.

**Why it matters:** Representing a material as a graph lets a neural network respect its
connectivity and local structure directly, rather than forcing it onto a grid.

**Common misunderstanding:** A graph need not be the chemical bond network; in practice edges are
usually drawn between all atoms within a cutoff distance, bonded or not.

### Node

**Plain meaning:** A single dot in a graph. In atomistic models a node is one atom, carrying a
small list of numbers describing it.

**In this handbook:** Part of the graph picture in
[Chapter 10 (graph neural networks)](../ch10-gnn/index.md).

**In equations:** Node $i$ has a feature vector $\mathbf{h}_i$ that gets updated as the network
runs; initially it might encode the element, e.g. a one-hot of $Z_i$.

**In code or input files:** A row of the node-feature matrix, e.g. `x[i]` of shape `(n_nodes,
n_features)`.

**Why it matters:** Per-atom predictions (forces, charges) live on nodes, and a total energy is
often a sum over node contributions.

**Common misunderstanding:** A node's feature vector is not fixed input data only — it is
repeatedly *updated* by message passing as information flows in from neighbours.

### Edge

**Plain meaning:** A line joining two nodes in a graph, representing a relationship between two
atoms — typically that they are close enough to interact.

**In this handbook:** Part of the graph picture in
[Chapter 10 (graph neural networks)](../ch10-gnn/index.md).

**In equations:** An edge between atoms $i$ and $j$ often carries a feature derived from their
separation, e.g. a function of the distance $\lVert \mathbf{r}_i - \mathbf{r}_j \rVert$ or the
vector $\mathbf{r}_{ij}$.

**In code or input files:** Listed in an `edge_index` of shape `(2, n_edges)`, with optional
`edge_attr` holding distances.

**Why it matters:** Edges define who talks to whom; the cutoff distance that creates them sets
how far an atom can "feel" its surroundings in the model.

**Common misunderstanding:** Edges are usually *directed* and stored as ordered pairs, so the
pair $(i, j)$ and $(j, i)$ may both appear — do not assume a single undirected line.

### Message passing

**Plain meaning:** The core step of a graph neural network: each atom gathers information from
its neighbours, combines it, and updates its own description. Repeat a few times and information
spreads across the structure.

**In this handbook:** The mechanism explained in
[Chapter 10 (graph neural networks)](../ch10-gnn/index.md).

**In equations:** Node features update as $\mathbf{h}_i' = U\big(\mathbf{h}_i,
\sum_{j \in \mathcal{N}(i)} M(\mathbf{h}_i, \mathbf{h}_j, \mathbf{e}_{ij})\big)$, where $M$ builds
a message and $U$ updates the node.

**In code or input files:** Implemented as a layer applied several times, e.g.
`for layer in self.layers: h = layer(h, edge_index)` (PyTorch; not browser-runnable).

**Why it matters:** Stacking message-passing steps lets an atom's representation reflect an ever
wider neighbourhood, which is how the network learns about extended structure.

**Common misunderstanding:** More rounds are not always better; too many can make all nodes look
alike ("over-smoothing"), and the effective range grows only one neighbour-shell per round.

### Uncertainty

**Plain meaning:** A model's own estimate of how much to trust each prediction — a way of saying
"I am confident here" versus "this input is unfamiliar, be careful".

**In this handbook:** Important for trustworthy MLIPs and for
[Chapter 11 (active learning)](../ch11-active/index.md).

**In equations:** Often a predicted variance or standard deviation $\sigma(\mathbf{x})$, for
example the spread across an ensemble of models, $\sigma^2 = \tfrac{1}{M}\sum_m (E_m - \bar{E})^2$.

**In code or input files:** Returned alongside the prediction, e.g.
`mean, std = model.predict(atoms, return_std=True)`.

**Why it matters:** Uncertainty tells you *when not to trust* a fast model, which is essential
when an MLIP is steering an expensive simulation.

**Common misunderstanding:** A confident (low-uncertainty) prediction is not guaranteed correct;
estimates can be overconfident, especially for inputs unlike anything in the training set.

### Active learning

**Plain meaning:** A loop where the model itself flags the cases it is unsure about, you compute
the true (DFT) answer only for those, add them to the training set, and retrain — focusing
expensive effort where it helps most.

**In this handbook:** The subject of [Chapter 11 (active learning)](../ch11-active/index.md).

**In equations:** Not a single formula but a cycle: train, predict with uncertainty
$\sigma(\mathbf{x})$, select inputs where $\sigma$ is largest, label them, and repeat. (This is a
procedure, so it is best described as steps rather than one equation.)

**In code or input files:** A loop alternating training and selective labelling, e.g.
`while not converged: train(); pick = high_uncertainty(pool); label(pick)`.

**Why it matters:** It builds an accurate potential with far fewer expensive DFT calculations
than blindly sampling everything, which is often the difference between feasible and not.

**Common misunderstanding:** Active learning does not just mean "add more data". The point is to
add the *most informative* data, guided by uncertainty, not random extra points.

### Foundation model

**Plain meaning:** A single large model pre-trained on a huge, diverse dataset so that it works
reasonably across many materials out of the box, rather than being trained from scratch for one
system.

**In this handbook:** The topic of [Chapter 12 (foundation models)](../ch12-foundation/index.md).

**In equations:** Like other MLIPs it is a parametrised $E_\theta(\{\mathbf{r}_i\})$, but with
very many parameters $\theta$ fitted to a broad dataset; downstream use may *fine-tune* $\theta$
on a smaller, specific set. (The novelty is in scale and breadth, not a new equation.)

**In code or input files:** Loaded from a released checkpoint, e.g.
`model = load_pretrained("mace-mp-0")` (PyTorch-based; not browser-runnable).

**Why it matters:** A good foundation model can give useful predictions for a new system
immediately, lowering the barrier to starting a project — sometimes with little or no extra
training.

**Common misunderstanding:** "Works across many materials" does not mean "works for everything".
Foundation models still have blind spots and, as a rule of thumb, benefit from checking against
DFT and from fine-tuning for demanding cases.

---

!!! tip "Still stuck on a word?"
    If a term here is not enough, try the terse one-liner in **Appendix C**, or follow the
    chapter link in the entry's *In this handbook* field to see it used in context. The
    [formula reading guide](formula-reading-guide.md) helps with the symbols, and the
    [code reading guide](code-reading-guide.md) helps with the snippets.
