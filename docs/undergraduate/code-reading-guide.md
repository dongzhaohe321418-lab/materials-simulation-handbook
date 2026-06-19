# How to Read the Code

This page teaches you how to *read* the scientific Python that appears throughout this
handbook, and — just as importantly — how a piece of code maps onto a mathematical or
physical equation. It is a companion to the [formula reading guide](formula-reading-guide.md):
that page helps you read the maths, this one helps you read the code that implements it.

!!! tip "Code is not magic"
    Code is not magic. A program is simply a list of instructions that carries out a
    mathematical or physical idea, one small step at a time. Anything you can write as an
    equation, you can (in principle) write as code — and anything written as code can be
    read back, line by line, into the idea it came from. If a block looks intimidating, slow
    down and read it the way you would read a proof: one line, one symbol, one operation at a
    time.

Throughout, the only libraries we use are **NumPy**, **Matplotlib** and **SciPy**, because
these run live in the browser version of the book (see [try in Jupyter](../try-in-jupyter.md)).
Code that uses PyTorch, ASE, pymatgen, LAMMPS or Quantum ESPRESSO does *not* run in the
browser — the book always says so where it appears.

---

## Reading `import` lines

Almost every snippet starts with `import` lines. They tell Python which extra toolboxes to
load, and what short name to call them by.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
```

- `import numpy as np` loads **NumPy**, the library for numerical arrays and linear algebra,
  and gives it the alias `np`. From then on `np.array`, `np.linalg.eigh` and so on all refer
  to NumPy. The alias `np` is a near-universal convention; you will see it everywhere.
- `import matplotlib.pyplot as plt` loads the plotting library, aliased `plt`. This is what
  draws graphs.
- `from scipy import linalg` pulls one specific submodule (`linalg`, for linear algebra) out
  of **SciPy**. SciPy builds on NumPy and adds more specialised routines.

When you meet an unfamiliar function such as `np.diag`, the `np.` prefix tells you it lives in
NumPy, so you can search "numpy diag" to find its documentation. The prefix is a signpost,
not decoration.

---

## Reading NumPy arrays, and why shape matters

A NumPy **array** is an ordered grid of numbers. Most scientific code is really a sequence of
operations on arrays.

```python
x = np.linspace(0.0, 1.0, 5)   # 5 evenly spaced points from 0 to 1
print(x)          # [0.   0.25 0.5  0.75 1.  ]
print(x.shape)    # (5,)
```

The single most useful habit when reading array code is to ask: **what is the shape of this
array?** The shape is a tuple of the sizes along each dimension. You can always check it with
`.shape`. Many bugs in numerical code are *shape mismatches* — you tried to combine arrays
whose shapes were not compatible — so knowing the shape at each step tells you what an
operation is really doing.

!!! note "Read shapes out loud"
    Reading `(5,)` as "a list of five numbers" or `(3, 4)` as "three rows of four" turns an
    abstract tuple into a picture in your head. Do this for every array you meet.

---

## Scalars, vectors, matrices, batches

These four words name arrays of increasing dimensionality. The shape tells them apart.

| Object | Maths | Example shape | Meaning |
|---|---|---|---|
| Scalar | $a$ | `()` | a single number |
| Vector (1-D array) | $\mathbf{v}$ | `(n,)` | a list of $n$ numbers |
| Matrix (2-D array) | $\mathbf{M}$ | `(n, n)` | a table of numbers |
| Batch (stack) | — | `(B, n)` | $B$ vectors stacked together |

```python
a = np.float64(2.0)            # scalar, shape ()
v = np.array([1.0, 2.0, 3.0])  # vector, shape (3,)
M = np.array([[1.0, 0.0],
              [0.0, 1.0]])     # matrix, shape (2, 2)
batch = np.zeros((10, 3))      # 10 position vectors, shape (10, 3)
```

A **batch** is a stack of objects of the same kind — for example the positions of $B=10$
atoms, each a 3-component vector, stored as a `(10, 3)` array. In materials simulation the
first axis is very often "which atom" or "which configuration", and the remaining axes
describe one item. When you see a leading dimension that looks like a count, that is usually
the batch axis.

---

## Reading a function definition

A function packages a calculation so it can be reused. Read the header first; the body second.

```python
def kinetic_energy(mass: float, velocity: np.ndarray) -> float:
    """Return the kinetic energy (1/2) m v^2 for a velocity vector."""
    speed_squared = np.dot(velocity, velocity)   # v . v
    return 0.5 * mass * speed_squared
```

Read it in this order:

- `def kinetic_energy(...)` — the keyword `def` *defines* a function named `kinetic_energy`.
- `mass: float, velocity: np.ndarray` — the **arguments** (inputs). The parts after the
  colons are **type hints**: `mass` is expected to be a single number, `velocity` a NumPy
  array. Type hints are advisory notes for the reader; Python does not enforce them.
- `-> float` — the **return type**: this function gives back a single number.
- The triple-quoted line is the **docstring**: a one-line statement of what the function does.
  Read this before the body — it tells you the *intent*.
- `return ...` — the **return value**, the result handed back to whoever called the function.

So the header alone tells you: "give me a mass and a velocity vector, and I will return one
number, the kinetic energy". You can often understand a function from its header and docstring
without reading the body at all.

---

## Reading a `for` loop

A `for` loop means **do this for each item in turn**.

```python
total = 0.0
for energy in [1.0, 2.5, 0.5]:   # for each energy in the list
    total = total + energy       # add it to the running total
print(total)                     # 4.0
```

Read `for energy in [...]:` as "for each `energy` in this list". The indented lines are the
**loop body**, run once per item. Here it accumulates a sum. In NumPy you can often replace a
loop with a single *vectorised* operation (e.g. `np.sum([...])`), which is faster and shorter —
but a plain loop is usually the easiest version to read first.

---

## Reading plotting code

Plotting code almost always follows the same skeleton: make data, draw it, label it, show it.

```python
x = np.linspace(0, 2 * np.pi, 200)
y = np.sin(x)
plt.plot(x, y)                 # draw y against x
plt.xlabel("x")                # label the horizontal axis
plt.ylabel("sin(x)")           # label the vertical axis
plt.title("A sine wave")
plt.show()                     # render the figure
```

- `plt.plot(x, y)` draws a curve of `y` against `x`; both must have the same length.
- `plt.xlabel`, `plt.ylabel`, `plt.title` add text — always read these to learn what the axes
  *mean*, including their units.
- `plt.show()` displays the finished figure. (In a notebook the figure often appears without
  it.)

When you meet a plot, read the axis labels first: they tell you what quantity is being shown
and in what units, which is exactly the information a graph is for.

---

## How code maps to equations

The central skill of this page is seeing the **one-to-one correspondence** between an equation
and the code that implements it. Take the root-mean-square displacement of $N$ atoms from
their reference positions,

$$ \mathrm{RMSD} = \sqrt{\frac{1}{N}\sum_{i=1}^{N} \lVert \mathbf{r}_i - \mathbf{r}_i^{0} \rVert^2 }. $$

A faithful translation reads almost like the formula:

```python
def rmsd(r: np.ndarray, r0: np.ndarray) -> float:
    diff = r - r0                       # r_i - r_i^0   (shape (N, 3))
    sq = np.sum(diff**2, axis=1)        # ||r_i - r_i^0||^2  (shape (N,))
    return np.sqrt(np.mean(sq))         # sqrt( mean over i )
```

Match each line to a piece of the maths:

- `diff = r - r0` is the subtraction $\mathbf{r}_i - \mathbf{r}_i^{0}$, done for every atom at
  once. With `r` of shape `(N, 3)`, `diff` is also `(N, 3)`.
- `np.sum(diff**2, axis=1)` squares each component and sums *along axis 1* (the 3 coordinates),
  giving the squared distance $\lVert\cdot\rVert^2$ for each atom — a vector of shape `(N,)`.
- `np.sqrt(np.mean(sq))` takes the mean over the $N$ atoms (the $\tfrac{1}{N}\sum_i$) and then
  the square root.

Reading this way — symbol by symbol — turns "I don't understand this code" into "this line is
the subtraction, that line is the sum". When the book shows an equation next to a snippet, do
this matching explicitly; it is the fastest route to understanding both.

---

## How to debug numerical code

When a calculation gives the wrong answer (or an error), do not stare at it — *investigate* it.

- **Print shapes.** `print(arr.shape)` after each step. Most numerical bugs are shape
  mismatches; the shape tells you whether you are summing over the axis you intended.
- **Print intermediate values.** Insert `print(...)` to inspect partial results. If the final
  number is wrong, the first intermediate value that looks wrong locates the bug.
- **Check units.** A result that is off by a clean factor (say $27.2$ or $0.529$) is often a
  unit conversion — hartree vs eV, bohr vs ångström. See the units subsection below.
- **Test on a known answer.** Run the code on a case whose result you can work out by hand or
  look up analytically, and confirm it agrees. The worked example below does exactly this: it
  checks a numerical eigenvalue against an exact formula.

!!! tip "Change one parameter at a time"
    When you experiment, change **one** parameter at a time and watch what happens. Halve the
    grid spacing; does the answer get closer to the exact value? Double the box length; does
    the energy fall as you expect? If you change three things at once and the result moves, you
    will not know which change caused it. This single habit will save you more time than any
    other.

---

## Recognising units and conventions in code

Computers do arithmetic on plain numbers; they do not know whether `1.0` means one electron-
volt or one joule. **Units live in the comments, variable names and constants — not in the
syntax.** Reading code therefore includes reading the human notes around it.

```python
HARTREE_TO_EV = 27.2114     # 1 hartree in electron-volts
BOHR_TO_ANG   = 0.529177    # 1 bohr in angstrom

energy_eV = energy_ha * HARTREE_TO_EV
```

Two conventions you will meet constantly in this book:

- **Atomic units.** Many quantum snippets set $\hbar = m_e = e = 1$ so the equations stay
  uncluttered. Lengths then come out in bohr ($a_0 = 0.529177$ Å) and energies in hartree
  ($1\,\mathrm{Ha} = 27.2114$ eV). The code looks unit-free, but a comment will say which
  units are in force.
- **eV and ångström.** Materials results are often quoted in eV (energy) and Å (length). A
  variable called `energy_eV` or a comment `# in eV` is telling you the unit; the number `0.5`
  by itself never can.

When you read a snippet, hunt for these signposts before trusting any number you compute.

---

## Worked example: a particle in a box by finite differences

This is one complete example in the spirit of the book. It solves a quantum-mechanics problem
that you can also solve exactly with pen and paper, so we can *check* the code against the
known answer — the debugging habit from above, built in.

### (a) The physics

A single particle confined to a one-dimensional box of length $L$, with infinitely high walls
(so $\psi = 0$ at both ends), obeys the time-independent Schrödinger equation

$$ \hat{H}\,\psi = E\,\psi, \qquad
   \hat{H} = -\frac{\hbar^2}{2m}\frac{\mathrm{d}^2}{\mathrm{d}x^2} + V(x), $$

with $V(x) = 0$ inside the box. The exact ground-state energy is

$$ E_1 = \frac{\pi^2 \hbar^2}{2 m L^2}. $$

We will **use simple units in which $\hbar = m = 1$** (atomic-unit style), so the target
becomes $E_1 = \pi^2 / (2 L^2)$. Saying so explicitly is part of reading the code honestly: the
numbers below only make sense in these units.

### (b) The code (NumPy only, runs in the browser)

```python
import numpy as np

# Units: hbar = m = 1.
L = 1.0                       # box length
N = 200                       # number of interior grid points
x = np.linspace(0.0, L, N + 2)   # grid including both walls
h = x[1] - x[0]               # grid spacing
xi = x[1:-1]                  # interior points only (psi = 0 at the walls)
n = xi.size                   # = N

# Finite-difference Hamiltonian H = -(1/2) d^2/dx^2  (since hbar = m = 1, V = 0)
main = (2.0 / h**2) * np.ones(n)        # main diagonal:  +2/h^2
off  = (-1.0 / h**2) * np.ones(n - 1)   # off-diagonals: -1/h^2
H = 0.5 * (np.diag(main)
           + np.diag(off, 1)
           + np.diag(off, -1))          # the prefactor hbar^2/(2m) = 1/2

E, psi = np.linalg.eigh(H)    # eigenvalues E (ascending) and eigenvectors
E1_numeric  = E[0]
E1_analytic = np.pi**2 / (2.0 * L**2)

print("numeric  E1 =", E1_numeric)    # ~ 4.9347
print("analytic E1 =", E1_analytic)   # ~ 4.9348
```

Running this prints two numbers that agree to about three significant figures. That agreement
*is* the test: the code reproduces a result we already know exactly, so we can trust it on
problems we do not.

### (c) Line-by-line, mapped to the maths

- `L = 1.0`, `N = 200` — the box length $L$ and how finely we chop it up. More points means a
  more accurate answer (try changing one at a time).
- `x = np.linspace(0.0, L, N + 2)` — a grid of $N+2$ points from one wall to the other. The two
  extra points are the walls themselves, where $\psi = 0$.
- `h = x[1] - x[0]` — the grid spacing $h$. The second derivative is approximated using $h$.
- `xi = x[1:-1]` — the *interior* points, dropping both walls. Because $\psi = 0$ there, the
  walls do not appear as unknowns; this enforces the boundary condition.
- The second derivative is approximated by the standard three-point formula
  $$ \frac{\mathrm{d}^2\psi}{\mathrm{d}x^2}\Big|_{x_i} \approx
     \frac{\psi_{i+1} - 2\psi_i + \psi_{i-1}}{h^2}. $$
  As a matrix this has $-2/h^2$ on the main diagonal and $+1/h^2$ on the neighbours. The
  Hamiltonian is $-\tfrac{1}{2}\,\mathrm{d}^2/\mathrm{d}x^2$, so the overall sign flips:
- `main = (2.0 / h**2) * np.ones(n)` — the main diagonal, $+2/h^2$ (the sign from the leading
  minus in $\hat H$).
- `off = (-1.0 / h**2) * np.ones(n - 1)` — the diagonals just above and below, $-1/h^2$.
- `H = 0.5 * (np.diag(main) + np.diag(off, 1) + np.diag(off, -1))` — assembles the matrix.
  `np.diag(v)` places `v` on the main diagonal; `np.diag(v, 1)` and `np.diag(v, -1)` place it
  on the first super- and sub-diagonal. The `0.5` is the prefactor $\hbar^2/(2m) = 1/2$ in our
  units. This matrix *is* $\hat H$ on the grid.
- `E, psi = np.linalg.eigh(H)` — `eigh` solves the eigenvalue problem $\hat H \psi = E\psi$ for
  a real symmetric (Hermitian) matrix, returning eigenvalues `E` in ascending order. The
  lowest, `E[0]`, is the ground-state energy.
- The final two lines compare the numerical $E_1$ with the exact $\pi^2/(2L^2)$.

!!! note "Why the diagonals are what they are"
    Double-check the construction whenever you build such a matrix: main diagonal $+2/h^2$,
    off-diagonals $-1/h^2$, the whole thing scaled by the prefactor $\hbar^2/(2m)$ (here
    $1/2$). Get a sign or a factor wrong and the energies come out negative, or too large by a
    constant factor — both caught immediately by the comparison with $E_1 = \pi^2/(2L^2)$.

The same recipe — build a grid, write the operator as a matrix, diagonalise — reappears in the
quantum and DFT chapters ([Chapter 4 (Quantum)](../ch04-quantum/index.md),
[Chapter 5 (DFT)](../ch05-dft/index.md)) with a non-zero potential $V(x)$ added to the main
diagonal. Once you can read this example, you can read those.

---

!!! question "Check yourself"
    1. An array has `.shape` equal to `(64, 3)`. In a materials context, what is each axis most
       likely to mean, and is this a scalar, vector, matrix or batch?
    2. In the function header `def f(x: np.ndarray) -> float:`, what does each of the three
       pieces — `x`, `: np.ndarray`, `-> float` — tell you, *before* you read the body?
    3. In the worked example, why do we keep only the *interior* grid points (`x[1:-1]`) when
       building the Hamiltonian?
    4. You halve the grid spacing `h` (by doubling `N`) and the numerical $E_1$ moves *closer*
       to $\pi^2/(2L^2)$. Why is changing only this one parameter a sensible test?

??? success "Answers"
    1. It is a **batch**: the first axis (length 64) is most likely "which atom" (or which
       configuration), and the second axis (length 3) is the $x,y,z$ components of one position
       vector. So it is a stack of 64 three-component vectors.
    2. `x` is the **argument** (the input you must supply); `: np.ndarray` is a **type hint**
       saying `x` is expected to be a NumPy array; `-> float` is the **return type**, telling
       you the function gives back a single number. Together they describe the function's
       inputs and output without your having to read its body.
    3. The walls impose the boundary condition $\psi = 0$ at both ends. Dropping the wall points
       means those zero values are not treated as unknowns, so the matrix represents only the
       free interior of the box — which is exactly the problem we are solving.
    4. Because the finite-difference formula is an approximation that becomes exact as
       $h \to 0$, a smaller $h$ *should* give a more accurate eigenvalue. Changing only `h` and
       seeing the error shrink confirms the method is converging as expected; if you had also
       changed `L` at the same time you could not tell which change moved the answer.

---

For the underlying mathematics of derivatives, matrices and eigenvalues, see
[Chapter 0 (Maths)](../ch00-math/index.md) and the [formula reading guide](formula-reading-guide.md);
for a gentler tour of Python itself, see [Chapter 1 (Python)](../ch01-python/index.md). When a
term is unfamiliar, the [beginner glossary](glossary-for-beginners.md) explains it slowly.
