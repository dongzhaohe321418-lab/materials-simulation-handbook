# 1.2 NumPy in thirty minutes

NumPy is the foundation of scientific Python. Every other library we use in this book — SciPy, Matplotlib, ASE, pymatgen, PyTorch — either is built on NumPy or borrows its array semantics. Thirty minutes of focused practice here will save you hundreds of hours later.

## Why arrays and not lists?

A Python list is a flexible container of pointers. Each element can be a different type; each element lives somewhere different in memory. That flexibility costs you on two fronts: memory (a list of a million floats uses about 28 bytes per element, plus the pointer table) and speed (iterating over a list and doing arithmetic involves type checks at every step).

A NumPy `ndarray` is the opposite: a single contiguous block of memory holding values of one fixed type. The same million floats fit in 8 MB exactly, and looping is delegated to compiled C code.

A quick demonstration. Save and run:

```python
import numpy as np
import time

N = 10_000_000
xs = list(range(N))
arr = np.arange(N)

t0 = time.perf_counter()
s = sum(x * x for x in xs)
print(f"list:  {time.perf_counter() - t0:.3f} s")

t0 = time.perf_counter()
s = (arr * arr).sum()
print(f"array: {time.perf_counter() - t0:.3f} s")
```

On a modern laptop, the array version is 50–100× faster. The reason is not that the Python loop is "slow"; the reason is that NumPy never enters the Python interpreter for the inner loop.

The convention everywhere in this book:

```python
import numpy as np
```

## Creating arrays

The fundamental constructor is `np.array`, which converts a Python list (or list of lists) into an `ndarray`:

```python
v = np.array([1.0, 2.0, 3.0])
M = np.array([[1, 2, 3],
              [4, 5, 6]])
print(v.shape, v.dtype)   # (3,) float64
print(M.shape, M.dtype)   # (2, 3) int64
```

Several helpers create arrays without typing every value:

```python
np.zeros((3, 3))                  # 3×3 of zeros, float64
np.ones(5)                        # length-5 of ones
np.full((2, 2), 7.0)              # 2×2 filled with 7.0
np.eye(4)                         # 4×4 identity
np.arange(0, 10, 2)               # [0 2 4 6 8] — like range()
np.linspace(0, 1, 11)             # 11 evenly spaced points in [0, 1]
np.random.default_rng(0).normal(size=(3, 3))  # 3×3 of N(0,1)
```

Use `arange` when you want a step size; use `linspace` when you want a number of points. `linspace(0, 1, 11)` gives 11 points and *includes* both endpoints by default, which is the most common gotcha.

### Shape and dtype

Every array has a `shape` (a tuple of dimension sizes) and a `dtype` (the element type).

```python
a = np.zeros((2, 3, 4))           # 3-dimensional array
print(a.shape)                     # (2, 3, 4)
print(a.ndim)                      # 3
print(a.size)                      # 24 = 2*3*4
print(a.dtype)                     # float64
```

The dtype matters. `float64` is the default in NumPy; PyTorch defaults to `float32`. Mixing them in materials code is a frequent silent bug. To be explicit:

```python
np.array([1, 2, 3], dtype=np.float32)
np.zeros(10, dtype=np.int64)
```

You can change dtype with `astype`:

```python
v = np.array([1, 2, 3])
v.astype(np.float64)
```

## Indexing and slicing

NumPy arrays look like Python lists when indexed in one dimension, but extend cleanly to many dimensions. A slice `a[i:j:k]` means "from index `i` up to (but not including) `j`, step `k`":

```
1D array (length 8):

   index    0   1   2   3   4   5   6   7
   value   10  20  30  40  50  60  70  80

a[2:6]      ->  [30 40 50 60]
a[::2]      ->  [10 30 50 70]   # every second element
a[::-1]     ->  [80 70 60 50 40 30 20 10]   # reversed
```

For 2D arrays, separate indices with commas:

```python
M = np.arange(12).reshape(3, 4)
# M is
# [[ 0  1  2  3]
#  [ 4  5  6  7]
#  [ 8  9 10 11]]

M[1, 2]        # 6
M[0]           # first row -> [0 1 2 3]
M[:, 0]        # first column -> [0 4 8]
M[1:, 1:3]     # last two rows, columns 1 and 2 -> [[5 6] [9 10]]
```

!!! tip "Slices are views, not copies"
    A slice shares memory with the parent array. Modifying it modifies the original. If you want an independent copy, write `a[2:6].copy()`. This is one of the most useful — and most surprising — properties of NumPy.

```python
a = np.arange(8)
b = a[2:6]
b[0] = 999
print(a)       # [  0   1 999   3   4   5   6   7]  — 'a' was mutated
```

## Boolean masks and fancy indexing

You can index an array with a boolean array of the same shape:

```python
energies = np.array([-1.2, 0.4, -0.8, 1.6, -0.1])
mask = energies < 0
print(mask)              # [ True False  True False  True]
print(energies[mask])    # [-1.2 -0.8 -0.1]
```

This is the idiom for filtering. To count, use `mask.sum()` (booleans coerce to 0/1):

```python
print(f"{mask.sum()} negative-energy structures")
```

Fancy indexing with an integer array picks out specific rows or elements:

```python
idx = np.array([0, 2, 4])
energies[idx]            # [-1.2 -0.8 -0.1]
```

## Broadcasting — the heart of NumPy

When you write `a + b` and `a` and `b` have different shapes, NumPy applies **broadcasting rules** to align them without copying data. Master this and you stop writing `for` loops.

The rule is simple: align shapes from the **right**. Two dimensions are compatible if they are equal *or* one of them is 1.

**Worked example.** A common task: compute the distance from one reference atom to every other atom.

```python
ref   = np.array([0.0, 0.0, 0.0])           # shape (3,)
atoms = np.array([[1.0, 0.0, 0.0],          # shape (4, 3)
                  [0.0, 1.0, 0.0],
                  [0.0, 0.0, 1.0],
                  [1.0, 1.0, 1.0]])

displacements = atoms - ref                  # shape (4, 3)
distances = np.linalg.norm(displacements, axis=1)
print(distances)        # [1. 1. 1. 1.732...]
```

What did broadcasting do step by step?

1. `atoms` has shape `(4, 3)`. `ref` has shape `(3,)`.
2. Aligning from the right: the trailing `3` matches the trailing `3`.
3. `ref` has no left-hand dimension; NumPy treats it as if it had shape `(1, 3)`.
4. NumPy "stretches" `ref` along axis 0 to produce a virtual `(4, 3)` array where every row is `[0, 0, 0]`.
5. The subtraction is then elementwise.

**Pairwise distances.** Now compute the full $N \times N$ distance matrix for a small Si cell:

```python
from ase.build import bulk

si = bulk("Si", cubic=True)           # 8-atom diamond cell
positions = si.get_positions()         # shape (8, 3)

# Add a new axis to make shapes (8, 1, 3) and (1, 8, 3).
# Broadcasting then produces an (8, 8, 3) array of pairwise displacements.
diff = positions[:, None, :] - positions[None, :, :]   # shape (8, 8, 3)
D = np.linalg.norm(diff, axis=-1)                       # shape (8, 8)

print("nearest-neighbour distance:", D[D > 0].min(), "Å")
```

Three lines, no Python loops, $O(N^2)$ work done in C. (For real materials we should also worry about periodic images — see Chapter 3 — but the broadcasting pattern is the same.)

A second canonical example: outer product.

```python
x = np.linspace(0, 1, 5)               # shape (5,)
y = np.linspace(0, 2, 3)               # shape (3,)
grid = x[:, None] * y[None, :]         # shape (5, 3)
```

When broadcasting confuses you, `print(a.shape, b.shape)` solves 90% of the problem.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. State the broadcasting rule: starting from which end are shapes aligned, and when are two dimensions compatible?
    2. Why does subtracting a shape-`(3,)` reference position from a shape-`(N, 3)` array of atoms produce the correct per-atom displacements without any explicit loop?
    3. Given two arrays of shapes `(8, 1, 3)` and `(1, 8, 3)`, what is the shape of their difference, and why is this the natural way to build a pairwise distance matrix?

    If any of these is shaky, re-read the preceding section before continuing.

## Elementwise versus matrix operations

NumPy operators are **elementwise by default**. Matrix multiplication uses `@`.

```python
A = np.array([[1, 2],
              [3, 4]], dtype=float)
B = np.array([[5, 6],
              [7, 8]], dtype=float)

A * B          # elementwise:    [[5 12] [21 32]]
A @ B          # matrix product: [[19 22] [43 50]]
A.T            # transpose:      [[1 3] [2 4]]
```

Forgetting this is the second most common bug, after forgetting that slices are views.

For 1D vectors, `@` is the dot product:

```python
u = np.array([1.0, 2.0, 3.0])
v = np.array([4.0, 5.0, 6.0])
u @ v          # 32.0
```

## Reductions

Most reductions take an `axis` argument. Without one, they collapse the whole array.

```python
M = np.arange(12).reshape(3, 4)        # 3 rows, 4 columns

M.sum()              # 66, scalar
M.sum(axis=0)        # [12 15 18 21], sum each column
M.sum(axis=1)        # [ 6 22 38],    sum each row
M.mean(axis=0)       # column means
M.std(axis=1)        # row standard deviations
M.argmax(axis=1)     # index of max in each row
```

A mnemonic: `axis=k` means "the dimension that disappears". `axis=0` collapses rows into columns; `axis=1` collapses columns into rows.

If you want the dimension *retained* (for broadcasting back), pass `keepdims=True`:

```python
row_means = M.mean(axis=1, keepdims=True)   # shape (3, 1)
centred = M - row_means                      # subtracts row mean from each row
```

## Linear algebra

The functions you will use most often live in `np.linalg`:

```python
import numpy as np

A = np.array([[4.0, 1.0],
              [1.0, 3.0]])
b = np.array([1.0, 2.0])

# Solve Ax = b for x. Always prefer solve() over inv(); it is faster and more accurate.
x = np.linalg.solve(A, b)

# Eigen-decomposition of a symmetric matrix
w, V = np.linalg.eigh(A)      # eigenvalues sorted ascending; V has eigenvectors as columns

# Norms
np.linalg.norm(b)             # Euclidean (L2) norm
np.linalg.norm(A, ord="fro")  # Frobenius norm

# Determinant and inverse (use sparingly)
np.linalg.det(A)
np.linalg.inv(A)
```

For a Hamiltonian matrix, you almost always want `eigh` (Hermitian) rather than `eig` — it is faster and returns real eigenvalues. You will use this in Chapter 4 when we diagonalise a particle-in-a-box and a tight-binding Hamiltonian.

A small worked example: find the eigenfrequencies of three masses on a chain of springs.

```python
import numpy as np

k = 1.0       # spring constant
m = 1.0       # mass
n = 3

# Construct the tridiagonal stiffness matrix K
K = (2.0 * np.eye(n)
     - np.eye(n, k=1)
     - np.eye(n, k=-1)) * (k / m)

omega2, modes = np.linalg.eigh(K)
omega = np.sqrt(omega2)
print("Mode frequencies:", omega)
```

You should see three real frequencies, the lowest corresponding to the "in phase" mode and the highest to the alternating mode. This is exactly the structure of a phonon calculation in one dimension; we extend it to three dimensions in Chapter 7.

## Reading and writing data

For small arrays, the simplest text format is:

```python
np.savetxt("energies.dat", energies)              # one number per line
loaded = np.loadtxt("energies.dat")

# With a header for clarity
np.savetxt("eos.dat",
           np.column_stack([volumes, energies]),
           header="V (Å^3)   E (eV)",
           fmt="%.6f")
```

For binary, faster, and more compact:

```python
np.save("structure.npy", positions)               # one array
arr = np.load("structure.npy")

np.savez("trajectory.npz", positions=positions, forces=forces)
data = np.load("trajectory.npz")
data["positions"], data["forces"]
```

For very large datasets we will switch to HDF5 (`h5py`) in later chapters; `np.savez` is fine up to ~1 GB.

## Putting it together: a worked materials example

Compute the radial distribution function of a small silicon cell.

```python
import numpy as np
from ase.build import bulk

# Build a 2x2x2 supercell of Si (64 atoms)
si = bulk("Si", cubic=True).repeat((2, 2, 2))
positions = si.get_positions()
cell = si.cell.array          # (3, 3) lattice vectors as rows

# Pairwise displacements with minimum-image convention
diff = positions[:, None, :] - positions[None, :, :]
fractional = np.linalg.solve(cell.T, diff.reshape(-1, 3).T).T
fractional = fractional.reshape(diff.shape)
fractional -= np.round(fractional)
diff_mic = fractional @ cell
distances = np.linalg.norm(diff_mic, axis=-1)

# Take the upper triangle (excluding diagonal) to avoid double counting
i, j = np.triu_indices(len(positions), k=1)
d = distances[i, j]

# Histogram into a g(r)
bins = np.linspace(0, 6.0, 121)
counts, edges = np.histogram(d, bins=bins)
centres = 0.5 * (edges[:-1] + edges[1:])
print("First peak at r ≈", centres[counts.argmax()], "Å")
```

You should find the first peak near 2.35 Å — the Si–Si bond length. Twelve lines of NumPy, no explicit loops, fully vectorised. This is the style we will use everywhere.

## Where to next

You now know enough NumPy to read every numerical example in this book without surprise. Most remaining NumPy idioms — `einsum`, structured arrays, masked arrays — we introduce when we first need them. In the next section we learn to *look* at the numbers we have just computed, with Matplotlib.
