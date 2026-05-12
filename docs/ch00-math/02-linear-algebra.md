# 0.2 Linear Algebra

In Section 0.1 we worked with single numbers. The leap to materials science is the leap to many numbers at once — the three coordinates of an atom, the $3N$ degrees of freedom of an $N$-atom system, the millions of plane-wave coefficients in a DFT calculation. Linear algebra is the systematic language for handling collections of numbers that transform together. It is the most important section in this chapter; spend time on it.

## Vectors as ordered tuples

A **vector** in $\mathbb{R}^n$ is an ordered list of $n$ real numbers:

$$
\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix} \in \mathbb{R}^n. \tag{0.2.1}
$$

We write vectors as columns by convention. The number $v_i$ is the $i$-th **component**. Two vectors are equal if and only if all their components are equal.

A vector is usually given a geometric interpretation as an arrow from the origin to the point with the same coordinates. In materials science the same algebra serves many purposes: a position, a velocity, a force, a lattice vector, a feature vector for a machine-learning model. The mathematical rules are identical.

### Addition and scalar multiplication

Two vectors of the same length add componentwise,

$$
\mathbf{u} + \mathbf{v} = \begin{pmatrix} u_1 + v_1 \\ \vdots \\ u_n + v_n \end{pmatrix}, \tag{0.2.2}
$$

and a real scalar $c$ multiplies each component,

$$
c\,\mathbf{v} = \begin{pmatrix} c v_1 \\ \vdots \\ c v_n \end{pmatrix}. \tag{0.2.3}
$$

Geometrically, addition is the parallelogram rule and scalar multiplication stretches or flips the arrow. These two operations together make $\mathbb{R}^n$ a **vector space**: closed under addition and scalar multiplication, with a zero vector and additive inverses.

```python
import numpy as np

u: np.ndarray = np.array([1.0, 2.0, 3.0])
v: np.ndarray = np.array([4.0, -1.0, 0.5])
print(u + v)        # [5. 1. 3.5]
print(2.0 * u - v)  # [-2.  5.  5.5]
```

## The dot product

The **dot product** (also: inner product, scalar product) of two vectors in $\mathbb{R}^n$ is the scalar

$$
\mathbf{u} \cdot \mathbf{v} = \sum_{i=1}^{n} u_i v_i. \tag{0.2.4}
$$

It distributes over addition, is symmetric ($\mathbf{u}\cdot\mathbf{v} = \mathbf{v}\cdot\mathbf{u}$), and is linear in each argument. The Euclidean **norm** of a vector is

$$
\lVert \mathbf{v} \rVert = \sqrt{\mathbf{v} \cdot \mathbf{v}} = \sqrt{\sum_i v_i^2}. \tag{0.2.5}
$$

### Geometric meaning

A celebrated identity, provable from the law of cosines, gives the dot product its geometric content:

$$
\mathbf{u} \cdot \mathbf{v} = \lVert \mathbf{u} \rVert\, \lVert \mathbf{v} \rVert \cos\theta, \tag{0.2.6}
$$

where $\theta$ is the angle between the two vectors. Consequences:

- If $\mathbf{u}\cdot\mathbf{v} > 0$ the vectors point into the same half-space; if $< 0$, opposite half-spaces.
- $\mathbf{u} \cdot \mathbf{v} = 0$ exactly when the vectors are **orthogonal**.
- $\mathbf{u} \cdot \mathbf{u} = \lVert\mathbf{u}\rVert^2 \ge 0$, with equality only for the zero vector.

```python
u = np.array([1.0, 0.0, 0.0])
v = np.array([1.0, 1.0, 0.0])
cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
print(np.degrees(np.arccos(cos_theta)))  # 45.0
```

In a crystallographic context, dot products show up in computing the angle between lattice vectors, projecting a force onto a bond direction, and constructing structural descriptors for machine learning.

## The cross product

In $\mathbb{R}^3$ specifically, the **cross product** $\mathbf{u} \times \mathbf{v}$ is the vector

$$
\mathbf{u} \times \mathbf{v} =
\begin{pmatrix}
u_2 v_3 - u_3 v_2 \\
u_3 v_1 - u_1 v_3 \\
u_1 v_2 - u_2 v_1
\end{pmatrix}. \tag{0.2.7}
$$

Its magnitude is $\lVert\mathbf{u}\rVert\,\lVert\mathbf{v}\rVert\,|\sin\theta|$, the area of the parallelogram spanned by $\mathbf{u}$ and $\mathbf{v}$. Its direction is perpendicular to both, chosen by the right-hand rule.

A key application: the volume of the unit cell spanned by lattice vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$ is

$$
V = \mathbf{a}_1 \cdot (\mathbf{a}_2 \times \mathbf{a}_3). \tag{0.2.8}
$$

This **scalar triple product** is positive for a right-handed basis and equals the determinant of the matrix whose columns are the lattice vectors.

```python
a1 = np.array([1.0, 0.0, 0.0])
a2 = np.array([0.5, np.sqrt(3) / 2, 0.0])
a3 = np.array([0.0, 0.0, 1.6])
volume = np.dot(a1, np.cross(a2, a3))
print(volume)  # ≈ 1.386, a hexagonal unit cell
```

## Matrices as linear maps

A real **matrix** $A$ of size $m \times n$ is a rectangular array of numbers with $m$ rows and $n$ columns. Its entries are $A_{ij}$, with $i$ indexing the row and $j$ the column.

Matrices act on vectors. The product $A\mathbf{v}$, defined when $A$ is $m \times n$ and $\mathbf{v} \in \mathbb{R}^n$, is the vector in $\mathbb{R}^m$ with components

$$
(A\mathbf{v})_i = \sum_{j=1}^{n} A_{ij}\, v_j. \tag{0.2.9}
$$

This operation is **linear**: $A(\mathbf{u} + \mathbf{v}) = A\mathbf{u} + A\mathbf{v}$ and $A(c\,\mathbf{v}) = c\,A\mathbf{v}$. The deep theorem of linear algebra is the converse: every linear map from $\mathbb{R}^n$ to $\mathbb{R}^m$ can be represented by some $m \times n$ matrix. Linear maps and matrices are essentially the same thing.

This is why matrices show up everywhere in physics. A rotation of three-dimensional space, a Lorentz boost, a quantum-mechanical observable acting on a finite basis, a graph-convolutional layer in a neural network — all are linear maps, hence all are matrices.

### Matrix multiplication

If $A$ is $m \times n$ and $B$ is $n \times p$, the product $C = AB$ is the $m \times p$ matrix with entries

$$
C_{ij} = \sum_{k=1}^{n} A_{ik}\, B_{kj}. \tag{0.2.10}
$$

Matrix multiplication is associative ($A(BC) = (AB)C$) and distributive ($A(B+C) = AB + AC$), but **not commutative** in general: $AB \neq BA$. The order in which you apply linear maps matters — rotate-then-stretch is not stretch-then-rotate.

!!! example "A 2x2 worked product"
    Let
    $$
    A = \begin{pmatrix} 1 & 2 \\ 0 & 3 \end{pmatrix}, \qquad
    B = \begin{pmatrix} 4 & -1 \\ 1 & 2 \end{pmatrix}.
    $$
    Then
    $$
    AB = \begin{pmatrix}
    1\cdot 4 + 2 \cdot 1 & 1\cdot(-1) + 2\cdot 2 \\
    0\cdot 4 + 3 \cdot 1 & 0\cdot(-1) + 3\cdot 2
    \end{pmatrix} = \begin{pmatrix} 6 & 3 \\ 3 & 6 \end{pmatrix}.
    $$
    Meanwhile
    $$
    BA = \begin{pmatrix}
    4\cdot 1 + (-1)\cdot 0 & 4\cdot 2 + (-1)\cdot 3 \\
    1\cdot 1 + 2 \cdot 0 & 1\cdot 2 + 2 \cdot 3
    \end{pmatrix} = \begin{pmatrix} 4 & 5 \\ 1 & 8 \end{pmatrix}.
    $$
    Clearly $AB \neq BA$.

```python
A = np.array([[1.0, 2.0], [0.0, 3.0]])
B = np.array([[4.0, -1.0], [1.0, 2.0]])
print(A @ B)
# [[6. 3.]
#  [3. 6.]]
print(B @ A)
# [[4. 5.]
#  [1. 8.]]
```

The `@` operator is Python's matrix-multiplication shorthand (PEP 465). Use it; never use the elementwise `*` when you mean matrix multiplication.

## The determinant

For a square $n \times n$ matrix $A$, the **determinant** $\det A$ is a single number that captures how the linear map $A$ scales volumes. For a $2 \times 2$ matrix,

$$
\det \begin{pmatrix} a & b \\ c & d \end{pmatrix} = ad - bc, \tag{0.2.11}
$$

and for a $3 \times 3$ matrix,

$$
\det A = a_{11}(a_{22}a_{33} - a_{23}a_{32}) - a_{12}(a_{21}a_{33} - a_{23}a_{31}) + a_{13}(a_{21}a_{32} - a_{22}a_{31}). \tag{0.2.12}
$$

Key properties:

- $\det(AB) = (\det A)(\det B)$.
- $\det A = 0$ iff $A$ is singular (no inverse, columns linearly dependent).
- $|\det A|$ equals the volume of the parallelepiped spanned by $A$'s columns; the sign tells you orientation.

The last point is what links determinants to the lattice-volume formula (0.2.8): stacking the lattice vectors as columns of a $3 \times 3$ matrix and taking the determinant gives the unit-cell volume.

## The inverse

If $A$ is square and $\det A \neq 0$, there is a unique matrix $A^{-1}$ satisfying

$$
A A^{-1} = A^{-1} A = I, \tag{0.2.13}
$$

where $I$ is the identity matrix (ones on the diagonal, zeros elsewhere). For $2 \times 2$,

$$
A^{-1} = \frac{1}{\det A} \begin{pmatrix} d & -b \\ -c & a \end{pmatrix}, \quad A = \begin{pmatrix} a & b \\ c & d \end{pmatrix}. \tag{0.2.14}
$$

For larger matrices, never compute $A^{-1}$ by hand; in code, prefer solving the linear system $A\mathbf{x} = \mathbf{b}$ via `numpy.linalg.solve`, which is faster and numerically more stable than forming the inverse explicitly.

```python
A = np.array([[2.0, 1.0], [1.0, 3.0]])
b = np.array([4.0, 5.0])
x = np.linalg.solve(A, b)
print(x)                # solution
print(A @ x - b)        # ≈ [0., 0.]
```

## Eigenvalues and eigenvectors

We arrive at the central object of linear algebra. An **eigenvector** of a square matrix $A$ is a non-zero vector $\mathbf{v}$ such that

$$
A \mathbf{v} = \lambda \mathbf{v} \tag{0.2.15}
$$

for some scalar $\lambda$, the corresponding **eigenvalue**. Geometrically, $\mathbf{v}$ is a direction along which $A$ acts purely as a stretch (or flip, if $\lambda < 0$): no rotation.

### Why eigenvectors matter

Most matrices encountered in physics are diagonalisable: there exists a basis of eigenvectors. In that basis the matrix becomes diagonal,

$$
A = V D V^{-1}, \qquad D = \mathrm{diag}(\lambda_1, \ldots, \lambda_n), \tag{0.2.16}
$$

with $V$ the matrix whose columns are the eigenvectors. Diagonal matrices are trivial to work with: they commute, their powers are componentwise powers, their exponentials are componentwise exponentials. The whole point of "diagonalising a Hamiltonian" in quantum mechanics is to find this basis.

### Worked example: full derivation for a 2x2

Let us solve the eigenvalue problem for

$$
A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}.
$$

**Step 1.** Rearrange (0.2.15) as $(A - \lambda I)\mathbf{v} = 0$. For a non-zero $\mathbf{v}$, the matrix $A - \lambda I$ must be singular, so

$$
\det(A - \lambda I) = 0. \tag{0.2.17}
$$

This is the **characteristic equation**.

**Step 2.** Compute the characteristic polynomial:

$$
\det \begin{pmatrix} 2 - \lambda & 1 \\ 1 & 2 - \lambda \end{pmatrix}
= (2 - \lambda)^2 - 1 = \lambda^2 - 4\lambda + 3.
$$

**Step 3.** Solve $\lambda^2 - 4\lambda + 3 = 0$:

$$
\lambda = \frac{4 \pm \sqrt{16 - 12}}{2} = \frac{4 \pm 2}{2} = 3 \text{ or } 1.
$$

**Step 4.** For each eigenvalue, find the eigenvector.

For $\lambda_1 = 3$, solve $(A - 3I)\mathbf{v} = 0$:

$$
\begin{pmatrix} -1 & 1 \\ 1 & -1 \end{pmatrix} \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = 0,
$$

giving $v_1 = v_2$. Choose $\mathbf{v}_1 = (1, 1)/\sqrt{2}$ (unit norm).

For $\lambda_2 = 1$, solve $(A - I)\mathbf{v} = 0$:

$$
\begin{pmatrix} 1 & 1 \\ 1 & 1 \end{pmatrix} \begin{pmatrix} v_1 \\ v_2 \end{pmatrix} = 0,
$$

giving $v_1 = -v_2$. Choose $\mathbf{v}_2 = (1, -1)/\sqrt{2}$.

Notice that $\mathbf{v}_1 \cdot \mathbf{v}_2 = 0$: the eigenvectors are orthogonal. This is no accident — $A$ is symmetric ($A^\top = A$), and symmetric real matrices always have real eigenvalues and orthogonal eigenvectors. This is the **spectral theorem**, and it is the mathematical foundation of why Hamiltonian eigenvalues (energies) are real.

```python
A = np.array([[2.0, 1.0], [1.0, 2.0]])
eigvals, eigvecs = np.linalg.eigh(A)  # 'eigh' for Hermitian matrices
print(eigvals)   # [1. 3.]
print(eigvecs)
# columns are eigenvectors; sign convention is arbitrary
```

### Why this matters in physics

Quantum mechanics is, at the level of computation, the eigenvalue problem

$$
\hat H \,\psi_n = E_n \, \psi_n,
$$

where $\hat H$ is the Hamiltonian operator, $E_n$ are the allowed energies, and $\psi_n$ are the stationary states. In a finite basis the operator $\hat H$ becomes a matrix and the equation becomes (0.2.15). Diagonalising the Hamiltonian, in Chapter 4 for the hydrogen atom and Chapter 5 for the Kohn–Sham equations of DFT, is *the* central numerical task.

Eigenvalue problems are equally pervasive in classical materials science. The vibrational frequencies of a crystal are eigenvalues of the dynamical matrix (Chapter 8). The principal axes of a stress or strain tensor are its eigenvectors. The leading principal components of a feature matrix in machine learning are again eigenvectors. Get comfortable with (0.2.15); it will pay back the effort hundreds of times.

## Basis sets and change of basis

A **basis** of $\mathbb{R}^n$ is a set of $n$ linearly independent vectors $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$. Every vector $\mathbf{v}$ admits a unique expansion

$$
\mathbf{v} = \sum_{i=1}^{n} c_i\, \mathbf{e}_i, \tag{0.2.18}
$$

with the $c_i$ called the **coefficients** of $\mathbf{v}$ in this basis. The standard basis is $\mathbf{e}_i = (0, \ldots, 1, \ldots, 0)$ with the $1$ in the $i$-th slot, and the components of a vector are its coefficients in this basis.

Changing basis is a change of coordinates. If $V$ is the matrix whose columns are a new basis (expressed in the old basis), then the coefficients $\mathbf{c}$ of a vector in the new basis relate to its old components $\mathbf{v}$ by

$$
\mathbf{v} = V \mathbf{c}, \qquad \mathbf{c} = V^{-1} \mathbf{v}. \tag{0.2.19}
$$

This is the formal content of "expressing a quantum state in the energy eigenbasis" or "rotating into the principal-axis frame".

## Orthonormality

A basis is **orthonormal** if its vectors are unit-norm and mutually orthogonal:

$$
\mathbf{e}_i \cdot \mathbf{e}_j = \delta_{ij}, \tag{0.2.20}
$$

where $\delta_{ij}$ is the Kronecker delta ($1$ if $i = j$, $0$ otherwise). In an orthonormal basis the coefficients can be read off by projection,

$$
c_i = \mathbf{e}_i \cdot \mathbf{v}, \tag{0.2.21}
$$

without ever computing an inverse matrix. Orthonormal bases are computationally privileged for exactly this reason. In quantum mechanics one almost always works in an orthonormal basis of states; the matrix $V$ in (0.2.19) becomes **orthogonal** ($V^\top V = I$), with the cleanest possible inverse $V^{-1} = V^\top$.

```python
A = np.array([[2.0, 1.0], [1.0, 2.0]])
eigvals, V = np.linalg.eigh(A)
print(V.T @ V)      # ≈ I, confirming orthonormal columns
print(V.T @ A @ V)  # ≈ diag(eigvals), confirming diagonalisation
```

## Where this is used

- Chapter 4 builds quantum mechanics on top of eigenvalue problems exactly of the form (0.2.15).
- Chapter 5 (DFT) is, computationally, repeated diagonalisation of a self-consistent Hamiltonian matrix.
- Chapter 10 (graph neural networks) layers linear maps with non-linearities; every "layer" is a matrix multiplication.

With vectors, matrices, and eigenproblems in hand, we now turn from algebra to analysis — the calculus of changing quantities.
