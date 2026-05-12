# 0.2 Linear Algebra

In Section 0.1 we worked with single numbers. The leap to materials science is the leap to many numbers at once — the three coordinates of an atom, the $3N$ degrees of freedom of an $N$-atom system, the millions of plane-wave coefficients in a DFT calculation. Linear algebra is the systematic language for handling collections of numbers that transform together. It is the most important section in this chapter; spend time on it.

## Vectors as ordered tuples

A **vector** in $\mathbb{R}^n$ is an ordered list of $n$ real numbers:

$$
\mathbf{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix} \in \mathbb{R}^n. \tag{0.2.1}
$$

We write vectors as columns by convention. The number $v_i$ is the $i$-th **component**. Two vectors are equal if and only if all their components are equal.

A vector is usually given a geometric interpretation as an arrow from the origin to the point with the same coordinates. In materials science the same algebra serves many purposes: a position, a velocity, a force, a lattice vector, a feature vector for a machine-learning model. The mathematical rules are identical.

!!! tip "What changes when $n$ grows large?"
    In three dimensions, you can draw vectors as arrows. In $3N$ dimensions — the configuration space of an $N$-atom system — drawing is impossible, but the algebra is the same. A common psychological trick is to reason about high-dimensional vectors *as if* they were three-dimensional, then trust the symbols to extend. This habit will serve you well in Chapter 4 (Hilbert space), Chapter 9 (machine-learning feature vectors), and especially Chapter 10 (graph-neural-network embeddings, which can live in $\mathbb{R}^{256}$ or higher).

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

### Algebraic and geometric faces of the dot product

The dot product wears two faces, and you must be fluent in both.

The **algebraic face** is the componentwise sum (0.2.4): a calculation any computer can do. It says nothing about geometry; given the components, you crank the handle.

The **geometric face** is (0.2.6): a statement about angles and lengths. It says nothing about coordinates; given two arrows, you measure their lengths and the angle between them.

The astonishing fact is that these two utterly different procedures yield the same number. The reconciliation goes through the law of cosines applied to the triangle with vertices at the origin, $\mathbf{u}$, and $\mathbf{v}$: the squared length of the side $\mathbf{u} - \mathbf{v}$ equals $\lVert\mathbf{u}\rVert^2 + \lVert\mathbf{v}\rVert^2 - 2\lVert\mathbf{u}\rVert\lVert\mathbf{v}\rVert\cos\theta$. Expanding the same quantity componentwise gives $\sum_i (u_i - v_i)^2 = \sum_i u_i^2 + \sum_i v_i^2 - 2 \sum_i u_i v_i$. Comparing yields (0.2.6).

When you are coding a force projection, you use the algebraic face. When you are reasoning about whether two basis functions are orthogonal, you use the geometric face. Real research demands both, often in the same line of code.

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

### Algebraic and geometric faces of the cross product

As with the dot product, the cross product has two complementary descriptions.

The **algebraic face** is formula (0.2.7), a bookkeeping device for the components. A useful mnemonic is to write the cross product as a $3 \times 3$ determinant with the standard basis vectors $\hat{\mathbf{x}}, \hat{\mathbf{y}}, \hat{\mathbf{z}}$ in the first row:

$$
\mathbf{u} \times \mathbf{v} = \det \begin{pmatrix} \hat{\mathbf{x}} & \hat{\mathbf{y}} & \hat{\mathbf{z}} \\ u_1 & u_2 & u_3 \\ v_1 & v_2 & v_3 \end{pmatrix}.
$$

Expanding along the first row reproduces (0.2.7).

The **geometric face** is the statement that $\mathbf{u} \times \mathbf{v}$ is the directed-area vector of the parallelogram with sides $\mathbf{u}$ and $\mathbf{v}$. Its magnitude is the parallelogram's area, its direction is normal to the parallelogram, and the sign is fixed by the right-hand rule (curl the fingers from $\mathbf{u}$ to $\mathbf{v}$; the thumb points along $\mathbf{u} \times \mathbf{v}$). Two consequences worth memorising: $\mathbf{u} \times \mathbf{u} = \mathbf{0}$ (a parallelogram with two parallel sides has zero area); $\mathbf{u} \times \mathbf{v} = -\mathbf{v} \times \mathbf{u}$ (reversing the orientation flips the sign).

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

### Special matrices to know by name

A few classes of matrices recur often enough that they deserve naming.

- **Identity matrix** $I$: ones on the diagonal, zeros elsewhere. Acts as a multiplicative identity: $IA = AI = A$.
- **Diagonal matrix** $D = \mathrm{diag}(d_1, \ldots, d_n)$: zeros off the diagonal. Multiplying a diagonal matrix by a vector simply scales each component.
- **Symmetric matrix**: $A^\top = A$. Examples: covariance matrices, Hessians, dynamical matrices in phonon calculations.
- **Orthogonal matrix**: $A^\top A = A A^\top = I$. Columns and rows form orthonormal bases. Geometrically: a rotation, possibly with a reflection. Determinant is $\pm 1$.
- **Unitary matrix**: complex analogue of orthogonal, $U^\dagger U = I$. Time-evolution operators in quantum mechanics are unitary.
- **Permutation matrix**: a row-permuted identity. Maps one basis vector to another.
- **Upper-triangular matrix**: zeros below the diagonal. Solving $A \mathbf{x} = \mathbf{b}$ when $A$ is triangular is trivial by back-substitution; this is why factorisations like LU and QR are valuable.

Identifying which class your matrix belongs to is the first step toward choosing an efficient algorithm. For example, `numpy.linalg.eigh` (for Hermitian matrices) is twice as fast and twice as accurate as `numpy.linalg.eig` because it exploits the spectral theorem.

### Algebraic and geometric faces of matrix multiplication

The algebraic formula (0.2.10) is "row of $A$ times column of $B$", repeated for each entry of $C$. This is what you implement in code. But there is a deeper geometric picture: matrix multiplication is **composition of linear maps**.

If $A$ is the matrix of a linear map $f$ and $B$ is the matrix of a linear map $g$, then $AB$ is the matrix of $f \circ g$ — the map that first applies $g$, then applies $f$. The non-commutativity $AB \ne BA$ is the matrix expression of the non-commutativity of function composition we met in Section 0.1. Try it physically: rotate a book by $90^\circ$ about the $x$-axis, then about the $z$-axis; now do them in the opposite order. The book ends up in different orientations.

A third, often-overlooked face of matrix multiplication is **column interpretation**: the $j$-th column of $AB$ is $A$ times the $j$-th column of $B$. This perspective is the right one for understanding how a change of basis transforms an entire family of vectors at once. We will use this in (0.2.19).

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

### The transpose, adjoint, and Hermitian property

Three closely-related operations rearrange matrix entries.

The **transpose** $A^\top$ has entries $(A^\top)_{ij} = A_{ji}$. Geometrically, transposition swaps the roles of rows and columns; algebraically, it satisfies $(AB)^\top = B^\top A^\top$ and $(A^\top)^\top = A$.

For complex matrices the more important operation is the **adjoint** (or conjugate transpose, or Hermitian conjugate) $A^\dagger$ with $(A^\dagger)_{ij} = (A_{ji})^*$. It combines transposition with complex conjugation. A matrix with $A^\dagger = A$ is **Hermitian**; in the real case this coincides with symmetric.

The adjoint encodes the most general inner-product compatibility: $\langle A \mathbf{u}, \mathbf{v} \rangle = \langle \mathbf{u}, A^\dagger \mathbf{v} \rangle$ for all vectors. This identity is the algebraic version of the integration-by-parts argument that moves a derivative from one wavefunction onto another in quantum mechanics. The fact that Hermitian operators act symmetrically in this inner-product sense is exactly what guarantees real eigenvalues — the topic of a worked proof below.

### The trace

A different scalar invariant is the **trace**, the sum of diagonal entries:
$$
\mathrm{tr}(A) = \sum_i A_{ii}.
$$
The trace has two remarkable properties: it is cyclic ($\mathrm{tr}(AB) = \mathrm{tr}(BA)$) and basis-independent ($\mathrm{tr}(V^{-1} A V) = \mathrm{tr}(A)$). Combining these with the diagonalisation (0.2.16),
$$
\mathrm{tr}(A) = \sum_i \lambda_i,
$$
the trace is the sum of eigenvalues. Likewise $\det A = \prod_i \lambda_i$.

In statistical mechanics (Chapter 8) you will meet traces of density matrices $\rho$, which encode thermal averages: $\langle O \rangle = \mathrm{tr}(\rho O)$. The cyclic property is what allows reshuffling of operator products inside such expressions.

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

### Worked example: the bonding/antibonding 2x2

Quantum mechanics in a two-state basis is the proving ground for every concept in this chapter. Consider two atomic orbitals $|a\rangle$ and $|b\rangle$ on neighbouring atoms, with on-site energies $H_{aa} = H_{bb} = -1\,\mathrm{eV}$ and hopping amplitude $H_{ab} = H_{ba} = t = -0.5\,\mathrm{eV}$. The Hamiltonian matrix in this basis is

$$
H = \begin{pmatrix} -1 & -0.5 \\ -0.5 & -1 \end{pmatrix} \,\mathrm{eV}.
$$

We solve $H\boldsymbol{\psi} = E \boldsymbol{\psi}$ step by step.

**Step (1).** Form $H - EI$:
$$
H - EI = \begin{pmatrix} -1 - E & -0.5 \\ -0.5 & -1 - E \end{pmatrix}.
$$

**Step (2).** Set $\det(H - EI) = 0$. Using the $2 \times 2$ formula (0.2.11):
$$
(-1 - E)(-1 - E) - (-0.5)(-0.5) = (1 + E)^2 - 0.25 = 0.
$$

!!! note "Why this step?"
    We multiplied the two diagonal entries and subtracted the product of the off-diagonal entries — exactly the $ad - bc$ pattern. The minus signs in both diagonal entries make their product $(1+E)^2$ positive; the two minus signs on the off-diagonals also multiply to a positive $0.25$.

**Step (3).** Solve $(1 + E)^2 = 0.25$, giving $1 + E = \pm 0.5$. Hence
$$
E_+ = -0.5\,\mathrm{eV}, \qquad E_- = -1.5\,\mathrm{eV}.
$$

The lower eigenvalue $E_-$ is the **bonding** energy; the upper $E_+$ is the **antibonding** energy. Their splitting is $|2t| = 1\,\mathrm{eV}$, as expected from the off-diagonal coupling.

**Step (4).** Find the eigenvector for $E_- = -1.5\,\mathrm{eV}$. Solve $(H - E_- I)\boldsymbol{\psi} = 0$:
$$
\begin{pmatrix} 0.5 & -0.5 \\ -0.5 & 0.5 \end{pmatrix} \begin{pmatrix} c_a \\ c_b \end{pmatrix} = 0,
$$
giving $c_a = c_b$. Normalising $c_a^2 + c_b^2 = 1$ yields
$$
\boldsymbol{\psi}_- = \frac{1}{\sqrt 2} \begin{pmatrix} 1 \\ 1 \end{pmatrix}.
$$

The bonding orbital is the **symmetric** combination — equal amplitude with the same sign on both atoms.

**Step (5).** For $E_+ = -0.5\,\mathrm{eV}$, $(H - E_+ I)\boldsymbol{\psi} = 0$ gives
$$
\begin{pmatrix} -0.5 & -0.5 \\ -0.5 & -0.5 \end{pmatrix} \begin{pmatrix} c_a \\ c_b \end{pmatrix} = 0,
$$
so $c_a = -c_b$, and
$$
\boldsymbol{\psi}_+ = \frac{1}{\sqrt 2} \begin{pmatrix} 1 \\ -1 \end{pmatrix}.
$$

The antibonding orbital is the **antisymmetric** combination, with a node between the atoms.

**Step (6).** Check $\boldsymbol{\psi}_- \cdot \boldsymbol{\psi}_+ = \tfrac{1}{2}(1\cdot 1 + 1\cdot(-1)) = 0$. Orthogonal, as the spectral theorem demands.

The squared amplitudes give probabilities: $|\boldsymbol{\psi}_-|^2 = (\tfrac{1}{2}, \tfrac{1}{2})$ — electron equally likely on either atom, with the wavefunction enhanced between them; $|\boldsymbol{\psi}_+|^2 = (\tfrac{1}{2}, \tfrac{1}{2})$ — same marginal probabilities, but with a node that reduces the electron density in the bonding region. The energy lowering of $\boldsymbol{\psi}_-$ relative to the isolated orbital energy $-1\,\mathrm{eV}$ is the *chemical bond*. We will revisit this exact two-state problem in Chapter 4 §3 when we build the Hückel and tight-binding models.

```python
import numpy as np
H = np.array([[-1.0, -0.5], [-0.5, -1.0]])
energies, states = np.linalg.eigh(H)
print(energies)   # [-1.5 -0.5]
print(states)     # columns: bonding (symm), antibonding (antisymm)
```

### Why eigenvalues are real for Hermitian matrices

In quantum mechanics, the Hamiltonian is *Hermitian*: $H^\dagger = H$, where the dagger denotes the conjugate transpose. The statement that energies are real-valued is equivalent to the statement that the eigenvalues of a Hermitian matrix are real. Let us prove this.

**Setup.** Let $H$ be a Hermitian matrix and suppose $H\boldsymbol{\psi} = \lambda \boldsymbol{\psi}$ with $\boldsymbol{\psi} \ne 0$. We want to show $\lambda \in \mathbb{R}$.

**Step (1).** Take the inner product of both sides with $\boldsymbol{\psi}$:
$$
\boldsymbol{\psi}^\dagger H \boldsymbol{\psi} = \lambda\, \boldsymbol{\psi}^\dagger \boldsymbol{\psi}.
$$
The right-hand side is $\lambda \lVert \boldsymbol{\psi}\rVert^2$, which is $\lambda$ times a positive real number.

**Step (2).** Take the complex conjugate of the entire equation. Since $(\boldsymbol{\psi}^\dagger H \boldsymbol{\psi})^* = \boldsymbol{\psi}^\dagger H^\dagger \boldsymbol{\psi}$ (taking the conjugate of a scalar is the same as taking the dagger of a $1 \times 1$ matrix), and using $H^\dagger = H$,
$$
\boldsymbol{\psi}^\dagger H \boldsymbol{\psi} = \lambda^* \, \lVert\boldsymbol{\psi}\rVert^2.
$$

!!! note "Why this step?"
    We used two facts: that scalars commute with dagger ($(c)^\dagger = c^*$), and that $(AB)^\dagger = B^\dagger A^\dagger$. Applying these to $\boldsymbol{\psi}^\dagger H \boldsymbol{\psi}$ gives $\boldsymbol{\psi}^\dagger H^\dagger \boldsymbol{\psi}$, and Hermiticity collapses $H^\dagger$ back to $H$.

**Step (3).** Compare the two expressions: $\lambda \lVert\boldsymbol{\psi}\rVert^2 = \lambda^* \lVert\boldsymbol{\psi}\rVert^2$. Since $\lVert\boldsymbol{\psi}\rVert^2 > 0$, we conclude $\lambda = \lambda^*$, i.e. $\lambda \in \mathbb{R}$. $\square$

The orthogonality of eigenvectors belonging to distinct eigenvalues follows by an analogous calculation: take the inner product of $H\boldsymbol{\psi}_i = \lambda_i \boldsymbol{\psi}_i$ with $\boldsymbol{\psi}_j$ and of $H\boldsymbol{\psi}_j = \lambda_j \boldsymbol{\psi}_j$ with $\boldsymbol{\psi}_i$, subtract, and observe that Hermiticity forces $(\lambda_i - \lambda_j) \boldsymbol{\psi}_j^\dagger \boldsymbol{\psi}_i = 0$. For distinct eigenvalues, the inner product must vanish. This is the spectral theorem at work.

### Why this matters in physics

Quantum mechanics is, at the level of computation, the eigenvalue problem

$$
\hat{H} \,\psi_n = E_n \, \psi_n,
$$

where $\hat{H}$ is the Hamiltonian operator, $E_n$ are the allowed energies, and $\psi_n$ are the stationary states. In a finite basis the operator $\hat{H}$ becomes a matrix and the equation becomes (0.2.15). Diagonalising the Hamiltonian, in Chapter 4 for the hydrogen atom and Chapter 5 for the Kohn–Sham equations of DFT, is *the* central numerical task.

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

## Gram–Schmidt orthogonalisation

Given an arbitrary linearly independent set $\{\mathbf{v}_1, \ldots, \mathbf{v}_n\}$, the **Gram–Schmidt** procedure constructs an orthonormal basis $\{\mathbf{e}_1, \ldots, \mathbf{e}_n\}$ that spans the same subspace. The recipe is geometric: take each new vector, subtract its projections onto all earlier orthonormal vectors, then normalise.

Formally, set $\mathbf{u}_1 = \mathbf{v}_1$ and $\mathbf{e}_1 = \mathbf{u}_1 / \lVert \mathbf{u}_1 \rVert$. For $k = 2, \ldots, n$,

$$
\mathbf{u}_k = \mathbf{v}_k - \sum_{j=1}^{k-1} (\mathbf{e}_j \cdot \mathbf{v}_k)\, \mathbf{e}_j, \qquad \mathbf{e}_k = \frac{\mathbf{u}_k}{\lVert \mathbf{u}_k \rVert}. \tag{0.2.22}
$$

!!! example "Gram–Schmidt on a 3D basis"
    Start with
    $$
    \mathbf{v}_1 = (1, 1, 0), \quad \mathbf{v}_2 = (1, 0, 1), \quad \mathbf{v}_3 = (0, 1, 1).
    $$

    **Step (1).** $\mathbf{u}_1 = \mathbf{v}_1$, $\lVert \mathbf{u}_1 \rVert = \sqrt 2$, so
    $$
    \mathbf{e}_1 = \frac{1}{\sqrt 2}(1, 1, 0).
    $$

    **Step (2).** Compute $\mathbf{e}_1 \cdot \mathbf{v}_2 = (1 \cdot 1 + 1 \cdot 0 + 0 \cdot 1)/\sqrt 2 = 1/\sqrt 2$. Then
    $$
    \mathbf{u}_2 = (1, 0, 1) - \frac{1}{\sqrt 2} \cdot \frac{1}{\sqrt 2}(1, 1, 0) = (1, 0, 1) - (\tfrac{1}{2}, \tfrac{1}{2}, 0) = (\tfrac{1}{2}, -\tfrac{1}{2}, 1).
    $$
    Its norm is $\sqrt{1/4 + 1/4 + 1} = \sqrt{3/2}$, so
    $$
    \mathbf{e}_2 = \sqrt{\tfrac{2}{3}}\,(\tfrac{1}{2}, -\tfrac{1}{2}, 1) = \frac{1}{\sqrt 6}(1, -1, 2).
    $$

    **Step (3).** Compute $\mathbf{e}_1 \cdot \mathbf{v}_3 = (0 + 1 + 0)/\sqrt 2 = 1/\sqrt 2$ and $\mathbf{e}_2 \cdot \mathbf{v}_3 = (0 - 1 + 2)/\sqrt 6 = 1/\sqrt 6$. Then
    $$
    \mathbf{u}_3 = (0, 1, 1) - \frac{1}{\sqrt 2}\,\mathbf{e}_1 - \frac{1}{\sqrt 6}\,\mathbf{e}_2.
    $$
    Substituting:
    $$
    \mathbf{u}_3 = (0,1,1) - (\tfrac{1}{2}, \tfrac{1}{2}, 0) - (\tfrac{1}{6}, -\tfrac{1}{6}, \tfrac{2}{6}) = (-\tfrac{2}{3}, \tfrac{2}{3}, \tfrac{2}{3}).
    $$
    Its norm is $\sqrt{4/9 + 4/9 + 4/9} = 2/\sqrt 3$. Hence
    $$
    \mathbf{e}_3 = \frac{1}{\sqrt 3}(-1, 1, 1).
    $$

    **Check.** All three pairs are orthogonal: $\mathbf{e}_1 \cdot \mathbf{e}_2 = (1 - 1 + 0)/\sqrt{12} = 0$, and similarly for the other pairs. All three have unit norm. $\square$

!!! warning "Numerical stability"
    Classical Gram–Schmidt as written is numerically unstable for large $n$ in finite precision: small floating-point errors accumulate and the output vectors lose orthogonality. The standard fix is **modified Gram–Schmidt**, which subtracts each projection immediately rather than all at once. In production code (Chapter 1 onward) use `numpy.linalg.qr`, which implements a Householder QR decomposition that is robust at any size.

The Gram–Schmidt procedure is the algebraic engine behind the **QR decomposition** of a matrix and is essential for the numerical solution of least-squares problems in Chapter 9.

## Positive-definite matrices

A symmetric (or Hermitian) matrix $A$ is **positive definite** if $\mathbf{v}^\top A \mathbf{v} > 0$ for every non-zero $\mathbf{v}$, and **positive semi-definite** if $\mathbf{v}^\top A \mathbf{v} \ge 0$. Equivalently — and this is the form that connects to Section 0.2.7 — all eigenvalues of $A$ are strictly positive (positive definite) or non-negative (positive semi-definite).

Three classes of positive-definite matrices appear constantly:

- **Hessians at local minima** are positive definite (Section 0.3): all curvatures are upward.
- **Covariance matrices** $\Sigma_{ij} = \mathrm{Cov}(X_i, X_j)$ from Section 0.5 are positive semi-definite by construction: variance of any linear combination is non-negative.
- **Gram matrices** $G_{ij} = \mathbf{v}_i \cdot \mathbf{v}_j$ of any collection of vectors are positive semi-definite.

Positive definiteness has the practical pay-off of admitting a **Cholesky decomposition** $A = L L^\top$ with $L$ lower-triangular. Solving $A \mathbf{x} = \mathbf{b}$ via Cholesky is twice as fast and numerically more stable than general LU decomposition. We will use this in Chapter 11 to invert covariance matrices in Gaussian-process regression.

## Vector and matrix norms

The Euclidean norm (0.2.5) is one of many ways to measure the "size" of a vector. Three families recur:

- $\ell^1$ **norm**: $\lVert \mathbf{v} \rVert_1 = \sum_i |v_i|$. Geometrically, the "taxicab" distance — the length of the shortest path along grid lines.
- $\ell^2$ **norm** (Euclidean, default in physics): $\lVert \mathbf{v} \rVert_2 = \sqrt{\sum_i v_i^2}$.
- $\ell^\infty$ **norm**: $\lVert \mathbf{v} \rVert_\infty = \max_i |v_i|$. The largest component in absolute value.

In machine learning (Chapter 9) you will meet **regularisation** terms that penalise the norm of model weights. $\ell^2$ regularisation (ridge regression) shrinks all weights uniformly. $\ell^1$ regularisation (lasso) drives many weights to exactly zero, producing sparse models — a useful property when interpreting which structural features matter.

Matrices also have norms. The most useful for our purposes is the **spectral norm**, the largest singular value: $\lVert A \rVert_2 = \sigma_1$. It quantifies how much $A$ stretches the worst-case input vector: $\lVert A \mathbf{v} \rVert_2 \le \lVert A \rVert_2 \cdot \lVert \mathbf{v} \rVert_2$. The **condition number** $\kappa(A) = \sigma_1 / \sigma_n$ measures how badly $A$ amplifies relative errors when solving $A \mathbf{x} = \mathbf{b}$; values $\gg 1$ flag ill-conditioned problems. A Hessian with condition number $10^{10}$, common in DFT optimisation, makes naive gradient descent crawl, motivating preconditioning techniques in Chapter 5.

## Singular value decomposition: a preview

The most general factorisation of a real matrix $A$ of any shape $m \times n$ is the **singular value decomposition** (SVD):
$$
A = U \Sigma V^\top, \tag{0.2.23}
$$
where $U$ is $m \times m$ orthogonal, $V$ is $n \times n$ orthogonal, and $\Sigma$ is $m \times n$ with non-negative entries on its main diagonal (the **singular values**) and zeros elsewhere. The singular values $\sigma_i$ are the square roots of the eigenvalues of $A^\top A$, and they measure how much $A$ stretches each principal direction.

The SVD is to general matrices what eigendecomposition is to square diagonalisable matrices: a complete structural description. It handles rank-deficient and rectangular matrices gracefully, gives the best low-rank approximation of $A$ (the **Eckart–Young theorem**), and underpins principal component analysis. Chapter 9 §4 uses the SVD to compress feature matrices; Chapter 10 uses it implicitly inside linear-algebra primitives for GNN weight updates.

## Projection operators

A square matrix $P$ is a **projection** if $P^2 = P$. Geometrically, $P$ maps a vector onto a subspace and is idempotent: projecting twice gives the same result as projecting once.

An **orthogonal projection** further satisfies $P^\top = P$ (or $P^\dagger = P$ in the complex case). It projects perpendicularly onto a subspace. If $\{\mathbf{e}_1, \ldots, \mathbf{e}_k\}$ is an orthonormal basis of the target subspace, the projector is
$$
P = \sum_{i=1}^{k} \mathbf{e}_i\, \mathbf{e}_i^\top, \tag{0.2.24}
$$
a sum of *outer products* (rank-1 matrices). Its eigenvalues are $1$ (with multiplicity $k$, eigenvectors spanning the target subspace) and $0$ (with multiplicity $n - k$, the orthogonal complement).

Projectors recur throughout quantum mechanics: $|\psi\rangle \langle \psi|$ is the projector onto the state $|\psi\rangle$, and any observable measurement is described by projecting the state onto an eigenspace of the observable. In machine learning, principal-component analysis amounts to projecting onto the top-$k$ eigenspace of a covariance matrix — exactly (0.2.24) with eigenvectors chosen by descending eigenvalue. The mathematical machinery is identical.

## A glimpse of tensors

A **tensor** is the natural generalisation of a vector and a matrix to higher rank. A scalar is a rank-0 tensor (one number); a vector is rank 1 (one index, $v_i$); a matrix is rank 2 (two indices, $A_{ij}$); a rank-3 tensor has three indices, $T_{ijk}$, and so on. The hallmark of a tensor is that its components transform predictably under a change of basis: if $V$ is the change-of-basis matrix from (0.2.19), the components of a rank-2 tensor transform as $T'_{ij} = \sum_{kl} V^{-1}_{ik} V^{-1}_{jl} T_{kl}$, with each index "rotated" individually.

Tensor algebra includes generalisations of the dot product and matrix multiplication. The two most common are **contraction** (summing over a shared index, generalising the trace and the dot product) and the **outer product** (forming a higher-rank object from lower-rank ones). The Einstein summation convention — repeated indices are summed over — is the universal shorthand.

In Chapter 9 (MACE and related equivariant neural networks) you will meet tensors of rank up to four or five, organised by their behaviour under three-dimensional rotations. The transformation rule above generalises to spherical tensors of any angular-momentum $\ell$: under a rotation $R$, the rank-$\ell$ components transform via the Wigner $D$-matrices $D^\ell(R)$. The deep reason this works is precisely that tensors are linear objects, and our entire Section 0.2 toolkit lifts to them. Get comfortable with rank-2 matrices first; the rest follows by analogy.

## Block-matrix arithmetic

A matrix can be partitioned into rectangular **blocks**, and matrix arithmetic respects this partitioning. For example, if
$$
M = \begin{pmatrix} A & B \\ C & D \end{pmatrix},
$$
with $A, B, C, D$ matrices of compatible sizes, then $M \mathbf{v}$ with $\mathbf{v} = (\mathbf{v}_1, \mathbf{v}_2)^\top$ gives $(A\mathbf{v}_1 + B\mathbf{v}_2,\, C\mathbf{v}_1 + D\mathbf{v}_2)$. Block decompositions are the basis of every divide-and-conquer linear-algebra algorithm and of the **Schur complement** identity used in Gaussian-process posteriors (Chapter 11).

A particularly useful form is the block-diagonal matrix, where off-diagonal blocks vanish: $M = \mathrm{diag}(A_1, A_2, \ldots)$. Such a matrix is invertible iff each block is, with $M^{-1} = \mathrm{diag}(A_1^{-1}, A_2^{-1}, \ldots)$. The eigenvalues of $M$ are the union of the eigenvalues of the blocks. **Symmetry-adapted bases** in quantum chemistry exploit this: by labelling basis functions by their transformation properties under the symmetry group, one block-diagonalises the Hamiltonian and reduces a large diagonalisation to several smaller ones.

## Sparse matrices, briefly

Most matrices in this book have very few non-zero entries. A real-space tight-binding Hamiltonian has hoppings only between neighbouring sites — say, $4$ to $20$ entries per row, regardless of system size $N$. Storing such a matrix densely as $N^2$ entries is wasteful; storing only the non-zero entries as a list of $(i, j, \text{value})$ triples is far cheaper, $O(N)$ instead of $O(N^2)$.

Sparse matrices come with their own algorithms. Sparse matrix-vector multiplication is $O(\text{nnz})$ in the number of non-zeros; sparse linear solves use iterative methods (conjugate gradient, GMRES) that touch the matrix only via its action on a vector. **SciPy** provides `scipy.sparse` for the data structures and `scipy.sparse.linalg` for the solvers. We will use these in Chapter 5 (large-scale Kohn–Sham diagonalisation) and Chapter 10 (graph-neural-network message passing).

## Where this is used

- Chapter 4 builds quantum mechanics on top of eigenvalue problems exactly of the form (0.2.15).
- Chapter 5 (DFT) is, computationally, repeated diagonalisation of a self-consistent Hamiltonian matrix.
- Chapter 10 (graph neural networks) layers linear maps with non-linearities; every "layer" is a matrix multiplication.

With vectors, matrices, and eigenproblems in hand, we now turn from algebra to analysis — the calculus of changing quantities.
