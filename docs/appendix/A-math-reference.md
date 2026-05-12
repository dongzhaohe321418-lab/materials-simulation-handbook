# Appendix A — Mathematical Reference

This appendix gathers the formulae, identities and unit conversions
that appear throughout the book. It is intended as a lookup resource;
none of the material is derived in detail here. Where a result is
non-trivial, a chapter cross-reference points to the derivation.

---

## A.1 Physical constants

In what follows, **SI** denotes the International System of Units and
**a.u.** denotes Hartree atomic units, in which
$\hbar = m_e = e = 1/(4\pi\epsilon_0) = 1$.

| Constant | Symbol | SI value | Atomic units |
|---|---|---|---|
| Speed of light in vacuum | $c$ | $2.998 \times 10^8$ m·s$^{-1}$ | $137.036\ (\alpha^{-1})$ |
| Planck constant | $h$ | $6.626 \times 10^{-34}$ J·s | $2\pi$ |
| Reduced Planck constant | $\hbar$ | $1.055 \times 10^{-34}$ J·s | $1$ |
| Elementary charge | $e$ | $1.602 \times 10^{-19}$ C | $1$ |
| Electron mass | $m_e$ | $9.109 \times 10^{-31}$ kg | $1$ |
| Proton mass | $m_p$ | $1.673 \times 10^{-27}$ kg | $1836.15$ |
| Boltzmann constant | $k_B$ | $1.381 \times 10^{-23}$ J·K$^{-1}$ | $3.167 \times 10^{-6}$ Ha·K$^{-1}$ |
| Vacuum permittivity | $\epsilon_0$ | $8.854 \times 10^{-12}$ F·m$^{-1}$ | $1/(4\pi)$ |
| Avogadro number | $N_A$ | $6.022 \times 10^{23}$ mol$^{-1}$ | — |
| Fine-structure constant | $\alpha$ | $7.297 \times 10^{-3}$ | $1/c$ |
| Bohr radius | $a_0$ | $5.292 \times 10^{-11}$ m | $1$ |
| Hartree energy | $E_h$ | $4.360 \times 10^{-18}$ J | $1$ |

---

## A.2 Unit conversions

### Energy

| From | To | Multiply by |
|---|---|---|
| Ha (Hartree) | eV | $27.2114$ |
| Ha | kcal·mol$^{-1}$ | $627.509$ |
| Ha | kJ·mol$^{-1}$ | $2625.50$ |
| Ha | J | $4.3597 \times 10^{-18}$ |
| Ha | K (via $k_B$) | $3.158 \times 10^{5}$ |
| Ha | cm$^{-1}$ | $2.1947 \times 10^{5}$ |
| eV | J | $1.602 \times 10^{-19}$ |
| eV | K | $1.160 \times 10^{4}$ |
| eV | cm$^{-1}$ | $8065.5$ |
| kcal·mol$^{-1}$ | meV | $43.36$ |
| kJ·mol$^{-1}$ | meV | $10.36$ |

The conversion energy $\leftrightarrow$ temperature is
$$
E = k_B T,
$$
so $1$ Ha corresponds to a temperature of $E_h / k_B
\approx 3.158 \times 10^5$ K, and *room temperature* ($T = 300$ K)
corresponds to $k_B T \approx 25.85$ meV $\approx 0.95 \times 10^{-3}$
Ha.

### Length

| From | To | Multiply by |
|---|---|---|
| bohr ($a_0$) | Å | $0.52918$ |
| bohr | nm | $0.052918$ |
| bohr | m | $5.2918 \times 10^{-11}$ |
| Å | bohr | $1.8897$ |
| Å | nm | $0.1$ |

### Time

| From | To | Multiply by |
|---|---|---|
| atomic time unit | s | $2.4189 \times 10^{-17}$ |
| atomic time unit | fs | $0.02419$ |
| ps | fs | $10^3$ |
| ns | ps | $10^3$ |

The atomic unit of time is $\hbar / E_h$, the time taken by an
electron in the ground state of hydrogen (Bohr orbit) to traverse
$1$ radian of phase.

### Force, pressure, dipole moment

| From | To | Multiply by |
|---|---|---|
| Ha / bohr | eV / Å | $51.422$ |
| Ha / bohr$^3$ | GPa | $2.942 \times 10^4$ |
| eV / Å$^3$ | GPa | $160.218$ |
| a.u. (dipole) | Debye | $2.5418$ |

---

## A.3 Linear algebra cheatsheet

Throughout this section, $A, B, C$ are square matrices of compatible
dimension; $\mathbf{x}, \mathbf{y}$ are column vectors.

### Matrix identities

$$
(AB)^\top = B^\top A^\top, \qquad
(AB)^{-1} = B^{-1} A^{-1},
$$
$$
\det(AB) = \det(A)\det(B), \qquad
\det(A^\top) = \det(A),
$$
$$
\det(A^{-1}) = \frac{1}{\det(A)}, \qquad
\det(cA) = c^n \det(A)\ \text{for}\ A \in \mathbb{R}^{n\times n},
$$
$$
\operatorname{tr}(AB) = \operatorname{tr}(BA), \qquad
\operatorname{tr}(A+B) = \operatorname{tr}(A) + \operatorname{tr}(B).
$$

### Eigenvalues

If $A\mathbf{v} = \lambda\mathbf{v}$:

$$
\det(A) = \prod_i \lambda_i, \qquad
\operatorname{tr}(A) = \sum_i \lambda_i.
$$

For a Hermitian (or symmetric real) matrix, all eigenvalues are
real and eigenvectors for distinct eigenvalues are orthogonal. Any
Hermitian $A$ admits a spectral decomposition
$$
A = U \Lambda U^\dagger,
$$
with $U$ unitary and $\Lambda$ diagonal.

The *Rayleigh quotient* $R_A(\mathbf{x}) = \mathbf{x}^\dagger A
\mathbf{x} / \mathbf{x}^\dagger \mathbf{x}$ is bounded between
$\lambda_\text{min}$ and $\lambda_\text{max}$, with equality on the
corresponding eigenvectors. This underlies the variational principle
of Chapter 4.

### Matrix derivatives

For $\mathbf{x} \in \mathbb{R}^n$ and constant $A$:

$$
\frac{\partial}{\partial \mathbf{x}}(\mathbf{a}^\top \mathbf{x})
= \mathbf{a}, \qquad
\frac{\partial}{\partial \mathbf{x}}(\mathbf{x}^\top A \mathbf{x})
= (A + A^\top)\mathbf{x}.
$$

For a scalar function $f$ of a matrix $X$:

$$
\frac{\partial}{\partial X}\operatorname{tr}(AX) = A^\top, \qquad
\frac{\partial}{\partial X}\operatorname{tr}(AX^\top) = A,
$$
$$
\frac{\partial}{\partial X}\log\det X = (X^{-1})^\top.
$$

### Block matrix inversion

$$
\begin{pmatrix} A & B \\ C & D \end{pmatrix}^{-1}
=
\begin{pmatrix}
A^{-1} + A^{-1} B S^{-1} C A^{-1} & -A^{-1} B S^{-1} \\
-S^{-1} C A^{-1} & S^{-1}
\end{pmatrix},
$$
where $S = D - C A^{-1} B$ is the Schur complement. Useful for
deriving the equations of constrained optimisation and for
manipulating Kohn–Sham systems with non-orthogonal bases.

### The Woodbury identity

$$
(A + UCV)^{-1} = A^{-1} - A^{-1} U (C^{-1} + V A^{-1} U)^{-1} V A^{-1}.
$$

Used whenever a low-rank update to a covariance or Fock matrix must
be inverted cheaply.

---

## A.4 Vector calculus

### Standard differential operators

In Cartesian coordinates, for a scalar field $f$ and a vector field
$\mathbf{A}$:

$$
\nabla f = \left(\partial_x f,\ \partial_y f,\ \partial_z f\right),
$$
$$
\nabla \cdot \mathbf{A} = \partial_x A_x + \partial_y A_y + \partial_z A_z,
$$
$$
\nabla \times \mathbf{A} =
\bigl(
\partial_y A_z - \partial_z A_y,\
\partial_z A_x - \partial_x A_z,\
\partial_x A_y - \partial_y A_x
\bigr),
$$
$$
\nabla^2 f = \partial_x^2 f + \partial_y^2 f + \partial_z^2 f.
$$

### Identities

For any smooth $f$ and $\mathbf{A}$:

$$
\nabla \times (\nabla f) = \mathbf{0}, \qquad
\nabla \cdot (\nabla \times \mathbf{A}) = 0.
$$

For vector fields $\mathbf{A}, \mathbf{B}$:

$$
\nabla \times (\nabla \times \mathbf{A})
= \nabla(\nabla \cdot \mathbf{A}) - \nabla^2 \mathbf{A},
$$
$$
\nabla \cdot (\mathbf{A} \times \mathbf{B})
= \mathbf{B}\cdot(\nabla\times\mathbf{A}) - \mathbf{A}\cdot(\nabla\times\mathbf{B}),
$$
$$
\nabla(\mathbf{A}\cdot\mathbf{B})
= (\mathbf{A}\cdot\nabla)\mathbf{B} + (\mathbf{B}\cdot\nabla)\mathbf{A}
+ \mathbf{A}\times(\nabla\times\mathbf{B})
+ \mathbf{B}\times(\nabla\times\mathbf{A}).
$$

### Integral theorems

**Divergence (Gauss) theorem.** For a vector field $\mathbf{F}$
defined on a region $V$ with boundary $\partial V$:
$$
\int_V \nabla \cdot \mathbf{F}\ \mathrm{d}V
= \oint_{\partial V} \mathbf{F} \cdot \mathbf{n}\ \mathrm{d}S.
$$

**Stokes' theorem.** For a vector field $\mathbf{F}$ on a surface
$S$ with boundary $\partial S$:
$$
\int_S (\nabla \times \mathbf{F}) \cdot \mathbf{n}\ \mathrm{d}S
= \oint_{\partial S} \mathbf{F} \cdot \mathrm{d}\boldsymbol{\ell}.
$$

**Green's theorem** (2D special case):
$$
\oint_{\partial D} (P\,\mathrm{d}x + Q\,\mathrm{d}y)
= \iint_D \left(\frac{\partial Q}{\partial x} - \frac{\partial P}{\partial y}\right) \mathrm{d}A.
$$

**Green's identities** (for scalar fields $\phi, \psi$):
$$
\int_V (\phi \nabla^2 \psi + \nabla\phi \cdot \nabla\psi)\ \mathrm{d}V
= \oint_{\partial V} \phi\,\nabla\psi \cdot \mathbf{n}\ \mathrm{d}S,
$$
$$
\int_V (\phi \nabla^2 \psi - \psi \nabla^2 \phi)\ \mathrm{d}V
= \oint_{\partial V} (\phi\nabla\psi - \psi\nabla\phi)\cdot\mathbf{n}\ \mathrm{d}S.
$$

The second of these underlies the derivation of the Hellmann–Feynman
theorem and the divergence-form Poisson solver of Chapter 6.

---

## A.5 Common integrals

### Gaussian integrals

$$
\int_{-\infty}^{\infty} e^{-\alpha x^2}\ \mathrm{d}x = \sqrt{\frac{\pi}{\alpha}}, \qquad \alpha > 0,
$$
$$
\int_{-\infty}^{\infty} x^2 e^{-\alpha x^2}\ \mathrm{d}x
= \frac{1}{2\alpha}\sqrt{\frac{\pi}{\alpha}},
$$
$$
\int_{-\infty}^{\infty} e^{-\alpha x^2 + \beta x}\ \mathrm{d}x
= \sqrt{\frac{\pi}{\alpha}}\,\exp\!\left(\frac{\beta^2}{4\alpha}\right).
$$

In $d$ dimensions, for a positive-definite matrix $A$:
$$
\int_{\mathbb{R}^d} e^{-\tfrac12 \mathbf{x}^\top A \mathbf{x} + \mathbf{b}^\top \mathbf{x}}\ \mathrm{d}^d x
= \frac{(2\pi)^{d/2}}{\sqrt{\det A}}\,
\exp\!\left(\tfrac12 \mathbf{b}^\top A^{-1} \mathbf{b}\right).
$$

### Exponential / gamma integrals

$$
\int_0^\infty x^{n-1} e^{-x}\,\mathrm{d}x = \Gamma(n),
\quad \Gamma(n+1) = n\,\Gamma(n),\
\Gamma(1) = 1,\ \Gamma(1/2) = \sqrt{\pi}.
$$
$$
\int_0^\infty x^n e^{-\alpha x}\,\mathrm{d}x = \frac{n!}{\alpha^{n+1}},
\qquad
\int_0^\infty \frac{x^{s-1}}{e^x - 1}\,\mathrm{d}x = \Gamma(s)\zeta(s).
$$

The last appears in the derivation of the Debye model in
Chapter 8.

### Useful trigonometric / spherical integrals

$$
\int_0^\pi \sin\theta\,\mathrm{d}\theta = 2, \qquad
\int_0^{2\pi}\mathrm{d}\phi = 2\pi,
$$
so the solid-angle integral over the sphere is $4\pi$. The
spherical harmonics satisfy
$$
\int_0^{2\pi}\!\!\int_0^\pi Y_\ell^{m\,*}(\theta,\phi) Y_{\ell'}^{m'}(\theta,\phi)
\sin\theta\,\mathrm{d}\theta\,\mathrm{d}\phi
= \delta_{\ell\ell'}\delta_{mm'}.
$$

### The Coulomb integral

$$
\int_{\mathbb{R}^3} \frac{e^{-2\zeta r}}{|\mathbf{r} - \mathbf{r}'|}
\,\mathrm{d}^3 r
= \frac{4\pi}{(2\zeta)^2} \cdot \frac{1}{r'}\left[1 - (1 + \zeta r')e^{-2\zeta r'}\right].
$$

The cornerstone evaluation underlying every Slater-type integral in
quantum chemistry. Used in Chapter 4.

---

## A.6 Taylor series

For a sufficiently smooth $f$ about $x_0$:

$$
f(x) = f(x_0) + f'(x_0)(x-x_0) + \tfrac12 f''(x_0)(x-x_0)^2 + \ldots
$$

Common one-variable series ($|x|$ small):

$$
e^x = 1 + x + \tfrac{x^2}{2!} + \tfrac{x^3}{3!} + \ldots,
$$
$$
\sin x = x - \tfrac{x^3}{3!} + \tfrac{x^5}{5!} - \ldots,
$$
$$
\cos x = 1 - \tfrac{x^2}{2!} + \tfrac{x^4}{4!} - \ldots,
$$
$$
\ln(1+x) = x - \tfrac{x^2}{2} + \tfrac{x^3}{3} - \ldots,
$$
$$
(1+x)^\alpha = 1 + \alpha x + \tfrac{\alpha(\alpha-1)}{2}x^2 + \ldots
$$

Multivariate Taylor expansion to second order about $\mathbf{x}_0$:
$$
f(\mathbf{x}) \approx f(\mathbf{x}_0) + \nabla f(\mathbf{x}_0)^\top
\delta\mathbf{x} + \tfrac12 \delta\mathbf{x}^\top H \delta\mathbf{x},
$$
where $H_{ij} = \partial_i\partial_j f(\mathbf{x}_0)$ is the
Hessian. The harmonic approximation to vibrations (Chapter 8) and
the second-order optimisation schemes (Chapter 6) are direct
applications.

---

## A.7 Probability and statistics

For a random variable $X$ with density $p(x)$:

$$
\mathbb{E}[X] = \int x\,p(x)\,\mathrm{d}x, \qquad
\mathrm{Var}(X) = \mathbb{E}[X^2] - \mathbb{E}[X]^2.
$$

For independent $X, Y$:
$$
\mathrm{Var}(X+Y) = \mathrm{Var}(X) + \mathrm{Var}(Y).
$$

**Gaussian density** in one dimension:
$$
p(x) = \frac{1}{\sqrt{2\pi\sigma^2}} \exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right).
$$

Multivariate Gaussian:
$$
p(\mathbf{x}) = \frac{1}{\sqrt{(2\pi)^d \det \Sigma}}
\exp\!\left(-\tfrac12 (\mathbf{x}-\boldsymbol{\mu})^\top \Sigma^{-1} (\mathbf{x}-\boldsymbol{\mu})\right).
$$

The **Boltzmann distribution** at temperature $T$ for a system with
energies $\{E_i\}$:
$$
p_i = \frac{e^{-\beta E_i}}{Z},
\qquad
\beta = 1/(k_B T),
\qquad
Z = \sum_i e^{-\beta E_i}.
$$

The connection to thermodynamic free energy is $F = -k_B T \ln Z$.

---

## A.8 Useful inequalities

**Cauchy–Schwarz.**
$$
|\langle \mathbf{x}, \mathbf{y}\rangle| \le \|\mathbf{x}\|\,\|\mathbf{y}\|.
$$

**Triangle.**
$$
\|\mathbf{x}+\mathbf{y}\| \le \|\mathbf{x}\| + \|\mathbf{y}\|.
$$

**Jensen.** For convex $\varphi$ and any random variable $X$,
$$
\varphi(\mathbb{E}[X]) \le \mathbb{E}[\varphi(X)].
$$

**Variational (Rayleigh–Ritz).** For any normalised trial state
$|\psi\rangle$ and Hermitian $H$,
$$
\langle\psi|H|\psi\rangle \ge E_0.
$$

This is the principle that underwrites *every* variational
method in this book — Hartree–Fock, DFT, configuration interaction,
quantum Monte Carlo.

---

The remainder of the appendix turns to operational matters: how to
run these methods on real hardware (Appendix B), what the terminology
means (Appendix C), where to read further (Appendix D), and a
consolidated bibliography (Appendix E).
