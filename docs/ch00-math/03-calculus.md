# 0.3 Calculus and Gradients

Section 0.2 gave us the algebra of constant quantities — vectors and matrices that simply *are*. Materials science, however, is dominated by quantities that *change*: energies as a function of atomic positions, densities as a function of position in space, loss functions as a function of model parameters. Calculus is the systematic language of change, and the gradient — its multivariable workhorse — is arguably the single most-used tool in computational physics and machine learning.

## The derivative

Given a real-valued function $f$ of a real variable $x$, the **derivative** at $x_0$ is the limit

$$
f'(x_0) \;=\; \lim_{h \to 0} \frac{f(x_0 + h) - f(x_0)}{h}, \tag{0.3.1}
$$

provided the limit exists. Geometrically, $f'(x_0)$ is the slope of the tangent line to the graph of $f$ at $x_0$. Physically, if $x$ is time and $f$ is position, $f'$ is velocity; if $f$ is energy and $x$ is a coordinate, $-f'$ is a force component.

We use several notations interchangeably:

$$
f'(x) \;=\; \frac{\mathrm{d} f}{\mathrm{d} x} \;=\; \frac{\mathrm{d}}{\mathrm{d} x} f(x).
$$

A function is **differentiable** at $x_0$ if (0.3.1) exists; it is differentiable on an interval if it is differentiable at every point of the interval. Differentiability implies continuity but not vice versa: $|x|$ is continuous at $0$ but not differentiable there.

### Rules of differentiation

The following rules are established by elementary manipulation of (0.3.1) and are worth memorising. Let $f$ and $g$ be differentiable.

**Linearity.** For constants $a, b$,
$$
(af + bg)'(x) = a f'(x) + b g'(x). \tag{0.3.2}
$$

**Product rule.**
$$
(fg)'(x) = f'(x)\, g(x) + f(x)\, g'(x). \tag{0.3.3}
$$

**Quotient rule.** Where $g(x) \neq 0$,
$$
\left( \frac{f}{g} \right)'(x) = \frac{f'(x) g(x) - f(x) g'(x)}{g(x)^2}. \tag{0.3.4}
$$

**Chain rule.** For the composition $h(x) = f(g(x))$,
$$
h'(x) = f'(g(x)) \cdot g'(x). \tag{0.3.5}
$$

The chain rule is the most important of the four. It is the engine of backpropagation in neural networks, of force evaluation through interatomic potentials, and of every reparameterisation trick in modern materials ML.

### A small catalogue

The following standard derivatives recur constantly:

| $f(x)$ | $f'(x)$ |
|---|---|
| $x^n$ | $n x^{n-1}$ |
| $e^{ax}$ | $a e^{ax}$ |
| $\ln x$ | $1/x$ |
| $\sin x$ | $\cos x$ |
| $\cos x$ | $-\sin x$ |

!!! example "Chain rule on the Boltzmann factor"
    Differentiate $p(T) = \exp\!\big(-E_\mathrm{a} / (k_\mathrm{B} T)\big)$ with respect to $T$.

    Write $p = e^{u}$ with $u = -E_\mathrm{a}/(k_\mathrm{B} T)$. Then $\mathrm{d}u/\mathrm{d}T = E_\mathrm{a}/(k_\mathrm{B} T^2)$, and the chain rule gives
    $$
    \frac{\mathrm{d} p}{\mathrm{d} T} = e^{u} \cdot \frac{\mathrm{d} u}{\mathrm{d} T} = \exp\!\left(-\frac{E_\mathrm{a}}{k_\mathrm{B} T}\right) \cdot \frac{E_\mathrm{a}}{k_\mathrm{B} T^2}.
    $$
    Positive, as expected: hotter systems explore higher-energy states more often.

## Taylor expansion

A smooth function can be approximated near a point $x_0$ by a polynomial in $(x - x_0)$. The $n$-th order **Taylor expansion** of $f$ about $x_0$ is

$$
f(x) \;\approx\; \sum_{k=0}^{n} \frac{f^{(k)}(x_0)}{k!} (x - x_0)^k, \tag{0.3.6}
$$

with $f^{(k)}$ the $k$-th derivative. The remainder shrinks as $|x - x_0|^{n+1}$ for sufficiently smooth $f$.

For $n = 1$ this is the **linearisation** $f(x) \approx f(x_0) + f'(x_0)(x - x_0)$ — the tangent line. For $n = 2$ we add a curvature term:

$$
f(x) \approx f(x_0) + f'(x_0)(x - x_0) + \tfrac{1}{2} f''(x_0) (x - x_0)^2. \tag{0.3.7}
$$

This second-order picture justifies the **harmonic approximation** in materials science: near an equilibrium structure, the energy as a function of displacement is locally quadratic, with curvature given by $f''$. The eigenvalues of the matrix of second derivatives (the *Hessian*) are the squared phonon frequencies.

### Worked example: Taylor expansion of $\cos$

Compute the second-order Taylor expansion of $\cos x$ about $x_0 = 0$.

Derivatives at zero: $\cos 0 = 1$, $-\sin 0 = 0$, $-\cos 0 = -1$. Substituting into (0.3.7) with $x_0 = 0$ gives

$$
\cos x \approx 1 - \tfrac{1}{2} x^2.
$$

A small numerical check: $\cos(0.1) = 0.995004\ldots$, and $1 - 0.5(0.1)^2 = 0.995$. The approximation is excellent for small $x$ — this is exactly the regime in which the small-oscillation approximation holds in classical mechanics.

```python
import numpy as np
x: np.ndarray = np.linspace(-0.5, 0.5, 11)
exact = np.cos(x)
approx = 1 - 0.5 * x**2
print(np.max(np.abs(exact - approx)))  # ~3e-3 across the interval
```

## Partial derivatives and the gradient

A function $f$ of several variables, $f(x_1, \ldots, x_n)$, has $n$ **partial derivatives**:

$$
\frac{\partial f}{\partial x_i} \;=\; \lim_{h \to 0} \frac{f(x_1, \ldots, x_i + h, \ldots, x_n) - f(x_1, \ldots, x_n)}{h}. \tag{0.3.8}
$$

The partial derivative with respect to $x_i$ treats every other variable as a constant. The collection of partial derivatives, arranged as a vector, is the **gradient**:

$$
\nabla f(\mathbf{x}) = \begin{pmatrix} \partial f / \partial x_1 \\ \vdots \\ \partial f / \partial x_n \end{pmatrix}. \tag{0.3.9}
$$

### What the gradient means

Two interpretations are essential.

**1. Direction of steepest ascent.** $\nabla f(\mathbf{x})$ points in the direction along which $f$ increases most rapidly from $\mathbf{x}$, and its magnitude is the rate of increase in that direction. Consequently, $-\nabla f$ is the direction of steepest descent — the basis of every gradient-descent optimiser.

**2. Linear approximation.** The multivariable analogue of the tangent-line approximation is
$$
f(\mathbf{x} + \mathbf{h}) \approx f(\mathbf{x}) + \nabla f(\mathbf{x}) \cdot \mathbf{h}, \tag{0.3.10}
$$
valid for small $\mathbf{h}$. The gradient is the unique vector that makes this first-order approximation accurate.

In materials science, the gradient of the potential energy with respect to atomic positions is — up to a sign — the force:

$$
\mathbf{F}_i = -\nabla_{\mathbf{r}_i} U(\mathbf{r}_1, \ldots, \mathbf{r}_N). \tag{0.3.11}
$$

Every molecular-dynamics integrator, every geometry optimiser, every machine-learning interatomic potential ultimately reduces to evaluating this expression cheaply and accurately.

### Directional derivative

The **directional derivative** of $f$ at $\mathbf{x}$ along a unit vector $\hat{\mathbf{u}}$ is

$$
D_{\hat{\mathbf{u}}} f(\mathbf{x}) \;=\; \lim_{h \to 0} \frac{f(\mathbf{x} + h \hat{\mathbf{u}}) - f(\mathbf{x})}{h} \;=\; \nabla f(\mathbf{x}) \cdot \hat{\mathbf{u}}. \tag{0.3.12}
$$

It measures the rate of change of $f$ along the chosen direction. The maximum value over unit directions, by the Cauchy–Schwarz inequality applied to (0.3.12), is $\lVert \nabla f \rVert$, achieved when $\hat{\mathbf{u}}$ is parallel to $\nabla f$. This is a one-line proof that the gradient is the direction of steepest ascent.

## The Laplacian

A second-order differential operator of great physical importance is the **Laplacian**,

$$
\nabla^2 f \;=\; \sum_{i=1}^{n} \frac{\partial^2 f}{\partial x_i^2}. \tag{0.3.13}
$$

In three dimensions $\nabla^2 = \partial_x^2 + \partial_y^2 + \partial_z^2$. It appears in the time-independent Schrödinger equation,

$$
\left[ -\frac{\hbar^2}{2m} \nabla^2 + V(\mathbf{r}) \right] \psi(\mathbf{r}) = E\, \psi(\mathbf{r}), \tag{0.3.14}
$$

and in the Poisson equation for the electrostatic potential of a charge density, $\nabla^2 \phi = -\rho/\varepsilon_0$. The Laplacian measures the deviation of $f$ from its local average: $\nabla^2 f > 0$ at a local minimum, $< 0$ at a local maximum.

## Integrals

A definite integral $\int_a^b f(x) \, \mathrm{d} x$ may be visualised as the signed area between the graph of $f$ and the $x$-axis over $[a, b]$. The rigorous definition via **Riemann sums** is: partition $[a, b]$ into $N$ subintervals of width $\Delta x_k$, choose a sample point $x_k^\star$ in each, and define

$$
\int_a^b f(x)\, \mathrm{d} x \;=\; \lim_{N \to \infty, \, \max \Delta x_k \to 0} \sum_{k=1}^{N} f(x_k^\star)\, \Delta x_k, \tag{0.3.15}
$$

whenever this limit exists and is independent of the partition. Functions for which this works are **Riemann integrable**; continuous functions on bounded intervals always are.

Numerical integration in this book reduces to evaluating sums like (0.3.15) on a finite grid. The simplest **rectangle rule** is

```python
import numpy as np
def integrate_rectangle(f, a: float, b: float, n: int) -> float:
    x = np.linspace(a, b, n, endpoint=False) + (b - a) / (2 * n)
    return float(np.sum(f(x)) * (b - a) / n)

print(integrate_rectangle(np.sin, 0.0, np.pi, 1000))  # ≈ 2.0
```

Better rules (trapezoidal, Simpson, Gauss quadrature) appear in Chapter 1.

## The fundamental theorem of calculus

Differentiation and integration are inverse operations. The **fundamental theorem of calculus** has two parts.

**Part I.** If $f$ is continuous on $[a, b]$, the function
$$
F(x) = \int_a^x f(t)\, \mathrm{d} t
$$
is differentiable on $(a, b)$ with $F'(x) = f(x)$.

**Part II.** If $F$ is any antiderivative of $f$ on $[a, b]$, then
$$
\int_a^b f(x)\, \mathrm{d} x = F(b) - F(a). \tag{0.3.16}
$$

This is why a table of derivatives is also a table of integrals: $\int e^{ax} \mathrm{d} x = e^{ax}/a$ because $\frac{\mathrm{d}}{\mathrm{d} x} e^{ax}/a = e^{ax}$.

## Integration by parts

From the product rule (0.3.3), $(uv)' = u'v + uv'$. Integrating both sides over $[a, b]$ and rearranging gives the **integration-by-parts** formula:

$$
\int_a^b u(x)\, v'(x)\, \mathrm{d} x = \big[ u(x)\, v(x) \big]_a^b - \int_a^b u'(x)\, v(x)\, \mathrm{d} x. \tag{0.3.17}
$$

This is the algebraic identity that allows derivatives in quantum-mechanical matrix elements to be moved from one wavefunction onto another. In Chapter 4 you will see expressions of the form

$$
\int \psi^*(x) \left( -\frac{\hbar^2}{2m} \frac{\mathrm{d}^2 \psi}{\mathrm{d} x^2} \right) \mathrm{d} x,
$$

and integration by parts, applied twice with the assumption that $\psi \to 0$ at infinity, rewrites this as
$$
\frac{\hbar^2}{2m} \int \left| \frac{\mathrm{d} \psi}{\mathrm{d} x} \right|^2 \mathrm{d} x,
$$
a manifestly positive quantity. That positivity is the mathematical reason kinetic energy expectation values are non-negative.

## Variational thinking

A great deal of modern materials physics is phrased as a **variational principle**: a physical state is the one that minimises (or makes stationary) some functional of a function. The quintessential example is the variational principle of quantum mechanics, which we will use in earnest in Chapter 4 and Chapter 5.

A **functional** $E[\psi]$ assigns a number to each function $\psi$. The textbook example is

$$
E[\psi] \;=\; \int \psi^*(\mathbf{r})\, \hat H\, \psi(\mathbf{r})\, \mathrm{d}\tau, \tag{0.3.18}
$$

the expectation value of the Hamiltonian $\hat H$ in the state $\psi$. The variational theorem states that the ground-state energy $E_0$ satisfies

$$
E_0 \;=\; \min_{\psi} \frac{E[\psi]}{\int |\psi|^2 \, \mathrm{d}\tau}, \tag{0.3.19}
$$

where the minimisation runs over all admissible wavefunctions. Equivalently, minimise $E[\psi]$ subject to the **normalisation constraint** $\int |\psi|^2 \mathrm{d}\tau = 1$.

The Lagrange-multiplier prescription handles the constraint. Form
$$
L[\psi, \lambda] = E[\psi] - \lambda \left( \int |\psi|^2 \mathrm{d}\tau - 1 \right),
$$
and require its **functional derivative** with respect to $\psi^*$ to vanish. The result is

$$
\hat H \psi = \lambda \, \psi, \tag{0.3.20}
$$

the time-independent Schrödinger equation with $\lambda$ identified as the eigenvalue $E$. The variational principle thus *derives* the eigenvalue problem of Section 0.2; the Lagrange multiplier *is* the energy. Pause to appreciate this: a question about the minimum of an integral becomes an eigenvalue problem for a linear operator. This is one of the deep bridges between calculus and linear algebra, and it underlies almost everything in Chapters 4 and 5.

We will not develop calculus of variations rigorously here; for our purposes the operational rule is that $\delta E / \delta \psi^* = 0$ behaves exactly like an ordinary derivative being set to zero, with $\psi$ and $\psi^*$ treated as independent variables.

## Numerical differentiation

In practice, derivatives of functions defined only by tables (or by expensive black-box codes) must be estimated numerically. The two standard schemes are derived directly from the definition (0.3.1).

**Forward difference.** With step size $h$,
$$
f'(x) \approx \frac{f(x + h) - f(x)}{h} \;+\; O(h). \tag{0.3.21}
$$

**Central difference.**
$$
f'(x) \approx \frac{f(x + h) - f(x - h)}{2h} \;+\; O(h^2). \tag{0.3.22}
$$

Central differences are usually preferable: same number of function evaluations, one order higher accuracy. The error analysis comes from Taylor-expanding $f(x \pm h)$ to third order and observing that the odd-order terms cancel.

```python
import numpy as np

def forward_diff(f, x: float, h: float = 1e-5) -> float:
    return (f(x + h) - f(x)) / h

def central_diff(f, x: float, h: float = 1e-5) -> float:
    return (f(x + h) - f(x - h)) / (2 * h)

x0: float = 1.0
print(forward_diff(np.sin, x0))   # ≈ cos(1) = 0.5403
print(central_diff(np.sin, x0))   # ≈ 0.5403, more accurate
print(np.cos(x0))                  # exact reference
```

!!! warning "Choosing $h$"
    Smaller $h$ reduces the truncation error in (0.3.21)–(0.3.22), but eventually floating-point round-off in the subtraction $f(x+h) - f(x)$ dominates and the answer gets worse. For double precision and central differences, $h \sim 10^{-5}$ is usually a sweet spot. This is one place where naive intuition ("smaller step is always better") fails badly.

## Where this is used

- Chapter 4 builds on (0.3.18)–(0.3.20) to introduce the variational method and the Schrödinger equation in earnest.
- Chapter 7 implements (0.3.11): forces are minus the gradient of the potential.
- Chapter 9 and Chapter 10 train neural networks by stochastic gradient descent on a loss function; every parameter update is an application of the chain rule.
- Numerical differentiation (0.3.21)–(0.3.22) reappears in Chapter 1 as a debugging tool: comparing analytic gradients against finite-difference reference values is the standard sanity check.

We have now mastered change in the real domain. The next section extends our analytical reach to the complex plane, where waves, oscillations, and reciprocal-space ideas live most naturally.
