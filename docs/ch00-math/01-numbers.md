# 0.1 Numbers, Sets, and Functions

The roadmap in the chapter overview promised that we would build mathematics from the ground up. The ground floor is the number systems themselves and the language of functions. Everything that follows — vectors, derivatives, Fourier transforms — is constructed out of these primitives. By the end of this section we will have arrived at the Arrhenius rate law and the Boltzmann factor, two expressions you will encounter in nearly every subsequent chapter.

## The number hierarchy

Mathematics is built up out of progressively richer number systems. Each one extends the previous to solve an equation that could not be solved before.

The **natural numbers** $\mathbb{N} = \{0, 1, 2, 3, \ldots\}$ are what you count with. Adding two natural numbers gives a natural number, but subtraction can take you outside the set: $3 - 5$ is not a natural number. To repair this we introduce the **integers** $\mathbb{Z} = \{\ldots, -2, -1, 0, 1, 2, \ldots\}$.

Integers are closed under addition, subtraction, and multiplication, but not division: $1 / 2$ is not an integer. The **rational numbers** $\mathbb{Q}$ are quotients $p / q$ of integers with $q \neq 0$. They are dense — between any two rationals there is another rational — but they are still incomplete. The diagonal of a unit square has length $\sqrt 2$, and Pythagoras' classic proof shows that no rational squares to $2$.

The **real numbers** $\mathbb{R}$ fill in the gaps. Intuitively, a real number is anything that can be approximated to arbitrary accuracy by a rational. Equivalently, $\mathbb{R}$ corresponds to the unbroken number line. Most quantities in physics — a bond length, a temperature, an energy — are modelled as real numbers.

Finally there are equations that have no real solution. The simplest is $x^2 = -1$. To handle these we introduce the imaginary unit $i$ with $i^2 = -1$ and the **complex numbers** $\mathbb{C} = \{a + bi : a, b \in \mathbb{R}\}$. Complex numbers turn out to be indispensable in quantum mechanics; we devote Section 0.4 to them.

The chain of containment is

$$
\mathbb{N} \subset \mathbb{Z} \subset \mathbb{Q} \subset \mathbb{R} \subset \mathbb{C}.
$$

!!! note "Floating-point arithmetic is not $\mathbb{R}$"
    A computer cannot store an arbitrary real number. The IEEE-754 double-precision format approximates $\mathbb{R}$ using $2^{64}$ discrete values, with relative precision around $10^{-16}$. This is usually invisible, but it explains why two mathematically equivalent expressions can give different numerical answers. We will return to this in Chapter 1.

## Sets and interval notation

A **set** is an unordered collection of distinct objects. We write $x \in S$ if $x$ is an element of $S$, and $S \subseteq T$ if every element of $S$ is also an element of $T$. The empty set is $\varnothing$.

For subsets of the real line we use **interval notation**:

- $[a, b] = \{x \in \mathbb{R} : a \le x \le b\}$ (closed interval, endpoints included);
- $(a, b) = \{x \in \mathbb{R} : a < x < b\}$ (open interval, endpoints excluded);
- $[a, b) = \{x \in \mathbb{R} : a \le x < b\}$ (half-open);
- $(-\infty, b]$ means $x \le b$; $(a, \infty)$ means $x > a$.

Infinity is never an element of an interval — it is a bound, not a number.

## Functions as mappings

A **function** $f \colon A \to B$ is a rule that assigns to each element $x$ of a set $A$ exactly one element $f(x)$ of a set $B$. The set $A$ is the **domain**; $B$ is the **codomain**. The set of values actually taken, $\{f(x) : x \in A\}$, is the **range** (or image).

Two pitfalls trip up beginners. First, a function is the rule together with its domain; the same formula can define different functions on different domains. Second, the codomain need not equal the range — only the range is what the function actually produces.

!!! example "Domain and range"
    The function $f(x) = \sqrt{x}$ has natural domain $[0, \infty)$ in $\mathbb{R}$ and range $[0, \infty)$. The function $g(x) = 1 / x$ has domain $\mathbb{R} \setminus \{0\}$ and range $\mathbb{R} \setminus \{0\}$. Trying to evaluate $g(0)$ is not "infinity"; it is undefined.

### Composition

If $f \colon A \to B$ and $g \colon B \to C$, the **composition** $g \circ f \colon A \to C$ is defined by

$$
(g \circ f)(x) = g(f(x)).
$$

Order matters: in general $g \circ f \neq f \circ g$. Composition is associative, $(h \circ g) \circ f = h \circ (g \circ f)$, which is why nested function calls in code can be written without parentheses-of-parentheses.

### Inverses

A function $f$ is **injective** (one-to-one) if $f(x_1) = f(x_2)$ implies $x_1 = x_2$, and **surjective** (onto) if its range equals its codomain. A function that is both is **bijective**, and only bijective functions have inverses.

The **inverse** $f^{-1}$ satisfies

$$
f^{-1}(f(x)) = x \quad \text{and} \quad f(f^{-1}(y)) = y.
$$

Graphically, the inverse is the reflection of $f$ across the line $y = x$.

!!! warning "Inverse versus reciprocal"
    $f^{-1}(x)$ denotes the inverse function, not $1/f(x)$. The notation is unfortunate but universal.

## Exponentials and logarithms

The exponential function $\exp(x) = e^x$, with $e \approx 2.71828$, is arguably the most important function in physics. It is defined by the series

$$
e^x = \sum_{n=0}^{\infty} \frac{x^n}{n!} = 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \cdots
\tag{0.1.1}
$$

and satisfies the multiplicative property

$$
e^{x + y} = e^x \, e^y. \tag{0.1.2}
$$

Its inverse is the **natural logarithm** $\ln \colon (0, \infty) \to \mathbb{R}$, with $\ln(e^x) = x$. From (0.1.2) we get the logarithmic identity

$$
\ln(xy) = \ln x + \ln y, \qquad \ln(x^n) = n \ln x. \tag{0.1.3}
$$

Bases other than $e$ are sometimes useful — $\log_{10}$ for orders of magnitude, $\log_2$ for information — but they are related by $\log_b x = \ln x / \ln b$, so $e$ is the only base we will need.

### Worked example: the Arrhenius rate law

In Chapter 7 (molecular dynamics) and Chapter 8 (statistical mechanics) you will meet a great many activated processes — atomic hopping, chemical reactions, defect formation. Their rates almost always follow the **Arrhenius law**:

$$
k(T) = A \exp\!\left(-\frac{E_\mathrm{a}}{k_\mathrm{B} T}\right), \tag{0.1.4}
$$

where $k$ is the rate constant, $T$ is absolute temperature, $E_\mathrm{a}$ is an activation energy, $k_\mathrm{B} \approx 1.381 \times 10^{-23}\,\mathrm{J/K}$ is Boltzmann's constant, and $A$ is a prefactor with the same units as $k$.

Taking the logarithm of both sides,

$$
\ln k = \ln A - \frac{E_\mathrm{a}}{k_\mathrm{B}} \cdot \frac{1}{T}. \tag{0.1.5}
$$

This is the equation of a straight line in the variables $\ln k$ versus $1/T$, with slope $-E_\mathrm{a} / k_\mathrm{B}$ and intercept $\ln A$. Plotting experimental rate data this way — an **Arrhenius plot** — lets you read off the activation energy from the slope. The trick of taking a logarithm to convert a multiplicative law into a straight line is one you will use repeatedly.

### Worked example: the Boltzmann factor

The probability that a system in thermal equilibrium at temperature $T$ occupies a microstate of energy $E$ is proportional to the **Boltzmann factor**

$$
p(E) \propto \exp\!\left(-\frac{E}{k_\mathrm{B} T}\right) = e^{-\beta E}, \qquad \beta \equiv \frac{1}{k_\mathrm{B} T}. \tag{0.1.6}
$$

We will derive this in Chapter 8. For now, note three properties that follow from elementary algebra:

1. Lower-energy states are exponentially more probable than higher-energy ones.
2. The ratio $p(E_1) / p(E_2) = e^{-\beta(E_1 - E_2)}$ depends only on the energy *difference*, not the absolute energy.
3. As $T \to \infty$ ($\beta \to 0$) all states become equally probable; as $T \to 0$ only the ground state survives.

These three statements alone explain a huge amount of materials physics, from defect concentrations to magnetic phase transitions.

## A short detour: proof by induction

Most of this book is computational, but every so often you will see an inductive argument. The principle of **mathematical induction** says that to prove a statement $P(n)$ for every natural number $n \ge n_0$, it suffices to prove:

- **Base case.** $P(n_0)$ is true.
- **Inductive step.** If $P(k)$ is true for some $k \ge n_0$, then $P(k+1)$ is also true.

These two pieces together imply $P(n)$ for all $n \ge n_0$, by a sort of mathematical domino effect.

!!! example "Sum of the first $n$ integers"
    Claim: $\sum_{i=1}^{n} i = \frac{n(n+1)}{2}$ for every $n \ge 1$.

    *Base case.* For $n = 1$, the left-hand side is $1$ and the right-hand side is $1 \cdot 2 / 2 = 1$. They agree.

    *Inductive step.* Assume the formula holds for $n = k$:
    $$
    \sum_{i=1}^{k} i = \frac{k(k+1)}{2}.
    $$
    Add $k+1$ to both sides:
    $$
    \sum_{i=1}^{k+1} i = \frac{k(k+1)}{2} + (k+1) = \frac{k(k+1) + 2(k+1)}{2} = \frac{(k+1)(k+2)}{2},
    $$
    which is the formula with $n = k + 1$. Done.

We will use induction sparingly — once to establish the binomial theorem, once when discussing recursion in Chapter 1, and never again in any serious way. But the *mode of reasoning* — establishing a base case and a step — recurs constantly in numerical methods, where you analyse the error per step and then propagate it.

## A small Python sanity check

The Arrhenius law makes a concrete numerical prediction. Let us verify with code that doubling the temperature near room temperature does *not* double the rate.

```python
import numpy as np

k_B: float = 8.617e-5  # eV / K, Boltzmann's constant in convenient units
E_a: float = 0.5       # eV, typical activation energy
A: float = 1.0e13      # 1 / s, attempt frequency

def arrhenius(T: float) -> float:
    """Return the Arrhenius rate constant k(T) in 1/s."""
    return A * np.exp(-E_a / (k_B * T))

T1: float = 300.0
T2: float = 600.0
print(f"k({T1} K) = {arrhenius(T1):.3e} 1/s")
print(f"k({T2} K) = {arrhenius(T2):.3e} 1/s")
print(f"ratio    = {arrhenius(T2) / arrhenius(T1):.3e}")
```

The ratio is roughly $10^4$, not $2$. Activated processes are exquisitely temperature-sensitive — a fact that will save you many hours of confused log-reading when you start running real molecular dynamics in Chapter 7.

## Where this is used

- The number-line intuition reappears in Section 0.2, where vectors generalise it to higher dimensions.
- Exponentials power Section 0.4 (Euler's formula $e^{i\theta}$) and Section 0.5 (the Gaussian distribution $e^{-x^2/2}$).
- The Arrhenius and Boltzmann expressions are the practical bridge to Chapters 7 and 8.

With numbers and functions in hand, we are ready to stack them into vectors and matrices.
