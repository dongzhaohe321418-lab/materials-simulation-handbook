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

### Why we obsess over set notation

A beginner might reasonably ask why we bother distinguishing $\mathbb{N}$ from $\mathbb{Z}$ from $\mathbb{R}$ from $\mathbb{C}$ when, in the end, every number we write down on a page or feed into a computer looks like a finite string of digits. The answer is that the *kind* of number a quantity takes encodes physical content.

Consider three concrete examples drawn from materials science. The number of atoms in a unit cell is an element of $\mathbb{N}$: it must be non-negative and it must be whole. A fractional coordinate inside the unit cell is an element of $[0, 1) \subset \mathbb{R}$: it is real-valued, bounded, and continuous. The phase of a Bloch wavefunction is an element of $\mathbb{C}$ with unit modulus, $e^{i\theta}$: real and imaginary parts together carry the wave's directional information.

Now the crucial point: **physical observables — quantities that can be measured with an experimental apparatus — must be real**. Energies, positions, temperatures, magnetic moments, all live in $\mathbb{R}$. Complex numbers enter our equations because they are the most economical book-keeping device for waves and rotations, but every measurement at the end of a calculation must produce a real number. This is why, in quantum mechanics (Chapter 4), an observable is represented by a *Hermitian* operator: such operators are guaranteed to have real eigenvalues (Section 0.2). The set-theoretic distinction $\mathbb{R} \subset \mathbb{C}$ is not pedantry; it is the mathematical statement that "what we measure" sits inside "what we compute".

!!! tip "A mental triage"
    When you read a new equation, ask yourself: which symbols here are integers (counts, quantum numbers, indices), which are real (energies, distances, temperatures), and which are complex (wavefunctions, phases, Fourier amplitudes)? Half of the bugs in scientific code come from a complex quantity being silently cast to a real one, or from a real index being computed as a float.

## Sets and interval notation

A **set** is an unordered collection of distinct objects. We write $x \in S$ if $x$ is an element of $S$, and $S \subseteq T$ if every element of $S$ is also an element of $T$. The empty set is $\varnothing$.

For subsets of the real line we use **interval notation**:

- $[a, b] = \{x \in \mathbb{R} : a \le x \le b\}$ (closed interval, endpoints included);
- $(a, b) = \{x \in \mathbb{R} : a < x < b\}$ (open interval, endpoints excluded);
- $[a, b) = \{x \in \mathbb{R} : a \le x < b\}$ (half-open);
- $(-\infty, b]$ means $x \le b$; $(a, \infty)$ means $x > a$.

Infinity is never an element of an interval — it is a bound, not a number.

### Set operations

Two sets can be combined in three standard ways. The **union** $S \cup T$ is the set of elements in $S$ or $T$ (or both); the **intersection** $S \cap T$ contains elements in both; the **set difference** $S \setminus T$ contains elements of $S$ that are not in $T$. The **Cartesian product** $S \times T$ is the set of ordered pairs $\{(s, t) : s \in S, t \in T\}$ — what gives us $\mathbb{R}^2$, $\mathbb{R}^3$, and so on as products of $\mathbb{R}$ with itself.

These operations underpin the language of probability events (Section 0.5) and the construction of configuration spaces in statistical mechanics (Chapter 8). When we say "the probability that an atom is in region $A$ *or* region $B$", we mean the probability of $A \cup B$. When we say "the joint probability of $A$ *and* $B$", we mean $P(A \cap B)$. When we describe the configuration space of $N$ atoms in a box, we mean $[0, L]^{3N}$ — a Cartesian product of intervals.

## Functions as mappings

A **function** $f \colon A \to B$ is a rule that assigns to each element $x$ of a set $A$ exactly one element $f(x)$ of a set $B$. The set $A$ is the **domain**; $B$ is the **codomain**. The set of values actually taken, $\{f(x) : x \in A\}$, is the **range** (or image).

Two pitfalls trip up beginners. First, a function is the rule together with its domain; the same formula can define different functions on different domains. Second, the codomain need not equal the range — only the range is what the function actually produces.

!!! example "Domain and range"
    The function $f(x) = \sqrt{x}$ has natural domain $[0, \infty)$ in $\mathbb{R}$ and range $[0, \infty)$. The function $g(x) = 1 / x$ has domain $\mathbb{R} \setminus \{0\}$ and range $\mathbb{R} \setminus \{0\}$. Trying to evaluate $g(0)$ is not "infinity"; it is undefined.

### Why domain bookkeeping matters in physics

Physical quantities come with implicit domains, and ignoring them is a classic source of bugs.

- A temperature $T$ in absolute units lives in $[0, \infty)$. Plugging negative temperatures into the Boltzmann factor $e^{-E/k_\mathrm{B} T}$ yields nonsense.
- A radial distance $r$ in spherical coordinates lives in $[0, \infty)$, while the polar angle $\theta$ lives in $[0, \pi]$ and the azimuth $\phi$ in $[0, 2\pi)$. Mixing these up produces wrong integration measures.
- A probability density is non-negative and integrates to one over its domain. Section 0.5 makes this precise.
- A wavefunction phase $\theta \in (-\pi, \pi]$ is a circular coordinate; subtracting two phases requires reducing the result mod $2\pi$.

The discipline of stating the domain of every function before manipulating it pays for itself many times over.

### Composition

If $f \colon A \to B$ and $g \colon B \to C$, the **composition** $g \circ f \colon A \to C$ is defined by

$$
(g \circ f)(x) = g(f(x)).
$$

Order matters: in general $g \circ f \neq f \circ g$. Composition is associative, $(h \circ g) \circ f = h \circ (g \circ f)$, which is why nested function calls in code can be written without parentheses-of-parentheses.

Composition is so fundamental that programming languages take it for granted. When you write `np.sqrt(np.abs(x))`, you are computing $(\sqrt{\cdot} \circ |\cdot|)(x)$. The implicit order — innermost function first — matches the mathematical convention. When debugging, mentally peel off the outermost function and ask "what does the inner expression produce?". A surprising fraction of bugs are caught by this single discipline. The whole field of *function-composition pipelines* in machine learning (PyTorch modules, JAX `vmap`/`grad` stacks, scikit-learn `Pipeline` objects) is built on this elementary observation.

### Inverses

A function $f$ is **injective** (one-to-one) if $f(x_1) = f(x_2)$ implies $x_1 = x_2$, and **surjective** (onto) if its range equals its codomain. A function that is both is **bijective**, and only bijective functions have inverses.

The **inverse** $f^{-1}$ satisfies

$$
f^{-1}(f(x)) = x \quad \text{and} \quad f(f^{-1}(y)) = y.
$$

Graphically, the inverse is the reflection of $f$ across the line $y = x$.

!!! warning "Inverse versus reciprocal"
    $f^{-1}(x)$ denotes the inverse function, not $1/f(x)$. The notation is unfortunate but universal.

### Why bijections matter in physics

In physics and chemistry, bijective maps appear under the name **diffeomorphisms** (smooth invertible maps) or, in the discrete case, **permutations**. They are the mathematical objects that translate one description of a system into an equally complete description.

For example, a change of coordinates from Cartesian $(x, y, z)$ to spherical $(r, \theta, \phi)$ is a bijection from $\mathbb{R}^3 \setminus \{0\}$ to $(0, \infty) \times (0, \pi) \times [0, 2\pi)$. Because the map is invertible, no information is lost; the description is just reorganised. The same holds for a change of basis in linear algebra (Section 0.2), a unitary transformation in quantum mechanics, and a canonical transformation in Hamiltonian mechanics. Each is a bijection on the relevant state space.

By contrast, *non*-bijective maps lose information. Projection onto a lower-dimensional subspace, integration over a coordinate, tracing out part of a quantum system — each compresses many states into fewer, and is irreversible. The distinction between bijective and non-bijective maps thus mirrors the distinction between *reversible* and *irreversible* operations in physics, an idea you will meet again in Chapter 8 in the discussion of entropy.

### Even, odd, and periodic functions

Three structural properties of functions are worth naming because they propagate into structure of Taylor series, Fourier series, and selection rules.

A function $f$ is **even** if $f(-x) = f(x)$ for all $x$ in its domain (the graph is symmetric about the $y$-axis). Examples: $x^2, \cos x, e^{-x^2}$. The Taylor series of an even function contains only even powers.

A function is **odd** if $f(-x) = -f(x)$ (the graph has $180^\circ$ rotational symmetry about the origin). Examples: $x, \sin x, \tan x$. The Taylor series of an odd function contains only odd powers, and $f(0) = 0$.

A function is **periodic** with period $L > 0$ if $f(x + L) = f(x)$ for all $x$. Examples: $\sin, \cos$, the electron density of a crystal. Periodic functions are the natural domain of Fourier series (Section 0.4), and crystalline solids — the central object of this book — are defined by periodicity of the atomic potential.

Selection rules in spectroscopy are nothing more than the observation that the integral $\int f(x)\, \mathrm{d} x$ over a symmetric interval vanishes when $f$ is odd. We will use this idea explicitly in Chapter 4 to argue why certain quantum transitions are forbidden.

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

### Intuition: why is $e$ special?

Among all bases, $e$ is the unique number for which the exponential function is its own derivative: $\frac{\mathrm{d}}{\mathrm{d} x} e^x = e^x$. Geometrically this means: at every point on the graph of $y = e^x$, the slope equals the height. Picture a curve so finely tuned that, no matter how high you climb, the steepness of the path is precisely your current altitude. This self-reproducing property is what makes $e^{-\beta E}$ the natural solution of the differential equation that governs equilibrium populations, and what makes $e^{i\omega t}$ the natural solution of the oscillator equation. Every other base, like $10$ or $2$, picks up an annoying factor of $\ln(\text{base})$ whenever you differentiate; with $e$ the factor is one.

This is not a coincidence but a *definition*: $e$ is the unique number that makes the exponential maximally simple under calculus. We will see in Section 0.3 that this simplicity propagates to every chain-rule computation in physics, and in Section 0.4 that it propagates to every Fourier-style decomposition.

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

!!! example "Arrhenius computation in convenient units"
    Pick $E_\mathrm{a} = 1\,\mathrm{eV}$, a typical activation energy for atomic diffusion in a solid. Use $k_\mathrm{B} = 8.617 \times 10^{-5}\,\mathrm{eV/K}$ — the same constant as before, just rescaled. We compute the dimensionless argument $E_\mathrm{a}/(k_\mathrm{B} T)$ at three temperatures:

    (1) At $T = 300\,\mathrm{K}$:
    $$
    \frac{E_\mathrm{a}}{k_\mathrm{B} T} = \frac{1}{8.617 \times 10^{-5} \times 300} = \frac{1}{0.02585} = 38.68.
    $$
    The Boltzmann factor is $e^{-38.68} \approx 1.6 \times 10^{-17}$.

    (2) At $T = 600\,\mathrm{K}$:
    $$
    \frac{E_\mathrm{a}}{k_\mathrm{B} T} = \frac{1}{8.617 \times 10^{-5} \times 600} = \frac{1}{0.0517} = 19.34.
    $$
    The Boltzmann factor is $e^{-19.34} \approx 4.0 \times 10^{-9}$.

    (3) At $T = 1000\,\mathrm{K}$:
    $$
    \frac{E_\mathrm{a}}{k_\mathrm{B} T} = \frac{1}{0.0862} = 11.60.
    $$
    The Boltzmann factor is $e^{-11.60} \approx 9.1 \times 10^{-6}$.

    Now the rate ratios with respect to $T = 300\,\mathrm{K}$:
    $$
    \frac{k(600)}{k(300)} = \frac{4.0 \times 10^{-9}}{1.6 \times 10^{-17}} \approx 2.5 \times 10^{8},
    $$
    $$
    \frac{k(1000)}{k(300)} = \frac{9.1 \times 10^{-6}}{1.6 \times 10^{-17}} \approx 5.7 \times 10^{11}.
    $$

    Doubling the temperature from $300\,\mathrm{K}$ to $600\,\mathrm{K}$ speeds the process up by **eight orders of magnitude**. Going to $1000\,\mathrm{K}$ adds another three orders. This is why annealing experiments work, why high-temperature simulations equilibrate fast, and why a $50\,\mathrm{K}$ drift in your thermostat can render an MD run physically meaningless.

!!! note "Why this step?"
    The logarithm in (0.1.5) is not magic; it is exploiting the algebraic property $\ln(e^x) = x$. Whenever a model predicts an *exponentially* multiplicative dependence on a parameter, taking the log linearises it. Linear models are statistically and visually easier to fit. Every modern materials informatics paper that plots $\ln \tau$ versus $1/T$, or $\log \sigma$ versus $\log E$, is using the same trick.

### Worked example: the Boltzmann factor

The probability that a system in thermal equilibrium at temperature $T$ occupies a microstate of energy $E$ is proportional to the **Boltzmann factor**

$$
p(E) \propto \exp\!\left(-\frac{E}{k_\mathrm{B} T}\right) = e^{-\beta E}, \qquad \beta \equiv \frac{1}{k_\mathrm{B} T}. \tag{0.1.6}
$$

We will derive this in Chapter 8 (and again in Section 0.5 as a maximum-entropy distribution). For now, note three properties that follow from elementary algebra:

1. Lower-energy states are exponentially more probable than higher-energy ones.
2. The ratio $p(E_1) / p(E_2) = e^{-\beta(E_1 - E_2)}$ depends only on the energy *difference*, not the absolute energy.
3. As $T \to \infty$ ($\beta \to 0$) all states become equally probable; as $T \to 0$ only the ground state survives.

These three statements alone explain a huge amount of materials physics, from defect concentrations to magnetic phase transitions.

!!! tip "Sanity check on the limits"
    The limit $T \to \infty$ giving uniform probability is intuitive: an enormous reservoir of thermal energy washes out all energetic preferences. The limit $T \to 0$ is the *third law of thermodynamics* in action: as the temperature drops to absolute zero, the system collapses into its ground state, and entropy approaches zero (for a non-degenerate ground state). Whenever you write down a probability that depends on $T$, evaluate these two limits as a free sanity check.

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

!!! example "Sum of the first $n$ squares"
    Claim: $\sum_{i=1}^{n} i^2 = \frac{n(n+1)(2n+1)}{6}$ for every $n \ge 1$.

    *Base case.* For $n = 1$, the left-hand side is $1^2 = 1$. The right-hand side is $1 \cdot 2 \cdot 3 / 6 = 1$. They agree.

    *Inductive step.* Assume $\sum_{i=1}^{k} i^2 = k(k+1)(2k+1)/6$. We compute
    $$
    \sum_{i=1}^{k+1} i^2 = \frac{k(k+1)(2k+1)}{6} + (k+1)^2.
    $$
    Factor out $(k+1)$:
    $$
    = \frac{(k+1)\big[ k(2k+1) + 6(k+1) \big]}{6}
    = \frac{(k+1)(2k^2 + 7k + 6)}{6}.
    $$
    The quadratic factors as $(k+2)(2k+3)$, giving
    $$
    \sum_{i=1}^{k+1} i^2 = \frac{(k+1)(k+2)(2k+3)}{6},
    $$
    which is the original formula with $n \to k+1$. Done.

    !!! note "Why this step?"
        Factoring $(k+1)$ out of two terms whose common factor is not visually obvious is the kind of algebraic step that gets glossed over in textbooks. The trick is to write both terms with the same denominator and inspect the numerator. We will use the same factoring pattern in Section 0.3 when collecting Taylor-series terms.

!!! example "Sum of a geometric series"
    Claim: for any $r \ne 1$ and $n \ge 0$,
    $$
    \sum_{i=0}^{n} r^i = \frac{1 - r^{n+1}}{1 - r}.
    $$

    *Base case.* For $n = 0$, the left-hand side is $r^0 = 1$. The right-hand side is $(1 - r)/(1 - r) = 1$. They agree.

    *Inductive step.* Assume $\sum_{i=0}^{k} r^i = (1 - r^{k+1})/(1 - r)$. Add $r^{k+1}$ to both sides:
    $$
    \sum_{i=0}^{k+1} r^i = \frac{1 - r^{k+1}}{1 - r} + r^{k+1}
    = \frac{1 - r^{k+1} + r^{k+1}(1 - r)}{1 - r}
    = \frac{1 - r^{k+2}}{1 - r}.
    $$
    This is the formula with $n \to k+1$. Done.

    The geometric series is the simplest non-trivial infinite sum and underlies the radius of convergence of every power series we will encounter — including the exponential, sine, and cosine series.

We will use induction sparingly — once to establish the binomial theorem, once when discussing recursion in Chapter 1, and never again in any serious way. But the *mode of reasoning* — establishing a base case and a step — recurs constantly in numerical methods, where you analyse the error per step and then propagate it. Convergence proofs for time-integration schemes (Chapter 7) and for stochastic gradient descent (Chapter 9) are recognisably inductive in spirit.

## Asymptotic notation

Throughout this book we describe how quantities scale with system size or precision using **big-$O$ notation**. We write $f(n) = O(g(n))$ as $n \to \infty$ if there exist constants $C$ and $n_0$ such that $|f(n)| \le C |g(n)|$ for all $n \ge n_0$. Informally: $f$ grows no faster than $g$ up to a constant factor.

The notation lets us compare algorithms without committing to constants that depend on hardware. Computing the Fourier transform of $N$ samples by the definition is $O(N^2)$; the FFT is $O(N \log N)$. For $N = 10^6$, this is the difference between $10^{12}$ operations (hours) and $10^7 \cdot 20 = 2 \times 10^8$ (seconds). In Chapter 1 we will see DFT scaling as $O(N^3)$ in the number of bands and Monte Carlo error scaling as $O(1/\sqrt{M})$ in the number of samples — both consequences of the same big-$O$ language.

A close relative is little-$o$: $f(n) = o(g(n))$ if $f(n)/g(n) \to 0$. This says $f$ grows *strictly slower* than $g$. We will see this in the Taylor-series remainders of Section 0.3.

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

## The hyperbolic functions

A pair of close relatives of the exponential, used heavily in statistical mechanics and machine learning, are the **hyperbolic** functions
$$
\cosh x = \frac{e^x + e^{-x}}{2}, \qquad \sinh x = \frac{e^x - e^{-x}}{2}, \qquad \tanh x = \frac{\sinh x}{\cosh x}.
$$
These are even, odd, and odd respectively. They satisfy the identity $\cosh^2 x - \sinh^2 x = 1$, a sign-flipped echo of the trigonometric identity. They arise naturally because they are the real and imaginary parts of $e^x$ extended to imaginary argument: $\cos\theta = \cosh(i\theta)$ and $i \sin\theta = \sinh(i\theta)$. We will see this connection in Section 0.4.

In statistical mechanics, $\tanh(\beta J)$ appears as the magnetisation of the one-dimensional Ising model. In machine learning, $\tanh$ and the related sigmoid $\sigma(x) = 1/(1 + e^{-x}) = (1 + \tanh(x/2))/2$ are the classical activation functions of neural networks. Get used to seeing exponentials in disguise.

!!! example "Limits of $\tanh$"
    The hyperbolic tangent saturates: as $x \to +\infty$, $\tanh x \to 1$; as $x \to -\infty$, $\tanh x \to -1$. Near zero, $\tanh x \approx x - x^3/3$. This S-shape is what makes it useful as an activation: it is approximately linear for small inputs and bounded for large ones, suppressing run-away gradients.

## Limits and continuity, briefly

We have used the limit symbol $\lim_{x \to a} f(x) = L$ informally; let us pin down the meaning. Intuitively, $L$ is the value $f(x)$ approaches as $x$ approaches $a$ — without ever necessarily reaching $a$. Formally, $\lim_{x \to a} f(x) = L$ if, for every $\varepsilon > 0$, there exists $\delta > 0$ such that $0 < |x - a| < \delta$ implies $|f(x) - L| < \varepsilon$. This is the famous $(\varepsilon, \delta)$ definition that haunts undergraduate analysis courses; we will rarely need it, but its conceptual content is worth understanding.

A function $f$ is **continuous** at $a$ if $\lim_{x \to a} f(x) = f(a)$. This requires three things: the limit exists, $f(a)$ is defined, and the two agree. Most "nice" functions (polynomials, exponentials, $\sin, \cos$, compositions of these) are continuous on their natural domains. Discontinuities occur at things like step functions $\Theta(x)$ and rational functions at their poles.

Continuity matters in numerical analysis. Differentiation amplifies discontinuities: even if $f$ is continuous, $f'$ may have jumps, and $f''$ may diverge. This is why we use smooth interpolation schemes (splines, smoothing kernels) for tabulated potentials in Chapter 7 rather than piecewise-linear ones; the integrator equations involve forces (gradients of energies), and discontinuous forces breed numerical chaos.

## Where this is used

- The number-line intuition reappears in Section 0.2, where vectors generalise it to higher dimensions.
- Exponentials power Section 0.4 (Euler's formula $e^{i\theta}$) and Section 0.5 (the Gaussian distribution $e^{-x^2/2}$).
- The Arrhenius and Boltzmann expressions are the practical bridge to Chapters 7 and 8.
- Set-theoretic distinctions $\mathbb{R} \subset \mathbb{C}$ underpin the proof in Section 0.2 that Hermitian operators have real eigenvalues.
- The induction template recurs in convergence arguments for Chapter 7 integrators and Chapter 9 optimisers.

With numbers and functions in hand, we are ready to stack them into vectors and matrices.
