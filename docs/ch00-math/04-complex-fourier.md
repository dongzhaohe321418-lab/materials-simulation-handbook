# 0.4 Complex Numbers and Fourier Intuition

Section 0.3 ended by promising that complex numbers would extend our analytic reach. They do more than that: complex numbers are the natural language for waves, oscillations, and — crucially — for the reciprocal-space picture of crystals. By the end of this section the symbol $e^{i\mathbf{k}\cdot\mathbf{r}}$, which will appear hundreds of times in Chapters 4 to 6, will read like ordinary algebra.

## The imaginary unit

The equation $x^2 = -1$ has no real solution. To proceed, we postulate a number $i$ satisfying

$$
i^2 = -1, \tag{0.4.1}
$$

and form the set of **complex numbers**

$$
\mathbb{C} = \{ a + b i : a, b \in \mathbb{R} \}. \tag{0.4.2}
$$

A complex number $z = a + b i$ has a **real part** $\operatorname{Re} z = a$ and an **imaginary part** $\operatorname{Im} z = b$. Addition and multiplication follow the usual algebraic rules with the substitution $i^2 = -1$:

$$
(a + b i)(c + d i) = (ac - bd) + (ad + bc) i.
$$

The **complex conjugate** of $z$ is $\bar z = a - b i$, and the **modulus** is

$$
|z| = \sqrt{\bar z z} = \sqrt{a^2 + b^2}. \tag{0.4.3}
$$

Geometrically, $\mathbb{C}$ is the plane: $z = a + b i$ corresponds to the point $(a, b)$. The modulus is then the Euclidean distance from the origin.

!!! tip "Complex numbers are 2D vectors with multiplication"
    The set $\mathbb{C}$ is, as an additive group, identical to $\mathbb{R}^2$. The novelty is that $\mathbb{C}$ comes equipped with a *multiplication* operation, with two-dimensional vectors do not have. Multiplying by $i$, in particular, rotates the plane by $90^\circ$, and multiplying by a unit-modulus number $e^{i\theta}$ rotates by $\theta$. This is what makes $\mathbb{C}$ the natural arena for rotations, oscillations, and waves. Vectors in $\mathbb{R}^2$ cannot rotate themselves; complex numbers can.

### Sanity check: $i^2 = -1$ as a rotation

Visualise $i = (0, 1)$ as the unit vector pointing along the imaginary axis. Multiplying $i$ by itself rotates that vector by $90^\circ$ a second time, landing on $(-1, 0) = -1$. Geometrically this is obvious; algebraically it required postulating $i^2 = -1$. The two descriptions agree, and the agreement is the whole content of "complex multiplication encodes rotations".

## Polar form and the argument

The same point $(a, b)$ in the plane can be written in polar coordinates,

$$
a = r \cos\theta, \qquad b = r \sin\theta,
$$

with $r = |z|$ and $\theta = \arg z$ the **argument**, conventionally in $(-\pi, \pi]$. Hence

$$
z = r(\cos\theta + i \sin\theta). \tag{0.4.4}
$$

Multiplication in polar form is gorgeously simple. If $z_1 = r_1(\cos\theta_1 + i\sin\theta_1)$ and $z_2 = r_2(\cos\theta_2 + i\sin\theta_2)$, then by the addition formulae,

$$
z_1 z_2 = r_1 r_2 \big[ \cos(\theta_1 + \theta_2) + i \sin(\theta_1 + \theta_2) \big]. \tag{0.4.5}
$$

Moduli multiply, arguments add. Multiplication by a complex number of unit modulus is a rotation of the plane.

## Euler's formula

The cleanest way to express (0.4.4) and (0.4.5) is via **Euler's formula**:

$$
e^{i\theta} = \cos\theta + i \sin\theta. \tag{0.4.6}
$$

This is not a mere abbreviation: it is a genuine identity that follows from the power-series definitions of the exponential, sine, and cosine.

### Derivation 1: via Taylor series

Recall the series (0.1.1) for the exponential:

$$
e^{x} = \sum_{n=0}^{\infty} \frac{x^n}{n!}.
$$

We extend it to complex arguments by substituting $x \to i\theta$:

$$
e^{i\theta} = \sum_{n=0}^{\infty} \frac{(i\theta)^n}{n!}.
$$

The powers of $i$ cycle: $i^0 = 1$, $i^1 = i$, $i^2 = -1$, $i^3 = -i$, $i^4 = 1$, and so on. Separating the sum by parity:

$$
e^{i\theta} = \sum_{k=0}^{\infty} \frac{(-1)^k \theta^{2k}}{(2k)!} \;+\; i \sum_{k=0}^{\infty} \frac{(-1)^k \theta^{2k+1}}{(2k+1)!}.
$$

The two series on the right are the Taylor series of $\cos\theta$ and $\sin\theta$ respectively. Done.

### Derivation 2: via the defining ODE

A second derivation, even quicker once you know calculus, treats $f(\theta) = e^{i\theta}$ as the unique solution of a differential equation.

**Step (1).** By the chain rule, if $f(\theta) = e^{i\theta}$ then $f'(\theta) = i\, e^{i\theta} = i\, f(\theta)$, with $f(0) = 1$.

**Step (2).** Now propose $g(\theta) = \cos\theta + i \sin\theta$ as a candidate solution to the same ODE. Differentiate: $g'(\theta) = -\sin\theta + i \cos\theta = i (\cos\theta + i \sin\theta) = i\, g(\theta)$. And $g(0) = \cos 0 + i \sin 0 = 1$.

**Step (3).** Both $f$ and $g$ satisfy the same first-order ODE with the same initial condition. By the uniqueness theorem for linear ODEs, $f = g$:
$$
e^{i\theta} = \cos\theta + i\sin\theta. \quad \square
$$

!!! note "Why this step?"
    Uniqueness of solutions to first-order linear ODEs is a theorem you may take on faith here (we will not prove it in this book). Operationally: if two functions both satisfy $y' = i y$ and both equal $1$ at $\theta = 0$, they must agree everywhere. This kind of "uniqueness pinning" is a workhorse technique in mathematical physics.

The two derivations complement one another: the Taylor-series approach reveals *why* the formula has the structure it does (the alternating powers of $i$ separate sines and cosines), while the ODE approach explains *what makes it natural* (it is the simplest solution of the simplest complex ODE).

### Euler's identity

A particularly celebrated consequence, obtained by setting $\theta = \pi$, is

$$
e^{i\pi} + 1 = 0,
$$

which combines five fundamental constants in a single equation: $0, 1, e, i, \pi$. This is sometimes called **Euler's identity** and is widely regarded as one of the most beautiful formulas in mathematics. Historically Euler published Euler's formula (0.4.6) in his *Introductio in analysin infinitorum* (1748), well before complex analysis was rigorously codified. He had no formal justification for substituting an imaginary argument into a real power series — his motivation was purely the structural elegance of the result. Two centuries later, we still teach his derivation almost unchanged because nothing simpler has ever been found.

### Consequences

Using (0.4.6), the polar form (0.4.4) becomes $z = r e^{i\theta}$, the multiplication rule (0.4.5) becomes the trivial $e^{i\theta_1} e^{i\theta_2} = e^{i(\theta_1 + \theta_2)}$, and the conjugate is $\bar z = r e^{-i\theta}$. Sines and cosines themselves are linear combinations of complex exponentials:

$$
\cos\theta = \frac{e^{i\theta} + e^{-i\theta}}{2}, \qquad \sin\theta = \frac{e^{i\theta} - e^{-i\theta}}{2i}. \tag{0.4.7}
$$

These identities turn trigonometric algebra — addition formulae, product-to-sum identities — into one-line manipulations of exponentials. We will use them constantly when discussing waves.

```python
import numpy as np

theta: float = np.pi / 3
z = np.exp(1j * theta)
print(z)                       # 0.5 + 0.866i
print(np.cos(theta) + 1j * np.sin(theta))  # the same
print(abs(z))                  # 1.0
```

### Complex roots of polynomials

The introduction of $\mathbb{C}$ pays an enormous algebraic dividend: the **fundamental theorem of algebra** states that every non-constant polynomial with complex coefficients has at least one complex root. By repeated application, a degree-$n$ polynomial has exactly $n$ complex roots (counted with multiplicity).

This is the algebraic reason a $2 \times 2$ matrix always has exactly two eigenvalues over $\mathbb{C}$ — they are the roots of a quadratic. The Hermitian matrices of physics are special in that their eigenvalues happen to all be real, but the surrounding mathematical scaffolding lives in $\mathbb{C}$. Whenever we say "diagonalise a matrix", we are implicitly working over $\mathbb{C}$ to guarantee the eigenvalues exist.

## Why complex numbers are natural in quantum mechanics

The state of a quantum system is a **complex-valued** wavefunction $\psi(\mathbf{r}, t)$, evolving according to the time-dependent Schrödinger equation

$$
i \hbar \frac{\partial \psi}{\partial t} = \hat{H} \psi.
$$

The $i$ on the left-hand side is not optional. If $\psi$ were forced to be real, $\hat{H} \psi$ would have to be purely imaginary at every instant, conflicting with the requirement that $\hat{H}$ — the energy — is a real-valued observable. The complex phase is what stores the dynamical information; the squared modulus $|\psi|^2$ is what we measure.

For a stationary state of energy $E$, the time dependence is a pure complex exponential,

$$
\psi(\mathbf{r}, t) = \psi(\mathbf{r})\, e^{-i E t / \hbar}, \tag{0.4.8}
$$

and the probability density $|\psi(\mathbf{r}, t)|^2 = |\psi(\mathbf{r})|^2$ is independent of time. The complex phase rotates uniformly; the observable quantities do not. This stationary-phase picture, together with the Bloch theorem we will see in Chapter 3.5, is the bridge from quantum mechanics to band structure.

A second reason: **plane waves**. A free particle of momentum $\hbar \mathbf{k}$ has wavefunction $e^{i\mathbf{k}\cdot\mathbf{r}}$. Real sines and cosines cannot represent a wave travelling in a definite direction; they always contain a counter-propagating component (see (0.4.7)). Plane-wave bases — the natural basis for DFT in periodic solids — are intrinsically complex.

!!! tip "The phase carries the physics"
    A common student error is to discard the complex phase as "just a normalisation" and work with real wavefunctions only. This works for the ground state of a real Hamiltonian (which can be chosen real) but fails the moment time evolution, currents, or magnetic fields enter the picture. Probability currents are bilinear in $\psi$ and $\psi^*$: $\mathbf{j} = (\hbar / m) \mathrm{Im}(\psi^* \nabla \psi)$. A real $\psi$ has $\mathbf{j} = \mathbf{0}$ identically. The phase is what makes a wavefunction move.

## Fourier series — a paragraph of motivation

A periodic function $f$ of period $L$ can be expanded as a sum of sines and cosines of frequencies that are integer multiples of $2\pi/L$. Writing this in complex form via (0.4.7),

$$
f(x) = \sum_{n=-\infty}^{\infty} c_n\, e^{i\, 2\pi n x / L}, \qquad c_n = \frac{1}{L} \int_0^L f(x)\, e^{-i\, 2\pi n x / L} \, \mathrm{d} x. \tag{0.4.9}
$$

The coefficients $c_n$ tell you "how much of each frequency" $f$ contains. In a crystal of lattice constant $L$, periodic functions — for instance the electron density — admit exactly this expansion, with $2\pi n / L$ the reciprocal-lattice wavevectors. This is the algebraic content of the reciprocal-lattice picture you will meet in Chapter 3.

### Worked example: Fourier series of a square wave

Take the square wave on $[-\pi, \pi]$ with $f(x) = +1$ for $x > 0$ and $f(x) = -1$ for $x < 0$, extended periodically with period $L = 2\pi$. By the symmetry $f(-x) = -f(x)$, the series contains only sine terms. A direct computation of the coefficients yields the real Fourier series

$$
f(x) = \frac{4}{\pi} \sum_{k=0}^{\infty} \frac{\sin\big((2k+1) x\big)}{2k+1}.
$$

The first four partial sums are

- $S_1(x) = \frac{4}{\pi} \sin x$,
- $S_2(x) = \frac{4}{\pi}\left[ \sin x + \frac{\sin 3x}{3} \right]$,
- $S_3(x) = \frac{4}{\pi}\left[ \sin x + \frac{\sin 3x}{3} + \frac{\sin 5x}{5} \right]$,
- $S_4(x) = \frac{4}{\pi}\left[ \sin x + \frac{\sin 3x}{3} + \frac{\sin 5x}{5} + \frac{\sin 7x}{7} \right]$.

At $x = \pi/2$, these evaluate to $4/\pi \approx 1.273$, $1.273 - 0.424 = 0.849$ wait — let's redo this. With $\sin(3\pi/2) = -1$: $S_2(\pi/2) = (4/\pi)(1 - 1/3) = (4/\pi)(2/3) \approx 0.849$. $S_3(\pi/2)$: $\sin(5\pi/2) = 1$, so $S_3 = (4/\pi)(1 - 1/3 + 1/5) \approx 1.103$. $S_4(\pi/2)$: $\sin(7\pi/2) = -1$, so $S_4 = (4/\pi)(1 - 1/3 + 1/5 - 1/7) \approx 0.922$.

The successive partial sums oscillate above and below the target value $1$, converging slowly toward it. Each new term knocks off a bit of error but introduces small ripples near the discontinuities at $x = 0, \pm\pi$. These ripples never vanish: their amplitude stays roughly $9\%$ of the jump height, and they merely shift closer to the discontinuity as more terms are added. This is the **Gibbs phenomenon**, and it is a generic feature of Fourier series of discontinuous functions. In materials science it appears whenever you try to represent a step-function potential (a sharp band edge, a hard wall) with a truncated plane-wave basis — the wiggles are real, and they motivate the use of smooth pseudopotentials in Chapter 5.

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-np.pi, np.pi, 2001)
plt.axhline(1, color="grey", lw=0.5)
plt.axhline(-1, color="grey", lw=0.5)
for K in [1, 2, 3, 4, 20]:
    series = sum(np.sin((2*k+1) * x) / (2*k+1) for k in range(K))
    plt.plot(x, (4/np.pi) * series, label=f"K={K}")
plt.legend()
plt.show()
```

Run this and watch the partial sums creep toward the square wave: each new harmonic adds a faster wiggle that flattens the body of the wave while leaving the overshoots intact at the jumps.

## The Fourier transform

Releasing the period to infinity converts the Fourier series into the **Fourier transform**. For a function $f$ on the real line, define

$$
\tilde f(k) \;=\; \int_{-\infty}^{\infty} f(x)\, e^{-i k x}\, \mathrm{d} x, \tag{0.4.10}
$$

with inverse

$$
f(x) \;=\; \frac{1}{2\pi} \int_{-\infty}^{\infty} \tilde f(k)\, e^{i k x}\, \mathrm{d} k. \tag{0.4.11}
$$

(Several sign and prefactor conventions exist; we follow the physics convention.)

The Fourier transform decomposes $f$ into plane waves of wavenumber $k$, with $\tilde f(k)$ the complex amplitude of the $k$-component. The variables $x$ and $k$ are **conjugate**: position and wavenumber, or in quantum mechanics, position and momentum (with $p = \hbar k$).

### Worked example: the Gaussian is its own Fourier transform

The Gaussian is uniquely placed in Fourier analysis as a fixed point (up to scaling) of the transform. Let us verify and quantify this with a calculation.

Take a real-space Gaussian centred at zero with width $\sigma_x$:
$$
f(x) = \frac{1}{\sigma_x \sqrt{2\pi}}\, e^{-x^2 / (2 \sigma_x^2)}.
$$

**Step (1).** Substitute into (0.4.10):
$$
\tilde f(k) = \int_{-\infty}^{\infty} \frac{1}{\sigma_x \sqrt{2\pi}}\, e^{-x^2 / (2\sigma_x^2)}\, e^{-i k x}\, \mathrm{d} x.
$$

**Step (2).** Complete the square in the exponent. Write
$$
-\frac{x^2}{2\sigma_x^2} - i k x = -\frac{1}{2\sigma_x^2}\left( x^2 + 2 i k \sigma_x^2 x \right) = -\frac{1}{2\sigma_x^2}\big( x + i k \sigma_x^2 \big)^2 - \frac{k^2 \sigma_x^2}{2}.
$$

!!! note "Why this step?"
    Completing the square is the standard tool for evaluating Gaussian integrals against complex exponentials. We pull the $k$-dependence out of the integral as a clean prefactor, leaving inside only a shifted Gaussian.

**Step (3).** Substitute and split:
$$
\tilde f(k) = e^{-k^2 \sigma_x^2 / 2} \cdot \frac{1}{\sigma_x \sqrt{2\pi}} \int_{-\infty}^{\infty} e^{-(x + i k \sigma_x^2)^2 / (2 \sigma_x^2)}\, \mathrm{d} x.
$$
By contour deformation (or by the fact that the integrand is analytic and decays at infinity), the shifted Gaussian integrates to the same value as the unshifted one, namely $\sigma_x \sqrt{2\pi}$. Hence
$$
\tilde f(k) = e^{-k^2 \sigma_x^2 / 2} = e^{-k^2 / (2 \sigma_k^2)} \quad \text{with} \quad \sigma_k = 1/\sigma_x.
$$

**Step (4).** Read off the duality. The Fourier transform of a Gaussian of width $\sigma_x$ is another Gaussian, with width $\sigma_k = 1/\sigma_x$. The product
$$
\sigma_x\, \sigma_k = 1
$$
is fixed. Narrow in real space ↔ broad in reciprocal space, and vice versa.

This is, up to the factor of $\hbar$ and a numerical prefactor, the **Heisenberg uncertainty relation**. A localised wavepacket of position spread $\sigma_x$ has momentum spread $\sigma_p = \hbar/\sigma_x$, and the product $\sigma_x \sigma_p = \hbar$ (a minimum, achieved by Gaussian wavepackets). The mathematical statement of the uncertainty relation is exactly the Fourier-duality statement above.

```python
import numpy as np
import matplotlib.pyplot as plt
# Numerical check via FFT
N = 2048
x = np.linspace(-20, 20, N, endpoint=False)
dx = x[1] - x[0]
sigma_x = 1.5
f = np.exp(-x**2 / (2*sigma_x**2)) / (sigma_x * np.sqrt(2*np.pi))
F = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(f))) * dx
k = np.fft.fftshift(np.fft.fftfreq(N, d=dx)) * 2 * np.pi
sigma_k = 1.0 / sigma_x
analytic = np.exp(-k**2 * sigma_x**2 / 2)
print(np.max(np.abs(F.real - analytic)))   # tiny
```

### The key intuition

The single most useful sentence in this section is:

> The Fourier transform exchanges real-space localisation for reciprocal-space localisation.

A sharp Gaussian peak in $x$-space has a broad Gaussian profile in $k$-space, and vice versa. A perfectly periodic function in $x$ (period $L$) has a Fourier transform supported only on the discrete set $k = 2\pi n / L$. A delta-function in $x$ has a flat, uniform Fourier transform; a flat function in $x$ has a delta-function transform. The two pictures are mathematically equivalent and physically complementary.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. State Euler's formula and sketch one of the two derivations given (Taylor series or the defining ODE).
    2. The Fourier transform of a Gaussian of real-space width $\sigma_x$ is another Gaussian — what is its width, and what fixed product do the two widths satisfy?
    3. Why is this Gaussian width relation the mathematical content of the Heisenberg uncertainty principle, and what does it imply about localising a wavepacket in both position and momentum?

    If any of these is shaky, re-read the preceding section before continuing.

This duality is **everywhere** in solid-state physics. Three illustrations:

- **Real space ↔ reciprocal space.** A crystal lattice $\{ \mathbf{R} \}$ has a dual reciprocal lattice $\{ \mathbf{G} \}$ with $e^{i \mathbf{G} \cdot \mathbf{R}} = 1$ for all lattice translations. Functions that are periodic in real space have Fourier components only at reciprocal-lattice vectors. The Brillouin zone, band structure, X-ray diffraction patterns, and plane-wave DFT bases all live in reciprocal space.
- **Position ↔ momentum.** In quantum mechanics the position-space wavefunction $\psi(\mathbf{r})$ and the momentum-space wavefunction $\tilde\psi(\mathbf{p})$ are Fourier partners. The Heisenberg uncertainty relation $\Delta x\, \Delta p \ge \hbar/2$ is, mathematically, a statement about the Fourier-transform pair: you cannot localise a function arbitrarily well in both spaces simultaneously.
- **Time ↔ frequency.** A signal $f(t)$ and its spectrum $\tilde f(\omega)$ are Fourier partners. A short pulse in time has a broad spectrum; a long monochromatic wave has a narrow spectrum. This is the operating principle behind every spectroscopy you will use.

## The discrete Fourier transform

In numerics we work with sampled data $f_0, f_1, \ldots, f_{N-1}$ at $N$ equally spaced points. The **discrete Fourier transform** (DFT) is

$$
\tilde f_m \;=\; \sum_{n=0}^{N-1} f_n\, e^{-i\, 2\pi m n / N}, \qquad m = 0, 1, \ldots, N-1, \tag{0.4.12}
$$

with the inverse

$$
f_n = \frac{1}{N} \sum_{m=0}^{N-1} \tilde f_m\, e^{i\, 2\pi m n / N}. \tag{0.4.13}
$$

A naive evaluation of (0.4.12) costs $O(N^2)$ operations, but the **Fast Fourier Transform** (FFT) algorithm reduces this to $O(N \log N)$, making large transforms tractable. NumPy's `np.fft.fft` is a high-quality implementation.

### How the FFT works (in two paragraphs)

The Cooley–Tukey FFT exploits a recursive factorisation. For $N$ even, split (0.4.12) into sums over even and odd indices:
$$
\tilde f_m = \sum_{n=0}^{N/2 - 1} f_{2n}\, e^{-i 2\pi m (2n)/N} + \sum_{n=0}^{N/2 - 1} f_{2n+1}\, e^{-i 2\pi m (2n+1)/N}.
$$
The first sum is a DFT of length $N/2$ on the even-indexed data; the second is a similar DFT on the odd-indexed data, multiplied by a *twiddle factor* $e^{-i 2\pi m / N}$. Recursing — splitting each $N/2$ transform into two $N/4$ transforms, and so on — produces $\log_2 N$ levels, each doing $O(N)$ work; the total cost is $O(N \log N)$.

The recursive splitting reorders the input array. After all the bisections, the data appears in a permutation called **bit-reversal**: index $n$, written in binary, ends up at the position obtained by reversing those bits. For $N = 8$, the binary indices $000, 001, 010, 011, 100, 101, 110, 111$ map under bit-reversal to $000, 100, 010, 110, 001, 101, 011, 111$, that is, $0, 4, 2, 6, 1, 5, 3, 7$. Picture the array being shuffled like a perfect riffle of a deck of cards: the even indices come first in their bit-reversed order, then the odd ones. Modern FFT libraries hide this permutation, but you may glimpse it if you ever step through an FFT in a debugger and wonder why the intermediate state looks scrambled. The pattern is essential: it is exactly the order in which the recursive even/odd splits read the data, and any in-place FFT must perform this permutation either at the start (decimation-in-time) or at the end (decimation-in-frequency).

### A worked example: FFT of a sinusoid

To make all of this concrete, let us sample a sinusoid, take its FFT, and confirm that the spectrum is what we expect.

```python
import numpy as np
import matplotlib.pyplot as plt

# Sampling parameters.
N: int = 512                       # number of samples
T: float = 1.0                     # total duration in seconds
dt: float = T / N
t: np.ndarray = np.arange(N) * dt

# Signal: a 50 Hz sinusoid plus a 120 Hz one.
f1: float = 50.0
f2: float = 120.0
signal: np.ndarray = np.sin(2 * np.pi * f1 * t) + 0.5 * np.sin(2 * np.pi * f2 * t)

# Discrete Fourier transform.
S: np.ndarray = np.fft.fft(signal)
freqs: np.ndarray = np.fft.fftfreq(N, d=dt)

# Plot magnitude spectrum on the positive-frequency half.
mask = freqs >= 0
fig, ax = plt.subplots(2, 1, figsize=(7, 5))
ax[0].plot(t, signal)
ax[0].set_xlabel("time / s")
ax[0].set_ylabel("signal")
ax[1].stem(freqs[mask], np.abs(S[mask]) / N * 2, basefmt=" ")
ax[1].set_xlim(0, 200)
ax[1].set_xlabel("frequency / Hz")
ax[1].set_ylabel("|amplitude|")
plt.tight_layout()
plt.show()
```

The amplitude spectrum should show two sharp peaks: one at 50 Hz with height 1, one at 120 Hz with height 0.5. Those are precisely the amplitudes of the two sinusoidal components we put in. The FFT has decomposed the time-domain signal into its frequency-domain content — exactly as (0.4.10)–(0.4.11) promise, just on a discrete grid.

!!! note "Aliasing"
    A sinusoid of frequency $f$ can only be reconstructed unambiguously if the sampling rate exceeds $2f$ — the **Nyquist criterion**. Higher-frequency components fold back to lower apparent frequencies, an effect called aliasing. In practice this means choosing a fine enough real-space grid for the highest reciprocal-space components in your problem; in DFT calculations this is what the plane-wave cutoff controls.

## The delta function

A device of enormous utility — and questionable rigour from a pure-mathematician's standpoint — is the **Dirac delta function** $\delta(x)$. It is defined operationally by the sifting property:
$$
\int_{-\infty}^{\infty} f(x)\, \delta(x - x_0)\, \mathrm{d} x = f(x_0)
$$
for any continuous $f$. Intuitively, $\delta(x)$ is a "spike" of infinite height and zero width at $x = 0$, with total integral one. Strictly speaking it is a *distribution* (a generalised function), but for our purposes the sifting property is enough.

Two key Fourier identities involving $\delta$:
$$
\delta(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} e^{i k x}\, \mathrm{d} k, \qquad \tilde\delta(k) = 1.
$$
The first reads: a delta function is the equal superposition of all plane waves. The second: the Fourier transform of a delta is a constant. These are the limiting cases of the Gaussian duality we worked out above — take $\sigma_x \to 0$, and the real-space Gaussian becomes a delta while its reciprocal-space counterpart becomes the constant function.

In crystallography (Chapter 3 §5), the Fourier transform of an infinite crystal lattice $\sum_\mathbf{R} \delta(\mathbf{r} - \mathbf{R})$ is a sum of deltas at reciprocal-lattice vectors $\sum_\mathbf{G} \delta(\mathbf{k} - \mathbf{G})$. This is what produces the sharp Bragg peaks of an X-ray diffraction pattern, and what makes the reciprocal lattice meaningful as a Fourier-space description of the crystal.

## A few useful properties

Three properties of the Fourier transform are used so often they are worth memorising.

**Linearity.** $\widetilde{\alpha f + \beta g} = \alpha \tilde f + \beta \tilde g$.

**Shift theorem.** A real-space translation becomes a phase factor in reciprocal space:
$$
\widetilde{f(x - x_0)}(k) = e^{-i k x_0} \tilde f(k). \tag{0.4.14}
$$

**Convolution theorem.** Define the convolution $(f \ast g)(x) = \int f(y)\, g(x - y)\, \mathrm{d} y$. Then
$$
\widetilde{f \ast g} = \tilde f \cdot \tilde g. \tag{0.4.15}
$$

The convolution theorem is the engine behind fast electrostatic solvers and the smooth-particle-mesh Ewald method used in classical MD: convolutions in real space cost $O(N^2)$ to evaluate directly, but $O(N \log N)$ via two FFTs, a multiplication, and an inverse FFT.

### Sanity check: convolution as smoothing

A picture for the convolution theorem: convolving a signal with a wide Gaussian *smooths* the signal. In real space, smoothing is an $O(N \times W)$ operation, with $W$ the kernel width. In reciprocal space, smoothing is *multiplication* by a function (the Fourier transform of the kernel) that is concentrated near $k = 0$ — i.e. it kills high-$k$ components, suppressing fine structure. The two descriptions are identical, but the second is computationally trivial.

This dichotomy — local in real space ↔ pointwise in reciprocal space — is omnipresent in solid-state codes. A kinetic-energy operator $-\hbar^2 \nabla^2 / 2m$ is a second-derivative in real space (costly to apply on a real-space grid) but a simple multiplication by $\hbar^2 k^2 / 2m$ in reciprocal space (essentially free). A local potential $V(\mathbf{r})$ is the opposite: pointwise multiplication in real space, but a convolution in reciprocal space. The plane-wave DFT code design philosophy is to keep wavefunctions in whichever representation makes the current operator cheap, and to FFT between them as needed.

### Parseval's theorem

A final identity that ties the duality to energy conservation: **Parseval's theorem** states
$$
\int_{-\infty}^{\infty} |f(x)|^2\, \mathrm{d} x = \frac{1}{2\pi} \int_{-\infty}^{\infty} |\tilde f(k)|^2\, \mathrm{d} k. \tag{0.4.16}
$$
The total "power" of the signal is the same in real space and reciprocal space. This guarantees that a wavefunction normalised in real space is also normalised, up to the conventional prefactor, in momentum space. We will use Parseval routinely when comparing real-space and reciprocal-space sums in Chapter 5.

## The Fourier transform in higher dimensions

The one-dimensional transform (0.4.10) extends to $\mathbb{R}^d$ in the obvious way:
$$
\tilde f(\mathbf{k}) = \int_{\mathbb{R}^d} f(\mathbf{r})\, e^{-i \mathbf{k} \cdot \mathbf{r}}\, \mathrm{d}^d \mathbf{r},
$$
with inverse $f(\mathbf{r}) = (2\pi)^{-d} \int \tilde f(\mathbf{k})\, e^{i \mathbf{k} \cdot \mathbf{r}}\, \mathrm{d}^d \mathbf{k}$. Each Cartesian direction is transformed independently, so a separable function $f(\mathbf{r}) = f_x(x) f_y(y) f_z(z)$ has separable transform $\tilde f_x(k_x) \tilde f_y(k_y) \tilde f_z(k_z)$.

A particularly elegant special case: the Fourier transform of a spherically symmetric function $f(r)$ is itself spherically symmetric and reduces to a one-dimensional integral known as the **Hankel transform**:
$$
\tilde f(k) = \frac{4\pi}{k} \int_0^\infty r\, f(r)\, \sin(k r)\, \mathrm{d} r.
$$
This is the calculation behind the X-ray scattering amplitude of an isotropic atom: take the spherical electron density $\rho(r)$, perform this integral, and out comes the atomic form factor $f(k)$ — Chapter 3 §6.

## Where this is used

- Chapter 3.5 (solid-state physics) introduces Bloch's theorem; the resulting plane-wave decomposition is a Fourier series of the periodic part of the wavefunction.
- Chapter 5 (DFT) and Chapter 6 (running plane-wave codes) live entirely in reciprocal space.
- Chapter 8 (statistical mechanics) uses the FFT for structure factors and pair correlation functions.

Equipped with both the real and complex sides of analysis, we move to the final mathematical pillar: probability.
