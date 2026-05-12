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

### Derivation via Taylor series

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

A particularly celebrated consequence, obtained by setting $\theta = \pi$, is

$$
e^{i\pi} + 1 = 0,
$$

which combines five fundamental constants in a single equation.

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

## Why complex numbers are natural in quantum mechanics

The state of a quantum system is a **complex-valued** wavefunction $\psi(\mathbf{r}, t)$, evolving according to the time-dependent Schrödinger equation

$$
i \hbar \frac{\partial \psi}{\partial t} = \hat H \psi.
$$

The $i$ on the left-hand side is not optional. If $\psi$ were forced to be real, $\hat H \psi$ would have to be purely imaginary at every instant, conflicting with the requirement that $\hat H$ — the energy — is a real-valued observable. The complex phase is what stores the dynamical information; the squared modulus $|\psi|^2$ is what we measure.

For a stationary state of energy $E$, the time dependence is a pure complex exponential,

$$
\psi(\mathbf{r}, t) = \psi(\mathbf{r})\, e^{-i E t / \hbar}, \tag{0.4.8}
$$

and the probability density $|\psi(\mathbf{r}, t)|^2 = |\psi(\mathbf{r})|^2$ is independent of time. The complex phase rotates uniformly; the observable quantities do not. This stationary-phase picture, together with the Bloch theorem we will see in Chapter 3.5, is the bridge from quantum mechanics to band structure.

A second reason: **plane waves**. A free particle of momentum $\hbar \mathbf{k}$ has wavefunction $e^{i\mathbf{k}\cdot\mathbf{r}}$. Real sines and cosines cannot represent a wave travelling in a definite direction; they always contain a counter-propagating component (see (0.4.7)). Plane-wave bases — the natural basis for DFT in periodic solids — are intrinsically complex.

## Fourier series — a paragraph of motivation

A periodic function $f$ of period $L$ can be expanded as a sum of sines and cosines of frequencies that are integer multiples of $2\pi/L$. Writing this in complex form via (0.4.7),

$$
f(x) = \sum_{n=-\infty}^{\infty} c_n\, e^{i\, 2\pi n x / L}, \qquad c_n = \frac{1}{L} \int_0^L f(x)\, e^{-i\, 2\pi n x / L} \, \mathrm{d} x. \tag{0.4.9}
$$

The coefficients $c_n$ tell you "how much of each frequency" $f$ contains. In a crystal of lattice constant $L$, periodic functions — for instance the electron density — admit exactly this expansion, with $2\pi n / L$ the reciprocal-lattice wavevectors. This is the algebraic content of the reciprocal-lattice picture you will meet in Chapter 3.

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

### The key intuition

The single most useful sentence in this section is:

> The Fourier transform exchanges real-space localisation for reciprocal-space localisation.

A sharp Gaussian peak in $x$-space has a broad Gaussian profile in $k$-space, and vice versa. A perfectly periodic function in $x$ (period $L$) has a Fourier transform supported only on the discrete set $k = 2\pi n / L$. A delta-function in $x$ has a flat, uniform Fourier transform; a flat function in $x$ has a delta-function transform. The two pictures are mathematically equivalent and physically complementary.

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

## Where this is used

- Chapter 3.5 (solid-state physics) introduces Bloch's theorem; the resulting plane-wave decomposition is a Fourier series of the periodic part of the wavefunction.
- Chapter 5 (DFT) and Chapter 6 (running plane-wave codes) live entirely in reciprocal space.
- Chapter 8 (statistical mechanics) uses the FFT for structure factors and pair correlation functions.

Equipped with both the real and complex sides of analysis, we move to the final mathematical pillar: probability.
