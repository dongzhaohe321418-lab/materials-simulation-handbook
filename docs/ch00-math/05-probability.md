# 0.5 Probability and Statistics

Section 0.4 closed our review of deterministic analysis. The remaining mathematical ingredient is probability — the language of uncertainty, fluctuation, and ensemble averaging. Probability enters computational materials science in two big ways. Internally, statistical mechanics describes thermal systems through probability distributions over microstates; the Boltzmann factor we met in Section 0.1 is its emblem. Externally, machine learning models are built and assessed in probabilistic terms; loss functions, uncertainty quantification, and Bayesian optimisation (Chapter 11) all require fluency with random variables.

## Random variables

Informally, a **random variable** $X$ is a quantity whose value is not known in advance but is drawn from some underlying distribution. We distinguish two species.

A **discrete** random variable takes values in a countable set $\{x_1, x_2, \ldots\}$. Its **probability mass function** $p$ assigns to each outcome a non-negative number with

$$
p(x_i) \ge 0, \qquad \sum_i p(x_i) = 1. \tag{0.5.1}
$$

Examples: the spin state of an Ising variable ($\pm 1$); the number of vacancies in a finite crystal; the predicted class label of a structure (metal vs. insulator).

A **continuous** random variable takes values in $\mathbb{R}$ (or some subset). It is described by a **probability density function** $\rho(x)$ with

$$
\rho(x) \ge 0, \qquad \int_{-\infty}^{\infty} \rho(x)\, \mathrm{d} x = 1. \tag{0.5.2}
$$

The probability that $X$ falls in an interval $[a, b]$ is

$$
P(a \le X \le b) = \int_a^b \rho(x)\, \mathrm{d} x. \tag{0.5.3}
$$

Examples: a thermally-fluctuating bond length; a force component on an atom in a Langevin thermostat; the predicted formation energy from a regression model.

### Density is not probability

A point of frequent confusion. For a continuous variable, $\rho(x)$ is a probability *per unit $x$*, not a probability. The probability of an exact value, $P(X = x)$, is zero; only the probability of being in an interval is meaningful. Concretely, $\rho(x)$ can exceed $1$ — that simply means the density is concentrated narrowly. What must integrate to one is $\rho(x)\, \mathrm{d} x$ over the entire domain.

This distinction is the source of countless bugs in numerical sampling code. If you discretise a continuous distribution onto a grid of spacing $\Delta x$, the *probability* of bin $i$ is approximately $\rho(x_i)\, \Delta x$, not $\rho(x_i)$.

## Expectation, mean, variance

The **expectation** of a random variable $X$ with density $\rho$ is

$$
\langle X \rangle \;\equiv\; \mathbb{E}[X] \;=\; \int x\, \rho(x)\, \mathrm{d} x, \tag{0.5.4}
$$

or, in the discrete case, $\sum_i x_i\, p(x_i)$. We write the expectation with angle brackets in physical contexts and as $\mathbb{E}$ in mathematical or ML contexts; both mean the same thing.

More generally, the expectation of any function $g$ of $X$ is

$$
\langle g(X) \rangle = \int g(x)\, \rho(x)\, \mathrm{d} x. \tag{0.5.5}
$$

Expectation is linear: $\langle aX + bY \rangle = a\langle X \rangle + b \langle Y \rangle$ for any random variables $X, Y$ and constants $a, b$.

The **variance** measures the spread of $X$ about its mean:

$$
\mathrm{Var}(X) \;=\; \langle (X - \langle X \rangle)^2 \rangle \;=\; \langle X^2 \rangle - \langle X \rangle^2. \tag{0.5.6}
$$

Its square root, the **standard deviation** $\sigma_X = \sqrt{\mathrm{Var}(X)}$, has the same units as $X$ and is the more interpretable summary of spread.

Two useful identities. If $X, Y$ are **independent** (their joint density factorises),

$$
\mathrm{Var}(X + Y) = \mathrm{Var}(X) + \mathrm{Var}(Y). \tag{0.5.7}
$$

For any constant $a$,

$$
\mathrm{Var}(aX) = a^2 \mathrm{Var}(X). \tag{0.5.8}
$$

```python
import numpy as np
rng = np.random.default_rng(seed=0)
x = rng.normal(loc=3.0, scale=2.0, size=100_000)
print(np.mean(x), np.var(x), np.std(x))  # ≈ 3.0, 4.0, 2.0
```

## The Gaussian (normal) distribution

The single most important continuous distribution is the **Gaussian**, with density

$$
\rho(x) \;=\; \frac{1}{\sqrt{2\pi}\, \sigma}\, \exp\!\left(-\frac{(x - \mu)^2}{2 \sigma^2}\right), \tag{0.5.9}
$$

parameterised by mean $\mu$ and standard deviation $\sigma$. The prefactor is the unique normalisation that makes $\int \rho \, \mathrm{d} x = 1$.

Three properties make the Gaussian central.

**1. Closure under linear combinations.** Any linear combination of independent Gaussians is Gaussian. If $X_i \sim \mathcal{N}(\mu_i, \sigma_i^2)$ independently, then $\sum_i a_i X_i \sim \mathcal{N}(\sum_i a_i \mu_i,\, \sum_i a_i^2 \sigma_i^2)$ by (0.5.7)–(0.5.8) and a small additional argument.

**2. Maximum entropy.** Among all distributions with a given mean and variance, the Gaussian maximises the entropy $-\int \rho \ln \rho \, \mathrm{d} x$. In the absence of additional information, the Gaussian is the "least biased" choice — the formal version of "innocent until proven guilty".

**3. The central limit theorem.** The sum of many independent random variables, suitably normalised, tends to a Gaussian regardless of the individual distributions. Precisely: if $X_1, \ldots, X_N$ are independent and identically distributed with finite mean $\mu$ and variance $\sigma^2$, then

$$
\frac{1}{\sqrt N} \sum_{i=1}^{N} \frac{X_i - \mu}{\sigma} \;\xrightarrow{d}\; \mathcal{N}(0, 1) \quad \text{as } N \to \infty. \tag{0.5.10}
$$

This theorem is why Gaussians appear in places that have nothing to do with Gaussian inputs. Measurement errors, thermal noise, financial returns, sums of many small independent contributions — they all tend toward normality. The Maxwell–Boltzmann velocity distribution of a classical ideal gas, which you will derive in Chapter 8, is exactly a Gaussian in each Cartesian component because each component is a sum of many small molecular collisions.

```python
import numpy as np
rng = np.random.default_rng(0)
# Sum many uniform[-1, 1] samples; the result is approximately Gaussian.
M, N = 50_000, 12
sums = rng.uniform(-1.0, 1.0, size=(M, N)).sum(axis=1)
print(np.mean(sums), np.std(sums))  # ≈ 0, sqrt(N/3) ≈ 2.0
```

## The Boltzmann distribution

The Boltzmann factor of Section 0.1 generalises to a full probability distribution over the microstates of a thermal system. If $\{s\}$ enumerates the microstates and $E(s)$ is the energy of state $s$, the probability of finding the system in state $s$ at temperature $T$ is

$$
p(s) \;=\; \frac{1}{Z}\, e^{-\beta E(s)}, \qquad \beta = \frac{1}{k_\mathrm{B} T}, \tag{0.5.11}
$$

with the **partition function** $Z = \sum_s e^{-\beta E(s)}$ ensuring $\sum_s p(s) = 1$. For continuous coordinates the sum becomes an integral over phase space, and $p(s)$ is interpreted as a density (subject to the cautions of Section 0.5 on density versus probability).

The Boltzmann distribution is the prototype example of a physically motivated probability distribution; it is **the** distribution of equilibrium statistical mechanics, and almost every observable in Chapter 8 is an expectation against it. Two of its features deserve attention now.

First, the partition function $Z$ encodes thermodynamics. The Helmholtz free energy is $F = -k_\mathrm{B} T \ln Z$. The mean energy is

$$
\langle E \rangle = -\frac{\partial \ln Z}{\partial \beta}. \tag{0.5.12}
$$

The heat capacity is the variance of the energy divided by $k_\mathrm{B} T^2$,

$$
C_V = \frac{1}{k_\mathrm{B} T^2}\, \mathrm{Var}(E). \tag{0.5.13}
$$

So thermodynamic response functions are statistical moments of the underlying microstate distribution. This is one of the most beautiful results in physics.

Second, Monte Carlo simulation works by drawing samples from (0.5.11) without ever computing $Z$. The Metropolis–Hastings algorithm, which you will implement in Chapter 8, only requires ratios $e^{-\beta \Delta E}$, sidestepping the typically intractable partition function entirely.

## Joint and conditional probability

For two random variables $X, Y$ with joint density $\rho(x, y)$, the **marginal** density of $X$ is

$$
\rho_X(x) = \int \rho(x, y)\, \mathrm{d} y, \tag{0.5.14}
$$

and the **conditional** density of $Y$ given $X = x$ is

$$
\rho(y \mid x) = \frac{\rho(x, y)}{\rho_X(x)}, \quad \rho_X(x) > 0. \tag{0.5.15}
$$

Two variables are **independent** if and only if $\rho(x, y) = \rho_X(x)\,\rho_Y(y)$, equivalently $\rho(y \mid x) = \rho_Y(y)$ — knowing $X$ tells you nothing about $Y$.

## Bayes' rule

From the definition (0.5.15) and its symmetric counterpart, $\rho(x, y) = \rho(y \mid x)\, \rho_X(x) = \rho(x \mid y)\, \rho_Y(y)$. Solving for $\rho(y \mid x)$:

$$
\rho(y \mid x) \;=\; \frac{\rho(x \mid y)\, \rho_Y(y)}{\rho_X(x)}. \tag{0.5.16}
$$

This is **Bayes' rule**. In the language of inference, with $\theta$ for parameters and $D$ for observed data,

$$
\underbrace{p(\theta \mid D)}_{\text{posterior}} \;=\; \frac{\overbrace{p(D \mid \theta)}^{\text{likelihood}}\; \overbrace{p(\theta)}^{\text{prior}}}{\underbrace{p(D)}_{\text{evidence}}}. \tag{0.5.17}
$$

In words: the posterior belief about parameters after seeing the data is proportional to the likelihood of the data under those parameters, times the prior belief.

Three remarks for later use.

**1.** The evidence $p(D)$ is just the normalising constant $\int p(D \mid \theta)\, p(\theta) \, \mathrm{d}\theta$. For point estimation we rarely need it explicitly — we just maximise the numerator. For model comparison it becomes important.

**2.** Bayesian updating is iterative. Yesterday's posterior is today's prior. After observing $D_1$ and $D_2$ in sequence, with the two observations independent given $\theta$, we have $p(\theta \mid D_1, D_2) \propto p(D_2 \mid \theta)\, p(D_1 \mid \theta)\, p(\theta)$. This is the mathematical engine of active learning: collect a measurement, update beliefs, choose the next measurement.

**3.** Gaussian processes — the workhorse surrogate model in Chapter 11's Bayesian optimisation — are Bayes' rule applied to function-valued random variables with Gaussian priors over functions. The posterior remains Gaussian, and acquisition functions like expected improvement are tractable integrals against this posterior. We will spell this out properly when the time comes; for now the formula (0.5.17) is the seed.

## A short numerical illustration: Bayesian coin flip

Suppose we are uncertain whether a coin is fair. Let $\theta \in [0, 1]$ be the unknown probability of heads, with a uniform prior $p(\theta) = 1$. We observe $n_H$ heads in $N$ flips. The likelihood is

$$
p(D \mid \theta) = \theta^{n_H} (1 - \theta)^{N - n_H},
$$

and the posterior is (up to a normalisation constant)

$$
p(\theta \mid D) \propto \theta^{n_H} (1 - \theta)^{N - n_H},
$$

a Beta distribution that peaks at $\theta = n_H / N$. As $N$ grows, the posterior narrows around the empirical frequency.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import beta

theta = np.linspace(0, 1, 400)
for nH, N in [(2, 5), (8, 20), (40, 100)]:
    plt.plot(theta, beta.pdf(theta, nH + 1, N - nH + 1),
             label=f"{nH}/{N} heads")
plt.xlabel(r"$\theta$ (probability of heads)")
plt.ylabel("posterior density")
plt.legend()
plt.show()
```

Three observations: more data sharpens the posterior; with uniform prior and Bernoulli likelihood, the posterior mode is the maximum-likelihood estimate; the spread of the posterior is the formal Bayesian quantification of remaining uncertainty.

!!! note "Frequentist versus Bayesian"
    Section 0.5 has implicitly adopted a Bayesian outlook — probabilities as degrees of belief, updated by data. The frequentist counterpart — probabilities as long-run frequencies, parameters as fixed unknowns — is equally valid and produces, in many cases, numerically similar answers. Modern materials ML mixes both: maximum-likelihood training of neural networks is frequentist in spirit; Bayesian optimisation and uncertainty quantification are Bayesian. You will need to be bilingual.

## Sampling and Monte Carlo, briefly

A great deal of computational materials science is, at heart, the problem of computing expectations like (0.5.5) when the integral is too high-dimensional for grid-based quadrature. The solution is **Monte Carlo**: draw samples $x_1, \ldots, x_M$ from $\rho$ and approximate

$$
\langle g(X) \rangle \;\approx\; \frac{1}{M} \sum_{i=1}^{M} g(x_i), \tag{0.5.18}
$$

with statistical error scaling as $1/\sqrt M$ by the central limit theorem (0.5.10). The error is independent of dimension — the Monte Carlo killer feature that makes statistical mechanics in $3N$ dimensions tractable.

When direct sampling from $\rho$ is impossible — as it is for the Boltzmann distribution at non-trivial energies — Markov-chain Monte Carlo methods construct a stochastic process whose stationary distribution is the target. This is the subject of Chapter 8.

## Where this is used

- Chapter 7 (MD) uses the Maxwell–Boltzmann distribution (a Gaussian) to assign initial velocities and a Langevin equation (Gaussian noise) for thermostatting.
- Chapter 8 (statistical mechanics) is built entirely around the Boltzmann distribution (0.5.11) and Monte Carlo sampling (0.5.18).
- Chapter 9 and Chapter 10 train ML models by maximising likelihoods or minimising loss functions that are negative log-likelihoods.
- Chapter 11 (active learning, Bayesian optimisation) uses Bayes' rule (0.5.17) to update Gaussian-process surrogates after each new measurement.

This concludes the mathematical groundwork. The eight exercises that follow will give you a chance to test the toolkit before we move on to the practical side: getting Python, NumPy, and the rest of the scientific stack working on your machine in Chapter 1.
