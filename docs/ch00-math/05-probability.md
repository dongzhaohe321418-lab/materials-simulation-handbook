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

### Covariance and correlation

When two random variables interact, their joint statistics are summarised by the **covariance**
$$
\mathrm{Cov}(X, Y) = \langle (X - \langle X \rangle)(Y - \langle Y \rangle) \rangle = \langle XY \rangle - \langle X \rangle \langle Y \rangle. \tag{0.5.8a}
$$
A positive covariance means $X$ and $Y$ tend to deviate from their means in the same direction; negative means opposite directions. The dimensionless analogue is the **Pearson correlation coefficient**
$$
\rho_{XY} = \frac{\mathrm{Cov}(X, Y)}{\sigma_X \sigma_Y} \in [-1, 1]. \tag{0.5.8b}
$$
Independence implies zero covariance, but zero covariance does *not* imply independence — a frequent source of confusion. Two variables can be uncorrelated yet strongly dependent through non-linear interactions. The variance of a sum generalises (0.5.7):
$$
\mathrm{Var}(X + Y) = \mathrm{Var}(X) + \mathrm{Var}(Y) + 2\, \mathrm{Cov}(X, Y).
$$
The covariance matrix $\Sigma_{ij} = \mathrm{Cov}(X_i, X_j)$ of a vector-valued random variable plays the same role for multivariate Gaussians that $\sigma^2$ plays for univariate ones. It is symmetric and positive semi-definite — once again, Section 0.2 linear algebra is the right toolkit.

## The Gaussian (normal) distribution

The single most important continuous distribution is the **Gaussian**, with density

$$
\rho(x) \;=\; \frac{1}{\sqrt{2\pi}\, \sigma}\, \exp\!\left(-\frac{(x - \mu)^2}{2 \sigma^2}\right), \tag{0.5.9}
$$

parameterised by mean $\mu$ and standard deviation $\sigma$. The prefactor is the unique normalisation that makes $\int \rho \, \mathrm{d} x = 1$.

Three properties make the Gaussian central.

**1. Closure under linear combinations.** Any linear combination of independent Gaussians is Gaussian. If $X_i \sim \mathcal{N}(\mu_i, \sigma_i^2)$ independently, then $\sum_i a_i X_i \sim \mathcal{N}(\sum_i a_i \mu_i,\, \sum_i a_i^2 \sigma_i^2)$ by (0.5.7)–(0.5.8) and a small additional argument.

**2. Maximum entropy.** Among all distributions with a given mean and variance, the Gaussian maximises the entropy $-\int \rho \ln \rho \, \mathrm{d} x$. In the absence of additional information, the Gaussian is the "least biased" choice — the formal version of "innocent until proven guilty".

### Deriving the Gaussian as a maximum-entropy distribution

The maximum-entropy claim deserves a full derivation; it is the cleanest way to motivate the Gaussian's appearance everywhere from thermal noise to neural-network priors.

The problem: maximise the entropy functional
$$
S[\rho] = -\int_{-\infty}^{\infty} \rho(x) \ln \rho(x)\, \mathrm{d} x
$$
subject to three constraints: normalisation $\int \rho\, \mathrm{d} x = 1$; fixed mean $\int x\,\rho\, \mathrm{d} x = \mu$; and fixed variance $\int (x - \mu)^2 \rho\, \mathrm{d} x = \sigma^2$.

**Step (1).** Form the Lagrangian with three multipliers $\alpha, \beta, \gamma$:
$$
L[\rho] = -\int \rho \ln \rho\, \mathrm{d} x - \alpha \left( \int \rho\, \mathrm{d} x - 1 \right) - \beta \left( \int x \rho\, \mathrm{d} x - \mu \right) - \gamma \left( \int (x-\mu)^2 \rho \, \mathrm{d} x - \sigma^2 \right).
$$

**Step (2).** Functional derivative with respect to $\rho(x)$, set to zero. The derivative of $-\rho \ln \rho$ at fixed $x$ is $-\ln \rho - 1$. The three constraints contribute $-\alpha$, $-\beta x$, $-\gamma (x - \mu)^2$ respectively. Setting the total to zero:
$$
-\ln \rho(x) - 1 - \alpha - \beta x - \gamma (x - \mu)^2 = 0.
$$

!!! note "Why this step?"
    The functional derivative of $-\int \rho \ln \rho\, \mathrm{d} x$ with respect to $\rho(x)$ is most easily found by treating $\rho \ln \rho$ as an ordinary function of the variable $\rho$, then differentiating: $\mathrm{d}(\rho \ln \rho)/\mathrm{d}\rho = \ln \rho + 1$. The functional formalism turns each constraint into an extra additive term that is linear in $\rho$.

**Step (3).** Solve for $\rho$:
$$
\rho(x) = \exp\!\Big( -1 - \alpha - \beta x - \gamma (x - \mu)^2 \Big).
$$
For the integral over $\mathbb{R}$ to converge we need $\gamma > 0$. The linear $\beta x$ term can be absorbed into the quadratic by completing the square (and by the mean constraint, $\beta = 0$ when the centre is already at $\mu$). What remains is
$$
\rho(x) = C\, e^{-\gamma (x - \mu)^2}.
$$

**Step (4).** Fix $C$ and $\gamma$ from the normalisation and variance constraints. Normalisation gives $C = \sqrt{\gamma / \pi}$ (the standard Gaussian integral). The variance constraint gives $\langle (x-\mu)^2 \rangle = 1/(2\gamma) = \sigma^2$, so $\gamma = 1/(2 \sigma^2)$. Substituting yields
$$
\rho(x) = \frac{1}{\sqrt{2\pi}\, \sigma}\, e^{-(x-\mu)^2 / (2\sigma^2)},
$$
which is exactly (0.5.9). $\square$

This derivation is more than a curiosity: it is the proper justification for using the Gaussian as a prior in Bayesian inference, as a noise model in regression, and as the energy distribution in the harmonic-oscillator partition function. Whenever we know nothing about a quantity except its mean and variance, the maximum-entropy principle says to model it as a Gaussian.

**3. The central limit theorem.** The sum of many independent random variables, suitably normalised, tends to a Gaussian regardless of the individual distributions. Precisely: if $X_1, \ldots, X_N$ are independent and identically distributed with finite mean $\mu$ and variance $\sigma^2$, then

$$
\frac{1}{\sqrt N} \sum_{i=1}^{N} \frac{X_i - \mu}{\sigma} \;\xrightarrow{d}\; \mathcal{N}(0, 1) \quad \text{as } N \to \infty. \tag{0.5.10}
$$

This theorem is why Gaussians appear in places that have nothing to do with Gaussian inputs. Measurement errors, thermal noise, financial returns, sums of many small independent contributions — they all tend toward normality. The Maxwell–Boltzmann velocity distribution of a classical ideal gas, which you will derive in Chapter 8, is exactly a Gaussian in each Cartesian component because each component is a sum of many small molecular collisions.

### Sketch of the CLT proof via moment-generating functions

The classical proof uses **characteristic functions** (Fourier transforms of densities). Define the moment-generating function $M_X(t) = \mathbb{E}[e^{t X}]$, assuming it exists in some neighbourhood of $t = 0$.

**Step (1).** Let $Y_i = (X_i - \mu)/\sigma$ be the centred and scaled samples, with $\mathbb{E}[Y_i] = 0$ and $\mathrm{Var}(Y_i) = 1$. By Taylor expansion,
$$
M_{Y_i}(t) = 1 + t \cdot 0 + \frac{t^2}{2} \cdot 1 + O(t^3) = 1 + \frac{t^2}{2} + O(t^3).
$$

**Step (2).** Consider the standardised sum $S_N = N^{-1/2} \sum_i Y_i$. By independence, moment-generating functions of independent variables multiply:
$$
M_{S_N}(t) = \prod_{i=1}^{N} M_{Y_i}(t/\sqrt N) = \left[ 1 + \frac{t^2}{2 N} + O\!\left(\frac{t^3}{N^{3/2}}\right) \right]^N.
$$

**Step (3).** Take $N \to \infty$. Using $\lim_{N \to \infty} (1 + a/N)^N = e^a$,
$$
\lim_{N \to \infty} M_{S_N}(t) = e^{t^2 / 2}.
$$

!!! note "Why this step?"
    The function $e^{t^2 / 2}$ is the moment-generating function of the standard normal $\mathcal{N}(0, 1)$. A theorem of Lévy guarantees that pointwise convergence of moment-generating functions in a neighbourhood of zero implies convergence in distribution. So $S_N \xrightarrow{d} \mathcal{N}(0, 1)$. $\square$

The proof reveals **why** the Gaussian is the universal limit: it is the unique distribution whose moment-generating function is the exponential of a quadratic in $t$. Any sum of independent variables, once we expand to second order and exponentiate, lands here. The shape of the limiting distribution is determined entirely by the first two moments — mean and variance — which is exactly the maximum-entropy framing from earlier.

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

### Deriving the Boltzmann distribution from counting

Where does (0.5.11) come from? The deepest derivation rests on counting microstates in the microcanonical ensemble. We sketch it here; the details belong to Chapter 8.

**Setup.** Imagine a small system $A$ in thermal contact with a vast reservoir $R$, the combined system being isolated. The total energy $E_\mathrm{tot}$ is fixed. Let $\Omega_R(E_R)$ count the microstates of $R$ at energy $E_R$. When $A$ is in microstate $s$ with energy $E(s)$, the reservoir has energy $E_\mathrm{tot} - E(s)$, and so the number of accessible joint microstates is $\Omega_R(E_\mathrm{tot} - E(s))$.

**Step (1).** By the postulate of equal a priori probability (Chapter 8), every microstate of the isolated total system is equally likely. Hence
$$
p(s) \propto \Omega_R(E_\mathrm{tot} - E(s)).
$$

**Step (2).** Take the logarithm and Taylor-expand around $E_\mathrm{tot}$:
$$
\ln \Omega_R(E_\mathrm{tot} - E(s)) = \ln \Omega_R(E_\mathrm{tot}) - E(s) \cdot \frac{\partial \ln \Omega_R}{\partial E_R}\bigg|_{E_\mathrm{tot}} + O(E(s)^2).
$$
The first term is a constant, irrelevant to relative probabilities. The coefficient of $-E(s)$ in the second term is identified with $\beta = 1/(k_\mathrm{B} T)$, the **inverse temperature** of the reservoir. (This is essentially the definition of temperature: how rapidly the reservoir's log-multiplicity grows with energy.)

!!! note "Why this step?"
    The Taylor expansion of $\ln \Omega_R$ rather than $\Omega_R$ itself is essential. $\Omega_R$ is astronomically large — easily $10^{10^{23}}$ — so its direct Taylor series converges glacially. Its logarithm, by contrast, is well-behaved and order $N$. The reservoir's energy scale $E_\mathrm{tot}$ is enormous compared with $E(s)$, so the linear term dominates.

**Step (3).** Exponentiating,
$$
p(s) \propto e^{-\beta E(s)}.
$$
This is the Boltzmann factor. Normalising by $Z = \sum_s e^{-\beta E(s)}$ recovers (0.5.11). $\square$

### Stirling's approximation

A workhorse formula used throughout the derivation above (and many in Chapter 8) is **Stirling's approximation** for factorials of large numbers:
$$
\ln N! \approx N \ln N - N \quad \text{for large } N, \tag{0.5.11a}
$$
with relative error of order $\ln N / N$, hence vanishing fast for thermodynamic $N \sim 10^{23}$. A quick derivation: write $\ln N! = \sum_{k=1}^{N} \ln k$, then approximate the sum by the integral $\int_1^N \ln x\, \mathrm{d} x = N \ln N - N + 1$. The constant terms are sub-leading.

A sanity check: for $N = 10$, $\ln(10!) = \ln(3628800) \approx 15.10$, and $N \ln N - N = 10 \ln 10 - 10 \approx 13.03$. The relative error is about $14\%$, falling to $1\%$ at $N = 100$ and to $0.1\%$ at $N = 1000$. For the thermodynamic limit $N \to \infty$, Stirling is exact.

Stirling's approximation is what turns counting problems (like the multinomial number of ways to distribute $N$ particles into energy levels) into smooth optimisation problems amenable to Lagrange multipliers — and ultimately what produces the Boltzmann exponential.

## A short catalogue of distributions

Beyond the Gaussian, a few discrete and continuous distributions appear so often that they deserve naming.

The **Bernoulli** distribution: $X \in \{0, 1\}$ with $P(X = 1) = p$. The two-outcome trial. Mean $p$, variance $p(1-p)$.

The **binomial** distribution: $X \sim \mathrm{Bin}(N, p)$ counts the number of successes in $N$ independent Bernoulli trials, with $P(X = k) = \binom{N}{k} p^k (1-p)^{N-k}$. Mean $Np$, variance $Np(1-p)$. By the CLT, $\mathrm{Bin}(N, p)$ tends to a Gaussian for large $N$.

The **Poisson** distribution: $X \sim \mathrm{Pois}(\lambda)$ has $P(X = k) = e^{-\lambda} \lambda^k / k!$. It arises as the limit of $\mathrm{Bin}(N, p)$ with $N \to \infty$, $p \to 0$, $Np = \lambda$ fixed. Mean and variance both equal $\lambda$. Used for radioactive decay counts, photon counts, defect counts in a crystal.

The **exponential** distribution: $X \sim \mathrm{Exp}(\lambda)$ with density $\rho(x) = \lambda e^{-\lambda x}$ for $x \ge 0$. Mean $1/\lambda$, variance $1/\lambda^2$. It is the waiting-time distribution between events of a Poisson process — for instance, between successive jumps in a Markov-chain Monte Carlo simulation. It is also the *only* continuous distribution that is **memoryless**: $P(X > s + t \mid X > s) = P(X > t)$.

The **uniform** distribution on $[a, b]$: constant density $1/(b - a)$. The default "I know nothing more than the range" distribution; the maximum-entropy distribution under a support constraint.

The **multivariate Gaussian** $\mathcal{N}(\boldsymbol{\mu}, \Sigma)$: density $\rho(\mathbf{x}) \propto \exp\!\big( -\tfrac{1}{2}(\mathbf{x} - \boldsymbol{\mu})^\top \Sigma^{-1} (\mathbf{x} - \boldsymbol{\mu}) \big)$ with mean vector $\boldsymbol{\mu}$ and covariance matrix $\Sigma$. The generalisation of the univariate Gaussian to vector-valued random variables, central to Gaussian processes in Chapter 11.

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

### Worked example: Bayes' rule with a medical test

A useful intuition pump: a disease has prevalence $0.1\%$ in the general population, so $p(\text{disease}) = 0.001$. A diagnostic test has sensitivity $99\%$ — $p(\text{positive} \mid \text{disease}) = 0.99$ — and specificity $99\%$ — $p(\text{positive} \mid \text{no disease}) = 0.01$. You test positive. What is the probability you have the disease?

**Step (1).** Apply Bayes' rule (0.5.17). With $\theta = \text{disease}$ and $D = \text{positive test}$,
$$
p(\text{disease} \mid +) = \frac{p(+ \mid \text{disease})\, p(\text{disease})}{p(+)}.
$$

**Step (2).** Compute the evidence $p(+)$ by the law of total probability:
$$
p(+) = p(+ \mid \text{disease})\, p(\text{disease}) + p(+ \mid \text{no disease})\, p(\text{no disease}).
$$
Substituting:
$$
p(+) = 0.99 \cdot 0.001 + 0.01 \cdot 0.999 = 0.00099 + 0.00999 = 0.01098.
$$

**Step (3).** Form the ratio:
$$
p(\text{disease} \mid +) = \frac{0.99 \cdot 0.001}{0.01098} = \frac{0.00099}{0.01098} \approx 0.0902.
$$

!!! note "Why this step?"
    Despite a $99\%$ accurate test, a positive result implies only a $\sim 9\%$ posterior probability of disease. The reason is the low base rate: the false positives from the $99.9\%$ of healthy people swamp the true positives from the $0.1\%$ sick. This is the **base-rate fallacy**, and it appears whenever a rare phenomenon is screened for with an imperfect test.

The same arithmetic governs anomaly detection in materials informatics: if your model flags $1\%$ of structures as "unstable" but the true rate of instability is $0.1\%$, the precision of your alarm — the fraction of flagged structures that really are unstable — depends critically on the base rate, not just the model's accuracy.

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

## Information and entropy

The entropy functional we maximised in deriving the Gaussian deserves a section of its own. For a discrete distribution $p$ on outcomes $\{x_i\}$, the **Shannon entropy** is
$$
S(p) = -\sum_i p(x_i) \ln p(x_i).
$$
It is non-negative, with $S = 0$ iff $p$ is concentrated on a single outcome (no uncertainty), and maximal when $p$ is uniform (complete uncertainty among the available outcomes). For continuous distributions the analogue is the **differential entropy** $-\int \rho \ln \rho\, \mathrm{d} x$.

Two related quantities appear in machine learning. The **Kullback–Leibler divergence**
$$
D_\mathrm{KL}(p \,\|\, q) = \sum_i p(x_i) \ln \frac{p(x_i)}{q(x_i)}
$$
measures how distinguishable distribution $p$ is from a reference $q$. It is non-negative, zero iff $p = q$, and not symmetric. The **cross-entropy** $H(p, q) = -\sum_i p(x_i) \ln q(x_i)$ is the loss function for classification: training a neural network by maximum likelihood is mathematically the same as minimising cross-entropy.

The physical interpretation: $S$ is the average number of nats of information needed to specify a sample drawn from $p$. The Boltzmann entropy of statistical mechanics, $S_\mathrm{thermo} = k_\mathrm{B} \ln \Omega$, agrees with the Shannon entropy of the microcanonical distribution up to the prefactor $k_\mathrm{B}$. Chapter 8 develops this connection in earnest.

## Sampling and Monte Carlo, briefly

A great deal of computational materials science is, at heart, the problem of computing expectations like (0.5.5) when the integral is too high-dimensional for grid-based quadrature. The solution is **Monte Carlo**: draw samples $x_1, \ldots, x_M$ from $\rho$ and approximate

$$
\langle g(X) \rangle \;\approx\; \frac{1}{M} \sum_{i=1}^{M} g(x_i), \tag{0.5.18}
$$

with statistical error scaling as $1/\sqrt M$ by the central limit theorem (0.5.10). The error is independent of dimension — the Monte Carlo killer feature that makes statistical mechanics in $3N$ dimensions tractable.

When direct sampling from $\rho$ is impossible — as it is for the Boltzmann distribution at non-trivial energies — Markov-chain Monte Carlo methods construct a stochastic process whose stationary distribution is the target. This is the subject of Chapter 8.

## Markov chains and detailed balance

A **Markov chain** is a stochastic process on a state space $\mathcal{S}$ defined by transition probabilities $T(s' \mid s)$ — the probability of moving from $s$ to $s'$ in one step. The key property is **memorylessness**: the probability of the next state depends only on the current state, not on the history. Formally, if $X_t$ denotes the state at step $t$,
$$
P(X_{t+1} = s' \mid X_t = s, X_{t-1}, \ldots, X_0) = T(s' \mid s).
$$

If the chain has been running long enough, it converges (under technical conditions: irreducibility and aperiodicity) to a **stationary distribution** $\pi(s)$ satisfying
$$
\pi(s') = \sum_s T(s' \mid s)\, \pi(s). \tag{0.5.19}
$$
The stationary distribution is the eigenvector of the transition operator with eigenvalue $1$ — note the connection to Section 0.2 eigenproblems.

### Detailed balance

A sufficient condition for $\pi$ to be the stationary distribution is **detailed balance**:
$$
\pi(s)\, T(s' \mid s) = \pi(s')\, T(s \mid s'), \quad \text{for all } s, s'. \tag{0.5.20}
$$
This says that, in equilibrium, the flow of probability from $s$ to $s'$ exactly cancels the reverse flow. Summing over $s$ on both sides and using $\sum_{s'} T(s' \mid s) = 1$ recovers (0.5.19), so detailed balance implies stationarity.

The Metropolis–Hastings algorithm constructs a chain that explicitly satisfies detailed balance with $\pi(s) \propto e^{-\beta E(s)}$. The acceptance probability $\min(1, e^{-\beta \Delta E})$ is precisely what is needed to balance forward and reverse moves. We will derive this in Chapter 8 §3 and use it for Ising-model simulations.

Detailed balance has a beautiful physical interpretation: it is the microscopic statement of **time-reversal symmetry**. In a chain that obeys detailed balance, no movie of the trajectory could give away whether it was running forward or backward in time. This is exactly the property a thermal equilibrium state must have. Chains that violate detailed balance correspond to driven, non-equilibrium systems — also useful, but outside the scope of this book.

## Where this is used

- Chapter 7 (MD) uses the Maxwell–Boltzmann distribution (a Gaussian) to assign initial velocities and a Langevin equation (Gaussian noise) for thermostatting.
- Chapter 8 (statistical mechanics) is built entirely around the Boltzmann distribution (0.5.11) and Monte Carlo sampling (0.5.18).
- Chapter 9 and Chapter 10 train ML models by maximising likelihoods or minimising loss functions that are negative log-likelihoods.
- Chapter 11 (active learning, Bayesian optimisation) uses Bayes' rule (0.5.17) to update Gaussian-process surrogates after each new measurement.

This concludes the mathematical groundwork. The eight exercises that follow will give you a chance to test the toolkit before we move on to the practical side: getting Python, NumPy, and the rest of the scientific stack working on your machine in Chapter 1.
