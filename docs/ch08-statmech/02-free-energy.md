# 8.2 Free Energy Methods

<figure markdown>
```mermaid
flowchart TD
    Q{"What kind of<br/>free-energy<br/>difference?"}
    Q -->|"between two<br/>well-defined states<br/>(close in phase space)"| FEP["<b>FEP</b><br/>Free-Energy Perturbation<br/>ΔF = −kT ln ⟨e^{−βΔU}⟩"]
    Q -->|"smooth path<br/>between two states"| TI["<b>Thermodynamic<br/>Integration</b><br/>ΔF = ∫ ⟨∂U/∂λ⟩ dλ"]
    Q -->|"free energy along<br/>a known reaction<br/>coordinate (rare event)"| US["<b>Umbrella Sampling</b><br/>biased windows +<br/>WHAM"]
    Q -->|"high-dim or<br/>unknown CV /<br/>need on-the-fly bias"| MD["<b>Metadynamics</b><br/>history-dependent<br/>Gaussian bias"]
```
<figcaption>Decision tree for choosing a free-energy method. The root question asks what kind of free-energy difference is needed: between two well-defined states close in phase space points to free-energy perturbation; a smooth path between two states points to thermodynamic integration; a free energy along a known reaction coordinate for a rare event points to umbrella sampling with WHAM; and a high-dimensional or unknown collective variable needing an on-the-fly bias points to metadynamics.</figcaption>
</figure>

## Why free energies are hard

The expectation value of any function $A(\mathbf{q})$ in the canonical ensemble is

$$
\langle A\rangle = \frac{1}{Z} \int d\mathbf{q}\, A(\mathbf{q})\, e^{-\beta U(\mathbf{q})}.
\tag{8.16}
$$

This is a phase-space integral with an exponentially peaked integrand. MD samples regions where $e^{-\beta U}$ is large in proportion to their weight, so $\langle A\rangle$ converges as $1/\sqrt{N_\mathrm{samples}}$ for $A$ that is reasonably smooth.

Free energies are different. $A = -k_B T \ln Z$ depends on the **value** of $Z$, not just on relative weights:

$$
A = -k_B T \ln \int d\mathbf{q}\, e^{-\beta U(\mathbf{q})}.
\tag{8.17}
$$

To compute $A$ directly we would need to know the integrand's value across all of phase space, including regions that MD never visits because $e^{-\beta U}$ is exponentially small there. This is hopeless — but the **difference** $A_2 - A_1$ between two systems is tractable, because it depends on the **ratio** $Z_2/Z_1$, which we can sample.

This section presents four classical recipes for computing free energy differences from MD: thermodynamic integration (TI), free energy perturbation (FEP), umbrella sampling with WHAM, and metadynamics. Each is the right tool somewhere.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. Why can an MD simulation estimate an ordinary ensemble average $\langle A \rangle$ well, but not the absolute free energy $A = -k_B T \ln Z$ directly?
    2. What makes a *difference* of free energies $A_2 - A_1$ tractable when the absolute values are not?
    3. Why does the zero-temperature energy difference alone give the wrong answer for phase stability?

    If any of these is shaky, re-read the preceding section before continuing.

## Why free energies matter

Three flagship applications:

- **Phase stability.** The phase with the lowest Gibbs free energy at $(T, P)$ is the equilibrium phase. To predict whether copper is fcc or hcp at $(T, P)$, compute $G_\mathrm{fcc}$ and $G_\mathrm{hcp}$; the smaller wins. The energy difference alone (zero-temperature) is wrong because it ignores entropy.
- **Reaction and activation barriers.** A catalyst's effectiveness is set by the free energy along the reaction coordinate. Activation free energies $\Delta G^\ddagger$ are the input to transition-state theory and to Arrhenius-rate predictions.
- **Solubility, partitioning, binding.** Drug binding affinities, solubilities, partition coefficients are all free energy differences. Computational drug discovery rests almost entirely on FEP-style calculations.

In each case, what you can extract from MD is a free energy **difference** between two well-defined states.

## Thermodynamic integration

Define a Hamiltonian that interpolates between two endpoints:

$$
H_\lambda(\mathbf{q}, \mathbf{p}) = (1 - \lambda)\, H_0(\mathbf{q}, \mathbf{p}) + \lambda\, H_1(\mathbf{q}, \mathbf{p}), \qquad \lambda \in [0, 1].
\tag{8.18}
$$

(In practice, only the potential is interpolated; the kinetic energy is the same for $H_0$ and $H_1$.) For each $\lambda$, $Z(\lambda) = \int e^{-\beta H_\lambda}$. The free energy difference between endpoints is

$$
\Delta A = A_1 - A_0 = -k_B T \ln \frac{Z_1}{Z_0} = -k_B T \int_0^1 \frac{d \ln Z(\lambda)}{d\lambda} d\lambda.
\tag{8.19}
$$

Differentiate $\ln Z(\lambda)$:

$$
\frac{d \ln Z(\lambda)}{d \lambda} = -\beta\, \frac{\int (\partial H_\lambda/\partial \lambda)\, e^{-\beta H_\lambda}}{\int e^{-\beta H_\lambda}} = -\beta\, \left\langle \frac{\partial H_\lambda}{\partial \lambda}\right\rangle_\lambda.
\tag{8.20}
$$

Substitute:

$$
\boxed{\Delta A = \int_0^1 \left\langle \frac{\partial H_\lambda}{\partial \lambda}\right\rangle_\lambda\, d\lambda.}
\tag{8.21}
$$

This is the **thermodynamic integration** formula. The strategy:

1. Choose a sequence of $\lambda$ values, $\lambda_1 < \lambda_2 < \ldots < \lambda_M$ (often 11 or 21 points equally spaced from 0 to 1).
2. At each $\lambda_k$, run an MD simulation with Hamiltonian $H_{\lambda_k}$, sufficiently long to equilibrate and to estimate $\langle \partial H/\partial \lambda\rangle_{\lambda_k}$ with small uncertainty.
3. Numerically integrate the values of $\langle \partial H/\partial \lambda\rangle$ over $\lambda$ (Simpson's rule or Gauss-Legendre quadrature).

The integrand $\langle \partial H/\partial \lambda\rangle_\lambda$ is just a force-field-style energy gradient; if $H_\lambda = U_0 + \lambda(U_1 - U_0)$, then $\partial H/\partial \lambda = U_1 - U_0$, the energy difference at the current configuration.

TI is robust and well-conditioned when the integrand varies slowly with $\lambda$ — typical case for smooth alchemical transformations between similar molecules. It struggles when the integrand has singularities (e.g., turning on a Lennard-Jones interaction from zero introduces a $\lambda^{-3}$ singularity in $\partial H/\partial \lambda$) which is handled by **soft-core potentials** that regularise the small-$\lambda$ behaviour.

## Free energy perturbation (Zwanzig)

Zwanzig's 1954 identity:

$$
\Delta A = A_1 - A_0 = -k_B T \ln \frac{Z_1}{Z_0} = -k_B T \ln \left\langle e^{-\beta(U_1 - U_0)}\right\rangle_0,
\tag{8.22}
$$

where $\langle \cdot\rangle_0$ is the average **in the ensemble of state 0**. That is, sample with $U_0$, accumulate $e^{-\beta \Delta U}$ at each step, log-average.

The derivation is one line:

$$
\frac{Z_1}{Z_0} = \frac{\int e^{-\beta U_1}}{\int e^{-\beta U_0}} = \frac{\int e^{-\beta(U_1 - U_0)} e^{-\beta U_0}}{\int e^{-\beta U_0}} = \langle e^{-\beta(U_1 - U_0)}\rangle_0.
$$

FEP is conceptually simpler than TI — one simulation, one expectation value — but suffers a severe practical problem: the estimator $e^{-\beta \Delta U}$ has poor statistics when $\Delta U$ is much larger than $k_B T$. The dominant configurations in state 1 may be exponentially rare in the sampled ensemble 0.

**Practical rule.** FEP works well when $|\beta \Delta U|$ is of order 2 or less — i.e., when the two states are similar enough that their canonical distributions overlap substantially. For larger differences, decompose the transformation into a chain of intermediate windows and FEP between adjacent windows (this is the multi-stage variant), which is essentially what TI does.

The Bennett Acceptance Ratio (BAR) and its extension to multiple states (MBAR, Shirts and Chodera) are the modern optimal estimators that improve on naive FEP by combining forward and backward perturbations between windows. They are now standard; `pymbar` is the canonical implementation.

### The Bennett acceptance ratio, in detail

BAR (Bennett, 1976) is the **minimum-variance** estimator of $\Delta A$ from two independent samples — one drawn from state $0$, one from state $1$ — both connected through the energy difference $\Delta U = U_1 - U_0$. Its derivation makes the optimality precise.

Let $\langle\cdot\rangle_0$ and $\langle\cdot\rangle_1$ denote canonical averages in states $0$ and $1$. For *any* function $w(\Delta U)$ with non-vanishing expectations in both ensembles, the ratio identity

$$
\frac{Z_1}{Z_0} = \frac{\langle w\, e^{-\beta U_1}\rangle_0 / \langle e^{-\beta U_1}\rangle_0}{\langle w\, e^{-\beta U_0}\rangle_1 / \langle e^{-\beta U_0}\rangle_1} = \frac{\langle w\, e^{-\beta U_1 + \beta U_0}\rangle_0}{\langle w\, e^{-\beta U_0 + \beta U_1}\rangle_1}
\cdot \frac{Z_0}{Z_0}
\tag{8.B1}
$$

holds because both numerator and denominator are computed by reweighting from their own ensemble. Bennett's insight was to optimise the choice of $w$ for minimum variance of $\ln(Z_1/Z_0)$. The variational calculation — minimise the asymptotic variance subject to (8.B1) — yields a *Fermi function* weight:

$$
w(\Delta U) \propto \frac{1}{1 + e^{\beta(\Delta U - C)}},
\tag{8.B2}
$$

with a constant $C$ to be determined self-consistently. Substituting (8.B2) into the ratio identity and rearranging produces the implicit equation

$$
\boxed{\sum_{n=1}^{N_0} \frac{1}{1 + e^{\beta(\Delta U_n^{(0)} - C)}} = \sum_{m=1}^{N_1} \frac{1}{1 + e^{-\beta(\Delta U_m^{(1)} - C)}},}
\tag{8.B3}
$$

where $\Delta U_n^{(0)}$ is the energy difference evaluated on configuration $n$ drawn from ensemble $0$, and analogously for ensemble $1$. The free energy is then

$$
\Delta A = C - k_B T \ln\frac{N_1}{N_0}.
\tag{8.B4}
$$

Solving (8.B3) is a one-dimensional root-find on $C$, trivial in practice. `pymbar.BAR` does this in one line.

### The maximum-likelihood interpretation

Equation (8.B3) is also the **maximum-likelihood estimator** of $\Delta A$ if one views the samples from $0$ and $1$ as realisations of a logistic-regression problem: given a configuration, predict whether it was drawn from $0$ or $1$. The log-likelihood of the logistic classifier with intercept $C$ is exactly the difference of the two sums in (8.B3), and its score equation is the BAR equation. This connection — re-derived by Shirts, Bair, Hooker and Pande in 2003 — places BAR on the firmest possible statistical footing: any other estimator that uses only the same two-sided samples has equal or worse asymptotic variance.

### Why BAR beats one-sided FEP

One-sided FEP (8.22) is the special case of BAR with $C \to -\infty$, throwing away the backward sample entirely. Two failure modes are then exposed:

- **Distribution overlap.** FEP from $0 \to 1$ relies on configurations sampled in $0$ being representative of $1$. When the canonical distributions of $0$ and $1$ overlap poorly — anywhere $\beta|\Delta U|$ exceeds a few — the estimator is dominated by exponentially rare configurations and converges only slowly, or biases low. The backward direction (FEP from $1 \to 0$) suffers symmetrically. BAR uses both directions and is well-conditioned whenever *either* direction has reasonable overlap.
- **Hysteresis diagnostic.** Forward and backward FEP estimates disagree when the simulations are too short. The discrepancy is itself a diagnostic; BAR collapses the two estimates into a single optimal value and reports a sensible uncertainty.

Empirically, BAR converges roughly four times faster than one-sided FEP for typical alchemical transformations and is essentially never worse. There is no reason to use one-sided FEP in production work.

### Practical use: two-state simulations

A canonical BAR setup interpolates between endpoints $0$ and $1$ through a sequence of $M$ intermediate windows $\lambda_0 = 0 < \lambda_1 < \ldots < \lambda_M = 1$, chosen so that consecutive windows have substantial canonical-distribution overlap. At each window $k$, run MD with potential $U_{\lambda_k}$ and accumulate two energy traces: $\Delta U_{k \to k+1}$ (evaluated on configurations from window $k$) and $\Delta U_{k+1 \to k}$ (from window $k+1$). For each pair of adjacent windows, solve (8.B3) for $\Delta A_{k \to k+1}$; then $\Delta A_{0 \to 1} = \sum_k \Delta A_{k \to k+1}$.

The optimal *spacing* of $\lambda$ values is the one that gives equal overlap (equal pair-BAR variance) between adjacent windows — typically a denser spacing near $\lambda = 0$ when activating Lennard-Jones interactions from zero, and uniform spacing for smooth chemical transformations. An adaptive scheme, where intermediate $\lambda$ values are added on the fly to equalise overlap, is implemented in modern alchemical packages.

### Connection to umbrella sampling and WHAM

The multi-state generalisation of BAR is the **Multistate Bennett Acceptance Ratio** (MBAR; Shirts and Chodera, 2008). Instead of pairwise window combinations, MBAR solves a system of coupled equations for free-energy offsets $f_k = -k_B T \ln Z_k$ across all windows simultaneously, using *every* sample evaluated at *every* window's energy function. The estimator is again maximum-likelihood, and again minimum-variance, but now the variance reduction comes from exploiting *non-adjacent* window overlap when it exists.

The Weighted Histogram Analysis Method (WHAM) of Kumar et al. — used in umbrella sampling as described above — is the discretised (binned) limit of MBAR. The conceptual identity is sharp: WHAM equations are the BAR equations of (8.B3) generalised to many windows, with histogram bins replacing the exact reweighting. Modern practice usually skips the binning step entirely and applies MBAR directly to the raw sample energies, which avoids histogram-bin-width sensitivity and is what `pymbar` does. WHAM remains widespread because of legacy and because for very long simulations the bin-averaged version is faster.

A single unifying picture: BAR, MBAR, WHAM, and the equilibrium reweighting underlying replica exchange are all instances of the same maximum-likelihood estimator applied to the same overlap structure between sampled ensembles. Once you have written one, you have understood them all.

## Umbrella sampling and WHAM

For free energies along a **collective variable** (CV) — a reaction coordinate, an order parameter — neither TI nor FEP is most natural. You want the potential of mean force (PMF)

$$
F(\xi) = -k_B T \ln P(\xi),
\qquad P(\xi) = \frac{1}{Z}\int d\mathbf{q}\, \delta(\xi(\mathbf{q}) - \xi)\, e^{-\beta U(\mathbf{q})},
\tag{8.23}
$$

as a function of the CV $\xi$. Naive sampling estimates $P(\xi)$ by histogramming MD output, but if $F(\xi)$ has a barrier of 10 $k_B T$, the system spends $e^{-10} \approx 5 \times 10^{-5}$ of its time near the barrier top, and the histogram is empty there.

**Umbrella sampling** (Torrie and Valleau, 1977) restores statistics in the barrier region by adding a bias potential $W_i(\xi)$ in each of several "windows" $i$:

$$
U_i^\mathrm{biased}(\mathbf{q}) = U(\mathbf{q}) + W_i(\xi(\mathbf{q})), \qquad W_i(\xi) = \tfrac{1}{2} k\, (\xi - \xi_i^0)^2.
\tag{8.24}
$$

The harmonic restraint keeps the system near $\xi_i^0$ regardless of how high the underlying barrier is. Each window samples a distribution $P_i^\mathrm{biased}(\xi)$, and the unbiased PMF is recovered window-by-window via

$$
P(\xi)\Big|_\mathrm{window\, i} \propto P_i^\mathrm{biased}(\xi)\, e^{+\beta W_i(\xi)}.
\tag{8.25}
$$

Combining the windows into a single PMF over the full range of $\xi$ requires resolving the unknown constant offsets between windows; the **Weighted Histogram Analysis Method** (WHAM, Kumar et al., 1992) does this by solving self-consistently for the offsets. Modern toolkits (`pymbar`, `wham`, `WESTPA`, PLUMED's `wham`) ship with this; you rarely write WHAM yourself.

US/WHAM is the workhorse for free energy profiles along a known reaction coordinate. Typical use cases:

- Dissociation curve of a complex (CV = centre-of-mass distance).
- Conformational change of a polymer (CV = end-to-end distance, or a dihedral).
- Translocation of an ion through a membrane (CV = $z$-position).

The cost is several tens of MD windows of order nanoseconds each, easily parallelisable.

## Metadynamics

Metadynamics (Laio and Parrinello, 2002) takes a different approach: rather than pre-defining windows, it **builds up** the bias adaptively. A Gaussian is deposited along the CV every $\tau$ steps centred at the current value $\xi(t)$, gradually filling up the wells of $F(\xi)$:

$$
V_\mathrm{bias}(\xi, t) = \sum_{t' \le t,\, t' = n\tau} w\, \exp\left[-\frac{(\xi - \xi(t'))^2}{2\sigma^2}\right].
\tag{8.26}
$$

Once the bias is large enough to flatten the underlying $F$, the system diffuses freely along $\xi$ and the bias itself becomes the negative of the free energy (up to a constant). Well-tempered metadynamics (Barducci et al., 2008) modifies the deposition rate adaptively, ensuring convergence to a smooth limit.

Compared to umbrella sampling:
- Metadynamics requires no choice of windows but does require sensible choices of $\sigma$, $w$, $\tau$, and bias factor.
- It extends naturally to two or three CVs (US becomes impractical past 1–2 CVs because the number of windows grows exponentially).
- It is the standard tool for chemistry-flavoured problems in materials simulation — nucleation, phase transitions in 2D order parameters, defect migration paths with multiple coordinates.

PLUMED is the universally adopted library; it interfaces with LAMMPS, GROMACS, NAMD, CP2K, QE, and most ASE-compatible codes.

## A worked example: chemical potential of an LJ liquid

Compute the **excess chemical potential** of a Lennard-Jones liquid at $T^* = 1.5, \rho^* = 0.8$ via thermodynamic integration from the ideal gas to the LJ liquid.

The chemical potential decomposes as $\mu = \mu^\mathrm{id}(\rho, T) + \mu^\mathrm{ex}(\rho, T)$, with $\mu^\mathrm{id}$ the ideal-gas value (analytical) and $\mu^\mathrm{ex}$ the excess due to interactions. Define a $\lambda$-coupled potential with $\lambda$ scaling the LJ interactions:

$$
U_\lambda(\{r_{ij}\}) = \lambda \sum_{i<j} U_\mathrm{LJ}(r_{ij}).
\tag{8.27}
$$

At $\lambda = 0$ atoms are an ideal gas; at $\lambda = 1$ they interact via the full LJ potential. By TI,

$$
\mu^\mathrm{ex} = \frac{1}{N}\int_0^1 \left\langle \sum_{i<j} U_\mathrm{LJ}(r_{ij})\right\rangle_\lambda\, d\lambda = \frac{1}{N}\int_0^1 \langle U\rangle_\lambda\, d\lambda.
\tag{8.28}
$$

Implementation (LAMMPS pseudocode):

```text
units           lj
atom_style      atomic
lattice         fcc 0.8
region          box block 0 6 0 6 0 6
create_box      1 box
create_atoms    1 box
mass            1 1.0

velocity        all create 1.5 12345
timestep        0.005
fix             1 all nvt temp 1.5 1.5 0.5

# Loop over lambda
variable lambda index 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
label loop
pair_style      lj/cut/soft 1.0 0.5 2.5     # soft-core LJ
pair_coeff      * * ${lambda} 1.0 1.0
run             50000                        # equilibration
fix             2 all ave/time 10 1000 10000 c_thermo_pe file U_${lambda}.dat
run             200000                       # production
unfix           2
next            lambda
jump            SELF loop
```

Run this; for each $\lambda$, read $\langle U\rangle/N$ from `U_*.dat`; numerically integrate via Simpson's rule:

```python
import numpy as np
lambdas = np.linspace(0, 1, 11)
U_per_atom = np.array([read_avg(f"U_{lam:.1f}.dat") for lam in lambdas]) / N
mu_ex = np.trapz(U_per_atom, lambdas)
print(f"mu_ex = {mu_ex:.4f} (LJ units)")
```

Reference value for LJ at $T^* = 1.5, \rho^* = 0.8$: $\mu^\mathrm{ex} \approx -1.6 \epsilon$. Your TI should reproduce this to within 1–2%.

A subtle point: the soft-core LJ (the `lj/cut/soft` style) is needed because turning on a $1/r^{12}$ repulsion from zero introduces a singularity in $\partial U/\partial \lambda$ at small $r$. Soft-core regularises the small-$r$ behaviour at low $\lambda$:

$$
U_\mathrm{SC}(r; \lambda) = 4\epsilon\lambda\, \left[\frac{1}{(\alpha(1-\lambda)^2 + (r/\sigma)^6)^2} - \frac{1}{\alpha(1-\lambda)^2 + (r/\sigma)^6}\right],
\tag{8.29}
$$

with $\alpha \approx 0.5$ the soft-core parameter. This is the standard fix and is built into every modern alchemical free energy code (LAMMPS, GROMACS, AMBER).

## Method selection

| Problem | Best method |
|---|---|
| Free energy difference between similar molecules (drug analogue) | FEP/BAR via alchemical mutation |
| Insertion free energy of a small molecule | Widom particle insertion + thermodynamic integration |
| Free energy of a phase transition | Two-phase coexistence (Chapter 8.4) + cross-validation by TI |
| Reaction barrier with known reaction coordinate | Umbrella sampling + WHAM |
| Multi-dimensional reaction surface | Metadynamics or path collective variables |
| Solubility / partitioning | Cycle of TI + appropriate thermodynamic cycle |

The most common rookie mistake is to use FEP across a transformation too large for it to converge ($|\Delta U| \gg k_B T$ in a single step), see a nice-looking number, and not realise the answer is wrong by tens of kJ/mol. Always cross-check forward vs reverse FEP estimates; they should agree to within the BAR error estimate.

!!! warning "Convergence of free energies is fragile"
    A free energy that "looks converged" — i.e., the running average has flattened — may still be wrong if the simulation has not visited all relevant states. Always run multiple independent replicates; if they agree to within their error bars, you have evidence of convergence. If they disagree, the longer-time-scale exploration is incomplete and your numbers are basin-trapped.

## What we have

The four canonical recipes for free energy differences — TI, FEP, US/WHAM, metadynamics — and a sense of which to deploy when. Combined with the ensemble framework of §8.1, they let us compute thermodynamic stability, reaction barriers, and binding affinities from atomistic simulation.

The next section turns from thermodynamics (which the ensemble averages give us) to **dynamics** — transport coefficients, accessed through equilibrium fluctuations via the Green-Kubo relations.
