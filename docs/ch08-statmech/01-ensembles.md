# 8.1 Ensembles and Partition Functions

## A recap of the four standard ensembles

Statistical mechanics describes macroscopic equilibrium as an average over a probability distribution on microstates. Which distribution depends on which macroscopic variables are held fixed.

**Microcanonical (NVE).** Fixed particle number $N$, volume $V$, total energy $E$. All microstates with $H(\mathbf{q}, \mathbf{p}) = E$ are equally likely:

$$
P_\mathrm{NVE}(\mathbf{q}, \mathbf{p}) = \frac{1}{\Omega(N, V, E)} \delta(H - E),
\tag{8.1}
$$

with $\Omega$ the density of states. The connection to thermodynamics is $S(N, V, E) = k_B \ln \Omega$.

**Canonical (NVT).** Fixed $N$, $V$, temperature $T$. The Boltzmann distribution

$$
P_\mathrm{NVT}(\mathbf{q}, \mathbf{p}) = \frac{e^{-\beta H(\mathbf{q}, \mathbf{p})}}{Z(N, V, T)},
\tag{8.2}
$$

with the canonical partition function

$$
Z(N, V, T) = \frac{1}{N!\, h^{3N}} \int d^{3N}q\, d^{3N}p\, e^{-\beta H(\mathbf{q}, \mathbf{p})}.
\tag{8.3}
$$

Connection to thermodynamics: Helmholtz free energy $A = -k_B T \ln Z$.

**Isobaric-isothermal (NPT).** Fixed $N$, pressure $P$, temperature $T$. Volume becomes a dynamical variable:

$$
P_\mathrm{NPT}(\mathbf{q}, \mathbf{p}, V) = \frac{e^{-\beta [H(\mathbf{q}, \mathbf{p}) + PV]}}{\Delta(N, P, T)},
\tag{8.4}
$$

with the isothermal-isobaric partition function

$$
\Delta(N, P, T) = \frac{1}{V_0}\int_0^\infty dV\, e^{-\beta PV}\, Z(N, V, T),
\tag{8.5}
$$

where $V_0$ is a reference volume that makes $\Delta$ dimensionless. Connection: Gibbs free energy $G = -k_B T \ln \Delta$.

**Grand canonical (µVT).** Fixed chemical potential $\mu$, $V$, $T$. Particle number fluctuates:

$$
P_{\mu VT}(N, \mathbf{q}, \mathbf{p}) = \frac{e^{\beta \mu N}\, e^{-\beta H_N}}{\Xi(\mu, V, T)},
\tag{8.6}
$$

with grand partition function

$$
\Xi(\mu, V, T) = \sum_{N=0}^\infty e^{\beta \mu N}\, Z(N, V, T).
\tag{8.7}
$$

Connection: grand potential $\Phi = -k_B T \ln \Xi = -PV$.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. For each of the four standard ensembles, name which macroscopic variables are held fixed and which thermodynamic potential it connects to.
    2. In the canonical ensemble, why does the partition function $Z$ carry the factor $1/(N!\,h^{3N})$ — what does each piece account for?
    3. Why does the choice of ensemble (NVE vs NVT vs NPT) not matter for averages in the thermodynamic limit, but can matter for fluctuations and for small systems?

    If any of these is shaky, re-read the preceding section before continuing.

## Free energies and what to compute them with

Every ensemble has a natural thermodynamic potential. From any free energy, derivatives recover the conjugate variables:

| Ensemble | Potential | First derivatives |
|---|---|---|
| NVE | $-T S(N, V, E)$ | $1/T = (\partial S/\partial E)_{N,V}$ |
| NVT | $A(N, V, T)$ | $S = -(\partial A/\partial T)_{N,V}$, $P = -(\partial A/\partial V)_{N,T}$ |
| NPT | $G(N, P, T)$ | $S = -(\partial G/\partial T)_{N,P}$, $V = (\partial G/\partial P)_{N,T}$ |
| µVT | $\Phi(\mu, V, T)$ | $N = -(\partial \Phi/\partial \mu)_{V,T}$, $P = -(\partial \Phi/\partial V)_{\mu,T}$ |

Computing the **value** of a free energy from simulation is hard — it is a free energy of activation, requiring integration over all microstates. Computing **differences** of free energy between states is the practical tool. The whole of [§8.2](02-free-energy.md) is devoted to methods for computing such differences.

## Which thermostat samples which ensemble

In Chapter 7 we introduced thermostats operationally; here we connect them to the ensembles they sample.

| MD configuration | Ensemble (formally) |
|---|---|
| Velocity-Verlet, no thermostat | NVE |
| Velocity-Verlet + velocity rescaling | None — does not sample any equilibrium distribution |
| Velocity-Verlet + Berendsen thermostat | None — wrong fluctuations |
| Velocity-Verlet + Nosé-Hoover chain | NVT |
| Velocity-Verlet + Langevin | NVT |
| Velocity-Verlet + CSVR (Bussi) | NVT |
| Velocity-Verlet + MTK barostat + NH chain | NPT |
| Hybrid MD/MC with particle insertion/deletion | µVT |

The last item — grand canonical MD — is unusual. Standard MD cannot insert or delete atoms during dynamics (Newton's equations are particle-conserving). To simulate µVT one alternates short MD blocks (for kinetic equilibration) with Monte Carlo insertion/deletion moves (for particle-number fluctuation). This is the **grand canonical molecular dynamics** (GCMD) protocol used to study, for instance, adsorption isotherms.

!!! warning "Sampling an ensemble takes more than the right thermostat"
    A Langevin thermostat will sample the canonical distribution **provided the simulation is ergodic and equilibrated**. A glassy system below its glass transition will not visit the full canonical distribution within accessible MD timescales; you will be sampling a particular metastable basin. The thermostat is necessary, not sufficient.

## Equivalence in the thermodynamic limit

For $N \to \infty$, all the ensembles give the same answer for **intensive** quantities. The reason is that fluctuations in extensive quantities (energy, volume, particle number) scale as $\sqrt{N}$, so their relative magnitude vanishes. Concretely, in NVT the fluctuations of the energy are

$$
\langle (E - \langle E\rangle)^2 \rangle = k_B T^2 C_V,
\tag{8.8}
$$

so $\sigma_E/E \sim N^{-1/2}$. For large $N$, the constraints "fix $E$" (NVE) and "fix $T$" (NVT) both pin the energy to within fluctuations of order $\sqrt{N}$ around the same mean.

In practice, when computing intensive properties at $N \sim 10^3$–$10^5$ on a workstation, the ensemble choice is dictated by convenience and by the property you want, not by which gives the correct answer. The same liquid density and the same radial distribution function will emerge from NVT (at the right $T$) or NPT (at the right $P$) within statistical noise.

For **extensive** quantities or **fluctuations**, the ensemble matters. Heat capacity in NVE differs in form from heat capacity in NVT (the variance of energy is zero in NVE!). Fluctuation formulas always carry the ensemble label.

## Fluctuation formulas

The deepest practical consequence of the ensemble structure is that thermodynamic response functions are equilibrium **fluctuations** in the appropriate ensemble. We derive two examples.

### Heat capacity from energy variance (NVT)

In NVT, the average energy is

$$
\langle E \rangle = \frac{1}{Z}\int d\mathbf{q}\, d\mathbf{p}\, H\, e^{-\beta H} = -\frac{\partial \ln Z}{\partial \beta}.
\tag{8.9}
$$

Differentiate again:

$$
\frac{\partial \langle E\rangle}{\partial \beta} = -\frac{1}{Z}\int d\mathbf{q}\, d\mathbf{p}\, H^2\, e^{-\beta H} + \frac{1}{Z^2}\left(\int H\, e^{-\beta H}\right)^2 = -\langle E^2\rangle + \langle E\rangle^2.
\tag{8.10}
$$

Since $\partial / \partial \beta = -k_B T^2 \partial / \partial T$:

$$
C_V = \frac{\partial \langle E\rangle}{\partial T} = \frac{\langle E^2\rangle - \langle E\rangle^2}{k_B T^2}.
\tag{8.11}
$$

The heat capacity is the equilibrium variance of the total energy divided by $k_B T^2$. To measure $C_V$, run an NVT simulation, collect $E(t)$, and compute its sample variance. No derivative of an average is needed.

This is a powerful method because variances converge faster than derivatives. A direct calculation of $C_V$ by computing $\langle E\rangle$ at $T$ and $T + \Delta T$ and finite-differencing requires statistical noise much smaller than $\Delta T \cdot C_V$ — a tight requirement. The fluctuation formula gets $C_V$ from a single run.

### Compressibility from volume variance (NPT)

Analogous derivation in NPT gives

$$
\kappa_T = -\frac{1}{V}\left(\frac{\partial V}{\partial P}\right)_T = \frac{\langle V^2\rangle - \langle V\rangle^2}{k_B T\, \langle V\rangle}.
\tag{8.12}
$$

The isothermal compressibility is the volume variance over $k_B T \langle V\rangle$. NPT simulations therefore deliver $\kappa_T$ for free, with no need to vary $P$.

### Specific heat at constant pressure

$$
C_P = \frac{\langle H^2\rangle_\mathrm{NPT} - \langle H\rangle_\mathrm{NPT}^2}{k_B T^2}, \quad H = E + PV.
\tag{8.13}
$$

The relevant variance is that of the **enthalpy** in NPT, not the energy.

### Pair correlations and thermodynamic derivatives

The radial distribution function $g(r)$ of [§7.6](../ch07-md/06-analysis.md) is itself an ensemble average. From $g(r)$ many thermodynamic quantities follow exactly via the **pressure equation**

$$
P = \rho k_B T - \frac{2\pi \rho^2}{3} \int_0^\infty r\, \frac{dU}{dr}\, g(r)\, r^2\, dr,
\tag{8.14}
$$

and the **energy equation**

$$
\frac{\langle U\rangle}{N} = \frac{1}{2}\, 4\pi \rho \int_0^\infty U(r)\, g(r)\, r^2\, dr.
\tag{8.15}
$$

Both are exact for pairwise potentials. The energy equation provides a cross-check on energy averages from the simulation; if the value of $\langle U \rangle / N$ from the LAMMPS log disagrees with the integral over $g(r)$, you have a bug.

For many-body potentials (EAM, ML potentials), the energy equation generalises in a more involved way; the pressure equation generalises to the virial expression already used in (7.50).

!!! note "Variance estimators are noisier than they look"
    Fluctuation formulas are unbiased but high-variance. For $C_V$ in NVT, the **variance of the variance** scales as $T^4 C_V^2 / N$ — so a 1% error on $C_V$ at $N = 1000$ requires aggregating roughly $10^4$ independent samples. With autocorrelation in the energy trace, that may mean nanoseconds of MD. Block averaging (Flyvbjerg-Petersen) is the standard tool for computing error bars on variances.

## Ergodicity and equilibration

All the formulas above assume the simulation is **ergodic** — that the time-average over a finite trajectory converges to the ensemble average. Equivalently, the trajectory must visit every microstate consistent with the constraints, with probabilities given by the ensemble weight.

Ergodicity is the silent assumption that breaks more MD simulations than any other. Symptoms:

- A protein that never escapes its starting conformation in 100 ns. The PMF along a CV (Chapter 8.2) shows a single basin even though there are obviously many.
- A glass below $T_g$ that maintains memory of its preparation cooling rate. Different cooling histories give different "equilibrium" properties.
- A defect cluster trapped in a metastable configuration; the simulation reports a free energy that is really the free energy of that one basin, not of the equilibrium ensemble.

For ergodicity-broken systems, the simulation samples a **constrained** Boltzmann distribution — restricted to one basin. This is sometimes what you want (the local free energy of a particular structural minimum) and sometimes catastrophic (the wrong absolute phase stability). Enhanced sampling methods (umbrella sampling, metadynamics, replica exchange) are responses to this problem.

A diagnostic: re-run the simulation with two different initial conditions and check whether the results agree. If a property differs by more than its statistical error bar between two independent runs, ergodicity has not been achieved on accessible timescales.

## What ensemble for what property?

A pragmatic guide:

| Property | Best ensemble | Why |
|---|---|---|
| Total energy at temperature $T$ | NVT | Directly samples $\langle E \rangle$ |
| Equilibrium volume at $(T, P)$ | NPT | Directly samples $\langle V\rangle$ |
| Pressure at $(T, V)$ | NVT | Directly samples $\langle P\rangle$ via virial |
| $C_V$ | NVT | Variance of $E$ |
| $C_P$ | NPT | Variance of $H = E + PV$ |
| $\kappa_T$ | NPT | Variance of $V$ |
| Adsorption isotherm | µVT (GCMD) | Particle number must fluctuate |
| Free energy difference | NVT or NPT with TI/FEP | See [§8.2](02-free-energy.md) |
| Diffusion coefficient | NVT (low $\gamma$) or NVE | Kinetics must not be distorted |
| Phonon spectrum | NVE | Microcanonical preserves vibrational modes |

The diffusion-coefficient row deserves emphasis: thermostats with large friction distort kinetics. For computing $D$ via MSD or Green-Kubo, either use NVE after equilibration (sampling NVE with a properly equilibrated initial condition) or use NVT with the gentlest possible thermostat (small Langevin friction, or large Nosé-Hoover time constant). LAMMPS users frequently equilibrate in NPT, then switch to NVE for production data collection; this is good practice.

## Metropolis Monte Carlo and the canonical ensemble

MD is not the only way to sample a Boltzmann distribution. The **Metropolis algorithm** (Metropolis, Rosenbluth, Rosenbluth, Teller and Teller, 1953) generates a Markov chain on configuration space whose stationary distribution is exactly $P_\mathrm{NVT}(\mathbf{q}) \propto e^{-\beta U(\mathbf{q})}$. The derivation is short, and worth doing once: it is the same logic underlying parallel tempering, replica exchange, the basin-hopping global optimisers of [Chapter 9](../ch09-mlip/index.md), and — pertinently for [Chapter 11](../ch11-active/index.md) — Metropolis-Hastings acquisition in Bayesian optimisation.

### Detailed balance from scratch

Consider a Markov chain with transition matrix $T_{i \to j}$ — the probability that, given the system is currently in microstate $i$, the next step moves it to $j$. If $P_i$ is the stationary distribution then by definition

$$
P_j = \sum_i P_i\, T_{i \to j},\qquad \sum_j T_{i \to j} = 1.
\tag{8.A1}
$$

A **sufficient** (not necessary) condition for $P_i$ to be stationary is **detailed balance**:

$$
P_i\, T_{i \to j} = P_j\, T_{j \to i}\qquad \forall i, j.
\tag{8.A2}
$$

Summing (8.A2) over $i$ reproduces (8.A1); the converse is not true. Detailed balance is the stronger statement that every microscopic transition is balanced by its time-reverse — equilibrium in the sense of zero net probability flow.

### The Metropolis acceptance criterion

Factor the transition as a *proposal* times an *acceptance*: $T_{i \to j} = g_{i \to j}\, A_{i \to j}$, where $g_{i \to j}$ is the probability of *proposing* move $i \to j$ (chosen by the simulator — e.g., displace a random atom by a uniform random vector) and $A_{i \to j} \in [0, 1]$ is the probability of accepting that proposal. Choose a *symmetric* proposal: $g_{i \to j} = g_{j \to i}$. Detailed balance then reduces to

$$
\frac{A_{i \to j}}{A_{j \to i}} = \frac{P_j}{P_i} = e^{-\beta (U_j - U_i)} = e^{-\beta \Delta U}.
\tag{8.A3}
$$

Any pair of acceptance probabilities satisfying (8.A3) gives a correct sampler. Metropolis chose the maximal-acceptance pair:

$$
\boxed{A_{i \to j} = \min\!\left(1,\, e^{-\beta \Delta U}\right).}
\tag{8.A4}
$$

When $\Delta U \le 0$ the move is accepted with probability $1$ (always go downhill); when $\Delta U > 0$ it is accepted with probability $e^{-\beta \Delta U}$, decaying exponentially with the cost in $k_B T$. Verifying (8.A3) is one line: if $\Delta U > 0$ then $A_{i \to j} = e^{-\beta\Delta U}$ and $A_{j \to i} = 1$, so the ratio is $e^{-\beta\Delta U}$ as required; the case $\Delta U < 0$ is symmetric.

The crucial observation is that the acceptance depends only on $\Delta U$, not on the absolute energies or on the partition function $Z$. That is why MC samples the canonical distribution at all: $Z$ cancels in every ratio.

### The algorithm

```text
choose initial configuration q_0
for step = 1, 2, ...:
    propose q' = q + dq          # symmetric random displacement
    DeltaU = U(q') - U(q)
    if DeltaU <= 0:
        accept: q <- q'
    else:
        draw u ~ Uniform(0, 1)
        if u < exp(-beta * DeltaU):
            accept: q <- q'
        else:
            reject: q unchanged
    record q  # for ensemble averages
```

Rejection is *not* idle: the rejected step contributes the same configuration to the running ensemble average, weighted by the time it spent there. Skipping the recording is a common bug.

A target acceptance rate of $30$–$50$% is the conventional sweet spot: a much higher rate means displacements are too small (slow decorrelation); a much lower rate means too many proposals are wasted. Tune the displacement step size during equilibration to land in this band.

### A 60-line 2D Ising demo

The two-dimensional Ising model on an $L \times L$ square lattice with periodic boundaries is the canonical test bed. Spins $s_i \in \{-1, +1\}$, energy $U = -J \sum_{\langle ij\rangle} s_i s_j$, exact critical temperature $T_c = 2J / [k_B \ln(1 + \sqrt{2})] \approx 2.269\, J/k_B$ (Onsager 1944). Below $T_c$ the system magnetises spontaneously; above $T_c$, $\langle m \rangle = 0$.

```python
from __future__ import annotations

import numpy as np


def ising_mc(
    L: int = 32,
    beta: float = 0.44,
    n_sweeps: int = 4000,
    n_burn: int = 1000,
    seed: int = 0,
) -> tuple[float, np.ndarray]:
    """Single-spin-flip Metropolis MC on the 2D Ising model.

    Returns the mean absolute magnetisation per spin (post-burn-in)
    and the magnetisation time series in units of [1].
    """
    rng = np.random.default_rng(seed)
    s = rng.choice([-1, 1], size=(L, L)).astype(np.int8)
    mag_trace = np.empty(n_sweeps, dtype=np.float64)

    # Precompute acceptance probabilities for the only possible
    # values of DeltaU in 2D Ising: -8, -4, 0, +4, +8 (units of J).
    accept = {dU: min(1.0, np.exp(-beta * dU)) for dU in (-8, -4, 0, 4, 8)}

    for sweep in range(n_sweeps + n_burn):
        # One sweep = L*L attempted flips.
        for _ in range(L * L):
            i, j = rng.integers(0, L, size=2)
            nb_sum = (
                s[(i + 1) % L, j] + s[(i - 1) % L, j]
                + s[i, (j + 1) % L] + s[i, (j - 1) % L]
            )
            dU = int(2 * s[i, j] * nb_sum)        # J = 1
            if rng.random() < accept[dU]:
                s[i, j] = -s[i, j]
        if sweep >= n_burn:
            mag_trace[sweep - n_burn] = s.mean()

    return float(np.abs(mag_trace).mean()), mag_trace


if __name__ == "__main__":
    # Sweep T from below to above the Onsager critical point.
    Ts = np.linspace(1.5, 3.5, 21)
    m_of_T = np.array([ising_mc(beta=1.0 / T)[0] for T in Ts])
    for T, m in zip(Ts, m_of_T):
        print(f"T = {T:.2f}   <|m|> = {m:.3f}")
```

The output shows $\langle |m|\rangle$ falling from near unity (ordered) to near $1/L$ (disordered fluctuations of a finite paramagnet) as $T$ crosses $T_c \approx 2.27$. The curve has the characteristic Onsager shape, and the finite-size rounding around $T_c$ is the textbook demonstration that critical singularities only become singular in the thermodynamic limit.

This sixty-line program contains everything: detailed-balance acceptance, the canonical-ensemble logic, and a real second-order phase transition. The same loop, with continuous-coordinate displacements and an inter-atomic potential, is the entire Metropolis MC method for fluids.

### When MC beats MD

For atomic systems Metropolis MC has three advantages over MD:

- **No forces.** Only $\Delta U$ is needed, which is often cheaper than the full gradient — useful for complicated potentials or DFT-on-the-fly methods.
- **Non-physical moves.** Particle insertions/deletions ($\mu VT$), atom swaps (alloy ordering), large rotational moves for molecules — none of these have natural MD analogues. Hybrid MD/MC samplers use MD for local equilibration and MC for moves that break the conservation laws Newton's equations enforce.
- **No timestep stability.** Hard cores, near-singular potentials, stiff bonds — MC handles them through rejections, where MD would explode.

The disadvantage is that MC gives no dynamics: there is no time, no diffusion coefficient, no transport. MC samples the equilibrium distribution; if you want kinetics, use MD.

!!! note "Metropolis-Hastings, the asymmetric generalisation"
    If the proposal is *not* symmetric — for instance, biased toward downhill moves, or drawn from a non-uniform distribution — the acceptance must include the proposal ratio:
    $$
    A_{i \to j} = \min\!\left(1,\, \frac{g_{j \to i}}{g_{i \to j}} e^{-\beta \Delta U}\right).
    $$
    This is the **Metropolis-Hastings** algorithm (Hastings, 1970). It is the engine underneath modern Markov-chain-Monte-Carlo samplers (HMC, NUTS in PyMC, parallel tempering), and — relevant for [Chapter 11](../ch11-active/index.md) — underneath the Thompson sampling and acquisition-function machinery used in Bayesian optimisation, where the "energy" $\beta U$ is replaced by a negative log-posterior.

## What we have

A precise vocabulary for what each MD ensemble samples, and for what quantities each ensemble gives most cleanly. The next two sections build on this: free energies ([§8.2](02-free-energy.md)) require sampling that **interpolates between** ensembles or potentials, and transport coefficients ([§8.3](03-transport.md)) come from equilibrium fluctuations of currents rather than of state variables.
