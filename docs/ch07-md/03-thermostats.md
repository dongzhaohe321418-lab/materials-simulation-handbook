# 7.3 Thermostats and Barostats

## The NVE ensemble

Pure velocity-Verlet on a closed system samples the **microcanonical ensemble**: fixed number of particles $N$, fixed volume $V$, fixed total energy $E$. The conserved quantity is the Hamiltonian itself,

$$
H(\{\mathbf{r}_i\}, \{\mathbf{p}_i\}) = \sum_i \frac{|\mathbf{p}_i|^2}{2 m_i} + U(\{\mathbf{r}_i\}).
\tag{7.25}
$$

Over a long ergodic trajectory the distribution of microstates that NVE samples is the constant-energy shell,

$$
P_\mathrm{NVE}(\{\mathbf{r}_i\}, \{\mathbf{p}_i\}) \propto \delta(H - E).
\tag{7.26}
$$

NVE is what Verlet gives you "for free". It is the right ensemble for studying energy flow and for benchmarking integrators. It is the wrong ensemble for almost any real material property, because experiments are done at constant $T$ or constant $T$ and $P$.

## Temperature from kinetic energy

In equilibrium the equipartition theorem gives, for any quadratic degree of freedom,

$$
\langle \tfrac{1}{2} m v_\alpha^2 \rangle = \tfrac{1}{2} k_B T.
\tag{7.27}
$$

For $N$ atoms in three dimensions with no external constraints, the total kinetic energy is

$$
\langle KE \rangle = \tfrac{3}{2} N k_B T,
\tag{7.28}
$$

from which the **instantaneous kinetic temperature** is

$$
T_\mathrm{kin} = \frac{2 \langle KE \rangle}{3 N k_B} = \frac{1}{3 N k_B} \sum_i m_i |\mathbf{v}_i|^2.
\tag{7.29}
$$

Strictly, the denominator should use the number of **dynamical** degrees of freedom $N_\mathrm{df}$. Subtract 3 for the centre-of-mass momentum if it is conserved at zero; subtract another 3 for total angular momentum in a non-periodic system; subtract one for each holonomic constraint (SHAKE). For large $N$ the difference is negligible; for $N \sim 10$ it matters.

The instantaneous $T_\mathrm{kin}$ fluctuates wildly — its standard deviation is $T\sqrt{2/(3N)}$. For $N = 1000$, the fluctuation is about 2.6% of $T$. A reported temperature is always a time-average.

## Velocity rescaling

The simplest temperature control is brute rescaling. Every $n_\mathrm{rescale}$ steps, measure $T_\mathrm{kin}$ and scale all velocities by

$$
\lambda = \sqrt{T_\mathrm{target}/T_\mathrm{kin}}.
\tag{7.30}
$$

This pins the kinetic energy to exactly $\tfrac{3}{2} N k_B T$ at the rescaling steps. It is convenient for **equilibration** — heating up an initial configuration in a controlled way. It is wrong for **sampling**: it produces a distribution with zero variance in $T_\mathrm{kin}$, whereas the canonical ensemble has $\sigma_T^2 = 2T^2/(3N)$. The phase-space distribution is not the canonical $e^{-\beta H}$.

Use velocity rescaling only during equilibration, never during production runs from which you compute averages.

## Berendsen thermostat

The Berendsen scheme softens the rescaling. Every step, rescale velocities by

$$
\lambda = \sqrt{1 + \frac{\Delta t}{\tau_T}\left(\frac{T_\mathrm{target}}{T_\mathrm{kin}} - 1\right)},
\tag{7.31}
$$

so the kinetic temperature relaxes exponentially to its target with time constant $\tau_T$:

$$
\frac{dT_\mathrm{kin}}{dt} = \frac{T_\mathrm{target} - T_\mathrm{kin}}{\tau_T}.
\tag{7.32}
$$

For $\tau_T \to 0$ Berendsen reduces to rescaling; for $\tau_T \to \infty$ it becomes NVE. A typical choice is $\tau_T = 0.1$–1 ps.

Berendsen is the most commonly misused thermostat in MD. It does the right thing macroscopically (the mean temperature is correct) but it does **not** generate the canonical distribution: it suppresses temperature fluctuations to about $\sqrt{1/N_\mathrm{df}}$ of their true value, and energy fluctuations correspondingly. Heat capacities computed from energy variance with a Berendsen thermostat are wrong by a factor that depends on $\tau_T$. The "flying ice cube" problem — anomalous accumulation of kinetic energy in centre-of-mass and rotational modes — is also a Berendsen artefact.

!!! warning "Berendsen is for equilibration only"
    Use Berendsen to bring a system to temperature quickly without crashing it. Then switch to a thermostat that samples the canonical ensemble correctly (Nosé-Hoover, Langevin, or stochastic velocity rescaling) for production.

A modern fix is the **stochastic velocity rescaling** thermostat (Bussi-Donadio-Parrinello, 2007), which adds a single stochastic kick at each rescaling step chosen so that the canonical kinetic-energy distribution is exactly preserved. It is essentially as cheap as Berendsen and is correct. Most LAMMPS users today reach for `fix temp/csvr` or `fix langevin` rather than Berendsen.

## Nosé-Hoover: extended Lagrangian for NVT

The right way to sample the canonical ensemble deterministically is to add a single extra degree of freedom that exchanges energy with the physical system. Nosé's 1984 construction introduces a fictitious variable $s$ with conjugate momentum $p_s$ and writes the extended Hamiltonian

$$
H_\mathrm{Nos\'e} = \sum_i \frac{|\mathbf{p}_i|^2}{2 m_i s^2} + U(\{\mathbf{r}_i\}) + \frac{p_s^2}{2 Q} + g\, k_B T \ln s,
\tag{7.33}
$$

where $Q$ is a "mass" for the new variable (with units of energy·time$^2$) and $g$ is chosen as $N_\mathrm{df} + 1$ so that integrating out $s$ gives the canonical distribution in the physical variables. The role of $s$ is to scale momenta: $\mathbf{p}_i^\mathrm{real} = \mathbf{p}_i / s$.

Hoover's reformulation (1985) eliminates $s$ in favour of a friction variable $\zeta = p_s s/Q$, giving the **Nosé-Hoover equations of motion** in real-time variables:

$$
\dot{\mathbf{r}}_i = \frac{\mathbf{p}_i}{m_i},
\tag{7.34}
$$

$$
\dot{\mathbf{p}}_i = -\nabla_i U - \zeta\, \mathbf{p}_i,
\tag{7.35}
$$

$$
\dot{\zeta} = \frac{1}{Q}\left(\sum_i \frac{|\mathbf{p}_i|^2}{m_i} - g k_B T\right) = \frac{2}{Q}\left(KE - \tfrac{g}{2} k_B T\right).
\tag{7.36}
$$

The interpretation is transparent. Equation (7.35) adds a friction term to Newton's equation, drawing kinetic energy out of the system when $\zeta > 0$. Equation (7.36) drives $\zeta$ towards positive values when the kinetic energy exceeds the target $\tfrac{1}{2} g k_B T$, and towards negative values when it falls short. The friction itself is a dynamical variable that fluctuates around zero in equilibrium.

The conserved quantity (not the Hamiltonian, but a related construct) is

$$
H' = H_\mathrm{phys} + \tfrac{1}{2} Q \zeta^2 + g k_B T \int_0^t \zeta(t')\, dt',
\tag{7.37}
$$

which provides a diagnostic: $H'$ should be flat in a correct NH simulation.

**Choice of $Q$.** The "thermostat mass" sets the timescale on which $\zeta$ responds. Too small a $Q$ gives stiff oscillations of $\zeta$ that the integrator cannot resolve; too large a $Q$ gives slow coupling and ergodicity problems. A common choice is

$$
Q = g k_B T \tau_T^2,
\tag{7.38}
$$

with $\tau_T \approx 0.1$–1 ps — same as Berendsen.

**Ergodicity and the chain extension.** A pure NH thermostat is **not ergodic** for small systems: a single harmonic oscillator coupled to a NH thermostat produces a non-canonical distribution because the conserved $H'$ over-constrains the dynamics. Martyna, Klein and Tuckerman (1992) solved this by chaining $M$ Nosé-Hoover thermostats — the first heats the physical system, the second heats the first, and so on. The chain shakes loose the spurious conservation law and restores ergodicity. In production MD, $M = 3$–5 is universal. LAMMPS' `fix nvt` uses a chain of 3 by default; you rarely need to touch this.

!!! note "Nosé-Hoover for an isolated oscillator"
    If you run a single 1D oscillator with NH and a single thermostat, the kinetic-energy distribution will not be canonical — you can verify this by histogram. Switch to a 3-chain, and the histogram converges to the correct $\chi^2$ form. This is the textbook test of an MD code's NVT implementation.

## Langevin thermostat

A different approach to canonical sampling is **stochastic**. Replace Newton's equation with a Langevin equation:

$$
m_i \ddot{\mathbf{r}}_i = -\nabla_i U - \gamma\, m_i \dot{\mathbf{r}}_i + \boldsymbol{\eta}_i(t),
\tag{7.39}
$$

with friction $\gamma$ and a Gaussian white-noise force

$$
\langle \eta_{i\alpha}(t) \rangle = 0,
\qquad
\langle \eta_{i\alpha}(t)\, \eta_{j\beta}(t') \rangle = 2 m_i \gamma k_B T \,\delta_{ij} \delta_{\alpha\beta} \delta(t - t').
\tag{7.40}
$$

The fluctuation-dissipation relation (7.40) ties the noise amplitude to the friction so that the stationary distribution is exactly canonical. Each Cartesian component of each atom is connected to its own thermal bath; ergodicity is built in.

Langevin dynamics samples NVT correctly for any $\gamma > 0$. Small $\gamma$ gives weak coupling and Newtonian-like dynamics on short timescales; large $\gamma$ overdamps the dynamics into Brownian motion. For most condensed-phase systems $\gamma \sim 1$–10 ps$^{-1}$ is a good range. Note that diffusion coefficients and other dynamical observables **are affected** by $\gamma$ — Langevin distorts kinetics — whereas static thermodynamic averages are correct for any $\gamma$.

If you need correct kinetics, use Nosé-Hoover; if you only need correct static averages and you want robust ergodicity, use Langevin.

## Parrinello-Rahman barostat: NPT

To control pressure as well as temperature, the simulation cell must also become dynamical. Parrinello and Rahman (1980, 1981) elevated the lattice matrix $\mathbf{H} = (\mathbf{a}_1\;\mathbf{a}_2\;\mathbf{a}_3)$ to a set of dynamical variables, with its own kinetic and potential energy contributions. The extended Lagrangian is

$$
\mathcal{L} = \sum_i \tfrac{1}{2} m_i \dot{\mathbf{s}}_i^\top \mathbf{H}^\top \mathbf{H}\, \dot{\mathbf{s}}_i - U(\mathbf{H}\, \mathbf{s}) + \tfrac{1}{2} W\, \mathrm{Tr}(\dot{\mathbf{H}}^\top \dot{\mathbf{H}}) - P_\mathrm{ext} V,
\tag{7.41}
$$

where $\mathbf{s}_i$ are fractional coordinates ($\mathbf{r}_i = \mathbf{H} \mathbf{s}_i$), $V = \det \mathbf{H}$, $W$ is a fictitious cell mass, and $P_\mathrm{ext}$ is the target pressure. The Euler-Lagrange equations give the cell-vector dynamics

$$
W \ddot{\mathbf{H}} = V\, (\boldsymbol{\sigma} - P_\mathrm{ext} \mathbf{1})\, \mathbf{H}^{-\top},
\tag{7.42}
$$

with $\boldsymbol{\sigma}$ the instantaneous stress tensor. Cell vectors respond to imbalances between internal stress and external pressure; the cell can shrink, grow, and shear in response to the underlying material.

Combined with a Nosé-Hoover thermostat on both the atomic and cell degrees of freedom (the Martyna-Tobias-Klein, "MTK", scheme is the standard reference), this gives a deterministic NPT integrator that correctly samples the isothermal-isobaric ensemble. LAMMPS' `fix npt` implements MTK.

**Choice of $W$.** The cell mass controls how fast the box responds. Too small $W$: the box rings at unphysical frequencies. Too large $W$: the box can't follow density changes during equilibration. A standard scaling is

$$
W = (N_\mathrm{df} + 3) k_B T\, \tau_P^2,
\tag{7.43}
$$

with $\tau_P \approx 0.5$–5 ps. As with $\tau_T$, defaults are usually fine.

**Isotropic vs anisotropic cells.** For liquids and cubic solids, isotropic pressure coupling (`fix npt iso`) constrains the cell to remain a scaled copy of itself, varying only the volume. For non-cubic crystals, anisotropic coupling (`fix npt aniso`) allows the three cell edges to fluctuate independently while keeping the angles fixed; full triclinic (`fix npt tri`) lets shear strains develop as well. Use the most restrictive setting consistent with your symmetry — extra cell degrees of freedom add noise without information.

!!! warning "NPT during phase transitions"
    Across a first-order transition the equilibrium volume jumps discontinuously. An NPT simulation crossing the transition will not give a sharp jump in the volume — it will give a slow drift over many nanoseconds, because the new phase has to nucleate. Don't mistake the slow drift for "not equilibrated yet" and abandon the simulation. The two-phase method ([§8.4](../ch08-statmech/04-phase-diagrams.md)) is the proper tool here.

## A summary table

| Thermostat | Ergodic? | Canonical distribution? | Affects kinetics? | Use case |
|---|---|---|---|---|
| Velocity rescaling | Yes | No | Yes | Equilibration only |
| Berendsen | Yes | No (suppressed fluctuations) | Mildly | Equilibration only |
| Stochastic velocity rescaling (CSVR) | Yes | Yes | Minimal | Production NVT |
| Nosé-Hoover (chained) | Yes (with chain) | Yes | No | Production NVT, correct kinetics |
| Langevin | Yes | Yes | Yes (depends on $\gamma$) | Production NVT, robust |

| Barostat | Ensemble | Notes |
|---|---|---|
| Berendsen | Approximate NPT | Equilibration only; cell relaxes exponentially |
| Parrinello-Rahman + MTK | NPT, correct | Full cell dynamics, standard for production |
| MC volume moves | NPT, correct | Used in hybrid MD/MC; less common |

## What we have

We can now sample NVE, NVT, and NPT ensembles, given a force model. We have not specified what produces those forces. That is [§7.4](04-force-fields.md): classical force fields, from the cartoon (Lennard-Jones) to the production (EAM, Tersoff, ReaxFF). The chapter on machine-learning interatomic potentials ([Chapter 9](../ch09-mlip/index.md)) will then take the next step, learning forces directly from DFT.
