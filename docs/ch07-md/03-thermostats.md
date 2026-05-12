# 7.3 Thermostats and Barostats

!!! tip "Why thermostats exist"
    Velocity-Verlet gives you NVE — constant number of particles, constant volume, constant total energy. The universe rarely cares about NVE. Experiments are done at a controlled temperature (a coffee cup at 300 K, a furnace at 1500 K) or a controlled temperature and pressure (atmospheric conditions). You cannot directly impose a temperature on Hamilton's equations of motion; "temperature" is a *statistical* property of an ensemble of trajectories.
    
    A thermostat is a *modification* of the dynamics that makes the resulting trajectory sample the canonical ensemble $P \propto e^{-\beta H}$ rather than the microcanonical $P \propto \delta(H-E)$. Different thermostats embody different philosophies:
    
    - **Velocity rescaling**: brute force — multiply velocities by a constant to pin $T$.
    - **Berendsen**: smooth rescaling — relax $T$ to target exponentially.
    - **Nosé-Hoover**: clever trick — invent a fictitious extra degree of freedom that couples to all velocities; the combined extended system has a Hamiltonian, and projecting onto the physical variables gives canonical NVT.
    - **Langevin**: stochastic — add random kicks that explicitly maintain Boltzmann statistics through fluctuation-dissipation.
    
    The deepest of these is Nosé's. His 1984 insight: if you cannot impose temperature on Hamilton's equations directly, *embed* the system in a larger one that still obeys Hamilton's equations but whose *marginal* distribution is canonical. The price is one extra "thermostat" variable; the payoff is exact NVT sampling with deterministic dynamics. The "aha" picture: a thermostat as a black box outside your system, exchanging energy with it on a controlled timescale, but the black box itself follows the rules of mechanics.

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

#### Why velocity rescaling is not canonical

The canonical (NVT) distribution for $N$ atoms in 3D is
$$P_\mathrm{NVT}(\{\mathbf{r}_i\},\{\mathbf{p}_i\}) \propto e^{-\beta H},$$
which factorises into a configurational part $e^{-\beta U}$ and a Maxwell-Boltzmann momentum part. The Maxwell-Boltzmann momentum distribution has *fluctuations* in the total kinetic energy:
$$P(KE) \propto KE^{3N/2 - 1}\, e^{-\beta KE},$$
with mean $\langle KE\rangle = \tfrac{3}{2}N k_BT$ and variance $\sigma_{KE}^2 = \tfrac{3}{2}N(k_BT)^2$. The relative fluctuation $\sigma_{KE}/\langle KE\rangle = \sqrt{2/(3N)}$ is small for large $N$ but non-zero.

Velocity rescaling kills these fluctuations at the rescaling steps: $KE$ is *exactly* $\tfrac{3}{2}Nk_BT$ then, with no variance. Even between rescaling steps the fluctuation is only what the (deterministic) dynamics can produce starting from a zero-variance initial condition, which is much less than canonical. Heat capacities computed from $\sigma_{KE}^2$ then come out wrong: $C_V^\mathrm{ideal} = \tfrac{3}{2}Nk_B = \sigma_{KE}^2/k_BT^2$ is severely underestimated.

The same critique applies in spirit to any thermostat that imposes a hard constraint on $KE$ rather than letting it fluctuate canonically. The way out is either (a) inject stochastic kicks that explicitly maintain the right $\chi^2$ distribution on $KE$ (the Bussi-Donadio-Parrinello "CSVR" thermostat), or (b) let the thermostat itself fluctuate as a dynamical variable (Nosé-Hoover).

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

!!! note "Deriving the Berendsen relaxation equation"
    Equation (7.32) is what the scaling (7.31) implements at first order in $\Delta t$. To see this, write the kinetic energy after rescaling: $KE^\mathrm{new} = \lambda^2 KE^\mathrm{old}$, so $T_\mathrm{kin}^\mathrm{new} = \lambda^2 T_\mathrm{kin}^\mathrm{old}$. Substituting (7.31),
    $$T_\mathrm{kin}^\mathrm{new} = \left[1 + \frac{\Delta t}{\tau_T}\!\left(\frac{T_\mathrm{target}}{T_\mathrm{kin}^\mathrm{old}} - 1\right)\right] T_\mathrm{kin}^\mathrm{old} = T_\mathrm{kin}^\mathrm{old} + \frac{\Delta t}{\tau_T}\!\left(T_\mathrm{target} - T_\mathrm{kin}^\mathrm{old}\right).$$
    Rearranging:
    $$\frac{T_\mathrm{kin}^\mathrm{new} - T_\mathrm{kin}^\mathrm{old}}{\Delta t} = \frac{T_\mathrm{target} - T_\mathrm{kin}}{\tau_T},$$
    which in the continuum limit is (7.32). The solution is $T_\mathrm{kin}(t) = T_\mathrm{target} + (T_0 - T_\mathrm{target})e^{-t/\tau_T}$ — first-order exponential approach.
    
    This is exactly the behaviour you want for equilibration: monotonic, no overshoot, controlled rate. It is exactly the behaviour you do *not* want for sampling: fluctuations are damped to zero on timescale $\tau_T$, not maintained at their canonical values.

For $\tau_T \to 0$ Berendsen reduces to rescaling; for $\tau_T \to \infty$ it becomes NVE. A typical choice is $\tau_T = 0.1$–1 ps.

Berendsen is the most commonly misused thermostat in MD. It does the right thing macroscopically (the mean temperature is correct) but it does **not** generate the canonical distribution: it suppresses temperature fluctuations to about $\sqrt{1/N_\mathrm{df}}$ of their true value, and energy fluctuations correspondingly. Heat capacities computed from energy variance with a Berendsen thermostat are wrong by a factor that depends on $\tau_T$. The "flying ice cube" problem — anomalous accumulation of kinetic energy in centre-of-mass and rotational modes — is also a Berendsen artefact.

!!! warning "Berendsen is for equilibration only"
    Use Berendsen to bring a system to temperature quickly without crashing it. Then switch to a thermostat that samples the canonical ensemble correctly (Nosé-Hoover, Langevin, or stochastic velocity rescaling) for production.

A modern fix is the **stochastic velocity rescaling** thermostat (Bussi-Donadio-Parrinello, 2007), which adds a single stochastic kick at each rescaling step chosen so that the canonical kinetic-energy distribution is exactly preserved. It is essentially as cheap as Berendsen and is correct. Most LAMMPS users today reach for `fix temp/csvr` or `fix langevin` rather than Berendsen.

## Nosé-Hoover: extended Lagrangian for NVT

!!! abstract "Key idea (Nosé's extended-system trick)"
    To sample NVT deterministically, embed the $N$-particle Hamiltonian system in an *extended* system with one additional degree of freedom $(s,p_s)$. The extended system is Hamiltonian (so admits a symplectic integrator and conserves an exact quantity). The marginal distribution on the physical variables, obtained by integrating out $(s,p_s)$, is *exactly* canonical at temperature $T$. The fictitious variable $s$ rescales atomic momenta and exchanges energy with the physical system; in equilibrium it fluctuates around 1.

**Definition 7.3.1 (Extended Nosé Hamiltonian).** Given a physical system $\{(\mathbf{r}_i,\mathbf{p}_i)\}_{i=1}^N$ with potential $U$, the *Nosé extended Hamiltonian* on the augmented phase space $\{(\mathbf{r}_i,\mathbf{p}_i)\}\cup\{(s,p_s)\}$ is given by (7.33):
$$H_\mathrm{Nosé} = \sum_i \frac{|\mathbf{p}_i|^2}{2m_i s^2} + U(\{\mathbf{r}_i\}) + \frac{p_s^2}{2Q} + g k_BT\,\ln s.$$
The parameter $Q>0$ is the thermostat mass; $g$ is set to $N_\mathrm{df}+1$ in the Nosé scheme (or $N_\mathrm{df}$ in Hoover's reformulation).

**Theorem 7.3.2 (Nosé's theorem).** *Microcanonical sampling of $H_\mathrm{Nosé}$ on the extended phase space induces a canonical distribution at temperature $T$ on the physical phase space.*

*Proof sketch.* The microcanonical density is $\delta(H_\mathrm{Nosé}-E)$. Change momentum variables to $\mathbf{p}^\mathrm{real}=\mathbf{p}/s$ (Jacobian $s^{-3N}$). Integrate the resulting expression over $s$ from $0$ to $\infty$ and over $p_s\in\mathbb{R}$. The $p_s$ integral is Gaussian; the $s$ integral, after the substitution $u = g k_BT\,\ln s$, gives the canonical factor $e^{-(H_\mathrm{phys}-E)/k_BT}$ provided the choice $g = N_\mathrm{df}+1$ in $d=3$. The result is a Boltzmann distribution over $(\mathbf{r},\mathbf{p}^\mathrm{real})$ at temperature $T$. $\square$

The right way to sample the canonical ensemble deterministically is to add a single extra degree of freedom that exchanges energy with the physical system. Nosé's 1984 construction introduces a fictitious variable $s$ with conjugate momentum $p_s$ and writes the extended Hamiltonian

$$
H_\mathrm{Nos\'e} = \sum_i \frac{|\mathbf{p}_i|^2}{2 m_i s^2} + U(\{\mathbf{r}_i\}) + \frac{p_s^2}{2 Q} + g\, k_B T \ln s,
\tag{7.33}
$$

where $Q$ is a "mass" for the new variable (with units of energy·time$^2$) and $g$ is chosen as $N_\mathrm{df} + 1$ so that integrating out $s$ gives the canonical distribution in the physical variables. The role of $s$ is to scale momenta: $\mathbf{p}_i^\mathrm{real} = \mathbf{p}_i / s$.

!!! note "Why this works — Nosé's microcanonical-to-canonical trick"
    The remarkable claim is that integrating microcanonically (NVE) over the *extended* phase space $(\mathbf{r}, \mathbf{p}, s, p_s)$ produces the *canonical* distribution in the physical subspace $(\mathbf{r}, \mathbf{p}^\mathrm{real})$. The proof is short. The microcanonical distribution on the extended system is $\delta(H_\mathrm{Nosé} - E)$. Change variables from $\mathbf{p}$ to $\mathbf{p}^\mathrm{real} = \mathbf{p}/s$; the Jacobian is $s^{-3N}$. Integrate over $s$ and $p_s$:
    $$\int ds\,dp_s\,s^{-3N}\,\delta(H_\mathrm{Nosé} - E).$$
    The Gaussian over $p_s$ gives $\sqrt{2\pi Q}$. The remaining $s$-integral, with the change of variables $u = g k_B T\,\ln s$ that converts the log term into $u$ and the $s^{-3N}$ Jacobian into $e^{-u\cdot 3N/(gk_BT)}$, evaluates to a constant times $e^{-(H_\mathrm{phys}-E)/k_BT}$ — provided $g = 3N + 1$ (in $d=3$). So the marginal distribution on the physical variables is exactly canonical.
    
    The role of $s$ is to *expand* the energy shell of the extended system into the full canonical ensemble of the physical system: $s>1$ regions correspond to lower-than-target physical kinetic energy, $s<1$ to higher. By exploring $s$ microcanonically, the system effectively visits all physical kinetic energies with the right Boltzmann weight.

#### Equations of motion in $s$-variables

Hamilton's equations for the Nosé Hamiltonian (7.33) are
$$\dot{\mathbf{r}}_i = \frac{\partial H_\mathrm{Nosé}}{\partial \mathbf{p}_i} = \frac{\mathbf{p}_i}{m_i s^2},$$
$$\dot{\mathbf{p}}_i = -\frac{\partial H_\mathrm{Nosé}}{\partial \mathbf{r}_i} = -\nabla_i U,$$
$$\dot{s} = \frac{\partial H_\mathrm{Nosé}}{\partial p_s} = \frac{p_s}{Q},$$
$$\dot{p}_s = -\frac{\partial H_\mathrm{Nosé}}{\partial s} = \sum_i \frac{|\mathbf{p}_i|^2}{m_i s^3} - \frac{g k_BT}{s}.$$
Note that in this $s$-time formulation, time is "scaled": the "physical" time is $t' = \int s^{-1}\,dt$. This is awkward to work with directly, hence Hoover's reformulation below.

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

#### Deriving Hoover's $\zeta$-equations from Nosé

The change of variables from $(s, p_s)$ to $(\zeta, \cdot)$ proceeds as follows. Define $\zeta = p_s/(Q)$ and switch to physical-time $t' = \int s^{-1}\,dt$ (so $dt'/dt = 1/s$). Under this transformation:

- Physical position: $\mathbf{r}_i' = \mathbf{r}_i$ (no change).
- Physical momentum: $\mathbf{p}_i' = \mathbf{p}_i/s$ (the scaling).
- Time derivative: $d/dt' = s\,d/dt$.

For positions:
$$\frac{d\mathbf{r}_i'}{dt'} = s\,\frac{d\mathbf{r}_i}{dt} = s\cdot\frac{\mathbf{p}_i}{m_i s^2} = \frac{\mathbf{p}_i'}{m_i},$$
which is just $\dot{\mathbf{r}}_i = \mathbf{p}_i/m_i$ in primed variables: equation (7.34).

For momenta:
$$\frac{d\mathbf{p}_i'}{dt'} = s\,\frac{d}{dt}\!\left(\frac{\mathbf{p}_i}{s}\right) = s\!\left[\frac{\dot{\mathbf{p}}_i}{s} - \frac{\mathbf{p}_i\,\dot{s}}{s^2}\right] = \dot{\mathbf{p}}_i - \frac{\dot{s}}{s}\mathbf{p}_i = -\nabla_i U - \zeta\,\mathbf{p}_i',$$
using $\dot{s}/s = p_s/Q = \zeta$ and $\dot{\mathbf{p}}_i = -\nabla_i U$ (and absorbing the $s$ rescaling of $\mathbf{p}_i$ into the prime). This is (7.35).

For the thermostat:
$$\frac{d\zeta}{dt'} = s\,\frac{d}{dt}\!\left(\frac{p_s}{Q}\right) = \frac{s\,\dot{p}_s}{Q} = \frac{1}{Q}\!\left[\sum_i\frac{|\mathbf{p}_i|^2}{m_i s^2} - g k_BT\right] = \frac{1}{Q}\!\left[\sum_i \frac{|\mathbf{p}_i'|^2}{m_i} - g k_BT\right],$$
which is (7.36). Hoover's reformulation is thus an exact rewriting of Nosé's original equations in physical time and physical momenta.

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

#### Nosé-Hoover chains (Martyna-Tuckerman) — equations

A chain of $M$ thermostats introduces friction variables $\zeta_1, \zeta_2, \ldots, \zeta_M$ with masses $Q_1, Q_2, \ldots, Q_M$. The equations of motion are
$$\dot{\mathbf{r}}_i = \mathbf{p}_i/m_i,$$
$$\dot{\mathbf{p}}_i = -\nabla_i U - \zeta_1\,\mathbf{p}_i,$$
$$\dot{\zeta}_1 = \frac{1}{Q_1}\!\left(\sum_i \frac{|\mathbf{p}_i|^2}{m_i} - gk_BT\right) - \zeta_2\,\zeta_1,$$
$$\dot{\zeta}_k = \frac{Q_{k-1}\zeta_{k-1}^2 - k_BT}{Q_k} - \zeta_{k+1}\,\zeta_k, \qquad k = 2, \ldots, M-1,$$
$$\dot{\zeta}_M = \frac{Q_{M-1}\zeta_{M-1}^2 - k_BT}{Q_M}.$$

The first thermostat $\zeta_1$ couples to the *physical* kinetic energy; each subsequent $\zeta_k$ acts as a heat bath for the previous one. The last thermostat $\zeta_M$ has no further bath and obeys a pure thermalisation equation.

**Why does this help?** A single NH thermostat ($M=1$) on a 1D harmonic oscillator has a *third* conserved quantity beyond energy and momentum (related to the Hamiltonian-like (7.37)), which restricts the trajectory to a 3-torus in 4-dimensional phase space rather than filling the energy shell. The system is not ergodic; histogramming $\mathbf{p}^2$ over a long trajectory does *not* converge to $\chi^2_1(T)$. Adding a second thermostat $\zeta_2$ couples $\zeta_1$ to its own bath and breaks the extra conservation law. With $M\ge 2$, ergodicity is restored and the canonical distribution emerges.

For a large many-atom system, $M=1$ is often "ergodic enough" because the atomic motions themselves provide effective thermal mixing; but the chain is cheap and robust, so production codes default to it. Typical mass choices: $Q_1 = gk_BT\,\tau_T^2$ as in (7.38), $Q_k = k_BT\,\tau_T^2$ for $k\ge 2$.

!!! note "Nosé-Hoover for an isolated oscillator"
    If you run a single 1D oscillator with NH and a single thermostat, the kinetic-energy distribution will not be canonical — you can verify this by histogram. Switch to a 3-chain, and the histogram converges to the correct $\chi^2$ form. This is the textbook test of an MD code's NVT implementation.

**Definition 7.3.3 (Nosé-Hoover chain).** A *Nosé-Hoover chain* of length $M$ is a set of $M$ friction variables $\zeta_1,\ldots,\zeta_M$ with masses $Q_1,\ldots,Q_M$ obeying the equations of motion given immediately after this definition (in the main text). The chain extension is due to Martyna, Klein, and Tuckerman (1992) and restores ergodicity for systems with few degrees of freedom.

**Remark 7.3.4 (Shuichi Nosé, 1984).** Nosé's original Hamiltonian formulation appeared in *Molecular Physics* (1984); Hoover's reformulation in physical-time variables in *Phys. Rev. A* (1985). The extension to chains by Martyna-Klein-Tuckerman (1992) addressed the long-standing ergodicity problem of small-system Nosé-Hoover. Together these three papers define the modern deterministic NVT methodology.

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

!!! note "Why does the Langevin equation sample canonical?"
    The proof is a one-line application of the Fokker-Planck equation associated with (7.39). The stochastic differential equation $m\ddot x = -\partial_x U - \gamma m\dot x + \eta$ with $\langle\eta(t)\eta(t')\rangle = 2m\gamma k_BT\,\delta(t-t')$ implies the Fokker-Planck equation
    $$\frac{\partial P}{\partial t} = -\frac{\partial}{\partial x}\!\left(\frac{p}{m}P\right) + \frac{\partial}{\partial p}\!\left[\partial_x U\,P + \gamma p\,P + \gamma m k_BT\,\frac{\partial P}{\partial p}\right]$$
    on the joint distribution $P(x,p,t)$. Substituting the canonical guess $P_\mathrm{can}(x,p) \propto e^{-\beta H(x,p)}$ with $H = p^2/(2m) + U$:
    $$\partial_x P_\mathrm{can} = -\beta\,\partial_x U\,P_\mathrm{can}, \quad \partial_p P_\mathrm{can} = -\beta\,(p/m)\,P_\mathrm{can}.$$
    The Liouville part vanishes (Hamiltonian flow is incompressible on the canonical density). The diffusion-drift part becomes
    $$\partial_p\!\left[\partial_x U\,P_\mathrm{can} + \gamma p\,P_\mathrm{can} - \gamma m k_BT\,\beta(p/m)\,P_\mathrm{can}\right] = \partial_p[\partial_x U\,P_\mathrm{can}] = 0$$
    because the $\gamma p$ and $-\gamma m k_BT\,\beta(p/m)$ pieces cancel exactly thanks to the fluctuation-dissipation relation (7.40). $P_\mathrm{can}$ is a stationary solution; ergodicity (provable for any $\gamma > 0$ and sufficiently regular $U$) ensures it is the *unique* stationary solution. Langevin samples canonical.
    
    The drift-fluctuation pairing in (7.40) is not arbitrary — it is the unique noise amplitude that produces the correct Boltzmann weight. Halving $\gamma$ without halving the noise variance gives a non-canonical stationary distribution; doubling them both gives the same canonical distribution but on a faster relaxation timescale.

Langevin dynamics samples NVT correctly for any $\gamma > 0$. Small $\gamma$ gives weak coupling and Newtonian-like dynamics on short timescales; large $\gamma$ overdamps the dynamics into Brownian motion. For most condensed-phase systems $\gamma \sim 1$–10 ps$^{-1}$ is a good range. Note that diffusion coefficients and other dynamical observables **are affected** by $\gamma$ — Langevin distorts kinetics — whereas static thermodynamic averages are correct for any $\gamma$.

If you need correct kinetics, use Nosé-Hoover; if you only need correct static averages and you want robust ergodicity, use Langevin.

### Numerical integration of Langevin: the BAOAB splitting

A subtle practical point: naive substitution of $\boldsymbol{\eta}$ into a Verlet step gives an integrator whose accuracy is only $O(\Delta t)$ — strong-order convergence rather than the $O(\Delta t^2)$ you would expect. The remedy is operator splitting. Leimkuhler-Matthews (2013) proposed **BAOAB**, in which one Langevin step decomposes into five sub-steps in the order: half-momentum-kick (B), half-position-update (A), Ornstein-Uhlenbeck thermalisation (O), half-position-update (A), half-momentum-kick (B). Concretely:

1. $\mathbf{p}^{n+1/4} = \mathbf{p}^n + \tfrac{\Delta t}{2}\,\mathbf{F}(\mathbf{r}^n)$. (B)
2. $\mathbf{r}^{n+1/2} = \mathbf{r}^n + \tfrac{\Delta t}{2}\,\mathbf{p}^{n+1/4}/m$. (A)
3. $\mathbf{p}^{n+3/4} = c\,\mathbf{p}^{n+1/4} + \sqrt{m k_BT(1-c^2)}\,\boldsymbol{\xi}^n$, with $c=e^{-\gamma\Delta t}$. (O)
4. $\mathbf{r}^{n+1} = \mathbf{r}^{n+1/2} + \tfrac{\Delta t}{2}\,\mathbf{p}^{n+3/4}/m$. (A)
5. $\mathbf{p}^{n+1} = \mathbf{p}^{n+3/4} + \tfrac{\Delta t}{2}\,\mathbf{F}(\mathbf{r}^{n+1})$. (B)

Here $\boldsymbol{\xi}^n \sim \mathcal{N}(0, 1)$. The "O" step is exact for the Ornstein-Uhlenbeck process $\dot{\mathbf{p}} = -\gamma\mathbf{p} + \eta$, which is why BAOAB recovers second-order weak accuracy and produces correct configurational averages even at moderately large $\Delta t$. LAMMPS' `fix langevin` and OpenMM's `LangevinMiddleIntegrator` use BAOAB or close relatives.

**Example 7.3.5 (Verifying NVT sampling on a 1D harmonic oscillator).**
*Problem.* Run a Langevin MD of a single 1D harmonic oscillator at temperature $T$. Verify that the kinetic-energy histogram matches the canonical distribution.

*Solution.* Integrate $\ddot x = -x - \gamma\dot x + \eta(t)$ for $10^6$ steps with $\Delta t = 0.01$, $\gamma = 1$, and noise $\eta$ scaled so $\langle\eta^2\rangle = 2\gamma k_BT/\Delta t$. Record $v(t)$ at each step. Histogram $\tfrac{1}{2}m v^2$; compare to the exponential $P(KE) \propto e^{-KE/k_BT}$.

*Discussion.* For Langevin, the agreement should be quantitative within shot-noise error bars. Repeating with Nosé-Hoover and $M=1$ chain length will give a *non*-canonical histogram for this single-DOF system — the kinetic energy gets stuck on certain values. Increasing to $M=3$ restores the canonical distribution. This is the standard verification test of an MD code's NVT implementation.

**Example 7.3.6 (Choosing the thermostat relaxation time $\tau_T$).**
*Problem.* For an MD of liquid water at 300 K with $N = 1000$ molecules, what $\tau_T$ should you choose for the Nosé-Hoover thermostat?

*Solution.* Two constraints:

1. $\tau_T \gg \Delta t$: the thermostat must be slow compared to integration. With $\Delta t = 1$ fs, $\tau_T > 10$ fs.
2. $\tau_T \gg \tau_\mathrm{coupling}$: the thermostat must be slow compared to the system's internal energy-transfer timescale, so it does not perturb dynamics. For water at 300 K, the dominant internal coupling is via the H-bond network with timescale $\sim 1$ ps.

A standard choice: $\tau_T = 0.1$ to $1$ ps. Too short and you suppress kinetic energy fluctuations relative to canonical. Too long and equilibration takes forever. The default of $\tau_T = 100\Delta t = 100$ fs in LAMMPS' `fix nvt` is a safe middle ground.

*Discussion.* For static thermodynamic averages, the choice of $\tau_T$ in the canonical-correct thermostats (NH, Langevin, CSVR) does not affect the final distribution — only the relaxation timescale. For *dynamic* properties (diffusion, viscosity, IR spectra), $\tau_T$ matters: too-strong coupling distorts the dynamics. If you care about kinetics, use NH with moderate $\tau_T$ or Langevin with small $\gamma$.

## Parrinello-Rahman barostat: NPT

!!! tip "Why we need a barostat as well"
    Constant-pressure rather than constant-volume ensembles are what nature gives us in most cases — a beaker of liquid is at atmospheric pressure, not fixed volume; a crystal in a heat bath expands as $T$ rises. To compute thermal expansion coefficients, melting points, equation-of-state behaviour, you need to let the cell volume *respond* to the internal stress. This requires upgrading the cell vectors from fixed parameters to dynamical variables — what Parrinello and Rahman did in 1980.
    
    The picture: the simulation cell is no longer a rigid box but a soft membrane with its own kinetic energy and stiffness. External pressure pushes inward; internal stress pushes outward; the cell oscillates until they balance. In addition, anisotropic stresses can shear or stretch the cell, allowing simulation of crystal phase transitions where the unit cell shape itself changes (martensite, ferroelectric switching).

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

#### Parrinello-Rahman: physical picture

Parrinello-Rahman dynamics treats the cell vectors $\mathbf{h}=(\mathbf{a}_1|\mathbf{a}_2|\mathbf{a}_3)$ as additional generalised coordinates. Each atom's position is decomposed as $\mathbf{r}_i = \mathbf{h}\,\mathbf{s}_i$ with $\mathbf{s}_i \in [0,1)^3$. The cell shape $\mathbf{h}$ has its own kinetic energy $\tfrac{1}{2}W\,\mathrm{Tr}(\dot{\mathbf{h}}^\top\dot{\mathbf{h}})$ and a potential cost $P_\mathrm{ext}V$ from the external pressure pushing on the boundary.

Two consequences:

1. **The cell can deform anisotropically.** Under an applied uniaxial stress, the cell stretches in the loading direction; under shear stress, it shears. This is essential for studying mechanical response, phase transitions with structural symmetry change, and martensitic transformations.
2. **Cell modes have their own dynamics.** The cell oscillates at frequencies $\omega_\mathrm{cell} \sim \sqrt{(B+G)/W}$ where $B,G$ are bulk and shear moduli. These oscillations should be slower than atomic vibrations but faster than the equilibration time; choosing $W$ such that $\tau_P = 2\pi/\omega_\mathrm{cell}$ falls in 1-10 ps is standard.

The Martyna-Tobias-Klein (MTK) scheme combines Parrinello-Rahman cell dynamics with a Nosé-Hoover chain on both atomic and cell degrees of freedom. The conserved quantity has eight or nine terms (atomic KE, atomic PE, cell KE, $PV$, two thermostat KEs, two thermostat potentials, and the cell-thermostat-thermostat coupling); its constancy is a diagnostic of correct MTK integration. The full equations of motion fill a page; you do not implement them by hand. LAMMPS' `fix npt` is the reference implementation.

**Choice of $W$.** The cell mass controls how fast the box responds. Too small $W$: the box rings at unphysical frequencies. Too large $W$: the box can't follow density changes during equilibration. A standard scaling is

$$
W = (N_\mathrm{df} + 3) k_B T\, \tau_P^2,
\tag{7.43}
$$

with $\tau_P \approx 0.5$–5 ps. As with $\tau_T$, defaults are usually fine.

**Isotropic vs anisotropic cells.** For liquids and cubic solids, isotropic pressure coupling (`fix npt iso`) constrains the cell to remain a scaled copy of itself, varying only the volume. For non-cubic crystals, anisotropic coupling (`fix npt aniso`) allows the three cell edges to fluctuate independently while keeping the angles fixed; full triclinic (`fix npt tri`) lets shear strains develop as well. Use the most restrictive setting consistent with your symmetry — extra cell degrees of freedom add noise without information.

!!! warning "NPT during phase transitions"
    Across a first-order transition the equilibrium volume jumps discontinuously. An NPT simulation crossing the transition will not give a sharp jump in the volume — it will give a slow drift over many nanoseconds, because the new phase has to nucleate. Don't mistake the slow drift for "not equilibrated yet" and abandon the simulation. The two-phase method ([§8.4](../ch08-statmech/04-phase-diagrams.md)) is the proper tool here.

**Definition 7.3.7 (Stress tensor).** The instantaneous stress tensor of an $N$-atom periodic system at volume $V$ is
$$\sigma_{\alpha\beta} = -\frac{1}{V}\left[\sum_i m_i v_{i\alpha} v_{i\beta} - \sum_{i<j} r_{ij,\alpha} F_{ij,\beta}\right],$$
where the first term is the kinetic contribution (ideal-gas pressure) and the second is the virial. The instantaneous pressure is $P = -\tfrac{1}{3}\mathrm{Tr}(\boldsymbol{\sigma})$.

**Remark 7.3.8 (MTK and the conserved energy).** The Martyna-Tobias-Klein NPT scheme combines Parrinello-Rahman cell dynamics with Nosé-Hoover chains on both atomic and barostat degrees of freedom. The conserved quantity has nine terms; its constancy to $10^{-5}$ of its initial value is the standard diagnostic that the integrator is correct. LAMMPS reports this as `c_thermo_pe + c_thermo_ke + cell-KE + PV + thermostat-terms`.

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

### Section summary

The thermostat hierarchy:

1. **NVE** is what Verlet gives you. Use for benchmarking integrators and for energy-flow studies; almost never the right ensemble for material properties.
2. **Velocity rescaling and Berendsen** are equilibration tools only. They do not sample canonical; they fail at fluctuation-derived properties like heat capacity.
3. **Stochastic velocity rescaling (CSVR)** is the simple-and-correct production thermostat: as cheap as Berendsen, samples canonical exactly, minimal kinetic distortion.
4. **Nosé-Hoover with chain length $\ge 2$** is the deterministic production thermostat. Best preservation of kinetics; chain extension restores ergodicity for small systems.
5. **Langevin** is the stochastic production thermostat. Robust ergodicity, canonical for any $\gamma > 0$, but distorts kinetics — choose $\gamma$ for the application.
6. **Parrinello-Rahman + MTK** is the NPT analogue, cell vectors as dynamical variables, full anisotropic cell response.

In LAMMPS:
- `fix nve` — NVE.
- `fix temp/csvr T T tau seed` — CSVR.
- `fix nvt temp T T tau` — Nosé-Hoover chain (3 by default).
- `fix langevin T T tau seed` — Langevin (combine with `fix nve` for the integrator).
- `fix npt temp T T tauT iso P P tauP` — MTK.

The right choice depends on the question: dynamics → NH or NVE; statics → any canonical thermostat; pressure control → NPT.

## What we have

**Remark 7.3.9 (Equilibration vs. production).** A common rookie mistake is to use a single thermostat throughout the simulation. The professional pattern is: run a short equilibration (10-100 ps) with a strongly-coupled thermostat (Berendsen, $\tau_T = 0.01-0.1$ ps) to quickly bring the system to target $T$, then switch to a canonical thermostat (NH, Langevin, CSVR) with moderate coupling for production. The equilibration phase samples nothing useful — its only job is to relax the system; the production phase is where averages are computed. Always discard the equilibration trajectory; never include it in any time-average.

We can now sample NVE, NVT, and NPT ensembles, given a force model. We have not specified what produces those forces. That is [§7.4](04-force-fields.md): classical force fields, from the cartoon (Lennard-Jones) to the production (EAM, Tersoff, ReaxFF). The chapter on machine-learning interatomic potentials ([Chapter 9](../ch09-mlip/index.md)) will then take the next step, learning forces directly from DFT.
