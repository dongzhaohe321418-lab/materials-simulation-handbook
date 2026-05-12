# Kinetic Monte Carlo and Phase Field

Atomistic MD struggles with two regimes: very slow processes (diffusion in
solids, where individual hops are rare events on the MD timescale) and very
large length scales (microstructure evolution in a polycrystal). Two distinct
techniques address these gaps without going all the way to continuum
mechanics: kinetic Monte Carlo (KMC) for rare-event dynamics on lattices, and
phase field for microstructure on the mesoscale.

Both methods consume *inputs from atomistic simulation* — barriers,
free-energy curves — and produce *outputs at mesoscopic scales* — diffusion
coefficients, grain morphologies. They are quintessentially multiscale.

## Kinetic Monte Carlo

The problem KMC solves: imagine a vacancy in a copper crystal. At room
temperature the vacancy jumps to a neighbouring site with rate $k = \nu_0
\exp(-E_a/k_B T)$. With $\nu_0 \sim 10^{13}$ Hz and $E_a \sim 0.7$ eV, the
hopping rate is $\sim 10^{1}$ s$^{-1}$. To see one hop in MD you would need
to integrate for $\sim 10^{12}$ MD steps. To see equilibration of a vacancy
concentration profile, ten orders of magnitude more.

MD cannot do this. But the *physics* is simple: a small set of discrete
events, each with a well-defined rate. KMC handles exactly this kind of
problem.

### The KMC algorithm

The setup:

1. **A catalogue of events.** Each event $i$ is a transition between two
   states of the system. For each event we know its rate $k_i$.
2. **An initial state.** Some configuration of atoms on a lattice.
3. **A time-stepping rule.** Advance the system stochastically by selecting
   one event at a time.

The core idea is that in a system with rates $k_1, k_2, \ldots, k_n$, the
time to the *next* event of any type is exponentially distributed with mean
$1/\sum_i k_i$, and the probability that the next event is event $i$ is
$k_i / \sum_j k_j$.

Concretely:

1. Enumerate all possible events from the current state. Compute their rates
   $\{k_i\}$.
2. Compute the total rate $R = \sum_i k_i$.
3. Draw $\xi_1, \xi_2 \sim \mathrm{Uniform}(0, 1)$.
4. Advance time by $\Delta t = -\ln(\xi_1) / R$.
5. Select event $j$ such that $\sum_{i < j} k_i / R \leq \xi_2 < \sum_{i \leq j} k_i / R$.
6. Execute event $j$: update the configuration.
7. Repeat.

This is sometimes called the BKL algorithm (Bortz-Kalos-Lebowitz, 1975) or
the *n*-fold way. It is statistically exact: the trajectory it produces is
a sample from the master equation for the discrete state space.

### Rates from DFT

Where do the rates come from? For a vacancy hop or an adatom diffusion event,
the standard approach is *harmonic transition state theory*:

$$
k = \nu \exp\left(-\frac{E_a}{k_B T}\right)
$$

where $E_a$ is the activation energy (the saddle-point energy minus the
initial-state energy) and $\nu$ is a prefactor related to the vibrational
modes at the initial and transition states.

$E_a$ is computed from DFT using the *nudged elastic band* (NEB) method or
the *climbing image NEB* (CI-NEB). $\nu$ is computed from the Hessian of the
energy at the initial and transition states.

The flow of information is a textbook sequential coupling:

DFT (one calculation per event type) → harmonic TST → KMC rate catalogue →
KMC simulation → diffusion coefficient.

For a simple system with a handful of distinct local environments, a few
dozen DFT calculations may suffice for the rate catalogue. For complex
systems with many local environments (alloys, surface reactions), you may
need hundreds or thousands of NEBs, and machine-learning the barriers
(*kinetic neural networks*) becomes attractive.

### Worked example: vacancy diffusion on a 2D lattice

Let us write a minimal KMC for vacancy diffusion. The system: a $L \times L$
square lattice with periodic boundaries, one vacancy. The vacancy hops to a
nearest neighbour with rate $k = \nu \exp(-E_a/k_B T)$. All four hops have
the same barrier (we are on a Bravais lattice), so the algorithm is
particularly simple.

```python
"""Vacancy diffusion on a 2D square lattice via kinetic Monte Carlo.

Computes the diffusion coefficient by mean-square displacement of a single
vacancy. Demonstrates the BKL algorithm in its simplest form.
"""
from __future__ import annotations

import numpy as np


def run_kmc(
    L: int = 50,
    n_steps: int = 200_000,
    e_a: float = 0.7,        # activation energy in eV
    T: float = 600.0,        # temperature in K
    nu: float = 1.0e13,      # attempt frequency in Hz
    seed: int = 0,
) -> tuple[float, float]:
    """Run KMC for a single vacancy on a 2D periodic lattice.

    Parameters
    ----------
    L : int
        Linear size of the lattice (lattice units).
    n_steps : int
        Number of KMC steps to take.
    e_a, T, nu : float
        Activation energy (eV), temperature (K), attempt frequency (Hz).
    seed : int
        RNG seed.

    Returns
    -------
    diffusion_coefficient : float
        Estimated D in units of lattice-units^2 / s.
    total_time : float
        Total simulated time in s.
    """
    kB = 8.617e-5  # eV / K
    k_hop = nu * np.exp(-e_a / (kB * T))   # rate per direction
    R = 4.0 * k_hop                         # all four hops equally likely

    rng = np.random.default_rng(seed)
    pos = np.array([0, 0], dtype=np.int64)  # vacancy position (unwrapped)
    total_time = 0.0
    msd_sum = 0.0
    n_msd_samples = 0

    moves = np.array([[1, 0], [-1, 0], [0, 1], [0, -1]])

    for _ in range(n_steps):
        # 1. Draw two uniform numbers.
        xi1 = rng.random()
        xi2 = rng.random()
        # 2. Advance time.
        dt = -np.log(xi1) / R
        total_time += dt
        # 3. Select next move (all equally likely here).
        idx = int(xi2 * 4)
        pos = pos + moves[idx]
        # 4. Accumulate squared displacement from origin.
        msd_sum += float(pos[0] ** 2 + pos[1] ** 2)
        n_msd_samples += 1

    # In 2D, <r^2> = 4 D t, so D = <r^2> / (4 t).
    mean_msd = msd_sum / n_msd_samples
    # Approximate mean time per sample by total_time / n_msd_samples.
    mean_t = total_time / n_msd_samples
    D = mean_msd / (4.0 * mean_t) if mean_t > 0 else float("nan")
    return D, total_time


if __name__ == "__main__":
    D, t = run_kmc()
    print(f"Total simulated time: {t:.3e} s")
    print(f"Estimated D:          {D:.3e} lattice^2 / s")
```

A few teaching points from this minimal example.

**The lattice $L$ never enters.** A single isolated vacancy on a periodic
lattice has the same dynamics for any $L$, because all sites are equivalent.
The reason to make $L$ large is to track the mean-square displacement without
unwrapping issues; in this code we keep an *unwrapped* position to avoid
periodic-boundary correction.

**Time really is real.** The output is in *seconds*. There is no "MD
timestep" here. The time step is whatever the exponential clock produces.
Long runs need many KMC steps.

**Generalising to multiple events.** When events have different rates (a
vacancy near a solute hops more slowly than a vacancy in the bulk), you
build a cumulative-rate array and use a binary search in step 5. For very
large rate catalogues, *rejection-free* and *list-based* variants of BKL
keep the time per step manageable.

**Generalising beyond the lattice.** *Off-lattice* KMC allows continuous
positions; the events are catalogued by local environment rather than by
lattice site. This is more general but requires a strategy for identifying
unique events as the system evolves.

### Limitations

KMC inherits the limitations of harmonic TST:

- **Only works for rare events.** If the barriers are low relative to
  $k_B T$, the assumptions of TST break down and you should just run MD.
- **Requires a complete catalogue.** If you forget an event type, it will
  not happen in the simulation, and you will not know.
- **Pre-tabulated rates.** If the barrier depends on the configuration
  (which it always does, weakly), you need either a many-rate catalogue or
  a machine-learning surrogate.

The first two are why KMC is often paired with *automatic event detection*
algorithms like ART-n or kART, which discover new event types as the
simulation runs.

## Phase field

Where KMC handles atomic-scale rare events, phase field handles the
*continuum* mesoscale. The protagonist is an *order parameter*
$\phi(\mathbf{r}, t)$, a smooth field defined over the simulation domain,
that distinguishes one phase from another. For solidification, $\phi$ might
be 1 in the solid and 0 in the liquid, with a smooth transition region of
finite width. For grain growth, multiple order parameters $\phi_\alpha$
distinguish grains of different orientation.

The phase-field method *does not* try to resolve atoms; it tries to resolve
*interfaces*. Where the interface is, where it moves, how it interacts with
other interfaces and with the bulk thermodynamics.

### Cahn-Hilliard: conserved order parameter

For a conserved order parameter — typically a composition variable, where
the total amount is fixed — the governing equation is *Cahn-Hilliard*:

$$
\frac{\partial \phi}{\partial t} = \nabla \cdot \left( M \nabla \frac{\delta F}{\delta \phi} \right)
$$

with

$$
F[\phi] = \int \left[ f(\phi) + \frac{\kappa}{2} |\nabla \phi|^2 \right] d^3 r
$$

where $f(\phi)$ is the bulk free energy density, $\kappa$ is the gradient
energy coefficient (penalising sharp interfaces), and $M$ is a mobility
parameter.

The Cahn-Hilliard equation is fourth-order in space (the $\nabla^2$ of a
$\nabla^2$), which makes it numerically more demanding than a standard
diffusion equation but well within reach of finite differences or spectral
methods.

It is the canonical equation for *spinodal decomposition*: a homogeneous
mixture that lowers its free energy by separating into phases conserves the
total composition, so the order parameter (local composition) is conserved.

### Allen-Cahn: non-conserved order parameter

For a non-conserved order parameter — typically a structural variable, such
as a measure of how solid-like a region is, where the total can change —
the governing equation is *Allen-Cahn*:

$$
\frac{\partial \phi}{\partial t} = -L \frac{\delta F}{\delta \phi}
$$

with the same form of free-energy functional $F[\phi]$ as above, and $L$
a kinetic coefficient.

This is second-order in space and much like a relaxational diffusion
equation. It describes *interface motion* under the thermodynamic driving
force; flat interfaces with no driving force are stationary, curved
interfaces move under capillarity (the Gibbs-Thomson effect), and
interfaces in the presence of a free-energy bias move with a velocity
proportional to that bias.

In practice, real materials simulations often use *coupled* Cahn-Hilliard
and Allen-Cahn equations: a conserved composition variable plus one or more
non-conserved structural variables (one per phase or per grain
orientation).

### Inputs from atomistic simulation

A phase-field simulation needs two kinds of input that are not internal to
the continuum theory:

1. **The bulk free energy $f(\phi)$.** For binary alloy decomposition,
   this comes from CALPHAD thermodynamic databases or from atomistic free
   energy calculations (thermodynamic integration, free-energy MD). For
   solidification, $f$ is parameterised to give the correct melting point
   and latent heat.

2. **The gradient energy coefficient $\kappa$.** This sets the interfacial
   energy $\gamma$ and the interface width $w$. In equilibrium:

   $$
   \gamma = \sqrt{2 \kappa\, f_\mathrm{barrier}}, \quad w \sim \sqrt{\kappa / f_\mathrm{barrier}}
   $$

   where $f_\mathrm{barrier}$ is the height of the double-well in $f(\phi)$.
   The interface energy can be computed from atomistic simulation by
   constructing a coherent two-phase interface and measuring the excess
   free energy.

The kinetic coefficients $M$ and $L$ are calibrated from atomistic mobility
data (diffusion coefficients, interface velocities under known driving
forces).

This is again a sequential pipeline: atomistic → free energies and
interfacial properties → phase field → microstructure evolution → mechanical
properties.

### What phase-field captures

Phase field is the workhorse of mesoscale microstructure simulation. It
handles:

- **Solidification of alloys** (dendrite morphology, microsegregation).
- **Solid-state phase transformations** (martensitic transformations,
  precipitation).
- **Grain growth and recrystallisation** in polycrystals.
- **Microstructure under stress and strain** (with coupling to elasticity).

It does *not* capture:

- Atomic-scale chemistry (defects, dislocations cores).
- Anything sharper than the imposed interface width.
- True nucleation events (these must be added explicitly, often by
  classical nucleation theory or by inserting nuclei stochastically).

### Connecting phase field upward to engineering

Phase-field outputs are typically *microstructures*: spatial fields of order
parameters, compositions, and stresses. To predict engineering properties
(strength, toughness, conductivity), you usually feed the microstructure
into:

- A *crystal-plasticity finite-element model* (CP-FEM) for mechanical
  response — see [Section 5](05-finite-element-bridge.md).
- An *effective-medium theory* or *full-field FFT solver* for transport
  properties.

This is the third sequential coupling in the chain: DFT → atomistic →
phase field → FEM → engineering.

## When to use which

A quick decision guide.

- **KMC** when you have a small number of discrete event types with known
  barriers, the dynamics is rare-event dominated, and you care about
  composition or population statistics over long times. Classic examples:
  diffusion in solids, surface catalysis, thin-film growth, defect kinetics
  in irradiated materials.

- **Phase field** when you have a continuum order parameter that varies
  smoothly in space, you care about morphology and interface dynamics on
  μm-mm scales over s-h timescales, and you can afford a free-energy
  functional. Classic examples: solidification microstructure, spinodal
  decomposition, grain growth.

- **Neither** when the system is so out of equilibrium that there is no
  well-defined free energy surface (driven systems far from equilibrium,
  glassy systems near $T_g$), or when atomistic detail is essential to
  the answer.

The shared feature of KMC and phase field is that they discard atomic-scale
information in favour of a coarser description that nonetheless retains the
right thermodynamics and kinetics. They are the *de facto* tools for the
mesoscale gap, and modern multiscale workflows often run several of them in
parallel, exchanging data with atomistic simulation.

!!! tip "Sanity-checking mesoscale outputs"
    Always look at the output of a KMC or phase-field simulation as a
    movie, not as a number. The morphology — dendrite shapes, grain
    distributions, segregation patterns — is the immediate sanity check. If
    your dendrites look like the wrong material, your inputs are wrong.
    Single-number outputs (mean grain size, diffusion coefficient) can hide
    a lot.

## Looking ahead

We have now seen most of the bridging methods that take atomistic data and
turn it into mesoscale predictions. The last step on the multiscale staircase
is the engineering scale — the finite element method, which takes mesoscale
material properties and predicts the behaviour of components and structures.
That is the subject of [the next section](05-finite-element-bridge.md).
