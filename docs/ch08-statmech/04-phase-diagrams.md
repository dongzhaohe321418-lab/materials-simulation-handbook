# 8.4 Phase Diagrams from Simulation

A phase diagram is a map: which phase is thermodynamically stable as a function of temperature, pressure, composition, and any other relevant intensive variables. For an unalloyed element like copper, the relevant axes are $T$ and $P$ and the phase diagram is a 2D plot with fcc, hcp, bcc, and liquid regions. For a binary alloy, composition adds a third axis. For an oxide under chemical potential control, the oxygen chemical potential adds another. Phase diagrams quickly become high-dimensional and computationally expensive.

Two complementary approaches to constructing them from atomistic simulation are presented here.

## Free energy minimisation across phases

At fixed $(T, P)$, the equilibrium phase is the one with the lowest Gibbs free energy. To locate a phase boundary between phases A and B in $(T, P)$ space, compute $G_A(T, P)$ and $G_B(T, P)$ and find the locus where they cross:

$$
G_A(T, P) = G_B(T, P).
\tag{8.41}
$$

This is the cleanest, most general approach. Its main difficulty is that **absolute** Gibbs free energies cannot be measured directly — only differences. Three classical decompositions help:

### Reference-state thermodynamic integration

Compute $G(T, P)$ along a path from a reference state where it is known analytically:

$$
G(T, P) = G_\mathrm{ref} + \int_\mathrm{path} dG.
\tag{8.42}
$$

Along an isobar at fixed $P$, $dG = -S\, dT$; along an isotherm at fixed $T$, $dG = V\, dP$. Both can be evaluated from MD averages: $S$ via the thermodynamic integration of $C_P/T$ from a reference $T_0$, $V$ directly from NPT.

For solids, a common reference is the **Einstein crystal**: a system of independent harmonic oscillators with the same equilibrium structure but no inter-atomic forces beyond a per-atom harmonic well. The Einstein crystal has an analytical free energy. The Frenkel-Ladd path turns on the LJ (or real) interactions gradually from this reference, accumulating $\Delta G$ via thermodynamic integration. The result is the **absolute** free energy of the crystalline phase.

For liquids, the analogous reference is an **ideal gas** at the same temperature; the path is the Widom particle insertion or a TI of (8.27)-style.

This sounds elaborate. It is — Frenkel-Ladd for a single solid phase at a single state point can require a few days of MD. For a small region of phase diagram and a few phases it is the rigorous gold standard.

### Equation-of-state fitting

Less rigorous but often sufficient: compute $E(V, T)$ on a grid and fit a thermodynamic model (Birch-Murnaghan, Vinet) to extract $V(P, T)$ and $G(T, P)$ semi-analytically. This requires only NVT simulations at multiple volumes and temperatures, no special TI paths.

## Coexistence simulations: the two-phase method

A more direct approach to locating a phase boundary, especially for melting: simulate the two phases **simultaneously** in a single box and watch which one grows.

**Two-phase setup.** Build a supercell with the two phases stacked along $z$ — top half liquid, bottom half solid — with a planar interface. Equilibrate at constant $T$ and $P$ (or constant $T$ and a fixed total volume that allows the interface to move).

**Three possible outcomes:**

1. $T > T_m$: the solid melts, and the simulation ends with the box entirely liquid.
2. $T < T_m$: the liquid crystallises, and the simulation ends entirely solid.
3. $T = T_m$: the interface is stationary (within fluctuations). Both phases coexist indefinitely.

Bracket $T_m$ by running at several temperatures and noting which way each interface moves.

The two-phase method avoids the nucleation-barrier problem that biases single-phase heating/cooling simulations: with an interface already present, no nucleation is needed for the phase transition to proceed. The result is the true thermodynamic $T_m$ for the chosen force field, free of kinetic hysteresis.

**Implementation in LAMMPS.** A typical sequence:

```text
# Build single crystal of Cu
lattice         fcc 3.615
region          box block 0 6 0 6 0 24
create_box      1 box
create_atoms    1 box

pair_style      eam/alloy
pair_coeff      * * Cu_mishin1.eam.alloy Cu

# Heat the top half far above Tm to melt it
region          top block INF INF INF INF 12 24
group           top region top
group           bot subtract all top

# Stage 1: melt the top half quickly
fix             1 top nvt temp 300 3000 0.05
fix             2 bot nvt temp 300 300 0.05
run             20000
unfix           1; unfix 2

# Stage 2: equilibrate both halves at target T (close to expected Tm)
fix             3 all nvt temp 1325 1325 0.1
run             50000

# Stage 3: production, watch which way the interface moves
fix             4 all npt temp 1325 1325 0.1 z 1.0 1.0 1.0
dump            1 all custom 500 cu_coex.lammpstrj id type x y z
thermo          200
thermo_style    custom step temp pe ke etotal press vol density
run             1000000
```

Monitor either the energy (slowly drifting up means melting won, down means crystallising won) or compute a structural order parameter — bond orientational order $Q_6$ separates liquid from crystal cleanly — and track the fraction of crystalline atoms over time.

After bracketing $T_m$ to within a few K by running at 1300, 1325, and 1350 K, you can refine further or trust the bracket. For Mishin Cu EAM, this gives $T_m \approx 1320$ K, consistent with rigorous TI calculations on the same potential.

## Melting from MD: subtleties

A few things to watch for:

**Premelting at surfaces.** Free surfaces of solids melt below their bulk melting temperature. If your supercell has free surfaces, you will measure a surface premelting temperature, not bulk $T_m$. Use periodic boundary conditions in all three directions, with the two-phase interface as the only inhomogeneity.

**Anisotropic systems.** For non-cubic crystals, the cell vectors are not equivalent, and the choice of interface orientation matters. Different crystal faces have different surface energies and may give slightly different $T_m$ values via the two-phase method (because of finite-size interfacial contributions). Use the lowest-energy face (typically the close-packed plane) for the most accurate result.

**Pressure equilibration.** The volume of the liquid is typically a few percent larger than the solid at the same conditions. NPT must allow this volume change. If you use isotropic NPT on an interface aligned along $z$, the box will both expand in $z$ (to accommodate the liquid) and shrink slightly in $x$ and $y$, distorting the crystal. Use anisotropic NPT (`fix npt aniso` with independent control over each cell length) or, better, only allow $z$ to fluctuate (`fix npt z 1.0 1.0 1.0`) while $x$ and $y$ are pinned at the equilibrium solid values.

**Finite-size effects.** Two-phase $T_m$ depends weakly on the supercell cross-section, with a $\sim 1/L^2$ scaling from the interface free energy. For Cu, a $6 \times 6$ cross-section gives $T_m$ accurate to $\sim 5$ K; smaller boxes overestimate.

## Solid-solid phase boundaries

For solid-solid transitions (fcc to hcp in cobalt, $\alpha$ to $\gamma$ in iron), the two-phase method is harder because solid-solid interfaces are slow to migrate. Free energy minimisation via reference-state TI is then the cleaner route. Alternatively, **non-equilibrium free energy methods** like the Jarzynski equality

$$
e^{-\beta \Delta A} = \langle e^{-\beta W}\rangle_\mathrm{neq},
\tag{8.43}
$$

let you compute free energy differences from finite-rate switching trajectories, with cost roughly comparable to TI but a different statistical noise profile.

For multicomponent systems with chemical-ordering transitions (e.g., the L1$_2$-disordered transition in Cu$_3$Au), specialised methods like the **semi-grand canonical** ensemble or **cluster expansion + Monte Carlo** are standard. These are outside the scope here; the Frenkel-Smit textbook is the reference.

## Connection to CALPHAD databases

The CALPHAD (Calculation of Phase Diagrams) framework, developed since the 1970s, builds phase diagrams from **assessed** thermodynamic functions — polynomial fits to experimental data, occasionally augmented by first-principles or MD inputs. CALPHAD databases like TCFE for steels, COST507 for aluminium alloys, and the SGTE pure-element data are the practical workhorses of industrial metallurgy.

When should you use simulation rather than CALPHAD?

| Use simulation when | Use CALPHAD when |
|---|---|
| The system is outside CALPHAD coverage (novel chemistry) | Mature alloy systems with extensive experimental data |
| You need atomistic mechanism, not just thermodynamics | Industrial design where only the equilibrium answer matters |
| Extreme conditions (high P, far from ambient) | Standard conditions, well-assessed databases |
| Sub-grid information for multiscale modelling | Macroscopic process simulation |
| Comparing trends across families of materials | Specific industrial alloy specifications |

Modern practice integrates both: high-throughput DFT + MLIP MD provides input to CALPHAD assessments, especially for finite-temperature solid free energies that were historically the weakest link. The **OQMD**, **Materials Project**, and **AFLOW** databases publish DFT-derived formation enthalpies that feed CALPHAD; the next generation of databases will include vibrational free energies from MLIP MD.

## Clausius-Clapeyron, nucleation barriers, and the hysteresis method

Two complementary ideas tie the previous subsections together: the *slope* of a phase boundary is fixed by a thermodynamic identity that any simulation must respect, and the *direct measurement* of melting from heating a single phase is systematically biased by a nucleation barrier that the two-phase method bypasses.

### Deriving Clausius-Clapeyron

Along a phase boundary between phases A and B, $G_A(T, P) = G_B(T, P)$. Differentiate this identity along the boundary, where both phases coexist:

$$
\left(\frac{\partial G_A}{\partial T}\right)_P\, dT + \left(\frac{\partial G_A}{\partial P}\right)_T\, dP = \left(\frac{\partial G_B}{\partial T}\right)_P\, dT + \left(\frac{\partial G_B}{\partial P}\right)_T\, dP.
\tag{8.44}
$$

The Gibbs derivatives are $(\partial G/\partial T)_P = -S$ and $(\partial G/\partial P)_T = V$. Substituting,

$$
-S_A\, dT + V_A\, dP = -S_B\, dT + V_B\, dP,
$$

and rearranging,

$$
\boxed{\frac{dP}{dT} = \frac{S_B - S_A}{V_B - V_A} = \frac{\Delta S}{\Delta V}.}
\tag{8.45}
$$

For a first-order transition the latent heat is $L = T \Delta S$, so (8.45) takes its more familiar form

$$
\frac{dP}{dT} = \frac{\Delta H}{T\, \Delta V},
\tag{8.46}
$$

with $\Delta H$ the enthalpy difference between phases at the boundary. Equation (8.46) is exact within equilibrium thermodynamics — it makes no model assumption beyond the first-order character of the transition.

For melting, $\Delta H_\mathrm{fus} > 0$ (heat must be supplied to melt) and typically $\Delta V > 0$ (liquid is less dense than solid), so $dP/dT > 0$: the melting line tilts up. The famous exception is water, with $\Delta V < 0$ (ice floats), so $dP/dT < 0$ and applying pressure melts ice. For solid-solid transitions, both signs are common and the boundary can be steep or shallow.

Clausius-Clapeyron is the basic cross-check on any simulated phase diagram. Compute $T_m(P)$ at two pressures, and also compute $\Delta H$ and $\Delta V$ directly from coexistence simulations. The numerical $dP/dT$ from the boundary trajectory and the right-hand side of (8.46) should agree to within statistical noise. If they disagree by a factor of two, the simulation is unequilibrated, the latent-heat estimate is wrong, or the coexistence has been mis-identified.

### Why direct heating overestimates $T_m$

A perfect single crystal in a periodic box, heated above its thermodynamic $T_m$, does not melt. To start melting it must *nucleate* a liquid embryo of size larger than the critical-nucleus radius $r^*$, against the bulk-liquid-solid interface free energy $\gamma_{sl}$. Classical nucleation theory gives the work of forming an embryo of radius $r$ as

$$
W(r) = -\frac{4\pi}{3} r^3\, \rho_l\, |\Delta\mu_{sl}(T)| + 4\pi r^2 \gamma_{sl},
\tag{8.47}
$$

with $\Delta\mu_{sl} = \mu_l - \mu_s$ the chemical-potential difference, vanishing at $T_m$ and growing approximately linearly with $T - T_m$. The barrier height is

$$
W^* = \frac{16\pi}{3} \frac{\gamma_{sl}^3}{(\rho_l\, |\Delta\mu_{sl}|)^2},
\tag{8.48}
$$

which diverges as $T \to T_m^+$. To melt on the simulation timescale, $T$ must be raised enough that $W^*$ drops to a few $k_B T$ — typically $10$–$30$% above $T_m$. The single-phase heating method therefore reports a **superheating limit** $T_+$ that overestimates $T_m$ systematically. The analogous problem on cooling — undercooling a liquid below $T_m$ before crystallisation nucleates — gives a $T_-$ that underestimates $T_m$.

For Cu with the Mishin EAM, single-phase heating typically reports $T_+ \approx 1500$–$1600$ K, against the true thermodynamic $T_m \approx 1325$ K. The error is $13$–$20$%. Reporting $T_+$ as "the melting point from MD" is one of the most common — and most embarrassing — errors in the literature.

### The hysteresis method

Luo, Strachan and Swift (2004) proposed a clever workaround. Heat a solid until it melts at $T_+$, then cool the resulting liquid until it crystallises at $T_-$. Both transitions are nucleation-limited and the gap is large. Their empirical observation, since validated on dozens of systems, is that

$$
T_m \approx \frac{T_+ + T_-}{2}
\tag{8.49}
$$

to within a few per cent, provided the heating and cooling rates are comparable. The qualitative argument is that the superheating limit and the undercooling limit lie roughly symmetrically about $T_m$ in $\Delta T$, because both barriers scale as $W^* \sim 1/(\Delta T)^2$.

This is cheaper than the two-phase method — one continuous heating-cooling trajectory rather than several distinct simulations at neighbouring temperatures — and is the natural first estimate when bracketing $T_m$. For publication-quality numbers the two-phase method remains the gold standard, but for a quick sanity check on an MLIP's predicted $T_m$ the hysteresis method gives results in a single overnight run.

A representative LAMMPS protocol:

```text
# Heating ramp
fix    1 all npt temp 300 2200 0.1 iso 0.0 0.0 1.0
run    400000
# Watch potential energy: a discontinuous jump marks T_+
# Cooling ramp
fix    1 all npt temp 2200 300 0.1 iso 0.0 0.0 1.0
run    400000
# Discontinuous jump on cooling marks T_-
```

Plot $\langle U\rangle$ versus $T$. The hysteresis loop is a classic double S-curve. Pick $T_+$ and $T_-$ at the inflection points; $T_m$ from (8.49).

### Connection to Project 2 (Cu melting via MLIP)

The Cu-melting project of [Chapter 13](../projects/02-melting-point-mlip/README.md) puts all of this together. The recipe is:

1. **Hysteresis estimate.** Drive a $4 \times 4 \times 4$ fcc Cu supercell with MACE-MP-0 through a heating-then-cooling cycle to bracket $T_m$ to within $\sim 100$ K.
2. **Two-phase refinement.** Build the half-solid-half-liquid supercell of [§8.4](04-phase-diagrams.md), equilibrate at the hysteresis midpoint, and run at three temperatures spanning $\pm 25$ K. Identify the temperature at which the interface is stationary.
3. **Clausius-Clapeyron cross-check.** Repeat steps 1-2 at $P = 1$ GPa and $P = 5$ GPa; verify $dT_m/dP$ is consistent with $\Delta V / \Delta S$ measured at the coexistence point.
4. **Comparison to experiment and to classical EAM.** Mishin EAM gives $T_m \approx 1325$ K; experiment is $1358$ K; a fine-tuned MACE-MP-0 typically lands in between, with the gap shrinking as the fine-tuning data approach the melting region. The size of the residual error tells you how much DFT-MD data on near-melting configurations is missing from the foundation model's training distribution.

The whole pipeline runs in a long weekend on a single A100. It also exposes the limits of every method along the way: the hysteresis estimate has the systematic bias just discussed, the two-phase result is sensitive to interface orientation and finite-size effects, and the Clausius-Clapeyron check fails if pressures are too low to resolve $\Delta V$ above the volume-fluctuation noise. Resolving each of these is the point.

## A worked example: argon's triple-point temperature

Lennard-Jones argon. Triple point known analytically: $T^* = 0.687, P^* = 0.0017, \rho^*_\mathrm{liq} = 0.847, \rho^*_\mathrm{sol} = 0.96$ (Hansen-Verlet 1969, refined since).

To estimate $T_m$ at $P = 0$ via two-phase coexistence, recall that $P^* \approx 0$ corresponds to the solid-liquid line just above the triple point. Build a 4x4x16 fcc supercell (1024 atoms), melt the top half, equilibrate at $T^* = 0.7$ in NPT with $P^* = 0$. Run for $\sim 10^6$ steps and observe:

- At $T^* = 0.75$: the liquid grows; the box ends up fully liquid.
- At $T^* = 0.70$: the interface is roughly stationary on the available simulation timescale.
- At $T^* = 0.65$: the solid grows.

The bracket $T^*_m \approx 0.70 \pm 0.03$ matches the literature triple-point value within statistical uncertainty.

In LAMMPS `lj` units, $T^* = 0.7$ is just `temp 0.7 0.7 0.5`; the rest of the script is exactly the Cu sketch above. A laptop completes the run overnight.

!!! note "Force-field-dependent phase diagrams"
    The phase diagram you compute is the phase diagram **of the force field**, not of the real material. Mishin Cu EAM gives $T_m \approx 1325$ K against experimental $T_m = 1358$ K — a small but non-zero error. Lennard-Jones argon's $T_m$ agrees better with experiment because argon is genuinely well-described by LJ. When MLIPs (Chapter 9) replace classical force fields, the phase-diagram agreement with experiment improves, sometimes dramatically (especially for non-trivial chemistries like ice polymorphs or oxide phase boundaries).

## What we have

The pipeline from atomistic simulation to a phase diagram is now complete: integrators (§7.1), thermostats (§7.3), force models (§7.4), ensembles (§8.1), free energies (§8.2), and the methods of this section. Each is an essential link.

The looming limitation is force-field accuracy. The classical force fields of §7.4 give phase boundaries good to a few percent — sometimes good enough, sometimes not. Closing that gap is what [Chapter 9](../ch09-mlip/index.md) is about: learning interatomic potentials directly from DFT, achieving near-DFT accuracy at near-classical cost. The downstream simulation machinery — everything in Chapters 7 and 8 — applies unchanged.

The exercises ([§8.5](exercises.md)) put these methods to work.
