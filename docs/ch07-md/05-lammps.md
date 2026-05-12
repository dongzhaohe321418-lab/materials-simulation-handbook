# 7.5 LAMMPS Workflows

LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator), developed at Sandia since the mid-1990s, is the de facto open-source molecular dynamics engine for materials science. It runs on anything from a laptop to a Top-500 machine, ships with hundreds of pair styles, fixes, and computes, and is the back end for almost every modern MLIP package (Allegro, NequIP, MACE, GAP, SNAP, ACE — all integrate via LAMMPS).

This section will get you to two working simulations: a Lennard-Jones argon equilibrium in NVE/NVT, and a copper melting calculation in NPT with EAM. Along the way you will learn what every standard LAMMPS input line does and when to reach for ASE instead.

## Installation

The simplest installation, sufficient for everything in this chapter:

```bash
conda install -c conda-forge lammps
```

This gives a serial executable `lmp` plus a parallel one launched with `mpirun lmp_mpi`. For MLIPs you generally need to build LAMMPS from source with the relevant package (`make yes-ml-iap` for SNAP/ACE; specific build instructions ship with each MLIP); this is covered in Chapter 9.

Test the install:

```bash
lmp -in - <<'EOF'
units lj
atom_style atomic
region box block 0 10 0 10 0 10
create_box 1 box
create_atoms 1 random 100 12345 box
mass 1 1.0
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0
velocity all create 1.0 87287
fix 1 all nve
run 100
EOF
```

If you see thermodynamic output ending in a `run 100` summary with no errors, LAMMPS is happy.

## Anatomy of a LAMMPS input script

A LAMMPS input script is read line by line and executed immediately. The main sections, in conventional order:

| Section | Commands | Purpose |
|---|---|---|
| Units & atom style | `units`, `atom_style`, `dimension`, `boundary` | Choose unit system and how atoms are represented |
| Domain | `region`, `create_box`, `lattice`, `read_data` | Define the simulation cell and lattice |
| Atoms | `create_atoms`, `mass`, `set` | Place atoms and set their properties |
| Interactions | `pair_style`, `pair_coeff`, `bond_style`, `kspace_style` | Define the potential |
| Initial conditions | `velocity`, `displace_atoms` | Set initial velocities and positions |
| Settings | `neighbor`, `timestep`, `thermo`, `thermo_style` | Numerical and output settings |
| Fixes | `fix nvt`, `fix npt`, `fix langevin` | Apply ensembles and constraints |
| Computes & dumps | `compute`, `dump` | Define and write observables |
| Run | `run`, `minimize` | Execute |

A few critical commands deserve detailed treatment.

**`units`.** The choice of unit system affects every numerical parameter in your script. Common choices:
- `lj`: dimensionless reduced units. Energies in $\epsilon$, lengths in $\sigma$, times in $\tau = \sqrt{m\sigma^2/\epsilon}$. Use for Lennard-Jones model studies.
- `metal`: eV, Å, ps, K, bar, g/mol. The standard for materials science.
- `real`: kcal/mol, Å, fs, K, atm. The standard for biomolecular MD inherited from CHARMM/AMBER.

**`atom_style`.** What information each atom carries.
- `atomic`: just position, velocity, type. Fine for metals, noble gases, simple solids.
- `charge`: also a charge per atom. Needed for ionic systems.
- `molecular`: also bonds, angles, dihedrals. For organic molecules with fixed topology.
- `full`: charges plus topology. The biomolecular workhorse.

**`boundary p p p`.** Periodic in all three directions (the default for bulk simulations). Use `p p f` for a slab geometry with a free surface in $z$. The `f` direction needs a `fix wall` or sufficient vacuum to prevent atoms escaping.

**`pair_style` and `pair_coeff`.** Choose the interatomic potential and its parameters. For LJ: `pair_style lj/cut 2.5` then `pair_coeff 1 1 epsilon sigma`. For EAM: `pair_style eam/alloy` then `pair_coeff * * Cu_mishin1.eam.alloy Cu`. For machine-learned potentials, a pair style like `mliap` or `allegro` is used.

**`fix`.** Applies a "fix" — a per-timestep modification of the dynamics. The standard production fixes are:
- `fix 1 all nve`: pure Newton, NVE.
- `fix 1 all nvt temp 300 300 0.1`: Nosé-Hoover thermostat at 300 K with $\tau_T = 0.1$ ps.
- `fix 1 all npt temp 300 300 0.1 iso 1.0 1.0 1.0`: MTK isotropic NPT at 300 K, 1 bar.
- `fix 1 all langevin 300 300 0.1 12345`: Langevin thermostat, seed 12345.

The `fix` ID (`1`) and group (`all`) let you stack multiple fixes — for example, Langevin on a thermostat region plus NVE on the bulk for non-equilibrium MD.

**`thermo` and `dump`.** `thermo 100` prints thermodynamic output every 100 steps; `dump 1 all atom 1000 traj.lammpstrj` writes a snapshot of positions every 1000 steps. The dump file is what you analyse later (§7.6).

**`run`.** Execute integration for $N$ steps. Multiple `run` commands can be chained, possibly with different fixes active.

## Workflow 1: Lennard-Jones argon

We will simulate 4000 argon atoms in a face-centred cubic crystal at low temperature, melt it, equilibrate the liquid in NVT, then compute mean potential energy and pressure. The script (`argon.lmp`):

```text
# === Lennard-Jones argon in metal units ===
units           metal
atom_style      atomic
boundary        p p p

# Argon parameters in metal units: epsilon = 0.0103 eV, sigma = 3.405 A
variable        eps equal 0.0103
variable        sig equal 3.405
variable        a   equal ${sig}*2^(1/6)*sqrt(2)  # fcc lattice constant

lattice         fcc ${a}
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 39.948

pair_style      lj/cut 8.5125          # 2.5 * sigma
pair_coeff      1 1 ${eps} ${sig}
pair_modify     shift yes

neighbor        2.0 bin
neigh_modify    every 10 delay 0 check yes

velocity        all create 100.0 4928459 mom yes rot yes dist gaussian

timestep        0.002                  # 2 fs

# === Stage 1: equilibration in NVT, slow heating ===
thermo          200
thermo_style    custom step temp pe ke etotal press vol density
fix             1 all nvt temp 100.0 100.0 0.5
run             20000
unfix           1

# === Stage 2: production in NVT ===
fix             1 all nvt temp 100.0 100.0 0.5
dump            1 all custom 500 argon.lammpstrj id type x y z vx vy vz
compute         msd_argon all msd
fix             msd_out all ave/time 100 1 100 c_msd_argon[4] file msd.dat
run             100000
```

What this does, step by step:

1. **Lines 1–3.** Choose `metal` units (eV/Å/ps/bar), atomistic atoms with no charges, periodic in all directions.
2. **Lines 6–7.** Define $\epsilon$ and $\sigma$ for argon. Both in `metal` units.
3. **Line 8.** Compute the fcc lattice constant from $\sigma$: the nearest-neighbour distance in fcc is the LJ minimum $2^{1/6}\sigma$, and the conventional cubic cell parameter is $\sqrt{2}$ times that.
4. **Lines 10–14.** Build a $10 \times 10 \times 10$ fcc supercell of argon: 4000 atoms.
5. **Lines 16–18.** Lennard-Jones with cutoff $2.5\sigma$, shifted to zero at the cutoff.
6. **Line 20.** Build a neighbour list with a 2 Å skin, rebuilt every 10 steps when needed.
7. **Line 23.** Initialise velocities from a Gaussian at 100 K, with zero centre-of-mass momentum and angular momentum.
8. **Line 25.** Time step of 2 fs — comfortable for argon (lightest atomic mass 40, no hydrogens).
9. **Lines 28–32.** Equilibration in NVT at 100 K for 40 ps, printing thermodynamic averages every 200 steps.
10. **Lines 34–40.** Production NVT for 200 ps, writing trajectory every 1 ps and time-averaged MSD to `msd.dat`.

Running this on a laptop takes a couple of minutes. The thermodynamic output should show a mean potential energy near $-0.066$ eV/atom (close to LJ liquid argon at 100 K and equilibrium density) and pressure near zero (since we initialised at the equilibrium lattice constant).

To compute the pressure as a thermodynamic observable: the LAMMPS variable `press` is the instantaneous virial pressure

$$
P = \frac{N k_B T}{V} + \frac{1}{3V}\sum_{i<j} \mathbf{r}_{ij} \cdot \mathbf{F}_{ij},
\tag{7.50}
$$

and `thermo_style custom press` writes it to the log; time-average it over the production run.

To compute total energy: the `etotal` thermo variable is $H = KE + PE$. In NVE this is constant by construction; in NVT it fluctuates and you average it.

## Workflow 2: copper melting in NPT with EAM

A more substantial run: melt fcc copper using the Mishin EAM potential, locating $T_m$ approximately by the volume jump in NPT. Sketch:

```text
units           metal
atom_style      atomic
boundary        p p p

lattice         fcc 3.615
region          box block 0 8 0 8 0 8
create_box      1 box
create_atoms    1 box
mass            1 63.546

pair_style      eam/alloy
pair_coeff      * * Cu_mishin1.eam.alloy Cu

velocity        all create 300.0 87287 mom yes rot yes dist gaussian
timestep        0.002

# Ramp temperature from 300 to 2000 K in NPT at 1 bar
fix             1 all npt temp 300.0 2000.0 0.1 iso 1.0 1.0 1.0
thermo          500
thermo_style    custom step temp pe ke etotal press vol density

dump            1 all custom 1000 cu_melt.lammpstrj id type x y z
run             500000     # 1 ns ramp
```

You will need the EAM file `Cu_mishin1.eam.alloy` from the LAMMPS `potentials/` directory (it ships with the conda install). Plot temperature vs. volume from the log: at the melting transition the volume rises sharply (about 5% jump in Cu) over a narrow temperature window. The midpoint of the jump is your estimate of $T_m$ for this potential.

Mishin EAM gives $T_m \approx 1340$ K for Cu, against an experimental value of 1358 K — about 1% accurate, which is excellent for a 1990s force field. Beware that NPT ramps are biased estimators: the real melting point is bracketed by the heating and cooling hysteresis, and the more reliable method is the two-phase coexistence simulation introduced in [§8.4](../ch08-statmech/04-phase-diagrams.md).

!!! tip "Avoid hysteresis with the two-phase method"
    Heating a perfect crystal in NPT will overshoot the true $T_m$ because nucleation of the liquid is slow. Cooling a liquid will undershoot $T_m$ because nucleation of the crystal is slower still. The two-phase coexistence method (start with half solid, half liquid) avoids nucleation barriers entirely.

## ASE driver alternative

For many production workflows you do not need to write a raw LAMMPS input. The ASE `LAMMPSlib` and `LAMMPSrun` calculators let you drive LAMMPS from Python:

```python
from ase.build import bulk
from ase.calculators.lammpsrun import LAMMPS
from ase.md.langevin import Langevin
from ase import units

# Build a copper supercell
atoms = bulk("Cu", "fcc", a=3.615, cubic=True).repeat((8, 8, 8))

# LAMMPS as a calculator: forces and energies come from LAMMPS
parameters = {
    "pair_style": "eam/alloy",
    "pair_coeff": ["* * Cu_mishin1.eam.alloy Cu"],
}
atoms.calc = LAMMPS(parameters=parameters, files=["Cu_mishin1.eam.alloy"])

# But run dynamics with ASE's integrator
dyn = Langevin(atoms, timestep=2 * units.fs, temperature_K=1200,
               friction=0.01)
dyn.run(50000)
```

The trade-off is clear:

- **Pure LAMMPS input**: faster, more efficient for large systems, gives access to all of LAMMPS' built-in fixes and computes, parallelises seamlessly.
- **ASE driving LAMMPS**: more flexible for unusual workflows, easy to mix in Python analysis on the fly, lets you use ASE's library of MD algorithms (NEB, dimer method, structure search) with LAMMPS forces. Some overhead per step.

A good rule: for production MD of fixed systems (thousands of atoms, nanosecond runs), use raw LAMMPS. For exploratory workflows, structure searches, NEB transition-state calculations, or pipelines that mix multiple calculators, use ASE. For MLIP work, the standard pattern is to write a LAMMPS input that loads the MLIP and dump trajectories that you analyse in Python — best of both worlds.

## Common LAMMPS pitfalls

A non-exhaustive list of footguns:

- **Units forgotten.** Specifying `pair_coeff 1 1 0.0103 3.405` in a script with `units real` will produce nonsense; in `units real`, energies are in kcal/mol, not eV. Match `units` to the parameters of your force field.
- **Velocity initialisation.** `velocity all create T seed mom no rot no` keeps the centre-of-mass and angular momentum, which then act as additional kinetic energy not associated with temperature, leading to an effective temperature lower than $T$. Always use `mom yes rot yes`.
- **NPT during phase transition.** A NPT simulation crossing a first-order phase boundary will not show a clean jump in volume but a slow drift while nucleation happens. Do not interpret this as poor equilibration.
- **Cutoff > L/2.** The minimum image convention breaks. LAMMPS will warn; do not ignore the warning.
- **`reset_timestep` resets only the integer step counter, not the dynamics.** Velocities and positions are unchanged.
- **`fix nvt` writes its own thermostat conserved quantity to `c_thermo_press` and friends — these are not the bare physical pressure but include the thermostat contribution.** When in doubt, define your own `compute pressure` outside the fix.

## What we have

LAMMPS scripts that run to completion and produce trajectory files. Those trajectories are the raw input to the analysis methods of [§7.6](06-analysis.md): radial distribution functions, mean squared displacements, velocity autocorrelations, structure factors. From those quantities we extract the actual physics: diffusion coefficients, phonon dispersions, structural order parameters.
