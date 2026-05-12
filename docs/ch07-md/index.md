# Chapter 7 — Molecular Dynamics

In [Chapter 6](../ch06-running-dft/index.md) we computed forces on atoms from DFT and stopped there. A force on a nucleus is the gradient of the Born-Oppenheimer energy surface; at a stationary point those forces are zero, and that is what a relaxation calculation hunts for. Everything we did in Chapters 5 and 6 was static. Time did not appear.

Real materials are not static. Atoms in a copper crystal at 300 K oscillate around their lattice sites at frequencies of a few THz, each carrying about $\tfrac{3}{2}k_B T \approx 39$ meV of kinetic energy. Hydrogen in palladium hops between octahedral sites every nanosecond or so. Water molecules in a biological pocket rearrange their hydrogen-bond network on timescales of picoseconds. None of this can be read off a single zero-temperature potential energy surface; it requires integrating Newton's equations.

This chapter is about how to do that integration well, and what physics you extract from the resulting trajectories. By the end you will be able to take a force model — whether from a classical potential, from DFT, or from one of the machine-learning potentials of [Chapter 9](../ch09-mlip/index.md) — and turn it into a numerical trajectory $\{\mathbf{r}_i(t)\}_{i=1\ldots N}$ that samples a well-defined statistical ensemble. You will also be able to test, after the fact, whether your sampling was adequate.

## What you will build

By the end of this chapter you will have:

1. Written a velocity-Verlet integrator from scratch in NumPy and verified that it conserves energy to one part in $10^6$ over a million steps.
2. Run a Lennard-Jones argon simulation in LAMMPS and computed its self-diffusion coefficient.
3. Driven a copper melting calculation using EAM at constant pressure, watched the volume jump at $T_m$, and located the (approximate) melting temperature.
4. Computed a radial distribution function and a mean squared displacement from a trajectory file using your own code, then cross-checked with MDAnalysis.

## The path

We proceed in six steps:

1. **Integration** ([§7.1](01-integration.md)). Newton's second law for atoms, the Verlet family of integrators, why they are symplectic and why that matters, and how to choose the time step.
2. **Periodic boundary conditions** ([§7.2](02-pbc.md)). Why we use them, the minimum image convention, the subtle distinction between wrapped and unwrapped coordinates, and how Ewald summation handles long-range Coulomb interactions.
3. **Thermostats and barostats** ([§7.3](03-thermostats.md)). How to drive an MD simulation to sample NVT, NPT, and other ensembles instead of just NVE. The Nosé-Hoover formalism, Langevin dynamics, and the Parrinello-Rahman barostat.
4. **Force fields** ([§7.4](04-force-fields.md)). Lennard-Jones, EAM for metals, Tersoff and bond-order potentials for covalent solids, and ReaxFF. The accuracy-cost gap that motivates everything in [Chapter 9](../ch09-mlip/index.md).
5. **LAMMPS workflows** ([§7.5](05-lammps.md)). The dominant open-source MD code. Installation, the anatomy of an input script, two end-to-end workflows, and when to drive LAMMPS from Python via ASE.
6. **Analysing trajectories** ([§7.6](06-analysis.md)). MSD and diffusion coefficients, radial distribution functions, velocity autocorrelations and vibrational densities of states, structure factors.

## What this chapter is not

This is not a survey of MD codes. GROMACS, NAMD, OpenMM, HOOMD-blue, DL_POLY, AMBER and ESPResSo each have devoted user communities; once you understand LAMMPS, switching is largely a matter of syntax. We also do not cover ab-initio MD in depth here — `cp.x` from Quantum ESPRESSO, CP2K's Born-Oppenheimer MD, and VASP's MD all use the integrators of [§7.1](01-integration.md) on top of DFT forces, and we have already met those forces in [Chapter 6](../ch06-running-dft/index.md). The integration is the same; the cost is $10^4$ to $10^6$ times higher per step.

Coarse-grained MD, dissipative particle dynamics, Brownian dynamics for polymers, and Monte Carlo sampling are out of scope. Monte Carlo for materials is covered briefly in [Chapter 8](../ch08-statmech/index.md) when we need it for free energies.

## Prerequisites

You should be comfortable with:

- Newton's laws and elementary classical mechanics.
- NumPy array manipulation at the level of [Chapter 1](../ch01-python/index.md).
- The notion of a potential energy surface from [Chapter 6](../ch06-running-dft/index.md).
- Some statistical mechanics — at least the Boltzmann distribution. We will recap the rest in [Chapter 8](../ch08-statmech/index.md).

You will need LAMMPS (`conda install -c conda-forge lammps`), ASE 3.23 or later, and MDAnalysis 2.6 or later. Total disk for the worked examples is about 1 GB.

## A note on units

MD codes use a confusion of unit systems. LAMMPS has `real`, `metal`, `lj`, and several others; ASE works in eV, Å, and amu; many force-field papers report energies in kcal/mol. We use eV/Å/amu (the ASE convention) in derivations, switch to LAMMPS' `metal` units (eV, Å, ps, bar) for production scripts, and quote final results in SI when comparing to experiment. The relevant unit conversions are tabulated in [Appendix A](../appendix/A-math-reference.md).

Let us begin with the equations we are trying to solve.
