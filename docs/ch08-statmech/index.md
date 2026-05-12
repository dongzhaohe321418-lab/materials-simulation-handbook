# Chapter 8 — Statistical Mechanics from Simulations

[Open in Jupyter (browser)](/materials-simulation-handbook/lite/lab/index.html?path=ch08-statmech.ipynb){ .md-button .md-button--primary }

In [Chapter 7](../ch07-md/index.md) we ran molecular dynamics simulations and pulled observable quantities — energies, pressures, mean squared displacements, radial distribution functions — out of trajectories. What we did not do was step back and ask what those quantities are estimators **of**. A pressure averaged over 100 ps of LAMMPS output is a sample mean. What population mean is it estimating? Under what conditions does the sample mean converge to the population mean? What if we want a free energy — a quantity that, unlike energy or pressure, has no instantaneous-snapshot value?

This chapter is about the statistical mechanics that gives those questions sharp answers. Most of it you will already have met in a physics undergraduate degree. The point here is not to teach statistical mechanics from scratch but to connect each piece of it explicitly to what MD does. We will go through ensembles (which one does each LAMMPS `fix` actually sample?), partition functions (how do we extract them from samples?), free energies (the methods that built modern computational drug discovery and computational materials design), transport coefficients (the Green-Kubo relations), and phase diagrams from simulation.

## What you will learn

By the end of this chapter you will be able to:

1. Match each MD ensemble to the appropriate thermodynamic potential and identify which fluctuations carry which information (heat capacities from energy variance, compressibilities from volume variance).
2. Compute the free energy difference between two states by thermodynamic integration, free energy perturbation, or umbrella sampling, and explain when each is the right tool.
3. Extract a diffusion coefficient or viscosity from a trajectory using the Green-Kubo integral, and recognise when that integral has not converged.
4. Locate phase boundaries by comparing free energies across phases and by two-phase coexistence simulations.

## The path

1. **Ensembles and partition functions** ([§8.1](01-ensembles.md)). The microcanonical, canonical, isobaric-isothermal, and grand canonical ensembles. Partition functions and free energies. Which thermostat samples which ensemble. Fluctuation formulas.
2. **Free energy methods** ([§8.2](02-free-energy.md)). Thermodynamic integration, free energy perturbation, umbrella sampling, metadynamics. A worked example computing the excess chemical potential of a Lennard-Jones liquid.
3. **Transport coefficients** ([§8.3](03-transport.md)). Green-Kubo relations for diffusion, viscosity, and thermal conductivity. The fluctuation-dissipation theorem. Comparison with the Einstein relation and with non-equilibrium MD.
4. **Phase diagrams** ([§8.4](04-phase-diagrams.md)). Free energies across phases, two-phase coexistence simulations, melting from MD, and when to defer to assessed CALPHAD databases.

## Why this matters

A great deal of computational materials work — phase stability, defect thermodynamics, adsorption energies, solubility, reaction barriers in solution — reduces to a free energy difference. You can compute energies trivially. Free energies you cannot: $F = -k_B T \ln Z$, and $Z$ is an integral over an exponentially large configuration space. The methods in this chapter are how that integral becomes tractable.

Similarly, transport properties — diffusion coefficients, viscosities, thermal conductivities — are not directly accessible from energy or force evaluations. They live in the **correlations** of equilibrium fluctuations. The Green-Kubo formalism turns those correlations into transport coefficients, but only if you understand which correlation function to compute and how long to integrate.

## What this chapter is not

This is not a textbook of statistical mechanics. We will not derive the canonical distribution from the postulate of equal a priori probabilities, nor revisit the saddle-point approximation that leads from microcanonical to canonical. Those are in any of Pathria, Reichl, or Tuckerman, and we will reference where appropriate. Our focus is on the link from formalism to simulation: each result is something you will compute from an MD trajectory.

We also do not cover Monte Carlo sampling in depth. MC is the natural partner to MD for sampling problems where dynamics is irrelevant — lattice models, polymer conformations, hard-sphere fluids. Where MC enters this chapter (umbrella sampling, parallel tempering, replica exchange) we will introduce it briefly; a fuller treatment lives in any of Frenkel & Smit or Landau & Binder.

## Prerequisites

You should be comfortable with:

- The Boltzmann distribution $P \propto e^{-\beta E}$ at the level of an introductory statistical mechanics course.
- Multivariable calculus through change-of-variables and partial derivatives of thermodynamic potentials.
- The molecular dynamics machinery of [Chapter 7](../ch07-md/index.md): integrators, thermostats, ensembles in their MD-operational sense.
- Numerical integration and basic stochastic-process intuition (autocorrelation, variance, central limit theorem).

The methods of this chapter generalise directly to any force model — classical force field, DFT, or MLIP. Where examples use a specific potential, it is LJ in `metal` units for concreteness; the algorithms are independent.

## A note on conventions

We will write $\beta = 1/(k_B T)$ throughout, free energies as $A$ (Helmholtz) and $G$ (Gibbs) following classical thermodynamics convention, and ensemble averages as $\langle \cdot \rangle$ without specifying the ensemble when context makes it clear. Energies in derivations are in units of $k_B T$ where it simplifies; in code examples we keep eV / kcal/mol / reduced LJ explicit.

Let us begin where every statistical-mechanical discussion must: with the choice of ensemble.
