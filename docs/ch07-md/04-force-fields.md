# 7.4 Force Fields

<figure markdown>
![Lennard-Jones potential](../assets/figures/ch07/fig_lj_potential.png){ width="600" }
<figcaption>Figure 7.4.1. The Lennard-Jones 12-6 potential: the canonical pair potential, with a steep repulsive wall, a single minimum at \(r_{\min} = 2^{1/6}\sigma\), and a smooth attractive tail. (Synthetic curve in reduced units.)</figcaption>
</figure>

The integrators of §7.1 and the thermostats of §7.3 are agnostic to where forces come from. In [Chapter 6](../ch06-running-dft/index.md) we computed forces from DFT — accurate but expensive, scaling as $O(N^3)$ for hybrid functionals and limiting routine MD to a few hundred atoms for picoseconds. To reach the thousands-of-atoms, nanoseconds-and-beyond regime that materials science actually needs, we replace the DFT force calculation by an analytic functional form fitted to reference data. This is a **classical force field**.

The art of force-field design is a half-century-old enterprise. We survey four families that span the chemistry of materials: simple pair potentials, density-dependent metal potentials, bond-order potentials for covalent solids, and reactive force fields. Each has a regime where it is excellent and a regime where it fails badly. Together they motivate the machine-learning potentials of [Chapter 9](../ch09-mlip/index.md).

## Lennard-Jones

The simplest non-trivial pair potential is

$$
U_\mathrm{LJ}(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right].
\tag{7.44}
$$

The $r^{-12}$ term is a phenomenological short-range repulsion (Pauli exclusion); the $r^{-6}$ term comes from the second-order multipole expansion of the dispersion interaction (Section 2.6.2 of any electronic-structure textbook gives the derivation from London's 1930 result). Two parameters, $\epsilon$ and $\sigma$, set the energy and length scales.

**Equilibrium spacing.** At the minimum of $U_\mathrm{LJ}$,

$$
\frac{dU_\mathrm{LJ}}{dr} = 0 \quad \Longrightarrow \quad
-12\,\frac{\sigma^{12}}{r^{13}} + 6\,\frac{\sigma^{6}}{r^{7}} = 0,
$$

giving $r_\mathrm{eq} = 2^{1/6}\sigma \approx 1.122\sigma$ and well depth $U_\mathrm{LJ}(r_\mathrm{eq}) = -\epsilon$.

**Reduced units.** Express lengths in $\sigma$, energies in $\epsilon$, times in $\tau = \sqrt{m\sigma^2/\epsilon}$, temperatures in $\epsilon/k_B$. The reduced LJ phase diagram is universal: triple point at $T^* \approx 0.69$, $\rho^* \approx 0.85$; critical point at $T^* \approx 1.31$, $\rho^* \approx 0.31$. These numbers are independent of the material — what makes "argon" different from "krypton" is the values of $\epsilon$ and $\sigma$.

**When LJ is enough.** Liquid noble gases (Ar, Kr, Xe) at moderate density are described almost quantitatively by LJ. Argon in particular has $\epsilon/k_B \approx 120$ K and $\sigma \approx 3.40$ Å, giving a predicted melting temperature near 84 K (experimental: 83.8 K) and a critical temperature of 158 K (experimental: 150.7 K) — agreement at the percent level. LJ is also the workhorse model for **conceptual** simulations of liquids: glass formation, nucleation, hydrodynamics, where universality matters more than chemical specificity.

LJ is the **wrong** model for: anything covalent, anything metallic, anything with hydrogen bonds, anything with significant polarisability. The bonding physics in those systems is not pairwise.

## Embedded Atom Method for metals

A copper atom in bulk copper has 12 nearest neighbours. Adding a 13th neighbour costs more energy than adding the first 12 — the cohesive energy per neighbour decreases as the coordination increases. This **bond order** effect cannot be captured by any pairwise potential, where the per-neighbour energy is constant.

The Embedded Atom Method (Daw and Baskes, 1984) handles this by writing the energy as

$$
E = \sum_i F_i(\bar\rho_i) + \tfrac{1}{2} \sum_{i \ne j} \phi_{ij}(r_{ij}),
\qquad
\bar\rho_i = \sum_{j \ne i} \rho_j(r_{ij}).
\tag{7.45}
$$

Atom $i$ feels a "host electron density" $\bar\rho_i$ contributed by all its neighbours through pairwise density functions $\rho_j(r)$. The energy to embed atom $i$ in that density is a non-linear function $F_i(\bar\rho)$ — the **embedding function**. On top of this sits a pairwise potential $\phi_{ij}(r_{ij})$ that captures core-core repulsion.

The non-linearity of $F$ is what allows EAM to encode coordination dependence. For a typical EAM functional, $F$ is concave: doubling the density does not double the embedding energy. This is exactly the right behaviour for metals.

The forces are slightly more involved than in a pair potential:

$$
\mathbf{F}_i = -\nabla_i E = -\sum_{j \ne i} \left[\phi'_{ij}(r_{ij}) + F'_i(\bar\rho_i) \rho'_j(r_{ij}) + F'_j(\bar\rho_j) \rho'_i(r_{ij})\right] \hat{\mathbf{r}}_{ij}.
\tag{7.46}
$$

Each force depends on the local densities at **both** atoms, which means EAM is not pairwise additive but its computational cost is only modestly higher than LJ — still $O(N)$ with neighbour lists.

EAM and its closely related cousins (Finnis-Sinclair, second-moment tight-binding, MEAM with angular terms) are the standard for simulating elemental and alloy metals. Public databases like the NIST Interatomic Potentials Repository host parameterisations for most metals; the original Foiles-Baskes-Daw set for Cu, Ag, Au, Ni, Pd, Pt is still in heavy use four decades on.

**Where EAM fails.** Surface reconstructions where directional bonding matters (Au(111) herringbone, Si(111) 7×7) are outside EAM's reach. Alloys with strong charge transfer (NiAl, FeAl) are marginal. Hydrogen in metals, where H sees a much sharper potential than the host atoms, often requires specialised parameterisations.

## Tersoff and bond-order potentials

For covalent solids — silicon, carbon, germanium, SiC, GaN — the bonding is highly directional. A silicon atom forms four tetrahedral bonds; bond angles around 109.5° are stabilised, others are penalised. A pair potential plus an EAM-style density cannot do this.

Tersoff's 1988 potential builds in angular information through a bond-order term. Write the energy as

$$
E = \tfrac{1}{2} \sum_{i \ne j} f_C(r_{ij}) \left[V_R(r_{ij}) - b_{ij} V_A(r_{ij})\right],
\tag{7.47}
$$

where $V_R$ is a repulsive Morse-like term, $V_A$ an attractive term, $f_C$ a smooth cutoff function, and $b_{ij}$ the **bond order** between atoms $i$ and $j$. The bond order depends on the local environment of $i$:

$$
b_{ij} = \left(1 + (\beta\, \zeta_{ij})^n\right)^{-1/(2n)},
\qquad
\zeta_{ij} = \sum_{k \ne i, j} f_C(r_{ik})\, g(\theta_{ijk}),
\tag{7.48}
$$

with $g(\theta)$ an angular function that favours specific bond angles. A bond between $i$ and $j$ is weakened (smaller $b_{ij}$) if there are many other neighbours $k$ of atom $i$ at "wrong" angles. The effect is to reproduce the saturated valence of covalent bonding: silicon's preference for four bonds, carbon's preference for three (sp$^2$) or four (sp$^3$) depending on environment.

The Tersoff family includes Brenner's REBO (Reactive Empirical Bond Order, for hydrocarbons), AIREBO (Adaptive Intermolecular REBO, adding non-bonded interactions), and parameterisations for III-V semiconductors. Bond-order potentials can describe bond-breaking and bond-formation events qualitatively, which is why they are popular for studies of mechanical damage in covalent solids — crack tips, indentation, irradiation cascades. Quantitative chemistry (activation barriers, reaction enthalpies) is mostly beyond them.

## ReaxFF

ReaxFF (van Duin, Goddard et al., 2001) takes the bond-order idea to its logical extreme. Every interaction in the system — bond, angle, torsion, van der Waals, Coulomb — is modulated by continuously variable bond orders that respond to the local environment. Charges are redistributed self-consistently at every step via the **electronegativity equalisation method** (EEM) or the charge-equilibration (QEq) scheme, solving a charge-flow equation at each step:

$$
\frac{\partial E}{\partial q_i} = \chi_i + J_{ii} q_i + \sum_{j \ne i} J_{ij}(r_{ij}) q_j = \mu \quad \forall i,
\tag{7.49}
$$

with $\chi_i$ an electronegativity, $J_{ij}$ a hardness matrix, $\mu$ a chemical potential, and the charges constrained to sum to the total charge.

The upshot is a force field that can describe bond rearrangement: combustion of hydrocarbons, oxidation of metals, polymer cross-linking, the silica-water interface. ReaxFF parameterisations exist for most of the periodic table. The cost is steep: charge equilibration at every step adds an $O(N^2)$ piece (or $O(N \log N)$ with FFT-based variants), and a ReaxFF simulation is typically 30–100$\times$ slower per step than EAM or Tersoff.

ReaxFF is the right tool when chemistry happens during the simulation and you cannot afford DFT. It is the wrong tool for high-precision thermodynamics: parameterisations are typically accurate to a few kcal/mol for energies of reaction, and bond-breaking events are described qualitatively rather than to chemical accuracy.

## The accuracy-cost gap

Pull together what each method costs and what it delivers, on a per-atom-per-step basis for a typical bulk simulation in 2026:

| Method | Cost/atom/step | Typical accuracy | Transferability |
|---|---|---|---|
| Lennard-Jones | $\sim 1\,\mu s$ | Cartoonish | Universal in reduced units |
| EAM | $\sim 5\,\mu s$ | $\sim 10\%$ for metals | Within one parameterisation |
| Tersoff/REBO | $\sim 10\,\mu s$ | $\sim 10\%$ for covalent | Within one parameterisation |
| ReaxFF | $\sim 300\,\mu s$ | $\sim 5\,\mathrm{kcal/mol}$ | Better, fitted to chemistry |
| DFT (PBE) | $\sim 10\,\mathrm{ms}$ — $1\,\mathrm{s}$ | $\sim 50$ meV/atom | True ab-initio |
| MLIP (Chapter 9) | $\sim 0.1$–$10$ ms | $\sim 5$ meV/atom | As good as training data |

The classical force fields are five to ten orders of magnitude faster than DFT but lose roughly an order of magnitude in accuracy. The energy resolution of EAM for a copper system, when compared to DFT, is around 50–100 meV/atom — large enough to get the wrong stacking-fault energy, miss vacancy-cluster binding energies, and place phase boundaries at the wrong temperature. For chemistry — bond breaking, reaction barriers, charge transfer — classical force fields lose at least another order of magnitude.

This gap is what machine-learning interatomic potentials (MLIPs) aim to close. A neural network or kernel method trained on a few thousand DFT energies and forces can deliver near-DFT accuracy at roughly $10^3$–$10^4 \times$ lower cost. The result is that a 10000-atom simulation that would have cost months on a supercomputer with DFT now runs in hours on a workstation. Chapter 9 is devoted to this transition; Chapter 10 to graph neural network architectures and Chapter 11 to active-learning protocols that build training sets efficiently.

!!! note "When to skip MLIPs"
    If you are simulating bulk argon, do not train a neural network. Lennard-Jones is fine. If you are simulating defect kinetics in tungsten on microsecond timescales for fusion-reactor materials work, classical EAM is the only thing that runs fast enough, and the missing physics often does not matter at the level of trends. MLIPs are the right tool for the **middle** regime: hundreds to tens of thousands of atoms, where you need DFT accuracy and classical speed simultaneously.

## A practical sequence

If you are setting up a new material in MD for the first time:

1. **Search for an existing force field** in the LAMMPS distribution, the NIST IPR, or the OpenKIM project. For metals and oxides, parameterisations almost always exist. For exotic compositions, they do not.
2. **Validate against DFT** on a small training set: pick a dozen configurations of interest (relaxed crystal, dilated crystal, simple defects), compute energies and forces with both your force field and DFT, and compare.
3. **Decide whether the agreement is good enough.** "Good enough" depends on the question. For trends in elastic moduli across a series, 10% error is fine. For melting temperatures, 1% is needed.
4. **If not good enough, go to MLIPs (Chapter 9).** Reparameterising a classical force field by hand is a project that ate the careers of people in the 1990s; today it is rarely the right move.

## What we have

A complete recipe: integrator (§7.1), periodic cell (§7.2), thermostat or barostat (§7.3), force field (this section). We can now run a real simulation. Time to meet LAMMPS.
