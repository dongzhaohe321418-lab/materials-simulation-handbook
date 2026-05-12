# Coarse-Graining

QM/MM zooms *in*: it puts a high-accuracy method on a small region of
interest. Coarse-graining (CG) does the opposite — it zooms *out*, throwing
away atomic-scale degrees of freedom in order to reach time and length scales
that fully atomistic MD cannot. A 100 ns all-atom simulation of a lipid bilayer
is hard; a 100 μs coarse-grained simulation is routine.

This section walks through:

- Why CG works at all (and when it does not).
- The bottom-up procedures for building a CG model from atomistic data:
  iterative Boltzmann inversion, force matching (MS-CG), and relative entropy.
- Top-down models, with MARTINI as the canonical example.
- Machine-learned coarse-graining — a young and active field.
- Practical failure modes.

## Why coarse-grain?

The number of pairwise interactions in MD scales as $\mathcal{O}(N \log N)$
with cutoffs, but the *useful* dynamics often live on a much smaller manifold.
A lipid molecule has of order 100 atoms and dozens of soft dihedrals, but for
the purposes of membrane mechanics it is essentially a polar head plus a
non-polar tail. A polymer chain has thousands of atoms but its conformational
statistics are well captured by a few coarse beads representing chain
segments.

The dream of coarse-graining is to replace many fast degrees of freedom with a
few slow ones, and to do this in such a way that the slow dynamics is
*reproduced correctly*. The slow dynamics is, by hypothesis, the part you
actually care about.

There are two distinct payoffs:

- **Reduced cost per step.** Fewer beads, fewer forces, faster.
- **Larger integration timestep.** The fastest oscillations have been
  removed. A CG model with 10× longer timestep is a 10× speed-up on top of
  the bead reduction.

The combination can easily yield $10^3$ – $10^5$ speed-ups over atomistic MD.
That is enough to turn a year of compute into a weekend.

The catch — and there is always a catch — is that CG models do not
automatically inherit the physics of the underlying atomistic system. They
inherit *whichever quantities the CG fitting procedure was tuned to
reproduce*. A CG model that reproduces the radial distribution function may
have completely wrong dynamics. A CG model that reproduces forces may have
the wrong free energy. Choosing the right target is the whole game.

## Mapping: from atoms to beads

The first step of CG is *mapping*: defining the function that takes an
atomistic configuration and returns a CG configuration. For each CG bead
$I$, we define a set of atoms that map onto it, and a weighted average

$$
\mathbf{R}_I = \sum_{i \in I} c_{Ii}\, \mathbf{r}_i, \qquad \sum_{i \in I} c_{Ii} = 1
$$

The standard choices for $c_{Ii}$ are:

- **Mass-weighted**: $c_{Ii} = m_i / \sum_{j \in I} m_j$. The CG bead position
  is the centre of mass of its constituent atoms. This is the most physical
  choice if you want a Newtonian CG equation of motion.
- **Geometric centre**: $c_{Ii} = 1/|I|$. Easier to interpret but does not
  preserve the centre of mass.
- **Charge-weighted**: $c_{Ii} = q_i / \sum_{j \in I} q_j$. Used when
  electrostatics dominates.

Some quantities aggregate naturally: the mass of a CG bead is the sum of the
masses of its atoms, and the total charge of a CG bead is the sum of its
atoms' partial charges. Others do not: the moment of inertia of a CG bead is
*not* the sum of its atoms' moments of inertia, and CG dynamics often ignores
internal rotational degrees of freedom entirely.

Choices of mapping that you will see in practice:

- **Lipids (MARTINI).** Four heavy atoms per bead, roughly. A typical
  phospholipid maps to about 12 beads.
- **Polymers (Kremer-Grest).** A whole monomer per bead, sometimes 2–3
  monomers per bead for very long chains.
- **Proteins.** One bead per amino acid (centre of mass of backbone +
  representative side-chain bead) is common in the MARTINI scheme and its
  successors.
- **Water.** In MARTINI, four water molecules per bead. The discreteness of
  water is sacrificed for speed, with consequences we will return to.

The mapping is not unique. Different mappings of the same atomistic system
will give different CG models with different transferabilities and different
failure modes. The CG modeller's first big design choice is the mapping.

## Top-down CG: MARTINI

The MARTINI force field, developed by Marrink and collaborators, is the most
widely used CG force field for biomolecules. It is *top-down*: the CG
interaction parameters are calibrated against experimental thermodynamic
quantities (partition coefficients, surface tensions, densities) rather than
against atomistic simulations.

MARTINI's design philosophy is worth describing because it illustrates a
general principle.

Bead types are limited — about 18 types in MARTINI 2, roughly classified by
polarity (P), nonpolar (C), apolar (N), and charged (Q), with sub-classes.
Interactions between bead types are encoded in a $\sim 18 \times 18$ matrix
of Lennard-Jones parameters, each calibrated to reproduce the experimental
water/oil partition coefficient of a small molecule that the bead type is
meant to represent.

The advantages: extreme transferability across biomolecular systems, a
single canonical parameter set, a large community of users, well-validated
behaviour for membrane systems.

The disadvantages: chemically specific effects (hydrogen bonding patterns,
specific solvation shells, secondary structure of proteins) are not well
captured; the "four heavy atoms per bead" rule forces awkward compromises
for small molecules; the water model is famously simplistic, leading to
issues with the dielectric constant.

For materials simulation, MARTINI proper is less relevant than for biology,
but the *philosophy* — bead types, thermodynamic-property calibration,
limited transferability tested across many systems — has been imported into
materials CG (e.g. dissipative particle dynamics force fields for polymer
melts).

## Bottom-up CG: iterative Boltzmann inversion

Bottom-up CG fits the CG force field to *atomistic simulation data* rather
than to experimental properties. The most-used method for fitting CG pair
potentials is iterative Boltzmann inversion (IBI).

The setup: you have run a long all-atom MD simulation. You apply your
chosen mapping to every frame, and compute the *target* radial distribution
function (RDF) between CG bead types:

$$
g_\mathrm{target}(r)
$$

You wish to find a CG pair potential $U(r)$ such that running CG-MD with
$U(r)$ produces the same RDF.

The initial guess is the *direct Boltzmann inversion* of the target:

$$
U^{(0)}(r) = -k_B T \ln g_\mathrm{target}(r)
$$

This is exact if there are no many-body correlations; in practice many-body
effects are not negligible, so the CG-MD using $U^{(0)}$ produces an RDF
$g^{(0)}(r)$ that does not quite match $g_\mathrm{target}$.

Iterative Boltzmann inversion updates the potential:

$$
U^{(n+1)}(r) = U^{(n)}(r) + \alpha\, k_B T \ln \frac{g^{(n)}(r)}{g_\mathrm{target}(r)}
$$

where $\alpha \in (0, 1]$ is a step size. The procedure is repeated until
$g^{(n)}$ converges to $g_\mathrm{target}$, which typically takes 10–50
iterations.

IBI's strengths: simple, physically motivated, reproduces structure by
construction.

Its weaknesses: it reproduces only the *quantity it was fitted to* (the RDF).
If you also want correct pressure, you need a separate pressure correction
(VOTCA has this built in). If you want correct angular distributions, you
need to do IBI for angles as well. There is no guarantee that the CG model
generalises beyond the state point at which the atomistic reference was
generated.

## Bottom-up CG: force matching (MS-CG)

A more principled approach is *force matching*, also called multiscale
coarse-graining (MS-CG). Rather than match a derived quantity like the RDF,
match the forces themselves.

For each CG bead $I$, the atomistic system exerts a total force that maps
(via a back-mapping) to a "mean force" $\mathbf{F}_I^\mathrm{atomistic}$.
The CG force field gives a CG force $\mathbf{F}_I^\mathrm{CG}$. Force
matching minimises

$$
\chi^2[U_\mathrm{CG}] = \sum_I \big\langle |\mathbf{F}_I^\mathrm{atomistic} - \mathbf{F}_I^\mathrm{CG}|^2 \big\rangle
$$

over the CG potential parameters.

It can be shown that the minimum-$\chi^2$ CG force field is the *potential
of mean force* of the CG variables, projected onto the chosen functional
form (pair potentials, three-body terms, etc.). This is a thermodynamic
quantity that, in principle, contains exactly the right physics for the
slow CG modes.

Practical force matching:

1. Run atomistic MD; collect a large number of snapshots.
2. For each snapshot, compute the atomistic forces and map them to CG
   beads.
3. Set up a linear (or nonlinear, depending on basis) regression problem:
   express $\mathbf{F}_I^\mathrm{CG}$ as a sum over basis functions of the
   CG geometry, and solve for the coefficients by least squares.
4. Validate by running CG-MD and comparing not only forces but also
   structural and dynamic observables.

Force matching is more elegant than IBI and avoids the iteration loop. Its
downside is that it requires forces at every snapshot — not just
configurations — which roughly triples the I/O burden of the atomistic
reference simulation. It also tends to be more sensitive to the choice of
CG basis functions.

## Bottom-up CG: relative entropy minimisation

The most rigorous CG framework, due to Shell, Voth, and others, is
*relative entropy minimisation*. The relative entropy between the atomistic
ensemble (projected onto the CG variables) and the CG ensemble is

$$
S_\mathrm{rel} = \int p_\mathrm{atomistic-mapped}(\mathbf{R}) \ln \frac{p_\mathrm{atomistic-mapped}(\mathbf{R})}{p_\mathrm{CG}(\mathbf{R}; \theta)}\, d\mathbf{R}
$$

where $\theta$ parameterises the CG potential.

Minimising $S_\mathrm{rel}$ over $\theta$ gives the CG model that best
reproduces the *full Boltzmann distribution* of CG variables, not just one
moment of it. It can be shown that:

- Minimising $S_\mathrm{rel}$ with the CG distribution restricted to a
  Boltzmann distribution of a CG potential equals force matching in a
  certain limit.
- IBI is recovered if the CG potential is restricted to pair-additive form
  and the target is the pair RDF.

So $S_\mathrm{rel}$ is the unifying framework; IBI and MS-CG are its
projections onto particular function spaces. In practice $S_\mathrm{rel}$
minimisation is computationally heavy, and most production CG fitting still
uses IBI or MS-CG.

## Machine-learned CG

In the last few years, the idea of using neural networks for the CG
potential has matured. The motivation is obvious: the true potential of mean
force is a complicated many-body function of CG coordinates, far from any
simple pairwise or three-body decomposition. Why force it into a small
analytic form when we have universal approximators?

Several distinct lines of work exist:

**CGnet** (Husic, Charron, Noé and collaborators) trains a neural network
to reproduce the mean forces obtained from atomistic MD, in the spirit of
force matching but with a flexible function class. The network is trained
on snapshots of CG coordinates and the mapped mean forces.

**CG-MACE** and equivariant CG potentials extend the message-passing
equivariant architectures (see [Chapter 10](../ch10-gnn/index.md)) to CG.
The architectural choices that make MACE work well at the atomistic scale —
$O(3)$-equivariance, message passing over local neighbourhoods, body-order
control — translate directly to CG. The bonds connect CG beads instead of
atoms, but the maths is the same.

**Score-matching CG.** Rather than fit a force field directly, fit the
score (gradient of log-density) of the CG ensemble using a denoising
diffusion network. The CG potential is then implicitly the integrated
score. This sidesteps the choice of basis functions entirely.

**Coarse-graining with foundation models** is starting to be explored.
The idea: pre-train a large equivariant network on diverse atomistic data,
then fine-tune on CG mapping for a specific system. This is at an early
stage but is plausibly the future direction.

The shared challenge of all ML-CG methods is *transferability*. A CG model
fit to a single state point (one temperature, one composition) tends to
fail when extrapolated. The trick is to train on a diverse set of state
points, with active learning to identify configurations where the CG model
diverges from the atomistic reference.

!!! note "ML-CG is not free physics"
    Even with a perfect machine-learned CG potential, you still face the
    fundamental limitations of any CG approach: chemically specific events
    that require atomic resolution cannot be recovered from the CG
    coordinates alone. Hydrogen bonding, charge transfer, and electronic
    polarisation must either be built into the CG mapping somehow, or
    accepted as out of scope.

## Dynamics in CG: a subtlety

A CG force field, however well-fitted, does not automatically give correct
*dynamics*. The reason is statistical mechanics, not numerics.

When you average away fast degrees of freedom to get a CG potential of mean
force, the resulting Hamiltonian dynamics on CG coordinates is in general
*not* the dynamics that the CG coordinates would have shown in the
atomistic system. The correct equation of motion is the *generalised
Langevin equation*:

$$
M\ddot{\mathbf{R}} = -\nabla U_\mathrm{PMF}(\mathbf{R}) - \int_0^t K(t-t') \dot{\mathbf{R}}(t')\, dt' + \boldsymbol{\eta}(t)
$$

with a memory kernel $K$ and noise $\boldsymbol{\eta}$ satisfying a
fluctuation-dissipation relation.

In practice, CG-MD is almost always run with a simple thermostat (Langevin
or Nosé-Hoover) and the memory kernel is ignored. This gives the right
*equilibrium* statistics — the CG configurations have the correct Boltzmann
distribution — but the *timescales* are wrong. CG dynamics is typically
faster than the underlying atomistic dynamics by a factor of 4 to 10. Some
authors recalibrate time by matching a single timescale (e.g. lipid lateral
diffusion); others quote results in "CG units" without claiming a mapping
to real time.

If you only care about thermodynamics, ignore this. If you care about
kinetics, be careful.

## When CG fails

Coarse-graining is a powerful technique, but it has well-defined failure
modes. You should know them.

### Chemically specific features

Hydrogen bonding has a directional, angular character that is hard to
capture with isotropic CG beads. Solvation shells around small ions have a
specific structure that is washed out when you replace four water molecules
with one bead. Reactive chemistry (bond making and breaking) is out of
reach by construction.

The cure is either to use a CG model that builds the relevant directionality
into the bead types (MARTINI's distinction between donor and acceptor beads
is a primitive version of this), or to give up and use atomistic MD.

### Anomalous diffusion

If your system has *glassy* or *sub-diffusive* dynamics — many materials
near a transition, polymers above $T_g$, ionic liquids — then the CG model
typically underestimates the dynamical heterogeneity. The fast modes that
you averaged away include precisely the rare events that gate slow dynamics.

### Pressure and density

CG models fit by IBI famously give the wrong pressure unless you add a
specific pressure correction. The reason is that the RDF does not fully
determine the equation of state; the same RDF can correspond to different
pressures for different potentials. If you care about the equation of state
(constant-pressure simulations, mechanical properties), use a method that
matches pressure as well, such as IBI with iterative pressure correction
or the conditional reversible work method.

### Transferability across temperature

A CG potential fit at 300 K may give wrong results at 400 K because the
free-energy surface of the underlying atomistic system depends on
temperature in a way the CG potential does not. The cure is to fit at
multiple temperatures simultaneously, or to make the CG potential
temperature-dependent explicitly.

### Cross-validation hygiene

The cardinal sin in CG, as in any data-driven modelling, is to validate
on the same configurations you trained on. The CG model will reproduce
the training RDF perfectly while failing badly on independent test
ensembles. Always set aside a held-out atomistic simulation as a test
set, ideally at a different state point.

## Summary

Coarse-graining is the natural complement to QM/MM: where QM/MM zooms in,
CG zooms out. Both are about choosing the right resolution for the
physics you care about.

For your own projects, the heuristic is:

- If you can afford fully atomistic MD with a good force field or MLIP,
  use it. CG is an optimisation; do not pay its costs unnecessarily.
- If you cannot, identify the slow modes you care about and design a
  mapping that preserves them.
- Always validate the CG model against atomistic data at a state point you
  did not train on.

We turn next to two methods that bridge the gap between atomistic and
continuum without explicitly representing every degree of freedom: kinetic
Monte Carlo and phase-field modelling.
