# The Finite-Element Bridge

The last rung of the multiscale ladder is the finite-element method (FEM),
which lives at the engineering scale: components, structures, devices. FEM
has been the workhorse of mechanical, civil and aerospace engineering for
fifty years. Its inputs — material parameters such as the Young's modulus,
the yield stress, the thermal conductivity — used to come from experimental
testing. Increasingly, they come from atomistic and DFT simulations.

This section is about the bridge from atomistic simulation to FEM:
*atomistically informed finite-element analysis*. The argument that this
matters is straightforward. Designing a turbine blade for the next generation
of aircraft engines requires a material that performs reliably at conditions
where exhaustive experimental characterisation is expensive and slow.
Predicting the material's behaviour from electronic structure, with FEM as
the engineering bridge, dramatically shortens the design cycle.

## What FEM does, very briefly

FEM solves boundary-value problems for partial differential equations on a
discretised geometry. For solid mechanics, the equation is the equilibrium
or momentum equation,

$$
\nabla \cdot \boldsymbol{\sigma} + \mathbf{f}_\mathrm{body} = \rho \ddot{\mathbf{u}}
$$

with $\boldsymbol{\sigma}$ the stress tensor, $\mathbf{u}$ the displacement
field, and $\mathbf{f}_\mathrm{body}$ a body force such as gravity. Closure is
provided by a *constitutive law* relating stress to strain. For linear
elasticity,

$$
\boldsymbol{\sigma} = \mathbf{C}\, \boldsymbol{\varepsilon}
$$

with $\mathbf{C}$ the fourth-rank elastic stiffness tensor and
$\boldsymbol{\varepsilon}$ the strain. For plasticity, more complicated.

The numerical method discretises the geometry into elements (tetrahedra,
hexahedra, etc.), interpolates the displacement field with shape functions,
and reduces the PDE to a large linear system

$$
\mathbf{K} \mathbf{u} = \mathbf{f}
$$

with $\mathbf{K}$ the global stiffness matrix assembled from element-level
contributions that depend on $\mathbf{C}$ and the mesh geometry. Solving
this system gives the displacement field; differentiating gives strains and
stresses everywhere in the component.

The detailed machinery of FEM is the subject of entire textbooks. For our
purposes the key observation is this: **FEM needs material parameters as
input.** The job of an atomistic simulation, in the context of FEM, is to
supply those parameters.

## The simplest bridge: elastic constants

The most established atomistic-to-FEM bridge is the computation of *elastic
constants*. The stiffness tensor $\mathbf{C}$ for a single crystal has, in
general, up to 21 independent components, reducing by symmetry to as few as
2 (isotropic), 3 (cubic) or 5 (hexagonal).

The atomistic recipe is:

1. Optimise the crystal at zero stress. Record the equilibrium volume
   $V_0$ and energy $E_0$.
2. Apply a small *strain* $\boldsymbol{\varepsilon}$ to the lattice
   vectors. For computing $C_{11}$, a uniaxial strain
   $\varepsilon_{11} = \pm 0.005$, with all other strains zero.
3. Re-relax the *internal* atomic positions at the strained cell. Compute
   the new energy $E(\boldsymbol{\varepsilon})$.
4. Fit the energy as a function of strain:

   $$
   E(\boldsymbol{\varepsilon}) = E_0 + V_0 \sum_{ij} \sigma_{ij}^{(0)} \varepsilon_{ij} + \frac{V_0}{2} \sum_{ijkl} C_{ijkl} \varepsilon_{ij} \varepsilon_{kl} + \ldots
   $$

   At a stress-free reference, the linear term vanishes; the elastic
   constants are read off from the quadratic term.

This is a textbook DFT calculation (or MLIP, or MD-with-fluctuations
calculation), and many automated packages handle it. AiiDA and Pymatgen
both have built-in elastic-constant workflows.

The result is the single-crystal stiffness $\mathbf{C}^\mathrm{sc}$. But
engineering parts are not single crystals; they are polycrystalline. The
next step is averaging.

## Polycrystalline averaging: Voigt, Reuss and Hill

A polycrystal is an aggregate of many crystallites with random orientations.
The effective elastic stiffness of the polycrystal depends on how stress and
strain are distributed among the crystallites — which depends on the
microstructure, which we generally cannot compute exactly.

Two simple bounds, the Voigt and Reuss averages, give an idea of the range.

The **Voigt average** assumes uniform strain across all crystallites. It is
an upper bound:

$$
C^\mathrm{V}_{ij} = \frac{1}{N} \sum_n C^{(n)}_{ij}
$$

(in the Voigt notation that flattens the fourth-rank tensor into a $6
\times 6$ matrix), with the sum over crystallite orientations.

The **Reuss average** assumes uniform stress and is a lower bound:

$$
S^\mathrm{R}_{ij} = \frac{1}{N} \sum_n S^{(n)}_{ij}
$$

with $\mathbf{S}$ the compliance ($\mathbf{S} = \mathbf{C}^{-1}$).

The **Hill (VRH) average** takes the arithmetic mean of the Voigt and Reuss
estimates:

$$
C^\mathrm{VRH} = \frac{1}{2}\left(C^\mathrm{V} + C^\mathrm{R}\right)
$$

The VRH average is what most "predicted bulk modulus" or "predicted shear
modulus" numbers in computational materials databases (Materials Project,
NOMAD, AFLOW) refer to. It is not exact — the true polycrystal modulus
depends on texture and microstructure — but it is a robust estimate that
takes only a single-crystal stiffness tensor as input.

From the VRH-averaged stiffness one extracts the engineering moduli:

$$
B = \frac{1}{9}\left(C_{11} + C_{22} + C_{33} + 2 C_{12} + 2 C_{13} + 2 C_{23}\right) \quad \text{(bulk modulus, Voigt form)}
$$

The shear modulus $G$ follows from a similar combination, and Young's
modulus and Poisson's ratio from

$$
E = \frac{9 B G}{3 B + G}, \quad \nu = \frac{3 B - 2 G}{2 (3 B + G)}.
$$

These two scalars, $E$ and $\nu$, are what an *isotropic* linear-elastic
FEM calculation needs.

## A case study: alloy strength

Let us walk through a complete bridge for an engineering question: predict
the elastic response of a steel component under load, starting from
electronic structure.

### Step 1 — DFT of the constituents

Choose a model alloy: bcc iron with 2% chromium. Run DFT on a 2x2x2
supercell of bcc iron with one Cr atom in place of one Fe. Optimise
geometry. Apply small strains in the six independent directions of the
Voigt notation. Compute $\mathbf{C}$. The result is a $6 \times 6$
single-crystal stiffness for this composition.

Practical notes: this is a metallic system, so smearing is essential
(see [Section 5 of the capstone](../ch14-capstone/05-common-pitfalls.md)).
The supercell must be large enough that the periodic images of the Cr atom
do not interact strongly. K-point convergence must be checked carefully.

### Step 2 — Voigt-Reuss-Hill averaging

Apply the VRH formula to the stiffness tensor of step 1. Extract bulk
modulus $B$, shear modulus $G$, Young's modulus $E$, Poisson's ratio $\nu$.

For 2% Cr in bcc Fe, you expect $E \approx 210$ GPa, $\nu \approx 0.29$,
similar to pure iron because Cr is similar in size and bonding to Fe.

### Step 3 — FEM of a tensile specimen

Build a 3D geometry of a standard tensile test specimen (the "dog-bone"
shape: a parallel gauge section of known cross-section, with grip regions
at each end). Mesh it into roughly $10^5$ tetrahedra.

Set the material model: linear elastic, isotropic, with $E = 210$ GPa,
$\nu = 0.29$.

Apply boundary conditions: fix one grip, apply a tensile displacement to
the other grip. Solve $\mathbf{K} \mathbf{u} = \mathbf{f}$.

The output is the stress and strain everywhere in the specimen. From the
gauge section, you can read the engineering stress-strain curve. The slope
gives Young's modulus (which we already know — it was our input). The
*non-trivial* outputs are stress concentrations at fillets, the spatial
distribution of stress at any given applied load, and (with plasticity
added) where yielding starts and how it propagates.

### Step 4 — extending to plasticity and to mixed loading

Linear elasticity is rarely enough for engineering. To predict yielding and
failure, you need a *plasticity model*. The standard choice is J2 plasticity
(von Mises) with a hardening law, which adds parameters: yield stress
$\sigma_y$, hardening exponent $n$, strain to failure $\varepsilon_f$.

These parameters come from a combination of:

- *Single-crystal plasticity* simulations using dislocation dynamics, which
  use DFT-computed core energies and Peierls barriers as inputs.
- *MLIP-MD* of polycrystalline aggregates under load, extracting effective
  yield surfaces.
- *Experimental calibration* — because at the end of the day, plasticity is
  hard and an experimental data point goes a long way.

This is where the multiscale chain gets long and the error budget gets
large. Predicting elastic constants from DFT is robust; predicting the
strain to failure of a polycrystalline alloy from DFT is, at present, an
active research problem.

## Implementation: getting hands on it

For an actual workflow, the common open-source tools are:

- **Atomistic side**: VASP, Quantum ESPRESSO, GPAW for DFT, with Pymatgen or
  AiiDA orchestrating the elastic-constant calculation. The Pymatgen
  `ElasticTensor` class implements the strain protocol and VRH averaging.
- **FEM side**: FEniCS (Python-native, great for prototyping), MOOSE
  (Idaho National Lab's framework, multi-physics, well-suited to coupling
  with materials simulation), Code_Aster (open-source, mature). Commercial
  tools: Abaqus, Ansys.
- **Bridge**: simply file-based. The FEM tool takes $E$, $\nu$ as input
  numbers, or a full $\mathbf{C}$ tensor if you set up an anisotropic
  material. The trick is reproducibility: track the DFT version, the
  k-point grid, the cutoff, the cell, the strain magnitudes — every
  setting that fed the FEM.

## Higher-fidelity bridges

The elastic-constant bridge is one-dimensional: scalar material parameters
fed into a standard continuum law. More sophisticated bridges pass
*structured data* between scales.

**Crystal plasticity FEM (CPFEM).** Instead of an isotropic continuum, the
FEM uses an anisotropic constitutive law that explicitly models slip on
crystallographic systems. Each integration point in the FEM has a local
crystal orientation; the constitutive law depends on the orientation. The
slip-system-level parameters — Peierls barriers, hardening coefficients —
come from atomistic simulation.

**Microstructure-aware FEM.** The mesh resolves grains explicitly, with each
grain assigned an orientation. The input is a synthetic microstructure
from a phase-field simulation (see [previous section](04-kmc-phase-field.md))
or from a CT scan of a real specimen. This is heavy but is the gold standard
for predicting how microstructure affects component behaviour.

**Concurrent atomistic-FEM (CADD, bridging domain).** A truly *concurrent*
coupling where a small region of the FEM mesh is replaced by atomistic
simulation. The handshake region is, again, where the difficulties live —
see [Section 1](01-coupling-schemes.md) for the general issues. Used for
crack-tip mechanics and dislocation-grain-boundary interactions; rarely used
in routine engineering analysis.

## Limitations: what FEM cannot see

FEM is the engineering scale. It excels at predicting how a *given material*
responds to *known loads in a known geometry*. It does not, by construction,
see:

- **Dislocation cores.** The constitutive law of plasticity is an averaged
  representation of dislocation behaviour; it does not know about cores.
- **Defect chemistry.** Hydrogen at a grain boundary, oxygen segregation,
  point-defect kinetics — all atomistic phenomena that FEM treats only
  through effective parameters in its constitutive law.
- **Phase transformations.** Unless explicitly built in as a phase-field
  coupled subproblem, FEM will not nucleate a new phase.
- **Surface and interface effects below the mesh resolution.** FEM with mm-
  scale elements cannot resolve a μm-thick interfacial layer.

This is not a flaw of FEM; it is a feature of *any* continuum description.
The job of atomistic and mesoscale simulation is to provide effective
parameters that capture these effects in the FEM's continuum language.

!!! warning "Don't fool yourself with FEM"
    A common failure mode is to take a sleek FEM simulation, with colourful
    von Mises stress contours, as more reliable than it is. The FEM is *only
    as good as its constitutive law*. Garbage parameters in, beautiful
    garbage out. When you see a multiscale paper that claims to predict
    component-level engineering properties from DFT alone, look hard at the
    constitutive law. Almost always, experimental calibration has crept in
    somewhere — sometimes acknowledged, sometimes not.

## Where this leaves us

The multiscale staircase, viewed in full:

| Scale | Method | Output | Time |
|---|---|---|---|
| Sub-Å, fs | DFT | Energy, forces, electronic structure | ps |
| Å-nm, ps-ns | MD with MLIP | Free energies, dynamics | ns-μs |
| nm-μm, μs-ms | KMC, CG-MD | Rates, microstructure | μs-s |
| μm-mm, ms-s | Phase field | Microstructure evolution | s-h |
| mm-m, s-y | FEM | Component response | as needed |

Information flows up the table: each method's output feeds the next.
Occasionally information flows back down — when a coarser-scale boundary
condition modifies a finer-scale problem — and that is when concurrent
coupling enters.

No single project will use all five rows. A practical PhD or undergraduate
thesis often touches two — atomistic + KMC, or DFT + FEM — and gets a great
deal of value out of the connection. The skill of a multiscale practitioner
is recognising *which two rows your problem actually requires*, and
resisting the temptation to add more for show.

We have now reached the end of the multiscale chapter. The
[exercises](exercises.md) give you a chance to think about coupling
choices, error budgets, and where the seams of a multiscale workflow are
likely to tear.

[Chapter 14](../ch14-capstone/index.md) turns from techniques to
*projects*: how to design your own undergraduate or master's thesis using
the methods of this handbook.
