# 3.2 Bonding

A solid exists because atoms stick to each other. *How* they stick determines almost everything about its macroscopic properties: whether it is hard or soft, whether it conducts electricity, whether it melts at 100 °C or 3000 °C, whether it dissolves in water or shrugs at concentrated acid. The chemical bond is the workhorse concept of materials science, and the qualitative bonding picture is a useful guide long before one writes down any equations.

This section walks through the five canonical bond types — metallic, covalent, ionic, van der Waals, and hydrogen — with a focus on what each predicts about properties and on how well each is captured by the different simulation methods.

## What bonding is, physically

At the quantum level, a chemical bond is simply the lowering of the total energy of a collection of atoms when they are brought close together. The lowering can arise in several distinct ways, and these ways are what distinguish bond types:

- **Sharing of electrons** between two or more atoms in delocalised states, lowering kinetic energy (covalent and metallic).
- **Charge transfer** from one atom to another, with the resulting ions held together by electrostatic attraction (ionic).
- **Fluctuation-induced dipole interactions** between neutral, closed-shell species (van der Waals / dispersion).
- **Electrostatic attraction between partial charges**, often involving a hydrogen sitting between two electronegative atoms (hydrogen bonding).

In a real solid these mechanisms coexist. A nominally ionic crystal like NaCl has a small covalent component; a metallic crystal like copper has tiny dispersion contributions at the surface. The taxonomy below is a *first* description: helpful for orientation, not the last word.

## Metallic bonding

In a metal, the valence electrons are not associated with any particular atom. They form a delocalised *sea* moving in the periodic potential of the positive ion cores. The system is held together because spreading the electrons over the whole crystal lowers their kinetic energy (the larger the box, the lower the zero-point kinetic energy of a quantum particle).

**Energy scale.** Cohesive energies of $1$–$10$ eV per atom. Strong bonds: 1 eV/atom corresponds to a melting point of order 1000 K.

**Length scale.** Nearest-neighbour distances of 2.4–3.2 Å for common metals. Coordination numbers of 8 (BCC) or 12 (FCC, HCP) — close packing dominates.

**What it predicts.**

- **High electrical conductivity.** Delocalised electrons respond to applied fields. Metals are good conductors.
- **High thermal conductivity.** Electrons carry heat as well as charge.
- **Ductility and malleability.** The bond is non-directional: planes of atoms can slide past each other without breaking individual bonds. Metals deform plastically.
- **Metallic lustre.** Free electrons reflect light in the visible range.
- **Moderate to high Young's modulus.** Stiff but not as stiff as the strongest covalent solids.

**Examples.**

- *Sodium* — soft, low melting point (98 °C), one valence electron per atom. The textbook *nearly-free electron* example.
- *Copper* — much harder, melts at 1085 °C, eleven valence electrons including the full $3d$ shell. Excellent conductor.
- *Tungsten* — refractory, melts at 3422 °C, strong $d$-electron contribution. The most refractory metal.

**Simulation considerations.** Metallic bonding is the easy case for both DFT and classical potentials. The embedded-atom method (EAM) is a successful classical scheme tailored to metals and gives reasonable elastic constants and surface energies for cheap. DFT handles metals well — but they require dense k-point sampling because the Fermi surface has structure (more on this in Chapter 6). The main subtlety is in *transition metals*, where $d$-electron correlation can make standard DFT subtly inaccurate; modern functionals (SCAN, r$^2$SCAN) help.

## Covalent bonding

In a covalent bond, electrons are shared between two atoms in localised molecular orbitals. The bond is highly directional, set by the geometry of the orbitals involved. Carbon's tetrahedral geometry — every bond at 109.5° — is the canonical example, set by the sp$^3$ hybridisation of the four valence orbitals.

**Energy scale.** Cohesive energies of $3$–$8$ eV per atom — comparable to metals or slightly higher. Individual bond energies (i.e., the energy required to break one specific bond) range from $2$ eV for a weak single bond to $9$ eV for a triple bond.

**Length scale.** Nearest-neighbour distances 1.0–2.5 Å, often shorter than metallic bonds because the shared electrons are concentrated between the atoms. Coordination numbers are *low* (typically 2 to 6) because each atom has a limited number of bonds it can form, set by its valence.

**What it predicts.**

- **High hardness and stiffness.** Diamond is the canonical example, with a Young's modulus of $\sim 1$ TPa, an order of magnitude above copper.
- **High melting points.** Diamond sublimes around 3550 °C (without ambient pressure stabilisation); silicon melts at 1414 °C.
- **Low electrical conductivity at low temperature** (in pure form). Covalent semiconductors and insulators have band gaps; carriers must be thermally or optically excited.
- **Brittleness.** Bond directionality means dislocations cannot easily glide; the crystal fractures rather than deforming plastically.

**Examples.**

- *Diamond* — covalent carbon, sp$^3$-hybridised, the prototypical hard insulator. Band gap $5.5$ eV.
- *Silicon* — same diamond structure but with a smaller gap (1.1 eV) and weaker bonds. The basis of microelectronics.
- *Silicon dioxide* — partially covalent SiO$_4$ tetrahedra; the structural basis of glass, sand, and zeolites.

**Simulation considerations.** Covalent bonding is well handled by DFT and by good classical potentials (Tersoff, Stillinger–Weber, ReaxFF for reactive systems). Modern MLIPs trained on DFT data handle covalent systems excellently. The main subtleties involve aromatic systems (where simple potentials miss conjugation) and bond-breaking reactions (where bond order changes and a non-reactive potential fails).

## Ionic bonding

When an electropositive atom (like a Group 1 metal) and an electronegative atom (like a Group 17 halogen) come together, the metal donates an electron to the non-metal. The result is two ions of opposite charge held together by electrostatic attraction.

**Energy scale.** Cohesive energies $3$–$8$ eV per ion pair. The electrostatic energy scales as $Q_1 Q_2 / r$, which is large and slowly decaying.

**Length scale.** Inter-ionic distances 2.0–3.5 Å, set by the sum of ionic radii. Coordination numbers vary widely with the radius ratio: 6 for NaCl, 8 for CsCl, 4 for ZnS.

**What it predicts.**

- **High melting points.** NaCl melts at 801 °C; MgO at 2852 °C.
- **Low electrical conductivity at low temperature**, high at high temperature (ions become mobile in the melt or via vacancies).
- **Brittleness.** Like covalent crystals, ionic crystals shatter; the ions cannot easily rearrange without bringing like charges into contact.
- **Solubility in polar solvents.** Water dissolves NaCl by solvating Na$^+$ and Cl$^-$ ions individually.

**Examples.**

- *Sodium chloride (NaCl)*, the canonical rock salt. Na$^+$ and Cl$^-$ in an FCC arrangement with 6:6 coordination.
- *Magnesium oxide (MgO)*. Same rock salt structure but with doubly charged ions, hence dramatically higher cohesive energy and melting point.
- *Caesium chloride (CsCl)*. Different structure because the radius ratio favours 8:8 coordination.

**Simulation considerations.** The long-range Coulomb interaction is the central technical challenge. Naive truncation of $1/r$ is *not* convergent; one uses Ewald summation or particle-mesh methods. Classical potentials for ionic crystals (e.g. the Buckingham potential plus point charges) are accurate for simple ionic solids but require careful treatment of polarisation effects (the shell model captures these). DFT handles ionic systems robustly, with the caveat that strong electron correlation in some oxides (NiO, FeO) requires beyond-GGA treatment.

## Van der Waals (dispersion) bonding

When two closed-shell atoms or molecules — argon atoms, two benzene rings — approach each other, they polarise each other's electron clouds. The instantaneous dipole on one induces a correlated dipole on the other, and the resulting interaction is weakly attractive. This is the *London dispersion force*, the universal background interaction between any two molecules.

**Energy scale.** Tiny: $1$–$10$ meV per atom pair, three orders of magnitude weaker than covalent bonds. The London dispersion potential scales as $-C_6 / r^6$ at large separations.

**Length scale.** Equilibrium distances of $3.5$–$5$ Å between atoms. Larger than chemical bonds because nothing pulls the electron clouds in.

**What it predicts.**

- **Very low melting and boiling points.** Argon melts at $-189$ °C; benzene at $5$ °C. Both are held together only by dispersion.
- **Low elastic moduli.** Soft crystals.
- **Layered structures, in some cases.** Graphite is held by strong covalent bonds within layers and weak dispersion between layers; the weak interlayer coupling is why graphite is a good lubricant.

**Examples.**

- *Argon crystal* — FCC, melts at 84 K. Pure dispersion, simplest possible test case.
- *Molecular crystals* — naphthalene, anthracene, the molecular pharmaceutical solids. Within-molecule bonds are covalent; between-molecule binding is dispersion.
- *Graphite, hexagonal BN, MoS$_2$* — layered materials with covalent or polar-covalent in-plane bonds and dispersion-dominated interlayer binding. The basis of two-dimensional materials.

**Simulation considerations.** *Standard DFT (LDA, PBE, even hybrid functionals) does not capture dispersion.* Local and semi-local exchange-correlation functionals see only the local electron density and miss the long-range, non-local correlation that gives rise to dispersion. The standard fix is to add an explicit dispersion correction:

- **DFT-D3** (Grimme) — adds a pairwise $-C_6 / r^6$ term with empirically tabulated coefficients. Cheap and routine.
- **DFT-D4** — refines D3 with charge-dependent $C_6$ coefficients.
- **TS** (Tkatchenko–Scheffler) and **MBD** (many-body dispersion) — derive $C_6$ from the actual electron density and, for MBD, include beyond-pairwise effects.
- **Non-local correlation functionals** (vdW-DF1, vdW-DF2, rVV10) — incorporate dispersion at the functional level.

For any system where layers, molecules, surfaces, or noble gases matter, dispersion corrections are essential. A DFT calculation of graphite without vdW correction predicts no interlayer binding at all. Modern MLIPs, when trained on dispersion-corrected DFT data, inherit the corrections automatically.

## Hydrogen bonding

Hydrogen bonding is a special case, important enough to its own category. When a hydrogen atom is covalently bound to an electronegative atom (N, O, F), the H becomes substantially positive (the electronegative partner pulls the bonding electrons away). This positive H can then form an electrostatic interaction with another electronegative atom nearby: $\mathrm{X{-}H} \cdots \mathrm{Y}$ where X, Y are electronegative.

**Energy scale.** $50$–$300$ meV per H-bond ($\sim 5$–$30$ kJ/mol), intermediate between van der Waals and covalent. Strong enough to determine structure, weak enough to be broken and reformed at room temperature.

**Length scale.** $1.5$–$2.0$ Å for the H$\cdots$Y distance (the *acceptor* distance), with the X–H covalent bond around 1.0 Å on the *donor* side.

**What it predicts.**

- **Anomalously high boiling points** for hydride compounds. Water boils at 100 °C; H$_2$S at $-60$ °C, despite sulfur being heavier than oxygen.
- **Open structures.** Ice has a tetrahedral hydrogen-bonded network that is less dense than liquid water — why ice floats.
- **Specificity.** The directionality of the H-bond is what makes DNA's base pairing specific.

**Examples.**

- *Ice and water* — the canonical hydrogen-bonded system. Each H$_2$O molecule donates two and accepts two H-bonds in ideal ice Ih.
- *DNA double helix* — base pairs are held together by H-bonds (two between A–T, three between G–C).
- *Carboxylic acid dimers* — two H-bonds between matching donor and acceptor groups.

**Simulation considerations.** Hydrogen bonds are partly electrostatic (well captured by ordinary DFT) and partly covalent in nature (a small charge-transfer component). Standard GGA-DFT does reasonably well on H-bond geometries but tends to overbind compared to experiment. Hybrid functionals do better. The most challenging aspect is *quantum nuclear effects*: the hydrogen is light enough that zero-point motion is comparable to the H-bond length, and proton tunnelling can be relevant. Path-integral molecular dynamics is the gold standard, used routinely in water research.

## Real systems mix bond types

Few real materials are purely one bond type. The pure cases are pedagogical idealisations.

- *Silicon carbide (SiC)* — partly covalent, partly ionic. The Si–C bond is covalent in character but has a small charge transfer from Si to C because C is more electronegative.
- *Transition-metal oxides (e.g., TiO$_2$)* — mixed ionic and covalent. The Ti–O bond has substantial ionic character but the $d$-electron interactions add covalency.
- *Metal–organic frameworks* — metallic clusters connected by covalent organic linkers, with the assembly held together by a combination of covalent and ionic motifs at the cluster–linker interface.
- *Proteins in water* — covalent backbones, hydrogen-bonded secondary structure, ionic interactions between charged residues, and dispersion-dominated hydrophobic effects.

A simulation method has to handle whichever bond types dominate in the system at hand. Choosing a force field for an MOF on the basis of its metal node alone, ignoring the covalent linker, will fail. The skill is in identifying which bond types are present *and* which of them carry the physics you want.

## How methods handle each bond type

Synthesising the preceding discussion into a table:

| Bond type | DFT (GGA) | DFT + vdW | Hybrid DFT | MLIP (foundation) | Classical (specialised) |
|---|---|---|---|---|---|
| Metallic | Excellent | Excellent | Excellent | Excellent | Excellent (EAM) |
| Covalent | Good | Good | Excellent | Excellent | Excellent (Tersoff, SW) |
| Ionic | Excellent | Excellent | Excellent | Excellent | Good (Buckingham + shell) |
| van der Waals | **Poor** | Excellent | **Poor** | Good (if trained) | Excellent (LJ) for vdW only |
| Hydrogen | Good | Good | Excellent | Excellent | Modest (TIP4P, etc.) |

The most important row is *van der Waals*: standard DFT misses it, and any system whose binding is dispersion-dominated requires either a correction (DFT-D3, MBD) or a non-local functional. Forgetting this is the single most common DFT mistake.

The most important column is *MLIP*: a foundation MLIP trained on dispersion-corrected DFT inherits the dispersion. The same MLIP, trained on vanilla PBE without correction, will not. The bond-type table is, in this sense, a guide to what your training data has to include.

## A computational view

A useful exercise is to compute approximate binding energies of small dimers, one per bond type, and compare. Here we will not run the calculations live, but the conceptual numbers (in eV per pair, equilibrium distance in Å) are:

- Na$_2$ (metallic, weak): $0.74$ eV, $3.08$ Å.
- C$_2$ (covalent, strong): $6.2$ eV (triple bond, in fact), $1.24$ Å.
- NaCl (ionic, gas-phase dimer): $4.3$ eV, $2.36$ Å.
- Ar$_2$ (vdW): $0.012$ eV, $3.76$ Å.
- (HF)$_2$ (hydrogen-bonded): $0.20$ eV, $1.83$ Å (F$\cdots$H), $2.79$ Å (F$\cdots$F).

These differ by *three orders of magnitude* in binding strength. A method that gets the metallic and covalent cases right but predicts zero binding for Ar$_2$ is not wrong about chemistry — it is just missing one of the channels. Knowing which channel is which is the point of this section.

!!! tip "Diagnostic question"
    Before starting any simulation, ask: *which bond types are present in my system, and which dominate the property I am trying to compute?* If van der Waals dominates, you need dispersion corrections. If covalent bond-breaking is involved, you need a reactive potential or DFT. If the property is hardness or stiffness, the strongest bonds in the system dominate; if the property is melting point, the weakest bonds often do.

## Looking forward

The bonding picture set up here will reappear at every level of the book. In Chapter 4, bonding emerges from the molecular-orbital picture of two-atom systems. In Chapter 5, DFT lets us compute it from first principles. In Chapter 7, classical force fields encode specific bond types in their functional forms. In Chapter 9 and beyond, machine-learning potentials *learn* the bonding from data — and the relevant question becomes whether the training set contains enough of each bond type.

With the bonding picture in hand, we can now ask how atoms organise themselves once there are many of them: crystals.
