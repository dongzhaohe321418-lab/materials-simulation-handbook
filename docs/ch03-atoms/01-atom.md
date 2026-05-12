# 3.1 The Atom — What We Keep and What We Throw Away

An atom is, in the model we will use throughout this book, a positively charged nucleus surrounded by negatively charged electrons. The nucleus is approximately $10^5$ times smaller than the electron cloud, weighs roughly $10^4$ times as much as all the electrons combined, and — for materials science purposes — moves slowly enough that we can usually treat its motion separately from that of the electrons. This separation, the Born–Oppenheimer approximation, is the first and most consequential simplification of computational materials science. We will return to it in Chapter 4. For now, take it as given: the nucleus has a position, a charge, a mass, and otherwise sits there.

The interesting structure of the atom — and the interesting variation between elements — lives in the electrons. This section is about what we do with them.

## What we keep

In a full first-principles description, every electron is a quantum particle. Its state is part of a many-electron wavefunction $\Psi(\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_N)$, antisymmetric under exchange of any two electrons because they are fermions, normalised to unity, and obeying the Schrödinger equation

$$
\hat{H} \Psi = E \Psi,
$$

where $\hat{H}$ contains the electrons' kinetic energy, their attraction to the nuclei, and their repulsion from each other. Solving this for more than two or three electrons is impossible analytically, and exact numerical solution is too expensive for systems larger than a small molecule. The history of electronic structure theory is largely the history of approximations to this equation that retain enough physics to be useful.

The single most important approximation we make is the **orbital picture**: we describe the electrons as occupying one-particle states (orbitals) and build approximate many-body wavefunctions from these orbitals. The orbital picture is exact for non-interacting electrons (it isn't fully right when electrons interact) but it survives, in modified form, as the language of essentially all practical electronic-structure methods. In density functional theory, the orbitals are the Kohn–Sham orbitals; in Hartree–Fock, they are the canonical molecular orbitals; in chemistry textbooks, they are simply *the orbitals*. The names differ but the concept is the same.

## The orbital picture

For an isolated atom, the orbitals are well-known special functions called *hydrogen-like orbitals* (or, with corrections, *Slater* and *Hartree–Fock orbitals*). Each orbital is labelled by three quantum numbers: the principal quantum number $n = 1, 2, 3, \ldots$, the orbital angular momentum quantum number $\ell = 0, 1, 2, \ldots, n-1$, and the magnetic quantum number $m = -\ell, \ldots, +\ell$. A fourth quantum number, the spin $s = \pm 1/2$, doubles the count.

Orbitals with $\ell = 0$ are called *s* orbitals. They are spherically symmetric — the wavefunction depends only on the distance from the nucleus, not on direction. The $1s$ orbital is the lowest-energy orbital and has the shape of a decaying exponential centred on the nucleus.

Orbitals with $\ell = 1$ are *p* orbitals. There are three of them ($m = -1, 0, +1$), each with a node at the nucleus and a dumb-bell shape pointing along the $x$, $y$, or $z$ axis (after a real-valued linear combination).

Orbitals with $\ell = 2$ are *d* orbitals: five of them, with more complex angular shapes (a four-lobed clover in the $xy$, $yz$, $xz$ planes, plus a doughnut-wrapped lobe along $z^2$).

Orbitals with $\ell = 3$ are *f* orbitals: seven of them, important for lanthanides and actinides, mostly irrelevant for the rest of the periodic table.

The Pauli exclusion principle says that no two electrons can occupy the same quantum state (same $n, \ell, m, s$). Accounting for the spin doubling, each orbital with given $(n, \ell, m)$ can hold *two* electrons of opposite spin. An $s$ subshell holds 2 electrons; a $p$ subshell holds 6 ($3 \times 2$); a $d$ subshell holds 10; an $f$ subshell holds 14.

## Electron configurations

The *electron configuration* of an atom lists which orbitals are occupied in its ground state. For an isolated atom, the orbitals are filled in approximately increasing order of energy according to the *aufbau* principle. The order, with some exceptions for transition metals, is:

$$
1s, \ 2s, \ 2p, \ 3s, \ 3p, \ 4s, \ 3d, \ 4p, \ 5s, \ 4d, \ 5p, \ 6s, \ 4f, \ 5d, \ 6p, \ \ldots
$$

Carbon, atomic number 6, has the ground-state configuration $1s^2 \, 2s^2 \, 2p^2$: two electrons in the $1s$ subshell, two in $2s$, two distributed among the three $2p$ orbitals. Silicon, atomic number 14, has $1s^2 \, 2s^2 \, 2p^6 \, 3s^2 \, 3p^2$, or equivalently $[\text{Ne}] \, 3s^2 \, 3p^2$, using the *core* shorthand. Copper, atomic number 29, has the slightly anomalous $[\text{Ar}] \, 3d^{10} \, 4s^1$ (the configuration with a full $3d$ shell is energetically favourable enough to break aufbau).

These configurations matter because the chemistry of an element is determined almost entirely by its outermost electrons.

## Valence versus core

Electrons in inner subshells are tightly bound to the nucleus and largely indifferent to the chemical environment. Their orbital shapes, energies, and densities are almost the same in an isolated atom as in a molecule or solid. Such electrons are called **core electrons**.

Electrons in the outermost subshell — the *valence* electrons — are the ones that form bonds, transfer charge between atoms, and respond to applied fields. They are weakly bound (their energies lie within a few eV of zero), their orbitals extend several ångströms from the nucleus, and their behaviour is sensitive to neighbouring atoms.

The distinction between core and valence is one of the most important practical simplifications in electronic-structure theory. For a heavy element like uranium, with 92 electrons, an explicit all-electron treatment requires representing the rapid spatial oscillations of inner-shell electrons close to the nucleus. The relevant length scales there are of order $a_0 / Z \approx 0.001$ Å, which would require an enormous plane-wave cutoff to resolve. The chemistry, however, is in the $5f$, $6d$, and $7s$ electrons — perhaps a dozen valence electrons out of 92. Nearly 90% of the electrons are along for the ride.

The standard practical move is to **freeze the core**: replace the inner-shell electrons by a fixed effective potential that the valence electrons feel, and solve only for the valence orbitals. This is the **pseudopotential** approximation. A pseudopotential is constructed for each element such that, beyond a chosen cutoff radius, the pseudo-orbital matches the true all-electron valence orbital; inside the cutoff, the pseudo-orbital is smoothed to avoid the rapid oscillations near the nucleus. Modern pseudopotentials — projector-augmented-wave (PAW) and ultrasoft (USPP) variants — are extremely accurate for total energy differences while drastically reducing computational cost.

There are subtleties. *Semi-core states* are inner-shell electrons that, in some chemical environments, are not quite frozen — the $3p$ electrons of early transition metals, for instance, can hybridise with valence states and must be treated explicitly. Modern pseudopotentials offer multiple versions per element (different valence configurations, with or without semi-core), and choosing the right one is part of a careful DFT setup. We will return to pseudopotentials in Chapter 6.

!!! note "When all-electron is necessary"
    For very heavy elements (lanthanides, actinides, post-transition metals where relativity matters), for core-level spectroscopies (X-ray absorption, XPS), and for benchmark calculations against experiment, an all-electron treatment is required. Codes such as FHI-aims and WIEN2k specialise in this. For routine total-energy and forces work on the elements up to Z $\sim 80$, pseudopotentials are the rule.

## The periodic table as a pattern of valence electrons

The genius of Mendeleev's periodic table is that it organises elements by valence-electron configuration. Each *row* corresponds to a value of the principal quantum number $n$ of the outermost shell. Each *block* corresponds to a value of $\ell$:

- the *s-block* (groups 1 and 2, plus helium) — outermost subshell is $ns$;
- the *p-block* (groups 13–18) — outermost subshell is $np$;
- the *d-block* (groups 3–12, the transition metals) — outermost is $(n{-}1)d$ with $ns$ also occupied;
- the *f-block* (lanthanides and actinides) — outermost is $(n{-}2)f$.

Elements in the same column (group) have the same valence-electron configuration up to principal quantum number, and therefore exhibit similar chemistry. Sodium and potassium both have one valence $s$ electron and behave similarly (both are alkali metals); the differences between them — sodium melts at 98 °C, potassium at 64 °C — come from the different size of their valence orbitals and the resulting differences in bond length and strength.

### Reading the periodic table for computational purposes

When you set up a DFT calculation, you need to know, for each element in your system:

1. **How many valence electrons** to include. For silicon, four ($3s^2 3p^2$). For copper, eleven ($3d^{10} 4s^1$). For iron, eight or sixteen depending on whether you include the $3p$ semi-core.

2. **What pseudopotentials are available** and which one is recommended for your accuracy and cost target. Pseudopotential libraries (PAW for VASP, SSSP for QE, GBRV) curate these choices.

3. **What basis-set parameters** to use. The plane-wave cutoff scales with the hardness of the pseudopotential, which depends on the element. Transition metals with $d$ electrons typically need higher cutoffs than main-group elements.

4. **Whether spin matters**. Elements with partially filled $d$ or $f$ subshells (Cr, Mn, Fe, Co, Ni, the lanthanides) are usually magnetic; you need spin-polarised DFT and an initial guess for the magnetic moment.

5. **Whether relativity matters**. For elements heavier than $Z \sim 50$, scalar-relativistic corrections (mass–velocity, Darwin terms) are essential. For $Z > 80$ or for spin–orbit-sensitive properties (Rashba splittings, topological phases, the band gap of lead halide perovskites), spin–orbit coupling must be included.

A computational materials scientist learns to read the periodic table for these signals at a glance. *Heavy transition metal, partially filled d* → magnetism plus modest relativity. *Lanthanide* → f electrons, strong correlation, all-electron preferable. *Lead or bismuth compound* → spin–orbit essential.

## Periodic trends

A few systematic trends across the periodic table affect bonding and properties:

**Atomic radius.** Atoms get larger as you go *down* a group (more shells) and smaller as you go *right* across a period (more nuclear charge pulling the same shell inward). Cs is the largest stable atom; H and He are among the smallest. Radii vary by a factor of $\sim 5$ across the table.

**Ionisation energy.** The energy to remove the outermost electron decreases down a group (less tightly bound) and increases across a period (more nuclear charge). The noble gases at the right edge have the highest ionisation energies; the alkali metals at the left edge have the lowest. Plot ionisation energy versus atomic number and you see periodic peaks at the noble gases.

**Electronegativity.** A composite measure of how strongly an atom attracts electrons in a bond. Fluorine is the most electronegative (Pauling scale: 4.0); caesium is the least (0.79). Electronegativity differences predict, roughly, whether a bond is ionic ($\Delta \chi$ large) or covalent ($\Delta \chi$ small).

**Polarisability.** The ease with which an atom's electron cloud distorts in response to an electric field. Polarisability increases down a group (larger atoms are more polarisable). Highly polarisable atoms participate in stronger van der Waals interactions and softer bonds.

These trends are taught in introductory chemistry but they remain useful day to day. Looking at a structure and thinking *which atom here has the loosest valence electrons?* is a habit worth cultivating.

## Special elements

A handful of elements deserve special mention because they appear repeatedly and have peculiarities.

**Hydrogen.** One electron, one proton. Sits awkwardly in the periodic table — it has the valence configuration of group 1 (one $s$ electron) but the electronegativity of group 17 (it forms covalent bonds). In simulation, hydrogen is light enough that quantum nuclear effects can matter at room temperature, and its lack of core electrons means standard pseudopotentials are sometimes replaced by an explicit treatment.

**Carbon.** Four valence electrons in a small atom. Capable of forming up to four bonds, and (uniquely among elements) of stable catenation in long chains. The basis of organic chemistry; a constant presence in materials simulation through polymers, biomolecules, graphene, carbides, and carbon-based catalysts.

**Transition metals.** Partially filled $d$ shells. Source of magnetism, of catalysis, of the rich chemistry of solids. Computationally awkward because $d$ electrons are localised enough that standard DFT often gives wrong magnetic and electronic answers. *Hubbard-U* corrections, hybrid functionals, or dynamical mean-field theory may be needed.

**Lanthanides and actinides.** Partially filled $f$ shells. Notoriously hard to simulate. Strong correlation, strong spin–orbit, multireference character. A research frontier; not for the new student.

**Noble gases.** Fully closed shells. Chemically inert to a first approximation; the rare gas crystals (solid argon, krypton) are bound only by van der Waals forces. The simplest possible test systems, and famously hard for standard DFT to describe accurately because of the vdW issue.

## What we throw away

The model we have arrived at — nuclei as classical point charges, electrons as quantum particles, core electrons replaced by pseudopotentials — discards a great deal:

- the internal structure of the nucleus (irrelevant for chemistry at the scales we care about);
- nuclear spins (relevant for NMR and hyperfine structure but absent from most simulations);
- quantum tunnelling of nuclei (relevant for light nuclei at low temperature; treated by path-integral methods when needed);
- the relativistic structure of inner-shell electrons (absorbed into pseudopotentials for heavy elements);
- subtle many-body correlations between electrons (the great approximation; addressed differently by every electronic-structure method);
- electron–photon coupling (relevant for optical properties; added back when needed).

Every approximation listed here can be relaxed at some computational cost. The art of materials simulation is in knowing which approximations matter for your question.

## A first computational view

ASE has a small chemical-data module that lets you query basic atomic information.

```python
from ase.data import chemical_symbols, atomic_numbers, atomic_masses, covalent_radii

for symbol in ("H", "C", "Si", "Cu", "Fe", "U"):
    Z = atomic_numbers[symbol]
    m = atomic_masses[Z]
    r = covalent_radii[Z]
    print(f"{symbol:>2s}  Z={Z:3d}  m={m:6.2f} u   r_cov={r:.2f} Å")
```

Running this prints a small table that already encodes a lot of what we will need: atomic number, mass, and a typical bonding radius. From here, the next section uses the bonding radius to discuss the qualitatively different ways atoms stick together.

!!! tip "Useful habit"
    When starting on a new system, write down on paper, for each element involved: the valence electron count, the typical oxidation state, whether magnetism is likely, and whether relativity matters. Five minutes spent doing this saves hours of confusion later.

With the atom established, we turn to bonding.
