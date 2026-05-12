# 3b.7 — Defects, Doping, and Band Engineering

!!! tip "Why do defects and band engineering matter?"
    **The puzzle:** if the entire chapter so far has assumed a perfect, infinite, periodic crystal, why does anyone care about *imperfections*? Because every technologically relevant material is *deliberately imperfect*. Silicon would be a useless insulator at room temperature without phosphorus or boron doping. The semiconductors in your phone, computer, and TV screen all depend on impurity atoms at the parts-per-million level changing the carrier concentration by *orders of magnitude*.

    **The historical "aha":** in the 1940s-50s, Shockley and his colleagues at Bell Labs realised that adding tiny amounts of group-V impurities (P, As) to silicon — at concentrations of one impurity per million silicon atoms — could turn the material from an effective insulator into an n-type semiconductor with appreciable conductivity. The reason: each impurity atom donates one extra electron to the conduction band, weakly bound by the residual Coulomb field. The binding energy is small (~30 meV) thanks to the large dielectric constant of silicon ($\epsilon_r \approx 11.7$) and the small effective mass of conduction electrons ($m^* \approx 0.3 m_e$). At room temperature ($k_B T \approx 26$ meV), most donor electrons are ionised and free to conduct.

    **The picture to keep:** a dopant in a semiconductor is a *giant hydrogen atom*. The bound electron orbits the impurity in an enormous orbit (~20 Å, encompassing ~700 host atoms), screened by the dielectric, with a tiny binding energy. This is the *hydrogenic donor model* and it explains 90% of semiconductor doping in one equation: $E_d = -13.6\,\text{eV}\cdot (m^*/m_e)/\epsilon_r^2$.

    **What it buys us:** the entire field of band engineering — alloying, strain, quantum confinement, doping — rests on this single conceptual framework. The capstone projects in this book ask you to design materials by tuning these levers; the chapter gives you the tools to make those designs from first principles.

> *"Crystals are like people; it is the defects in them which tend to make them interesting."* — Sir Charles Frank

Every model so far has assumed a perfect, infinite, periodic crystal. Real materials are not like that. They contain *defects*: vacancies where an atom is missing, interstitials where an extra atom is squeezed in, substitutional impurities where one species replaces another, dislocations, grain boundaries, surfaces. Far from being a nuisance, defects are *the entire point* of most technological materials: dopants in silicon make transistors work; oxygen vacancies in oxides make memristors; lithium intercalates into graphite to store charge. This section is the briefest of introductions — enough to make sense of the band engineering language used throughout Tier 1 and the capstone projects.

The section covers four topics: point defect taxonomy, the shallow-impurity (hydrogenic) donor/acceptor picture, the effective mass approximation, and the language of band engineering — alloying, strain, quantum confinement.

!!! tip "Intuition: a giant hydrogen atom inside a crystal"
    Take a crystal of silicon and replace one silicon atom with a phosphorus atom. Silicon has four valence electrons; phosphorus has five. Four of phosphorus's electrons take part in the four covalent bonds with neighbouring silicon atoms — just like silicon's own electrons would. The fifth electron is left over, weakly bound to the phosphorus nucleus by the residual Coulomb attraction. But it does not see the bare Coulomb potential: the surrounding silicon lattice *screens* the charge by a factor of $\epsilon_r \approx 11.7$, dramatically weakening the binding. Moreover, the loosely bound electron does not move like a free electron — it moves like a conduction-band electron with effective mass $m^*\approx 0.3 m_e$, lighter than the bare electron mass.

    Both effects — small effective mass, large dielectric screening — push the binding energy *far below* the hydrogen-atom value of 13.6 eV. Specifically, the binding scales as $m^*/\epsilon_r^2 \approx 0.3/137 \approx 0.002$, giving a binding energy of order $\sim 30$ meV. That is the right order of magnitude for shallow donors in silicon, and it is comparable to $k_B T$ at room temperature, so these electrons are essentially fully ionised and contribute to the conduction band. The same logic in reverse explains shallow acceptors: a missing electron (hole) bound to an extra negative ion (boron in silicon).

    The picture extends to *every* form of band engineering. Confinement gives particle-in-a-box energies; strain shifts band edges; alloying interpolates band gaps. In each case the only material parameters that matter are the band-edge effective masses and the dielectric constant of the host — the same handful of numbers that come out of a single DFT calculation of the perfect crystal.

## 3b.7.1 Point defect taxonomy

Point defects are zero-dimensional: they perturb the lattice at a single site (or a few neighbouring sites). Four canonical types:

**Vacancy.** A missing atom. Denoted $V_X$ (sometimes $V_X^{q}$ if charged, with $q$ the integer charge state). A vacancy in silicon, $V_\text{Si}$, leaves four dangling bonds; the relaxed geometry typically distorts (Jahn–Teller) and the four levels split into deep and shallow states inside the gap.

**Interstitial.** An atom inserted at a non-lattice site. Denoted $X_i$. Silicon self-interstitials and oxygen interstitials in oxides are important examples. Interstitials tend to be high-energy defects that anneal out at moderate temperatures, but they dominate the kinetics during ion implantation.

**Substitutional impurity.** An impurity atom $Y$ replacing a host atom $X$. Denoted $Y_X$. The textbook example is phosphorus on a silicon site, $\text P_\text{Si}$, which contributes one extra electron compared to silicon and acts as a donor.

**Antisite defect.** In a compound $AB$, an atom of $A$ sitting on a $B$ site (or vice versa). Denoted $A_B$. In GaAs, the $\text{As}_\text{Ga}$ antisite is the famous EL2 deep level that pins the Fermi level near mid-gap in undoped GaAs and makes semi-insulating GaAs possible.

In each case, the *defect formation energy* — the free energy cost of creating one defect — is the central thermodynamic quantity. From it the equilibrium concentration of the defect is

$$c_\text{def} = N_\text{sites}\, e^{-E_\text{form}/k_B T}, \tag{3b.7.1}$$

where $N_\text{sites}$ is the number of available sites per unit volume. Typical $E_\text{form}$ values are 1–5 eV, giving negligible concentrations at room temperature ($e^{-E_\text{form}/k_B T} \sim 10^{-17}$ for $E_\text{form} = 1$ eV at 300 K) but appreciable concentrations at growth temperatures (1000 K).

Computing $E_\text{form}$ from DFT is the standard workflow you will execute in §6.5. The recipe: build a large supercell containing one defect; relax the geometry; compute the total energy; compute the total energy of the *pristine* supercell of the same size; take the difference, with chemical-potential corrections for the missing/added atoms. The "supercell finite-size error" — interactions between a defect and its periodic images — is the principal source of uncertainty, and is corrected by various schemes (Freysoldt, Makov–Payne) that you will learn in Ch 6.

## 3b.7.2 The shallow-impurity (hydrogenic) picture

Take silicon (group IV, four valence electrons per atom) and replace one silicon with phosphorus (group V, five valence electrons). The phosphorus atom donates four electrons to bonds with its silicon neighbours, leaving one extra electron loosely associated with the positively charged phosphorus nucleus. We can model this extra electron as a hydrogen-like atom: an electron of effective mass $m^*$ orbiting a positive point charge in a medium of dielectric constant $\epsilon_r$.

The Schrödinger equation for the electron, in atomic units rescaled by the medium, is

$$\left[-\frac{\hbar^2}{2m^*}\nabla^2 - \frac{e^2}{4\pi\epsilon_0\epsilon_r r}\right]\psi = E\psi. \tag{3b.7.2}$$

This is *exactly* the hydrogen atom, with two replacements:

- Bare electron mass $m_e \to m^*$.
- Vacuum permittivity $\epsilon_0 \to \epsilon_0\epsilon_r$.

!!! note "Why this step? — deriving the scaling from the hydrogen Rydberg"
    Recall the hydrogen-atom Rydberg energy:
    $$E_1 = -\frac{m_e e^4}{2(4\pi\epsilon_0)^2\hbar^2} = -13.6 \text{ eV},$$
    and the Bohr radius:
    $$a_0 = \frac{4\pi\epsilon_0\hbar^2}{m_e e^2} = 0.529 \text{ Å}.$$
    Now replace $m_e \to m^*$ and $\epsilon_0 \to \epsilon_0\epsilon_r$. The energy involves $m_e\cdot 1/\epsilon_0^2$, so it scales by $(m^*/m_e)/\epsilon_r^2$. The Bohr radius involves $\epsilon_0/m_e$, so it scales by $\epsilon_r/(m^*/m_e) = \epsilon_r m_e/m^*$. Putting it together:

The hydrogen-atom ground state binding energy and Bohr radius are $E_1 = -13.6$ eV, $a_0 = 0.529$ Å. With the substitutions,

$$\boxed{\; E_d = -13.6 \text{ eV}\cdot\frac{m^*/m_e}{\epsilon_r^2}, \qquad a_d = 0.529\text{ Å}\cdot\frac{\epsilon_r}{m^*/m_e}. \;} \tag{3b.7.3}$$

The donor level $E_d$ is measured *below* the conduction band edge $E_c$, so the binding energy in the formula above means $E_c - E_d$ for a donor, equivalently $E_d - E_v$ for an acceptor.

!!! example "Worked example: phosphorus donor in silicon"
    For silicon: $m^*/m_e \approx 0.26$ (averaged conduction-band effective mass, somewhat reduced from $0.3$ because of valley anisotropy averaging), $\epsilon_r \approx 11.7$. So

    $$E_d \approx 13.6 \cdot \frac{0.26}{(11.7)^2} \approx 0.026 \text{ eV} = 26 \text{ meV}, \tag{3b.7.4}$$

    $$a_d \approx 0.529 \cdot \frac{11.7}{0.26} \approx 23.8 \text{ Å}. \tag{3b.7.5}$$

    Using $m^*/m_e = 0.3$ instead (an alternative commonly tabulated value) gives $E_d \approx 30$ meV — the textbook figure for the phosphorus donor in silicon. The experimental binding energy of P in Si is $\approx 45$ meV — the discrepancy with the EMA value is the *central-cell correction*: the EMA assumes the donor electron sees only the long-range Coulomb tail of the impurity, but in reality the wavefunction has some weight near the impurity nucleus where the bare (unscreened) nuclear charge dominates. The correction is roughly a constant 10–20 meV shift downwards in binding energy for the various group-V donors in Si (P at 45 meV, As at 54 meV, Sb at 43 meV) — all near $k_B T$ at room temperature.

The donor electron sits ~26 meV *below* the conduction band edge and is delocalised over a sphere of radius ~24 Å — a region containing $(24/2.7)^3 \approx 700$ silicon atoms (with silicon's $\sim 5.4$ Å lattice). The state is "shallow": the binding energy is only about $k_B T$ at room temperature, so the donor electrons are fully ionised at 300 K and contribute to the conduction-band electron density.

Acceptors are the symmetric story for the *hole*: boron in silicon (group III) is missing one valence electron compared to silicon, leaving a positively charged hole bound to a negatively charged $\text{B}_\text{Si}^-$ centre. The same effective Bohr formula gives a similar shallow level ~45 meV above the valence band edge. The asymmetry in $E_d$ between donors (P, As, Sb) and acceptors (B, Al, Ga) comes from the difference in $m^*$ between conduction and valence bands and the band-structure anisotropy.

!!! note "When the hydrogenic picture breaks"
    The hydrogenic model assumes (i) the electron is weakly bound (so the dielectric constant of the bulk medium is the right screening parameter), and (ii) the effective mass approximation holds (so the electron's wavefunction is built from states near a single band edge). For *deep* defects — vacancies, transition-metal impurities, complex defects — neither assumption holds. The electron is localised on the defect site, sees the bare (unscreened) potential, and mixes states from multiple bands. Deep defects must be computed with full DFT, and they are responsible for most non-radiative recombination losses in solar cells.

### Acceptor levels: the symmetric story

A boron atom (group III, three valence electrons) substituted on a silicon site has *one too few* electrons to complete the four-coordinated covalent bonding. In a picture where silicon has fully bonded $sp^3$ hybrids, the missing electron leaves a *hole* — a vacant valence-band state — weakly bound to the negatively charged $\text B_\text{Si}^-$ centre (negative because the boron core has fewer protons than silicon). The hole is delocalised across many silicon atoms by the same hydrogenic mechanism, just with the *valence-band* effective mass $m_h^*$ instead of $m_e^*$.

For silicon, $m_h^*/m_e \approx 0.49$ (the heavy hole mass; light hole is $\approx 0.16$, and the actual binding depends on a Luttinger-Kohn calculation of the warped valence band, not a single $m_h^*$). Using $0.49$:

$$E_a \approx 13.6\cdot \frac{0.49}{(11.7)^2} \approx 0.049 \text{ eV} = 49 \text{ meV}.$$

The experimental boron acceptor level in silicon is $\approx 45$ meV above the valence band edge, in surprisingly close agreement. As with donors, "above the valence band" means: the hole sits in a bound state $\approx 45$ meV above $E_v$, energetically slightly below the bare valence band edge by the same amount (energies of holes are measured downward from $E_v$). At room temperature ($k_B T \approx 26$ meV) the bound hole is largely ionised and contributes to the hole density in the valence band.

The asymmetry $E_d^\text{P} = 45$ meV vs $E_a^\text{B} = 45$ meV is coincidental for these particular dopants. Across other group V donors (As, Sb) and group III acceptors (Al, Ga, In) the values cluster around 30–70 meV — all in the shallow regime, all available as room-temperature carriers.

## 3b.7.3 The effective mass approximation

The hydrogenic argument above was a special case of the **effective mass approximation** (EMA), which is the workhorse method for treating *any* slowly varying perturbation of a perfect crystal. The setup: take a perfect crystal with bands $E_n(\mathbf k)$ and Bloch states $\psi_{n\mathbf k}$. Add a slowly varying potential $U(\mathbf r)$ — slow meaning that $U$ does not vary appreciably on the scale of a unit cell. We ask: what are the new eigenstates and energies?

In the EMA one expands around a band extremum at $\mathbf k_0$, say the conduction-band minimum. To quadratic order in $\mathbf k - \mathbf k_0$,

$$E_n(\mathbf k) \approx E_n(\mathbf k_0) + \sum_{\alpha\beta} \frac{\hbar^2}{2 m^*_{\alpha\beta}}(k_\alpha - k_{0\alpha})(k_\beta - k_{0\beta}), \tag{3b.7.6}$$

with the **effective mass tensor**

$$\frac{1}{m^*_{\alpha\beta}} = \frac{1}{\hbar^2}\frac{\partial^2 E_n}{\partial k_\alpha\, \partial k_\beta}\bigg|_{\mathbf k_0}. \tag{3b.7.7}$$

In isotropic bands (silicon's conduction band is approximately isotropic if averaged over the six equivalent valleys near $X$), $m^*_{\alpha\beta} = m^*\delta_{\alpha\beta}$ and equation (3b.7.6) reduces to a single scalar effective mass.

The wavefunction of the slowly perturbed crystal state is approximately

$$\Psi(\mathbf r) \approx F(\mathbf r)\, \psi_{n\mathbf k_0}(\mathbf r), \tag{3b.7.8}$$

where $F(\mathbf r)$ is an *envelope function* that varies slowly on the scale of a unit cell, and $\psi_{n\mathbf k_0}$ is the Bloch state at the band edge. The envelope satisfies an effective Schrödinger equation:

$$\left[-\frac{\hbar^2}{2 m^*}\nabla^2 + U(\mathbf r)\right] F(\mathbf r) = (E - E_n(\mathbf k_0))\, F(\mathbf r). \tag{3b.7.9}$$

This is a clean piece of physics: the *envelope* of the wavefunction obeys a Schrödinger equation with the bare electron mass replaced by the effective mass, and with the potential being the *extra* perturbation $U$ (not the periodic crystal potential, which has been absorbed into the band structure). The Bloch state at the band edge serves as a "wallpaper" pattern that the envelope modulates.

The EMA is the language in which essentially all device physics is written: the semiconductor equations in a transistor, the band-bending profiles in a heterojunction, the bound states of an exciton in a quantum well — all are envelope-function calculations with appropriate $m^*$ and $\epsilon_r$ as inputs from the underlying ab initio band structure.

!!! note "Why this step? — the EMA is exact in a specific limit"
    The EMA can be made rigorous: one expands the wavefunction in Bloch states at the band edge $\mathbf k_0$, $\psi(\mathbf r) = \int d^3 q\, F(\mathbf q)\, \psi_{n,\mathbf k_0 + \mathbf q}(\mathbf r)$, and recognises that for slowly varying $U(\mathbf r)$ the envelope $F(\mathbf q)$ is sharply peaked near $\mathbf q = 0$. In this limit the dispersion can be replaced by its quadratic approximation (3b.7.6), and the resulting Schrödinger-like equation for $F$ is exactly (3b.7.9). The approximation is controlled by the smallness of $\nabla U/k_F$ — i.e. how rapidly the perturbation varies compared to the Fermi wavelength. For donor wavefunctions extending over $\sim 20$ Å in silicon, $U$ varies on the scale of $20$ Å while the Fermi wavelength is $\sim 1$ Å, so the small parameter is $\sim 1/20$ and the EMA is good to a few per cent.

    The EMA fails dramatically when the binding potential is *sharp* on the scale of the unit cell (deep defects), or when several bands lie close in energy and mix (heavy-hole / light-hole coupling in p-type semiconductors). The remedy in the second case is multi-band EMA, also called the Kane or Luttinger-Kohn model; in the first case, one must abandon EMA entirely and do a full DFT supercell calculation.

## 3b.7.4 Band engineering: tuning bands by composition, strain, and confinement

Three principal levers for tuning the band structure of a material:

**Alloying / composition.** Mix two semiconductors $A$ and $B$ continuously, e.g. $\text{Si}_{1-x}\text{Ge}_x$ or $\text{Al}_{x}\text{Ga}_{1-x}\text{As}$. To first approximation, the band gap interpolates linearly with $x$ between the gaps of the endpoints, with a quadratic correction known as the *bowing parameter*:

$$E_g(x) = (1-x)\, E_g(A) + x\, E_g(B) - b\, x(1-x), \tag{3b.7.\text{alloy}}$$

where $b$ is the bowing parameter (typically 0–1 eV; for $\text{In}_x\text{Ga}_{1-x}\text{As}$, $b \approx 0.4$ eV). For $\text{Al}_x\text{Ga}_{1-x}\text{As}$ in the direct-gap regime ($x < 0.45$), $E_g(0) = 1.42$ eV (GaAs), $E_g(0.45) \approx 2.0$ eV, with a small bowing. This is how we tune the bandgap of III-V semiconductors continuously from 0.3 eV (InSb) to 6 eV (AlN). Modern band engineering at the alloy level is high-throughput DFT screening across composition space — exactly what you will do in capstone project 2 for a halide perovskite alloy library.

**Strain.** Apply a uniaxial or biaxial strain. The strain breaks the cubic symmetry of the bands and splits degeneracies — for instance, the heavy-hole and light-hole bands of silicon, degenerate at $\Gamma$ in the unstrained crystal, split by $\sim 100$ meV under 1% uniaxial strain. This is the basis of *strained silicon* technology in modern CMOS: a thin Si layer grown on a slightly larger-lattice SiGe substrate is biaxially strained and has higher carrier mobility.

To leading order, the *hydrostatic* component of strain (uniform expansion or compression) shifts the band edges through the **deformation potential**:

$$\boxed{\; \Delta E_g = a_v\, (\varepsilon_{xx} + \varepsilon_{yy} + \varepsilon_{zz}) = a_v\, \text{Tr}(\varepsilon). \;} \tag{3b.7.\text{strain}}$$

Here $\varepsilon_{\alpha\beta}$ is the strain tensor (with $\varepsilon_{\alpha\alpha} = (L_\alpha - L_\alpha^0)/L_\alpha^0$ being the fractional change in length along axis $\alpha$), and $a_v$ is the *volume deformation potential* — a material-specific quantity, usually of order 1–10 eV. The deformation potential captures the physically intuitive fact that compressing a crystal increases overlaps between neighbouring orbitals, broadening bands and typically *reducing* the gap.

!!! note "Why this step? — derivation of the deformation potential to lowest order"
    Treat the strain $\varepsilon_{\alpha\beta}$ as a small perturbation that modifies the crystal Hamiltonian. To linear order in $\varepsilon$, the change in the band energy at any fixed $\mathbf k$ is
    $$\Delta E_n(\mathbf k) = \sum_{\alpha\beta} D_{n,\alpha\beta}(\mathbf k)\, \varepsilon_{\alpha\beta},$$
    where $D_{n,\alpha\beta}(\mathbf k) := \partial E_n(\mathbf k)/\partial\varepsilon_{\alpha\beta}|_{\varepsilon=0}$ is the deformation potential tensor for band $n$. For a cubic crystal at a high-symmetry point like $\Gamma$, the only allowed component (by cubic symmetry) is the trace: $D_{n,\alpha\beta} = a_n\delta_{\alpha\beta}$, so $\Delta E_n = a_n\, \text{Tr}(\varepsilon)$. The volume deformation potential of the gap is $a_v = a_c - a_v^{(v)}$, where $a_c$ is for the conduction band edge and $a_v^{(v)}$ for the valence band edge. For silicon, $a_v \approx -1.5$ eV: the gap *shrinks* by 15 meV per 1% hydrostatic compression.

    Beyond hydrostatic strain, shear components $\varepsilon_{xy}$ etc. *split* degeneracies (e.g. heavy-hole vs light-hole) without shifting the average, and require additional deformation-potential constants $b$ and $d$. The full theory (the Bir–Pikus Hamiltonian) is in any solid-state textbook.

**Quantum confinement.** Restrict the electron's motion in one or more dimensions, by sandwiching a thin layer of low-bandgap material between two high-bandgap barriers (a *quantum well*), or by reducing all three dimensions (a *quantum dot*). The confinement raises the ground-state energy by roughly $\hbar^2\pi^2/(2m^* L^2)$ — the particle-in-a-box energy from Chapter 4, but with effective mass. By controlling $L$ (the well width or dot diameter) one can tune the optical absorption edge of the structure continuously across the visible spectrum. This is the principle of CdSe quantum-dot displays and InGaN blue LEDs.

!!! note "Why this step? — deriving the confinement shift in detail"
    For a quantum well of width $L$ along $z$ with infinite barriers, the envelope function satisfies
    $$-\frac{\hbar^2}{2m^*}\frac{d^2 F}{dz^2} = (E - E_c)\, F, \qquad F(0) = F(L) = 0.$$
    The standard particle-in-a-box solution gives
    $$F_n(z) = \sqrt{2/L}\, \sin(n\pi z/L), \qquad E_n - E_c = \frac{\hbar^2 \pi^2 n^2}{2 m^* L^2}.$$
    The ground state ($n=1$) is shifted by $\Delta E = \hbar^2\pi^2/(2 m^* L^2)$. For the valence-band hole, the same calculation with $m^*_h$ gives a downward shift of the hole ground state by $\hbar^2\pi^2/(2 m_h^* L^2)$. The optical band gap of the quantum well is therefore
    $$E_g^\text{QW}(L) = E_g^\text{bulk} + \frac{\hbar^2\pi^2}{2 L^2}\left(\frac{1}{m_e^*} + \frac{1}{m_h^*}\right) = E_g^\text{bulk} + \frac{\hbar^2\pi^2}{2\mu^* L^2},$$
    with $1/\mu^* = 1/m_e^* + 1/m_h^*$. The confinement raises the gap *quadratically* with $1/L$. For a 5 nm CdSe quantum dot with $\mu^* \approx 0.1\, m_e$:
    $$\Delta E_g = \frac{(1.055\times 10^{-34})^2 \pi^2}{2\cdot 0.1\cdot 9.11\times 10^{-31}\cdot (5\times 10^{-9})^2} \approx 2.4\times 10^{-20}\text{ J} \approx 0.15\text{ eV}.$$
    Adding to the bulk CdSe gap of $\approx 1.74$ eV gives $\approx 1.89$ eV — corresponding to orange emission. Shrink the dot to 2 nm: $\Delta E_g \approx 0.95$ eV, total gap $\approx 2.7$ eV, blue. This continuous size-tunable colour is the central design lever for quantum-dot displays and solar cells.

In each case, the EMA equation (3b.7.9) governs the bound-state energies. The band-edge effective mass and the dielectric constant of the host material are the only material parameters you need. Both come from DFT band-structure calculations of the bulk host — which is the practical reason that band-structure calculations are *the* central computational task in semiconductor device design.

### Worked example: thin quantum well from GaAs/AlGaAs

Consider a 5 nm GaAs quantum well sandwiched between two AlGaAs barriers (a typical structure for near-infrared lasers). Parameters: $m^*_e = 0.067 m_e$ for the GaAs conduction band, $m^*_h \approx 0.5 m_e$ for the heavy-hole valence band (weighted average), bulk gap of GaAs $E_g^\text{bulk} = 1.42$ eV.

**Reduced mass:**
$$\frac{1}{\mu^*} = \frac{1}{0.067} + \frac{1}{0.5} \approx 16.93 \quad\Rightarrow\quad \mu^* \approx 0.059\, m_e.$$

**Confinement shift:**
$$\Delta E_g = \frac{\hbar^2\pi^2}{2\mu^* L^2} = \frac{(1.055\times 10^{-34})^2 \pi^2}{2 \cdot 0.059 \cdot 9.11\times 10^{-31} \cdot (5\times 10^{-9})^2} \approx 4.07\times 10^{-20}\text{ J} \approx 0.25\text{ eV}.$$

**Effective optical gap:**
$$E_g^\text{QW} = 1.42 + 0.25 = 1.67 \text{ eV}, \quad \lambda \approx 742 \text{ nm} \text{ (near-infrared)}.$$

Shrink the well to 2 nm: $\Delta E_g \propto 1/L^2$ gives a factor of $(5/2)^2 = 6.25$ enhancement, so $\Delta E_g \approx 1.56$ eV. Total $E_g^\text{QW} \approx 2.98$ eV, $\lambda \approx 416$ nm — deep violet. This is the basic design principle of size-tunable laser diodes.

In a *quantum dot* (confined in all three dimensions), the shift is roughly 3x larger because confinement enters in each direction:

$$\Delta E_g^\text{dot} = \frac{\hbar^2\pi^2}{2\mu^* L^2}\cdot 3 \quad\text{(approximately, for cubic confinement)}.$$

For a 5 nm GaAs dot, $\Delta E_g^\text{dot} \approx 0.75$ eV, and the optical gap moves into the visible regime even from a starting material with a near-IR bulk gap. This is exactly the principle of CdSe/CdS quantum-dot LEDs used in modern televisions: colour-tunable emission from a single chemical compound by changing the dot size during synthesis.

## 3b.7.5 Defects in MLIPs and graph neural networks

A word on how defects show up in the machine-learning side of the book.

**MLIPs.** A typical MLIP is trained on pristine bulk configurations plus a careful sprinkling of defect configurations. The transferability of the potential — whether it predicts the right energetics for an *unseen* defect — depends entirely on whether the local environments around the new defect are spanned by the training set. Modern active-learning workflows (§11) iteratively augment training data with high-uncertainty defect configurations until the MLIP converges. The success metric: defect formation energies within $\sim 50$ meV of DFT, across vacancy + interstitial + antisite + substitutional types.

**Graph neural networks.** GNNs naturally accommodate defects: a defect is just a graph in which one node has been removed (vacancy) or has different atomic features (substitutional). The same model architecture handles both pristine and defective crystals. This is one of the principal pragmatic advantages of GNN-based property prediction over traditional descriptor methods, which often require feature engineering specific to a defect type.

**Band engineering with ML.** The capstone project 1 will use a trained band-gap GNN to screen $10^4$ candidate dopants in a host semiconductor, predicting which substitutional impurities will produce shallow (useful) versus deep (recombination-active) levels. The shallow-vs-deep distinction is, microscopically, the question of whether the defect level is well-described by the EMA — and that distinction is exactly what a band-aware GNN can learn to predict.

!!! warning "What GNNs cannot do today"
    Predicting defect formation *energies* with chemical accuracy ($\le 0.1$ eV) remains an open problem. Most published GNNs achieve $\sim 0.2 - 0.5$ eV mean absolute error. The reason is that defects involve charge-transfer states whose energies depend sensitively on the dielectric environment — a long-ranged property that local message-passing struggles to capture. Active research in 2025–2026 focuses on equivariant GNNs with explicit charge prediction and screened-Coulomb message terms.

## 3b.7.6 Section summary

!!! abstract "Key ideas"
    1. **Defect taxonomy.** Vacancies, interstitials, substitutional impurities, antisites. Concentrations $c \propto N_\text{sites}\, e^{-E_\text{form}/k_BT}$. Computing $E_\text{form}$ from DFT supercells is a Ch 6 workflow.
    2. **Hydrogenic donors/acceptors.** Shallow impurity bound state behaves like a hydrogen atom with $m_e\to m^*$ and $\epsilon_0\to\epsilon_0\epsilon_r$. Binding energy $E_d = 13.6\,\text{eV}\cdot(m^*/m_e)/\epsilon_r^2 \sim 30$ meV in Si. Orbital radius $\sim 20$ Å.
    3. **Effective mass approximation.** Slowly-varying perturbations of a perfect crystal can be treated by replacing $\hat p^2/2m$ with $\hat p^2/(2m^*)$ for an envelope function. The Bloch state at the band edge is the wallpaper.
    4. **Band engineering.** Three levers: composition (alloying with bowing), strain (deformation potential $\Delta E_g = a_v\text{Tr}(\varepsilon)$), and quantum confinement (particle-in-a-box energies $\sim \hbar^2\pi^2/(2m^*L^2)$).
    5. **ML on defects.** GNNs accommodate defects naturally; current state-of-the-art has $\sim 0.2$–$0.5$ eV MAE on formation energies; charge-transfer levels remain open problem.

What to remember three months from now: *"A dopant in a semiconductor is a giant hydrogen atom; its binding energy is $13.6\,\text{eV}\cdot(m^*/m_e)/\epsilon_r^2$, typically $\sim 30$ meV, so dopants are mostly ionised at room temperature. Band engineering — alloying, strain, confinement — works in the effective-mass framework with these few material parameters."*

## Where this is used later

- **Tier 1.** §6.5 (defect formation energies in DFT supercells, with Freysoldt/Kumagai corrections), §6.7 (charge-transition levels and the configuration-coordinate diagram), §6.8 (point defects in 2D materials).
- **Tier 2.** §8.6 (defect-mediated diffusion and the Vineyard formula), §9.7 (training MLIPs for defective materials), §10.5 (GNNs for defect property prediction), §11.2 (active learning across the defect-configuration space).
- **Capstone Project 1.** Screening dopants for a target band-gap semiconductor — the entire workflow you will assemble rests on the EMA, hydrogenic donor model, and high-throughput band-edge effective-mass extraction.

This is the last content section. Proceed to §3b.8 for exercises that draw on all seven sections and to consolidate the picture before launching into quantum mechanics in Chapter 4 and DFT proper in Chapter 5.
