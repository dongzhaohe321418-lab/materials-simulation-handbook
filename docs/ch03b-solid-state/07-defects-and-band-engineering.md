# 3b.7 — Defects, Doping, and Band Engineering

> *"Crystals are like people; it is the defects in them which tend to make them interesting."* — Sir Charles Frank

Every model so far has assumed a perfect, infinite, periodic crystal. Real materials are not like that. They contain *defects*: vacancies where an atom is missing, interstitials where an extra atom is squeezed in, substitutional impurities where one species replaces another, dislocations, grain boundaries, surfaces. Far from being a nuisance, defects are *the entire point* of most technological materials: dopants in silicon make transistors work; oxygen vacancies in oxides make memristors; lithium intercalates into graphite to store charge. This section is the briefest of introductions — enough to make sense of the band engineering language used throughout Tier 1 and the capstone projects.

The section covers four topics: point defect taxonomy, the shallow-impurity (hydrogenic) donor/acceptor picture, the effective mass approximation, and the language of band engineering — alloying, strain, quantum confinement.

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

The hydrogen-atom ground state binding energy and Bohr radius are $E_1 = -13.6$ eV, $a_0 = 0.529$ Å. With the substitutions,

$$\boxed{\; E_d = -13.6 \text{ eV}\cdot\frac{m^*/m_e}{\epsilon_r^2}, \qquad a_d = 0.529\text{ Å}\cdot\frac{\epsilon_r}{m^*/m_e}. \;} \tag{3b.7.3}$$

For silicon: $m^*/m_e \approx 0.26$ (averaged conduction-band effective mass), $\epsilon_r \approx 11.7$. So

$$E_d \approx 13.6 \cdot \frac{0.26}{(11.7)^2} \approx 0.026 \text{ eV} = 26 \text{ meV}, \tag{3b.7.4}$$

$$a_d \approx 0.529 \cdot \frac{11.7}{0.26} \approx 23.8 \text{ Å}. \tag{3b.7.5}$$

The donor electron sits ~26 meV *below* the conduction band edge and is delocalised over a sphere of radius ~24 Å — a region containing $(24/2.7)^3 \approx 700$ silicon atoms (with silicon's $\sim 5.4$ Å lattice). The state is "shallow": the binding energy is only about $k_B T$ at room temperature, so the donor electrons are fully ionised at 300 K and contribute to the conduction-band electron density.

Acceptors are the symmetric story for the *hole*: boron in silicon (group III) is missing one valence electron compared to silicon, leaving a positively charged hole bound to a negatively charged $\text{B}_\text{Si}^-$ centre. The same effective Bohr formula gives a similar shallow level ~45 meV above the valence band edge. The asymmetry in $E_d$ between donors (P, As, Sb) and acceptors (B, Al, Ga) comes from the difference in $m^*$ between conduction and valence bands and the band-structure anisotropy.

!!! note "When the hydrogenic picture breaks"
    The hydrogenic model assumes (i) the electron is weakly bound (so the dielectric constant of the bulk medium is the right screening parameter), and (ii) the effective mass approximation holds (so the electron's wavefunction is built from states near a single band edge). For *deep* defects — vacancies, transition-metal impurities, complex defects — neither assumption holds. The electron is localised on the defect site, sees the bare (unscreened) potential, and mixes states from multiple bands. Deep defects must be computed with full DFT, and they are responsible for most non-radiative recombination losses in solar cells.

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

## 3b.7.4 Band engineering: tuning bands by composition, strain, and confinement

Three principal levers for tuning the band structure of a material:

**Alloying / composition.** Mix two semiconductors $A$ and $B$ continuously, e.g. $\text{Si}_{1-x}\text{Ge}_x$ or $\text{Al}_{x}\text{Ga}_{1-x}\text{As}$. To first approximation, the band gap interpolates linearly with $x$ (Vegard's law) between the gaps of the endpoints, with a quadratic correction known as the *bowing parameter*. This is how we tune the bandgap of III-V semiconductors continuously from 0.3 eV (InSb) to 6 eV (AlN). Modern band engineering at the alloy level is high-throughput DFT screening across composition space — exactly what you will do in capstone project 2 for a halide perovskite alloy library.

**Strain.** Apply a uniaxial or biaxial strain. The strain breaks the cubic symmetry of the bands and splits degeneracies — for instance, the heavy-hole and light-hole bands of silicon, degenerate at $\Gamma$ in the unstrained crystal, split by $\sim 100$ meV under 1% uniaxial strain. This is the basis of *strained silicon* technology in modern CMOS: a thin Si layer grown on a slightly larger-lattice SiGe substrate is biaxially strained and has higher carrier mobility.

**Quantum confinement.** Restrict the electron's motion in one or more dimensions, by sandwiching a thin layer of low-bandgap material between two high-bandgap barriers (a *quantum well*), or by reducing all three dimensions (a *quantum dot*). The confinement raises the ground-state energy by roughly $\hbar^2\pi^2/(2m^* L^2)$ — the particle-in-a-box energy from Chapter 4, but with effective mass. By controlling $L$ (the well width or dot diameter) one can tune the optical absorption edge of the structure continuously across the visible spectrum. This is the principle of CdSe quantum-dot displays and InGaN blue LEDs.

In each case, the EMA equation (3b.7.9) governs the bound-state energies. The band-edge effective mass and the dielectric constant of the host material are the only material parameters you need. Both come from DFT band-structure calculations of the bulk host — which is the practical reason that band-structure calculations are *the* central computational task in semiconductor device design.

## 3b.7.5 Defects in MLIPs and graph neural networks

A word on how defects show up in the machine-learning side of the book.

**MLIPs.** A typical MLIP is trained on pristine bulk configurations plus a careful sprinkling of defect configurations. The transferability of the potential — whether it predicts the right energetics for an *unseen* defect — depends entirely on whether the local environments around the new defect are spanned by the training set. Modern active-learning workflows (§11) iteratively augment training data with high-uncertainty defect configurations until the MLIP converges. The success metric: defect formation energies within $\sim 50$ meV of DFT, across vacancy + interstitial + antisite + substitutional types.

**Graph neural networks.** GNNs naturally accommodate defects: a defect is just a graph in which one node has been removed (vacancy) or has different atomic features (substitutional). The same model architecture handles both pristine and defective crystals. This is one of the principal pragmatic advantages of GNN-based property prediction over traditional descriptor methods, which often require feature engineering specific to a defect type.

**Band engineering with ML.** The capstone project 1 will use a trained band-gap GNN to screen $10^4$ candidate dopants in a host semiconductor, predicting which substitutional impurities will produce shallow (useful) versus deep (recombination-active) levels. The shallow-vs-deep distinction is, microscopically, the question of whether the defect level is well-described by the EMA — and that distinction is exactly what a band-aware GNN can learn to predict.

!!! warning "What GNNs cannot do today"
    Predicting defect formation *energies* with chemical accuracy ($\le 0.1$ eV) remains an open problem. Most published GNNs achieve $\sim 0.2 - 0.5$ eV mean absolute error. The reason is that defects involve charge-transfer states whose energies depend sensitively on the dielectric environment — a long-ranged property that local message-passing struggles to capture. Active research in 2025–2026 focuses on equivariant GNNs with explicit charge prediction and screened-Coulomb message terms.

## Where this is used later

- **Tier 1.** §6.5 (defect formation energies in DFT supercells, with Freysoldt/Kumagai corrections), §6.7 (charge-transition levels and the configuration-coordinate diagram), §6.8 (point defects in 2D materials).
- **Tier 2.** §8.6 (defect-mediated diffusion and the Vineyard formula), §9.7 (training MLIPs for defective materials), §10.5 (GNNs for defect property prediction), §11.2 (active learning across the defect-configuration space).
- **Capstone Project 1.** Screening dopants for a target band-gap semiconductor — the entire workflow you will assemble rests on the EMA, hydrogenic donor model, and high-throughput band-edge effective-mass extraction.

This is the last content section. Proceed to §3b.8 for exercises that draw on all seven sections and to consolidate the picture before launching into quantum mechanics in Chapter 4 and DFT proper in Chapter 5.
