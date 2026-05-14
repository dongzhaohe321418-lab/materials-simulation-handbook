# 5.4 Exchange–Correlation Functionals

!!! note "Why does this chapter exist?"
    Kohn–Sham theory (§5.3) is exact *in principle* — if we knew the exchange–correlation functional $E_{xc}[n]$ in closed form, we could compute the energy of any molecule or crystal to arbitrary accuracy. We do not know it in closed form. What we have is *approximations* — and over sixty years, hundreds of them. Choosing the right one is the central practical decision in any DFT calculation.
    
    A useful analogy. Imagine you are a chef, and the perfect recipe for a sauce is locked in a vault. You know the sauce exists, you know it would taste exactly right, but you cannot open the vault. You have to *guess* the recipe based on principles (it should contain salt, fat, acid), constraints (no ingredient may be negative, the total weight must be conserved), and benchmark dishes (the right sauce should taste like *this* on chicken, like *that* on fish). Different chefs make different guesses; we call them *functionals*. LDA is the simplest plausible guess; PBE adds a sensible gradient correction; SCAN adds more constraints; HSE06 mixes in some exact (but expensive) ingredient. Each climbs one rung of John Perdew's metaphorical *Jacob's ladder* toward the unreachable exact functional in the heavens.
    
    Why does LDA work as well as it does, given how crude it is? The answer is one of the great minor miracles of modern physics. LDA gets exchange right in regions where the density is uniform (which is *no* region of any real atom!). The non-uniform errors should be huge — and they are, individually. But they *cancel* between exchange and correlation, in ways that the so-called "adiabatic-connection" picture (§5.3) makes transparent. So LDA gives binding energies accurate to ~30 kcal/mol (off by enough to be useless for chemistry, but in the right ballpark for atomic energies) — better than its crude derivation deserves. This cancellation is the deep reason DFT works.

The Kohn–Sham construction (§5.3) is exact in principle. It becomes an approximation only because we do not know the exchange–correlation energy functional $E_{xc}[n]$ in closed form. Every DFT calculation in the world today involves a choice of approximate $E_{xc}$: a *functional*. The choice matters. The same molecule can have its bond length reproduced within 1 pm by one functional and missed by 10 pm by another; a band gap can come out qualitatively right or qualitatively wrong; a magnetic ground state can flip.

The functional zoo is large — thousands have been proposed; perhaps fifty are in common use. To navigate it, John Perdew suggested a marvellous metaphor: **Jacob's ladder**. Each rung adds an ingredient and, on average, climbs toward chemical accuracy ($\sim 1$ kcal/mol $\approx 0.04$ eV). Each rung also costs more. We climb the ladder in turn.

!!! abstract "Key idea (Chapter 5.4)"
    The exchange–correlation functional $E_{xc}[n]$ is the *one* thing standing between Kohn–Sham DFT and exact results. We do not know its closed form, but we can build approximations: LDA (uniform-gas limit), GGA (gradient corrections), meta-GGA (kinetic-energy density), hybrid (mix in exact exchange), double-hybrid (add MP2-like correlation). Each rung up *Jacob's ladder* adds ingredients and accuracy at increasing cost. The choice of functional governs the accuracy of every downstream prediction.

## 5.4.1 Jacob's ladder

From bottom to top:

1. **LDA**: uses $n(\mathbf r)$ only.
2. **GGA**: uses $n$ and $|\nabla n|$.
3. **meta-GGA**: adds the kinetic energy density $\tau(\mathbf r)$ or $\nabla^{2}n$.
4. **Hybrid**: mixes in a fraction of exact (Hartree–Fock) exchange.
5. **Double hybrid** / RPA / wavefunction methods: include unoccupied orbitals.

The first three rungs are purely *semi-local* — the value of $\epsilon_{xc}$ at $\mathbf r$ depends only on quantities at $\mathbf r$ (or its immediate gradients). The fourth introduces non-locality through exact exchange and is roughly an order of magnitude more expensive. The fifth introduces virtual orbitals and another order of magnitude.

!!! note "Why this step?"
    Each rung adds a new physical ingredient — and a new opportunity to break exact constraints if one is not careful. The ladder is *not* monotonically more accurate: a poorly-constructed meta-GGA can be worse than a well-constructed GGA. The metaphor is aspirational, not guaranteed. What is true is that *more ingredients* allow *more constraints* to be satisfied simultaneously; whether they actually are depends on the functional designer's craft. PBE is a clean GGA that has aged well precisely because its constants were determined by exact constraints, not fits to data; SCAN is a meta-GGA that satisfies seventeen exact constraints simultaneously. Both stand the test of time better than their many empirically-fit competitors.

!!! example "Cost ratio summary"
    For a benchmark of organic molecules on the same hardware:
    
    | Rung | Functional | Time / GGA |
    |---|---|---|
    | 1 LDA | PW92 | 0.9 |
    | 2 GGA | PBE | 1.0 (reference) |
    | 3 meta-GGA | SCAN | 1.5 |
    | 4 hybrid (screened) | HSE06 | 10–30 |
    | 4 hybrid (global) | PBE0 | 15–50 |
    | 5 double-hybrid | XYG3 | 100–500 |
    | wavefunction | CCSD(T) | $10^{4}$–$10^{6}$ |
    
    The cost increase with rung is dominated by the cost of evaluating exact exchange (non-local, scales as $\mathcal O(N_\mathrm{occ}^{2})$ rather than $\mathcal O(N)$). Modern range-separation tricks reduce the prefactor but not the scaling.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. What single quantity stands between Kohn–Sham DFT and exact results, and why is it only ever known approximately?
    2. Name the five rungs of Jacob's ladder in order, and state the new ingredient each one adds.
    3. Why is climbing the ladder not guaranteed to improve accuracy, and why are the first three rungs much cheaper than the fourth?

    If any of these is shaky, re-read the preceding section before continuing.

## 5.4.2 LDA: the local density approximation

The simplest approximation: pretend that, locally, the electron gas is uniform. Define an exchange–correlation energy density per particle, $\epsilon_{xc}^\mathrm{unif}(n)$, for a uniform electron gas of density $n$. Then

$$
\boxed{\;\;E_{xc}^\mathrm{LDA}[n] \;=\; \int n(\mathbf r)\,\epsilon_{xc}^\mathrm{unif}\!\big(n(\mathbf r)\big)\,\mathrm d\mathbf r.\;\;}
\tag{5.34}
$$

Write $\epsilon_{xc}^\mathrm{unif} = \epsilon_{x}^\mathrm{unif} + \epsilon_{c}^\mathrm{unif}$: the exchange part can be computed in closed form; the correlation part is known very accurately from quantum Monte Carlo (Ceperley and Alder, 1980) and fitted to convenient analytic forms (VWN, Perdew–Zunger, Perdew–Wang).

### Derivation of LDA exchange

For a uniform electron gas of density $n$, the exchange energy per unit volume is

$$
\epsilon_x^\mathrm{vol}(n) = -\frac{1}{2}\int\frac{\rho_x(\mathbf r,\mathbf r')\,n}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r',
$$

where $\rho_x$ is the exchange hole. For a single Slater determinant of plane waves, the exchange hole is computable analytically. We take a more direct route via the Fock energy of the Hartree–Fock ground state of the uniform gas.

The Hartree–Fock exchange energy of $N$ plane-wave electrons (two per $\mathbf k$ up to $k_F$) is

$$
E_x = -\frac{1}{2}\sum_{\mathbf k,\mathbf k'}^\mathrm{occ}\int\frac{e^{-i(\mathbf k-\mathbf k')\cdot\mathbf r}\,e^{i(\mathbf k-\mathbf k')\cdot\mathbf r'}}{L^{6}\,|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r'.
$$

Converting sums to integrals and doing the spatial integral (Fourier transform of $1/|\mathbf r|$ is $4\pi/q^{2}$), one obtains for the exchange energy per unit volume

$$
\frac{E_x}{L^{3}} = -\frac{1}{L^{3}}\cdot 2 \cdot \left(\frac{L^{3}}{(2\pi)^{3}}\right)^{2}\int_{k<k_F}\int_{k'<k_F}\frac{4\pi}{|\mathbf k-\mathbf k'|^{2}}\,\mathrm d\mathbf k\,\mathrm d\mathbf k'.
$$

The double integral evaluates (the calculation is in many textbooks; the trick is the substitution $\mathbf q = \mathbf k - \mathbf k'$):

$$
\int_{k<k_F}\int_{k'<k_F}\frac{1}{|\mathbf k - \mathbf k'|^{2}}\,\mathrm d\mathbf k\,\mathrm d\mathbf k' = 2\pi k_F^{4}.
$$

Substituting and dividing out $L^{3}$ to obtain the energy per particle (using $n = k_F^{3}/(3\pi^{2})$),

$$
\epsilon_x^\mathrm{unif}(n) = -\frac{3}{4\pi}k_F = -\frac{3}{4}\Big(\frac{3}{\pi}\Big)^{1/3}\,n^{1/3}.
\tag{5.35}
$$

Equation (5.35) is Dirac's 1930 result. The LDA exchange functional is therefore

$$
\boxed{\;\;E_{x}^\mathrm{LDA}[n] = -\frac{3}{4}\Big(\frac{3}{\pi}\Big)^{1/3}\int n(\mathbf r)^{4/3}\,\mathrm d\mathbf r.\;\;}
\tag{5.36}
$$

The corresponding potential is $v_x^\mathrm{LDA}(\mathbf r) = -(3/\pi)^{1/3}n(\mathbf r)^{1/3}$. The LDA exchange–correlation potential $v_{xc}^\mathrm{LDA}$ is the sum of (5.36)'s functional derivative and the (numerical) $v_c^\mathrm{LDA}$ from the parametrised correlation energy.

!!! note "Why this step?"
    A useful reparametrisation: define the Wigner–Seitz radius $r_s$ by $\tfrac{4}{3}\pi r_s^{3}\,n = 1$, so that $r_s = (3/(4\pi n))^{1/3}$ is the radius of the sphere containing one electron on average. In terms of $r_s$, the LDA exchange energy per particle is
    $$
    \epsilon_x^\mathrm{unif} = -\frac{3}{4}\Big(\frac{9}{4\pi^{2}}\Big)^{1/3}\frac{1}{r_s} \approx -\frac{0.4582}{r_s}\;\text{Ha},
    $$
    which makes manifest that exchange scales as $1/r_s$. For metallic densities, $r_s\sim 2$–$5$, giving $\epsilon_x\sim -0.1$ to $-0.2\;\text{Ha}$ per electron, i.e. a few eV. This is the right ballpark for the exchange contribution to atomic and molecular binding.

!!! example "Worked example: $H_2$ binding energy at the LDA level"
    The hydrogen molecule has an experimental binding energy of $D_e = 4.748\;\text{eV}$ at $R_e = 0.741\;\text{Å}$. An LDA calculation gives $D_e^\mathrm{LDA}\approx 4.91\;\text{eV}$ and $R_e^\mathrm{LDA}\approx 0.766\;\text{Å}$: overbinding by $\sim 3\%$, bond too short by $\sim 3\%$. This is the prototypical "LDA overbinding" pattern — the bond is too strong, the lattice constant too short, the cohesive energy too large. The error is roughly half attributable to exchange (LDA-X is too soft in the bond region) and half to correlation (LDA-C overestimates the magnitude of correlation in inhomogeneous systems).

### Where LDA works, where it fails

**Strengths.**

- *Free-electron-like systems*: bulk metals (Na, Al, Cu) — bond lengths and bulk moduli within a few per cent. The uniform-gas reference is a reasonable starting point when the density is genuinely slowly varying.
- *Total energies*: the cohesive energy of a metal comes out plausibly, often within 0.5 eV/atom.
- *Geometries*: equilibrium bond lengths in simple solids are usually within 1–2% of experiment.

**Weaknesses.**

- *Overbinding*: LDA systematically overestimates binding energies, often by tens of per cent. The H$_2$ binding energy is 4.75 eV experimentally; LDA gives roughly 4.9 eV with a too-short bond, and the error grows for larger molecules.
- *Lattice constants*: LDA gives lattice constants about 1–3% *too small* (the famous "LDA overbinding").
- *Band gaps*: LDA underestimates band gaps by 30–100%. (This has two distinct causes — the derivative discontinuity and self-interaction — both discussed in §5.6.)
- *Strongly correlated electrons*: LDA misses Mott insulating gaps entirely; predicts FeO, CoO, and many other transition metal oxides to be metals when they are antiferromagnetic insulators.
- *Van der Waals*: no dispersion at all (no $-C_6/R^{6}$ tail).

LDA is the bottom rung. It is rarely the right choice today, but it is the reference against which all other functionals are calibrated.

## 5.4.3 GGA: gradient corrections

The next rung uses the local density gradient $|\nabla n|$. The reasoning: real systems are not uniform. Adding sensitivity to how rapidly $n$ varies should help, especially in regions of bond formation and at surfaces where the density changes quickly.

A general GGA has the form

$$
E_{xc}^\mathrm{GGA}[n] = \int n(\mathbf r)\,\epsilon_{xc}^\mathrm{unif}\!\big(n(\mathbf r)\big)\,F_{xc}\!\big(n, |\nabla n|\big)\,\mathrm d\mathbf r,
\tag{5.37}
$$

where $F_{xc}$ is a dimensionless **enhancement factor** that depends on the local density and its gradient through the dimensionless **reduced gradient**

$$
s = \frac{|\nabla n|}{2k_F(n)\,n} = \frac{|\nabla n|}{2(3\pi^{2})^{1/3}\,n^{4/3}}.
\tag{5.38}
$$

At $s = 0$ we recover LDA, $F_{xc}(s=0) = 1$. For larger $s$ — bond regions, surface tails — $F_{xc}$ deviates from unity.

### PBE: the workhorse

The Perdew–Burke–Ernzerhof functional (PBE, 1996) is the most widely used GGA in materials science. Its construction is principled: PBE is built from a small set of exact constraints on the exchange–correlation energy, with no fits to experimental data.

The PBE exchange enhancement factor is

$$
F_x^\mathrm{PBE}(s) = 1 + \kappa - \frac{\kappa}{1 + \mu s^{2}/\kappa},
\tag{5.39}
$$

with constants $\mu = 0.21951$ (chosen to recover the linear-response of the uniform gas in the small-$s$ limit, equivalent to second-order gradient expansion) and $\kappa = 0.804$ (chosen to satisfy the **Lieb–Oxford bound** $E_x \geq -1.679\int n^{4/3}$). At small $s$, $F_x \approx 1 + \mu s^{2}$, and at large $s$, $F_x \to 1 + \kappa \approx 1.804$. The correlation part of PBE is similarly built from exact constraints; we shall not reproduce its full form here (the reader can find it in Perdew, Burke, and Ernzerhof, Phys. Rev. Lett. **77**, 3865 (1996)).

!!! note "Why this step?"
    What does the reduced gradient $s = |\nabla n|/(2k_F n)$ measure physically? It is the change in $n$ over a distance $\sim 1/k_F$ (the Fermi wavelength) divided by $n$ itself — a dimensionless measure of how rapidly the density varies on the scale of the local "quantum mechanical wiggle length" of the Fermi sea. In a homogeneous gas, $s=0$. In the tail of an atomic density (where $n\sim \mathrm e^{-\alpha r}$), $s$ grows without bound. In the bond region between two atoms, $s$ is moderate. PBE's enhancement factor saturates at $F_x = 1 + \kappa\approx 1.804$ as $s\to\infty$: an explicit imposition of the Lieb–Oxford bound that prevents the exchange energy from becoming too negative in low-density regions.

!!! example "Worked example: Si band gap with PBE"
    Crystalline silicon has an indirect experimental gap of $E_g^\mathrm{exp} = 1.17\;\text{eV}$. PBE predicts $E_g^\mathrm{PBE}\approx 0.6$–$0.7\;\text{eV}$ — a $40$–$50\%$ underestimate. The error has two roughly equal sources: (i) the missing derivative discontinuity (§5.6), worth $\sim 0.3\;\text{eV}$, and (ii) self-interaction error pushing the HOMO too high and LUMO too low. The PBE *band structure* is qualitatively correct (right ordering of bands, right symmetry character), but the gap is quantitatively wrong. HSE06 below fixes this.

The total PBE functional is

$$
E_{xc}^\mathrm{PBE}[n] = \int n\,\epsilon_x^\mathrm{unif}(n)\,F_x^\mathrm{PBE}(s)\,\mathrm d\mathbf r + E_c^\mathrm{PBE}[n].
$$

Variants — PBEsol (tuned for solids), revPBE, RPBE — adjust the constants for better lattice constants or surface energies. PBE itself is the default in many materials-science DFT codes and the workhorse of large-scale databases like the Materials Project.

!!! note "Why this step?"
    PBE's exchange enhancement factor is *constructed*, not fitted. The constants $\mu$ and $\kappa$ are chosen so that (i) the small-$s$ limit recovers the gradient expansion of the slowly-varying electron gas (constraint on $\mu$), and (ii) the large-$s$ limit respects the Lieb–Oxford bound, the rigorous lower bound on the exchange energy for any density (constraint on $\kappa$). This *constraint-based* construction philosophy — pioneered by Perdew and colleagues — is one of the most enduring contributions to functional design. It produces functionals that, by construction, get certain limits exactly right; what they fail at is often illuminating in the diagnosis sense (it points to a constraint they cannot simultaneously satisfy).

### When GGA helps and when it hurts

GGAs cure the worst of LDA's pathologies:

- **Atomisation energies of molecules**: LDA errors of $\sim 30$ kcal/mol drop to $\sim 8$ kcal/mol with PBE.
- **Lattice constants**: PBE typically gives lattice constants slightly *over* the experimental value, in contrast to LDA's under-estimate. PBEsol corrects this for solids.
- **Surface energies**: GGAs are an improvement, but PBE has a known small *underestimate* — PBEsol again does better.

GGAs do *not* cure:

- **Band gap underestimation**: GGAs lower the LDA band gap further, or leave it essentially unchanged.
- **vdW**: still no dispersion tail.
- **Strongly correlated systems**: still wrong.
- **Self-interaction error**: still present.
- **Barrier heights**: GGAs systematically *underestimate* reaction barriers by 5–10 kcal/mol — a problem for chemistry that hybrids partly fix.

For routine materials calculations, PBE is the modern minimum.

## 5.4.4 meta-GGA: SCAN

The next rung adds the **kinetic energy density**

$$
\tau(\mathbf r) = \tfrac{1}{2}\sum_i^\mathrm{occ}|\nabla\phi_i(\mathbf r)|^{2}.
\tag{5.40}
$$

(Some meta-GGAs use the Laplacian $\nabla^{2}n$ instead; SCAN uses $\tau$.) The new ingredient distinguishes single-orbital regions (where $\tau$ equals the von Weizsäcker bound $\tau_W = |\nabla n|^{2}/(8n)$) from regions of overlapping orbitals (where $\tau$ exceeds $\tau_W$).

The **strongly constrained and appropriately normed (SCAN)** functional of Sun, Ruzsinszky, and Perdew (2015) is built to satisfy all 17 known exact constraints on $E_{xc}$ that can be obeyed by a semi-local functional. It often outperforms PBE on diverse benchmarks — atomic energies, molecular binding, hydrogen-bonded systems, even some weakly bound systems through the implicit treatment of intermediate-range correlation. The cost is roughly the same as a GGA (no exact exchange to evaluate), though convergence can be more delicate due to the more complex functional dependence.

SCAN is increasingly the default for high-accuracy materials calculations where hybrid cost is prohibitive. Variants like r$^{2}$SCAN (Furness et al., 2020) improve numerical stability for solids.

!!! example "TiO$_2$ rutile gap across rungs"
    The optical gap of rutile TiO$_2$ is $E_g^\mathrm{exp}\approx 3.03\;\text{eV}$. DFT values:
    
    | Functional | $E_g$ (eV) | Error |
    |---|---|---|
    | LDA | $\sim 1.7$ | $-1.3$ |
    | PBE | $\sim 2.0$ | $-1.0$ |
    | SCAN | $\sim 2.3$ | $-0.7$ |
    | HSE06 | $\sim 2.7$–$3.0$ | $-0.3$ to $0$ |
    | $G_0W_0$@PBE | $\sim 3.3$ | $+0.3$ |
    
    The progression up Jacob's ladder steadily reduces the gap error from $\sim 1\;\text{eV}$ (LDA/PBE) to $\sim 0.3\;\text{eV}$ (HSE06), at the cost of an order-of-magnitude increase in compute. For materials screening at the gap-prediction level, HSE06 is the workhorse.

## 5.4.5 Hybrid functionals: mixing in exact exchange

The next major rung adds a fraction of **exact (Hartree–Fock) exchange**:

$$
E_x^\mathrm{HF} = -\tfrac{1}{2}\sum_{ij}\iint\frac{\phi_i^{*}(\mathbf r)\phi_j(\mathbf r)\phi_j^{*}(\mathbf r')\phi_i(\mathbf r')}{|\mathbf r - \mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r'.
\tag{5.41}
$$

This is non-local in the orbitals — it depends on $\phi_i(\mathbf r)\phi_i(\mathbf r')$ at *two* points — and is significantly more expensive to compute than any semi-local functional.

### B3LYP

The first widely successful hybrid was Becke's three-parameter mix (1993), most often used with the Lee–Yang–Parr correlation functional:

$$
E_{xc}^\mathrm{B3LYP} = (1-a)E_x^\mathrm{LDA} + a\,E_x^\mathrm{HF} + b\,\Delta E_x^\mathrm{B88} + (1-c)E_c^\mathrm{LDA} + c\,E_c^\mathrm{LYP},
$$

with $a=0.20$, $b=0.72$, $c=0.81$, fit to atomisation energies. B3LYP became the workhorse of computational chemistry — for molecules, organic systems, and biomolecules it gives near-chemical accuracy.

For *solids*, however, B3LYP is problematic: the LYP correlation does not recover the uniform electron gas correlation in the high-density limit, so it misbehaves for metals. In solid-state physics one rarely uses B3LYP.

### HSE06: the screened hybrid for solids

For solids, the dominant hybrid is **HSE** (Heyd–Scuseria–Ernzerhof). The idea is to apply exact exchange *only at short range*, where it matters most for chemical bonding, and use the much cheaper PBE exchange at long range. The Coulomb operator is split via a screening parameter $\omega$:

$$
\frac{1}{r_{12}} = \underbrace{\frac{\mathrm{erfc}(\omega r_{12})}{r_{12}}}_\mathrm{short\;range} + \underbrace{\frac{\mathrm{erf}(\omega r_{12})}{r_{12}}}_\mathrm{long\;range}.
$$

HSE06 mixes 25% exact exchange into the short-range part and uses pure PBE for everything else:

$$
E_{xc}^\mathrm{HSE} = 0.25\,E_x^\mathrm{HF,SR}(\omega) + 0.75\,E_x^\mathrm{PBE,SR}(\omega) + E_x^\mathrm{PBE,LR}(\omega) + E_c^\mathrm{PBE},
$$

with $\omega = 0.11$ a.u.$^{-1}$. HSE06 dramatically improves band gaps relative to PBE — typical errors drop from $\sim 1$ eV to $\sim 0.3$ eV — while the screened exchange makes it tractable in metallic and small-gap systems where pure global hybrids (PBE0, B3LYP) develop convergence pathologies.

### Cost of hybrids

Computing $E_x^\mathrm{HF}$ requires evaluating four-centre integrals (or, equivalently, doing Fourier transforms with all pairs of occupied orbitals). For plane-wave codes, the dominant cost scales as $\mathcal O(N_\mathrm{occ}^{2}N_\mathrm{plane-wave}\log N_\mathrm{plane-wave})$ per SCF step — typically 10–30 times more expensive than a GGA calculation on the same system. For very large systems or molecular dynamics, this can be prohibitive; HSE06's range-separation softens but does not eliminate the cost.

### Double hybrids: rung 5

Double-hybrid functionals (B2PLYP, XYG3, DSD-BLYP) add not only exact exchange but also a fraction of *second-order perturbation* (MP2-style) correlation, computed from the virtual orbitals of a hybrid KS calculation:

$$
E_{xc}^\mathrm{DH} = (1-a)E_x^\mathrm{GGA} + a\,E_x^\mathrm{HF} + (1-b)E_c^\mathrm{GGA} + b\,E_c^\mathrm{MP2},
$$

with typical mixing fractions $a\sim 0.5$–$0.7$, $b\sim 0.3$–$0.5$. The MP2 correlation costs $\mathcal O(N^{5})$ to evaluate, an order of magnitude more than pure hybrids. The reward is near-chemical accuracy ($\sim 1\;\text{kcal/mol}$) for thermochemistry, reaction barriers, and noncovalent interactions. For very accurate quantum chemistry on systems of $\sim 30$–$50$ atoms, double hybrids are now the workhorse alongside DLPNO-CCSD(T).

!!! note "Why this step?"
    The fifth rung makes contact with traditional wavefunction-based quantum chemistry: it uses information about *unoccupied* (virtual) KS orbitals, not just occupied ones. This is a fundamentally different ingredient — it depends on the entire spectrum of the KS Hamiltonian, not just the ground-state density. The price is that one loses the strict density-functional pedigree (the functional depends explicitly on orbitals, not just $n$). The reward is access to the missing long-range correlation (dispersion, $-C_6/R^{6}$) that semi-local functionals cannot reach.

## 5.4.6 Van der Waals corrections

A defining failure of all semi-local (LDA, GGA, meta-GGA) functionals — and of hybrids that mix only short-range exact exchange — is that they have *no* $-C_6/R^{6}$ dispersion attraction between non-overlapping fragments. London dispersion is a fundamentally non-local correlation effect: instantaneous dipole fluctuations on one fragment induce dipoles on the other, with energy $\sim -\alpha_A\alpha_B/R^{6}$ for polarisabilities $\alpha_A,\alpha_B$. Where the densities do not overlap, the local functionals see nothing.

For a long list of important systems — molecular crystals, layered materials (graphite, transition metal dichalcogenides), surface adsorption of organics, biological molecules — this matters quantitatively.

Several pragmatic fixes are available.

**D3 / D4 (Grimme).** Add an explicit pairwise correction:

$$
E_\mathrm{disp} = -\sum_{A<B}\Big[s_6\frac{C_6^{AB}}{R_{AB}^{6}}f_6(R_{AB}) + s_8\frac{C_8^{AB}}{R_{AB}^{8}}f_8(R_{AB})\Big].
$$

The $C_n^{AB}$ coefficients are pre-tabulated (D3) or made geometry-dependent through fractional coordination numbers (D4); the damping functions $f_n$ kill the divergence at short range; the scaling factors $s_n$ are fit per functional. D3 and D4 are essentially free to compute and improve binding of dispersion-bound systems by orders of magnitude.

**Tkatchenko–Scheffler (TS).** Like D3 but with $C_6$ coefficients computed *self-consistently* from the actual electron density via Hirshfeld partitioning. Captures environment dependence of $C_6$ better than D3's tabulated values; modest extra cost.

**Non-local vdW functionals (vdW-DF1, vdW-DF2, rVV10, MBD).** Add a non-local correlation kernel directly to the functional:

$$
E_c^\mathrm{nl}[n] = \tfrac{1}{2}\iint n(\mathbf r)\,\Phi(\mathbf r,\mathbf r')\,n(\mathbf r')\,\mathrm d\mathbf r\,\mathrm d\mathbf r',
$$

with $\Phi$ a kernel encoding the dispersion physics. These can be implemented efficiently via FFT and are routinely available in plane-wave codes.

For materials applications today, *not* applying some kind of dispersion correction when a system has non-bonded fragments is a methodological error.

## 5.4.7 Self-interaction error

The Hartree energy (5.26) integrates $n(\mathbf r)n(\mathbf r')/|\mathbf r-\mathbf r'|$ over the *entire* density — including, for a single electron, an electron's interaction with its own charge distribution. Exact exchange (5.41) for a single electron exactly cancels this spurious self-Hartree, but approximate exchange functionals do *not*: the LDA or GGA exchange of a hydrogen atom does not fully cancel its self-Hartree. The leftover is **self-interaction error** (SIE).

### Explicit demonstration: hydrogen atom

For the exact ground state of hydrogen ($N=1$), the only electron cannot interact with itself, so the *exact* Hartree-plus-exchange-correlation energy must vanish:

$$
U_H[n_\mathrm H] + E_{xc}^\mathrm{exact}[n_\mathrm H] \;=\; 0.
\tag{5.43}
$$

The hydrogenic density is $n_\mathrm H(\mathbf r) = (1/\pi)\,\mathrm e^{-2r}$ in atomic units. The Hartree integral can be done in closed form:

$$
U_H[n_\mathrm H] = \tfrac{1}{2}\iint\frac{n_\mathrm H(\mathbf r)n_\mathrm H(\mathbf r')}{|\mathbf r - \mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r' = \tfrac{5}{16}\;\text{Ha} = 0.3125\;\text{Ha}.
$$

For (5.43) to hold, the exact $E_{xc}[n_\mathrm H] = -0.3125\;\text{Ha}$. The LDA exchange is

$$
E_x^\mathrm{LDA}[n_\mathrm H] = -\tfrac{3}{4}(3/\pi)^{1/3}\int n_\mathrm H^{4/3}\,\mathrm d\mathbf r \approx -0.260\;\text{Ha},
$$

leaving an LDA correlation contribution of $\approx -0.04\;\text{Ha}$ and a total $U_H + E_{xc}^\mathrm{LDA}\approx 0.013\;\text{Ha}\approx 0.34\;\text{eV}$ — *not* zero. The leftover $\sim 0.34\;\text{eV}$ is the LDA self-interaction error for hydrogen, and it is what causes the LDA HOMO of H to be too high in energy by several eV (the binding of the electron to the nucleus is artificially weakened by this residual self-repulsion).

The corresponding calculation for PBE gives a slightly smaller residue ($\approx 0.005\;\text{Ha}\approx 0.14\;\text{eV}$); for HSE06 it is even smaller. For pure Hartree–Fock, the residue is *exactly* zero by construction — but at the cost of missing all correlation, which gives Hartree–Fock its own characteristic errors.

!!! note "Why this step?"
    The single-electron self-interaction is the *cleanest* diagnostic of approximate functionals because we know exactly what the answer should be (zero). Many-electron self-interaction is harder to define rigorously but has the same flavour: in a generic many-electron system, the Hartree integral over the density of any *single occupied orbital* should be exactly cancelled by an exchange contribution involving that same orbital. LDA/GGA exchange does this approximately at best.

SIE has well-known consequences:

- *Over-delocalisation*: electrons artificially spread out (e.g., LDA breaks H$_2^{+}$ symmetry wrongly at large bond distance; gives fractional charges on dissociating molecules).
- *Underestimated band gaps* in semiconductors (the HOMO is too high; the LUMO is too low).
- *Bad treatment of polarons and small radicals*.

Hybrid functionals partially cure SIE because their exact-exchange fraction cancels part of the self-Hartree. Range-separated hybrids, the optimised effective potential method, self-interaction-corrected (Perdew–Zunger SIC) functionals, and DFT+U each tackle SIE from a different angle. The problem is fundamental to local and semi-local exchange and is the deepest reason for the band gap problem (§5.6).

### The derivative discontinuity and the band gap problem

A closely related issue is the **derivative discontinuity** of the exchange–correlation potential. The exact $v_{xc}$ jumps by a uniform constant $\Delta_{xc}$ as the electron number passes through an integer. The fundamental gap of an $N$-electron system equals

$$
E_g = (\varepsilon_\mathrm{LUMO} - \varepsilon_\mathrm{HOMO}) + \Delta_{xc},
$$

with the *KS gap* (the eigenvalue difference) generally smaller than the fundamental gap. For LDA, GGA, and SCAN, $\Delta_{xc}\equiv 0$ — the local potential is a smooth function of $n$ at integer occupations. So the reported band gap equals the KS gap, missing the (typically $\sim 1\;\text{eV}$) derivative-discontinuity contribution. Hybrid functionals partially restore $\Delta_{xc}$ through their non-local exchange piece, which is why HSE06 band gaps are quantitatively much closer to experiment.

The deep connection between SIE and the derivative discontinuity is that *both* are signatures of the failure of approximate functionals to obey the *piecewise-linearity* of the exact $E(N)$ between integer electron numbers (§5.3.4). Indeed, the magnitude of the derivative discontinuity at $N$ equals the discontinuity of $\partial E/\partial N$ at $N$:

$$
\Delta_{xc} = \lim_{\eta\to 0^{+}}\big[\partial E/\partial N|_{N+\eta} - \partial E/\partial N|_{N-\eta}\big] - (\varepsilon_\mathrm{LUMO} - \varepsilon_\mathrm{HOMO}).
$$

For LDA/GGA, $E(N)$ is smooth (no jump), so $\Delta_{xc}=0$. Restoring piecewise linearity — via exact exchange, optimal tuning, or many-body perturbation theory — also restores the derivative discontinuity and fixes band gaps.

## 5.4.7a Decision flowchart for choosing a functional

<figure markdown>
```mermaid
flowchart TD
    A[Start: what system?] --> B{Periodic solid?}
    B -- Yes --> C{Strongly correlated?<br/>3d/4f transition metal?}
    B -- No --> M{Molecule}
    C -- Yes --> D[DFT+U / HSE06 / DMFT]
    C -- No --> E{Need band gap?}
    E -- Yes --> F[HSE06 + dispersion]
    E -- No --> G{Lattice constants critical?}
    G -- Yes --> H[PBEsol or SCAN]
    G -- No --> I[PBE + D3]
    M --> N{Reaction barriers / thermochemistry?}
    N -- Yes --> O[B3LYP-D3, ωB97X-D, or M06-2X]
    N -- No --> P{Charge transfer / Rydberg?}
    P -- Yes --> Q[Range-separated: CAM-B3LYP, ωB97X]
    P -- No --> R[B3LYP + D3]
```
<figcaption>A practical decision tree for choosing an exchange–correlation functional. The first branch asks whether the system is a periodic solid or a molecule; periodic solids then branch on strong correlation (pointing to DFT+U, HSE06 or DMFT), on whether a band gap is needed (HSE06 plus dispersion), and on whether lattice constants are critical (PBEsol or SCAN, otherwise PBE+D3); molecules branch on reaction barriers and thermochemistry (B3LYP-D3, ωB97X-D or M06-2X) and on charge-transfer or Rydberg character (range-separated functionals, otherwise B3LYP+D3). Always include a dispersion correction such as D3, D4 or vdW-DF when non-bonded fragments are present.</figcaption>
</figure>

## 5.4.7b Benchmarks: a tour across the ladder

To make the rungs concrete, here are typical performance numbers for three representative quantities — atomic ionisation potentials, molecular atomisation energies, and semiconductor band gaps — averaged over standard benchmark sets:

| Functional | IP error (eV) | Atomisation error (kcal/mol) | Gap error (eV) | Cost (rel. PBE) |
|---|---|---|---|---|
| LDA (PW92) | $\sim 0.5$ | $\sim 35$ (overbind) | $-0.9$ | 0.9 |
| PBE | $\sim 0.4$ | $\sim 8$ (overbind) | $-0.8$ | 1.0 |
| PBEsol | $\sim 0.5$ | $\sim 10$ | $-0.8$ | 1.0 |
| SCAN | $\sim 0.2$ | $\sim 5$ | $-0.6$ | 1.5 |
| B3LYP | $\sim 0.2$ | $\sim 3$ | n/a (mol.) | 15 |
| HSE06 | $\sim 0.2$ | $\sim 4$ | $-0.3$ | 20 |
| $\omega$B97X-D | $\sim 0.15$ | $\sim 2$ | n/a | 25 |
| double-hybrid (B2PLYP) | $\sim 0.1$ | $\sim 2$ | n/a | 100 |
| $G_0W_0$@PBE | n/a | n/a | $+0.2$ | 100 |

Two patterns emerge. First, the *atomisation energy* of molecules — the gold-standard chemical accuracy benchmark — improves monotonically up the ladder, from $\sim 35\;\text{kcal/mol}$ (LDA) to $\sim 2\;\text{kcal/mol}$ (double-hybrid). Second, the *band gap* of solids requires the non-local exchange of a hybrid before substantial improvement is seen. For high-accuracy gap predictions, $G_0W_0$ on a hybrid starting point is the practical state of the art.

!!! example "Graphite interlayer binding"
    The experimental interlayer binding energy of graphite is $\sim 50\;\text{meV/atom}$ and the interlayer spacing is $3.35\;\text{Å}$. Computational predictions:
    
    | Method | Binding (meV/atom) | Spacing (Å) |
    |---|---|---|
    | PBE | $\sim 1$ | $\sim 4.4$ (way off) |
    | PBE+D3 | $\sim 70$ (overbind) | $3.30$ |
    | SCAN | $\sim 50$ | $3.32$ |
    | optB88-vdW | $\sim 55$ | $3.34$ |
    | RPA | $\sim 50$ | $3.35$ |
    | QMC | $\sim 50$ | $3.35$ |
    
    Pure PBE is *catastrophically wrong* — graphite is barely bound and the layers float to twice the experimental spacing. Any dispersion correction or vdW functional fixes the binding; SCAN does so without any explicit dispersion correction thanks to its built-in intermediate-range correlation.

## 5.4.8 Which functional should I use?

There is no single answer. Match the tool to the question.

| Question / system | Sensible default | Cost |
|---|---|---|
| Bulk metals, simple oxides, geometries | PBE (or PBEsol for lattice constants) | low |
| General materials screening | PBE + D3 | low |
| Organic molecules, gas-phase chemistry | B3LYP + D3 (or $\omega$B97X-D) | medium |
| Semiconductor band gaps | HSE06 | high |
| Magnetic transition-metal oxides | DFT+U (PBE+U) or HSE06 | low / high |
| Layered materials, molecular crystals | PBE + D3, optB88-vdW, or SCAN | low |
| Surface adsorption (chemical) | RPBE / BEEF-vdW | low |
| Hydrogen bonding, water | SCAN, revPBE+D3 | low |
| Strongly correlated (Mott) | DFT+U, hybrid, or DMFT (Ch. 5.6) | low / high |
| Excited states, optical absorption | TD-DFT with a hybrid, BSE | high |
| Reaction barriers in chemistry | Hybrid (B3LYP, M06-2X, $\omega$B97X) | high |

!!! example "How much does the functional cost in absolute terms?"
    For a 100-atom supercell on a 16-core CPU node:
    
    - PBE: $\sim 30\;\text{min}$ per SCF.
    - SCAN: $\sim 45\;\text{min}$.
    - HSE06: $\sim 10\;\text{h}$.
    - PBE + D3: PBE time + $<1\;\text{s}$ for D3.
    
    For high-throughput screening of $10^{4}$–$10^{5}$ structures (typical for materials databases), only PBE and PBEsol are economically feasible at scale. Hybrids are reserved for "second-pass" refinement on shortlisted candidates. This is why every major materials database (Materials Project, AFLOW, OQMD) is built on PBE: a deliberate accuracy-throughput trade-off.

A few rules of thumb:

- *Always include a dispersion correction* (D3 or vdW-DF) for any system with non-bonded fragments. The cost is negligible.
- *Try at least two functionals.* If your conclusion changes between PBE and PBE+D3, or between PBE and HSE06, your result is functional-sensitive and you should report both.
- *For high-throughput databases*, the answer is usually PBE — partly because the data was generated with PBE, and consistency matters more than absolute accuracy.
- *For band gaps, never trust LDA or GGA at face value.* Use HSE06 or correct with a GW calculation (Chapter 5.6).

In Chapter 9 we shall see that machine-learning interatomic potentials inherit, in a precise sense, the errors of the functional they are trained on. A model trained on PBE energies will reproduce PBE bond lengths, including PBE's slight systematic overestimate. Awareness of the functional is therefore not just a methodological nicety; it propagates into every downstream tool that consumes DFT data.

### Functional choice and machine-learning surrogates

Modern materials informatics increasingly relies on machine-learning interatomic potentials (MLIPs) trained on DFT data (Chapter 9). The choice of functional in the training data has subtle but important consequences for the trained model. A model trained on PBE energies will reproduce PBE bond lengths (slightly long), PBE band gaps (too small), PBE-overbound dispersion (typically absent). Switching the training-set functional from PBE to SCAN typically improves the accuracy of derived properties — but at $\sim 1.5\times$ the DFT cost during data generation.

A few practical recommendations for ML-potential training-set design:

- *Be consistent.* Mix-and-match between PBE and HSE06 energies in the same training set produces incoherent models. Pick one functional and stick with it.
- *Include dispersion.* If the system has any non-bonded fragments, use PBE+D3 or SCAN throughout.
- *Match the deployment context.* If your downstream task involves predicting band gaps, train on HSE06. If it involves predicting lattice constants for a wide range of materials, PBEsol or SCAN is a good choice.
- *Document the functional explicitly.* The functional choice is part of the model card.

In Chapter 9 we shall return to this point in the context of training set construction for crystal-property prediction, where the choice of reference functional shapes the achievable accuracy of any downstream ML model.

### Summary of §5.4 — what to remember in 3 months

- **Jacob's ladder**: LDA → GGA → meta-GGA → hybrid → double-hybrid. Each rung adds an ingredient and ~10× cost.
- **LDA exchange (Dirac)**: $E_x^\mathrm{LDA} = -\tfrac{3}{4}(3/\pi)^{1/3}\int n^{4/3}\,\mathrm d\mathbf r$. Derived from the uniform electron gas.
- **PBE**: GGA built from exact constraints (uniform-gas limit, Lieb–Oxford bound). The workhorse of materials science.
- **SCAN**: meta-GGA satisfying 17 exact constraints. Often competitive with hybrids at GGA cost.
- **HSE06**: screened hybrid with 25% short-range exact exchange. Standard for solid-state band gaps.
- **B3LYP**: hybrid empirically fit; chemistry standard but unreliable for solids.
- **Self-interaction error (SIE)**: spurious self-Hartree not cancelled by approximate exchange. The deepest pathology of LDA/GGA, root cause of band-gap underestimation and over-delocalisation.
- **Derivative discontinuity**: exact $v_{xc}$ jumps by $\Delta_{xc}$ at integer $N$; LDA/GGA have $\Delta_{xc}=0$.
- **vdW dispersion**: $-C_6/R^{6}$ tail is missing from all semi-local functionals. Add D3/D4 or use vdW-DF/SCAN+rVV10.
- **No universal best functional**: choose based on system and property of interest (see decision table).

!!! note "Remark: how many functionals are there?"
    The Libxc library, the standard open-source functional library used by many DFT codes, currently contains over 600 distinct functionals. In practice, fewer than 20 are in routine use; the rest are historical curiosities, narrow-purpose variants, or recent proposals awaiting community testing. The bewildering number is one reason functional choice has become a methodological topic in its own right.

We have, in this section:

- Derived LDA exchange from the uniform electron gas.
- Stated PBE and explained its enhancement factor.
- Surveyed meta-GGAs (SCAN), hybrids (B3LYP, HSE06), van der Waals corrections, and self-interaction.
- Distilled a practical decision table.

We have not solved any actual equations. Choosing a functional gives us $v_{xc}[n](\mathbf r)$; what we now need is an algorithm for solving the Kohn–Sham equations self-consistently with that $v_{xc}$. That algorithm — the self-consistent field loop — and a complete Python implementation are the subject of §5.5.
