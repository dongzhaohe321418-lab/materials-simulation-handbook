# 5.4 Exchange–Correlation Functionals

The Kohn–Sham construction (§5.3) is exact in principle. It becomes an approximation only because we do not know the exchange–correlation energy functional $E_{xc}[n]$ in closed form. Every DFT calculation in the world today involves a choice of approximate $E_{xc}$: a *functional*. The choice matters. The same molecule can have its bond length reproduced within 1 pm by one functional and missed by 10 pm by another; a band gap can come out qualitatively right or qualitatively wrong; a magnetic ground state can flip.

The functional zoo is large — thousands have been proposed; perhaps fifty are in common use. To navigate it, John Perdew suggested a marvellous metaphor: **Jacob's ladder**. Each rung adds an ingredient and, on average, climbs toward chemical accuracy ($\sim 1$ kcal/mol $\approx 0.04$ eV). Each rung also costs more. We climb the ladder in turn.

## 5.4.1 Jacob's ladder

From bottom to top:

1. **LDA**: uses $n(\mathbf r)$ only.
2. **GGA**: uses $n$ and $|\nabla n|$.
3. **meta-GGA**: adds the kinetic energy density $\tau(\mathbf r)$ or $\nabla^{2}n$.
4. **Hybrid**: mixes in a fraction of exact (Hartree–Fock) exchange.
5. **Double hybrid** / RPA / wavefunction methods: include unoccupied orbitals.

The first three rungs are purely *semi-local* — the value of $\epsilon_{xc}$ at $\mathbf r$ depends only on quantities at $\mathbf r$ (or its immediate gradients). The fourth introduces non-locality through exact exchange and is roughly an order of magnitude more expensive. The fifth introduces virtual orbitals and another order of magnitude.

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

### Where LDA works, where it fails

**Strengths.**

- *Free-electron-like systems*: bulk metals (Na, Al, Cu) — bond lengths and bulk moduli within a few per cent. The uniform-gas reference is a reasonable starting point when the density is genuinely slowly varying.
- *Total energies*: the cohesive energy of a metal comes out plausibly, often within 0.5 eV/atom.
- *Geometries*: equilibrium bond lengths in simple solids are usually within 1–2% of experiment.

**Weaknesses.**

- *Overbinding*: LDA systematically overestimates binding energies, often by tens of per cent. The H$_2$ binding energy is 4.91 eV experimentally; LDA gives 4.79 eV with a too-short bond, then is unreliable for everything bigger.
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

The total PBE functional is

$$
E_{xc}^\mathrm{PBE}[n] = \int n\,\epsilon_x^\mathrm{unif}(n)\,F_x^\mathrm{PBE}(s)\,\mathrm d\mathbf r + E_c^\mathrm{PBE}[n].
$$

Variants — PBEsol (tuned for solids), revPBE, RPBE — adjust the constants for better lattice constants or surface energies. PBE itself is the default in many materials-science DFT codes and the workhorse of large-scale databases like the Materials Project.

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

SIE has well-known consequences:

- *Over-delocalisation*: electrons artificially spread out (e.g., LDA breaks H$_2^{+}$ symmetry wrongly at large bond distance; gives fractional charges on dissociating molecules).
- *Underestimated band gaps* in semiconductors (the HOMO is too high; the LUMO is too low).
- *Bad treatment of polarons and small radicals*.

Hybrid functionals partially cure SIE because their exact-exchange fraction cancels part of the self-Hartree. Range-separated hybrids, the optimised effective potential method, self-interaction-corrected (Perdew–Zunger SIC) functionals, and DFT+U each tackle SIE from a different angle. The problem is fundamental to local and semi-local exchange and is the deepest reason for the band gap problem (§5.6).

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

A few rules of thumb:

- *Always include a dispersion correction* (D3 or vdW-DF) for any system with non-bonded fragments. The cost is negligible.
- *Try at least two functionals.* If your conclusion changes between PBE and PBE+D3, or between PBE and HSE06, your result is functional-sensitive and you should report both.
- *For high-throughput databases*, the answer is usually PBE — partly because the data was generated with PBE, and consistency matters more than absolute accuracy.
- *For band gaps, never trust LDA or GGA at face value.* Use HSE06 or correct with a GW calculation (Chapter 5.6).

In Chapter 9 we shall see that machine-learning interatomic potentials inherit, in a precise sense, the errors of the functional they are trained on. A model trained on PBE energies will reproduce PBE bond lengths, including PBE's slight systematic overestimate. Awareness of the functional is therefore not just a methodological nicety; it propagates into every downstream tool that consumes DFT data.

## 5.4.9 Looking ahead

We have, in this section:

- Derived LDA exchange from the uniform electron gas.
- Stated PBE and explained its enhancement factor.
- Surveyed meta-GGAs (SCAN), hybrids (B3LYP, HSE06), van der Waals corrections, and self-interaction.
- Distilled a practical decision table.

We have not solved any actual equations. Choosing a functional gives us $v_{xc}[n](\mathbf r)$; what we now need is an algorithm for solving the Kohn–Sham equations self-consistently with that $v_{xc}$. That algorithm — the self-consistent field loop — and a complete Python implementation are the subject of §5.5.
