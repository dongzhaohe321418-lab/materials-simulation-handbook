# 5.6 Where DFT Fails

!!! note "Why does this chapter exist?"
    Every powerful tool has a comfort zone and a failure zone. A hammer drives nails magnificently but is hopeless for screws; a microscope reveals cells but cannot see atoms. Kohn–Sham DFT is the same: extraordinary inside its comfort zone, dangerous outside. The comfort zone is *most of materials science* — bond lengths, lattice parameters, vibrational spectra, cohesive energies of typical materials. The failure zone is a long, dispiriting list: van der Waals interactions, band gaps, Mott insulators, charge transfer, excited states, dissociation curves of stretched bonds.
    
    What makes this chapter important is that the failures are not random. They cluster around well-understood physical mechanisms — *self-interaction error*, *missing non-locality*, *derivative discontinuity* of the exact potential, *static correlation*. If you can recognise the symptoms (PBE predicting graphite to float, PBE predicting NiO to be a metal, PBE collapsing the dissociation of H$_2^{+}$), you can diagnose the cause and reach for the right remedy: DFT+U for Mott systems, hybrids for band gaps, vdW corrections for dispersion, GW for quasiparticle energies, DMFT for strong correlation.
    
    A useful analogy. A general practitioner can treat 90% of medical complaints; for the other 10% you need a specialist. DFT is the general practitioner of materials science; this section is your reference for when to call the specialist (CCSD(T), QMC, GW+BSE, DMFT) and what they cost. The crucial skill is knowing what your tool cannot do — and being honest about it in your published work.

DFT is the workhorse of computational materials science because it is, on the whole, *good enough* — it gets bond lengths within a few per cent, cohesive energies within tens of kJ/mol, vibrational frequencies within 10%, and the qualitative ground-state physics right for the vast majority of systems. But "the vast majority" leaves a long tail of systems and properties where standard Kohn–Sham DFT, in any practical approximation, is simply wrong. Some of these failures are technical (chosen functional too crude); others are fundamental (no semi-local functional can fix them).

This section gives an honest tour. Knowing where DFT breaks is the difference between a trustworthy calculation and a published mistake. For each failure mode we identify the symptom, the underlying physics, and the higher-level methods one reaches for instead.

**What problem are we solving?** DFT is the workhorse of the field — but a workhorse is only safe in the hands of someone who knows where it stumbles. A calculation that has converged tidily and printed a clean number can still be *systematically* wrong, and it will not warn you. The purpose of this section is to make you a responsible user: someone who recognises the handful of situations where standard DFT gives a confident but incorrect answer, so you are never fooled into trusting a number just because the code ran without complaint.

!!! note "In plain language"
    DFT fails in a small number of recognisable *families*, and it is worth carrying them in your head:

    - **Band gaps come out too small.** Plain LDA/GGA functionals underestimate semiconductor and insulator gaps, sometimes badly enough to predict a metal where there is none.
    - **Van der Waals (dispersion) forces go missing.** The gentle "stickiness" between non-bonded fragments — graphite layers, molecular crystals, noble-gas dimers — is a long-range correlation effect that a local functional simply cannot see.
    - **Strong correlation breaks the picture.** In some transition-metal oxides (NiO, CoO, FeO) and $f$-electron systems, electrons localise so strongly that the single-orbital Kohn–Sham picture is the wrong starting point, and DFT can wrongly predict a metal.
    - **Self-interaction (delocalisation) error.** An electron spuriously feels its own charge, so functionals tend to over-spread electrons — visible most starkly when a stretched bond dissociates into unphysical fractional charges.
    - **Excited states are out of reach.** Kohn–Sham DFT is a *ground-state* theory by construction; optical absorption, fluorescence and the like need a different framework (TD-DFT, GW+BSE).

    None of these are bugs. They are the known limits of an approximate ground-state theory, and each has a known remedy described below.

!!! abstract "Key idea (Chapter 5.6)"
    Approximate Kohn–Sham DFT fails for a small but important set of systems and properties: band gaps (LDA/GGA underestimate by 30–100%); van der Waals binding (no $-C_6/R^{6}$ tail); strongly correlated electrons (Mott insulators predicted as metals); charge-transfer states (fractional charges at dissociation); excited states (a fundamentally ground-state theory). Each failure has a known physical mechanism — self-interaction error, missing derivative discontinuity, missing non-local correlation — and a known remedy: hybrids, DFT+U, GW, DMFT, CCSD(T). Knowing the failure modes is as important as knowing the theory.

## 5.6.1 Band gap underestimation

**The symptom.** Compute the band structure of silicon with PBE. The fundamental gap comes out around 0.6 eV. The experimental value is 1.17 eV. For germanium, PBE predicts a *metal* with no gap at all; the experimental gap is 0.66 eV. Across the semiconductors, LDA and GGA underestimate band gaps by 30–100%, with the error worse for narrow-gap and small-gap systems.

This is sometimes called *the band gap problem*. It has two distinct contributors, and conflating them has caused decades of confusion.

### The derivative discontinuity

The fundamental gap of an $N$-electron system is

$$
E_g = I - A = [E(N-1) - E(N)] - [E(N) - E(N+1)],
\tag{5.45}
$$

the difference between the ionisation potential and the electron affinity. In *exact* Kohn–Sham theory, the gap relates to the KS eigenvalues by

$$
E_g = \varepsilon_\mathrm{LUMO}^{N} - \varepsilon_\mathrm{HOMO}^{N} + \Delta_{xc},
\tag{5.46}
$$

where $\Delta_{xc}$ is the **derivative discontinuity** of the exchange–correlation potential: as the total electron number passes through an integer $N$, the exact $v_{xc}(\mathbf r)$ jumps by a *uniform constant* $\Delta_{xc}$. The KS eigenvalue difference $\varepsilon_\mathrm{LUMO}-\varepsilon_\mathrm{HOMO}$ is called the *KS gap*; the true *fundamental gap* exceeds it by $\Delta_{xc}$.

For LDA, GGA, and meta-GGA functionals, $v_{xc}$ is a smooth function of $n$ at integer occupation — there is *no* derivative discontinuity. $\Delta_{xc}^\mathrm{LDA} = \Delta_{xc}^\mathrm{GGA} = 0$. So the KS gap is reported as the band gap, and is missing a structural piece that, for real materials, is of order 0.5–2 eV.

### Self-interaction error

Even setting aside the derivative discontinuity, LDA/GGA *Kohn–Sham gaps themselves* are too small, because of self-interaction error (§5.4). Self-interaction artificially raises the HOMO (electrons see their own Coulomb repulsion) and lowers the LUMO. Both effects shrink the apparent gap.

### What to do

- **HSE06** and other range-separated hybrids partially restore the derivative discontinuity through their exact-exchange fraction. Typical gap errors drop to 0.3 eV.
- **GW** (Green's function method, named after the product of the Green's function $G$ and the screened interaction $W$ in Hedin's equations) is the next step up: a many-body perturbation theory correction to the KS quasiparticle energies. $G_0 W_0$ on top of a PBE calculation typically gives gaps within 0.1–0.3 eV of experiment. Cost is $\mathcal O(N^{4})$.
- **Δ-SCF** for small molecules: separately compute the $(N\!-\!1)$ and $(N\!+\!1)$ systems and take the energy difference. Cheap; surprisingly accurate.

!!! example "Worked example: TiO$_2$ rutile gap"
    Experimental optical gap: $E_g\approx 3.03\;\text{eV}$.
    
    | Method | $E_g$ (eV) | Cost (rel. PBE) |
    |---|---|---|
    | LDA | $\sim 1.7$ | 0.9 |
    | PBE | $\sim 2.0$ | 1 |
    | SCAN | $\sim 2.3$ | 1.5 |
    | HSE06 | $\sim 2.7$ | 20 |
    | $G_0W_0$@PBE | $\sim 3.3$ | 100 |
    | scGW@PBE | $\sim 3.1$ | 500 |
    
    The PBE gap error of $\sim 1\;\text{eV}$ is roughly half attributable to the derivative discontinuity (which PBE sets to zero) and half to self-interaction error in the Ti $3d$ states. HSE06 fixes about $70\%$ of the gap error through its 25% exact exchange; $G_0W_0$ corrects the remainder via dynamical screening. The price is a $\sim 100\times$ cost increase for $G_0W_0$ relative to PBE.

!!! warning "Common misunderstanding: the underestimated gap is not a bug"
    A natural first reaction to a too-small gap is to assume something went wrong: a loose convergence threshold, too few $k$-points, an unconverged SCF. **It is none of these.** Tightening every numerical setting will not move a PBE gap towards experiment. The underestimate is *intrinsic* to the theory, for two reasons. First, the quantity you read off is a difference of **Kohn–Sham eigenvalues**, which belong to the *auxiliary non-interacting system*, not to the real interacting electrons. Second, the exact fundamental gap contains the **derivative-discontinuity** contribution $\Delta_{xc}$ of Eq. (5.46), and any semi-local functional (LDA, GGA, meta-GGA) has $\Delta_{xc}=0$ — so it omits a structural piece of the gap of order $0.5$–$2$ eV. The cure is therefore a *better physical model*, not better numerics: a hybrid such as HSE06 restores part of the discontinuity, and **GW** (named above) supplies the proper quasiparticle correction.

!!! warning "Do not over-interpret PBE band gaps"
    A PBE band gap is not "the band gap". It is the Kohn–Sham gap of a particular approximate functional. For predictions of optical or transport gaps, use HSE06 or GW; for ordering of mid-band features, PBE often suffices. Always state the functional alongside the gap.

!!! note "Kohn–Sham eigenvalues are interpretive, not physical excitation energies"
    The eigenvalues $\varepsilon_i$ of the Kohn–Sham equations are genuinely useful for interpretation — they order the states, sketch the band structure, and locate features — but they are properties of the fictitious non-interacting system (Section 5.3), not, in general, true excitation or removal energies. The one rigorous exception is the highest occupied eigenvalue, which in *exact* KS theory equals minus the ionisation potential. Treat the rest as a guide, and obtain real excitation energies from a method designed for them (GW, $\Delta$-SCF, TD-DFT).

## 5.6.2 Van der Waals dispersion

**The symptom.** Stack two graphene sheets at 3.35 Å, the experimental interlayer spacing of graphite. Compute the binding energy with PBE: about 1 meV/atom, essentially zero. The experimental value is around 50 meV/atom. PBE predicts graphite to be barely bound, when it is a robust layered solid.

Or: try a benzene dimer. PBE gives no binding. Or: rare-gas dimers — argon, krypton, xenon — all repulsive with semi-local DFT.

**The physics.** London dispersion forces arise from instantaneous quantum fluctuations of the charge density on one fragment polarising another. The induced dipole pair gives the famous $-C_6/R^{6}$ attraction at large separation. This is a long-range correlation effect: the densities of the two fragments do not overlap, so any *local* functional sees nothing happening between them. Semi-local exchange-correlation, by construction, cannot reproduce $-C_6/R^{6}$.

**What to do.**

- **DFT-D3 / D4** (Grimme): add an empirical pairwise correction. Cheap and effective for most systems.
- **vdW-DF / vdW-DF2 / rVV10**: a non-local correlation kernel built into the functional. Computationally tractable via FFT; available in most plane-wave codes.
- **Tkatchenko–Scheffler / MBD**: density-dependent dispersion coefficients, including many-body screening effects (MBD: many-body dispersion). Best-in-class for systems where polarisability matters.

For materials with non-bonded fragments — molecular crystals, layered materials, surface adsorption, polymers, biomolecules — *not* including a vdW correction in DFT is a methodological error. Modern best practice always includes one.

!!! example "Worked example: graphite interlayer spacing"
    Experimental: $c$-axis layer spacing $\approx 3.35\;\text{Å}$, binding energy $\approx 50\;\text{meV/atom}$.
    
    | Method | $c$ (Å) | $E_b$ (meV/atom) |
    |---|---|---|
    | PBE | $4.4$ (no minimum) | $\sim 1$ |
    | LDA | $3.32$ | $20$ (accidental binding) |
    | PBE+D3 | $3.30$ | $70$ (slight overbind) |
    | optB88-vdW | $3.34$ | $55$ |
    | SCAN+rVV10 | $3.33$ | $52$ |
    | RPA | $3.35$ | $50$ |
    | DMC | $3.35$ | $48\pm 2$ |
    
    Pure PBE is *qualitatively wrong*: it predicts essentially no binding and a layer spacing $\sim 30\%$ too large. LDA accidentally binds graphite at the right spacing because its overbinding bias compensates for its lack of dispersion — a coincidence, not a virtue. Modern vdW-corrected functionals (D3, vdW-DF, rVV10) give the right answer; the gold-standard RPA and DMC confirm them.

### Q&A: when is a PBE gap "good enough"?

*Q: My collaborator runs PBE on a series of perovskites and orders them by band gap. The absolute values are wrong, but the trend looks sensible. Can I trust the ordering?*

A: Often, but with caveats. PBE's systematic underestimate is *roughly constant* (often $\sim 1\;\text{eV}$ across a chemical family) so the *relative* ordering of gaps within a family is usually preserved. But two cautions: (i) when the underestimate flips the qualitative result (PBE metal where experiment is insulator, or PBE direct gap where experiment is indirect), the ordering breaks. (ii) For materials at the metal–insulator transition, the absolute error can flip across the threshold. Best practice: confirm a few key compounds with HSE06 to anchor the trend.

*Q: I see papers reporting "PBE band gap of 3 eV". Is this OK?*

A: It depends on what is claimed. Stating "the PBE Kohn–Sham gap is X eV" is fine and reproducible. Claiming "the band gap is X eV" without specifying the method is sloppy. The fundamental gap, the optical gap, and the KS gap are three different numbers in general; conflating them is one of the most common mistakes in DFT papers.

## 5.6.3 Strongly correlated electrons

**The symptom.** Apply PBE to FeO, CoO, NiO. PBE predicts all three to be metals. Experimentally, all three are antiferromagnetic insulators with gaps of 2–4 eV. Apply PBE to cerium oxide: the famous Ce $4f$ electrons come out delocalised, when in CeO$_2$ they are localised on cerium sites.

**The physics.** In these systems, the dominant energy scale is the on-site Coulomb repulsion $U$ between electrons in the same localised orbital (typically a $3d$ or $4f$ shell). When $U$ exceeds the hopping integral $t$, electrons localise on individual atoms and the system is a Mott insulator. The KS density of a Mott insulator is not the density of any non-interacting system in any reasonable potential: the single-Slater-determinant ansatz of KS theory is not a good starting point.

**What to do.**

- **DFT+U** (Anisimov–Liechtenstein–Zaanen): add a Hubbard-$U$ correction term to the energy functional, penalising fractional occupation of the localised shell. Choice of $U$ is empirical (3–8 eV typical), or computable via linear response. Cheap; often dramatically improves gaps and magnetic order in transition metal oxides.
- **DMFT** (dynamical mean-field theory): treat the local correlations on the correlated site exactly using an impurity solver (e.g., continuous-time quantum Monte Carlo), embedded in a DFT bath. The state of the art for strongly correlated materials. Cost: orders of magnitude beyond DFT.
- **Multireference quantum chemistry** (CASSCF, CASPT2, NEVPT2): for small clusters where the active space is manageable. Cost: prohibitive beyond ~20 active orbitals.
- **Hybrid functionals** sometimes fix Mott gaps via the exact-exchange fraction, but the result is functional- and parameter-dependent.

Strong correlation is the area where DFT is most likely to be qualitatively wrong, and where one most needs a higher-level method. The 2010s and 2020s have seen rapid development of DFT+DMFT codes (TRIQS, EDMFT) that automate the process.

!!! example "Worked example: NiO band gap and magnetism"
    NiO is a classic late transition-metal antiferromagnetic insulator. Experiment: AFM-II ground state, gap $\sim 4.0\;\text{eV}$, local moment $\sim 1.9\;\mu_B$ per Ni. DFT predictions:
    
    | Method | Gap (eV) | Moment ($\mu_B$) | Ground state |
    |---|---|---|---|
    | PBE (NM) | 0 | 0 | metal (wrong) |
    | PBE+SP | $\sim 0.5$ | $1.4$ | AFM metal (wrong) |
    | PBE+U ($U=6$) | $\sim 3.2$ | $1.7$ | AFM insulator |
    | HSE06 | $\sim 4.1$ | $1.8$ | AFM insulator |
    | DFT+DMFT | $\sim 4.0$ | $1.9$ | AFM insulator |
    
    Without any Hubbard correction or hybrid mixing, PBE predicts NiO to be a *metal* — a qualitative failure. Adding a Hubbard $U$ on the Ni $3d$ states penalises double occupation and opens the gap; HSE06 achieves the same effect through its 25% exact exchange. DMFT captures the full local correlation physics including spectral weight transfer to the upper Hubbard band, at the cost of an impurity solver.

## 5.6.4 Self-interaction error and charge transfer

We met self-interaction in §5.4: approximate exchange-correlation functionals do not cancel the spurious self-Hartree term, with the result that electrons artificially delocalise. Two consequences are worth singling out.

**Fractional charges in dissociation.** Take H$_2^{+}$, one electron, two protons. Stretch the bond to infinity. The correct answer is one electron localised on *one* of the protons (the other proton is a bare H$^{+}$). PBE instead delocalises the electron equally over both protons, giving a fractionally charged H$^{0.5+}$ — H$^{0.5+}$ configuration at infinite separation. The total energy is too low by tens of kcal/mol.

!!! example "Worked example: H$_2^{+}$ at infinite separation"
    For H$_2^{+}$ with one electron, the exact dissociation energy is $E(\mathrm H_2^{+}\to \mathrm H + \mathrm H^{+}) = -0.5\;\text{Ha}$ (the bound 1s electron on the remaining H atom). DFT calculations at $R = 10\;a_0$ (effectively infinite separation):
    
    | Method | $E$ (Ha) | Charge on each H |
    |---|---|---|
    | Exact / FCI | $-0.5000$ | 1.0 / 0.0 (symmetry-broken) |
    | Hartree–Fock | $-0.5000$ | 0.5 / 0.5 (incorrect delocalisation, but exchange exact) |
    | LDA | $-0.4523$ | 0.5 / 0.5 |
    | PBE | $-0.4567$ | 0.5 / 0.5 |
    | HSE06 (25% HF) | $-0.4823$ | 0.5 / 0.5 |
    | LC-$\omega$PBE (100% HF at LR) | $-0.4998$ | symmetry-broken |
    
    All semi-local functionals delocalise the electron equally over both protons because of self-interaction error — they prefer to spread the density to lower the (artificially included) self-Hartree. Only functionals with 100% exact exchange at long range cure this, by exactly cancelling the spurious self-repulsion at any separation.

**Charge transfer excitations.** Time-dependent DFT (TDDFT) with semi-local functionals notoriously fails for excited states involving long-range charge transfer (e.g., between a donor and an acceptor in a complex). The TDDFT excitation energy collapses to nearly the KS HOMO–LUMO gap — far below the true excitation energy, which should include the Coulomb attraction $-1/R$ of the resulting electron–hole pair.

**What to do.**

- **Range-separated hybrids** (CAM-B3LYP, $\omega$B97X) include 100% exact exchange at long range and cure long-range SIE.
- **Self-interaction correction (SIC)** functionals explicitly subtract the orbital self-interaction (Perdew–Zunger SIC). Computationally awkward (orbital-dependent potentials, non-Janak occupations) but available.
- **Constrained DFT (cDFT)** lets you fix charge configurations by hand for specific applications (electron transfer rate calculations).
- For chemical-accuracy needs: post-Hartree–Fock or coupled cluster.

## 5.6.5 Excited states

Kohn–Sham DFT is, by construction, a **ground-state theory**. The Hohenberg–Kohn theorems (§5.2) prove that the ground-state density determines everything; they say nothing about excited states. The KS eigenvalues are mathematical objects (§5.3), not excitation energies.

This is a serious limitation. Many of the most interesting properties of materials — optical absorption, fluorescence, photochemistry, photovoltaics — are excited-state phenomena.

### Time-dependent DFT (TD-DFT)

The cleanest extension is **time-dependent DFT**, based on the Runge–Gross theorem (1984): for fixed initial state, the time-dependent density $n(\mathbf r,t)$ determines the time-dependent external potential $v_\mathrm{ext}(\mathbf r,t)$ up to a purely time-dependent constant. This justifies a time-dependent KS scheme,

$$
i\hbar\frac{\partial\phi_i(\mathbf r,t)}{\partial t} = \Big[-\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}[n](\mathbf r,t)\Big]\phi_i(\mathbf r,t),
$$

with $v_\mathrm{KS}$ a time-dependent functional of the time-dependent density. Linearising around the ground state and Fourier transforming gives **linear-response TD-DFT**, the standard method for computing vertical excitation energies in molecules.

TD-DFT has its own pathologies. The required exchange-correlation kernel $f_{xc}(\mathbf r,\mathbf r',\omega) = \delta v_{xc}/\delta n$ is approximated in the **adiabatic** limit (frequency-independent kernel, taken from the ground-state functional). Adiabatic TD-DFT systematically fails for:

- Charge-transfer excitations (above).
- Rydberg states (need long-range exact exchange).
- Double excitations (require frequency dependence of the kernel; adiabatic TD-DFT misses them entirely).
- Conical intersections (the topology is wrong).

For routine vertical singlet excitations in organic molecules, TD-DFT with a hybrid functional (CAM-B3LYP, $\omega$B97X-D) often gives errors of 0.2–0.4 eV. For anything outside this comfort zone, more sophisticated methods are needed.

### GW and Bethe–Salpeter

For solids, the modern gold standard for optical absorption is **GW+BSE**. First compute quasiparticle energies via GW (corrected band structure). Then solve the **Bethe–Salpeter equation** for the electron–hole interaction kernel, giving exciton binding energies and optical spectra including excitonic effects. Cost: $\mathcal O(N^{4})$–$\mathcal O(N^{6})$. Accuracy for absorption peaks: 0.1–0.2 eV.

### Quantum chemistry methods for excited states

- **EOM-CCSD** (equation-of-motion coupled cluster): chemical accuracy for low-lying excitations in molecules. Scales $\mathcal O(N^{6})$.
- **CASPT2 / NEVPT2** (complete active space perturbation theory): handles multireference systems including double excitations and conical intersections. Bespoke selection of active space required.
- **ADC(2), ADC(3)**: algebraic diagrammatic construction, intermediate accuracy and cost.

## 5.6.6 When to reach for higher methods

A practical decision tree. If your DFT calculation is suspicious or if your science demands more than DFT can deliver, consider:

| Failure mode | First-line remedy | Higher-level method |
|---|---|---|
| Band gap | HSE06 hybrid | GW ($G_0W_0$, scGW) |
| Optical absorption | TD-DFT with hybrid | GW+BSE; EOM-CCSD |
| van der Waals binding | PBE+D3, optB88-vdW | RPA; QMC |
| Reaction barriers | M06-2X, $\omega$B97X-D | CCSD(T) (gold standard) |
| Mott insulator | DFT+U; HSE06 | DFT+DMFT |
| Multireference (transition states, biradicals) | broken-symmetry DFT | CASSCF/CASPT2; MRCI |
| Charge transfer excitations | range-separated hybrid | EOM-CCSD; ADC(2) |
| Photochemistry / conical intersections | not safely DFT | CASSCF; MS-CASPT2 |
| Quantitative formation energies of small molecules | hybrid + D3 | composite methods (G4, W1); CCSD(T) |
| Thermochemistry of large organics | hybrid + D3 | DLPNO-CCSD(T) (local correlation) |

A few rules:

1. **CCSD(T)** — coupled cluster with singles, doubles, and perturbative triples — is the "gold standard" of quantum chemistry. For systems where it is tractable (~30 atoms), its accuracy is essentially benchmark-quality. Beyond ~50 atoms, even local approximations to CCSD(T) become expensive.
2. **Quantum Monte Carlo (QMC)** — variational and diffusion Monte Carlo — gives benchmark accuracy with $\mathcal O(N^{3})$ scaling but a large prefactor. For carefully chosen problems (cohesive energies of solids, vdW binding) it is unrivalled.
3. **GW+BSE** is the standard for optical properties of solids; pair with a good DFT starting point.
4. **DMFT** is the standard for strongly correlated materials; needs an impurity solver and considerable expertise.
5. **Machine learning interatomic potentials** (Chapter 9) cannot save you from a bad reference: a model trained on DFT data inherits DFT's errors. If DFT is wrong for your system, an ML potential trained on DFT will be wrong in the same way.

### Escalation table: cost ratios for higher-level methods

For a system of $N$ atoms (and a fixed basis), approximate cost scalings:

| Method | Scaling | Cost vs. PBE | When to use |
|---|---|---|---|
| LDA, PBE | $\mathcal O(N^{3})$ | 1 | High-throughput, screening |
| SCAN | $\mathcal O(N^{3})$ | 1.5 | Improved binding, no exact exchange |
| HSE06 | $\mathcal O(N^{3})$ (large prefactor) | 20 | Band gaps, defects in semiconductors |
| DFT+U | $\mathcal O(N^{3})$ | 1.2 | Transition metal oxides, $f$-electrons |
| $G_0W_0$ | $\mathcal O(N^{4})$ | 100 | Quasiparticle band structures |
| Self-consistent GW | $\mathcal O(N^{4})$ | 500 | When $G_0W_0$ starting-point sensitivity matters |
| GW+BSE | $\mathcal O(N^{4})$–$\mathcal O(N^{6})$ | 200–2000 | Optical absorption with excitons |
| RPA correlation | $\mathcal O(N^{5})$ | 200 | Total energies of vdW systems |
| MP2 | $\mathcal O(N^{5})$ | 50 | Reasonable cost weakly-correlated chemistry |
| CCSD | $\mathcal O(N^{6})$ | 500 | Single-reference systems, small molecules |
| CCSD(T) | $\mathcal O(N^{7})$ | 5000 | "Gold standard" for small molecules |
| DLPNO-CCSD(T) | $\mathcal O(N)$ (with constant) | 100–1000 | Local correlation for larger systems |
| FCI / DMRG | exponential | $10^{6}$+ | Tiny systems, benchmarks only |
| DMFT (CT-QMC impurity) | depends | 100–10000 | Strong correlation |
| QMC (DMC) | $\mathcal O(N^{3})$ large prefactor | $10^{4}$–$10^{6}$ | Benchmarks; cohesive energies |

A few orientation points: $G_0W_0$ on a PBE starting point is the most-used post-DFT method for quasiparticle gaps and adds about two orders of magnitude to the cost. CCSD(T) is the gold standard for molecules up to ~50 atoms; beyond that, local approximations (DLPNO) are essential. QMC is used as a benchmark — accurate to ~1 meV/atom for cohesive energies — but is rarely the first-line method because of its prefactor.

## 5.6.7 An honest assessment

DFT is, for an extraordinary range of systems, the right tool: fast enough for high-throughput screening, accurate enough for materials prediction, and based on a rigorous theoretical foundation. It is the engine behind essentially every materials database, every ML-potential training set, every large-scale electronic-structure calculation done in industry. None of that is going to change soon.

But it is not magic. There is no single functional that is best for everything; there are systems where any practical functional is qualitatively wrong; there are properties (excited states, fundamental gaps) where the framework itself is not designed for the question being asked. A good computational materials scientist knows:

- **What** their functional gets right and wrong for the class of system they study.
- **Why** — at the level of physics, not just empirics.
- **When** to escalate to a higher-level method.

Chapter 6 turns to the practical business of running DFT calculations: plane waves, pseudopotentials, $k$-point sampling, convergence testing, and the choice of code. Chapter 7 covers the post-DFT methods touched on here — GW, BSE, DMFT — in more depth. The Hohenberg–Kohn–Kohn–Sham theorem is, in the end, an existence proof; the practical art begins with knowing how to use it well, and when to put it down.

!!! question "Check yourself"
    1. Name **three** of DFT's systematic failure modes and, for each, the underlying cause in one phrase.
    2. A colleague says "my silicon gap is too small, I'll just tighten the convergence and add more $k$-points." Will that fix it? Why or why not?
    3. Why are excited states fundamentally hard for Kohn–Sham DFT?

    ??? success "Answers"
        1. Any three of, for example: **band-gap underestimation** (missing derivative discontinuity $\Delta_{xc}$ plus self-interaction error); **missing van der Waals binding** (a local functional cannot capture the long-range $-C_6/R^{6}$ dispersion correlation); **strong-correlation failures** in Mott insulators (the on-site repulsion $U$ dominates, so the single-determinant KS picture breaks down); **self-interaction / delocalisation error** (an electron feels its own charge, spreading density and giving fractional charges at dissociation); **excited states** (the framework is ground-state by construction).
        2. **No.** The gap error is intrinsic, not numerical. Semi-local functionals set $\Delta_{xc}=0$ and the quantity read off is a Kohn–Sham eigenvalue difference for the auxiliary non-interacting system. No amount of tighter SCF convergence or denser $k$-mesh changes that; you need a better functional (HSE06) or a quasiparticle method (GW).
        3. Because the Hohenberg–Kohn theorems (and the Kohn–Sham construction) establish that the *ground-state* density determines the system; they say nothing about excited states, and the KS eigenvalues are auxiliary objects (Section 5.3), not excitation energies. Accessing excited states requires an extended framework — TD-DFT for molecules, GW+BSE for solids, or wavefunction methods such as EOM-CCSD.

### Summary of §5.6 — what to remember in 3 months

- **Band gaps**: LDA/GGA underestimate by 30–100%, due to missing derivative discontinuity + SIE. Use HSE06 (cheap fix) or GW (expensive correct).
- **vdW**: semi-local functionals miss $-C_6/R^{6}$. Always include D3/D4 or use vdW-DF/SCAN+rVV10 when non-bonded fragments are present.
- **Mott insulators**: PBE predicts metal for NiO/CoO/FeO. Use DFT+U, hybrid, or DMFT.
- **Fractional charges**: stretched H$_2^{+}$ and similar; due to SIE. Use range-separated hybrids with 100% LR exact exchange.
- **Excited states**: KS-DFT is ground-state only. Use TD-DFT (molecules) or GW+BSE (solids).
- **Escalation methods**: HSE06 → GW → GW+BSE → CCSD(T) → DMFT/QMC. Each is roughly 10–100× more expensive than the previous.
- **The rule of thumb**: never trust a single DFT calculation; compare functionals, compare with experiment, escalate when stakes are high.

!!! note "Remark: ML potentials and DFT errors"
    Machine-learning interatomic potentials inherit the errors of the DFT functional they are trained on, *exactly*. If you train a GAP or MACE model on PBE forces, it will reproduce PBE-overbound vdW interactions, PBE band gaps, PBE bond lengths. The model is at best a surrogate for the functional; it cannot exceed the accuracy of the training labels. This is the *garbage-in-garbage-out* principle of ML, made specific. In Chapter 9 we shall see how to choose training functionals appropriate to downstream tasks.

!!! success "What to remember from Chapter 5"
    Stepping back over the whole chapter, the logical chain of density-functional theory is:

    - The **electron density** $n(\mathbf r)$ replaces the many-body wavefunction as the basic variable — three coordinates instead of $3N$.
    - The **Hohenberg–Kohn theorems** make this legitimate: the ground-state density determines the external potential, and hence everything, and a variational principle holds for the density.
    - The **Kohn–Sham construction** makes it computable, mapping the interacting problem onto a fictitious system of non-interacting orbitals that reproduce the true density.
    - The **exchange–correlation energy** $E_{xc}$ holds everything we do not know exactly and *must be approximated* (LDA, GGA, meta-GGA, hybrids) — this is where all the modelling choices, and most of the errors, live.
    - The **self-consistent field (SCF) loop** solves the resulting coupled equations, iterating density and potential to self-consistency.
    - The **Kohn–Sham eigenvalues** are auxiliary quantities — excellent for interpretation, but not in general physical excitation energies.
    - And you should **know the failure modes** — band gaps, van der Waals, strong correlation, self-interaction, excited states — so you can recognise when DFT is being quietly wrong.

    With the theory in hand, [Chapter 6](../ch06-running-dft/index.md) turns to running real DFT: plane waves, pseudopotentials, $k$-point sampling and convergence testing.
