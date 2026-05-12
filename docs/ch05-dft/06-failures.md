# 5.6 Where DFT Fails

DFT is the workhorse of computational materials science because it is, on the whole, *good enough* — it gets bond lengths within a few per cent, cohesive energies within tens of kJ/mol, vibrational frequencies within 10%, and the qualitative ground-state physics right for the vast majority of systems. But "the vast majority" leaves a long tail of systems and properties where standard Kohn–Sham DFT, in any practical approximation, is simply wrong. Some of these failures are technical (chosen functional too crude); others are fundamental (no semi-local functional can fix them).

This section gives an honest tour. Knowing where DFT breaks is the difference between a trustworthy calculation and a published mistake. For each failure mode we identify the symptom, the underlying physics, and the higher-level methods one reaches for instead.

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

!!! warning "Do not over-interpret PBE band gaps"
    A PBE band gap is not "the band gap". It is the Kohn–Sham gap of a particular approximate functional. For predictions of optical or transport gaps, use HSE06 or GW; for ordering of mid-band features, PBE often suffices. Always state the functional alongside the gap.

## 5.6.2 Van der Waals dispersion

**The symptom.** Stack two graphene sheets at 3.35 Å, the experimental interlayer spacing of graphite. Compute the binding energy with PBE: about 1 meV/atom, essentially zero. The experimental value is around 50 meV/atom. PBE predicts graphite to be barely bound, when it is a robust layered solid.

Or: try a benzene dimer. PBE gives no binding. Or: rare-gas dimers — argon, krypton, xenon — all repulsive with semi-local DFT.

**The physics.** London dispersion forces arise from instantaneous quantum fluctuations of the charge density on one fragment polarising another. The induced dipole pair gives the famous $-C_6/R^{6}$ attraction at large separation. This is a long-range correlation effect: the densities of the two fragments do not overlap, so any *local* functional sees nothing happening between them. Semi-local exchange-correlation, by construction, cannot reproduce $-C_6/R^{6}$.

**What to do.**

- **DFT-D3 / D4** (Grimme): add an empirical pairwise correction. Cheap and effective for most systems.
- **vdW-DF / vdW-DF2 / rVV10**: a non-local correlation kernel built into the functional. Computationally tractable via FFT; available in most plane-wave codes.
- **Tkatchenko–Scheffler / MBD**: density-dependent dispersion coefficients, including many-body screening effects (MBD: many-body dispersion). Best-in-class for systems where polarisability matters.

For materials with non-bonded fragments — molecular crystals, layered materials, surface adsorption, polymers, biomolecules — *not* including a vdW correction in DFT is a methodological error. Modern best practice always includes one.

## 5.6.3 Strongly correlated electrons

**The symptom.** Apply PBE to FeO, CoO, NiO. PBE predicts all three to be metals. Experimentally, all three are antiferromagnetic insulators with gaps of 2–4 eV. Apply PBE to cerium oxide: the famous Ce $4f$ electrons come out delocalised, when in CeO$_2$ they are localised on cerium sites.

**The physics.** In these systems, the dominant energy scale is the on-site Coulomb repulsion $U$ between electrons in the same localised orbital (typically a $3d$ or $4f$ shell). When $U$ exceeds the hopping integral $t$, electrons localise on individual atoms and the system is a Mott insulator. The KS density of a Mott insulator is not the density of any non-interacting system in any reasonable potential: the single-Slater-determinant ansatz of KS theory is not a good starting point.

**What to do.**

- **DFT+U** (Anisimov–Liechtenstein–Zaanen): add a Hubbard-$U$ correction term to the energy functional, penalising fractional occupation of the localised shell. Choice of $U$ is empirical (3–8 eV typical), or computable via linear response. Cheap; often dramatically improves gaps and magnetic order in transition metal oxides.
- **DMFT** (dynamical mean-field theory): treat the local correlations on the correlated site exactly using an impurity solver (e.g., continuous-time quantum Monte Carlo), embedded in a DFT bath. The state of the art for strongly correlated materials. Cost: orders of magnitude beyond DFT.
- **Multireference quantum chemistry** (CASSCF, CASPT2, NEVPT2): for small clusters where the active space is manageable. Cost: prohibitive beyond ~20 active orbitals.
- **Hybrid functionals** sometimes fix Mott gaps via the exact-exchange fraction, but the result is functional- and parameter-dependent.

Strong correlation is the area where DFT is most likely to be qualitatively wrong, and where one most needs a higher-level method. The 2010s and 2020s have seen rapid development of DFT+DMFT codes (TRIQS, EDMFT) that automate the process.

## 5.6.4 Self-interaction error and charge transfer

We met self-interaction in §5.4: approximate exchange-correlation functionals do not cancel the spurious self-Hartree term, with the result that electrons artificially delocalise. Two consequences are worth singling out.

**Fractional charges in dissociation.** Take H$_2^{+}$, one electron, two protons. Stretch the bond to infinity. The correct answer is one electron localised on *one* of the protons (the other proton is a bare H$^{+}$). PBE instead delocalises the electron equally over both protons, giving a fractionally charged H$^{0.5+}$ — H$^{0.5+}$ configuration at infinite separation. The total energy is too low by tens of kcal/mol.

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

## 5.6.7 An honest assessment

DFT is, for an extraordinary range of systems, the right tool: fast enough for high-throughput screening, accurate enough for materials prediction, and based on a rigorous theoretical foundation. It is the engine behind essentially every materials database, every ML-potential training set, every large-scale electronic-structure calculation done in industry. None of that is going to change soon.

But it is not magic. There is no single functional that is best for everything; there are systems where any practical functional is qualitatively wrong; there are properties (excited states, fundamental gaps) where the framework itself is not designed for the question being asked. A good computational materials scientist knows:

- **What** their functional gets right and wrong for the class of system they study.
- **Why** — at the level of physics, not just empirics.
- **When** to escalate to a higher-level method.

Chapter 6 turns to the practical business of running DFT calculations: plane waves, pseudopotentials, $k$-point sampling, convergence testing, and the choice of code. Chapter 7 covers the post-DFT methods touched on here — GW, BSE, DMFT — in more depth. The Hohenberg–Kohn–Kohn–Sham theorem is, in the end, an existence proof; the practical art begins with knowing how to use it well, and when to put it down.
