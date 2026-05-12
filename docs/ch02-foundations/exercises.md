# Chapter 2 — Exercises

Six problems, weighted toward conceptual fluency. Difficulty: <span class="diff-easy">★</span> easy, <span class="diff-medium">★★</span> medium, <span class="diff-hard">★★★</span> hard. Solutions are collapsed; attempt the problem before expanding.

---

## Exercise 2.1 — Pick a rung <span class="diff-easy">★</span>

For each of the following research questions, identify the most appropriate rung of the scale ladder (electronic, atomistic, mesoscale, continuum) and a representative method or code.

(a) What is the equilibrium lattice constant of FCC nickel at $T = 0$?

(b) How does the grain-size distribution of an aluminium alloy evolve over 10 hours of annealing at 400 °C?

(c) What is the diffusion coefficient of liquid water at 300 K?

(d) Will a 5-metre titanium beam buckle under a load of 100 kN?

(e) What is the activation energy for vacancy migration in copper?

(f) How does a polymer melt flow through a 1 mm channel?

??? success "Solution"
    (a) **Electronic.** A simple DFT relaxation. Code: Quantum ESPRESSO, VASP, GPAW. Expected error: $\sim 1\%$.

    (b) **Mesoscale.** Phase-field simulation of grain growth, or kinetic Monte Carlo. Code: MOOSE, PRISMS-PF, SPPARKS. The mesoscale is set by the requirement to resolve grain shapes (µm scale) over hours.

    (c) **Atomistic.** Molecular dynamics with an appropriate water model (TIP4P, MB-pol, or a foundation MLIP like MACE-MP-0). Code: GROMACS, LAMMPS, or ASE. A nanosecond trajectory yields a converged diffusion coefficient.

    (d) **Continuum.** Finite-element elasticity calculation. Code: ABAQUS, FEniCS. Atomic-scale physics is irrelevant; the input is the elastic modulus, which could itself come from DFT.

    (e) **Electronic.** Nudged elastic band (NEB) calculation in DFT, computing the energy along the vacancy hop pathway. Code: VASP or QE plus the NEB extension. Output feeds upward into KMC simulations of vacancy diffusion.

    (f) **Continuum.** Computational fluid dynamics, perhaps coupled to a viscoelastic constitutive model whose parameters come from coarse-grained MD. Code: OpenFOAM, COMSOL.

    The general pattern: identify the *resolution* of the question. The smallest length and time scales mentioned in the question are typically the dominant constraint.

---

## Exercise 2.2 — Capabilities audit <span class="diff-easy">★</span>

Classify each statement below as **trustworthy** (a routine 2026 capability), **trustworthy with caveats** (works, but requires care), or **not yet** (an open problem). Briefly justify each.

(a) Predict the band gap of a new lead halide perovskite to within $0.3$ eV.

(b) Compute the elastic constants of a known FCC metal to within 10%.

(c) Predict the superconducting transition temperature of a new cuprate.

(d) Screen 50000 candidate alloys by formation energy to identify the most stable.

(e) Predict the melting point of UO$_2$ to within 50 K.

(f) Predict from first principles that a new compound, never synthesised, will form by mixing two precursors at 1000 °C for two hours.

??? success "Solution"
    (a) **Trustworthy with caveats.** Hybrid functionals (HSE06) or GW reach $\sim 0.3$ eV accuracy for typical semiconductors. For lead halide perovskites, spin–orbit coupling is essential and dramatically changes the gap; standard GGA gives the wrong answer for the wrong reason. With SOC + HSE the prediction is in range. Expect $\sim 0.3$ eV error rather than $0.1$ eV.

    (b) **Trustworthy.** Elastic constants of metals from DFT typically agree with experiment to 5–10%. Standard workflow with energy-versus-strain fits.

    (c) **Not yet.** Cuprates are strongly correlated and unconventional. No reliable workflow exists.

    (d) **Trustworthy with caveats.** Routine for ranking. The caveat is that *stable on the convex hull of structures we considered* is not the same as *stable in nature*. If you missed a competing phase, ranking is wrong.

    (e) **Trustworthy with caveats.** Possible with a good potential (modern MLIP, or specialised classical) and two-phase coexistence or interface-pinning methods. UO$_2$ is specifically challenging because it is a Mott insulator and DFT+U is required for accurate energetics. With those caveats, 50 K accuracy is plausible.

    (f) **Not yet.** Synthesis prediction — including precursor choice, temperature schedule, and yield — remains beyond routine simulation. Active-learning and large-language-model approaches are nascent (see Chapter 12), but no robust workflow exists.

---

## Exercise 2.3 — Reading a phase diagram <span class="diff-medium">★★</span>

Consider a hypothetical binary system A–B with the following features (sketch it for yourself):

- A melts at 1200 K; B melts at 800 K.
- A eutectic at $x_B = 0.6$, $T_\mathrm{eut} = 600$ K.
- Maximum solubility of B in solid A (the $\alpha$ phase) is $0.10$ at $T_\mathrm{eut}$, falling to $0.02$ at room temperature.
- Maximum solubility of A in solid B (the $\beta$ phase) is $0.08$ at $T_\mathrm{eut}$, negligible at room temperature.

Answer the following.

(a) Describe the phases present, at equilibrium, in an alloy of $x_B = 0.4$ at $T = 700$ K. Apply the lever rule to estimate the phase fractions.

(b) Describe what happens, qualitatively, as the same alloy is cooled slowly from 1300 K to 300 K. Identify each phase boundary it crosses.

(c) Why does the solubility of B in $\alpha$ decrease with temperature? What microstructural feature would you expect to see in a room-temperature sample of an alloy with $x_B = 0.08$ that has been cooled slowly from above the eutectic?

??? success "Solution"
    (a) At $T = 700$ K and $x_B = 0.4$, we are above the eutectic temperature, so we are in the two-phase region between $\alpha$ and liquid. The $\alpha$ phase boundary at 700 K lies somewhere near $x_B \approx 0.08$ (the solubility limit decreases roughly linearly from 0.10 at 600 K toward 0.02 at room temperature, but is still close to 0.10 just above the eutectic). The liquidus at 700 K can be read off the diagram; for a steep liquidus typical of dilute eutectics, it is around $x_B \approx 0.5$.

    Lever rule: fraction of liquid $= (0.40 - 0.08) / (0.50 - 0.08) \approx 0.76$. Fraction of $\alpha$ solid is $\approx 0.24$. The $\alpha$ has composition $x_B \approx 0.08$; the liquid has composition $x_B \approx 0.50$.

    (b) Cooling from 1300 K: completely liquid. Around 900–1000 K (somewhere on the liquidus), $\alpha$ phase begins to crystallise — we enter the $\alpha + L$ two-phase region. As $T$ falls, more $\alpha$ forms; the liquid composition runs along the liquidus toward the eutectic at $x_B = 0.6$. At $T = T_\mathrm{eut} = 600$ K, the remaining liquid solidifies as a eutectic mixture of $\alpha$ ($x_B = 0.10$) and $\beta$ ($x_B = 0.92$). Below 600 K we are in the $\alpha + \beta$ two-phase region. As $T$ falls further, the solubility lines move inward: $\alpha$ rejects B (precipitating fine $\beta$) and $\beta$ rejects A.

    (c) Solubility decreases with temperature because the solid solution has an unfavourable enthalpy of mixing relative to phase separation; entropy stabilises the solution at high $T$ but loses out at low $T$. The room-temperature microstructure of an alloy starting at $x_B = 0.08$, cooled slowly, would show a matrix of $\alpha$ (with $x_B \approx 0.02$) containing precipitates of $\beta$. This is precipitation hardening — the same mechanism that gives age-hardened aluminium alloys their strength.

---

## Exercise 2.4 — Reading a band structure <span class="diff-medium">★★</span>

You are shown a calculated band structure for an unknown crystalline material. You observe:

- the valence-band maximum lies at $\Gamma$, with energy 0;
- the conduction-band minimum lies at $\Gamma$, with energy $+2.4$ eV;
- the valence band near $\Gamma$ has two branches with curvatures corresponding to effective masses of $0.5 \, m_e$ (light hole) and $1.2 \, m_e$ (heavy hole);
- the conduction band near $\Gamma$ has curvature $0.2 \, m_e$;
- the DOS is approximately zero in a window of $2.4$ eV around $E_\mathrm{F} = 0$.

(a) Classify the material: metal, semiconductor, or insulator? Direct or indirect gap?

(b) Will this material absorb visible light efficiently?

(c) Would you expect this material to be a good electron or hole conductor (which carrier has higher mobility)?

(d) Suggest a class of materials this might belong to.

??? success "Solution"
    (a) **Semiconductor with a direct gap.** A 2.4-eV gap is too large to be a metal but small enough that, with appropriate doping, the material conducts; the conduction-band minimum and valence-band maximum sit at the same $\mathbf{k}$, so the gap is direct.

    (b) Visible light spans roughly 1.6–3.1 eV. A 2.4-eV direct gap absorbs visible light from $\sim 520$ nm and shorter (green, blue, violet, UV). The material will appear yellow-orange in transmission. Crucially, the gap being *direct* means absorption is dipole-allowed and strong; this material would be a good photovoltaic absorber and a good light emitter.

    (c) **Electrons.** The conduction-band effective mass ($0.2 \, m_e$) is much smaller than the heavy-hole mass ($1.2 \, m_e$). Mobility scales inversely with effective mass at fixed scattering rate, so electrons will be more mobile.

    (d) The combination — direct gap of $\sim 2.4$ eV, light conduction electrons, heavy and light holes at $\Gamma$ — is characteristic of **direct-gap III–V or II–VI semiconductors** (GaP at 2.26 eV indirect, ZnSe at 2.7 eV direct, CdS at 2.4 eV direct). Among these, CdS or a similar wide-gap II–VI compound is the closest fit.

---

## Exercise 2.5 — Interpreting an RDF <span class="diff-medium">★★</span>

You receive an RDF $g(r)$ from a collaborator's MD simulation of a metal at unknown conditions. The features:

- sharp first peak at $r_1 = 2.55$ Å with peak height $g(r_1) \approx 3.5$;
- a clear trough near $r = 3.0$ Å, dropping to $g(r) \approx 0.4$;
- a second peak at $r_2 \approx 4.4$ Å, broader, with $g(r_2) \approx 1.4$;
- a third weaker peak near $r_3 \approx 5.6$ Å;
- $g(r) \to 1$ smoothly by $r \approx 8$ Å.

(a) Is the system crystalline or liquid? Justify.

(b) Estimate the nearest-neighbour coordination number, assuming the density is $\rho = 0.060$ Å$^{-3}$. Use $N_\mathrm{nn} \approx 4\pi \rho \int_0^{r_\mathrm{min}} r^2 g(r) \, \mathrm{d}r$ where $r_\mathrm{min}$ is the first minimum. Make a reasonable approximation.

(c) Guess the structure that the liquid is likely sitting near.

??? success "Solution"
    (a) **Liquid (or amorphous).** The signature is the broad second peak with washed-out structure, followed by smooth decay to unity by $\sim 8$ Å. A crystalline RDF would have well-separated sharp peaks at every coordination-shell distance, with $g(r) \approx 0$ in between.

    (b) Approximate the first peak as a triangle (or just use the fact that the integral up to the first minimum picks up roughly the area under the first peak):

    $$
    N_\mathrm{nn} = 4\pi \rho \int_0^{r_\mathrm{min}} r^2 g(r) \, \mathrm{d}r.
    $$

    Replace $g(r) r^2$ by its peak value $g(r_1) r_1^2 \approx 3.5 \times (2.55)^2 \approx 22.8$ Å$^2$, and let the effective width be the full width of the peak, perhaps $\Delta r \approx 0.7$ Å:

    $$
    N_\mathrm{nn} \approx 4\pi \times 0.060 \times 22.8 \times 0.7 \approx 12.
    $$

    A coordination number of $\sim 12$ is consistent with a close-packed liquid metal (the FCC and HCP coordination number at $T = 0$ is 12; liquids near close packing keep $\sim 11$–12).

    (c) The system is most likely a **simple metallic liquid near close packing** — copper, aluminium, or another FCC metal just above its melting point. The peak distance of 2.55 Å is in the right range for liquid copper (whose nearest-neighbour distance in FCC is 2.56 Å). Confirming this would require the temperature and the density at that temperature, but the RDF alone is suggestive.

---

## Exercise 2.6 — Designing a workflow <span class="diff-hard">★★★</span>

You are tasked with finding a new oxide catalyst for the oxygen evolution reaction (OER) at the anode of a water electrolyser. The catalyst should be more active than IrO$_2$ (the current state of the art) and contain only earth-abundant elements.

Design a multi-stage computational workflow that uses methods from each rung of the scale ladder where appropriate. For each stage, specify (i) the input, (ii) the method and code, (iii) the output, (iv) the role of the stage in the overall search. Justify your design choices. Address explicitly: where DFT is essential, where an MLIP might be used to accelerate sampling, where you would consult a database, and where you would invoke a generative or screening model. Acknowledge realistic limitations.

??? success "Solution"
    A defensible workflow has roughly the following stages.

    **Stage 1 — Composition screening (database query).**
    *Input:* periodic table restricted to earth-abundant elements.
    *Method:* query Materials Project and/or Alexandria for known oxides containing transition metals (Fe, Co, Ni, Mn, Cu, ...) and one or more secondary elements.
    *Output:* a list of $\sim 10^4$ candidate compositions with known DFT-relaxed structures and formation energies.
    *Role:* eliminate compositions that are not thermodynamically accessible.

    **Stage 2 — Property pre-screening (ML descriptor).**
    *Input:* the $\sim 10^4$ candidates.
    *Method:* a trained property-prediction model (graph neural network on Materials Project data, or a foundation model) that predicts a known OER-relevant descriptor — for instance, the OH adsorption energy on a typical surface, or simpler bulk descriptors known to correlate with activity (band-centre proxies).
    *Output:* a ranked list, top 100 candidates.
    *Role:* reduce the candidate list cheaply by orders of magnitude.

    **Stage 3 — Surface enumeration (atomistic, ASE/pymatgen).**
    *Input:* top 100 bulk candidates.
    *Method:* for each, enumerate plausible low-index surface terminations using pymatgen's slab-builder.
    *Output:* $\sim 1000$ surface slabs.
    *Role:* the catalysis happens on a surface, not in the bulk; the surface choice matters.

    **Stage 4 — Adsorption energy calculation (DFT).**
    *Input:* the $\sim 1000$ surfaces.
    *Method:* DFT with a hybrid functional (or PBE+U where needed). Compute the binding energies of OER intermediates — *OH, *O, *OOH — on relevant surface sites. Use the standard OER scaling-relation analysis and the Sabatier overpotential.
    *Output:* predicted theoretical overpotential $\eta$ for each surface.
    *Code:* VASP or QE, automated with atomate2.
    *Role:* this is where the real chemistry enters; it is also the most expensive step. Aim to spend $\sim 80\%$ of compute here.

    **Stage 5 — Stability under operating conditions (DFT thermodynamics).**
    *Input:* the top 10 candidates by predicted overpotential.
    *Method:* Pourbaix-diagram analysis using DFT bulk and surface energies plus aqueous-ion experimental data. Determine whether the surface dissolves, oxidises further, or remains stable at the relevant pH and potential.
    *Output:* stability map and a shortlist of candidates that are both active and stable.
    *Role:* a catalyst that dissolves in the electrolyser is useless regardless of its activity.

    **Stage 6 — Kinetic refinement (DFT + MLIP).**
    *Input:* the shortlist.
    *Method:* NEB calculations for proton-coupled electron-transfer steps; if necessary, train an MLIP on the DFT data and run reactive MD to sample disordered surfaces.
    *Output:* refined activation energies and a final ranking.
    *Role:* close the gap between thermodynamic descriptors and actual reaction kinetics.

    **Stage 7 — Experimental hand-off.**
    *Output:* a list of 3–5 compositions with synthesis-relevant data (target structure, expected stability, suggested precursors based on Materials Project precursor analysis or related tools).
    *Role:* commit to experimental tests. No simulation pipeline replaces this step.

    **Realistic limitations to acknowledge.**

    - DFT errors on adsorption energies are typically 0.1–0.2 eV. Differences within the OER scaling relations smaller than this are not significant.
    - The Sabatier framework collapses a real catalytic cycle to two or three descriptors; real catalysis often violates the assumptions.
    - Property predictors at Stage 2 work well *in distribution*. If your target chemistry is far from the training set, expect surprises.
    - Synthesis-route prediction at Stage 7 is currently weak. Expect a high rate of *predicted but not made* failures.
    - Real OER catalysts in industrial use are often nanostructured, hydrated, and amorphous; an idealised crystalline-slab calculation may miss what the real catalyst is doing.

    A defensible workflow does not promise to *discover* a new catalyst. It promises to *narrow the experimental search space* by one to two orders of magnitude. That is what computational catalysis has actually delivered in 2026, and a sober proposal reflects this.

---

When you can do all six of these confidently, you are ready for Chapter 3, where we descend into the atomic and crystalline detail that everything in Chapter 2 quietly assumed.
