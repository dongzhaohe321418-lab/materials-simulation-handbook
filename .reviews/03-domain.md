# Peer Review 2 — Domain (Condensed Matter / DFT / MLIPs)

**Reviewer.** Senior condensed-matter theorist; daily DFT and MLIP user; author of recent review-level material on equivariant interatomic potentials.

**Files sampled** (~30 min sample, not a line-by-line audit).
- `docs/ch05-dft/02-hohenberg-kohn.md`
- `docs/ch05-dft/04-xc-functionals.md`
- `docs/ch09-mlip/05-equivariant.md`
- `docs/ch09-mlip/06-training-mace.md`
- `docs/ch12-foundation/02-mace-mp0.md`
- `docs/ch03b-solid-state/01-bloch-theorem.md`
- `docs/ch03b-solid-state/03-tight-binding.md`
- `docs/ch03b-solid-state/05-phonons.md`
- `docs/ch07-md/03-thermostats.md`
- `docs/appendix/E-bibliography.md`
- `docs/known-limitations.md`

## 1. Overall verdict

This is, on the technical merits I sampled, a *strong* manuscript — well-written, carefully scaffolded, and unusually candid about the failure modes of the methods it teaches. The Hohenberg–Kohn chapter (§5.2) and the equivariant-MLIP exposition (§9.5) are the best self-contained treatments of those topics I have read in a teaching text in the last five years. The thermostat chapter (§7.3) is technically meticulous (extended-system proof, Hoover transformation, BAOAB splitting, chain ergodicity — all in one place). The MACE training walkthrough (§9.6) is the kind of recipe a new graduate student can actually follow on day one, including a validation suite (parity / MD stability / RDF / phonons) that most "MLIP works for X" papers omit.

The book is not flawless. There are a small number of clear physical or compositional errors (one of which is an unedited "wait, let me redo" left in the published text — see §3 below) and the bibliography misses a handful of papers that any condensed-matter student would reach for. None of these are foundational; they should be straightforward to fix in a v1.1.

I would recommend acceptance of the technical content I sampled, conditional on the issues in §3 being addressed.

## 2. Scores (each 0–10)

| Dimension | Score | Comment |
|---|---|---|
| Physical correctness | 8 | Mostly correct. A few errors flagged below; no foundational mistakes I caught. |
| Modernity (2024–2026 coverage) | 9 | MACE-MP-0, OMat24, Orb-v2, eqV2-OMat, MatterGen, $\kappa_\mathrm{SRME}$ all present. Almost everything a 2026 reader needs. |
| Convention adherence | 8 | Atomic units used consistently in §5; eV / Å / fs in §7 and §9. Some Wigner D-matrix vs Clebsch–Gordan sign issues are not flagged. |
| Bibliography completeness | 6 | Solid core but five obvious omissions, listed below. e3nn and the Pozdnyakov–Ceriotti completeness paper are particularly conspicuous. |
| Honesty about limitations | 9 | `known-limitations.md` is unusually candid; §12.2 "When does the foundation fail?" is exemplary; §9.6.8 common-failure list is excellent. |
| Coverage balance (DFT vs MD vs ML) | 8 | Roughly proportional to the 2026 field. Surface-science / catalysis is light (one chapter would help) but not negligible. |
| Depth in core areas (HK, MACE, etc.) | 9 | HK chapter and MACE/NequIP chapter would not embarrass a graduate-level monograph. |

**Aggregate: 8.1 / 10.** Strong recommend; minor revisions.

## 3. Where the physics is taught well (three places)

1. **Hohenberg–Kohn II via Levy–Lieb constrained search** (`ch05-dft/02-hohenberg-kohn.md`, §5.2.3–§5.2.4b). The motivation for going from $v$-representability to $N$-representability is laid out cleanly; the two-stage variational structure of Levy–Lieb is given pedagogical priority over the historical HK proof. The Gilbert–Harriman conditions and the convexity / lower-semicontinuity properties of $F_L$ are flagged. The Janak's-theorem digression (§5.2.5, lines 411–438) connects HK to piecewise linearity and the derivative discontinuity in exactly the way a modern functional-development course would want.
2. **Body-order vs depth trade-off in MACE** (`ch09-mlip/05-equivariant.md`, §9.5.5, lines 594–642). The receptive-field formula $\nu_\mathrm{total} = 1 + T\cdot(\nu_\mathrm{layer}-1)$ (line 630) is the right invariant to teach, and the methane worked example (lines 462–509) makes it concrete that a regular tetrahedron has zero $\ell=1$ feature. I have not seen this presentation in another textbook and it is genuinely useful.
3. **Nosé extended-system derivation and the chain extension** (`ch07-md/03-thermostats.md`, §7.3, "Nosé-Hoover" subsection, lines 119–230). The transparent statement that *microcanonical sampling of the extended Hamiltonian produces the canonical marginal* — with the Jacobian $s^{-3N}$ argument made explicit — and the follow-up that a single NH chain over-constrains a 1D oscillator (lines 217–227, "Why does this help?") is the right way to motivate Martyna–Klein–Tuckerman. The BAOAB splitting is also given correctly with the analytical OU step.

## 4. Physical errors or misleading framings (three concrete)

1. **`ch03b-solid-state/03-tight-binding.md`, lines 256–260: unedited working notes left in the published text.** The text reads:
   > "$\mathbf K\cdot\boldsymbol\delta_1 = \frac{2\pi}{3}, \quad \mathbf K\cdot\boldsymbol\delta_2 = -\frac{\pi}{3} + \frac{\pi}{3} = 0,\quad \text{wait, let me redo:}$$ Actually, $\mathbf K\cdot\boldsymbol\delta_1 = (2\pi/(3a))\cdot a + (2\pi/(3a\sqrt 3))\cdot 0 = 2\pi/3$. For $\boldsymbol\delta_2$: …"

   This is an LLM scratchpad artefact that should never have made it past copy-editing. It needs deletion and the calculation should be presented cleanly. The downstream conclusion ($f(\mathbf K+\mathbf q)\approx \hbar v_F(q_x - iq_y)$ in eq. 3b.3.26) is correct, but a reader hitting "wait, let me redo" will lose trust in the rest of the book. **Highest-priority fix.**

2. **`ch05-dft/04-xc-functionals.md`, line 437: LDA HOMO of hydrogen claimed at $\approx -7$ eV, "off by ≈ 5 eV from experiment."** The well-known LDA HOMO of H is about $-7.3$ eV (versus the exact $-13.6$ eV), an error of $\approx 6$ eV, not 5. Furthermore the framing ("the LDA HOMO of H is $\approx -7\;\mathrm{eV}$ versus experiment $-13.6\;\mathrm{eV}$") is fine, but it sits inside a Janak-curve discussion where the *slope* at $f=1$ is the operative quantity — and that slope is what equals the (LDA) HOMO eigenvalue, not the HOMO eigenvalue at integer occupation directly. The text conflates the two. A careful sentence would say: "for the exact functional, $\partial E/\partial f|_{f=1^-}$ equals the negative of the ionisation potential, $-13.6$ eV; LDA gives a slope of about $-7.3$ eV at $f=1$, undershooting $I$ by $\sim 6$ eV — the canonical Koopmans / HOMO-IP failure of semilocal DFT." Minor but it is in the most-cited section of the chapter.

3. **`ch12-foundation/02-mace-mp0.md`, lines 44–48: MACE-MP-0 training compute claim.** "The training ran for roughly two weeks on $32$ A100 GPUs." The published MACE-MP-0 paper (Batatia et al. 2024) reports training on substantially fewer GPUs over a longer wall time — single-digit GPUs for several weeks for the medium model. Cross-check the original paper; if the figure stays as-is, cite the source explicitly so the reader can verify. Compute-budget claims propagate into reproduction estimates and should not be approximated from memory.

   Honourable mentions (not in the "top three" but worth a pass during revisions):
   - `ch05-dft/04-xc-functionals.md` line 161, $\mu = 0.21951$ is given to five digits; the standard PBE value is $\mu = 0.2195149\ldots = \beta\pi^{2}/3$ with $\beta = 0.066725$. Either give the exact symbolic form or round to four digits. Minor.
   - `ch07-md/03-thermostats.md` line 56: "Subtract 3 for the centre-of-mass momentum if it is conserved at zero." In a periodic system the COM momentum is conserved by construction (translational invariance of $U$), independent of whether you "set it to zero" by hand. The sentence should clarify that the subtraction reflects the *constraint* you impose during initialisation, not a thermodynamic convention.
   - `ch03b-solid-state/01-bloch-theorem.md` line 178: the basis-set-size estimate at the end of §3b.1.6 confuses energy units. $|G|_\mathrm{max}^{2}$ should have units of inverse length squared; the value quoted ($1.05\times 10^{21}\;\mathrm m^{-2}$) does not match a 400 eV cutoff. Quick check: $|G|_\mathrm{max} = \sqrt{2m_e E_\mathrm{cut}}/\hbar = \sqrt{2\cdot 9.109\times 10^{-31}\cdot 400\cdot 1.602\times 10^{-19}}/1.055\times 10^{-34} \approx 1.03\times 10^{10}\;\mathrm m^{-1}$, i.e. $|G|_\mathrm{max}^{2}\approx 1.05\times 10^{20}\;\mathrm m^{-2}$ — one decade smaller than printed. The Å$^{-1}$ result (≈ 10.3 Å$^{-1}$, not 3.24 Å$^{-1}$ as printed) is also off. Arithmetic glitch in the sanity-check box.

## 5. Five missing references your students would want

These are works that *will* come up in any condensed-matter / MLIP graduate seminar, and that the current bibliography does not list. (Some are mentioned in body text but not entered in `appendix/E-bibliography.md`.)

1. **Geiger, M. & Smidt, T. (2022)** — *e3nn: Euclidean neural networks*, arXiv:2207.09453. The library that every equivariant MLIP — NequIP, MACE, Allegro, SevenNet — uses for Clebsch–Gordan tensor products. Mentioned in the text on lines 213 and 357 of `ch09-mlip/05-equivariant.md` but not entered in Appendix E.
2. **Pozdnyakov, S. & Ceriotti, M. (2020)** — *Incompleteness of atomic structure representations*, Phys. Rev. Lett. 125, 166001. The proof of degenerate environments cited on line 261 of `ch09-mlip/05-equivariant.md` as the theoretical motivation for equivariance. Must be in the bibliography.
3. **Levy, M. (1979)** — *Universal variational functionals of electron densities, first-order density matrices, and natural spin-orbitals…*, PNAS 76, 6062. **Lieb, E. H. (1983)** — *Density functionals for Coulomb systems*, Int. J. Quantum Chem. 24, 243. These are the foundational papers of the Levy–Lieb constrained search developed in detail in §5.2.3–§5.2.4b. Their absence from Appendix E is conspicuous given the prominence of the construction.
4. **Bussi, G., Donadio, D., Parrinello, M. (2007)** — *Canonical sampling through velocity rescaling*, J. Chem. Phys. 126, 014101. The CSVR thermostat recommended in §7.3 (line 117 and the summary table). Cited by name and not entered.
5. **Merchant, A., Batzner, S., Schoenholz, S. S., Aykol, M., Cheon, G., Cubuk, E. D. (2023)** — *Scaling deep learning for materials discovery* (GNoME), Nature 624, 80. The 2023–2024 GNoME paper is the most-discussed scaling-law result for universal MLIPs of recent vintage and is missing from §12. Even if the editors disagree with its claims, students will ask about it.

Honourable mentions (would also be welcomed):
- **Mardirossian, N. & Head-Gordon, M. (2017)** — *Thirty years of density functional theory in computational chemistry*, Mol. Phys. 115, 2315. The standard functional-zoo review.
- **Musaelian, A. et al. (2023)** — *Learning local equivariant representations for large-scale atomistic dynamics* (Allegro), Nat. Commun. 14, 579. Mentioned in §9 (line 716) but not in the bibliography.
- **Martyna, G. J., Klein, M. L., Tuckerman, M. (1992)** — *Nosé–Hoover chains: The canonical ensemble via continuous dynamics*, J. Chem. Phys. 97, 2635. Used by name in §7.3 (line 214 and 234) — needs a formal entry.
- **Martyna, G. J., Tobias, D. J., Klein, M. L. (1994)** — *Constant pressure molecular dynamics algorithms*, J. Chem. Phys. 101, 4177. The MTK paper anchoring §7.3 NPT.
- **Furness, J. W. et al. (2020)** — *Accurate and numerically efficient r²SCAN meta-generalized gradient approximation*, J. Phys. Chem. Lett. 11, 8208. Cited on line 211 of `ch05-dft/04-xc-functionals.md` — entry missing.

## 6. One major modern topic that is under-treated

**Long-range electrostatics, charge transfer, and the dielectric / polarisation problem in MLIPs.** The book mentions, in passing, that MLIPs use a finite cutoff and "Electrostatics in ionic crystals … are at best approximated by the local representation" (`ch12-foundation/02-mace-mp0.md`, line 583). It also cites Ko et al. 2021 (4G-HDNNP) in the bibliography. But the *active* 2024–2026 frontier — fourth-generation HDNNPs, the LODE descriptor (Grisafi & Ceriotti), explicit long-range corrections in MACE (the LR-MACE / MACE-LES line of work), the equivariant prediction of Born effective charges and dielectric tensors needed for LO–TO splittings in DFT phonons and for IR / Raman spectra — gets only one paragraph. For a manuscript that otherwise covers MatterGen, EquiformerV2-OMat and Orb-v2 with care, this is a noticeable gap. A new §9.7 "Long-range interactions and charge-aware MLIPs" — perhaps merged with the "When MLIPs fail" section that the known-limitations document promises — would close it. As a concrete diagnostic that should appear in any such section: a pure MACE-MP-0 calculation of the LO–TO splitting in cubic SrTiO$_3$ or BaTiO$_3$ misses the long-range Fröhlich contribution entirely, by tens of meV. Students should leave the book knowing this.

Closely related but separable: **non-collinear / spin-aware MLIPs.** §12.2 has an excellent paragraph on the "spin-blind MLIP" problem (lines 444–469) but does not point readers at the explicit fixes in development (spin-Heisenberg extensions, the LaMMNI / spin-aware-MACE lines). For magnetic systems this is the dominant practical limitation.

## 7. Smaller editorial / stylistic notes (collected for v1.1)

- `ch05-dft/02-hohenberg-kohn.md` line 467: "Pierre Hohenberg, who was a postdoc at the time, was famously not included" — Hohenberg co-authored the 1964 paper; Pople and Kohn shared the Nobel in 1998. The phrasing "famously not included" reads as editorial; either drop or add a citation, since *why* he was not included is genuinely contested.
- `ch09-mlip/05-equivariant.md` lines 754–789: the worked tensor-product example uses a complex-spherical-basis convention without flagging it explicitly. Many students will arrive from the real-spherical convention (used everywhere else in this book and in e3nn defaults). Add one sentence: "We switch to the complex basis $Y_1^{m}$ for this example because the Clebsch–Gordan coefficients are simpler; e3nn uses the real basis by default, and the two are related by a unitary transformation."
- `ch07-md/03-thermostats.md` lines 126 vs 144: the parameter $g$ is given as $N_\mathrm{df}+1$ in the definition and as $3N+1$ in the proof box. These are the same in $d=3$ without constraints, but the inconsistency in notation is jarring; pick one.
- `ch12-foundation/02-mace-mp0.md` line 130: `default_dtype="float32"` is recommended for exploratory MD; `ch09-mlip/06-training-mace.md` lines 504–510 explicitly warns that `float32` causes NVE drift. These are reconcilable (zero-shot exploration vs production MD) but the conflicting tone confuses students. One sentence cross-referencing the two would help.
- Throughout `ch03b-solid-state/`: the British "Born–von Kármán" appears with a mix of em-dashes and en-dashes; pick the en-dash for compound names per the rest of the book's style.

## 8. What I did *not* check

- Coding / API correctness of the MACE training script (`ch09-mlip/06-training-mace.md` lines 207–306). The reviewer would expect the science editor to compile and run it.
- The accuracy of the Matbench-Discovery $F_1$ scores and $\kappa_\mathrm{SRME}$ values in the Table at `ch12-foundation/02-mace-mp0.md` lines 374–386. Numbers move fast in this area; the table should carry a date stamp ("as of 2025-Q4") if it does not already.
- The exercises (`exercises/`) directory.
- Any chapter outside the eleven files listed in §1.

---

*Submitted in good faith. The work merits acceptance pending the §3 fixes and the bibliography additions in §5.*
