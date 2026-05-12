# Evaluation rubric — Bayesian optimisation for catalyst composition

Total: 100 marks. Five categories:

| Category | Weight |
| --- | --- |
| A. Oracle definition and validation | 20 |
| B. BO methodology | 25 |
| C. Baselines and seed averaging | 15 |
| D. Result, regret curve, and posterior | 25 |
| E. Code, write-up, and literature | 15 |

Grade boundaries:

| Total | Grade | Descriptor |
| --- | --- | --- |
| 80–100 | Excellent | A BO loop that demonstrably beats random search, with a validated oracle, a posterior heatmap, and a verified top-5. |
| 60–79 | Good | A working BO loop with reasonable seeds and a documented comparison to random search. |
| 40–59 | Pass | A working BO loop but with weak validation or no statistical comparison to baseline. |
| < 40 | Fail | BO loop does not run, or the oracle is broken. |

---

## A. Oracle definition and validation (20 marks)

| Marks | Criterion |
| --- | --- |
| 18–20 | Oracle is fully defined (slab geometry, site, ZPE+entropy correction, averaging procedure). Validated against literature at all three pure endpoints (Pt, Pd, Ag) with values reported in a table. If MLIP oracle is used, the validation includes DFT cross-checks at the endpoints and at one alloyed composition. |
| 13–17 | Oracle defined and validated at 2 endpoints. |
| 7–12 | Oracle defined but only validated at one endpoint. |
| 0–6 | Oracle is undocumented or unvalidated. |

Common deductions: forgetting the ZPE + entropy correction so that
reported $\Delta G_\mathrm{H}$ are off by 0.24 eV (−4); using a
different adsorption site for different compositions (−3); not
averaging over random arrangements (−2); using a slab thinner than
4 layers (−2).

---

## B. BO methodology (25 marks)

| Marks | Criterion |
| --- | --- |
| 22–25 | GP with Matérn kernel (or another defensible choice). Acquisition function explicitly chosen and motivated. Discrete projection handled cleanly. Initial design via space-filling sampling (LHS or Dirichlet) rather than ad hoc. Duplicates handled. Bounds and constraints explicit. |
| 16–21 | GP with reasonable kernel. Acquisition function used but barely motivated. Discrete projection done. |
| 10–15 | GP set up but kernel hyperparameters fixed (not learnt via MLL). Duplicates not handled. |
| 4–9 | Acquisition function ill-defined; the loop is essentially "pick high-mean prediction" with no exploration term. |
| 0–3 | No GP; "BO" is just sequential random with cherry-picking. |

Common deductions: maximising $\Delta G_\mathrm{H}$ instead of
$-|\Delta G_\mathrm{H}|$ (−6; the objective is the wrong sign); not
fitting GP hyperparameters (−3); not handling the simplex constraint
(−2); using an RBF kernel with very long lengthscales because the
data are too sparse (−1; flag and discuss instead).

---

## C. Baselines and seed averaging (15 marks)

| Marks | Criterion |
| --- | --- |
| 13–15 | At least 3 BO seeds and 3 random-search seeds. Regret curve shows mean ± standard error. Compositions sampled uniformly from the simplex (Dirichlet) for the random baseline, not from $[0, 1]^2$ then rejected. |
| 9–12 | 2 BO + 2 random seeds, or 3 of each but the random baseline samples from the wrong domain. |
| 5–8 | 1 BO + 1 random seed. No statistical comparison. |
| 0–4 | No random-search baseline. |

Common deductions: random search drawn from $[0, 1]^2$ with rejection
(halves effective density; biases against high-fraction endpoints)
(−3); not reporting error bars on the regret curve (−2).

---

## D. Result, regret curve, and posterior (25 marks)

| Marks | Criterion |
| --- | --- |
| 22–25 | Regret curve shows BO strictly outperforming random over most of the iteration range (or, if it does not, the student investigates and explains why). Posterior heatmap rendered on the simplex (ternary or triangular plot, not a square). Top 5 compositions reported with oracle values, and at least 2 verified independently (e.g., with the DFT oracle if MLIP was primary). |
| 16–21 | Regret curve and posterior plotted. Top 5 reported but unverified. |
| 10–15 | Regret curve only; no posterior or top-5 table. |
| 4–9 | A single best-composition value; no curves. |
| 0–3 | No results, or the reported "best" is not the minimum-$|\Delta G_\mathrm{H}|$ point in the history. |

Common deductions: regret curve shown but not on $|\Delta G_\mathrm{H}|$
(−4; this is the right axis, not signed $\Delta G_\mathrm{H}$); BO
indistinguishable from random and no discussion (−3); top-5
verification missing for MLIP-oracle projects (−4).

A *negative* result — BO does *not* beat random for this problem and
budget — is acceptable if (a) the analysis is honest, (b) plausible
reasons are discussed (e.g., the descriptor is too smooth, the
simplex is too small to need BO, GP hyperparameter learning failed
with only 5 initial points), and (c) the student suggests
constructive next steps.

---

## E. Code, write-up, and literature (15 marks)

| Marks | Criterion |
| --- | --- |
| 13–15 | Reproducible pipeline; type hints throughout; configuration centralised. Report cites at least three of Nørskov, Greeley, Frazier, Balandat (BoTorch), Tran–Ulissi. Figures publication-quality. Limitations discussed. |
| 9–12 | Reproducible pipeline. Some literature engagement. Figures present. |
| 4–8 | Pipeline partially reproducible. Sparse citations. |
| 0–3 | No reproducible pipeline; no critical discussion. |

---

## Stretch goals (up to +10, capped at 100)

1. **Multi-fidelity BO.** Implement a two-fidelity BO with MLIP as
   the cheap surrogate and DFT as the expensive oracle; compare
   convergence to single-fidelity BO. (+6)
2. **Categorical extension.** Add a 4th element (Au, Cu, or Ru) and
   run BO on the quaternary simplex with the same total budget. (+4)
3. **Acquisition-function comparison.** Run BO with qEI *and* UCB *and*
   probability-of-improvement on the same problem; compare their
   regret curves. (+3)
4. **GP kernel-comparison study.** Try Matérn 5/2, Matérn 3/2, RBF, and
   a periodic kernel. Quantify which gives the best validation
   log-likelihood. (+3)

---

## Self-assessment checklist

- [ ] My oracle is validated at all three pure endpoints.
- [ ] I report values for $E_\mathrm{slab}$, $E_\mathrm{slab+H}$,
      $E_\mathrm{H_2}$ separately for the validation table, not just
      the final $\Delta G_\mathrm{H}$.
- [ ] I optimise $-|\Delta G_\mathrm{H}|$ in the GP (the right sign).
- [ ] My initial design is space-filling (LHS / Dirichlet), not the
      first 5 random points.
- [ ] I ran ≥ 3 BO seeds and ≥ 3 random-search seeds.
- [ ] My regret curve has mean and standard-error bands.
- [ ] My posterior is plotted on a ternary diagram (or with explicit
      acknowledgement of the simplex shape).
- [ ] My top-5 candidates have been spot-verified.
- [ ] My report includes a discussion of why BO did (or did not) win.

A complete checklist puts you in the Good–Excellent range.
