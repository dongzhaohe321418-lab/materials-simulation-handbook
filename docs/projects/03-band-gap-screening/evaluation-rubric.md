# Evaluation rubric — High-throughput band-gap screening

Total: 100 marks. Five categories:

| Category | Weight |
| --- | --- |
| A. Data acquisition and cleaning | 15 |
| B. Model training and validation | 25 |
| C. Screening pipeline | 20 |
| D. DFT verification | 20 |
| E. Discussion, code, and write-up | 20 |

Grade boundaries:

| Total | Grade | Descriptor |
| --- | --- | --- |
| 80–100 | Excellent | A defensible shortlist of 10 candidates with DFT-verified top 3 and honest novelty/literature analysis. |
| 60–79 | Good | A trained model with MAE < 0.5 eV and a usable shortlist; partial verification. |
| 40–59 | Pass | A trained model with MAE > 0.5 eV but the pipeline produces *some* candidates. |
| < 40 | Fail | No working ensemble; the screen never runs. |

---

## A. Data acquisition and cleaning (15 marks)

| Marks | Criterion |
| --- | --- |
| 13–15 | MP query cached locally with the snapshot date recorded. All cleaning filters documented and applied in a reproducible script. Per-formula train/val/test split. The cleaning script reports counts at each step (`pandas` "before / after" log). |
| 9–12 | Cleaning applied but the split is per-record (with potential polymorph leakage). |
| 5–8 | Cleaning incomplete (e.g., no stability filter, or duplicate formulas not deduplicated). |
| 0–4 | The dataset is downloaded each run, without caching; no documented filters. |

Common deductions: failing to record the MP snapshot date (−2); using
the default random split that does not respect formula groups (−2);
no logging of cleaning step counts (−1).

---

## B. Model training and validation (25 marks)

| Marks | Criterion |
| --- | --- |
| 22–25 | Ensemble of ≥ 5 CGCNNs with different seeds. Overall test-set MAE < 0.5 eV. Stratified MAE reported in at least 5 gap bins, with the [1.5, 2.5] eV bin called out. Parity plot is publication-quality. Calibration plot (predicted σ vs absolute error) included. |
| 16–21 | Ensemble of 3–5 CGCNNs. Overall MAE in [0.5, 0.7] eV. Stratified MAE reported. |
| 10–15 | Single model. Overall MAE in [0.7, 1.0] eV. No stratification. |
| 4–9 | Single model with MAE > 1.0 eV, or model fails to converge. |
| 0–3 | No working model, or test-set leakage detected. |

Common deductions: reporting train MAE as if it were test MAE (−6);
not reporting a stratified MAE (−4); using a CGCNN with one
convolution layer (architecturally crippled) (−3); claiming the
ensemble standard deviation is "uncertainty" without showing a
calibration plot (−2).

---

## C. Screening pipeline (20 marks)

| Marks | Criterion |
| --- | --- |
| 18–20 | Screen set of ≥ 15 000 candidates. Filters applied in a documented order (gap window → uncertainty → stability → element availability). Ranking function explicit and motivated (e.g., $(\mu - 2)^2 + \sigma^2$). Top 10 saved as a CSV with all relevant metadata. Each candidate manually inspected and flagged with a notes column. |
| 13–17 | Screen set of 5 000–15 000. Filters applied but ranking is "by predicted gap distance only". Top 10 saved without manual inspection. |
| 7–12 | Screen set < 5 000 (e.g., just the held-out test set used as a "screen"). Top 10 produced mechanically. |
| 0–6 | No screen step, or top-10 is just the lowest-MAE test points (a methodological error). |

Common deductions: forgetting an uncertainty filter so that
high-σ candidates dominate the top 10 (−3); ranking by raw μ
without distance to the target gap (−2); not recording the screen
source (MP-not-in-train vs prototype substitution) (−2).

---

## D. DFT verification (20 marks)

| Marks | Criterion |
| --- | --- |
| 18–20 | At least 3 candidates relaxed with PBE; SCF on dense **k**-mesh; band gap extracted and compared to CGCNN prediction. Discrepancy discussed in terms of CGCNN uncertainty. (Optional) One HSE06 calculation. |
| 13–17 | 2 candidates verified with PBE; comparison reported. |
| 7–12 | 1 candidate verified; or 3 candidates with single-point SCF only (no relaxation). |
| 0–6 | No DFT verification, or DFT inputs prepared but never run. |

Common deductions: comparing CGCNN prediction to a single-point SCF
on an unrelaxed structure (−4); not stating the **k**-mesh used for
the gap evaluation (−2); using the MP-relaxed structure unchanged
and treating that as "verification" (−5; this is just re-reading
the MP value the model was trained on).

---

## E. Discussion, code, and write-up (20 marks)

| Marks | Criterion |
| --- | --- |
| 17–20 | Reproducible pipeline. Report cites at least three of Xie–Grossman, Jain–MP, Castelli, Tran–Ulissi, Pilania. Honest discussion of (i) PBE band-gap underestimation and what it implies for the 1.5–2.5 eV target, (ii) how many shortlisted candidates are already known photocatalysts (a brief literature check), and (iii) the limitations of stability filtering by `energy_above_hull`. Figures publication-quality. |
| 12–16 | Reproducible pipeline. Some literature engagement. Honest discussion of one of the three points above. |
| 6–11 | Pipeline partially reproducible. Sparse citations. |
| 0–5 | No reproducible pipeline; no critical discussion. |

Common deductions: claiming "novel candidates discovered" without
a literature search (−3); writing "we discovered" when the shortlist
includes BiVO$_4$ or TiO$_2$ (the model is recovering known
photocatalysts — celebrate this as a sanity check, not a discovery)
(−2).

---

## Stretch goals (up to +10, capped at 100)

1. Swap CGCNN for ALIGNN or MEGNet and re-report the stratified MAE.
   Discuss what improves. (+5)
2. Apply a learned PBE → HSE correction (e.g., a per-element-class
   linear shift) and re-rank. Discuss how the shortlist changes. (+3)
3. Run one HSE06 verification on a top candidate. Report the
   HSE/PBE gap ratio and compare with the empirical 1.3× rule. (+5)
4. Include a "false positive" analysis: re-predict the gap of three
   *known* materials with HSE gaps outside the 1.5–2.5 eV window;
   if the CGCNN places them inside the window, your model has bias
   you should discuss. (+3)

---

## Self-assessment checklist

- [ ] My MP query is cached and the snapshot date is recorded.
- [ ] My train/val/test split is per-formula (no polymorph leakage).
- [ ] I trained ≥ 3 CGCNNs with different random seeds.
- [ ] I report a stratified MAE table with at least 5 bins.
- [ ] My parity plot has equal axes.
- [ ] My screen-set predictions include both mean and standard
      deviation per candidate.
- [ ] My top-10 CSV has at least seven columns of metadata.
- [ ] I relaxed at least 3 candidates with PBE and report the
      relaxed band gap.
- [ ] My report discusses PBE band-gap underestimation explicitly.
- [ ] My report includes a literature check of every shortlisted
      candidate.

A complete checklist puts you in the Good–Excellent range.
