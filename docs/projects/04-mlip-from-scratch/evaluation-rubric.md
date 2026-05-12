# Evaluation rubric — Training an MLIP from scratch

Total: 100 marks. Five categories:

| Category | Weight |
| --- | --- |
| A. Problem framing and baseline | 15 |
| B. Sampling strategy and training set | 20 |
| C. Model training and validation suite | 30 |
| D. Active learning | 15 |
| E. Production MD and write-up | 20 |

Grade boundaries:

| Total | Grade | Descriptor |
| --- | --- | --- |
| 80–100 | Excellent | A novel-system MLIP that passes the full validation suite with documented AL iteration and a 100-ps production trajectory. |
| 60–79 | Good | A trained model that passes most validation tests but with some unresolved issues. |
| 40–59 | Pass | A model that trains and produces parity plots, but MD is unstable or validation is incomplete. |
| < 40 | Fail | No working MLIP, or no validation beyond parity. |

---

## A. Problem framing and baseline (15 marks)

| Marks | Criterion |
| --- | --- |
| 13–15 | The system is clearly motivated and shown to be inadequately covered by MACE-MP-0 (or another foundation model). Baseline E and F MAE are reported with numerical values. The choice between "train from scratch" and "fine-tune" is explicitly justified. |
| 9–12 | The system is described, but the foundation-model baseline is qualitative ("MACE-MP-0 didn't work well") rather than numerical. |
| 5–8 | The system is described but no baseline check was performed. |
| 0–4 | The choice of system is not motivated; the project is "MLIP for whatever I picked first". |

Common deductions: claiming a system is novel without an MP search to
confirm (−3); not reporting the baseline numerically (−3); choosing a
target system covered well by an existing MLIP (−5).

---

## B. Sampling strategy and training set (20 marks)

| Marks | Criterion |
| --- | --- |
| 18–20 | Training set covers equilibrium, two or more thermal regimes (low and high $T$), strained, and rattled configurations. Subsampling addresses AIMD correlation. Held-out test set is curated *before* training and is honestly distinct from the training distribution. Total ≥ 400 frames. Configuration type is recorded for every frame. |
| 13–17 | Training set covers three of the four categories. Test set is held out. Total 250–400 frames. |
| 7–12 | Training set covers only equilibrium + one AIMD trajectory. Test set is a random split. ≤ 250 frames. |
| 0–6 | A single AIMD run treated as the entire training set; no diversity. |

Common deductions: subsampling AIMD at < 5 fs (correlation) (−3);
test-set frames drawn from the same trajectory as training (−5);
no record of `config_type` per frame, making stratified analysis
impossible (−2).

---

## C. Model training and validation suite (30 marks)

This is the most heavily weighted category. The validation hierarchy
is:

1. Parity (energy + force).
2. MD stability (NVE drift + minimum interatomic distance).
3. RDF agreement with AIMD reference.
4. Phonon spectrum agreement with DFT (when applicable).
5. Optional: stress / elastic constants.

| Marks | Criterion |
| --- | --- |
| 26–30 | Model trained with documented hyperparameters. All four validation tests performed, with quantitative comparisons. Force MAE on held-out test set is within 5 % of typical force magnitude. NVE drift < 0.5 meV/atom over 50 ps. RDF first-peak position agrees to ≤ 0.02 Å. Phonon spectrum within 10 % over the dispersion. |
| 19–25 | Parity + MD stability + RDF performed. Force MAE within 10 %. NVE drift < 5 meV/atom. RDF qualitatively reasonable. |
| 12–18 | Parity + one other validation test (typically MD stability). |
| 4–11 | Parity only. |
| 0–3 | No validation; only training-loss curves. |

Common deductions: reporting parity on the training set as if it were
the test set (−6); MD stability test < 10 ps (−3); RDF computed only
on the equilibrium snapshot, not on a thermal trajectory (−3); phonon
calculation done with too-tight thresholds and reports "many small
imaginary modes near $\Gamma$" without recognising these as numerical
artefacts (−1).

---

## D. Active learning (15 marks)

| Marks | Criterion |
| --- | --- |
| 13–15 | At least one AL round performed. Committee of ≥ 2 MACE models. Uncertainty trigger and selection criteria documented. Before/after validation metrics show *quantitative* improvement (e.g., MD stability improves from drift 8 meV/atom to 0.3 meV/atom). |
| 9–12 | One AL round but the before/after improvement is qualitative ("seems better"). |
| 5–8 | AL attempted but no clear improvement reported, *or* committee is just one member with random initialisation (not a real committee). |
| 0–4 | No AL attempted. The student notes "I had no time"; this is honest but does not earn marks. |

Common deductions: AL candidates selected only by max-sigma without
diversity (−2); failing to DFT-relabel the candidates (−5; this is
not active learning); claiming "AL improved my model" without a
controlled comparison on the same metric (−3).

---

## E. Production MD and write-up (20 marks)

| Marks | Criterion |
| --- | --- |
| 17–20 | Final production MD ≥ 100 ps at two temperatures, with one quantitative structural observable (RDF, MSD, ADF). Report is 4–6 pages, cites at least three of Behler, Batatia, Vandermause, Smith, Stocker. Pipeline is reproducible. Limitations honestly stated. |
| 12–16 | Production MD of 50–100 ps at one temperature; one structural observable. Reasonable write-up. |
| 6–11 | Production MD < 50 ps; sparse write-up. |
| 0–5 | No production MD or no write-up. |

Common deductions: production MD that crashes within the run (−5);
not reporting the simulation temperature (−1); using floating-point
seconds where picoseconds were meant (−1); citing reference for the
foundation model but not for MACE itself (−2).

---

## Stretch goals (up to +10, capped at 100)

1. **Free-energy calculation.** Use the final MLIP to compute a free
   energy (e.g., point-defect formation free energy via
   thermodynamic integration or via the quasi-harmonic
   approximation). (+5)
2. **Cross-architecture comparison.** Train an alternative MLIP
   (e.g., ALIGNN, NequIP, or Allegro) on the same dataset and compare
   the validation suite. (+5)
3. **Foundation-model fine-tuning ablation.** Even after training
   from scratch, fine-tune MACE-MP-0 on your dataset and report the
   metrics side-by-side. Discuss the trade-off. (+5)
4. **Transferability test.** Apply your final MLIP to a related but
   unseen system (e.g., a similar Janus dichalcogenide; a different
   peptide sequence) and report how it fails. (+3)

---

## Self-assessment checklist

- [ ] My system is justified with a quantitative MACE-MP-0 baseline.
- [ ] My training set contains ≥ 400 frames across ≥ 4 configuration
      types.
- [ ] My test set is curated before training and not drawn from the
      same trajectories.
- [ ] I trained at least 2 MACE models with different seeds.
- [ ] My parity plot shows per-element force MAE.
- [ ] My MD stability test reports both energy drift and minimum
      interatomic distance.
- [ ] My RDF comparison plots MLIP and AIMD on the same axes.
- [ ] I ran at least one AL round with quantified before/after
      metrics.
- [ ] My production MD ran for ≥ 100 ps without crashing.
- [ ] My report names the residual limitations of the final MLIP.

A complete checklist puts you in the Good–Excellent range.
