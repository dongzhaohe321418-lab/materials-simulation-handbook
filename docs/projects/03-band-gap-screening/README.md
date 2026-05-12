# Project 3 — High-throughput band-gap screening of oxides

## Research question

Can you train a Crystal Graph Convolutional Neural Network (CGCNN) on
roughly 5000 oxides from the Materials Project, use it to screen
≈ 20 000 hypothetical or untested oxide candidates for band gaps in the
photocatalytically interesting range (1.5–2.5 eV), and then verify
your top 10 with DFT? You will compare your shortlist against the
known photocatalyst literature.

This project is a small-scale but methodologically complete instance
of *in-silico* materials discovery: data acquisition, surrogate model
training, large-scale prediction, top-K selection, and high-fidelity
verification of the shortlist.

## Why this project

Several characteristics make oxide band-gap screening a friendly first
high-throughput project:

- Materials Project provides clean, queryable oxide data with PBE
  (and, increasingly, HSE06) gaps. The API is well-documented and
  rate-limited but free.
- The target — a band gap in a narrow window — is a simple scalar
  with a known physical meaning, not a derived quantity that requires
  multiple post-processing steps.
- CGCNN is a battle-tested architecture with abundant reference
  implementations; you will not be debugging novel architectures.
- The full pipeline (query → train → predict → verify) is short
  enough to fit in eight weeks, yet exposes you to every step that a
  real screening effort needs.

The honest caveat: PBE systematically underestimates band gaps by
about 40 %. We work around this either by training on HSE06 data
where available, or — more commonly — by applying a learned
correction. Both approaches are described in the methods.

## Expected outcomes

By the end of the project you will deliver:

- A trained CGCNN with held-out test-set MAE < 0.5 eV on PBE band
  gaps. Excellent submissions also report a *stratified* MAE in the
  1.5–2.5 eV target window.
- A shortlist of 10 candidate oxides predicted to fall in the
  1.5–2.5 eV photocatalyst window, ranked by predicted gap, with
  uncertainty estimates from a CGCNN ensemble.
- DFT verification (PBE relaxation + band-gap calculation) of at
  least 3 of the top 10 candidates, with the discrepancy between
  CGCNN and DFT documented honestly.
- A discussion of how many of your shortlisted materials are *novel*
  (not in MP), how many were already studied as photocatalysts (you
  will have to do a literature check), and how many are spurious
  (decompose, are unstable, or contain rare elements).
- A short report (4–6 pages) describing the pipeline and presenting
  the shortlist.

## Time estimate

Eight weeks, part-time.

| Week | Activity |
| --- | --- |
| 1 | Background reading. Register a Materials Project API key. Build an oxide query, cache locally. |
| 2 | Cleanup: deduplicate by reduced formula + spacegroup, remove magnetic compounds where MP's calculation is unstable, partition train/val/test. |
| 3 | Implement / install CGCNN. Train a first model. Inspect parity plot. |
| 4 | Ensemble of 5 CGCNNs (different random seeds) for uncertainty. Re-validate. |
| 5 | Construct the 20 000-candidate screen set. Predict gap + uncertainty for each. |
| 6 | Rank candidates; filter by stability (energy_above_hull < 0.1 eV) and by uncertainty. |
| 7 | Pick top 10; run PBE relaxation + band-gap calculation in QE for top 3. |
| 8 | Comparison, literature lookup of the shortlist, write-up. |

## Compute budget

- Modest GPU (e.g., RTX 3090, RTX 4090, A100, or even a V100). Each
  CGCNN trains in 2–4 hours; an ensemble of 5 is 10–20 GPU-hours.
- ≈ 200 CPU-hours for the DFT verification step. A typical oxide
  relaxation is a 20–40-atom cell, 5–10 CPU-hours; for 3 candidates
  including a non-self-consistent dense **k**-mesh band-structure
  step, plan 100 CPU-hours, with margin.

## Prerequisites

- [Chapter 5 — Density functional theory](../../ch05-dft/index.md) and
  [Chapter 6 — Running DFT in practice](../../ch06-running-dft/index.md)
  for the underlying DFT band-gap concepts (the gap as the difference
  between band edges, the Kohn–Sham gap problem).
- [Chapter 10 — Graph neural networks for crystals](../../ch10-gnn/index.md)
  for the CGCNN architecture, message passing, and the role of
  per-node and per-edge features.

You should be comfortable with PyTorch (or PyTorch Geometric), with
`pandas`, and with at least the basics of crystal symmetry and
band-structure terminology.

## What "photocatalyst gap" actually means

The 1.5–2.5 eV window is the rough range where:

- The gap is large enough to straddle the water-oxidation and
  water-reduction potentials (O$_2$/H$_2$O at +1.23 V and H$_2$/H$^+$
  at 0 V, separated by 1.23 V; plus an overpotential margin of
  ≈ 0.5 V on each side).
- The gap is small enough that the material absorbs a useful fraction
  of solar photons (the AM1.5 spectrum peaks around 1.6–2.0 eV).

PBE-underestimated gaps in this range will most likely correspond to
HSE06 gaps of 2.0–3.0 eV — which is precisely where the most-studied
photocatalysts (anatase TiO$_2$, BiVO$_4$, $\alpha$-Fe$_2$O$_3$) lie.
This is no accident: PBE shifts the entire shortlist consistently.
You should account for this in your interpretation.

## Pitfalls flagged up front

1. **Caching the MP query.** The Materials Project API rate-limits
   you to a few thousand requests per day. Cache your query results
   locally (a Parquet or pickled file) and never re-query
   unnecessarily.
2. **PBE-zero-gap entries.** Metallic compounds have band gap = 0 in
   MP. Decide whether to keep them (you can frame the model as
   "predict any gap, including zero") or drop them (you may want a
   gap > 0.1 eV filter). Most students keep them; the dataset is more
   varied.
3. **Stratified evaluation.** A CGCNN with MAE 0.45 eV averaged over
   the whole dataset may have MAE 0.8 eV in the 1.5–2.5 eV slice —
   precisely where you care. Always report stratified errors.
4. **Stability filter.** A predicted-low-gap candidate that has
   `energy_above_hull` > 0.1 eV/atom is almost certainly thermodynamically
   unstable. Filter on stability before ranking.
5. **Atomic-fraction leakage.** If your train and test sets share the
   same chemical formula but different polymorphs, you may be
   inflating apparent performance. Consider a strict per-formula
   split for one of your held-out sets.
6. **DFT verification trap.** A PBE relaxation of a downloaded
   structure can change its band gap by 0.3 eV. Always relax before
   reporting; do not just take the MP relaxed structure and run a
   single-point SCF.

## Deliverables checklist

- [ ] `data/mp_oxides.parquet` — cached MP query results.
- [ ] `data/train.json`, `data/val.json`, `data/test.json` — splits.
- [ ] `models/cgcnn_seed{0..4}.pt` — five trained models for the
      ensemble.
- [ ] `analysis/parity.png` — predicted vs MP gap on the test set.
- [ ] `analysis/stratified_mae.csv` — MAE in 0.5-eV-wide gap bins.
- [ ] `screen/candidates.parquet` — predictions for the 20 000-candidate
      set with ensemble means and standard deviations.
- [ ] `shortlist/top10.csv` — the final shortlist with predictions,
      uncertainties, and stability filters.
- [ ] `dft/<formula>/` — PBE inputs and outputs for at least 3 verified
      candidates.
- [ ] `report.pdf`.
