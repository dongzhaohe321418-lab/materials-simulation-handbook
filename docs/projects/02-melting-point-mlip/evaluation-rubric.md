# Evaluation rubric — Melting point of copper via MLIP-driven MD

Total: 100 marks. The five categories are:

| Category | Weight |
| --- | --- |
| A. Training-set quality | 20 |
| B. MLIP validation | 20 |
| C. Coexistence methodology | 25 |
| D. Result and uncertainty | 20 |
| E. Code, write-up, and literature | 15 |

Grade boundaries:

| Total | Grade | Descriptor |
| --- | --- | --- |
| 80–100 | Excellent | A defensible $T_m$ within 50 K of experiment, with documented uncertainty. |
| 60–79 | Good | A working pipeline with a $T_m$ estimate within 200 K. |
| 40–59 | Pass | The MLIP is trained and equilibrates a single-phase cell, but $T_m$ is unreliable. |
| < 40 | Fail | The pipeline does not run end-to-end. |

---

## A. Training-set quality (20 marks)

| Marks | Criterion |
| --- | --- |
| 18–20 | Training data span at least three thermodynamic states (cold solid, hot solid, liquid) plus volumetric strain. Each state is sampled long enough (≥ 1.5 ps) to be statistically independent. Pre-melt protocol used for the liquid. Total ≥ 400 frames. The sampling rationale is documented. |
| 13–17 | Three states sampled but with shorter trajectories (< 1 ps) or with the liquid initiated from the solid without pre-melting. ≥ 300 frames. |
| 8–12 | Only two thermodynamic states (e.g., 300 K and 1500 K). ≤ 250 frames. |
| 0–7 | A single thermodynamic state, or no AIMD generation at all (e.g., the student tried to reuse an off-the-shelf foundation potential without retraining). |

Common deductions: failing to pre-melt the liquid sample (the
trajectory does not represent a true equilibrium liquid) (−4); no
strained cells (will fail elastic-constant test later) (−3).

---

## B. MLIP validation (20 marks)

| Marks | Criterion |
| --- | --- |
| 18–20 | Held-out test set with energies and forces from configurations not seen in training. Parity plots present. Reported metrics: E MAE ≤ 5 meV/atom and F MAE ≤ 50 meV/Å. A separate validation on the *liquid only* shows that liquid-state errors are quantified explicitly (typically 1.5× higher than solid). |
| 13–17 | Parity plots and metrics present. E MAE ≤ 10 meV/atom and F MAE ≤ 80 meV/Å. |
| 8–12 | Parity plots present but with single metric only (no liquid-vs-solid breakdown). Errors exceed the target by ≤ 2×. |
| 0–7 | No held-out test set, or test set leakage (frames from the same trajectory appear in train and test). |

Common deductions: reporting train-set metrics as if they were test-set
metrics (−5); using `mae` of *total* energy instead of *per-atom*
energy (−2); using float64 in MD where the MACE model is float32, so
that drift is masked (−1).

---

## C. Coexistence methodology (25 marks)

| Marks | Criterion |
| --- | --- |
| 22–25 | Elongated cell ($L_z \gtrsim 100$ Å, ≥ 6000 atoms). Half-melt protocol applied correctly. NPT at 1 bar with appropriate thermostat/barostat times. ≥ 4 trial temperatures bracketing $T_m$, each ≥ 200 ps. Interface position is identified via a defensible order parameter ($q_6$, coordination, or local-structure index). Linear regression of $z_\mathrm{interface}(t)$ is shown with its residuals. |
| 16–21 | The above but with three temperatures or 100-ps runs. Interface position identified but the order parameter is somewhat noisy. |
| 10–15 | Two temperatures or short (≤ 50 ps) runs; the interface velocity is reported as a single number without an uncertainty estimate. |
| 4–9 | Coexistence cell is set up but never properly equilibrated, or analyses use the bulk volume rather than an explicit interface velocity. |
| 0–3 | No coexistence simulation; only hysteresis is attempted. |

Common deductions: NVT instead of NPT (the cell volume cannot adapt
to the changing solid:liquid fraction) (−6); $L_z < 50$ Å (periodic
self-interaction) (−4); using a *coordination-number cutoff* without
checking that the chosen cutoff cleanly separates solid from liquid
(−2).

---

## D. Result and uncertainty (20 marks)

| Marks | Criterion |
| --- | --- |
| 18–20 | $T_m$ reported with a documented standard error from at least two sources: the linear fit of $v(T) = 0$ and a repeat coexistence run with a different random seed at the bracketed $T_m$. Reported $T_m$ within 100 K of experiment (1358 K). |
| 13–17 | $T_m$ within 150 K of experiment, with an uncertainty estimate from only one source. |
| 8–12 | $T_m$ within 250 K of experiment; no uncertainty estimate. |
| 0–7 | $T_m$ outside 250 K of experiment, or only a hysteresis bracket reported (no coexistence estimate). |

Common deductions: reporting more significant figures than the
uncertainty supports (−2); not stating the pressure at which $T_m$ is
defined (−1); not comparing against the EAM-Foiles baseline (−2).

A *result better than experiment* is suspicious. PBE typically
overestimates $T_m$ for Cu by 50–100 K; reporting a value within ±10 K
of 1358 K may reflect a serendipitous cancellation rather than a
genuinely accurate calculation. The student should at least
acknowledge this in the discussion.

---

## E. Code, write-up, and literature (15 marks)

| Marks | Criterion |
| --- | --- |
| 13–15 | Reproducible pipeline (each phase regenerable from scripts). Type hints throughout. Report cites at least three of Morris et al., Belonoshko et al., Foiles et al., Batatia et al., Vega et al. Figures are publication-quality. EAM-Foiles baseline included. |
| 9–12 | Reproducible pipeline. At least one literature comparison. Figures present but unpolished. |
| 4–8 | Pipeline partially reproducible. Citations sparse. |
| 0–3 | No reproducible pipeline; ad-hoc shell scripts only. |

---

## Stretch goals (up to +10 marks, capped at 100)

1. Compute the melting curve $T_m(p)$ at one additional pressure (e.g.,
   10 GPa). Compare the slope with the Clausius–Clapeyron prediction
   $dT_m/dp = T_m \Delta V / \Delta H$, using the volume and enthalpy
   jumps from your own simulations. (+5)
2. Run an *active learning* loop: identify high-uncertainty
   configurations from the production MD (using a committee of two
   MACE models trained with different random seeds), re-label with
   DFT, and retrain. Report the improvement in liquid-state F MAE
   and the corresponding refinement in $T_m$. (+5)
3. Include a third comparison: the **Foiles EAM** and a **foundation
   model** like MACE-MP-0 (off-the-shelf, no retraining). Tabulate
   $T_m$ for all three approaches. (+5)

---

## Self-assessment checklist

- [ ] My training set spans cold solid, hot solid, and equilibrium
      liquid.
- [ ] My parity-plot script writes both `parity_energy.png` and
      `parity_force.png` with axes labelled in meV/atom and eV/Å.
- [ ] My coexistence cell is at least 6000 atoms and 100 Å long.
- [ ] I ran at least four temperatures, each for ≥ 200 ps, in NPT at
      1 bar.
- [ ] My interface-position time series is plotted as a separate
      figure, not buried in the report.
- [ ] My final $T_m$ has an uncertainty estimate derived from two
      independent sources.
- [ ] I compare with experiment (1358 K) and with at least one
      published MLIP or EAM value.

A complete checklist puts you in the Good–Excellent range.
