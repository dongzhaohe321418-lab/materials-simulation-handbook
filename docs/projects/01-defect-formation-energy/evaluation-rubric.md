# Evaluation rubric — Defect formation energy in silicon

This rubric is intended both for the student (so you know what
"finished" looks like) and for the assessor. Total: 100 marks.

Marks are awarded in five categories:

| Category | Weight |
| --- | --- |
| A. Numerical convergence | 20 |
| B. Defect physics correctness | 25 |
| C. Finite-size extrapolation | 20 |
| D. Code quality and reproducibility | 15 |
| E. Write-up and engagement with literature | 20 |

Letter grades are mapped as follows:

| Total mark | Grade | Descriptor |
| --- | --- | --- |
| 80–100 | Excellent | Goes meaningfully beyond a textbook reproduction. |
| 60–79 | Good | A solid, defensible answer for the neutral defect. |
| 40–59 | Pass | The calculation runs but is not fully converged. |
| < 40 | Fail | Either no working pipeline or fundamental physics errors. |

---

## A. Numerical convergence (20 marks)

What is being assessed: do you know how converged your numbers are,
and can you defend the values you chose?

| Marks | Criterion |
| --- | --- |
| 20 | Plane-wave cut-off and **k**-mesh both converged to ≤ 1 meV/atom, demonstrated by explicit sweeps with tabulated data. Equivalent-**k**-mesh prescription explicitly applied across cells. Smearing chosen and justified. |
| 14–19 | Convergence demonstrated for one of cut-off or **k**-mesh; the other is justified by citing SSSP recommendations. |
| 8–13 | A single point computed at each cut-off / mesh value; the choice is asserted rather than defended. |
| 0–7 | No convergence study; parameters copied without justification. |

Common deductions: using different cut-offs for bulk and defect (−5);
forgetting to use the equivalent **k**-mesh (−5); reporting energies
without an estimate of numerical noise (−3).

---

## B. Defect physics correctness (25 marks)

What is being assessed: does your final number correspond to the
right physical object?

| Marks | Criterion |
| --- | --- |
| 22–25 | Neutral vacancy relaxed *and* shown to exhibit the $D_{2d}$ Jahn–Teller distortion. Bond-length pairing reported numerically. The formation-energy formula is stated explicitly with each term identified. At least one charge state ($q = +1$ or $-1$) also computed, with a Freysoldt-style correction applied and the correction term reported separately. |
| 16–21 | Neutral vacancy relaxed and JT distortion observed. Formation-energy formula stated correctly. No charged-defect work. |
| 10–15 | Neutral vacancy computed but $T_d$-locked (no JT). Reported $E_f$ is therefore ≈ 0.4 eV too high but the workflow is mechanically correct. |
| 4–9 | Sign or atom-count error in the formation-energy formula, *or* inconsistent pseudopotentials between bulk and defect runs. |
| 0–3 | No defect calculation, or the reported "formation energy" has no defensible physical interpretation. |

Common deductions: not perturbing the initial positions and so missing
the JT distortion (−4); using `degauss = 0` and seeing SCF
non-convergence (−2); relaxing the cell parameter of the defect cell
(−5).

---

## C. Finite-size extrapolation (20 marks)

What is being assessed: do you take seriously the limit of an
isolated defect?

| Marks | Criterion |
| --- | --- |
| 18–20 | At least three supercell sizes computed (e.g., 64, 216, 512 atoms). $E_f$ vs $1/N$ plotted with a fit, intercept reported with a standard error ≤ 0.05 eV. A complementary $L^{-3}$ fit is shown for cross-checking. A qualitative discussion of the source of the dominant correction (elastic vs electronic). |
| 13–17 | Two supercell sizes computed; linear extrapolation performed but with only one degree of freedom in the fit. The student acknowledges this limitation. |
| 7–12 | Only one supercell size; no extrapolation. The reported number is presented as if it were converged. |
| 0–6 | No discussion of supercell scaling at all. |

Common deductions: extrapolating from two unequally-converged cells
without weighting (−3); reporting an uncertainty smaller than the
fit's standard error (−4).

---

## D. Code quality and reproducibility (15 marks)

What is being assessed: could a second student reproduce your numbers
from your repository?

| Marks | Criterion |
| --- | --- |
| 13–15 | The repository contains a single configuration source, all scripts are runnable with `python -m`, random seeds for the perturbation are recorded, every `pw.x` input is regenerable from the scripts, and there is a `README` describing the run order. The code passes a linter (e.g., `ruff`). Type hints throughout. |
| 9–12 | The repository is reproducible in principle, but some parameters are hard-coded in multiple files; a few `pw.x` inputs are committed by hand. |
| 4–8 | The repository contains output files but no scripts that regenerate them, or several scripts duplicate parameter definitions. |
| 0–3 | Output files only; no scripts. |

A common positive signal: a small `tests/` directory with a single
fast unit test (e.g., that `make_vacancy` removes exactly one atom and
that the formation-energy formula reduces to zero for a defect-free
cell). Add three marks if present.

---

## E. Write-up and engagement with literature (20 marks)

What is being assessed: can you place your numbers in context?

| Marks | Criterion |
| --- | --- |
| 17–20 | The report cites at least two of Corsetti–Mostofi, Probert–Payne, Wright, Puska et al., and compares the extrapolated $E_f^\infty$ value quantitatively. Plot is publication-quality with clear axes, units, error bars, and a fit line. The Jahn–Teller distortion is discussed qualitatively. Limitations are stated honestly. |
| 12–16 | At least one literature comparison; figures are reasonable but not publication-quality. Some discussion of limitations. |
| 6–11 | A bare statement of the result with one literature citation; figures lack uncertainty estimates. |
| 0–5 | No literature comparison or no figures. |

Common deductions: citing a literature value without checking which
functional was used (−2); reporting a number with more significant
figures than the convergence supports (−2).

---

## Boundary cases

### What separates "Excellent" from "Good"

Two things in practice push a report from Good (60–79) to Excellent (≥80):

1. At least one *charge-state* calculation with a proper electrostatic
   correction. This is the natural next step after the neutral defect
   and demonstrates that you have understood the full formation-energy
   formalism.
2. A genuinely careful uncertainty analysis. Most undergraduate reports
   quote a single decimal place with no error bar; the excellent report
   identifies the dominant source of error (often the **k**-mesh) and
   estimates its contribution.

### What separates "Pass" from "Fail"

A failing report typically exhibits one of:

- The defect cell and the bulk cell were computed with different
  pseudopotentials, cut-offs, or smearing settings, making the energy
  difference meaningless.
- The reported formation energy is negative, or wildly different from
  the literature (e.g., 1 eV or 8 eV), with no discussion of why.
- The Python code does not run, or there is no Python code at all.

A passing report has a working pipeline, a sensible (if unconverged)
number, and acknowledges what is not yet converged.

---

## Stretch goals — only attempt after the core deliverables

Each adds up to 5 marks (and may push you over 100; the cap is 100):

1. Compute the migration barrier between two equivalent vacancy sites
   using the NEB method. Compare with Watkins's experimental value
   (≈ 0.4 eV for $V^0$).
2. Run a spin-polarised calculation on $V^{-1}$ and verify the
   negative-U behaviour (the $-1/+1$ thermodynamic transition level
   sits below the $0/-1$ level).
3. Compute the formation energy with a hybrid functional (HSE06) on
   the 64-atom cell. Report the difference against PBE and discuss the
   gap-error contribution.

---

## Self-assessment checklist (use *before* submitting)

- [ ] My `analysis/scaling.png` plot has axes labelled with units.
- [ ] My report explicitly states the pseudopotential file used and
      cites the SSSP paper.
- [ ] My extrapolated value is within 0.1 eV of at least one published
      PBE value, and I cite that publication.
- [ ] My code runs from a clean checkout in under one command per
      stage.
- [ ] My figures are vector or ≥ 200 dpi raster, not low-resolution
      screenshots.
- [ ] I have stated which charge state and which spin state my
      calculation corresponds to.
- [ ] My finite-size fit reports a standard error, not just an
      intercept.

If you can tick every box, you are in the Good–Excellent range. If you
miss two or more, expect to land in the Pass band.
