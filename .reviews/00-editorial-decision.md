# Editorial Decision Letter

**Manuscript**: Materials Simulation Handbook (v1.3 release candidate)
**Date**: 2026-05-12
**Reviewers**: EIC (commissioning editor), R1 Methodology, R2 Domain, R3 Perspective, Devil's Advocate

## Aggregate scores

| Reviewer | Score | Recommendation |
|---|---|---|
| EIC | 7 / 10 | Minor Revision |
| Methodology | 8.3 / 10 | Accept with minor revisions |
| Domain | 8.1 / 10 | Strong recommend, minor revisions |
| Perspective | 6.6 / 10 | Major Revision |
| Devil's Advocate | — | 8 CRITICAL + 7 MAJOR issues |

**Iron Rule applied**: Devil's Advocate found CRITICAL issues → Decision cannot be Accept. **Decision: Minor Revision**, with all P0 items below as mandatory.

## P0 — MUST FIX before v1.3 push (consensus of ≥2 reviewers)

| # | File | Issue | Reviewers flagging |
|---|---|---|---|
| 1 | `README.md` line 18 | Three typos: "ideas starts", "mist importantly", "thhe" | EIC, DA |
| 2 | `docs/ch03b-solid-state/03-tight-binding.md` lines 256-260 | **LLM scratchpad "wait, let me redo:" leaked into published text** | Domain |
| 3 | `docs/ch03b-solid-state/01-bloch-theorem.md` line 178 | ENCUT 400 eV worked example off by factor of 10 | Domain |
| 4 | Capstone project metadata (3 places) | Durations 4/6/6/8/10 vs 6/8/8/8-10/6 vs 6/8/8/10/6; folder names disagree; bare github.com links in `docs/ch14-capstone/07-the-five-projects.md` | EIC, DA |
| 5 | `docs/ch00-math/01-numbers.md` lines 288, 294 | Promises Ch 1 covers recursion / DFT scaling / MC error — Ch 1 is Python install, none of those topics appear | DA |
| 6 | Ch 3.5 / Ch 4 ordering contradiction across `learning-path.md` and `docs/ch03b-solid-state/index.md` | Path A puts Ch 3.5 before Ch 4; Ch 3.5 index says "read after Ch 4" with Ch 4 content listed as prereqs | Perspective |
| 7 | `docs/index.md` line 13 | "you do not need a physics or chemistry degree" contradicts the A-level prereq floor stated in README and prerequisites-checker | EIC, Perspective |

## P1 — should fix in v1.3 (single-reviewer high-severity)

| # | File | Issue | Reviewer |
|---|---|---|---|
| 8 | `docs/ch05-dft/04-xc-functionals.md` line 437 | Janak slope conflated with HOMO eigenvalue at integer occupation; LDA HOMO of H stated ~5 eV (should be ~6 eV) | Domain |
| 9 | `docs/ch12-foundation/02-mace-mp0.md` lines 44-48 | MACE-MP-0 training compute claim "32 A100s × 2 weeks" appears inflated vs original paper | Domain |
| 10 | `code/tier1/ch11/gp.py` lines 53-55, 68-70 | Silent auto-transpose of `X` when dimensions mismatch — bug magnet | Methodology |
| 11 | `code/tier1/ch05/scf_1d.py` `kinetic_matrix` | Uses dense storage; prose explicitly recommends sparse | Methodology |
| 12 | `scripts/figures/fig_rdf_real.py` lines 246-249 | Same density 0.030 Å⁻³ for solid / liquid / gas; the 600 K "gas" is really a supercritical fluid at this density | Methodology |
| 13 | `.github/workflows/deploy.yml` | Deploy job doesn't depend on test job — red tests can ship green site | Methodology |

## P2 — defer to v1.4 (improvements, not defects)

- 5 missing references the field expects: e3nn (Geiger & Smidt 2022), Pozdnyakov & Ceriotti 2020, Levy 1979 / Lieb 1983 (constrained search originals), Bussi-Donadio-Parrinello 2007 (CSVR), GNoME (Merchant 2023)
- Long-range electrostatics / 4G-HDNNP / LODE / LR-MACE — currently one paragraph in Ch 12, should be a whole subsection in Ch 12.4
- "Pause and recall" retrieval boxes (≥3 per chapter) per Perspective reviewer
- Ch 9 MLIP failure modes section expansion (12 lines is asymmetric vs Ch 5.6 at 280 lines)
- Reorder Ch 6.1: defer NC/USPP/PAW taxonomy to AFTER §6.2 first calculation
- Replace prerequisites-checker self-graded short-answer with explicit decision rubrics
- Ch 11.2: move intuition box before formal GP definition
- Promote `known-limitations.md` to Start-Here block with explicit version + timestamp

## Synthesizer arbitration notes

- Methodology and Domain reviewers agree the *physics* and *code* are largely sound; the issues are polish-level, not foundational.
- Perspective reviewer is harshest (6.6) because pedagogical inconsistencies (the 3.5/4 ordering, the "zero-foundation" claim, the prereq-checker grading) compound for a real new reader.
- Devil's Advocate is hardest on broken promises and project-link integrity. These are the highest-priority fixes by ROI (each takes minutes; impact is high credibility).
- No reviewer recommends Reject. All agree this is a v1.0 publishable open educational resource pending the P0 fixes.

## Path forward

Editorial Synthesizer (me) will now:
1. Apply all 7 P0 fixes inline
2. Apply 5 of 6 P1 fixes (skip the MACE-MP-0 compute number — needs the original paper to verify)
3. Commit P0/P1 batch with explicit reference back to this Decision Letter
4. Re-run `mkdocs build --strict`, pytest, and a final reviewer pass (Devil's Advocate only — verify CRITICAL items now resolved)
5. Push to GitHub Pages
6. Move P2 to a tracked GitHub Issues list for the v1.4 milestone
