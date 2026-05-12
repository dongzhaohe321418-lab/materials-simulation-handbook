# EIC Review

**Reviewer persona**: Senior commissioning editor, computational science textbooks (CUP / Springer open educational resources programme)
**Date**: 2026-05-12

## Summary recommendation

This is an unusually ambitious and, on the whole, unusually coherent piece
of work for an open textbook at v1.0. The three-tier architecture (`docs/learning-path.md`),
the explicit prerequisite floor at A-level / first-year university
(`README.md` §"Who is this for?"), the dependency graph in
`docs/learning-path.md`, the 20-question self-test
(`docs/prerequisites-checker.md`), the deliberately bounded compromise of
Ch 3.5 between Ashcroft–Mermin and Sholl–Steckel, and the explicit
"Known Limitations" page (`docs/known-limitations.md`) together amount to
an editorial concept the field does not yet have a competitor for. The
book sits in the genuinely empty quadrant of "zero-floor /
foundation-model-ceiling" pedagogy. My overall verdict is **Minor
Revision**: the manuscript is publishable in essentially its current form
as a v1.0 open educational resource, but it cannot ship to a paying
imprint until a small number of production defects — internal numerical
inconsistencies, several typographical errors in the front matter, and
two structural over-claims — are corrected. Confidence: high on
structure and originality; medium on production polish (sampled, not
exhaustively audited).

## Strengths (3)

1. **The three-tier promise is genuinely delivered, and the seams are
   visible by design.** The dependency graph in `docs/learning-path.md`
   lines 13–67, the three named reading paths (A linear, B deep-core, C
   project-driven, lines 103–151), the cross-reference matrix (lines
   159–177), and the explicit "what you can skip on a first pass" list
   (lines 184–191) are exactly the editorial infrastructure that
   distinguishes a *taught* text from a *recorded* one. Most competing
   open textbooks (e.g. Marder, Sholl–Steckel) assume the reader
   reconstructs this scaffolding themselves; this book hands it over.

2. **Honest, dated, opinionated voice on a moving field.** Chapter 12
   opens (line 21) with "as of 2026, [the answer is] a qualified *yes*"
   on whether foundation models are the GPT moment for materials science,
   and §12.4 (lines 478–500) closes with a deliberately measured
   summary that explicitly lists the open frontiers — long-range
   interactions, charge transfer, synthesisability, autonomous labs.
   The willingness to date claims is rare in textbooks and is precisely
   what an open resource — updatable on push to `main` — can do that a
   print edition cannot. Similarly the `docs/known-limitations.md` page
   (in particular the candid §1 admission that Tier 0 linear algebra
   does not cover Hilbert-space technique demanded later, and §3 on the
   cluster-access gap) is the kind of honesty that earns reader trust.

3. **The five-project capstone scaffold is a publishable contribution
   in its own right.** Each project README — sampled here in
   `docs/projects/01-defect-formation-energy/README.md` — supplies a
   research question, an *uncertainty target* (≤ 0.05 eV on the
   vacancy formation energy), a week-by-week schedule, a CPU-hour
   budget table, a prerequisites list keyed to specific chapters, a
   named pitfalls list, and a deliverables checklist. This is more
   structure than most MSc thesis briefs currently provide in the
   author's home institution, never mind in a free textbook.

## Concerns (3)

1. **Production defects in the front matter that any commissioning
   editor would flag on first read.** The root `README.md` contains
   several typos and one ungrammatical sentence in the project's
   *opening paragraph*: line 18, "The ideas starts from a normal
   undergraduate student, who feels it hard and mist importantly, no
   guide… detailed contents are filled under human supervision using
   Claude. Minor mistakes might happen for thhe relevant graphs and
   plots…". A casual reader who reaches "mist importantly" and "thhe"
   in the *first paragraph after the badges* will, fairly or not,
   downgrade their expectation of the rest. This is a thirty-second
   fix and must happen before any v1.0 announcement.

2. **Internal numerical inconsistency in the project durations.** The
   root `README.md` lines 46–50 advertise Projects 1–5 as ~4 / ~6 / ~6
   / ~8 / ~10 weeks. The actual project READMEs report 6 / 8 / 8 / 8–10
   / 6 weeks respectively. The comparative table in
   `docs/ch14-capstone/07-the-five-projects.md` line 293 reports 6 / 8
   / 8 / 10 / 6. No two of these three sources agree. For a textbook
   whose Chapter 14 explicitly teaches students that "quoting MLIP
   accuracy on the training set" (line 53) is a hallmark of careless
   research, internal numerical drift of this magnitude in the load-
   bearing planning tables is awkward.

3. **A handful of unsupported scope and audience claims in the
   marketing surface.** `docs/index.md` line 9 promises a "zero-foundation"
   text from "high-school algebra all the way to foundation models",
   and the README's `mkdocs.yml` site-description (line 2) echoes
   this; but `README.md` line 24 sets the actual prerequisite floor at
   A-level mathematics, single-variable calculus, and "you know what a
   wavefunction is even if you have never solved Schrödinger's
   equation by yourself", and the Tier-0 self-test (`docs/prerequisites-checker.md`)
   contains a question (Q5) on unitary matrices that a literal "high-
   school algebra" reader will not answer. The "zero-foundation" claim
   is a marketing over-shoot. Similarly the homepage table at
   `docs/index.md` lines 35–41 lists "five parts" running to Chapter
   12; but the nav in `mkdocs.yml` (lines 254–272) has six
   sections (adding "VI · Frontier" with Chs 13–14), and `README.md`
   describes fifteen chapters in three tiers reaching Ch 14. The
   reader is asked to hold three slightly different table-of-contents
   in their head before they have read a single section.

## Scores (0-10)

| Dimension | Score | Justification |
|---|---|---|
| Audience clarity | 7 | Three reading paths (Path A/B/C in `docs/learning-path.md` 103–151) and the prerequisites checker are excellent. Loses two points for the "zero-foundation" vs "A-level" contradiction between `docs/index.md` line 9 and `README.md` line 24, and one for the homepage that addresses no fewer than four reader profiles in two paragraphs (`docs/index.md` lines 48–53) without committing to a primary one. |
| Originality / unique contribution | 9 | I cannot name a competitor that goes from Ch 0 "what is a vector" to Ch 12 "fine-tune MACE-MP-0 on a perovskite" in a single coherent voice. Marder skips the ML; Frenkel–Smit predates it; Sholl–Steckel stops at hybrid functionals; the Schmidt et al. and Choudhary et al. reviews cited in `docs/ch12-foundation/04-frontiers.md` lines 447–460 are reviews, not curricula. The book occupies an empty quadrant. |
| Scope completeness | 8 | The arc Ch 0 → Ch 14 is real and the seams are deliberate (Ch 3.5 is explicitly inserted between Ch 3 and Ch 4 to repair the bridge, as called out in `docs/ch03b-solid-state/index.md` line 9). One point lost for `docs/known-limitations.md` §1: Tier 0 linear algebra does not in fact cover Hilbert-space and bra–ket technique required by Chs 4–5. Another lost for the Tier-2 cluster-access gap (Known Limitations §3) — for a self-study learner only Project 3 is genuinely accessible. |
| Pedagogical voice | 8 | The voice is consistent across the chapters I sampled: friendly, opinionated, willing to flag what *not* to read on a first pass (`docs/learning-path.md` 184–191). The §3b.7 worked example on a 5 nm GaAs quantum well (`docs/ch03b-solid-state/07-defects-and-band-engineering.md` 170–183) is the kind of small, concrete, numerical anchor that good textbooks live or die on. Loses two points for occasional drift into self-marketing tone (Ch 14 §"Why we wrote this chapter" 46–61 is honest but borderline preachy). |
| Editorial polish | 5 | The biggest weakness. README typos ("mist importantly", "thhe", "material sicence", a stray space "The structure"), inconsistent project durations (Concern 2), inconsistent ToC counts ("five parts" vs six sections in nav), and a handful of dead-looking cross-references (e.g. `docs/ch14-capstone/index.md` line 7 links to the bare URL `https://github.com/` rather than the projects folder). These are individually small; together they cumulate into a v0.9 feel. |
| Production readiness | 6 | A reader can pick this up *today* and learn from it — the navigation, the search, the MathJax rendering, the in-browser Pyodide runner (`docs/how-to-use.md` 50–62), and the JupyterLite buttons (e.g. `docs/ch00-math/index.md` line 3) all work. The compute story for projects 1, 2, 4, 5 (Known Limitations §3) is the production hole that prevents a clean v1.0 ship to a self-study audience. |
| Overall readiness for v1.0 | 7 | A confident "ship as v0.9.5, fix concerns, ship as v1.0 in six weeks". The structural work is done; what remains is polish and a small list of fact-checks. |

## Specific suggestions for v1.1

1. **Single source of truth for project metadata.** Maintain project
   duration, compute budget, prerequisites, and difficulty in *one*
   YAML or JSON file per project, and render the README, the Ch 14
   comparative table (`docs/ch14-capstone/07-the-five-projects.md`
   line 293), and the root README capstone list (lines 46–50) from
   that file at build time. The current three-way disagreement on
   durations (Concern 2) cannot survive a second motivated reader.

2. **Reconcile the "zero-foundation" promise with the actual A-level
   floor, in the reader's favour.** Either lift `docs/index.md` line 9
   to read "A-level-foundation handbook…" — the honest statement that
   matches `README.md` line 24 — or commission a short pre-Chapter-0
   "What you absolutely need" page that delivers genuine high-school-
   to-A-level remediation (one page on functions, one on matrices, one
   on derivatives). The current gap between the marketing surface and
   the prerequisite checker rewards the wrong reader.

3. **Add laptop-tier variants to projects 1, 2, 4, 5.** Known
   Limitations §3 already concedes this. Make it concrete:
   `docs/projects/01-defect-formation-energy/README.md` could add a
   "laptop-only path" section that caps the supercell at 64 atoms,
   skips the 512-atom point, and substitutes a pre-computed reference
   value for the missing data points. Similar treatments for Projects
   2, 4, 5 would close the largest accessibility gap the book
   currently has, and would let the *Open Educational Resources*
   imprint genuinely fulfil the "open" promise on its title page.

4. **Tighten the front matter to a single voice.** The homepage
   (`docs/index.md`), the root README, the `mkdocs.yml` site
   description, and the "How to use this book" page each address the
   reader in a slightly different voice with a slightly different ToC
   structure ("five parts" vs three tiers vs six nav sections vs
   fifteen chapters). Choose the three-tier framing — it is the book's
   genuine architectural innovation — and propagate it consistently.
   A reader who lands on `docs/index.md` should see the same Tier 0 /
   Tier 1 / Tier 2 vocabulary they will use for the next 200 hours.

5. **Move "Known Limitations" out of the back-of-book and into the
   start-here block, and timestamp it.** `docs/known-limitations.md`
   is already excellent; it is also currently exposed only in the
   "Start Here" nav section, which a linear reader hits before they
   know what any of its three sections are talking about. Either
   gate it ("read this *after* you have finished Tier 1") or, better,
   split it: a short "What this book is honest about" panel at the top
   of `docs/index.md`, plus the detailed three-section version where
   it lives now. Add a "last updated" date — for a 2026 book that
   discusses MACE-MP-0 and MatterGen by name, datedness is a feature
   and the book should flaunt it.

## Editorial Decision

**Minor Revision.** The manuscript is structurally complete, the
pedagogical concept is original and well-executed, and the three-tier
promise is genuinely delivered. The required revisions are confined to:
(i) front-matter copy-editing (Concern 1), (ii) reconciliation of the
three-way disagreement on project durations and the two-way disagreement
on table-of-contents structure (Concerns 2 and 3), and (iii) the five
suggestions above, of which only suggestion 1 is genuinely mechanical
work. None of these require re-writing any chapter. I would expect
v1.0 in six to eight weeks, and I would commission this title for the
imprint's open educational resources line on the strength of the
material I have read.
