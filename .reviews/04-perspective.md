# Peer Review 3 — Pedagogical Perspective

**Reviewer.** Physics-education / learning-science perspective. I read the
handbook as a sceptical PER researcher would: paying attention to how a
new reader would actually move through the material, where their working
memory would saturate, and whether the exercises and self-tests do real
diagnostic work. I tried, throughout, to imagine a strong Cambridge IB
Natural-Sciences student in the first long vacation — competent at
calculus and linear algebra, fluent in Python, but with no DFT, almost
no statistical mechanics, and only the vaguest contact with reciprocal
space.

**Scope of reading (≈ 35 min).**

- `docs/index.md`, `docs/how-to-use.md`, `docs/learning-path.md`,
  `docs/prerequisites-checker.md`.
- Pedagogical transitions: end of Ch 0 → Ch 4 (via
  `docs/ch00-math/05-probability.md`, `docs/ch04-quantum/01-why.md`,
  `docs/ch04-quantum/index.md`); end of Ch 3.5 → Ch 4 (via
  `docs/ch03b-solid-state/index.md`, `docs/ch03b-solid-state/01-bloch-theorem.md`);
  end of Ch 5 → Ch 6 (via `docs/ch05-dft/06-failures.md`,
  `docs/ch06-running-dft/index.md`, `docs/ch06-running-dft/01-pseudo-basis.md`).
- Exercise sets: `docs/ch00-math/exercises.md`,
  `docs/ch05-dft/exercises.md`, `docs/ch11-active/exercises.md`.
- Project template: `docs/projects/01-defect-formation-energy/README.md`.
- Random spot reads for cognitive-load assessment:
  `docs/ch05-dft/03-kohn-sham.md`, `docs/ch11-active/02-gp.md`,
  `docs/ch00-math/index.md`, `docs/ch00-math/03-calculus.md`.

---

## Scores (0–10)

| Dimension | Score | One-line basis |
|---|---|---|
| Difficulty progression smoothness | 6.5 | Within-chapter pacing is excellent; **between** Ch 3.5–4 and Ch 5–6 there are real cliffs, and Path A in `learning-path.md` orders Ch 3.5 before Ch 4 in direct contradiction to Ch 3.5's own stated prerequisites. |
| Misconception coverage | 7.5 | The "density is not probability" box (`05-probability.md`, §"Density is not probability"), the "KS orbitals are not the real wavefunction" insistence (`03-kohn-sham.md`, end of §5.3.1), and the explicit photoelectric misconception list (`01-why.md`, §4.1.2) are model PER moves. Coverage is uneven though: nothing flags the "average over k means k-grid" confusion in Ch 6.1, nothing flags the universal "GP uncertainty is calibrated" misconception in Ch 11.2. |
| Cognitive load per page | 6 | Several pages cross the seven-items-in-working-memory line in a single paragraph. Ch 6.1.2 dumps NC / USPP / PAW with new equations, new operators ($\hat{S}$, $q^{(a)}_{ij}$), new notation, all before the reader has run a single calculation. Ch 11.2.1 opens with a definition involving "any finite collection… joint Gaussian" — formally correct, pedagogically opaque to anyone meeting GPs for the first time. |
| Multi-modal reinforcement (text + code + figure) | 7 | The "▶ Run" promise in `how-to-use.md` is genuinely innovative if delivered; the code snippets in Ch 0 and Ch 5 exercises tie text to executable verification well. Figures are mostly absent from the prose I sampled (the GP figure in `ch11/02-gp.md` is a positive exception). Analogies are present but inconsistent — the "swarm of bees" analogy for KS (`03-kohn-sham.md`) is brilliant; nothing comparable appears for self-interaction error or Bloch's theorem's "snake-on-corrugation" depends entirely on the reader having a mental image. |
| Self-assessment quality | 6.5 | Exercise solutions are excellent **answers** but mediocre **diagnostics**: they tell the student the right answer, but rarely tell them what their wrong path implied about which misconception was active. The prerequisites checker (`prerequisites-checker.md`) attempts remediation pointers per question — this is the right idea, but only the multiple-choice questions are realistically self-gradable; the short-answer ones (Q3, Q9, Q13, Q16) ask the student to grade their own prose, which IB-level students reliably over-rate. |
| Vocabulary/jargon management | 7 | The "What you need / What you do not need" pattern at the top of each chapter index is genuinely good practice. Bold-on-first-use is mostly observed. However, jargon does occasionally appear before definition: `ch04-quantum/index.md` uses "Hermitian", "eigenvalue problem", and "Slater determinant" in the roadmap before §4.2 defines any of them; the chapter goal box on the same page deploys $\hat{H}\Psi(\mathbf{r}_1, \ldots, \mathbf{r}_N) = E\Psi$ before any reader of Path A has seen what a multi-electron $\hat{H}$ looks like. |
| Accessibility for a true zero-foundation reader | 5.5 | The book *says* "zero foundation, high-school algebra" (`docs/index.md` line 9; `ch00-math/index.md` line 6). It mostly delivers on this in Ch 0. But the operational reality from Ch 3.5 onwards is that the implicit prerequisite is closer to second-year university physics — degenerate perturbation theory, the dynamical matrix, Wirtinger calculus (mentioned offhand in `ch00-math/exercises.md` Exercise 8 solution) are all assumed once needed. The promise and the delivery diverge. |

**Overall pedagogical score: 6.6 / 10.** A genuinely excellent project for
its intended Path B / Path C audience (a finalist with solid maths). It
overpromises for the true beginner and underprovides the diagnostic
scaffolding that would let such a reader recover from a wrong turn.

---

## Three pedagogical strengths

### 1. "Why this step?" boxes are PER-grade scaffolding

These boxes appear throughout Ch 4 and Ch 5 (e.g.
`ch04-quantum/01-why.md` lines 69–70; `ch05-dft/03-kohn-sham.md` lines
18–19, 65–66). They do something rare: they explicitly model the
**cognitive move** the author is making, not just the algebraic step.
The "Why translations commute" note in `ch03b-solid-state/01-bloch-theorem.md`
(lines 49–50) is a textbook example — it pre-empts the standard
misconception that abelian-vs-non-abelian is a technicality, and shows
the student that the entire structure of Bloch labels is downstream of
this fact. Wieman would approve.

### 2. The roadmap-then-content-then-recap structure of each chapter

Every chapter `index.md` I read (`ch00-math/index.md`,
`ch04-quantum/index.md`, `ch05-dft/index.md`, `ch06-running-dft/index.md`)
follows the same three-part pattern: motivation, section-by-section
roadmap, and an explicit "What you need / What you do not need" box. The
practical advice on `learning-path.md` (line 201) — "Re-read the index
pages. Every chapter's `index.md` is a one-page map of the chapter.
Read it before starting, and again after finishing; the gap between
the two is what you have learnt." — is exactly the metacognitive prompt
that PER literature shows produces real learning gains. This deserves
to stay.

### 3. The "honest assessment" / "when does this method fail" sections

`ch05-dft/06-failures.md` is one of the strongest pieces of writing in
the handbook. It does what almost no DFT textbook does: tells the
student, in advance, which classes of problem the tool they are about
to learn will get systematically wrong. The escalation table
(LDA → SCAN → HSE06 → GW → CCSD(T)) is a fantastic learning artefact
because it gives the student a *map* of the methodology landscape, not
just a recipe for one square of it. The summary block "what to
remember in 3 months" anticipates the forgetting curve and gives the
reader a deliberate, narrow set of takeaways to consolidate. This is
pedagogically sophisticated; many a senior textbook does worse.

---

## Three walls a real IB student will hit

### Wall A — The Ch 3.5 / Ch 4 ordering contradiction

`docs/learning-path.md` Path A orders the chapters thus (lines 109–118):
"Weeks 3–4: Ch 2, Ch 3, Ch 3.5. Weeks 5–7: Ch 4, Ch 5." So a Path A
reader meets Ch 3.5 *before* Ch 4.

But `docs/ch03b-solid-state/index.md` line 9 explicitly says: "We are
deliberately inserting this material between Chapter 3 (descriptive
crystallography) and Chapter 4 (quantum mechanics in a box) in spirit,
but you should read it after Chapter 4". The same chapter's "What you
need" (line 47) lists, from Chapter 4: "the time-independent Schrödinger
equation, eigenvalues of Hermitian operators, bra-ket notation, the
harmonic oscillator".

A real Path-A reader who follows the calendar will arrive at
`ch03b-solid-state/01-bloch-theorem.md` and immediately encounter
$\hat{H} = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf{r})$ (eq. 3b.1.3)
applied to "eigenstates of $\hat{H}$" — a phrase the Path-A reader has
never been formally introduced to. The cross-reference matrix on
`learning-path.md` (line 167) confirms that Ch 4 depends on Ch 0
only, and Ch 5 depends on Ch 3.5 *and* Ch 4 — implying Ch 3.5 and Ch 4
are siblings, neither dependent on the other. That is wrong: 3.5
genuinely needs 4.

For a confident finalist this is a 20-minute confusion that gets
resolved by trusting the chapter's own preface over the path. For a
true zero-foundation reader on Path A, it is a wall: they will read
Bloch's theorem on a Hamiltonian whose definition the book formally
postpones to next week.

### Wall B — The Ch 5 → Ch 6 detail cliff

Ch 5 is mostly theoretical (functionals, HK, KS, SCF, $E_{xc}$ on
Jacob's ladder). Its closing in `ch05-dft/06-failures.md` is
beautifully reflective. The student then opens
`ch06-running-dft/01-pseudo-basis.md` and within 1500 words meets:
norm-conserving vs ultrasoft vs PAW pseudopotentials, the augmentation
operator $\hat{S}$ (eq. on line 59), augmentation integrals
$q^{(a)}_{ij}$, semi-core states, projector indices $i,j$, the SSSP
library, cutoff energies in Ry, and `pw.x` input files (in §6.2).

This is a classic cognitive-load failure: too many novel concepts,
each justified only by appeal to its place in the practical pipeline.
The chapter index (`ch06-running-dft/index.md` line 13) says "Two
numbers — the cutoff $E_\mathrm{cut}$ and the k-point grid — control
the entire calculation." That sentence is the pedagogical anchor; but
the chapter then *front-loads* every other detail before the student
has held those two knobs in their hand. The right move pedagogically
is to defer NC / USPP / PAW until §6.3 (after the student has run one
calculation), but it is currently §6.1.2.

A real Cambridge IB student opening §6.1 cold will produce a page of
notes, none of which they will be able to use until they have the
muscle memory of `pw.x silicon.in > silicon.out`. The chapter's own
table of contents puts that experience at §6.2 — the chapter would
work better with §6.1 and §6.2 swapped or aggressively pruned.

### Wall C — The "this is the hardest section" problem in Ch 11.2

`ch11-active/02-gp.md` opens with a one-line definition
("any finite collection… joint Gaussian", lines 41–48) that is
mathematically the right definition, but pedagogically inert: there is
no path from the definition to an intuition for what a sample from a
GP looks like. The chapter does include a `!!! note "Why this is the
hardest section in Chapter 11"` (line 23) — which is honest, and the
"Intuition: a GP as an infinite-dimensional Gaussian" box (line 71)
helps — but the *order* is still: formal definition first, intuition
second. PER work on threshold concepts (Meyer & Land) is unambiguous
that intuition has to come first; the formal object is what the
intuition resolves into, not the seed from which it grows.

A reader who survived Ch 4–6 and is now arriving in Ch 11 will hit
this section, recognise they "don't get it", look for a remediation
pointer, and find only "rehearse the linear algebra in §11.2.4". That
is not enough scaffold. The exercises (Exercise 11.4, "RBF kernel and
length scale" — `ch11-active/exercises.md` line 113) would actually be
a much better *opening* than a closing for the section.

---

## Five specific pedagogical improvements

### 1. Fix the Ch 3.5 / Ch 4 ordering and dependency declarations

Three pieces of text are mutually inconsistent:

- `learning-path.md` line 111 ("Weeks 3–4: Ch 2, Ch 3, Ch 3.5. Weeks 5–7: Ch 4, Ch 5.").
- `learning-path.md` cross-reference table line 167 (Ch 5 depends on both Ch 3.5 and Ch 4; Ch 3.5 is listed as depending only on Ch 0 and Ch 3, **not on Ch 4**).
- `ch03b-solid-state/index.md` lines 9 and 47 (you should read 3.5
  after Ch 4; "what you need" includes time-independent Schrödinger
  equation, Hermitian eigenvalues, bra-ket, harmonic oscillator).

Pick one. The pedagogically right answer is: Ch 4 before Ch 3.5.
Reorder Path A to "Weeks 3: Ch 2, Ch 3. Weeks 4–6: Ch 4. Week 7: Ch 3.5.
Weeks 8–9: Ch 5." Then update the cross-reference matrix to list Ch 4
as a hard prerequisite for Ch 3.5, and renumber Ch 3.5 → Ch 4b (or
keep the "3.5" name as historical, but make its position in the linear
path unambiguous).

### 2. Add a "what would going wrong here look like?" block to every
exercise solution

Currently every solution gives the right answer. Add, in collapsed
form below the solution, an explicit "common wrong answers and what
they reveal":

- For `ch00-math/exercises.md` Ex 3 (dot product, angle), the common
  wrong answer is to compute the angle as $\arctan$ rather than
  $\arccos$ — diagnostic for confusion of $\tan\theta$ and $\cos\theta$
  in 2D.
- For `ch05-dft/exercises.md` Ex 5.1 (HK I), common wrong answers
  forget the non-degeneracy hypothesis or use the wrong direction of
  the variational inequality — both are PhD-thesis-grade-mistakes that
  the student should be warned about.
- For `ch11-active/exercises.md` Ex 11.3 (UCB exploration), the
  common wrong intuition is that $\kappa \to 0$ "kills exploration" —
  the solution should explicitly say "no, it does not; the exploration
  bonus has the same *sign*, only smaller *magnitude*, and this is why
  you cannot get pure greedy behaviour from UCB by lowering $\kappa$
  alone".

This is the PER-standard "elicit-confront-resolve" pattern, and it is
exactly what turns an answer sheet into a diagnostic.

### 3. Move the Ch 6.1 pseudopotential taxonomy after §6.2

`ch06-running-dft/01-pseudo-basis.md` currently dumps NC / USPP / PAW
in §6.1.2 before the student has run any calculation. The cognitive
load argument is overwhelming: NC vs USPP vs PAW is **only** a
meaningful distinction for someone who has felt the pain of an
unconverged cutoff. Restructure as:

- §6.1: "The frozen core: why we replace 10 of silicon's 14
  electrons by an effective potential." Single picture, one
  pseudopotential file path, no taxonomy.
- §6.2: First calculation. `pw.x`, the input file, the output, the
  total energy.
- §6.3: Now that you have a number, what would happen if you changed
  the pseudo? Here is the NC/USPP/PAW taxonomy — and here is what each
  change does to your number.

This sequencing matches Sanjoy Mahajan's "concrete-before-abstract"
heuristic and respects the reader's working-memory budget.

### 4. Add real diagnostic scoring to the prerequisites checker

`prerequisites-checker.md` has the right idea (per-question remediation
pointers, scoring rubric, sub-category triage at the end). What it
lacks is:

- **Calibration on short-answer questions.** Q3, Q9, Q13, Q16 ask the
  reader to grade their own prose. IB students reliably over-rate
  themselves. Either replace these with multiple-choice (e.g. for
  Q9, offer four candidate one-line Schrödinger equations and ask
  which is correct) or give a rubric ("did your answer contain
  $\hat{H}$, $\psi$, $E$, and the eigenvalue structure? If any of
  these is missing, mark wrong.").
- **An explicit warning about the Q14 vs Q19 asymmetry.** A reader
  scoring 12 by getting all the Python right and missing four QM
  questions is in a fundamentally different place from a reader
  scoring 12 by getting all the QM right and missing four Python ones.
  The end-of-page sub-category triage (lines 320–325) is good but
  buried; promote it to a prominent decision tree, not a paragraph.
- **A re-test pointer with a date stamp.** "Re-take this test after
  every Tier completes" (line 327) is good advice that no reader will
  follow without a prompt; embed a self-set re-test reminder in the
  rubric.

### 5. Slot in three "pause-and-recall" prompts per chapter

The PER literature on retrieval practice is unambiguous: forcing a
student to recall before they re-read multiplies retention. The
handbook currently relies on the exercise set at the end of each
chapter to serve this purpose. That is too sparse and too late.

Add, at the natural section breaks (typically after every 2000–3000
words), a small `!!! question "Pause and recall"` box with one or two
questions whose answers can be checked in 30 seconds against a
single-line answer that follows. Examples:

- After Ch 4.1 (the historical motivation), before Ch 4.2: "Without
  looking back: name three experiments that classical physics fails
  to explain, and the order-of-magnitude discrepancy in each."
- After Ch 5.3 (KS construction): "What is the *only* physical
  content of a Kohn–Sham orbital?" (Answer: it contributes
  $|\phi_i|^2$ to the density; it is not, by itself, a real
  single-particle wavefunction.)
- After Ch 6.2 (first calculation): "List the four numbers in your
  `pw.x` input that, if changed by 10%, would change your total
  energy by more than 1 meV/atom."

These cost the author thirty minutes per chapter and raise retention
by a measurable fraction. They also create explicit checkpoints where
a reader can decide "I'm not ready — go back" without losing face,
which is what zero-foundation readers most need.

---

## Notes on specific files

- `docs/index.md` line 13 ("You do **not** need a physics or chemistry
  degree to start"): I would soften this to "You do not need a degree
  in physics or chemistry, but a comfort with university-level
  calculus and a genuine willingness to invest 200 hours will both be
  required by Chapter 5." Honesty here protects the reader from
  abandoning the book at Ch 6 thinking they have failed.
- `docs/how-to-use.md` line 56 (in-browser execution and Colab
  pointer): this is excellent, but assumes a reader who can install
  WSL2 (line 32 of `ch06-running-dft/index.md`). For the "true
  beginner", explicitly flag the Colab path as the default for Ch 6+
  and only mention local installation as an enthusiast option.
- `docs/projects/01-defect-formation-energy/README.md`: the pitfalls
  list (lines 122–145) is gold-standard pedagogy — it pre-empts the
  five common errors before the student makes them. Promote this
  pattern to every chapter, not just the project templates.

---

## Concluding remark

This handbook is, at the level of individual sections, one of the best
pedagogical artefacts I have read in computational materials science.
"Why this step?" boxes, honest "where this fails" sections,
explicit roadmaps, and runnable code in the browser put it ahead of
every textbook I would otherwise recommend. The weaknesses are
structural rather than local: an over-optimistic ordering for the
zero-foundation path, an over-loaded first section in Ch 6, and
insufficient diagnostic feedback in the exercises and the
prerequisites checker. All five of the improvements above are
mechanical to implement and would, together, lift the book from
"excellent if you survive it" to "excellent for the audience it
actually claims".
