# Devil's Advocate Review

*Reviewer brief: find every reason not to recommend this book. Read as a
hostile referee. Today's date: 2026-05-12. Target reader: a Cambridge IB
undergraduate.*

## Strongest counter-argument

The strongest case against recommending the *Materials Simulation
Handbook* to a Cambridge IB student is that **the book sells breadth at
the price of depth in exactly the places where depth is the only thing
that matters**. The IB student already has Mott, Ashcroft–Mermin, and a
shelf of supervisor-curated DFT notes; what they actually need is (a) a
working DFT calculation on their first day of Part III project work, and
(b) one — *one* — fully reproducible foundation-model pipeline that
they could submit as a thesis chapter. This book gives them neither
quickly. It gives them a 16-chapter tour. By the time they have read
through Chapter 0 (maths from numbers up) and Chapter 1 (Miniforge
install — fine, but they already had Anaconda), Chapter 2 (three
diagrams), Chapter 3a + 3b (atoms and solid state — competently covered
in every physics undergraduate course at Cambridge already), and four
chapters of DFT/MD theory, they could have *finished* the official
MACE-MP-0 Colab tutorial and the Quantum ESPRESSO foundation tutorial,
done a vacancy calculation, and read the original PBE paper.

Worse, the book's distinctive material — the foundation-model chapter
(12), the multiscale chapter (13), and the capstone (14) — is the part
that is most freely available elsewhere as well-maintained Colab
notebooks, blog posts, and the Pó-et-al. Matbench-Discovery write-up.
The early chapters, which the book has comparative advantage to write,
are precisely those the IB student has the *least* need of. A
Coursera/edX DFT course plus the MACE-MP-0 tutorial plus *Modern
Quantum Chemistry* (Szabo–Ostlund) and *Electronic Structure* (Martin)
would get the same student further faster, at zero cost, with stronger
authority. The handbook's main value is editorial cohesion — but that
cohesion is undercut by the unkept cross-references and self-admitted
gaps listed below.

## CRITICAL issues (must fix before release)

| Where | What | Why critical | Fix |
|---|---|---|---|
| `docs/ch00-math/01-numbers.md:294` | Claims "In **Chapter 1** we will see DFT scaling as $O(N^3)$ in the number of bands and Monte Carlo error scaling as $O(1/\sqrt{M})$". Chapter 1 is the Python/install chapter; it contains no DFT scaling discussion and no Monte Carlo content. DFT scaling is in Ch 5.3 (Kohn–Sham). | The book's opening chapter promises content that does not exist where it points to. A first-time reader following the promise finds Miniforge installation instructions. The `.cross-check-report.md` Round 3 explicitly claimed "Unkept promises: 0" — that claim is false. | Re-point to Chapter 5 (§5.3 for DFT scaling, Ch 7 / Ch 8 for Monte Carlo error). |
| `docs/ch00-math/01-numbers.md:288` | "We will use induction sparingly — once to establish the binomial theorem, **once when discussing recursion in Chapter 1**". `grep -i recursion` in `ch01-python/*.md` returns zero hits. | Same as above — a promise unmet at the very first cross-reference. Undermines reader's trust on page 5. | Either add a brief recursion subsection to Ch 1 (e.g. when motivating numpy vectorisation against Python recursion) or remove the promise. |
| `docs/ch14-capstone/07-the-five-projects.md:25,77,129,183,238` | All five "Project folder" links point to bare `https://github.com/` — i.e. nothing. Quotes: `**Project folder:** [\`projects/01-silicon-defects/README.md\`](https://github.com/)` etc. | The capstone chapter is the book's killer feature ("doable in 4–10 weeks, produce something publishable"). Every link to the actual project material is dead in v1.0. | Replace with relative MkDocs links `../projects/01-defect-formation-energy/README.md` (note the folder names also differ — `01-silicon-defects` vs `01-defect-formation-energy`, `02-copper-melting` vs `02-melting-point-mlip`, etc.). |
| `docs/ch14-capstone/07-the-five-projects.md:25 ff.` | The displayed *path* in each link's link text disagrees with the actual on-disk folder names. E.g. the chapter calls the folder `projects/01-silicon-defects/` whereas disk has `docs/projects/01-defect-formation-energy/`. All five are wrong. | Path-shaped strings in a textbook are commitments. Even if the hyperlink is later fixed, copy-paste of the visible path text fails. | Rename folders to match the chapter text, or update the chapter text to match disk. The latter is cheaper. |
| `README.md:18` | "The ideas starts from a normal undergraduate student, who feels it hard and mist importantly, no guide into the computational simulation of materials. … Minor mistakes might happen for thhe relevant graphs and plots". Three typos in one paragraph of the public-facing README: "ideas starts", "mist importantly", "thhe". | This is the first English the reader sees. It is the single line judging the book's editorial standards. The student is being asked to trust a multi-chapter physics textbook produced by an author who has not proofread their own README. | Rewrite paragraph; run a spell checker on README, `docs/index.md`, and `learning-path.md`. |
| `docs/known-limitations.md:20-31` | The author admits in writing that "Chapter 9 introduces machine-learning interatomic potentials enthusiastically … but contains relatively little material on the documented failure modes of MLIPs — instabilities …, the long-range / charge-transfer limitations …, the silent-failure problem". | A textbook chapter on MLIPs is selling MLIPs without honest failure-mode coverage *and the author knows it*. The honesty in `known-limitations.md` does not absolve the chapter; it indicts it. The frontier honesty in Ch 12.4 does not reach a reader who only reads Ch 9. | Write the promised §9.7 *before* declaring v1.0. Include the Stocker et al. high-T stability tests, the long-range failure on polar phonons, the silent-extrapolation example. |
| `docs/known-limitations.md:33-43` | The author admits the capstone projects "assume cluster access that not every reader will have. … a reader without an institutional allocation effectively has only Project 3 (which can be done on a small GPU) accessible." | This contradicts the book's marketed audience ("self-study learners"; learning-path's "Path C — Project-driven"). Four of five flagship projects are unavailable to the unfunded reader. | Add a laptop-tier variant per project, *before* claiming v1.0 — not in a future revision. |
| `docs/ch12-foundation/02-mace-mp0.md:323-339` | The headline universal-MLIP comparison table is captioned "late 2025" and quotes "Orb-v2 — Top of Matbench-Discovery as of 2025". The publication date implied throughout is mid-2026, and the leaderboard has demonstrably moved (eqV2-OMat large, $F_1 = 0.83$, is listed at line 385 but the body text at line 354–360 still says "MACE-MP-0 medium remains, as of mid-2026, the most balanced choice"). MatterSim (Microsoft, 2024) is mentioned in a single drive-by sentence in Ch 2 with no entry in the comparison table at all. | A book about foundation models that does not include the FAIR universal model line (UMA, eqV2 family treated only partially), MatterSim, or Orb-v3 — and prints stale leaderboard numbers as if current — is anti-current the day it ships. | (a) Either drop dated leaderboard numbers and refer the reader to the live Matbench-Discovery site; (b) add a one-row "What's missing here" disclaimer; (c) update the zoo to include MatterSim and any v3 releases extant on 2026-05-12. |

## MAJOR issues

| Where | What | Why major | Fix |
|---|---|---|---|
| `docs/ch00-math/04-complex-fourier.md:164, 435` and chapter-naming throughout | The reader is told "Bloch theorem we will see in Chapter 3.5" but the directory is `ch03b-solid-state`, the section header is `# 3b.1 — Bloch's Theorem`, the MkDocs nav says `3.5.1 Bloch's theorem`, and the URL slug is `ch03b-solid-state/01-bloch-theorem/`. Three different numbering systems for the same chapter. | A textbook with three names for one chapter signals editorial instability. The cross-check report claims this was fixed; the substance is unchanged — only the prose was updated, not the internal section numbering. | Pick one canonical scheme (recommend keeping "3.5" externally as in nav and updating internal `## 3b.1.x` headings to `## 3.5.x`). |
| `docs/ch00-math/01-numbers.md:24` | Says "We will return to this in Chapter 1" for IEEE-754 floating point. Ch 1 is the install + numpy chapter; it does *not* return to IEEE-754 in any explicit treatment. | Reinforces the unkept-promise pattern at the book's start. | Either redirect to a specific Ch 1 section that actually treats IEEE-754, or strip the promise. |
| `docs/ch00-math/03-calculus.md:280, 596` | Promises Ch 1 will contain "Better rules (trapezoidal, Simpson, Gauss quadrature)" and "Numerical differentiation … reappears in Chapter 1 as a debugging tool". `grep -i 'trapezoid\|simpson\|gauss quadrature'` in `ch01-python/*.md`: zero hits. | More Ch0→Ch1 promises that Ch1 does not deliver. | Add the quadrature material to Ch 1 (it would *belong* there) or rewrite the promises. |
| Ch 9 (MLIP) chapter as a whole | Reads as an enthusiastic tour of descriptors → BP/GAP → equivariant → MACE training, with the "What MLIPs are not" section (9.1.5) running 12 lines. Compare with Ch 5.6 (DFT failures), which is a stand-alone 280-line section. | The book is structurally honest about DFT failures and structurally optimistic about MLIPs. The asymmetry trains the reader to under-discount the latter. | Promote 9.1.5 to a full §9.7 with the depth of 5.6 and the worked stability/extrapolation case studies the author flagged in `known-limitations.md`. |
| `docs/ch12-foundation/02-mace-mp0.md:62-67` | Treats MPtrj's "equilibrium-heavy" and "chemistry-imbalanced" biases as a known issue then waves them off ("Section 12.4 will discuss the OMat24 follow-on dataset, which was designed in part to address them") — but never quantifies the bias in residual force error on, say, high-T MD configurations. | This is the single most consequential fact about every model in the chapter. The book under-quantifies it. | One bar chart: zero-shot force MAE on (a) equilibrium MPtrj, (b) WBM substitutional, (c) 1500 K MD configurations, for MACE-MP-0 medium. Without this figure, the chapter's caveats are decorative. |
| `docs/ch14-capstone/07-the-five-projects.md:11-21` | "Tractable: each fits in 6–10 weeks of part-time work, or 3–5 weeks full-time". Cross-check with project READMEs: Project 4 (MLIP from scratch) and Project 5 (Bayesian catalyst search) both compress months of skill acquisition into the same window. | The "tractable" claim is cherry-picked from the easiest project. An IB student attempting Project 4 (implement a BP potential end-to-end) in 6 weeks is being misled. | Either widen the time band per project, or add a "this project is hard" warning to Projects 4 and 5. |
| `docs/ch10-gnn/03-cgcnn.md:582-610` | Uses UMAP on a CGCNN element embedding to claim "alkalis cluster together" etc. This is the canonical cherry-picked figure of GNN papers. | The same plot can be produced from one-hot encodings + a random projection. Without a baseline (random embedding, or a non-trained CGCNN), it is not evidence of anything. | Add a control: same UMAP plot from an untrained CGCNN. If the structure persists, say so; if it changes, that *is* the point. |
| `docs/ch12-foundation/01-paradigm.md:206-289` | The "MACE-MP-0 on lithium battery materials" subsection picks lithium batteries — a domain where MPtrj has rich training coverage. Reads as a sales pitch. | The book's own §12.4 admits OOD chemistry is a real problem. Yet the empirical case study chooses the friendliest possible chemistry. | Show one in-distribution success (current) and one out-of-distribution failure (e.g., a high-pressure hydride, an organic crystal, or an actinide oxide) side by side. |

## MINOR issues

| Where | What | Why minor | Fix |
|---|---|---|---|
| `README.md:114` | "If you cite a specific chapter, please give the chapter title and section number alongside the repository URL, e.g. *'Materials Simulation Handbook, Ch 9 §5 (equivariant networks).'*" — but the internal section numbering for Ch 9 (e.g. `## 9.1.5`) is multi-level; the example "§5" is ambiguous. | Citing this book is awkward. | Use unambiguous example, e.g. `Ch 9 §9.5 (equivariant networks)`. |
| `docs/ch01-python/01-install.md:18` | "Anaconda's commercial terms have also restricted its use in organisations with more than 200 employees." Correct at time of writing but worth a date stamp; the policy has moved before. | Stale policy claims age badly. | Add "(as of 2024)" so a 2027 reader knows to verify. |
| `docs/ch12-foundation/02-mace-mp0.md:80` | Says zero-shot MACE-MP-0 is "often within a factor of two of a hand-fitted potential trained on 1,000 system-specific configurations" — no citation. | Empirically plausible but uncited. | Cite Batatia et al. 2024 supplementary, or the Pó et al. Matbench-Discovery paper. |
| `docs/ch07-md/index.md` (and ten others) | Cross-check report admits "ch07-md/index.md — link to non-existent `../appendix/units.md`; pointed to `../appendix/A-math-reference.md`." Many such silent redirects exist; readers expecting "units" find "math reference". | Document hygiene; not blocking. | Either add a §units to Appendix A or rename the link text. |
| `docs/ch11-active/02-gp.md` (983 lines) | Single section is 983 lines. That is longer than the entire Ch 13 multiscale chapter. | Pedagogical pacing imbalance. | Split into §11.2a (theory of GPs) and §11.2b (computational details). |
| `docs/ch14-capstone/01-how-research-works.md:179-220` | The anonymised cautionary tales ("the student who never converged k-points", "the student who over-claimed MLIP transferability") read as composite anecdotes presented as data. | Pedagogically reasonable but methodologically weak in a book that elsewhere demands citations. | Mark explicitly as illustrative, or replace with citable case studies (e.g., the published retractions in DFT functional benchmarks). |
| `docs/ch12-foundation/02-mace-mp0.md` (whole section) | Mentions `mace-torch` PyPI package and pre-trained weights "downloaded on first use" — no checksum, no pinned version. | A 2027 reader hits a different model under the same name. | Pin: e.g., `mace-torch==0.3.x`, and record the model hash. |
| `docs/ch02-foundations/exercises.md:67` | "(f) **Not yet.** Synthesis prediction — including precursor choice, temperature schedule, and yield — remains beyond routine simulation." — this is true *and* in tension with the optimistic A-Lab framing in Ch 12. | Internal tone mismatch. | Reconcile: A-Lab gives autonomous *attempts*, not *prediction*; sharpen wording. |

## Ignored alternative paths

A book like this could have been radically more useful by *narrowing*. Three uninvestigated paths:

1. **"DFT in one weekend"**: a 60-page minibook that takes the reader from zero install to a converged vacancy energy in two days. The handbook spends Chapters 0–6 to get there. The vast majority of IB students need exactly the minibook; they will not read 16 chapters before running anything.
2. **"MACE-MP-0 only"**: skip GAP, skip Behler–Parrinello as a historical curiosity, treat MACE-MP-0 as the unit and teach fine-tuning, active learning, validation, and uncertainty quantification *on top of it*. Ch 9 currently builds the historical scaffold and only reaches MACE at §9.5 — a pedagogical inversion if the student's actual research artifact will be an MLIP foundation model.
3. **"Reproduce one published paper"**: an entire textbook structured around faithfully reproducing one foundational paper (e.g., the original PBE paper, or the CGCNN paper) end-to-end. This is the format used by *Deep Learning from Scratch*; it forces honesty because every claim must compile.

The handbook chose breadth. Breadth has costs that the book does not pay for: editorial coherence across 16 chapters is hard, and the patches show.

## Missing stakeholder perspectives

- **Self-study learners without cluster access.** Author admits Project 1, 2, 4, 5 require resources most don't have. So 80% of capstones serve only the funded reader.
- **Industry users.** The book trains the reader on Materials Project + VASP-style PBE+U workflows. Industry computational materials groups (e.g., Toyota Research Institute, Microsoft Quantum, BASF) increasingly use proprietary MLIPs and licensed codes the book does not mention.
- **Software engineers writing the next MLIP.** No discussion of inference deployment (ONNX export, TorchScript pinning, mixed-precision inference, batched evaluation throughput on consumer GPUs). The book trains *users* of MLIPs, not *builders*.
- **Experimentalists asked to interpret a computational paper.** A 50-page "computational results for experimentalists" chapter would be more impactful than half the existing chapters. None exists.
- **Educators.** No instructor manual, no exam questions with rubrics, no slide decks. The book is harder to teach *from* than to learn *from*.

## So-what test result

**Half pass.**

What the handbook *uniquely* offers, that I cannot find freely elsewhere in one place:

- The Ch 5.6 honest-failures-of-DFT section is genuinely strong and not easy to assemble from Martin / Burke / Cohen alone.
- The Ch 14 "how research works" + "choosing a problem" framework is decent thesis-mentor material and rare in textbooks.
- The Ch 11 + Ch 12 combined treatment of active learning *layered on top of* universal MLIPs is current and not in the established texts.

What the handbook offers that is **freely available elsewhere** in equal or better form:

- Ch 0 (maths): every undergraduate has these notes; Strang and 3Blue1Brown are better.
- Ch 1 (Python): Software Carpentry, the SciPy lecture notes.
- Ch 3, 3b (atoms/solid-state): Ashcroft–Mermin, Kittel, Simon (*Solid State Basics*).
- Ch 4 (quantum): Griffiths.
- Ch 5 (DFT): Martin *Electronic Structure*, Sholl & Steckel *DFT: a Practical Introduction*, Burke's online ABC of DFT.
- Ch 6 (running DFT): the QE tutorials are excellent and current.
- Ch 7 (MD): Frenkel & Smit, Tuckerman.
- Ch 9, 10 (MLIPs/GNNs): Csányi-group materials, the MACE tutorial, the DimeNet/SchNet papers.
- Ch 12 (foundation MLIPs): the Matbench-Discovery paper and the MACE-MP-0 Colab are state-of-the-art and updated.

The honest verdict: a Cambridge IB student would benefit most from **Chapters 5, 11, 12, 14 only** — about 50,000 words — bound separately as a "foundation-model project handbook". The other 100,000 words duplicate material that is more authoritative elsewhere and that the IB student has already met. As a single 16-chapter recommendation, the book competes against its own better halves and loses.

The recommendation: read it for Ch 5.6, Ch 11, Ch 12 (with the caveats above), and Ch 14. Skip the rest. Use Coursera + the MACE-MP-0 tutorial + Martin's *Electronic Structure* for the prerequisites.
