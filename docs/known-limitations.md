# Known Limitations (v1.0)

This page is the honest companion to the rest of the handbook. The book is
intentionally broad, and breadth costs depth. Below are the three structural
critiques we would expect a careful journal referee to raise at v1.0. They
are flagged here for the reader's awareness and for future revisions.

## 1. The Tier 0 prerequisites do not fully cover the linear algebra demanded later

The prerequisite checklist (Tier 0, Problem 4) assumes a level of comfort
with bra–ket notation, the spectral theorem in finite dimensions, and
projection operators that Chapter 0.2 (Linear Algebra) does not explicitly
teach. A reader who arrives with only A-level / first-year matrix algebra
can pass the self-test but will find Chapter 4 (Quantum Mechanics) and
Chapter 5 (DFT) harder than the difficulty grading implies. Future editions
should either expand Chapter 0.2 to include a half-section on Hilbert spaces
and Dirac notation, or weaken the Tier 0 self-test to flag the gap
explicitly.

## 2. The MLIP chapter under-represents failure modes

Chapter 9 introduces machine-learning interatomic potentials enthusiastically
and reaches MACE training by §6, but contains relatively little material on
the documented failure modes of MLIPs — instabilities under aggressive
sampling, the long-range / charge-transfer limitations of strictly local
message passing, the silent-failure problem when an MLIP extrapolates outside
its training distribution. The frontier section in Chapter 12 acknowledges
these issues, but a reader who studies only Chapter 9 leaves with an
optimistic picture. A future revision should add a §9.7 "When MLIPs Fail"
with concrete worked examples (e.g., the Stocker et al. high-temperature
stability tests).

## 3. The capstone projects assume cluster access that not every reader will have

Projects 1, 2, 4 and 5 specify CPU-hour or GPU-hour budgets in the hundreds.
While Chapter 14 discusses laptop-only alternatives at a high level, the
project READMEs themselves do not offer scaled-down variants — a reader
without an institutional allocation effectively has only Project 3 (which can
be done on a small GPU) accessible. A future revision should add an explicit
"laptop-tier" pathway to each project, with reduced supercell sizes,
single-defect cases, or surrogate-oracle substitutes, so that the project
infrastructure works for self-study learners as well as for graduate students
on a managed cluster.

---

These are issues the editors and authors are aware of and intend to address
in v1.1. Reports of additional structural problems are welcome via the
repository issue tracker.
