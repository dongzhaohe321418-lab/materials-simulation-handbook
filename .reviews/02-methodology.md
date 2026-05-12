# Peer review 2 — Methodology

**Reviewer:** Reviewer 1 (research software engineer, computational physics)
**Date:** 2026-05-12
**Scope:** Numerical methods, code reproducibility, derivation rigour, computational physics best practices.

## Summary

The Materials Simulation Handbook is, by computational-physics standards, an
unusually careful piece of work. The Tier-1 code modules under `code/tier1/`
are extracted from the book text rather than copied separately, the pytest
suite (48 tests + 1 skip) passes cleanly in 2.7 s, CI is wired up on GitHub
Actions, derivations are dense but mathematically correct in every sample I
checked, and the Slater–Koster and LJ-MD figure scripts are written with the
care of someone who has built these things before. The book reads as a piece
of working scientific software rather than a tutorial.

There are, however, several concrete methodological weak spots: a CI matrix
that tests only one Python version on one OS while `environment.yml` declares
support for `python>=3.11,<3.13`; a small but real physics issue in the
LJ-MD `fig_rdf_real.py` (a single density used for all three phases); the GP
implementation in `code/tier1/ch11/gp.py` carries an internal API
inconsistency around input shapes that will eventually surface as a bug; and
the figure-script collection has no pinned package versions or a
"reproduce-all-figures" Make target. None of these are show-stoppers; they
are the items I would expect a v1.1 to clean up before a JOSS submission.

## Strengths

1. **Code is genuinely extracted, not just transcribed.** The Tier-1
   modules under `code/tier1/` (e.g. `ch04/schrodinger_1d.py`,
   `ch07/velocity_verlet.py`, `ch11/gp.py`) are byte-for-byte the same
   algorithms shown in the docs. The tests in `tests/` exercise these
   modules — *not* the docs separately — which means the in-text examples
   are runnable by construction. `tests/conftest.py` even injects
   `code/` onto `sys.path` so the import surface mirrors what a reader
   would type. This is the right architectural choice and is rare.

2. **Symplecticity, shadow-Hamiltonian, and reversibility are derived,
   not asserted.** `docs/ch07-md/01-integration.md` §7.1 contains a full
   Taylor-cancellation derivation of position-Verlet (lines 107–124),
   the explicit shear-decomposition proof of symplecticity (lines
   239–245), a quantitative stability analysis with eigenvalues of the
   Verlet update map on a harmonic oscillator (lines 312–314), and a
   correct factor-of-four diagnostic derived from σ_H ∝ Δt^p (lines
   299–301). The chapter also explicitly equates leapfrog and
   velocity-Verlet trajectories at floating-point precision (lines
   449–461), which is the right level of care. Tests in
   `tests/test_md_integrator.py` then check energy spread, secular
   drift, and oscillation period — exactly the three diagnostics named
   in the text.

3. **The Hohenberg–Kohn proof is sharp.** `docs/ch05-dft/02-hohenberg-kohn.md`
   lines 67–136 gives the three-step proof in the form
   one would actually want to read: the strict-inequality argument
   includes the inline justification that |Ψ⁽¹⁾⟩ ≠ |Ψ⁽²⁾⟩ via the
   one-body difference operator (line 90), and the
   non-degeneracy assumption is foregrounded. The corollaries (lines
   141–145) draw out the Kato cusp / inverse-DFT picture without
   over-claiming. This is graduate-level material done at graduate-level
   rigour.

## Concerns

### 1. CI matrix is too narrow for the declared support window

`environment.yml` line 6 declares `python>=3.11,<3.13`, but
`.github/workflows/test.yml` line 15 sets `python-version: ["3.11"]` —
a single value, not a matrix. There is no `os` axis either; only
`ubuntu-latest`. Yet `tests/` runs fine on macOS Python 3.9.6 in my
session (48 passed, 1 skipped, 2.7 s), so the de-facto support is
broader than what CI verifies. Either tighten `environment.yml` to
`python==3.11` and document that as the supported version, or
broaden the matrix to `["3.11", "3.12"]` × `[ubuntu-latest, macos-latest]`
so the declared support is checked.

The `deploy.yml` workflow runs `mkdocs build --strict --verbose` (line
45), which is good, but it does *not* run pytest as a precondition for
deployment. A red test run can ship a green site. The two workflows
should be chained, or `deploy.yml` should depend on the `test` job.

### 2. `fig_rdf_real.py` uses a single density for solid/liquid/gas

`scripts/figures/fig_rdf_real.py` lines 246–249 sets `rho = 0.030 Å⁻³`
once and reuses the same FCC box for all three temperatures. For
argon, 0.030 Å⁻³ is approximately the liquid density at the triple
point but is *too dense* for the gas phase (real gas-phase Ar at 600 K
and ambient P has rho ≈ 10⁻⁵ Å⁻³) and slightly too low for the FCC
solid at 50 K (rho_solid ≈ 0.026 Å⁻³ if we are honest about lattice
parameter, but compressed to 0.030 it is over-bound). The resulting
g(r) curves are still pedagogically correct in *shape* — sharp peaks
for the cold case, smeared peaks for the warm case, decay to 1 for the
hot case — but the caption ("characteristic of … a gas (T = 600 K)") is
quantitatively misleading: this is a *supercritical fluid* at 600 K,
600× critical density, not a gas. The fix is one of:
 (a) tune ρ to phase-appropriate values for each run,
 (b) re-label the third curve "hot fluid" or "supercritical",
 (c) keep the density fixed and add one sentence to the docstring and
     caption to that effect.

### 3. GP API has a load-bearing transpose heuristic

`code/tier1/ch11/gp.py` lines 53–55 contain
```python
if X.shape[0] != y.shape[0]:
    X = X.T
```
and the same pattern recurs in `predict` (lines 68–70). This silently
transposes the input when the row/column ordering does not match the
target, but it makes the contract ambiguous: a user passing `X` of
shape `(d, n)` with `n != d` gets auto-corrected; a user passing
`X` of shape `(d, d)` does not. The correct posture is to require
`X` shape `(n, d)` (as the docstring already implies in `rbf_kernel`)
and `raise ValueError` if the dimensions disagree. Heuristic transposes
are the source of half of all GP bugs I have seen in undergraduate
projects.

A second, smaller issue: `optimise_hyperparameters` parameterises the
three hyperparameters with `np.exp(params)` (line 87) but seeds the
restart with `rng.normal(0.0, 1.0, size=3)` (line 114). The initial
guess therefore samples sigma_f, ℓ, sigma_n uniformly in [exp(-3),
exp(3)] ≈ [0.05, 20], which is reasonable but undocumented. State this
in the docstring or, better, parameterise the prior explicitly.

### 4. Minor: figure scripts are not pinned and have no batch runner

`scripts/figures/` has 22+ scripts. They mostly compile cleanly (I
spot-checked `fig_rdf_real.py`, `fig_si_bands_real.py`,
`fig_bo_convergence.py`, `fig_ei_acquisition.py`,
`fig_msd_diffusion.py` — all pass `python3 -m py_compile`). All three
scripts using randomness seed an explicit `np.random.default_rng(...)`
(see `fig_bo_convergence.py` line 10, `fig_gaussian_kernel.py` line 10,
`fig_rdf_real.py` line 239). That's already better than most repos.
But there is no `scripts/figures/Makefile`, no `make figures` target,
no requirements lock-file, and no per-figure metadata file that says
"this PNG was produced from this script at this commit". For a v1.0
release that wants to be cited in papers, this is the gap.

### 5. Minor: `kinetic_matrix` in scf_1d uses dense matrices for periodic FD

`code/tier1/ch05/scf_1d.py` `kinetic_matrix` (lines 27–34) constructs
the periodic-BC kinetic matrix densely. For the 1D hydrogen-chain
demo (a few dozen grid points) this is fine. But the text in
`docs/ch04-quantum/03-particle-in-box.md` lines 466–467 already flags
sparse storage as essential for production work. The tier-1 SCF
module would be a better example of the right approach if it
demonstrated `scipy.sparse.diags` + `eigsh` for the same problem;
right now it teaches the wrong habit by example, even though it
explains the right habit in prose.

## Scores (0–10)

| Dimension | Score | Note |
|---|---|---|
| Code runtime correctness | 9 | 48/48 tests pass; sample figure scripts compile; force formulae in `fig_rdf_real.py` correct after manual sign-check. |
| Test coverage and CI quality | 7 | Test suite is small but pointed; CI matrix too narrow; deploy job does not require tests. |
| Derivation rigour | 9 | Verlet, HK I, Schur–GP and SCF Jacobian derivations all check out; one minor inline argument in HK I (line 90) is correct but could be promoted out of a parenthesis. |
| Numerical method appropriateness | 8 | Right choices throughout (FD second-order, velocity-Verlet, RBF GP with Cholesky, sp3d5s* TB for Si). Dense storage in `scf_1d.kinetic_matrix` is the only weak spot. |
| Reproducibility (seeds, versions, citations) | 7 | All three random scripts seed RNG; conftest seeds rng=12345; environment.yml not pinned at minor version; no figure-rebuild target. |
| Convention consistency (units) | 9 | SI in ch04, atomic units announced from ch05 onwards (consistent with docs/ch04 line 469), eV/Å/ps in ch07 MD. The conversion constant `AMU_TO_EV_PS2_PER_A2 = 1.0364e-4` checks out manually. |
| Best-practice adherence (type hints, docstrings, dep mgmt) | 9 | `from __future__ import annotations` throughout, full type hints with `numpy.typing.NDArray`, dataclasses for state containers, NumPy-style docstrings. Above average for a research codebase. |

**Mean score:** 8.3 / 10.

## Suggested fixes for v1.1

1. **Broaden CI matrix** in `.github/workflows/test.yml` to
   `python-version: ["3.11", "3.12"]` and `os: [ubuntu-latest, macos-latest]`.
   Add a `needs: test` clause to the `build` job in `deploy.yml` so the
   site cannot ship when tests are red.

2. **Fix `fig_rdf_real.py` phase densities.** Either parametrise each
   phase with its own (rho, T) pair or relabel the 600 K curve and add
   a docstring/caption sentence acknowledging the fixed-density
   simplification.

3. **Tighten the GP input contract** in `code/tier1/ch11/gp.py`:
   require `X` of shape `(n_samples, n_features)`, raise on mismatch
   instead of auto-transposing, and document the
   hyperparameter-restart prior in `optimise_hyperparameters`.

4. **Add `scripts/figures/Makefile`** with one target per figure (or a
   single `make figures` target that re-runs all 22+ scripts in
   parallel) and a per-figure docstring header recording inputs and
   expected runtime. Also pin `requirements-lite.txt` and
   `environment.yml` to specific minor versions so figure regeneration
   reproduces bit-for-bit on CI.

5. **Replace dense FD kinetic in `code/tier1/ch05/scf_1d.py`** with the
   sparse construction (`scipy.sparse.diags`, `scipy.sparse.linalg.eigsh`).
   This both teaches the right pattern and matches the prose in
   `docs/ch04-quantum/03-particle-in-box.md` §4.3.5.

## Recommendation

**Accept with minor revisions.**

The methodology is solid, derivations are correct, code is genuinely
runnable and tested, and the few methodological flaws I found are all
fixable in a single afternoon of focused work. With the five v1.1
items addressed this would clear a JOSS review and most
Nature-Reviews-Physics computational-methods bars.

Signed,
Reviewer 1 (Methodology)
