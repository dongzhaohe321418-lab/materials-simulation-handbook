# Methods — Bayesian optimisation for catalyst composition

The protocol has six stages:

1. Define the oracle.
2. Validate the oracle on three reference compositions.
3. Featurise the composition.
4. Build the GP surrogate and choose an acquisition.
5. Run BO + random-search baselines.
6. Analyse regret and posterior.

---

## Stage 0 — Setup

```bash
conda create -n bo-catalyst python=3.11
conda activate bo-catalyst
pip install torch gpytorch botorch ase pymatgen matplotlib pandas
# For the oracle:
conda install -c conda-forge qe   # if DFT oracle
pip install mace-torch            # if MLIP oracle
```

Pin BoTorch (e.g., `botorch==0.10.0`); the API moves between
versions.

---

## Stage 1 — Define the oracle

### 1.1 Slab geometry

Pt–Pd–Ag are all fcc. Use a (111) slab of 4 layers with a fixed
$3 \times 3$ surface unit cell, giving 36 atoms total (9 per layer).
The top 2 layers are relaxed; the bottom 2 are held at the bulk
positions of a virtual-crystal-approximation (VCA) lattice
parameter $a_\mathrm{VCA}$ that interpolates linearly between the
three pure metals.

Set the vacuum layer to $\geq 15$ Å. Include the dipole correction
along **z** (`dipfield = .true.`, `tefield = .true.` in QE).

### 1.2 Adsorption-site convention

Place a single H atom at the fcc hollow site above the centre of
the top layer. This is the canonical site for Pt(111) and remains
defensible across mild alloying. Document the convention.

### 1.3 Composition parametrisation

A composition $(x_\mathrm{Pt}, x_\mathrm{Pd}, x_\mathrm{Ag})$ with
$x_\mathrm{Pt} + x_\mathrm{Pd} + x_\mathrm{Ag} = 1$ defines the
expected number of each atom on a 9-atom surface layer:

$$
n_\mathrm{Pt} = \mathrm{round}(9 \cdot x_\mathrm{Pt}),\quad n_\mathrm{Pd} = \mathrm{round}(9 \cdot x_\mathrm{Pd}),\quad n_\mathrm{Ag} = 9 - n_\mathrm{Pt} - n_\mathrm{Pd}.
$$

For the *surface* layer, sample 3 random arrangements with the
given counts; for the subsurface layers, use the same counts but
fix them as a single arrangement to keep the oracle deterministic.
Average the oracle over the 3 surface arrangements.

(This is the "SQS-light" approach. A proper SQS construction is
more rigorous; for an undergraduate project the 3-arrangement
average is acceptable and is faster.)

### 1.4 The oracle function

```python
def oracle(x: tuple[float, float, float]) -> float:
    """Return Delta G_H in eV for composition (x_Pt, x_Pd, x_Ag)."""
```

The function builds the slab, runs DFT (or MLIP) for `slab` and
`slab+H`, looks up the (precomputed) $E_\mathrm{H_2}$ in vacuum,
and combines.

$$
\Delta G_\mathrm{H} = E_\mathrm{slab+H} - E_\mathrm{slab} - \tfrac{1}{2} E_\mathrm{H_2} + 0.24\,\mathrm{eV}.
$$

Average over the 3 surface arrangements.

### 1.5 Choice of oracle

For most undergraduate projects, oracle B (MLIP) is recommended:

- Use MACE-MP-0 (medium) as the surrogate.
- Relax the slab and the slab+H configuration.
- Compute $\Delta G_\mathrm{H}$.

Validate the MLIP oracle against the DFT oracle on 3 reference
compositions (pure Pt, pure Pd, pure Ag). If the MLIP $\Delta G_\mathrm{H}$
agrees with DFT to within 0.1 eV, the MLIP is fit for purpose at
this descriptor level; otherwise fall back to DFT oracle for the
production run.

---

## Stage 2 — Reference validation

Compute $\Delta G_\mathrm{H}$ for:

- Pt(111): expected $\approx -0.07$ eV (literature).
- Pd(111): expected $\approx -0.20$ eV (literature; H likes Pd).
- Ag(111): expected $\approx +0.40$ eV (Ag is a poor HER catalyst).

If your oracle gives values within ≈ 0.1 eV of these, proceed. If
not, the oracle is broken — likely an incorrect H-reference energy,
a wrong site choice, or an unconverged slab.

Document the validation table in your report.

---

## Stage 3 — Featurisation

Two options:

**Option 1 (recommended): raw composition fractions.** Two-dimensional
feature space $(x_\mathrm{Pt}, x_\mathrm{Pd})$, with $x_\mathrm{Ag} = 1 - x_\mathrm{Pt} - x_\mathrm{Pd}$.
Domain: the unit simplex $\{x_\mathrm{Pt}, x_\mathrm{Pd} \geq 0,\,
x_\mathrm{Pt} + x_\mathrm{Pd} \leq 1\}$.

**Option 2: Magpie compositional features.** A higher-dimensional
feature vector that may help generalisation but adds modelling
complexity. Skip unless you have time.

Use option 1 by default.

---

## Stage 4 — GP surrogate and acquisition

### 4.1 GP setup

In BoTorch:

```python
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from gpytorch.mlls import ExactMarginalLogLikelihood
from gpytorch.kernels import MaternKernel, ScaleKernel

train_x = ...  # shape (n, 2)
train_y = ...  # shape (n, 1), the *negative* |Delta G_H| (we maximise)

model = SingleTaskGP(train_x, train_y,
                     covar_module=ScaleKernel(MaternKernel(nu=2.5, ard_num_dims=2)))
mll = ExactMarginalLogLikelihood(model.likelihood, model)
fit_gpytorch_mll(mll)
```

Normalise inputs to $[0, 1]^2$ and outcomes to zero mean/unit
variance for kernel-hyperparameter learning.

Note the objective transformation: BO maximises by convention; we
want $\Delta G_\mathrm{H}$ close to zero, so we maximise the
*negative absolute value*: $y = -|\Delta G_\mathrm{H}|$. Document
this.

### 4.2 Acquisition

Use Expected Improvement (qEI for batch consistency):

```python
from botorch.acquisition import qExpectedImprovement
from botorch.optim import optimize_acqf

acq = qExpectedImprovement(model, best_f=train_y.max())
candidate, _ = optimize_acqf(
    acq, bounds=torch.tensor([[0.0, 0.0], [1.0, 1.0]]),
    q=1, num_restarts=20, raw_samples=512,
)
```

Reject candidates that violate the simplex constraint $x_\mathrm{Pt}
+ x_\mathrm{Pd} > 1$; re-sample if so.

For comparison, run an additional BO with UCB acquisition:

```python
from botorch.acquisition import UpperConfidenceBound
acq_ucb = UpperConfidenceBound(model, beta=2.0)
```

Use UCB at a beta value of 2.0 for a moderate-exploration setting.

### 4.3 Discretisation

Round each chosen continuous composition to the nearest discrete
one (achievable on the 9-atom surface). Reject duplicates of
already-queried points.

### 4.4 Initial design

Use Latin Hypercube Sampling for 5 initial points inside the simplex.
These bootstrap the GP before the acquisition function takes over.

---

## Stage 5 — Run the BO loop

```python
def bo_loop(n_init: int = 5, n_iter: int = 45,
            acquisition: str = "qEI", seed: int = 0) -> list[dict]:
    history: list[dict] = []
    # 1. Initial design (LHS in simplex)
    init = lhs_simplex(n_init, seed=seed)
    for x in init:
        y = oracle(tuple(x))
        history.append({"x": x.tolist(), "y": y, "stage": "init"})
    # 2. BO iterations
    for it in range(n_iter):
        model = train_gp(history)
        cand = pick_next(model, acquisition)
        cand_d = discretise(cand)
        if any(np.allclose(cand_d, h["x"]) for h in history):
            cand_d = perturb(cand_d)   # avoid exact duplicates
        y = oracle(tuple(cand_d))
        history.append({"x": cand_d.tolist(), "y": y, "stage": f"bo_{it}"})
    return history
```

Run with 3 different `seed` values for statistical robustness.

### 5.2 Random-search baseline

Same total budget (50 oracle calls), but each composition drawn
uniformly from the simplex (Dirichlet(1, 1, 1) is equivalent). Three
seeds.

---

## Stage 6 — Analyse

### 6.1 Regret curve

For each method (BO-qEI, BO-UCB, random) and each seed, compute the
*running best* $|\Delta G_\mathrm{H}|$ after each oracle call:

$$
r(t) = \min_{i \le t} |\Delta G_\mathrm{H}^{(i)}|.
$$

Plot $r(t)$ vs $t$ as a mean ± standard-error band over the 3 seeds.
BO should beat random search; the gap is your headline result.

### 6.2 Posterior heatmap

At the end of the BO loop (with the most data), fit one final GP
and plot the posterior mean over a fine grid of the simplex.
Overlay the queried points. This is your "discovery map".

### 6.3 Top 5 table

Take the 5 lowest-$|\Delta G_\mathrm{H}|$ compositions queried.
Report their actual $\Delta G_\mathrm{H}$, their GP-posterior
uncertainty, and (if MLIP oracle) verify the top 2 with DFT.

---

## Pitfalls

1. **Forgetting to negate.** BO maximises; $|\Delta G_\mathrm{H}|$
   should be minimised. If you forget to feed $-|\Delta G_\mathrm{H}|$
   into the GP, your acquisition will steer toward the worst
   catalysts.
2. **Single seed.** A single BO run can luck out (or fail) on a
   particular initial design. Three seeds is the minimum for a
   meaningful comparison.
3. **Wrong domain.** The simplex has area 1/2 in the
   $(x_\mathrm{Pt}, x_\mathrm{Pd})$ plane, not 1. If you draw random
   compositions uniformly from $[0, 1]^2$ and reject those with sum >
   1, you have effectively halved your random-search density. Use
   `np.random.dirichlet` for uniform sampling on the simplex.
4. **Site fixed to fcc hollow but compositions change.** On a pure-Ag
   slab the H atom may prefer a top site; on a pure-Pd slab it may
   prefer a bridge. To avoid the messiness, *constrain* the site
   choice to fcc hollow throughout. Document this and acknowledge
   that the descriptor you optimise is "fcc-hollow $\Delta G_\mathrm{H}$",
   which is a reasonable proxy but not the true optimum site.
5. **Pre-relaxation differences.** A slab with H at the hollow may
   relax to a slightly different geometry than a slab without H; the
   *energy difference* is what enters $\Delta G_\mathrm{H}$, so the
   absolute geometries do not need to be identical, only consistent.
6. **MLIP-oracle bias.** MACE-MP-0 was not trained specifically on
   Pt–Pd–Ag surfaces with adsorbed H, and its $\Delta G_\mathrm{H}$
   predictions can be biased by ≈ 0.1–0.2 eV. Always validate at
   pure-element endpoints.
7. **Forgetting the discrete projection.** A BO candidate of
   $(0.347, 0.512, 0.141)$ maps to discrete $(3, 5, 1)$ atoms on a
   9-atom surface, equivalent to $(0.333, 0.556, 0.111)$. Always
   record the *discrete* composition that was actually queried, not
   the continuous candidate.

---

## What "done" looks like

You have:

- A working, validated oracle.
- 3 BO runs and 3 random-search runs.
- A regret-curve plot showing BO ≤ random for most of the iteration
  range (i.e., BO actually wins).
- A posterior heatmap of the surrogate.
- A top-5 table with verification.
- A discussion of why BO won, by how much, and where it failed.

If BO did *not* beat random search, that is also a valid result —
investigate why (poor kernel choice, too few initial points,
constrained simplex). Report it honestly.
