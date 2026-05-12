# 11.4 BO for Materials Discovery

<figure markdown>
![BO vs random search regret](../assets/figures/ch11/fig_bo_convergence.png){ width="650" }
<figcaption>Figure 11.4.1. Simple-regret curves for Bayesian optimisation (EI) versus random search on a synthetic benchmark, averaged over 30 random seeds. BO converges roughly an order of magnitude faster on this problem — the same advantage typically seen on real materials-discovery tasks where each evaluation is expensive.</figcaption>
</figure>

Sections 11.1–11.3 built the Bayesian optimisation machinery in the
abstract. This section connects it to materials. The questions we must
answer to apply BO to a real campaign are practical and unglamorous.
What does a candidate material look like as a numerical vector? What
counts as the oracle, and how expensive is it? How do we wire all of
this together into a campaign that interleaves prediction, querying
and update? Two case studies — perovskite band-gap optimisation and
catalyst screening with BoTorch — provide the structure; the broader
discussion treats featurisation, multi-fidelity coupling, and
autonomous experimentation.

## 11.4.1 Featurising a material for BO

Bayesian optimisation operates on a continuous input space — the GP
kernel needs to compute distances between inputs. A material, however,
is naturally described by a discrete chemical formula and an
arrangement of atoms in a unit cell. Bridging these two views is the
featurisation problem.

Three approaches dominate.

**Compositional descriptors.** Encode each material by a fixed-length
vector summarising its composition. The *Magpie* set of Ward et al.
(2016) is the workhorse: 132 features computed by taking weighted
averages and ranges of elemental properties (atomic number,
electronegativity, group, period, atomic radius, etc.) over the
constituent elements. For a binary oxide A$_x$B$_y$O$_z$, the Magpie
features summarise the elemental properties of A and B weighted by
their stoichiometry. Magpie features are cheap, deterministic, and
respect the natural smoothness of composition space — small changes
in stoichiometry produce small changes in the feature vector.

**Structural descriptors.** When the structure (not just composition)
matters, we need descriptors like SOAP, ACSF, or pre-trained GNN
embeddings (Chapter 9, Chapter 10). These capture atomic
neighbourhoods directly but can be high-dimensional. For BO with a
GP, high-dimensional features are problematic — the curse of
dimensionality bites — so one typically projects to 50–100 dimensions
via PCA or learns a low-dimensional embedding.

**Learned latent features.** Train a property-prediction GNN
(Chapter 10) on a large database, then use the last hidden layer's
activations as a fixed feature for BO. The features are tailored to
the property of interest, dense, and typically 64–256-dimensional. The
cost is training the GNN once upfront; the benefit is featurisation
quality that no hand-crafted descriptor can match. This is the
modern default.

For one-off campaigns with small candidate sets ($\lesssim 10^3$),
Magpie features and an RBF-kernel GP work well. For larger campaigns
or for problems where structure-property relationships dominate, GNN
embeddings are worth the upfront training cost.

## 11.4.2 Case study 1: perovskite band-gap optimisation

The first case study is a stylised but realistic problem. We have a
library of cubic perovskite candidates ABX$_3$ with $A \in
\{$Cs, Rb, K, Na, MA$\}$ (MA = methylammonium), $B \in
\{$Pb, Sn, Ge$\}$, $X \in \{$I, Br, Cl$\}$, giving $5 \times 3 \times 3 = 45$
candidates. Each candidate's "true" band gap is, in our simulation,
a known function returned by a hypothetical DFT calculation. The goal
is to find the candidate with band gap closest to 1.4 eV (the
Shockley–Queisser optimum for solar cells) while running as few DFT
calculations as possible.

We featurise each candidate via Magpie:

```python
from __future__ import annotations

import itertools

import numpy as np
from matminer.featurizers.composition import ElementProperty
from pymatgen.core import Composition


A_set = ["Cs", "Rb", "K", "Na", "C N H6"]  # methylammonium proxy
B_set = ["Pb", "Sn", "Ge"]
X_set = ["I", "Br", "Cl"]

featuriser = ElementProperty.from_preset("magpie")

candidates: list[dict] = []
for A, B, X in itertools.product(A_set, B_set, X_set):
    formula = f"({A})({B})({X})3"
    comp = Composition(formula)
    features = featuriser.featurize(comp)
    candidates.append({"A": A, "B": B, "X": X, "features": np.array(features)})

X_all = np.vstack([c["features"] for c in candidates])  # (45, 132)

# Standardise.
mean = X_all.mean(axis=0)
std = X_all.std(axis=0) + 1e-8
X_all = (X_all - mean) / std
```

For the simulated objective we use a smooth function of three
chemical proxies — average ionic radius, Pauling electronegativity
difference and electron affinity — designed so the minimum-gap-error
candidate lies near MA-Pb-I (the real-world champion lead-halide
perovskite, which has gap $\sim 1.55$ eV).

```python
def simulated_band_gap(features: np.ndarray) -> float:
    # Some smooth function of the standardised features.
    # In reality this would call DFT or an MLIP.
    rng_global = np.random.default_rng(13)
    coefs = rng_global.normal(0, 1, size=features.shape[0])
    return float(0.7 * np.sin(features @ coefs / 8) + 1.5)


def gap_error(features: np.ndarray) -> float:
    """Distance from the target band gap of 1.4 eV."""
    return -abs(simulated_band_gap(features) - 1.4)
```

We negate the gap error so the BO objective is to *maximise*
$-|\mathrm{gap} - 1.4|$, which puts the optimum at gap exactly 1.4 eV.

The BO loop, using our GP from §11.2 and EI from §11.3:

```python
# Initialise with three random candidates.
rng = np.random.default_rng(42)
seen_idx = list(rng.choice(45, size=3, replace=False))
remaining_idx = [i for i in range(45) if i not in seen_idx]

X_train = X_all[seen_idx]
y_train = np.array([gap_error(X_all[i]) for i in seen_idx])

history = list(y_train.copy())

for iteration in range(15):
    gp = GP()
    gp.optimise_hyperparameters(X_train, y_train)
    X_remaining = X_all[remaining_idx]
    ei = expected_improvement(gp, X_remaining, float(np.max(y_train)), xi=0.01)
    next_local_idx = int(np.argmax(ei))
    next_global_idx = remaining_idx[next_local_idx]
    next_y = gap_error(X_all[next_global_idx])

    seen_idx.append(next_global_idx)
    remaining_idx.pop(next_local_idx)
    X_train = np.vstack([X_train, X_all[next_global_idx]])
    y_train = np.append(y_train, next_y)
    history.append(float(np.max(y_train)))

    print(
        f"iter {iteration:2d}: queried "
        f"({candidates[next_global_idx]['A']}, "
        f"{candidates[next_global_idx]['B']}, "
        f"{candidates[next_global_idx]['X']}) "
        f"gap_error = {next_y:+.4f}, "
        f"best so far = {np.max(y_train):+.4f}"
    )
```

A typical run identifies the global optimum within about ten
iterations, versus 23 for random sampling (which on average needs
half the candidates). The BO has roughly halved the number of DFT
calculations required to find the best candidate, at the cost of
some bookkeeping.

The point of this small example is not the absolute numbers — the
simulated objective is a toy — but the workflow. The same loop, with
a real DFT oracle in place of `simulated_band_gap`, drives real
campaigns. The user types the loop once, hits go, and the campaign
runs for as many iterations as the DFT budget permits, returning the
best candidate found at the end.

## 11.4.3 Case study 2: catalyst screening with BoTorch

For larger campaigns, hand-rolling GPs becomes burdensome. The
production tool is *BoTorch*, the Bayesian optimisation library built
on PyTorch and GPyTorch. BoTorch handles GP fitting, acquisition
function optimisation, batch acquisitions and constraints with a
unified API.

The second case study is catalyst screening: a library of 200
candidate bimetallic catalysts, each featurised by 32-dimensional GNN
embeddings, with the objective being CO$_2$ reduction selectivity
(higher is better). DFT-based microkinetic models give the labels.

```python
import torch
from botorch.models import SingleTaskGP
from botorch.fit import fit_gpytorch_mll
from botorch.acquisition import qExpectedImprovement
from botorch.optim import optimize_acqf_discrete
from gpytorch.mlls import ExactMarginalLogLikelihood


# Suppose features and labels are precomputed.
# features: (200, 32) GNN embeddings, normalised to unit variance.
# labels:   (200,)  selectivities.

features = torch.tensor(features_np, dtype=torch.float64)
labels = torch.tensor(labels_np, dtype=torch.float64)

# Pick three initial candidates at random.
torch.manual_seed(0)
init_idx = torch.randperm(200)[:3]
seen_idx = set(init_idx.tolist())

X_train = features[init_idx]
y_train = labels[init_idx].unsqueeze(-1)

for iteration in range(20):
    gp = SingleTaskGP(X_train, y_train)
    mll = ExactMarginalLogLikelihood(gp.likelihood, gp)
    fit_gpytorch_mll(mll)

    # Acquisition function: batch EI with batch size 1.
    ei = qExpectedImprovement(model=gp, best_f=y_train.max())

    # Restrict optimisation to the unseen candidates.
    candidates = torch.stack([
        features[i] for i in range(200) if i not in seen_idx
    ])
    candidate_indices = [i for i in range(200) if i not in seen_idx]

    # Discrete optimisation: evaluate EI at every candidate and pick the best.
    with torch.no_grad():
        ei_values = ei(candidates.unsqueeze(1)).squeeze(-1)
    best_local = int(torch.argmax(ei_values))
    best_global = candidate_indices[best_local]

    seen_idx.add(best_global)
    X_train = torch.cat([X_train, features[best_global].unsqueeze(0)])
    y_train = torch.cat([y_train, labels[best_global].view(1, 1)])

    print(
        f"iter {iteration:2d}: queried candidate {best_global}, "
        f"y = {labels[best_global].item():.4f}, "
        f"best so far = {y_train.max().item():.4f}"
    )
```

A few features of BoTorch worth highlighting:

- `SingleTaskGP` handles GP construction, including reasonable kernel
  defaults (Matérn-5/2), automatic relevance determination (a separate
  length scale per input dimension), and noise estimation.
- `fit_gpytorch_mll` runs L-BFGS to optimise the marginal log
  likelihood.
- `qExpectedImprovement` is the *batch* extension of EI; with batch
  size 1 it reduces to the standard EI we derived in §11.3. It uses
  Monte Carlo sampling to estimate the expectation, which is necessary
  in the batch case where closed-form EI no longer applies.
- `optimize_acqf` or `optimize_acqf_discrete` finds the maximiser of
  the acquisition over either a continuous box or a discrete candidate
  set, respectively.

For continuous search spaces — say, optimising a continuous
composition $x \in [0, 1]$ in (Pb$_{1-x}$Sn$_x$)I$_2$ — replace
`optimize_acqf_discrete` with `optimize_acqf` and supply box bounds.
BoTorch will then run multi-start gradient ascent on the acquisition
surface.

## 11.4.4 Coupling BO with DFT and MLIPs as oracles

In the examples above, the *oracle* — the function we are querying —
was hidden behind a function call. In a real campaign the oracle is
either DFT (slow, accurate, deterministic) or an MLIP from Chapter 9
(faster, approximately accurate, ideally trained on the relevant
chemistry).

Two patterns dominate.

**Single-fidelity DFT.** Each BO query triggers a full DFT calculation
of the candidate's target property. Calculation time is hours; the BO
loop runs for tens of iterations. This is the right pattern when DFT
accuracy is non-negotiable and the candidate count is modest. A
trained MLIP can be used to *initialise* candidate geometries before
DFT relaxation, dramatically reducing the DFT cost per query.

**Multi-fidelity DFT + MLIP.** Each BO query chooses both *what* to
evaluate and *at which fidelity*. The MLIP is the cheap fidelity, DFT
the expensive one. The GP models the correlation between fidelities;
many MLIP evaluations refine the surrogate's view of the candidate
space, with sparse DFT evaluations correcting any MLIP bias. BoTorch's
`MultiFidelityGP` and `qMultiFidelityKnowledgeGradient` make this
workflow accessible.

The conceptual gain is large. Suppose DFT costs 1 hour per query and
MLIP costs 1 second; a 100-query DFT-only campaign takes 100 hours,
but a campaign with 10 DFT queries and 1000 MLIP queries can locate
the optimum at perhaps 10 hours of DFT plus 17 minutes of MLIP — an
order of magnitude saving in DFT cost. The catch is that the MLIP
must be reasonably calibrated; if its predictions are systematically
biased, the multi-fidelity model will inherit the bias.

## 11.4.5 Closed-loop autonomous experimentation: A-Lab and beyond

The single most cited recent demonstration of materials BO is the
A-Lab from Lawrence Berkeley National Laboratory (Szymanski et al.,
*Nature* 2023). A-Lab combines:

- A generative model to propose candidate compositions.
- A trained GNN to predict synthesis feasibility and stability.
- An active-learning loop that selects which candidates to try next.
- A robotic synthesis platform that physically synthesises the
  candidates.
- Automated characterisation (XRD, electron microscopy) to evaluate
  the results.

The system ran continuously for 17 days, attempted synthesis of 58
target compositions and successfully made 41 of them — over 70%
success rate, where the baseline manual rate was estimated at 20–30%.
The BO loop adapted on the fly: failed syntheses updated the
feasibility model and reshaped subsequent proposals.

The lessons for the practitioner — even one not building a robotic
lab — are general.

First, *the loop must be tight*. Long iteration cycles (weeks between
candidate proposal and the result coming back) starve the BO of
information; the algorithm wastes budget on queries that overlap with
each other. A-Lab's tight cycle of synthesis-characterise-update
within hours is what enables rapid adaptation.

Second, *failure is informative*. A failed synthesis is not a
*missing* data point; it is a labelled negative example. The BO
should incorporate failures into the surrogate, ideally via a
classifier that predicts feasibility alongside the regression of the
target.

Third, *uncertainty calibration matters more than mean accuracy*. A
surrogate with biased means but well-calibrated uncertainties allows
the BO to recognise its blind spots and explore them. A surrogate
with low mean error but miscalibrated uncertainties leads to
overconfident queries in regions where the model is silently wrong.

Most academic BO campaigns will not have a robotic lab attached. The
*workflow* lessons translate anyway: write a loop that closes the
data-update cycle as quickly as your infrastructure permits, log
failures as labelled data, and choose acquisition functions whose
behaviour you understand.

## 11.4.6 Practical workflow checklist

To run BO on a new materials problem, the checklist is:

1. **Choose a featurisation.** Magpie for composition-only problems,
   GNN embeddings for structure-dependent problems. Standardise to
   unit variance.
2. **Define the oracle.** DFT, MLIP, experiment, or simulation. Note
   the cost and the noise level.
3. **Initialise.** A small initial sample (5–20 candidates) via
   random or Latin hypercube sampling. The initial sample's coverage
   matters disproportionately for the rest of the campaign.
4. **Choose a surrogate.** Single-task GP for low-dimensional inputs,
   sparse or deep-kernel GP for high-dimensional, multi-fidelity GP if
   you have multiple oracle levels.
5. **Choose an acquisition.** EI is the safe default; UCB if you want
   explicit exploration control; EHVI if multi-objective.
6. **Run the loop.** Fit GP, optimise acquisition, query oracle,
   append, repeat. Log every step.
7. **Monitor.** Track best-found value versus iteration. Plot the
   GP posterior periodically. Inspect the candidates being chosen for
   plausibility — a BO that consistently queries similar candidates
   needs more exploration; one that ignores promising candidates
   needs less.
8. **Stop.** When the best-found value has plateaued for, say, 10
   iterations, or when the GP posterior maximum's uncertainty drops
   below a threshold.

The workflow is mechanical once the components are in place. The hard
part — and the part that varies most by problem — is choosing the
right featurisation and oracle. The BoTorch documentation and the
materials-BO survey of Lookman et al. (2019) are the right next
readings for a real campaign.

## 11.4.7 Where this leaves us

Chapter 11 has built, from first principles, the tools that turn a
fast surrogate (Chapter 10) into a budget-aware discovery loop. The
mathematical content was the Gaussian process and its acquisition
functions; the materials content was featurisation and oracle
selection; the practical content was BoTorch and the autonomous-lab
exemplar.

Chapter 12 takes the next step: *foundation models* for materials —
universally pre-trained networks whose representations make the
featurisation problem nearly trivial and whose fine-tuning makes the
surrogate fitting nearly automatic. In the chapter-12 framing, the BO
loop becomes one ingredient in a much larger machine that ingests any
materials question and returns the next experiment to run. The
conceptual content of Chapter 11 — the exploration-exploitation
trade-off and the acquisition functions that formalise it — persists
unchanged into that larger machine. It is the part of materials ML
that no amount of foundation-model scale displaces.
