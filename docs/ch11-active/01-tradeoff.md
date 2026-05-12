# 11.1 The Exploration–Exploitation Trade-off

A scientist has a budget — of synthesis runs, of DFT-core-hours, of
graduate-student months — and a list of candidates. Each candidate, if
investigated, returns some scalar measure of how good it is: a band
gap, a formation energy, a catalytic turnover frequency. The goal is to
identify the best candidate, or the best handful, while spending no
more of the budget than necessary.

This sounds straightforward; it is not. The candidates the scientist
ends up evaluating depend on choices made along the way, and those
choices interact with the budget in a way that has a precise
mathematical structure. The structure is the *exploration–exploitation
trade-off*, and the language for thinking about it crystallised in the
study of a deceptively simple problem: the multi-armed bandit.

## 11.1.1 The bandit framing

A multi-armed bandit is a row of slot machines, each with an unknown
expected payout. At each pull you choose one machine, observe a noisy
sample of its payout, and update your beliefs. After $T$ pulls the
total payout — the sum of what you collected — is what you want to
maximise.

The trade-off is starkest with two machines. Machine A has paid out an
average of 1.0 over your first ten pulls. Machine B has paid out an
average of 0.8 over your first ten pulls. With ninety pulls remaining,
should you pull A every time (exploit your best estimate) or
occasionally pull B (explore in case your estimate of B was unlucky)?

A purely exploitative agent always pulls A, locks in the suboptimal
expected payout if B is in fact better, and never finds out. A purely
explorative agent splits time evenly between the two and converges on
the correct ranking but spends half its budget on the worse machine.
Neither extreme is right; the optimum is some adaptive mixture.

The Gittins index theorem (1979) showed that for the discounted-reward
bandit, there is an exact optimal policy: compute an index for each
machine that depends on its posterior payoff distribution and the
discount factor, and pull whichever machine has the highest index. The
index inflates by an amount that grows with the uncertainty of the
machine's estimated payoff — the explore bonus. Bayesian optimisation
acquisition functions are continuous-space generalisations of this
idea.

For our purposes, the bandit gives the right *language*: every
candidate has an estimated value and an uncertainty about that
estimate; rational decisions weigh both.

## 11.1.2 The materials version

Two recurring materials examples make the abstract trade-off concrete.

**Which alloy to synthesise next?** A materials scientist screening
high-entropy alloys for hardness has thirty candidates ranked by a
neural-network predictor. The top-ranked candidate has predicted
hardness 12 GPa with model uncertainty $\pm 1$ GPa. A candidate ranked
fifteenth has predicted hardness 10 GPa with uncertainty $\pm 4$ GPa.
The first is the safer bet; the second has a meaningful probability of
exceeding the first if the model is wrong. With one synthesis run
available, pick the first. With ten runs available, allocating some to
the high-uncertainty candidates is the better strategy — they could be
better, and learning that they are or are not will inform the
remaining campaign.

**Which candidate to compute DFT on?** A graduate student has trained a
CGCNN (Chapter 10) on 5000 oxides and now has predictions for 100 000
candidates from a generative model. DFT can verify perhaps a hundred
before the conference deadline. The naïve strategy is to compute the
top hundred by predicted formation energy. The active-learning strategy
is to compute the hundred candidates whose verification will most
improve future CGCNN training — typically a mix of high-confidence top
candidates (to confirm them) and high-uncertainty edge cases (to teach
the model). The two strategies frequently produce overlapping but not
identical selections; the active-learning strategy generally yields
both better final predictions and a better-calibrated CGCNN.

**A more domain-specific framing.** When the objective is to *discover*
materials in a region of composition space that has not been studied —
say, exploring the quaternary Mg–Al–O–N system — exploration dominates.
When the objective is to *refine* an already-promising candidate via
small composition variations — say, tuning the Ti/Zr ratio in BaTiO$_3$
to maximise piezoelectric response — exploitation dominates. The
trade-off is not abstract; it is set by the question being asked.

## 11.1.3 Failure modes of pure strategies

It is worth spelling out, in pictures, what goes wrong at each
extreme.

**Pure exploitation.** Every iteration, pick the candidate with the
highest predicted value. The model has no incentive to update. If the
initial training data contained, by chance, no examples in a region
that *actually* contains the optimum, the model will never explore
there. Practically, this manifests as a campaign converging quickly to
a local optimum and never finding the global one. The worse the
initial model, the more egregious the trap.

A real example. A 2018 search for organic photovoltaic candidates used
a random-forest surrogate to suggest the next synthesis from a 30 000-
candidate library. The strategy was pure exploitation. After fifty
synthesis runs, the best efficiency had plateaued. A retrospective
analysis showed that the true optimum lay in a region of chemical
space that contained no initial training examples and which the random
forest extrapolated wildly into negative numbers; the candidates from
that region therefore looked unpromising and were never selected.

**Pure exploration.** Every iteration, pick the candidate the model is
most uncertain about. The model rapidly improves its uncertainty
calibration everywhere — but the campaign's actual best-found value
improves slowly. The model knows a lot about the function; it has not
found the optimum. This is the right strategy if your objective is to
*build a model* (active learning), not to *find an optimum* (BO). The
two objectives are different, and the right strategy depends on which
you have.

**Random sampling.** Pick the next candidate at random. Surprisingly,
this often works better than pure exploitation, because at least it
explores. It is the conservative baseline against which any
acquisition strategy should be measured: if your fancy BO scheme does
not beat random sampling on your problem, you have a calibration bug
or a wrong-prior bug somewhere.

## 11.1.4 Cost-aware acquisition

Real experiments cost money, and they cost different amounts. Synthesising
a sample at the Diamond Light Source costs thousands of pounds in beam
time; running a CGCNN prediction costs milliseconds of GPU time. A
sensible acquisition function should weight each candidate by what it
costs to evaluate. Two natural ways to do this.

**Cost-normalised acquisition.** Replace the acquisition $\alpha(\mathbf{x})$
with $\alpha(\mathbf{x}) / c(\mathbf{x})$, where $c(\mathbf{x})$ is the
cost of evaluating candidate $\mathbf{x}$. The optimum then becomes the
candidate with the best expected improvement *per pound spent*, not
the best absolute improvement. This is the right form when cost is
roughly proportional to time and the budget is in time units.

**Multi-fidelity BO.** Distinguish a *cheap* oracle (a CGCNN
prediction, $c \sim 1$ ms) from an *expensive* oracle (a DFT
calculation, $c \sim$ hours) from an *extremely expensive* oracle
(synthesis and measurement, $c \sim$ days). Build a probabilistic
surrogate over all fidelity levels. At each iteration the algorithm
chooses both which input $\mathbf{x}$ to query *and* at which fidelity
level. Querying at the cheap fidelity many times to localise the
optimum, then verifying at expensive fidelity only at the most
promising point, is often the right strategy.

Multi-fidelity is mathematically cleaner than it sounds: model the
function as a co-kriging GP with a correlation parameter between
fidelities. BoTorch supports this directly via `MultiFidelityGP` and
the multi-fidelity knowledge-gradient acquisition. The case study in
§11.4 will use a two-fidelity GP–CGCNN combination on a perovskite
band-gap problem.

A practical observation about cost. The cost of a single DFT run is
predictable; the cost of a single synthesis run is not — many
syntheses fail outright. *Failure-aware* acquisition functions add a
classifier on top of the regressor that predicts the probability of
successful evaluation, and weight the acquisition by that probability.
This is an active research area; for our purposes we treat costs as
known and deterministic.

## 11.1.5 How the trade-off enters Bayesian optimisation

The remainder of the chapter is mathematics designed to formalise the
trade-off. Section 11.2 builds a probabilistic surrogate — a Gaussian
process — that gives every candidate both a predicted value (the
posterior mean) and an uncertainty (the posterior standard deviation).
Section 11.3 introduces acquisition functions that combine these two
numbers into a single scalar to optimise. The specific form of each
acquisition function corresponds to a specific stance on the trade-off:

- The *upper confidence bound* $\mathrm{UCB}(\mathbf{x}) = \mu(\mathbf{x}) + \kappa \sigma(\mathbf{x})$
  with large $\kappa$ favours exploration; with small $\kappa$ it favours
  exploitation. The parameter $\kappa$ is the trade-off dial, set by the
  user or by a decreasing schedule that explores early and exploits late.
- *Expected improvement* contains its own implicit trade-off — uncertain
  candidates with predicted values close to the current best score
  highly because the chance of a large improvement is nonzero even
  if their mean prediction is modest.
- *Thompson sampling* trades off by *sampling*: it draws a random
  function from the posterior and picks its maximum, which produces a
  natural mixture of exploration and exploitation by virtue of the
  sampled function's randomness.

Different acquisition functions parametrise the trade-off differently;
they all instantiate it.

## 11.1.6 A frame for the rest of the chapter

We will return repeatedly to a simple visual mental model. Plot the
unknown function as a smooth curve over a one-dimensional input. Plot
the GP posterior as a mean line with a shaded uncertainty band. Mark
the observed data points. The acquisition function is a second curve
below; its maximum is the next query.

In that picture, exploitation means the acquisition peaks near the
maximum of the posterior mean; exploration means it peaks where the
uncertainty band is widest. A well-designed acquisition smoothly blends
the two: it peaks slightly above the posterior maximum, where the mean
is high *and* the uncertainty is non-negligible. Sections 11.3 and 11.4
make this concrete; in particular, Section 11.3's worked code plots
exactly this picture for a 1D regression example.

The trade-off is the conceptual backbone. Whenever a step of a BO
algorithm seems counterintuitive — why did it pick the third-ranked
candidate rather than the first? — the answer is almost always that the
algorithm was paying for information rather than for immediate gain.
Whether to applaud or correct that decision depends on whether
information acquisition is part of the campaign's goal.

In the next section we develop the Gaussian process machinery that
makes "uncertainty" a precise quantity. After that, the rest of the
chapter is the systematic study of how to act on it.
