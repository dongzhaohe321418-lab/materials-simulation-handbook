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

!!! note "Why this section?"
    Before any equation, we have to agree on what we are optimising and
    over what. The bandit picture is the cleanest formal toy that
    captures the exploration-exploitation trade-off, and every Bayesian
    optimisation algorithm in §11.3 is a continuous generalisation of a
    bandit strategy. Get the bandit intuition right and the rest of the
    chapter is just careful bookkeeping.

!!! note "Key Idea (Box 11.1.A)"
    Every BO algorithm trades two desiderata against each other:
    *exploitation* (pick the candidate the surrogate says is best)
    and *exploration* (pick a candidate whose evaluation most
    reduces the surrogate's uncertainty). Acquisition functions in
    §11.3 are precise ways of scalarising this trade-off, but the
    underlying intuition is fully captured by the multi-armed bandit
    of this section.

## 11.1.1 The bandit framing

A multi-armed bandit is a row of slot machines, each with an unknown
expected payout. At each pull you choose one machine, observe a noisy
sample of its payout, and update your beliefs. After $T$ pulls the
total payout — the sum of what you collected — is what you want to
maximise.

!!! note "Formal definition"
    A $K$-armed bandit is specified by $K$ unknown distributions
    $\mathcal{P}_1, \ldots, \mathcal{P}_K$ over $\mathbb{R}$ with means
    $\mu_1, \ldots, \mu_K$. At each round $t = 1, \ldots, T$, the
    agent selects an arm $A_t \in \{1, \ldots, K\}$ and observes a
    reward $X_t \sim \mathcal{P}_{A_t}$ independently of all past
    rewards. The agent's *cumulative regret* after $T$ rounds is
    $$
    R_T = T \mu^* - \mathbb{E}\!\left[\sum_{t=1}^T X_t\right],
    \quad \mu^* = \max_k \mu_k.
    $$
    A *sublinear* regret algorithm satisfies $R_T = o(T)$, meaning
    the average per-round suboptimality $R_T / T \to 0$ as
    $T \to \infty$. Algorithms that achieve this — UCB1, Thompson
    sampling, $\epsilon_t$-greedy with $\epsilon_t \to 0$ — are
    *consistent*; pure exploitation is not.

!!! note "Derivation: $\epsilon$-greedy regret"
    The $\epsilon_t$-greedy algorithm pulls the empirically best arm
    with probability $1 - \epsilon_t$ and a uniform random arm with
    probability $\epsilon_t$. Auer, Cesa-Bianchi and Fischer (2002)
    show that with $\epsilon_t = \min(1, cK / (d^2 t))$ for suitable
    constants $c$ and $d$, the regret is bounded by
    $R_T \leq \frac{cK}{d^2} \log T + O(K)$, i.e. $R_T = O(\log T)$.
    The intuition is that exploration costs are $\sum_t \epsilon_t / K$
    per arm, which is $O(\log T)$ if $\epsilon_t \propto 1/t$; this
    explores enough to identify the best arm but no more than necessary.

    A simpler, slightly looser argument gives $R_T = O(\sqrt{T \log T})$
    for fixed $\epsilon = T^{-1/3}$, which is the textbook bound. The
    point is: even simple strategies achieve regret growing much more
    slowly than $T$, provided they balance exploration and
    exploitation with care.

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

!!! example "Bandit as materials screening"
    Map the bandit to a discrete materials problem. The arms are
    $K = 30$ candidate alloys; each arm's reward distribution
    $\mathcal{P}_k$ is the (noisy) hardness of alloy $k$ when measured
    in the lab. Each pull = one synthesis-and-measurement cycle, taking
    a week and costing £2 000. The budget is $T = 50$ syntheses. The
    goal is not strictly cumulative reward — we usually care about
    finding the best alloy by the end — but the framework adapts: this
    is the *pure exploration* or *best-arm identification* variant of
    the bandit, where the relevant guarantee is on the *simple regret*
    $\mu^* - \mu_{\hat k_T}$, with $\hat k_T$ the arm recommended after
    the budget is spent. Algorithms like *successive elimination* and
    *UCB-E* are designed exactly for this scenario.

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

!!! note "Concrete numbers for cost-aware acquisition"
    Suppose two candidates have predicted improvements 0.5 and 0.3 in
    arbitrary units, with respective evaluation costs £100 (CGCNN
    prediction, milliseconds of compute) and £1 (CGCNN with
    standardised cost) and £10 000 (one synthesis run). Cost-normalised
    EI gives 0.5/100 = 0.005 and 0.3/1 = 0.3. The cheaper candidate is
    chosen even though its predicted improvement is lower, because the
    bang-for-buck is sixty times better. This is the right way to use
    a fast surrogate alongside expensive oracles: cheap evaluations
    can be wasted; expensive ones must justify themselves.

    For a multi-fidelity GP with two fidelities, the correlation
    coefficient $\rho \in [0, 1]$ between fidelities determines how
    much the cheap oracle informs the expensive one. $\rho = 0.9$ (a
    well-calibrated MLIP) means cheap evaluations almost entirely
    pin down expensive ones; $\rho = 0.3$ (a poorly-calibrated
    surrogate) means the cheap oracle informs but does not replace.
    Measuring $\rho$ empirically — on a small set of paired
    cheap/expensive evaluations — is a prerequisite for sensible
    multi-fidelity BO.

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

## 11.1.5a A worked toy problem: two arms, ten pulls

To make the trade-off arithmetic concrete, consider a two-armed
Bernoulli bandit. Arm A pays out 1 with unknown probability $p_A$ and 0
otherwise; same for arm B with probability $p_B$. The true values are
$p_A = 0.6, p_B = 0.5$, but the agent does not know these.

Suppose after 4 pulls each the agent has empirical means
$\hat p_A = 3/4 = 0.75$ (three successes), $\hat p_B = 1/4 = 0.25$
(one success). The 95% confidence intervals (using Hoeffding's
inequality, $\pm \sqrt{\log(2/0.05) / (2 \cdot 4)} \approx \pm 0.68$)
are $[0.07, 1.0]$ for A and $[-0.43, 0.93]$ for B (clipped to $[0, 1]$).

*Pure exploitation* picks A every remaining round; expected total
reward over the remaining 2 rounds is $2 p_A = 1.2$. But because the
intervals overlap heavily, this is a high-variance commitment.

*UCB1* computes $\hat p_k + \sqrt{2 \log t / n_k}$ for each arm and
picks the larger; this favours A initially but shifts towards B as the
exploration bonus for B (with smaller $n_B$ later) accumulates. The
expected regret over 100 rounds with this strategy is roughly
$8 \log 100 \cdot \Delta / \Delta^2 = 8 \cdot 4.6 \cdot 0.1 / 0.01 \approx 370$
unit pulls in the worst case — a sublinear bound that translates to an
average regret per round of $370/100 = 3.7$ scaled units. In practice
the empirical regret is much smaller, but the upper bound illustrates
the form.

The lesson generalises to materials: identifying the better of two
candidate compositions requires a number of pulls that scales like
$1/\Delta^2$, where $\Delta$ is the difference in their true values.
Twin candidates with similar predicted properties take many
experiments to discriminate; obvious favourites are confirmed quickly.

## 11.1.5b When randomised baselines beat clever schemes

A surprisingly stubborn empirical observation: on many real materials
problems, a *random sampling* baseline performs nearly as well as a
sophisticated BO scheme. Three reasons.

First, *small datasets*. With $\leq 20$ queries and high noise, the
surrogate is so uncertain that its mean predictions are barely better
than random. The BO acquisition function is then optimising over noise.

Second, *flat objective landscapes*. If most candidates are roughly
equivalent and only a few are dramatically better, random sampling
will find one of the few by sheer coincidence with reasonable
probability. BO's advantage shows when the landscape is structured.

Third, *miscalibrated uncertainty*. If $\sigma(\mathbf{x})$ from the
surrogate is systematically wrong — too confident, or too uncertain in
the wrong places — the acquisition function is built on a faulty
foundation. Random sampling does not depend on uncertainty estimates
and is robust to this failure mode.

The corollary: always include a random baseline in any BO experiment.
If your fancy scheme does not beat random, the right response is to
diagnose the surrogate (uncertainty calibration, kernel choice) before
blaming the acquisition function.

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
is high *and* the uncertainty is non-negligible.

!!! tip "Cross-reference"
    The CGCNN of Chapter 10 will appear in §11.4 as a *featuriser*: its
    final-layer embedding becomes the input vector to the GP we are
    about to build. The trade-off discussed here governs how that
    integrated CGCNN-GP pipeline allocates experimental budget.
    Chapter 5's lattice-symmetry analysis enters indirectly through
    constraints on candidate compositions: only certain symmetries are
    physically accessible, and a candidate space restricted by these
    constraints is the right domain for BO.

### Section summary

- The multi-armed bandit formalises exploration-exploitation; the
  cumulative-regret framework (Theorem-style bounds for
  $\epsilon$-greedy, UCB1) tells us pure exploitation has linear
  regret while balanced strategies have sublinear regret.
- Materials versions of the bandit: which composition to synthesise
  next, which candidate to compute DFT on; the cost of "pulls" is
  highly non-uniform.
- Cost-aware acquisition functions divide $\alpha(\mathbf{x})$ by
  evaluation cost; multi-fidelity BO couples cheap surrogates (MLIPs)
  with expensive oracles (DFT, experiment).
- Active learning and Bayesian optimisation share machinery but
  differ in objective: build a better model vs. find a better
  candidate. Sections 11.3 and 11.4
make this concrete; in particular, Section 11.3's worked code plots
exactly this picture for a 1D regression example.

The trade-off is the conceptual backbone. Whenever a step of a BO
algorithm seems counterintuitive — why did it pick the third-ranked
candidate rather than the first? — the answer is almost always that the
algorithm was paying for information rather than for immediate gain.
Whether to applaud or correct that decision depends on whether
information acquisition is part of the campaign's goal.

!!! tip "Forward reference"
    Chapter 12 will reinterpret the same trade-off in the context of
    *foundation models for materials*: a pre-trained model arrives with
    its own prior over function values, and the user's role becomes
    choosing which queries best update that prior in the direction of
    the user's specific goal. The mathematics — acquisition functions
    over a Bayesian surrogate — is unchanged; what changes is the
    quality of the prior.

!!! note "Active learning versus Bayesian optimisation"
    A semantic distinction worth pinning down. *Active learning* and
    *Bayesian optimisation* share the surrogate-and-acquisition
    machinery but differ in their objective.

    *Active learning*: minimise the model's uncertainty over the input
    space. The acquisition is some uncertainty score (variance,
    information gain). The output is a *better model*. Used when the
    goal is to build a property predictor with as few labels as
    possible.

    *Bayesian optimisation*: maximise (or find the maximum of) the
    underlying function. The acquisition trades off mean and
    uncertainty. The output is *the best candidate found*. Used when
    the goal is to identify a high-performing material, not to predict
    every material's property.

    The two are confused often enough to deserve clarification. The
    surrogate construction (§11.2) is identical; the acquisition
    function differs. Most of Chapter 11 is BO; Exercise 11.3 contains
    an active-learning variant.

In the next section we develop the Gaussian process machinery that
makes "uncertainty" a precise quantity. After that, the rest of the
chapter is the systematic study of how to act on it.
