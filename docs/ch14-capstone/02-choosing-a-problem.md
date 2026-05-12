# Choosing a Problem

The single most consequential decision in a thesis is the first one: what
problem you decide to work on. Get this right and the rest of the project
will be hard but tractable. Get it wrong and no amount of cleverness in
the methods section will save you.

Most theses succeed or fail based on this initial scoping, and most
students under-invest in it. They are anxious to start *doing* something
— writing code, running calculations — and they accept the first problem
their supervisor proposes without examining whether it is the right
problem for them, at this stage of their training, with the resources
available.

This section is a framework for thinking about problem choice. The
framework has three parts (the *three lenses*), each of which must hold
for the project to be worth attempting.

## The three lenses

A good thesis problem must pass all three of:

**Lens 1 — Interesting physics or chemistry.** The problem must engage
you scientifically. You will spend six months on it; if the question is
boring on day 1, it will be unbearable by month 4.

**Lens 2 — Computationally feasible.** The problem must be solvable with
the resources you actually have. Not the resources you wish you had. Not
the resources your supervisor had ten years ago. The resources today, in
your current group, with your current account on the cluster.

**Lens 3 — Writable.** When the work is done, you must be able to write
it up as a coherent story. The result, whether positive or negative,
should be a defensible claim that fits in a thesis chapter.

A problem that passes only one or two lenses is a trap. Let us look at
each in detail.

## Lens 1: Interesting physics or chemistry

By "interesting" we do not mean "famous" or "trendy". We mean: there is a
specific question whose answer would teach you something you did not
already know. The question can be very modest in scope — predicting one
property of one material — and still be interesting in this sense.

What is *not* interesting:

- "Reproduce the literature value of the bulk modulus of copper." If the
  literature value is well established, your reproduction is not a
  result; it is a sanity check. Sanity checks belong *inside* a thesis,
  not as the thesis itself.

- "Compute the band gap of silicon with DFT." This has been done
  approximately $10^4$ times. Adding one more calculation contributes
  nothing.

- "Apply the latest deep-learning method to a problem in materials
  science." Trendy and vague. Without a specific question, the project
  has no shape.

What *is* interesting:

- "Quantify how the vacancy formation energy in iron changes under
  hydrostatic compression at the conditions of the Earth's outer core."
  Specific question, specific system, useful answer.

- "Train an MLIP for the Al-Cu-Mg ternary system and use it to predict
  the metastable phase boundary at 600 K." Specific question, but
  ambitious on the methodology side.

- "Reproduce the published prediction of a new metastable polymorph of
  TiO$_2$ using a different exchange-correlation functional, and assess
  whether the prediction is functional-robust." A clean test of a
  literature claim, well-bounded.

The pattern in the "interesting" examples: each has a *specific
question*, a *specific system*, and an answer that would change someone's
view (yours, your supervisor's, the field's) by some specific amount.

### How to find a question

If your supervisor has handed you one, your job is to understand it
deeply: take the supervisor's question, restate it in your own words,
identify what is implicit, and pose your version back to them. If their
question was "study the diffusion of Li in spinel-structured cathodes",
your version might be: "compute the migration barrier of Li between
tetrahedral and octahedral sites in LiMn$_2$O$_4$, at three Li
concentrations, using DFT-NEB, and validate against published
electrochemical impedance data".

If you have a list of possible projects (as in this handbook), the
question is which to pick. Use the criteria below.

If you have no question at all and your supervisor is letting you find
your own, the literature is your friend. Read recent reviews in a
sub-field that interests you. The questions a review identifies as *open*
— "the role of X is still poorly understood", "more work is needed on
Y" — are candidate thesis questions. Pick one and run it past your
supervisor.

## Lens 2: Computationally feasible

This is where students most often deceive themselves. The instinct is to
choose a problem that *sounds* important and to assume the calculations
will somehow fit in the available time.

Here are concrete rules of thumb. These are not laws; they vary by
group, by hardware, by software, by your own efficiency. But they
represent realistic ranges for an undergraduate or master's student
working with shared facilities.

### DFT-driven projects

For a DFT-driven thesis on a typical academic cluster (≤ 64 cores
allocated to you on average):

- **Maximum useful supercell size**: 100-200 atoms for routine
  calculations; up to ~500 for a small number of headline calculations.
- **Total number of converged calculations**: 100-500 over the project.
- **Heaviest single calculation**: hours to a day or two.
- **Hybrid functionals (HSE, etc.)**: only on small systems (~30
  atoms) and only when essential; a hybrid-DFT-only thesis is hard at
  this scale.
- **GW or BSE**: not feasible without dedicated computing time.

A laptop-only DFT thesis is possible if you stick to ≤ 50 atom cells, use
a fast plane-wave code (Quantum ESPRESSO) or, more realistically, a
linear-scaling DFT code (ONETEP, BigDFT) or all-electron localised-basis
code (FHI-aims).

### MD-driven projects

For an MD thesis:

- **If a well-validated classical force field exists for your system**:
  great, run MD on $10^3$-$10^5$ atom systems for $10^7$-$10^9$
  timesteps.
- **If you need a force field that does not exist**: training it
  becomes the thesis. Allocate at least three months of the timeline to
  training and validation, with the actual application as the last
  third.
- **If you plan to use an MLIP**: training one from scratch is a thesis
  in itself ([Project 4](07-the-five-projects.md)). Using an existing
  pre-trained foundation model (MACE-MP-0, CHGNet, etc.) is much faster
  but requires careful validation on your system.

### High-throughput projects

For high-throughput or database-driven projects:

- **Query an existing database (Materials Project, OQMD, AFLOW)**: easy,
  scriptable in days. The science is in the analysis.
- **Generate new high-throughput DFT data**: hard to do well in six
  months; expect to spend more time on workflow automation than on
  science.
- **Train an ML model on database data**: feasible, especially with
  established models like MEGNet, ALIGNN, M3GNet.

### Active learning / Bayesian optimisation

For an AL project:

- **Use an existing simulator as the oracle**: feasible. The AL loop
  itself is light; the simulator dominates compute.
- **Use real experiments as the oracle**: outside the scope of a
  computational thesis (unless you are co-supervised by an experimental
  group with this set up already).

### A sanity check before you commit

The single most useful check: estimate the *total CPU-hours* the project
will require and compare to what you actually have.

If your project requires 500 DFT calculations and each takes 4 CPU-hours
on average, that is 2000 CPU-hours of production. Add 50% for convergence
studies, 30% for failed runs, and you are at 3500 CPU-hours. Divide by
the cores you typically get and convert to wall-clock time. If the
answer is "five months on the cluster I share with twelve other people",
your project is too big.

This estimate takes thirty minutes to make and saves months of regret.

!!! warning "The cluster availability tax"
    Shared cluster time is rarely all yours. Queue waits, maintenance
    windows, holidays, and other users' jobs all eat into your real
    wall-clock budget. As a rule of thumb, assume you will get 50-70% of
    the cluster time on paper. Plan accordingly.

## Lens 3: Writable

The third lens is the one students think about least. A "writable"
project is one that, when done, has a clear story.

A story has:

- **A setup**: what was the question, why does it matter, what was
  unknown?
- **A method**: how did you attempt to answer the question?
- **A result**: what did you find?
- **A discussion**: what does the result tell us, and what does it not
  tell us?

If your project, *no matter how it turns out*, does not allow you to
write a story with these four parts, it is not writable.

The most common way for a project to be unwritable: the outcome is
ambiguous. You computed a number, but you have no benchmark to compare
it to; the result is just "the number is what it is". Or: you trained a
model, but the MAE is "okay" without context.

The fix: ensure that at planning time, you can articulate what the
*possible outcomes* are, and that each is writable.

For example, suppose your project is "compute the migration barrier of Li
in a candidate solid electrolyte". Possible outcomes:

- Barrier is *low* (< 0.3 eV): the material is a candidate fast-ion
  conductor; we can compare to published values for known electrolytes
  and discuss in that context. **Writable.**
- Barrier is *high* (> 0.6 eV): the material is not a good electrolyte,
  but we can discuss why (geometry, coordination, charge state) and
  compare to predictions from descriptor models. **Writable.**
- Barrier is *intermediate*: we can frame the discussion around the
  trade-offs that intermediate barriers imply. **Writable.**
- We cannot converge the calculation: this needs to be flagged at the
  *convergence study* stage, not after months of production. **Unwritable
  unless caught early.**

The key habit: think about the discussion section before you start.

## The minimum publishable unit

In the research world, there is an old idea called the "minimum
publishable unit" — the smallest piece of work that can be published as a
standalone paper. It is partly a joke about salami-slicing, but it has a
serious core: there is a minimum unit of *complete result* below which a
piece of work is not a result at all, just a fragment.

For a thesis, the minimum unit is slightly different. A thesis chapter
does not have to be a complete paper, but it does have to *tell a
complete story*. That is, in addition to the work, it has to have:

- An explicit question,
- A defended methodology,
- Results with uncertainty quantified,
- Comparison to prior work or experiment where possible,
- Honest discussion of limitations.

A thesis with three such chapters is much stronger than one with seven
half-finished ones. **Build for completeness, not coverage.**

## Risk register

Once you have a candidate problem, draft a *risk register*: a list of the
things that could go wrong, what would happen if they did, and what your
mitigation is. Examples:

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| K-point convergence not achievable in available wall-time | Low | High | Test on day 1; switch to smaller cell if needed |
| MLIP fails to extrapolate to high pressure | Medium | High | Active-learning loop with held-out high-pressure validation |
| Cluster downtime > 1 month | Low | Medium | Local laptop fallback for small-cell runs |
| The phenomenon I want to study is not present in my simulation cell | Medium | Catastrophic | Use literature precedent to set cell size; test for system-size effects |
| My supervisor goes on sabbatical for two months | High | Medium | Identify a second contact in the group, schedule weekly meetings in advance |

Risks rated *high* probability and *high* impact must have a serious
mitigation strategy or you should reconsider the project. Risks rated
*low/low* are fine to note and move on. The point of the exercise is to
make the failure modes explicit rather than implicit.

## Two more rules of thumb

**Pick a problem your supervisor cares about.** A thesis project where
the supervisor has ongoing engagement gets weekly attention. A thesis
project the supervisor is mildly indifferent to gets attention twice a
semester. The difference, in finished thesis quality, is enormous. If
your supervisor has a preferred direction, *strongly* consider taking
it.

**Pick a problem where you can build on the group's existing
infrastructure.** Every group has scripts, conventions, calculator
setups, established workflows. If your project sits within these, you
get to inherit a year of accumulated work. If your project requires the
group to learn a new piece of software (a different DFT code, a
different MD package), you are simultaneously doing your thesis *and*
porting infrastructure. This is usually a mistake at the undergraduate
or master's level.

## What good scoping looks like, in one paragraph

A well-scoped thesis project, in a single paragraph, might read:

> *I will compute the vacancy formation energy in bcc iron as a function
> of hydrostatic pressure from 0 to 360 GPa, using PBE-DFT in 128-atom
> supercells, with full geometric relaxation. I will compare to two
> published DFT calculations at 0 GPa and to one shock-compression
> experimental dataset at higher pressures. The expected outcome is a
> formation energy curve that can be fit to a simple Murnaghan-style
> form, allowing extrapolation to Earth-core pressures. The main risk is
> that magnetism, which I will not treat explicitly, becomes important
> below the alpha→epsilon transition; I will treat results in that
> region with caution and flag this in the discussion. The total cost is
> approximately 60 production calculations at 6 hours each on 32 cores,
> plus a 2-week convergence study, comfortably fitting in 4 months of
> shared cluster access.*

Note what is in this paragraph: specific question, specific method,
specific comparison points, explicit risk, explicit cost estimate. Note
what is *not*: vague aspirations, "novel applications of", "we will
explore".

The next time you are about to commit to a project, write one of these
paragraphs first.

[Section 3](03-literature-search.md) covers how to fill in the
literature side of this paragraph.
