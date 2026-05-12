# Background reading — Bayesian optimisation for catalyst composition

Read in order. Annotations explain how each paper feeds your project.

---

## 1. Nørskov, Bligaard, Logadottir, Kitchin, Chen, Pandelov, and Stimming (2005) — *Trends in the exchange current for hydrogen evolution*

> Nørskov, J. K. et al. *J. Electrochem. Soc.* **152**, J23 (2005).

The foundational HER descriptor paper. Establishes that the
exchange-current density correlates with the Gibbs free energy of
H adsorption, with the optimum at $\Delta G_\mathrm{H} \approx 0$.
Introduces the "computational hydrogen electrode" reference.

**Why read this.** This is your objective function. You must
understand what $\Delta G_\mathrm{H}$ means physically, why it
predicts catalytic activity, and what the volcano plot looks like.

**Extract.** The definition of $\Delta G_\mathrm{H}$ (including the
ZPE + entropy correction of +0.24 eV); the position of Pt, Pd, Ni
on the volcano; the role of the computational hydrogen electrode.

---

## 2. Greeley, Jaramillo, Bonde, Chorkendorff, and Nørskov (2006) — *Computational high-throughput screening of electrocatalytic materials for hydrogen evolution*

> Greeley, J. et al. *Nat. Mater.* **5**, 909 (2006).

The original HT screening of HER catalysts using $\Delta G_\mathrm{H}$.
Surveys many bimetallic surfaces and identifies BiPt and other
candidates.

**Why read this.** Gives you a concrete benchmark — what counts as
"interesting" $\Delta G_\mathrm{H}$ values, what the dynamic range
of the descriptor across alloys is.

**Extract.** The empirical range of $\Delta G_\mathrm{H}$ for
bimetallic surfaces (roughly $-0.5$ to $+0.5$ eV); the typical
slab thickness (3–5 layers); the standard site-coverage convention
(1/4 monolayer is the usual reference).

---

## 3. Lookman, Balachandran, Xue, and Yuan (2019) — *Active learning in materials science with emphasis on adaptive sampling using uncertainties for targeted design*

> Lookman, T. et al. *npj Comput. Mater.* **5**, 21 (2019).

A pedagogical review of active learning / BO for materials. Frames
the loop you are implementing.

**Why read this.** Best single review for the methodology you are
about to implement. Sections on EGO (efficient global optimisation),
acquisition functions, and the difference between exploration and
exploitation are essential.

**Extract.** The acquisition-function zoo (EI, UCB, PI, EGO); the
notion of *regret*; the trade-off between exploration and
exploitation; the recommendation that one should report multiple
seeds.

---

## 4. Frazier (2018) — *A tutorial on Bayesian optimization*

> Frazier, P. I. arXiv:1807.02811 (2018).

The most readable tutorial on BO. Mathematical but accessible.

**Why read this.** Required reading. You should be able to derive
Expected Improvement from first principles after this paper.

**Extract.** The EI formula in closed form; the qEI extension for
batch acquisition; the optimisation of the acquisition function
itself (a nested optimisation problem); the discussion of GP
hyperparameter learning via maximum likelihood.

---

## 5. Balandat, Karrer, Jiang et al. (2020) — *BoTorch: A Framework for Efficient Monte-Carlo Bayesian Optimization*

> Balandat, M. et al. *NeurIPS* 2020.

The BoTorch paper. You are using this library; this paper is the
reference for it.

**Why read this.** Explains the Monte Carlo formulation of
acquisition functions (qEI, qUCB) and the role of the
sampling-based gradient. Helps you read BoTorch's API documentation.

**Extract.** The MC formulation of qEI; the recommended optimiser
for the acquisition function (`optimize_acqf` with multi-start
L-BFGS); the standard BO loop boilerplate.

---

## 6. Ueno, Rhone, Hou, Mizoguchi, and Tsuda (2016) — *COMBO: An efficient Bayesian optimization library for materials science*

> Ueno, T. et al. *Mater. Discov.* **4**, 18 (2016).

An older but still-relevant materials-focused BO library. Useful for
seeing how BO is framed *outside* the BoTorch ecosystem.

**Why read this.** Helps you avoid the trap of thinking BO = BoTorch.
The methodology is independent of the library; you should be able to
implement it in a few hundred lines.

**Extract.** The use of a Thompson-sampling acquisition for batch
queries; the role of categorical features (atomic species).

---

## 7. Tran and Ulissi (2018) — *Active learning across intermetallics to guide discovery of electrocatalysts for the hydrogen evolution reaction*

> Tran, K. and Ulissi, Z. W. *Nat. Catal.* **1**, 696 (2018).

A real BO/AL study of HER intermetallics. Methodologically the
closest analogue to your project.

**Why read this.** This is what your project is "trying to be". Read
their description of the oracle (a single-site $\Delta G_\mathrm{H}$
DFT calculation on a slab), their featurisation, their acquisition
function, and their regret curves.

**Extract.** The featurisation choices (Magpie features + structural
descriptors); the budget (≈ 800 calculations across many alloys);
the role of uncertainty in the candidate ranking.

---

## 8. Ward, Agrawal, Choudhary, and Wolverton (2016) — *A general-purpose machine learning framework for predicting properties of inorganic materials*

> Ward, L. et al. *npj Comput. Mater.* **2**, 16028 (2016).

The Magpie composition-featurisation paper.

**Why read this.** If you use Magpie features in your GP kernel,
this is the reference. Even if you use simple Pt/Pd/Ag fractions,
the Magpie paper is a useful baseline for what compositional features
*can* encode.

**Extract.** The feature categories (statistics over electronegativity,
atomic number, radii, ...); the standard ML pipeline that Magpie
underlies.

---

## 9. Zhang, Apley, and Chen (2020) — *Bayesian optimization for materials design with mixed quantitative and qualitative variables*

> Zhang, Y. et al. *Sci. Rep.* **10**, 4924 (2020).

Methodological paper on mixed-integer BO — exactly your problem,
because composition counts are integers but BoTorch's acquisition
function optimiser prefers continuous variables.

**Why read this.** Helps you make the round-to-discrete decision
intelligently. Their alternative ("latent variable" embedding of
categorical variables) is more sophisticated than rounding but
overkill for an undergraduate project.

**Extract.** The rounding strategies (nearest, projection,
multi-start with rounding) and their pitfalls.

---

## 10. Hennig and Schuler (2012) — *Entropy search for information-efficient global optimization*

> Hennig, P. and Schuler, C. J.
> *J. Mach. Learn. Res.* **13**, 1809 (2012).

An information-theoretic acquisition function that contrasts nicely
with EI/UCB. Not required for the project, but good intellectual
context.

**Why read this.** Helps you appreciate that EI is *one* choice
among many, and that an information-theoretic objective leads to
different exploration behaviour.

**Extract.** The notion of entropy reduction over the optimum
distribution; the relationship to the predictive variance.

---

## Optional eleventh — for the multi-fidelity extension

If you want to compare DFT and MLIP oracles in a single multi-fidelity
BO loop, see Kandasamy, Dasarathy, Schneider, and Póczos (2017),
*Multi-fidelity Bayesian Optimisation with Continuous Approximations*,
*ICML* 2017.

---

## Synthesis exercise

After reading, write a one-page memo:

1. What is the precise definition of $\Delta G_\mathrm{H}$ you will
   compute, including ZPE + entropy correction?
2. What is your composition-feature representation, and what kernel
   will you use?
3. What is your acquisition function, and why?
4. How will you handle the discrete-composition constraint?
5. How many BO seeds and random-search seeds will you run?

If you cannot answer these, re-read papers 1, 3, and 7 before
implementing.
