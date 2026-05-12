# Chapter 12 — Exercises

Six exercises follow, mixing conceptual questions and hands-on work.
Exercises that require a GPU are marked **[GPU]**. Solutions are at
the end of the file.

---

## Exercise 12.1 — The vocabulary argument (conceptual)

Section 12.1 argued that the periodic table's "small alphabet" is one
of the reasons the foundation-model paradigm transfers from NLP to
materials.

(a) Suppose we tokenise an atomistic configuration by treating each
$(Z_i, \text{local environment})$ pair as a token. Estimate the
*effective vocabulary size* of this tokenisation for a typical
inorganic compound, taking into account the number of distinct
elements that appear and the diversity of local environments. How
does this compare to the vocabulary of a modern BPE-tokenised
language model ($\sim 5 \times 10^4$)?

(b) The argument also relied on *shared local physics*. Give a
counterexample — a setting in which the chemistry of an element in
one compound is genuinely not a useful predictor of its chemistry
in another. (Hint: think about oxidation states or about elements
appearing in radically different coordination environments.)

---

## Exercise 12.2 — Cost–accuracy curves (conceptual / numerical)

The fine-tuning table in Section 12.2 reported the following errors
on a perovskite oxide:

| $N_\text{ft}$ | $\text{MAE}(F)$ / meV·Å$^{-1}$ |
|---|---|
| $25$ | $52$ |
| $100$ | $28$ |
| $500$ | $14$ |
| $2{,}000$ | $9$ |

(a) Assume the empirical scaling $\text{MAE}(F) \approx A N^{-\alpha} + B$.
Estimate $A$, $\alpha$ and $B$ by fitting in log-log space (you may
ignore $B$ to a first approximation; treat the four points as a
log-log line). What is the implied exponent $\alpha$?

(b) Using your fit, predict the fine-tuning size required to reach
$\text{MAE}(F) = 5$ meV/Å. Comment on whether the prediction is
trustworthy.

(c) The asymptote $B$ corresponds to an *irreducible error floor*.
What physical phenomena might set this floor for a system the
foundation model was not trained on?

---

## Exercise 12.3 — Zero-shot MD diagnostic **[GPU]**

Using the script in Section 12.2, run a $100$-step Langevin MD
trajectory for an MgO supercell at $600$ K with MACE-MP-0 (medium).

(a) Plot the potential energy as a function of step. Comment on
whether the trajectory has equilibrated.

(b) Sample five configurations from the trajectory (e.g., steps
$20, 40, 60, 80, 100$). For each, run a single-point DFT calculation
in your favourite code (Quantum ESPRESSO, VASP, GPAW, …). Compute the
mean absolute energy difference per atom between MACE and DFT, and
the mean absolute force difference per atom.

(c) Repeat with the *large* MACE-MP-0 variant. Does the error
improve enough to justify the increased inference cost?

---

## Exercise 12.4 — Fine-tuning a foundation MLIP **[GPU]**

Take a system of your choice for which you have $\sim 100$
DFT-relaxed configurations with energies and forces. (Materials
Project relaxation trajectories for a single material, downloaded
through `mp-api`, are a convenient choice. So is any reaction-path
dataset from a previous Chapter 9 exercise.)

(a) Split your dataset into $80$%–$20$% train/validation. Evaluate
zero-shot MACE-MP-0 on the validation set; record the energy and
force MAEs.

(b) Fine-tune MACE-MP-0 using the recipe of Section 12.2: freeze
the lower layers, use a forces-weighted loss, train for $\sim 200$
epochs. Plot validation loss vs epoch.

(c) Repeat with $25$, $50$ and the full $\sim 80$ training
configurations. Reproduce the cost–accuracy curve for your system
and compare with Exercise 12.2.

(d) Pick one validation configuration and inspect, visually, where
the largest force errors occur (e.g., on which atoms, in which
local environment). What does this tell you about the model's
remaining weaknesses?

---

## Exercise 12.5 — Generative validation pipeline (conceptual)

Section 12.3 sketched a six-stage validation pipeline for outputs
from a generative model: generate, sanity-filter, MLIP-relax, DFT-
verify, synthesisability screen, experiment.

(a) For each of the first five stages, identify the *false negative*
mode (i.e., the way in which a genuinely interesting candidate could
be incorrectly rejected). Suggest a remedy.

(b) For each of the first five stages, identify the *false positive*
mode (a candidate that passes but should not). Suggest a remedy.

(c) Suppose the generative model is conditioned on a target band gap
of $2.5$ eV, and the model's training predicted band gaps from
PBE-DFT (known to underestimate gaps by a factor of $\sim 1.5$). What
target should you actually condition on if the goal is to find a real
material with a measured gap of $2.5$ eV?

---

## Exercise 12.6 — A reading exercise (conceptual)

Choose one paper from the reading list in Section 12.4. Read it
carefully and write a one-page summary covering:

(a) the central technical contribution (one paragraph);

(b) the experimental or numerical evidence supporting it (one
paragraph);

(c) the limitations the authors acknowledge, and any further
limitations *you* believe the paper does not adequately address (one
paragraph);

(d) one concrete follow-up experiment or extension you would propose.

This exercise has no "solution" below: it is for your own use.

---

## Solutions

### 12.1

(a) An inorganic compound typically contains $\sim 3$–$5$ distinct
elements. Each element appears in a small number of distinct local
environments (say, $\sim 5$–$10$ coordination motifs across all
crystals the model might encounter). The product gives an effective
vocabulary of perhaps $96 \times 10 \approx 10^3$ "tokens", which is
indeed two orders of magnitude smaller than the BPE vocabulary of a
language model.

The right interpretation is not that one *literally* tokenises this
way (one does not; the model learns continuous local descriptors).
The conceptual point is that the *underlying entropy of the data* is
much lower than for natural language, and a model with comparable
capacity is therefore in a much better position relative to the size
of the space it must cover.

(b) Lanthanides provide the canonical counterexample. The chemistry
of Eu$^{2+}$ (large divalent cation, ionic, $7s$-like behaviour
through screening) is essentially unrelated to that of Eu$^{3+}$
(strong magnetic character, $4f$ states sharply localised, important
in optical applications). A foundation model that has only seen
Eu$^{2+}$ environments will not perform well on Eu$^{3+}$. Similar
remarks apply to early transition metals in unusual oxidation states
(Ti$^{2+}$ vs Ti$^{4+}$), and to actinides generally.

### 12.2

(a) Working in log-log space and ignoring $B$ for the moment, the
four points sit approximately on a line. Taking the slope from
$(N, F) = (25, 52)$ to $(2000, 9)$:
$$
\alpha = -\frac{\log(9/52)}{\log(2000/25)}
       = \frac{\log 5.78}{\log 80}
       \approx \frac{0.762}{1.903}
       \approx 0.40.
$$
The prefactor follows from $A = 52 \times 25^{0.40} \approx 192$
meV/Å. So $\text{MAE}(F) \approx 190 \cdot N^{-0.40}$ meV/Å.

(b) Setting $\text{MAE}(F) = 5$:
$$
5 = 190 \cdot N^{-0.40}
\quad\Rightarrow\quad
N = (190/5)^{1/0.40} = 38^{2.5} \approx 9 \times 10^3.
$$
The prediction should be treated with caution. Empirical scaling
laws in this regime typically *break* before the predicted target is
reached, because the irreducible floor $B$ becomes dominant. A more
realistic estimate is that $\text{MAE}(F) = 5$ meV/Å is not
reachable by fine-tuning alone for this system, and approaching it
would require either training-distribution-matched data or moving to
a larger model variant.

(c) The floor can be set by: (i) the gap between the level of theory
of the foundation training data (PBE+U) and the level used to label
your fine-tuning data (e.g., a hybrid functional); (ii) intrinsic
inability of the model architecture to represent some long-range or
many-body effect relevant to your system; (iii) noise in the
fine-tuning labels themselves, especially in forces, where DFT
convergence at $\sim 1$ meV/Å is not always achieved.

### 12.3

(a) On a $108$-atom MgO supercell, the trajectory typically takes
$\sim 30$ steps to lose memory of the initial conditions and reach
an approximately stationary distribution of $E$ and $T$.
Equilibration is indicated by stabilisation of the running mean of
the potential energy and of the temperature around the target. The
chosen friction ($0.01$ fs$^{-1}$) is light; longer trajectories
than $100$ steps will produce cleaner statistics.

(b) For MgO at $600$ K, zero-shot MACE-MP-0 medium typically
delivers $\text{MAE}(E) \sim 3$–$5$ meV/atom and $\text{MAE}(F)
\sim 30$–$60$ meV/Å against PBE single points, depending on which
configurations are sampled. The largest errors usually appear on
oxygen forces near transient distortions of the rocksalt structure.

(c) The large model typically reduces force errors by $\sim 30$–$40$%
on this system, at a cost of roughly $2.5\times$ inference time.
For production MD where the bottleneck is wallclock time, the medium
model is usually the better choice; the large model is preferable
when one is post-hoc evaluating snapshots from another trajectory.

### 12.4

This is a hands-on exercise; the precise numbers will depend on the
chosen system. Two general patterns to expect.

(b) The validation loss should drop rapidly during the first $\sim 50$
epochs and then plateau. If it continues to drop after $200$ epochs,
the training is probably underfit; if it drops then rises, the model
is overfitting and the patience parameter should be reduced.

(c) The cost–accuracy curve should produce an exponent $\alpha$ in
the range $0.3$–$0.5$, consistent with Exercise 12.2.

(d) Common pathologies include: large errors on atoms near defects
or surfaces; systematic underestimation of forces in
high-coordination environments not well represented in MPtrj;
unexpectedly poor performance for low-frequency vibrational modes
where the relevant forces are small in absolute magnitude but
important for thermodynamics.

### 12.5

(a) False negatives by stage:

| Stage | False-negative mode | Remedy |
|---|---|---|
| Generate | Diversity collapse, target not represented in latent | Multiple guidance scales; ensemble of models |
| Sanity filter | Aggressive bond-length cutoff rejects unusual but real chemistry | Use compound-specific covalent-radii ranges; soft thresholds |
| MLIP relax | Foundation MLIP unfamiliar with the candidate's chemistry, predicts a spurious instability | Fine-tune the MLIP on the closest known chemistry first |
| DFT verify | Local-minimum convergence on a saddle point or a higher-energy polymorph | Use multiple starting configurations; check phonons |
| Synthesisability | Heuristic too tied to known precursor routes; rejects genuine novelties | Combine with literature-mined precedent searches |

(b) False positives by stage:

| Stage | False-positive mode | Remedy |
|---|---|---|
| Generate | Conditioning satisfied in expectation, not in any particular sample | Sample-by-sample re-evaluation of target property |
| Sanity filter | Heuristic too generous, lets through valence-violating structures | Add explicit oxidation-state checks |
| MLIP relax | MLIP misjudges stability of an exotic structure | Cross-check with multiple universal MLIPs |
| DFT verify | PBE under- or over-estimates the target property | Use higher-level theory ($GW$, hybrid) for the property of interest |
| Synthesisability | Match to a known precursor route does not guarantee successful synthesis under realistic conditions | Send to actual experiment, accept some failures |

(c) PBE underestimates gaps by roughly $30$–$50$% in most cases. A
measured gap of $2.5$ eV typically corresponds to a PBE gap of
$\sim 1.5$–$1.7$ eV. If the generative model was conditioned on a
PBE-trained property predictor, one should target a PBE gap of about
$1.6$ eV to land near a measured $2.5$ eV. (More carefully: use a
trained PBE-to-experiment correction, of which several are published
for oxides and chalcogenides.)
