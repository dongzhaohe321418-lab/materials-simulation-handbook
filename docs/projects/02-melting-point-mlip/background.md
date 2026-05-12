# Background reading — Melting point of copper via MLIP-driven MD

Read these in roughly the order given. The annotation explains how
each paper feeds into your project, and what to record before moving
on.

---

## 1. Foiles, Baskes, and Daw (1986) — *Embedded-atom method functions for the fcc metals*

> Foiles, S. M., Baskes, M. I., and Daw, M. S. *Phys. Rev. B* **33**,
> 7983 (1986).

The classical EAM potential for Cu (and other fcc metals) that
remains, after forty years, the workhorse baseline for any classical
MD on copper.

**Why read this.** You will use the Foiles-Cu potential as a control:
both as a sanity check on your LAMMPS workflow before MACE is ready,
and as a baseline against which to compare the MLIP's $T_m$ prediction
(EAM-Foiles for Cu gives $T_m \approx 1325$ K — close to experiment
by virtue of fitting).

**Extract.** The form of the EAM functional ($E = \sum_i F(\bar\rho_i)
+ \tfrac{1}{2} \sum_{ij} \phi(r_{ij})$), the fitting targets (lattice
parameter, cohesive energy, bulk modulus, vacancy formation energy),
and the published melting point.

---

## 2. Morris, Wang, Ho, and Chan (1994) — *Melting line of aluminum from simulations of coexisting phases*

> Morris, J. R., Wang, C. Z., Ho, K. M., and Chan, C. T.
> *Phys. Rev. B* **49**, 3109 (1994).

The classic paper introducing the two-phase coexistence method for
extracting melting points from MD. They studied Al, but the
methodology is identical for Cu.

**Why read this.** This is *the* paper that describes your primary
$T_m$ estimator. Read it carefully; the cell geometry, the equilibration
protocol, and the interface-velocity diagnostic all come from here.

**Extract.** The construction of the half-solid–half-liquid cell with
periodic boundary conditions, the role of NPT vs NPH ensembles, the
definition of interface velocity, and the importance of running at
multiple temperatures and bracketing the $T_m$ root.

---

## 3. Belonoshko, Ahuja, and Johansson (2000) — *Quasi-Ab Initio Molecular Dynamic Study of Fe Melting*

> Belonoshko, A. B., Ahuja, R., and Johansson, B.
> *Phys. Rev. Lett.* **84**, 3638 (2000).

Demonstrates the two-phase method for a transition metal at finite
pressure. Establishes that ≈ 1000–10000-atom cells are needed and that
nucleation barriers in pure-phase cells can yield hysteresis of 30 % in
$T_m$.

**Why read this.** Justifies your choice of cell size and your
decision to use coexistence rather than hysteresis as the primary
estimator.

**Extract.** The cell-size convergence study (their figure 2), the
quantitative magnitude of superheating in pure-solid simulations, and
the pressure-corrected melting curves.

---

## 4. Batatia, Kovács, Simm, Ortner, and Csányi (2022) — *MACE: Higher Order Equivariant Message Passing Neural Networks*

> Batatia, I., Kovács, D. P., Simm, G., Ortner, C., and Csányi, G.
> *NeurIPS* 2022.

The MACE architecture paper. Introduces the equivariant message-passing
formalism with higher-body-order features (typically up to body order
4) and the spherical-harmonic-based irreducible representations.

**Why read this.** You are using this exact model. Understand at
minimum: the meaning of `max_L` (the maximum angular momentum of the
representation), the role of body order, the difference between hidden
features and node features, and the rationale for choosing energy +
force training over force-only training.

**Extract.** The architectural overview (the message-passing diagram),
the typical hyperparameter ranges (e.g., `num_interactions=2`, `max_L=2`,
`hidden_irreps='128x0e+128x1o'`), and the training-loss decomposition.

---

## 5. Behler (2015) — *Constructing high-dimensional neural network potentials: A tutorial review*

> Behler, J. *Int. J. Quantum Chem.* **115**, 1032 (2015).

A practical guide to constructing MLIP training sets, with emphasis on
sampling, descriptor choice, and validation. Predates MACE but the
philosophy applies verbatim.

**Why read this.** Best single source for "how do I build a good
training set". The discussion of active learning, of force-vs-energy
weighting, and of out-of-distribution detection are all relevant.

**Extract.** The recommendation to sample several thermodynamic
states explicitly; the rule of thumb that the training-set force MAE
should be of order 5 % of the typical force magnitude; the warning
against under-sampling unusual configurations (e.g., the liquid).

---

## 6. Bartók, Kondor, and Csányi (2013) — *On representing chemical environments*

> Bartók, A. P., Kondor, R., and Csányi, G.
> *Phys. Rev. B* **87**, 184115 (2013).

The SOAP descriptor paper. Even if you do not use SOAP explicitly,
MACE's invariants can be thought of as a structured generalisation,
and this paper teaches you the language of "atomic environments".

**Why read this.** Develops your intuition for why nearest-neighbour
environments encode most of the relevant local physics. Useful when
you later need to decide on the receptive-field radius (typically
4–5 Å for fcc metals).

**Extract.** The notion of an atomic environment as a smooth density;
the cutoff function; the invariance under permutation, rotation, and
translation.

---

## 7. Vega, Sanz, Abascal, and Noya (2008) — *Determination of phase diagrams via computer simulation*

> Vega, C., Sanz, E., Abascal, J. L. F., and Noya, E. G.
> *J. Phys.: Condens. Matter* **20**, 153101 (2008).

A pedagogical review of free-energy and direct-coexistence methods
for solid–liquid equilibria. Compares Gibbs–Duhem integration,
hysteresis, thermodynamic integration, and direct coexistence.

**Why read this.** Lets you justify the choice of two-phase
coexistence as the most robust method for a single-component system,
and gives you the analytic forms for the Clausius–Clapeyron relation
should you wish to estimate $dT_m/dp$.

**Extract.** The hierarchy of methods (free energies > coexistence >
hysteresis); the sample-size requirements; the role of finite-size
corrections.

---

## 8. Pozzo and Alfè (2013) — *Melting curve of face-centred-cubic nickel from first-principles calculations*

> Pozzo, M. and Alfè, D. *Phys. Rev. B* **88**, 024111 (2013).

A close-relative study (Ni instead of Cu) using DFT-MD with
direct coexistence. Establishes the practical workflow at the DFT
level — important because you are doing the same workflow at the MLIP
level.

**Why read this.** Quantifies the convergence of the DFT $T_m$ with
respect to cell size, **k**-mesh, and pseudopotential. Numbers are
similar for Cu.

**Extract.** Recommended cell sizes for direct coexistence; typical
runtime of 30–50 ps per temperature; the magnitude of finite-size
$T_m$ corrections.

---

## 9. Bonati and Parrinello (2018) — *Silicon Liquid Structure and Crystallization from an Ab Initio Deep Potential*

> Bonati, L. and Parrinello, M. *Phys. Rev. Lett.* **121**, 265701
> (2018).

An MLIP-driven coexistence simulation for Si. Methodologically the
closest precedent to what you are doing, although for a different
element and with a different MLIP architecture (Deep Potential).

**Why read this.** Best published example of an MLIP-derived
melting-point prediction. Their final value matched experiment to
≈ 50 K with PBE training data.

**Extract.** The training-set composition (≈ 800 frames spanning
solid, liquid, and intermediate temperatures), the parity-error
targets, and the comparison against experiment.

---

## 10. Kovács, Batatia, Arany, and Csányi (2023) — *Evaluation of the MACE Force Field Architecture: From Medicinal Chemistry to Materials Science*

> Kovács, D. P., Batatia, I., Arany, E. S., and Csányi, G.
> *J. Chem. Phys.* **159**, 044118 (2023).

A practical "how to use MACE" paper with extensive benchmarks and
default-hyperparameter recommendations. Most useful section is on
training-data size and hyperparameter sensitivity.

**Why read this.** Tells you what hyperparameter values *not* to
spend two weeks tuning. The defaults are sane.

**Extract.** The sensitivity of force MAE to `max_L` and
`hidden_irreps` (modest); the sensitivity to training-set size
(strong below 500 frames, weak above 2000); the typical training-curve
shape.

---

## Optional eleventh — for the EAM control run

Mishin, Mehl, Papaconstantopoulos, Voter, and Kress (2001) on
modified-EAM Cu (*Phys. Rev. B* **63**, 224106) is a more modern EAM
that beats Foiles–Baskes–Daw on most properties. If you want a
fairer EAM baseline, use this one.

---

## Synthesis exercise

After reading, write a one-page memo answering:

1. What is the difference between hysteresis and coexistence as
   melting-point estimators, and which gives the *true* thermodynamic
   $T_m$?
2. What thermodynamic states must your DFT training set cover, and
   why? What would go wrong if you only included 300 K configurations?
3. What MACE hyperparameters will you fix at defaults, and which one
   (if any) will you tune?
4. What is the expected force MAE on a held-out test set for a
   well-trained Cu MACE, and what does a 50 meV/Å force MAE mean
   physically?

If you cannot write the memo confidently, re-read papers 2, 4, and 5
before generating data.
