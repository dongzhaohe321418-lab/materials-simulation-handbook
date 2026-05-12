# Background reading — High-throughput band-gap screening

Each paper below is the minimum context for a section of your project.
Read in order; the annotations explain what to extract.

---

## 1. Xie and Grossman (2018) — *Crystal Graph Convolutional Neural Networks for an Accurate and Interpretable Prediction of Material Properties*

> Xie, T. and Grossman, J. C. *Phys. Rev. Lett.* **120**, 145301 (2018).

The original CGCNN paper. Introduces the crystal-graph representation
(nodes = atoms, edges = bonds within a cutoff) and the graph
convolution that produces a fixed-length crystal embedding.

**Why read this.** This is your model. Understand the message-passing
equation, the role of the radial bond expansion, and the pooling
operation that produces a per-crystal embedding.

**Extract.** The architecture diagram; the form of the per-node
update; the choice of feature dimensions (typically 64); the
hyperparameter values reported (number of conv layers, hidden
dimension, learning rate); the reported MAE on the Materials Project
band-gap task (≈ 0.39 eV at the time of publication).

---

## 2. Jain, Ong, Hautier et al. (2013) — *Commentary: The Materials Project: A materials genome approach to accelerating materials innovation*

> Jain, A. et al. *APL Materials* **1**, 011002 (2013).

The Materials Project announcement. Describes the calculation
workflow, the database schema, and the design principles.

**Why read this.** You are downloading thousands of records from MP;
you should know what they are. In particular, the difference between
GGA, GGA+U, and HSE entries; the meaning of `energy_above_hull`; and
the role of the materials taxonomy.

**Extract.** The standard PBE/PAW calculation parameters (`ENCUT`
520 eV, `EDIFF` 1e-5, default Monkhorst–Pack); the documented
limitations (no van der Waals correction, no spin-orbit unless
flagged); the distinction between `band_gap` (Kohn–Sham gap from a
band-structure calculation) and `formation_energy_per_atom`.

---

## 3. Zhuo, Mansouri Tehrani, and Brgoch (2018) — *Predicting the Band Gaps of Inorganic Solids by Machine Learning*

> Zhuo, Y., Mansouri Tehrani, A., and Brgoch, J. *J. Phys. Chem. Lett.*
> **9**, 1668 (2018).

A composition-only (Magpie-style) baseline for band-gap prediction.
Less accurate than a structure-aware GNN, but methodologically
informative.

**Why read this.** Gives you a baseline number to beat. Their reported
MAE was around 0.45–0.50 eV with a Magpie + Random Forest pipeline.
Your CGCNN should be similar or better.

**Extract.** The composition-only feature set; the importance ranking
of features (electronegativity, mean atomic number, ...); the
discussion of which materials classes are hardest to predict.

---

## 4. Castelli, Olsen, Datta et al. (2012) — *Computational screening of perovskite metal oxides for optimal solar light capture*

> Castelli, I. E. et al. *Energy Environ. Sci.* **5**, 5814 (2012).

The original computational photocatalyst screen for cubic perovskites.
Used DFT directly on a few thousand candidates to predict band gaps
and band-edge positions relative to the water-splitting potentials.

**Why read this.** Frames the photocatalyst-screening problem at the
level of *thermodynamics* — not just "is the gap in 1.5–2.5 eV" but
"are the band edges correctly aligned for water splitting". You should
at least mention these complications in your report, even if you do
not address them.

**Extract.** The 1.5–3.0 eV "useful" gap range; the band-edge
alignment criteria (VBM above O$_2$/H$_2$O, CBM below H$^+$/H$_2$);
the typical fraction of stable candidates (only ≈ 10 % of nominally
stable perovskites end up satisfying all criteria).

---

## 5. Chen, Ye, Zuo, Zheng, and Ong (2019) — *Graph Networks as a Universal Machine Learning Framework for Molecules and Crystals*

> Chen, C., Ye, W., Zuo, Y., Zheng, C., and Ong, S. P.
> *Chem. Mater.* **31**, 3564 (2019).

The MEGNet paper. Improves on CGCNN by adding a global-state node
and explicit edge updates. Often beats CGCNN by 5–10 % on the same
tasks.

**Why read this.** A natural alternative architecture, and a good
point of comparison. You will not be required to implement MEGNet,
but knowing it exists prevents you from over-claiming CGCNN's
performance.

**Extract.** The role of the global state; the recommended
hyperparameters; the comparison table against CGCNN on Materials
Project band gaps.

---

## 6. Choudhary and DeCost (2021) — *Atomistic Line Graph Neural Network for improved materials property predictions*

> Choudhary, K. and DeCost, B. *npj Comput. Mater.* **7**, 185 (2021).

ALIGNN, which adds an explicit line-graph (bond-angle) representation
to the GNN. Beats CGCNN and MEGNet on most tasks by another 10–20 %.

**Why read this.** Same reason as MEGNet — a more accurate
alternative. If you have time you may swap in ALIGNN; otherwise cite
it as a point of comparison.

**Extract.** The line-graph construction; the typical training set
size where ALIGNN starts to outperform CGCNN (above ~ 10 000 samples).

---

## 7. Tran, Stein, Ulissi (2018) — *Active learning across intermetallics to guide discovery of electrocatalysts for the hydrogen evolution reaction*

> Tran, K. and Ulissi, Z. W. *Nat. Catal.* **1**, 696 (2018).

Although focused on catalysts rather than photocatalysts, this paper
demonstrates the *uncertainty-aware screening* approach you will use:
prediction + uncertainty, then DFT verification only of the
most-promising-with-low-uncertainty candidates.

**Why read this.** You will use an ensemble of CGCNNs to estimate
uncertainty. This paper formalises why that is the right thing to do.

**Extract.** The interplay between prediction mean and standard
deviation in candidate ranking; the discussion of why "predicted
best" candidates with high uncertainty are *not* the best DFT
candidates.

---

## 8. Stanev, Oses, Kusne, Rodriguez, Paglione, Curtarolo, and Takeuchi (2018) — *Machine learning modeling of superconducting critical temperature*

> Stanev, V. et al. *npj Comput. Mater.* **4**, 29 (2018).

A different target ($T_c$ of superconductors) but methodologically
similar: train on a large clean dataset, predict on a candidate set,
verify a subset. Useful for the *stratified-evaluation* methodology
the authors use.

**Why read this.** Their honest discussion of how a model with
0.3 K MAE on $T_c$ in the global average has 5 K MAE in the
high-$T_c$ regime is exactly the lesson you must take to your own
project: average errors lie.

**Extract.** The stratified-MAE table; the discussion of selection
bias in the candidate set.

---

## 9. Schmidt, Hoffmann, Wang, Borlido, Carriço, Cerqueira, Botti, and Marques (2019) — *Predicting the Thermodynamic Stability of Solids Combining Density Functional Theory and Machine Learning*

> Schmidt, J. et al. *Chem. Mater.* **29**, 5090 (2019). [Note: 2017 paper.]

ML prediction of formation energies / stability, which you will use
as a stability filter on your candidates.

**Why read this.** Stability is not optional — a predicted 2-eV-gap
oxide that decomposes spontaneously is not a candidate, it is a
mistake. This paper teaches you the standard ways to estimate
stability without an explicit hull calculation.

**Extract.** The use of `energy_above_hull` from MP as a stability
filter; the typical fraction of "stable enough" entries (≈ 60 % of MP
is below 0.1 eV/atom); the discussion of why hull computation is
nontrivial.

---

## 10. Pilania, Mannodi-Kanakkithodi, Uberuaga, Ramprasad, Gubernatis, and Lookman (2016) — *Machine learning bandgaps of double perovskites*

> Pilania, G. et al. *Sci. Rep.* **6**, 19375 (2016).

A focused, well-controlled screening study: train on 600 known double
perovskites, predict on a much larger combinatorial space,
verify with DFT. Identifies several previously unknown candidates.

**Why read this.** A real example of the project workflow you are
attempting. Their final shortlist was small (10–20 candidates) and
their verification rate was high. Read their honest discussion of
how many predicted candidates turned out to be metastable but
synthesisable.

**Extract.** The screening pipeline; the rate of successful DFT
verification; the literature follow-up they performed.

---

## Optional eleventh — for the band-gap-correction step

If you want to apply a simple PBE → HSE correction to your CGCNN
predictions, see Morales-García, Valero, and Illas (2017),
*J. Phys. Chem. C* **121**, 18862. A linear correction of the form
$E_g^\mathrm{HSE} \approx 1.3\, E_g^\mathrm{PBE} + 0.4$ eV is a
crude but useful first approximation; better corrections from the
literature improve it modestly.

---

## Synthesis exercise

After reading, write a one-page memo answering:

1. What is the typical accuracy of a CGCNN on PBE band gaps from
   Materials Project? On HSE gaps?
2. Why is *stratified MAE* important when screening for a specific
   gap window?
3. How does PBE band-gap underestimation affect the apparent
   "photocatalyst window" you should target?
4. What stability filter will you apply, and why?
5. How will you estimate uncertainty for your predictions, and how
   will it enter the candidate-ranking step?

If you cannot answer these confidently, re-read papers 1, 4, and 7
before training anything.
