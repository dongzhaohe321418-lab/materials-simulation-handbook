# Chapter 10 — Graph Neural Networks for Materials

[Open in Colab](https://colab.research.google.com/github/dongzhaohe321418-lab/materials-simulation-handbook/blob/main/notebooks/ch10-gnn-colab.ipynb){ .md-button .md-button--primary }

Chapter 9 closed with a quiet observation. A machine-learning interatomic
potential, stripped of its symmetry plumbing, is a regression model that
takes a local atomic neighbourhood and returns a scalar — the atomic
energy. Sum over atoms, and you have the total energy. Differentiate, and
you have the forces. The architecture choices that distinguish SchNet from
MACE are details of *how* the neighbourhood is encoded; the conceptual
content is the same.

That observation, taken seriously, leads naturally to this chapter. If we
already think of a crystal as a collection of atoms each surrounded by
neighbours within a cutoff, we are halfway to a graph. Promote the atoms
to nodes, the neighbour relations to edges, and the whole apparatus of
*graph neural networks* (GNNs) becomes available. We can then ask the same
network to predict any scalar property of the crystal — formation energy,
band gap, bulk modulus, even thermodynamic stability — not just the
potential energy surface that drives molecular dynamics.

This shift in viewpoint matters for two reasons.

First, *property prediction* is the bread-and-butter task of high-throughput
materials informatics. We have databases — Materials Project, OQMD, AFLOW,
JARVIS — each containing hundreds of thousands of DFT-computed structures
labelled with formation energies, gaps, magnetic moments and elastic
constants. We would like a model that, given a new structure, predicts
these properties in milliseconds rather than the hours a DFT calculation
would take. GNNs are the dominant architecture for this task, and a single
trained CGCNN or ALIGNN can screen a million candidate compounds overnight
on a single GPU.

Second, the graph viewpoint generalises. Sections of later chapters will
treat heterogeneous graphs (atoms plus voids, atoms plus surface sites)
and line graphs (edges become nodes, bond angles become edges). The same
message-passing machinery handles them all. Understanding GNNs deeply
once means you can read most of the modern materials-ML literature.

The chapter proceeds as follows. Section 10.1 establishes the dictionary
between crystals and graphs. We define $G = (V, E)$ carefully, attach
features to nodes (one-hot or learned element embeddings) and edges
(interatomic distances expanded in a Gaussian basis), and confront the
subtle point that periodic boundary conditions force a single atom to
connect to many of its own periodic images. Code snippets construct a
graph from a `pymatgen.Structure` using either `CrystalNN` or a plain
cutoff.

Section 10.2 introduces the message-passing neural network (MPNN)
framework that underpins essentially every modern GNN. We write down the
abstract update equations
$$
m_v^{(t+1)} = \sum_{u \in \mathcal{N}(v)} M_t(h_v^{(t)}, h_u^{(t)}, e_{uv}),
\qquad
h_v^{(t+1)} = U_t(h_v^{(t)}, m_v^{(t+1)}),
$$
explain why summing over neighbours yields *permutation invariance for
free*, and confront the practical pathologies — finite receptive field,
over-smoothing — that limit how deep these networks can be.

Section 10.3 builds CGCNN, the Crystal Graph Convolutional Neural Network
of Xie and Grossman (2018), from scratch. It is a deliberately simple
model — one of the first to treat crystals as graphs — and that simplicity
makes it the right pedagogical target. We implement the gated convolution,
wire it into a full PyTorch Geometric pipeline complete with dataset
class, training loop and evaluation, and verify it on fifty structures
pulled from the Materials Project.

Section 10.4 surveys what came next. MEGNet added a global *state*
attribute, permitting temperature- or pressure-conditioned predictions.
ALIGNN added a *line graph* so that bond angles enter explicitly, and it
remains the top performer on most static-property benchmarks. M3GNet
took the same line-graph idea and trained on the full Materials Project
relaxation trajectory set, producing the first credible universal
interatomic potential. We discuss honestly what each model is good for —
ALIGNN for property regression, M3GNet (and its successors CHGNet, MACE-MP)
for dynamics — rather than promoting a single architecture.

Section 10.5 walks through a full Materials Project pipeline. You will
register for an API key, query five thousand oxide structures via
`mp-api`, deduplicate them, split into in- and out-of-domain test sets,
train CGCNN, and produce a parity plot. The exercise highlights a
mistake that pervades the literature: random train/test splits over
databases that contain many near-duplicate polymorphs systematically
overestimate model accuracy.

By the end of the chapter you will have written, trained, and evaluated
a complete graph neural network for materials property prediction, and
you will understand the design space well enough to read — and modify —
any of the modern architectures. Chapter 11 will then turn the question
around: given such a fast surrogate model, how do we use it intelligently
to guide expensive experiments? That is the territory of Bayesian
optimisation and active learning.

A final pointer. Chapter 9 treated the MLIP problem as one of fitting an
energy surface; this chapter treats the more general problem of fitting
any scalar property of a crystal. The two are not separate worlds. An
MLIP is exactly a GNN whose output is the energy and whose training loss
includes a gradient term for the forces. Chapter 12 will return to this
unification when we examine *foundation models* for materials, where one
network is pre-trained on all available data and fine-tuned for whatever
property the user happens to care about.

## Chapter-end summary

By the end of Chapter 10 the reader should be able to:

1. **State the formal definition** of a crystal graph $G = (V, E)$ as a
   directed multigraph with self-loops, with node features encoding
   chemistry and edge features encoding Gaussian-expanded distances
   (§10.1).
2. **Derive permutation invariance** of the MPNN framework from the
   commutativity of the neighbour-sum, and **diagnose over-smoothing**
   via the spectral-contraction argument on the normalised adjacency
   (§10.2).
3. **Implement CGCNN** end-to-end in PyTorch Geometric — gated message
   function, residual update with batch normalisation, mean-pool
   readout, MLP head — and reproduce published Materials Project MAE
   numbers within a factor of two on a held-out test set (§10.3).
4. **Choose between architectures** (CGCNN / MEGNet / ALIGNN / M3GNet)
   using the §10.4.5a decision tree, with a clear sense of when angular
   information, state attributes or forces are needed.
5. **Run an honest Materials Project pipeline**: API setup, query
   construction, deduplication with `StructureMatcher`, composition-
   disjoint splitting, training, parity-plot evaluation, and ensemble
   uncertainty quantification (§10.5).

The single most important takeaway: every GNN architecture in the
materials literature is a specific point in the five-dimensional design
space of (message function, aggregation, update, edge update, readout).
Reading any new GNN paper, locate it in this space and the rest
follows.
