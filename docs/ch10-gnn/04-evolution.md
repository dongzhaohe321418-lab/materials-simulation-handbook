# 10.4 Evolution: MEGNet, ALIGNN, M3GNet

CGCNN appeared in 2018 and demonstrated the basic point: graph neural
networks reach DFT-level accuracy on bulk property prediction at orders
of magnitude lower inference cost. Five years of refinements then
followed, each motivated by a specific architectural shortcoming. This
section traces the line through three landmark architectures —
MEGNet (2019), ALIGNN (2021) and M3GNet (2022) — and is candid about
what each is and is not good for.

!!! note "Key Idea (Box 10.4.A)"
    Each post-CGCNN architecture identifies a specific shortcoming and
    fixes it. MEGNet adds a *global state* for conditioning on external
    variables (T, P). ALIGNN adds a *line graph* to inject bond
    angles. M3GNet adds *three-body terms* and trains on the full MP
    relaxation trajectory set, becoming the first credible universal
    MLIP. Each fix costs compute and parameters; choose by the
    physics of your target.

## 10.4.1 MEGNet: state attributes

CGCNN as we built it has no mechanism to condition predictions on
external variables. The temperature, pressure or applied field at which
a property is measured is simply absent from the model. For static
ground-state properties — formation energy, band gap at zero
temperature, equilibrium lattice constants — this is fine. For
properties that depend on thermodynamic state — finite-temperature heat
capacity, thermal conductivity at 300 K, pressure-induced phase
transitions — it is a fatal omission.

Chen, Ye, Zuo, Zheng and Ong (2019) addressed this with MEGNet
(MatErials Graph Network). The architecture is a CGCNN with an added
*global state* vector $s$ that participates in every message-passing
operation. The state vector is, in the simplest case, a one-dimensional
container for "the temperature at which this property was measured" or
"the magnetisation of this configuration"; it can be multidimensional if
several external variables are relevant.

The message-passing update in MEGNet has three stages, run in sequence
at each layer.

**Edge update.** Each edge feature is refined using the current node
states and the global state:
$$
e_{uv}^{(t+1)} = \phi_e\!\left( e_{uv}^{(t)},\, h_u^{(t)},\, h_v^{(t)},\, s^{(t)} \right).
$$

**Node update.** The classical message passing, but with the global
state appended:
$$
h_v^{(t+1)} = \phi_h\!\left( h_v^{(t)},\, \sum_{u \in \mathcal{N}(v)} e_{uv}^{(t+1)},\, s^{(t)} \right).
$$

**State update.** The state itself evolves by aggregating over the
whole graph:
$$
s^{(t+1)} = \phi_s\!\left( s^{(t)},\, \frac{1}{|E|}\sum_{(u,v) \in E} e_{uv}^{(t+1)},\, \frac{1}{|V|}\sum_{v \in V} h_v^{(t+1)} \right).
$$
The functions $\phi_e, \phi_h, \phi_s$ are MLPs.

The architectural lesson is general: every level — edge, node, graph —
can carry its own state and be updated. CGCNN updates only nodes;
MEGNet updates all three. The cost is more parameters and slower
forward passes; the benefit is the ability to make state-conditioned
predictions.

!!! note "Derivation: how the state vector enters a single message"
    Spelling out the MEGNet edge update at the level of inputs and
    outputs makes the difference with CGCNN concrete. CGCNN's message
    function is $M(h_v, h_u, e_{uv})$ — three inputs. MEGNet's edge
    update is
    $$
    e_{uv}^{(t+1)} = \phi_e\!\left([\,h_u^{(t)};\, h_v^{(t)};\, e_{uv}^{(t)};\, s^{(t)}\,]\right),
    $$
    where the concatenated input has dimension $2d + d_e + d_s$ rather
    than $2d + d_e$ (CGCNN). The additional $d_s$ components carry the
    global state into every edge computation. If $s = (T)$ is just
    temperature, $d_s = 1$ and the model sees the same temperature on
    every edge — equivalent to passing $T$ as an extra input feature
    to a CGCNN that has been retrained per temperature.

    The trick that pays off is that the *same* trained MEGNet works for
    all temperatures: training data at $T = 100, 300, 500$ K is pooled,
    and at test time we set $s = T_\text{test}$ and predict. CGCNN
    would require three separate models or temperature as a node
    feature (with the same effective consequence but more awkward
    bookkeeping).

!!! example "Numerical example: temperature-dependent heat capacity"
    Suppose we train MEGNet on Materials Project heat-capacity data at
    100, 200, 300 K and want to predict at 250 K. After training, set
    $s^{(0)} = 250$ (after any normalisation the model used during
    training) and run the forward pass. Because every layer's update
    uses $s$, the prediction smoothly interpolates between the
    network's behaviour at the training temperatures. With CGCNN we
    would either need separate networks or to retrain with $T$ as a
    seventh node feature — equivalent but less elegant.

MEGNet's headline result was a single network trained on
$\sim 64\,000$ Materials Project crystals reaching $0.028$ eV/atom MAE
on formation energies, comfortably beating CGCNN's $\sim 0.04$ eV/atom
on the same benchmark. More striking was that MEGNet matched or beat
specialised SchNet variants on QM9 molecules using the same architecture
with no modifications, suggesting that the inductive bias was right.

When to use MEGNet rather than CGCNN: when your target depends on a
continuous external variable (temperature, pressure, doping), or when
you want a single model to predict several properties at once and you
need a place to inject the property-identity tag. For one-property,
zero-temperature regression, CGCNN gets you most of the way there.

## 10.4.2 ALIGNN: bond angles via the line graph

A persistent weakness of all the architectures so far — CGCNN, SchNet,
MEGNet — is that they see only *distances* between atoms. A pair of
crystals with identical sets of interatomic distances but different
*bond angles* would be assigned identical embeddings. In practice such
ambiguous pairs are rare, but the missing angular information leaves
accuracy on the table: bond angles encode hybridisation and local
coordination geometry, both highly informative.

The Atomistic Line Graph Neural Network (ALIGNN) of Choudhary and
DeCost (2021) fixes this with an elegant construction. They observe
that a graph's *edges* can themselves be promoted to nodes of a new
graph, the *line graph*. In the line graph, each node corresponds to a
bond in the original; two line-graph nodes are connected by a
line-graph edge if the two original bonds share an atom. The line-graph
edge then naturally carries the *bond angle* between the two bonds.

Formally, given $G = (V, E)$, the line graph $L(G) = (E, E')$ has node
set $E$ and edge set $E' = \{(e_1, e_2) : e_1, e_2 \in E, e_1 \cap e_2 \neq \emptyset\}$.
Each line-graph edge $(e_1, e_2)$ — where $e_1 = (u, v)$ and
$e_2 = (v, w)$ share atom $v$ — carries the angle $\theta = \angle uvw$
expanded in a Gaussian basis.

!!! example "Worked construction: line graph of a triangle"
    Take $G$ as a triangle with nodes $\{1, 2, 3\}$ and edges
    $e_{12}, e_{23}, e_{13}$ (treated as undirected). The line graph
    $L(G)$ has three nodes (one per original edge). Edges between
    them: $e_{12}$ and $e_{23}$ share atom 2; $e_{23}$ and $e_{13}$
    share atom 3; $e_{12}$ and $e_{13}$ share atom 1. So $L(G)$ is
    itself a triangle — every pair of original edges shares an atom.

    Each $L(G)$-edge carries a bond angle. Suppose the original
    triangle has angles $60°, 60°, 60°$ at atoms 1, 2, 3 respectively.
    Then the line-graph edges between (12,13), (12,23), (13,23) carry
    those three angles, expanded in a basis.

    ALIGNN then runs MPNN on $G$ to update edge features, then MPNN on
    $L(G)$ (with edge features as nodes) to incorporate angle
    information, alternating. After several alternations, the central
    atom's embedding has absorbed information about both its bond
    lengths *and* its bond angles — the missing chemistry of
    coordination geometry.

!!! note "Why angles distinguish polymorphs"
    Anatase and rutile TiO$_2$ share the same composition and similar
    Ti–O bond lengths (around 1.94 Å in both). They differ in the
    *bond angles*: rutile's TiO$_6$ octahedra are nearly regular
    (angles near 90°), while anatase's are distorted (angles around
    78° and 102°). A distance-only GNN like CGCNN cannot distinguish
    them on geometry alone — only on the slightly different
    distributions of bond lengths and the topology of the edge graph,
    which is a much weaker signal. ALIGNN, by including angles
    explicitly, makes the distinction sharp. This is why ALIGNN's
    advantage over CGCNN is largest on polymorph-rich datasets.

ALIGNN then runs message passing *alternately* on the original crystal
graph and on its line graph: each block updates the bond
representations using both adjacent atoms and adjacent angles, and the
atom representations using updated bonds. After several alternating
blocks, the atom embeddings have absorbed angular information through
the intermediate bond representations.

The accuracy gain is real and consistent. On the Materials Project
formation-energy benchmark, ALIGNN reaches $0.022$ eV/atom MAE versus
MEGNet's $0.028$ and CGCNN's $0.039$. On bulk modulus, shear modulus
and many other elastic properties the gap is larger — angles matter
more for mechanical stiffness than for energetics. As of writing,
ALIGNN remains the strongest published architecture on the full
Matbench suite of property-regression tasks.

The cost is computational. The line graph has up to $|E| \times \bar{k}$
edges (where $\bar{k}$ is the average node degree), which can be ten
times larger than the original edge count. ALIGNN training takes
roughly four times longer than CGCNN training for the same number of
epochs, and the memory footprint is higher. For high-throughput
screening this matters; for one-off model training it does not.

When to use ALIGNN: when angular information is plausibly important
(elastic properties, polymorph energetics, anything involving
hybridisation changes) and you have the compute. For applications
where inference speed matters — embedding millions of candidates — a
faster model like CGCNN with the trade-off explicitly accepted may be
the right call.

## 10.4.3 M3GNet: three-body terms and a universal MLIP

Chen and Ong's 2022 M3GNet paper makes a different bet. The authors set
out not to write the best property-regression GNN — ALIGNN already
existed — but to write the first credible *universal* machine-learning
interatomic potential, capable of predicting energies, forces and
stresses for any element in the periodic table.

To do this they took the line-graph idea, refined it as the inclusion
of explicit three-body terms in the message-passing scheme, and trained
on the full Materials Project *relaxation trajectory* set — every
intermediate geometry produced during the structural optimisation of
every entry in the database.

!!! note "What 'explicit three-body terms' actually means"
    A standard two-body MPNN computes messages that depend on
    $(h_v, h_u, r_{uv})$ — pairwise. A three-body term depends on
    *three* atoms simultaneously: $(h_v, h_u, h_w, r_{uv}, r_{vw}, \theta_{uvw})$.
    The simplest way to inject this into an MPNN is to compute, for
    each pair of edges $(uv, vw)$ sharing atom $v$, a three-body
    message
    $$
    m_{uw|v}^{(3b)} = \phi_3\!\left( h_v, h_u, h_w, e_{uv}, e_{vw}, e_{uvw}^{\text{angle}} \right),
    $$
    and aggregate over all such triples into the central atom.
    Formally this is a tensor-field-network-style construction: the
    three-body message transforms as a scalar under rotation if
    $\phi_3$ uses only invariants ($r_{uv}, r_{vw}, \theta_{uvw}$), or
    as a higher-rank tensor if $\phi_3$ involves spherical harmonics
    of the unit vectors $\hat{\mathbf{r}}_{uv}$ and
    $\hat{\mathbf{r}}_{vw}$ (Chapter 9, §9.2). The training set is about 187 000
relaxation trajectories with $\sim 1.6$ million single-point
calculations. The model is trained jointly on energy, force and stress
losses (Chapter 9, §9.4).

The result is a single model with about $\sim 0.7$ million parameters
that runs on a laptop at hundreds of structures per second and reaches
formation-energy MAE of $0.035$ eV/atom (below the DFT error itself for
many materials), force MAE of about $70$ meV/Å, and supports
geometry relaxations of arbitrary crystals out of the box.

M3GNet's significance is conceptual rather than architectural. It is the
first piece of evidence that one can have a *foundation MLIP* — a
single network trained once, fine-tuned for specific systems, and
reused as a drop-in replacement for DFT in geometry optimisations and
short MD runs. Chapter 12 will revisit this when we discuss CHGNet,
MACE-MP-0 and the broader foundation-model story.

When to use M3GNet (or its successors): whenever you need *forces*, not
just energies, and you are willing to accept errors of around 100 meV/Å
in those forces. M3GNet is the go-to for high-throughput relaxation —
take a database of CIF files, run M3GNet relaxation on each, and you
have approximately-DFT-converged geometries in seconds per structure.

When *not* to use it: when you need DFT-level accuracy on a specific
material class. M3GNet's universal training distributes its capacity
across the periodic table, and a specialised model fine-tuned on, say,
iron-based superconductors will beat M3GNet on iron-based
superconductors. The relationship is exactly the foundation-model-and-
fine-tune pattern of modern machine learning.

## 10.4.4 An honest comparison

| Property | CGCNN | MEGNet | ALIGNN | M3GNet |
| --- | --- | --- | --- | --- |
| Formation energy (eV/atom) | 0.039 | 0.028 | **0.022** | 0.035 |
| Band gap (eV) | 0.388 | 0.330 | **0.218** | 0.330 |
| Bulk modulus ($\log_{10}$ GPa) | 0.071 | 0.060 | **0.051** | 0.068 |
| State-conditioned predictions | no | **yes** | no | no |
| Forces / stresses available | no | no | no | **yes** |
| Inference speed | **fast** | fast | slow | medium |
| Training speed | **fast** | fast | slow | medium |
| Code complexity | **low** | medium | high | high |
| Pre-trained universal weights | no | partial | partial | **yes** |

Numbers above are from the published papers on the Materials Project
test split as of late 2024; details and exact protocols vary, and the
relative ordering is more meaningful than the absolute values.

!!! note "Matbench formation-energy MAE — exact comparison"
    The Matbench `matbench_mp_e_form` task (132 000 stable inorganic
    crystals, 5-fold cross-validation, fixed splits, formation energy
    per atom) gives a cleaner comparison than the Materials Project
    test split. Approximate published MAE values (eV/atom) at the time
    of writing:

    | Architecture | MAE | Year |
    | --- | --- | --- |
    | CGCNN | 0.039 | 2018 |
    | MEGNet | 0.028 | 2019 |
    | SchNet | 0.034 | 2018 |
    | ALIGNN | 0.022 | 2021 |
    | M3GNet | 0.025 | 2022 |
    | MACE-MP-0 | 0.020 | 2023 |

    The trend is clear: each generation has dropped MAE by 5–25%, and
    we are now within a factor of two of the DFT noise floor itself
    (different functionals disagree by 30–50 meV/atom on the same
    structure). Further architectural gains are likely small; the
    leverage has shifted to *data* — larger, cleaner, more diverse
    training sets — and to *fine-tuning* foundation models on
    specialised subsets, the subject of Chapter 12.

The reader who has followed §10.3 to a working CGCNN now has the right
foundation to read the MEGNet, ALIGNN and M3GNet papers. Each is a
specific instantiation of the abstract MPNN framework of §10.2; each
trades complexity for accuracy or generality. The right choice depends
on the question being asked.

## 10.4.5 What the field is currently doing

Beyond the four architectures above, three threads dominate recent
work. The first is *equivariance*: NequIP and MACE (covered in Chapter
9) extend M3GNet's universal-MLIP ambition with explicit
$\mathrm{SO}(3)$-equivariant features, and the resulting models
(MACE-MP-0, MACE-OFF) are competitive with M3GNet at lower parameter
counts. The second is *scale*: the Open Catalyst Project's OC20 and
OC22 datasets, with twenty million single-point calculations, have
become a benchmark for how large materials GNNs can be trained.
GemNet-OC and equivariant successors push these scales further. The
third is *foundation models for graphs*: pre-trained universal
backbones (CHGNet, MACE-MP-0, ORB, MatterSim) that fine-tune on small
target datasets. Chapter 12 takes up this story in detail.

The mood of the field has shifted accordingly. The question is no
longer "what architecture should I train from scratch on my $10^3$
crystals?" — it is "which pre-trained foundation model should I
fine-tune?" CGCNN, MEGNet and ALIGNN remain valuable both as
pedagogical examples and as efficient task-specific models when you
*do* train from scratch. But the production answer in 2026 is
increasingly to start from a foundation model and adapt it. Section
10.5 considers what that adaptation looks like in practice on the
Materials Project, and what subtleties — particularly around
train/test splitting — make the difference between honest accuracy
estimates and self-flattering ones.

## 10.4.5a A decision tree for choosing an architecture

The reader confronting a new property-regression task will rationally
ask: which of these four should I train? A pragmatic decision tree.

*Do I need forces or stresses?* If yes (geometry optimisation,
molecular dynamics, phonons), use M3GNet or MACE. If no, continue.

*Do I have a state variable to condition on (temperature, doping
fraction)?* If yes, use MEGNet. If no, continue.

*Are angles likely to matter? (polymorph energetics, elastic
properties, anything involving hybridisation changes)?* If yes, use
ALIGNN. If no, use CGCNN.

*Do I have $< 10^4$ training labels?* In all cases, *prefer fine-tuning
a foundation model* (CHGNet, MACE-MP-0; Chapter 12) over training from
scratch. The pre-trained representations carry most of the chemistry
your training labels would otherwise have to teach.

This decision tree is not exhaustive — in practice you should also
consider training time, hardware constraints, and the availability of
reference implementations — but it covers 80% of real choices. For the
remaining 20%, ablation studies on your specific data are the only
honest answer.

## 10.4.5b What changes if I want forces?

Returning briefly to the forces question. A CGCNN trained on formation
energies produces $\hat{E}(\text{structure})$. To get forces one would
naively autodifferentiate: $\hat{\mathbf{F}}_i = -\nabla_{\mathbf{r}_i}
\hat{E}$. The technical problem is that CGCNN's edge features depend on
distances $r_{uv}$, not on coordinates directly, and the graph itself
(neighbour assignments) depends on coordinates via the cutoff. The
latter dependence is *discontinuous*: as an atom crosses the cutoff
sphere of a neighbour, the graph topology changes abruptly, and the
forces blow up.

The fix used by M3GNet and successors is to *smoothly attenuate*
messages near the cutoff using a polynomial envelope
$f_{\text{env}}(r) = (1 - r/r_c)^p$ for $r < r_c$, zero beyond.
Messages then vanish smoothly as atoms approach the cutoff, and the
gradient of $\hat E$ with respect to coordinates is well-defined and
continuous. This is one of several technical reasons that a
force-predicting GNN is meaningfully harder to build than an
energy-predicting one. Chapter 9 develops the details.

!!! tip "Cross-reference"
    Chapter 9's MLIP discussion covers the smoothness-of-cutoff
    question, equivariance, and force-loss design. The MPNN
    architectures here are the *structural* basis for those MLIPs; the
    *training* differences are covered there.

## 10.4.6 A note on benchmarking

Every number in this section depends on the test split, the data
preprocessing and the hyperparameter budget. Comparisons in the
literature are not always apples-to-apples. The Matbench benchmark
(Dunn et al., 2020) is the closest the field has to a standardised
playing field for property-regression GNNs; it specifies fixed
train/test splits and evaluation protocols across thirteen tasks and
publishes a public leaderboard. If you are choosing between
architectures for a real project, look at Matbench rankings filtered to
the property class you care about — not the headline number from a
single paper.

We will return to the benchmarking question in §10.5, where the same
pitfall reappears in a different guise: random splits of a polymorph-
rich database systematically overstate model accuracy. The fix is the
*structurally disjoint* split, and the difference can be a factor of
two or three in reported MAE.

!!! note "A reproducible benchmarking workflow"
    To compare two architectures honestly on your own data:

    1. Use exactly the same train/val/test split for both.
    2. Use the same featurisation (cutoff, basis size, etc.).
    3. Allocate the same wall-clock training budget — do not let a
       slow architecture run for longer at the expense of a faster
       one.
    4. Report a mean and standard deviation over $\geq 3$ random
       seeds. A single-seed comparison routinely produces 10% spread
       between runs of the same architecture.
    5. Specify the exact MP database version (or commit hash of your
       data file). The same query can return different structures
       weeks apart.

    This is more bookkeeping than most papers do. Doing it
    nevertheless will reveal that on small datasets the four
    architectures are often within statistical noise of each other,
    and on large datasets the gaps are real but smaller than the
    advertised headline numbers.

### Section summary

- MEGNet adds a global state vector that participates in every
  message-passing update; the right choice for state-conditioned
  prediction (T, P, doping).
- ALIGNN promotes bonds to nodes of a line graph carrying bond
  angles, gaining accuracy on polymorphs and elastic properties at
  $\sim 4\times$ training cost.
- M3GNet adds three-body terms and trains on relaxation trajectories,
  yielding the first credible universal MLIP for forces and stresses.
- Choose the architecture by the physics of the target — angles,
  forces, state-dependence — not by headline benchmark numbers.
