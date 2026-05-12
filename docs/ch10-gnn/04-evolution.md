# 10.4 Evolution: MEGNet, ALIGNN, M3GNet

CGCNN appeared in 2018 and demonstrated the basic point: graph neural
networks reach DFT-level accuracy on bulk property prediction at orders
of magnitude lower inference cost. Five years of refinements then
followed, each motivated by a specific architectural shortcoming. This
section traces the line through three landmark architectures —
MEGNet (2019), ALIGNN (2021) and M3GNet (2022) — and is candid about
what each is and is not good for.

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
every entry in the database. The training set is about 187 000
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
