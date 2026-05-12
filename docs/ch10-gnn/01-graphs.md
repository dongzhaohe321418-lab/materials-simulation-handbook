# 10.1 Crystals as Graphs

The central abstraction of this chapter is small enough to write in a
single line and rich enough to fill a textbook: a crystal is a graph. To
make that statement operational we need to settle three questions. What
is a graph, formally? How do we attach the physical content of a crystal
— elements, distances, angles — to that graph? And how do we cope with
the fact that a crystal is periodic, so any cutoff-based neighbourhood
must respect the unit-cell translations?

## 10.1.1 The mathematical object

!!! note "Why this step?"
    Before we attach any physics, we need a precise vocabulary for the
    combinatorial object underneath. Sloppy definitions here propagate
    silently into bugs four hundred lines deeper in the dataloader.

A graph $G = (V, E)$ is an ordered pair consisting of a finite set $V$ of
*nodes* (also called vertices) and a set $E \subseteq V \times V$ of
*edges*. An edge is an unordered pair $\{u, v\}$ of distinct nodes, or,
in a *directed* graph, an ordered pair $(u, v)$. In a *multigraph* we
permit several distinct edges between the same pair of nodes; in a graph
with *self-loops* a node may be connected to itself.

Stated more formally still: in a simple undirected graph $E$ is a subset
of the unordered pairs $\binom{V}{2}$ and $|E| \leq \binom{|V|}{2}$. In a
directed graph $E \subseteq V \times V$ and $|E| \leq |V|^2$. In a
directed *multigraph*, $E$ is a multiset rather than a set, so the same
pair $(u, v)$ may appear several times — each occurrence a separate
edge with potentially distinct features. This last case is exactly the
one we need for crystals, for reasons that will become clear once
periodicity enters in §10.1.5.

!!! example "Intuition: graph as a labelled diagram"
    Picture five dots on a sheet of paper, with arrows between some of
    them. Each dot carries a small ticket of numbers (a feature vector);
    each arrow carries its own small ticket. That picture, formalised, is
    a directed graph with node features and edge features. The crystal
    graphs of this chapter are precisely such pictures, drawn in three
    dimensions and tiled periodically.

We will need two derived quantities throughout the chapter. The
*degree* of a node $v$ in a directed graph is the number of edges with
$v$ as their source, $\deg_+(v) = |\{u : (v, u) \in E\}|$, plus the
number with $v$ as their target, $\deg_-(v) = |\{u : (u, v) \in E\}|$.
For an undirected graph the two are equal. The *neighbourhood*
$\mathcal{N}(v) = \{u : (u, v) \in E\}$ collects the source nodes of all
edges incoming to $v$ — exactly the set we will sum over when computing
messages in §10.2.

For our purposes the natural object is a *directed multigraph with
self-loops*. The directionality matters because the message from atom
$u$ to atom $v$ (it sees a neighbour at displacement $+\mathbf{r}_{uv}$)
is generally different from the message from $v$ to $u$ (which sees a
displacement $-\mathbf{r}_{uv}$). The multiplicity matters because, under
periodic boundary conditions, a single atom may be connected to several
periodic images of another atom, each at a different displacement. The
self-loops matter because an atom may be connected to its own periodic
images when the unit cell is small enough.

!!! tip "Cross-reference"
    Compare with Chapter 5, where we wrote crystals as Bravais lattices
    plus a basis; the graph representation is *complementary*. It throws
    away the global lattice information after computing displacements,
    keeping only local neighbourhood structure. The two descriptions are
    formally equivalent for finite cutoffs but emphasise different
    invariants — lattice symmetry versus local chemistry.

Every node carries a feature vector $h_v \in \mathbb{R}^{d_V}$ and every
edge carries a feature vector $e_{uv} \in \mathbb{R}^{d_E}$. These
features are the only information that subsequent neural-network layers
will see. The graph topology and its features together constitute the
*input representation*.

## 10.1.2 Node features: encoding elements

A node represents an atom. The minimum information it must carry is the
chemical element. There are two routes.

**One-hot encoding.** Reserve one component per element in the periodic
table — usually 100, occasionally 118 — and set the entry corresponding
to the actual element to one, all others to zero. The advantage is
extreme simplicity; the disadvantage is that the network has no prior
knowledge that carbon resembles silicon more than sodium. It must learn
the entire periodic system from data.

**Hand-crafted descriptors.** CGCNN, in its original incarnation,
concatenates nine properties: group number, period, electronegativity,
covalent radius, valence electrons, first ionisation energy, electron
affinity, block (s/p/d/f) and atomic volume. These are categorised into
bins, each bin one-hot encoded, and the resulting 92-dimensional vector
forms the initial $h_v^{(0)}$. The handcrafted prior accelerates training
on small datasets.

**Learned embeddings.** Most modern architectures simply allocate an
embedding table — a learnable matrix $E \in \mathbb{R}^{N_\text{elem}
\times d}$ — and look up the embedding for each atom by atomic number. On
large datasets these learned embeddings cluster the periodic table into
recognisable patterns (transition metals together, halogens together)
without any hand-coded chemistry, which is reassuring.

For the rest of this chapter we use a learned embedding with $d = 64$
unless otherwise stated.

!!! note "Tradeoffs in featurisation"
    The choice between one-hot, hand-crafted and learned embeddings is
    not purely aesthetic. Three axes matter.

    *Sparse versus dense.* One-hot is sparse (one nonzero entry per
    node); learned embeddings are dense. Dense inputs train faster on
    GPUs because the first linear layer becomes a matrix multiplication
    over the whole batch rather than an indexed look-up plus a near-zero
    matmul.

    *Prior knowledge versus data hunger.* Hand-crafted descriptors
    inject chemistry into the model from line one. They pay off when
    the training set is small (a few thousand crystals): the network
    does not waste capacity learning that O has six valence electrons.
    On the other hand, a hand-crafted descriptor risks encoding the
    practitioner's bias, and on large datasets ($10^6$ crystals) a
    learned embedding generally surpasses any hand-crafted set.

    *Atomic number versus full embedding.* Treating an atom as just its
    atomic number $Z$ is the most extreme reduction: it deliberately
    discards isotope, oxidation state, and spin information. For
    formation-energy regression this is fine — the property does not
    depend on isotope and oxidation state is implicit in the
    coordination. For magnetic property prediction one would extend the
    node feature with spin or magnetic moment; for nuclear properties
    (NMR shifts, hyperfine couplings) one would add isotopic mass and
    nuclear spin.

!!! example "What CGCNN's nine-feature vector actually looks like"
    For an oxygen atom, the nine features (each binned into roughly ten
    one-hot categories) encode: group = 16, period = 2, electronegativity
    $\approx 3.44$, covalent radius $\approx 66$ pm, valence electrons
    = 6, first ionisation energy $\approx 13.6$ eV, electron affinity
    $\approx 1.46$ eV, block = p, atomic volume $\approx 14$ cm$^3$/mol.
    The 92-dimensional one-hot expansion of these nine values is the
    initial node feature. For silicon, the same nine slots take
    different values, and the encoded vector overlaps with oxygen's by
    perhaps two of the ninety-two bins (mostly through proximity of
    covalent radius and atomic volume). The model uses this overlap as
    an inductive bias: similar Magpie patterns should produce similar
    contributions to the total property.

## 10.1.3 Edge features: encoding distances

An edge represents a neighbour relationship. Its central physical content
is the interatomic distance $r_{uv} = \|\mathbf{r}_v - \mathbf{r}_u\|$.

A single scalar is too brittle a feature for a neural network. The
network must compute, for example, "is this distance close to 2.4 Å or
to 2.6 Å?" — and a raw scalar input requires the early layers to learn a
sharp threshold function before any chemistry can be extracted. The fix
is universal: expand the distance in a fixed *radial basis*. The
standard choice is a Gaussian basis,
$$
\phi_k(r) = \exp\!\left[ -\frac{(r - \mu_k)^2}{2 \sigma^2} \right],
\qquad k = 1, \ldots, K,
$$
with means $\mu_k$ uniformly spaced between $r_\text{min} = 0$ and
$r_\text{max} = 8$ Å in steps of roughly $\Delta = 0.2$ Å, and width
$\sigma = \Delta$. The result is a $K$-dimensional vector that varies
smoothly with $r$ and whose components have natural interpretations as
"closeness to $\mu_k$".

More sophisticated bases exist — Bessel functions with cutoff envelopes
yield a discrete basis with the right asymptotic behaviour — but the
Gaussian basis is the workhorse and we will use it throughout the CGCNN
implementation.

!!! note "Why expand the distance at all?"
    A single scalar input to a ReLU network can only express piecewise
    linear functions of $r$. To approximate a chemically sensible
    distance dependence — high response near a covalent bond length,
    rapid decay beyond — the network needs to combine many ReLU pieces.
    Expanding into a basis of $K$ smooth functions $\phi_k(r)$ converts
    the problem to a linear combination over $K$ inputs and lets the
    network express the right shape after a *single* linear layer. The
    Gaussian basis is the radial analogue of the Fourier basis on the
    interval $[0, r_\text{max}]$: it provides smooth, localised features
    that interact cleanly with subsequent dense layers.

!!! example "Numerical check: smoothness of the expansion"
    With $K = 41$ centres uniformly spaced on $[0, 8]$ Å and
    $\sigma = 0.2$ Å, each $\phi_k(r)$ has effective support of roughly
    $\pm 0.4$ Å. At any given distance the expansion is dominated by
    three or four nearby basis functions. Two distances $r$ and
    $r + \Delta r$ with $\Delta r \ll \sigma$ produce expansion vectors
    that differ by no more than $\|\partial \phi / \partial r\| \cdot
    \Delta r$ in each component, ensuring continuity. The corresponding
    crystal property — formation energy, say — is therefore predicted
    as a smooth function of $r$, with no surprise jumps at basis
    boundaries.

Beyond distances, additional edge features encode richer chemistry.
Line graphs (§10.4.2) promote each edge to a node carrying a *bond angle*
between two edges sharing an atom, $\theta_{uvw} = \angle(uvw)$, also
expanded in a Gaussian (or cosine) basis on $[0, \pi]$. Some
architectures attach *bond-length variance* over the neighbourhood of
$v$ as a node-level feature, capturing local strain. For most
property-regression tasks the distance expansion alone suffices.

Edges may additionally carry the displacement vector
$\mathbf{r}_{uv}$ (not merely its norm), bond-order indicators, or
chemical labels from a bonding analysis. We will see in §10.1.5 below
that retaining the full displacement is necessary if downstream layers
are to be *equivariant* (Chapter 9, §9.2); CGCNN itself uses only the
scalar distance and is therefore strictly invariant.

## 10.1.4 Constructing the neighbour list

Given a `pymatgen.Structure` or `ase.Atoms` object, how do we enumerate
the edges? Two strategies dominate.

**Fixed cutoff.** Pick a radius $r_\text{cut}$ — usually 5 to 8 Å — and
connect every pair of atoms whose distance is less than $r_\text{cut}$.
Simple, deterministic, and well suited to neural networks because the
neighbour count is bounded and roughly homogeneous across structures.

!!! note "Derivation: average degree as a function of cutoff"
    In a crystal with atomic number density $\rho$ (atoms per Å$^3$), the
    expected number of neighbours within a sphere of radius $r_\text{cut}$
    around any given atom equals $\rho$ multiplied by the volume of that
    sphere, minus one for the central atom itself:
    $$
    \langle \deg \rangle \approx \frac{4}{3}\pi r_\text{cut}^3 \rho - 1.
    $$
    A typical oxide has $\rho \approx 0.07$ Å$^{-3}$ (around 5 atoms per
    70 Å$^3$ unit cell). At $r_\text{cut} = 5$ Å this gives
    $\langle\deg\rangle \approx 36$; at $r_\text{cut} = 8$ Å it gives
    $\langle\deg\rangle \approx 150$. The graph cost — and the memory
    footprint of the edge tensor — therefore scales as the *cube* of the
    cutoff, which is why doubling the cutoff is a much bigger decision
    than it appears. For dense intermetallics
    ($\rho \approx 0.09$ Å$^{-3}$) the degree at 5 Å is closer to 47;
    for porous metal-organic frameworks ($\rho \approx 0.02$ Å$^{-3}$)
    it can be as low as 10. Inspect $\langle\deg\rangle$ before settling
    on a cutoff: it sets both the GPU memory budget and, indirectly, the
    receptive-field size of the GNN (§10.2.5).

!!! example "Worked number"
    For rutile TiO$_2$ ($a = b = 4.59$ Å, $c = 2.96$ Å, $Z = 2$
    formula units per cell, six atoms per cell), the cell volume is
    $V \approx 62.3$ Å$^3$, giving $\rho = 6/62.3 \approx 0.096$ Å$^{-3}$.
    A 5 Å cutoff therefore yields about $\frac{4}{3}\pi \cdot 125 \cdot
    0.096 - 1 \approx 49$ neighbours per atom, in good agreement with
    the empirical 40–55 found by enumerating with `pymatgen`'s
    `get_neighbor_list` (the discrepancy is because the formula
    overestimates very slightly: the sphere extends just beyond the
    unit-cell boundary).

**Coordination geometry.** Use `pymatgen.analysis.local_env.CrystalNN`
or the older `VoronoiNN` to identify only chemically meaningful
neighbours. Each algorithm computes a weighted Voronoi tessellation and
discards weak contacts. The resulting graph is sparser and arguably more
physical, but the cutoff is data-dependent and varies from atom to atom.

The two strategies need not be exclusive. CGCNN originally used a
$k$-nearest-neighbour list with $k = 12$ and a maximum cutoff of 8 Å,
which is essentially a hybrid.

Here is a minimal example using a fixed cutoff with `pymatgen`:

```python
from __future__ import annotations

import numpy as np
from pymatgen.core import Structure


def build_graph_cutoff(
    structure: Structure,
    cutoff: float = 5.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build edge list and edge features for a periodic structure.

    Returns:
        Z:          (N,) atomic numbers.
        edge_index: (2, E) source/target atom indices.
        distances:  (E,) edge distances in Angstrom.
        offsets:    (E, 3) integer image vectors n such that
                    r_uv = r_v + n . lattice - r_u.
    """
    Z = np.array(structure.atomic_numbers, dtype=np.int64)
    src: list[int] = []
    dst: list[int] = []
    dists: list[float] = []
    offsets: list[np.ndarray] = []

    # get_neighbor_list returns one entry per directed edge already.
    centres, points, images, distances = structure.get_neighbor_list(
        r=cutoff, exclude_self=True
    )
    for c, p, image, d in zip(centres, points, images, distances):
        src.append(int(c))
        dst.append(int(p))
        dists.append(float(d))
        offsets.append(np.asarray(image, dtype=np.int64))

    edge_index = np.stack([np.array(src), np.array(dst)], axis=0)
    return Z, edge_index, np.array(dists), np.stack(offsets)
```

And here is the alternative using `CrystalNN`, which prefers chemically
sensible neighbours but is roughly two orders of magnitude slower:

```python
from pymatgen.analysis.local_env import CrystalNN


def build_graph_crystalnn(
    structure: Structure,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build a graph using CrystalNN's chemically informed neighbours."""
    Z = np.array(structure.atomic_numbers, dtype=np.int64)
    nn = CrystalNN()

    src: list[int] = []
    dst: list[int] = []
    dists: list[float] = []
    offsets: list[np.ndarray] = []
    for i in range(len(structure)):
        for entry in nn.get_nn_info(structure, i):
            j = entry["site_index"]
            image = np.array(entry["image"], dtype=np.int64)
            site = structure[i]
            neighbour = entry["site"]
            d = float(site.distance(neighbour))
            src.append(i)
            dst.append(j)
            dists.append(d)
            offsets.append(image)

    edge_index = np.stack([np.array(src), np.array(dst)], axis=0)
    return Z, edge_index, np.array(dists), np.stack(offsets)
```

Use the first for high-throughput training and the second when the
chemistry of the bonding is the object of study.

## 10.1.5 Periodicity and image vectors

The single subtlest point in graph construction for crystals is
periodicity. A periodic crystal is not really a finite collection of $N$
atoms; it is an infinite collection generated by translating the unit
cell by all integer combinations of the lattice vectors. When we draw
edges with a cutoff we must include edges that *leave the unit cell* and
arrive at periodic images of atoms.

!!! note "Why bother with image offsets at all?"
    A naïve graph constructor would say: atom $j$ is a neighbour of atom
    $i$ if $\|\mathbf{r}_j - \mathbf{r}_i\| < r_\text{cut}$. That works
    for an isolated molecule but is catastrophically wrong for a
    crystal. In NaCl with $a = 5.64$ Å, the nearest Na–Cl distance is
    $a/2 = 2.82$ Å — but the *home cell* Na and Cl coordinates differ
    by $(0.5, 0.5, 0.5)$ in fractional units, giving a Cartesian
    distance of $\sqrt{3}/2 \cdot a \approx 4.88$ Å. The actual nearest
    neighbour sits in an *adjacent* image. Without image offsets, the
    code would miss the Na–Cl bond entirely.

Concretely: let $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$ be the lattice
vectors and let $\mathbf{r}_i$ denote the fractional position of atom
$i$. The displacement from atom $u$ in the home cell to atom $v$ in the
image cell labelled by integer vector $\mathbf{n} = (n_1, n_2, n_3)$ is
$$
\mathbf{r}_{uv}^{(\mathbf{n})}
= \mathbf{r}_v - \mathbf{r}_u + n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + n_3 \mathbf{a}_3.
$$
Each combination $(v, \mathbf{n})$ within $r_\text{cut}$ becomes its own
edge in the graph. So a single bond in real space — say the Si–O
contact in $\alpha$-quartz — typically corresponds to several edges in
the graph, one per equivalent image. For small cells, an atom may even
be connected to its own periodic images, giving a self-loop.

The image vector $\mathbf{n}$ is therefore an additional piece of edge
data. CGCNN, since it operates only on scalar distances, can discard
$\mathbf{n}$ after computing $r_{uv}$. Equivariant networks like NequIP,
MACE and M3GNet must retain it because they pass the displacement
$\mathbf{r}_{uv}^{(\mathbf{n})}$ into spherical-harmonic features.

A common bug: practitioners building their own graph constructors
forget the image offsets and use the minimum-image convention naïvely,
which silently truncates the neighbour list when the cutoff exceeds half
the shortest lattice vector. The `pymatgen` and `ase` neighbour-list
utilities handle this correctly; rolling your own is rarely worth the
risk.

!!! note "Derivation: how many image cells must we search?"
    Given a cutoff $r_\text{cut}$ and lattice vectors
    $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$, we need only consider
    image cells $\mathbf{n}$ whose centre lies within $r_\text{cut} +
    d_\text{max}$ of the home cell, where $d_\text{max}$ is the maximum
    intra-cell atom separation. A safe upper bound on the number of
    images to enumerate along axis $i$ is
    $n_i^{\max} = \lceil r_\text{cut} / h_i \rceil$, where
    $h_i$ is the perpendicular distance from the origin to the plane
    spanned by the other two lattice vectors. For a cubic cell with
    $a = 4$ Å and $r_\text{cut} = 8$ Å, $n_i^{\max} = 2$, so we sweep
    $\mathbf{n} \in \{-2, \ldots, 2\}^3$ — 125 image cells. For a
    plate-like cell with $c = 1$ Å in the third direction, we sweep
    $\mathbf{n}_3 \in \{-8, \ldots, 8\}$, and the cost rises sharply.
    Production code uses smarter spatial-hashing schemes that skip
    irrelevant image cells, but the *upper bound* above is what guards
    against forgotten images.

!!! tip "Forward reference"
    The displacement vector $\mathbf{r}_{uv}^{(\mathbf{n})}$ is the key
    object that connects this chapter to Chapter 9: equivariant networks
    pass it through spherical-harmonic projectors $Y_\ell^m(\hat{\mathbf{r}})$
    to obtain features that transform as $\mathrm{SO}(3)$ irreps. CGCNN,
    being invariant, throws away the direction and keeps only the norm.

## 10.1.6 Bipartite and heterogeneous graphs

So far every node is an atom and every edge a neighbour contact — a
*homogeneous* graph. Materials science offers natural generalisations.

A *bipartite* graph has two disjoint node sets and edges only between
them. A surface adsorption study might give one set to substrate atoms
and another to adsorbate atoms, with edges only across the interface.
The model can then maintain two separate embedding tables and two sets
of message-passing weights, capturing the asymmetry of the problem.

A *heterogeneous* graph allows multiple node and edge *types*. A grain
boundary might contain bulk-like atoms, interface atoms and defect
atoms; a metal–organic framework might distinguish framework metals,
organic linker atoms and guest molecules. Each type has its own
embedding, and message-passing operations are defined per edge type.
The PyTorch Geometric library has first-class support for these via
`HeteroData` objects.

These extensions are conceptually straightforward — the message-passing
framework of §10.2 handles them with only minor notational changes — but
they require care with data. Most materials databases provide only the
plain structure, and the user is responsible for assigning types
manually or via a clustering algorithm. We will not use heterogeneous
graphs in the CGCNN implementation that follows, but the reader should
be aware that they exist and that recent work on defects and interfaces
relies heavily on them.

!!! note "Key Idea (Box 10.1.A)"
    A crystal becomes a *directed multigraph with self-loops* in three
    steps: (i) atoms become nodes carrying element-based features,
    (ii) pairs within cutoff become edges carrying Gaussian-expanded
    distances, (iii) periodic boundary conditions are respected by
    recording an integer image-vector $\mathbf{n}$ per edge. The
    resulting object is the universal input to every GNN in this
    chapter.

## 10.1.7 Where we are

We now have a precise pipeline. Take a `Structure`. Build a directed
multigraph with image offsets respecting periodic boundary conditions.
Attach learned element embeddings to nodes, Gaussian-expanded distances
to edges. The output is a tuple `(h_V, e_E, edge_index, offsets)` ready
for a neural network.

Section 10.2 introduces the neural network. We will see that nearly all
modern GNNs are special cases of a single abstraction — message passing
— and that designing a new architecture amounts to choosing two
functions, $M_t$ and $U_t$.

### Section summary

- A crystal graph is a *directed multigraph with self-loops*; the
  multiplicity and self-loops are required by periodic boundary
  conditions.
- Node features encode chemistry (one-hot / Magpie / learned
  embedding); edge features encode distances expanded in a Gaussian
  basis.
- Average node degree scales as $\frac{4}{3}\pi r_c^3 \rho - 1$, so
  cutoff choice has a *cubic* effect on graph cost.
- Image vectors $\mathbf{n}$ record which periodic image of a
  neighbour an edge connects to; CGCNN can discard them after
  computing distances, equivariant networks cannot.
