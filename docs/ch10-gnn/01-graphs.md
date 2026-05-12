# 10.1 Crystals as Graphs

The central abstraction of this chapter is small enough to write in a
single line and rich enough to fill a textbook: a crystal is a graph. To
make that statement operational we need to settle three questions. What
is a graph, formally? How do we attach the physical content of a crystal
— elements, distances, angles — to that graph? And how do we cope with
the fact that a crystal is periodic, so any cutoff-based neighbourhood
must respect the unit-cell translations?

## 10.1.1 The mathematical object

A graph $G = (V, E)$ is an ordered pair consisting of a finite set $V$ of
*nodes* (also called vertices) and a set $E$ of *edges*. An edge is an
unordered pair $\{u, v\}$ of distinct nodes, or, in a *directed* graph,
an ordered pair $(u, v)$. In a *multigraph* we permit several distinct
edges between the same pair of nodes; in a graph with *self-loops* a node
may be connected to itself.

For our purposes the natural object is a *directed multigraph with
self-loops*. The directionality matters because the message from atom
$u$ to atom $v$ (it sees a neighbour at displacement $+\mathbf{r}_{uv}$)
is generally different from the message from $v$ to $u$ (which sees a
displacement $-\mathbf{r}_{uv}$). The multiplicity matters because, under
periodic boundary conditions, a single atom may be connected to several
periodic images of another atom, each at a different displacement. The
self-loops matter because an atom may be connected to its own periodic
images when the unit cell is small enough.

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
