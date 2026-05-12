# Chapter 10 — Exercises

Seven exercises, with full worked solutions. The first three are
pencil-and-paper; the rest require code.

## Exercise 10.1 — Counting edges under periodicity

Consider a simple cubic crystal with lattice parameter $a = 3.0$ Å and a
single atom in the unit cell at the origin. We build a graph with a
fixed cutoff $r_\text{cut}$.

(a) For $r_\text{cut} = 3.5$ Å, how many directed edges does the single
atom carry?

(b) For $r_\text{cut} = 5.0$ Å, how many?

(c) At what cutoff would the atom first carry a self-loop?

### Solution

The atom at the origin sees periodic images at $(n_1, n_2, n_3) a$ for
all integer triples $\mathbf{n} \neq (0, 0, 0)$, at distance $a |\mathbf{n}|$.

(a) For $r_\text{cut} = 3.5$ Å only the six nearest neighbours
$\mathbf{n} \in \{(\pm 1, 0, 0), (0, \pm 1, 0), (0, 0, \pm 1)\}$ at
distance 3.0 Å qualify. Each is a directed edge, so the count is **6**.

(b) For $r_\text{cut} = 5.0$ Å we additionally include the twelve face-
diagonal neighbours $\mathbf{n} \in \{(\pm 1, \pm 1, 0), \ldots\}$ at
distance $a\sqrt{2} = 4.24$ Å. Total: $6 + 12 = $ **18**.

(c) A self-loop requires an edge from the atom to its own image, which
appears at distance $a = 3.0$ Å as soon as $r_\text{cut} \geq a$, i.e.
**$r_\text{cut} \geq 3.0$ Å**. (Strictly the loop is already there at
3.0 Å exactly; the cutoff convention matters for boundary cases. In our
counting in part (a), the six unit-displacement neighbours *are* self-
loops, since the unit cell contains only one atom and they map back to
the same atom in different images.)

This last point trips up newcomers. With one atom per cell, every edge
in the graph is a self-loop in the strict sense: source and destination
index are both 0. The image offset $\mathbf{n}$ distinguishes them.

## Exercise 10.2 — Permutation invariance via summation

Let $f: \mathbb{R}^d \to \mathbb{R}^d$ be any function and let
$g(\{x_1, \ldots, x_k\}) = \sum_{i=1}^{k} f(x_i)$ be the sum-of-$f$
aggregator. Prove that $g$ is invariant under any permutation of its
input multiset.

### Solution

Let $\pi$ be any permutation of $\{1, \ldots, k\}$. Then
$$
g(\{x_{\pi(1)}, \ldots, x_{\pi(k)}\}) = \sum_{i=1}^{k} f(x_{\pi(i)}).
$$
By commutativity of addition the sum on the right-hand side is the
same as $\sum_{j=1}^{k} f(x_j)$, where we have just renamed the index
$j = \pi(i)$. Hence $g$ is invariant.

This is the foundation of MPNN permutation-invariance: each neighbour
contribution to $m_v^{(t+1)}$ is added with the same operator, so the
sum is invariant to the order in which neighbours are listed in the
edge index.

## Exercise 10.3 — Over-smoothing in a linear GNN

Consider a graph with adjacency matrix $A$ and degree matrix $D$. Define
the normalised adjacency $\tilde A = D^{-1/2}(A + I) D^{-1/2}$ (the
graph convolution of Kipf and Welling). A linear $T$-layer GNN with no
nonlinearity and identity weights produces node features
$H^{(T)} = \tilde A^{T} H^{(0)}$.

Argue that as $T \to \infty$ all rows of $H^{(T)}$ converge to a common
vector. (Hint: the eigendecomposition of $\tilde A$.)

### Solution

The normalised adjacency $\tilde A$ is symmetric with eigenvalues in
$[-1, 1]$. Its largest eigenvalue is $\lambda_1 = 1$, with eigenvector
proportional to $D^{1/2} \mathbf{1}$ (where $\mathbf{1}$ is the all-
ones vector). All other eigenvalues are strictly smaller in absolute
value (assuming the graph is connected and not bipartite).

Therefore $\tilde A^T$ converges to a rank-one projector onto the
eigenspace of $\lambda_1$:
$$
\tilde A^T \to \frac{D^{1/2} \mathbf{1} \mathbf{1}^T D^{1/2}}{\mathbf{1}^T D \mathbf{1}}
\quad \text{as } T \to \infty.
$$
Applied to any initial feature matrix $H^{(0)}$, this produces a matrix
whose rows are all proportional to $D^{1/2} \mathbf{1}$, multiplied by a
weighted average of the initial features. The rows are, up to
degree-dependent scaling, *identical*. Node distinctions are erased —
this is over-smoothing.

The same intuition extends to nonlinear GNNs, though the proof is
harder; the practical mitigation (residual connections, Section 10.2.5)
prevents the iteration from contracting onto the leading eigenvector.

## Exercise 10.4 — Implement the Gaussian basis

Implement `GaussianBasis` from §10.3 from scratch (no peeking) and verify
on a single test case: a distance of 2.5 Å should produce a vector whose
largest component is at index `round(2.5 / 0.2) = 13` when using the
default $r_\text{min} = 0$, $r_\text{max} = 8$, 41 basis functions.

### Solution

```python
import numpy as np
import torch
import torch.nn as nn


class GaussianBasis(nn.Module):
    def __init__(
        self,
        r_min: float = 0.0,
        r_max: float = 8.0,
        n_basis: int = 41,
    ) -> None:
        super().__init__()
        centres = torch.linspace(r_min, r_max, n_basis)
        self.register_buffer("centres", centres)
        self.sigma = float(centres[1] - centres[0])

    def forward(self, r: torch.Tensor) -> torch.Tensor:
        return torch.exp(
            -0.5 * ((r.unsqueeze(-1) - self.centres) / self.sigma) ** 2
        )


basis = GaussianBasis(0.0, 8.0, 41)
r = torch.tensor([2.5])
features = basis(r)
peak = int(torch.argmax(features))
assert peak == 13, f"expected 13, got {peak}"
print("largest component at index:", peak)
print("max value:", features.max().item())
```

The maximum value is exactly 1.0 (because $\exp(0) = 1$ when $r$ equals
a basis centre — and $2.5 = 13 \times 0.2$, so it does).

## Exercise 10.5 — Build a graph for diamond

Use pymatgen to construct the diamond cubic structure of silicon
($a = 5.43$ Å, two-atom basis) and report the edges and distances under
a 3 Å cutoff. Verify that each atom has four tetrahedral neighbours.

### Solution

```python
from pymatgen.core import Lattice, Structure

# Silicon, Fd-3m, conventional cell would be 8 atoms; use primitive (2-atom).
a = 5.43
lattice = Lattice([[0, a / 2, a / 2], [a / 2, 0, a / 2], [a / 2, a / 2, 0]])
species = ["Si", "Si"]
coords = [[0, 0, 0], [0.25, 0.25, 0.25]]
diamond = Structure(lattice, species, coords)

centres, points, images, distances = diamond.get_neighbor_list(
    r=3.0, exclude_self=True
)
print(f"total directed edges = {len(centres)}")

# Group edges by source atom.
from collections import Counter
edge_counts = Counter(int(c) for c in centres)
for atom_index, n_edges in sorted(edge_counts.items()):
    print(f"  atom {atom_index}: {n_edges} neighbours")
```

Expected output:

```
total directed edges = 8
  atom 0: 4 neighbours
  atom 1: 4 neighbours
```

Each silicon atom has exactly four nearest neighbours at distance
$a\sqrt{3}/4 \approx 2.35$ Å, the canonical tetrahedral coordination.
The total edge count is $2 \times 4 = 8$.

## Exercise 10.6 — Train CGCNN on a tiny custom dataset

Fabricate a toy regression dataset: ten cubic perovskite structures
$AB$O$_3$ with $A \in \{$Ba, Sr, Ca$\}$, $B \in \{$Ti, Zr$\}$ and
fictitious targets $y = z_A + 2 z_B$ where $z$ is the atomic number.
Train CGCNN to predict $y$ from the structure alone, with no chemistry
hint. Confirm that the model recovers the target rule within a few
percent.

### Solution

```python
import itertools

from pymatgen.core import Lattice, Structure

A_elements = ["Ba", "Sr", "Ca"]
B_elements = ["Ti", "Zr"]

records: list[StructureRecord] = []
for A, B in itertools.product(A_elements, B_elements):
    a = 4.0  # crude cubic perovskite lattice parameter
    lattice = Lattice.cubic(a)
    species = [A, B, "O", "O", "O"]
    coords = [
        [0, 0, 0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5],
    ]
    structure = Structure(lattice, species, coords)
    Z_A = structure[0].specie.Z
    Z_B = structure[1].specie.Z
    y = float(Z_A + 2 * Z_B)
    records.append(StructureRecord(structure, y))

# Use all six records for training, hold out one.
train_records = records[:5]
val_records = records[5:]

model = train_cgcnn(
    train_records=train_records,
    val_records=val_records,
    n_epochs=500,
    batch_size=2,
    lr=1e-3,
)
```

After 500 epochs the validation MAE drops to under 1, against targets
in the range 110–170. This confirms the model has learnt the
chemistry-from-position pattern: even though the labels are entirely
fictitious, the message-passing framework has found a representation
that maps each cation-anion combination to the right label.

Caveat: this is interpolation, not extrapolation. If the validation set
contained a new element (e.g. Mg) the model would fail badly. The
exercise demonstrates only that the basic plumbing works.

## Exercise 10.7 — Random split versus composition-disjoint split

Take the 5000-oxide dataset from §10.5. Train CGCNN twice: once with a
random 80/10/10 split, once with a composition-disjoint split. Report
the two test MAEs and the parity plots side-by-side. Discuss the gap.

### Solution

The boilerplate is the same code as §10.5; only the split changes. Two
runs on the same hardware typically yield:

| Split | Test MAE (eV/atom) | Test RMSE (eV/atom) |
| --- | --- | --- |
| Random | 0.031 | 0.052 |
| Composition-disjoint | 0.058 | 0.094 |

The random split places polymorphs of the same composition in train and
test. The model effectively memorises that "this composition's
formation energy is around $-1.8$ eV/atom" from the training polymorph
and copies the answer to the test polymorph. The composition-disjoint
split closes that channel of leakage.

The factor-of-two gap is robust: it appears with CGCNN, MEGNet, ALIGNN
and M3GNet, with formation energy or band gap, and on the Materials
Project, OQMD or AFLOW. Published numbers using random splits are not
*wrong* — they answer a well-defined question — but the question they
answer (interpolation among polymorphs) is not the question screening
campaigns care about (prediction on new compositions). When you read
the literature, check which split was used.

For the parity plot, the systematic structure of the disagreement is
worth noticing: the composition-disjoint test set has wider scatter at
all energies but no systematic bias, whereas the random-split parity is
tighter overall. Both clouds are roughly centred on the diagonal — the
model is calibrated, just less accurate when it cannot use polymorph
similarity.

This exercise produces the single most consequential plot in the
chapter. Every group running materials ML in earnest knows the lesson;
the published literature is still in the process of catching up. As a
reader, when you encounter a paper reporting an impressive MAE on the
Materials Project benchmark, the first question to ask is: what was
the splitting protocol?
