# 9 — Exercises

Eight exercises, mixing pencil-and-paper derivations with PyTorch and
`mace-torch`. Difficulty markers: (E) easy, (M) medium, (H) hard.
Solutions follow each exercise.

---

## Exercise 9.1 (E) — Force from a Behler–Parrinello energy

A BPNN computes the total energy as $U = \sum_i E_i(\mathbf{G}_i)$,
where $\mathbf{G}_i$ is the symmetry-function descriptor of atom $i$.
Show that the force on atom $k$ is

$$
\mathbf{F}_k = - \sum_i \frac{\partial E_i}{\partial \mathbf{G}_i}
               \cdot \frac{\partial \mathbf{G}_i}{\partial \mathbf{r}_k},
$$

and explain why the sum runs over *all* atoms $i$ rather than only
$i = k$.

### Solution

By definition $\mathbf{F}_k = -\nabla_{\mathbf{r}_k} U$. Apply the chain
rule to each term in the sum:

$$
\nabla_{\mathbf{r}_k} E_i(\mathbf{G}_i)
  = \frac{\partial E_i}{\partial \mathbf{G}_i}
    \cdot \nabla_{\mathbf{r}_k} \mathbf{G}_i.
$$

The descriptor $\mathbf{G}_i$ depends on $\mathbf{r}_k$ whenever atom
$k$ is a neighbour of atom $i$ — i.e. whenever $r_{ik} < r_\mathrm{c}$.
Summing over $i$ gives the stated expression. The reason the sum is
not restricted to $i = k$ is that moving atom $k$ changes not only
$E_k$ (because $\mathbf{G}_k$ depends on $k$'s neighbours) but also the
energies of every atom for which $k$ is a neighbour. This is the
familiar action–reaction structure of pairwise interactions.

---

## Exercise 9.2 (M) — Implement and test a Behler $G^2$ descriptor

Using the code in §9.3.1 as a starting point, implement
`g2_descriptor` and verify the following properties on a randomly
generated 8-atom periodic cell:

(a) **Translation invariance.** Translate every atom by a random
vector; the descriptor should be unchanged to floating-point tolerance.

(b) **Rotation invariance.** Rotate every atom and the lattice by a
random $\mathrm{SO}(3)$ matrix; the descriptor should be unchanged.

(c) **Permutation invariance.** Reorder the atom indices; the
descriptor of the (new) atom $i$ should equal the descriptor of the
(old) atom that maps to position $i$.

(d) **Force consistency via finite differences.** With the descriptor
$\mathbf{G}_i$, write a toy energy $E_i = \sum_k G_{ik}$ and compute
the force on atom $0$ both analytically (by differentiating the
descriptor) and by finite differences. Confirm agreement to at least
six significant figures.

### Solution

```python
import numpy as np
from scipy.spatial.transform import Rotation

# Set up a random 8-atom cubic cell.
rng = np.random.default_rng(42)
cell = np.eye(3) * 8.0
positions = rng.uniform(0, 8, size=(8, 3))
eta = np.array([1.0, 2.0, 4.0])
r_s = np.array([1.5, 2.5, 3.5])

G0 = g2_descriptor(positions, cell, eta, r_s)

# (a) Translation
t = rng.uniform(-5, 5, size=3)
G_t = g2_descriptor((positions + t) % 8.0, cell, eta, r_s)
assert np.allclose(G0, G_t, atol=1e-10)

# (b) Rotation (must rotate both positions and the cell)
R = Rotation.random(random_state=rng).as_matrix()
G_r = g2_descriptor(positions @ R.T, cell @ R.T, eta, r_s)
assert np.allclose(G0, G_r, atol=1e-8)

# (c) Permutation
perm = rng.permutation(8)
G_p = g2_descriptor(positions[perm], cell, eta, r_s)
assert np.allclose(G0[perm], G_p, atol=1e-10)
```

For force consistency, the simplest route is `torch.autograd` on a
PyTorch reimplementation; finite differences with step
$h = 10^{-5}\,\text{\AA}$ should agree with autograd to roughly
$10^{-7}$. The full implementation is left to the reader; the key
check is

```python
E_plus  = float(g2_descriptor(positions + h*e0, cell, eta, r_s).sum())
E_minus = float(g2_descriptor(positions - h*e0, cell, eta, r_s).sum())
F_num   = -(E_plus - E_minus) / (2*h)
assert np.allclose(F_num, F_analytic, rtol=1e-5)
```

with `e0` the unit vector along $x$ applied only to atom $0$.

---

## Exercise 9.3 (M) — Prove SOAP power-spectrum rotation invariance

Show that the SOAP power spectrum

$$
p_{nn'\ell} = \frac{1}{\sqrt{2\ell+1}}
              \sum_m c_{n\ell m} c_{n'\ell m}^*
$$

is invariant under rotations, using only the rotation rule
$c_{n\ell m} \mapsto \sum_{m'} D^{(\ell)}_{m m'}(R) c_{n\ell m'}$ and
the unitarity of the Wigner D-matrix.

### Solution

Under rotation,
$c_{n\ell m} \mapsto \tilde c_{n\ell m} = \sum_{m_1} D^{(\ell)}_{m m_1} c_{n\ell m_1}$
and likewise $c_{n'\ell m}^* \mapsto \sum_{m_2} D^{(\ell)\,*}_{m m_2} c_{n'\ell m_2}^*$.
Substitute and exchange the order of summation:

$$
\tilde p_{nn'\ell}
  = \frac{1}{\sqrt{2\ell+1}}
    \sum_{m, m_1, m_2}
    D^{(\ell)}_{m m_1}\, D^{(\ell)\,*}_{m m_2}\,
    c_{n\ell m_1}\, c_{n'\ell m_2}^*.
$$

The Wigner D-matrix is unitary:
$\sum_m D^{(\ell)}_{m m_1} D^{(\ell)\,*}_{m m_2} = \delta_{m_1 m_2}$.
Applying this collapses the sum over $m$,

$$
\tilde p_{nn'\ell}
  = \frac{1}{\sqrt{2\ell+1}}
    \sum_{m_1} c_{n\ell m_1}\, c_{n'\ell m_1}^*
  = p_{nn'\ell}.\qquad\square
$$

---

## Exercise 9.4 (M) — Cutoff continuity

The Behler cosine cutoff is
$f_\mathrm{c}(r) = \tfrac12[\cos(\pi r / r_\mathrm{c}) + 1]$ for
$r < r_\mathrm{c}$ and $0$ otherwise. Show that $f_\mathrm{c}$ is
$C^1$ but not $C^2$ at $r = r_\mathrm{c}$. What goes wrong in a
molecular dynamics simulation if the cutoff is only $C^0$?

### Solution

At $r = r_\mathrm{c}$, $f_\mathrm{c}(r_\mathrm{c}) = \tfrac12(\cos\pi + 1) = 0$;
$f_\mathrm{c}'(r) = -\tfrac{\pi}{2 r_\mathrm{c}} \sin(\pi r/r_\mathrm{c})$,
so $f_\mathrm{c}'(r_\mathrm{c}) = 0$. The function and its first
derivative both vanish at the boundary, so $f_\mathrm{c}$ is $C^1$.
The second derivative is
$f_\mathrm{c}''(r) = -\tfrac{\pi^2}{2 r_\mathrm{c}^2} \cos(\pi r/r_\mathrm{c})$,
so $f_\mathrm{c}''(r_\mathrm{c}) = \tfrac{\pi^2}{2 r_\mathrm{c}^2} \neq 0$
on the inside but $0$ on the outside: a jump discontinuity. Hence not
$C^2$.

If the cutoff were only $C^0$ (function continuous but derivative
jumps at $r_\mathrm{c}$), the force — which is a derivative of the
energy — would be discontinuous whenever a neighbour crossed the
cutoff. In NVE molecular dynamics, the energy receives a delta-function
impulse each crossing, and the total energy drifts at a rate
proportional to the crossing rate (typically a few percent per
picosecond for a $C^0$ cutoff). In thermostatted ensembles the drift
is masked but shows up as spurious heating from the thermostat.

---

## Exercise 9.5 (E) — Atomic energy decomposition

Argue that the decomposition $U = \sum_i E_i$ with strictly local
$E_i$ (depending only on neighbours within $r_\mathrm{c}$ of atom
$i$) cannot exactly reproduce long-range Coulomb interactions. What
strategies exist to recover long-range physics within an MLIP
framework?

### Solution

A strictly local decomposition truncates *all* contributions to $E_i$
at $r_\mathrm{c}$. Coulomb interactions decay only as $1/r$, so the
sum over neighbours outside $r_\mathrm{c}$ is conditionally convergent
in three dimensions and contributes finite tails that cannot be
captured by local sums. Equivalently, a charged defect has an energy
contribution from screening that extends over many atomic spacings.

Standard strategies to recover long-range physics:

- **Charge-equilibration models.** Predict per-atom partial charges
  with a small auxiliary network, then add an Ewald sum analytically.
  Examples: PhysNet, BpopNN, the charge-equilibration variant of MACE.
- **Long-range message passing.** Use a second network with a larger
  cutoff to capture slowly decaying interactions, at higher
  computational cost.
- **Hybrid embedding.** Run DFT or a polarisable model on a small QM
  region embedded in an MLIP-treated bulk.

For most condensed-phase materials problems, a $5$ to $7\,\text{\AA}$
cutoff captures enough of the electrostatic screening that explicit
long-range terms are unnecessary; for ionic crystals and electrochemistry
they become essential.

---

## Exercise 9.6 (H) — Tensor-product dimensions

Using the Clebsch–Gordan decomposition
$\ell_1 \otimes \ell_2 = (\ell_1 + \ell_2) \oplus \cdots \oplus |\ell_1 - \ell_2|$,
compute:

(a) The decomposition of $1 \otimes 1$, $1 \otimes 2$, and
$2 \otimes 2$.

(b) The total dimension of each tensor product, and verify it equals
$(2\ell_1 + 1)(2\ell_2 + 1)$.

(c) In a NequIP layer with $\ell_\mathrm{max} = 2$, what tensor-product
pairs $(\ell_1, \ell_2)$ contribute to the output $\ell = 1$ channel?

### Solution

(a) Apply the formula:

- $1 \otimes 1 = 2 \oplus 1 \oplus 0$. Dimensions: $5 + 3 + 1 = 9$.
- $1 \otimes 2 = 3 \oplus 2 \oplus 1$. Dimensions: $7 + 5 + 3 = 15$.
- $2 \otimes 2 = 4 \oplus 3 \oplus 2 \oplus 1 \oplus 0$. Dimensions:
  $9 + 7 + 5 + 3 + 1 = 25$.

(b) Check $(2 \cdot 1 + 1)(2 \cdot 1 + 1) = 9$, $3 \cdot 5 = 15$,
$5 \cdot 5 = 25$. All match.

(c) The output $\ell = 1$ appears in $\ell_1 \otimes \ell_2$ whenever
$|\ell_1 - \ell_2| \le 1 \le \ell_1 + \ell_2$. With
$\ell_\mathrm{max} = 2$ this gives the pairs
$(0, 1), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2)$. The
NequIP layer will combine features from each of these channels
(weighted by their respective radial MLP outputs) when constructing the
new $\ell = 1$ messages.

---

## Exercise 9.7 (H) — Train a tiny MACE potential

Take 50 configurations from the water dataset of §9.6 (or download a
small benchmark such as the rMD17 aspirin set). Train a small MACE
model (hidden_irreps `"32x0e + 32x1o"`, two layers, correlation $3$)
for $50$ epochs. Compute energy and force MAE on a held-out 10-frame
test set. Compare to (i) the rule-of-thumb of $\sim 25\,\mathrm{meV}/\text{\AA}$
force MAE on $1000$ training configurations and (ii) the value you
might predict by scaling as $1/\sqrt{N_\mathrm{train}}$.

### Solution

A 50-frame MACE training run on water typically yields:

| Set size | Force MAE | Energy MAE |
|---|---|---|
| 50 | 120 meV/Å | 4 meV/atom |
| 1000 | 25 meV/Å | 0.5 meV/atom |

The naive $1/\sqrt{N}$ scaling predicts
$25 \times \sqrt{1000/50} \approx 110\,\mathrm{meV}/\text{\AA}$ — close to
the observed $120$. The match is partly coincidental: at very small
training sizes the model is underdetermined and the error is dominated
by model bias rather than by data noise, so the actual scaling can
flatten or steepen. Nonetheless, $1/\sqrt{N}$ is a useful
order-of-magnitude prior when planning training budgets.

If you cannot get below $\sim 200\,\mathrm{meV}/\text{\AA}$ on the 50-frame
training set, the most likely cause is overfitting (training loss
much lower than validation loss). Try reducing hidden_irreps to
`"16x0e + 16x1o"` and raising weight decay to $10^{-5}$.

---

## Exercise 9.8 (M) — Diagnose a bad parity plot

A colleague hands you a trained MACE potential for a Cu-Zn alloy.
The energy parity plot looks tight, but the force parity plot has
*two distinct diagonal clouds*, one for each element, offset from
the true diagonal in opposite directions: predicted Cu forces are
systematically $5\%$ too small, predicted Zn forces $5\%$ too large.

(a) What is the most likely cause?

(b) How would you verify your diagnosis?

(c) Suggest two fixes, ranked by effort.

### Solution

(a) The most likely cause is *mismatched isolated-atom energies
($E_0$ values)*. If $E_0(\mathrm{Cu})$ is too high and
$E_0(\mathrm{Zn})$ is too low, the model compensates by tilting the
per-atom-energy contribution: it learns to under-credit Cu's
contribution to the cohesive energy and over-credit Zn's. The
gradients of these biased atomic energies — i.e. the forces — inherit
the tilt. Alternatively, the issue could be a bug in unit conversion
specific to one element (e.g. a misread of the DFT output file that
swapped Cu and Zn columns), but the symmetric $\pm 5\%$ pattern is
the fingerprint of $E_0$ imbalance.

(b) Verify by (i) recomputing the isolated-atom energies with the
same DFT settings used for the bulk training data, and comparing to
the values used during training; (ii) checking that the *mean* energy
prediction error is zero when both Cu and Zn are present in equal
proportion, but nonzero in pure-Cu or pure-Zn configurations; (iii)
inspecting `model.state_dict()['atomic_energies']` to confirm the
stored values match what the user intended.

(c) Two fixes:

1. **Cheap.** Recompute $E_0$ values correctly and retrain (or
   fine-tune the existing checkpoint for $5$–$10$ epochs with the
   corrected $E_0$). The model's per-atom-energy bias will absorb the
   correction quickly. Total cost: a few hours.

2. **Robust.** Add a per-element scalar offset as a *learnable*
   parameter of the model, with weak regularisation, so that the
   model can absorb small $E_0$ errors automatically. This is what
   `mace-torch` does when `atomic_energies` is set to `'average'`
   instead of an explicit table; the offsets are fitted by least
   squares to the training set. Total cost: refactor and retrain.
   The drawback is that the absolute energies are no longer
   physically interpretable.

For an immediate paper deadline, fix (1) is the right answer; for a
production code intended to ingest community datasets with varying
DFT conventions, fix (2) is more durable.

---

## Closing comment

You now have all the components needed to read the current MLIP
literature and to fit potentials for your own systems. Chapter 10
generalises the message-passing structure introduced here to *graph
neural networks for materials properties*: instead of forces, the
target is band gap, formation energy, elastic tensor — but the
representation-theoretic machinery is the same. Chapter 11 shows how
to close the active-learning loop, and Chapter 12 explores foundation
models such as MACE-MP-0 that amortise the training cost across all
of materials chemistry.
