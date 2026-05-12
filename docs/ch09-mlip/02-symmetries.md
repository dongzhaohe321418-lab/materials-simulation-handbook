# 9.2 Required symmetries

Before we choose any regression model we must decide what function we
are trying to fit. The potential energy of a configuration of atoms is
not an arbitrary function of $3N$ Cartesian coordinates: it inherits a
collection of symmetries from the physics of the problem. Every
successful MLIP architecture enforces these symmetries *exactly*, in
the structure of its representation, rather than hoping the regressor
will learn them from data. This section catalogues the five symmetries
that matter — translation, rotation, permutation, smoothness,
compactness — and introduces the distinction between *invariant* and
*equivariant* features that organises the modern literature.

## 9.2.0 Invariance and equivariance — the central distinction

Before walking through the individual symmetries, we make sharp the
single most important distinction in the modern MLIP literature: that
between *invariance* and *equivariance*. Many beginners conflate the
two and then find Section 9.5 mysterious. The conceptual investment
here pays back many times over.

!!! note "Definitions"
    Let $G$ be a group (e.g. $\mathrm{O}(3)$, rotations and
    reflections in 3D) acting on the input space of a function $f$.
    Let $\rho_\mathrm{in}$ denote that action on the input and
    $\rho_\mathrm{out}$ its action on the output. The function $f$ is

    - **invariant** under $G$ if
      $f(\rho_\mathrm{in}(g) \cdot \mathbf{x}) = f(\mathbf{x})$ for
      every group element $g \in G$. The output never changes.
    - **equivariant** under $G$ if
      $f(\rho_\mathrm{in}(g) \cdot \mathbf{x}) = \rho_\mathrm{out}(g) \cdot
      f(\mathbf{x})$ for every $g \in G$. The output transforms in a
      *predictable, group-consistent way*.

    Invariance is the special case of equivariance in which
    $\rho_\mathrm{out}$ is the *trivial* representation, i.e. the
    output is a scalar.

The energy $U(\{\mathbf{r}_i\})$ is invariant under translations,
rotations, reflections, and same-species permutations. Forces and
stresses are *not* invariant: they transform predictably under those
operations. Specifically, if we apply a rotation $R$ to every position
in the system, then

$$
U(\{R\mathbf{r}_i\}) = U(\{\mathbf{r}_i\}), \qquad
\mathbf{F}_i(\{R\mathbf{r}_j\}) = R\,\mathbf{F}_i(\{\mathbf{r}_j\}),
$$

and the stress tensor transforms as
$\sigma'_{\alpha\beta} = R_{\alpha\gamma} R_{\beta\delta} \sigma_{\gamma\delta}$.
Energy is scalar, forces are vectors ($\ell = 1$), stress is a
symmetric rank-2 tensor ($\ell = 0$ plus $\ell = 2$).

!!! example "A concrete check"
    Take three atoms at $\mathbf{r}_1 = (0, 0, 0)$,
    $\mathbf{r}_2 = (1, 0, 0)$, $\mathbf{r}_3 = (0, 1, 0)$ and any
    rotation, say a $90^\circ$ rotation $R$ about the $z$-axis,
    $R: (x, y, z) \mapsto (-y, x, z)$. The rotated configuration is
    $R\mathbf{r}_1 = (0, 0, 0)$, $R\mathbf{r}_2 = (0, 1, 0)$,
    $R\mathbf{r}_3 = (-1, 0, 0)$.

    The pairwise distances are unchanged: $|r_{12}| = |R r_{12}| = 1$,
    $|r_{13}| = |R r_{13}| = 1$, $|r_{23}| = |R r_{23}| = \sqrt{2}$.
    The energy, being a function of these distances, is unchanged
    (invariance). If atom $1$ originally felt force $\mathbf{F}_1 =
    (F_x, F_y, F_z)$, in the rotated frame it feels $R\mathbf{F}_1 =
    (-F_y, F_x, F_z)$ — the force has rotated alongside the
    configuration (equivariance).

The whole point of this chapter is that one can choose to build an
MLIP whose *internal features* are either invariant or equivariant.
Both yield invariant energies (just take a final inner product to
collapse equivariant features to scalars). But equivariant features
preserve more information per dimension and therefore yield more
data-efficient potentials — the empirical lesson of §9.5.

The slogan to memorise: **invariance throws away information,
equivariance keeps it.** Both are correct; one is wasteful.

## 9.2.0a Why differentiability matters

The forces that drive molecular dynamics are not measured; they are
*computed* as gradients of the energy with respect to positions:

$$
\mathbf{F}_i = -\nabla_{\mathbf{r}_i} U(\{\mathbf{r}_j\}).
$$

This identity, banal in classical mechanics, has sharp consequences
for MLIP design.

**Consequence 1: continuity.** The energy must be a *continuous*
function of every atomic coordinate. A discontinuity at any
configuration produces a delta-function in the force, which integrates
to a finite-momentum kick in MD. The kick injects energy on every
neighbour-list update, and energy conservation in NVE fails
catastrophically.

**Consequence 2: differentiability.** The energy must be
*differentiable*. A continuous but kinked function (think
$|r - r_0|$ as a stand-in bond term) has finite forces on either
side but no well-defined force at the kink, and the integrator
oscillates pathologically near it.

**Consequence 3: smooth derivatives for thermodynamic integration
and free-energy methods.** Many post-MD analyses — free-energy
perturbation, thermodynamic integration, response-function evaluation —
require *second* derivatives of the energy (the Hessian). Models with
$C^1$ but not $C^2$ smoothness will produce noisy phonon spectra and
unreliable enhanced-sampling free-energy estimates.

The implication for design: we will require not just *invariance and
equivariance* of the energy, but *smooth invariance and equivariance*.
This rules out hard-cutoff descriptors, requires smooth cutoff
envelopes (§9.2.4), and forbids $\mathrm{ReLU}$-like non-smooth
activations in the network (we use $\tanh$, SiLU, softplus instead).

A useful heuristic: the energy should be at least $C^2$ for honest
phonon work and $C^1$ at the very least for stable MD. Best-of-breed
MACE potentials are $C^\infty$ inside the cutoff sphere and $C^p$ at
the cutoff with $p \approx 5$.

## 9.2.0b On the body-order hierarchy

The atomic energy $E_i$ depends on the geometry of all neighbours of
$i$ inside the cutoff. With $n$ neighbours one can in principle
construct features of *body order* $1, 2, 3, \dots, n + 1$ (the $+1$
counts atom $i$ itself). A two-body feature depends on $i$ and one
neighbour; a three-body feature depends on $i$ and two neighbours
(and thus involves an angle); a four-body feature involves three
neighbours; and so on.

Why truncate? Two reasons:

1. **Combinatorial explosion.** With $n = 12$ neighbours (a typical
   first coordination shell), there are $n$ two-body terms, $\binom{n}{2}
   = 66$ three-body terms, $\binom{n}{3} = 220$ four-body terms,
   $\binom{n}{4} = 495$ five-body terms. The number of distinct
   feature parameters grows correspondingly.
2. **Diminishing physical relevance.** Most chemistries are well
   described by terms up to four-body order: covalent bonds
   (two-body), bond angles (three-body), dihedrals (four-body). The
   five-body and higher contributions to the energy are small
   corrections.

The practical consequence is the body-order truncation: Behler
descriptors capture two- and three-body, ACE and MACE go to four- or
five-body, and that is essentially sufficient. Where it is not — in
some metallic systems with strong electronic-structure non-locality,
in conjugated $\pi$-systems with delocalised electrons — one needs
deeper message-passing networks whose effective receptive field grows
beyond a single cutoff sphere.

We will return to body order constantly in §9.3 (it is the organising
principle of the ACE expansion) and §9.5 (it sets MACE apart from
NequIP).

## 9.2.1 Translation invariance

The first symmetry is the simplest. Translate every atom by the same
vector $\mathbf{t}$ and nothing physical changes:

$$
U(\{\mathbf{r}_i + \mathbf{t}\}) = U(\{\mathbf{r}_i\})
\quad \text{for all } \mathbf{t} \in \mathbb{R}^3.
$$

This rules out functional forms that depend on absolute positions.
Concretely, a network that takes $\mathbf{r}_i$ as input and learns
biases on the position coordinates will memorise the simulation box
and fail the moment you translate. The fix is also simple: depend only
on *relative* coordinates,

$$
\mathbf{r}_{ij} \equiv \mathbf{r}_j - \mathbf{r}_i,
$$

or on functions of relative coordinates such as scalar distances
$r_{ij} = \|\mathbf{r}_{ij}\|$ and angles
$\cos\theta_{ijk} = \hat{\mathbf{r}}_{ij}\cdot\hat{\mathbf{r}}_{ik}$.

!!! note "Why distances are automatically translation- and rotation-invariant"
    The pairwise distance $r_{ij} = \|\mathbf{r}_j - \mathbf{r}_i\|$ is
    invariant under any global translation $\mathbf{r}_i \mapsto
    \mathbf{r}_i + \mathbf{t}$, because the translation cancels in
    the difference. It is also invariant under any rotation $R \in
    \mathrm{O}(3)$, because

    $$
    \|R\mathbf{r}_j - R\mathbf{r}_i\| =
    \|R(\mathbf{r}_j - \mathbf{r}_i)\| =
    \|\mathbf{r}_j - \mathbf{r}_i\|,
    $$

    using the fact that rotations preserve the Euclidean norm. So
    *any* function $f$ built from pairwise distances alone satisfies
    $f(\{R\mathbf{r}_i + \mathbf{t}\}) = f(\{\mathbf{r}_i\})$ for
    every $(R, \mathbf{t})$ — translation and rotation invariance
    come for free.

    The corresponding statement for angles, which are inner products of
    unit vectors $\hat{\mathbf{r}}_{ij} \cdot \hat{\mathbf{r}}_{ik}$, is
    the same: inner products are rotation-invariant scalars. Two- and
    three-body invariants built from distances and angles are
    automatically $\mathrm{O}(3) \times \mathbb{R}^3$-invariant.

    This is why §9.3.1's Behler symmetry functions, built only from
    $r_{ij}$ and $\cos\theta_{ijk}$, are invariant by inspection — no
    extra machinery is needed.

In a periodic cell with lattice vectors $\mathbf{L}_a$ the relative
coordinate must be taken modulo the lattice: we use the *minimum-image*
convention or, more generally, list all neighbours $j$ such that
$\|\mathbf{r}_j + \mathbf{L} - \mathbf{r}_i\| < r_\mathrm{c}$ for any
lattice translation $\mathbf{L}$. Every MLIP architecture in this
chapter operates on a neighbour list of relative vectors, never on
absolute coordinates.

!!! tip "Common mistake"
    Beginners sometimes encode the periodic cell into the model by
    concatenating the lattice matrix into the input. This is unnecessary
    and counter-productive: it breaks the local decomposition of
    energy and gives the network a route to overfit on box size. The
    correct treatment is to build the neighbour list using the periodic
    boundary conditions and feed only relative neighbour vectors to the
    model.

## 9.2.2 Rotation invariance — and rotation equivariance

Rotate the entire configuration about a point: again nothing physical
changes. Letting $R \in \mathrm{SO}(3)$ act on each position,

$$
U(\{R \mathbf{r}_i\}) = U(\{\mathbf{r}_i\})
\quad \text{for all } R \in \mathrm{SO}(3).
$$

The same statement extends to $\mathrm{O}(3)$ if the energy is parity
invariant, as it is for non-chiral systems. Energy is a *scalar*; it
transforms trivially under rotation.

Forces, however, are *vectors*. Under a rotation they pick up the
rotation:

$$
\mathbf{F}_i(\{R\mathbf{r}_j\}) = R\, \mathbf{F}_i(\{\mathbf{r}_j\}).
$$

A function that transforms in this predictable way is called
*equivariant*: rotate the input, the output rotates correspondingly.
Equivariance generalises invariance — an invariant is an equivariant
that happens to transform as the trivial (scalar) representation.

Internally, an MLIP can carry either invariant or equivariant features.
Both choices yield an invariant energy (you can always project an
equivariant feature down to a scalar at the end), but they differ in
how much information they preserve about the local geometry.

Consider two carbon atoms in a tetrahedral environment, one at the
centre and another at the apex of a methylene bridge. The scalar
description of the central atom — interatomic distances and bond angles
to its four neighbours — does not by itself reveal where any one of
those neighbours sits in space; that information is collapsed the
moment you take inner products. An equivariant feature, by contrast,
keeps the direction information explicitly, as a vector indexed by
$\mathrm{O}(3)$ irreducible-representation labels. Later layers can
combine those vectors, taking dot products to form invariants when they
are needed and keeping the directionality when they are not.

The mathematics of equivariance lives in representation theory. The
irreducible representations of $\mathrm{O}(3)$ are labelled by an
integer $\ell = 0, 1, 2, \dots$ and a parity. The $\ell = 0$ irrep is
the scalar, the $\ell = 1$ irrep is the vector (three components
transforming as $\mathbf{v} \mapsto R\mathbf{v}$), the $\ell = 2$ irrep
is the symmetric traceless rank-2 tensor (five components transforming
under the rotation matrices $D^{(2)}(R)$), and so on. A general
equivariant feature is a list of vectors $\mathbf{x}^{(\ell)}$ each of
which transforms as
$\mathbf{x}^{(\ell)} \mapsto D^{(\ell)}(R)\,\mathbf{x}^{(\ell)}$.
Sections 9.3 and 9.5 will make this concrete via the spherical
harmonics $Y_\ell^m$, which are the canonical basis for the $\ell$-th
irrep on the unit sphere.

The empirical lesson, which we will revisit, is that throwing away
direction information at every layer is wasteful. Equivariant networks
need far less data to reach a given accuracy because their inductive
bias matches the symmetry of the underlying physics. We will see in
§9.5 that MACE on rMD17 achieves SchNet-level accuracy with roughly
$1/20$ the training data, a difference attributable almost entirely
to the choice of equivariant features.

!!! tip
    The slogan: *Energy is invariant. Forces are equivariant. Internal
    features can be either, but equivariant features waste less
    information.*

## 9.2.3 Permutation invariance

If atoms $i$ and $j$ are of the same chemical species, swapping their
labels cannot change the energy. Letting $\pi$ be any permutation of
atom indices that respects element identity,

$$
U(\{\mathbf{r}_{\pi(i)}\}) = U(\{\mathbf{r}_i\}).
$$

This rules out architectures that index atoms by their position in a
list. It is the reason MLIPs are built as *atom-centred sums*: the
energy is

$$
U = \sum_i E_i(\text{environment of } i),
$$

where each atomic contribution $E_i$ is a function of the *unordered
set* of neighbours of atom $i$, partitioned by element. A function of
an unordered set can be implemented in two ways. The first is to sum a
function of each element of the set:

$$
\phi(\{x_j\}) = \sum_j f(x_j).
$$

This *deep set* construction is permutation-invariant by inspection.
Most descriptors in §9.3 — Behler symmetry functions, SOAP power
spectra, ACE basis functions — are sums over neighbours and inherit
permutation invariance for free.

The second route is message passing. A graph neural network defines
features on each atom and updates them by aggregating messages from
neighbours,

$$
h_i^{(t+1)} = \mathrm{update}\!\left(h_i^{(t)},\;
  \mathop{\mathrm{aggregate}}_{j \in \mathcal{N}(i)} m(h_i^{(t)}, h_j^{(t)}, \mathbf{r}_{ij})\right),
$$

where the aggregation is sum, mean, or max — all permutation-invariant.
NequIP and MACE are message-passing networks of this kind, with
equivariant features playing the role of $h$.

!!! example "Permutation invariance — three atoms by hand"
    Let us verify permutation invariance for three identical atoms
    explicitly. Take three hydrogens at $\mathbf{r}_1 = (0, 0, 0)$,
    $\mathbf{r}_2 = (1, 0, 0)$, $\mathbf{r}_3 = (0, 1, 0)$, and a
    Behler radial descriptor

    $$
    G_i = \sum_{j \neq i} g(r_{ij}), \qquad g(r) = e^{-r^2}.
    $$

    With the original labelling,
    $G_1 = g(1) + g(1) = 2e^{-1}$;
    $G_2 = g(1) + g(\sqrt{2}) = e^{-1} + e^{-2}$;
    $G_3 = g(1) + g(\sqrt{2}) = e^{-1} + e^{-2}$.

    Now swap labels $1 \leftrightarrow 2$, so the same three atoms are
    relabelled as $\mathbf{r}'_1 = (1, 0, 0)$, $\mathbf{r}'_2 =
    (0, 0, 0)$, $\mathbf{r}'_3 = (0, 1, 0)$.
    $G'_1 = g(1) + g(\sqrt{2}) = e^{-1} + e^{-2}$;
    $G'_2 = g(1) + g(1) = 2e^{-1}$;
    $G'_3 = g(1) + g(\sqrt{2}) = e^{-1} + e^{-2}$.

    The descriptor values are *permuted along with the labels*, but
    the multiset $\{G_1, G_2, G_3\}$ is identical. The sum
    $E = E_1 + E_2 + E_3 = f(G_1) + f(G_2) + f(G_3)$ — for any choice
    of regression function $f$ — is therefore invariant under the
    relabelling. The atom-centred sum structure is what makes
    permutation invariance automatic.

    Now consider what would go wrong with a non-atom-centred model.
    Suppose we tried $E = \mathrm{NN}(\mathbf{r}_1, \mathbf{r}_2,
    \mathbf{r}_3)$ as a single function of nine ordered coordinates.
    The neural network would treat $(0, 0, 0; 1, 0, 0; 0, 1, 0)$ and
    $(1, 0, 0; 0, 0, 0; 0, 1, 0)$ as different inputs, and would
    return different energies unless trained to symmetrise — which it
    will only learn approximately, never exactly. Atom-centred sums
    are the architectural commitment that makes permutation invariance
    exact rather than approximate.

## 9.2.4 Smoothness and the cutoff function

Molecular dynamics integrators rely on forces that are continuous and
differentiable functions of position. A potential whose forces jump or
diverge will break energy conservation in the NVE ensemble and produce
artefacts in any thermostatted ensemble. Smoothness is therefore not a
luxury but a hard requirement.

The challenge arises at the cutoff. To keep the neighbour list finite
we discard atoms beyond a distance $r_\mathrm{c}$. If we simply
truncate the sum at $r_\mathrm{c}$, every time an atom drifts across
the cutoff boundary the energy jumps by a finite amount and the force
acquires a delta-function spike. The remedy is a *smooth cutoff
function* $f_\mathrm{c}(r)$ that decays smoothly to zero at
$r_\mathrm{c}$:

$$
f_\mathrm{c}(r) = \begin{cases}
\tfrac12\!\left[\cos\!\left(\pi r / r_\mathrm{c}\right) + 1\right]
  & r < r_\mathrm{c},\\
0 & r \ge r_\mathrm{c}.
\end{cases}
$$

This is the standard Behler cutoff. It is $C^1$: continuous with
continuous first derivative at $r = r_\mathrm{c}$, which is what
molecular dynamics requires. Higher-order alternatives — polynomial
$(1 - r/r_\mathrm{c})^p$ envelopes with $p \ge 4$, or the
$1/r^p$-style envelopes used in MACE — buy additional smoothness at
the cost of slightly more arithmetic.

!!! note "Why this step? Derivative of the cosine cutoff"
    Let us check that the cosine cutoff vanishes smoothly. Set $u =
    \pi r / r_\mathrm{c}$ for brevity. Then $f_\mathrm{c}(r) =
    \tfrac{1}{2}(\cos u + 1)$ and

    $$
    f_\mathrm{c}'(r) = -\frac{\pi}{2 r_\mathrm{c}} \sin\!\left(\frac{\pi r}{r_\mathrm{c}}\right).
    $$

    At $r = r_\mathrm{c}$: $\cos(\pi) = -1$, so $f_\mathrm{c}(r_\mathrm{c}) =
    \tfrac{1}{2}(-1 + 1) = 0$. Good — the value goes to zero.

    Also at $r = r_\mathrm{c}$: $\sin(\pi) = 0$, so
    $f_\mathrm{c}'(r_\mathrm{c}) = 0$. The first derivative also vanishes
    at the cutoff. This is the $C^1$ property: both $f_\mathrm{c}$ and
    its first derivative are continuous at $r = r_\mathrm{c}$ (with the
    convention $f_\mathrm{c}(r) = 0$ for $r > r_\mathrm{c}$). Forces,
    which involve $\partial r_{ij}/\partial \mathbf{r}_i$ multiplied
    by $f_\mathrm{c}'(r_{ij})$ (chain rule through a descriptor),
    therefore go smoothly to zero at the cutoff and there is no jump.

    Sanity check: $f_\mathrm{c}(0) = \tfrac{1}{2}(\cos 0 + 1) = 1$. The
    cutoff is unity at the origin and decays smoothly to zero at
    $r_\mathrm{c}$ — exactly the envelope we want.

    The second derivative does *not* vanish at $r_\mathrm{c}$:
    $f_\mathrm{c}''(r_\mathrm{c}) = -\pi^2/(2 r_\mathrm{c}^2) \cos(\pi)
    = \pi^2 / (2 r_\mathrm{c}^2) \neq 0$. So the cosine cutoff is $C^1$
    but not $C^2$ at the boundary. For phonon work and for sensitive
    free-energy estimates one prefers a polynomial cutoff
    $(1 - r/r_\mathrm{c})^p$ with $p \ge 4$, which is $C^{p-1}$ at the
    boundary.

Every descriptor in this chapter applies $f_\mathrm{c}$ at every place
where a neighbour enters a sum:

$$
G_i = \sum_{j \in \mathcal{N}(i)} g(r_{ij}) f_\mathrm{c}(r_{ij}),
$$

so that adding or removing a neighbour at the cutoff has vanishing
effect.

!!! tip "Common mistake"
    A symmetric variant — multiplying by $f_\mathrm{c}(r_{ij})$ inside
    the radial part *and* inside the angular part of a three-body
    descriptor — is required, not optional. If you forget to apply
    $f_\mathrm{c}$ to the angular triple $\{ij, ik\}$ you will see
    spurious oscillations in the radial distribution function and the
    NVE drift will become visible after a few picoseconds.

## 9.2.5 Compactness

The descriptor of atom $i$ should be a *fixed-length vector*,
independent of how many neighbours $i$ has. A copper atom in the bulk
has twelve neighbours within a typical cutoff; the same atom on a
surface might have nine; an interstitial copper near a dislocation
might have fifteen. The regression model downstream — a feed-forward
neural network, a Gaussian process — expects inputs of a fixed
dimension and cannot accommodate a variable-length list directly.

The deep-set construction discussed above solves this. By writing the
descriptor as a sum of one-neighbour or two-neighbour contributions,
$\sum_j f(\mathbf{r}_{ij})$ or $\sum_{j,k} g(\mathbf{r}_{ij}, \mathbf{r}_{ik})$,
the dimension of the result is determined by the parameterisation of
$f$ or $g$ and not by the number of neighbours. Bartók's SOAP power
spectrum (§9.3.2), for example, is a vector of dimension
$N_\mathrm{rad}^2 (\ell_\mathrm{max} + 1)/2$, independent of
neighbour count.

Compactness sits alongside locality: the radial cutoff $r_\mathrm{c}$
bounds the number of neighbours physically (by the average atomic
density times $\tfrac43 \pi r_\mathrm{c}^3$), and the descriptor
parameterisation bounds the dimension mathematically. Together they
make $E_i$ a function whose computational cost per atom is *constant
in $N$*, which is what gives MLIPs their favourable scaling.

## 9.2.6 Invariant versus equivariant features: a worked example

To make the invariant/equivariant distinction concrete, consider a
single atom with two neighbours at relative positions
$\mathbf{r}_1, \mathbf{r}_2$. We will build two two-body features and
ask what each tells us about the geometry.

The *invariant* feature is the pair of distances and the cosine of the
included angle:

$$
\phi_\mathrm{inv}(\mathbf{r}_1, \mathbf{r}_2)
  = \big(\;\|\mathbf{r}_1\|,\; \|\mathbf{r}_2\|,\;
          \hat{\mathbf{r}}_1\cdot \hat{\mathbf{r}}_2\;\big)
  \in \mathbb{R}^3.
$$

Under any rotation $R$, these three numbers are unchanged. Good: that
is what invariance means. But the information lost is also clear: the
absolute orientation of the pair in space is gone.

The *equivariant* feature, in the spirit of NequIP, keeps the vectors
themselves but groups them by irreducible representation:

$$
\phi_\mathrm{eq}(\mathbf{r}_1, \mathbf{r}_2)
  = \Big(\;\|\mathbf{r}_1\|,\; \|\mathbf{r}_2\|\;\Big)^{\ell=0}
  \oplus \Big(\;\hat{\mathbf{r}}_1,\; \hat{\mathbf{r}}_2\;\Big)^{\ell=1}.
$$

The $\ell = 0$ part has two components and is scalar; the $\ell = 1$
part has six components (two vectors of three components each) and
transforms as $\mathbf{x} \mapsto R\mathbf{x}$ under rotation. The
total dimensionality is larger, eight numbers rather than three, but
no geometric information has been thrown away. From the equivariant
feature we can reconstruct the original $\mathbf{r}_1, \mathbf{r}_2$
up to a single global rotation; from the invariant feature we cannot.

This loss matters. Consider three neighbours at the vertices of an
equilateral triangle versus three neighbours at the vertices of an
isoceles triangle of the same edge length sum. Their distance
multisets are different — three equal versus two equal and one
different — and an invariant descriptor distinguishes them. But there
are pairs of configurations that have identical distance and angle
multisets while being *not related by rotation*. These are the so-called
*degenerate environments* of Pozdnyakov et al. (2020). Strictly
two-body and three-body invariant descriptors cannot distinguish
them; one needs either higher body-order or equivariant features
that propagate direction information through the network.

This is the deep reason equivariant networks outperform invariant
ones at small training-set sizes. They have a representation rich
enough to distinguish geometries that invariant descriptors collapse,
so they need fewer examples to learn the right mapping.

## 9.2.7 Putting the symmetries together

A correct MLIP architecture, then, must:

1. Operate on a neighbour list of relative vectors $\mathbf{r}_{ij}$
   (translation invariance).
2. Construct features that either are invariant under $\mathrm{O}(3)$
   or transform as irreps of $\mathrm{O}(3)$ (rotation
   invariance/equivariance), with the final energy a scalar.
3. Build atomic contributions $E_i$ as functions of an unordered set
   of neighbour features (permutation invariance), and sum atomic
   contributions to obtain the total energy.
4. Multiply every neighbour contribution by a smooth cutoff
   $f_\mathrm{c}(r_{ij})$ (smoothness).
5. Produce descriptors of fixed dimension independent of the number
   of neighbours (compactness).

The remainder of the chapter is a catalogue of architectures that
satisfy these constraints in different ways, and an exploration of
the trade-offs between them. We begin with the oldest and simplest
of the modern descriptors — Behler–Parrinello symmetry functions —
which makes every constraint visible in a few lines of code.

## 9.2.8 Locality and the cutoff radius

A symmetry not yet stated but implicit throughout is **locality**: the
atomic energy $E_i$ depends only on neighbours within a finite cutoff
$r_\mathrm{c}$. This is not a symmetry of nature — Coulomb
interactions are long-ranged — but a *modelling decision*. Three
considerations justify it.

**Physical.** In condensed matter at non-trivial density, screening
suppresses long-range correlations rapidly. The bonded and
short-range non-bonded contributions to the energy decay over a few
$\text{\AA}$; the residual long-range dispersion and Coulomb
contributions are small enough to be absorbed into smooth, slowly
varying corrections, or, in the case of strongly ionic systems, treated
by an explicit long-range term (a learnable charge model or an Ewald
sum on top of the local MLIP).

**Algorithmic.** Locality makes the per-atom evaluation cost
*constant in system size*. Without locality the cost would scale as
$O(N)$ per atom and $O(N^2)$ for the whole configuration, the same
unfavourable scaling that dooms naive DFT. The factor that turns MLIPs
into practical tools is exactly the linear total cost in $N$.

**Statistical.** Locality is a strong regulariser. By forcing $E_i$ to
be a function of only $\mathcal{O}(50)$ inputs (the neighbour list
within $r_\mathrm{c}$), we shrink the hypothesis class enormously and
require far less training data. A non-local architecture would have to
learn that distant atoms do not influence the energy — wasteful, when
we already know this is true to good approximation.

The choice of $r_\mathrm{c}$ is the main hyperparameter that
controls the locality assumption. Typical values are:

- $r_\mathrm{c} = 3.5$ to $4.0\,\text{\AA}$ for tightly bound
  molecular systems (organic chemistry, biomolecules) — captures
  first-shell covalent bonding and immediate hydrogen bonds.
- $r_\mathrm{c} = 5$ to $6\,\text{\AA}$ for general condensed
  matter — captures the first two coordination shells in most metals
  and oxides.
- $r_\mathrm{c} = 6$ to $8\,\text{\AA}$ for ionic and polar systems
  where electrostatic screening is partial.

For a message-passing network with $T$ layers, the *effective*
receptive field is $T \times r_\mathrm{c}$; a two-layer MACE network
with $r_\mathrm{c} = 5\,\text{\AA}$ sees information up to
$10\,\text{\AA}$ away from each atom, even though no single layer
extends beyond $5\,\text{\AA}$. This is the trick that lets modest
cutoffs cover surprisingly long-range correlations.

!!! tip "How to choose $r_\mathrm{c}$ in practice"
    A good empirical rule: pick the smallest cutoff such that the
    radial distribution function $g(r)$ of your system is essentially
    unity ($g(r) \approx 1.0 \pm 0.05$) at $r = r_\mathrm{c}$. This
    ensures the cutoff lies past the structured part of the local
    environment, in the bulk-like regime where neighbours look like a
    uniform sea. Pre-computing $g(r)$ from a short classical-force-field
    or DFT-MD trajectory costs almost nothing and saves considerable
    later debugging.

## 9.2.9 Why the symmetry constraints come *first*

A pedagogical note before we move on. Many introductions to MLIPs
present the architecture (Behler network, GAP, MACE) and then add
"oh, and it respects rotation/translation/permutation symmetry" as a
side remark. The presentation in this chapter inverts that emphasis
deliberately. We have stated the symmetry constraints first because:

1. **The constraints are universal.** Any correct MLIP — past, present,
   future — must satisfy them. Architectures come and go; the
   symmetries persist.
2. **The constraints determine the design space.** Once you have
   accepted that you must be translation-, rotation-, and
   permutation-invariant with smooth cutoffs and bounded body order,
   the set of possible architectures is greatly narrowed. The
   architectures we will study are *natural answers* to the
   constraints, not arbitrary inventions.
3. **The constraints are how you debug.** If your trained potential
   fails — energy drifts in NVE, the radial distribution function
   has spurious sharp features, MD blows up in five picoseconds —
   the first questions are: does it respect the symmetries exactly?
   Is the cutoff smooth to the required order? Are the descriptor
   gradients correct at boundary? Almost every MLIP bug traces back
   to a violated symmetry constraint.

Reader's compass: as you read §9.3 through §9.5, ask yourself for
each architecture which of the five symmetries it satisfies *exactly*,
which it satisfies *approximately*, and which it relies on the data
to teach. The honest architectures (modern equivariant networks)
satisfy all five exactly; older architectures cut corners on one or
two (e.g. older message-passing nets approximate rotation invariance
rather than enforcing it). The trend is clear: enforce more, learn
less, generalise better.
