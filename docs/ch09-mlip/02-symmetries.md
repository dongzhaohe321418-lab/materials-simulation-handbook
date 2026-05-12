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
