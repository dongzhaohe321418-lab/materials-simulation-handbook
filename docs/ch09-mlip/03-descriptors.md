# 9.3 Descriptors

A *descriptor* is a function that maps the local environment of an atom
to a fixed-length vector. It is the input to whatever regression model
predicts the atomic energy. The job of a good descriptor is to satisfy
the symmetries of §9.2 while losing as little information as possible
about the geometry. This section develops the three descriptor
families that organise the literature: Behler–Parrinello symmetry
functions, the Smooth Overlap of Atomic Positions (SOAP), and the
Atomic Cluster Expansion (ACE).

## 9.3.1 Behler–Parrinello symmetry functions

The simplest invariant descriptor consists of sums of radial and
angular functions over neighbours, evaluated on a fixed parameter
grid. Behler and Parrinello introduced two families.

### Why descriptor design matters

It is tempting to think that once we have a large enough neural
network, the descriptor is a mere bookkeeping detail. This view is
wrong, and the reason is informative.

A neural network trained on raw atomic Cartesians can in principle
learn the symmetries, but in practice it never does so exactly. The
symmetry group $\mathrm{O}(3) \times S_N \times \mathbb{R}^3$ has
infinite cardinality (the rotation group is continuous), and no finite
training set can densely sample its orbits. The network will achieve
some approximate symmetry on the training distribution and fail at
test time the moment a rotation takes it slightly out of distribution.

The proper engineering response is to *enforce* the symmetries in the
descriptor, so that the network never sees the same configuration in
two different rotational guises. This is the prior knowledge that the
descriptor encodes; the better the prior, the less data we need.

A second reason descriptors matter is *injectivity*. Two distinct
physical configurations should map to two distinct descriptors. If
the descriptor collapses inequivalent geometries to the same vector,
no regressor — however expressive — can recover the energy difference
between them. We will see in §9.2.6 and again at the end of §9.3.2
that some old descriptors fail injectivity in subtle ways (the
"degenerate environment" problem). The SOAP power spectrum, in
contrast, is provably complete up to inversion, which is one reason
GAP achieves the data-efficiency frontier it does.

The slogan: **a bad descriptor cannot be saved by a large network;
a good descriptor lets a small network suffice.**

### Radial symmetry functions

The radial function $G^2$ is a Gaussian centred at distance $r_s$:

$$
G^2_i(\eta, r_s) \;=\; \sum_{j \neq i}
  \exp\!\left[-\eta (r_{ij} - r_s)^2\right] f_\mathrm{c}(r_{ij}).
$$

Each choice of $(\eta, r_s)$ gives one scalar. By varying $r_s$ over
a grid of distances inside the cutoff, and $\eta$ over a grid of
widths, we obtain a vector of $N_\mathrm{rad}$ radial descriptors per
atom that together resolve the radial distribution of neighbours.

!!! note "Why this step? Completeness of the radial $G^2$ family"
    Why is the $(\eta, r_s)$ grid able to represent *any* radial
    distribution of neighbours? Because the family of Gaussians
    $\{e^{-\eta(r - r_s)^2}\}$ centred on a dense grid of $r_s$ values
    is a *complete basis* for sufficiently smooth functions on
    $[0, r_\mathrm{c}]$. Concretely, by varying $r_s$ continuously,
    summed Gaussians can approximate any smooth function (by the
    standard Gaussian-radial-basis-function approximation theorem). On
    a finite grid, the family approximates smooth functions to a
    resolution set by the Gaussian width $\sigma = 1/\sqrt{2\eta}$ and
    the grid spacing.

    The implication: the *fingerprint* of the radial density of
    neighbours is recovered, in principle, up to a resolution
    controlled by $\eta$ and the $r_s$ grid density. Choose them too
    coarse and you blur shells together; choose them too fine and you
    overfit on a single MD snapshot's noise. A useful default is to
    set the $r_s$ spacing equal to $1/\sqrt{2\eta}$, so that adjacent
    Gaussians overlap at the half-maximum and the basis tiles the
    radial axis without large gaps.

**Worked example.** Take $r_\mathrm{c} = 6.0\,\text{\AA}$, place
$r_s$ on a uniform grid of eight values $\{0.9, 1.6, 2.3, 3.0, 3.7,
4.4, 5.1, 5.8\}\,\text{\AA}$, and fix $\eta = 4.0\,\text{\AA}^{-2}$.
The width of each Gaussian is roughly
$\sigma = 1/\sqrt{2\eta} \approx 0.35\,\text{\AA}$, so neighbouring
shells in $r_s$ overlap mildly. For a copper atom in the fcc bulk, the
twelve first-shell neighbours sit at $r \approx 2.56\,\text{\AA}$ and
will contribute strongly to the $r_s = 2.3$ and $r_s = 3.0$
descriptors, weakly to the $r_s = 1.6$ descriptor, and negligibly to
the others. The second shell at $r \approx 3.62\,\text{\AA}$ contributes
to the $r_s = 3.7$ descriptor. The eight-component vector thus
fingerprints the radial structure of the environment.

Permutation invariance is automatic — the sum over $j$ does not depend
on neighbour ordering. Translation invariance is automatic — we use
only $r_{ij}$. Rotation invariance is automatic — distance is a scalar.
Smoothness follows from $f_\mathrm{c}$, which we choose to be the
Behler cosine envelope of §9.2.4. Compactness follows from the fixed
size of the $(\eta, r_s)$ grid.

### Angular symmetry functions

Two-body radial functions cannot distinguish environments with the
same radial distribution but different angular structure — square
planar versus tetrahedral, for instance. To capture three-body
information Behler introduced an angular function. The most widely
used form is

$$
G^4_i(\eta, \zeta, \lambda) = 2^{1-\zeta}
  \!\!\sum_{j,k \neq i,\, j \neq k}\!\!
  (1 + \lambda \cos\theta_{ijk})^\zeta \,
  e^{-\eta (r_{ij}^2 + r_{ik}^2 + r_{jk}^2)}
  f_\mathrm{c}(r_{ij})\, f_\mathrm{c}(r_{ik})\, f_\mathrm{c}(r_{jk}),
$$

where $\theta_{ijk}$ is the angle at atom $i$ subtended by $j$ and $k$,
and $\lambda \in \{+1, -1\}$ selects either small angles
($\lambda = +1$, peak at $\theta = 0$) or large angles
($\lambda = -1$, peak at $\theta = \pi$). The exponent $\zeta$ controls
the angular sharpness. Each choice of $(\eta, \zeta, \lambda)$
yields one scalar; a typical Behler grid uses $\sim 20$ angular
parameters per element pair.

!!! note "Anatomy of every parameter in $G^4$"
    Let us unpack the expression piece by piece, because the design
    choices in $G^4$ recur in many later descriptors.

    - $\lambda \in \{+1, -1\}$: the *angular handedness*. With
      $\lambda = +1$, the factor $(1 + \cos\theta)$ is maximised at
      $\theta = 0$ (collinear configurations) and vanishes at $\theta =
      \pi$. With $\lambda = -1$, the factor $(1 - \cos\theta)$ is
      maximised at $\theta = \pi$ (back-to-back). Both choices are
      used; both are needed to span the angular axis.
    - $\zeta \in \{1, 2, 4, 8, 16\}$: the *angular sharpness*. The
      $(1 + \lambda \cos\theta)^\zeta$ factor becomes more peaked as
      $\zeta$ increases — at $\zeta = 1$ it is a broad lobe, at
      $\zeta = 16$ it is a narrow spike. Varying $\zeta$ resolves
      angular features at different angular scales.
    - $\eta$: the *radial scale*, governing how rapidly the triple
      contribution decays with the squared distances $r_{ij}^2 +
      r_{ik}^2 + r_{jk}^2$. Large $\eta$ confines the triple to small
      bond-length triangles; small $\eta$ lets larger triangles
      contribute.
    - $2^{1-\zeta}$: a *normalisation prefactor* so that the maximum
      possible value of $(1 + \lambda \cos\theta)^\zeta /\, 2^{\zeta-1}$
      is unity — keeping different $\zeta$ values on a comparable
      scale.
    - $f_\mathrm{c}(r_{ij}) f_\mathrm{c}(r_{ik}) f_\mathrm{c}(r_{jk})$:
      *all three* edges of the triple must be cut off, not just the
      ones from atom $i$. If we omitted $f_\mathrm{c}(r_{jk})$, the
      descriptor would have a discontinuity when neighbours $j$ and
      $k$ drift past their mutual cutoff while remaining inside the
      cutoff sphere of $i$ — a subtle bug that produces NVE drift
      after a few picoseconds. Beginners frequently forget this third
      cutoff.

    A typical Behler angular grid for a single-element system uses
    $\eta \in \{0.001, 0.01, 0.05, 0.1\}\,\text{\AA}^{-2}$,
    $\zeta \in \{1, 2, 4, 16\}$, $\lambda \in \{\pm 1\}$, giving
    $4 \times 4 \times 2 = 32$ angular features per element pair. The
    full grid is a balance between coverage of the $(\theta, r)$
    angular-radial space and the cost of evaluating, storing, and
    training on a high-dimensional descriptor.

The $G^4$ descriptor is rotation-invariant (angles and distances are
scalars), permutation-invariant (the double sum runs over all pairs),
translation-invariant (we use relative quantities), smooth (every
distance is enveloped by $f_\mathrm{c}$), and compact (fixed grid size).

The full Behler–Parrinello descriptor is the concatenation of all
radial and angular components, partitioned by element type of the
neighbour. For a binary system AB with $N_\mathrm{rad} = 8$ radial
parameters and $N_\mathrm{ang} = 20$ angular parameters, the
descriptor of an A atom has $2 \times 8 + 3 \times 20 = 76$
components (8 for A–A radial, 8 for A–B radial, and 20 each for the
three element-pair channels AA, AB, BB on angles), where the factor
of three is the number of distinct element pairs that can occupy the
roles of $j$ and $k$.

### Implementation

A clean from-scratch implementation of the radial part follows.

```python
import numpy as np
from numpy.typing import NDArray

def cosine_cutoff(r: NDArray[np.float64], r_c: float) -> NDArray[np.float64]:
    """Behler cosine cutoff: smooth to zero at r = r_c."""
    fc = 0.5 * (np.cos(np.pi * r / r_c) + 1.0)
    return np.where(r < r_c, fc, 0.0)

def g2_descriptor(
    positions: NDArray[np.float64],
    cell: NDArray[np.float64],
    eta: NDArray[np.float64],
    r_s: NDArray[np.float64],
    r_c: float = 6.0,
) -> NDArray[np.float64]:
    """
    Compute Behler G^2 radial symmetry functions for every atom.

    Parameters
    ----------
    positions : (N, 3) array of atom Cartesian positions in angstrom.
    cell      : (3, 3) array of lattice vectors (rows).
    eta       : (K,) array of Gaussian widths in inverse angstrom squared.
    r_s       : (K,) array of Gaussian centres in angstrom.
    r_c       : radial cutoff in angstrom.

    Returns
    -------
    G : (N, K) array of descriptor values.
    """
    n_atoms = positions.shape[0]
    n_param = eta.shape[0]
    G = np.zeros((n_atoms, n_param))

    # Build minimum-image displacement matrix r_ij (vectorised).
    # For larger cells use a proper neighbour list with periodic images.
    inv_cell = np.linalg.inv(cell)
    frac = positions @ inv_cell             # (N, 3) fractional
    dfrac = frac[:, None, :] - frac[None, :, :]   # (N, N, 3)
    dfrac -= np.round(dfrac)                # minimum image
    dxyz = dfrac @ cell                     # (N, N, 3)
    r = np.linalg.norm(dxyz, axis=-1)       # (N, N)
    np.fill_diagonal(r, np.inf)             # exclude self

    mask = r < r_c                          # (N, N)
    fc = cosine_cutoff(r, r_c)

    # G_i^k = sum_{j != i} exp(-eta_k (r_ij - rs_k)^2) f_c(r_ij)
    # Loop over parameters K (small) to keep memory at O(N^2).
    for k in range(n_param):
        contrib = np.exp(-eta[k] * (r - r_s[k]) ** 2) * fc * mask
        G[:, k] = contrib.sum(axis=1)
    return G
```

A few features of this implementation are worth pointing out. The
neighbour search uses the minimum-image convention, which is correct
only when $r_\mathrm{c}$ is smaller than half the shortest cell
dimension. Production codes use a Verlet list with full periodic
replication; the principles are unchanged. The inner loop is over
descriptor parameters $K \approx 8$, not over atom pairs, so memory
scales as $O(N^2)$ and time scales as $O(K N^2)$. For cells larger
than a few hundred atoms one should switch to a linked-cell neighbour
list and replace the dense double loop with a sparse iteration.

Exercise 9.2 asks you to extend this to $G^4$ and verify
rotation invariance numerically.

### Strengths and limitations

Behler–Parrinello symmetry functions are interpretable, easy to
implement, and good enough for many production potentials. Their
limitations are:

1. The parameter grid $(\eta, r_s, \zeta, \lambda)$ must be chosen by
   hand. Bad choices leave gaps in coverage.
2. The descriptor size grows as $N_\mathrm{species}^2$ for radial
   functions and $N_\mathrm{species}^3$ for angular functions —
   awkward for systems with many element types.
3. Only two- and three-body information is captured. Distinguishing
   higher-order patterns requires additional ad hoc features.
4. The descriptor is invariant, not equivariant — the limitations of
   §9.2.6 apply.

SOAP and ACE address points (1) and (3) by switching to a *complete
basis* on the local environment, which we develop next.

## 9.3.2 SOAP — the Smooth Overlap of Atomic Positions

SOAP, introduced by Bartók, Kondor, and Csányi in 2010, replaces
hand-tuned symmetry functions with a systematic spherical-harmonic
expansion of the local atomic density.

### SOAP from intuition

Before equations, the intuition. Imagine spraying a thin coat of
Gaussian "paint" on every neighbour of atom $i$, with intensity
weighted by the smooth cutoff. The result is a continuous density
field $\rho_i(\mathbf{r})$ in the cutoff sphere — a smeared map of
where atoms are. Two environments that differ only by an overall
rotation produce $\rho_i$ fields that are rotated copies of each
other; if we extract any *rotation-invariant property* of the field,
those two environments will give the same descriptor by construction.

How do we extract rotation-invariant properties of a function on the
sphere? Expand in spherical harmonics. The coefficients $c_{n\ell m}$
of the expansion transform predictably under rotation — specifically,
as the $\ell$-th irrep at each fixed $(n, \ell)$. Inner products over
the $m$ index annihilate the rotation (because the rotation matrices
$D^{(\ell)}$ are unitary), leaving a rotation-invariant scalar.

The resulting descriptor — the *power spectrum* — captures the
two-point correlation of the atomic density on the cutoff sphere.
Higher correlations (three-point, four-point) give more discriminating
descriptors at the cost of more arithmetic; these are the *bispectrum*
and beyond.

### The local atomic density

Define the atomic density seen by atom $i$ as a sum of Gaussians
centred on each neighbour, weighted by the cutoff:

$$
\rho_i(\mathbf{r})
 = \sum_{j \in \mathcal{N}(i)}
   \exp\!\left[-\frac{(\mathbf{r} - \mathbf{r}_{ij})^2}{2\sigma_\mathrm{at}^2}\right]
   f_\mathrm{c}(r_{ij}).
$$

The width $\sigma_\mathrm{at}$ controls how sharply each atom is
resolved. SOAP treats $\rho_i$ as a function on $\mathbb{R}^3$ and
expands it in a basis of radial functions $g_n(r)$ (orthonormal on
$[0, r_\mathrm{c}]$) and spherical harmonics $Y_\ell^m(\hat{\mathbf{r}})$:

$$
\rho_i(\mathbf{r}) = \sum_{n,\ell,m} c^i_{n\ell m}\, g_n(r)\, Y_\ell^m(\hat{\mathbf{r}}),
\qquad
c^i_{n\ell m} = \int d^3\mathbf{r}\; \rho_i(\mathbf{r})\, g_n(r)\, Y_\ell^{m\,*}(\hat{\mathbf{r}}).
$$

The expansion coefficients $c^i_{n\ell m}$ are translation-invariant
(we shifted the density to atom $i$'s frame) and permutation-invariant
(the density is a sum over neighbours). They are *not*
rotation-invariant: under a rotation $R$ the spherical-harmonic index
$m$ rotates within each $\ell$ block via the Wigner D-matrix,

$$
c^i_{n\ell m} \mapsto \sum_{m'} D^{(\ell)}_{m m'}(R)\, c^i_{n\ell m'}.
$$

In other words, $(c^i_{n\ell m})_m$ transforms as the $\ell$-th irrep
of $\mathrm{O}(3)$ at each fixed $n$. The coefficients are themselves
*equivariant* features (this will be the entry point for NequIP in
§9.5).

### The rotation-invariant power spectrum

To make a rotation-invariant descriptor, take the inner product over
$m$:

$$
p^i_{n n' \ell} = \frac{1}{\sqrt{2\ell+1}}
  \sum_m c^i_{n\ell m}\, c^{i\,*}_{n'\ell m}.
$$

This is the SOAP *power spectrum*. To see that it is rotation-invariant,
apply a rotation and use the unitarity of the Wigner D-matrices:

$$
\sum_m \!\!\bigg[\sum_{m_1} D^{(\ell)}_{m m_1} c^i_{n\ell m_1}\bigg]
       \bigg[\sum_{m_2} D^{(\ell)\,*}_{m m_2} c^{i\,*}_{n'\ell m_2}\bigg]
 = \sum_{m_1 m_2} \delta_{m_1 m_2} c^i_{n\ell m_1} c^{i\,*}_{n'\ell m_2}
 = \sum_{m_1} c^i_{n\ell m_1} c^{i\,*}_{n' \ell m_1},
$$

which is the original $p^i_{nn'\ell}$ (up to the normalisation). The
detailed proof appears as Exercise 9.3.

The dimension of the power spectrum is
$N_\mathrm{rad}(N_\mathrm{rad}+1)/2 \times (\ell_\mathrm{max}+1)$
once one accounts for the symmetry $p_{nn'\ell} = p_{n'n\ell}^*$.
A typical choice is $N_\mathrm{rad} = 8$ and $\ell_\mathrm{max} = 6$,
giving $36 \times 7 = 252$ scalar descriptors per atom.

!!! note "Why this step? Unitarity of Wigner D-matrices"
    The proof above relied on the unitarity identity
    $\sum_m D^{(\ell)}_{m m_1}(R) D^{(\ell)*}_{m m_2}(R) =
    \delta_{m_1 m_2}$. This is not magic; it follows because $D^{(\ell)}$
    is the *unitary matrix* of the $\ell$-th irrep of the rotation
    group acting on a $(2\ell+1)$-dimensional vector space. Unitarity
    means $D^{(\ell)\dagger} D^{(\ell)} = I$, which when written
    component-wise is exactly the identity above.

    More concretely, for the $\ell = 1$ irrep, $D^{(1)}(R)$ is the
    $3 \times 3$ rotation matrix itself (in a particular spherical
    basis); unitarity is the standard $R^\top R = I$ for rotation
    matrices. For higher $\ell$, $D^{(\ell)}$ is a $(2\ell+1) \times
    (2\ell+1)$ matrix whose explicit form can be derived from the
    Euler-angle parameterisation of $R$ (the "Wigner d-matrix"
    formulae). For our purposes we only need that they are unitary,
    which is a property of *any* matrix irrep of a compact group.

    The take-home: rotation invariance of inner products of irrep
    components is *automatic*. We will use this fact again in §9.5
    to argue that the tensor-product MACE layer preserves
    equivariance — same unitarity, same algebraic identity.

### A 40-line SOAP implementation

To make the derivation concrete, here is a from-scratch SOAP
implementation for a single atom with two neighbours, in NumPy only.
It computes the $\ell$-only power spectrum (summed over radial
indices) for compactness and verifies invariance under rotation by
applying a random rotation to the coordinates and recomputing.

```python
import numpy as np
from numpy.typing import NDArray
from scipy.special import sph_harm     # Y_l^m(theta, phi)

def soap_power_spectrum(
    neighbours: NDArray[np.float64],   # (M, 3) relative positions
    r_c: float = 4.0,
    sigma: float = 0.5,
    n_rad: int = 4,
    l_max: int = 4,
) -> NDArray[np.float64]:
    """Compute SOAP power spectrum p_{nn'l} summed over m.

    Atomic density: rho(r) = sum_j G(r - r_j; sigma) f_c(r_j)
    Expand: c_{nlm} = int rho(r) g_n(r) Y_l^m*(r) d^3 r
    Power spectrum: p_{nn'l} = sum_m c_{nlm}^* c_{n'lm}
    """
    M = neighbours.shape[0]
    # Cosine cutoff per neighbour.
    r = np.linalg.norm(neighbours, axis=1)
    fc = 0.5 * (np.cos(np.pi * r / r_c) + 1.0) * (r < r_c)
    # Spherical coordinates of each neighbour.
    theta = np.arccos(neighbours[:, 2] / np.where(r > 0, r, 1.0))
    phi = np.arctan2(neighbours[:, 1], neighbours[:, 0])
    # Gaussian radial basis: g_n(r) = exp(-(r - r_n)^2 / (2 sigma^2))
    r_n = np.linspace(0.5, r_c - 0.5, n_rad)
    g = np.exp(-(r[None, :] - r_n[:, None]) ** 2 / (2 * sigma**2))
    # c_{nlm} = sum_j g_n(r_j) Y_l^m*(theta_j, phi_j) f_c(r_j)
    p = np.zeros((n_rad, n_rad, l_max + 1))
    for l in range(l_max + 1):
        c_l = np.zeros((n_rad, 2 * l + 1), dtype=complex)
        for m_idx, m in enumerate(range(-l, l + 1)):
            Y = sph_harm(m, l, phi, theta)                # (M,)
            c_l[:, m_idx] = (g * fc[None, :] * np.conj(Y)[None, :]).sum(axis=1)
        # p_{nn'l} = sum_m c_{nlm}^* c_{n'lm}
        p[:, :, l] = np.real(c_l @ c_l.conj().T)
    return p

# Demonstration: invariance under random rotation.
rng = np.random.default_rng(0)
neigh = rng.normal(size=(5, 3)) * 1.5   # 5 neighbours
p0 = soap_power_spectrum(neigh)
# Apply a random rotation.
from scipy.spatial.transform import Rotation
R = Rotation.random(random_state=rng).as_matrix()
neigh_rot = neigh @ R.T
p1 = soap_power_spectrum(neigh_rot)
assert np.allclose(p0, p1, atol=1e-10), "SOAP failed rotation invariance!"
print("SOAP power spectrum is rotation invariant up to 1e-10.")
```

A few comments on the implementation. The `scipy.special.sph_harm`
convention has `(m, l, phi, theta)` argument order — the *azimuthal*
angle first, then the polar angle — which surprises many readers; the
complex-conjugate-then-multiply pattern in `c_l @ c_l.conj().T`
implements the $m$-contraction that gives rotation invariance. The
test under a random rotation confirms invariance to machine precision
($10^{-10}$), as the algebra promises.

For production work one uses an optimised SOAP routine — `dscribe`,
`librascal`, or `quippy` — that vectorises over atoms and handles
periodic boundary conditions, multi-element systems, and orthonormal
radial bases. But the 40-line core captures every essential idea:
neighbour density on the sphere, spherical-harmonic expansion,
$m$-contraction.

### From the power spectrum to higher-order invariants

The power spectrum is the *two-point correlation* of the atomic
density on the unit sphere — informally, "how similar is the density
at two random directions, averaged over all directions". Two
configurations that share the same two-point correlation can still
differ in their three-point correlation (the *bispectrum*),
four-point correlation, and so on.

The bispectrum is the triple contraction of expansion coefficients
with a Clebsch–Gordan symbol:

$$
b^i_{n_1 n_2 n_3, \ell_1 \ell_2 \ell_3}
   = \sum_{m_1 m_2 m_3} C^{\ell_1 \ell_2 \ell_3}_{m_1 m_2 m_3}
     c^i_{n_1 \ell_1 m_1}\,
     c^i_{n_2 \ell_2 m_2}\,
     c^i_{n_3 \ell_3 m_3},
$$

invariant under rotation because the Clebsch–Gordan symbol contracts
three irrep indices to a scalar. The bispectrum is the basis of the
*SNAP* (Spectral Neighbour Analysis Potential) descriptor.

As body order grows, more invariants become possible, and at each
body order the *complete* set of invariants (those generated by
contracting coefficients with the appropriate $n$-point Clebsch–Gordan
generalisation) gives an injective fingerprint of the local
environment. This is the point of view that ACE makes systematic in
§9.3.3.

### Why SOAP is complete

The key property of SOAP, proven by Bartók et al., is *completeness up
to inversion*: two environments with the same SOAP power spectrum are
identical up to an overall rotation and inversion (and the
permutation-by-element of neighbours), provided $N_\mathrm{rad}$ and
$\ell_\mathrm{max}$ are large enough. This contrasts with
Behler–Parrinello descriptors, which can collapse genuinely distinct
environments (the degenerate-environment issue of §9.2.6).

Completeness explains SOAP's success on small datasets. When paired
with a Gaussian process (§9.4), a SOAP-based GAP potential typically
reaches the data-efficiency frontier among invariant methods. The
power spectrum is, however, only the second-order term in a series:
the *bispectrum*, obtained by triple-contracting coefficients with a
Clebsch–Gordan symbol, captures three-body information; higher
contractions give four-body and beyond.

### Per-element channels

A practical wrinkle is that the atomic density above conflates all
neighbour species. In a heterogeneous system one keeps a separate
density per element,

$$
\rho_i^{(\alpha)}(\mathbf{r}) =
 \sum_{j \in \mathcal{N}(i),\, Z_j = \alpha}
 \exp[\cdots] f_\mathrm{c}(r_{ij}),
$$

and the power spectrum carries an extra pair of element indices:
$p^i_{(\alpha,\alpha'),nn'\ell}$. The descriptor dimension grows as
$N_\mathrm{species}^2$, the same as in BP but with no choice of
hand-tuned parameters.

## 9.3.3 ACE — the Atomic Cluster Expansion

The Atomic Cluster Expansion, introduced by Drautz in 2019, recasts
the SOAP-style construction as a *systematic body-order expansion*
of the atomic energy. ACE is the mathematical foundation of MACE
(§9.5) and provides a clean language for comparing descriptors.

### ACE in one paragraph of intuition

Imagine writing the energy as a polynomial in *neighbour density
features*. The single-particle features $A^i_v$ are the SOAP
coefficients $c^i_{n\ell m}$: a sum over neighbours of an angular
spherical harmonic times a radial function. Linear combinations of
$A^i_v$ give two-body terms. Products of two $A^i_v$'s give
three-body terms (each $A$ is a sum over one neighbour, so the
product is a double sum over two neighbours). Products of three
$A^i_v$'s give four-body terms. And so on. By taking polynomials of
the $A^i_v$, we generate the entire many-body hierarchy from a single
elegant construction.

The subtlety is that products of equivariant features
($A^i_{n\ell m}$ has an $m$ index that rotates) must be contracted
with appropriate Clebsch–Gordan symbols to obtain rotation
invariants. This is where the machinery of §9.5 enters: the same
contractions that make ACE features rotation-invariant also let MACE
build equivariant features of any chosen output irrep.

### The body-order expansion

Write the atomic energy as a sum of contributions from increasing
numbers of neighbours:

$$
E_i = V_0
    + \sum_j V_1(\mathbf{r}_{ij})
    + \sum_{j<k} V_2(\mathbf{r}_{ij}, \mathbf{r}_{ik})
    + \sum_{j<k<l} V_3(\mathbf{r}_{ij}, \mathbf{r}_{ik}, \mathbf{r}_{il})
    + \cdots
$$

This is the *cluster expansion*. The two-body term $V_1$ depends on
one neighbour, the three-body term $V_2$ on two, the four-body term
$V_3$ on three. The expansion converges rapidly for most chemistries:
truncating at four- or five-body order suffices.

A direct fit of $V_2, V_3$ as multivariate functions is intractable —
the dimension grows combinatorially. ACE makes the fit tractable in
two steps.

### Step 1: a one-particle basis

Expand each one-body function on a basis of the form

$$
\phi_v(\mathbf{r}) = R_{nl}(r)\, Y_\ell^m(\hat{\mathbf{r}}),
$$

with $v = (n, \ell, m)$ a composite index. The radial part $R_{nl}$
is a smooth function on $[0, r_\mathrm{c}]$ vanishing at the cutoff,
typically a polynomial or a sum of Bessel functions. The angular part
is the spherical harmonic, exactly as in SOAP.

The atomic basis function for atom $i$ is the sum over neighbours,

$$
A^i_v = \sum_{j \in \mathcal{N}(i)} \phi_v(\mathbf{r}_{ij}) f_\mathrm{c}(r_{ij}).
$$

These $A^i_v$ are precisely the SOAP coefficients $c^i_{n\ell m}$, up
to choices of normalisation and radial basis. So far the construction
is identical to SOAP.

### Step 2: products and contraction

The body-order-$N$ basis function is a *product* of $N-1$ one-particle
functions, summed over distinct neighbours, then contracted with a
generalised Clebsch–Gordan symbol $\mathcal{C}$ to produce a
rotation-invariant scalar:

$$
B^i_{\mathbf{v}, L}
 = \sum_{\mathbf{m}} \mathcal{C}^{\,L}_{\mathbf{v}, \mathbf{m}}
   \prod_{\alpha=1}^{N-1} A^i_{v_\alpha m_\alpha},
$$

where $\mathbf{v} = (v_1, \dots, v_{N-1})$ collects the composite
indices and $L$ is the total angular-momentum label of the product
representation. Choosing $L = 0$ gives an invariant; non-zero $L$
gives equivariants of the corresponding order, which is the route
MACE exploits.

The crucial observation is that **products of densities give
many-body invariants automatically**. The $N$-body ACE feature
$B^i_{\mathbf{v}, 0}$ depends on $N - 1$ neighbours of atom $i$, and
the full set $\{B^i_{\mathbf{v}, 0}\}_{|\mathbf{v}| = N-1}$ forms a
complete basis for the $N$-body cluster term $V_{N-1}$ in the body-order
expansion.

### Why products of densities give many-body terms

A useful sanity check: how do products of two-body density coefficients
$A^i_v = \sum_j \phi_v(\mathbf{r}_{ij}) f_\mathrm{c}(r_{ij})$ become
three- and four-body terms?

Consider a single $A$, which is a sum over one neighbour: $A^i_v =
\phi_v(\mathbf{r}_{i1}) + \phi_v(\mathbf{r}_{i2}) + \cdots$. The
product of two such sums is

$$
A^i_{v_1} A^i_{v_2}
   = \sum_{j_1, j_2} \phi_{v_1}(\mathbf{r}_{i j_1})\,
                      \phi_{v_2}(\mathbf{r}_{i j_2}).
$$

The double sum runs over *all pairs* of neighbours $(j_1, j_2)$,
including the diagonal $j_1 = j_2$. The off-diagonal terms $j_1 \neq
j_2$ are three-body contributions (atom $i$ plus two distinct
neighbours); the diagonal terms are two-body "self-products" that
can be absorbed into the two-body channel. Similarly, the product
$A^i_{v_1} A^i_{v_2} A^i_{v_3}$ contains four-body terms (atom $i$
plus three distinct neighbours), three-body terms (when two indices
coincide), and so on. By taking polynomials of arbitrary degree, we
generate the full body-order hierarchy.

The contraction with a generalised Clebsch–Gordan symbol is what
projects this product onto a specific output irrep — typically the
scalar irrep $L = 0$ for invariant energy contributions. The full
many-body basis is therefore parameterised by the multi-index
$\mathbf{v} = (v_1, \dots, v_{\nu})$ specifying which $\phi_{v_\alpha}$
single-particle functions are multiplied, plus the choice of coupling
scheme implied by $\mathcal{C}$.

### Polynomial-in-density expansion

The ACE atomic energy is then a *polynomial* in the basis functions:

$$
E_i = \sum_{\mathbf{v}} c_{\mathbf{v}}\, B^i_{\mathbf{v}, 0}.
$$

The coefficients $c_{\mathbf{v}}$ are the parameters of the model.
Linear ACE — that is, linear regression on the basis $B^i_{\mathbf{v}, 0}$
— is already a competitive interatomic potential and is the basis of
performant codes such as PACEMAKER.

### A worked example of the ACE basis at body order 3

To see the body-order construction explicitly, let us write down the
$\nu = 2$ (three-body, atom $i$ plus two neighbours) ACE features for
a single carbon environment.

The two-body density coefficients are $A^i_{n\ell m}$ for some chosen
$(n, \ell, m)$. We pick $n = 0$ (a single radial channel), and let
$\ell$ run over $\{0, 1, 2\}$. So we have $1 + 3 + 5 = 9$ single-particle
features.

A three-body ACE feature is a product $A^i_{n_1 \ell_1 m_1} A^i_{n_2
\ell_2 m_2}$, contracted with a Clebsch–Gordan symbol to produce a
scalar (output irrep $L = 0$):

$$
B^i_{(n_1 \ell_1, n_2 \ell_2), L=0}
   = \sum_{m_1 m_2} C^{0\,0}_{\ell_1 m_1; \ell_2 m_2}
     A^i_{n_1 \ell_1 m_1} A^i_{n_2 \ell_2 m_2}.
$$

The selection rule for the Clebsch–Gordan symbol coupling to $L = 0$
is $\ell_1 = \ell_2$. So the surviving features have the form

$$
B^i_{\ell, \ell} = \sum_m (-1)^m A^i_{0 \ell m} A^i_{0 \ell, -m}
\quad \text{(up to normalisation)},
$$

which is the SOAP $\ell$-component power spectrum. The ACE three-body
basis *is* the SOAP power spectrum at this body order; the body-order
expansion makes the equivalence transparent.

At $\nu = 3$ (four-body), we get triple products $A_{n_1 \ell_1 m_1}
A_{n_2 \ell_2 m_2} A_{n_3 \ell_3 m_3}$ contracted to $L = 0$. The
selection rules now involve triangle inequalities on $(\ell_1,
\ell_2, \ell_3)$ — the same algebra that arises in atomic physics
when coupling three angular momenta. This is the bispectrum.

The hierarchy continues. At each body order $\nu$, ACE produces a
complete (over-)basis of $\nu$-body invariants; the redundancy
between different couplings of the same physical content is removed
by a Gram-Schmidt-like orthogonalisation (or, in practice, by L1
regularisation that prunes unused features). What MACE accomplishes is
to *learn* this orthogonalisation jointly with the regression, using
small message-passing networks to produce a sparse low-dimensional
basis that captures the most useful many-body features.

### Connection to MACE

Linear ACE has one weakness: the number of basis functions grows
combinatorially with body order and angular cutoff. To capture
four-body interactions with $\ell_\mathrm{max} = 3$ one needs tens of
thousands of features. MACE solves this by *learning* the basis
hierarchically with a small message-passing network. The starting
features are the equivariant atomic densities $A^i_{v}$ (the SOAP
coefficients). A tensor-product layer (see §9.5) combines features of
irrep order $\ell_1$ and $\ell_2$ via Clebsch–Gordan coupling to
produce features of irrep order $\ell$. After $K$ such layers the
features encode body-order up to $\nu = K(\nu_\mathrm{layer} - 1) + 1$,
where $\nu_\mathrm{layer}$ is the in-layer correlation order. With
$\nu_\mathrm{layer} = 3$ and $K = 2$ layers, MACE reaches four-body
correlations.

The point of view we should carry forward is:

- Behler symmetry functions are a *fixed*, hand-designed set of
  two- and three-body invariants.
- SOAP is a *complete* set of two-body invariants on a chosen
  spherical-harmonic basis.
- ACE generalises SOAP to arbitrary body order via products of
  density coefficients.
- MACE *learns* a sparse, low-dimensional version of the ACE basis
  via message-passing tensor-product layers.

Each successive level pays additional implementation complexity to
buy either richer geometry (more body order) or richer information
flow (equivariance and message passing). The pay-off is data
efficiency: MACE reaches a given accuracy with roughly twenty times
fewer DFT calculations than a Behler–Parrinello network.

## 9.3.4 Choosing a descriptor in practice

For the practitioner choosing among descriptors, the rough heuristics
are:

- **Behler symmetry functions** are appropriate for one- or two-element
  systems with modest accuracy requirements ($\sim 10\,\mathrm{meV}/\text{atom}$),
  large training sets, and a need for a simple, transparent
  implementation. The descriptor itself is one screenful of code.
- **SOAP** is the right choice for Bayesian potentials with built-in
  uncertainty (GAP, §9.4), and for systems where you want a
  parameter-free, complete invariant descriptor with strong
  theoretical guarantees.
- **ACE / MACE** is the right choice when accuracy and data
  efficiency matter most, particularly for chemistries with many
  elements or for foundation-model fine-tuning. This is what we
  train in §9.6.

The rest of the chapter develops the regression models that consume
these descriptors. Section 9.4 covers Behler–Parrinello networks and
GAP; §9.5 develops the equivariant message-passing networks that have
come to dominate the field.

## 9.3.5 A worked numerical comparison

To anchor the descriptors in numerical reality, consider a single
methane CH$_4$ molecule at equilibrium geometry: a central carbon at
the origin, four hydrogens at the corners of a regular tetrahedron at
distance $1.09\,\text{\AA}$. We compute three descriptors for the
central carbon and report a representative scalar from each.

**Behler $G^2$ at $r_s = 1.1\,\text{\AA}, \eta = 4.0\,\text{\AA}^{-2}$.**

$$
G^2_C = \sum_{j=1}^4 e^{-4 (1.09 - 1.1)^2} f_\mathrm{c}(1.09)
      = 4 \times e^{-0.0004} \times f_\mathrm{c}(1.09).
$$

With $r_\mathrm{c} = 4.0\,\text{\AA}$, $f_\mathrm{c}(1.09) = 0.5(\cos(\pi
\times 1.09 / 4.0) + 1) \approx 0.815$. So $G^2_C \approx 4 \times
1.00 \times 0.815 = 3.26$.

**Behler $G^4$ at $\eta = 0.01, \zeta = 4, \lambda = -1$.**

A tetrahedron has six $\binom{4}{2} = 6$ pairs of H–C–H angles, each
equal to $\theta_\mathrm{tet} = \arccos(-1/3) \approx 109.47^\circ$,
$\cos\theta = -1/3$. So $(1 + \lambda \cos\theta)^\zeta = (1 - (-1/3))^4
= (4/3)^4 = 256/81 \approx 3.16$. The radial exponential is
$\exp(-0.01 \times (1.09^2 + 1.09^2 + 1.78^2)) \approx \exp(-0.056)
\approx 0.945$. Three cutoff factors at $r = 1.09$ ($f_\mathrm{c}
\approx 0.815$) and $r = 1.78$ for the H–H distance ($f_\mathrm{c}
\approx 0.616$) give $0.815 \times 0.815 \times 0.616 \approx 0.409$.

Combining: $G^4_C \approx 2^{1-4} \times 6 \times 2 \times 3.16 \times
0.945 \times 0.409 \approx 1.83$. (The factor of $2$ is the
double-counting in the unrestricted $\sum_{j, k \neq i, j \neq k}$ sum
of ordered pairs.)

**SOAP power spectrum, $n = n' = 0, \ell = 0$.**

With $\sigma_\mathrm{at} = 0.5\,\text{\AA}$ and $g_0(r) = 1$ (a
constant radial basis for illustration), and $Y_0^0 = 1/\sqrt{4\pi}$,

$$
c^C_{000} = \sum_{j=1}^4 \int d^3 \mathbf{r}\;
   e^{-(\mathbf{r} - \mathbf{r}_{Cj})^2 / 2\sigma^2}
   \cdot 1 \cdot \frac{1}{\sqrt{4\pi}}\, f_\mathrm{c}(r_{Cj}).
$$

The Gaussian integral $\int d^3\mathbf{r}\, e^{-(\mathbf{r} -
\mathbf{r}_0)^2/2\sigma^2} = (2\pi\sigma^2)^{3/2}$, so $c^C_{000} =
4 \times (2\pi \times 0.25)^{3/2} / \sqrt{4\pi} \times 0.815 \approx
4 \times 1.97 \times 0.282 \times 0.815 \approx 1.81$. Then
$p^C_{000} = |c^C_{000}|^2 / \sqrt{1} \approx 3.28$.

The Behler $G^2$ and the SOAP $p_{000}$ are both rough measures of
"how many neighbours are within a Gaussian-broadened shell". They
return numerically similar values (3.26 vs 3.28) because they are
both, in essence, integrating the neighbour density against a smooth
radial kernel. The Behler $G^4$ captures the angular structure of the
tetrahedron explicitly.

The whole point of these descriptor frameworks is that we never have
to do this arithmetic by hand at scale — `dscribe`, `librascal`,
`mace-torch`, and similar packages compute them vectorised over
configurations on a GPU in milliseconds. But knowing what the numbers
*mean* is what lets you diagnose problems when descriptor-level bugs
arise.

## 9.3.6 Forward references

In §9.4 we will plug the Behler symmetry functions and SOAP
descriptors into specific regressors:

- **Behler symmetry functions $+$ small neural network** = the
  original Behler–Parrinello architecture (§9.4.1).
- **SOAP power spectrum $+$ Gaussian process** = the GAP framework
  (§9.4.4).

In §9.5 we will revisit the ACE formalism in its full glory: the
single-particle basis $A^i_v$ becomes the message-passing
"two-body message" of a NequIP/MACE layer, the symmetric tensor
product replaces the explicit products of $A$'s, and the entire
construction becomes end-to-end trainable. The conceptual continuity
from Behler to MACE — increasingly powerful invariants, increasingly
deep regressors — should be evident by the time we reach §9.5.4.

A final point on terminology that has confused many readers. In the
SOAP literature, "power spectrum" often refers to the *full*
$p_{nn'\ell}$ tensor (a vector of $\sim 250$ scalars), and sometimes
just to its dimension-summed scalar projections. In the ACE
literature, $B^i_\mathbf{v}$ are called "basis functions" and the
specific contraction conventions matter. When reading papers, always
check the dimensionality of the descriptors quoted, not just the name.
