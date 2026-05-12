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

### Polynomial-in-density expansion

The ACE atomic energy is then a *polynomial* in the basis functions:

$$
E_i = \sum_{\mathbf{v}} c_{\mathbf{v}}\, B^i_{\mathbf{v}, 0}.
$$

The coefficients $c_{\mathbf{v}}$ are the parameters of the model.
Linear ACE — that is, linear regression on the basis $B^i_{\mathbf{v}, 0}$
— is already a competitive interatomic potential and is the basis of
performant codes such as PACEMAKER.

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
