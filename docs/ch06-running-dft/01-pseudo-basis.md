# 6.1 Pseudopotentials and basis sets

A solid is, atomistically, a lattice of nuclei surrounded by electrons. To solve the Kohn-Sham equations of [Chapter 5](../ch05-dft/index.md) we need two things: a description of the potential each electron feels from the nuclei (and from itself, mean-field), and a basis in which to expand the Kohn-Sham orbitals. Choosing both is the first practical decision in any DFT calculation, and the rest of the chapter depends on it.

!!! tip "Why this whole chapter exists"
    The Kohn-Sham equations of Chapter 5 are exact in principle. In practice you cannot solve them — they are infinitely-many coupled three-dimensional partial differential equations in real space. To make progress you must (a) reduce the number of electrons you actually track (pseudopotentials), (b) replace continuous wavefunctions by a finite set of basis functions (plane waves), and (c) replace the BZ integral by a discrete sum (k-point grid). Each step trades exactness for tractability. The art of solid-state DFT is making each trade-off small enough to ignore while keeping the calculation runnable.
    
    The picture in your head: an atom has a tiny core of fast electrons and a fluffy halo of slow valence electrons. We freeze the core, smooth the valence, and replace the whole atom by an "effective pseudo-atom" that valence electrons see. Then we tile space with plane waves — Fourier modes of the cell — and pick k-points by symmetry. By the end of this section you have three knobs (pseudopotential, $E_\mathrm{cut}$, k-grid) and a recipe to set them.

## 6.1.1 Why we do not need the core electrons

Consider silicon. It has 14 electrons in the configuration $1s^2\, 2s^2\, 2p^6\, 3s^2\, 3p^2$. When silicon atoms come together to form crystalline Si, only the four $3s$/$3p$ valence electrons participate in bonding. The ten core electrons sit tightly bound around the nucleus, do not rearrange when neighbours change, and contribute essentially nothing to the chemistry.

Yet they cost an enormous amount to compute. The $1s$ orbital of silicon is highly localised: it has a maximum at about $0.05\,a_0$ from the nucleus and oscillates rapidly. To represent it on a plane-wave basis would require Fourier components up to wavevectors $|\mathbf{G}| \sim 2\pi / 0.05\,a_0 \approx 100\,a_0^{-1}$, corresponding to plane-wave cutoffs of order $10^4$ Ry. For comparison, valence physics is converged by 30-50 Ry.

The **frozen-core approximation** is the observation that we can fix the core electrons at their atomic-state density and treat them as part of the ionic environment. The valence electrons then see an *effective* potential: the bare nuclear potential, plus the Hartree and exchange-correlation contributions from the core. This effective potential is the **pseudopotential**.

A pseudopotential must do two things:

1. **Replace the singular $-Z/r$ Coulomb potential** by something smooth inside a cutoff radius $r_c$, so that valence wavefunctions (the *pseudo-wavefunctions*) become nodeless and smooth there too.
2. **Reproduce, outside $r_c$, the scattering of valence electrons from the true all-electron atom** at the energies relevant to bonding.

Outside $r_c$ the pseudo-wavefunction equals the all-electron wavefunction; inside $r_c$ it is smooth and easy to expand in plane waves. The price is that the pseudopotential is non-local (it acts differently on different angular momenta $\ell$, since the all-electron $\ell = 0, 1, 2, \ldots$ channels scatter differently) and energy-dependent in principle (in practice we evaluate at a reference energy).

!!! note "Semi-core states"
    For some elements the boundary between "core" and "valence" is fuzzy. Transition metals such as Ti or Zn often need their $3s$ and $3p$ "semi-core" electrons treated as valence, because they overlap with neighbours and respond to chemistry. Pseudopotential libraries provide multiple variants; reading the metadata is not optional.

That is the whole idea: freeze the core, smooth the valence, and hand the valence electrons an effective pseudo-atom. In practice three families of pseudopotential — norm-conserving, ultrasoft, and PAW — implement this idea with different trade-offs between smoothness and bookkeeping. You do *not* need those distinctions to run your first calculation: a curated library (§6.1.5) picks a good pseudopotential and its cutoffs for you. The detailed NC/USPP/PAW comparison is therefore deferred to the end of this section (§6.1.6); on a first pass you may skip ahead to [§6.2](02-first-qe.md) and come back to it once you have a running calculation to anchor the theory.

## 6.1.2 Choosing a basis: plane waves vs everything else

The Kohn-Sham orbital must be expanded in some basis $\{\chi_\alpha(\mathbf{r})\}$:

$$ \psi_{n\mathbf{k}}(\mathbf{r}) = \sum_\alpha c_{n\mathbf{k},\alpha}\, \chi_\alpha(\mathbf{r}). $$

Three families are common.

### Plane waves

!!! note "Why this step? — Plane waves diagonalise translation"
    The reason plane waves are the *natural* basis for a periodic solid is not aesthetic; it is operational. The single-electron Hamiltonian $\hat{H} = \hat{T} + \hat{V}_\mathrm{ext} + \hat{V}_\mathrm{H} + \hat{V}_\mathrm{xc}$ commutes with every lattice translation $\hat{T}_\mathbf{R}$ (where $\mathbf{R}$ is a Bravais lattice vector). Two commuting Hermitian operators share an eigenbasis. The eigenfunctions of $\hat{T}_\mathbf{R}$ are exactly the Bloch functions: $\hat{T}_\mathbf{R}\,\psi_\mathbf{k}(\mathbf{r}) = e^{-i\mathbf{k}\cdot\mathbf{R}}\,\psi_\mathbf{k}(\mathbf{r})$. Plane waves $e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}$ are the canonical realisation of these eigenstates. Choosing plane waves *automatically* block-diagonalises $\hat{H}$ in the crystal momentum $\mathbf{k}$, so the original $\infty$-dimensional problem decomposes into one finite problem per $\mathbf{k}$-point.

#### A derivation of Bloch's theorem from translational symmetry

Let $\hat{T}_\mathbf{R}$ be the operator that shifts a wavefunction by a Bravais lattice vector $\mathbf{R}$: $(\hat{T}_\mathbf{R}\psi)(\mathbf{r}) = \psi(\mathbf{r}+\mathbf{R})$. Lattice translations form an abelian group: $\hat{T}_\mathbf{R}\hat{T}_{\mathbf{R}'} = \hat{T}_{\mathbf{R}+\mathbf{R}'}$. Because the Kohn-Sham potential is lattice-periodic, $[\hat{H},\hat{T}_\mathbf{R}] = 0$ for all $\mathbf{R}$.

The eigenvalues of any unitary operator have unit modulus, so we can write
$$\hat{T}_\mathbf{R}\,\psi(\mathbf{r}) = e^{-i\mathbf{k}\cdot\mathbf{R}}\,\psi(\mathbf{r})$$
for some real vector $\mathbf{k}$. The group structure $\hat{T}_\mathbf{R}\hat{T}_{\mathbf{R}'} = \hat{T}_{\mathbf{R}+\mathbf{R}'}$ forces the dependence on $\mathbf{R}$ to be linear in the exponent; the value of $\mathbf{k}$ labels the simultaneous eigenstate.

Now define $u_\mathbf{k}(\mathbf{r}) \equiv e^{-i\mathbf{k}\cdot\mathbf{r}}\psi_\mathbf{k}(\mathbf{r})$. Apply $\hat{T}_\mathbf{R}$:
$$u_\mathbf{k}(\mathbf{r}+\mathbf{R}) = e^{-i\mathbf{k}\cdot(\mathbf{r}+\mathbf{R})}\psi_\mathbf{k}(\mathbf{r}+\mathbf{R}) = e^{-i\mathbf{k}\cdot(\mathbf{r}+\mathbf{R})}\,e^{i\mathbf{k}\cdot\mathbf{R}}\psi_\mathbf{k}(\mathbf{r}) = u_\mathbf{k}(\mathbf{r}).$$
The function $u_\mathbf{k}$ is lattice-periodic. Bloch's theorem follows: $\psi_\mathbf{k}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}}u_\mathbf{k}(\mathbf{r})$.

In a crystal with lattice vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$, Bloch's theorem says

$$ \psi_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}} u_{n\mathbf{k}}(\mathbf{r}), \qquad u_{n\mathbf{k}}(\mathbf{r}+\mathbf{R}) = u_{n\mathbf{k}}(\mathbf{r}), $$

and the cell-periodic part $u_{n\mathbf{k}}$ has a natural Fourier expansion over reciprocal lattice vectors $\mathbf{G}$:

$$ u_{n\mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt{\Omega}} \sum_\mathbf{G} c_{n\mathbf{k},\mathbf{G}}\, e^{i\mathbf{G}\cdot\mathbf{r}}, $$

so

$$ \psi_{n\mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt{\Omega}} \sum_\mathbf{G} c_{n\mathbf{k},\mathbf{G}}\, e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}. $$

#### Matrix elements of $\hat{T}$ and $\hat{V}_\mathrm{ext}$

Plane waves are not merely convenient labels: every operator that appears in Kohn-Sham theory has a closed-form expression in this basis. Two examples make the point.

The **kinetic-energy operator** $\hat{T} = -\tfrac{1}{2}\nabla^2$ is diagonal. For two plane waves $|\mathbf{k}+\mathbf{G}\rangle = \Omega^{-1/2}e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}$,
$$
\langle\mathbf{k}+\mathbf{G}|\hat{T}|\mathbf{k}+\mathbf{G}'\rangle
= -\tfrac{1}{2}\int_\Omega \Omega^{-1}\,e^{-i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}\,\nabla^2 e^{i(\mathbf{k}+\mathbf{G}')\cdot\mathbf{r}}\,d^3r
= \tfrac{1}{2}|\mathbf{k}+\mathbf{G}'|^2\,\delta_{\mathbf{G},\mathbf{G}'}.
$$
This is the only completely diagonal piece of the Hamiltonian and the reason the plane-wave cutoff in Eq. (6.1) reads as a kinetic-energy bound.

The **external (ionic, lattice-periodic) potential** $\hat{V}_\mathrm{ext}(\mathbf{r})$ has Fourier components $V_\mathrm{ext}(\mathbf{G})$, and its plane-wave matrix elements are simply convolutions in reciprocal space:
$$
\langle\mathbf{k}+\mathbf{G}|\hat{V}_\mathrm{ext}|\mathbf{k}+\mathbf{G}'\rangle
= \frac{1}{\Omega}\int_\Omega e^{-i(\mathbf{G}-\mathbf{G}')\cdot\mathbf{r}}\,V_\mathrm{ext}(\mathbf{r})\,d^3r
= V_\mathrm{ext}(\mathbf{G}-\mathbf{G}').
$$
The Hartree term $\hat{V}_\mathrm{H}$ is similarly diagonal in the *density* Fourier components: $V_\mathrm{H}(\mathbf{G}) = 4\pi\rho(\mathbf{G})/|\mathbf{G}|^2$ for $\mathbf{G}\neq 0$. The exchange-correlation operator $\hat{V}_\mathrm{xc}$ is multiplicative in $\mathbf{r}$ and computed by FFT round-trip (real space $\to$ apply $V_\mathrm{xc}[\rho(\mathbf{r})]$ pointwise $\to$ back to $\mathbf{G}$-space). Every matrix element is either trivially diagonal or efficiently FFT-able. This algorithmic transparency is why plane-wave codes parallelise so cleanly.

#### Completeness and systematic convergence

The set $\{e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}/\sqrt{\Omega}\}_{\mathbf{G}\in\Lambda^*}$ is a **complete orthonormal basis** for any function on the torus $\mathbb{R}^3/\Lambda$. As $E_\mathrm{cut}\to\infty$ we retain every $\mathbf{G}$ and recover the full Hilbert space. The variational theorem then tells us:

1. The truncated $E^\mathrm{trunc}_\mathrm{tot}(E_\mathrm{cut})$ is a **monotonically decreasing** function of $E_\mathrm{cut}$.
2. $E^\mathrm{trunc}_\mathrm{tot}(E_\mathrm{cut}) \to E^\mathrm{exact}_\mathrm{tot}$ from above as $E_\mathrm{cut}\to\infty$.
3. Adding *any* plane wave with $\tfrac{1}{2}|\mathbf{k}+\mathbf{G}|^2 \le E_\mathrm{cut}'$ for $E_\mathrm{cut}' > E_\mathrm{cut}$ can only lower (or leave unchanged) the energy.

This **systematic convergence** is what people mean by "one knob". Compare to Gaussian basis sets, where moving from cc-pVDZ to cc-pVTZ is a discrete jump that can in principle raise the energy if the larger set lacks a function from the smaller set (it usually does not, but the guarantee is not automatic). Plane waves give you a single continuous lever you turn until the answer stops moving.

The basis is fully specified by one number: the **kinetic-energy cutoff** $E_\mathrm{cut}$. Only plane waves with

$$ \tfrac{1}{2}|\mathbf{k} + \mathbf{G}|^2 \leq E_\mathrm{cut} $$

are kept. (In Hartree atomic units; in Rydberg units, divide by 2 — Ry uses $\frac{1}{2}|\mathbf{k}+\mathbf{G}|^2$ in Hartree but quotes it in Ry. QE asks for `ecutwfc` in Ry, so do not double-convert.)

Plane waves are the natural basis for periodic solids. They have several decisive advantages:

- **One convergence parameter.** Increase $E_\mathrm{cut}$ monotonically; the answer improves monotonically (variational principle). No basis-set superposition error.
- **Orthonormal and unbiased.** Every region of the cell is treated equally; no preference for one atom over another.
- **Fast Fourier transforms** make the Hamiltonian-on-vector operation $O(N \log N)$ rather than $O(N^2)$.
- **Forces and stresses come almost for free** via the Hellmann-Feynman theorem, with no Pulay corrections (basis functions do not move with atoms because they are pinned to the cell).

Their disadvantages: they describe vacuum equally well as electron-rich regions, so isolated molecules in large boxes are wasteful; and core electrons would require absurd cutoffs, so plane waves *require* pseudopotentials.

### Gaussian-type orbitals (GTO)

Quantum chemistry's choice. Each atom carries a set of Gaussians $\chi_\alpha(\mathbf{r}) \propto (x-X_a)^l(y-Y_a)^m(z-Z_a)^n e^{-\zeta |\mathbf{r}-\mathbf{R}_a|^2}$. Compact for molecules; integrals over four centres have closed-form expressions. Used by Gaussian, NWChem, CP2K, FHI-aims. Drawback for solids: Pulay forces (basis moves with atoms), basis-set superposition error, and no monotonic convergence — you choose between cc-pVDZ, cc-pVTZ, cc-pVQZ in discrete steps.

### Numerical atomic orbitals (NAO)

Tabulated radial functions $R_{n\ell}(r)$ times spherical harmonics $Y_{\ell m}$, anchored on atoms. Used by SIESTA, FHI-aims, OpenMX. Very compact, scales linearly with system size for sparse-matrix algorithms, and excellent for large systems (thousands of atoms). Convergence is non-trivial: you pick a "minimal", "double-$\zeta$", "double-$\zeta$ polarised" basis from a catalogue.

### Why solid-state DFT prefers plane waves

For periodic, dense solids, plane waves give you variational convergence, no Pulay forces, and trivial parallelisation over $\mathbf{G}$-vectors. Almost every major solid-state code (QE, VASP, ABINIT, Castep) is plane-wave based. For surfaces with thick vacuum, large unit cells, or large systems, NAO-based codes such as FHI-aims become competitive or superior, but plane waves remain the default for first calculations and for any cross-code comparison.

## 6.1.3 The cutoff energy

The condition

$$ \frac{1}{2}|\mathbf{k} + \mathbf{G}|^2 \leq E_\mathrm{cut} \tag{6.1} $$

selects a finite sphere of $\mathbf{G}$-vectors. The number of plane waves in the basis is

$$ N_\mathrm{PW} \approx \frac{\Omega}{6\pi^2}\, (2 E_\mathrm{cut})^{3/2}, $$

where $\Omega$ is the cell volume in Bohr$^3$. For Si in its 2-atom primitive cell ($\Omega \approx 270\,a_0^3$) at $E_\mathrm{cut} = 40$ Ry, $N_\mathrm{PW} \approx 1100$. Doubling the cutoff multiplies the basis size by $2^{3/2} \approx 2.8$, and the cost of the FFT and the iterative diagonalisation scales similarly.

The right cutoff is whatever your pseudopotential needs. SSSP-PBE-efficiency lists, for each element, a recommended `ecutwfc` (wavefunction cutoff) and a recommended `ecutrho` (density cutoff, used because the density $\rho(\mathbf{r}) = \sum_n |\psi_n(\mathbf{r})|^2$ has Fourier components up to $2|\mathbf{G}|_\mathrm{max}$). For norm-conserving pseudopotentials, $E_\mathrm{cut}^\rho = 4 E_\mathrm{cut}^\psi$ is exact; for USPP and PAW you typically need $E_\mathrm{cut}^\rho = 8$-$12 \times E_\mathrm{cut}^\psi$ to handle augmentation charges.

!!! warning "ecutrho is not optional"
    A common mistake is to converge `ecutwfc` (e.g., 40 Ry) and leave `ecutrho` at its default (4 × ecutwfc, i.e. 160 Ry). For USPP/PAW this is too low — typical values are 320-480 Ry. Symptoms: spurious total-energy oscillations as you change the cell, bad stresses, charge sloshing. Always set `ecutrho` explicitly, following the pseudopotential's recommendation.

For our silicon example using SSSP-PBE-efficiency, the recommended values are `ecutwfc = 30` Ry and `ecutrho = 240` Ry. We will use 40 Ry / 320 Ry to be a little conservative, and verify convergence in [§6.3](03-convergence.md).

## 6.1.4 k-point sampling

!!! tip "Why k-points exist — the picture"
    A crystal extends to infinity. The Schrödinger equation in such a system has uncountably many solutions, indexed by a continuous label $\mathbf{k}$ in the Brillouin zone. To compute *anything* — the total energy, the density, a force — you must integrate over $\mathbf{k}$. The integral is infinite-dimensional in principle (one Hamiltonian per $\mathbf{k}$) but only $\sim$ tens of "characteristic" $\mathbf{k}$-values actually matter, because the Bloch eigenvalues $\epsilon_{n\mathbf{k}}$ vary smoothly with $\mathbf{k}$ (gradient ≈ group velocity). The job of "k-point sampling" is to pick a finite set of $\mathbf{k}$ that integrates this smooth function accurately.
    
    The Monkhorst-Pack idea: place $\mathbf{k}$'s on a uniform grid in the BZ, like sampling a smooth function on a Cartesian grid. The denser the grid, the better the integral. The clever choice of grid offset (cell centres, not corners) makes the integration error die fast in grid spacing — it is exactly the same idea as the midpoint rule in 1D quadrature, dressed up for 3D.

The total energy of a crystal is an integral over the Brillouin zone (BZ):

$$ E_\mathrm{tot} = \sum_n \int_\mathrm{BZ} \frac{d^3\mathbf{k}}{(2\pi)^3} f_{n\mathbf{k}}\, \epsilon_{n\mathbf{k}} + \text{(double-counting)}. $$

We can only afford to evaluate the Hamiltonian at a discrete set of $\mathbf{k}$-points and replace the integral by a weighted sum:

$$ \int_\mathrm{BZ} \frac{d^3\mathbf{k}}{(2\pi)^3}\, F(\mathbf{k}) \to \sum_i w_i\, F(\mathbf{k}_i). $$

The question is how to choose the $\mathbf{k}_i$ and $w_i$ to converge the integral fastest.

### Monkhorst-Pack grids

Monkhorst and Pack (1976) gave the standard answer for periodic functions: a uniform grid

$$ \mathbf{k}_{n_1 n_2 n_3} = \sum_{\alpha=1}^3 \frac{2 n_\alpha - N_\alpha - 1}{2 N_\alpha}\, \mathbf{b}_\alpha, \qquad n_\alpha = 1, \ldots, N_\alpha, $$

where $\mathbf{b}_\alpha$ are the reciprocal lattice vectors. We specify the grid by three integers $N_1 \times N_2 \times N_3$. In QE this is `K_POINTS automatic` followed by `N1 N2 N3 s1 s2 s3`, where `s_alpha = 0` or `1` shifts the grid by half a step.

!!! note "Why this formula? — Deriving the MP special points"
    The starting observation is that any lattice-periodic function on the BZ can be expanded in a Fourier series in the *direct* lattice vectors $\mathbf{R}$: $F(\mathbf{k}) = \sum_\mathbf{R} F_\mathbf{R}\,e^{i\mathbf{k}\cdot\mathbf{R}}$. The BZ-average is $\langle F\rangle = F_{\mathbf{0}}$. To estimate $F_\mathbf{0}$ from a finite sample, one picks a grid that makes the "noise" contributions $\sum_{\mathbf{R}\neq\mathbf{0}} F_\mathbf{R}\,e^{i\mathbf{k}_n\cdot\mathbf{R}}$ vanish for as many short $\mathbf{R}$ as possible.
    
    Consider one dimension first. The choice $k_n = (2n-N-1)/(2N)\cdot b$ for $n=1,\ldots,N$ places the grid points at the *midpoints* of $N$ equal sub-intervals of $(0,b)$, equivalently at $k_n/b = (n-\tfrac{1}{2})/N - \tfrac{1}{2}$. The discrete sum $\sum_n e^{i k_n R} = \sum_n e^{i(2n-N-1)\pi R / N}$ is a geometric series; it equals $N$ if $R$ is an integer multiple of $N$ (times $a$) and zero otherwise. So all Fourier components with $|R| < Na$ are killed exactly; the smallest surviving "leakage" comes from $|R|=Na$, controlling convergence.
    
    Generalising to three dimensions and shifting the origin to half-spacing as in the formula above yields the Monkhorst-Pack grid. The reason the points sit at the *centres* (not corners) of cells is that the centre placement avoids putting weight on high-symmetry zone-boundary points like $L$ or $X$, which lie on a measure-zero set in the integral and would have to be weighted by hand.

The general MP formula above reduces, when written in the equivalent shifted form
$$\mathbf{k}_{n_1n_2n_3} = \sum_\alpha \frac{n_\alpha - \tfrac{1}{2}(N_\alpha+1) + \tfrac{1}{2}s_\alpha}{N_\alpha}\,\mathbf{b}_\alpha,$$
to two conventional choices: unshifted ($s_\alpha=0$, gives the original Monkhorst-Pack offset) and $\Gamma$-shifted ($s_\alpha=1$, moves the grid by half a step so that $\Gamma$ falls on the grid when $N_\alpha$ is even).

!!! example "8-atom Si — counting irreducible k-points"
    A $4\times 4\times 4$ unshifted MP grid on the FCC primitive lattice generates 64 k-points. The point group $O_h$ has 48 operations; folding by the symmetry, only **10 irreducible points** remain. With shift $(1,1,1)$ that number is also 10 — the same density, different points. With a smaller $2\times 2\times 2$ unshifted grid the count is just **3 irreducible**: $\Gamma$, $L$, and the midpoint $W$. This is the count `pw.x` prints near the top of the output.

### $\Gamma$-centred vs shifted

The simplest choice (no shift) places one $\mathbf{k}$-point at $\Gamma = (0,0,0)$ when $N_\alpha$ is odd, and avoids $\Gamma$ when $N_\alpha$ is even. Putting the shift `s = 1 1 1` moves a point to $\Gamma$ when $N$ is even.

- **Use $\Gamma$-centred** (no shift, or shift such that $\Gamma$ is included) for hexagonal lattices (the standard MP rule is wrong for hexagonal symmetry — see Pack and Monkhorst 1977 erratum); for cells with low symmetry; whenever a band crossing sits at $\Gamma$; and for any calculation where you want consistency between scf and band-path runs.
- **Use shifted grids** (`s = 1 1 1`) for cubic systems with $N$ even to get a denser effective sampling: a shifted $4\times4\times4$ grid samples 8 inequivalent points versus 8 for the unshifted version, but the shifted points avoid the high-symmetry corners and converge integrals over smooth quantities faster.

For silicon (FCC, cubic), a shifted $8\times8\times8$ or unshifted $4\times4\times4$ is standard. We use unshifted $4\times4\times4$ in [§6.2](02-first-qe.md) for simplicity, then converge it in [§6.3](03-convergence.md).

### How dense is dense enough

For insulators, the integrand is smooth (the occupations $f_{n\mathbf{k}}$ are either 0 or 1 everywhere, with a gap between), and even coarse grids ($4\times4\times4$ for Si's 2-atom cell) suffice. For metals, the occupations have a step at the Fermi surface that requires either dense sampling or smearing (see [§6.3](03-convergence.md)); typical metallic k-grids are $12\times12\times12$ or denser per primitive cell.

Rule of thumb: $N_i \cdot |\mathbf{a}_i| \gtrsim 30$ Å is a reasonable starting density. A 10 Å unit cell wants $N_i = 3$ or 4; a 30 Å supercell wants $N_i = 1$ (just $\Gamma$).

!!! warning "k-grids and supercells"
    When you make a $2\times2\times2$ supercell of a primitive cell, the BZ shrinks by 8 in volume. A k-grid that was $8\times8\times8$ for the primitive cell becomes $4\times4\times4$ for the supercell — the *density* of k-points in absolute reciprocal-space units is preserved, not the *number*. Forgetting this rule is the most common bug in supercell calculations.

## 6.1.4b Smearing functions for metals — derivations

For metals the integrand in the BZ sum has a step at the Fermi level $\mu$ and the integral converges as $1/N_k$ — almost useless. The cure is to replace the sharp Heaviside occupation $\theta(\mu-\epsilon)$ by a smooth function $f((\epsilon-\mu)/\sigma)$ of width $\sigma$. Three choices dominate; each is best understood as a derivation from a different approximation strategy.

### Gaussian smearing

The simplest smooth replacement for the step is its derivative-of-a-Gaussian primitive:
$$f_\mathrm{G}(x) = \tfrac{1}{2}\,\mathrm{erfc}(x), \qquad \frac{df_\mathrm{G}}{dx} = -\frac{1}{\sqrt{\pi}}\,e^{-x^2},$$
with $x = (\epsilon-\mu)/\sigma$. The "free-energy" entropy contribution per state is
$$S_\mathrm{G}(x) = -\tfrac{1}{2}\sigma^{-1}\,e^{-x^2}/\sqrt{\pi},$$
and the corrected total energy at $T=0$ is $E_\mathrm{tot} = E^\mathrm{KS} - TS \approx E^\mathrm{KS}_\mathrm{Mermin}$. Convergence is *linear in $\sigma^2$*: doubling $\sigma$ quadruples the smearing error. Use only if MP or MV is unavailable.

### Methfessel-Paxton

Methfessel and Paxton (1989) noticed that you can construct a smoothing function as a finite-order Hermite expansion such that all moments up to order $2N+1$ of the original step function are reproduced. The first-order MP function is
$$f_\mathrm{MP}^{(1)}(x) = \tfrac{1}{2}\,\mathrm{erfc}(x) + \tfrac{1}{2\sqrt{\pi}}\,x\,e^{-x^2}.$$
The polynomial correction subtracts the leading error of Gaussian smearing. The result is that integrals of smooth functions over the smeared Fermi sea differ from their $\sigma\to 0$ limit by $O(\sigma^4)$ rather than $O(\sigma^2)$ — converges *four times* faster in $\sigma$.

### Marzari-Vanderbilt cold smearing

Marzari and Vanderbilt (1999) noted that MP can produce *negative* occupations $f^{(1)}_\mathrm{MP} < 0$ for $|x|\gtrsim 1$, which is unphysical for some downstream uses (orbital-dependent functionals, Wannier construction). They proposed an alternative
$$f_\mathrm{MV}(x) = \tfrac{1}{2}\,\mathrm{erfc}\!\left(x + \tfrac{1}{\sqrt{2}}\right) + \tfrac{1}{\sqrt{2\pi}}\,e^{-(x+1/\sqrt{2})^2}$$
based on a *non-symmetric* combination that remains in $[0,1]$ and still produces $O(\sigma^4)$ convergence. The entropy contribution is small even at moderately large $\sigma$ — hence "cold smearing" — making MV the default for production metallic DFT in QE today.

### A practical comparison

| Scheme | Order of error | Stays in $[0,1]$? | QE keyword |
|---|---|---|---|
| Fermi-Dirac | $\sigma^2$ | Yes | `smearing='fd'` |
| Gaussian | $\sigma^2$ | Yes | `smearing='gauss'` |
| Methfessel-Paxton 1 | $\sigma^4$ | No (can be slightly negative) | `smearing='mp'` |
| Marzari-Vanderbilt | $\sigma^4$ | Yes | `smearing='mv'` |

For everything in this book that involves a metal, we use MV with $\sigma = 0.01$-$0.02$ Ry; the reasoning follows the convergence study of [§6.3](03-convergence.md).

## 6.1.5 Pseudopotential libraries

Three curated libraries are the practical choices for plane-wave DFT.

### SSSP — Standard Solid-State Pseudopotentials

Hosted by [Materials Cloud](https://www.materialscloud.org/discover/sssp). For each element they pick the best-performing pseudopotential (from several libraries) by benchmarking against all-electron references on a curated test set (delta-factor, phonon frequencies, stress). They publish two libraries:

- **SSSP-PBE-efficiency**: chosen for low cutoff energies; typical `ecutwfc` 30-50 Ry; ideal for high-throughput screening and first calculations.
- **SSSP-PBE-precision**: chosen for accuracy; typical `ecutwfc` 50-80 Ry; use when you need converged equations of state or vibrational frequencies.

Both come with metadata: recommended cutoffs, the original library each pseudopotential is from, and validation data. **For this chapter we use SSSP-PBE-efficiency 1.3.0.** Download from [materialscloud.org](https://www.materialscloud.org/discover/sssp).

### PseudoDojo

[pseudo-dojo.org](http://www.pseudo-dojo.org). A library of norm-conserving (ONCV) pseudopotentials with extensive testing. Multiple accuracy tiers (standard / stringent) and multiple XC functionals (PBE, PBEsol, LDA, hybrid-ready). Use when you need NC for hybrids, GW, or DFPT phonons, or for cross-checking SSSP results.

### GBRV

The Garrity-Bennett-Rabe-Vanderbilt library of USPP. Very low cutoffs (around 25-30 Ry), well-tested for high-throughput, slightly older. SSSP draws on GBRV for many elements.

### Recommended starter setup

For everything in this chapter:

| Item | Choice |
|---|---|
| Library | SSSP-PBE-efficiency 1.3.0 |
| Pseudopotential format | UPF (the QE-native format) |
| For Si | `Si.pbe-n-rrkjus_psl.1.0.0.UPF` (the file SSSP-efficiency points to) |
| Cutoff | `ecutwfc = 40` Ry, `ecutrho = 320` Ry |
| k-grid (2-atom primitive Si) | $4\times4\times4$ unshifted (starter) — converge in [§6.3](03-convergence.md) |

Set a single environment variable so QE finds the files:

```bash
export ESPRESSO_PSEUDO=$HOME/pseudo/SSSP_1.3.0_PBE_efficiency
```

and `pw.x` will look there for any pseudopotential named in your input. We will use this in [§6.2](02-first-qe.md).

### Section summary

A solid-state DFT calculation needs three numerical ingredients:

1. **A pseudopotential** for each species — a smooth replacement for the bare nucleus-plus-core potential, justified by the frozen-core approximation. Norm-conserving, ultrasoft, or PAW; the choice trades off cutoff energy against complexity of the basis-set bookkeeping. PAW is the modern default.
2. **A plane-wave basis** characterised by the cutoff $E_\mathrm{cut}$: keep all $\mathbf{G}$ with $\tfrac{1}{2}|\mathbf{k}+\mathbf{G}|^2 \le E_\mathrm{cut}$. The basis is complete in the limit $E_\mathrm{cut}\to\infty$ and converges systematically. The density cutoff `ecutrho` is a separate parameter for USPP/PAW.
3. **A Monkhorst-Pack k-grid** $N_1\times N_2\times N_3$ — the discrete sampling of the BZ integral. Density chosen for the BZ-integration tolerance you need; metals need denser grids than insulators because of Fermi-surface steps, which require smearing functions (MV cold smearing is the production default).

The combination $(E_\mathrm{cut}, N_k, \sigma)$ is what you converge in §6.3.

## 6.1.6 Three flavours of pseudopotential — the detailed comparison

With a calculation now running, we can return to the question deferred from §6.1.1: *which* pseudopotential, and why. A curated library has already made a sound choice for you, but understanding the trade-offs is what lets you read pseudopotential metadata, diagnose convergence trouble, and pick the right family for properties beyond total energy. Three families dominate modern DFT codes. They differ in how aggressively they smooth the pseudo-wavefunction, and therefore in what cutoff energy you need.

### Norm-conserving (NC)

Introduced by Hamann, Schlüter and Chiang in 1979. The pseudo-wavefunction $\tilde\psi_{n\ell}$ is required to have the same norm inside $r_c$ as the all-electron wavefunction $\psi_{n\ell}^\mathrm{AE}$:

$$ \int_0^{r_c} |\tilde\psi_{n\ell}(r)|^2 r^2\, dr = \int_0^{r_c} |\psi_{n\ell}^\mathrm{AE}(r)|^2 r^2\, dr. $$

This norm condition guarantees correct scattering at the reference energy and, crucially, correct *energy derivative* of the scattering — the pseudopotential is transferable across chemical environments. Modern "ONCV" (optimised norm-conserving Vanderbilt) pseudopotentials are tight, accurate, and the gold standard for high-precision work. They require cutoffs of 50-100 Ry for typical elements.

### Ultrasoft (USPP)

Vanderbilt (1990) relaxed the norm condition. The pseudo-wavefunction can be much smoother — *ultrasoft* — at the price that the orthonormality condition changes from $\langle \tilde\psi_i | \tilde\psi_j \rangle = \delta_{ij}$ to a generalised form $\langle \tilde\psi_i | S | \tilde\psi_j \rangle = \delta_{ij}$ with an overlap operator $S$ that depends on augmentation charges localised near each atom. The eigenvalue problem becomes generalised. USPP cutoffs are typically 25-40 Ry — a factor of 2-4 cheaper than NC.

### Projector augmented wave (PAW)

Blöchl (1994). A linear transformation $|\psi\rangle = |\tilde\psi\rangle + \sum_i (|\phi_i\rangle - |\tilde\phi_i\rangle)\langle \tilde p_i | \tilde\psi\rangle$ takes a smooth pseudo-wavefunction $|\tilde\psi\rangle$ to the true all-electron wavefunction $|\psi\rangle$ by adding back, atom by atom, the difference between all-electron partial waves $\phi_i$ and pseudo partial waves $\tilde\phi_i$. PAW reproduces all-electron results to within a few meV/atom for most properties, is as cheap as USPP, and gives you access to the true wavefunction near the nucleus (useful for hyperfine fields, NMR, EFGs). It is the default in VASP and GPAW, and widely available in QE.

#### Valence wavefunction reconstruction in one picture

A schematic comparison helps. Inside the core radius $r_c$, the all-electron valence orbital $\psi^\mathrm{AE}_{n\ell}(r)$ oscillates rapidly because it must remain orthogonal to the (likewise rapidly oscillating) core orbitals; outside, $\psi^\mathrm{AE}$ is smooth. Norm-conserving, ultrasoft, and PAW pseudopotentials all introduce a smooth pseudo-orbital $\tilde\psi_{n\ell}$ that *equals* $\psi^\mathrm{AE}$ for $r>r_c$ and *replaces* the oscillations inside. The flavours differ in **what they keep track of**:

- **NC**: only $\tilde\psi$ is stored; the all-electron $\psi^\mathrm{AE}$ is forgotten, but the norm condition guarantees correct scattering. Orthonormality $\langle\tilde\psi_i|\tilde\psi_j\rangle=\delta_{ij}$ is preserved.
- **USPP**: $\tilde\psi$ is *more* aggressively smoothed (no norm conservation); a localised augmentation charge $Q_{ij}(\mathbf{r})$ is added to the density to restore the correct total charge.
- **PAW**: $\tilde\psi$ is stored *and* the transformation back to $\psi^\mathrm{AE}$ is stored. PAW thus retains, formally, all the information of an all-electron calculation; properties involving the wavefunction near the nucleus (NMR shielding, hyperfine fields, EFGs) are accessible.

#### The augmentation operator $\hat{S}$

For USPP and PAW the orthogonality of the all-electron orbitals translates into a *generalised* orthonormality for the pseudo-orbitals:
$$\langle\tilde\psi_i|\hat{S}|\tilde\psi_j\rangle = \delta_{ij}, \qquad \hat{S} = \mathbb{1} + \sum_{a,ij} q^{(a)}_{ij}\,|\tilde{p}^{(a)}_i\rangle\langle\tilde{p}^{(a)}_j|.$$
Here the sum is over atoms $a$ and projector indices $i,j$, and $q^{(a)}_{ij} = \int_{r<r_c}[\phi^{(a)}_i(\mathbf{r})\phi^{(a)}_j(\mathbf{r}) - \tilde\phi^{(a)}_i(\mathbf{r})\tilde\phi^{(a)}_j(\mathbf{r})]\,d^3r$ are augmentation integrals. The Kohn-Sham equation becomes a generalised eigenvalue problem
$$\hat{H}|\tilde\psi_n\rangle = \epsilon_n\,\hat{S}|\tilde\psi_n\rangle,$$
which costs roughly the same as the standard problem to solve. Forces inherit a Pulay-like correction from the implicit $\mathbf{r}$-dependence of $\hat{S}$ through the projectors, which production codes handle automatically; you just have to know it is there when you read the QE output and see a term called "augmentation".

#### Pseudopotential construction recipe

The half-page algorithm used in every modern NC/PAW generator (Hamann, Vanderbilt, Blöchl) is:

1. **All-electron atomic reference.** Solve the radial Kohn-Sham equation for the isolated atom; obtain $\psi^\mathrm{AE}_{n\ell}(r)$ and eigenvalues $\epsilon^\mathrm{AE}_{n\ell}$ for each valence channel.
2. **Choose a core radius $r_c^{(\ell)}$.** Typically $r_c \approx 1.0$-$2.0\,a_0$, large enough that $\tilde\psi$ can be smooth but small enough that valence chemistry happens entirely outside.
3. **Construct a smooth $\tilde\psi$.** For $r > r_c$, set $\tilde\psi = \psi^\mathrm{AE}$. For $r < r_c$, replace by a polynomial (or spherical Bessel sum, or whatever the scheme prefers) that is nodeless and matches $\psi^\mathrm{AE}$ in value and derivative at $r_c$.
4. **Invert the radial Schrödinger equation.** Given $\tilde\psi$ and the reference energy $\epsilon^\mathrm{AE}_{n\ell}$, the radial Schrödinger equation
   $$\left[-\tfrac{1}{2r^2}\partial_r(r^2\partial_r) + \tfrac{\ell(\ell+1)}{2r^2} + V^{(\ell)}_\mathrm{ps}(r)\right]\tilde\psi(r) = \epsilon^\mathrm{AE}_{n\ell}\,\tilde\psi(r)$$
   is solved for the *channel-dependent* potential $V^{(\ell)}_\mathrm{ps}(r)$. Different $\ell$'s see different pseudopotentials: the pseudopotential is **non-local**.
5. **Log-derivative test.** The accuracy of a pseudopotential is the agreement between $\psi^\mathrm{AE}$ and $\tilde\psi$ scattering log-derivatives,
   $$L^{(\ell)}(\epsilon, r_c) \equiv \frac{r\,\partial_r \psi(r,\epsilon)}{\psi(r,\epsilon)}\bigg|_{r=r_c},$$
   over a range of energies $\epsilon$. The norm condition (Hamann-Schlüter-Chiang) guarantees that the first energy derivative $\partial_\epsilon L^{(\ell)}$ agrees at the reference energy; ONCV-style optimisation also matches higher derivatives.
6. **Unscreen.** Subtract the Hartree and XC contributions from the *atomic* valence density to leave a pure ionic pseudopotential, transferable to any chemical environment.

This recipe is implemented in `ld1.x` (QE), `oncvpsp` (Hamann's code, the engine behind PseudoDojo), and `atompaw`. Modern users download finished UPF/PSP8 files; understanding the recipe explains why the "PBE pseudopotential" is not unique and why curated libraries matter.

### Which to use

For a starter calculation: **PAW or efficient USPP via SSSP-PBE-efficiency** (see §6.1.5). Cutoffs around 40-50 Ry, runs in seconds for small cells. For high-precision properties (especially absolute energies and pressures), prefer **norm-conserving** (PseudoDojo) with cutoffs around 80 Ry. For NMR shielding tensors, hyperfine, or anything that probes the wavefunction near the nucleus, use **PAW**.

!!! warning "Mix-and-match pseudopotentials at your peril"
    A pseudopotential is generated together with a particular exchange-correlation functional. A PBE pseudopotential is not interchangeable with an LDA one, and even between PBE pseudopotentials from different libraries the all-electron reference may differ in subtle ways (relativistic treatment, frozen-core boundary, choice of valence configuration). Within one calculation, use pseudopotentials from a *single, consistent set*. Mixing SSSP for one element and PseudoDojo for another is a common bug.

!!! tip "Why pseudopotentials exist — the historical moment"
    Before the 1970s, all-electron DFT of solids was the only option, and it was painfully expensive: a $1s$ orbital of Si oscillates at the Bohr-radius scale ($\sim 0.05$ Å) while bonding happens at the bond-length scale ($\sim 2.4$ Å). Resolving both requires a basis 50× finer than needed for chemistry. The pseudopotential idea — already implicit in the orthogonalised-plane-wave method of Herring (1940) and the empirical pseudopotential method of Phillips and Kleinman (1959) — is to *replace* the core's oscillation by a smooth function whose scattering of valence electrons is unchanged. Hamann, Schlüter, and Chiang's 1979 paper on *ab initio* norm-conserving pseudopotentials was the breakthrough: a recipe for constructing the smooth replacement directly from atomic DFT, with provable accuracy. Modern PAW and ONCV are descendants. The "aha" is that *chemistry only sees the valence shell*; the core is geometrically inert. Once you accept that, the smooth replacement is the obvious move.

## 6.1.7 Summary

You now have the knobs you need:

- **Pseudopotential** — frozen-core, smooth representation of nuclear plus core potential. Choose NC, USPP, or PAW based on accuracy/cost trade-off. SSSP-PBE-efficiency is a safe default.
- **Plane-wave basis** — single number $E_\mathrm{cut}$ controls completeness; chosen large enough for the pseudopotential to be accurate.
- **k-point grid** — Monkhorst-Pack $N_1 \times N_2 \times N_3$, $\Gamma$-centred or shifted; chosen dense enough for the BZ integration tolerance you need.

**Remark 6.1.5 (The cubic root in basis count).** The plane-wave basis size $N_\mathrm{PW} \propto \Omega\,E_\mathrm{cut}^{3/2}$ grows as the cell volume times the 3/2 power of the cutoff. This explains a recurring observation: doubling the supercell linear size (×8 in volume) makes a calculation 8× harder; doubling the cutoff (×2.8 in $N_\mathrm{PW}$) makes it ~3× harder per k-point. Combined, an 8-atom supercell at 50 Ry is roughly 25× more expensive than a 2-atom primitive cell at 30 Ry — exactly what you see in wall-time benchmarks.

In the next section we put numbers into these knobs and run a calculation.
