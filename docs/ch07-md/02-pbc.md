# 7.2 Periodic Boundary Conditions

!!! tip "Why periodic boundary conditions exist"
    A real crystal has $\sim 10^{23}$ atoms; an MD simulation can manage $\sim 10^3$ to $10^8$. The atoms at the surface of a small simulation box are not surrounded by their proper environment — they see vacuum. For a box of $10^3$ argon atoms, more than half of them are within two atomic diameters of a face. The "simulation" then describes an unphysically small cluster, not a bulk material.
    
    The cure is a kind of mathematical magic: identify opposite faces of the box, turning the cube into a topological 3-torus. An atom that walks off the right face re-enters from the left at the same height. Periodic copies of the cell tile space; each atom is surrounded by an infinite homogeneous environment.
    
    Picture in your head: a video game world that wraps. You walk east and arrive at your starting point from the west. Physically, PBC are how we make a small box pretend to be a chunk of bulk material. The price: artefacts when interactions are longer than the box (the cause is that one atom can interact with multiple periodic copies of another), and missing physics for interfaces (you cannot have a free surface in a torus).

## Why we need PBC

Take 1000 argon atoms in a cubic box of side 30 Å. If the box has hard walls, the atoms within one or two atomic diameters of a wall experience an environment qualitatively different from those in the bulk: fewer neighbours, asymmetric forces, anomalous density. Of our 1000 atoms, perhaps 600 are "surface" atoms in this sense. The simulation is then a study of an extremely small cluster, not of bulk argon.

This is a problem for materials science, where we typically want bulk properties of a notional infinite system. The cure is to identify opposite faces of the simulation cell, turning a cube into a 3-torus. An atom that exits the right face re-enters from the left at the same height and depth; its periodic images tile space. Every atom now lives in an infinite homogeneous medium, and surface artefacts vanish.

Formally: given a lattice with primitive vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$, the potential satisfies

$$
U(\mathbf{r}_1, \ldots, \mathbf{r}_i + \mathbf{n}, \ldots, \mathbf{r}_N) = U(\mathbf{r}_1, \ldots, \mathbf{r}_i, \ldots, \mathbf{r}_N),
\qquad
\mathbf{n} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + n_3 \mathbf{a}_3,
\quad n_k \in \mathbb{Z}.
\tag{7.14}
$$

The simulation cell holds $N$ atoms; the infinite lattice holds their images. We integrate Newton's equations only for the $N$ in-cell atoms; image positions are determined by the lattice translation.

**Definition 7.2.1 (Periodic system).** A system of $N$ atoms in a simulation cell with lattice matrix $\mathbf{H} = (\mathbf{a}_1\,|\,\mathbf{a}_2\,|\,\mathbf{a}_3)$ is *periodic* if the potential satisfies (7.14) for all integer triples $(n_1,n_2,n_3)\in\mathbb{Z}^3$.

**Definition 7.2.2 (Minimum image).** Given two atoms at positions $\mathbf{r}_i, \mathbf{r}_j$ in a periodic cell, the *minimum image* of $j$ relative to $i$ is the periodic image $\mathbf{r}_j + \mathbf{n}\cdot\mathbf{H}$ minimising $|\mathbf{r}_j + \mathbf{n}\cdot\mathbf{H} - \mathbf{r}_i|$. For an orthorhombic cell with sides $L_\alpha$ and cutoff $r_c < L_\alpha/2$, the minimum image is unique.

!!! abstract "Key idea (Half-box rule)"
    Periodic boundary conditions and a finite-range force coexist *only* if the force range is strictly less than half the cell's smallest perpendicular width. Beyond that limit, multiple periodic images come into range and the minimum-image convention breaks. This is the single most important geometric constraint in MD.

## The minimum image convention

The pair potential between atom $i$ at $\mathbf{r}_i$ and atom $j$ at $\mathbf{r}_j$ is now ambiguous: which image of $j$ do we use? For a short-ranged potential whose range is less than half the smallest box dimension, the answer is unique — only one image of $j$ lies within range of $i$. This is the **minimum image convention**.

For a cubic cell of side $L$ with $r_\mathrm{cut} < L/2$, the minimum image displacement is

$$
\Delta x_{ij} = x_j - x_i - L \cdot \mathrm{round}\!\left(\frac{x_j - x_i}{L}\right),
\tag{7.15}
$$

and similarly for $y$ and $z$. The pair distance is then $r_{ij} = \|\Delta \mathbf{r}_{ij}\|$. In Python:

```python
import numpy as np

def minimum_image(dr: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Apply minimum image convention to displacements.

    Parameters
    ----------
    dr : (N, 3) or (3,) array of displacement vectors.
    box : (3,) array of orthogonal box lengths.

    Returns
    -------
    Wrapped displacement(s) in the central image.
    """
    return dr - box * np.round(dr / box)
```

For a non-orthogonal (triclinic) cell with lattice matrix $\mathbf{H} = (\mathbf{a}_1\;\mathbf{a}_2\;\mathbf{a}_3)$, the convention generalises to fractional coordinates $\mathbf{s} = \mathbf{H}^{-1}\mathbf{r}$:

```python
def minimum_image_triclinic(dr: np.ndarray, H: np.ndarray) -> np.ndarray:
    """Minimum image for a triclinic cell with column-vector lattice matrix H."""
    s = np.linalg.solve(H, dr.T).T          # fractional
    s -= np.round(s)                         # wrap into (-0.5, 0.5]
    return (H @ s.T).T                       # back to Cartesian
```

This works whenever the cell is sufficiently "boxy" — formally when the cutoff is less than half the perpendicular width of the cell, which for triclinic geometries is more restrictive than half the longest edge. LAMMPS' `neighbor` command warns you when this is violated.

!!! warning "Cell sizes that are too small"
    A frequent mistake is to set up a 5x5x5 unit cell supercell, decide the cutoff is 12 Å "because the literature says so", and then discover the cell is only 20 Å across. The minimum image convention breaks and you are now simulating an unphysical periodic structure. Either grow the supercell or shrink the cutoff.

## Wrapped versus unwrapped coordinates

This distinction trips up nearly every new MD user.

- **Wrapped** coordinates lie inside the primary simulation cell: $0 \le x < L$ (or $-L/2 \le x < L/2$, depending on convention). They are what visualisation tools want.
- **Unwrapped** coordinates ignore the wrapping. As an atom diffuses through the cell, its $x$ coordinate grows monotonically beyond $L$; jumps are smooth.

For computing forces, wrapping is fine — only displacements within the minimum image matter. For computing **mean squared displacement** (§7.6), wrapping is fatal: the MSD will saturate at a value of order $L^2$ rather than growing linearly with time. Always compute MSD from unwrapped coordinates.

The relationship is

$$
\mathbf{r}^\mathrm{unwrapped}_i(t) = \mathbf{r}^\mathrm{wrapped}_i(t) + \mathbf{n}_i(t) \cdot \mathbf{H},
\tag{7.16}
$$

where $\mathbf{n}_i(t) \in \mathbb{Z}^3$ counts how many times atom $i$ has crossed each face. LAMMPS' `dump` command stores the image flags $\mathbf{n}_i$ alongside the wrapped positions in standard output. ASE's `Atoms.get_positions(wrap=False)` returns the wrapped positions; you have to reconstruct unwrapped trajectories by tracking jumps:

```python
def unwrap(positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Unwrap a trajectory (T, N, 3) wrapped in a box (3,)."""
    unwrapped = positions.copy()
    for t in range(1, len(positions)):
        dr = positions[t] - positions[t - 1]
        jumps = np.round(dr / box)
        unwrapped[t:] -= jumps * box  # propagate correction to all future frames
    return unwrapped
```

A subtler issue arises with molecular systems. A water molecule straddling the cell boundary may have its oxygen on one side and its hydrogens on the other. The geometry is still correct under minimum-image (the O–H bond is short), but if you forget and compute the centre of mass with naive averaging, it will land in the middle of the cell — nowhere near either atom. Molecular trajectories should be "unwrapped within molecule, wrapped between molecules" for visualisation.

## Cutoffs and shifted potentials

For short-ranged pair potentials we never sum over all image pairs; we truncate at a cutoff $r_\mathrm{cut}$. The potential becomes

$$
U_\mathrm{trunc}(r) = \begin{cases} U(r) & r < r_\mathrm{cut} \\ 0 & r \ge r_\mathrm{cut}. \end{cases}
\tag{7.17}
$$

This introduces a discontinuity at $r_\mathrm{cut}$: a delta-function force as a pair crosses the cutoff. The discontinuity is small in magnitude (LJ at $r_\mathrm{cut} = 2.5\sigma$ has $U \approx -0.016\,\epsilon$) but feeds the integrator a non-conservative perturbation, producing a slow energy drift in NVE simulations.

The standard fix is to **shift** the potential so it vanishes at the cutoff:

$$
U_\mathrm{shift}(r) = U(r) - U(r_\mathrm{cut}),
\qquad r < r_\mathrm{cut}.
\tag{7.18}
$$

The force is unchanged (the shift is a constant) so equilibrium properties are unaffected, and there is no jump in $U$. Forces still have a small discontinuity at $r_\mathrm{cut}$. To kill that too, one can shift the force as well:

$$
F_\mathrm{shift}(r) = F(r) - F(r_\mathrm{cut}),
$$

but this changes the equation of state slightly. LAMMPS supports both shifting modes via the `pair_modify shift yes` and `pair_style lj/cut/smooth` family.

For the Lennard-Jones potential, the standard choice $r_\mathrm{cut} = 2.5\sigma$ truncates at about $-0.016\,\epsilon$. Production calculations of bulk thermodynamic properties either go further (typically $5\sigma$) or apply analytic **tail corrections** that integrate the potential and pressure contributions of pairs beyond the cutoff assuming a uniform $g(r) = 1$:

$$
U_\mathrm{tail} = 2\pi N \rho \int_{r_\mathrm{cut}}^\infty U(r)\, r^2\, dr,
\tag{7.19}
$$

with an analogous expression for the pressure. Tail corrections recover the long-range thermodynamics at zero additional cost; they fail near interfaces or whenever the assumption of uniform density breaks.

## A worked minimum-image example

Take a 2D snapshot for clarity. Two atoms at $\mathbf{r}_1 = (0.1, 0.1)$ and $\mathbf{r}_2 = (0.9, 0.9)$ in a unit-side square cell. The naive displacement is $\mathbf{r}_2 - \mathbf{r}_1 = (0.8, 0.8)$, with magnitude $\sqrt{1.28}\approx 1.13$. But $r_\mathrm{cut} < L/2 = 0.5$, so this distance is meaningless — there is a periodic image of atom 2 closer to atom 1.

Apply (7.15) with $L=1$:
$$\Delta x = 0.8 - 1\cdot \mathrm{round}(0.8) = 0.8 - 1 = -0.2.$$
Similarly $\Delta y = -0.2$. The minimum-image displacement is $(-0.2, -0.2)$, magnitude $\sqrt{0.08}\approx 0.28$. This is the image of atom 2 sitting at $(-0.1, -0.1)$ — across the boundary in both directions.

A unit test of the minimum-image function should include this 2D case: it catches off-by-one and sign errors that pass silently in symmetric 3D tests.

```python
def test_minimum_image() -> None:
    box = np.array([1.0, 1.0, 1.0])
    # Cross-corner pair
    dr = np.array([0.8, 0.8, 0.0])
    np.testing.assert_allclose(minimum_image(dr, box),
                                np.array([-0.2, -0.2, 0.0]))
    # Exactly half-box (ambiguous; numpy's round-to-even kicks in)
    dr = np.array([0.5, 0.0, 0.0])
    result = minimum_image(dr, box)
    assert abs(result[0]) < 1e-10 or abs(abs(result[0]) - 0.5) < 1e-10
    # Pair within range
    dr = np.array([0.1, 0.0, 0.0])
    np.testing.assert_allclose(minimum_image(dr, box), dr)
```

The "exactly half-box" case is an edge case: the closer image is mathematically tied. Most implementations let `numpy.round` (banker's rounding) decide; some break ties by always rounding down. For simulations near the half-box limit the choice does not matter because $r_\mathrm{cut}$ is strictly less than $L/2$ — atoms at the half-box separation are *outside* the cutoff and contribute zero force in either case.

### The half-box rule, formally

The minimum-image convention is well-defined as long as $r_\mathrm{cut} < L_\mathrm{perp}/2$, where $L_\mathrm{perp}$ is the *perpendicular* distance between opposite faces of the cell (not the longest cell vector). For an orthorhombic cell $L_\mathrm{perp,\alpha} = L_\alpha$; for a triclinic cell with one long axis and a small tilt, $L_\mathrm{perp}$ can be much smaller than the longest edge. LAMMPS computes these perpendicular widths and warns if you cross the threshold.

Why "half"? Because if $r_\mathrm{cut} = L/2$, there are *two* images of atom $j$ at distance exactly $L/2$ from atom $i$ — one in each direction — and the minimum-image convention cannot decide between them. For any $r_\mathrm{cut} > L/2$ multiple images come into range; the convention breaks completely. The strict inequality $r_\mathrm{cut} < L/2$ guarantees uniqueness.

## Long-range Coulomb: Ewald summation

Truncating Coulomb interactions does not work. The $1/r$ tail is too long-ranged: the contribution of pairs beyond any cutoff is not small. Naive truncation produces simulations whose dielectric response, ion mobility and thermodynamics are all wrong.

!!! note "Why splitting is the right idea"
    The Coulomb kernel $1/r$ has two awkward features simultaneously: it is *singular* at $r=0$ and *long-ranged* as $r\to\infty$. No single summation strategy handles both. Real-space sums converge slowly (the tail does not decay fast enough); reciprocal-space sums converge slowly at small $|\mathbf{k}|$ (the singularity at the origin). The Ewald trick is to algebraically split $1/r$ into two pieces, one of which is short-ranged (handle in real space) and the other of which is smooth at the origin (handle in reciprocal space). The standard identity used is the partition of unity
    $$\frac{1}{r} = \frac{\mathrm{erfc}(\alpha r)}{r} + \frac{\mathrm{erf}(\alpha r)}{r},$$
    where $\alpha$ is a tunable parameter with units of inverse length. The first piece decays as $e^{-\alpha^2 r^2}/r$ for large $r$ — short-ranged. The second piece is finite at $r=0$ (since $\mathrm{erf}(\alpha r)/r \to 2\alpha/\sqrt{\pi}$) and represents the smooth long-range part — sum in reciprocal space via Fourier transform.

Ewald's 1921 insight: split the Coulomb interaction into a short-ranged part that sums in real space, and a long-ranged smooth part that sums in reciprocal space. Around each point charge $q_i$ at $\mathbf{r}_i$ place a compensating Gaussian charge density of opposite sign:

$$
\rho_i^\mathrm{Gauss}(\mathbf{r}) = -q_i \left(\frac{\alpha^2}{\pi}\right)^{3/2} e^{-\alpha^2 |\mathbf{r} - \mathbf{r}_i|^2}.
\tag{7.20}
$$

The sum of point charges plus screening Gaussians is short-ranged and can be summed in real space with a small cutoff. The screening Gaussians themselves form a smooth periodic density that is summed in reciprocal space by Fourier transform. The total Coulomb energy becomes

$$
U_\mathrm{Coul} = U_\mathrm{real} + U_\mathrm{recip} - U_\mathrm{self},
\tag{7.21}
$$

where

$$
U_\mathrm{real} = \tfrac{1}{2} \sum_{i \ne j}\, \sum_\mathbf{n}\, q_i q_j \frac{\mathrm{erfc}(\alpha |\mathbf{r}_{ij} + \mathbf{n}|)}{|\mathbf{r}_{ij} + \mathbf{n}|},
\tag{7.22}
$$

$$
U_\mathrm{recip} = \frac{2\pi}{V} \sum_{\mathbf{k} \ne 0} \frac{e^{-k^2/4\alpha^2}}{k^2} |S(\mathbf{k})|^2,
\qquad
S(\mathbf{k}) = \sum_i q_i e^{i\mathbf{k}\cdot \mathbf{r}_i},
\tag{7.23}
$$

and $U_\mathrm{self} = (\alpha/\sqrt{\pi}) \sum_i q_i^2$ corrects for each charge interacting with its own Gaussian. The parameter $\alpha$ controls the split: larger $\alpha$ shifts work to reciprocal space, smaller $\alpha$ to real space. Tuning $\alpha$ optimally gives $O(N^{3/2})$ scaling.

For production work the **particle-mesh Ewald (PME)** method replaces the explicit reciprocal sum (7.23) by evaluation on a regular mesh via FFT, achieving $O(N \log N)$ scaling. PME is the default Coulomb method in GROMACS, AMBER, NAMD, OpenMM and LAMMPS' `kspace_style pppm`. You will rarely write Ewald code yourself; you will tune two numbers: the real-space cutoff and the PME grid spacing. Typical settings for bulk water are a real-space cutoff of 10 Å and a grid spacing of about 1 Å, giving forces accurate to $10^{-5}\,e/$Å.

### Particle-mesh Ewald in three steps

The PME algorithm reduces the long-range Coulomb cost from $O(N^{3/2})$ (best naive Ewald) to $O(N\log N)$. The recipe:

1. **Charge assignment.** For each point charge $q_i$ at $\mathbf{r}_i$, distribute the charge to a small set of nearby grid points using a cardinal B-spline interpolation of order $p$ (typically $p=4$ or 6). This produces a gridded charge density $\rho(\mathbf{R}_g)$ on a regular mesh of $N_g \sim L/h$ points per side, where $h$ is the grid spacing.
2. **FFT the gridded density** to reciprocal space, $\hat\rho(\mathbf{k}) = \sum_g e^{-i\mathbf{k}\cdot\mathbf{R}_g}\rho(\mathbf{R}_g)$. Multiply by the Ewald reciprocal kernel $(4\pi/k^2)e^{-k^2/(4\alpha^2)}$ to get the Fourier-space potential $\hat V(\mathbf{k})$. FFT back to real space, $V(\mathbf{R}_g)$.
3. **Force interpolation.** For each charge, gather the gradient of $V$ at the charge position by differentiating the B-spline interpolant, $\mathbf{F}_i = -q_i \sum_g \nabla M_p(\mathbf{r}_i - \mathbf{R}_g)\,V(\mathbf{R}_g)$.

The FFT step dominates: $O(N_g^3\log N_g) = O(N\log N)$ since $N_g\propto L\propto N^{1/3}$. Charge assignment and force interpolation are $O(N p^3)$ with the small prefactor $p^3 = 64$ for $p=4$. Real-space short-ranged part scales as $O(N)$ with the neighbour list.

**Two tunable parameters.** The Ewald split parameter $\alpha$ and the grid spacing $h$ trade off real-space and reciprocal-space accuracy:

- Larger $\alpha$: real-space erfc decays faster (smaller real-space cutoff possible), but the reciprocal-space Gaussian width $1/\alpha$ shrinks and a finer grid is needed.
- Smaller $\alpha$: coarser grid suffices but real-space cutoff must grow.

In practice you fix the desired force accuracy (typically $10^{-4}$-$10^{-5}\,e$/Å), then solve for the cheapest $\alpha$ given hardware FFT performance. Codes do this automatically: in LAMMPS, `kspace_style pppm 1e-4` asks for a relative force error of $10^{-4}$ and the code picks $\alpha$ and grid.

!!! note "PME is not optional for ionic systems"
    A simulation of NaCl in water with a plain Coulomb cutoff will give wrong radial distribution functions, wrong ion-pair lifetimes, and possibly the wrong sign of the Soret coefficient. Use PME (or PPPM, or its variants) for any system with explicit point charges.

## When PBC bites you

Periodic boundary conditions assume the property of interest is bulk-like and that the unit cell is large enough to suppress spurious interactions across the periodic boundary. Several common situations break this assumption.

### Defects in supercells

A single vacancy in a $4\times 4\times 4$ Si supercell (256 atoms) interacts elastically with all its periodic images, separated by only 22 Å. The strain field of a vacancy decays as $1/r^3$, so the spurious self-interaction is small but not negligible — vacancy formation energies have a finite-size error of order $0.05$–$0.1$ eV at this size. For charged defects (a charged vacancy in a wide-gap insulator), the spurious interaction is Coulomb, $1/r$, and dramatic: $V^\mathrm{q+}$ in MgO at $4\times 4\times 4$ has a self-energy of order eV that must be removed by image-charge corrections (Freysoldt-Neugebauer-Van de Walle, Lany-Zunger).

### Surfaces and slabs

To model a surface you build a slab — a finite-thickness film periodic in two directions, with vacuum above and below. The vacuum needs to be wide enough that the slabs do not interact across it; for metals 10–15 Å is usually adequate. For polar surfaces with non-zero dipole moments perpendicular to the slab, periodic images create a macroscopic electric field through the vacuum, and you must apply a **dipole correction** (the Bengtsson scheme, `dipole_corr` in QE, `LDIPOL` in VASP).

### Permanent dipoles in bulk

Ferroelectric phases have a macroscopic dipole, which under PBC produces an unphysical depolarising field unless you handle the boundary conditions explicitly. The Berry-phase formalism (King-Smith-Vanderbilt) is the modern fix for static dipoles; for MD of strongly polar liquids like water, PME implicitly handles the bulk limit correctly via the "tinfoil" boundary condition.

### Finite-size scaling

For thermodynamic properties of bulk liquids and solids, finite-size effects scale as $1/N$. A diffusion coefficient computed in a 256-atom box and a 2048-atom box should not be wildly different, but they will not be identical either. The Yeh-Hummer correction

$$
D_\infty = D_\mathrm{PBC} + \frac{k_B T\, \xi}{6\pi \eta L},
\qquad \xi \approx 2.837297,
\tag{7.24}
$$

removes the leading hydrodynamic finite-size error from MD self-diffusion coefficients. Use it whenever quoting a diffusion coefficient to more than one significant figure.

## Self-Coulomb correction and tinfoil boundary

A subtle point: Ewald's reciprocal sum (7.23) excludes $\mathbf{k}=0$. The omitted term, sometimes called the "tinfoil" or "metallic" boundary, is finite if the cell has a net dipole moment but undefined if it does not. The physical interpretation is that the system is embedded in a conducting medium (tinfoil), so any net dipole is screened. For most condensed-phase simulations this is the appropriate choice and is the default in LAMMPS' `kspace_style pppm`.

If you simulate something with intrinsic dipole order — a ferroelectric, a polar slab — the tinfoil boundary condition has consequences: the long-range depolarisation field induced by periodic copies of the dipole is removed, which may or may not be what you want for the experimental setup you are modelling.

For polar slabs perpendicular to a periodic direction, the standard fix is a *dipole correction*: a triangular electric field added to cancel the macroscopic dipole. In LAMMPS this is `kspace_modify slab 3.0` (the 3.0 sets the vacuum-to-slab ratio); in VASP `LDIPOL` and `IDIPOL`; in QE `dipfield=.true.` in `&CONTROL`.

## Neighbour lists — how PBC interacts with the inner loop

The minimum-image convention is the geometric foundation; the *algorithmic* expression is the **neighbour list**. For each atom $i$, we maintain a list of all atoms $j$ within $r_\mathrm{cut} + r_\mathrm{skin}$ of $i$, including periodic images. The list is built periodically (every 10-20 steps) and used at each timestep to skip distant pairs.

Two algorithms dominate:

**Cell list (linked-cell)**. Partition the simulation box into cells of side $\ge r_\mathrm{cut}$. Each atom is assigned to one cell. To find neighbours of atom $i$, scan only the 27 cells surrounding $i$'s cell (3D). Cost: $O(N)$ to build, $O(N\cdot M_\mathrm{neighbours\_per\_cell})$ to use.

**Verlet list**. For each atom, store explicitly the indices of nearby atoms. Built from a cell list. Cost: $O(N)$ to build, $O(N\cdot Z)$ to use where $Z$ is the average coordination.

LAMMPS uses both: cell list to find candidates, Verlet list to store the result. The `neighbor 2.0 bin` command sets a 2.0 Å skin: atoms within $r_\mathrm{cut} + 2.0$ Å are kept in the list. The skin is rebuilt only when an atom has moved more than `skin/2` since the last rebuild — checked every `neigh_modify every N` steps.

A practical consequence: the *cell size* $r_\mathrm{cut} + r_\mathrm{skin}$ must be less than $L_\mathrm{perp}/2$ for the minimum-image convention to apply within neighbour search. If your cell is small, LAMMPS uses Newton's third law to fold pair forces back to a single computation per pair, but the geometry still requires the cell-size constraint.

## What we have

We can now integrate atoms in time (§7.1) inside a topologically toroidal cell that mimics an infinite material (this section). What we cannot yet do is set the temperature or pressure — the equations are still pure Newton, which conserves energy and volume. Adding control over thermodynamic variables is [§7.3](03-thermostats.md).
