# 7.2 Periodic Boundary Conditions

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

## Long-range Coulomb: Ewald summation

Truncating Coulomb interactions does not work. The $1/r$ tail is too long-ranged: the contribution of pairs beyond any cutoff is not small. Naive truncation produces simulations whose dielectric response, ion mobility and thermodynamics are all wrong.

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

## What we have

We can now integrate atoms in time (§7.1) inside a topologically toroidal cell that mimics an infinite material (this section). What we cannot yet do is set the temperature or pressure — the equations are still pure Newton, which conserves energy and volume. Adding control over thermodynamic variables is [§7.3](03-thermostats.md).
