# Methods — Defect formation energy in silicon

This is a twelve-step protocol. The intention is that an undergraduate
who has worked through chapters 4–6 can follow it end-to-end without
external guidance. Where a step looks short, that is because it is
short; where it looks long, the length is genuinely warranted by the
number of pitfalls.

---

## Step 0 — Notation and conventions

- Energies are reported in eV unless otherwise stated. Quantum
  ESPRESSO uses Ry internally (1 Ry = 13.6057 eV); always convert at
  the boundary, never in the middle of an analysis.
- Lengths are in Å. QE uses bohr (1 bohr = 0.52918 Å).
- $N$ is the number of atoms in the *defect-free* supercell. The
  vacancy cell therefore has $N - 1$ atoms.
- $\mu_\mathrm{Si} = E_\mathrm{tot}(\text{bulk Si}) / N_\mathrm{atoms,bulk}$.
- The formation energy of the neutral vacancy is
  $E_f^{V^0} = E_\mathrm{tot}(V_\mathrm{Si}^0, N{-}1) - \tfrac{N-1}{N}\,E_\mathrm{tot}(\text{bulk}, N).$

  Use the second form — fractional bulk — to remove a common
  bookkeeping error.

---

## Step 1 — Install the toolchain

Install Quantum ESPRESSO from source (preferred on a cluster) or via
`conda` for a laptop:

```bash
conda create -n defects python=3.11 ase pymatgen matplotlib numpy
conda activate defects
conda install -c conda-forge qe
which pw.x  # verify
```

Download the SSSP "efficiency" PBE pseudopotential library and unpack
it into `$HOME/pseudo/sssp_efficiency_PBE/`. Set the environment
variable:

```bash
export ESPRESSO_PSEUDO=$HOME/pseudo/sssp_efficiency_PBE
```

Verify the toolchain by running a hydrogen-atom SCF (the smallest
non-trivial test). If `pw.x` finishes with `JOB DONE.`, proceed.

---

## Step 2 — Bulk silicon: lattice parameter and pseudopotential check

Use the SSSP-recommended `Si.pbe-n-rrkjus_psl.1.0.0.UPF` pseudopotential
and the SSSP-recommended cut-off (40 Ry wavefunction, 320 Ry charge
density).

Compute the equilibrium lattice parameter $a_0$ by an equation-of-state
fit: run SCF at seven lattice parameters in the range 5.35–5.55 Å on
the 2-atom primitive cell, fit a Birch–Murnaghan EoS via ASE's
`ase.eos.EquationOfState`, and report $a_0$. You should get
5.469 ± 0.005 Å (the PBE value).

This $a_0$ is the lattice parameter you must use for *all* subsequent
calculations. Do not relax the cell of the defect supercell.

The `pw.x` input for one EoS point looks like:

```text
&CONTROL
  calculation = 'scf'
  prefix      = 'Si_bulk'
  pseudo_dir  = '/home/user/pseudo/sssp_efficiency_PBE'
  outdir      = './tmp/'
  verbosity   = 'low'
/
&SYSTEM
  ibrav      = 2          ! face-centred cubic
  celldm(1)  = 10.336     ! lattice parameter in bohr; convert from your value
  nat        = 2
  ntyp       = 1
  ecutwfc    = 40.0
  ecutrho    = 320.0
  occupations = 'fixed'
/
&ELECTRONS
  conv_thr    = 1.0d-9
  mixing_beta = 0.5
/
ATOMIC_SPECIES
  Si  28.0855  Si.pbe-n-rrkjus_psl.1.0.0.UPF
ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25
K_POINTS automatic
  8 8 8 0 0 0
```

The `celldm(1)` value of 10.336 bohr is 5.469 Å — convert with
`celldm = a / 0.529177`.

---

## Step 3 — Plane-wave cut-off convergence

On the 2-atom primitive cell, run SCF at $E_\mathrm{cut}^\mathrm{wfc} =
30, 35, 40, 45, 50, 60$ Ry, keeping
$E_\mathrm{cut}^\mathrm{rho} = 8 E_\mathrm{cut}^\mathrm{wfc}$. Plot the
total energy per atom against the cut-off and select the value at
which successive points differ by less than 1 meV per atom. With the
SSSP-recommended pseudopotential this is reached at 40 Ry.

Document this in `convergence/ecut.png` and `convergence/ecut.csv`.

---

## Step 4 — **k**-mesh convergence

Hold $E_\mathrm{cut}^\mathrm{wfc}$ at the converged value. Sweep the
**k**-mesh on the 2-atom cell over $4^3, 6^3, 8^3, 10^3, 12^3$. The
energy should be converged to 1 meV per atom by $8^3$ for silicon.

The relevant principle for the defect work is the *equivalent
**k**-mesh*: the BZ of an $n \times n \times n$ supercell is $1/n^3$
that of the primitive cell. Therefore an $8 \times 8 \times 8$
primitive mesh is equivalent to:

| Supercell | Conventional cell multiple | Equivalent k-mesh |
| --- | --- | --- |
| 64-atom (2×2×2 of conventional 8-atom) | $n = 2$ | $4 \times 4 \times 4$ |
| 216-atom (3×3×3) | $n = 3$ | $3 \times 3 \times 3$ (use 2×2×2 in practice — Γ-centred) |
| 512-atom (4×4×4) | $n = 4$ | $2 \times 2 \times 2$ |

Always include Γ for the defect cells (Γ-centred mesh, no shift).

---

## Step 5 — Build the perfect supercells

Use ASE rather than writing `pw.x` input by hand:

```python
from ase.build import bulk
from ase.io import write

si = bulk("Si", crystalstructure="diamond", a=5.469, cubic=True)
# 'cubic=True' returns the 8-atom conventional cell.
super64 = si.repeat((2, 2, 2))   # 64 atoms
super216 = si.repeat((3, 3, 3))  # 216 atoms
super512 = si.repeat((4, 4, 4))  # 512 atoms
write("Si_N064.xyz", super64)
write("Si_N216.xyz", super216)
write("Si_N512.xyz", super512)
```

For each cell, generate a `pw.x` SCF input with the appropriate
equivalent **k**-mesh from step 4. Run SCF and record the total energy
$E_\mathrm{tot}(\text{bulk}, N)$.

---

## Step 6 — Make the vacancy

Take the perfect supercell, remove atom 0 (or any single atom — the
choice is arbitrary by translational symmetry), and write the
result. Crucially, *perturb* the four neighbours of the missing atom
by 0.05 Å in random directions to allow the Jahn–Teller distortion to
break $T_d$:

```python
import numpy as np
from ase import Atoms

def make_vacancy(perfect: Atoms, index: int = 0,
                 perturb: float = 0.05, rng_seed: int = 42) -> Atoms:
    rng = np.random.default_rng(rng_seed)
    vac = perfect.copy()
    centre = vac.positions[index]
    # Identify four nearest neighbours by distance (diamond NN distance is a*sqrt(3)/4).
    d = np.linalg.norm(vac.positions - centre, axis=1)
    nn = np.argsort(d)[1:5]  # exclude the atom itself
    del vac[index]
    # The indices in 'vac' are now one less above the removed atom; recompute.
    centre = perfect.positions[index]
    d_new = np.linalg.norm(vac.positions - centre, axis=1)
    nn_new = np.argsort(d_new)[:4]
    displacements = perturb * rng.standard_normal(size=(4, 3))
    vac.positions[nn_new] += displacements
    return vac
```

Record the random seed so the geometry is reproducible.

---

## Step 7 — Write the defect `pw.x` input

Two crucial flags:

```text
&SYSTEM
  ...
  nosym = .true.       ! do not symmetrise: allow J-T distortion
  occupations = 'smearing'
  smearing = 'gaussian'
  degauss  = 0.005     ! 0.07 eV — soft enough to converge SCF, not so soft as to wash out the gap
/
&IONS
  ion_dynamics = 'bfgs'
/
```

For the *perfect* cell you must use the **same** `smearing` setting so
that energies are comparable. Re-run the perfect bulk SCF with these
flags if you ran it earlier without smearing.

Force the calculation to be a `vc-relax = .false.`, i.e. a `relax`
calculation, *not* a variable-cell relaxation: never relax the cell
parameter of a defect supercell.

---

## Step 8 — Geometry relaxation

Submit the vacancy relaxation. Useful convergence thresholds:

- `etot_conv_thr = 1.0d-5` (Ry, total-energy change between ionic steps)
- `forc_conv_thr = 1.0d-4` (Ry/bohr, ≈ 2.6 meV/Å)

Expect 20–40 ionic steps to converge from the perturbed start. Each
SCF inside that is one call to `pw.x`. The 64-atom relaxation finishes
in roughly 30 CPU-hours on modern Xeon nodes; the 216-atom case takes
≈ 200; the 512-atom case ≈ 600. Use checkpoint restart
(`startingwfc = 'file'`, `startingpot = 'file'`) liberally so that a
wall-clock timeout does not lose your progress.

---

## Step 9 — Verify the Jahn–Teller distortion

Read the relaxed structure with ASE and compute the four
nearest-neighbour distances from the vacancy site (which is just the
removed atom's original position). The four NN distances should pair
into two short (~2.30 Å) and two long (~2.42 Å) values, characteristic
of the $D_{2d}$ distortion. Use `pymatgen.symmetry.analyzer.PointGroupAnalyzer`
to confirm the local symmetry of a cluster cut around the vacancy.

If all four distances are equal, your calculation has stayed in $T_d$.
Re-run with a larger perturbation or check that `nosym = .true.` was
honoured.

---

## Step 10 — Compute $E_f$ for each cell

```python
def formation_energy(e_defect: float, e_bulk: float, n_bulk: int) -> float:
    """Neutral monovacancy formation energy.

    All energies in eV; n_bulk is the number of atoms in the
    defect-free supercell."""
    mu_si = e_bulk / n_bulk
    return e_defect - (n_bulk - 1) * mu_si
```

Tabulate $E_f$ vs $N$ in `analysis/summary.csv`.

---

## Step 11 — Finite-size extrapolation

The leading correction to a neutral defect's formation energy in a
periodic supercell scales as $1/N$ (equivalently, $1/L^3$ where $L$
is the cubic supercell parameter), arising from the elastic image
interaction. Fit

$$
E_f(N) = E_f^\infty + \alpha\, N^{-1}.
$$

Use weighted least squares with weights ∝ $N$ (the smaller cells are
less informative). Report $E_f^\infty$ and the standard error of the
fit. Also fit against $L^{-3}$ as a cross-check; the intercept should
agree within the uncertainty.

A representative result with three points (64, 216, 512) should give
$E_f^\infty \approx 3.55$–$3.65$ eV with a standard error well below
0.05 eV.

---

## Step 12 — Write up and compare

Compare with at least two literature values. Reasonable choices:

- Corsetti and Mostofi (2011): $3.56 \pm 0.05$ eV (PBE, asymptotic).
- Wright (2006): $3.52$ eV (PBE, 216-atom).
- Probert and Payne (2003): $\approx 3.6$ eV (denser **k**-mesh
  variant).

State plainly whether your number lies within their reported scatter.
If it does not, attempt to identify the cause: most often it is a
pseudopotential difference, an inconsistent **k**-mesh, or a residual
symmetry constraint preventing Jahn–Teller relaxation.

---

## Pitfalls

1. **Forgetting to relax.** Reporting the unrelaxed vacancy energy
   gives a value of ≈ 4.0 eV. The 0.4 eV gain on relaxation is the
   Jahn–Teller energy. Always relax.
2. **Inconsistent smearing.** As stated in step 7, perfect and defect
   cells must use identical smearing.
3. **Inconsistent **k**-meshes.** As stated in step 4, use the
   equivalent-mesh prescription.
4. **Mistaking the chemical potential.** For elemental Si, the
   reservoir is bulk Si itself; do not introduce silane or any other
   reference.
5. **Re-using `outdir`.** If you run multiple calculations from the
   same directory, give each a unique `prefix` and `outdir`, or QE
   will quietly overwrite charge-density files and produce subtly
   wrong restart behaviour.
6. **Comparing different functionals.** If you used PBE for the
   bulk lattice constant, use PBE for the vacancy. Do not mix.
7. **Negative formation energies.** If you see one, you have almost
   certainly made a sign or atom-count error in the formation-energy
   formula — re-check $\mu_\mathrm{Si} = E_\mathrm{bulk} / N$.

If at any step a calculation crashes or produces nonsensical numbers,
inspect the QE output for the lines `total energy`, `convergence has
been achieved`, and `Forces acting on atoms`. The output is verbose;
read it.
