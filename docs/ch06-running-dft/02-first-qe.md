# 6.2 Your first Quantum ESPRESSO calculation

```mermaid
flowchart LR
    A["<b>Input file</b><br/>silicon.scf.in<br/>(geometry, pseudos,<br/>k-points, E_cut)"]
    B["<b>pw.x</b><br/>plane-wave SCF<br/>solver"]
    C["<b>Output files</b><br/>silicon.scf.out<br/>+ outdir/*.save"]
    D["<b>Analysis</b><br/>energies, forces,<br/>charge density"]
    P["Pseudopotentials<br/>(.UPF)"]
    A --> B
    P --> B
    B --> C --> D
```
*The standard Quantum ESPRESSO workflow: input + pseudopotentials feed `pw.x`, which writes both a human-readable log and a binary save directory for downstream analysis.*

We compute the ground state of silicon. By the end of this section you will have run `pw.x` on a real input file, inspected the output, and driven the same calculation from Python via ASE.

## 6.2.1 Installing Quantum ESPRESSO

QE is a suite of Fortran/MPI codes. The main executable for self-consistent DFT is `pw.x`. There are platform-specific shortcuts:

=== "macOS (Homebrew)"

    ```bash
    brew install quantum-espresso
    pw.x --version
    ```

=== "Linux (apt)"

    ```bash
    sudo apt install quantum-espresso
    pw.x --version
    ```

=== "Conda (any platform)"

    ```bash
    conda install -c conda-forge qe
    pw.x --version
    ```

=== "From source"

    ```bash
    git clone https://gitlab.com/QEF/q-e.git
    cd q-e
    ./configure
    make pw pp -j 8
    export PATH=$PWD/bin:$PATH
    pw.x --version
    ```

Expected output is something like:

```text
     Program PWSCF v.7.3.1 starts on ...
```

If `pw.x --version` reports anything older than v.7.2, upgrade. Older versions have different defaults for some flags and will give slightly different numbers from those quoted here.

### Pseudopotentials

Download SSSP-PBE-efficiency 1.3.0 from [materialscloud.org](https://www.materialscloud.org/discover/sssp). Unpack to a directory of your choice:

```bash
mkdir -p ~/pseudo
cd ~/pseudo
# unpack SSSP_1.3.0_PBE_efficiency.tar.gz here
tar -xzf SSSP_1.3.0_PBE_efficiency.tar.gz
export ESPRESSO_PSEUDO=$HOME/pseudo/SSSP_1.3.0_PBE_efficiency
```

Add the `export` line to your `~/.zshrc` or `~/.bashrc` so it persists.

Verify silicon is there:

```bash
ls $ESPRESSO_PSEUDO | grep -i ^si
# Si.pbe-n-rrkjus_psl.1.0.0.UPF
```

The exact file name might drift between SSSP versions; we will refer to it as `Si.UPF` symbolically and use the SSSP filename in inputs.

### Scratch directory

QE writes large temporary files during SCF (wavefunctions, charge density). Point it somewhere with room:

```bash
mkdir -p ~/qe-scratch
export ESPRESSO_TMPDIR=$HOME/qe-scratch
```

## 6.2.2 Anatomy of a `pw.x` input file

A `pw.x` input has five required sections (`&CONTROL`, `&SYSTEM`, `&ELECTRONS`, `ATOMIC_SPECIES`, `ATOMIC_POSITIONS`) and one nearly always required section (`K_POINTS`), plus optional ones (`&IONS`, `&CELL`, `CELL_PARAMETERS`, etc.).

A minimal annotated skeleton:

```fortran
&CONTROL
  calculation = 'scf'           ! 'scf','nscf','bands','relax','vc-relax','md'
  prefix      = 'si'            ! tag for output files
  pseudo_dir  = '/path/to/upf'  ! optional; overrides $ESPRESSO_PSEUDO
  outdir      = './tmp'         ! scratch dir; overrides $ESPRESSO_TMPDIR
  verbosity   = 'low'           ! 'low' or 'high'
  tprnfor     = .true.          ! print forces
  tstress     = .true.          ! print stress tensor
/
&SYSTEM
  ibrav      = 0                ! 0 = lattice given by CELL_PARAMETERS card
  nat        = 2                ! number of atoms in cell
  ntyp       = 1                ! number of distinct species
  ecutwfc    = 40.0             ! plane-wave cutoff for wavefunctions (Ry)
  ecutrho    = 320.0            ! plane-wave cutoff for density       (Ry)
  occupations= 'fixed'          ! 'fixed' (insulators) or 'smearing' (metals)
/
&ELECTRONS
  conv_thr   = 1.0d-8           ! SCF energy threshold (Ry)
  mixing_beta= 0.4              ! charge mixing parameter
/
ATOMIC_SPECIES
  Si  28.085  Si.pbe-n-rrkjus_psl.1.0.0.UPF

CELL_PARAMETERS angstrom
  0.0000  2.7150  2.7150
  2.7150  0.0000  2.7150
  2.7150  2.7150  0.0000

ATOMIC_POSITIONS crystal
  Si  0.00  0.00  0.00
  Si  0.25  0.25  0.25

K_POINTS automatic
  4 4 4   0 0 0
```

A few things worth knowing.

- **`ibrav`**. Either `0` (lattice in `CELL_PARAMETERS`) or one of QE's 14 built-in Bravais lattice codes. For example `ibrav=2` is FCC; you then specify only `celldm(1) = a0` (the conventional cubic lattice parameter in Bohr). We use `ibrav=0` because it is explicit and unambiguous.
- **`celldm` vs `CELL_PARAMETERS`**. If `ibrav /= 0`, the lattice is built from `celldm(1)`-`celldm(6)`. If `ibrav = 0`, you write the three primitive lattice vectors as rows in `CELL_PARAMETERS`.
- **Coordinate types in `ATOMIC_POSITIONS`**: `alat` (units of `celldm(1)`), `bohr`, `angstrom`, or `crystal` (fractional). `crystal` is safest for symmetry.
- **`K_POINTS automatic`** triggers a Monkhorst-Pack grid; the six integers are $N_1\,N_2\,N_3\,s_1\,s_2\,s_3$.

!!! warning "Forget `occupations` for a metal and the SCF will not converge"
    By default `pw.x` assumes insulating occupations (`fixed`). For a metal you must add `occupations = 'smearing'`, `smearing = 'mv'` (Marzari-Vanderbilt cold smearing is a safe default), `degauss = 0.01` (Ry, roughly $k_B T / 7$ for $T = 300$ K). Forgetting this for a metallic system leads to charge sloshing, oscillating energies and SCF failure. Silicon is an insulator so it does not bite us here, but Al, Cu, Fe will.

!!! warning "k-grid for spin-polarised systems"
    Adding `nspin = 2` doubles the number of orbital channels but does *not* automatically densify the k-grid. Magnetic systems often need denser sampling than their non-magnetic counterparts because the Fermi surface differs between spin channels. Convergence-test separately.

## 6.2.3 A complete silicon input

We will use the **8-atom conventional FCC cell** rather than the 2-atom primitive cell, partly because it is conceptually simpler (a cube), partly because the same file will be useful for defect calculations later. Silicon's conventional cubic lattice parameter is $a = 5.43$ Å, and the cell contains 8 Si atoms in the diamond structure (two interpenetrating FCC sublattices, the second offset by $(\tfrac14,\tfrac14,\tfrac14)a$).

Save this as `si.scf.in`:

```fortran
&CONTROL
  calculation  = 'scf'
  prefix       = 'si'
  outdir       = './tmp/'
  pseudo_dir   = './pseudo/'
  verbosity    = 'high'
  tprnfor      = .true.
  tstress      = .true.
  wf_collect   = .true.
/
&SYSTEM
  ibrav        = 0
  nat          = 8
  ntyp         = 1
  ecutwfc      = 40.0
  ecutrho      = 320.0
  occupations  = 'fixed'
  nbnd         = 16
/
&ELECTRONS
  electron_maxstep = 100
  conv_thr         = 1.0d-8
  mixing_mode      = 'plain'
  mixing_beta      = 0.4
  diagonalization  = 'david'
/
ATOMIC_SPECIES
  Si  28.085  Si.pbe-n-rrkjus_psl.1.0.0.UPF

CELL_PARAMETERS angstrom
  5.43000000  0.00000000  0.00000000
  0.00000000  5.43000000  0.00000000
  0.00000000  0.00000000  5.43000000

ATOMIC_POSITIONS crystal
  Si  0.000  0.000  0.000
  Si  0.000  0.500  0.500
  Si  0.500  0.000  0.500
  Si  0.500  0.500  0.000
  Si  0.250  0.250  0.250
  Si  0.250  0.750  0.750
  Si  0.750  0.250  0.750
  Si  0.750  0.750  0.250

K_POINTS automatic
  4 4 4   0 0 0
```

A few decisions explained:

- **`nbnd = 16`** asks for 16 Kohn-Sham bands. The eight Si atoms have 8 × 4 = 32 valence electrons, occupying 16 bands (with spin degeneracy). For an scf this is enough; for plotting bands or DOS we will request more.
- **`mixing_mode = 'plain'`, `mixing_beta = 0.4`** is the simplest charge mixer with a moderately damped update. For tricky cases (transition metals, magnetic insulators) try `'local-TF'` or lower `mixing_beta`.
- **`diagonalization = 'david'`** is the Davidson iterative diagonaliser. Switch to `'cg'` (conjugate gradient) if Davidson fails — slower per step but more robust.
- **No symmetry flags**: QE auto-detects symmetry; the diamond structure has $Fd\bar{3}m$ symmetry, and the 4×4×4 unshifted MP grid reduces to **10 irreducible k-points**, which is what you should see in the output.

Make sure the pseudopotential file is reachable. Either copy it into `./pseudo/`:

```bash
mkdir -p ./pseudo
cp $ESPRESSO_PSEUDO/Si.pbe-n-rrkjus_psl.1.0.0.UPF ./pseudo/
```

or set `pseudo_dir = '/absolute/path/to/SSSP_1.3.0_PBE_efficiency/'` in the input.

## 6.2.4 Running it

From a terminal in the directory containing `si.scf.in`:

```bash
mkdir -p tmp
pw.x -inp si.scf.in > si.scf.out
```

On a 2023 MacBook (M2), this finishes in about 20 seconds with 4 threads. For parallel runs over MPI:

```bash
mpirun -np 4 pw.x -inp si.scf.in > si.scf.out
```

QE supports several parallel layers (`-nk`, `-nb`, `-nt`, `-nd`) — see the QE manual for `-nk` (k-point parallelisation), which is the most useful one. For our small Si run any value of `-nk` from 1 to 4 is reasonable.

## 6.2.5 Reading the output

Open `si.scf.out`. The interesting passages, in order:

**Header — version and parallelisation:**

```text
     Program PWSCF v.7.3.1 starts on  ...
     Parallel version (MPI), running on     1 processors
```

**Symmetry and k-points:**

```text
     number of Bravais lattice symmetry operations =   48
     number of inequivalent operations =   48
     ...
     number of k points=    10
```

The diamond structure has 48 point-group operations; the unshifted 4×4×4 MP grid reduces by symmetry to 10 irreducible k-points.

**Basis size:**

```text
     Sum of charges of pseudopotential ZVAL = 4 (Si)
     ...
     Plane waves found:    8771 (for kinetic energy < ecutwfc)
     Density plane waves: 35073 (for kinetic energy < ecutrho)
```

**SCF iterations:**

```text
     iteration #  1     ecut=    40.00 Ry      beta= 0.40
     Davidson diagonalization with overlap
     ethr =  1.00E-02,  avg # of iterations =  3.0
     total cpu time spent up to now is        1.2 secs
     total energy              =     -93.45412318 Ry
     estimated scf accuracy    <       0.07810122 Ry
     ...
     convergence has been achieved in   8 iterations
```

A typical Si SCF needs 7-12 iterations to hit `conv_thr = 1e-8`.

**Final results:**

```text
!    total energy              =     -93.45698312 Ry
     Harris-Foulkes estimate   =     -93.45698312 Ry
     estimated scf accuracy    <          1.0E-09 Ry

     The total energy is the sum of the following terms:
     one-electron contribution =     ...
     hartree contribution      =     ...
     xc contribution           =     ...
     ewald contribution        =     ...

     Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =      0.0000  0.0000  0.0000
     ...

     total   stress  (Ry/bohr**3)         (kbar)
       -0.00021    0.0   0.0          -30.5   0.0   0.0
        0.0  -0.00021  0.0             0.0 -30.5   0.0
        0.0   0.0   -0.00021           0.0   0.0 -30.5
```

The line beginning `!    total energy` is the SCF total energy: $E_\mathrm{tot} \approx -93.457$ Ry $\approx -1271.5$ eV. Divide by 8 atoms to get $-158.94$ eV/atom.

Forces are exactly zero by symmetry — silicon's diamond positions are stationary. The stress is approximately $-30$ kbar (about $-3$ GPa) at this lattice parameter and cutoff, meaning the cell wants to expand slightly. To find the true equilibrium $a$, you would run `calculation = 'vc-relax'`. We will not do that here; we accept the experimental $a = 5.43$ Å.

**Eigenvalues at each k-point:**

```text
     k =-0.3750-0.3750 0.1250 (   1097 PWs)   bands (ev):
       -5.6011  -2.5683  -1.0394  -1.0394   3.1234   3.1234   6.4218   6.4218
        7.1234   7.8923   ...
```

Band 1-8 are the occupied valence bands (since 8 atoms × 4 electrons / 2 spin = 16 occupied per spin, but the conventional cell has half the BZ of the primitive, so 8 occupied bands at each k-point in the conventional cell). The Fermi energy:

```text
     highest occupied, lowest unoccupied level (ev):     5.6798    6.3134
```

The PBE-predicted gap is $6.3134 - 5.6798 = 0.63$ eV. Experimental Si gap is 1.17 eV. This is the famous "**PBE band gap underestimation**" — PBE consistently gives ~50% of the experimental gap for sp-bonded semiconductors. We will see how to do better in later chapters (hybrid functionals, GW).

## 6.2.6 Driving QE from Python with ASE

For anything beyond a one-off calculation — convergence tests, sweeps, defect supercells, equation-of-state fits — you want a programmatic driver. ASE (the Atomic Simulation Environment) provides this. Install:

```bash
pip install "ase>=3.23"
```

ASE represents a system as an `Atoms` object and dispatches the actual energy/force evaluation to a *calculator*. The `Espresso` calculator writes a `pw.x` input, runs the binary, and parses the output.

### Profile setup

ASE 3.23+ uses `EspressoProfile` to encapsulate the QE binary and pseudopotential paths. A one-time setup in your script:

```python
from pathlib import Path
from ase.calculators.espresso import EspressoProfile

profile = EspressoProfile(
    command="pw.x",                                  # or "mpirun -np 4 pw.x"
    pseudo_dir=Path.home() / "pseudo/SSSP_1.3.0_PBE_efficiency",
)
```

### Single-point silicon

```python
from __future__ import annotations
from pathlib import Path
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile


def make_si_calculator(workdir: Path, kpts: tuple[int, int, int] = (4, 4, 4),
                       ecutwfc: float = 40.0, ecutrho: float = 320.0) -> Espresso:
    """Return a configured Quantum ESPRESSO calculator for silicon."""
    profile = EspressoProfile(
        command="pw.x",
        pseudo_dir=Path.home() / "pseudo/SSSP_1.3.0_PBE_efficiency",
    )
    pseudopotentials = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}
    input_data: dict = {
        "control": {
            "calculation": "scf",
            "verbosity": "high",
            "tprnfor": True,
            "tstress": True,
        },
        "system": {
            "ecutwfc": ecutwfc,
            "ecutrho": ecutrho,
            "occupations": "fixed",
        },
        "electrons": {
            "conv_thr": 1.0e-8,
            "mixing_beta": 0.4,
        },
    }
    return Espresso(
        profile=profile,
        directory=str(workdir),
        pseudopotentials=pseudopotentials,
        input_data=input_data,
        kpts=kpts,
    )


def main() -> None:
    si = bulk("Si", crystalstructure="diamond", a=5.43)  # 2-atom primitive cell
    workdir = Path("./si_scf")
    workdir.mkdir(exist_ok=True)
    si.calc = make_si_calculator(workdir)
    energy_eV = si.get_potential_energy()
    forces_eV_per_A = si.get_forces()
    stress_eV_per_A3 = si.get_stress()
    print(f"E_tot = {energy_eV:.6f} eV  ({energy_eV / len(si):.6f} eV/atom)")
    print(f"max |F| = {abs(forces_eV_per_A).max():.3e} eV/Å")
    print(f"stress (xx,yy,zz,yz,xz,xy) eV/Å³: {stress_eV_per_A3}")


if __name__ == "__main__":
    main()
```

Run it:

```bash
python si_scf.py
```

ASE creates `./si_scf/espresso.pwi` (input) and `./si_scf/espresso.pwo` (output), invokes `pw.x`, parses the output, and gives you:

```text
E_tot = -317.876214 eV  (-158.938107 eV/atom)
max |F| = 5.241e-05 eV/Å
stress (xx,yy,zz,yz,xz,xy) eV/Å³: [-0.001895 -0.001895 -0.001895 0 0 0]
```

Note that ASE converts everything to eV and eV/Å — its convention — regardless of QE's internal Ry. The energy per atom $-158.94$ eV/atom matches the 8-atom cell calculation up to round-off (the 2-atom primitive and 8-atom conventional cells are equivalent in the thermodynamic limit; small differences come from the k-grid being defined per cell, not per atom — a $4\times4\times4$ grid on the primitive cell is denser than on the conventional cell).

### Saving and reloading

ASE writes the parsed output to `espresso.pwo`. You can re-read it later:

```python
from ase.io import read
si = read("si_scf/espresso.pwo")
print(si.get_potential_energy())  # cached from the file, no re-run
```

This is the foundation of every workflow in the rest of the chapter: build an `Atoms`, attach a calculator, get energy/forces/stress, repeat with different parameters.

### Where ASE helps and where it does not

ASE excels at scripted scans: convergence tests, equations of state, defect supercells, surface slabs, NEB calculations. It hides the QE input syntax behind Python dicts.

It does *not* hide the physics. You still have to know that you need `occupations = 'smearing'` for a metal, that `nspin = 2` requires a starting magnetisation, that `ecutrho = 4*ecutwfc` is wrong for USPP. ASE is a driver, not an oracle.

## 6.2.7 What you have

You have run a complete DFT calculation. You have a total energy, forces (zero by symmetry), and a stress tensor. You can drive the same calculation from Python and you have a basis for everything that follows.

What you do not have is a *converged* calculation. The values 40 Ry and $4\times4\times4$ were guesses, informed by the SSSP recommendation but not verified for this particular cell. In the next section we test convergence — and find out which of those numbers we got away with.
