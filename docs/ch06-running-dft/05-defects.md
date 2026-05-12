# 6.5 Defects and formation energies

!!! tip "Why defect calculations exist"
    Real materials are defective: every silicon wafer has $10^{14}$-$10^{16}$ point defects per cm$^3$ at thermal equilibrium, and many more during processing. Defects control nearly all the practical properties of a material: dopant atoms enable transistors; vacancies and interstitials govern atomic diffusion (and therefore semiconductor processing, battery degradation, metal corrosion); deep levels destroy minority-carrier lifetime in solar cells. To predict any of these from first principles you need the formation energy of each defect — the thermodynamic cost of creating one, which sets the equilibrium concentration via $c \propto e^{-E_f/k_BT}$.
    
    The picture: take a chunk of perfect crystal and remove one atom. The total energy changes by some amount. Add an atom back to a chemical reservoir (e.g. metallic Si). The defect formation energy is the leftover. Doing this in a small box (a "supercell") with periodic boundary conditions creates an *array* of defects, not an isolated one; the calculation has artefacts proportional to the inverse cell size. The art is making the cell big enough that the artefacts are smaller than the answer.

A perfect crystal exists nowhere. Every real silicon wafer contains $10^{14}$-$10^{16}$ defects per cm$^3$ at thermal equilibrium — vacancies, interstitials, substitutional impurities, dislocations — and every one of them matters for some practical purpose.

This section explains how to compute the **formation energy** of a point defect from DFT, using the silicon vacancy as our example. The formation energy is the thermodynamic price you pay to create the defect; it controls the equilibrium concentration via the Boltzmann factor

$$ c_\mathrm{defect} = N_\mathrm{sites}\, e^{-E_f / k_B T}. $$

A defect with $E_f = 2$ eV is present at concentrations near melting; a defect with $E_f = 4$ eV is rare except under non-equilibrium conditions like irradiation.

## 6.5.1 Why defects matter

Three application domains drive most defect calculations.

### Doping

Silicon for transistors is doped with boron (p-type) or phosphorus (n-type). The dopant substitutes for a Si atom and contributes one extra hole or electron. Whether the dopant prefers a substitutional or interstitial site, and whether it ionises completely at room temperature, are first-principles questions answered by computing $E_f$ for each configuration and each charge state.

### Diffusion

Atomic transport in solids — and therefore semiconductor processing, battery degradation, metal corrosion — proceeds by point-defect motion. A vacancy diffuses by hopping into a neighbour site; an interstitial squeezes between lattice sites. The **migration barrier** is computed by NEB (nudged elastic band) calculations between defect configurations. The defect concentration times the hop rate gives the macroscopic diffusion coefficient.

### Recombination centres

In photovoltaics and LEDs, "deep" defect levels in the band gap act as Shockley-Read-Hall recombination centres, killing minority-carrier lifetime. The defect's *charge transition level* — the Fermi level at which one charge state becomes more stable than another — is a key DFT-computable quantity.

For a starter calculation we will compute the formation energy of the *neutral* silicon vacancy, $V_\mathrm{Si}^0$. We will not attempt charged defects in detail; a brief mention follows.

## 6.5.2 The supercell approach

DFT calculations are periodic. To simulate an *isolated* defect we use a **supercell**: a large cell, periodically repeated, with one defect per cell. If the cell is large enough, neighbouring defects in adjacent cells do not interact and we obtain the dilute-limit formation energy.

For silicon, the standard starting supercell is **a $2\times2\times2$ tiling of the 8-atom conventional cubic cell**, giving 64 atoms in a cubic box of side $2a \approx 10.86$ Å. Remove one atom and you have a 63-atom cell with one vacancy. The vacancy-vacancy distance across the periodic boundary is then $2a \approx 10.86$ Å, which is large enough that the elastic and electronic interactions are small for the neutral vacancy. For charged defects you need bigger cells — 216 or 512 atoms, see §6.5.4.

## 6.5.3 Formation energy: the formula

For a *neutral* defect in a single-component crystal, the formation energy is:

$$ E_f[V_\mathrm{Si}] = E[\mathrm{defective}] - E[\mathrm{perfect}] + \mu_\mathrm{Si} \tag{6.2} $$

where:

- $E[\mathrm{defective}]$ is the total DFT energy of the supercell with the vacancy.
- $E[\mathrm{perfect}]$ is the total DFT energy of the same supercell without the vacancy.
- $\mu_\mathrm{Si}$ is the chemical potential of silicon. We added a Si atom back to the *reservoir* when we created the vacancy, so we add its chemical potential to the energy.

For elemental silicon at thermal equilibrium with bulk silicon, $\mu_\mathrm{Si}$ is simply the total energy per atom of bulk Si:

$$ \mu_\mathrm{Si} = \frac{E[\mathrm{perfect}]}{N_\mathrm{perfect}}. $$

So with $N_\mathrm{perfect} = 64$:

$$ E_f[V_\mathrm{Si}] = E[\mathrm{defective, 63 atoms}] - \frac{63}{64}\, E[\mathrm{perfect, 64 atoms}]. \tag{6.3} $$

The "$63/64$" is the bookkeeping: we are comparing 63 atoms in a defective cell against 63 atoms' worth of bulk silicon.

For binary compounds (say, GaAs), the chemical potentials of the two species are not independent — they are constrained by the formation energy of the compound and are bounded above and below by elemental and impurity phases. The general expression is

$$ E_f^q[D] = E_\mathrm{tot}[D^q] - E_\mathrm{tot}[\mathrm{bulk}] - \sum_i n_i \mu_i + q (E_F + E_\mathrm{VBM} + \Delta V) + E_\mathrm{corr}^q, \tag{6.4} $$

where $n_i$ is the change in the number of species $i$ (negative for atoms removed, positive for atoms added), $q$ is the charge state, $E_F$ is the Fermi level relative to the host VBM, $\Delta V$ is a potential-alignment correction, and $E_\mathrm{corr}^q$ is the image-charge correction (§6.5.4). For our neutral elemental case the last three "$q$" and "$E_\mathrm{corr}$" terms vanish and we recover (6.3).

#### Term-by-term reading of (6.4)

- **$E_\mathrm{tot}[D^q] - E_\mathrm{tot}[\mathrm{bulk}]$**: the raw DFT energy difference between the supercell with the defect (in charge state $q$) and the same-size supercell of perfect bulk. This is what your two `pw.x` runs give you.
- **$-\sum_i n_i \mu_i$**: chemical-potential bookkeeping. Each atom removed from the supercell goes back into a *reservoir* at chemical potential $\mu_i$. The reservoir is whatever phase the atom equilibrates with — bulk Si for a Si vacancy in pure silicon, or one of several phases for an impurity in a multi-component compound. The sign convention is that $n_i > 0$ means *adding* species $i$ to the perfect cell.
- **$q(E_F + E_\mathrm{VBM})$**: an electron reservoir term. To make a $q$-charged defect, you have transferred $q$ electrons to or from a reservoir at electrochemical potential $E_F$. Inside DFT we report energies relative to the VBM, so we write $E_F$ as a number between 0 and the gap; $E_\mathrm{VBM}$ is the absolute VBM eigenvalue. Defect formation energies are then *functions* of $E_F$, and the plot of $E_f^q(E_F)$ for several $q$ is what defines the **charge transition levels** $\epsilon(q/q')$ — the crossings.
- **$\Delta V$**: a potential-alignment correction. The reference "zero" of the electrostatic potential differs between the perfect and defective supercells (because the charged defect changes the long-range monopole). $\Delta V$ aligns them by comparing a bulk-like region of the defective cell to the perfect cell; without it, the VBM you wrote down has been moved by an arbitrary constant.
- **$E_\mathrm{corr}^q$**: removes the spurious electrostatic interaction between the periodic images of the charged defect — the **image-charge** correction. The leading term scales as $q^2/(\epsilon_\infty L)$ where $L$ is the supercell linear size and $\epsilon_\infty$ is the host high-frequency dielectric constant.

For the neutral vacancy of Si in pure Si, $n_\mathrm{Si} = -1$, $\mu_\mathrm{Si} = E_\mathrm{tot}[\mathrm{bulk}]/N_\mathrm{bulk}$, $q=0$, and $E_\mathrm{corr}=0$. Equation (6.4) collapses to (6.3).

## 6.5.4 Charged defects — what the formula leaves out

If you create a charged defect ($q \neq 0$), three additional complications arise:

1. **A reservoir for electrons**: the term $q(E_F + E_\mathrm{VBM})$ accounts for the electrons added or removed, taken from a reservoir at chemical potential $E_F$. The formation energy then becomes a function of $E_F$, and you typically plot it from $E_F = 0$ (VBM) to $E_F = E_g$ (CBM).
2. **A compensating background charge**: a charged supercell with PBC has infinite electrostatic energy; QE adds a uniform neutralising jellium background to make it finite, but the resulting energy is *not* the energy of an isolated charged defect. It is the energy of a periodic array of charges in a jellium.
3. **Image-charge correction**: the **Freysoldt-Neugebauer-Van de Walle (FNV)** scheme corrects for the spurious electrostatic interaction between periodic images of the charged defect. Lany-Zunger and Kumagai-Oba are alternative formulations. The correction scales as $q^2/(\epsilon L)$ where $L$ is the supercell linear size and $\epsilon$ is the host dielectric constant.

Implementations of these corrections live in post-processing tools — `sxdefectalign` (Freysoldt), `pylada-defects`, `pyCDT`, `pymatgen.analysis.defects` — and require the dielectric tensor of the host, itself computed by a separate DFPT or finite-difference calculation. The full charged-defect workflow is the subject of a chapter of its own; we mention it here so you know the formula above is only the first step.

#### Freysoldt-Neugebauer-Van de Walle correction: a 1D derivation

To see where the FNV correction comes from, consider the simplest model: a 1D infinite array of point charges $q$ separated by length $L$, embedded in a uniform dielectric with constant $\epsilon$. The Coulomb energy of one image-cell sum, with a neutralising jellium background to make the lattice sum convergent, gives a Madelung-like result
$$E_\mathrm{Madelung}^{1\mathrm{D}} = -\frac{q^2 \ln 2}{\epsilon\,L}.$$
This is the *spurious* interaction we want to remove: in the limit $L\to\infty$ the formation energy should approach the isolated-defect value, but a finite-cell DFT calculation contains the Madelung energy as a finite-size error.

In 3D, the analogous quantity is
$$E_\mathrm{Madelung}^{3\mathrm{D}} = -\frac{q^2\,\alpha_\mathrm{M}}{2\epsilon\,L},$$
with $\alpha_\mathrm{M} \approx 2.84$ the Madelung constant for a simple cubic lattice of point charges in jellium. FNV (Freysoldt, Neugebauer, Van de Walle 2009) refines this by:

1. **Decomposing the defect charge** into a localised Gaussian model charge $\rho_\mathrm{model}(\mathbf{r})$ plus a long-range monopole.
2. **Computing the model's Madelung energy analytically** in the periodic supercell — this is the *image-charge* term proper:
   $$E_\mathrm{img} = E_\mathrm{lat}[\rho_\mathrm{model}] - E_\mathrm{iso}[\rho_\mathrm{model}]$$
   where $E_\mathrm{lat}$ is the periodic-array energy and $E_\mathrm{iso}$ the isolated-charge energy.
3. **Aligning the potential** to handle the fact that the average potential of the periodic monopole differs from that of an isolated charge:
   $$\Delta V = \langle V_\mathrm{defective}^\mathrm{far} - V_\mathrm{perfect}^\mathrm{far}\rangle_\mathrm{cell\,average},$$
   measured in a bulk-like region "far" from the defect.

The total FNV correction is then $E_\mathrm{corr} = E_\mathrm{img} - q\,\Delta V$. For a $4\times 4\times 4$ Si supercell with $q=+1$, FNV typically corrects the raw formation energy by 0.3-0.5 eV — large enough that ignoring it gives qualitatively wrong charge transition levels.

The take-home: even for the *neutral* Si vacancy you should be able to compute $E_f$ by (6.3) alone, but the moment you go to $V^+_\mathrm{Si}$ or $V^-_\mathrm{Si}$ you need the full machinery of (6.4) including FNV. Don't extrapolate from your neutral result to charge transition levels without it.

For the rest of this section we focus on the neutral vacancy, where (6.3) is exact (modulo supercell-size convergence).

**Definition 6.5.1 (Defect formation energy).** Given a perfect host crystal with total energy $E_\mathrm{bulk}$ and chemical potentials $\{\mu_i\}$ for each species, the formation energy of a defect $D$ in charge state $q$ is
$$E_f^q[D] = E_\mathrm{tot}[D^q] - E_\mathrm{tot}[\mathrm{bulk}] - \sum_i n_i \mu_i + q(E_F + E_\mathrm{VBM} + \Delta V) + E_\mathrm{corr}^q,$$
where $n_i$ is the change in species count, $E_F$ the Fermi level relative to host VBM, $\Delta V$ a potential alignment, and $E_\mathrm{corr}^q$ the image-charge correction.

**Example 6.5.2 (Neutral Si vacancy formation energy).**
*Problem.* Compute $E_f[V_\mathrm{Si}^0]$ in a 64-atom supercell.

*Solution.* Build the 8-atom conventional cell × $2\times 2\times 2$ supercell (64 atoms). Run SCF for the perfect cell ($E_\mathrm{bulk}$) and the cell with one atom deleted ($E_\mathrm{defective}$, 63 atoms). The chemical potential is $\mu_\mathrm{Si} = E_\mathrm{bulk}/64$ for bulk Si. Apply (6.3):
$$E_f = E_\mathrm{defective} - (63/64)\,E_\mathrm{bulk}.$$
Allowing the 63 atoms to relax (Jahn-Teller distortion), the result is $E_f \approx 3.3$-$3.6$ eV depending on the exact JT minimum reached.

*Discussion.* Experimental estimates from Si self-diffusion give $E_f \sim 3.6$ eV; PBE in a 64-atom cell is consistent within typical accuracy. The 50-meV finite-size error from the 64-atom cell is the dominant systematic; tightening to a 216-atom cell brings $E_f$ closer to the infinite-cell limit at the cost of $\sim 3\times$ wall time.

## 6.5.5 An ASE workflow for the Si vacancy

The workflow:

1. Build a 64-atom Si supercell (perfect).
2. Build a 63-atom Si supercell with one Si removed (vacancy).
3. Run SCF on each, with consistent parameters.
4. Apply (6.3) to get $E_f$.

A complete script:

```python
"""
si_vacancy.py — Compute the neutral Si vacancy formation energy by DFT.

Uses ASE to build cells, drive Quantum ESPRESSO, and post-process.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ase import Atoms
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile
from ase.optimize import BFGS


PSEUDO_DIR = Path.home() / "pseudo/SSSP_1.3.0_PBE_efficiency"
PSEUDOS = {"Si": "Si.pbe-n-rrkjus_psl.1.0.0.UPF"}


@dataclass
class DefectResult:
    """Container for one cell's energy result."""
    label: str
    natoms: int
    energy_eV: float


def make_calc(workdir: Path, *, kpts: tuple[int, int, int],
              relax: bool = False) -> Espresso:
    """Espresso calculator with our converged Si parameters."""
    profile = EspressoProfile(command="pw.x", pseudo_dir=PSEUDO_DIR)
    input_data: dict = {
        "control": {
            "calculation": "relax" if relax else "scf",
            "verbosity":  "low",
            "tprnfor":    True,
            "tstress":    False,
            "etot_conv_thr": 1.0e-5,    # Ry, for ionic relaxation
            "forc_conv_thr": 1.0e-4,    # Ry/Bohr
        },
        "system": {
            "ecutwfc":      50.0,
            "ecutrho":     400.0,
            "occupations": "fixed",
        },
        "electrons": {
            "conv_thr":    1.0e-9,
            "mixing_beta": 0.4,
        },
        "ions": {"ion_dynamics": "bfgs"} if relax else {},
    }
    return Espresso(
        profile=profile,
        directory=str(workdir),
        pseudopotentials=PSEUDOS,
        input_data=input_data,
        kpts=kpts,
    )


def build_perfect_supercell(a: float = 5.43) -> Atoms:
    """8-atom conventional cell × 2×2×2 = 64-atom supercell."""
    si_conv = bulk("Si", crystalstructure="diamond", a=a, cubic=True)  # 8 atoms
    assert len(si_conv) == 8
    super_cell = si_conv.repeat((2, 2, 2))                              # 64 atoms
    assert len(super_cell) == 64
    return super_cell


def build_vacancy_supercell(perfect: Atoms, remove_index: int = 0) -> Atoms:
    """Delete one atom from a copy of the perfect cell."""
    cell = perfect.copy()
    del cell[remove_index]
    assert len(cell) == len(perfect) - 1
    return cell


def run_scf(atoms: Atoms, label: str, kpts: tuple[int, int, int],
            relax: bool) -> DefectResult:
    """Run one calculation, caching by directory."""
    workdir = Path("defects") / label
    workdir.mkdir(parents=True, exist_ok=True)
    out = workdir / "espresso.pwo"
    if out.exists():
        from ase.io import read
        E = read(out).get_potential_energy()
    else:
        a = atoms.copy()
        a.calc = make_calc(workdir, kpts=kpts, relax=relax)
        E = a.get_potential_energy()
    return DefectResult(label=label, natoms=len(atoms), energy_eV=E)


def formation_energy(defective: DefectResult,
                      perfect: DefectResult) -> float:
    """E_f[V_Si] = E_defective - (N_def/N_perf) * E_perfect."""
    return (defective.energy_eV
            - (defective.natoms / perfect.natoms) * perfect.energy_eV)


def main() -> None:
    # ---- Build cells -------------------------------------------------
    perfect = build_perfect_supercell(a=5.43)
    vacancy = build_vacancy_supercell(perfect, remove_index=0)

    # ---- k-grid for a 2×2×2 conventional supercell -------------------
    # The supercell BZ is 8x smaller than the primitive BZ, so a 2×2×2
    # k-grid here is equivalent in absolute density to 4×4×4 on the
    # 8-atom conventional cell (which we showed converges in §6.3).
    kpts = (2, 2, 2)

    # ---- Run SCFs ----------------------------------------------------
    # For the vacancy, atoms around the void will want to relax.
    # We do a static SCF first (unrelaxed), then a relax run for the
    # true formation energy.
    print("=== Unrelaxed (static) ===")
    perfect_static = run_scf(perfect, "perfect_static", kpts, relax=False)
    vacancy_static = run_scf(vacancy, "vacancy_static", kpts, relax=False)
    Ef_static = formation_energy(vacancy_static, perfect_static)
    print(f"  E[perfect] = {perfect_static.energy_eV:.4f} eV "
          f"({perfect_static.energy_eV/64:.4f} eV/atom)")
    print(f"  E[vacancy, unrelaxed] = {vacancy_static.energy_eV:.4f} eV")
    print(f"  E_f[V_Si, unrelaxed]   = {Ef_static:.3f} eV")

    print("\n=== Relaxed ===")
    perfect_relax = run_scf(perfect, "perfect_relax", kpts, relax=True)
    vacancy_relax = run_scf(vacancy, "vacancy_relax", kpts, relax=True)
    Ef_relax = formation_energy(vacancy_relax, perfect_relax)
    print(f"  E[perfect, relaxed] = {perfect_relax.energy_eV:.4f} eV")
    print(f"  E[vacancy, relaxed] = {vacancy_relax.energy_eV:.4f} eV")
    print(f"  E_f[V_Si, relaxed]   = {Ef_relax:.3f} eV")
    print(f"  Relaxation lowers E_f by {Ef_static - Ef_relax:.3f} eV")


if __name__ == "__main__":
    main()
```

Run:

```bash
python si_vacancy.py
```

Wall time on a 2023 MacBook (M2):

- 64-atom static SCF: ~8 minutes.
- 63-atom static SCF: ~8 minutes (same cost; the vacancy doesn't change the basis size).
- Relaxation: 20-40 BFGS steps, ~30-60 minutes per cell.

Total: about 1.5-2 hours. Run overnight if you are starting fresh.

## 6.5.6 What you should see

Expected results for the neutral Si vacancy at the PBE level in a 64-atom cell:

- $E[\mathrm{perfect}, 64\text{-atom}] \approx -10171.8$ eV (= 64 × $-158.94$ eV/atom).
- $E[\mathrm{vacancy}, 63\text{-atom}, \text{unrelaxed}] \approx -10009.1$ eV.
- $E_f[V_\mathrm{Si}, \text{unrelaxed}] \approx 3.7$ eV.

After relaxation, the four atoms around the vacancy displace inward by ~0.2-0.3 Å (the **Jahn-Teller distortion** of the Si vacancy — its degenerate t$_2$-symmetric dangling bonds break degeneracy by lowering symmetry from $T_d$ to $D_{2d}$). The relaxed formation energy drops to:

- $E_f[V_\mathrm{Si}, \text{relaxed}] \approx 3.3$-3.6 eV depending on exactly which JT distortion converges first.

Compare to experiment: the experimental formation energy of $V_\mathrm{Si}$ is hard to measure directly (you can't make a sample of isolated vacancies) but indirect estimates from self-diffusion give ~3.6 eV. PBE in a 64-atom cell agrees within typical accuracy.

!!! warning "Symmetry-breaking and SCF convergence"
    The Si vacancy has a degenerate electronic ground state at the high-symmetry $T_d$ configuration. SCF starting from the symmetric geometry can lock into a high-symmetry but wrong minimum. To find the true Jahn-Teller distorted minimum: (i) break the symmetry by hand — displace one of the four nearest neighbours by 0.05 Å along (110); (ii) set `nspin = 2` even if you expect a singlet (the vacancy is a famously open-shell defect); (iii) try a spin-polarised starting guess with non-zero `starting_magnetization(1)`. Convergence to the correct ground state is harder than just running the script.

## 6.5.7 Sanity checks

Three things to verify before trusting any defect formation energy.

### Charge neutrality

For the neutral vacancy, the total number of electrons in the defective cell must equal $4 \times 63 = 252$. QE prints this near the start of the output as `Number of electrons`. If it does not match, you have made an input mistake.

### Force convergence

After relaxation, the maximum force on any atom should be below `forc_conv_thr` (we set $10^{-4}$ Ry/Bohr ≈ 0.003 eV/Å). Print it with `atoms.get_forces()` and check.

### Supercell-size convergence

Repeat the calculation with a $3\times3\times3$ supercell (216 atoms). If $E_f$ changes by more than your target tolerance, the 64-atom cell is not converged for *this* defect, and you should use the larger cell. For the neutral Si vacancy, the change from 64 to 216 atoms is typically below 50 meV; for charged vacancies it can be 0.5 eV and dominated by image-charge effects.

#### Worked example: $E_f$ vs. supercell size

A small extension to the `si_vacancy.py` script above sweeps the supercell linear size:

```python
def sweep_supercell_size() -> dict[int, float]:
    """Compute E_f[V_Si] in 2×2×2, 3×3×3, 4×4×4 supercells of the conventional cell."""
    results: dict[int, float] = {}
    for n in (2, 3, 4):
        perfect = bulk("Si", crystalstructure="diamond",
                       a=5.43, cubic=True).repeat((n, n, n))
        vacancy = perfect.copy()
        del vacancy[0]
        natoms = n ** 3 * 8
        # k-grid density is preserved by scaling: 4/n in each direction (rounded up)
        kgrid = max(1, 8 // n)
        kpts = (kgrid, kgrid, kgrid)
        Ep = run_scf(perfect, f"perfect_{n}", kpts, relax=False)
        Ev = run_scf(vacancy, f"vacancy_{n}", kpts, relax=False)
        results[natoms] = formation_energy(Ev, Ep)
        print(f"  N={natoms}  E_f = {results[natoms]:.3f} eV")
    return results
```

For the *unrelaxed* neutral Si vacancy at the PBE level, a typical sweep gives:

| Supercell | $N_\mathrm{atoms}$ | $E_f$ (eV) | $\Delta E_f$ vs. 512 (meV) |
|---|---|---|---|
| $2\times 2\times 2$ | 64 | 3.71 | $-70$ |
| $3\times 3\times 3$ | 216 | 3.78 | $\,\,0$ |
| $4\times 4\times 4$ | 512 | 3.78 | $\,\,0$ |

The 64-atom cell is converged to ~50-70 meV — typically adequate for trends but tight for absolute accuracy. The 216-atom cell is essentially converged. (For *relaxed* energies the picture is similar but absolute values drop by 0.2-0.3 eV due to the Jahn-Teller relaxation; the *convergence behaviour* in cell size is similar.) Plot this and you get a monotonic curve that flattens by $L=3$, the hallmark of a well-localised neutral defect. Charged vacancies show much slower convergence dominated by the FNV image-charge contribution, which is why we left them out.

### Chemical-potential reservoirs — multi-component systems

For elemental Si the chemical-potential bookkeeping is trivial: $\mu_\mathrm{Si}$ is set by bulk Si. For binary or ternary compounds, the chemical potentials are not independent and must be chosen consistently.

Consider GaAs. The formation energy of a Ga vacancy ($n_\mathrm{Ga} = -1$, $n_\mathrm{As} = 0$) is
$$E_f[V_\mathrm{Ga}] = E_\mathrm{tot}[\mathrm{cell\ with\ V_{Ga}}] - E_\mathrm{tot}[\mathrm{GaAs\ bulk}] + \mu_\mathrm{Ga}.$$
The chemical potential $\mu_\mathrm{Ga}$ is *not* simply $E_\mathrm{tot}[\mathrm{Ga\ metal}]/N$ — that would be the value in equilibrium with metallic Ga. In a GaAs sample equilibrated with some other phase (As-rich vapour, intermediate compound), $\mu_\mathrm{Ga}$ can take a range of values.

The thermodynamic bounds are:

1. **Ga-rich condition**: $\mu_\mathrm{Ga} = \mu_\mathrm{Ga}^\mathrm{bulk\ metal}$. Then the formation energy of GaAs is $\Delta H_f^\mathrm{GaAs} = \mu_\mathrm{GaAs}^\mathrm{bulk} - \mu_\mathrm{Ga}^\mathrm{bulk} - \mu_\mathrm{As}^\mathrm{bulk}$, so $\mu_\mathrm{As} = \mu_\mathrm{GaAs}^\mathrm{bulk} - \mu_\mathrm{Ga}^\mathrm{bulk\ metal} = \mu_\mathrm{As}^\mathrm{bulk} + \Delta H_f^\mathrm{GaAs}$.
2. **As-rich condition**: $\mu_\mathrm{As} = \mu_\mathrm{As}^\mathrm{bulk}$. Then $\mu_\mathrm{Ga} = \mu_\mathrm{Ga}^\mathrm{bulk\ metal} + \Delta H_f^\mathrm{GaAs}$.

The two conditions bracket the range of physical chemical potentials. For GaAs at growth temperatures, the actual condition is set by the partial pressure of As$_4$ in the growth chamber; reporting "Ga-rich" and "As-rich" formation energies gives the band you can expect across experimental conditions.

For a *defect formation energy phase diagram* — the formation energies of all defects as a function of $\mu$ and $E_F$ — this analysis multiplies out to $2^{n}$ corner conditions for $n$ chemical components, with intermediate phases (e.g. Ga$_2$O$_3$ if oxygen is present) further restricting the allowed range.

## 6.5.8 From formation energy to materials science

Given $E_f$ you can compute several derived quantities.

### Equilibrium concentration

$$ c(T) = N_\mathrm{sites}\, e^{-E_f / k_B T} $$

with $N_\mathrm{sites} = 5 \times 10^{22}$ cm$^{-3}$ for Si. At the silicon melting point (1687 K), with $E_f = 3.6$ eV:

$$ c \sim 5 \times 10^{22} \times e^{-3.6 / (8.617 \times 10^{-5} \cdot 1687)} \approx 4 \times 10^{11} \text{ cm}^{-3}. $$

Very dilute — but vacancy diffusion in Si nonetheless governs dopant diffusion at growth temperatures.

### Self-diffusion coefficient

Combining $E_f$ with the vacancy migration barrier $E_m$ (computed by NEB between adjacent vacancy configurations) gives the activation energy for self-diffusion:

$$ Q = E_f + E_m. $$

Experimental $Q$ for Si self-diffusion is ~4.7 eV, consistent with $E_f \sim 3.6$ eV plus $E_m \sim 1$ eV from NEB.

### Charge transition levels

For each defect, computing $E_f$ as a function of $E_F$ at each charge state $q$ gives a set of lines $E_f^q(E_F)$. The crossings of these lines are the **charge transition levels** $\epsilon(q/q')$ — the Fermi energies at which the defect prefers to change charge. These levels are what spectroscopy (DLTS, photoluminescence) measures.

### A worked $E_f$ vs $E_F$ plot for charged defects

The standard "defect formation energy diagram" shows $E_f^q(E_F)$ as lines for each charge state $q$, parametrised by $E_F$ between $E_\mathrm{VBM}=0$ and $E_\mathrm{CBM}=E_g$. Crossings between lines define charge transition levels.

```python
def plot_defect_diagram(defects: dict[str, dict[int, float]],
                         E_g: float, save: Path) -> None:
    """Plot E_f^q vs E_F for several defects.

    defects: nested dict {defect_name: {q: E_f_at_EF=0}}
    """
    EF = np.linspace(0, E_g, 200)
    fig, ax = plt.subplots(figsize=(6, 5))
    for name, qdict in defects.items():
        for q, Ef0 in qdict.items():
            line = Ef0 + q * EF
            ax.plot(EF, line, label=f"{name}^{q:+d}")
    ax.set_xlim(0, E_g)
    ax.set_xlabel(r"$E_F - E_\mathrm{VBM}$ (eV)")
    ax.set_ylabel(r"$E_f^q$ (eV)")
    ax.legend()
    fig.savefig(save, dpi=180)
```

Each line is straight: $E_f^q(E_F) = E_f^q(0) + q\,E_F$. The slope is the charge state. The *lower envelope* of all lines at a given $E_F$ is the stable charge state at that Fermi level; the *crossings* of the envelope are the charge transition levels $\epsilon(q/q')$, which determine the defect levels measured by DLTS or photoluminescence.

For the Si vacancy with PBE and FNV corrections, typical results:

| Charge $q$ | $E_f^q(E_F=0)$ (eV) |
|---|---|
| $V^{++}_\mathrm{Si}$ | 3.0 |
| $V^{+}_\mathrm{Si}$ | 3.2 |
| $V^{0}_\mathrm{Si}$ | 3.5 |
| $V^{-}_\mathrm{Si}$ | 3.9 |
| $V^{--}_\mathrm{Si}$ | 4.5 |

Crossings: $\epsilon(++/+)$ at $E_F=0.2$ eV, $\epsilon(+/0)$ at $E_F=0.3$ eV, $\epsilon(0/-)$ at $E_F=0.4$ eV. These are the predicted in-gap defect levels. PBE underestimates the gap (Si gap 0.6 eV instead of 1.17 eV), so the calculated levels span only the first half of the experimental gap — a known systematic error of semilocal DFT for charged defects, partially fixed by HSE06 or DFT+U.

### Section summary

The defect formation energy workflow:

1. **Build a supercell** large enough that the defect does not interact with its periodic images. 64 atoms for neutral Si vacancy; 216-512 for charged defects.
2. **Run SCF on the perfect and defective supercells** with matched parameters (cutoff, k-grid scaled by inverse supercell side).
3. **Apply the formation-energy formula** (6.3) for neutral defects, (6.4) with chemical-potential bookkeeping, Fermi-level term, potential alignment, and image-charge correction for charged defects.
4. **Allow relaxation** to capture Jahn-Teller and similar lattice distortions, which can lower $E_f$ by hundreds of meV.
5. **Check convergence in cell size** by comparing two supercells; for charged defects, apply the FNV correction.

Given $E_f$, you compute equilibrium concentration via $c \propto e^{-E_f/k_BT}$, self-diffusion coefficients via $E_f + E_m$ (combined with NEB-computed migration barriers), and charge transition levels via $E_f^q$-vs-$E_F$ crossings.

## 6.5.9 Where this scales up

Real defect studies look at dozens of defects (vacancies, antisites, interstitials, common impurity substitutions) in several charge states, all in supercells of 200-500 atoms. Each calculation is hours; the full study is weeks of compute. High-throughput defect databases — the [DEFAP database](https://www.materialscloud.org), [pymatgen-analysis-defects](https://github.com/materialsproject/pymatgen-analysis-defects), and Open Quantum Materials Database extensions — automate the workflow we sketched above and apply it to thousands of host materials.

For our purposes you have built the foundation: a script that constructs a supercell, runs DFT on perfect and defective copies, applies the formation-energy formula, and respects relaxation. Switching from silicon to your own host requires changing two lines (the species name and the pseudopotential file). The physics — the periodic-image issue, the image-charge correction for charged defects, the chemical-potential bookkeeping for multi-component hosts — is the same everywhere.

**Remark 6.5.3 (PBE and the gap problem for charged defects).** Charged-defect formation energies depend on the host's band gap, which PBE underestimates by 50% for sp-bonded semiconductors. As a result, PBE charge transition levels are squeezed into a smaller gap than they should occupy. Two practical workarounds: (i) hybrid functionals (HSE06) at $\sim 5\times$ PBE cost, which restore most of the gap; (ii) GW corrections to the perfect-host band structure followed by PBE for the defective cells. Both fixes are routine in production but beyond the scope of this introductory chapter.

In the next chapter we make the atoms move: ab initio molecular dynamics, where DFT forces propagate the system through time, and from which we can compute diffusion coefficients, phonons, and finite-temperature free energies.
