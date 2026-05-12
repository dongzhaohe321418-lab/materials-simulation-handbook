# 6.5 Defects and formation energies

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

$$ E_f^q[D] = E[D, q] - E[\mathrm{perfect}] - \sum_i n_i \mu_i + q (E_F + E_\mathrm{VBM}) + E_\mathrm{corr}^q, $$

where $n_i$ is the change in the number of species $i$ (negative for atoms removed, positive for atoms added), $q$ is the charge state, $E_F$ is the Fermi level relative to the host VBM, and $E_\mathrm{corr}^q$ is the image-charge correction (§6.5.4). For our neutral elemental case the last two terms vanish and we recover (6.3).

## 6.5.4 Charged defects — what the formula leaves out

If you create a charged defect ($q \neq 0$), three additional complications arise:

1. **A reservoir for electrons**: the term $q(E_F + E_\mathrm{VBM})$ accounts for the electrons added or removed, taken from a reservoir at chemical potential $E_F$. The formation energy then becomes a function of $E_F$, and you typically plot it from $E_F = 0$ (VBM) to $E_F = E_g$ (CBM).
2. **A compensating background charge**: a charged supercell with PBC has infinite electrostatic energy; QE adds a uniform neutralising jellium background to make it finite, but the resulting energy is *not* the energy of an isolated charged defect. It is the energy of a periodic array of charges in a jellium.
3. **Image-charge correction**: the **Freysoldt-Neugebauer-Van de Walle (FNV)** scheme corrects for the spurious electrostatic interaction between periodic images of the charged defect. Lany-Zunger and Kumagai-Oba are alternative formulations. The correction scales as $q^2/(\epsilon L)$ where $L$ is the supercell linear size and $\epsilon$ is the host dielectric constant.

Implementations of these corrections live in post-processing tools — `sxdefectalign` (Freysoldt), `pylada-defects`, `pyCDT`, `pymatgen.analysis.defects` — and require the dielectric tensor of the host, itself computed by a separate DFPT or finite-difference calculation. The full charged-defect workflow is the subject of a chapter of its own; we mention it here so you know the formula above is only the first step.

For the rest of this section we focus on the neutral vacancy, where (6.3) is exact (modulo supercell-size convergence).

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

## 6.5.9 Where this scales up

Real defect studies look at dozens of defects (vacancies, antisites, interstitials, common impurity substitutions) in several charge states, all in supercells of 200-500 atoms. Each calculation is hours; the full study is weeks of compute. High-throughput defect databases — the [DEFAP database](https://www.materialscloud.org), [pymatgen-analysis-defects](https://github.com/materialsproject/pymatgen-analysis-defects), and Open Quantum Materials Database extensions — automate the workflow we sketched above and apply it to thousands of host materials.

For our purposes you have built the foundation: a script that constructs a supercell, runs DFT on perfect and defective copies, applies the formation-energy formula, and respects relaxation. Switching from silicon to your own host requires changing two lines (the species name and the pseudopotential file). The physics — the periodic-image issue, the image-charge correction for charged defects, the chemical-potential bookkeeping for multi-component hosts — is the same everywhere.

In the next chapter we make the atoms move: ab initio molecular dynamics, where DFT forces propagate the system through time, and from which we can compute diffusion coefficients, phonons, and finite-temperature free energies.
