# 6.6 Exercises

Difficulty: **(E)** easy / **(M)** medium / **(H)** hard. Hands-on exercises marked **[QE]** assume Quantum ESPRESSO is installed; those marked **[paper]** use toy data and can be done with pen and paper or a notebook.

---

## Exercise 6.1 — Plane-wave basis size **(E, paper)**

Estimate the number of plane waves in the basis for a cubic silicon cell of volume $\Omega = 270\,a_0^3$ (primitive cell, $a = 5.43$ Å) at cutoffs of (a) 30 Ry, (b) 60 Ry, (c) 120 Ry. Comment on the scaling.

**Solution.**

The number of plane waves inside the cutoff sphere is

$$ N_\mathrm{PW} \approx \frac{\Omega}{6\pi^2}(2 E_\mathrm{cut})^{3/2}, $$

with $E_\mathrm{cut}$ in Hartree (so divide a Rydberg value by 2 before substituting).

(a) 30 Ry = 15 Hartree. $N_\mathrm{PW} \approx (270/6\pi^2)(30)^{3/2} = 4.56 \times 164.3 \approx 750$.

(b) 60 Ry = 30 Hartree. $N_\mathrm{PW} \approx 4.56 \times 464.8 \approx 2120$.

(c) 120 Ry = 60 Hartree. $N_\mathrm{PW} \approx 4.56 \times 1315 \approx 6000$.

Doubling the cutoff multiplies the basis by $2^{3/2} \approx 2.83$ — consistent with our cases (750 → 2120 → 6000). The compute cost of iterative diagonalisation scales as $N_\mathrm{PW} \log N_\mathrm{PW}$ per Hamiltonian application and as $N_\mathrm{bnd}^2 N_\mathrm{PW}$ overall; multiplying the cutoff by 4 thus multiplies cost by roughly $4^{3/2} = 8$.

---

## Exercise 6.2 — k-grid density for a slab **(E, paper)**

A graphene unit cell has $|\mathbf{a}_1| = |\mathbf{a}_2| = 2.46$ Å with the third axis $|\mathbf{c}| = 20$ Å (vacuum). Suggest a Monkhorst-Pack grid that achieves a k-point density of at least 30 Å in each periodic direction.

**Solution.**

Aim for $N_i \cdot |\mathbf{a}_i| \geq 30$ Å.

- $N_1 = N_2 \geq \lceil 30/2.46 \rceil = 13$. Round to 14 for an even (shifted) grid, or 13 for Γ-centred.
- $N_3$ is the vacuum direction; one k-point ($N_3 = 1$) is sufficient because the system is periodic in 2D only. Mathematically the dispersion along $k_z$ is flat (no coupling between slabs through the vacuum).

A standard choice: $13 \times 13 \times 1$ Γ-centred for graphene. Because graphene is a semimetal with Dirac cones at K and K', also consider including those points in the grid (use a $12 \times 12 \times 1$ or $18 \times 18 \times 1$ shifted grid).

---

## Exercise 6.3 — Diagnosing a failed SCF **(M, paper)**

A student runs `pw.x` on bcc iron with `ecutwfc = 40`, `occupations = 'fixed'`, `nspin = 1`, and a $4\times4\times4$ k-grid. The SCF oscillates and never converges:

```text
   iteration #  1    total energy =  -329.45 Ry
   iteration #  2    total energy =  -329.21 Ry
   iteration #  3    total energy =  -329.65 Ry
   iteration #  4    total energy =  -329.30 Ry
   ...
```

List four things wrong with the input and explain what to change.

**Solution.**

1. **`occupations = 'fixed'`** is wrong for a metal. Iron has bands crossing the Fermi level, so fixed integer occupations on a 4×4×4 grid will see different bands ordered differently at each k-point and oscillate. Change to `occupations = 'smearing'`, `smearing = 'mv'`, `degauss = 0.01`.
2. **`nspin = 1`** forces a non-magnetic solution. bcc Fe is ferromagnetic at room temperature with a moment of ~2.2 $\mu_B$ per atom; constraining $n_\uparrow = n_\downarrow$ misses this and lands on a higher-energy paramagnetic state with bad SCF. Set `nspin = 2`, `starting_magnetization(1) = 0.5` (positive, fractional initial guess).
3. **k-grid too coarse for a metal**. A $4\times4\times4$ grid samples 8 inequivalent points after symmetry — insufficient for the BCC Fermi surface. Use $12\times12\times12$ or denser.
4. **`ecutwfc = 40` Ry is probably too low** for Fe with the SSSP pseudopotential, which recommends ~50-90 Ry and `ecutrho` of ~500-800 Ry. The student also did not set `ecutrho` and is implicitly getting `4*ecutwfc = 160 Ry`, far too low for USPP/PAW iron.

Fixed input:

```fortran
&SYSTEM
  ibrav = 3,  celldm(1) = 5.42,  nat = 1, ntyp = 1
  ecutwfc     = 90.0
  ecutrho     = 800.0
  occupations = 'smearing'
  smearing    = 'mv'
  degauss     = 0.01
  nspin       = 2
  starting_magnetization(1) = 0.5
/
K_POINTS automatic
  12 12 12   1 1 1
```

---

## Exercise 6.4 — Convergence interpretation **(M, paper)**

A convergence sweep for the total energy per atom of bulk silver gives:

| `ecutwfc` (Ry) | E (eV/atom) |
|---|---|
| 30 | $-2906.4521$ |
| 40 | $-2906.5012$ |
| 50 | $-2906.5184$ |
| 60 | $-2906.5226$ |
| 70 | $-2906.5234$ |
| 80 | $-2906.5236$ |

(a) Plot $|E - E(80)|$ in meV/atom against `ecutwfc` on a log-y axis. (b) What is the smallest cutoff that converges $E$ to within 1 meV/atom? (c) Why might you need a higher cutoff for some other property?

**Solution.**

(a) Differences from the 80 Ry value:

| `ecutwfc` | $\Delta E$ (meV/atom) |
|---|---|
| 30 | 71.5 |
| 40 | 22.4 |
| 50 | 5.2 |
| 60 | 1.0 |
| 70 | 0.2 |
| 80 | 0 (ref) |

On a log-y plot this is a nearly linear decline (exponential convergence).

(b) The 1 meV/atom target is reached at **60 Ry**.

(c) Higher cutoffs may be needed for:

- **Stresses and elastic moduli** — derivatives of energy with respect to strain are more sensitive to basis-set incompleteness.
- **Phonon frequencies** — second derivatives of energy with respect to displacement, also derivative quantities.
- **Magnetic moments and orbital character at the Fermi level** — sensitive to small changes in the wavefunction near the nucleus.
- **Equation-of-state fits** — small fitting errors in $E(V)$ propagate to large errors in $K_0 = -V dP/dV$.

Quote convergence for the property you report, not just for $E$.

---

## Exercise 6.5 — Reproducing the silicon band gap **(M, [QE])**

Run the silicon band-structure pipeline from §6.4 yourself. Report:

(a) The PBE gap and the k-point of the conduction band minimum.

(b) The valence-band width (from the bottom of the s-derived band to the VBM at $\Gamma$).

(c) The energy difference between the heavy-hole top at $\Gamma$ and the conduction-band minimum (the *indirect* gap), and between the top at $\Gamma$ and the conduction band at $\Gamma$ (the *direct* gap at $\Gamma$).

**Solution sketch (expected values for PBE / SSSP-efficiency / 50 Ry / 6×6×6).**

(a) **PBE indirect gap**: $E_g \approx 0.60-0.65$ eV. The conduction band minimum lies along $\Gamma \to X$ at about $0.85 \times \mathbf{k}_X$, with momentum $\approx (0, 0.85, 0)$ in units of $2\pi/a$.

(b) **Valence-band width** (s-band bottom at $\Gamma$ to VBM at $\Gamma$): approximately 12.0 eV. The s-derived band sits around $-12$ to $-9$ eV, and the p-derived bands span $-9$ eV up to the VBM at 0 eV.

(c) **Indirect gap**: $E_g^\mathrm{indirect} \approx 0.63$ eV (VBM at $\Gamma$ to CBM along $\Gamma$-$X$).

**Direct gap at $\Gamma$**: $E_g^\mathrm{direct,\Gamma} \approx 2.5$ eV (VBM to CB at $\Gamma$). Both are underestimates of experiment ($1.17$ eV and $3.4$ eV respectively), but the *ratio* of direct to indirect is close to the experimental ratio. PBE shifts everything down by roughly a constant.

A correct plot also shows the triple degeneracy at the VBM ($\Gamma$, ignoring SOC) and the parabolic conduction-band valley along $\Gamma$-$X$.

---

## Exercise 6.6 — DOS integration check **(E, [QE] or paper from a saved dataset)**

You compute the silicon DOS on a $24\times24\times24$ grid and integrate it from $-\infty$ to the VBM. The result is $7.92$ electrons. Is this correct, what does it tell you about the calculation, and what would you do?

**Solution.**

The 2-atom primitive cell of silicon has $4 \times 2 = 8$ valence electrons. The integrated DOS up to the VBM should equal exactly 8.00 (the count of fully occupied states).

A value of 7.92 means the integration is missing 0.08 electrons (about 1%). The most likely causes:

1. **The energy grid is too coarse**. `DeltaE = 0.01` eV is fine, but if you use `DeltaE = 0.1` eV the trapezoidal integration underestimates the area near peaks.
2. **The k-grid is too coarse**. The DOS itself is built from eigenvalues at the nscf k-points; insufficiently dense k-points miss spectral weight.
3. **Smearing/tetrahedron not active**. If you used Gaussian smearing with too small a $\sigma$, the delta functions on a discrete k-grid leave gaps in the DOS.

Fixes, in order of cost: (i) refine `DeltaE` to 0.005 eV (cheap); (ii) check the VBM energy is correctly identified; (iii) increase k-grid to $32\times32\times32$ or use tetrahedron method (`occupations = 'tetrahedra_opt'`).

The lesson: integrated DOS at the VBM is a free, hard sanity check — always compute it and verify it matches the electron count to better than 0.5%.

---

## Exercise 6.7 — Designing a defect study **(H, paper)**

You are asked to predict whether substitutional nitrogen in diamond (N$_C$) is a deep or shallow donor. Outline the calculations you would perform, including:

(a) Reference calculations needed.

(b) Supercell choice and convergence tests.

(c) Charge states to compute.

(d) Corrections needed for each charge state.

(e) How you would determine "deep vs shallow".

**Solution.**

(a) **Reference calculations.**

1. **Bulk diamond**: total energy, equilibrium lattice constant, band gap, VBM and CBM positions (in absolute energy if you use an electrostatic alignment scheme, or in relative energy to a reference plane).
2. **Bulk nitrogen reference**: total energy per N atom in N$_2$ molecule (the conventional N reservoir), or the $\mu_N$ at chosen growth conditions.
3. **Dielectric constant** of diamond (for image-charge corrections), computed by DFPT (`ph.x` in QE) or finite-difference Berry-phase.

(b) **Supercell convergence.**

- Start with a $3\times3\times3$ conventional diamond supercell = 216 atoms; the C-C bond length is short and you need a few coordination shells between the impurity and its periodic image.
- Test $4\times4\times4$ (512 atoms) for at least one charge state to verify convergence.
- k-grid for the $3\times3\times3$ cell: $2\times2\times2$ shifted (equivalent to $6\times6\times6$ on the primitive cell).
- Plane-wave cutoffs: use the SSSP-recommended values for both C and N, with `ecutrho` set generously.

(c) **Charge states.** For a donor:

- $q = 0$ (neutral N$_C^0$): the dangling electron either sits on a localised level or is delocalised into the CBM.
- $q = +1$ (ionised N$_C^+$): one electron removed, expected to be the stable state when $E_F$ is below the donor level.

Run separate SCF + relaxation for each charge state in the same supercell.

(d) **Corrections** for $q = +1$:

- A neutralising jellium background is added automatically by `pw.x` (`tot_charge = 1.0` in input).
- **Image-charge correction**: Freysoldt-Neugebauer-Van de Walle or Kumagai-Oba. Required because $q^2 \alpha_M / (\epsilon L)$ scales as $\sim 0.3-0.5$ eV at 216 atoms for diamond.
- **Potential alignment**: align the average electrostatic potential of the defective cell to the bulk reference far from the defect.
- **Band-edge correction**: PBE underestimates the diamond gap (predicted ~4.2 eV vs experimental 5.5 eV); to predict the transition level *relative to the experimental CBM*, you might use a hybrid functional (HSE06) or apply a scissor shift.

(e) **Deep vs shallow.**

Compute the **charge transition level** $\epsilon(+/0)$:

$$ \epsilon(+/0) = E_f^0(E_F = 0) - E_f^{+1}(E_F = 0) $$

(i.e. the Fermi level at which $E_f^0 = E_f^{+}$). Then:

- If $\epsilon(+/0)$ sits within ~$k_B T \sim 25$ meV of the CBM → **shallow donor**. The dopant ionises at room temperature, the electron lives in the conduction band, conductivity is high.
- If $\epsilon(+/0)$ lies more than a few hundred meV below the CBM → **deep donor**. The dopant traps its electron, ionisation requires significant thermal activation, the defect is electrically inactive at room temperature.

The experimental answer for N$_C$ in diamond is: **deep donor**, with $\epsilon(+/0) \approx E_\mathrm{CBM} - 1.7$ eV. (This is why nitrogen-doped diamond is brown but not conductive at room temperature.) A correctly executed DFT study with HSE06 reproduces this number to ~0.2 eV.

This last exercise illustrates the gap between "running a DFT calculation" — which the rest of this chapter has taught you — and "doing a defect physics study" — which requires the workflow above for *every* charge state of *every* defect of interest, repeated for several supercell sizes, with electrostatic corrections applied carefully. It is several months of work and the subject of an entire research community. With this chapter behind you, you have the tools to start.
