# 6.1 Pseudopotentials and basis sets

A solid is, atomistically, a lattice of nuclei surrounded by electrons. To solve the Kohn-Sham equations of [Chapter 5](../ch05-dft/index.md) we need two things: a description of the potential each electron feels from the nuclei (and from itself, mean-field), and a basis in which to expand the Kohn-Sham orbitals. Choosing both is the first practical decision in any DFT calculation, and the rest of the chapter depends on it.

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

## 6.1.2 Three flavours of pseudopotential

Three families dominate modern DFT codes. They differ in how aggressively they smooth the pseudo-wavefunction, and therefore in what cutoff energy you need.

### Norm-conserving (NC)

Introduced by Hamann, Schlüter and Chiang in 1979. The pseudo-wavefunction $\tilde\psi_{n\ell}$ is required to have the same norm inside $r_c$ as the all-electron wavefunction $\psi_{n\ell}^\mathrm{AE}$:

$$ \int_0^{r_c} |\tilde\psi_{n\ell}(r)|^2 r^2\, dr = \int_0^{r_c} |\psi_{n\ell}^\mathrm{AE}(r)|^2 r^2\, dr. $$

This norm condition guarantees correct scattering at the reference energy and, crucially, correct *energy derivative* of the scattering — the pseudopotential is transferable across chemical environments. Modern "ONCV" (optimised norm-conserving Vanderbilt) pseudopotentials are tight, accurate, and the gold standard for high-precision work. They require cutoffs of 50-100 Ry for typical elements.

### Ultrasoft (USPP)

Vanderbilt (1990) relaxed the norm condition. The pseudo-wavefunction can be much smoother — *ultrasoft* — at the price that the orthonormality condition changes from $\langle \tilde\psi_i | \tilde\psi_j \rangle = \delta_{ij}$ to a generalised form $\langle \tilde\psi_i | S | \tilde\psi_j \rangle = \delta_{ij}$ with an overlap operator $S$ that depends on augmentation charges localised near each atom. The eigenvalue problem becomes generalised. USPP cutoffs are typically 25-40 Ry — a factor of 2-4 cheaper than NC.

### Projector augmented wave (PAW)

Blöchl (1994). A linear transformation $|\psi\rangle = |\tilde\psi\rangle + \sum_i (|\phi_i\rangle - |\tilde\phi_i\rangle)\langle \tilde p_i | \tilde\psi\rangle$ takes a smooth pseudo-wavefunction $|\tilde\psi\rangle$ to the true all-electron wavefunction $|\psi\rangle$ by adding back, atom by atom, the difference between all-electron partial waves $\phi_i$ and pseudo partial waves $\tilde\phi_i$. PAW reproduces all-electron results to within a few meV/atom for most properties, is as cheap as USPP, and gives you access to the true wavefunction near the nucleus (useful for hyperfine fields, NMR, EFGs). It is the default in VASP and GPAW, and widely available in QE.

### Which to use

For a starter calculation: **PAW or efficient USPP via SSSP-PBE-efficiency** (see below). Cutoffs around 40-50 Ry, runs in seconds for small cells. For high-precision properties (especially absolute energies and pressures), prefer **norm-conserving** (PseudoDojo) with cutoffs around 80 Ry. For NMR shielding tensors, hyperfine, or anything that probes the wavefunction near the nucleus, use **PAW**.

!!! warning "Mix-and-match pseudopotentials at your peril"
    A pseudopotential is generated together with a particular exchange-correlation functional. A PBE pseudopotential is not interchangeable with an LDA one, and even between PBE pseudopotentials from different libraries the all-electron reference may differ in subtle ways (relativistic treatment, frozen-core boundary, choice of valence configuration). Within one calculation, use pseudopotentials from a *single, consistent set*. Mixing SSSP for one element and PseudoDojo for another is a common bug.

## 6.1.3 Choosing a basis: plane waves vs everything else

The Kohn-Sham orbital must be expanded in some basis $\{\chi_\alpha(\mathbf{r})\}$:

$$ \psi_{n\mathbf{k}}(\mathbf{r}) = \sum_\alpha c_{n\mathbf{k},\alpha}\, \chi_\alpha(\mathbf{r}). $$

Three families are common.

### Plane waves

In a crystal with lattice vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$, Bloch's theorem says

$$ \psi_{n\mathbf{k}}(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}} u_{n\mathbf{k}}(\mathbf{r}), \qquad u_{n\mathbf{k}}(\mathbf{r}+\mathbf{R}) = u_{n\mathbf{k}}(\mathbf{r}), $$

and the cell-periodic part $u_{n\mathbf{k}}$ has a natural Fourier expansion over reciprocal lattice vectors $\mathbf{G}$:

$$ u_{n\mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt{\Omega}} \sum_\mathbf{G} c_{n\mathbf{k},\mathbf{G}}\, e^{i\mathbf{G}\cdot\mathbf{r}}, $$

so

$$ \psi_{n\mathbf{k}}(\mathbf{r}) = \frac{1}{\sqrt{\Omega}} \sum_\mathbf{G} c_{n\mathbf{k},\mathbf{G}}\, e^{i(\mathbf{k}+\mathbf{G})\cdot\mathbf{r}}. $$

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

## 6.1.4 The cutoff energy

The condition

$$ \frac{1}{2}|\mathbf{k} + \mathbf{G}|^2 \leq E_\mathrm{cut} \tag{6.1} $$

selects a finite sphere of $\mathbf{G}$-vectors. The number of plane waves in the basis is

$$ N_\mathrm{PW} \approx \frac{\Omega}{6\pi^2}\, (2 E_\mathrm{cut})^{3/2}, $$

where $\Omega$ is the cell volume in Bohr$^3$. For Si in its 2-atom primitive cell ($\Omega \approx 270\,a_0^3$) at $E_\mathrm{cut} = 40$ Ry, $N_\mathrm{PW} \approx 1100$. Doubling the cutoff multiplies the basis size by $2^{3/2} \approx 2.8$, and the cost of the FFT and the iterative diagonalisation scales similarly.

The right cutoff is whatever your pseudopotential needs. SSSP-PBE-efficiency lists, for each element, a recommended `ecutwfc` (wavefunction cutoff) and a recommended `ecutrho` (density cutoff, used because the density $\rho(\mathbf{r}) = \sum_n |\psi_n(\mathbf{r})|^2$ has Fourier components up to $2|\mathbf{G}|_\mathrm{max}$). For norm-conserving pseudopotentials, $E_\mathrm{cut}^\rho = 4 E_\mathrm{cut}^\psi$ is exact; for USPP and PAW you typically need $E_\mathrm{cut}^\rho = 8$-$12 \times E_\mathrm{cut}^\psi$ to handle augmentation charges.

!!! warning "ecutrho is not optional"
    A common mistake is to converge `ecutwfc` (e.g., 40 Ry) and leave `ecutrho` at its default (4 × ecutwfc, i.e. 160 Ry). For USPP/PAW this is too low — typical values are 320-480 Ry. Symptoms: spurious total-energy oscillations as you change the cell, bad stresses, charge sloshing. Always set `ecutrho` explicitly, following the pseudopotential's recommendation.

For our silicon example using SSSP-PBE-efficiency, the recommended values are `ecutwfc = 30` Ry and `ecutrho = 240` Ry. We will use 40 Ry / 320 Ry to be a little conservative, and verify convergence in [§6.3](03-convergence.md).

## 6.1.5 k-point sampling

The total energy of a crystal is an integral over the Brillouin zone (BZ):

$$ E_\mathrm{tot} = \sum_n \int_\mathrm{BZ} \frac{d^3\mathbf{k}}{(2\pi)^3} f_{n\mathbf{k}}\, \epsilon_{n\mathbf{k}} + \text{(double-counting)}. $$

We can only afford to evaluate the Hamiltonian at a discrete set of $\mathbf{k}$-points and replace the integral by a weighted sum:

$$ \int_\mathrm{BZ} \frac{d^3\mathbf{k}}{(2\pi)^3}\, F(\mathbf{k}) \to \sum_i w_i\, F(\mathbf{k}_i). $$

The question is how to choose the $\mathbf{k}_i$ and $w_i$ to converge the integral fastest.

### Monkhorst-Pack grids

Monkhorst and Pack (1976) gave the standard answer for periodic functions: a uniform grid

$$ \mathbf{k}_{n_1 n_2 n_3} = \sum_{\alpha=1}^3 \frac{2 n_\alpha - N_\alpha - 1}{2 N_\alpha}\, \mathbf{b}_\alpha, \qquad n_\alpha = 1, \ldots, N_\alpha, $$

where $\mathbf{b}_\alpha$ are the reciprocal lattice vectors. We specify the grid by three integers $N_1 \times N_2 \times N_3$. In QE this is `K_POINTS automatic` followed by `N1 N2 N3 s1 s2 s3`, where `s_alpha = 0` or `1` shifts the grid by half a step.

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

## 6.1.6 Pseudopotential libraries

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

## 6.1.7 Summary

You now have the knobs you need:

- **Pseudopotential** — frozen-core, smooth representation of nuclear plus core potential. Choose NC, USPP, or PAW based on accuracy/cost trade-off. SSSP-PBE-efficiency is a safe default.
- **Plane-wave basis** — single number $E_\mathrm{cut}$ controls completeness; chosen large enough for the pseudopotential to be accurate.
- **k-point grid** — Monkhorst-Pack $N_1 \times N_2 \times N_3$, $\Gamma$-centred or shifted; chosen dense enough for the BZ integration tolerance you need.

In the next section we put numbers into these knobs and run a calculation.
