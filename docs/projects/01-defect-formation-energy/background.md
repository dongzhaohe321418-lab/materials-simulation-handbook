# Background reading — Defect formation energy in silicon

Read these in roughly the order given. For each paper, the annotation
tells you (a) why it matters for *this* project and (b) what you should
extract before moving on. You do **not** need to read every reference
in every paper, but you should be able to write a paragraph in your
own words summarising the contribution.

---

## 1. Zhang and Northrup (1991) — *Chemical potential dependence of defect formation energies in GaAs*

> Zhang, S. B. and Northrup, J. E. *Phys. Rev. Lett.* **67**, 2339 (1991).

This is the paper that formalised the modern way of writing a defect
formation energy as an explicit function of the chemical potentials of
the species added or removed. Although the example is GaAs, the
formalism transfers verbatim to a monovacancy in Si:

$$
E_f(V_\mathrm{Si}^{q}) = E_\mathrm{tot}(V_\mathrm{Si}^{q}) - E_\mathrm{tot}(\text{bulk}) + \mu_\mathrm{Si} + q\,(E_\mathrm{VBM} + E_F) + E_\mathrm{corr}.
$$

**Why read this.** You must understand each term in this expression
before you write a single QE input file. In particular, the role of the
electron chemical potential $E_F$ for charged defects, and the
reservoir argument for $\mu_\mathrm{Si}$, are made explicit here.

**Extract.** A clean derivation of the formula above, the meaning of
"Si-rich" and "Si-poor" conditions (degenerate for an elemental crystal,
trivially $\mu_\mathrm{Si} = E_\mathrm{tot}(\text{bulk Si})/N_\mathrm{atoms}$),
and the conceptual separation between the supercell energies and the
"reservoir" chemical potentials.

---

## 2. Freysoldt, Neugebauer, and Van de Walle (2009) — *Fully ab initio finite-size corrections for charged-defect supercell calculations*

> Freysoldt, C., Neugebauer, J., and Van de Walle, C. G. *Phys. Rev. Lett.*
> **102**, 016402 (2009).

The "FNV" correction is the cleanest method available for removing the
spurious electrostatic interaction between a charged defect and its
periodic images. Although the project as scoped here focuses on the
neutral vacancy (where the FNV term vanishes), the scheme is essential
once you extend to charge states $\pm 1$.

**Why read this.** Even for the neutral defect you should understand
*why* the correction is zero — it is not "because the defect is neutral"
but because the leading multipole moment vanishes. This subtlety
matters when you later move to charged states.

**Extract.** The decomposition of the correction into a Madelung-like
term and a potential-alignment term; the role of the host dielectric
constant ($\varepsilon_\infty \approx 11.7$ for Si); the practical
recipe for evaluating the alignment from the planar-averaged
electrostatic potential.

---

## 3. Makov and Payne (1995) — *Periodic boundary conditions in ab initio calculations*

> Makov, G. and Payne, M. C. *Phys. Rev. B* **51**, 4014 (1995).

The historical precursor to FNV. It gives the asymptotic $L^{-1}$ and
$L^{-3}$ scaling of the image-charge energy and provides the simplest
intuition for why supercells must be extrapolated.

**Why read this.** The $L^{-3}$ term in your finite-size scaling plot
is *exactly* the Makov–Payne quadrupole-monopole interaction in
disguise (even for neutral defects, the elastic image dipole produces
a $L^{-3}$ tail).

**Extract.** The multipole expansion of the spurious interaction; the
distinction between charge-monopole, charge-dipole, and elastic
contributions; the empirical observation that a single Madelung term
captures most — but not all — of the artefact.

---

## 4. Watkins (1986) — *The lattice vacancy in silicon*

> Watkins, G. D. in *Deep Centers in Semiconductors* (ed. Pantelides),
> Gordon & Breach, 1986. Chapter 3.

The definitive experimental review. EPR on the four charge states of
the Si vacancy, the negative-U behaviour, and the Jahn–Teller
distortion — all from the experimentalist who actually measured them.

**Why read this.** Without this paper you do not know what your
calculation is meant to reproduce. The relaxed geometry of the neutral
vacancy is $D_{2d}$ (pair-wise bonding of the four neighbours into
two reconstructed bonds), with an energy gain of roughly 0.4 eV
relative to the $T_d$-symmetric configuration. If your relaxation
does not show this, you have either locked the symmetry by mistake or
not perturbed the initial positions.

**Extract.** The qualitative defect-level structure ($a_1$ singlet
below the gap, $t_2$ triplet near the conduction band edge); the
Jahn–Teller distortion modes; the experimental migration energy
(roughly 0.4 eV for $V^0$).

---

## 5. Puska, Pöykkö, Pesola, and Nieminen (1998) — *Convergence of supercell calculations for point defects in semiconductors: vacancy in silicon*

> Puska, M. J. et al. *Phys. Rev. B* **58**, 1318 (1998).

The original, careful convergence study of the Si vacancy. They
compared 32, 64, 128, and 216-atom cells and demonstrated the
$N^{-1}$-like scaling of $E_f$.

**Why read this.** Your project is, in essence, a modern repeat of
this study. Compare your numbers to theirs once you have them.

**Extract.** Numerical values of $E_f^{V^0}$ at each supercell size;
the residual extrapolation; the role of the **k**-mesh; the
discussion of breathing-mode vs Jahn–Teller relaxation contributions.

---

## 6. Probert and Payne (2003) — *Improving the convergence of defect calculations in supercells: an ab initio study of the neutral silicon vacancy*

> Probert, M. I. J. and Payne, M. C. *Phys. Rev. B* **67**, 075204 (2003).

A more sophisticated convergence analysis arguing that some of the
apparent $N^{-1}$ scaling is in fact slow convergence of the
**k**-mesh, not a true elastic-image effect. They obtain
$E_f^{V^0} \approx 3.6$ eV with denser **k**-sampling on smaller
cells.

**Why read this.** You will see in your own data that the **k**-mesh
choice and the cell size are entangled. This paper teaches the trick
of treating each cell with an equivalent BZ sampling.

**Extract.** The "equivalent-k-mesh" prescription; the magnitude of
**k**-error vs size-error; the recommendation for production runs.

---

## 7. Wright (2006) — *Density-functional-theory calculations for the silicon vacancy*

> Wright, A. F. *Phys. Rev. B* **74**, 165116 (2006).

A systematic study comparing LDA and GGA, with and without spin
polarisation, across charge states. Includes a useful tabulation
that you should compare against.

**Why read this.** It is the cleanest comparison table you will find.
Your single-functional, single-charge-state number must lie within
this scatter.

**Extract.** Table I: $E_f$ values from LDA, PBE, and several other
functionals. The spin state of the neutral vacancy (singlet, $S = 0$).

---

## 8. Corsetti and Mostofi (2011) — *System-size convergence of point-defect properties: the case of the silicon vacancy*

> Corsetti, F. and Mostofi, A. A. *Phys. Rev. B* **84**, 035209 (2011).

A modern, large-scale study going up to 1000 atoms with linear-scaling
DFT. They show that the asymptotic extrapolation is genuinely linear
in $N^{-1}$ all the way to the dilute limit and that the converged
$E_f^{V^0}$ is $3.56 \pm 0.05$ eV in PBE.

**Why read this.** This is the number you should be aiming at. If
your extrapolated value falls in the band 3.5–3.7 eV you are doing
well.

**Extract.** The asymptotic value; the convincing demonstration of
linearity; the discussion of what "converged" means in practice.

---

## 9. Mathew et al. (2016) — *Finite-size errors in continuum electrostatics of point defects: a case study with charged defects in silicon*

> Mathew, K. et al. *Phys. Rev. B* **94**, 044307 (2016).

A practical comparison of correction schemes (FNV, Lany–Zunger,
Kumagai) for various charged defects in Si.

**Why read this.** When you extend to charged states later, this is
your decision tree: which correction to use, when each fails, what
the residual scatter looks like.

**Extract.** The recommendation that FNV is the default choice; the
typical magnitudes of corrections (≈ 0.1–0.3 eV at 216 atoms for
$q = \pm 1$); the alignment-term recipe.

---

## 10. Standard Solid State Pseudopotentials (SSSP) library — Prandini, Marrazzo, Castelli, Mounet, Marzari (2018)

> Prandini, G. et al. *npj Comput. Mater.* **4**, 72 (2018).

The validation paper for the SSSP library. Explains how the
"efficiency" and "precision" sets were curated and how the recommended
cut-offs were derived.

**Why read this.** You need a defensible justification for your
choice of pseudopotential and cut-off energy. SSSP gives both.

**Extract.** The recommended Si pseudopotential for the "efficiency"
set, its $E_\mathrm{cut}$ for wavefunctions and density, and the
phonon and equation-of-state errors quoted for it. Cite these in your
report when justifying your numerical parameters.

---

## Optional eleventh — for context

If you have time, skim Allerdt, Hauschild, and Lany (2020) on
embedded-cluster vs supercell methods (*Phys. Rev. B* **102**, 195203).
It is a useful reminder that the supercell approach is one option
among several; for very deep defects an alternative DMET / embedded
approach is sometimes superior. Not required for the project, but
broadens your perspective.

---

## Synthesis exercise

After you finish reading, write a one-page memo that answers:

1. What is the formal definition of $E_f^{V_\mathrm{Si}^{0}}$?
2. What does the experimental literature say about its value, and
   about the relaxed symmetry of the defect?
3. What is the dominant source of finite-size error for a *neutral*
   defect, and how does it scale with cell size?
4. Which pseudopotential, which functional, and which **k**-mesh
   strategy will you use, and what is your justification?

If you cannot write the memo confidently, re-read papers 1, 5, and 8
before starting calculations.
