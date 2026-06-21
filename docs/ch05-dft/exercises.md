# 5.7 Exercises

Exercises organised in five levels (A–E), with worked solutions inline. Difficulty levels:

- **(★)** routine — checks understanding.
- **(★★)** moderate — requires non-trivial work or coding.
- **(★★★)** challenging — could be the basis of a project.

Answer in your own words, derive every step (where applicable; Levels A, B and E are conceptual/open-ended), and run the code yourself.

!!! tip "How to use these exercises"

    These exercises are now organised into five levels, from gentlest to hardest:

    - **A. Recall** — short vocabulary and definition checks (one line each).
    - **B. Explain in words** — conceptual short answers, no algebra required.
    - **C. Work through the mathematics** — the existing derivations and pen-and-paper problems (Exercises 5.1–5.5).
    - **D. Code and algorithms** — modify and run the SCF code (Exercises 5.6–5.8).
    - **E. Apply and critique** — judgement calls: when to trust DFT, how to diagnose a wrong result.

    If you are studying on your own for the first time, **start at A and B**. They build the vocabulary and intuition you need before the existing ★/★★ problems in C–E will make sense. Do not jump straight to Exercise 5.1 if the words "functional", "self-consistent", or "exchange–correlation" are still hazy.

---

## A. Recall

Short questions to check you know the vocabulary. Each should take one line.

<span class="diff-easy">★ easy</span> **A1.** What variable does DFT use as its fundamental quantity, instead of the many-electron wavefunction?

<span class="diff-easy">★ easy</span> **A2.** What two contributions does the universal functional $F[n]$ contain?

<span class="diff-easy">★ easy</span> **A3.** What is the exchange–correlation energy $E_{xc}$, in words?

<span class="diff-easy">★ easy</span> **A4.** In the Kohn–Sham scheme, what is computed from orbitals rather than directly from the density?

<span class="diff-easy">★ easy</span> **A5.** What quantity does the SCF loop iterate until it stops changing (self-consistency)?

<span class="diff-easy">★ easy</span> **A6.** What does "LDA" stand for, and what is its defining assumption?

??? success "Answers"

    **A1.** The electron density $n(\mathbf r)$, a function of three spatial coordinates — replacing the $3N$-coordinate wavefunction.

    **A2.** The (interacting) kinetic energy and the electron–electron interaction energy. In the Kohn–Sham split it is written as the non-interacting kinetic energy $T_s[n]$, the Hartree energy $U_H[n]$, and the exchange–correlation energy $E_{xc}[n]$.

    **A3.** Everything not captured by the non-interacting kinetic energy and the classical Hartree term: exchange, correlation, and the kinetic-energy correction. It is the small "everything we don't know exactly" term that must be approximated.

    **A4.** The non-interacting kinetic energy $T_s$, evaluated from the Kohn–Sham orbitals $\phi_i$. This is the whole point of the Kohn–Sham construction.

    **A5.** The electron density (equivalently the Kohn–Sham potential $v_{ks}[n]$ built from it). The loop stops when the density that comes out matches the density that went in.

    **A6.** Local Density Approximation. It assumes that at each point the exchange–correlation energy density is that of a uniform electron gas with the local density $n(\mathbf r)$.

---

## B. Explain in words

Conceptual short answers. No algebra — a paragraph each.

<span class="diff-easy">★ easy</span> **B1.** In one paragraph, why did Kohn and Sham choose to compute the kinetic energy from orbitals instead of from an explicit functional of the density?

<span class="diff-easy">★ easy</span> **B2.** Explain why a converged SCF calculation is not necessarily a *correct* result.

<span class="diff-easy">★ easy</span> **B3.** Explain in words why semi-local DFT (LDA/GGA) systematically underestimates band gaps.

<span class="diff-easy">★ easy</span> **B4.** A colleague says "DFT is exact, so any disagreement with experiment must be a bug." In a sentence or two, say what is right and what is wrong about this.

??? success "Hints and answers"

    **B1.** A pure density functional for the kinetic energy (Thomas–Fermi style) is very inaccurate — kinetic energy is large and sensitive to the shell structure of the density, which a local functional misses (you see this quantitatively in Exercise 5.3). By reintroducing a set of single-particle orbitals, the dominant part of the kinetic energy, $T_s$, can be evaluated essentially exactly, leaving only a small unknown remainder ($E_{xc}$) to approximate. The trade is more computational cost (solving for orbitals) in exchange for far higher accuracy.

    **B2.** Convergence only means the loop reached a fixed point: the output density equals the input density for the *chosen* functional. The result is only as good as that functional. A converged LDA number can still be wrong because LDA is approximate (self-interaction, delocalisation, gap errors). Convergence is necessary, not sufficient — see Exercises 5.4 and 5.8.

    **B3.** The Kohn–Sham gap (difference of the LUMO and HOMO eigenvalues) is not the true fundamental gap; they differ by the derivative discontinuity of $E_{xc}$, which semi-local functionals lack. On top of that, self-interaction and delocalisation errors push occupied and unoccupied levels closer together. The result is a systematic underestimate, often by 30–100%.

    **B4.** Right: the *theory* (Hohenberg–Kohn plus Kohn–Sham) is exact in principle. Wrong: in practice we always use an *approximate* exchange–correlation functional, so disagreement with experiment is usually the functional's error, not a code bug — though convergence settings and pseudopotentials can also be at fault.

---

## C. Work through the mathematics

*Level C (work through the mathematics)*

## Exercise 5.1 — Hohenberg–Kohn I in your own words **(★)**

State and prove Hohenberg–Kohn Theorem I without consulting §5.2. Your proof should:

(a) State the precise hypothesis: two non-degenerate ground states from different external potentials yielding the same density.
(b) Apply the variational principle to obtain two strict inequalities.
(c) Add them and derive a contradiction.
(d) State the conclusion: $v_\mathrm{ext}^{(1)} - v_\mathrm{ext}^{(2)} = \mathrm{const}$.

**Solution.** Suppose $v^{(1)}$ and $v^{(2)}$ differ by more than a constant but produce the same ground-state density $n_0$. Their non-degenerate ground states $|\Psi^{(1)}\rangle, |\Psi^{(2)}\rangle$ are distinct: were they equal, applying $\hat{H}^{(1)} - \hat{H}^{(2)} = \hat V^{(1)} - \hat V^{(2)}$ would give $v^{(1)} - v^{(2)} = E^{(1)} - E^{(2)} = \mathrm{const}$, contradicting the hypothesis.

Use $|\Psi^{(2)}\rangle$ as a trial in $\hat{H}^{(1)}$: by the variational principle,

$$
E^{(1)} < \langle\Psi^{(2)}|\hat{H}^{(1)}|\Psi^{(2)}\rangle = E^{(2)} + \int n_0(v^{(1)} - v^{(2)})\,\mathrm d\mathbf r.
$$

Symmetrically, $E^{(2)} < E^{(1)} + \int n_0(v^{(2)} - v^{(1)})\,\mathrm d\mathbf r$. Adding gives $E^{(1)} + E^{(2)} < E^{(1)} + E^{(2)}$ — contradiction. Hence $v^{(1)} - v^{(2)} = \mathrm{const}$. $\blacksquare$

---

*Level C*

## Exercise 5.2 — LDA exchange energy of a uniform density **(★)**

Consider a cubic box of side $L = 10\,a_0$ containing a uniform electron density $n = 0.05\,a_0^{-3}$.

(a) Compute the LDA exchange energy using equation (5.36), the local-density exchange functional $E_x^\mathrm{LDA}[n] = -\tfrac{3}{4}(3/\pi)^{1/3}\int n^{4/3}\,\mathrm d\mathbf r$ (equivalently $E_x^\mathrm{LDA} = -C_x\int n^{4/3}\,\mathrm d\mathbf r$ with $C_x = \tfrac{3}{4}(3/\pi)^{1/3} \approx 0.7386$).
(b) Compute the corresponding Fermi wavevector $k_F$ and Fermi energy, using the uniform-gas relation (5.3), $k_F = (3\pi^2 n)^{1/3}$, and $\varepsilon_F = k_F^2/2$.
(c) The Thomas–Fermi kinetic energy from (5.5), $T_\mathrm{TF}[n] = C_F\int n^{5/3}\,\mathrm d\mathbf r$ with $C_F = \tfrac{3}{10}(3\pi^2)^{2/3} \approx 2.871$.

**Solution.**

(a) From (5.36),
$$
E_x^\mathrm{LDA} = -\tfrac{3}{4}(3/\pi)^{1/3}\int n^{4/3}\,\mathrm d\mathbf r = -\tfrac{3}{4}(3/\pi)^{1/3}\,n^{4/3}\,L^{3}.
$$
With $n = 0.05$, $n^{4/3} = 0.05^{4/3} \approx 0.01842\,a_0^{-4}$. $(3/\pi)^{1/3} \approx 0.9847$. So
$$
E_x^\mathrm{LDA} \approx -0.75 \times 0.9847 \times 0.01842 \times 1000 \approx -13.61\;\mathrm{Ha}.
$$

(b) From (5.3), $k_F = (3\pi^{2}n)^{1/3} = (3\pi^{2}\times 0.05)^{1/3} \approx (1.4804)^{1/3} \approx 1.140\,a_0^{-1}$. Fermi energy $\varepsilon_F = k_F^{2}/2 \approx 0.650\;\mathrm{Ha}$.

(c) From (5.5) with $C_F \approx 2.871$,
$$
T_\mathrm{TF} = C_F\,n^{5/3}\,L^{3} = 2.871\times 0.05^{5/3}\times 1000.
$$
$0.05^{5/3} \approx 9.21\times 10^{-4}$, so $T_\mathrm{TF} \approx 2.871\times 0.921 \approx 2.64\;\mathrm{Ha}$. Note exchange is roughly $-5$ times the kinetic energy at this low density — characteristic of the low-density regime where exchange dominates.

---

*Level C*

## Exercise 5.3 — Constrained search for two electrons in a 1D box **(★★)**

Consider two non-interacting spin-paired electrons in a 1D box $[0,L]$ with hard walls. The single-particle eigenstates are $\phi_n(x) = \sqrt{2/L}\sin(n\pi x/L)$ with $\varepsilon_n = n^{2}\pi^{2}/(2L^{2})$.

(a) Write the ground-state density $n_0(x)$.
(b) Compute the non-interacting kinetic energy $T_s$ directly from the orbitals.
(c) Apply the Thomas–Fermi approximation: $T_\mathrm{TF}[n_0] = C_F\int n_0^{5/3}\mathrm dx$ with the 3D constant $C_F = \tfrac{3}{10}(3\pi^2)^{2/3} \approx 2.871$ (used here as an approximation; the genuine 1D constant differs). Compare $T_\mathrm{TF}$ with the exact $T_s$.

**Solution.**

(a) Both electrons in $\phi_1$: $n_0(x) = 2\,|\phi_1(x)|^{2} = (4/L)\sin^{2}(\pi x/L)$.

(b) $T_s = 2\,\varepsilon_1 = \pi^{2}/L^{2}$. For $L = 5\,a_0$, $T_s = \pi^{2}/25 \approx 0.395\;\mathrm{Ha}$.

(c) Compute numerically. Take $L = 5$.

```python
import numpy as np
L = 5.0
x = np.linspace(1e-6, L - 1e-6, 1000)
n = (4 / L) * np.sin(np.pi * x / L) ** 2
CF = 0.3 * (3 * np.pi ** 2) ** (2 / 3)  # 2.871 in 3D atomic units
T_TF = CF * np.trapezoid(n ** (5 / 3), x)
T_s = np.pi ** 2 / L ** 2
print(T_TF, T_s)
# Output: ~0.249, ~0.395
```

!!! note "NumPy version"

    `np.trapezoid` is the spelling from NumPy 2.0 onwards. On NumPy 1.x the same function is called `np.trapz` (now deprecated); substitute it if you are on an older install.

The Thomas–Fermi value (~0.25 Ha) underestimates the exact non-interacting kinetic energy (~0.40 Ha) by about 35%. The local kinetic functional is insufficient — a known fact, and the precise reason Kohn–Sham theory chose to compute $T_s$ from orbitals rather than as an explicit density functional.

---

*Level C*

## Exercise 5.4 — Self-interaction of one electron **(★★)**

A single electron in a hydrogen-like 1s orbital has density $n(r) = (Z^{3}/\pi)e^{-2Zr}$.

(a) Compute the Hartree self-energy $U_H[n] = \tfrac{1}{2}\iint n(r)n(r')/|r-r'|\,\mathrm d^{3}r\,\mathrm d^{3}r'$.
(b) Compute the LDA exchange energy $E_x^\mathrm{LDA}[n]$ from (5.36), namely $E_x^\mathrm{LDA}[n] = -\tfrac{3}{4}(3/\pi)^{1/3}\int n^{4/3}\,\mathrm d\mathbf r$.
(c) For exact DFT, $E_x[n] + U_H[n] = 0$ for a one-electron system (the self-Coulomb is exactly cancelled by exact exchange). Compute the LDA self-interaction error $U_H + E_x^\mathrm{LDA}$ for $Z = 1$ (hydrogen).

**Solution.**

(a) Standard atomic integral: $U_H = (5/16)Z = 0.3125$ Ha for $Z=1$.

(b) $E_x^\mathrm{LDA} = -\tfrac{3}{4}(3/\pi)^{1/3}\int n^{4/3}\,\mathrm d\mathbf r$. With $n = (Z^{3}/\pi)e^{-2Zr}$,

$$
\int n^{4/3}\,\mathrm d\mathbf r = (Z^{3}/\pi)^{4/3}\cdot 4\pi\int_0^{\infty}r^{2}e^{-8Zr/3}\,\mathrm dr = (Z^{3}/\pi)^{4/3}\cdot 4\pi\cdot\frac{2}{(8Z/3)^{3}} = \frac{27}{64\pi^{1/3}}Z.
$$

So $E_x^\mathrm{LDA} = -(3/4)(3/\pi)^{1/3}\cdot(27/(64\pi^{1/3}))\cdot Z = -(81/(256))\cdot(3/\pi)^{1/3}/\pi^{1/3}\cdot Z$.

Numerically for $Z = 1$: $(3/\pi)^{1/3}/\pi^{1/3} = (3)^{1/3}/\pi^{2/3} \approx 1.4422/2.1450 \approx 0.6725$, so $E_x^\mathrm{LDA} \approx -(81/256)\cdot 0.6725 \approx -0.2128\;\mathrm{Ha}$.

(c) $U_H + E_x^\mathrm{LDA} \approx 0.3125 - 0.2128 \approx +0.100\;\mathrm{Ha}$. *This is the self-interaction error of LDA for one electron in hydrogen*: about 0.1 Ha (~2.7 eV) of spurious self-repulsion that exact DFT would cancel but LDA does not. This is why LDA total energies for hydrogen are too high and why one-electron systems are pathological in LDA.

---

*Level C*

## Exercise 5.5 — Derive the Hartree potential **(★)**

Starting from $U_H[n] = \tfrac{1}{2}\iint n(\mathbf r)n(\mathbf r')/|\mathbf r-\mathbf r'|\,\mathrm d\mathbf r\,\mathrm d\mathbf r'$, derive the Hartree potential $v_H(\mathbf r) = \delta U_H/\delta n(\mathbf r)$.

**Solution.** Vary $n(\mathbf r) \to n(\mathbf r) + \delta n(\mathbf r)$. To first order,

$$
\delta U_H = \tfrac{1}{2}\iint \frac{\delta n(\mathbf r)\,n(\mathbf r') + n(\mathbf r)\,\delta n(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r'.
$$

By symmetry of $|\mathbf r-\mathbf r'|^{-1}$, both terms are equal:

$$
\delta U_H = \iint\frac{\delta n(\mathbf r)\,n(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r' = \int\delta n(\mathbf r)\left[\int\frac{n(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r'\right]\mathrm d\mathbf r.
$$

Reading off,

$$
v_H(\mathbf r) = \frac{\delta U_H}{\delta n(\mathbf r)} = \int\frac{n(\mathbf r')}{|\mathbf r - \mathbf r'|}\,\mathrm d\mathbf r'.
$$

This is the classical electrostatic potential of the charge distribution $n$.

---

## D. Code and algorithms

*Level D (code and algorithms)*

## Exercise 5.6 — Modify the SCF code: change the XC functional **(★★)**

In the SCF code from §5.5, the exchange–correlation potential is purely LDA-X. Modify the code to add a simple LDA correlation in the Perdew–Wang (high-density) form,

$$
\epsilon_c(n) \approx -A\ln(1 + B/r_s),\qquad r_s = (3/(4\pi n))^{1/3},
$$

with $A = 0.0311$ and $B = 1.0$ (very crude). Add the corresponding $v_c$ to `v_ks`. Because $\epsilon_c$ is the energy *per particle*, the correlation potential is the functional derivative of the energy density $n\epsilon_c$:

$$
v_c = \frac{\delta(n\epsilon_c)}{\delta n} = \epsilon_c + n\,\frac{\partial \epsilon_c}{\partial n},
\qquad
\frac{\partial \epsilon_c}{\partial n} = \frac{\partial \epsilon_c}{\partial r_s}\,\frac{\partial r_s}{\partial n}.
$$

For the form above, $\partial\epsilon_c/\partial r_s = -A\,(-B/r_s^2)/(1 + B/r_s)$ and, from $r_s = (3/(4\pi n))^{1/3}$, $\partial r_s/\partial n = -r_s/(3n)$. Substituting these into the chain rule gives the exact expression coded below — derive it yourself rather than copying the code. Re-run for the H$_4$ chain. Report:

(a) The converged total energy.
(b) The number of SCF iterations.
(c) Whether the converged density differs visually from the LDA-X only result.

**Solution sketch.** Add a function:

```python
def lda_correlation_potential(n: NDArray[np.float64]) -> tuple:
    A, B = 0.0311, 1.0
    rs = (3.0 / (4.0 * np.pi * np.maximum(n, 1e-12))) ** (1.0/3.0)
    eps_c = -A * np.log(1.0 + B / rs)
    # v_c = eps_c + n * d(eps_c)/dn
    drs_dn = -rs / (3.0 * np.maximum(n, 1e-12))
    deps_drs = -A * (-B / rs ** 2) / (1.0 + B / rs)
    v_c = eps_c + n * deps_drs * drs_dn
    return v_c, eps_c
```

Add `v_c, _ = lda_correlation_potential(n)` and `v_ks += v_c` inside the loop; update `total_energy()` to include $\int n\epsilon_c\,\mathrm dx$ and subtract $\int n v_c\,\mathrm dx$. Expected: total energy a few percent lower; density qualitatively unchanged; convergence in a similar number of iterations.

---

*Level D (code and algorithms)*

## Exercise 5.7 — Initial guess sensitivity **(★★)**

Re-run the SCF code with three different initial densities:

(a) Uniform: $n^{(0)} = N/L$ (current default).
(b) Superposition of Gaussians centred on the nuclei: $n^{(0)}(x) = \sum_\alpha N_\alpha\,e^{-(x-x_\alpha)^{2}/(2\sigma^{2})}/\sqrt{2\pi}\sigma$ with $\sigma = 1\,a_0$, normalised to $N$ total.
(c) A poor guess: a single Gaussian at $x = 0$ containing all the charge.

For each, report the number of SCF iterations to convergence and the converged total energy. Do they all converge to the same density? Why?

**Solution sketch.** All three should converge to the same density and total energy (the ground state of a convex variational problem with one minimiser, modulo Hohenberg–Kohn). The number of iterations will differ: (b) typically fastest (~6–8 iterations), (a) intermediate (~15), (c) slowest and possibly requires smaller $\alpha$ to avoid divergence. This is the practical reason production codes use atomic-density superpositions as initial guesses.

Caveat: for strongly correlated systems (which our toy LDA-X H chain is *not*) different initial guesses can converge to different local minima — broken-symmetry solutions, different magnetic orderings — and the choice of initial guess becomes a physically meaningful decision.

---

*Level D (code and algorithms)*

## Exercise 5.8 — A Mott-like failure **(★★★)**

In the SCF code, modify the geometry: place two protons close together at $x = 9, 11$ a.u. (an H$_2$ molecule) and *also* set the box length to $L = 40$ a.u. with $n=512$ grid points. (Here $a$ denotes the softening parameter of the 1D soft-Coulomb electron–nucleus potential $v(x) = -Z/\sqrt{(x-x_\alpha)^2 + a^2}$ used in the §5.5 code; it regularises the singularity at the proton and its value shifts the numerical bond length slightly.)

(a) Run with $n_\mathrm{electrons} = 2$. Report the bond length: vary the proton positions and find the energy minimum. Compare with the experimental H$_2$ bond length of 1.40 a.u.

(b) Now run with $n_\mathrm{electrons} = 1$ (one-electron H$_2^{+}$ ion). Stretch the bond to 10 a.u. and look at the converged density. Is the electron localised on one proton (as it should be in the exact ground state at large separation, where the molecule dissociates to H + H$^{+}$) or symmetrically delocalised between both?

(c) For the stretched H$_2^{+}$, compute the total energy with LDA-X and compare with the exact energy of a hydrogen atom (-0.5 Ha). What is the self-interaction error in this configuration?

**Solution sketch.**

(a) LDA-X (with our pedagogical 3D form on 1D) gives an H$_2$ bond length close to the experimental value, perhaps 1.5–1.7 a.u. (the exact number depends on softening parameter $a$). LDA in real 3D famously gives 1.45 a.u., a small overbind.

(b) The converged density will be *symmetric* between the two protons — half an electron on each. This is the famous *fractional-charge failure* of approximate functionals (§5.6). The true ground state at large separation is one full electron on one proton and a bare $\mathrm{H}^{+}$ on the other.

(c) For widely separated H$_2^{+}$, the exact energy is $-0.5$ Ha (one neutral H atom). LDA-X with the symmetric delocalised solution will give a *lower* (more negative) energy due to spurious self-interaction stabilisation — typically 0.1–0.2 Ha too low. This is the *delocalisation error* of semi-local functionals, the underlying cause of many DFT pathologies discussed in §5.6.

The lesson: even a small toy SCF code, faithfully implemented, exhibits the same systematic failure modes as production DFT. Self-interaction error is not a numerical artefact; it is intrinsic to the choice of approximate functional and survives any amount of numerical care.

??? note "Hint"

    The qualitative failure in parts (b) and (c) is the *many-electron* face of the *one-electron* self-interaction error you already quantified analytically in Exercise 5.4. There, for a single electron, $U_H + E_x^\mathrm{LDA} \approx +0.1$ Ha did not cancel as it must for an exact functional. Stretched H$_2^+$ is the same defect in disguise: the lone electron sees a spurious self-repulsion, and LDA-X lowers it by smearing the electron over both protons (fractional charge) instead of localising it. Use the Exercise 5.4 number as your analytic anchor — the per-electron self-interaction it gives is the order of magnitude of the energy error you should expect to see numerically in part (c).

---

## E. Apply and critique

Higher-level questions about choosing DFT, reading its output, and recognising when it fails. These are open-ended — there is no single numerical answer.

!!! note "Level E is open-ended by design"

    Unlike Levels A–D, the Level E items below are open-ended judgement calls. They come with *hints only* — there is deliberately no answer key, because a good answer depends on your system, your constraints, and your reading of the literature. Treat the hints as a checklist of considerations, not a marking scheme.

**E1.** You compute a semiconductor's band gap and get a value substantially smaller than the experimental gap. List three plausible causes and, for each, how you would check it.

??? note "Hint"

    Think about (i) the functional itself — the well-known semi-local gap underestimate and the missing derivative discontinuity (see B3); (ii) the calculation setup — k-point sampling, plane-wave cutoff, pseudopotential, or an under-converged structure; and (iii) what "the gap" means — are you comparing the Kohn–Sham eigenvalue gap to an optical or fundamental gap? To check each: try a hybrid or $GW$ calculation for (i); run convergence tests for (ii); confirm which gap the experiment measured for (iii).

**E2.** Give a situation where you would *not* trust standard semi-local DFT, and name the method you would reach for instead and why.

??? note "Hint"

    Strong correlation (transition-metal oxides, $f$-electron systems), fractional-charge or dissociation problems, and van der Waals–bound systems are classic failure modes. Depending on the case you might reach for DFT+U, a hybrid functional, a dispersion correction, or a higher-rung wavefunction method. Tie your choice to the specific error you are trying to fix (self-interaction, missing dispersion, static correlation).

**E3.** A converged DFT relaxation gives a structure and a set of energies. Before quoting them, what three checks would you run to convince yourself the numbers are meaningful rather than artefacts?

??? note "Hint"

    Convergence with respect to the numerical parameters (k-points, cutoff, smearing); sensible comparison of *energy differences* rather than absolute energies; and a sanity check against a known reference (a related system, experiment, or a different functional). Recognising a failure mode early — as in Exercises 5.4 and 5.8 — is part of the job.

---

## Beyond these exercises

To go further, consider any of the following as small projects:

- Implement Pulay/DIIS for a real 3D system (e.g., helium atom in a spherical box).
- Implement a simple GGA (PBE-style) by extending the code to use $|\nabla n|$ in $v_{xc}$.
- Implement DFT+U: add a Hubbard-$U$ correction term that penalises fractional occupation of an "atomic" orbital localised on each proton site.
- Plot the convergence rates of linear, Anderson, and Pulay mixing on a metallic test system (long uniform 1D chain).
- Add fractional Fermi–Dirac occupations and study the convergence for a half-filled chain (a 1D metal).

These connect the algebra and code of this chapter directly to the practical concerns of Chapter 6.
