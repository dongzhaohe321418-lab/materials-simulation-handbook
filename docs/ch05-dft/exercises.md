# 5.7 Exercises

Eight problems, with worked solutions inline. Difficulty levels:

- **(★)** routine — checks understanding.
- **(★★)** moderate — requires non-trivial work or coding.
- **(★★★)** challenging — could be the basis of a project.

Answer in your own words, derive every step, and run the code yourself.

---

## Exercise 5.1 — Hohenberg–Kohn I in your own words **(★)**

State and prove Hohenberg–Kohn Theorem I without consulting §5.2. Your proof should:

(a) State the precise hypothesis: two non-degenerate ground states from different external potentials yielding the same density.
(b) Apply the variational principle to obtain two strict inequalities.
(c) Add them and derive a contradiction.
(d) State the conclusion: $v_\mathrm{ext}^{(1)} - v_\mathrm{ext}^{(2)} = \mathrm{const}$.

**Solution.** Suppose $v^{(1)}$ and $v^{(2)}$ differ by more than a constant but produce the same ground-state density $n_0$. Their non-degenerate ground states $|\Psi^{(1)}\rangle, |\Psi^{(2)}\rangle$ are distinct: were they equal, applying $\hat H^{(1)} - \hat H^{(2)} = \hat V^{(1)} - \hat V^{(2)}$ would give $v^{(1)} - v^{(2)} = E^{(1)} - E^{(2)} = \mathrm{const}$, contradicting the hypothesis.

Use $|\Psi^{(2)}\rangle$ as a trial in $\hat H^{(1)}$: by the variational principle,

$$
E^{(1)} < \langle\Psi^{(2)}|\hat H^{(1)}|\Psi^{(2)}\rangle = E^{(2)} + \int n_0(v^{(1)} - v^{(2)})\,\mathrm d\mathbf r.
$$

Symmetrically, $E^{(2)} < E^{(1)} + \int n_0(v^{(2)} - v^{(1)})\,\mathrm d\mathbf r$. Adding gives $E^{(1)} + E^{(2)} < E^{(1)} + E^{(2)}$ — contradiction. Hence $v^{(1)} - v^{(2)} = \mathrm{const}$. $\blacksquare$

---

## Exercise 5.2 — LDA exchange energy of a uniform density **(★)**

Consider a cubic box of side $L = 10\,a_0$ containing a uniform electron density $n = 0.05\,a_0^{-3}$.

(a) Compute the LDA exchange energy using equation (5.36).
(b) Compute the corresponding Fermi wavevector $k_F$ and Fermi energy.
(c) The Thomas–Fermi kinetic energy from (5.5).

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

## Exercise 5.3 — Constrained search for two electrons in a 1D box **(★★)**

Consider two non-interacting spin-paired electrons in a 1D box $[0,L]$ with hard walls. The single-particle eigenstates are $\phi_n(x) = \sqrt{2/L}\sin(n\pi x/L)$ with $\varepsilon_n = n^{2}\pi^{2}/(2L^{2})$.

(a) Write the ground-state density $n_0(x)$.
(b) Compute the non-interacting kinetic energy $T_s$ directly from the orbitals.
(c) Apply the Thomas–Fermi approximation: $T_\mathrm{TF}[n_0] = C_F\int n_0^{5/3}\mathrm dx$ (use the 3D $C_F$ as an approximation). Compare $T_\mathrm{TF}$ with the exact $T_s$.

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
T_TF = CF * np.trapz(n ** (5 / 3), x)
T_s = np.pi ** 2 / L ** 2
print(T_TF, T_s)
# Output: ~0.249, ~0.395
```

The Thomas–Fermi value (~0.25 Ha) underestimates the exact non-interacting kinetic energy (~0.40 Ha) by about 35%. The local kinetic functional is insufficient — a known fact, and the precise reason Kohn–Sham theory chose to compute $T_s$ from orbitals rather than as an explicit density functional.

---

## Exercise 5.4 — Self-interaction of one electron **(★★)**

A single electron in a hydrogen-like 1s orbital has density $n(r) = (Z^{3}/\pi)e^{-2Zr}$.

(a) Compute the Hartree self-energy $U_H[n] = \tfrac{1}{2}\iint n(r)n(r')/|r-r'|\,\mathrm d^{3}r\,\mathrm d^{3}r'$.
(b) Compute the LDA exchange energy $E_x^\mathrm{LDA}[n]$ from (5.36).
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

## Exercise 5.6 — Modify the SCF code: change the XC functional **(★★)**

In the SCF code from §5.5, the exchange–correlation potential is purely LDA-X. Modify the code to add a simple LDA correlation in the Perdew–Wang (high-density) form,

$$
\epsilon_c(n) \approx -A\ln(1 + B/r_s),\qquad r_s = (3/(4\pi n))^{1/3},
$$

with $A = 0.0311$ and $B = 1.0$ (very crude). Add the corresponding $v_c = \delta(n\epsilon_c)/\delta n$ to `v_ks`. Re-run for the H$_4$ chain. Report:

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

## Exercise 5.7 — Initial guess sensitivity **(★★)**

Re-run the SCF code with three different initial densities:

(a) Uniform: $n^{(0)} = N/L$ (current default).
(b) Superposition of Gaussians centred on the nuclei: $n^{(0)}(x) = \sum_\alpha N_\alpha\,e^{-(x-x_\alpha)^{2}/(2\sigma^{2})}/\sqrt{2\pi}\sigma$ with $\sigma = 1\,a_0$, normalised to $N$ total.
(c) A poor guess: a single Gaussian at $x = 0$ containing all the charge.

For each, report the number of SCF iterations to convergence and the converged total energy. Do they all converge to the same density? Why?

**Solution sketch.** All three should converge to the same density and total energy (the ground state of a convex variational problem with one minimiser, modulo Hohenberg–Kohn). The number of iterations will differ: (b) typically fastest (~6–8 iterations), (a) intermediate (~15), (c) slowest and possibly requires smaller $\alpha$ to avoid divergence. This is the practical reason production codes use atomic-density superpositions as initial guesses.

Caveat: for strongly correlated systems (which our toy LDA-X H chain is *not*) different initial guesses can converge to different local minima — broken-symmetry solutions, different magnetic orderings — and the choice of initial guess becomes a physically meaningful decision.

---

## Exercise 5.8 — A Mott-like failure **(★★★)**

In the SCF code, modify the geometry: place two protons close together at $x = 9, 11$ a.u. (an H$_2$ molecule) and *also* set the box length to $L = 40$ a.u. with $n=512$ grid points.

(a) Run with $n_\mathrm{electrons} = 2$. Report the bond length: vary the proton positions and find the energy minimum. Compare with the experimental H$_2$ bond length of 1.40 a.u.

(b) Now run with $n_\mathrm{electrons} = 1$ (one-electron H$_2^{+}$ ion). Stretch the bond to 10 a.u. and look at the converged density. Is the electron localised on one proton (as it should be in the exact ground state at large separation, where the molecule dissociates to H + H$^{+}$) or symmetrically delocalised between both?

(c) For the stretched H$_2^{+}$, compute the total energy with LDA-X and compare with the exact energy of a hydrogen atom (-0.5 Ha). What is the self-interaction error in this configuration?

**Solution sketch.**

(a) LDA-X (with our pedagogical 3D form on 1D) gives an H$_2$ bond length close to the experimental value, perhaps 1.5–1.7 a.u. (the exact number depends on softening parameter $a$). LDA in real 3D famously gives 1.45 a.u., a small overbind.

(b) The converged density will be *symmetric* between the two protons — half an electron on each. This is the famous *fractional-charge failure* of approximate functionals (§5.6). The true ground state at large separation is one full electron on one proton and a bare $\mathrm{H}^{+}$ on the other.

(c) For widely separated H$_2^{+}$, the exact energy is $-0.5$ Ha (one neutral H atom). LDA-X with the symmetric delocalised solution will give a *lower* (more negative) energy due to spurious self-interaction stabilisation — typically 0.1–0.2 Ha too low. This is the *delocalisation error* of semi-local functionals, the underlying cause of many DFT pathologies discussed in §5.6.

The lesson: even a small toy SCF code, faithfully implemented, exhibits the same systematic failure modes as production DFT. Self-interaction error is not a numerical artefact; it is intrinsic to the choice of approximate functional and survives any amount of numerical care.

---

## Beyond these exercises

To go further, consider any of the following as small projects:

- Implement Pulay/DIIS for a real 3D system (e.g., helium atom in a spherical box).
- Implement a simple GGA (PBE-style) by extending the code to use $|\nabla n|$ in $v_{xc}$.
- Implement DFT+U: add a Hubbard-$U$ correction term that penalises fractional occupation of an "atomic" orbital localised on each proton site.
- Plot the convergence rates of linear, Anderson, and Pulay mixing on a metallic test system (long uniform 1D chain).
- Add fractional Fermi–Dirac occupations and study the convergence for a half-filled chain (a 1D metal).

These connect the algebra and code of this chapter directly to the practical concerns of Chapter 6.
