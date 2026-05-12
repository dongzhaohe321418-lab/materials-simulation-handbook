# 5.3 The Kohn–Sham Construction

The Hohenberg–Kohn theorems (§5.2) prove that an exact energy functional of the density exists and is variational. They leave us with a problem: the universal functional $F[n] = T[n] + V_{ee}[n]$ is defined by a constrained search and is, in practice, unavailable. Thomas–Fermi tried to write $T[n]$ as a local functional and failed (§5.1). Without a good approximation to the kinetic energy, density functional theory is stuck.

In 1965 Walter Kohn and Lu Jeu Sham, in a paper titled "Self-Consistent Equations Including Exchange and Correlation Effects", broke the impasse with a beautiful manoeuvre. They sidestepped the problem of writing $T[n]$ as an explicit functional of $n$ by reintroducing single-particle orbitals — not as approximations to the wavefunction, but as a *bookkeeping device for computing the kinetic energy exactly* for a fictitious system. The price is that one still has to approximate a small piece — the exchange–correlation energy — but the kinetic energy is no longer the bottleneck.

This is the Kohn–Sham (KS) construction. It is the basis of every practical DFT calculation done today.

## 5.3.1 The auxiliary non-interacting system

Consider a fictitious system of $N$ **non-interacting** electrons moving in a one-body potential $v_s(\mathbf r)$. The Hamiltonian is

$$
\hat H_s = -\tfrac{1}{2}\sum_i \nabla_i^{2} \;+\; \sum_i v_s(\mathbf r_i).
\tag{5.20}
$$

Because the electrons do not interact, the eigenstates of $\hat H_s$ are antisymmetrised products of single-particle orbitals $\phi_i(\mathbf r)$ — Slater determinants. The orbitals satisfy

$$
\Big[-\tfrac{1}{2}\nabla^{2} + v_s(\mathbf r)\Big]\phi_i(\mathbf r) = \varepsilon_i\,\phi_i(\mathbf r).
\tag{5.21}
$$

The ground-state density of this non-interacting system is

$$
n_s(\mathbf r) = \sum_{i=1}^{N}|\phi_i(\mathbf r)|^{2},
\tag{5.22}
$$

where the sum runs over the $N$ orbitals of lowest single-particle energy (with two electrons per spatial orbital in the spin-unpolarised case). The kinetic energy of this non-interacting system is *exactly*

$$
T_s[n] = -\tfrac{1}{2}\sum_{i=1}^{N}\int \phi_i^{*}(\mathbf r)\nabla^{2}\phi_i(\mathbf r)\,\mathrm d\mathbf r,
\tag{5.23}
$$

a number that can be computed once the orbitals are known. There is nothing approximate about (5.23) — it is the exact kinetic energy of the non-interacting system.

The Kohn–Sham postulate is the following: **assume there exists a one-body potential $v_s(\mathbf r)$ such that the ground-state density of the non-interacting system equals the ground-state density of the real, interacting system.**

$$
n_s(\mathbf r) \;=\; n_0(\mathbf r).
\tag{5.24}
$$

This is called *non-interacting $v$-representability* of $n_0$. For most physical systems it is believed (and in many cases proven) to hold; there are pathological counter-examples but they do not arise in typical electronic structure problems.

Given (5.24), the Kohn–Sham orbitals $\{\phi_i\}$ — eigenstates of (5.21) with a *specific* potential $v_s$ — are mathematical objects whose only physical content is that they reproduce the true density via (5.22). They are not the wavefunction of the real system. They are not approximations to the real wavefunction. They are auxiliary quantities, related to the real system *only* through the density they produce.

## 5.3.2 Decomposing the energy functional

The exact total energy of the interacting system is, by HK,

$$
E[n] = T[n] + V_{ee}[n] + \int v_\mathrm{ext}(\mathbf r)\,n(\mathbf r)\,\mathrm d\mathbf r.
\tag{5.25}
$$

We do not know $T[n]$ or $V_{ee}[n]$. But we can identify two pieces we *do* know how to compute and bundle the unknowns into a small remainder.

**The non-interacting kinetic energy $T_s[n]$.** Given by (5.23) once the KS orbitals are known.

**The Hartree energy $U_H[n]$.** The classical electrostatic self-energy of the density,

$$
U_H[n] = \tfrac{1}{2}\iint \frac{n(\mathbf r)\,n(\mathbf r')}{|\mathbf r - \mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r'.
\tag{5.26}
$$

Both $T_s[n]$ and $U_H[n]$ are well-defined functionals of the density (given the orbitals, in the case of $T_s$). Now define the **exchange–correlation energy** by what is left:

$$
\boxed{\;\;E_{xc}[n] \;\equiv\; \big(T[n] - T_s[n]\big) \;+\; \big(V_{ee}[n] - U_H[n]\big).\;\;}
\tag{5.27}
$$

$E_{xc}$ collects everything we do not know:

- $T[n] - T_s[n]$ is the difference between the true kinetic energy and the kinetic energy of the auxiliary non-interacting system. It is small — typically a few per cent of $T$ — because the interacting and non-interacting systems share the same density, hence the same gross spatial extent of the electrons.
- $V_{ee}[n] - U_H[n]$ contains the non-classical part of the electron–electron interaction: the **exchange** energy, which comes from antisymmetry of the wavefunction (Pauli exclusion), and the **correlation** energy, which encodes the dynamical avoidance of electrons.

By construction, the total energy decomposes exactly as

$$
\boxed{\;\;E[n] \;=\; T_s[n] \;+\; U_H[n] \;+\; E_{xc}[n] \;+\; \int v_\mathrm{ext} n\,\mathrm d\mathbf r.\;\;}
\tag{5.28}
$$

This is *exact*: every approximation in DFT now lives in $E_{xc}[n]$ alone. The kinetic energy is treated exactly for the non-interacting system; the dominant Coulomb energy (classical Hartree) is treated exactly; only the small remainder — typically 1–10% of the total energy in solids — needs to be approximated.

That is the central trick of Kohn–Sham theory.

## 5.3.3 Deriving the Kohn–Sham equations

We now minimise (5.28) over densities subject to the constraint $\int n\,\mathrm d\mathbf r = N$, taking advantage of the parameterisation $n = \sum_i |\phi_i|^{2}$ in terms of orthonormal KS orbitals.

Treat the orbitals $\{\phi_i\}$ as the variational degrees of freedom, with Lagrange multipliers $\varepsilon_i$ enforcing orthonormality $\langle\phi_i|\phi_j\rangle = \delta_{ij}$. The Lagrangian is

$$
\mathcal L = E[n] - \sum_{ij}\varepsilon_{ij}\Big(\langle\phi_i|\phi_j\rangle - \delta_{ij}\Big),
$$

with $n = \sum_i |\phi_i|^{2}$. Varying with respect to $\phi_i^{*}(\mathbf r)$ and using the chain rule $\delta n/\delta\phi_i^{*}(\mathbf r) = \phi_i(\mathbf r)$,

$$
\frac{\delta E}{\delta\phi_i^{*}(\mathbf r)} = \frac{\delta E}{\delta n(\mathbf r)}\,\phi_i(\mathbf r) + \frac{\delta T_s}{\delta\phi_i^{*}(\mathbf r)}.
$$

The kinetic term varies directly:

$$
\frac{\delta T_s}{\delta\phi_i^{*}(\mathbf r)} = -\tfrac{1}{2}\nabla^{2}\phi_i(\mathbf r).
$$

The remaining pieces of $E$ are explicit functionals of $n$, with

$$
\frac{\delta E}{\delta n(\mathbf r)} = v_\mathrm{ext}(\mathbf r) + \frac{\delta U_H}{\delta n(\mathbf r)} + \frac{\delta E_{xc}}{\delta n(\mathbf r)}.
$$

Evaluating each:

$$
\frac{\delta U_H}{\delta n(\mathbf r)} = \int\frac{n(\mathbf r')}{|\mathbf r - \mathbf r'|}\,\mathrm d\mathbf r' \;\equiv\; v_H(\mathbf r),
\qquad
\frac{\delta E_{xc}}{\delta n(\mathbf r)} \;\equiv\; v_{xc}(\mathbf r).
$$

Putting everything together and diagonalising $\varepsilon_{ij}$ into its eigenbasis $\varepsilon_i \delta_{ij}$, the variational condition $\delta\mathcal L/\delta\phi_i^{*} = 0$ gives

$$
\boxed{\;\;\Big[-\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}(\mathbf r)\Big]\phi_i(\mathbf r) \;=\; \varepsilon_i\,\phi_i(\mathbf r),\;\;}
\tag{5.29}
$$

with the **Kohn–Sham potential**

$$
\boxed{\;\;v_\mathrm{KS}[n](\mathbf r) \;=\; v_\mathrm{ext}(\mathbf r) \;+\; v_H[n](\mathbf r) \;+\; v_{xc}[n](\mathbf r),\;\;}
\tag{5.30}
$$

and the density

$$
n(\mathbf r) = \sum_{i=1}^{N_\mathrm{occ}} f_i\,|\phi_i(\mathbf r)|^{2}.
\tag{5.31}
$$

The occupations $f_i$ are 1 or 2 (with or without spin) for the lowest $N$ (or $N/2$) eigenstates. Equations (5.29)–(5.31) are the **Kohn–Sham equations**. They look formally like Hartree–Fock equations (Chapter 4), with the non-local Fock exchange operator replaced by the *local, multiplicative* potential $v_{xc}(\mathbf r)$. The non-locality of exchange is hidden inside $v_{xc}$ as a functional of $n$.

Note the self-consistency: $v_\mathrm{KS}$ depends on $n$, and $n$ depends on the orbitals, which are determined by $v_\mathrm{KS}$. The equations are nonlinear and must be solved iteratively — the self-consistent field loop, which is the subject of §5.5.

## 5.3.4 The KS potential as a tool, not a physical field

It is essential to be clear about the status of the various quantities in (5.29).

- **$v_\mathrm{KS}(\mathbf r)$** is the effective one-body potential of the auxiliary non-interacting system. It is *defined* by the requirement that the non-interacting density match the interacting density. It is not "the potential felt by an electron" in any direct physical sense; the real electrons feel each other through Coulomb interactions, not through $v_\mathrm{KS}$.
- **$\phi_i(\mathbf r)$** are eigenstates of $\hat h_s = -\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}$. They are *mathematical objects* whose sole physical content is that their squared moduli sum to the true density.
- **$\varepsilon_i$** are eigenvalues of $\hat h_s$. They are *Lagrange multipliers* for the orthonormality constraint in the variational problem, with one important exception detailed below.

This is a recurring source of misconceptions. We highlight the two most damaging.

!!! warning "KS orbitals are not the real wavefunction"
    There is no physical electron in the orbital $\phi_3$. The interacting many-body wavefunction $\Psi$ is *not* a Slater determinant of $\{\phi_i\}$, and there is no claim that $\Psi \approx \mathrm{det}[\phi_i]$. The KS orbitals provide a way to compute $T_s$ and $n$; that is all they are. Visualising "the highest occupied molecular orbital" as a charge cloud, or interpreting bond patterns from individual $\phi_i$, is using KS orbitals beyond what the theory guarantees. They often *look* reasonable — they tend to resemble Hartree–Fock orbitals — but this is a happy coincidence, not a theorem.

!!! warning "KS eigenvalues are not ionisation energies"
    Hartree–Fock obeys Koopmans' theorem: the negative of the highest occupied orbital energy approximates the first ionisation potential. Kohn–Sham *does not*. The KS eigenvalues are Lagrange multipliers in a variational problem and have no general physical meaning, except as differences (e.g., $\varepsilon_a - \varepsilon_i$ as a crude estimate of an excitation energy, with corrections from TD-DFT). Calculated KS band gaps systematically *underestimate* experimental band gaps by 30–100%, a phenomenon that has nothing to do with the quality of the functional and everything to do with what KS eigenvalues mean (or do not). See §5.6.

### The HOMO exception: Janak's theorem and exact DFT

There is exactly one KS eigenvalue with a guaranteed physical interpretation. In *exact* Kohn–Sham theory (with the exact, unknown $E_{xc}$), the highest occupied KS eigenvalue $\varepsilon_\mathrm{HOMO}$ equals minus the first ionisation potential of the system:

$$
\varepsilon_\mathrm{HOMO} = -I.
\tag{5.32}
$$

This follows from the asymptotic decay of the density: the density of a finite system decays as $n(\mathbf r) \sim e^{-2\sqrt{2I}\,r}$ for large $r$, and since the KS orbitals must reproduce this decay, and the most slowly-decaying occupied orbital governs it, $\varepsilon_\mathrm{HOMO}$ must equal $-I$ exactly. (The derivation is due to Almbladh and von Barth.)

In *approximate* KS theory — i.e., every calculation you will ever do — (5.32) holds only approximately, and badly so for LDA and GGA. Functionals with the correct asymptotic behaviour (some range-separated hybrids, the optimised effective potential method, etc.) do better. The other KS eigenvalues do *not* have such an exact interpretation, even with the exact functional.

### Janak's theorem, more generally

Janak (1978) proved that, in KS theory with fractional occupations $f_i \in [0,1]$,

$$
\frac{\partial E}{\partial f_i} = \varepsilon_i.
\tag{5.33}
$$

This identifies $\varepsilon_i$ as the derivative of the total energy with respect to occupation — a useful identity for some manipulations, but again *not* an ionisation-energy interpretation except for the HOMO at integer occupations.

## 5.3.5 The kinetic energy: $T_s$ versus $T$

A subtle and beautiful aspect of the KS construction is what it does to the kinetic energy.

The exact interacting kinetic energy is

$$
T[n] = -\tfrac{1}{2}\langle\Psi_0[n]|\sum_i\nabla_i^{2}|\Psi_0[n]\rangle,
$$

a functional defined through the *interacting* ground-state wavefunction $\Psi_0[n]$. The KS non-interacting kinetic energy is

$$
T_s[n] = -\tfrac{1}{2}\sum_i\langle\phi_i|\nabla^{2}|\phi_i\rangle,
$$

computed from the KS orbitals of the *non-interacting* auxiliary system that shares the density.

Since both systems have the same density $n$ but $T$ is a minimum over wavefunctions yielding $n$ with the interaction $\hat V_{ee}$ active, while $T_s$ is a minimum over wavefunctions yielding $n$ without interaction, we have

$$
T_s[n] \;\leq\; T[n].
$$

The difference $T_c[n] \equiv T[n] - T_s[n] \geq 0$ is the **correlation kinetic energy** and is part of $E_{xc}$. For most chemical systems $T_c$ is of order 10–50 millihartree per electron, while $T_s$ is of order 1 Hartree per electron: a few per cent correction. By computing $T_s$ exactly via the orbitals, KS theory captures the dominant kinetic energy without approximation; only the small remainder $T_c$ is bundled into $E_{xc}$ and approximated.

This is precisely the failure mode that doomed Thomas–Fermi (§5.1). TF tried to write the entire $T$ as an explicit functional of $n$ — $C_F\int n^{5/3}$ — and got the magnitude roughly right but the spatial dependence so wrong that no molecule could bind. KS bypasses the problem by introducing orbitals to compute $T_s$ directly, leaving $T_c$ (which is small and smooth) to be approximated as a functional of $n$. The orbitals are the price; the prize is chemistry.

!!! note "KS is *not* a free lunch"
    The KS construction costs us the linear scaling with system size that orbital-free DFT (Thomas–Fermi-like methods) enjoys. Diagonalising the KS Hamiltonian scales formally as $\mathcal O(N^{3})$ in the number of electrons; the Hartree term as $\mathcal O(N^{2})$ (or $\mathcal O(N\log N)$ with FFT). Modern algorithms reduce these to near-linear for sparse systems, but the constant in front is large compared with orbital-free methods. The trade is worth it because the result is *chemically accurate*; we discuss this scaling in detail in Chapter 6.

## 5.3.6 Spin-Kohn–Sham

For systems with magnetic order, electron pairing in open shells, or any source of net spin polarisation, one promotes the density to a spin density $(n_\uparrow, n_\downarrow)$ and writes separate KS equations for each spin channel:

$$
\Big[-\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}^{\sigma}[n_\uparrow,n_\downarrow](\mathbf r)\Big]\phi_{i\sigma}(\mathbf r) = \varepsilon_{i\sigma}\,\phi_{i\sigma}(\mathbf r),
$$

$$
n_\sigma(\mathbf r) = \sum_{i\;\mathrm{occ}}|\phi_{i\sigma}(\mathbf r)|^{2},
\qquad
n = n_\uparrow + n_\downarrow.
$$

The exchange–correlation potential is now spin-dependent: $v_{xc}^{\sigma} = \delta E_{xc}[n_\uparrow,n_\downarrow]/\delta n_\sigma$. This is **spin-polarised DFT** or **collinear-spin DFT**. For non-collinear magnetism (spin-orbit coupling, spin spirals), the density becomes a $2\times 2$ matrix in spin space and the KS Hamiltonian is correspondingly generalised.

In what follows, unless stated otherwise, we work with the spin-unpolarised theory; the generalisation is straightforward but notationally heavier.

## 5.3.7 What we have, and what remains

The Kohn–Sham construction gives us an algorithm:

1. Guess an initial density $n^{(0)}$.
2. Build the KS potential $v_\mathrm{KS}[n^{(0)}](\mathbf r) = v_\mathrm{ext} + v_H[n^{(0)}] + v_{xc}[n^{(0)}]$.
3. Solve the eigenvalue problem (5.29) for the orbitals $\{\phi_i\}$.
4. Build a new density $n^{(1)} = \sum_i|\phi_i|^{2}$.
5. Mix old and new densities and iterate to self-consistency.

What is still missing? Two things.

**The exchange–correlation functional $E_{xc}[n]$.** We have defined it by (5.27) but written nothing explicit. Every approximate $E_{xc}$ defines a different "flavour" of DFT — LDA, GGA, hybrid, meta-GGA, and so on. Choosing one is a science in itself, and is the subject of §5.4.

**The self-consistent iteration.** Naive fixed-point iteration on $n$ generally fails to converge — it oscillates between charge-rich and charge-poor solutions. Robust mixing schemes (Pulay/DIIS, Anderson, Broyden) are required, and we work through them in §5.5 with a complete Python implementation.

The KS framework gives chemistry, materials science, and condensed-matter physics a working electronic-structure theory at the cost of solving $N$ coupled non-interacting Schrödinger equations self-consistently. That is a stunning achievement of the 1960s; the next sections build on it.
