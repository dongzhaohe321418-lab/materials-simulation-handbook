# 5.3 The Kohn–Sham Construction

**What problem are we solving?** The Hohenberg–Kohn theorems (Section 5.2) guarantee that the ground-state energy *is* a functional of the density $n(\mathbf r)$ — but the hardest piece of that functional is the kinetic energy, and nobody knows how to write it accurately in terms of $n$ alone. Thomas–Fermi (Section 5.1) tried a purely local kinetic functional and it was too crude to bind a single molecule. Kohn and Sham's trick is to invent an *auxiliary* system of **non-interacting** electrons, deliberately chosen so that it has the *same density* as the real interacting system. For non-interacting electrons we have orbitals, and from orbitals we can compute almost all of the kinetic energy *exactly* — which is precisely the part that defeated every earlier attempt.

!!! note "In plain language"
    We swap one hard problem — many electrons all pushing on each other — for a set of easy ones: a handful of single-particle equations, each like an ordinary one-electron Schrödinger equation, solved over and over until they agree with themselves (self-consistency). The catch is that everything we *cannot* compute exactly gets swept into a single leftover term, the **exchange–correlation energy** $E_{xc}$. That one term is small, but it is where all the approximation lives.

!!! note "Symbol guide"
    | Symbol | Meaning | Units |
    |---|---|---|
    | $\psi_i$ (also $\phi_i$) | Kohn–Sham orbital $i$ of the auxiliary non-interacting system | $a_0^{-3/2}$ |
    | $\varepsilon_i$ | Kohn–Sham eigenvalue (Lagrange multiplier) of orbital $i$ | hartree (Ha) |
    | $n(\mathbf r)$ | electron density (same for auxiliary and real systems by design) | $a_0^{-3}$ |
    | $v_\mathrm{KS}$ | Kohn–Sham (effective) potential of the auxiliary system | hartree |
    | $v_H$ | Hartree potential — classical electrostatics of $n$ | hartree |
    | $v_{xc}$ | exchange–correlation potential, $\delta E_{xc}/\delta n$ | hartree |
    | $T_s$ | kinetic energy of the non-interacting auxiliary system (computed exactly from orbitals) | hartree |

!!! note "Why does this chapter exist?"
    HK said: there exists an exact density functional. Wonderful — but they did not tell us *what it is*. We have the existence proof and no recipe. Without a recipe, DFT is a beautiful theorem with no calculations.
    
    Kohn and Sham's idea is a magic trick. They said: forget about writing the energy directly in terms of the density. Instead, *pretend the electrons do not interact*, and pretend there is an effective potential $v_\mathrm{KS}(\mathbf r)$ — clever, fictitious, to be determined — such that the *density* of the non-interacting electrons in this potential equals the *density* of the real interacting electrons. If you can find such a $v_\mathrm{KS}$, you have a hugely simpler problem: $N$ independent one-electron Schrödinger equations. You solve those, get the density, compute the energy.
    
    A useful analogy. Suppose you want to predict where a swarm of $N$ bees will fly, taking into account their pheromone interactions. Hard problem. But suppose you can find a "mean pheromone field" — an effective potential — such that *non-interacting* bees flying in this field would form exactly the same swarm density as the real interacting bees. Then you can study one bee at a time. The trick is choosing the effective field. The Kohn–Sham construction tells us how.
    
    The price is one small unknown: the *exchange–correlation energy* $E_{xc}[n]$, which captures the difference between the real interacting system and our fictitious non-interacting one. Everything we cannot compute exactly is pushed into $E_{xc}$ — and $E_{xc}$ is small (typically $\sim 10\%$ of the total energy). The whole industry of DFT (§5.4) is the search for accurate approximations to $E_{xc}$.

The Hohenberg–Kohn theorems (§5.2) prove that an exact energy functional of the density exists and is variational. They leave us with a problem: the universal functional $F[n] = T[n] + V_{ee}[n]$ is defined by a constrained search and is, in practice, unavailable. Thomas–Fermi tried to write $T[n]$ as a local functional and failed (§5.1). Without a good approximation to the kinetic energy, density functional theory is stuck.

In 1965 Walter Kohn and Lu Jeu Sham, in a paper titled "Self-Consistent Equations Including Exchange and Correlation Effects", broke the impasse with a beautiful manoeuvre. They sidestepped the problem of writing $T[n]$ as an explicit functional of $n$ by reintroducing single-particle orbitals — not as approximations to the wavefunction, but as a *bookkeeping device for computing the kinetic energy exactly* for a fictitious system. The price is that one still has to approximate a small piece — the exchange–correlation energy — but the kinetic energy is no longer the bottleneck.

This is the Kohn–Sham (KS) construction. It is the basis of every practical DFT calculation done today.

!!! note "Why this step?"
    The intuition behind the KS construction is that we *partially* abandon the goal of a pure density theory in exchange for tractability. HK proved that an exact $E[n]$ exists, but the kinetic-energy piece $T[n]$ defined by a constrained search over wavefunctions is, in practice, no easier than the original Schrödinger equation. Kohn and Sham gave up on writing $T[n]$ explicitly *as a functional of $n$* and instead introduced an auxiliary set of orbitals from which $T$ (well, an approximation to it called $T_s$) can be computed directly. This is *not* a step back from HK; it is a strategic concession that lets us compute the *biggest* piece of the energy exactly, while approximating only a *small* remainder. The trade is one of the most successful in theoretical physics.

!!! abstract "Key idea (Chapter 5.3)"
    The Kohn–Sham construction introduces a fictitious non-interacting system that, by design, has the same density as the real interacting system. Solving $N$ one-electron equations $[-\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}]\phi_i = \varepsilon_i\phi_i$ gives the density $n = \sum_i|\phi_i|^{2}$ and the *non-interacting* kinetic energy $T_s = -\tfrac{1}{2}\sum_i\langle\phi_i|\nabla^{2}|\phi_i\rangle$. The remaining unknown, the *exchange–correlation* energy $E_{xc}[n] = (T-T_s) + (V_{ee}-U_H)$, is small ($\sim 10\%$ of total) and must be approximated. This single move turns DFT from an existence theorem (§5.2) into an algorithm.

## 5.3.1 The auxiliary non-interacting system

Consider a fictitious system of $N$ **non-interacting** electrons moving in a one-body potential $v_s(\mathbf r)$. The Hamiltonian is

$$
\hat{H}_s = -\tfrac{1}{2}\sum_i \nabla_i^{2} \;+\; \sum_i v_s(\mathbf r_i).
\tag{5.20}
$$

Because the electrons do not interact, the eigenstates of $\hat{H}_s$ are antisymmetrised products of single-particle orbitals $\phi_i(\mathbf r)$ — Slater determinants. The orbitals satisfy

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

!!! note "Why this step?"
    Why is it reasonable to expect such a $v_s$ to exist? Intuitively, the interacting density is a smooth, positive function with appropriate decay at infinity. The non-interacting (Slater-determinant) wavefunction has enough freedom — through the choice of $N$ orbitals — to reproduce *any* such density (this is the content of Harriman's construction mentioned in §5.2.4a). The non-trivial step is that the orbitals so constructed should *also* be eigenstates of a *single* local potential. For non-degenerate ground states of typical many-body systems, the existence and uniqueness of such a potential is the subject of considerable mathematical work; for the systems treated in materials science it is taken for granted. The pathological counter-examples (Levy, 1982; Lieb, 1983) involve carefully constructed densities that are *not* the ground-state density of any local potential — but they do not arise from any physical Hamiltonian.

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

!!! example "How small is $E_{xc}$ for a real system?"
    For a silicon atom, the total energy is $E\approx -288.9\;\text{Ha}$, comprising $T_s\approx +288.9\;\text{Ha}$, $V_\mathrm{ext}\approx -697.0\;\text{Ha}$, $U_H\approx +145.2\;\text{Ha}$, $E_x\approx -24.6\;\text{Ha}$ (LDA), and $E_c\approx -1.4\;\text{Ha}$ (LDA). These balance: $288.9 - 697.0 + 145.2 - 24.6 - 1.4 = -288.9\;\text{Ha}$, and $T_s\approx -E$ as the virial theorem requires for a Coulombic system. Exchange is roughly $8.5\%$ of the total kinetic energy and correlation roughly $0.5\%$. In a typical chemical bonding context, the energies *we care about* (binding energies, barriers) are $\sim 0.1\;\text{Ha}$, which is also the typical *error* in $E_{xc}$. So the bond-energy accuracy of a DFT calculation is governed almost entirely by the quality of the XC approximation — exactly the point that motivates the entire functional-development industry of §5.4.

    !!! note "Numbers are representative"
        These are representative LDA self-consistent all-electron results for the neutral Si atom (14 electrons), rounded so the ledger balances exactly ($E = T_s + V_\mathrm{ext} + U_H + E_x + E_c$, with $T_s\approx -E$ from the virial theorem); precise values depend on the code, basis, and functional variant. The qualitative conclusion — exchange a few per cent of $T_s$, correlation a fraction of a per cent — is what matters.

That is the central trick of Kohn–Sham theory.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. What is the Kohn–Sham postulate — what does the fictitious non-interacting system have in common with the real interacting one?
    2. The exact energy decomposes as $E = T_s + U_H + E_{xc} + \int v_\mathrm{ext}\,n$. Which of these terms is computed exactly, and which one contains every approximation in DFT?
    3. Write down the definition of $E_{xc}[n]$ — which two physical contributions does it bundle together, and why is it a relatively small fraction of the total energy?

    If any of these is shaky, re-read the preceding section before continuing.

## 5.3.3 Deriving the Kohn–Sham equations

We now minimise (5.28) over densities subject to the constraint $\int n\,\mathrm d\mathbf r = N$, taking advantage of the parameterisation $n = \sum_i |\phi_i|^{2}$ in terms of orthonormal KS orbitals. We walk through every step of the derivation, since the result is the operational heart of practical DFT.

**Step (i): identify the variational degrees of freedom.** The density $n$ enters the energy functional (5.28) only through $T_s$, $U_H$, $E_{xc}$, and $\int v_\mathrm{ext}\,n$. Of these, $T_s$ is a functional of the *orbitals* directly (it is the kinetic energy of the non-interacting wavefunction $\det[\phi_i]$), while $U_H$ and $E_{xc}$ depend on $n$ which itself depends on $\{\phi_i\}$ through (5.22). So we may treat $\{\phi_i\}$ as the independent variables.

**Step (ii): impose orthonormality.** The KS orbitals must satisfy $\langle\phi_i|\phi_j\rangle = \delta_{ij}$; otherwise the density $n = \sum_i|\phi_i|^{2}$ over-counts. Introduce Lagrange multipliers $\varepsilon_{ij}$ enforcing each constraint.

**Step (iii): write the Lagrangian.**

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

!!! note "Why this step?"
    The chain rule $\delta E/\delta\phi_i^{*} = (\delta E/\delta n)\,\delta n/\delta\phi_i^{*}$ exploits the fact that $U_H$, $E_{xc}$ and $\int v_\mathrm{ext} n$ depend on $\phi_i^{*}$ *only* through $n$. The kinetic functional $T_s$, in contrast, depends on $\phi_i^{*}$ directly through its gradient and so contributes the $-\tfrac{1}{2}\nabla^{2}\phi_i$ term separately. This split is the technical reason the KS Hamiltonian comes out as a *local* operator plus a kinetic operator — and is also why hybrid functionals, which depend on the *orbitals* through Fock exchange, give a Hamiltonian with a *non-local* exchange operator (the generalised Kohn–Sham scheme).

??? note "Full derivation: from the energy functional to the KS equations"
    The four steps above are each compressed to their result. Here we carry out every algebraic line, in the order (a) chain rule, (b) kinetic variation, (c) Hartree variation, (d) the Lagrange-multiplier constraint and the rotation to canonical form.

    **(a) The functional chain rule and the delta function.** The Hartree, exchange–correlation and external pieces of $E$ depend on $\phi_i^{*}(\mathbf r)$ *only* through the density $n(\mathbf r') = \sum_j |\phi_j(\mathbf r')|^{2}$. For a functional $E[n]$ the chain rule for functional derivatives reads

    $$
    \frac{\delta E}{\delta\phi_i^{*}(\mathbf r)} = \int \frac{\delta E}{\delta n(\mathbf r')}\,\frac{\delta n(\mathbf r')}{\delta\phi_i^{*}(\mathbf r)}\,\mathrm d\mathbf r'.
    \tag{5.29a}
    $$

    The inner factor is obtained by differentiating $n(\mathbf r') = \sum_j \phi_j^{*}(\mathbf r')\phi_j(\mathbf r')$ with respect to $\phi_i^{*}(\mathbf r)$. Only the $j=i$ term survives, and $\delta\phi_i^{*}(\mathbf r')/\delta\phi_i^{*}(\mathbf r) = \delta(\mathbf r-\mathbf r')$, so

    $$
    \frac{\delta n(\mathbf r')}{\delta\phi_i^{*}(\mathbf r)} = \phi_i(\mathbf r')\,\delta(\mathbf r-\mathbf r') = \phi_i(\mathbf r)\,\delta(\mathbf r-\mathbf r').
    \tag{5.29b}
    $$

    Inserting (5.29b) into (5.29a), the delta function collapses the integral, evaluating $\delta E/\delta n$ at $\mathbf r' = \mathbf r$:

    $$
    \frac{\delta E}{\delta\phi_i^{*}(\mathbf r)}\bigg|_{n\text{-part}} = \frac{\delta E}{\delta n(\mathbf r)}\,\phi_i(\mathbf r).
    \tag{5.29c}
    $$

    This is the $(\delta E/\delta n)\,\phi_i$ term quoted in the main text.

    **(b) The kinetic variation.** $T_s$ is *not* a functional of $n$ alone; it depends on the orbitals directly. From (5.23),

    $$
    T_s = -\tfrac12 \sum_j \int \phi_j^{*}(\mathbf r')\,\nabla'^{2}\phi_j(\mathbf r')\,\mathrm d\mathbf r'.
    $$

    Differentiating with respect to $\phi_i^{*}(\mathbf r)$, only the $j=i$ term contributes, and $\delta\phi_i^{*}(\mathbf r')/\delta\phi_i^{*}(\mathbf r) = \delta(\mathbf r-\mathbf r')$ again collapses the integral:

    $$
    \frac{\delta T_s}{\delta\phi_i^{*}(\mathbf r)} = -\tfrac12 \int \delta(\mathbf r-\mathbf r')\,\nabla'^{2}\phi_i(\mathbf r')\,\mathrm d\mathbf r' = -\tfrac12\nabla^{2}\phi_i(\mathbf r).
    \tag{5.29d}
    $$

    **(c) The Hartree variation and the factor-of-2 cancellation.** Write the Hartree energy (5.26) with the symmetric kernel $K(\mathbf r_1,\mathbf r_2) = |\mathbf r_1 - \mathbf r_2|^{-1} = K(\mathbf r_2,\mathbf r_1)$:

    $$
    U_H[n] = \tfrac12 \iint n(\mathbf r_1)\,K(\mathbf r_1,\mathbf r_2)\,n(\mathbf r_2)\,\mathrm d\mathbf r_1\,\mathrm d\mathbf r_2.
    $$

    A first-order change $n \to n + \delta n$ produces two linear-in-$\delta n$ terms, one from varying the first $n$ and one from varying the second:

    $$
    \delta U_H = \tfrac12 \iint \delta n(\mathbf r_1)\,K(\mathbf r_1,\mathbf r_2)\,n(\mathbf r_2)\,\mathrm d\mathbf r_1\,\mathrm d\mathbf r_2
    + \tfrac12 \iint n(\mathbf r_1)\,K(\mathbf r_1,\mathbf r_2)\,\delta n(\mathbf r_2)\,\mathrm d\mathbf r_1\,\mathrm d\mathbf r_2.
    $$

    In the second term relabel $\mathbf r_1 \leftrightarrow \mathbf r_2$ and use $K(\mathbf r_2,\mathbf r_1) = K(\mathbf r_1,\mathbf r_2)$; it becomes *identical* to the first. The two equal halves add, cancelling the $\tfrac12$:

    $$
    \delta U_H = \iint \delta n(\mathbf r_1)\,K(\mathbf r_1,\mathbf r_2)\,n(\mathbf r_2)\,\mathrm d\mathbf r_1\,\mathrm d\mathbf r_2
    = \int \delta n(\mathbf r_1)\left[\int \frac{n(\mathbf r_2)}{|\mathbf r_1-\mathbf r_2|}\,\mathrm d\mathbf r_2\right]\mathrm d\mathbf r_1.
    $$

    Reading off the coefficient of $\delta n$ gives the functional derivative

    $$
    \frac{\delta U_H}{\delta n(\mathbf r)} = \int \frac{n(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r' \equiv v_H(\mathbf r).
    \tag{5.29e}
    $$

    The factor of 2 from the two symmetric terms is exactly what cancels the $\tfrac12$ in the definition of $U_H$ — this is why $v_H$ carries no $\tfrac12$. (By contrast $\delta E_{xc}/\delta n \equiv v_{xc}$ is *defined* as the XC functional derivative; it has no closed form because $E_{xc}$ has none.)

    **(d) The constraint variation and the rotation to canonical form.** Collecting (5.29c), (5.29d), (5.29e) and $\delta(\int v_\mathrm{ext} n)/\delta n = v_\mathrm{ext}$, the variation of the energy is

    $$
    \frac{\delta E}{\delta\phi_i^{*}(\mathbf r)} = \Big[-\tfrac12\nabla^{2} + v_\mathrm{ext}(\mathbf r) + v_H(\mathbf r) + v_{xc}(\mathbf r)\Big]\phi_i(\mathbf r) = \Big[-\tfrac12\nabla^{2} + v_\mathrm{KS}(\mathbf r)\Big]\phi_i(\mathbf r).
    $$

    The orthonormality constraints $\langle\phi_i|\phi_j\rangle - \delta_{ij} = 0$ are enforced by the full multiplier matrix $\varepsilon_{ij}$ (one multiplier per constraint; the constraints are not independent of each other, so in general the matrix is *not* diagonal). Varying the constraint term,

    $$
    \frac{\delta}{\delta\phi_i^{*}(\mathbf r)}\sum_{kl}\varepsilon_{kl}\Big(\langle\phi_k|\phi_l\rangle - \delta_{kl}\Big) = \sum_{l}\varepsilon_{il}\,\phi_l(\mathbf r),
    $$

    since $\delta\langle\phi_k|\phi_l\rangle/\delta\phi_i^{*}(\mathbf r) = \delta_{ki}\,\phi_l(\mathbf r)$. Setting $\delta\mathcal L/\delta\phi_i^{*} = 0$ gives the **off-diagonal Euler–Lagrange equation**

    $$
    \Big[-\tfrac12\nabla^{2} + v_\mathrm{KS}(\mathbf r)\Big]\phi_i(\mathbf r) = \sum_{j}\varepsilon_{ij}\,\phi_j(\mathbf r).
    \tag{5.29f}
    $$

    Here $\varepsilon_{ij}$ is the full Hermitian multiplier *matrix* (Hermiticity follows because $\langle\phi_i|\phi_j\rangle = \langle\phi_j|\phi_i\rangle^{*}$, so the multipliers satisfy $\varepsilon_{ij} = \varepsilon_{ji}^{*}$). It is *not yet* the diagonal $\varepsilon_i\delta_{ij}$ of the canonical equation. To diagonalise it, perform a unitary rotation of the occupied orbitals, $\tilde\phi_i = \sum_j U_{ji}\phi_j$ with $U^{\dagger}U = \mathbb 1$. Because both $-\tfrac12\nabla^{2}+v_\mathrm{KS}$ (which depends only on $n$, and $n = \sum_i|\phi_i|^{2}$ is invariant under any unitary mixing of occupied orbitals) and the constraint structure are unchanged, (5.29f) transforms covariantly into

    $$
    \Big[-\tfrac12\nabla^{2} + v_\mathrm{KS}\Big]\tilde\phi_i = \sum_j (U^{\dagger}\varepsilon\,U)_{ij}\,\tilde\phi_j.
    $$

    Choosing $U$ as the matrix of eigenvectors of the Hermitian matrix $\varepsilon$ makes $U^{\dagger}\varepsilon\,U = \mathrm{diag}(\varepsilon_1,\dots,\varepsilon_N)$. In this **canonical basis** the right-hand side collapses to $\varepsilon_i\tilde\phi_i$, and dropping the tilde recovers the diagonal form (5.29). Thus $\varepsilon_{ij}$ is the full multiplier matrix and $\varepsilon_i$ are its real eigenvalues — the canonical KS eigenvalues. The density and total energy are identical in any choice of $U$; only the canonical basis yields the clean one-electron eigenvalue equation.

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

!!! note "Why this step?"
    The off-diagonal Lagrange multipliers $\varepsilon_{ij}$ form a Hermitian matrix; we can always rotate the orbital basis to diagonalise it without changing the Slater determinant (because $\det[U\phi] = \det U\cdot\det[\phi]$ and a unitary $U$ has $|\det U|=1$). So replacing $\varepsilon_{ij}$ by its diagonal form $\varepsilon_i\delta_{ij}$ costs nothing physically but yields the clean one-electron eigenvalue equation (5.29). The eigenvalues $\varepsilon_i$ are the *canonical* KS eigenvalues. Any unitary rotation of occupied orbitals leaves the density and the total energy unchanged, but the eigenvalues are only defined in the canonical basis.

Note the self-consistency: $v_\mathrm{KS}$ depends on $n$, and $n$ depends on the orbitals, which are determined by $v_\mathrm{KS}$. The equations are nonlinear and must be solved iteratively — the self-consistent field loop, which is the subject of §5.5.

## 5.3.4 The KS potential as a tool, not a physical field

It is essential to be clear about the status of the various quantities in (5.29).

- **$v_\mathrm{KS}(\mathbf r)$** is the effective one-body potential of the auxiliary non-interacting system. It is *defined* by the requirement that the non-interacting density match the interacting density. It is not "the potential felt by an electron" in any direct physical sense; the real electrons feel each other through Coulomb interactions, not through $v_\mathrm{KS}$.
- **$\phi_i(\mathbf r)$** are eigenstates of $\hat h_s = -\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}$. They are *mathematical objects* whose sole physical content is that their squared moduli sum to the true density.
- **$\varepsilon_i$** are eigenvalues of $\hat h_s$. They are *Lagrange multipliers* for the orthonormality constraint in the variational problem, with one important exception detailed below.

This is a recurring source of misconceptions. We highlight the two most damaging.

!!! warning "Common misunderstanding"
    The Kohn–Sham orbitals $\psi_i$ and eigenvalues $\varepsilon_i$ belong to the *fictitious* non-interacting system. They are **not** the real electrons, and the $\varepsilon_i$ are **not** the true energy levels (the energies you would measure when adding or removing an electron). Treating $-\varepsilon_i$ as an ionisation energy, or reading $\varepsilon_a - \varepsilon_i$ off as an exact excitation energy, is reading more into the auxiliary system than the theory promises.

    There are two important partial exceptions. (i) In *exact* DFT the highest occupied eigenvalue does have meaning: $\varepsilon_\mathrm{HOMO} = -I$, the negative of the first ionisation energy (see (5.32) below). (ii) Kohn–Sham gaps are not meaningless, but they *systematically underestimate* true fundamental gaps — a structural feature of the theory, not just a flaw of any one functional (Section 5.6).

!!! warning "KS orbitals are not the real wavefunction"
    There is no physical electron in the orbital $\phi_3$. The interacting many-body wavefunction $\Psi$ is *not* a Slater determinant of $\{\phi_i\}$, and there is no claim that $\Psi \approx \mathrm{det}[\phi_i]$. The KS orbitals provide a way to compute $T_s$ and $n$; that is all they are. Visualising "the highest occupied molecular orbital" as a charge cloud, or interpreting bond patterns from individual $\phi_i$, is using KS orbitals beyond what the theory guarantees. They often *look* reasonable — they tend to resemble Hartree–Fock orbitals — but this is a happy coincidence, not a theorem.

!!! warning "KS eigenvalues are not ionisation energies"
    Hartree–Fock obeys Koopmans' theorem: the negative of the highest occupied orbital energy approximates the first ionisation potential. Kohn–Sham *does not*. The KS eigenvalues are Lagrange multipliers in a variational problem and have no general physical meaning, except as differences (e.g., $\varepsilon_a - \varepsilon_i$ as a crude estimate of an excitation energy, with corrections from TD-DFT). Calculated KS band gaps systematically *underestimate* experimental band gaps by 30–100%, a phenomenon that has nothing to do with the quality of the functional and everything to do with what KS eigenvalues mean (or do not). See §5.6.

### Worked example: KS eigenvalues of the helium atom

Consider the helium atom (two electrons, both in the 1s orbital). The exact ionisation potential is $I = 24.59\;\text{eV} = 0.904\;\text{Ha}$. Different methods give very different orbital energies:

| Method | $\varepsilon_\mathrm{HOMO}$ (Ha) | $-\varepsilon_\mathrm{HOMO}$ (eV) | Error vs. $I$ |
|---|---|---|---|
| Hartree–Fock | $-0.918$ | $24.98$ | $+0.4\;\text{eV}$ (Koopmans good) |
| LDA | $-0.570$ | $15.51$ | $-9.1\;\text{eV}$ |
| PBE | $-0.585$ | $15.92$ | $-8.7\;\text{eV}$ |
| Exact KS (van Leeuwen reconstruction) | $-0.904$ | $24.59$ | $0$ exact |

The LDA and PBE HOMOs are off by 8–9 eV. The reason is self-interaction error and the wrong long-range behaviour of $v_{xc}$ (§5.6). The exact KS construction — found by inverting an exact reference density (e.g., from quantum Monte Carlo) — recovers the correct HOMO eigenvalue, confirming (5.32). The Hartree–Fock orbital energy happens to be a good approximation here because helium has only two electrons and self-interaction in HF exchange is exactly zero by construction.

!!! note "The exact-KS HOMO in the table is not hand-reproducible"
    The "Exact KS (van Leeuwen reconstruction)" row is *not* something you can derive with pencil and paper. It is obtained by numerically *inverting* an essentially exact reference density: one starts from a highly accurate interacting density $n_0(\mathbf r)$ (from a near-exact wavefunction calculation, e.g. quantum Monte Carlo or a large configuration-interaction expansion) and iteratively searches for the local potential $v_s(\mathbf r)$ whose non-interacting ground-state density reproduces $n_0$ (van Leeuwen and Baerends, 1994). The resulting $v_s$, and hence $\varepsilon_\mathrm{HOMO}$, come out of a numerical optimisation, not a closed-form expression. The value $-0.904\;\text{Ha}$ is quoted to *illustrate* that the exact KS theory satisfies (5.32) — it is a confirmation of the theorem, not a calculation you should expect to repeat by hand.

### The HOMO exception: Janak's theorem and exact DFT

There is exactly one KS eigenvalue with a guaranteed physical interpretation. In *exact* Kohn–Sham theory (with the exact, unknown $E_{xc}$), the highest occupied KS eigenvalue $\varepsilon_\mathrm{HOMO}$ equals minus the first ionisation potential of the system:

$$
\varepsilon_\mathrm{HOMO} = -I.
\tag{5.32}
$$

This follows from the asymptotic decay of the density: the density of a finite system decays as $n(\mathbf r) \sim e^{-2\sqrt{2I}\,r}$ for large $r$, and since the KS orbitals must reproduce this decay, and the most slowly-decaying occupied orbital governs it, $\varepsilon_\mathrm{HOMO}$ must equal $-I$ exactly. (The derivation is due to Almbladh and von Barth.)

??? note "Full derivation: equating the asymptotic exponents"
    **The KS-orbital tail.** A finite neutral system has $v_s(\mathbf r) \to 0$ as $r\to\infty$ (the effective potential decays to zero far from the molecule). In that region the canonical KS equation (5.29) for an occupied orbital of eigenvalue $\varepsilon < 0$ reduces to

    $$
    -\tfrac12\nabla^{2}\phi(\mathbf r) = \varepsilon\,\phi(\mathbf r), \qquad r\to\infty.
    $$

    Seeking a spherically-symmetric decaying solution $\phi \sim e^{-\kappa r}$, the dominant large-$r$ term of the radial Laplacian is $\nabla^{2}\phi \to \kappa^{2}\phi$ (the $-2\kappa/r$ correction is subleading). Substituting,

    $$
    -\tfrac12\,\kappa^{2}\,\phi = \varepsilon\,\phi \;\;\Longrightarrow\;\; \kappa = \sqrt{2|\varepsilon|}, \qquad \phi(\mathbf r) \sim e^{-\sqrt{2|\varepsilon|}\,r}.
    \tag{5.32a}
    $$

    A more negative $\varepsilon$ gives a faster decay; therefore the *slowest*-decaying occupied orbital is the one with the *least* negative eigenvalue — the HOMO. At large $r$ the density $n = \sum_i f_i|\phi_i|^{2}$ is dominated by this slowest tail, so

    $$
    n(\mathbf r) \sim |\phi_\mathrm{HOMO}|^{2} \sim e^{-2\sqrt{2|\varepsilon_\mathrm{HOMO}|}\,r}.
    \tag{5.32b}
    $$

    **The exact many-body tail.** Independently, the exact interacting density of a finite system decays as

    $$
    n(\mathbf r) \sim e^{-2\sqrt{2I}\,r}, \qquad r\to\infty,
    \tag{5.32c}
    $$

    where $I = E(N-1) - E(N)$ is the first ionisation potential. This is a rigorous many-body result (Morrell, Parr and Levy, 1975; Almbladh and von Barth, 1985): the rate of decay of the density is set by the energy required to remove the most weakly-bound electron.

    **Equate the exponents.** By construction the KS density equals the exact density, so the two tails (5.32b) and (5.32c) must have identical exponents:

    $$
    2\sqrt{2|\varepsilon_\mathrm{HOMO}|}\,r = 2\sqrt{2I}\,r \;\;\Longrightarrow\;\; |\varepsilon_\mathrm{HOMO}| = I.
    $$

    Since the HOMO of a bound system has $\varepsilon_\mathrm{HOMO} < 0$, this is $\varepsilon_\mathrm{HOMO} = -I$, which is (5.32).

In *approximate* KS theory — i.e., every calculation you will ever do — (5.32) holds only approximately, and badly so for LDA and GGA. Functionals with the correct asymptotic behaviour (some range-separated hybrids, the optimised effective potential method, etc.) do better. The other KS eigenvalues do *not* have such an exact interpretation, even with the exact functional.

### Janak's theorem, more generally

Janak (1978) proved that, in KS theory with fractional occupations $f_i \in [0,1]$,

$$
\frac{\partial E}{\partial f_i} = \varepsilon_i.
\tag{5.33}
$$

This identifies $\varepsilon_i$ as the derivative of the total energy with respect to occupation — a useful identity for some manipulations, but again *not* an ionisation-energy interpretation except for the HOMO at integer occupations.

!!! note "Why this step?"
    Janak's theorem (5.33) is the deepest justification for the *thermal smearing* of occupations in SCF calculations (§5.5). It says: changing the occupation of orbital $i$ by an infinitesimal amount $\mathrm d f_i$ changes the total energy by $\varepsilon_i\,\mathrm d f_i$. In a metal where many states cluster near the Fermi level $\varepsilon_F$, the energy difference between adjacent integer-occupation states is small but the *which* states are occupied changes abruptly from one SCF iteration to the next. Janak smoothes this out: by allowing each near-Fermi state to be fractionally occupied with a Fermi–Dirac factor at some artificial temperature, the total energy becomes a smooth function of the orbital energies and the SCF converges. The penalty for this trick is that one must extrapolate to $T\to 0$ to obtain the true ground-state energy — typically done by the entropy correction $E_{T=0} = E - \tfrac{1}{2}T S$ for Gaussian smearing (the so-called "thermal smearing" correction).

    The coefficient $\tfrac{1}{2}$ here is *not* universal: it is specific to Gaussian (and Methfessel–Paxton) smearing, where the free energy $F = E - TS$ and the internal energy are related such that the $T\to0$ extrapolated energy is the average $\tfrac12(E + F) = E - \tfrac12 TS$. Other smearing schemes (Fermi–Dirac, cold smearing) carry different coefficients; the schemes and their corrections are treated in full in Section 5.5.

### Connection to fractional charges

Combining Janak's theorem with the variational principle gives a powerful diagnostic of approximate functionals. For the *exact* functional, the total energy $E(N)$ as a function of (continuous) electron number is piecewise *linear* between integers, with kinks at each integer. The slope between $N$ and $N+1$ is $-A$ (negative of the electron affinity); the slope between $N-1$ and $N$ is $-I$. The kink at integer $N$ is the *derivative discontinuity* $\Delta = I - A - (\varepsilon_\mathrm{LUMO}-\varepsilon_\mathrm{HOMO})$ discussed in §5.6.

!!! note "Where $\Delta$ comes from (defined in Section 5.6, quoted here)"
    The quantity $\Delta$ is the **derivative discontinuity**, defined and analysed in Section 5.6; we only quote it here. The one-line argument is a bookkeeping of two different gaps. The **fundamental gap** is the true charged-excitation gap, $E_g = I - A$, obtained from total-energy differences of the $N$- and $(N\pm1)$-electron systems. The **Kohn–Sham gap** is the eigenvalue difference $E_g^\mathrm{KS} = \varepsilon_\mathrm{LUMO} - \varepsilon_\mathrm{HOMO}$ of the *fixed-$N$* auxiliary system. In exact KS theory $\varepsilon_\mathrm{HOMO} = -I$ (equation (5.32)), but $\varepsilon_\mathrm{LUMO} \neq -A$ in general; the two differ precisely because the exact $v_{xc}$ jumps by a spatially-constant amount $\Delta$ as the electron number crosses an integer. Subtracting,

    $$
    \Delta = \underbrace{(I - A)}_{\text{fundamental gap}} - \underbrace{(\varepsilon_\mathrm{LUMO} - \varepsilon_\mathrm{HOMO})}_{\text{KS gap}},
    $$

    i.e. the derivative discontinuity is exactly the amount by which the KS gap falls short of the true gap. Even with the exact functional, the KS gap underestimates $E_g$ by $\Delta$; with LDA/GGA, where $\Delta$ is effectively missing, the underestimate is compounded.

Approximate functionals (LDA, GGA) violate piecewise linearity: their $E(N)$ is *concave-up* between integers, with no kink. This violation is the formal expression of *self-interaction error*: electrons artificially want to delocalise to non-integer numbers because the approximate functional underpenalises fractional charges. The integral of this curvature error over an integer interval $\int_{N-1}^{N}[E_\mathrm{approx}(N')-E_\mathrm{linear}(N')]\,\mathrm dN' \sim 0.1$–$1\;\text{eV}$ is the diagnostic, and modern functional development (e.g., $\omega$B97X-V, optimally-tuned range-separated hybrids) explicitly minimises it.

!!! question "Check yourself"
    1. What *is* the auxiliary system in the Kohn–Sham construction, and what single property does it share — by design — with the real interacting system?
    2. In the decomposition $E = T_s + U_H + E_{xc} + \int v_\mathrm{ext}\,n$, which pieces are computed exactly and which one holds all the approximation?
    3. Are the Kohn–Sham eigenvalues $\varepsilon_i$ excitation energies? Is there any exception?

    ??? success "Answers"
        1. A fictitious system of **non-interacting** electrons moving in an effective potential $v_\mathrm{KS}$. It is chosen so that its ground-state density $n(\mathbf r) = \sum_i|\psi_i|^2$ equals the density of the real interacting system — that shared density is the only thing tying the two together.
        2. $T_s$ (exact, from the orbitals), $U_H$ (exact, classical electrostatics) and $\int v_\mathrm{ext}\,n$ (exact) are all computed without approximation; every approximation in a real DFT calculation lives in $E_{xc}$ alone.
        3. No — in general the $\varepsilon_i$ are Lagrange multipliers of the auxiliary system, not excitation or removal energies. The one guaranteed exception is the HOMO in *exact* DFT, where $\varepsilon_\mathrm{HOMO} = -I$. Differences such as $\varepsilon_a - \varepsilon_i$ are at best crude estimates, and KS gaps underestimate true gaps (Section 5.6).

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

!!! note "The clean Levy-constrained-search statement"
    The phrasing above ("$T$ is a minimum over wavefunctions ... while $T_s$ is a minimum ...") is correct but easy to misread, because $T$ and $T_s$ are not minima of the *same* object. The Levy constrained search makes the inequality transparent by searching both over the **same density class** — all antisymmetric wavefunctions $\Psi$ that yield the given density $n$:

    $$
    T_s[n] = \min_{\Psi \to n} \langle\Psi|\hat T|\Psi\rangle,
    \tag{5.33c}
    $$

    where $\hat T = -\tfrac12\sum_i\nabla_i^{2}$ is the bare kinetic operator (no interaction appears in this search — that is what makes the minimiser a single Slater determinant, the KS state). The interacting ground state $\Psi_0$ also belongs to the search set, because by the KS postulate $\Psi_0 \to n$ reproduces the same density. It is therefore merely *one candidate* in the minimisation (5.33c), not generally the minimiser. Since a minimum is no larger than the value at any candidate,

    $$
    T_s[n] = \min_{\Psi \to n}\langle\Psi|\hat T|\Psi\rangle \;\le\; \langle\Psi_0|\hat T|\Psi_0\rangle = T[n].
    $$

    The inequality is *not* a statement that interaction raises kinetic energy in some dynamical sense; it is the elementary fact that constraining the search to non-interacting trial states can only lower (or equal) the minimum of $\langle\hat T\rangle$ over that fixed density.

The difference $T_c[n] \equiv T[n] - T_s[n] \geq 0$ is the **correlation kinetic energy** and is part of $E_{xc}$. For most chemical systems $T_c$ is of order 10–50 millihartree per electron, while $T_s$ is of order 1 Hartree per electron: a few per cent correction. By computing $T_s$ exactly via the orbitals, KS theory captures the dominant kinetic energy without approximation; only the small remainder $T_c$ is bundled into $E_{xc}$ and approximated.

This is precisely the failure mode that doomed Thomas–Fermi (§5.1). TF tried to write the entire $T$ as an explicit functional of $n$ — $C_F\int n^{5/3}$ — and got the magnitude roughly right but the spatial dependence so wrong that no molecule could bind. KS bypasses the problem by introducing orbitals to compute $T_s$ directly, leaving $T_c$ (which is small and smooth) to be approximated as a functional of $n$. The orbitals are the price; the prize is chemistry.

### The exchange–correlation hole picture

A complementary way to understand $E_{xc}$ is through the *exchange–correlation hole* $n_{xc}(\mathbf r, \mathbf r')$, defined as the difference between the conditional probability of finding a second electron at $\mathbf r'$ given one at $\mathbf r$ and the unconditional density $n(\mathbf r')$. The hole encodes Pauli exclusion (no two same-spin electrons in the same place) and dynamical correlation (electrons avoid each other beyond the Pauli requirement). The exact exchange–correlation energy is

$$
E_{xc}[n] = \tfrac{1}{2}\iint\frac{n(\mathbf r)\,\bar n_{xc}(\mathbf r,\mathbf r')}{|\mathbf r - \mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r',
\tag{5.33b}
$$

where $\bar n_{xc}$ is the coupling-strength-averaged hole. Two exact sum rules constrain the hole: $\int\bar n_{xc}(\mathbf r,\mathbf r')\,\mathrm d\mathbf r' = -1$ (the hole removes exactly one electron's worth of charge), and the hole is non-positive in the same-spin channel (Pauli). Modern functional construction proceeds largely by *modelling the hole* — designing $\bar n_{xc}$ to satisfy as many exact constraints as possible — and then integrating (5.33b). This is the underlying philosophy of PBE (built from a real-space hole model) and SCAN (built from constraint satisfaction).

!!! note "KS is *not* a free lunch"
    The KS construction costs us the linear scaling with system size that orbital-free DFT (Thomas–Fermi-like methods) enjoys. Diagonalising the KS Hamiltonian scales formally as $\mathcal O(N^{3})$ in the number of electrons; the Hartree term as $\mathcal O(N^{2})$ (or $\mathcal O(N\log N)$ with FFT). Modern algorithms reduce these to near-linear for sparse systems, but the constant in front is large compared with orbital-free methods. The trade is worth it because the result is *chemically accurate*; we discuss this scaling in detail in Chapter 6.

### Fractional occupations and Mermin's theorem

The KS construction generalises naturally to *finite electronic temperature*. Mermin (1965) proved that the Hohenberg–Kohn theorems extend to the grand-canonical ensemble: at temperature $T$ and chemical potential $\mu$, the equilibrium density $n_\mathrm{eq}(\mathbf r; T, \mu)$ determines $v_\mathrm{ext}$ up to a constant, and a free-energy functional $\Omega[n; T, \mu]$ is minimised at the equilibrium density. The KS implementation replaces the integer step-function occupations with the Fermi–Dirac distribution

$$
f_i = \frac{1}{1 + \mathrm e^{(\varepsilon_i - \mu)/k_B T}},
\tag{5.31a}
$$

so that the density becomes $n(\mathbf r) = \sum_i f_i |\phi_i(\mathbf r)|^{2}$ where the sum is now over *all* eigenstates (not just the lowest $N$). The free energy adds an entropy contribution

$$
S = -k_B\sum_i\big[f_i\ln f_i + (1-f_i)\ln(1-f_i)\big],
$$

and the variational object is $\Omega = E - T S - \mu N$. Setting $\partial\Omega/\partial f_i = 0$ gives back (5.31a) — Janak's theorem $\partial E/\partial f_i = \varepsilon_i$ (which we derived in §5.2.4b) plus the entropy derivative.

??? note "Full derivation: minimising $\Omega$ recovers Fermi–Dirac"
    Differentiate $\Omega = E - TS - \mu N$ with respect to the occupation $f_i$, treating each $f_i$ as an independent variable. There are three contributions.

    *Energy term.* By Janak's theorem (5.33), $\dfrac{\partial E}{\partial f_i} = \varepsilon_i$.

    *Particle-number term.* With $N = \sum_j f_j$, we have $\dfrac{\partial N}{\partial f_i} = 1$, so $-\mu\,\partial N/\partial f_i = -\mu$.

    *Entropy term.* From $S = -k_B\sum_j[f_j\ln f_j + (1-f_j)\ln(1-f_j)]$, differentiate the single $j=i$ term:

    $$
    \frac{\partial}{\partial f_i}\big[f_i\ln f_i\big] = \ln f_i + 1,
    \qquad
    \frac{\partial}{\partial f_i}\big[(1-f_i)\ln(1-f_i)\big] = -\ln(1-f_i) - 1.
    $$

    Adding these, the $+1$ and $-1$ cancel:

    $$
    \frac{\partial S}{\partial f_i} = -k_B\big[\ln f_i - \ln(1-f_i)\big] = -k_B\ln\frac{f_i}{1-f_i}.
    \tag{5.31d}
    $$

    *Assemble.* The stationarity condition $\partial\Omega/\partial f_i = 0$ is

    $$
    \varepsilon_i - T\left(-k_B\ln\frac{f_i}{1-f_i}\right) - \mu = 0
    \;\;\Longrightarrow\;\;
    \varepsilon_i + k_B T\ln\frac{f_i}{1-f_i} - \mu = 0.
    $$

    Solve for $f_i$. Isolate the logarithm:

    $$
    \ln\frac{f_i}{1-f_i} = -\frac{\varepsilon_i - \mu}{k_B T}
    \;\;\Longrightarrow\;\;
    \frac{f_i}{1-f_i} = \mathrm e^{-(\varepsilon_i-\mu)/k_B T}.
    $$

    Write $x = \mathrm e^{(\varepsilon_i-\mu)/k_B T}$, so $f_i/(1-f_i) = 1/x$, i.e. $x f_i = 1 - f_i$, hence $f_i(1+x) = 1$ and

    $$
    f_i = \frac{1}{1+x} = \frac{1}{1 + \mathrm e^{(\varepsilon_i-\mu)/k_B T}},
    $$

    which is exactly the Fermi–Dirac distribution (5.31a). The entropy term is therefore precisely what turns the integer step-function occupations into smooth fractional ones.

!!! note "Why this step?"
    The mathematical content of finite-$T$ DFT is that smoothing the occupations from discontinuous step functions to smooth Fermi–Dirac distributions makes the energy functional differentiable in *every* variable, including the orbital occupations. For metals, where states near the Fermi level are nearly degenerate and tiny changes in the KS potential can flip which states are occupied, this smoothing is *essential* for SCF convergence. We return to this in §5.5. Even for $T=0$ insulator calculations, the underlying mathematical framework requires the freedom of fractional occupations because the constrained-search functional $F_L[n]$ is defined on the full $N$-representable class, which includes mixed-state densities corresponding to fractional occupations.

### Pseudo-eigenvalue equation in matrix form

In a finite basis $\{\chi_\mu\}_{\mu=1}^{M}$ (Gaussian-type orbitals, plane waves, atom-centred numerical orbitals), each KS orbital is expanded as $\phi_i(\mathbf r) = \sum_\mu c_{\mu i}\chi_\mu(\mathbf r)$, and the KS equations become a generalised matrix eigenvalue problem

$$
\mathbf H\mathbf c_i = \varepsilon_i\,\mathbf S\mathbf c_i,
\tag{5.31b}
$$

where $H_{\mu\nu} = \langle\chi_\mu|\hat h_\mathrm{KS}|\chi_\nu\rangle$ is the KS Hamiltonian matrix and $S_{\mu\nu} = \langle\chi_\mu|\chi_\nu\rangle$ is the overlap matrix (which equals the identity in an orthonormal basis such as plane waves). The dimension of the matrix is $M$, the basis size — typically $10^{2}$–$10^{6}$ depending on system size and basis choice. The cost of one diagonalisation step is $\mathcal O(M^{3})$, which dominates the total SCF cost for moderate-size systems and motivates the development of $\mathcal O(N)$ methods (Chapter 6).

### The adiabatic connection: from non-interacting to interacting

A central conceptual tool for understanding the Kohn–Sham construction is the *adiabatic connection*. Consider a family of Hamiltonians

$$
\hat H_\lambda \;=\; \hat T \;+\; \lambda\,\hat V_{ee} \;+\; \hat V_\mathrm{ext}^{\lambda},
$$

where $\lambda\in[0,1]$ scales the electron–electron interaction continuously from zero (non-interacting KS system) to the physical value (real interacting system). At each $\lambda$, choose $\hat V_\mathrm{ext}^{\lambda}$ so that the ground-state density equals the *physical* $n_0$. The exchange–correlation energy can then be written as the integral

$$
E_{xc}[n] \;=\; \int_0^{1}\!\mathrm d\lambda\,\langle\Psi_\lambda|\hat V_{ee}|\Psi_\lambda\rangle - U_H[n].
\tag{5.33a}
$$

??? note "Full derivation: the adiabatic-connection formula via Hellmann–Feynman"
    Let $E_\lambda = \langle\Psi_\lambda|\hat H_\lambda|\Psi_\lambda\rangle$ be the ground-state energy of $\hat H_\lambda = \hat T + \lambda\hat V_{ee} + \hat V_\mathrm{ext}^{\lambda}$, where at each $\lambda$ the external potential $\hat V_\mathrm{ext}^{\lambda}$ is *adjusted* so that the ground-state density is held fixed at the physical $n_0$ for all $\lambda$.

    **Hellmann–Feynman along the path.** Differentiating $E_\lambda$ with respect to $\lambda$, the terms where $\partial/\partial\lambda$ hits the wavefunction vanish because $\Psi_\lambda$ is an eigenstate (normalised), leaving only the explicit derivative of the Hamiltonian:

    $$
    \frac{\mathrm dE_\lambda}{\mathrm d\lambda} = \Big\langle\Psi_\lambda\Big|\frac{\partial\hat H_\lambda}{\partial\lambda}\Big|\Psi_\lambda\Big\rangle
    = \langle\Psi_\lambda|\hat V_{ee}|\Psi_\lambda\rangle + \Big\langle\Psi_\lambda\Big|\frac{\partial\hat V_\mathrm{ext}^{\lambda}}{\partial\lambda}\Big|\Psi_\lambda\Big\rangle.
    $$

    Because $\hat V_\mathrm{ext}^{\lambda} = \int v_\mathrm{ext}^{\lambda}(\mathbf r)\,\hat n(\mathbf r)\,\mathrm d\mathbf r$ is a one-body operator and the density is pinned to $n_0$ independent of $\lambda$,

    $$
    \Big\langle\Psi_\lambda\Big|\frac{\partial\hat V_\mathrm{ext}^{\lambda}}{\partial\lambda}\Big|\Psi_\lambda\Big\rangle = \int \frac{\partial v_\mathrm{ext}^{\lambda}(\mathbf r)}{\partial\lambda}\,n_0(\mathbf r)\,\mathrm d\mathbf r.
    $$

    **Integrate from 0 to 1.** Integrating the total derivative recovers the endpoint difference:

    $$
    E_{\lambda=1} - E_{\lambda=0} = \int_0^1 \langle\Psi_\lambda|\hat V_{ee}|\Psi_\lambda\rangle\,\mathrm d\lambda + \int_0^1\!\!\int \frac{\partial v_\mathrm{ext}^{\lambda}}{\partial\lambda}\,n_0\,\mathrm d\mathbf r\,\mathrm d\lambda.
    $$

    The last term integrates trivially in $\lambda$ to $\int (v_\mathrm{ext}^{1} - v_\mathrm{ext}^{0})\,n_0\,\mathrm d\mathbf r$.

    **Identify the endpoints.** At $\lambda=1$ the system is fully interacting: $\Psi_1 = \Psi_0$, $v_\mathrm{ext}^{1} = v_\mathrm{ext}$ (the physical potential), and $E_1 = T + V_{ee} + \int v_\mathrm{ext} n_0$ is the exact energy (5.25). At $\lambda=0$ there is no interaction: $\Psi_0^{\lambda=0}$ is the KS determinant, $v_\mathrm{ext}^{0} = v_\mathrm{KS}$, and $E_0 = T_s + \int v_\mathrm{KS} n_0$. Subtracting,

    $$
    E_1 - E_0 = \big(T - T_s\big) + V_{ee} + \int (v_\mathrm{ext} - v_\mathrm{KS})\,n_0\,\mathrm d\mathbf r,
    $$

    while the right-hand side of the integrated Hellmann–Feynman relation gives

    $$
    E_1 - E_0 = \int_0^1 \langle\Psi_\lambda|\hat V_{ee}|\Psi_\lambda\rangle\,\mathrm d\lambda + \int (v_\mathrm{ext} - v_\mathrm{KS})\,n_0\,\mathrm d\mathbf r.
    $$

    The external-potential pieces are identical on both sides and cancel, leaving

    $$
    \big(T - T_s\big) + V_{ee} = \int_0^1 \langle\Psi_\lambda|\hat V_{ee}|\Psi_\lambda\rangle\,\mathrm d\lambda.
    \tag{5.33d}
    $$

    **Assemble $E_{xc}$.** By the definition (5.27), $E_{xc} = (T - T_s) + (V_{ee} - U_H)$. Substituting (5.33d) for $(T-T_s)+V_{ee}$ and subtracting the Hartree energy $U_H$ (which is $\lambda$-independent — it depends only on the fixed density $n_0$):

    $$
    E_{xc}[n] = \int_0^1 \langle\Psi_\lambda|\hat V_{ee}|\Psi_\lambda\rangle\,\mathrm d\lambda - U_H[n],
    $$

    which is (5.33a). The coupling-strength integral thus assembles the kinetic correlation $T-T_s$ together with the non-classical part of $V_{ee}$ into $E_{xc}$, with the classical Hartree piece removed.

This *Levy–Langreth adiabatic connection formula* (independently derived by Harris–Jones and Gunnarsson–Lundqvist) gives an exact, explicit expression for $E_{xc}$ as an integral of the interaction energy along the adiabatic path. It is the foundation of the modern view of exchange–correlation as an "exchange–correlation hole" averaged over coupling strength, and it underpins virtually every modern functional construction including hybrids (which can be derived as the trapezoidal-rule approximation to (5.33a)) and double-hybrids (which add information from the $\lambda=0$ endpoint via MP2-style correlation).

!!! note "Why this step?"
    The adiabatic-connection identity (5.33a) is *exact*: no approximation has been made. It tells us that the unknown $E_{xc}[n]$ can in principle be computed by solving the many-body problem at each $\lambda$ — which is hard, but at least is a well-defined recipe. More importantly, it gives us a way to design *approximate* functionals: we can interpolate the integrand between known endpoints ($\lambda=0$, where the wavefunction is the KS determinant; $\lambda=1$, where the system is fully interacting) using a sensible model. The simplest interpolation — a straight line — gives a 50/50 mix of exact and non-interacting exchange and is one rigorous justification for the hybrid functional ansatz.

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

A more convenient pair of variables is the *total density* $n = n_\uparrow + n_\downarrow$ and the *spin polarisation*

$$
\zeta(\mathbf r) \;\equiv\; \frac{n_\uparrow(\mathbf r) - n_\downarrow(\mathbf r)}{n_\uparrow(\mathbf r) + n_\downarrow(\mathbf r)}\;\in\;[-1, +1].
\tag{5.31c}
$$

A non-magnetic system has $\zeta\equiv 0$ everywhere; a fully spin-polarised electron gas ($n_\downarrow = 0$) has $\zeta = +1$. The exchange energy of a spin-polarised uniform gas obeys the elegant *spin-interpolation formula*

$$
\epsilon_x(n,\zeta) = \epsilon_x(n,0)\cdot\Big[\tfrac{1}{2}(1+\zeta)^{4/3} + \tfrac{1}{2}(1-\zeta)^{4/3}\Big],
$$

which interpolates smoothly between the unpolarised and fully polarised limits. PBE and similar GGAs use analogous interpolations for the spin-dependent correlation energy. The two limits $\zeta=0$ and $\zeta=\pm 1$ are exactly known from many-body theory; the interpolation in between is constrained but has some freedom.

??? note "Full derivation: the $\tfrac12[(1+\zeta)^{4/3}+(1-\zeta)^{4/3}]$ spin-scaling factor"
    The bracket follows from two exact facts about exchange.

    **Fact 1 — the spin-scaling relation.** Exchange acts only between *same-spin* electrons (opposite-spin electrons have no exchange interaction). Hence the total exchange energy splits cleanly into independent spin channels, and the exact relation

    $$
    E_x[n_\uparrow, n_\downarrow] = \tfrac12\Big(E_x[2n_\uparrow] + E_x[2n_\downarrow]\Big)
    \tag{5.34a}
    $$

    holds, where $E_x[n]$ on the right is the *spin-unpolarised* functional. The factor of 2 inside the brackets and the $\tfrac12$ outside come from treating each spin channel as a fully polarised system of density $2n_\sigma$ and weighting it by one-half.

    **Fact 2 — the density scaling of the exchange energy density.** For the uniform electron gas (Dirac exchange, derived in §5.1), the exchange energy per particle scales as the cube root of the density,

    $$
    \epsilon_x(n) = -C_x\,n^{1/3}, \qquad C_x = \tfrac34\Big(\tfrac{3}{\pi}\Big)^{1/3},
    $$

    so the exchange energy *density* (energy per unit volume) is $n\,\epsilon_x(n) = -C_x\,n^{4/3}$, scaling as $n^{4/3}$.

    **Combine.** Apply (5.34a) at the level of energy densities. Writing the spin densities in terms of total density and polarisation,

    $$
    n_\uparrow = \tfrac12 n(1+\zeta), \qquad n_\downarrow = \tfrac12 n(1-\zeta),
    $$

    so that $2n_\uparrow = n(1+\zeta)$ and $2n_\downarrow = n(1-\zeta)$. The spin-resolved exchange energy density is

    $$
    n\,\epsilon_x(n,\zeta) = \tfrac12\Big[(2n_\uparrow)\,\epsilon_x(2n_\uparrow) + (2n_\downarrow)\,\epsilon_x(2n_\downarrow)\Big]
    = \tfrac12\Big[-C_x(2n_\uparrow)^{4/3} - C_x(2n_\downarrow)^{4/3}\Big].
    $$

    Substitute $2n_\sigma = n(1\pm\zeta)$ and factor out $n^{4/3}$:

    $$
    n\,\epsilon_x(n,\zeta) = -C_x\,n^{4/3}\cdot\tfrac12\Big[(1+\zeta)^{4/3} + (1-\zeta)^{4/3}\Big].
    $$

    Divide both sides by $n$ and recognise $-C_x n^{1/3} = \epsilon_x(n,0)$, the unpolarised energy per particle:

    $$
    \epsilon_x(n,\zeta) = \epsilon_x(n,0)\cdot\tfrac12\Big[(1+\zeta)^{4/3} + (1-\zeta)^{4/3}\Big].
    $$

    **Check the limits.** At $\zeta=0$ the bracket is $\tfrac12(1+1)=1$, recovering $\epsilon_x(n,0)$. At $\zeta=1$ (full polarisation) the bracket is $\tfrac12(2^{4/3}+0) = 2^{1/3}\approx 1.26$, the exact enhancement of a fully spin-polarised gas. Both endpoints match the many-body results, as claimed in the main text.

In what follows, unless stated otherwise, we work with the spin-unpolarised theory; the generalisation is straightforward but notationally heavier.

### Non-collinear magnetism

For systems where the local magnetisation vector $\mathbf m(\mathbf r) = (n_{\uparrow\downarrow} + n_{\downarrow\uparrow},\;-i(n_{\uparrow\downarrow}-n_{\downarrow\uparrow}),\;n_\uparrow - n_\downarrow)$ rotates in space, the density becomes a $2\times 2$ Hermitian matrix in spin space and the KS Hamiltonian is a $2\times 2$ matrix-valued operator. Spin–orbit coupling, in particular, couples the spin and spatial degrees of freedom and demands this full non-collinear treatment. Practical implementations exist in essentially every modern plane-wave code (VASP, Quantum ESPRESSO, ABINIT), at roughly twice the cost of collinear-spin DFT.

## 5.3.6a Reading the KS orbitals: what they do and do not mean

Because KS orbitals look like wavefunctions (one orbital per electron, occupation numbers 0 or 2, eigenvalues that order themselves into shells and bands), it is enormously tempting to interpret them as physical objects: a "HOMO orbital" that an electron occupies, a "LUMO" that it could jump to, a "band structure" that describes electron motion. Some of this interpretation is justified, some of it is approximate, and some of it is wrong.

**What KS orbitals are guaranteed to do.** Reproduce the exact density via $n = \sum_i|\phi_i|^{2}$. Give the exact non-interacting kinetic energy $T_s$. Provide a Slater-determinant approximation to the wavefunction that can be used as a starting point for post-DFT methods (MP2, coupled cluster, GW).

**What KS orbitals approximately do.** Look qualitatively like Hartree–Fock orbitals (with similar nodal structure and spatial extent). Provide ballpark molecular-orbital pictures useful for chemical intuition: where electrons "live", the symmetries of bonds, the orientation of lone pairs. For most molecules, the qualitative HOMO–LUMO picture from a PBE calculation is correct.

**What KS orbitals do NOT do.** They are not the wavefunction of any physical state — the interacting wavefunction $\Psi_0$ is not a Slater determinant of $\{\phi_i\}$ in general, and the energy of $\Psi_0$ is not the sum of orbital energies. The KS eigenvalues are not ionisation potentials or electron affinities, except for the HOMO in the exact-functional limit (see §5.3.4). And: orbital energies *do not transfer* between similar systems. The HOMO energy of a hydrogen atom calculated with PBE is *not* directly comparable to the HOMO energy of methane calculated with PBE; the absolute reference (which is set by the asymptotic behaviour of $v_{xc}$) varies.

!!! warning "Pitfall: 'KS band structure'"
    Plotting $\varepsilon_i(\mathbf k)$ for a solid and calling it the "band structure" is a common practice. The result is a useful diagnostic of the KS eigenspectrum, often qualitatively similar to the quasiparticle band structure measured in ARPES, but quantitatively *not* the same thing. The KS band gap is *not* the fundamental gap, the band widths are typically too narrow (LDA/GGA) or too broad (HSE06), and the orbital character at any $\mathbf k$-point is approximate. For semiquantitative band-structure work, GW corrections on top of a KS starting point are the standard remedy.

### Summary of §5.3 — what to remember in 3 months

- **The trick**: introduce a fictitious non-interacting system with the *same density* as the real one. Solve $N$ one-electron Schrödinger equations.
- **KS equations (Theorem 5.3 of the chapter)**: $[-\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}(\mathbf r)]\phi_i = \varepsilon_i\phi_i$, with $v_\mathrm{KS} = v_\mathrm{ext} + v_H + v_{xc}$ and $n = \sum_i|\phi_i|^{2}$.
- **Energy decomposition**: $E[n] = T_s[n] + U_H[n] + E_{xc}[n] + \int v_\mathrm{ext} n$. Exact for the right $E_{xc}$.
- **KS eigenvalues are not energies**: they are Lagrange multipliers. Only the HOMO has a guaranteed meaning ($\varepsilon_\mathrm{HOMO} = -I$ with the exact functional); the rest are diagnostic but not physical.
- **Janak's theorem**: $\partial E/\partial f_i = \varepsilon_i$, justifying fractional occupations / thermal smearing.
- **Adiabatic connection (5.33a)**: $E_{xc} = \int_0^{1}\langle\Psi_\lambda|\hat V_{ee}|\Psi_\lambda\rangle\mathrm d\lambda - U_H$ — gives an exact route to $E_{xc}$ via the exchange–correlation hole.
- **Generalised KS**: hybrid functionals use a non-local Fock operator instead of a local $v_x$; the variational structure is preserved.
- **Spin and non-collinearity**: replace $n$ by $(n_\uparrow,n_\downarrow)$ or by a $2\times 2$ density matrix for magnets.

!!! note "Remark: orbital interpretation in the literature"
    Despite the warnings here, almost every chemistry paper that discusses "the HOMO–LUMO gap" or "frontier orbital theory" implicitly treats KS orbitals as physical. This is mostly harmless — KS orbitals look like Hartree–Fock orbitals, and chemical intuition built on them tends to be qualitatively right. But it is wrong in principle, and it is wrong in practice for band gaps, ionisation potentials of inner orbitals, and excited-state ordering. Use KS orbitals as a heuristic, not as a source of ground truth.

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

### A note on generalised Kohn–Sham

The standard Kohn–Sham construction uses a *local* exchange–correlation potential $v_{xc}(\mathbf r)$. Hybrid functionals (§5.4.5) introduce non-local *Fock* exchange, whose corresponding "potential" is an integral operator $\hat\Sigma_x$ that takes an orbital to another orbital. The variational principle in this case gives a *generalised Kohn–Sham* (gKS) scheme, in which the one-electron Hamiltonian is

$$
\hat h_\mathrm{gKS} \;=\; -\tfrac{1}{2}\nabla^{2} + v_\mathrm{ext}(\mathbf r) + v_H[n](\mathbf r) + (1-a)v_x^\mathrm{LDA/GGA}[n](\mathbf r) + a\,\hat\Sigma_x^\mathrm{HF} + v_c[n](\mathbf r),
$$

with a fraction $a$ of exact exchange. This still defines a one-electron eigenvalue problem, still gives orbitals that reproduce the gKS density, but the auxiliary system is now *not* purely non-interacting — it interacts through the Fock operator. The HK theorems extend to gKS (Seidl, Görling, Vogl, Majewski, and Levy, 1996), and most modern functionals (HSE06, B3LYP, $\omega$B97X) are gKS schemes in this sense. The practical implication is that hybrid-functional KS eigenvalues are *closer* to quasiparticle energies than pure KS eigenvalues, because they include some of the exact-exchange physics that pure KS hides inside $v_{xc}$.


