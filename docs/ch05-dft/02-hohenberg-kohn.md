# 5.2 The Hohenberg–Kohn Theorems

!!! info "Why this section exists"
    Here is the question this whole section answers: can the ground-state *density* alone — a single function $n(\mathbf r)$ of three coordinates — in principle determine *everything* about a many-electron system? The full wavefunction $\Psi$ lives in $3N$-dimensional space and is impossibly complicated; the density is a humble object you could plot on a sheet of paper. It seems too cheap to be enough. Hohenberg and Kohn proved that, for the ground state, it *is* enough. That result is what makes "density functional theory" a legitimate, exact reformulation of quantum mechanics — not merely a convenient approximation that happens to work. Everything later in the chapter (Kohn–Sham, exchange–correlation functionals, real calculations) stands on this guarantee. If you are meeting functionals for the first time, the notation guide in `../undergraduate/formula-reading-guide.md` may help.

!!! note "Why does this chapter exist?"
    Thomas and Fermi *guessed* that you could compute everything from the density alone (§5.1). It was a beautiful guess, but a guess. For thirty-five years, nobody knew whether the guess was *exactly* right or only approximately right. Was there some hidden information in the wavefunction that the density did not capture? Or did the density really determine everything — bond lengths, magnetism, the colour of a crystal — uniquely?
    
    In 1964 Hohenberg and Kohn proved: *yes, exactly right*. Two theorems in a four-page paper. Theorem I says: if you know the ground-state density of an electron system, you know *everything*. There is no hidden information; the density is sufficient. Theorem II says: among all candidate densities, the *correct* one is the one that minimises a certain (unknown but exact) energy functional. Together they justify, retrospectively, why Thomas and Fermi's idea was the right idea — and they tell us what we should be approximating.
    
    A useful analogy. Imagine you are looking at a photograph of a crowd. The photograph (a 2D projection) clearly throws away information about where each person was in depth. Yet for certain *summary statistics* (the average height, the colour of the most common jacket), the photograph is enough. HK is the statement that, for the *ground-state energy and every ground-state observable* of an electron system, the density is enough. We give the proof in detail because the argument is short, beautiful, and rests on nothing more than the textbook variational principle of quantum mechanics applied twice.

In 1964, Pierre Hohenberg and Walter Kohn published a short paper titled "Inhomogeneous Electron Gas". In four pages, they put the density-as-fundamental-variable idea on a rigorous footing. The paper contains two theorems — what we now call **Hohenberg–Kohn I** and **Hohenberg–Kohn II** — that, taken together, justify the entire programme of density functional theory. This section states and proves both, in detail.

!!! abstract "Key idea (Chapter 5.2)"
    *Hohenberg–Kohn I*: the ground-state density $n_0(\mathbf r)$ uniquely determines the external potential $v_\mathrm{ext}(\mathbf r)$ (up to a constant), and therefore the Hamiltonian, the wavefunction, and every ground-state observable. *Hohenberg–Kohn II*: there exists a universal functional $F[n]$ such that $E_v[n] = F[n] + \int v_\mathrm{ext} n\,\mathrm d\mathbf r$ is minimised at $n = n_0$ with minimum value $E_0$. Together: the density is a complete, variational basic variable. The Levy–Lieb constrained search makes this rigorous over a clean function space.

## 5.2.1 The setting

Throughout this section, fix the number of electrons $N$ and consider the Hamiltonian

$$
\hat{H} = \hat T + \hat V_{ee} + \hat V_\mathrm{ext},
\tag{5.10}
$$

where $\hat T = -\tfrac{1}{2}\sum_i \nabla_i^{2}$ is the kinetic energy operator, $\hat V_{ee} = \sum_{i<j}1/|\mathbf r_i - \mathbf r_j|$ is the electron–electron interaction, and

$$
\hat V_\mathrm{ext} = \sum_{i=1}^{N} v_\mathrm{ext}(\mathbf r_i)
$$

couples each electron to a one-body external potential. The first two operators are *universal* — they have the same form for every $N$-electron problem in atomic units. The third is what makes a hydrogen atom different from a diamond lattice: it specifies which external potential, that is, which arrangement of nuclei or applied field, we are studying.

Assume the ground state $|\Psi_0\rangle$ is non-degenerate. (The degenerate case requires a refinement we discuss in §5.2.5.) The ground-state density is

$$
n_0(\mathbf r) = N\int |\Psi_0(\mathbf r,\mathbf r_2,\dots,\mathbf r_N)|^{2}\,\mathrm d\mathbf r_2\cdots\mathrm d\mathbf r_N,
\tag{5.11}
$$

normalised so that $\int n_0\,\mathrm d\mathbf r = N$.

We will use the variational principle for the ground-state energy: for any normalised trial state $|\Psi\rangle$,

$$
\langle\Psi|\hat{H}|\Psi\rangle \;\geq\; E_0,
\tag{5.12}
$$

with equality iff $|\Psi\rangle = |\Psi_0\rangle$.

### The role of the variational principle

Equation (5.12) is the *only* non-trivial physics we will use in the proof of HK I. Everything else is algebraic bookkeeping. It is worth pausing to appreciate this: the deepest theorem of density functional theory follows from the textbook quantum-mechanical variational principle on wavefunctions, applied twice and added together. No new physical input is required — only a clever choice of trial states.

### Where do degeneracies sit in this picture?

If the ground state of $\hat H^{(1)}$ is $g$-fold degenerate, the inequality in Step 1 can become an equality when $|\Psi^{(2)}\rangle$ happens to lie in the degenerate ground-state manifold of $\hat H^{(1)}$ (which can occur even though it is the ground state of a *different* Hamiltonian). The clean fix is to work with ensembles or with the constrained-search formulation; we return to this in §5.2.5.

!!! note "In plain language"
    Strip away the formalism and the two theorems say something almost cheeky.

    **Theorem I — the density is a secret fingerprint.** Two different physical systems are different because they have different external potentials $v_\mathrm{ext}$ (different nuclei in different places). Theorem I says the ground-state density quietly *encodes* that potential: no two distinct potentials can share the same ground-state density. So once you know $n_0(\mathbf r)$, you know $v_\mathrm{ext}$, which fixes the full Hamiltonian $\hat H$, which fixes the wavefunction and therefore *every* property of the system. The cheap object turns out to carry all the information of the expensive one.

    **Theorem II — finding the ground state becomes a minimisation over densities.** There exists an energy functional $E[n]$ — a machine that eats a candidate density and returns a number — whose lowest value is the true ground-state energy $E_0$, reached exactly when you feed it the true ground-state density $n_0$. Any wrong density gives a *higher* energy. So instead of solving a $3N$-dimensional Schrödinger equation, we "merely" have to hunt for the density that minimises $E[n]$. This is the variational principle of DFT, now phrased in terms of the density rather than the wavefunction.

    The full statements and proofs follow. The plain-language versions above are exactly correct in spirit; the formal ones add the precise hypotheses (non-degeneracy, "up to a constant", admissibility) that make them airtight.

!!! note "Symbol guide"
    | Symbol | Meaning | Units (atomic units) |
    |---|---|---|
    | $n(\mathbf r)$ | electron (number) density at point $\mathbf r$; $n_0$ is the *true* ground-state density | electrons per $a_0^{3}$ |
    | $v_\mathrm{ext}(\mathbf r)$ | external one-body potential (the nuclei / applied field the electrons feel) | hartree |
    | $E[n]$ | total-energy functional: a number assigned to each candidate density | hartree |
    | $F[n]$ | universal functional $T[n]+V_{ee}[n]$ — same for every system, but not known in closed form | hartree |
    | $\Psi$ | the $N$-electron wavefunction; $\Psi_0$ is the ground state | (normalised; dimensionless after integration) |

## 5.2.2 Theorem I: the density determines the potential

!!! abstract "Theorem 5.2.1 (Hohenberg–Kohn I)"
    Let $v_\mathrm{ext}^{(1)}(\mathbf r)$ and $v_\mathrm{ext}^{(2)}(\mathbf r)$ be two external potentials, each giving rise to a non-degenerate $N$-electron ground state, that produce the same ground-state density $n_0(\mathbf r)$. Then $v_\mathrm{ext}^{(1)} - v_\mathrm{ext}^{(2)} = \mathrm{const}$.

In other (older) words, repeated for clarity: *let $v_\mathrm{ext}^{(1)}(\mathbf r)$ and $v_\mathrm{ext}^{(2)}(\mathbf r)$ be two external potentials that give the same ground-state density $n_0(\mathbf r)$. Then $v_\mathrm{ext}^{(1)} - v_\mathrm{ext}^{(2)} = \mathrm{const}$.*

In other words: the ground-state density determines the external potential **uniquely up to an additive constant**. Combined with the fact that $\hat T$ and $\hat V_{ee}$ are universal, this means that $n_0$ determines the entire Hamiltonian, hence every eigenstate, hence every property of the system.

### Proof by contradiction

The proof is short, elementary, and one of the most beautiful arguments in many-body physics. We give it in detail. The key inputs are the variational principle (5.12) and the fact that the *only* term in the Hamiltonian that distinguishes one external potential from another is the one-body integral $\int v\, n$ — which depends on $\Psi$ only through its density.

Suppose, for the sake of contradiction, that two external potentials $v^{(1)}$ and $v^{(2)}$ differ by more than a constant — $v^{(1)} - v^{(2)} \neq \mathrm{const}$ — yet produce the same ground-state density $n_0$. Let the corresponding Hamiltonians, ground states, and ground-state energies be

$$
\hat{H}^{(1)} = \hat T + \hat V_{ee} + \hat V^{(1)},
\qquad |\Psi^{(1)}\rangle,\;\;E^{(1)},
$$
$$
\hat{H}^{(2)} = \hat T + \hat V_{ee} + \hat V^{(2)},
\qquad |\Psi^{(2)}\rangle,\;\;E^{(2)}.
$$

The two Hamiltonians differ only in their one-body external potential. Their ground states are different — they are eigenstates of different operators — but by assumption they yield the same density.

**Step 1.** Use $|\Psi^{(2)}\rangle$ as a trial state for $\hat{H}^{(1)}$. By the variational principle (5.12),

$$
\langle\Psi^{(2)}|\hat{H}^{(1)}|\Psi^{(2)}\rangle > E^{(1)},
$$

with strict inequality because $|\Psi^{(2)}\rangle \neq |\Psi^{(1)}\rangle$ — they are non-degenerate ground states of different Hamiltonians, so they cannot be the same wavefunction (if they were, applying $\hat{H}^{(1)}-\hat{H}^{(2)} = \hat V^{(1)} - \hat V^{(2)}$ would give $(E^{(1)}-E^{(2)})|\Psi^{(2)}\rangle = (\hat V^{(1)}-\hat V^{(2)})|\Psi^{(2)}\rangle$, so $v^{(1)}-v^{(2)}$ would equal the constant $E^{(1)}-E^{(2)}$, contradicting our hypothesis).

!!! note "Why this step?"
    The variational principle (5.12) is the load-bearing wall of the entire HK proof. It says: among all normalised wavefunctions, the ground state achieves the minimum of $\langle\Psi|\hat H|\Psi\rangle$, with strict inequality for any other (non-degenerate) state. By using $|\Psi^{(2)}\rangle$ — the ground state of a *different* Hamiltonian — as a trial state for $\hat H^{(1)}$, we get a *strict* upper bound on $E^{(1)}$ from below, $\langle\Psi^{(2)}|\hat H^{(1)}|\Psi^{(2)}\rangle > E^{(1)}$. We then express this trial-state expectation value in a particularly useful form: it equals $E^{(2)}$ (the ground-state energy of the *other* Hamiltonian) plus a one-body correction that depends only on the density. This is the crux: the density-dependence of the correction is what lets the cancellation in Step 3 happen.

Now write

$$
\langle\Psi^{(2)}|\hat{H}^{(1)}|\Psi^{(2)}\rangle
= \langle\Psi^{(2)}|\hat{H}^{(2)}|\Psi^{(2)}\rangle
+ \langle\Psi^{(2)}|\hat V^{(1)} - \hat V^{(2)}|\Psi^{(2)}\rangle.
$$

The first term is $E^{(2)}$ since $|\Psi^{(2)}\rangle$ is the ground state of $\hat{H}^{(2)}$. The second is the expectation of a one-body operator, which depends only on the density of $|\Psi^{(2)}\rangle$, which is $n_0$ by assumption:

$$
\langle\Psi^{(2)}|\hat V^{(1)} - \hat V^{(2)}|\Psi^{(2)}\rangle
= \int n_0(\mathbf r)\big[v^{(1)}(\mathbf r) - v^{(2)}(\mathbf r)\big]\,\mathrm d\mathbf r.
$$

So

$$
E^{(2)} + \int n_0\big[v^{(1)} - v^{(2)}\big]\,\mathrm d\mathbf r \;>\; E^{(1)}.
\tag{5.13}
$$

**Step 2.** By exactly the symmetric argument, using $|\Psi^{(1)}\rangle$ as a trial state for $\hat{H}^{(2)}$,

$$
E^{(1)} + \int n_0\big[v^{(2)} - v^{(1)}\big]\,\mathrm d\mathbf r \;>\; E^{(2)}.
\tag{5.14}
$$

**Step 3.** Add (5.13) and (5.14). The density integrals cancel exactly:

$$
\int n_0\big[v^{(1)} - v^{(2)}\big]\,\mathrm d\mathbf r + \int n_0\big[v^{(2)} - v^{(1)}\big]\,\mathrm d\mathbf r = 0.
$$

What is left is

$$
E^{(1)} + E^{(2)} \;>\; E^{(1)} + E^{(2)},
$$

a contradiction. The only escape is to abandon the hypothesis that $v^{(1)} - v^{(2)}$ is non-constant. So $v^{(1)} - v^{(2)} = \mathrm{const}$. $\blacksquare$

!!! example "Sanity check: what if $v^{(1)}=v^{(2)}+\mathrm{const}$?"
    The proof is internally consistent: if the two potentials differ only by a constant $c$, then $\hat H^{(1)} = \hat H^{(2)} + c\,\hat N$, both Hamiltonians share the same eigenstates (up to a phase) and the same density, and the energies differ by $c\cdot N$. Steps 1 and 2 then degenerate to equalities (the strict inequality becomes equality since the wavefunctions coincide), and the contradiction in Step 3 evaporates. The "up to a constant" caveat in the theorem statement is therefore both necessary and trivial.

### Two corollaries worth stating explicitly

**Corollary 1.** *Every ground-state expectation value is a unique functional of the ground-state density.* This follows because $n_0 \to v_\mathrm{ext} \to \hat H \to |\Psi_0\rangle$, and any operator $\hat O$ has a unique ground-state expectation value $\langle\Psi_0|\hat O|\Psi_0\rangle = O[n_0]$. For example, the kinetic energy functional $T[n_0]$, the electron–electron interaction functional $V_{ee}[n_0]$, the pair correlation function $g(\mathbf r, \mathbf r')[n_0]$, and the dipole moment all become functionals of the density.

**Corollary 2.** *The map $n_0 \mapsto v_\mathrm{ext}$ is in general highly non-trivial.* Although it is mathematically unique (modulo a constant), there is no closed-form expression for it. Inverting the map — given $n_0$, find $v_\mathrm{ext}$ — is the subject of *inverse DFT*, an active research area used to construct exact reference data for benchmarking functionals. The Levy–Lieb procedure in §5.2.3 gives one constructive route, but not a practical one.

### What the theorem really says

There is a chain of one-to-one maps:

$$
v_\mathrm{ext}(\mathbf r) \;\xleftrightarrow{\;1\!:\!1\;}\; |\Psi_0\rangle \;\xleftrightarrow{\;1\!:\!1\;}\; n_0(\mathbf r).
\tag{5.15}
$$

The first map is "up to a constant in $v$ and a phase in $|\Psi\rangle$": adding a constant to $v$ shifts every eigenvalue but does not change eigenstates. The second map is what HK Theorem I proves. The density carries all the information of the wavefunction — for the ground state — even though it lives in a much smaller space.

!!! warning "Excited states are not determined"
    HK I applies to the *ground state* density. Two different excited eigenstates of two different Hamiltonians can in principle share the same density without contradiction. The fundamental theorem of DFT is a ground-state theorem. Excited-state extensions (TD-DFT, ensemble DFT) require additional machinery — see §5.6.

### A note on the inner product and the density picture

Several subtle points often confuse first-time readers of HK I. We list them.

1. *Why the density and not the wavefunction?* Both $|\Psi_0|^{2}$ in 3$N$-dimensional configuration space and $n_0$ in 3-dimensional space carry the same information *given the Hamiltonian*. The HK theorem says that, for ground states, even the *huge* compression $|\Psi_0|^{2}\to n_0$ is information-preserving — provided we know we are looking at a ground state of some Hamiltonian of the universal form $\hat T + \hat V_{ee} + \hat V_\mathrm{ext}$.
2. *Where is information stored in $n_0$?* The information lives in the spatial profile of $n_0$, in particular in its asymptotic decay (which encodes the ionisation potential, §5.3.4) and in the cusps at the nuclei (which encode $Z_\alpha$ via Kato's cusp condition $-\tfrac{1}{2n}\,\partial n/\partial r|_{\mathbf r=\mathbf R_\alpha} = Z_\alpha$). Roughly, the high-$r$ tail tells you the HOMO; the cusp at each nucleus tells you the nuclear charge and position; everything in between encodes the bonding pattern.
3. *Why is the theorem hard to use in practice?* Because the map $n_0 \to v_\mathrm{ext}$ is not given by any explicit formula. It exists, is unique, but is not computable in closed form. Inverse DFT — given $n$, find $v$ — uses iterative procedures (e.g., the van Leeuwen–Baerends construction) that are themselves tricky to converge.

### Consequence: every observable is a functional of $n_0$

Since $n_0$ determines $v_\mathrm{ext}$, hence $\hat{H}$, hence $|\Psi_0\rangle$, every ground-state expectation value is a functional of $n_0$. In particular the kinetic energy $T[n_0] = \langle\Psi_0[n_0]|\hat T|\Psi_0[n_0]\rangle$ and the electron–electron energy $V_{ee}[n_0] = \langle\Psi_0[n_0]|\hat V_{ee}|\Psi_0[n_0]\rangle$ are *exact* density functionals. They are universal — defined by the operators $\hat T$ and $\hat V_{ee}$, which do not depend on the system — but they are also unknown: the existence of $\Psi_0[n]$ does not give us a way to compute it.

## 5.2.3 Theorem II: the variational principle for the density

!!! abstract "Theorem 5.2.2 (Hohenberg–Kohn II)"
    Let $v_\mathrm{ext}(\mathbf r)$ be a given external potential with ground-state density $n_0(\mathbf r)$ and ground-state energy $E_0$. Define the total-energy functional $E_{v_\mathrm{ext}}[n] = F[n] + \int v_\mathrm{ext}(\mathbf r)\,n(\mathbf r)\,\mathrm d\mathbf r$ with $F[n] = T[n] + V_{ee}[n]$ the universal functional. Then for every admissible trial density $n$, $E_{v_\mathrm{ext}}[n] \geq E_0$, with equality iff $n = n_0$.

For the longer historical statement and proof, read on. *Define the total energy functional*

$$
E_{v_\mathrm{ext}}[n] \;=\; F[n] \;+\; \int v_\mathrm{ext}(\mathbf r)\,n(\mathbf r)\,\mathrm d\mathbf r,
\tag{5.16}
$$

*where $F[n] = T[n] + V_{ee}[n]$ is the universal functional defined above. Then for any admissible trial density $n(\mathbf r)$,*

$$
E_{v_\mathrm{ext}}[n] \;\geq\; E_0,
$$

*with equality iff $n = n_0$, the true ground-state density.*

The functional $E_{v_\mathrm{ext}}[n]$ is minimised by the true ground-state density and the minimum equals the ground-state energy. This is the variational principle that justifies *minimising over densities* — exactly the kind of programme Thomas–Fermi tried to implement, but now with the assurance that, if we knew $F[n]$, the answer would be exact.

### The original proof, and what "admissible" means

HK's original proof restricted attention to densities that arise as ground-state densities of *some* external potential — so-called **$v$-representable** densities. For such a density $n$ there exists a unique $v$ (by HK I) and hence a unique $|\Psi[n]\rangle$ and a unique value $F[n] = \langle\Psi[n]|\hat T + \hat V_{ee}|\Psi[n]\rangle$. The variational principle on wavefunctions then gives, for the system with potential $v_\mathrm{ext}$,

$$
\langle\Psi[n]|\hat{H}|\Psi[n]\rangle = F[n] + \int v_\mathrm{ext} n\,\mathrm d\mathbf r \;\geq\; E_0,
$$

with equality only when $|\Psi[n]\rangle = |\Psi_0\rangle$, i.e., when $n = n_0$.

The trouble is that the space of $v$-representable densities is awkward: not every reasonable density (a smooth, positive function integrating to $N$) is the ground-state density of some $v$. In fact pathological counter-examples exist. This makes the variational principle restricted to $v$-representable densities theoretically inelegant and practically unusable, because in any numerical scheme we want to vary $n$ over a larger, simpler space.

### The constrained-search reformulation (Levy, Lieb)

Mel Levy in 1979, and independently Elliott Lieb in 1983, gave a cleaner statement that bypasses $v$-representability entirely. Define the universal functional by a *constrained search* over all antisymmetric $N$-electron wavefunctions $\Psi$ that yield the density $n$:

$$
F_L[n] \;\equiv\; \min_{\Psi \to n}\;\langle\Psi|\hat T + \hat V_{ee}|\Psi\rangle.
\tag{5.17}
$$

That is, among all $N$-electron antisymmetric wavefunctions whose density is the prescribed $n$, find the one with the lowest kinetic-plus-interaction expectation value. The minimiser exists for any **$N$-representable** density — any non-negative function $n(\mathbf r)$ with $\int n = N$, $\int|\nabla\sqrt n|^{2} < \infty$. This class is much larger and easier to characterise than the $v$-representable class.

!!! note "Why this step?"
    Restricting to $v$-representable densities is the *original* HK strategy: it guarantees that for each admissible $n$, there is a well-defined wavefunction $\Psi[n]$ from which to compute $F[n]$. The restriction is natural in the sense that the physical answer (the true ground-state density) is automatically $v$-representable. The cost is mathematical: minimising over a strangely-shaped set is bad for variational calculus. Levy and Lieb's reformulation enlarges the search domain and pays the price by an inner minimisation.

### Levy–Lieb proof of HK II

Let $n_0$ be the ground-state density of the system with potential $v_\mathrm{ext}$ and ground-state wavefunction $\Psi_0$. We must show

$$
F_L[n] + \int v_\mathrm{ext} n\,\mathrm d\mathbf r \;\geq\; E_0 \quad\text{for all }n,\text{ with equality at }n = n_0.
$$

**Step 1 (inequality).** Fix any $N$-representable density $n$. By definition (5.17), there is a wavefunction $\Psi_n$ with density $n$ achieving the constrained minimum:

$$
F_L[n] = \langle\Psi_n|\hat T + \hat V_{ee}|\Psi_n\rangle.
$$

Compute the full energy expectation value of $\Psi_n$ in the system of interest:

$$
\langle\Psi_n|\hat{H}|\Psi_n\rangle = \langle\Psi_n|\hat T + \hat V_{ee}|\Psi_n\rangle + \langle\Psi_n|\hat V_\mathrm{ext}|\Psi_n\rangle = F_L[n] + \int v_\mathrm{ext} n\,\mathrm d\mathbf r.
$$

By the variational principle on wavefunctions, $\langle\Psi_n|\hat{H}|\Psi_n\rangle \geq E_0$. Therefore

$$
F_L[n] + \int v_\mathrm{ext} n\,\mathrm d\mathbf r \;\geq\; E_0. \tag{5.18}
$$

!!! note "Why this step?"
    The constrained-search minimum $\Psi_n$ may not be the ground state of any external potential — it is whatever wavefunction makes $\langle\hat T + \hat V_{ee}\rangle$ smallest while reproducing $n$. So $\Psi_n$ is in general *not* an eigenstate of $\hat H$, but the variational principle still gives an upper bound on its energy expectation value, $\langle\Psi_n|\hat H|\Psi_n\rangle \geq E_0$. This is the key inequality.

**Step 2 (equality at $n=n_0$).** The true ground-state wavefunction $\Psi_0$ has density $n_0$, so it is a candidate in the constrained search (5.17) for $F_L[n_0]$:

$$
F_L[n_0] \;\leq\; \langle\Psi_0|\hat T + \hat V_{ee}|\Psi_0\rangle = E_0 - \int v_\mathrm{ext} n_0\,\mathrm d\mathbf r.
$$

Rearranging, $F_L[n_0] + \int v_\mathrm{ext} n_0\,\mathrm d\mathbf r \leq E_0$. Combined with the opposite inequality from Step 1, applied at $n = n_0$,

$$
F_L[n_0] + \int v_\mathrm{ext} n_0\,\mathrm d\mathbf r = E_0.
$$

This proves both the variational inequality and equality at the true ground-state density. $\blacksquare$

### A worked illustration: HK II for two electrons in a harmonic trap

To make the abstract proof concrete, consider two non-interacting electrons in a 1D harmonic trap $v_\mathrm{ext}(x) = \tfrac{1}{2}\omega^{2}x^{2}$. The exact ground state is a Slater determinant of the lowest two single-particle eigenstates of the harmonic oscillator, $\phi_0(x)$ and $\phi_1(x)$ with opposite spins (in fact, both spin-up and spin-down filling $\phi_0$ — for two electrons of opposite spin both occupy the same spatial orbital). The density is $n_0(x) = 2|\phi_0(x)|^{2} = 2\sqrt{\omega/\pi}\,\mathrm e^{-\omega x^{2}}$, normalised so $\int n_0\,\mathrm dx = 2$.

Now suppose we hand HK II a *different* trial density: $n_\alpha(x) = 2\sqrt{\alpha/\pi}\,\mathrm e^{-\alpha x^{2}}$, with $\alpha\neq \omega$. For non-interacting electrons, the Levy–Lieb minimum is achieved by the Slater determinant of orbitals $\phi^{(\alpha)}_0$ with the corresponding Gaussian — this gives $T_s[n_\alpha] = \tfrac{1}{2}\alpha$ and $\int v_\mathrm{ext}\,n_\alpha\,\mathrm dx = \omega^{2}/(2\alpha)$. The trial energy is

$$
E_v[n_\alpha] = \tfrac{1}{2}\alpha + \frac{\omega^{2}}{2\alpha}.
$$

Minimising over $\alpha$ recovers $\alpha = \omega$ as expected, with minimum value $E_0 = \omega$. For any $\alpha\neq\omega$, $E_v[n_\alpha]>\omega$, illustrating the strict HK II inequality. This toy example also shows that the *minimum* over trial densities is what gives the ground-state energy — exactly the variational structure HK II guarantees.

!!! warning "Common misunderstanding"
    Three traps catch almost every newcomer to the Hohenberg–Kohn theorems.

    - **"In principle" is doing enormous work.** HK guarantees that the universal functional $F[n]$ *exists* and is unique. It does **not** tell us how to write it down. You cannot open the HK paper, find a formula for $F[n]$, and start computing. The existence is exact; the recipe is missing — and supplying useful approximations to the missing piece has occupied the field for sixty years (this is §5.4).
    - **The variational principle is over densities, not wavefunctions.** Ordinary quantum-mechanical variation searches over trial *wavefunctions* $\Psi$. HK II searches over trial *densities* $n(\mathbf r)$. That is the whole point — the density lives in 3 dimensions, the wavefunction in $3N$. Confusing the two collapses the entire advantage of DFT.
    - **The theorems are about the GROUND state only.** Both theorems are statements about the ground state. They do *not* say that an excited-state density determines its potential, nor that excited-state energies follow from $E[n]$. Excited states need extra machinery (TD-DFT, ensemble DFT; see §5.6).

## 5.2.4 The universal functional $F[n]$ — exact, and unknowable

We have proved that the functional

$$
F[n] \;=\; T[n] \;+\; V_{ee}[n]
\tag{5.19}
$$

exists, is *universal* (its definition does not refer to any external potential — the same $F$ is used for hydrogen, for copper, for water), and combined with $\int v_\mathrm{ext} n$ yields a variational principle whose minimum is the exact ground-state energy.

Why, then, do we not just minimise $E_{v_\mathrm{ext}}[n] = F[n] + \int v_\mathrm{ext} n$ over densities and report the answer? Because we do not know $F[n]$ explicitly. The constrained-search definition (5.17) requires that, for each candidate $n$, we *minimise over all $N$-electron wavefunctions yielding $n$* — a problem at least as hard as the original Schrödinger equation we were trying to escape.

The pieces of $F[n]$ are no easier. The kinetic energy is

$$
T[n] = \min_{\Psi\to n}\langle\Psi|\hat T|\Psi\rangle,
$$

a non-local, non-trivial functional. We saw in §5.1 that even the leading-order local approximation $T_\mathrm{TF} \propto \int n^{5/3}$, while qualitatively reasonable, is too poor to bind a molecule. The electron–electron functional

$$
V_{ee}[n] = \langle\Psi[n]|\hat V_{ee}|\Psi[n]\rangle
$$

contains the classical Hartree piece (5.7) plus an **exchange–correlation** piece encoding the antisymmetry of $\Psi$ and the correlation hole around each electron. We cannot write either of these in closed form.

!!! warning "HK does not give us $F[n]$"
    A common misconception is that the Hohenberg–Kohn theorems supply a usable energy functional. They do not. They prove that an exact $F[n]$ *exists*; they say nothing about its form. Constructing accurate approximations to $F[n]$ — or rather to its still-mysterious piece, the exchange–correlation energy — has been the entire research enterprise of DFT for sixty years (§5.4).

## 5.2.4a $v$-representability and $N$-representability

A pair of technical issues lurks beneath the surface of the Hohenberg–Kohn theorems. They are usually glossed over in introductory treatments but matter for both rigorous mathematical statements and certain practical questions (orbital-free DFT, density-functional embedding, exact-exchange methods).

### $v$-representability

A density $n(\mathbf r)$ is **$v$-representable** if there exists *some* external potential $v(\mathbf r)$ for which $n$ is the ground-state density of the interacting $N$-electron Hamiltonian $\hat T + \hat V_{ee} + \sum_i v(\mathbf r_i)$. The HK proof uses this implicitly: when we wrote "let the corresponding ground state be $|\Psi^{(1)}\rangle$ with density $n_0$", we assumed $n_0$ comes from *some* well-defined ground state of *some* potential.

The trouble is that not every nice-looking density is $v$-representable. Known counter-examples include:

- *Densities with cusps too sharp or too smooth* relative to the Kato cusp condition $\partial n/\partial r|_{r=0} = -2Z\,n(0)$ at a nucleus of charge $Z$.
- *Degenerate ground-state densities* of certain pathological Hamiltonians.
- *Non-integer densities*: $\int n\,\mathrm d\mathbf r$ must be an integer for an interacting ground state of a fixed-$N$ Hamiltonian — yet variational schemes naturally allow fractional total charges (mixed states of $N$ and $N\pm 1$).
- *Densities with isolated zeros*: a density vanishing at some interior point of space cannot be the ground state of a finite local potential.

Englisch and Englisch (1983) exhibited explicit non-$v$-representable smooth positive densities integrating to an integer. The set of $v$-representable densities is therefore neither convex nor closed in any natural topology, which makes the original HK variational principle awkward as a rigorous foundation.

!!! note "Why this step?"
    Without $v$-representability, the original HK derivation of HK II goes through but only on a restricted (and ill-characterised) set of "admissible" densities. To do variational calculus we need a domain — a space of allowed inputs — and that domain should be large, simple, and closed under the operations we want to perform. The fix, due to Levy and Lieb, is to enlarge the domain to the class of $N$-representable densities, which is much better behaved.

### $N$-representability

A density $n(\mathbf r)$ is **$N$-representable** if there exists some antisymmetric $N$-electron wavefunction $\Psi(\mathbf r_1,\dots,\mathbf r_N)$ whose marginal density is $n$. The characterisation is much cleaner than for $v$-representability. Gilbert (1975) and Harriman (1981) proved:

!!! note "Gilbert–Harriman conditions"
    A non-negative function $n(\mathbf r)$ is $N$-representable by an antisymmetric pure-state $N$-electron wavefunction if and only if
    1. $n(\mathbf r) \geq 0$ everywhere,
    2. $\int n(\mathbf r)\,\mathrm d\mathbf r = N$ (an integer),
    3. $\int |\nabla\sqrt{n(\mathbf r)}|^{2}\,\mathrm d\mathbf r < \infty$ (the so-called *kinetic-energy bound*).

Condition (3) ensures the von Weizsäcker kinetic energy is finite and is the only non-trivial constraint; it rules out densities with too-sharp features. Harriman gave an explicit construction of an $N$-orbital Slater determinant that exactly reproduces any density satisfying (1)–(3), which proves the conditions are sufficient.

The $N$-representable class is convex and closed under reasonable topologies, which makes it the natural domain for variational principles. The Levy–Lieb constrained search below works on this set.

## 5.2.4b The Levy–Lieb constrained search, in detail

We restate the Levy–Lieb construction (already sketched above) with care, because it is the conceptual foundation of *every* modern formulation of DFT — including the Kohn–Sham construction of §5.3 and the generalised Kohn–Sham used by hybrid functionals.

For each $N$-representable density $n$, define

$$
F_L[n] \;\equiv\; \min_{\Psi \to n}\;\langle\Psi|\hat T + \hat V_{ee}|\Psi\rangle.
\tag{5.19a}
$$

The minimum is over all antisymmetric $N$-electron wavefunctions $\Psi$ with $|\Psi|^{2}$-marginal equal to $n$. The total-energy functional is

$$
E_v[n] \;=\; F_L[n] \;+\; \int v(\mathbf r)\,n(\mathbf r)\,\mathrm d\mathbf r,
\tag{5.19b}
$$

and the ground-state energy is

$$
E_0 \;=\; \min_{n}\,E_v[n] \;=\; \min_{n}\Big\{ F_L[n] + \int v n\,\mathrm d\mathbf r\Big\}.
\tag{5.19c}
$$

!!! note "Why this step?"
    The genius of the Levy–Lieb formulation is its *two-stage* variational structure. The inner minimisation in (5.19a) handles the wavefunction part — the hard, $3N$-dimensional problem — but produces a number depending only on $n$. The outer minimisation in (5.19c) is over the much simpler $3$-dimensional density. By splitting the problem this way, Levy and Lieb show that the HK variational principle is *equivalent* to the standard Ritz variational principle on wavefunctions, but with the optimisation over $\Psi$ broken into pieces classified by the density they produce.

### Properties of $F_L$

The Levy–Lieb functional has several non-obvious properties:

- *Universality.* $F_L[n]$ does not depend on the external potential. The same $F_L$ applies to a hydrogen atom and to copper metal; it depends only on $\hat T + \hat V_{ee}$, the universal operators of an $N$-electron problem.
- *Convexity.* $F_L$ is convex: $F_L[\tfrac{1}{2}(n_1 + n_2)] \leq \tfrac{1}{2}(F_L[n_1] + F_L[n_2])$. This follows from the linearity of the constraint $\Psi\to n$ and the convexity of the kinetic energy in $\Psi$.
- *Lower semi-continuity.* The functional is lower semi-continuous in appropriate topologies, ensuring the outer minimum in (5.19c) is attained.

The Lieb functional $F_L$ thus realises the dream of a clean variational principle on densities, defined on a clean function space, with all the abstract pathology of $v$-representability moved out of the way. Of course $F_L$ remains unknown in closed form — its constrained search is at least as hard as the original Schrödinger equation — but it gives a *well-posed* problem to approximate.

## 5.2.5 Refinements and subtleties

### Degenerate ground states

If the ground state is $g$-fold degenerate, the proof of HK I needs amending: the conclusion "$|\Psi^{(2)}\rangle \neq |\Psi^{(1)}\rangle$" can fail when both are members of a degenerate manifold of the same Hamiltonian. The correct statement is that *the set of ground-state densities* determines the external potential up to a constant — equivalently, one works with *ensembles* of degenerate ground states and the corresponding ensemble density. Levy's constrained search is naturally extended to mixed states $\hat\rho \to n$, giving an analogous variational principle.

### The Hohenberg–Kohn–Sham bridge

A natural objection at this point is: if $F[n]$ is exact but unknowable, how does DFT do anything in practice? The answer, supplied by Kohn and Sham one year after HK, is that we do *not* need an explicit $F[n]$. We need only an *approximation* to the small remainder $E_{xc}[n]$ defined as the difference between $F[n]$ and a *known*, exactly-computable piece (the non-interacting kinetic energy plus the classical Hartree energy). This is the subject of §5.3.

In particular, the HK theorems do not by themselves tell us *how accurate* an approximate functional needs to be. They tell us only that an exact one exists. The empirical fact — discovered through fifty years of trying — is that even crude approximations to $E_{xc}$ (like LDA) give qualitatively correct chemistry, and modest improvements (like PBE or SCAN) give quantitative accuracy. Why this should be so is itself a deep question: it is essentially the statement that the exchange–correlation energy is *near-local* in the density, i.e. it does not depend strongly on the density far from the point of interest. That fact is the foundation of all semi-local DFT.

### Spin-DFT

When external magnetic fields or spin polarisation are relevant, one promotes the density to the pair $(n_\uparrow,n_\downarrow)$ — or equivalently $(n, m_z)$ where $n = n_\uparrow + n_\downarrow$ and $m_z = n_\uparrow - n_\downarrow$. The Hohenberg–Kohn theorems generalise: the pair $(n_\uparrow,n_\downarrow)$ determines the pair $(v_\mathrm{ext}, B_z)$ up to constants. For non-magnetic systems the unpolarised theory suffices, but most modern functionals are written in their spin-polarised form because that is the more general expression.

### Finite temperature

Mermin (1965) extended the theorems to finite temperature: at temperature $T$, the equilibrium density determines the external potential, and a grand-canonical free-energy functional $\Omega[n]$ is minimised at the equilibrium density. This is the basis of "finite-$T$ DFT" used for warm dense matter.

### Non-collinear and current-density DFT

When spin–orbit coupling is important or magnetic moments rotate in space (skyrmions, helical magnets, frustrated antiferromagnets), the spin density becomes a $2\times 2$ Hermitian matrix $\hat n_{\alpha\beta}(\mathbf r)$ in spin space, equivalent to the pair $(n(\mathbf r), \mathbf m(\mathbf r))$ where $\mathbf m$ is the vector magnetisation. The corresponding HK-style theorem (Vignale–Rasolt, 1987) replaces the external potential by the pair $(v_\mathrm{ext}, \mathbf B_\mathrm{ext})$ and the density by the pair $(n, \mathbf m)$. For systems with orbital currents (induced magnetic fields, the integer quantum Hall effect), one needs **current-density functional theory** (CDFT), where the basic variable is enriched to include the paramagnetic current density $\mathbf j_p$. These extensions are technically demanding but follow the same logical pattern as HK I: a constructive uniqueness proof by variational contradiction.

### What about excited states?

HK theorems are emphatically ground-state theorems. The cleanest generalisation to excited states is **time-dependent DFT** (TD-DFT), based on the Runge–Gross theorem (1984): for fixed initial state $\Psi(0)$, the time-dependent density $n(\mathbf r,t)$ determines the time-dependent external potential $v_\mathrm{ext}(\mathbf r,t)$ up to a purely time-dependent constant. We discuss TD-DFT briefly in §5.6.

### Extended-system DFT and the Hohenberg–Kohn theorem with PBCs

For solids — the bread and butter of materials simulation — the system is infinite and periodic. The Hohenberg–Kohn theorems extend to this setting but require careful definition. The "external potential" is now a periodic function $v_\mathrm{ext}(\mathbf r)$ with the lattice periodicity, and the ground-state density is similarly periodic. The total energy per unit cell becomes the relevant variational object, and the uniqueness statement is "the periodic density determines the periodic potential up to a constant *per unit cell*". The proof goes through almost verbatim. Practically, this is why DFT codes for crystals can simply impose Born–von Kármán periodic boundary conditions and run the same machinery developed for finite systems.

### Janak's theorem and fractional occupations

Janak (1978) proved an identity of considerable practical and conceptual importance, which we shall use repeatedly in §5.3. Consider the Kohn–Sham total energy as a function of orbital occupations $f_i \in [0,1]$ (the spin-summed convention has $f_i\in[0,2]$, but to avoid factors of two we work with $f_i\in[0,1]$ and one spin channel here). Then

$$
\boxed{\;\;\frac{\partial E_\mathrm{KS}[\{f_i\}]}{\partial f_i} \;=\; \varepsilon_i,\;\;}
\tag{5.19d}
$$

where $\varepsilon_i$ is the $i$-th KS eigenvalue. The derivation is short and instructive. The KS density is $n(\mathbf r) = \sum_i f_i |\phi_i(\mathbf r)|^{2}$, and the KS total energy decomposes as

$$
E_\mathrm{KS} = \sum_i f_i\langle\phi_i| -\tfrac{1}{2}\nabla^{2}|\phi_i\rangle + \int v_\mathrm{ext} n\,\mathrm d\mathbf r + U_H[n] + E_{xc}[n].
$$

Differentiating with respect to $f_i$ and using $\partial n/\partial f_i = |\phi_i|^{2}$,

$$
\frac{\partial E_\mathrm{KS}}{\partial f_i} = \langle\phi_i| -\tfrac{1}{2}\nabla^{2}|\phi_i\rangle + \int (v_\mathrm{ext} + v_H + v_{xc})\,|\phi_i|^{2}\,\mathrm d\mathbf r = \langle\phi_i|\hat h_\mathrm{KS}|\phi_i\rangle = \varepsilon_i,
$$

since $\hat h_\mathrm{KS}\phi_i = \varepsilon_i\phi_i$. The orbitals themselves depend on $f_i$ through the self-consistent potential, but their contribution to $\partial E/\partial f_i$ vanishes by the variational property of the KS orbitals (this is the analogue of the Hellmann–Feynman theorem).

!!! note "Why this step?"
    Janak's theorem is the deepest justification for the practice of *partial occupations* (Fermi smearing) for metals — see §5.5. It also implies the *piecewise-linearity* of $E[N]$ for the exact functional: if we allow the total electron number $N$ to vary continuously between integers, the exact energy must be a *straight line* between consecutive integers, because $\partial E/\partial N$ equals the chemical potential, which jumps by $\Delta_{xc}$ (the derivative discontinuity, §5.6) at each integer. Approximate functionals violate this piecewise linearity, and the resulting *curvature* of $E[N]$ is one of the most diagnostic measures of self-interaction error in modern functional development.

!!! example "Janak in action: ionisation potential of hydrogen"
    For the hydrogen atom with exact $E_{xc}$, sweeping the occupation $f$ of the 1s orbital from 1 to 0 traces out the curve $E(f) = E(1) - I\cdot(1-f)$ exactly, where $I = 0.5\;\text{Ha} = 13.6\;\text{eV}$ is the ionisation potential, and the slope at any $f\in(0,1)$ equals $-I$. With approximate functionals like LDA the curve is concave-up — $E(f)$ bows below the straight line — so the slope at $f=1$, which is the LDA HOMO eigenvalue, underestimates $I$. The LDA HOMO of hydrogen sits at roughly $-7\;\text{eV}$, against the exact $-13.6\;\text{eV}$: an error of about $6$–$7\;\text{eV}$, almost half the binding energy, and a direct fingerprint of the self-interaction error discussed in §5.4.

### Asymptotic behaviour of the exact KS potential

A surprisingly informative consequence of HK is the *asymptotic decay* of the exact KS potential. For a neutral finite system, the exact $v_\mathrm{KS}(\mathbf r)$ must satisfy $v_\mathrm{KS}(\mathbf r) \to -1/|\mathbf r|$ as $|\mathbf r|\to\infty$, because the outermost electron sees the unscreened nuclear charge (the inner $N-1$ electrons fully screen the other $N-1$ protons). The Hartree potential gives $-1/|\mathbf r|\cdot(N-1)/N$ for a system of $N$ electrons asymptotically; the exchange–correlation potential must therefore contribute the missing $-1/(N|\mathbf r|)$ to leading order. LDA and GGA potentials, however, decay *exponentially* with $n(\mathbf r)$ — they violate the $-1/|\mathbf r|$ asymptote by orders of magnitude at large $r$, with the result that Rydberg states, the HOMO, and the long-range tails of orbitals are all wrong. This is a structural failure that no semi-local functional can repair, and one of the motivations for range-separated hybrids (§5.4) and for asymptotically-corrected potentials.

!!! question "Check yourself"
    1. Theorem I says the ground-state density $n_0$ determines the external potential up to a constant. Starting from $n_0$, list the chain of objects you can then reconstruct, ending in "every ground-state property".
    2. In Theorem II, *what* quantity is minimised, and *over what set* is the minimisation performed? What is the value at the minimum?
    3. If the Hohenberg–Kohn theorems prove an exact energy functional exists, why can we still not compute molecular energies exactly from them?
    4. True or false: HK II lets you variationally compute the energy of the first excited state by minimising $E[n]$ over excited-state densities.

??? success "Answers"
    1. $n_0 \;\rightarrow\; v_\mathrm{ext}$ (Theorem I, up to a constant) $\;\rightarrow\; \hat H = \hat T + \hat V_{ee} + \hat V_\mathrm{ext}$ (since $\hat T$ and $\hat V_{ee}$ are universal) $\;\rightarrow\; \Psi_0$ (solve the Schrödinger equation) $\;\rightarrow\;$ every ground-state observable $\langle\Psi_0|\hat O|\Psi_0\rangle$.
    2. The total-energy functional $E_{v_\mathrm{ext}}[n] = F[n] + \int v_\mathrm{ext}\,n\,\mathrm d\mathbf r$ is minimised. The minimisation is over admissible trial *densities* $n(\mathbf r)$ (those that are $N$-representable: $n\geq 0$, $\int n = N$, finite $\int|\nabla\sqrt n|^{2}$). The value at the minimum is the true ground-state energy $E_0$, attained at $n = n_0$.
    3. Because HK proves the universal functional $F[n]$ *exists* but gives no closed-form expression for it. Its constrained-search definition is as hard as the original many-body Schrödinger problem. We can only minimise $E[n]$ if we have a usable $F[n]$, and we do not — we have approximations (§5.4).
    4. False. HK II is a *ground-state* variational principle; minimising $E[n]$ returns the ground-state energy and density, not excited states. Excited states require separate theory (§5.6).

## 5.2.6 Summary

The Hohenberg–Kohn theorems are an existence proof. They tell us:

1. The ground-state density $n_0(\mathbf r)$ uniquely determines the external potential $v_\mathrm{ext}(\mathbf r)$ (up to a constant), and hence every property of the system.
2. There is a universal functional $F[n] = T[n] + V_{ee}[n]$, and the total energy functional $E_{v_\mathrm{ext}}[n] = F[n] + \int v_\mathrm{ext} n$ is minimised by the true ground-state density.

What we do *not* yet have:

- An explicit form for $F[n]$. Both $T[n]$ and $V_{ee}[n]$ are defined by constrained searches no more tractable than the original many-body problem.
- A practical algorithm.

Kohn and Sham, one year after Hohenberg and Kohn, supplied the second. They split $F[n]$ into a part we can compute exactly via auxiliary single-particle orbitals (the non-interacting kinetic energy $T_s$ and the classical Hartree energy), and a small remainder — the exchange–correlation energy $E_{xc}[n]$ — to be approximated. That is §5.3.

### Summary of §5.2 — what to remember in 3 months

- **HK I (Theorem 5.2.1)**: ground-state density uniquely determines the external potential, hence the Hamiltonian, hence every ground-state observable. Proof: variational principle applied twice + contradiction.
- **HK II (Theorem 5.2.2)**: there is a universal functional $F[n] = T[n] + V_{ee}[n]$, and the total-energy functional $E_v[n] = F[n] + \int v\,n$ is minimised by the true ground-state density.
- **Levy–Lieb constrained search**: $F_L[n] = \min_{\Psi\to n}\langle\Psi|\hat T + \hat V_{ee}|\Psi\rangle$ extends $F[n]$ to all $N$-representable densities, avoiding the $v$-representability headaches of the original HK proof.
- **$N$-representability**: a density $n\geq 0$ with $\int n = N$ and $\int|\nabla\sqrt n|^{2}<\infty$ is $N$-representable (Gilbert–Harriman conditions). Much cleaner than $v$-representability.
- **Janak's theorem**: $\partial E/\partial f_i = \varepsilon_i$ — the $i$-th KS eigenvalue equals the derivative of total energy with respect to occupation.
- **What HK does NOT give**: an explicit form for $F[n]$, or for $E_{xc}[n]$. Existence proof only. The recipe is in §5.3 and §5.4.

!!! note "Remark: Nobel Prize"
    Walter Kohn shared the 1998 Nobel Prize in Chemistry with John Pople "for his development of the density-functional theory". The citation specifically called out the Hohenberg–Kohn paper (1964) and the Kohn–Sham paper (1965). Pierre Hohenberg, who was a postdoc at the time, was famously not included; he later remarked that being absent from the prize was less important than the work being recognised.
