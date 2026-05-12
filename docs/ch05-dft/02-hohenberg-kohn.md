# 5.2 The Hohenberg–Kohn Theorems

In 1964, Pierre Hohenberg and Walter Kohn published a short paper titled "Inhomogeneous Electron Gas". In four pages, they put the density-as-fundamental-variable idea on a rigorous footing. The paper contains two theorems — what we now call **Hohenberg–Kohn I** and **Hohenberg–Kohn II** — that, taken together, justify the entire programme of density functional theory. This section states and proves both, in detail.

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

## 5.2.2 Theorem I: the density determines the potential

**Hohenberg–Kohn Theorem I.** *Let $v_\mathrm{ext}^{(1)}(\mathbf r)$ and $v_\mathrm{ext}^{(2)}(\mathbf r)$ be two external potentials that give the same ground-state density $n_0(\mathbf r)$. Then $v_\mathrm{ext}^{(1)} - v_\mathrm{ext}^{(2)} = \mathrm{const}$.*

In other words: the ground-state density determines the external potential **uniquely up to an additive constant**. Combined with the fact that $\hat T$ and $\hat V_{ee}$ are universal, this means that $n_0$ determines the entire Hamiltonian, hence every eigenstate, hence every property of the system.

### Proof by contradiction

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

### What the theorem really says

There is a chain of one-to-one maps:

$$
v_\mathrm{ext}(\mathbf r) \;\xleftrightarrow{\;1\!:\!1\;}\; |\Psi_0\rangle \;\xleftrightarrow{\;1\!:\!1\;}\; n_0(\mathbf r).
\tag{5.15}
$$

The first map is "up to a constant in $v$ and a phase in $|\Psi\rangle$": adding a constant to $v$ shifts every eigenvalue but does not change eigenstates. The second map is what HK Theorem I proves. The density carries all the information of the wavefunction — for the ground state — even though it lives in a much smaller space.

!!! warning "Excited states are not determined"
    HK I applies to the *ground state* density. Two different excited eigenstates of two different Hamiltonians can in principle share the same density without contradiction. The fundamental theorem of DFT is a ground-state theorem. Excited-state extensions (TD-DFT, ensemble DFT) require additional machinery — see §5.6.

### Consequence: every observable is a functional of $n_0$

Since $n_0$ determines $v_\mathrm{ext}$, hence $\hat{H}$, hence $|\Psi_0\rangle$, every ground-state expectation value is a functional of $n_0$. In particular the kinetic energy $T[n_0] = \langle\Psi_0[n_0]|\hat T|\Psi_0[n_0]\rangle$ and the electron–electron energy $V_{ee}[n_0] = \langle\Psi_0[n_0]|\hat V_{ee}|\Psi_0[n_0]\rangle$ are *exact* density functionals. They are universal — defined by the operators $\hat T$ and $\hat V_{ee}$, which do not depend on the system — but they are also unknown: the existence of $\Psi_0[n]$ does not give us a way to compute it.

## 5.2.3 Theorem II: the variational principle for the density

**Hohenberg–Kohn Theorem II.** *Define the total energy functional*

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

**Step 2 (equality at $n=n_0$).** The true ground-state wavefunction $\Psi_0$ has density $n_0$, so it is a candidate in the constrained search (5.17) for $F_L[n_0]$:

$$
F_L[n_0] \;\leq\; \langle\Psi_0|\hat T + \hat V_{ee}|\Psi_0\rangle = E_0 - \int v_\mathrm{ext} n_0\,\mathrm d\mathbf r.
$$

Rearranging, $F_L[n_0] + \int v_\mathrm{ext} n_0\,\mathrm d\mathbf r \leq E_0$. Combined with the opposite inequality from Step 1, applied at $n = n_0$,

$$
F_L[n_0] + \int v_\mathrm{ext} n_0\,\mathrm d\mathbf r = E_0.
$$

This proves both the variational inequality and equality at the true ground-state density. $\blacksquare$

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

## 5.2.5 Refinements and subtleties

### Degenerate ground states

If the ground state is $g$-fold degenerate, the proof of HK I needs amending: the conclusion "$|\Psi^{(2)}\rangle \neq |\Psi^{(1)}\rangle$" can fail when both are members of a degenerate manifold of the same Hamiltonian. The correct statement is that *the set of ground-state densities* determines the external potential up to a constant — equivalently, one works with *ensembles* of degenerate ground states and the corresponding ensemble density. Levy's constrained search is naturally extended to mixed states $\hat\rho \to n$, giving an analogous variational principle.

### Spin-DFT

When external magnetic fields or spin polarisation are relevant, one promotes the density to the pair $(n_\uparrow,n_\downarrow)$ — or equivalently $(n, m_z)$ where $n = n_\uparrow + n_\downarrow$ and $m_z = n_\uparrow - n_\downarrow$. The Hohenberg–Kohn theorems generalise: the pair $(n_\uparrow,n_\downarrow)$ determines the pair $(v_\mathrm{ext}, B_z)$ up to constants. For non-magnetic systems the unpolarised theory suffices, but most modern functionals are written in their spin-polarised form because that is the more general expression.

### Finite temperature

Mermin (1965) extended the theorems to finite temperature: at temperature $T$, the equilibrium density determines the external potential, and a grand-canonical free-energy functional $\Omega[n]$ is minimised at the equilibrium density. This is the basis of "finite-$T$ DFT" used for warm dense matter.

### What about excited states?

HK theorems are emphatically ground-state theorems. The cleanest generalisation to excited states is **time-dependent DFT** (TD-DFT), based on the Runge–Gross theorem (1984): for fixed initial state $\Psi(0)$, the time-dependent density $n(\mathbf r,t)$ determines the time-dependent external potential $v_\mathrm{ext}(\mathbf r,t)$ up to a purely time-dependent constant. We discuss TD-DFT briefly in §5.6.

## 5.2.6 Summary

The Hohenberg–Kohn theorems are an existence proof. They tell us:

1. The ground-state density $n_0(\mathbf r)$ uniquely determines the external potential $v_\mathrm{ext}(\mathbf r)$ (up to a constant), and hence every property of the system.
2. There is a universal functional $F[n] = T[n] + V_{ee}[n]$, and the total energy functional $E_{v_\mathrm{ext}}[n] = F[n] + \int v_\mathrm{ext} n$ is minimised by the true ground-state density.

What we do *not* yet have:

- An explicit form for $F[n]$. Both $T[n]$ and $V_{ee}[n]$ are defined by constrained searches no more tractable than the original many-body problem.
- A practical algorithm.

Kohn and Sham, one year after Hohenberg and Kohn, supplied the second. They split $F[n]$ into a part we can compute exactly via auxiliary single-particle orbitals (the non-interacting kinetic energy $T_s$ and the classical Hartree energy), and a small remainder — the exchange–correlation energy $E_{xc}[n]$ — to be approximated. That is §5.3.
