# 5.1 The Thomas–Fermi Idea

## 5.1.1 A radical economy of variables

The year is 1927. Schrödinger's wave mechanics is two years old. Hartree has just proposed his self-consistent field method, in which each electron moves in the average potential of all the others. Llewellyn Thomas, a 24-year-old at Cambridge, and independently Enrico Fermi in Rome, ask a more daring question. The many-electron wavefunction $\Psi(\mathbf r_1,\dots,\mathbf r_N)$ is a function of $3N$ spatial coordinates. The electron density

$$
n(\mathbf r) = N\int |\Psi(\mathbf r,\mathbf r_2,\dots,\mathbf r_N)|^{2}\,\mathrm d\mathbf r_2\cdots \mathrm d\mathbf r_N
\tag{5.1}
$$

is a function of just three. *Can we work directly with $n(\mathbf r)$ and never write down $\Psi$?*

The trade is dazzling. For $N=100$ electrons, the wavefunction on a $10^3$ grid would require $10^{300}$ floating-point numbers. The density on the same grid requires $10^3$. If we could express the total energy as a functional of $n$ alone,

$$
E[n] = T[n] + V_\mathrm{ext}[n] + V_{ee}[n],
\tag{5.2}
$$

and minimise it subject to $\int n\,\mathrm d\mathbf r = N$, we would have an electronic structure theory of stupendous economy. This is the founding fantasy of DFT. Thomas and Fermi gave the first concrete attempt to realise it, and it is well worth understanding both for what it gets right and for the very specific way in which it fails.

## 5.1.2 The kinetic energy of a uniform electron gas

The hard part of (5.2) is $T[n]$. We know the kinetic energy operator acts on the wavefunction:

$$
T = -\frac{\hbar^{2}}{2m}\sum_i \int \Psi^{*}\nabla_i^{2}\Psi\,\mathrm d\mathbf r_1\cdots\mathrm d\mathbf r_N,
$$

so writing $T$ as a functional of $n$ alone is not obviously possible. Thomas and Fermi made a *local* approximation. Imagine the electron gas were uniform with density $n$ everywhere; compute its kinetic energy per unit volume exactly; then pretend that, even in a real inhomogeneous system, the kinetic energy density at the point $\mathbf r$ is the same as that of a uniform gas with density $n(\mathbf r)$. We derive the uniform result.

### A uniform non-interacting electron gas in a box

Consider $N$ non-interacting electrons in a cubic box of side $L$ with periodic boundary conditions. The single-particle eigenstates are plane waves

$$
\phi_{\mathbf k}(\mathbf r) = \frac{1}{L^{3/2}}e^{i\mathbf k\cdot\mathbf r},
\qquad
\mathbf k = \frac{2\pi}{L}(n_x,n_y,n_z),\;\; n_i\in\mathbb Z,
$$

with single-particle energies $\varepsilon_{\mathbf k} = \hbar^{2}k^{2}/(2m)$. At zero temperature electrons fill all states up to the **Fermi wavevector** $k_F$, two per $\mathbf k$ for spin. The number of allowed $\mathbf k$-points inside a sphere of radius $k_F$ is the volume of the sphere divided by the volume of one $\mathbf k$-cell:

$$
\#\mathbf k = \frac{\tfrac{4}{3}\pi k_F^{3}}{(2\pi/L)^{3}} = \frac{L^{3} k_F^{3}}{6\pi^{2}}.
$$

Doubling for spin, the total number of electrons is $N = L^{3}k_F^{3}/(3\pi^{2})$, so the density and Fermi wavevector are linked by

$$
n = \frac{N}{L^{3}} = \frac{k_F^{3}}{3\pi^{2}}
\;\;\Longleftrightarrow\;\;
k_F(n) = (3\pi^{2}n)^{1/3}.
\tag{5.3}
$$

The total kinetic energy is the sum of $\varepsilon_{\mathbf k}$ over occupied states. Converting the sum to an integral, $\sum_{\mathbf k} \to \frac{L^3}{(2\pi)^3}\int\mathrm d^{3}k$, and including the factor of two for spin,

$$
T = 2\cdot\frac{L^{3}}{(2\pi)^{3}}\int_{|\mathbf k|<k_F}\frac{\hbar^{2}k^{2}}{2m}\,\mathrm d^{3}k
= \frac{L^{3}}{(2\pi)^{3}}\cdot\frac{\hbar^{2}}{m}\cdot 4\pi\int_0^{k_F}k^{4}\,\mathrm dk
= \frac{L^{3}\hbar^{2}}{10\pi^{2}m}k_F^{5}.
$$

The kinetic energy per unit volume is therefore

$$
t(n) \equiv \frac{T}{L^{3}} = \frac{\hbar^{2}}{10\pi^{2}m}\,k_F^{5}
= \frac{\hbar^{2}}{10\pi^{2}m}\big(3\pi^{2}n\big)^{5/3}
= \frac{3\hbar^{2}}{10m}(3\pi^{2})^{2/3}\,n^{5/3}.
$$

Defining the **Thomas–Fermi constant**

$$
C_F \equiv \frac{3\hbar^{2}}{10m}(3\pi^{2})^{2/3} \approx 2.871\;\text{(Hartree atomic units)},
\tag{5.4}
$$

we obtain the kinetic energy density of the uniform electron gas as a function of its density:

$$
t(n) = C_F\, n^{5/3}.
$$

!!! note "Atomic units"
    Throughout the chapter we work in Hartree atomic units: $\hbar = m_e = e = 4\pi\varepsilon_0 = 1$. Energies are in Hartrees ($\approx 27.211\;\text{eV}$), lengths in Bohr radii ($\approx 0.529\;\text{Å}$). In these units $C_F = \tfrac{3}{10}(3\pi^{2})^{2/3}\approx 2.871$.

### The local-density step

For a uniform gas the kinetic energy density is *exactly* $C_F n^{5/3}$. For an inhomogeneous system this is not true; the kinetic energy depends on how the density varies in space (it can have terms involving $\nabla n$). The Thomas–Fermi *ansatz* is to ignore those terms and integrate the uniform result over the inhomogeneous density:

$$
\boxed{\;\;T_\mathrm{TF}[n] \;=\; C_F\int n(\mathbf r)^{5/3}\,\mathrm d\mathbf r.\;\;}
\tag{5.5}
$$

This is the first **local density approximation** in the history of electronic structure theory. It is not exact; it is the leading term in a gradient expansion. We will meet that expansion again, with a vengeance, in §5.4.

## 5.1.3 The other pieces

### External potential energy

The interaction of the density with the external potential — for atoms and molecules, the nuclei — is just

$$
V_\mathrm{ext}[n] = \int v_\mathrm{ext}(\mathbf r)\,n(\mathbf r)\,\mathrm d\mathbf r.
\tag{5.6}
$$

There is nothing approximate about this expression: it follows from the diagonal-in-$\mathbf r$ nature of the external potential. For nuclei of charge $Z_\alpha$ at positions $\mathbf R_\alpha$,

$$
v_\mathrm{ext}(\mathbf r) = -\sum_\alpha \frac{Z_\alpha}{|\mathbf r - \mathbf R_\alpha|}.
$$

### Electron–electron energy: the Hartree term

The full electron–electron operator is a pair operator and depends on the *two*-body density. Thomas and Fermi replaced it with its mean-field, classical electrostatic part — the **Hartree energy**:

$$
U_H[n] \;=\; \frac{1}{2}\iint \frac{n(\mathbf r)\,n(\mathbf r')}{|\mathbf r - \mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r'.
\tag{5.7}
$$

This is the electrostatic self-energy of the charge cloud $n(\mathbf r)$. The factor of one half avoids double-counting the pair $(i,j)$ and $(j,i)$. It misses (i) the **exchange** energy demanded by antisymmetry of the fermionic wavefunction, and (ii) **correlation**, which encodes the avoidance dance electrons perform beyond the Pauli principle. Both are absent in the original Thomas–Fermi functional; Dirac later added an LDA-style exchange term, giving **Thomas–Fermi–Dirac** theory. We meet Dirac's exchange in §5.4.

It also includes an unphysical *self-interaction*: the field $\mathbf r'$ in the integral runs over the entire density, including the bit corresponding to electron at $\mathbf r$ itself. We return to this self-interaction error in §5.6.

## 5.1.4 The Thomas–Fermi energy functional

Putting the pieces together,

$$
E_\mathrm{TF}[n] = C_F\int n^{5/3}\,\mathrm d\mathbf r
\;+\; \int v_\mathrm{ext}(\mathbf r)\,n(\mathbf r)\,\mathrm d\mathbf r
\;+\; \frac{1}{2}\iint \frac{n(\mathbf r)\,n(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r'.
\tag{5.8}
$$

We minimise (5.8) over densities with the constraint $\int n\,\mathrm d\mathbf r = N$. Introducing a Lagrange multiplier $\mu$ for the particle-number constraint, the stationarity condition $\delta(E_\mathrm{TF} - \mu N)/\delta n(\mathbf r) = 0$ gives, term by term:

- $\delta T_\mathrm{TF}/\delta n = \tfrac{5}{3}C_F n^{2/3}$,
- $\delta V_\mathrm{ext}/\delta n = v_\mathrm{ext}(\mathbf r)$,
- $\delta U_H/\delta n = \int n(\mathbf r')/|\mathbf r-\mathbf r'|\,\mathrm d\mathbf r' \equiv v_H(\mathbf r)$,

so the **Thomas–Fermi equation** reads

$$
\boxed{\;\;\tfrac{5}{3}C_F\, n(\mathbf r)^{2/3} \;+\; v_\mathrm{ext}(\mathbf r) \;+\; v_H(\mathbf r) \;=\; \mu.\;\;}
\tag{5.9}
$$

The Lagrange multiplier $\mu$ plays the role of a chemical potential. Equation (5.9) is a single nonlinear equation in the three-dimensional density $n(\mathbf r)$. For a spherical atom it reduces to a one-dimensional ODE (the celebrated Thomas–Fermi equation $\phi''=\phi^{3/2}/x^{1/2}$ in scaled variables), which can be solved numerically with ease.

## 5.1.5 What Thomas–Fermi gets right

Thomas–Fermi is not a useless theory. It captures the correct overall trends of atomic binding energies: scaling arguments on (5.8) give the total energy of a neutral atom going as $-Z^{7/3}$ for large $Z$, which is in fact the leading term of the exact non-relativistic atomic energy. The total kinetic and Coulomb energies of atoms come out within ten or twenty per cent. As an order-of-magnitude theory for the gross structure of heavy atoms, it works.

It also has an important honourable mention: by minimising over a class of normalised densities (the *variational* formulation of Levy and Lieb, which we revisit in §5.2), one can make Thomas–Fermi theory mathematically rigorous as a *lower bound* on the true ground-state energy when the kinetic functional is replaced by the Lieb–Thirring inequality. This is a foundational result in mathematical physics; for our purposes, what matters is that the *variational structure* of (5.8) is the right idea. Hohenberg and Kohn will keep it; they will only replace the approximate functional with an exact one.

## 5.1.6 Why Thomas–Fermi fails for chemistry

For chemistry — which is to say, for everything we care about in materials science: bond lengths, lattice constants, surface energies, reaction barriers — Thomas–Fermi is hopeless.

**No shell structure.** The 5/3 power in the kinetic functional treats the density in a hydrogen atom as if it were a slab of uniform gas of the same local density. Real atoms have shells, sharp radial features arising from the orthogonality of $\phi_{1s}, \phi_{2s}, \phi_{2p},\dots$ — quantum interference between orbitals. The local-uniform kinetic functional cannot represent this; the predicted radial density for, say, argon is a smooth monotone decay, with no hint of the K, L, M shells.

**No covalent bonding: Teller's theorem.** Edward Teller proved in 1962 a striking result: *within Thomas–Fermi theory, no molecule is stable*. That is, for any arrangement of nuclei, the energy of the molecule as a function of internuclear separation has its minimum at infinite separation — atoms always prefer to dissociate. The proof uses scaling arguments on (5.8) and is not difficult; we sketch it.

Consider two atoms with nuclear positions $\mathbf R_A,\mathbf R_B$ and combined density $n(\mathbf r;\mathbf R_A,\mathbf R_B)$ minimising $E_\mathrm{TF}$. Take the limit $|\mathbf R_A-\mathbf R_B|\to\infty$: the optimal density factorises into two separated atomic densities and $E_\mathrm{TF} \to E_\mathrm{TF}[n_A] + E_\mathrm{TF}[n_B] = E_A + E_B$. Teller shows by a clever rearrangement inequality on the kinetic and Coulomb terms that the molecular energy is *never below* this separated-atom limit. Hence no chemical bond.

This is fatal. Materials science is the science of bonds — covalent, metallic, ionic, hydrogen, van der Waals. A theory that cannot predict the existence of $\mathrm H_2$ is not a theory of materials.

**The kinetic functional is too soft.** The deeper reason behind Teller's theorem is that $C_F\int n^{5/3}$ is a poor approximation to the true kinetic energy for rapidly varying densities — and the inter-nuclear region, with its build-up of charge characteristic of a bond, is exactly where the density varies rapidly. Gradient corrections to (5.5) help (von Weizsäcker added a $\frac{1}{8}|\nabla n|^{2}/n$ term that improves matters near nuclei), but the resulting **orbital-free DFT** has remained, decades later, a research field rather than a workhorse.

!!! warning "Why this matters"
    Modern *orbital-free DFT* is essentially Thomas–Fermi with better kinetic functionals, and it is an active research area: it has the unbeatable property of scaling linearly with system size and being trivially parallelisable. But for general systems no kinetic functional accurate enough to compete with Kohn–Sham has yet been found. The kinetic energy is, in a precise sense, the hardest part of the energy to write as an explicit functional of $n$.

## 5.1.7 What we have learnt, and what comes next

Thomas and Fermi taught us three things:

1. *In principle*, the density alone can be enough — at least for total energies — and minimising an energy functional gives a well-defined variational scheme.
2. The hard part is the kinetic energy. Approximating it as a *local* functional of the density misses physics — shell structure, bonding — that no chemistry can do without.
3. Without exchange and correlation, the electrostatic Hartree energy is not enough either.

Thirty-five years later, in 1964, Hohenberg and Kohn turned (1) into a theorem: there *is* an exact density functional. In 1965 Kohn and Sham resolved (2) by re-introducing single-particle orbitals not to represent the wavefunction (we have given that up) but to *compute the kinetic energy*, leaving only a small unknown — the exchange–correlation energy — to be approximated. The Thomas–Fermi failure of (3) became the problem of choosing an exchange–correlation functional, which is what §5.4 is about.

A practical density functional theory exists. The next sections build it.
