# 5.1 The Thomas–Fermi Idea

!!! info "Why this section exists"
    Here is the problem we are solving. In 1927, before there were any computers, could you do electronic structure — work out the energy of atoms and molecules — using only the electron density $n(\mathbf r)$ instead of the full wavefunction $\Psi$? The density is a function of just three coordinates; the wavefunction is a function of $3N$. Thomas and Fermi made the first serious attempt at this. It only half-works: it gets the gross size and energy of heavy atoms roughly right, but it cannot bind a single molecule. Studying *why* it half-works is the sharpest possible way to see what a genuine density functional has to capture, and it is the natural on-ramp to the rest of this chapter.

!!! note "In plain language"
    The Thomas–Fermi idea, in one breath: approximate the kinetic energy of the electrons with a *local* formula that depends only on how dense the electrons are at each point — a formula borrowed from the uniform electron gas (electrons spread out evenly, like a featureless jelly) — add the obvious electrostatic energies, and then slide the density around until the total energy is as low as possible while keeping the right number of electrons. That "slide until lowest" step is a variational principle: the best density is the one that minimises the energy.

!!! note "New words in this section"
    - **Functional** — a rule that takes a *whole function* as its input and returns a single number. Ordinary functions eat a number and give a number; a functional eats a function (here the density $n(\mathbf r)$, defined at every point in space) and gives a number (here an energy). We write the input in square brackets, $E[n]$, to signal "functional of $n$".
    - **Kinetic energy functional** — the particular functional $T[n]$ that returns the electrons' kinetic energy given their density. Writing kinetic energy in terms of $n$ alone is the hard, central problem of this whole approach; Thomas–Fermi's guess for it, $T_\mathrm{TF}[n] = C_F\int n^{5/3}$, is what makes their theory tick — and what makes it fail.

    New to this vocabulary? See the [beginner glossary](../undergraduate/glossary-for-beginners.md).

!!! note "Why does this chapter exist?"
    Imagine you are a physicist in 1927. The Schrödinger equation has just been written down. You know in principle how to compute the properties of every atom and molecule in the universe — you just have to solve a partial differential equation in $3N$ variables, where $N$ is the number of electrons. For helium ($N=2$) that is six variables; for benzene ($N=42$) it is one hundred and twenty-six. Nobody has any hope of doing this on paper.
    
    Thomas and Fermi noticed something liberating. The electron *density* $n(\mathbf r)$ — how many electrons are at each point in space — lives in only three dimensions, no matter how many electrons there are. If you could write the total energy as a recipe involving only $n(\mathbf r)$, you would have escaped the curse of dimensionality entirely. This chapter is the story of their first attempt. It fails for molecules (you cannot bind hydrogen!) but it succeeds in showing that the *idea* is right. Forty years later Hohenberg, Kohn and Sham will turn this guess into a theorem and an algorithm. Everything in modern DFT begins here.

!!! abstract "Key idea (Chapter 5.1)"
    Thomas–Fermi theory replaces the $3N$-dimensional many-electron wavefunction $\Psi$ with the 3-dimensional electron density $n(\mathbf r)$. The total energy is approximated as a *local* functional of $n$: kinetic energy $C_F\int n^{5/3}$, classical Coulomb attraction $\int v_\mathrm{ext} n$, and electrostatic self-repulsion $\tfrac{1}{2}\iint n n'/|\mathbf r - \mathbf r'|$. This is the first density-functional theory, the prototype for everything that follows, and a beautiful illustration of why pure locality is not enough — Teller proved no molecule binds within TF.

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

!!! note "Why this step?"
    The local-density philosophy is an *educated guess*. It assumes that, on the length scale over which the electron density varies in a real atom, the gas looks locally uniform — much like saying that the local temperature inside a kettle obeys the bulk thermodynamics of water at that temperature. The guess is rigorously controllable only when $n$ varies slowly on the scale of the Fermi wavelength $\lambda_F \sim 1/k_F$. Inside an atom this is *not* the case (the density falls off exponentially on the scale of the Bohr radius, faster than $\lambda_F$), which already foreshadows trouble.

### A uniform non-interacting electron gas in a box

We will derive, step by step, the density of states of the homogeneous electron gas and from it the kinetic energy density as a function of $n$. Every algebraic line is included so that nothing in (5.4)–(5.5) has to be taken on faith.

Consider $N$ non-interacting electrons in a cubic box of side $L$ with periodic boundary conditions. The single-particle eigenstates are plane waves

$$
\phi_{\mathbf k}(\mathbf r) = \frac{1}{L^{3/2}}e^{i\mathbf k\cdot\mathbf r},
\qquad
\mathbf k = \frac{2\pi}{L}(n_x,n_y,n_z),\;\; n_i\in\mathbb Z,
$$

with single-particle energies $\varepsilon_{\mathbf k} = \hbar^{2}k^{2}/(2m)$. At zero temperature electrons fill all states up to the **Fermi wavevector** $k_F$, two per $\mathbf k$ for spin.

!!! note "Why this step?"
    The periodic boundary conditions $\phi_{\mathbf k}(\mathbf r+L\hat{\mathbf x}_i) = \phi_{\mathbf k}(\mathbf r)$ quantise the allowed momenta. Each $\mathbf k$-point occupies a cell of volume $(2\pi/L)^{3}$ in reciprocal space. In the thermodynamic limit $L\to\infty$ this cell shrinks and the discrete sum becomes a Riemann integral, which is the source of the standard substitution $\sum_{\mathbf k} \to (L^{3}/(2\pi)^{3})\int\mathrm d^{3}k$ used below.

The number of allowed $\mathbf k$-points inside a sphere of radius $k_F$ is the volume of the sphere divided by the volume of one $\mathbf k$-cell:

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

!!! note "Why this step?"
    The factor $3/5$ that appears implicitly through $\int_0^{k_F} k^{4}\,\mathrm dk = k_F^{5}/5$ is precisely the average kinetic energy of the filled Fermi sea expressed as a fraction of the Fermi energy $\varepsilon_F = \hbar^{2}k_F^{2}/(2m)$. One can rewrite the result more transparently as
    $$
    t(n) = \tfrac{3}{5}\,\varepsilon_F(n)\,n,
    $$
    i.e. *three-fifths of the Fermi energy times the density*. This is a useful sanity check: the per-particle kinetic energy of the Fermi sea is $\tfrac{3}{5}\varepsilon_F$.

!!! note "Atomic units"
    Throughout the chapter we work in Hartree atomic units: $\hbar = m_e = e = 4\pi\varepsilon_0 = 1$. Energies are in Hartrees ($1\;\text{Ha} = 27.21138\;\text{eV} = 4.35974\times 10^{-18}\;\text{J}$), lengths in Bohr radii ($1\;a_0 = 0.52918\;\text{Å} = 5.2918\times 10^{-11}\;\text{m}$). In these units $C_F = \tfrac{3}{10}(3\pi^{2})^{2/3}\approx 2.871$. Some references quote $C_F = \tfrac{3}{10}(3\pi^{2})^{2/3}\,\hbar^{2}/m_e$ in SI units, numerically $\approx 4.785\times 10^{-37}\;\text{J m}^{2}$ when $n$ is in $\text{m}^{-3}$. Always check the unit system before plugging into a formula.

!!! example "Numerical sanity check: copper at metallic density"
    Take metallic copper, $n \approx 8.47\times 10^{28}\;\text{m}^{-3}$, equivalently $n \approx 0.0126\;a_0^{-3}$ in atomic units. Then
    $$
    k_F = (3\pi^{2}\cdot 0.0126)^{1/3} \approx 0.722\;a_0^{-1}
    \;\;\Longrightarrow\;\; \varepsilon_F = \tfrac{1}{2}k_F^{2}\approx 0.261\;\text{Ha}\approx 7.10\;\text{eV},
    $$
    in line with the textbook value of $\sim 7.0\;\text{eV}$ for the Fermi energy of copper. The TF kinetic energy density is then $t(n) = C_F n^{5/3}\approx 2.871\times 0.0126^{5/3}\approx 2.0\times 10^{-3}\;\text{Ha}/a_0^{3}$, and the total kinetic energy per electron $t/n \approx 0.157\;\text{Ha}\approx 4.3\;\text{eV}$, which is $\tfrac{3}{5}\varepsilon_F$ as required.

### The local-density step

For a uniform gas the kinetic energy density is *exactly* $C_F n^{5/3}$. For an inhomogeneous system this is not true; the kinetic energy depends on how the density varies in space (it can have terms involving $\nabla n$). The Thomas–Fermi *ansatz* is to ignore those terms and integrate the uniform result over the inhomogeneous density:

$$
\boxed{\;\;T_\mathrm{TF}[n] \;=\; C_F\int n(\mathbf r)^{5/3}\,\mathrm d\mathbf r.\;\;}
\tag{5.5}
$$

This is the first **local density approximation** in the history of electronic structure theory. It is not exact; it is the leading term in a gradient expansion. We will meet that expansion again, with a vengeance, in §5.4.

!!! note "Symbol guide"
    The symbols introduced so far in this section, gathered in one place (atomic units, where $\hbar = m_e = e = 4\pi\varepsilon_0 = 1$):

    | Symbol | Meaning | Units (atomic) |
    | --- | --- | --- |
    | $n(\mathbf r)$ | electron density: electrons per unit volume at point $\mathbf r$ | $a_0^{-3}$ |
    | $N$ | total number of electrons, $N = \int n\,\mathrm d\mathbf r$ | dimensionless |
    | $\Psi$ | many-electron wavefunction (the thing we are trying to avoid) | $a_0^{-3N/2}$ |
    | $k_F$ | Fermi wavevector, $k_F = (3\pi^2 n)^{1/3}$ | $a_0^{-1}$ |
    | $\varepsilon_F$ | Fermi energy, $\varepsilon_F = \tfrac{1}{2}k_F^2$ | Ha |
    | $C_F$ | Thomas–Fermi constant, $\tfrac{3}{10}(3\pi^2)^{2/3}\approx 2.871$ | Ha$\cdot a_0^{2}$ |
    | $t(n)$ | kinetic energy *density* of the uniform gas, $C_F n^{5/3}$ | Ha$\cdot a_0^{-3}$ |
    | $T_\mathrm{TF}[n]$ | Thomas–Fermi kinetic energy *functional*, $C_F\int n^{5/3}\,\mathrm d\mathbf r$ | Ha |

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

!!! note "Why this step?"
    The derivative $\delta T_\mathrm{TF}/\delta n$ is straightforward but instructive. Writing $T_\mathrm{TF} = C_F\int f(n)\,\mathrm d\mathbf r$ with $f(n) = n^{5/3}$, the chain rule gives $\delta T_\mathrm{TF}/\delta n(\mathbf r) = C_F f'(n(\mathbf r)) = \tfrac{5}{3}C_F n^{2/3}$. There is no $\nabla n$-dependent piece because $T_\mathrm{TF}$ is purely local. This *locality* is what TF buys algebraic simplicity at the cost of physical fidelity. Compare with the Kohn–Sham kinetic functional in §5.3, where $T_s$ depends on the *orbitals*, hence on the *non-local* response of $n$ to the potential.

!!! example "The TF atomic ODE in scaled variables"
    For a spherical atom of nuclear charge $Z$, write $n(r) = (Z/4\pi a_0^{3})\,\chi(x)^{3/2}/x^{3/2}$ with $r = a_0 x$ and $a_0 = (3\pi/4)^{2/3}\cdot Z^{-1/3}/2$. The TF equation (5.9) then reduces to the dimensionless **Thomas–Fermi equation**
    $$
    \frac{\mathrm d^{2}\chi}{\mathrm d x^{2}} = \frac{\chi^{3/2}}{x^{1/2}},
    \qquad \chi(0) = 1, \;\;\chi(\infty) = 0\;\;\text{(neutral atom)}.
    $$
    This is one ODE for *every* atom — the variable $Z$ has been completely scaled out. Numerical solution gives $\chi'(0) \approx -1.588$, from which one extracts the total TF energy $E_\mathrm{TF}(Z) = -0.7687\,Z^{7/3}\;\text{Ha}$. Compare with helium ($Z=2$): TF predicts $-2.43\;\text{Ha}$; the exact non-relativistic energy is $-2.904\;\text{Ha}$. TF underbinds by about 16%, a remarkably good number for a theory with no orbital structure.

## 5.1.4a Dirac's exchange correction: Thomas–Fermi–Dirac

The original Thomas–Fermi functional contains no exchange. In 1930 Paul Dirac added one in the same spirit: compute the exchange energy of a uniform electron gas, then apply it locally. We rederive Dirac's constant here because we shall meet the same calculation again, dressed up as "LDA exchange", in §5.4.

For a single Slater determinant of plane waves filled up to $k_F$, the Hartree–Fock exchange energy per unit volume of the uniform gas is

$$
\epsilon_x^\mathrm{unif}(n) \cdot n \;=\; -\frac{1}{4\pi^{3}}\int_{k<k_F}\int_{k'<k_F}\frac{1}{|\mathbf k-\mathbf k'|^{2}}\,\mathrm d\mathbf k\,\mathrm d\mathbf k' \;\cdot\;\frac{1}{L^{0}}.
$$

The double integral over two Fermi spheres is a standard exercise; substituting $\mathbf q = \mathbf k - \mathbf k'$ and doing the angular integral gives, after algebra (we postpone the explicit steps to §5.4),

$$
\epsilon_x^\mathrm{unif}(n) = -\frac{3}{4\pi}\,k_F(n) = -\frac{3}{4}\Big(\frac{3}{\pi}\Big)^{1/3}\, n^{1/3}.
\tag{5.7a}
$$

!!! note "Why this step?"
    The combination $(3/\pi)^{1/3}$ is, up to a factor of order unity, just $k_F/n^{1/3}\sim (3\pi^{2})^{1/3}/\pi$. The proportionality $\epsilon_x \propto k_F \propto n^{1/3}$ is dictated by dimensions: the exchange energy density must have dimensions of energy per volume, and the only combination of $\hbar, m, e^{2}, n$ giving that is $\propto n^{4/3}$ in atomic units.

Applying the same local-density logic that gave (5.5), the **Thomas–Fermi–Dirac** exchange functional is

$$
\boxed{\;\;E_x^\mathrm{TFD}[n] \;=\; -\frac{3}{4}\Big(\frac{3}{\pi}\Big)^{1/3}\int n(\mathbf r)^{4/3}\,\mathrm d\mathbf r \;\;\equiv\;\; -C_x\int n^{4/3}\,\mathrm d\mathbf r,\;\;}
\tag{5.7b}
$$

with $C_x = \tfrac{3}{4}(3/\pi)^{1/3} \approx 0.7386$ in atomic units. The full Thomas–Fermi–Dirac energy functional is then

$$
E_\mathrm{TFD}[n] = C_F\int n^{5/3}\,\mathrm d\mathbf r \;-\; C_x\int n^{4/3}\,\mathrm d\mathbf r \;+\; \int v_\mathrm{ext} n\,\mathrm d\mathbf r \;+\; U_H[n].
\tag{5.7c}
$$

!!! example "Sanity check: exchange energy of helium"
    The electron density of helium has $\int n\,\mathrm d\mathbf r = 2$ and a peak near the nucleus of order $n\sim 3\;a_0^{-3}$. Plugging the converged Hartree–Fock density into (5.7b) gives $E_x^\mathrm{TFD}\approx -0.88\;\text{Ha}$ compared with the exact HF exchange of helium $E_x^\mathrm{HF}\approx -1.026\;\text{Ha}$. The LDA-style exchange underestimates the magnitude by roughly 14% — a typical error for LDA exchange in atoms, and one of the reasons §5.4 climbs further up Jacob's ladder.

The Dirac correction reduces the spurious self-interaction in the Hartree term — but only partially. For a hydrogen atom ($N=1$), the exact answer is that the Hartree-plus-exchange energy must vanish identically (the single electron cannot interact with itself). One can check by direct calculation that with the exact hydrogenic density $n_\mathrm{H}(\mathbf r) = |\psi_{1s}|^{2}$, the TFD combination $U_H[n_\mathrm{H}] + E_x^\mathrm{TFD}[n_\mathrm{H}]$ does *not* vanish: it leaves behind a residual self-interaction of about $0.08\;\text{Ha}\approx 2\;\text{eV}$. This residue is the **self-interaction error** that has haunted DFT ever since; we return to it in §5.4.7 and §5.6.4.

## 5.1.5 What Thomas–Fermi gets right

Thomas–Fermi is not a useless theory. It captures the correct overall trends of atomic binding energies: scaling arguments on (5.8) give the total energy of a neutral atom going as $-Z^{7/3}$ for large $Z$, which is in fact the leading term of the exact non-relativistic atomic energy. The total kinetic and Coulomb energies of atoms come out within ten or twenty per cent. As an order-of-magnitude theory for the gross structure of heavy atoms, it works.

!!! example "The $Z^{7/3}$ scaling law: where does it come from?"
    Minimise (5.8) for a neutral atom of nuclear charge $Z$ and ignore the Hartree term momentarily for orientation. The kinetic term scales as $\int n^{5/3}\,\mathrm d\mathbf r \sim n_\mathrm{typ}^{5/3}\,V_\mathrm{typ}$ where $V_\mathrm{typ}\sim a^{3}$ for a characteristic atomic radius $a$ and the electron count is $Z \sim n_\mathrm{typ}\,a^{3}$, so $n_\mathrm{typ}\sim Z/a^{3}$ and $T_\mathrm{TF}\sim Z^{5/3} a^{-2}$. The external potential energy scales as $V_\mathrm{ext}\sim -Z\cdot Z/a = -Z^{2}/a$. Minimising the sum over $a$ gives $a_\mathrm{opt}\sim Z^{-1/3}$ (atoms get smaller as $Z$ grows!) and $E_\mathrm{TF}\sim -Z^{7/3}$. For $Z=82$ (lead), this predicts $E\approx -(82)^{7/3}\approx -2\times 10^{4}\;\text{Ha}$, in the right ballpark; the exact non-relativistic value is $\approx -2.0\times 10^{4}\;\text{Ha}$. Not bad for the wrong theory.

It also has an important honourable mention: by minimising over a class of normalised densities (the *variational* formulation of Levy and Lieb, which we revisit in §5.2), one can make Thomas–Fermi theory mathematically rigorous as a *lower bound* on the true ground-state energy when the kinetic functional is replaced by the Lieb–Thirring inequality. This is a foundational result in mathematical physics; for our purposes, what matters is that the *variational structure* of (5.8) is the right idea. Hohenberg and Kohn will keep it; they will only replace the approximate functional with an exact one.

## 5.1.6 Why Thomas–Fermi fails for chemistry

For chemistry — which is to say, for everything we care about in materials science: bond lengths, lattice constants, surface energies, reaction barriers — Thomas–Fermi is hopeless.

!!! warning "Common misunderstanding"
    It is tempting to conclude from this section that "density-based methods cannot work". That is the wrong lesson. Thomas–Fermi fails not because the density is an inadequate variable, but because one *specific* ingredient — a purely **local** kinetic-energy functional $C_F\int n^{5/3}$ — is too crude. It cannot see shell structure, and (Teller's theorem) it binds no molecule. The density itself is fine: Hohenberg and Kohn later prove it carries *all* the information. The fix, supplied by Kohn and Sham (§5.3), is not to abandon the density but to compute the kinetic energy from orbitals instead of from a local formula. So: blame the local kinetic functional, not the density.

**No shell structure.** The 5/3 power in the kinetic functional treats the density in a hydrogen atom as if it were a slab of uniform gas of the same local density. Real atoms have shells, sharp radial features arising from the orthogonality of $\phi_{1s}, \phi_{2s}, \phi_{2p},\dots$ — quantum interference between orbitals. The local-uniform kinetic functional cannot represent this; the predicted radial density for, say, argon is a smooth monotone decay, with no hint of the K, L, M shells.

**No covalent bonding: Teller's no-binding theorem.** Edward Teller proved in 1962 a striking and discouraging result: *within Thomas–Fermi theory, no molecule is stable*. That is, for any arrangement of nuclei, the Thomas–Fermi energy of the molecule as a function of internuclear separation has its minimum at infinite separation — atoms always prefer to dissociate. We state the theorem carefully and sketch the structure of the proof; for the full argument see Lieb and Simon's mathematically rigorous treatment (Adv. Math. **23**, 22 (1977)).

!!! note "Teller's theorem (statement)"
    Let $E_\mathrm{TF}(\{\mathbf R_\alpha\}, \{Z_\alpha\})$ denote the Thomas–Fermi ground-state energy of a system of nuclei of charges $Z_\alpha$ at positions $\mathbf R_\alpha$ with a total of $N \leq \sum_\alpha Z_\alpha$ electrons. Then
    $$
    E_\mathrm{TF}(\{\mathbf R_\alpha\}, \{Z_\alpha\}) \;\geq\; \sum_\alpha E_\mathrm{TF}^\mathrm{atom}(Z_\alpha),
    $$
    with strict inequality unless all nuclei coincide. *No molecular Thomas–Fermi system has a binding energy minimum at finite separation.*

The proof relies on three key ingredients:

1. *Scaling.* The TF functional (5.8) has a specific scaling under uniform rescaling of $n$ and lengths. Defining $n_\lambda(\mathbf r) = \lambda^{3}n(\lambda \mathbf r)$ one finds $T_\mathrm{TF}[n_\lambda] = \lambda^{2}T_\mathrm{TF}[n]$, $V_\mathrm{ext}[n_\lambda] = \lambda V_\mathrm{ext}[n]$, $U_H[n_\lambda] = \lambda U_H[n]$. The virial theorem $2T_\mathrm{TF} + V_\mathrm{ext} + 2U_H = 0$ at the minimum follows, fixing the equilibrium scale.
2. *Concavity of the energy in nuclear positions.* A clever rearrangement inequality (used by Lieb, Simon, and Thirring) shows that the TF energy is *concave* in $-Z_\alpha/|\mathbf r - \mathbf R_\alpha|$, which means the energy is at least as large for a "spread out" arrangement as for the coincident one.
3. *No barrier.* Combining 1 and 2 with the obvious limit $E_\mathrm{TF}\to\sum_\alpha E_\mathrm{TF}^\mathrm{atom}$ at infinite separation shows that pulling nuclei apart never raises the energy above the dissociation limit.

!!! note "Why this step?"
    The intuition behind Teller's theorem is that the Thomas–Fermi kinetic functional $\int n^{5/3}$ is *convex* in $n$. When two atomic densities overlap, the combined density at the overlap point is *larger* than either one alone, so the kinetic cost grows superlinearly. The Coulomb attraction gained by bringing the nuclei together cannot compensate for this kinetic price within the local approximation. Real quantum mechanics evades this trap because the *exact* kinetic energy depends on the gradients of $n$, not just its local value: a smoother density in the bond region lowers $T$ even as it raises $\int n^{5/3}$.

This is fatal. Materials science is the science of bonds — covalent, metallic, ionic, hydrogen, van der Waals. A theory that cannot predict the existence of $\mathrm H_2$ is not a theory of materials. Adding Dirac's exchange (5.7b) does not save the day: TFD also fails to bind molecules, though the proof is more involved. One needs a fundamentally non-local treatment of the kinetic energy — which is exactly what Kohn and Sham will provide in §5.3 by re-introducing orbitals.

**The kinetic functional is too soft.** The deeper reason behind Teller's theorem is that $C_F\int n^{5/3}$ is a poor approximation to the true kinetic energy for rapidly varying densities — and the inter-nuclear region, with its build-up of charge characteristic of a bond, is exactly where the density varies rapidly. Gradient corrections to (5.5) help (von Weizsäcker added a $\frac{1}{8}|\nabla n|^{2}/n$ term that improves matters near nuclei), but the resulting **orbital-free DFT** has remained, decades later, a research field rather than a workhorse.

!!! warning "Why this matters"
    Modern *orbital-free DFT* is essentially Thomas–Fermi with better kinetic functionals, and it is an active research area: it has the unbeatable property of scaling linearly with system size and being trivially parallelisable. But for general systems no kinetic functional accurate enough to compete with Kohn–Sham has yet been found. The kinetic energy is, in a precise sense, the hardest part of the energy to write as an explicit functional of $n$.

**No negative ions.** A subtle but striking failure: within TF and TFD, no atom can bind more electrons than its nuclear charge. The chemical potential $\mu$ in (5.9) vanishes at neutrality, and adding one more electron forces $\mu>0$, which has no normalisable solution. So Cl$^{-}$, F$^{-}$, H$^{-}$ — all of chemistry's anions — are unbound in Thomas–Fermi. The proof is again a scaling argument: a system can only support the extra electron through correlation effects entirely absent in TF. Lieb proved in 1981 that the maximum number of electrons bound by a nucleus of charge $Z$ in *any* TF-type theory is at most $Z + 1$ (in fact $Z$ exactly for plain TF), whereas real atoms bind H$^{-}$ ($Z=1, N=2$), O$^{2-}$ ($Z=8, N=10$), and so on. This is another piece of chemistry the local-density picture cannot reach.

**The smooth density problem.** Inside any real atom, the radial density has visible kinks where shells start and stop: at the boundary between K and L shells in argon, for instance, $\mathrm d n/\mathrm d r$ changes abruptly. The TF density, governed by the smooth ODE $\chi''=\chi^{3/2}/x^{1/2}$, is analytic everywhere except at $r=0$. It cannot resolve shell boundaries by construction. This is the same defect that, in a different guise, causes TF to miss the periodicity of the periodic table: ionisation potentials computed from TF vary monotonically with $Z$, not in the famous zigzag pattern of experiment.

### Beyond TF: the von Weizsäcker correction

Carl Friedrich von Weizsäcker (1935) noted that the TF functional misses the kinetic energy associated with *spatial variation* of $n$. By expanding the exact kinetic energy in gradients of $n$, he proposed the correction

$$
T_W[n] \;=\; \frac{1}{8}\int\frac{|\nabla n(\mathbf r)|^{2}}{n(\mathbf r)}\,\mathrm d\mathbf r,
\tag{5.9a}
$$

which is *exact* for a single orbital (one-electron system) and gives the leading gradient correction to the uniform-gas result. Modern orbital-free DFT typically uses a combination

$$
T_\mathrm{OF}[n] \;=\; T_\mathrm{TF}[n] \;+\; \lambda\,T_W[n] \;+\; \cdots,
$$

with $\lambda$ between 1/9 (the value from the formal gradient expansion at small $s$) and 1 (the von Weizsäcker upper bound). The full gradient expansion is asymptotic, not convergent — one cannot simply keep adding higher-order terms. Even with the best modern kinetic-energy functionals, orbital-free DFT only matches Kohn–Sham accuracy for nearly-free-electron systems (alkali metals, aluminium); for transition metals, semiconductors, and molecules it remains stubbornly off.

!!! note "Why this step?"
    The von Weizsäcker term has a transparent physical meaning. For a single orbital $\phi$ with density $n = |\phi|^{2}$ and a real wavefunction, $T_W = \tfrac{1}{2}\int|\nabla\phi|^{2}$ exactly — it *is* the kinetic energy. For an inhomogeneous many-electron system, $T_W$ is a *lower bound* on the true kinetic energy (a result due to Sears, Parr, and Dinur). The TF term, by contrast, undershoots in regions of rapid variation and overshoots in regions of slow variation; combining the two with $\lambda \sim 1/9$ partially heals both errors.

## 5.1.7 What we have learnt, and what comes next

Thomas and Fermi taught us three things:

1. *In principle*, the density alone can be enough — at least for total energies — and minimising an energy functional gives a well-defined variational scheme.
2. The hard part is the kinetic energy. Approximating it as a *local* functional of the density misses physics — shell structure, bonding — that no chemistry can do without.
3. Without exchange and correlation, the electrostatic Hartree energy is not enough either.

Thirty-five years later, in 1964, Hohenberg and Kohn turned (1) into a theorem: there *is* an exact density functional. In 1965 Kohn and Sham resolved (2) by re-introducing single-particle orbitals not to represent the wavefunction (we have given that up) but to *compute the kinetic energy*, leaving only a small unknown — the exchange–correlation energy — to be approximated. The Thomas–Fermi failure of (3) became the problem of choosing an exchange–correlation functional, which is what §5.4 is about.

A practical density functional theory exists. The next sections build it.

### Summary of §5.1

- *Density-functional ansatz*: replace $\Psi(\mathbf r_1,\dots,\mathbf r_N)$ with $n(\mathbf r)$ as the basic variable. Trade: $3N$ dimensions for $3$.
- *Thomas–Fermi kinetic functional*: $T_\mathrm{TF}[n] = C_F\int n^{5/3}\,\mathrm d\mathbf r$ with $C_F = \tfrac{3}{10}(3\pi^{2})^{2/3}\approx 2.871\;\text{Ha}\cdot a_0^{2}$, derived from filling free-electron states to the Fermi sea.
- *Dirac exchange*: $E_x^\mathrm{TFD}[n] = -C_x\int n^{4/3}\,\mathrm d\mathbf r$ with $C_x \approx 0.7386$, the first local exchange functional.
- *Successes*: correct $-Z^{7/3}$ scaling of atomic energies; reasonable kinetic and Coulomb energies for heavy atoms.
- *Fatal failures*: no shell structure (the density is monotonic), no molecular binding (Teller's theorem), no negative ions.
- *Historical context*: Thomas (1927), Fermi (1927), Dirac (1930), Teller's no-binding theorem (1962), Lieb–Simon rigorous treatment (1977).

!!! note "Remark: historical attribution"
    Thomas published his analysis in *Proc. Camb. Phil. Soc.* **23**, 542 (1927), independently of Fermi who published in *Rend. Accad. Naz. Lincei* **6**, 602 (1927). Both papers proposed essentially the same functional, and the field has used "Thomas–Fermi" symmetrically ever since. Dirac's contribution (1930) added the exchange term and originated the name *exchange energy* in the DFT context, though the exchange operator itself goes back to Fock's 1930 paper on Hartree–Fock theory.

!!! question "Check yourself"
    1. What basic variable does Thomas–Fermi theory use *instead of* the many-electron wavefunction $\Psi$, and why is that such a large saving?
    2. Why is Thomas–Fermi theory called *variational*? What quantity is minimised, and under what constraint?
    3. The kinetic energy functional $T_\mathrm{TF}[n] = C_F\int n^{5/3}$ is "local". What does *local* mean here, and which physical feature does that locality fail to capture inside a real atom?
    4. Name one concrete thing Thomas–Fermi gets wrong, and say in one sentence why.

??? success "Answers"
    1. It uses the **electron density** $n(\mathbf r)$, a function of just $3$ coordinates, in place of $\Psi(\mathbf r_1,\dots,\mathbf r_N)$, a function of $3N$ coordinates. For $N$ electrons the wavefunction needs of order (grid points)$^{N}$ numbers while the density needs only (grid points) numbers — the curse of dimensionality is removed.
    2. Because the working principle is to **minimise the total energy functional** $E_\mathrm{TF}[n]$ over all candidate densities, subject to the constraint that they integrate to the correct electron number, $\int n\,\mathrm d\mathbf r = N$. The ground-state density is the minimiser; the constraint is enforced by the Lagrange multiplier $\mu$ in equation (5.9).
    3. *Local* means the kinetic energy density at a point $\mathbf r$ is fixed by the value of $n$ *at that same point alone* — it uses the uniform-gas formula $C_F n(\mathbf r)^{5/3}$ and ignores how $n$ varies nearby (no $\nabla n$). Inside a real atom this misses **shell structure** (K, L, M shells), which comes from orbital orthogonality and shows up only through the spatial variation of the density.
    4. Any one of: **no molecular binding** (Teller's theorem — the energy is lowest with the atoms infinitely apart, because the local kinetic functional overcharges for the density build-up in a bond); **no shell structure** (the predicted radial density is a smooth monotone decay); or **no negative ions** (TF cannot bind more electrons than the nuclear charge, since correlation, absent in TF, is what stabilises anions).
