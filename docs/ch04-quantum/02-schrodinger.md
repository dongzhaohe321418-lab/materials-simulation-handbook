# 4.2 The Schrödinger equation

We now state the central equation of non-relativistic quantum mechanics. Everything that follows in this book — the band structure of silicon, the binding energy of a benzene molecule, the vibrational spectrum of a zeolite — is a consequence, or an approximation to a consequence, of this single postulate. We will *not* derive it. Schrödinger himself arrived at it by guesswork guided by analogy with classical wave optics, and there is no logical sense in which it can be deduced from earlier physics. It is a postulate, justified solely by the overwhelming agreement of its predictions with experiment.

## 4.2.1 The time-dependent Schrödinger equation

For a single non-relativistic particle of mass $m$ moving in a potential $V(\mathbf r, t)$, the state is described by a complex-valued wavefunction $\psi(\mathbf r, t)$, and the wavefunction evolves according to the **time-dependent Schrödinger equation** (TDSE):

$$\boxed{\; i\hbar \frac{\partial \psi(\mathbf r, t)}{\partial t} = \hat H\, \psi(\mathbf r, t) \;} \tag{4.2.1}$$

where $\hat H$ is the *Hamiltonian operator*,

$$\hat H = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf r, t). \tag{4.2.2}$$

Several features deserve immediate comment.

- **It is first order in time.** Unlike Newton's equations (second order) or the classical wave equation (also second order), the TDSE needs only an initial condition $\psi(\mathbf r, 0)$ to determine all future evolution. There is no need for an initial "velocity of the wavefunction".

- **It is complex.** The factor of $i$ on the left is essential. A complex-valued wavefunction is not a calculational convenience that can be replaced by a real two-component vector at the end of the day; the phase relationships between different parts of $\psi$ encode interference effects that have been observed in countless experiments.

- **It is linear.** If $\psi_1$ and $\psi_2$ both satisfy (4.2.1), so does $\alpha\psi_1 + \beta\psi_2$ for any complex constants $\alpha, \beta$. This is the principle of superposition, and it is the formal source of every characteristically quantum phenomenon.

- **The Hamiltonian is the energy operator.** Notice that $-\frac{\hbar^2}{2m}\nabla^2$ is what you get if you take the classical kinetic energy $p^2/2m$ and substitute $\mathbf p \to -i\hbar \nabla$. This substitution — momentum becomes a differential operator — is the heuristic device by which Schrödinger guessed his equation, but again, it is no derivation.

!!! note "A plausibility argument, not a derivation"
    For a free particle ($V = 0$) one can check that the plane wave
    $\psi(\mathbf r, t) = \exp[i(\mathbf k\cdot\mathbf r - \omega t)]$
    solves (4.2.1) provided $\hbar\omega = \hbar^2 k^2/2m$, i.e.\ $E = p^2/2m$ with de Broglie's $p = \hbar k$ and the Planck relation $E = \hbar\omega$. This is the standard motivational pattern: the Schrödinger equation is the simplest linear PDE whose plane-wave solutions reproduce the non-relativistic dispersion relation that de Broglie and Einstein together demand. Schrödinger generalised the kinetic-energy term to include $V(\mathbf r)$ by analogy with the classical Hamiltonian.

## 4.2.2 The Born rule

A complex wavefunction is not, by itself, an observable. We have to specify how $\psi$ relates to measurement outcomes. Max Born's 1926 proposal, for which he received the 1954 Nobel Prize, is now universal:

$$\rho(\mathbf r, t) \equiv |\psi(\mathbf r, t)|^2 = \psi^*(\mathbf r, t)\, \psi(\mathbf r, t) \tag{4.2.3}$$

is the **probability density** of finding the particle at position $\mathbf r$ at time $t$. The probability of finding it in a small volume $d^3 r$ around $\mathbf r$ is $\rho(\mathbf r, t)\, d^3 r$.

This is a probabilistic, not deterministic, theory. The Schrödinger equation evolves $\psi$ deterministically, but $\psi$ only tells you the *odds* of various measurement outcomes. A single experiment yields a single result; the predictions of quantum mechanics are statistical and only become sharp upon averaging over many identically-prepared systems.

Because $\rho$ is a probability density, it must integrate to one:

$$\int |\psi(\mathbf r, t)|^2\, d^3 r = 1. \tag{4.2.4}$$

This is the **normalisation condition**. A wavefunction satisfying (4.2.4) is called *normalised*. Notice that any solution of (4.2.1) can be rescaled by a constant without breaking the equation, so we always choose the constant so that (4.2.4) holds.

!!! warning "Conservation of probability"
    For (4.2.4) to remain valid at all times, the *total* probability must be conserved. Differentiating (4.2.4) and using (4.2.1) one finds
    $$\frac{\partial \rho}{\partial t} + \nabla\cdot \mathbf j = 0, \qquad \mathbf j = \frac{\hbar}{2mi}(\psi^* \nabla\psi - \psi \nabla\psi^*).$$
    This continuity equation has the same form as the conservation of electric charge or fluid mass: $\rho$ is a density, $\mathbf j$ is the **probability current**. We will use $\mathbf j$ implicitly when discussing transport in Chapter 10.

## 4.2.3 Stationary states and the time-independent Schrödinger equation

A vast amount of practical quantum mechanics — almost everything we do in materials physics — boils down to looking for *stationary states*: solutions whose probability density does not change with time. We try a separable ansatz,

$$\psi(\mathbf r, t) = \phi(\mathbf r)\, f(t),$$

and substitute into (4.2.1) with a time-independent potential $V(\mathbf r)$. The TDSE becomes

$$i\hbar\, \phi(\mathbf r)\, \dot f(t) = f(t)\, \hat H \phi(\mathbf r),$$

and dividing by $\phi f$,

$$i\hbar\, \frac{\dot f(t)}{f(t)} = \frac{\hat H \phi(\mathbf r)}{\phi(\mathbf r)}.$$

The left side depends only on $t$, the right only on $\mathbf r$; both must therefore equal a common constant, which we call $E$. This gives two equations. The time part is solved immediately,

$$f(t) = e^{-iEt/\hbar}, \tag{4.2.5}$$

and the spatial part is the **time-independent Schrödinger equation** (TISE):

$$\boxed{\; \hat H\, \phi(\mathbf r) = E\, \phi(\mathbf r) \;} \tag{4.2.6}$$

with $\hat H = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf r)$.

Equation (4.2.6) is an **eigenvalue equation**: we seek functions $\phi$ that are mapped, by the Hamiltonian, into multiples of themselves. The multipliers $E$ are the allowed *energies* of the system. For typical bound-state problems they form a discrete set $\{E_0, E_1, E_2, \ldots\}$, the spectrum of the Hamiltonian.

The probability density of a stationary state is

$$|\psi(\mathbf r, t)|^2 = |\phi(\mathbf r)\, e^{-iEt/\hbar}|^2 = |\phi(\mathbf r)|^2,$$

independent of time — which justifies the name. Stationary states are the closest quantum analogues of classical bound orbits: the electron has a definite energy and an unchanging spatial distribution.

!!! note "Why we care about stationary states"
    The Hamiltonian eigenstates form a complete basis for the space of physically allowed wavefunctions (under the conditions of the spectral theorem, which hold for the operators we meet). Any wavefunction can be expanded as
    $$\psi(\mathbf r, t) = \sum_n c_n\, \phi_n(\mathbf r)\, e^{-iE_n t/\hbar},$$
    so solving the TISE *is* solving the TDSE for any initial condition. The whole problem of quantum dynamics reduces to the spectral problem (4.2.6).

## 4.2.4 Expectation values and operators

In quantum mechanics each physical observable $A$ — position, momentum, energy, angular momentum — is represented by a linear operator $\hat A$ acting on wavefunctions. The rule for predicting the average outcome of measuring $A$ on a system in state $\psi$ is

$$\boxed{\; \langle A \rangle = \int \psi^*(\mathbf r)\, \hat A\, \psi(\mathbf r)\, d^3 r. \;} \tag{4.2.7}$$

The integral is over all space, and $\psi$ is assumed normalised. This is the **expectation value** of $\hat A$ in the state $\psi$. Examples:

- Position: $\hat x = x$ (multiplication by $x$), so $\langle x\rangle = \int x\,|\psi|^2\,d^3r$ — the centre of mass of the probability density.
- Momentum: $\hat{\mathbf p} = -i\hbar \nabla$, so $\langle p_x\rangle = -i\hbar\int \psi^* \partial_x \psi \,d^3r$.
- Kinetic energy: $\hat T = \hat p^2/2m = -\frac{\hbar^2}{2m}\nabla^2$.
- Potential energy: $\hat V = V(\mathbf r)$ (multiplication).
- Total energy: $\hat H = \hat T + \hat V$.

!!! warning "Operators do not commute"
    Unlike numbers, operators in general satisfy $\hat A \hat B \neq \hat B \hat A$. The canonical example is position and momentum:
    $$[\hat x, \hat p_x] \equiv \hat x \hat p_x - \hat p_x \hat x = i\hbar. \tag{4.2.8}$$
    This non-commutativity is the algebraic root of the Heisenberg uncertainty principle, $\Delta x\, \Delta p_x \geq \hbar/2$.

## 4.2.5 Hermitian operators have real eigenvalues

For observable quantities (position, energy, …) measurement outcomes must be real numbers. The operators corresponding to observables therefore cannot be arbitrary linear operators; they must have *real* eigenvalues. The relevant condition is that they be **Hermitian** (also called *self-adjoint*).

An operator $\hat A$ is Hermitian if, for all square-integrable $\phi$ and $\psi$,

$$\int \phi^* (\hat A \psi)\, d^3 r = \int (\hat A \phi)^* \psi\, d^3 r. \tag{4.2.9}$$

Equivalently, in bra-ket notation (which we introduce below), $\langle \phi | \hat A | \psi \rangle = \langle \psi | \hat A | \phi \rangle^*$.

**Theorem.** *The eigenvalues of a Hermitian operator are real.*

**Proof.** Let $\hat A \phi = a\phi$ with $\phi \neq 0$. Take the inner product of both sides with $\phi$:

$$\int \phi^*\, \hat A \phi\, d^3 r = a \int \phi^* \phi\, d^3 r = a\, \|\phi\|^2.$$

By Hermiticity (4.2.9) with $\psi = \phi$,

$$\int \phi^*\, \hat A \phi\, d^3 r = \int (\hat A \phi)^*\, \phi\, d^3 r = (a)^* \int \phi^* \phi\, d^3 r = a^*\, \|\phi\|^2.$$

Equating the two expressions and dividing by $\|\phi\|^2 \neq 0$ gives $a = a^*$, hence $a \in \mathbb R$. $\blacksquare$

The Hamiltonian (4.2.2) is Hermitian: $V$ is real (so $V$ acts as a Hermitian multiplication operator), and $-\frac{\hbar^2}{2m}\nabla^2$ is Hermitian under integration by parts, provided the wavefunctions decay at infinity (which they do for bound states). Hence the energy eigenvalues $E_n$ in (4.2.6) are real, as required.

## 4.2.6 Orthogonality of eigenfunctions

A second crucial property: eigenfunctions of a Hermitian operator belonging to *different* eigenvalues are automatically orthogonal.

**Theorem.** *If $\hat A \phi_m = a_m \phi_m$ and $\hat A \phi_n = a_n \phi_n$ with $a_m \neq a_n$, then $\int \phi_m^* \phi_n\, d^3 r = 0$.*

**Proof.** By Hermiticity,

$$\int \phi_m^*\, \hat A \phi_n\, d^3 r = \int (\hat A \phi_m)^*\, \phi_n\, d^3 r.$$

The left side is $a_n \int \phi_m^* \phi_n \, d^3 r$ (since $\hat A \phi_n = a_n \phi_n$). The right side is $a_m^* \int \phi_m^* \phi_n \, d^3 r = a_m \int \phi_m^* \phi_n\, d^3 r$ (using the previous theorem to drop the conjugate on $a_m$). Hence

$$(a_n - a_m) \int \phi_m^* \phi_n \, d^3 r = 0,$$

and since $a_n \neq a_m$ the integral vanishes. $\blacksquare$

Within a degenerate subspace (multiple eigenfunctions sharing the same eigenvalue), one can always *choose* an orthogonal basis by Gram–Schmidt. Combined with normalisation, this gives an **orthonormal** set:

$$\int \phi_m^*(\mathbf r)\, \phi_n(\mathbf r)\, d^3 r = \delta_{mn}. \tag{4.2.10}$$

The completeness statement — that this orthonormal set spans the space of admissible wavefunctions — is the spectral theorem for self-adjoint operators on $L^2$. We will treat it as given.

## 4.2.7 Bra-ket notation

The integrals we have been writing soon become cumbersome, and Paul Dirac's bra-ket notation is much tidier. Define:

- A **ket** $|\psi\rangle$ is an abstract state vector — informally, the wavefunction $\psi$ without committing to a particular coordinate representation.
- A **bra** $\langle\phi|$ is the corresponding linear functional. In the position representation, $\langle\phi|$ acts by $\int \phi^*(\mathbf r) \cdot\, d^3 r$.
- The **inner product** is $\langle\phi|\psi\rangle \equiv \int \phi^*(\mathbf r)\, \psi(\mathbf r)\, d^3 r$.
- An **operator** acts on a ket from the left: $\hat A |\psi\rangle$.
- **Matrix elements** of an operator: $\langle\phi|\hat A|\psi\rangle \equiv \int \phi^*(\mathbf r)\, \hat A\, \psi(\mathbf r)\, d^3 r$.

In this notation, our results compress neatly. Normalisation:

$$\langle\psi|\psi\rangle = 1. \tag{4.2.11}$$

Expectation value of $\hat A$ in state $|\psi\rangle$:

$$\langle A\rangle = \langle \psi|\hat A|\psi\rangle. \tag{4.2.12}$$

Hermiticity:

$$\langle\phi|\hat A|\psi\rangle = \langle\psi|\hat A|\phi\rangle^*. \tag{4.2.13}$$

Orthonormality of eigenstates:

$$\langle \phi_m|\phi_n\rangle = \delta_{mn}. \tag{4.2.14}$$

Completeness:

$$\sum_n |\phi_n\rangle\langle\phi_n| = \hat 1, \tag{4.2.15}$$

where $\hat 1$ is the identity operator; this is sometimes called the *resolution of the identity*. Expanding an arbitrary state,

$$|\psi\rangle = \sum_n c_n |\phi_n\rangle, \quad c_n = \langle\phi_n|\psi\rangle,$$

is then a one-line consequence of (4.2.15).

The eigenvalue equation $\hat H \phi_n = E_n \phi_n$ becomes

$$\hat H |n\rangle = E_n |n\rangle, \tag{4.2.16}$$

where we have shortened $|\phi_n\rangle$ to $|n\rangle$.

## 4.2.8 What we have built

In the space of a few pages we have constructed the entire algebraic framework of non-relativistic quantum mechanics.

- States live in a Hilbert space; in position representation they are square-integrable complex functions $\psi(\mathbf r)$.
- Time evolution is governed by the Schrödinger equation (4.2.1); for time-independent $V$, the dynamics is determined by the eigenstates and eigenvalues of $\hat H$.
- Observables correspond to Hermitian operators; their measurement statistics are encoded by (4.2.7), or in bra-ket form (4.2.12).
- The eigenvalues of a Hermitian operator are real, and the eigenstates belonging to distinct eigenvalues are orthogonal — proved above by direct calculation.
- The mathematical fact that the eigenstates of $\hat H$ form a complete basis reduces the entire problem of quantum dynamics to a *spectral problem*: find the eigenvalues and eigenvectors of $\hat H$.

The rest of the chapter is the systematic exploration of that spectral problem in increasingly realistic settings. In §4.3 we solve it for a single particle in a 1D box — analytically and on the computer. In §4.4 we do the same for the harmonic oscillator. In §4.5 we write down the Hamiltonian for a real solid, and discover that the spectral problem is, in practice, hopeless. The remaining sections of the chapter — and indeed the rest of the book — are devoted to the approximations that put it within reach.

A computational footnote before moving on. Equation (4.2.6) is a *linear* eigenvalue problem — the same kind of problem you met in Chapter 0.3 when diagonalising a 3×3 symmetric matrix. The only essential difference is that the operator $\hat H$ acts on an infinite-dimensional function space rather than a finite-dimensional vector space. In §4.3 we will *discretise* the position coordinate onto a finite grid, at which point $\hat H$ becomes a literal matrix and the Schrödinger equation becomes a problem for `scipy.linalg.eigh`. That is the bridge from quantum mechanics to computational quantum mechanics, and it is shorter than you might think.
