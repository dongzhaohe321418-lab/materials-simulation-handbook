# 4.2 The Schrödinger equation

We now state the central equation of non-relativistic quantum mechanics. Everything that follows in this book — the band structure of silicon, the binding energy of a benzene molecule, the vibrational spectrum of a zeolite — is a consequence, or an approximation to a consequence, of this single postulate. We will *not* derive it. Schrödinger himself arrived at it by guesswork guided by analogy with classical wave optics, and there is no logical sense in which it can be deduced from earlier physics. It is a postulate, justified solely by the overwhelming agreement of its predictions with experiment.

## 4.2.0 The four postulates of quantum mechanics

Before writing down the equation itself it is useful to lay out the conceptual scaffolding on which all of non-relativistic quantum mechanics rests. There is no universally agreed list of axioms — different books split the postulates differently — but the following four are reasonably standard and span everything we will need for this book.

**Postulate 1 (states).** The state of an isolated physical system is fully described by a unit vector $|\psi\rangle$ in a complex Hilbert space $\mathcal H$. Two vectors that differ only by an overall phase $|\psi\rangle$ and $e^{i\alpha}|\psi\rangle$ describe the same physical state. In the position representation $|\psi\rangle$ becomes the wavefunction $\psi(\mathbf r) = \langle \mathbf r | \psi\rangle$, a square-integrable function of position; in the momentum representation it becomes $\tilde\psi(\mathbf p) = \langle \mathbf p|\psi\rangle$, its Fourier transform. The two representations carry the same information.

**Postulate 2 (observables).** Every measurable physical quantity $A$ is represented by a Hermitian linear operator $\hat A$ acting on $\mathcal H$. Position is $\hat x$ (multiplication by $x$ in the position representation), momentum is $\hat p = -i\hbar\,\partial_x$, energy is the Hamiltonian $\hat H$, angular momentum is $\hat L = \hat r\times\hat p$, and so on. The spectrum of $\hat A$ — its eigenvalues — is the set of possible outcomes of a measurement of $A$.

**Postulate 3 (measurement).** A measurement of $A$ on a system in state $|\psi\rangle$ yields one of the eigenvalues $a_n$ of $\hat A$ with probability $|\langle a_n | \psi\rangle|^2$, where $|a_n\rangle$ is the corresponding eigenvector. Immediately after the measurement the system is in the state $|a_n\rangle$ (the "projection postulate", or "collapse"). The average outcome over many identically prepared systems is the expectation value $\langle A\rangle = \langle\psi|\hat A|\psi\rangle$.

**Postulate 4 (time evolution).** Between measurements the state vector evolves deterministically and unitarily according to the Schrödinger equation,
$$i\hbar\,\partial_t |\psi(t)\rangle = \hat H |\psi(t)\rangle,$$
where $\hat H$ is the Hamiltonian operator of the system.

These four postulates are the entirety of the theory. They are not derived from classical physics; classical physics is a *limit* of them. The whole of this chapter is the unpacking of Postulates 1, 2, and 4 in concrete cases (Postulate 3 we shall use silently, when we compute expectation values and transition probabilities). When you read DFT in Chapter 5 or coupled-cluster theory in the chemistry literature, you are reading consequences of these four statements.

!!! warning "What about measurement?"
    Postulate 3 — the collapse postulate — is famously contentious. The "Copenhagen interpretation", "many-worlds", "decoherent histories", and the Bohmian "pilot wave" all agree on the predictive content of quantum mechanics (the Born rule, Postulate 3) but disagree about *what is really happening*. For materials simulation, none of this matters: we predict expectation values, never the outcomes of individual measurements. Whatever one's interpretive preference, the equations are the same.

## 4.2.1 The time-dependent Schrödinger equation

For a single non-relativistic particle of mass $m$ moving in a potential $V(\mathbf r, t)$, the state is described by a complex-valued wavefunction $\psi(\mathbf r, t)$, and the wavefunction evolves according to the **time-dependent Schrödinger equation** (TDSE):

$$\boxed{\; i\hbar \frac{\partial \psi(\mathbf r, t)}{\partial t} = \hat{H}\, \psi(\mathbf r, t) \;} \tag{4.2.1}$$

where $\hat{H}$ is the *Hamiltonian operator*,

$$\hat{H} = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf r, t). \tag{4.2.2}$$

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

### Derivation of the continuity equation

The continuity equation is sufficiently central that we work it out step by step. Start from the TDSE and its complex conjugate,

$$i\hbar\, \partial_t \psi = -\frac{\hbar^2}{2m}\nabla^2 \psi + V\psi,$$
$$-i\hbar\, \partial_t \psi^* = -\frac{\hbar^2}{2m}\nabla^2 \psi^* + V\psi^*.$$

Multiply the first by $\psi^*$, the second by $\psi$, and subtract:

$$i\hbar\, (\psi^* \partial_t \psi + \psi\, \partial_t\psi^*) = -\frac{\hbar^2}{2m}\left(\psi^*\nabla^2\psi - \psi\nabla^2\psi^*\right).$$

!!! note "Why this step?"
    The $V\psi$ and $V\psi^*$ terms cancel because $V$ is real. This is why the continuity equation requires Hermiticity of the Hamiltonian: a non-Hermitian $V$ (an "optical potential" with imaginary part) would describe absorption or creation of probability, and (4.2.4) would not be preserved.

The left side is $i\hbar\, \partial_t (\psi^* \psi) = i\hbar\, \partial_t \rho$. The right side rearranges, using the identity $\psi^*\nabla^2\psi - \psi\nabla^2\psi^* = \nabla\cdot(\psi^*\nabla\psi - \psi\nabla\psi^*)$, to $-(\hbar^2/2m)\,\nabla\cdot(\psi^*\nabla\psi - \psi\nabla\psi^*)$. Dividing by $i\hbar$ and rearranging,

$$\partial_t \rho + \nabla\cdot \mathbf j = 0, \qquad \mathbf j = \frac{\hbar}{2mi}\left(\psi^* \nabla\psi - \psi \nabla\psi^*\right) = \frac{1}{m}\,\text{Re}\!\left[\psi^*(-i\hbar\nabla)\psi\right].$$

The second form has a transparent reading: $\mathbf j = \text{Re}\,\langle \psi|\hat{\mathbf v}|\psi\rangle$ where $\hat{\mathbf v} = \hat{\mathbf p}/m$ is the velocity operator. The probability current is "density times velocity", as in classical fluid mechanics.

!!! example "A free-particle plane wave"
    For $\psi(\mathbf r, t) = e^{i(\mathbf k\cdot\mathbf r - \omega t)}$, $|\psi|^2 = 1$ is uniform and $\mathbf j = \hbar\mathbf k/m = \mathbf p/m$ is constant: probability flows at the group velocity. This is the quantum-mechanical analogue of a uniform classical beam.

## 4.2.3 Stationary states and the time-independent Schrödinger equation

A vast amount of practical quantum mechanics — almost everything we do in materials physics — boils down to looking for *stationary states*: solutions whose probability density does not change with time. We try a separable ansatz,

$$\psi(\mathbf r, t) = \phi(\mathbf r)\, f(t),$$

and substitute into (4.2.1) with a time-independent potential $V(\mathbf r)$. The TDSE becomes

$$i\hbar\, \phi(\mathbf r)\, \dot f(t) = f(t)\, \hat{H} \phi(\mathbf r),$$

!!! note "Why this step?"
    On the left, $\partial_t [\phi(\mathbf r) f(t)] = \phi(\mathbf r)\, \dot f(t)$ because $\phi$ is time-independent. On the right, $\hat H$ contains only spatial derivatives and the time-independent potential $V(\mathbf r)$, so it commutes with $f(t)$ and we can pull $f(t)$ out: $\hat H[\phi f] = f(t)\,\hat H \phi(\mathbf r)$. The separable ansatz is consistent only if $V$ has no explicit time-dependence — which is *almost always* the case in materials science, since static lattices and frozen nuclei (BO approximation, §4.6) make the Hamiltonian time-independent.

and dividing by $\phi f$,

$$i\hbar\, \frac{\dot f(t)}{f(t)} = \frac{\hat{H} \phi(\mathbf r)}{\phi(\mathbf r)}.$$

The left side depends only on $t$, the right only on $\mathbf r$; both must therefore equal a common constant, which we call $E$.

!!! note "Why this step? — the separation-of-variables argument"
    This is a standard PDE trick worth pausing on. If $g(t) = h(\mathbf r)$ holds for *every* $\mathbf r$ and *every* $t$ then, pick any fixed $\mathbf r_0$ and vary $t$: the left side changes but the right side is constant. Conclude that the left side is independent of $t$. Pick any fixed $t_0$ and vary $\mathbf r$: now the right side is independent of $\mathbf r$. Both sides are therefore constant. We name this constant $E$ in anticipation of its physical interpretation as the energy eigenvalue. This gives two equations. The time part is

$$i\hbar\, \dot f(t) = E\, f(t),$$

a first-order linear ODE with the immediate solution

$$f(t) = e^{-iEt/\hbar}, \tag{4.2.5}$$

and the spatial part is the **time-independent Schrödinger equation** (TISE):

$$\boxed{\; \hat{H}\, \phi(\mathbf r) = E\, \phi(\mathbf r) \;} \tag{4.2.6}$$

with $\hat{H} = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf r)$.

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
- Total energy: $\hat{H} = \hat T + \hat V$.

!!! warning "Operators do not commute"
    Unlike numbers, operators in general satisfy $\hat A \hat B \neq \hat B \hat A$. The canonical example is position and momentum:
    $$[\hat x, \hat p_x] \equiv \hat x \hat p_x - \hat p_x \hat x = i\hbar. \tag{4.2.8}$$
    This non-commutativity is the algebraic root of the Heisenberg uncertainty principle, $\Delta x\, \Delta p_x \geq \hbar/2$.

### Derivation of the canonical commutator

It is worth seeing why $[\hat x, \hat p_x] = i\hbar$ explicitly. Act with the commutator on an arbitrary smooth test function $\psi(x)$:

$$[\hat x, \hat p_x]\,\psi = \hat x\,(-i\hbar\, \partial_x \psi) - (-i\hbar\, \partial_x)\,(\hat x \psi).$$

Use the product rule on the second term: $\partial_x(x\psi) = \psi + x\,\partial_x\psi$. Then

$$[\hat x, \hat p_x]\,\psi = -i\hbar\, x\,\partial_x\psi + i\hbar\,(\psi + x\partial_x\psi) = i\hbar\,\psi.$$

Since this holds for every $\psi$, the operator identity $[\hat x, \hat p_x] = i\hbar\,\hat I$ follows.

!!! note "Why this step?"
    The non-trivial commutator is forced by representing $\hat p$ as $-i\hbar\,\partial_x$. The product rule is exactly the leakage: applying $\partial_x$ to $x\psi$ produces both $x\partial_x\psi$ (the "expected" term) and $\psi$ (the residue, the source of the $i\hbar$). This single fact — the position–momentum commutator — propagates through every quantum calculation. It is, in a meaningful sense, the *whole content* of canonical quantisation.

The general uncertainty relation now follows by a few lines of inequality manipulation (the **Robertson uncertainty relation**): for any two Hermitian operators $\hat A, \hat B$ and any state $|\psi\rangle$,
$$\Delta A\,\Delta B \geq \tfrac12 |\langle[\hat A, \hat B]\rangle|, \qquad (\Delta A)^2 \equiv \langle\hat A^2\rangle - \langle\hat A\rangle^2.$$
For position and momentum, $\langle [\hat x, \hat p_x]\rangle = i\hbar$, so $\Delta x\,\Delta p_x \geq \hbar/2$. The uncertainty principle is not a statement about experimental clumsiness; it is a theorem in linear algebra.

## 4.2.5 Hermitian operators have real eigenvalues

For observable quantities (position, energy, …) measurement outcomes must be real numbers. The operators corresponding to observables therefore cannot be arbitrary linear operators; they must have *real* eigenvalues. The relevant condition is that they be **Hermitian** (also called *self-adjoint*).

An operator $\hat A$ is Hermitian if, for all square-integrable $\phi$ and $\psi$,

$$\int \phi^* (\hat A \psi)\, d^3 r = \int (\hat A \phi)^* \psi\, d^3 r. \tag{4.2.9}$$

Equivalently, in bra-ket notation (which we introduce below), $\langle \phi | \hat A | \psi \rangle = \langle \psi | \hat A | \phi \rangle^*$.

!!! note "Why this step? — what Hermiticity really demands"
    Condition (4.2.9) says that the operator $\hat A$ is the same whether we let it act to the right (on $\psi$) or to the left (on $\phi$). It is the operator generalisation of "$A$ equals its conjugate transpose" for matrices. For position-representation operators, Hermiticity often follows from integration by parts. For example, $-i\hbar\,\partial_x$ is Hermitian because
    $$\int \phi^*(-i\hbar\partial_x \psi)\,dx = -i\hbar [\phi^*\psi]_{-\infty}^{\infty} + \int (i\hbar\partial_x\phi)^* \psi\,dx,$$
    and the boundary term vanishes for square-integrable $\phi,\psi$. The kinetic-energy operator $-\hbar^2 \nabla^2/(2m)$ is Hermitian by a second integration by parts. Note that the *boundary conditions* matter: on a finite interval with rigid walls (particle in a box, §4.3) Hermiticity is preserved because $\psi = 0$ at the walls; with the wrong boundary conditions Hermiticity fails and energies cease to be real.

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

A worked example will help. Consider a three-state system (a $\Lambda$-system, common in quantum optics) with orthonormal basis $\{|1\rangle, |2\rangle, |3\rangle\}$. The most general state is

$$|\psi\rangle = c_1|1\rangle + c_2|2\rangle + c_3|3\rangle,$$

with $\sum_n |c_n|^2 = 1$. The expansion coefficients are recovered by projection: $c_n = \langle n|\psi\rangle$. This is the *defining* operation of Dirac notation — apply a bra to a ket and you get a complex number, the amplitude. Apply a ket to a bra (the other order, $|n\rangle\langle m|$) and you get an operator, called an *outer product*. The completeness relation (4.2.15) is the statement that the sum of outer products $|n\rangle\langle n|$ over a complete orthonormal set acts as the identity:

$$\hat I |\psi\rangle = \sum_n |n\rangle\langle n|\psi\rangle = \sum_n c_n |n\rangle = |\psi\rangle. \quad\checkmark$$

This is how we *insert* a basis. Faced with a matrix element $\langle \phi|\hat A|\psi\rangle$, sandwich completeness in the middle:

$$\langle\phi|\hat A|\psi\rangle = \sum_{m,n} \langle\phi|m\rangle\langle m|\hat A|n\rangle\langle n|\psi\rangle = \sum_{m,n} \phi_m^* A_{mn} \psi_n,$$

which is just $\boldsymbol\phi^\dagger \mathbf A \boldsymbol\psi$ in matrix notation. Quantum mechanics is linear algebra; Dirac notation is the typography that makes that fact visible.

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

The eigenvalue equation $\hat{H} \phi_n = E_n \phi_n$ becomes

$$\hat{H} |n\rangle = E_n |n\rangle, \tag{4.2.16}$$

where we have shortened $|\phi_n\rangle$ to $|n\rangle$.

### Worked example: a spin-1/2 particle

To exercise the Dirac notation on the simplest possible Hilbert space, consider a spin-1/2 system: $\dim \mathcal H = 2$, basis vectors $|\!\uparrow\rangle$ and $|\!\downarrow\rangle$. An arbitrary state is

$$|\psi\rangle = \alpha\,|\!\uparrow\rangle + \beta\,|\!\downarrow\rangle, \qquad |\alpha|^2 + |\beta|^2 = 1,$$

and operators are $2\times 2$ Hermitian matrices. The three components of the spin operator in this basis are $\hat S_a = (\hbar/2)\,\hat\sigma_a$ with the **Pauli matrices**

$$\hat\sigma_x = \begin{pmatrix} 0 & 1 \\ 1 & 0\end{pmatrix}, \quad \hat\sigma_y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}, \quad \hat\sigma_z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}.$$

These satisfy $\hat\sigma_a^2 = \hat I$, $[\hat\sigma_a, \hat\sigma_b] = 2i\,\varepsilon_{abc}\,\hat\sigma_c$ and $\{\hat\sigma_a, \hat\sigma_b\} = 2\delta_{ab}\,\hat I$.

Take the specific state $|\psi\rangle = (|\!\uparrow\rangle + |\!\downarrow\rangle)/\sqrt 2$ — spin pointing along $+x$. Compute the expectation values:

$$\langle\hat\sigma_z\rangle = \langle\psi|\hat\sigma_z|\psi\rangle = \tfrac12 (1, 1)\begin{pmatrix}1 & 0\\ 0 & -1\end{pmatrix}\begin{pmatrix} 1\\ 1\end{pmatrix} = \tfrac12(1 - 1) = 0,$$

$$\langle\hat\sigma_x\rangle = \tfrac12(1, 1)\begin{pmatrix} 0 & 1\\ 1 & 0\end{pmatrix}\begin{pmatrix} 1\\ 1\end{pmatrix} = \tfrac12(1 + 1) = 1,$$

$$\langle\hat\sigma_y\rangle = \tfrac12(1, 1)\begin{pmatrix} 0 & -i\\ i & 0\end{pmatrix}\begin{pmatrix} 1\\ 1\end{pmatrix} = \tfrac12(-i + i) = 0.$$

So $\langle\hat{\mathbf S}\rangle = (\hbar/2)\hat{\mathbf x}$, consistent with the state pointing along $+x$. The variance in $\sigma_z$ is $\langle\hat\sigma_z^2\rangle - \langle\hat\sigma_z\rangle^2 = 1 - 0 = 1$: measuring $\sigma_z$ on this state yields $+1$ or $-1$ with equal probability, and the spread $(\Delta\sigma_z)^2 = 1$ is maximal — entirely consistent with the uncertainty relation $\Delta\sigma_y\, \Delta\sigma_z \geq |\langle\hat\sigma_x\rangle|$.

!!! tip "Why spin-1/2 is the cleanest pedagogical example"
    In a two-dimensional Hilbert space every concept of quantum mechanics — superposition, measurement statistics, non-commuting observables, expectation values, time evolution — appears in its sharpest possible form. The integrals of position-space wave mechanics are replaced by finite $2\times 2$ matrix products. The Bloch sphere makes the geometry visible. The same algebra runs the spin part of every electronic-structure code, the qubit of every quantum computer, and the NMR experiments we will simulate in Chapter 11.

## 4.2.8 What we have built

In the space of a few pages we have constructed the entire algebraic framework of non-relativistic quantum mechanics.

- States live in a Hilbert space; in position representation they are square-integrable complex functions $\psi(\mathbf r)$.
- Time evolution is governed by the Schrödinger equation (4.2.1); for time-independent $V$, the dynamics is determined by the eigenstates and eigenvalues of $\hat{H}$.
- Observables correspond to Hermitian operators; their measurement statistics are encoded by (4.2.7), or in bra-ket form (4.2.12).
- The eigenvalues of a Hermitian operator are real, and the eigenstates belonging to distinct eigenvalues are orthogonal — proved above by direct calculation.
- The mathematical fact that the eigenstates of $\hat{H}$ form a complete basis reduces the entire problem of quantum dynamics to a *spectral problem*: find the eigenvalues and eigenvectors of $\hat{H}$.

The rest of the chapter is the systematic exploration of that spectral problem in increasingly realistic settings. In §4.3 we solve it for a single particle in a 1D box — analytically and on the computer. In §4.4 we do the same for the harmonic oscillator. In §4.5 we write down the Hamiltonian for a real solid, and discover that the spectral problem is, in practice, hopeless. The remaining sections of the chapter — and indeed the rest of the book — are devoted to the approximations that put it within reach.

A computational footnote before moving on. Equation (4.2.6) is a *linear* eigenvalue problem — the same kind of problem you met in Chapter 0.3 when diagonalising a 3×3 symmetric matrix. The only essential difference is that the operator $\hat{H}$ acts on an infinite-dimensional function space rather than a finite-dimensional vector space. In §4.3 we will *discretise* the position coordinate onto a finite grid, at which point $\hat{H}$ becomes a literal matrix and the Schrödinger equation becomes a problem for `scipy.linalg.eigh`. That is the bridge from quantum mechanics to computational quantum mechanics, and it is shorter than you might think.

## 4.2.9 Ehrenfest's theorem and the classical limit

Before turning to specific solvable problems, one further consequence of the formalism deserves note: the expectation values of position and momentum obey *exactly* Newton's equations. This is **Ehrenfest's theorem**, and it provides a rigorous bridge from quantum to classical mechanics.

Differentiate $\langle \hat x\rangle = \langle\psi(t)|\hat x|\psi(t)\rangle$ with respect to time using the TDSE. After standard manipulation (writing out the time derivatives of $|\psi\rangle$ and $\langle\psi|$, applying the Hamiltonian),

$$\frac{d\langle\hat A\rangle}{dt} = \frac{1}{i\hbar}\,\langle[\hat A, \hat H]\rangle + \left\langle\frac{\partial \hat A}{\partial t}\right\rangle. \tag{4.2.17}$$

For $\hat A = \hat x$ with no explicit time dependence,

$$\frac{d\langle\hat x\rangle}{dt} = \frac{1}{i\hbar}\langle[\hat x, \hat H]\rangle = \frac{1}{i\hbar}\,\frac{1}{2m}\langle[\hat x, \hat p^2]\rangle = \frac{\langle\hat p\rangle}{m},$$

using $[\hat x, \hat p^2] = 2i\hbar\,\hat p$. Similarly $d\langle\hat p\rangle/dt = -\langle V'(\hat x)\rangle$. Together,

$$m\frac{d^2\langle\hat x\rangle}{dt^2} = -\langle V'(\hat x)\rangle. \tag{4.2.18}$$

For a narrow wavepacket where $\langle V'(\hat x)\rangle \approx V'(\langle\hat x\rangle)$, this is Newton's second law applied to the wavepacket centroid. **The centre of mass of a quantum particle obeys classical mechanics**, to the extent that the wavepacket is narrow compared to the length scale over which the force varies. This is the rigorous content of the correspondence principle introduced in §4.1.

!!! example "When Ehrenfest fails"
    The approximation $\langle V'(\hat x)\rangle \approx V'(\langle\hat x\rangle)$ is exact only for $V$ linear or quadratic in $x$ — i.e.\ for free particles, uniform fields, and harmonic oscillators. For *any other* potential the wavepacket eventually spreads, $\langle V'(\hat x)\rangle$ deviates from $V'(\langle\hat x\rangle)$, and the classical trajectory of the centroid diverges from the true quantum evolution. This is the deep reason why a literal "classical limit" of the Schrödinger equation is subtle: pointwise convergence requires careful control of wavepacket dynamics, formalised in semiclassical analysis (WKB, stationary phase, geometric optics).

## 4.2.10 Looking ahead

The next three sections solve the TISE in increasingly realistic settings:

- §4.3 (particle in a box) — the simplest bound-state problem, solved twice (analytically and on a grid).
- §4.4 (harmonic oscillator) — the universal model for any potential near a minimum, including molecular vibrations and phonons.
- §4.5 onwards — many electrons, where the formalism still applies but the dimensionality explodes.

Throughout, the structure is the same: write the Hamiltonian, identify the boundary conditions, diagonalise. The framework of this section — wavefunctions, Hermitian operators, Dirac notation — is what makes that procedure mechanical.
