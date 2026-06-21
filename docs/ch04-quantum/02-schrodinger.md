# 4.2 The Schrödinger equation

!!! info "What problem are we solving?"
    In classical mechanics, if you know where a particle is and how fast it
    is moving, Newton's laws tell you everything it will ever do. For an
    electron in an atom that recipe fails completely: electrons show
    interference, tunnelling and discrete energy levels that no force law can
    reproduce. We need a *different* rulebook. This section states that
    rulebook. It introduces the single object that holds everything we are
    allowed to know about a quantum particle — the **wavefunction** $\psi$ —
    and the single equation that tells us how $\psi$ changes in time — the
    **Schrödinger equation**. Once we have those, the rest of the chapter (and
    most of the book) is just *solving* that equation in harder and harder
    settings.

!!! note "Plain-language version — what a wavefunction *is*"
    Forget, for a moment, the words "Hilbert space" and "operator" in the
    section below; we will build up to them slowly. The core idea is small
    enough to say in one breath:

    > A quantum particle does not have a definite position. Instead it carries
    > a **cloud of possibility**, written $\psi(\mathbf r)$, that is spread out
    > over space. Where the cloud is large, the particle is likely to be
    > found; where the cloud is zero, it will never be found.

    The wavefunction $\psi$ *is* that cloud, written as a number attached to
    every point in space. There is one twist that makes quantum mechanics
    strange: the number $\psi(\mathbf r)$ is a **complex number** (it has a
    size and a phase, like an arrow with a length and a direction), and it can
    be positive, negative, or genuinely complex. Only its *size-squared*,
    $|\psi(\mathbf r)|^2$, is a probability. The phase — the "direction of the
    arrow" — never shows up in a single measurement, but it controls how two
    clouds add up when they overlap, and *that* is where interference comes
    from. Hold on to this picture: everything formal below is machinery for
    computing with this cloud.

!!! note "Physical picture"
    Picture a hydrogen atom. Classically you might draw the electron as a dot
    on a circular orbit. Quantum mechanically you should instead picture a
    faint, fuzzy ball of "electron-stuff" centred on the proton — densest near
    the nucleus, fading smoothly outward, with no sharp edge and no dot
    anywhere. That fuzzy ball is $|\psi|^2$. It does not flicker or jiggle: for
    a state of definite energy (a *stationary state*, §4.2.3) the ball just
    sits there, unchanging, even though the underlying complex $\psi$ is
    quietly rotating its phase. When we later compute the energy levels of a
    box, an oscillator, or a whole crystal, we are computing the *shapes* these
    fuzzy clouds are allowed to take and the energy that goes with each shape.

!!! tip "New vocabulary"
    Four words below will carry the whole section. Skim them now; each is
    defined properly where it first does real work, and each links to the
    [beginner glossary](../undergraduate/glossary-for-beginners.md) for a
    slower treatment.

    - **[Wavefunction](../undergraduate/glossary-for-beginners.md)** $\psi$ —
      the complex "cloud of possibility" just described.
    - **[Operator](../undergraduate/glossary-for-beginners.md)** $\hat A$ — a
      *rule that acts on a wavefunction and returns another function*; the
      quantum stand-in for a measurable quantity. "Multiply by $x$" and
      "differentiate" are operators.
    - **[Eigenvalue / eigenstate](../undergraduate/glossary-for-beginners.md)**
      — a special function that an operator merely *rescales* (an eigenstate),
      and the number it is rescaled by (the eigenvalue). The allowed energies
      are the eigenvalues of the energy operator.
    - **Hermitian** — the precise mathematical condition (defined in §4.2.5)
      that guarantees an operator's eigenvalues are *real* numbers — which they
      must be, because a real experiment returns a real number.

We now state the central equation of non-relativistic quantum mechanics. Everything that follows in this book — the band structure of silicon, the binding energy of a benzene molecule, the vibrational spectrum of a zeolite — is a consequence, or an approximation to a consequence, of this single postulate. We will *not* derive it. Schrödinger himself arrived at it by guesswork guided by analogy with classical wave optics, and there is no logical sense in which it can be deduced from earlier physics. It is a postulate, justified solely by the overwhelming agreement of its predictions with experiment.

## 4.2.0 The four postulates of quantum mechanics

Before writing down the equation itself it is useful to lay out the conceptual scaffolding on which all of non-relativistic quantum mechanics rests. There is no universally agreed list of axioms — different books split the postulates differently — but the following four are reasonably standard and span everything we will need for this book.

!!! note "Plain-language version of the four postulates"
    The formal list below uses notation — $|\psi\rangle$, $\hat A$,
    $\langle a_n|\psi\rangle$ — that we only unpack later in the section. On a
    first reading, do not get stuck on the symbols; read the four postulates as
    four plain statements, and come back once §4.2.7 has explained the
    bracket notation.

    1. **A state is a wavefunction.** Everything knowable about the system is
       packed into one object $\psi$ (the cloud).
    2. **Each measurable quantity is an operator.** Position, momentum and
       energy are not numbers attached to the particle; they are *rules*
       ($\hat x$, $\hat p$, $\hat H$) that act on $\psi$.
    3. **Measurement is probabilistic (the Born rule).** You cannot predict a
       single outcome, only the *odds*. The possible outcomes are the
       operator's eigenvalues, and the average over many runs is the
       expectation value (§4.2.4).
    4. **Between measurements $\psi$ evolves smoothly** according to the
       Schrödinger equation — deterministically, with no randomness at all.

    The randomness in quantum mechanics lives entirely in Postulate 3
    (measurement); Postulate 4 (evolution) is as deterministic as Newton.

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

Before unpacking these equations, here is every symbol that appears in them
and in the rest of the section, so that nothing is left undefined. (For a
slower walk through how to read an equation symbol by symbol, see the
[formula reading guide](../undergraduate/formula-reading-guide.md).)

| Symbol | Read it as | Meaning | Units (SI) |
|---|---|---|---|
| $\psi(\mathbf r, t)$ | "psi of r and t" | the wavefunction: a complex number at each point $\mathbf r$ and time $t$ | $\mathrm{m}^{-3/2}$ (in 3D) |
| $i$ | "i" | the imaginary unit, $i^2 = -1$ | — |
| $\hbar$ | "h-bar" | reduced Planck constant, $h/2\pi \approx 1.055\times10^{-34}\ \mathrm{J\,s}$ | J s |
| $m$ | "m" | mass of the particle | kg |
| $\partial/\partial t$ | "partial dee by dee t" | rate of change in time, holding position fixed | $\mathrm{s}^{-1}$ |
| $\nabla^2$ | "del-squared" / Laplacian | $\partial_x^2 + \partial_y^2 + \partial_z^2$, a sum of second spatial derivatives | $\mathrm{m}^{-2}$ |
| $V(\mathbf r, t)$ | "the potential" | potential energy as a function of position (and possibly time) | J |
| $\hat H$ | "H-hat" | the Hamiltonian *operator* — the rule "kinetic + potential energy" | J (acts to give energy) |
| $\hat p = -i\hbar\nabla$ | "p-hat" | the momentum operator | kg m s$^{-1}$ |
| $\hat x$ | "x-hat" | the position operator (multiply by $x$) | m |
| $\langle A\rangle$ | "expectation value of $A$" | the average of measuring $A$ over many identical systems | units of $A$ |
| $E$ | "E" | an energy eigenvalue (allowed energy) | J |

!!! warning "Common misunderstandings — the hat and the cloud"
    - A **hat** ($\hat H$, $\hat p$) marks an *operator*, not a number. $\hat p$
      is the instruction "differentiate and multiply by $-i\hbar$"; it only
      becomes a number after you sandwich it in an expectation value (§4.2.4).
      Writing "$\hat p = mv$" as if it were a velocity times a mass is a
      category error.
    - $\nabla^2$ is **not** a vector and not "$\nabla$ squared as a dot
      product of two arrows you can choose freely" — it is one fixed scalar
      operator, the sum of three second derivatives.
    - The $i$ in front is **not** decorative. Drop it and you get the heat
      (diffusion) equation, whose solutions decay rather than oscillate; the
      $i$ is exactly what turns spreading into wave-like interference.

Several features deserve immediate comment.

- **It is first order in time.** Unlike Newton's equations (second order) or the classical wave equation (also second order), the TDSE needs only an initial condition $\psi(\mathbf r, 0)$ to determine all future evolution. There is no need for an initial "velocity of the wavefunction".

- **It is complex.** The factor of $i$ on the left is essential. A complex-valued wavefunction is not a calculational convenience that can be replaced by a real two-component vector at the end of the day; the phase relationships between different parts of $\psi$ encode interference effects that have been observed in countless experiments.

- **It is linear.** If $\psi_1$ and $\psi_2$ both satisfy (4.2.1), so does $\alpha\psi_1 + \beta\psi_2$ for any complex constants $\alpha, \beta$. This is the principle of superposition, and it is the formal source of every characteristically quantum phenomenon.

- **The Hamiltonian is the energy operator.** Notice that $-\frac{\hbar^2}{2m}\nabla^2$ is what you get if you take the classical kinetic energy $p^2/2m$ and substitute $\mathbf p \to -i\hbar \nabla$. This substitution — momentum becomes a differential operator — is the heuristic device by which Schrödinger guessed his equation, but again, it is no derivation.

!!! note "A plausibility argument, not a derivation"
    For a free particle ($V = 0$) one can check that the plane wave
    $\psi(\mathbf r, t) = \exp[i(\mathbf k\cdot\mathbf r - \omega t)]$
    solves (4.2.1) provided $\hbar\omega = \hbar^2 k^2/2m$, i.e.\ $E = p^2/2m$ with de Broglie's $p = \hbar k$ and the Planck relation $E = \hbar\omega$. This is the standard motivational pattern: the Schrödinger equation is the simplest linear PDE whose plane-wave solutions reproduce the non-relativistic dispersion relation that de Broglie and Einstein together demand. Schrödinger generalised the kinetic-energy term to include $V(\mathbf r)$ by analogy with the classical Hamiltonian.

??? note "Full derivation: the plane wave and the dispersion relation"
    We verify the claim line by line, so that the constant $\hbar^2 k^2/2m$ does
    not appear from nowhere. Take the free-particle TDSE ($V=0$) and the trial
    wave
    $$\psi(\mathbf r,t)=e^{i(\mathbf k\cdot\mathbf r-\omega t)},\qquad
      \mathbf k\cdot\mathbf r = k_x x + k_y y + k_z z.$$

    **Left-hand side.** Differentiating with respect to $t$ pulls down a factor
    $-i\omega$ from the exponent:
    $$i\hbar\,\frac{\partial\psi}{\partial t}
      = i\hbar\,(-i\omega)\,\psi
      = \hbar\omega\,\psi,$$
    using $i\cdot(-i)=-i^2=1$.

    **Right-hand side.** Each spatial derivative pulls down a factor of the
    corresponding $ik$:
    $$\frac{\partial\psi}{\partial x}=ik_x\,\psi,\qquad
      \frac{\partial^2\psi}{\partial x^2}=(ik_x)^2\,\psi=-k_x^2\,\psi,$$
    and likewise for $y$ and $z$. Adding the three,
    $$\nabla^2\psi=-(k_x^2+k_y^2+k_z^2)\,\psi=-k^2\,\psi,\qquad
      k^2\equiv|\mathbf k|^2.$$
    Hence the Hamiltonian acting on $\psi$ gives
    $$\hat H\psi=-\frac{\hbar^2}{2m}\nabla^2\psi
      =-\frac{\hbar^2}{2m}\,(-k^2)\,\psi
      =\frac{\hbar^2 k^2}{2m}\,\psi.$$

    **Match the two sides.** Equation (4.2.1) demands
    $i\hbar\,\partial_t\psi=\hat H\psi$, i.e.
    $$\hbar\omega\,\psi=\frac{\hbar^2 k^2}{2m}\,\psi
      \;\;\Longrightarrow\;\;
      \boxed{\;\hbar\omega=\frac{\hbar^2 k^2}{2m}\;}\tag{4.2.1a}$$
    after cancelling the common non-zero factor $\psi$. This is the
    **dispersion relation**.

    **Read off the physics.** Insert the Planck relation $E=\hbar\omega$ and the
    de Broglie relation $p=\hbar k$ (so $p^2=\hbar^2k^2$):
    $$E=\hbar\omega=\frac{\hbar^2 k^2}{2m}=\frac{(\hbar k)^2}{2m}=\frac{p^2}{2m}.$$
    The plane wave reproduces *exactly* the non-relativistic free-particle
    energy $E=p^2/2m$. That agreement is the whole reason for the particular
    coefficient $-\hbar^2/2m$ in front of $\nabla^2$: it is chosen so that this
    line comes out right.

## 4.2.2 The Born rule

A complex wavefunction is not, by itself, an observable. We have to specify how $\psi$ relates to measurement outcomes. Max Born's 1926 proposal, for which he received the 1954 Nobel Prize, is now universal:

!!! note "Plain-language version — from amplitude to probability"
    The wavefunction $\psi$ is a complex number at each point: it has a *size*
    and a *phase*. You can never measure the phase directly. What an experiment
    can report is a *probability*, and the rule connecting the two is the
    simplest one that could possibly work: take the size of $\psi$ and square
    it. Big $\psi$ (in magnitude) means "likely to be found here"; $\psi=0$
    means "never found here". The technical notation $|\psi|^2=\psi^*\psi$ just
    spells out "size squared": for a complex number $z=a+ib$, $z^*z=a^2+b^2$,
    the squared length of its arrow. The phase cancels out of $\psi^*\psi$
    entirely — which is exactly why a single measurement can never reveal it.

$$\rho(\mathbf r, t) \equiv |\psi(\mathbf r, t)|^2 = \psi^*(\mathbf r, t)\, \psi(\mathbf r, t) \tag{4.2.3}$$

is the **probability density** of finding the particle at position $\mathbf r$ at time $t$. The probability of finding it in a small volume $d^3 r$ around $\mathbf r$ is $\rho(\mathbf r, t)\, d^3 r$.

This is a probabilistic, not deterministic, theory. The Schrödinger equation evolves $\psi$ deterministically, but $\psi$ only tells you the *odds* of various measurement outcomes. A single experiment yields a single result; the predictions of quantum mechanics are statistical and only become sharp upon averaging over many identically-prepared systems.

Because $\rho$ is a probability density, it must integrate to one:

$$\int |\psi(\mathbf r, t)|^2\, d^3 r = 1. \tag{4.2.4}$$

This is the **normalisation condition**. A wavefunction satisfying (4.2.4) is called *normalised*. Notice that any solution of (4.2.1) can be rescaled by a constant without breaking the equation, so we always choose the constant so that (4.2.4) holds.

!!! warning "Common misunderstandings about $\psi$ and $|\psi|^2$"
    These three errors trip up almost everyone on a first pass. They are worth
    reading slowly.

    - **"$\psi$ is a physical, measurable wave, like a ripple on a pond."** No.
      A water wave's height is a real, directly measurable displacement.
      $\psi$ is complex and cannot be measured at all; only $|\psi|^2$ — a
      probability density — connects to experiment. The phase of $\psi$ is
      physically real in its *consequences* (interference), but it is not a
      thing a single detector ever reads off.
    - **"The electron is literally smeared out like jelly, filling the
      cloud."** No. $|\psi|^2$ is a probability density, not a density of
      smeared-out charge that you would find spread thin if you looked. Every
      time you *measure* the position you find the electron at one point. The
      cloud describes the *statistics* of where that one point lands over many
      identical experiments, not a physically dilute electron. (In §4.2.3 we
      will see that for a single electron $-e|\psi|^2$ does behave as an
      effective charge density when it interacts with other particles — but the
      electron is still found whole, never in pieces.)
    - **"Normalising changes the physics."** No. Multiplying $\psi$ by a
      constant changes neither the energy nor any probability *ratio*; (4.2.4)
      just fixes the overall scale so the total probability is exactly $1$.
      Likewise, multiplying by a pure phase $e^{i\alpha}$ leaves $|\psi|^2$
      untouched and so describes the very same physical state.

!!! question "Check yourself — the Born rule"
    Pause here before reading on. Answer these from the Born rule alone.

    1. A particle in 1D has $\psi(x)=A\,e^{-x^2/2a^2}$ (a Gaussian). At which
       single point is the particle *most* likely to be found, and where is it
       *never* found?
    2. Two students prepare the same state but write it as $\psi$ and
       $-\psi$ respectively. Will any experiment ever distinguish their
       predictions?
    3. Someone hands you $\psi(x)$ with $\int|\psi|^2\,dx = 4$. What number
       must you divide $\psi$ by to normalise it, and why is that allowed?
    4. True or false: "Between measurements, the particle is really at a
       definite but unknown place, and $|\psi|^2$ is just our ignorance of it."

    ??? note "Hint"
        For (1), $|\psi|^2 \propto e^{-x^2/a^2}$ — where is *that* largest?
        For (3), normalisation needs $\int |\psi/c|^2\,dx = 1$, so
        $\int |\psi|^2\,dx = |c|^2$. For (4), recall that the phase of $\psi$
        produces interference; a "hidden definite position" cannot.

    ??? success "Answer"
        1. Most likely at $x=0$, where $|\psi|^2=|A|^2$ is largest. The Gaussian
           is never exactly zero for finite $x$, so strictly there is no point
           where it *cannot* be found; it only $\to 0$ as $x\to\pm\infty$. (If
           the question were a particle-in-a-box ground state $\sin(\pi x/L)$,
           it would be most likely at the centre and *never* at the walls
           $x=0,L$.)
        2. No. $|-\psi|^2 = |\psi|^2$, so every probability — and hence every
           prediction — is identical. An overall sign (or phase) is physically
           invisible.
        3. Divide by $c=2$ (i.e. by $\sqrt{4}$). Normalising needs
           $\int|\psi/c|^2 = \tfrac{1}{|c|^2}\int|\psi|^2 = \tfrac{4}{|c|^2}=1$,
           so $|c|=2$. It is allowed because rescaling a solution of the linear
           equation (4.2.1) gives another valid solution describing the same
           physics.
        4. **False** — and this is the deep point. The interference experiments
           that the phase of $\psi$ predicts cannot be reproduced by any
           "definite but unknown position" picture. Quantum indeterminacy is
           not mere ignorance; it is a property of nature. (This is the content
           of Bell's theorem, beyond our scope here.)

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

!!! note "Why this step? — solving the time ODE"
    The equation $i\hbar\,\dot f = E f$ rearranges to $\dot f/f = E/(i\hbar) =
    -iE/\hbar$ (using $1/i = -i$). Any equation of the form "rate of change
    divided by value equals a constant $\lambda$" has the exponential solution
    $f(t)=f(0)\,e^{\lambda t}$, here with $\lambda = -iE/\hbar$. We absorb the
    constant $f(0)$ into the spatial part $\phi$, so we may set $f(0)=1$ and
    write $f(t)=e^{-iEt/\hbar}$. You can check it by substituting back:
    $\dot f = (-iE/\hbar)e^{-iEt/\hbar}$, so $i\hbar\,\dot f =
    i\hbar(-iE/\hbar)f = E f$. ✓ Because $|e^{-iEt/\hbar}|=1$ for real $E$, this
    time factor is a pure *phase* that rotates but never grows or decays — the
    hallmark of a state that conserves probability.

and the spatial part is the **time-independent Schrödinger equation** (TISE):

$$\boxed{\; \hat{H}\, \phi(\mathbf r) = E\, \phi(\mathbf r) \;} \tag{4.2.6}$$

with $\hat{H} = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf r)$.

Equation (4.2.6) is an **eigenvalue equation**: we seek functions $\phi$ that are mapped, by the Hamiltonian, into multiples of themselves. The multipliers $E$ are the allowed *energies* of the system. For typical bound-state problems they form a discrete set $\{E_0, E_1, E_2, \ldots\}$, the spectrum of the Hamiltonian.

!!! warning "Common misunderstanding — eigenvalue vs eigenfunction"
    The two halves of $\hat H\phi = E\phi$ are different *kinds* of object and
    must never be confused.

    - $\phi(\mathbf r)$ is the **eigenfunction** (or eigenstate): a whole
      *function* of position — the shape of the cloud.
    - $E$ is the **eigenvalue**: a single *number* — the energy that goes with
      that shape.

    The equation says: "applying $\hat H$ to this particular function gives the
    *same* function back, merely scaled by the number $E$." A function is not a
    number; "the eigenvalue $\phi$" or "the energy is $\phi(\mathbf r)$" are
    both nonsense. They come in matched *pairs* $(E_n, \phi_n)$: the lowest
    energy $E_0$ has its own ground-state function $\phi_0$, the next energy
    $E_1$ its own $\phi_1$, and so on. When `scipy.linalg.eigh` returns
    `evals, evecs` (§4.3), `evals[n]` is the eigenvalue $E_n$ and the column
    `evecs[:, n]` is the (discretised) eigenfunction $\phi_n$ — two outputs,
    two meanings.

The probability density of a stationary state is

$$|\psi(\mathbf r, t)|^2 = |\phi(\mathbf r)\, e^{-iEt/\hbar}|^2 = |\phi(\mathbf r)|^2,$$

independent of time — which justifies the name. Stationary states are the closest quantum analogues of classical bound orbits: the electron has a definite energy and an unchanging spatial distribution.

!!! note "Why we care about stationary states"
    The Hamiltonian eigenstates form a complete basis for the space of physically allowed wavefunctions (under the conditions of the spectral theorem, which hold for the operators we meet). Any wavefunction can be expanded as
    $$\psi(\mathbf r, t) = \sum_n c_n\, \phi_n(\mathbf r)\, e^{-iE_n t/\hbar},$$
    so solving the TISE *is* solving the TDSE for any initial condition. The whole problem of quantum dynamics reduces to the spectral problem (4.2.6).

## 4.2.4 Expectation values and operators

In quantum mechanics each physical observable $A$ — position, momentum, energy, angular momentum — is represented by a linear operator $\hat A$ acting on wavefunctions. The rule for predicting the average outcome of measuring $A$ on a system in state $\psi$ is

!!! note "Plain-language version — what an expectation value is"
    An **expectation value** $\langle A\rangle$ is *not* the value you get in
    one measurement (that would be a single eigenvalue, with a random spread).
    It is the *average* you would get if you measured $A$ on a great many
    systems all prepared in the same state $\psi$, and then took the mean. The
    word "expectation" is borrowed from probability theory, where it means
    exactly "average over the distribution". Equation (4.2.7) is the recipe for
    that average, written for quantum mechanics.

!!! tip "New vocabulary"
    - **Expectation value** $\langle A\rangle$ — the statistical mean of $A$
      over many identical measurements. See the
      [beginner glossary](../undergraduate/glossary-for-beginners.md).
    - **Observable** — any physical quantity you can measure, hence (Postulate
      2) any Hermitian operator: position, momentum, energy, spin.

The formula (4.2.7) is not a new postulate; it follows from the Born rule.

??? note "Full derivation: why $\langle A\rangle = \int\psi^*\hat A\,\psi$"
    **Start with position, where the Born rule gives the answer directly.** The
    average position is the ordinary statistical mean of $x$ weighted by the
    probability density $|\psi(x)|^2$ (Born rule):
    $$\langle x\rangle = \int x\,\underbrace{|\psi(x)|^2}_{\text{prob.\ density}}\,dx
      = \int \psi^*(x)\,x\,\psi(x)\,dx,$$
    where in the last step we wrote $|\psi|^2=\psi^*\psi$ and slid the (real)
    number $x$ between them. Since the position operator $\hat x$ *is*
    "multiply by $x$", this is exactly $\int\psi^*\,\hat x\,\psi\,dx$.

    **Generalise to functions of position.** The same logic gives, for any
    function $V(x)$,
    $$\langle V(x)\rangle = \int V(x)\,|\psi|^2\,dx = \int\psi^*\,V(x)\,\psi\,dx,$$
    because $V(x)$ is also just multiplication.

    **Momentum needs the operator form.** For momentum there is no "$|\psi|^2$
    in $x$" that does the job, because momentum is not a function of position.
    The clean way is to use the momentum-space wavefunction
    $\tilde\psi(p)=\langle p|\psi\rangle$, whose $|\tilde\psi(p)|^2$ *is* the
    probability density in momentum (Postulate 1). The average momentum is then
    $$\langle p\rangle = \int p\,|\tilde\psi(p)|^2\,dp.$$
    Transforming this single line back to position space (a standard Fourier
    calculation: $p$ in momentum space becomes $-i\hbar\,\partial_x$ in position
    space) turns it into
    $$\langle p\rangle = \int \psi^*(x)\,(-i\hbar\,\partial_x)\,\psi(x)\,dx
      = \int \psi^*\,\hat p\,\psi\,dx.$$

    **Read off the general pattern.** In every case the average of the
    observable equals "$\psi^*$ times *the operator acting on* $\psi$,
    integrated". Promoting this observed pattern to a rule for *every*
    observable gives (4.2.7). This is consistent with — indeed equivalent to —
    Postulate 3: if you expand $\psi$ in the eigenstates of $\hat A$, equation
    (4.2.7) reproduces "$\sum_n a_n\times(\text{probability }|c_n|^2)$", the
    textbook definition of a weighted average (we prove this once the spectral
    machinery is in place, §4.2.6–4.2.7).

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

!!! info "What problem are we solving? — why energies come out real"
    We have built a theory in which the wavefunction is *complex* and operators
    can be complex too ($\hat p = -i\hbar\,\partial_x$ has an explicit $i$). Yet
    every real measurement — an energy on a meter, a position on a screen —
    returns an ordinary *real* number. Something must guarantee that the
    numbers the theory predicts (the eigenvalues, and the expectation values)
    are real, never $3+2i$ joules. This subsection identifies that something:
    the operator must be **Hermitian**, and we *prove* that Hermiticity forces
    real eigenvalues and real averages. This is the precise reason energies are
    real.

!!! note "Plain-language version"
    "Hermitian" is the quantum version of "symmetric" for a matrix. A real
    symmetric matrix has real eigenvalues — you may have seen this in linear
    algebra (Chapter 0). Hermitian operators are the infinite-dimensional
    generalisation, and the same conclusion holds: real eigenvalues. The two
    short proofs below are the whole story; everything physical (real energies,
    real averages, orthogonal states) flows from them.

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

**Corollary.** *The expectation value of a Hermitian operator is real, in any state — even a state that is not an eigenstate.* This is the second half of "why energies are real": not only are the *possible outcomes* (eigenvalues) real, but the *average* of any measurement is real too.

**Proof.** Let $\hat A$ be Hermitian and $|\psi\rangle$ any normalised state. Apply the Hermiticity condition (4.2.9) with $\phi=\psi$:

$$\langle A\rangle = \int \psi^*\,(\hat A\psi)\,d^3r
   = \int (\hat A\psi)^*\,\psi\,d^3r
   = \left(\int \psi^*\,(\hat A\psi)\,d^3r\right)^{\!*}
   = \langle A\rangle^*,$$

where the middle equality is Hermiticity and the third uses
$\int (\hat A\psi)^*\psi = \big(\int \psi^*(\hat A\psi)\big)^*$ (taking the
complex conjugate of a number swaps which factor wears the star). A number
equal to its own complex conjugate is real, so $\langle A\rangle\in\mathbb R$.
$\blacksquare$

In particular $\langle H\rangle$ — the average energy of *any* wavefunction, eigenstate or not — is real. This is exactly the property that makes the variational principle of §4.7 work: it minimises a real energy.

!!! question "Check yourself — Hermiticity and reality"
    1. The matrix $\hat A=\begin{pmatrix}0&-i\\ i&0\end{pmatrix}$ (this is the
       Pauli matrix $\sigma_y$). Is it Hermitian? What does that tell you about
       its eigenvalues *before* you compute them?
    2. Why does the proof that eigenvalues are real *not* work for a
       non-Hermitian operator — which exact line fails?
    3. The corollary says $\langle A\rangle$ is real even when $\psi$ is not an
       eigenstate of $\hat A$. Reconcile this with the fact that a single
       measurement of $A$ on such a $\psi$ gives a random eigenvalue.

    ??? note "Hint"
        For (1), a matrix is Hermitian iff it equals its conjugate-transpose:
        take the transpose, then conjugate every entry, and compare. For (2),
        look for the step that used "$a=a^*$". For (3), "real" and "definite"
        are different properties — an average can be sharp even when individual
        outcomes scatter.

    ??? success "Answer"
        1. Transpose gives $\begin{pmatrix}0&i\\ -i&0\end{pmatrix}$; conjugating
           gives back $\begin{pmatrix}0&-i\\ i&0\end{pmatrix}=\hat A$. So $\hat
           A$ *is* Hermitian, and its eigenvalues are guaranteed real without
           any calculation (they are in fact $\pm1$).
        2. The proof equated $a\|\phi\|^2$ with $a^*\|\phi\|^2$ *using
           Hermiticity*. Without Hermiticity the two integrals
           $\int\phi^*\hat A\phi$ and $\int(\hat A\phi)^*\phi$ need not be
           equal, so the conclusion $a=a^*$ does not follow and eigenvalues may
           be complex.
        3. The expectation value is an *average* over many runs; it can be a
           sharp real number even though each individual run returns one of
           several (real) eigenvalues at random. For example, the $+x$ spin
           state has $\langle\sigma_z\rangle=0$ exactly, yet every single
           measurement of $\sigma_z$ yields $+1$ or $-1$, never $0$ — the
           average of equally likely $\pm1$ is $0$.

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

!!! tip "Where this appears later"
    Two pillars of this section reappear as the foundation of nearly everything
    that follows.

    - The fact that $\langle H\rangle$ is **real** and bounded below (proved in
      §4.2.5) is what makes the **variational principle** of §4.7 legitimate:
      one can *minimise* the energy of a trial wavefunction and never fall
      below the true ground state. That principle is, in turn, the engine of
      [Chapter 5 (DFT)](../ch05-dft/index.md) — density functional theory
      recasts the whole many-electron problem as a variational minimisation of
      an energy functional, and the self-consistent loop you will meet there is
      exactly "minimise a real $\langle H\rangle$".
    - The **eigenvalue equation** $\hat H\phi=E\phi$ becomes a literal matrix
      diagonalisation in §4.3, the same diagonalisation that produces the
      orbitals and band energies of [Chapter 5 (DFT)](../ch05-dft/index.md) and
      the normal-mode frequencies of [Chapter 7 (molecular dynamics)](../ch07-md/index.md).

## 4.2.9 Ehrenfest's theorem and the classical limit

Before turning to specific solvable problems, one further consequence of the formalism deserves note: the expectation values of position and momentum obey *exactly* Newton's equations. This is **Ehrenfest's theorem**, and it provides a rigorous bridge from quantum to classical mechanics.

Differentiate $\langle \hat x\rangle = \langle\psi(t)|\hat x|\psi(t)\rangle$ with respect to time using the TDSE. After standard manipulation (writing out the time derivatives of $|\psi\rangle$ and $\langle\psi|$, applying the Hamiltonian),

$$\frac{d\langle\hat A\rangle}{dt} = \frac{1}{i\hbar}\,\langle[\hat A, \hat H]\rangle + \left\langle\frac{\partial \hat A}{\partial t}\right\rangle. \tag{4.2.17}$$

??? note "Full derivation: the Ehrenfest equation (4.2.17) and the two commutators"
    **The general rate equation.** Differentiate the expectation value
    $\langle A\rangle=\langle\psi|\hat A|\psi\rangle$ by the product rule. Three
    terms appear — one for the bra, one for the operator, one for the ket:
    $$\frac{d}{dt}\langle\psi|\hat A|\psi\rangle
      = \Big(\tfrac{d}{dt}\langle\psi|\Big)\hat A|\psi\rangle
      + \langle\psi|\Big(\tfrac{\partial\hat A}{\partial t}\Big)|\psi\rangle
      + \langle\psi|\hat A\Big(\tfrac{d}{dt}|\psi\rangle\Big).$$
    From the TDSE, $\tfrac{d}{dt}|\psi\rangle=\tfrac{1}{i\hbar}\hat H|\psi\rangle$.
    Taking the complex conjugate (and using that $\hat H$ is Hermitian) gives the
    bra version, $\tfrac{d}{dt}\langle\psi|=-\tfrac{1}{i\hbar}\langle\psi|\hat H$.
    Substitute both:
    $$\frac{d\langle A\rangle}{dt}
      = -\tfrac{1}{i\hbar}\langle\psi|\hat H\hat A|\psi\rangle
      + \Big\langle\tfrac{\partial\hat A}{\partial t}\Big\rangle
      + \tfrac{1}{i\hbar}\langle\psi|\hat A\hat H|\psi\rangle.$$
    Collect the first and third terms:
    $$\frac{d\langle A\rangle}{dt}
      = \tfrac{1}{i\hbar}\langle\psi|(\hat A\hat H-\hat H\hat A)|\psi\rangle
      + \Big\langle\tfrac{\partial\hat A}{\partial t}\Big\rangle
      = \tfrac{1}{i\hbar}\langle[\hat A,\hat H]\rangle
      + \Big\langle\tfrac{\partial\hat A}{\partial t}\Big\rangle,$$
    which is (4.2.17). The whole time-dependence of any average is carried by a
    single commutator with $\hat H$.

    **The commutator $[\hat x,\hat p^2]=2i\hbar\,\hat p$.** Use the algebraic
    identity $[\hat x,\hat B\hat C]=[\hat x,\hat B]\hat C+\hat B[\hat x,\hat C]$
    with $\hat B=\hat C=\hat p$, together with the canonical relation
    $[\hat x,\hat p]=i\hbar$ from (4.2.8):
    $$[\hat x,\hat p^2]=[\hat x,\hat p]\hat p+\hat p[\hat x,\hat p]
      = i\hbar\,\hat p+\hat p\,(i\hbar)=2i\hbar\,\hat p.$$
    Then, since $\hat H=\hat p^2/2m+V$ and $\hat x$ commutes with $V(\hat x)$,
    $$\frac{d\langle\hat x\rangle}{dt}
      =\frac{1}{i\hbar}\langle[\hat x,\hat H]\rangle
      =\frac{1}{i\hbar}\frac{1}{2m}\langle[\hat x,\hat p^2]\rangle
      =\frac{1}{i\hbar}\frac{1}{2m}\,2i\hbar\,\langle\hat p\rangle
      =\frac{\langle\hat p\rangle}{m}.$$

    **The commutator $[\hat p,V(\hat x)]=-i\hbar\,V'(\hat x)$.** Act on a test
    function $\psi$ and use the product rule, exactly as in the canonical
    commutator derivation:
    $$[\hat p,V]\psi
      =-i\hbar\,\partial_x(V\psi)-V(-i\hbar\,\partial_x\psi)
      =-i\hbar\big(V'\psi+V\,\partial_x\psi\big)+i\hbar\,V\,\partial_x\psi
      =-i\hbar\,V'\,\psi.$$
    Since $\hat p$ commutes with $\hat p^2$, only the potential contributes to
    $[\hat p,\hat H]$, giving
    $$\frac{d\langle\hat p\rangle}{dt}
      =\frac{1}{i\hbar}\langle[\hat p,\hat H]\rangle
      =\frac{1}{i\hbar}\langle[\hat p,V]\rangle
      =\frac{1}{i\hbar}\,(-i\hbar)\langle V'(\hat x)\rangle
      =-\langle V'(\hat x)\rangle.$$

    **Combine.** Differentiate the first result once more and insert the second:
    $$m\frac{d^2\langle\hat x\rangle}{dt^2}
      =m\frac{d}{dt}\frac{\langle\hat p\rangle}{m}
      =\frac{d\langle\hat p\rangle}{dt}
      =-\langle V'(\hat x)\rangle,$$
    which is (4.2.18) — the quantum centroid obeying Newton's second law.

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

!!! question "Check yourself — the whole section in five questions"
    If you can answer these, you have the working content of Section 4.2.

    1. State, in one sentence each, what the **time-dependent** and the
       **time-independent** Schrödinger equations are *for* — what does each one
       compute?
    2. Where does the time-independent equation *come from*? Name the
       mathematical trick and the single assumption it requires about $V$.
    3. Why must observables be represented by **Hermitian** operators? Give the
       two consequences this guarantees.
    4. A wavefunction is $\psi(x)=N\sin(\pi x/L)$ on $0\le x\le L$ and zero
       outside. Without computing $N$, where is the particle most likely to be
       found, and what is $\langle x\rangle$ by symmetry?
    5. What single computational object does the entire section reduce the
       problem of quantum dynamics to, and which `scipy` routine will solve it
       in §4.3?

    ??? note "Hint"
        (2) revisit §4.2.3. (3) revisit §4.2.5 — reality of eigenvalues *and*
        reality of $\langle A\rangle$. (4) $|\psi|^2\propto\sin^2(\pi x/L)$ is
        symmetric about the midpoint $x=L/2$.

    ??? success "Answer"
        1. The **TDSE** (4.2.1) evolves a given wavefunction forward in time:
           given $\psi(\mathbf r,0)$ it predicts $\psi(\mathbf r,t)$. The
           **TISE** (4.2.6) finds the allowed *energies* and stationary-state
           *shapes* — the eigenvalues $E_n$ and eigenfunctions $\phi_n$ of
           $\hat H$.
        2. **Separation of variables**: write $\psi(\mathbf r,t)=\phi(\mathbf
           r)f(t)$, substitute, and divide. It works only if $V$ is
           **time-independent**, so that the two sides can depend on $\mathbf r$
           and $t$ separately and equal a common constant $E$.
        3. Because a real experiment returns a real number. Hermiticity
           guarantees (i) the eigenvalues (possible single-measurement outcomes)
           are real, and (ii) the expectation value $\langle A\rangle$ (the
           average) is real in *every* state.
        4. $|\psi|^2\propto\sin^2(\pi x/L)$ peaks at the centre $x=L/2$, so the
           particle is most likely found there; by the symmetry of $|\psi|^2$
           about $x=L/2$, $\langle x\rangle=L/2$.
        5. A **matrix eigenvalue problem** $\hat H\phi=E\phi$ — solved once
           $\hat H$ is discretised on a grid (§4.3) by
           `scipy.linalg.eigh`.
