# 4.6 Born–Oppenheimer separation

The full Hamiltonian (4.5.1) couples electrons and nuclei: every electron interacts with every nucleus through the Coulomb term $\hat V_{\mathrm{en}}$, and every nucleus interacts with every other nucleus through $\hat V_{\mathrm{nn}}$. Solving this jointly is hopeless, as the previous section established. The first — and arguably the most important — approximation in all of atomistic simulation is to *decouple* the electronic and nuclear motions on the basis of the enormous mass disparity between them. This is the **Born–Oppenheimer (BO) approximation**, introduced by Max Born and Robert Oppenheimer in 1927.

The BO separation has two pay-offs. First, it shrinks the electronic problem to a well-defined sub-problem: solve the electronic Schrödinger equation with nuclei held fixed. Second — and this is the dramatic one — it produces the **potential energy surface** $E_{\mathrm{BO}}(\mathbf R_1, \ldots, \mathbf R_{N_{\mathrm n}})$, the central object on which every classical molecular dynamics simulation, every geometry optimisation, every transition-state search, and every machine-learning interatomic potential is built. Without BO there is no PES; without a PES there is no atomistic simulation of the kind we will be doing for the rest of the book.

## 4.6.1 Mass and timescale disparity

A proton is approximately 1836 times heavier than an electron. A typical nucleus (carbon, silicon, …) is 10$^4$ – 10$^5$ times heavier. From this single fact a great deal follows.

Suppose we have an electron and a nucleus, both with comparable kinetic energies — as is the case in any bound state where the virial theorem applies. Equipartition gives roughly $\tfrac12 m v^2 \sim k_B T$ or, more relevantly for chemistry, comparable contributions to the total energy. Comparing velocities at equal kinetic energy:

$$\frac{v_{\mathrm{n}}}{v_{\mathrm{e}}} = \sqrt{\frac{m_{\mathrm e}}{M}}\sim \sqrt{\frac{1}{1836}} \approx 0.023. \tag{4.6.1}$$

Nuclei move roughly a hundredth of the speed of electrons. Equivalently, the characteristic *timescale* of nuclear motion is two orders of magnitude longer than that of electronic motion. Typical molecular vibrations are at $\sim 10^{13}$–$10^{14}$ Hz (period $\sim 10$ fs); typical electronic transitions are at $\sim 10^{15}$–$10^{16}$ Hz (period $\sim 0.1$ fs).

!!! example "A concrete timescale calculation"
    The C–H stretching vibration in methane has wavenumber $\tilde\nu \approx 3000$ cm$^{-1}$. Its period is
    $$T_{\mathrm n} = \frac{1}{c\,\tilde\nu} = \frac{1}{(3\times 10^{10}\,\text{cm/s})(3000\,\text{cm}^{-1})} \approx 1.1\times 10^{-14}\ \mathrm{s} = 11\ \mathrm{fs}.$$
    A core 1s electron in carbon has orbital energy $\approx -290$ eV, giving a characteristic period
    $$T_{\mathrm e} = h/|E| = (4.14\times 10^{-15}\,\text{eV s})/(290\,\text{eV}) \approx 1.4\times 10^{-17}\ \mathrm{s} = 14\ \mathrm{as}.$$
    The ratio is $T_{\mathrm n}/T_{\mathrm e} \approx 800$. Valence electrons are slower than core electrons (smaller $|E|$, longer period); using a 1 eV valence-band timescale gives $T_{\mathrm e} \approx 4$ fs and a ratio of only $\sim 3$. The BO approximation is most accurate when the relevant *electronic* states are well-separated in energy from the ground state (so the electronic dynamics is fast); it begins to fail near level crossings, which we discuss in §4.6.6.

The physical picture is clear. Electrons, fast and light, follow the nuclei *instantaneously*. Every time a nucleus moves a small amount, the electron cloud adjusts to its new equilibrium configuration almost immediately, on a timescale invisibly short compared to the nuclear motion. From the electrons' perspective, the nuclei are essentially stationary external charges. From the nuclei's perspective, the electrons provide an averaged force.

This is the same logic that lets us treat a slow-moving ship as instantaneously surrounded by an equilibrium pattern of water waves: as long as the ship moves much more slowly than the wave propagation, the wave field has time to relax to its quasi-static configuration. The BO approximation makes this precise.

## 4.6.2 Setting up the separation

We write the full molecular/solid Hamiltonian (4.5.1) as

$$\hat{H} = \hat T_{\mathrm n}(\mathbf R) + \hat{H}_{\mathrm e}(\mathbf r; \mathbf R), \tag{4.6.2}$$

where $\hat T_{\mathrm n}(\mathbf R) = -\sum_I \hbar^2 \nabla_I^2/(2M_I)$ is the nuclear kinetic energy and

$$\hat{H}_{\mathrm e}(\mathbf r; \mathbf R) = \hat T_{\mathrm e} + \hat V_{\mathrm{ee}} + \hat V_{\mathrm{en}}(\mathbf r, \mathbf R) + \hat V_{\mathrm{nn}}(\mathbf R) \tag{4.6.3}$$

is everything else. We have grouped $\hat V_{\mathrm{nn}}$ with the electronic part because it depends on $\mathbf R$ only (not on $\mathbf r$); it contributes a constant for any fixed nuclear configuration, simply shifting the electronic eigenvalues.

**Step 1: the electronic problem at fixed nuclei.** For each nuclear configuration $\mathbf R \equiv (\mathbf R_1, \ldots, \mathbf R_{N_{\mathrm n}})$, solve the *electronic* eigenvalue equation

$$\boxed{\; \hat{H}_{\mathrm e}(\mathbf r; \mathbf R)\, \psi_k(\mathbf r; \mathbf R) = E_k(\mathbf R)\, \psi_k(\mathbf r; \mathbf R). \;} \tag{4.6.4}$$

Here the nuclear coordinates $\mathbf R$ appear as parameters (note the semicolon): $\hat{H}_{\mathrm e}$ depends on them but does not differentiate with respect to them. For each $\mathbf R$ this is a many-electron problem of the kind discussed in §4.5, with $\mathbf R$-dependent external potential and Coulomb repulsion. The eigenfunctions $\psi_k(\mathbf r; \mathbf R)$ are called the **electronic states** at nuclear configuration $\mathbf R$.

For each $k$, the eigenvalue $E_k(\mathbf R)$, considered as a function of $\mathbf R$, is the **Born–Oppenheimer potential energy surface** of the $k$th electronic state. The ground-state PES $E_0(\mathbf R)$ is what we normally mean when we say "the energy of the molecule as a function of geometry"; it is what gets optimised, sampled by molecular dynamics, and learned by interatomic potentials.

**Step 2: the BO ansatz.** Suppose the electrons stay on a single electronic surface (typically the ground state) throughout the nuclear motion. Then we can write the *total* wavefunction as a product,

$$\boxed{\; \Psi(\mathbf r, \mathbf R) \approx \chi(\mathbf R)\, \psi_0(\mathbf r; \mathbf R), \;} \tag{4.6.5}$$

where $\chi(\mathbf R)$ is a nuclear wavefunction yet to be determined.

!!! note "Why this step? — the ansatz is *not* exact"
    The exact total wavefunction has the form $\Psi(\mathbf r, \mathbf R) = \sum_n \chi_n(\mathbf R)\,\psi_n(\mathbf r; \mathbf R)$, a sum over *all* electronic states $n$, each with its own nuclear amplitude $\chi_n(\mathbf R)$ (this is the **Born expansion**). The ansatz (4.6.5) truncates this sum after one term. Truncation is justified provided the electronic ground state is well-separated from excited states (a finite gap) and the nuclear motion does not "kick" the system into excited states (low temperature, or equivalently slow nuclei). We will identify exactly when this fails in §4.6.6.

The Born expansion would be exact if we kept all terms, but it does not help — we are back to the full coupled electron–nucleus problem. The pay-off of BO is that *one term suffices*, by a factor of $\sim m_{\mathrm e}/M$, for most chemistry. The next correction, the diagonal Born–Oppenheimer correction ($B_I$ kept, $\mathbf A_I^{(0n)}$ neglected for $n \neq 0$), is computable and is sometimes included in high-accuracy spectroscopy of light molecules.

## 4.6.3 Deriving the nuclear equation

Substitute (4.6.5) into the full time-independent Schrödinger equation $\hat{H} \Psi = E \Psi$:

$$\left[\hat T_{\mathrm n} + \hat{H}_{\mathrm e}\right] \chi(\mathbf R) \psi_0(\mathbf r; \mathbf R) = E\, \chi(\mathbf R) \psi_0(\mathbf r; \mathbf R). \tag{4.6.6}$$

The electronic part is straightforward, since $\psi_0$ is an eigenstate of $\hat{H}_{\mathrm e}$:

$$\hat{H}_{\mathrm e}\, \chi(\mathbf R)\, \psi_0(\mathbf r; \mathbf R) = E_0(\mathbf R)\, \chi(\mathbf R)\, \psi_0(\mathbf r; \mathbf R). \tag{4.6.7}$$

The nuclear kinetic term needs care. The nuclear Laplacian acts on *both* factors of (4.6.5):

$$-\frac{\hbar^2}{2M_I} \nabla_I^2 \!\left[\chi(\mathbf R)\, \psi_0(\mathbf r; \mathbf R)\right] = -\frac{\hbar^2}{2M_I}\!\left[\psi_0\, \nabla_I^2 \chi + 2(\nabla_I \chi)\cdot(\nabla_I \psi_0) + \chi\, \nabla_I^2 \psi_0\right]. \tag{4.6.8}$$

Project (4.6.6) onto $\psi_0$ by multiplying both sides by $\psi_0^*(\mathbf r; \mathbf R)$ and integrating over the electronic coordinates $\mathbf r$:

$$-\sum_I \frac{\hbar^2}{2M_I}\left[\nabla_I^2 \chi + 2 \mathbf A_I(\mathbf R) \cdot \nabla_I \chi + B_I(\mathbf R)\, \chi\right] + E_0(\mathbf R)\, \chi = E\, \chi, \tag{4.6.9}$$

with the electronic matrix elements

$$\mathbf A_I(\mathbf R) \equiv \int \psi_0^* \nabla_I \psi_0\, d\mathbf r, \qquad B_I(\mathbf R) \equiv \int \psi_0^* \nabla_I^2 \psi_0\, d\mathbf r. \tag{4.6.10}$$

These extra terms — $\mathbf A_I$ (a "geometric vector potential") and $B_I$ (a "diagonal correction") — encode how the electronic wavefunction *changes* as the nuclei move. They are present because $\psi_0(\mathbf r; \mathbf R)$ depends on $\mathbf R$.

The Born–Oppenheimer approximation, in its strict form, is the statement that these terms are negligible compared to the leading $\nabla_I^2 \chi$:

$$\mathbf A_I \approx 0, \qquad B_I \approx 0. \tag{4.6.11}$$

Why is this justified? Both $\mathbf A_I$ and $B_I$ involve derivatives of the electronic wavefunction with respect to the *slow* nuclear coordinate, with no factor of $M_I^{-1/2}$ to enhance them. The electronic wavefunction changes smoothly with $\mathbf R$ except at special points (see §4.6.6 below), and the change is of order unity over a typical nuclear length scale $\sim 1$ Å. Meanwhile $\nabla_I^2 \chi$ is large: the nuclear wavefunction varies on the typical nuclear de Broglie wavelength, $\sim 0.1$ Å, so $\nabla_I^2 \chi/\chi \sim 100\ \mathrm{\AA}^{-2}$. The neglected terms are smaller than the kept term by typically $\sim m_{\mathrm e}/M$ — exactly the small parameter we identified at the start of the section.

With (4.6.11) the nuclear equation simplifies dramatically:

$$\boxed{\; \left[-\sum_I \frac{\hbar^2}{2M_I}\nabla_I^2 + E_0(\mathbf R)\right] \chi(\mathbf R) = E\, \chi(\mathbf R). \;} \tag{4.6.12}$$

This is a Schrödinger equation for the *nuclei* alone, with the BO surface $E_0(\mathbf R)$ playing the role of the potential. The electrons have been integrated out, leaving behind their average influence as a function of $\mathbf R$.

## 4.6.4 What the BO approximation has done for us

Three statements summarise the achievement.

1. **The full $\Psi(\mathbf r, \mathbf R)$ has been factorised** into an electronic part $\psi_0(\mathbf r; \mathbf R)$ (the ground-state electron cloud at fixed nuclei) and a nuclear part $\chi(\mathbf R)$ (the wavefunction of the nuclei moving on the BO surface).

2. **The electronic problem becomes parametric.** We solve $\hat{H}_{\mathrm e}\psi_0 = E_0 \psi_0$ once for each nuclear configuration. We never need to track time-dependence of the electronic state; it follows the nuclei adiabatically.

3. **The PES $E_0(\mathbf R)$ becomes a function of the $3N_{\mathrm n}$ nuclear coordinates only.** This is the central object of atomistic simulation. Geometry optimisation = find a local minimum of $E_0(\mathbf R)$. MD = integrate Newton's equations $M_I \ddot{\mathbf R}_I = -\nabla_I E_0(\mathbf R)$. Reaction pathways = trace minimum-energy paths on $E_0$. Vibrational analysis = diagonalise $\nabla\nabla E_0$ at a minimum (recall §4.4.5). Machine-learning potentials = learn $E_0(\mathbf R)$ from training data.

In short, the BO approximation cleanly separates the problem of *electronic structure* (Chapter 5: DFT) from the problem of *nuclear motion* (Chapters 7–9: MD, lattice dynamics, ML potentials). Modern computational materials science is largely a matter of computing $E_0$ accurately enough by electronic-structure methods and then using it efficiently in some nuclear-dynamics scheme.

## 4.6.5 Classical limit and the force theorem

If nuclei are sufficiently heavy and the temperature sufficiently high, we can replace the quantum nuclear equation (4.6.12) by its classical limit. Take $\chi(\mathbf R, t) \approx A(\mathbf R, t) \exp[i S(\mathbf R, t)/\hbar]$ and apply the WKB argument: in the limit $\hbar \to 0$, $S$ satisfies the classical Hamilton–Jacobi equation with potential $E_0(\mathbf R)$, and the nuclear trajectories obey

$$M_I \ddot{\mathbf R}_I = -\nabla_I E_0(\mathbf R). \tag{4.6.13}$$

This is the Born–Oppenheimer molecular dynamics equation — Newton's second law with forces derived from the BO surface. It is the foundation of *ab initio* molecular dynamics (Chapter 7).

The forces $\mathbf F_I = -\nabla_I E_0(\mathbf R)$ can be computed efficiently using the **Hellmann–Feynman theorem**: for a normalised electronic eigenstate $\psi_0$,

$$\nabla_I E_0(\mathbf R) = \int \psi_0^*\, (\nabla_I \hat{H}_{\mathrm e})\, \psi_0\, d\mathbf r. \tag{4.6.14}$$

The gradient acts only on the explicit $\mathbf R$-dependence of $\hat{H}_{\mathrm e}$ (the electron–nucleus and nucleus–nucleus terms), not on the wavefunction.

### Proof of the Hellmann–Feynman theorem

Starting from $E_0(\mathbf R) = \langle\psi_0(\mathbf R)|\hat H_{\mathrm e}(\mathbf R)|\psi_0(\mathbf R)\rangle$ with $\langle\psi_0|\psi_0\rangle = 1$, differentiate with respect to a nuclear coordinate $R_I^\alpha$ (one component of $\mathbf R_I$):

$$\partial_{R_I^\alpha} E_0 = \langle\partial_{R_I^\alpha}\psi_0|\hat H_{\mathrm e}|\psi_0\rangle + \langle\psi_0|\partial_{R_I^\alpha}\hat H_{\mathrm e}|\psi_0\rangle + \langle\psi_0|\hat H_{\mathrm e}|\partial_{R_I^\alpha}\psi_0\rangle.$$

In the first and third terms, $\hat H_{\mathrm e}|\psi_0\rangle = E_0|\psi_0\rangle$ (and conjugate), so they combine to

$$E_0 \bigl[\langle\partial_{R_I^\alpha}\psi_0|\psi_0\rangle + \langle\psi_0|\partial_{R_I^\alpha}\psi_0\rangle\bigr] = E_0 \,\partial_{R_I^\alpha}\langle\psi_0|\psi_0\rangle = E_0 \cdot 0 = 0,$$

because the normalisation $\langle\psi_0|\psi_0\rangle = 1$ is a constant. Only the second term — the explicit derivative of the Hamiltonian — survives:

$$\boxed{\;\partial_{R_I^\alpha} E_0 = \langle\psi_0|\partial_{R_I^\alpha}\hat H_{\mathrm e}|\psi_0\rangle.\;}$$

!!! note "Why this step?"
    The vanishing of the "wavefunction-derivative" terms is the content of the theorem and the source of its computational power. It relies on three facts: (i) $\psi_0$ is normalised; (ii) $\psi_0$ is an *exact* eigenstate of $\hat H_{\mathrm e}$ (a variational stationary point); (iii) $\hat H_{\mathrm e}$ is Hermitian. The theorem fails when $\psi_0$ is only approximate (e.g.\ Hartree–Fock with a finite basis), in which case one must add the so-called **Pulay forces** to compensate. Modern plane-wave DFT codes use complete basis sets per *k*-point and the Pulay correction is zero; localised-orbital codes (Gaussian-type orbitals) require it.

The pay-off is enormous. The right-hand side is an integral over the *explicit* $\mathbf R$-dependence of $\hat H_{\mathrm e}$, which (looking back at (4.6.3)) consists only of the electron–nucleus and nucleus–nucleus Coulomb terms — both of them simple to differentiate analytically. Forces on the nuclei are then computed as expectation values of these simple analytical derivatives, never by finite-differencing the (expensive) total energy. Every DFT-based molecular-dynamics calculation, every geometry optimisation, every transition-state search relies on this theorem. We revisit Hellmann–Feynman in §5.5.

### Adiabatic vs diabatic representations

A subtle point worth mentioning: the eigenstate basis $\{\psi_k(\mathbf r; \mathbf R)\}$ we have used is called the **adiabatic basis** — each $\psi_k$ is, at each $\mathbf R$, an exact eigenstate of $\hat H_{\mathrm e}(\mathbf R)$. The adiabatic basis is uniquely defined (up to overall phase) and is what one computes by diagonalising the electronic Hamiltonian. Its drawback is that the basis *rotates* as $\mathbf R$ changes, producing the non-adiabatic coupling terms $\mathbf A_I, B_I$ of (4.6.10).

An alternative is the **diabatic basis**: a smoothly $\mathbf R$-varying set of states $\{\tilde\psi_k\}$ chosen so that the non-adiabatic couplings *vanish*, at the cost that the diabatic states are not (in general) eigenstates of $\hat H_{\mathrm e}$. The off-diagonal matrix elements $\langle\tilde\psi_m|\hat H_{\mathrm e}|\tilde\psi_n\rangle$ are now non-zero and describe transitions between diabatic states. Diabatic representations are convenient for surface-hopping dynamics and for problems involving conical intersections, where the adiabatic basis is singular. The two representations carry the same physics, differently partitioned. For computational materials science the adiabatic representation is overwhelmingly the standard.

For light nuclei (hydrogen, deuterium) the classical approximation can fail measurably, and one must keep the nuclear quantum problem (4.6.12). Quantum effects on nuclear motion show up as zero-point energies, tunnelling, and isotope effects; methods such as path-integral molecular dynamics (PIMD) preserve nuclear quantum effects while still using the BO surface.

!!! example "When nuclei are quantum"
    Three signatures of nuclear quantum effects you should recognise:

    - **Hydrogen-bond geometries.** Water and ice exhibit unusual proton-position distributions even at 0 K because of zero-point motion of the H atom. Classical MD simulations underestimate proton delocalisation.
    - **Hydrogen diffusion in metals.** H atoms tunnel through energy barriers in palladium and niobium at temperatures where classical thermal hopping is negligible. The "tunnelling cross-over temperature" $T_c \sim \hbar\omega_b/(2\pi k_B)$ (where $\omega_b$ is the barrier-top imaginary frequency) is typically 100–200 K for H — well above room temperature for hydrogen in many metals.
    - **Kinetic isotope effects.** Replacing H by D in a chemical reaction changes the rate by a factor of 2–10 even though the BO surface is *identical*. The difference is entirely in the zero-point energy and tunnelling, both of which depend on nuclear mass.

    Path-integral methods imaginary-time-discretise the nuclear partition function and represent each quantum nucleus as a "ring polymer" of $P$ classical replicas. The replicas are coupled by harmonic springs of stiffness $P m k_B^2 T^2/\hbar^2$; in the $P \to \infty$ limit one recovers exact quantum statistical mechanics. We will not need PIMD for most of this book — classical nuclei suffice for the structure and energetics of heavy elements — but it is the right tool for hydrogen-dominated and low-temperature systems.

## 4.6.6 When Born–Oppenheimer breaks down

The BO approximation is excellent for the vast majority of materials simulations, but it does break down — and where it breaks, interesting physics happens. There are essentially two failure modes.

**Non-adiabatic coupling.** If two electronic states $\psi_0$ and $\psi_1$ have nearly degenerate energies at some nuclear configuration, the matrix elements

$$\mathbf A_I^{(01)} = \int \psi_0^*(\mathbf r; \mathbf R)\, \nabla_I \psi_1(\mathbf r; \mathbf R)\, d\mathbf r$$

that we dropped become very large. A useful identity helps see why. Take the gradient $\nabla_I$ of the electronic eigenvalue equation $\hat H_{\mathrm e}\psi_n = E_n\psi_n$, project onto another eigenstate $\psi_m$ with $m \neq n$, and use the Hellmann–Feynman manipulation above:

$$\mathbf A_I^{(mn)} = \frac{\langle\psi_m|\nabla_I \hat H_{\mathrm e}|\psi_n\rangle}{E_n - E_m}, \qquad m \neq n. \tag{4.6.NA}$$

The denominator is the electronic energy gap. When the gap shrinks toward zero, the coupling diverges. The neglect of $\mathbf A_I, B_I$ is justified precisely as long as the gap is large compared to the typical "kick" delivered by nuclear motion, $\hbar\omega_{\mathrm n} \sim 10$ meV. For typical insulators with gaps of several eV, the suppression is by a factor of $\sim 100$, and BO is excellent. For metals (zero gap at the Fermi level) and for photo-excited molecules near a level crossing, BO can fail.

The electronic state of the system can no longer be assumed to be the ground state at all times; nuclear motion can "kick" the system from $\psi_0$ to $\psi_1$. This is the regime of photochemistry, internal conversion in molecules, and electronic stopping in radiation damage.

**Conical intersections.** Two BO surfaces $E_0(\mathbf R)$ and $E_1(\mathbf R)$ can touch at specific points $\mathbf R^*$ where they meet in a double cone (locally $E_\pm \approx E^* \pm |\mathbf k\cdot \delta \mathbf R|$). At such a point the BO approximation fails catastrophically: the electronic wavefunction is undefined, the geometric phase becomes singular, and proper treatment requires *coupled* dynamics on both surfaces. Conical intersections are the dominant decay channels in photophysics (e.g.\ the cis–trans photoisomerisation of retinal in your eye proceeds through one). They cannot be captured by single-surface BO dynamics.

The von Neumann–Wigner *non-crossing rule* says that in a 1D system (only one nuclear coordinate) two states of the same symmetry cannot cross — they form an avoided crossing instead, and BO is locally valid although strained. In higher dimensions exact crossings are generic, occurring on subspaces of co-dimension 2. For a diatomic (1 nuclear coord) this means crossings are absent for states of the same symmetry; for any polyatomic ($\geq 2$ coords) conical intersections are unavoidable wherever excited states are involved. They are not pathologies; they are the *typical* topology of multi-state electronic structure, and modern photochemistry codes (e.g.\ NEWTON-X, SHARC) handle them explicitly.

Methods that go beyond BO — surface-hopping schemes (Tully), Ehrenfest dynamics, exact factorisation, multi-configuration time-dependent Hartree — are an active research area. They are beyond the scope of this book; we will treat BO as exact throughout the remaining chapters.

!!! warning "A practical note"
    Even within ground-state BO dynamics, you should be aware that *most* DFT codes assume the BO approximation silently. When you read in Chapter 6 about "static" DFT, you are computing a single point on $E_0(\mathbf R)$; when you read about *ab initio* MD, you are integrating (4.6.13). Non-adiabatic effects are typically negligible for ground-state thermodynamics and structure, but become important for electron-phonon coupling, metal–insulator transitions, and any phenomenon involving electronic excited states.

## 4.6.6a A worked example: the diatomic potential energy curve

To make BO concrete it is useful to walk through the simplest non-trivial example, the diatomic molecule. Two nuclei of mass $M_A, M_B$ with charges $Z_A, Z_B$ are separated by a distance $R$. Step 1 of the BO programme: at each fixed $R$, solve the electronic Schrödinger equation. The result is a family of curves $E_k(R)$ — the BO potential energy curves of the molecule.

For H$_2$ in its ground state, $E_0(R)$ has the shape sketched in Fig. 4.6.1 (not reproduced here): a steep repulsive wall at small $R$ (nuclear–nuclear repulsion overwhelms electronic attraction), a minimum at $R = R_e \approx 0.74$ Å (the bond length) with depth $D_e \approx 4.75$ eV (the bond dissociation energy), and a flat asymptote at large $R$ (two separated hydrogen atoms). This single curve contains an enormous amount of physics:

- The minimum gives the equilibrium geometry.
- The curvature at the minimum, $V''(R_e)$, gives the harmonic vibrational frequency via $\omega = \sqrt{V''/\mu}$ where $\mu = M_A M_B/(M_A + M_B)$ is the reduced mass. For H$_2$, $\omega \approx 4400$ cm$^{-1}$.
- The depth gives the dissociation energy.
- The shape away from the minimum gives the anharmonicities, the centrifugal distortion constants, and the vibrational overtones.

Step 2: solve the nuclear Schrödinger equation (4.6.12) on this potential. For a diatomic, after separating centre-of-mass motion, this reduces to a 1D problem in $R$ for the radial wavefunction. The vibrational levels $\{E_n^{\mathrm{vib}}\}$ are approximately those of the harmonic oscillator near the minimum, with corrections from anharmonicity (cf.\ §4.4.6).

This two-step procedure — compute $E(R)$ point by point, then solve nuclear dynamics on the resulting curve — is exactly what happens in *every* atomistic simulation. Geometry optimisation finds the minimum of $E(R)$ over the full $3N_{\mathrm n}$-dimensional nuclear coordinate space. Molecular dynamics integrates Newton's equations on $E(R)$. Vibrational analysis computes second derivatives of $E(R)$ at the minimum. The BO surface $E(R)$ is the *interface* between electronic structure (the expensive part) and nuclear motion (the cheap part), and the BO approximation is what makes the interface clean.

!!! tip "Why this matters for machine-learning potentials"
    Modern ML interatomic potentials (Chapter 9) are functions that take nuclear coordinates as input and return the BO energy and forces. They are entirely conceptually downstream of the BO approximation: their *output* is $E_{\mathrm{BO}}(\mathbf R)$ and $-\nabla E_{\mathrm{BO}}(\mathbf R)$, learned from DFT or higher-level reference calculations. An ML potential does not "know" about electrons; it knows about the BO surface that electrons produce. Pre-BO, there is no such thing as an interatomic potential.

## 4.6.7 Looking ahead

The BO approximation has done two things for us. It has reduced the problem of materials simulation to two coupled but separable sub-problems:

1. **Electronic structure.** Given fixed nuclei $\mathbf R$, compute $E_0(\mathbf R)$ and (optionally) the forces $-\nabla_I E_0(\mathbf R)$. This is the domain of Hartree–Fock (§4.7) and DFT (Chapter 5).

2. **Nuclear dynamics.** Given $E_0(\mathbf R)$, propagate the nuclei — classically by (4.6.13) or quantum-mechanically by (4.6.12). This is the domain of geometry optimisation (Chapter 6), molecular dynamics (Chapter 7), phonon analysis, and ML potentials (Chapter 9).

For the electronic problem, the bare wavefunction approach is still exponentially intractable for any non-trivial system; we must approximate further. In §4.7 we introduce the simplest such approximation, Hartree–Fock, and identify what it misses (correlation) — setting the stage for the density-functional revolution of Chapter 5.

!!! tip "The BO approximation in one sentence"
    Because nuclei are heavy, the electrons solve their own Schrödinger equation at every fixed nuclear configuration, producing an effective potential $E_0(\mathbf R)$ on which the nuclei then move according to Newton's (or Schrödinger's) equations. The full coupled $(\mathbf r, \mathbf R)$ problem reduces to two sequential, decoupled problems — and the entire field of atomistic simulation is the systematic exploitation of this fact.

A historical note. Born and Oppenheimer's original 1927 paper treated the corrections as a power series in $\kappa = (m_{\mathrm e}/M)^{1/4}$ and showed that the leading effects on the electronic spectrum scale as $\kappa^2$ (vibrational levels), $\kappa^4$ (rotational levels), and $\kappa^6$ (non-adiabatic mixing). For typical molecules these correspond to electronic transitions in the visible, vibrational transitions in the infrared, and rotational transitions in the microwave — a clean hierarchy that explains why molecular spectroscopy is taught as three separate subjects. The BO expansion is so successful that violations of it (e.g.\ in NO$_2$ or in ozone, where vibronic coupling is strong) are flagged as "anomalies".
