# 4.6 Born–Oppenheimer separation

The full Hamiltonian (4.5.1) couples electrons and nuclei: every electron interacts with every nucleus through the Coulomb term $\hat V_{\mathrm{en}}$, and every nucleus interacts with every other nucleus through $\hat V_{\mathrm{nn}}$. Solving this jointly is hopeless, as the previous section established. The first — and arguably the most important — approximation in all of atomistic simulation is to *decouple* the electronic and nuclear motions on the basis of the enormous mass disparity between them. This is the **Born–Oppenheimer (BO) approximation**, introduced by Max Born and Robert Oppenheimer in 1927.

The BO separation has two pay-offs. First, it shrinks the electronic problem to a well-defined sub-problem: solve the electronic Schrödinger equation with nuclei held fixed. Second — and this is the dramatic one — it produces the **potential energy surface** $E_{\mathrm{BO}}(\mathbf R_1, \ldots, \mathbf R_{N_{\mathrm n}})$, the central object on which every classical molecular dynamics simulation, every geometry optimisation, every transition-state search, and every machine-learning interatomic potential is built. Without BO there is no PES; without a PES there is no atomistic simulation of the kind we will be doing for the rest of the book.

## 4.6.1 Mass and timescale disparity

A proton is approximately 1836 times heavier than an electron. A typical nucleus (carbon, silicon, …) is 10$^4$ – 10$^5$ times heavier. From this single fact a great deal follows.

Suppose we have an electron and a nucleus, both with comparable kinetic energies — as is the case in any bound state where the virial theorem applies. Equipartition gives roughly $\tfrac12 m v^2 \sim k_B T$ or, more relevantly for chemistry, comparable contributions to the total energy. Comparing velocities at equal kinetic energy:

$$\frac{v_{\mathrm{n}}}{v_{\mathrm{e}}} = \sqrt{\frac{m_{\mathrm e}}{M}}\sim \sqrt{\frac{1}{1836}} \approx 0.023. \tag{4.6.1}$$

Nuclei move roughly a hundredth of the speed of electrons. Equivalently, the characteristic *timescale* of nuclear motion is two orders of magnitude longer than that of electronic motion. Typical molecular vibrations are at $\sim 10^{13}$–$10^{14}$ Hz (period $\sim 10$ fs); typical electronic transitions are at $\sim 10^{15}$–$10^{16}$ Hz (period $\sim 0.1$ fs).

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

The gradient acts only on the explicit $\mathbf R$-dependence of $\hat{H}_{\mathrm e}$ (the electron–nucleus and nucleus–nucleus terms), not on the wavefunction. This is what makes BOMD tractable: forces are computed as expectation values, not by numerical differentiation of the energy. We revisit Hellmann–Feynman in §5.5.

For light nuclei (hydrogen, deuterium) the classical approximation can fail measurably, and one must keep the nuclear quantum problem (4.6.12). Quantum effects on nuclear motion show up as zero-point energies, tunnelling, and isotope effects; methods such as path-integral molecular dynamics (PIMD) preserve nuclear quantum effects while still using the BO surface.

## 4.6.6 When Born–Oppenheimer breaks down

The BO approximation is excellent for the vast majority of materials simulations, but it does break down — and where it breaks, interesting physics happens. There are essentially two failure modes.

**Non-adiabatic coupling.** If two electronic states $\psi_0$ and $\psi_1$ have nearly degenerate energies at some nuclear configuration, the matrix elements $\mathbf A_I^{(01)} = \int \psi_0^* \nabla_I \psi_1\, d\mathbf r$ that we dropped become very large (they diverge at exact degeneracy). The electronic state of the system can no longer be assumed to be the ground state at all times; nuclear motion can "kick" the system from $\psi_0$ to $\psi_1$. This is the regime of photochemistry, internal conversion in molecules, and electronic stopping in radiation damage.

**Conical intersections.** Two BO surfaces $E_0(\mathbf R)$ and $E_1(\mathbf R)$ can touch at specific points $\mathbf R^*$ where they meet in a double cone (locally $E_\pm \approx E^* \pm |\mathbf k\cdot \delta \mathbf R|$). At such a point the BO approximation fails catastrophically: the electronic wavefunction is undefined, the geometric phase becomes singular, and proper treatment requires *coupled* dynamics on both surfaces. Conical intersections are the dominant decay channels in photophysics (e.g.\ the cis–trans photoisomerisation of retinal in your eye proceeds through one). They cannot be captured by single-surface BO dynamics.

Methods that go beyond BO — surface-hopping schemes (Tully), Ehrenfest dynamics, exact factorisation, multi-configuration time-dependent Hartree — are an active research area. They are beyond the scope of this book; we will treat BO as exact throughout the remaining chapters.

!!! warning "A practical note"
    Even within ground-state BO dynamics, you should be aware that *most* DFT codes assume the BO approximation silently. When you read in Chapter 6 about "static" DFT, you are computing a single point on $E_0(\mathbf R)$; when you read about *ab initio* MD, you are integrating (4.6.13). Non-adiabatic effects are typically negligible for ground-state thermodynamics and structure, but become important for electron-phonon coupling, metal–insulator transitions, and any phenomenon involving electronic excited states.

## 4.6.7 Looking ahead

The BO approximation has done two things for us. It has reduced the problem of materials simulation to two coupled but separable sub-problems:

1. **Electronic structure.** Given fixed nuclei $\mathbf R$, compute $E_0(\mathbf R)$ and (optionally) the forces $-\nabla_I E_0(\mathbf R)$. This is the domain of Hartree–Fock (§4.7) and DFT (Chapter 5).

2. **Nuclear dynamics.** Given $E_0(\mathbf R)$, propagate the nuclei — classically by (4.6.13) or quantum-mechanically by (4.6.12). This is the domain of geometry optimisation (Chapter 6), molecular dynamics (Chapter 7), phonon analysis, and ML potentials (Chapter 9).

For the electronic problem, the bare wavefunction approach is still exponentially intractable for any non-trivial system; we must approximate further. In §4.7 we introduce the simplest such approximation, Hartree–Fock, and identify what it misses (correlation) — setting the stage for the density-functional revolution of Chapter 5.
