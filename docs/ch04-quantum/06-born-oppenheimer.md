# 4.6 Born–Oppenheimer separation

The full Hamiltonian (4.5.1) couples electrons and nuclei: every electron interacts with every nucleus through the Coulomb term $\hat V_{\mathrm{en}}$, and every nucleus interacts with every other nucleus through $\hat V_{\mathrm{nn}}$. Solving this jointly is hopeless, as the previous section established. The first — and arguably the most important — approximation in all of atomistic simulation is to *decouple* the electronic and nuclear motions on the basis of the enormous mass disparity between them. This is the **Born–Oppenheimer (BO) approximation**, introduced by Max Born and Robert Oppenheimer in 1927.

The BO separation has two pay-offs. First, it shrinks the electronic problem to a well-defined sub-problem: solve the electronic Schrödinger equation with nuclei held fixed. Second — and this is the dramatic one — it produces the **potential energy surface** $E_{\mathrm{BO}}(\mathbf R_1, \ldots, \mathbf R_{N_{\mathrm n}})$, the central object on which every classical molecular dynamics simulation, every geometry optimisation, every transition-state search, and every machine-learning interatomic potential is built. Without BO there is no PES; without a PES there is no atomistic simulation of the kind we will be doing for the rest of the book.

!!! info "What problem are we solving?"
    The exact Schrödinger equation for a molecule or solid contains
    *both* the electrons and the nuclei as quantum particles, all moving
    at once and all pulling on one another through the Coulomb force. The
    wavefunction $\Psi(\mathbf r, \mathbf R)$ then depends on the
    positions of *every* electron ($\mathbf r$) **and** *every* nucleus
    ($\mathbf R$) simultaneously — far too many coupled variables to
    solve. We want to split this single impossible problem into two
    smaller, sequential problems: first work out what the electrons do
    while the nuclei sit still, then let the nuclei move under the
    averaged influence of the electrons. This section shows exactly when
    that split is allowed and what we throw away to make it.

!!! note "Plain-language version"
    The trick rests on one fact: electrons are roughly two thousand times
    lighter than nuclei, so they move roughly a hundred times faster. From
    the electrons' point of view the nuclei are almost frozen; from the
    nuclei's point of view the electrons are a fast-moving blur that has
    already settled into place wherever the nuclei happen to be. So we
    *freeze the nuclei*, solve for the fast electrons, and read off their
    energy as a function of where we froze the nuclei. That energy,
    plotted against nuclear positions, is the potential energy surface the
    nuclei then roll around on. "Fast electrons, slow nuclei" is the whole
    idea in four words.

!!! note "Physical picture"
    Imagine pulling two bonded atoms slowly apart. At every separation the
    electron cloud between them re-shapes itself *instantly* to the
    lowest-energy arrangement for that separation — it never lags behind.
    Because the cloud is always in its electronic ground state for the
    current geometry, its energy is a well-defined number $E_0(\mathbf R)$
    that depends only on where the nuclei are, not on how fast they got
    there. The nuclei then feel a force equal to minus the slope of this
    energy landscape, exactly as a ball feels a force down a hill. The
    word "adiabatic" means precisely this: the electrons follow the slow
    change without ever jumping to an excited state.

!!! tip "New vocabulary"
    - **Born–Oppenheimer approximation** — the assumption that electrons
      adjust instantaneously to the nuclear positions, letting us solve
      the electronic and nuclear problems one after the other.
    - **Potential energy surface (PES)** — the electronic ground-state
      energy as a function of the nuclear positions, $E_0(\mathbf R)$. See
      the [beginner glossary](../undergraduate/glossary-for-beginners.md).
    - **Adiabatic** — a change slow enough that a system stays in the same
      (instantaneous) eigenstate throughout; here, the electrons stay in
      their ground state as the nuclei creep along.
    - **Parametric dependence** — when a quantity depends on a variable
      that is held fixed rather than differentiated. We write
      $\psi_0(\mathbf r; \mathbf R)$ with a semicolon to flag that
      $\mathbf R$ is a *parameter*, not a coordinate the electronic
      operator acts on. Terms like *Hamiltonian*, *operator*,
      *eigenvalue* and *wavefunction* are in the
      [beginner glossary](../undergraduate/glossary-for-beginners.md).

Before the algebra begins, here is every symbol used in this section in one place.

| Symbol | Meaning | Units (SI) |
|---|---|---|
| $\mathbf r$ | collective coordinates of *all* electrons, $(\mathbf r_1, \ldots)$ | m |
| $\mathbf R$ | collective coordinates of *all* nuclei, $(\mathbf R_1, \ldots, \mathbf R_{N_{\mathrm n}})$ | m |
| $\Psi(\mathbf r, \mathbf R)$ | full electron-plus-nucleus wavefunction | — |
| $\psi_k(\mathbf r; \mathbf R)$ | $k$-th electronic eigenstate at *fixed* nuclei $\mathbf R$ | — |
| $\chi(\mathbf R)$ | nuclear wavefunction | — |
| $\hat T_{\mathrm n}$ | nuclear kinetic-energy operator | J |
| $\hat T_{\mathrm e}$ | electronic kinetic-energy operator | J |
| $\hat H_{\mathrm e}(\mathbf r; \mathbf R)$ | electronic Hamiltonian at fixed nuclei (includes $\hat V_{\mathrm{nn}}$) | J |
| $E_k(\mathbf R)$ | $k$-th electronic eigenvalue; $E_0$ is the ground-state PES | J |
| $E$ | total energy of the whole system | J |
| $\nabla_I$ | gradient with respect to the position of nucleus $I$ | m$^{-1}$ |
| $m_{\mathrm e}$ | electron mass, $9.11\times10^{-31}$ kg | kg |
| $M_I$ | mass of nucleus $I$ (proton $\approx 1836\,m_{\mathrm e}$) | kg |
| $\mathbf A_I, B_I$ | non-adiabatic coupling matrix elements (4.6.10) | m$^{-1}$, m$^{-2}$ |

## 4.6.1 Mass and timescale disparity

A proton is approximately 1836 times heavier than an electron. A typical nucleus (carbon, silicon, …) is 10$^4$ – 10$^5$ times heavier. From this single fact a great deal follows.

Suppose we have an electron and a nucleus, both with comparable kinetic energies — as is the case in any bound state where the virial theorem applies. Equipartition gives roughly $\tfrac12 m v^2 \sim k_B T$ or, more relevantly for chemistry, comparable contributions to the total energy. Comparing velocities at equal kinetic energy:

$$\frac{v_{\mathrm{n}}}{v_{\mathrm{e}}} = \sqrt{\frac{m_{\mathrm e}}{M}}\sim \sqrt{\frac{1}{1836}} \approx 0.023. \tag{4.6.1}$$

Nuclei move roughly a hundredth of the speed of electrons. Equivalently, the characteristic *timescale* of nuclear motion is two orders of magnitude longer than that of electronic motion. Typical molecular vibrations are at $\sim 10^{13}$–$10^{14}$ Hz (period $\sim 10$ fs); typical electronic transitions are at $\sim 10^{15}$–$10^{16}$ Hz (period $\sim 0.1$ fs).

??? note "Full derivation: why nuclear energies are $\sqrt{m_{\mathrm e}/M}$ times electronic ones"
    The single small parameter $m_{\mathrm e}/M \sim 1/1836$ controls *every* energy scale in the problem. Here is the dimensional argument Born and Oppenheimer made precise; it tells us how big vibrational and rotational energies are compared with electronic ones, and is the reason the approximation works.

    **Electronic scale.** An electron is confined by the molecule to a region of size $a$ (an atomic bond length, $\sim 1$ Å). By the uncertainty principle its momentum is at least $p_{\mathrm e}\sim \hbar/a$, so its kinetic — and hence its characteristic total — energy is

    $$E_{\mathrm e}\sim \frac{p_{\mathrm e}^2}{2m_{\mathrm e}} \sim \frac{\hbar^2}{2 m_{\mathrm e} a^2}. \tag{4.6.1a}$$

    This sets the spacing of *electronic* energy levels (a few eV for valence electrons).

    **Vibrational scale.** Now consider a nucleus vibrating in the bottom of the PES. Near the minimum the surface is harmonic, $E_0(\mathbf R)\approx E_0^{\min} + \tfrac12 M\omega_{\mathrm{vib}}^2 (\Delta R)^2$. What is the curvature $M\omega_{\mathrm{vib}}^2$? It is set by the *electronic* energy: if you displace a nucleus by the full bond length $a$ you change the molecular energy by roughly one electronic quantum $E_{\mathrm e}$ — that is what "the bond breaks" means. So the spring constant is

    $$k = M\omega_{\mathrm{vib}}^2 \sim \frac{E_{\mathrm e}}{a^2} \sim \frac{\hbar^2}{2 m_{\mathrm e} a^4}, \tag{4.6.1b}$$

    using (4.6.1a). The crucial point is that $k$ contains the *electron* mass, not the nuclear mass — the stiffness comes from the electron cloud, which is the glue. Solving for the vibrational frequency,

    $$\omega_{\mathrm{vib}} = \sqrt{\frac{k}{M}} \sim \sqrt{\frac{\hbar^2}{2 m_{\mathrm e} a^4 M}}, \tag{4.6.1c}$$

    so the vibrational quantum is

    $$E_{\mathrm{vib}} = \hbar\omega_{\mathrm{vib}} \sim \frac{\hbar^2}{2 m_{\mathrm e} a^2}\sqrt{\frac{m_{\mathrm e}}{M}} = E_{\mathrm e}\,\sqrt{\frac{m_{\mathrm e}}{M}}. \tag{4.6.1d}$$

    There it is: **a vibrational quantum is smaller than an electronic one by the factor $\sqrt{m_{\mathrm e}/M} = (m_{\mathrm e}/M)^{1/2}\approx 0.023$.** With $E_{\mathrm e}\sim$ a few eV this gives $E_{\mathrm{vib}}\sim$ a few $\times 10^{-2}$ eV $\sim 0.1$ eV $\sim$ a few hundred cm$^{-1}$ to a few thousand cm$^{-1}$, exactly the infrared range observed.

    **Rotational scale.** A rotating molecule of moment of inertia $I\sim Ma^2$ has rotational quantum $E_{\mathrm{rot}}\sim \hbar^2/(2I)\sim \hbar^2/(2Ma^2)$. Comparing with (4.6.1a),

    $$E_{\mathrm{rot}} \sim E_{\mathrm e}\,\frac{m_{\mathrm e}}{M} = E_{\mathrm e}\,(m_{\mathrm e}/M)^{1}. \tag{4.6.1e}$$

    **The hierarchy.** Writing $\kappa \equiv (m_{\mathrm e}/M)^{1/4}\approx 0.15$ (so $\kappa^2 = (m_{\mathrm e}/M)^{1/2}\approx 0.023$, $\kappa^4 = m_{\mathrm e}/M \approx 5\times10^{-4}$), the three energy scales separate cleanly as powers of $\kappa$:

    $$E_{\mathrm e} : E_{\mathrm{vib}} : E_{\mathrm{rot}} \;\sim\; 1 : \kappa^2 : \kappa^4. \tag{4.6.1f}$$

    Each is smaller than the last by two powers of $\kappa$, i.e. by roughly a factor of $0.15^2 \approx 0.02$. This is precisely the clean separation — electronic (eV, visible), vibrational ($10^{-2}$–$10^{-1}$ eV, infrared), rotational ($10^{-3}$ eV, microwave) — referred to in the historical note at the end of the section, and it is *why* the single-surface ansatz works: the energy a slow nucleus carries is far too small to bridge the electronic gap and excite the electrons.

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

!!! example "Step-by-step: substituting the ansatz and projecting onto $\psi_0$"
    The jump from (4.6.6) to (4.6.9) packs three operations into one line. Here it is broken out.

    1. **Write out the full Hamiltonian acting on the product.** With $\hat H = \hat T_{\mathrm n} + \hat H_{\mathrm e}$ and $\Psi = \chi\,\psi_0$, the Schrödinger equation $\hat H\Psi = E\Psi$ reads
       $$\hat T_{\mathrm n}\bigl[\chi\,\psi_0\bigr] + \hat H_{\mathrm e}\bigl[\chi\,\psi_0\bigr] = E\,\chi\,\psi_0.$$
    2. **The electronic term is easy.** $\hat H_{\mathrm e}$ acts only on electronic coordinates, and $\chi(\mathbf R)$ is a constant as far as it is concerned, so $\hat H_{\mathrm e}[\chi\,\psi_0] = \chi\,\hat H_{\mathrm e}\psi_0 = \chi\, E_0(\mathbf R)\,\psi_0$, using (4.6.4). This is equation (4.6.7).
    3. **The nuclear-kinetic term needs the product rule.** $\hat T_{\mathrm n}=-\sum_I \tfrac{\hbar^2}{2M_I}\nabla_I^2$ differentiates with respect to the nuclei, and *both* $\chi(\mathbf R)$ and $\psi_0(\mathbf r;\mathbf R)$ depend on $\mathbf R$. Applying $\nabla_I^2$ to a product of two $\mathbf R$-dependent factors gives three pieces — this is equation (4.6.8):
       $$\nabla_I^2(\chi\,\psi_0) = \psi_0\,\nabla_I^2\chi + 2(\nabla_I\chi)\cdot(\nabla_I\psi_0) + \chi\,\nabla_I^2\psi_0.$$
       (Recall $\nabla^2(fg)=g\nabla^2 f + 2\nabla f\cdot\nabla g + f\nabla^2 g$ — the operator version of $(fg)''=f''g+2f'g'+fg''$.)
    4. **Project onto $\psi_0$.** Multiply the whole equation from the left by $\psi_0^*(\mathbf r;\mathbf R)$ and integrate over *electronic* coordinates only, $\int(\cdots)\,d\mathbf r$. The nuclear factor $\chi$ and the operators $\nabla_I$ acting on it pass straight through this integral, because they do not touch $\mathbf r$. Term by term:
       - From $\psi_0\nabla_I^2\chi$: $\;\int\psi_0^*\psi_0\,d\mathbf r\;\nabla_I^2\chi = \nabla_I^2\chi$, since $\psi_0$ is **normalised**, $\int|\psi_0|^2\,d\mathbf r = 1$.
       - From $2(\nabla_I\chi)\cdot(\nabla_I\psi_0)$: $\;2\Bigl(\int\psi_0^*\nabla_I\psi_0\,d\mathbf r\Bigr)\cdot\nabla_I\chi = 2\,\mathbf A_I\cdot\nabla_I\chi$, defining $\mathbf A_I$ as in (4.6.10).
       - From $\chi\,\nabla_I^2\psi_0$: $\;\Bigl(\int\psi_0^*\nabla_I^2\psi_0\,d\mathbf r\Bigr)\chi = B_I\,\chi$, defining $B_I$.
    5. **Collect.** The electronic term contributes $E_0(\mathbf R)\,\chi\int|\psi_0|^2 d\mathbf r = E_0(\mathbf R)\,\chi$, and the right-hand side gives $E\,\chi$ the same way. Putting the kinetic pieces together yields exactly (4.6.9). The terms $2\mathbf A_I\cdot\nabla_I\chi$ and $B_I\chi$ — the ones that survive *only because $\psi_0$ depends on $\mathbf R$* — are the non-adiabatic couplings we are about to neglect.

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

!!! warning "Common misunderstandings"
    - **"Born–Oppenheimer is exact."** It is not. It is an *approximation*
      — we deliberately dropped the terms $\mathbf A_I$ and $B_I$ in
      (4.6.11). The full wavefunction is the Born expansion
      $\Psi=\sum_n\chi_n\psi_n$ over *all* electronic states; BO keeps one
      term. The approximation is superbly accurate when the electronic gap
      is large, but it is never identically exact. (Many later chapters
      *treat* it as exact — that is a working convenience, not a claim
      about the underlying physics.)
    - **"BO works everywhere."** It fails wherever two electronic surfaces
      come close: at level crossings, in metals (zero gap), and most
      dramatically at **conical intersections**, where the dropped
      couplings (4.6.NA) diverge. We deal with this in §4.6.6 — do not
      read the single-surface picture as universal.
    - **"Fixing the nuclei means the nuclei have no kinetic energy."** No.
      We fix the nuclei only to *solve the electronic problem* (Step 1).
      The nuclei still move — their kinetic energy reappears in the
      nuclear equation (4.6.12). "Clamped nuclei" is a calculational
      stage, not a physical claim that the nuclei are frozen.
    - **"$E_0(\mathbf R)$ depends on the nuclear masses."** It does not.
      The PES comes from the electronic Hamiltonian (4.6.3), which
      contains no nuclear masses. Swapping H for D leaves the surface
      identical; only the nuclear motion *on* it changes (hence isotope
      effects, §4.6.5).
    - **The semicolon in $\psi_0(\mathbf r; \mathbf R)$ is not decoration.**
      It means $\mathbf R$ is a *parameter*: the electronic operator never
      differentiates with respect to it. The *nuclear* operator
      $\nabla_I$, by contrast, does — which is exactly why the coupling
      terms appear.

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

??? note "Full derivation: the off-diagonal coupling and where the gap comes from"
    We want $\mathbf A_I^{(mn)}=\langle\psi_m|\nabla_I\psi_n\rangle$ for two *different* electronic states, $m\neq n$. (The diagonal case $m=n$ behaves quite differently — see the remark at the end.)

    1. **Start from the electronic eigenvalue equation** for state $n$ at every $\mathbf R$:
       $$\hat H_{\mathrm e}\,|\psi_n\rangle = E_n\,|\psi_n\rangle.$$
    2. **Differentiate both sides with respect to a nuclear coordinate** $R_I^\alpha$ (write $\nabla_I$ for short). Use the product rule on each side:
       $$(\nabla_I\hat H_{\mathrm e})|\psi_n\rangle + \hat H_{\mathrm e}|\nabla_I\psi_n\rangle = (\nabla_I E_n)|\psi_n\rangle + E_n|\nabla_I\psi_n\rangle.$$
    3. **Project onto a different eigenstate $\langle\psi_m|$ with $m\neq n$.** Take the inner product of the whole equation with $\langle\psi_m|$:
       $$\langle\psi_m|\nabla_I\hat H_{\mathrm e}|\psi_n\rangle + \langle\psi_m|\hat H_{\mathrm e}|\nabla_I\psi_n\rangle = (\nabla_I E_n)\underbrace{\langle\psi_m|\psi_n\rangle}_{=\,0} + E_n\langle\psi_m|\nabla_I\psi_n\rangle.$$
       The first right-hand term dies because eigenstates of a Hermitian operator with different labels are **orthogonal**, $\langle\psi_m|\psi_n\rangle=0$ for $m\neq n$.
    4. **Move the Hamiltonian onto the bra.** $\hat H_{\mathrm e}$ is Hermitian, so $\langle\psi_m|\hat H_{\mathrm e} = E_m\langle\psi_m|$. The second left-hand term becomes $E_m\langle\psi_m|\nabla_I\psi_n\rangle$. The equation is now
       $$\langle\psi_m|\nabla_I\hat H_{\mathrm e}|\psi_n\rangle + E_m\langle\psi_m|\nabla_I\psi_n\rangle = E_n\langle\psi_m|\nabla_I\psi_n\rangle.$$
    5. **Solve for the coupling.** Collect the $\langle\psi_m|\nabla_I\psi_n\rangle$ terms:
       $$\langle\psi_m|\nabla_I\hat H_{\mathrm e}|\psi_n\rangle = (E_n - E_m)\,\langle\psi_m|\nabla_I\psi_n\rangle,$$
       and divide by the gap $(E_n-E_m)$, which is non-zero precisely while the states are distinct:
       $$\mathbf A_I^{(mn)} = \langle\psi_m|\nabla_I\psi_n\rangle = \frac{\langle\psi_m|\nabla_I\hat H_{\mathrm e}|\psi_n\rangle}{E_n - E_m},$$
       which is (4.6.NA). The numerator is a smooth, finite quantity (the same kind of object that gave us forces in the Hellmann–Feynman theorem); all the danger sits in the denominator.

    **Why the diagonal case is different.** Setting $m=n$ in step 3 would put $\langle\psi_n|\psi_n\rangle=1$, not $0$, and the $(E_n-E_m)$ factor in step 5 would vanish — division by zero. So this derivation says *nothing* about $\mathbf A_I^{(nn)}$; that diagonal term is instead fixed (to a pure phase) by differentiating the normalisation $\langle\psi_n|\psi_n\rangle=1$, exactly as in the Hellmann–Feynman proof above, and can be made zero by a suitable choice of phase for $\psi_n$.

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

!!! question "Check yourself"
    1. In the ansatz $\Psi(\mathbf r,\mathbf R)\approx\chi(\mathbf R)\,\psi_0(\mathbf r;\mathbf R)$, which factor describes the electrons and which the nuclei? What does the semicolon in $\psi_0(\mathbf r;\mathbf R)$ tell you?
    2. The full nuclear equation (4.6.9) has three terms inside the brackets: $\nabla_I^2\chi$, $2\mathbf A_I\cdot\nabla_I\chi$ and $B_I\chi$. Which two does the Born–Oppenheimer approximation drop, and what physical situation makes them small?
    3. Using $E_{\mathrm e}\sim 5$ eV and the scaling $E_{\mathrm{vib}}\sim E_{\mathrm e}\sqrt{m_{\mathrm e}/M}$, estimate the size of a vibrational quantum (take $\sqrt{m_{\mathrm e}/M}\approx 0.023$). Does your answer sit in the infrared, as claimed?
    4. What is a potential energy surface, and why is it the *same* curve for an H$_2$ molecule and a D$_2$ (deuterium) molecule even though their vibrational frequencies differ?
    5. Name the two ways the BO approximation breaks down, and state in one phrase what the two electronic states do to each other in each case.

    ??? note "Hint"
        For 2, look at (4.6.10) and (4.6.11): the dropped terms are the ones containing *derivatives of the electronic wavefunction* with respect to the nuclei. For 3, just multiply. For 4, ask which Hamiltonian — electronic or nuclear — the masses appear in. For 5, re-read §4.6.6: one mode involves a small but finite gap, the other a gap that closes to zero.

    ??? success "Answer"
        1. $\chi(\mathbf R)$ is the **nuclear** wavefunction; $\psi_0(\mathbf r;\mathbf R)$ is the **electronic** ground state. The semicolon marks $\mathbf R$ as a *parameter*: the electronic Hamiltonian depends on the nuclear positions but never differentiates with respect to them — the nuclei are clamped while we solve for the electrons.
        2. BO drops $2\mathbf A_I\cdot\nabla_I\chi$ and $B_I\chi$ (the non-adiabatic couplings of (4.6.10)), keeping only $\nabla_I^2\chi$. They are small when the electronic ground state is **well separated in energy** from the excited states (a large gap) and the nuclei move slowly, so the electronic wavefunction changes only gently with $\mathbf R$. Their size relative to the kept term is of order $m_{\mathrm e}/M$.
        3. $E_{\mathrm{vib}}\sim 5\ \text{eV}\times 0.023 \approx 0.12\ \text{eV}$. Converting, $0.12\ \text{eV}\approx 0.12/(1.24\times10^{-4})\ \text{cm}^{-1}\approx 9\times10^{2}\ \text{cm}^{-1}$ — a few hundred to a thousand cm$^{-1}$, squarely in the infrared. Yes.
        4. A PES is the electronic ground-state energy as a function of nuclear positions, $E_0(\mathbf R)$. It is produced by the **electronic** Hamiltonian (4.6.3), which contains *no nuclear masses*. Replacing H by D changes only the nuclear mass $M$, which enters the *nuclear* equation (4.6.12), so the surface is unchanged but the vibrational levels — set by $\omega=\sqrt{V''/\mu}$ with the heavier reduced mass $\mu$ — shift down. This is the origin of isotope effects.
        5. **Non-adiabatic coupling**: two states with a small-but-finite gap mix, so nuclear motion can kick the system from $\psi_0$ into $\psi_1$. **Conical intersection**: two surfaces actually *touch* (gap $\to 0$), the coupling (4.6.NA) diverges, and the electronic wavefunction becomes undefined — single-surface BO fails completely.
