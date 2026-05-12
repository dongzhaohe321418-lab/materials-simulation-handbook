# 4.1 Why we need quantum mechanics

By the end of the nineteenth century physics looked, to its practitioners, very nearly finished. Newton's laws governed mechanics. Maxwell's equations described the electromagnetic field. Boltzmann had reduced thermodynamics to the statistical motion of atoms. It was widely believed that within a generation every observed phenomenon would be expressible as a corollary of these three pillars. Lord Kelvin, addressing the Royal Institution in 1900, famously remarked that there were only "two small clouds" remaining on the horizon — the failure of the aether-drift experiments, and the puzzle of blackbody radiation. Within twenty-five years those two clouds had grown into the two great storms of modern physics: relativity and quantum mechanics.

This section tells the story of the second storm. We will not develop the theory yet — that is the work of §4.2. Our aim here is more modest, and more important: to make the reader *feel* why a wave description of matter is not optional. Classical mechanics does not just give slightly wrong answers at small scales; it gives qualitatively wrong answers, and in some cases predicts catastrophes that simply do not occur. The arguments below are the same ones that convinced Planck, Einstein, Bohr and de Broglie — they are worth working through carefully even though every physicist alive accepts the conclusion.

## 4.1.1 The ultraviolet catastrophe

Heat any object — a poker, a star, the tungsten filament of a light bulb — and it glows. The spectrum of the emitted light depends only on the temperature: a "blackbody", an idealised perfect absorber and emitter, radiates with a universal spectral curve $u(\nu, T)$ that gives the energy per unit volume per unit frequency interval.

Classical electromagnetism makes a very definite prediction for $u(\nu, T)$. The argument, due to Rayleigh and refined by Jeans, runs as follows. Consider a cubical cavity of side $L$ in thermal equilibrium. The electromagnetic field inside can be decomposed into standing-wave modes, each labelled by a wavevector $\mathbf k$. The number of modes per unit volume with frequency between $\nu$ and $\nu + d\nu$ is

$$g(\nu)\, d\nu = \frac{8\pi \nu^2}{c^3}\, d\nu.$$

By the classical equipartition theorem, each mode is an independent harmonic oscillator and carries an average energy $k_{\mathrm B}T$ in thermal equilibrium. Multiplying gives the Rayleigh–Jeans law,

$$u_{\mathrm{RJ}}(\nu, T) = \frac{8\pi \nu^2}{c^3}\, k_{\mathrm B} T. \tag{4.1.1}$$

At low frequencies — radio waves, microwaves, the red end of the visible spectrum — this formula matches experiment beautifully. But it has a fatal feature: it grows without bound as $\nu \to \infty$. The total radiated energy per unit volume,

$$U = \int_0^\infty u_{\mathrm{RJ}}(\nu, T)\, d\nu = \infty,$$

diverges. This is the "ultraviolet catastrophe": classical physics predicts that *any* warm object should emit an infinite amount of high-frequency radiation. A glowing coal should incinerate the room. It does not.

In December 1900 Max Planck produced a fix. He postulated — at first as a purely mathematical trick — that the energy of each electromagnetic mode of frequency $\nu$ is not continuously variable but comes in discrete lumps of size

$$\varepsilon = h\nu, \tag{4.1.2}$$

where $h \approx 6.626 \times 10^{-34}$ J s. Repeating the equipartition calculation with quantised energies replaces $k_{\mathrm B} T$ by $h\nu /(e^{h\nu/k_{\mathrm B}T} - 1)$, and the resulting Planck law,

$$u(\nu, T) = \frac{8\pi h \nu^3}{c^3}\, \frac{1}{e^{h\nu / k_{\mathrm B} T} - 1}, \tag{4.1.3}$$

reduces to Rayleigh–Jeans at low frequency and cuts off exponentially at high frequency. The total energy is finite and agrees with the Stefan–Boltzmann law. The data fit is perfect.

!!! note "Why this matters for materials"
    Equation (4.1.2) is the first crack in classical physics. Whatever else is true, energy at the atomic scale comes in *discrete amounts*. This is the seed of every band gap, every vibrational quantum, every laser line. When we compute the electronic levels of a solid in Chapter 5 we are computing the descendants of Planck's quanta.

Planck himself was uncomfortable with his own postulate and spent years trying to derive it from classical physics. He failed, because there is no such derivation: nature is not classical.

!!! tip "Reading the Planck spectrum"
    There are two useful limiting forms of (4.1.3). For $h\nu \ll k_{\mathrm B}T$ we may expand the exponential as $e^{h\nu/k_{\mathrm B}T} \approx 1 + h\nu/k_{\mathrm B}T$ and recover Rayleigh–Jeans:
    $$u(\nu, T) \approx \frac{8\pi\nu^2}{c^3}\, k_{\mathrm B}T, \qquad (h\nu \ll k_{\mathrm B}T).$$
    For $h\nu \gg k_{\mathrm B}T$ the denominator $e^{h\nu/k_{\mathrm B}T} - 1 \approx e^{h\nu/k_{\mathrm B}T}$ and the spectrum decays exponentially: this is **Wien's law**, $u(\nu,T) \propto \nu^3 e^{-h\nu/k_{\mathrm B}T}$. The crossover happens at $h\nu \sim k_{\mathrm B}T$, which for $T = 300$ K is $\nu \sim 6\times 10^{12}$ Hz — in the far infrared. This is why room-temperature objects glow in the IR, not the visible: their thermal energy $k_{\mathrm B}T \approx 0.026$ eV is well below visible photon energies $\sim 2$ eV.

The numerical check is instructive. Differentiating $u(\nu, T)$ with respect to $\nu$ and setting the derivative to zero yields **Wien's displacement law**, $h\nu_{\max} \approx 2.82\, k_{\mathrm B}T$. For the surface of the Sun ($T \approx 5800$ K), $\nu_{\max} \approx 3.4\times 10^{14}$ Hz, i.e.\ $\lambda_{\max} \approx 880$ nm — in the near-infrared, with the peak energy density spilling into the visible. Human vision and the solar spectrum are matched, and the match is set by a single equation that classical physics could not have produced.

## 4.1.2 The photoelectric effect

The next crack appeared in 1905, in one of the four miraculous papers Einstein published that year. Shine light onto a clean metal surface and electrons are ejected — the photoelectric effect, discovered by Hertz in 1887. Classical electromagnetism makes three predictions:

1. The kinetic energy of the ejected electrons should depend on the *intensity* of the light (brighter light, more energy per electron).
2. There should be a measurable delay between switching the light on and observing the first electron — time enough for the wave to deposit sufficient energy.
3. Electrons should be ejected at any frequency, provided one waits long enough.

Every one of these predictions is wrong. Experiment shows:

1. The maximum kinetic energy depends on the *frequency*, not the intensity. Doubling the intensity doubles the number of electrons emitted per second but leaves their energies unchanged.
2. There is no detectable delay, even with light so dim that classical estimates would require minutes of accumulation.
3. Below a threshold frequency $\nu_0$, no electrons emerge at all, no matter how bright the light.

Einstein's resolution was breathtakingly simple. Take Planck's quanta literally. Light of frequency $\nu$ consists of particles — *photons* — each carrying energy $h\nu$. An electron in the metal is bound by a "work function" $\phi$; absorbing a single photon either liberates it (if $h\nu > \phi$, with leftover kinetic energy $h\nu - \phi$) or does not (if $h\nu < \phi$):

$$E_{\mathrm{kin}}^{\mathrm{max}} = h\nu - \phi. \tag{4.1.4}$$

This explains the threshold (the existence of $\nu_0 = \phi/h$), the linearity in $\nu$, the independence from intensity, and the absence of delay all at once.

Einstein received the 1921 Nobel Prize for this work — *not* for relativity, which the Swedish Academy still considered too speculative. Equation (4.1.4) is the moment light became dual: a wave (interference, diffraction, Maxwell's equations) *and* a particle (photoelectric effect, Compton scattering).

!!! note "Why this step? — Why classical waves cannot explain the threshold"
    In a classical wave the energy flux is $\langle S\rangle = \tfrac12 c\varepsilon_0 E_0^2$, proportional to the *intensity*, not the frequency. The energy is spread continuously over the wavefront, and a bound electron should be able to accumulate $\phi$ worth of energy from a sufficiently long exposure to *any* frequency. Concretely: a $1$ mW HeNe laser ($\lambda = 633$ nm) striking a $1$ cm$^2$ caesium photocathode delivers about $10^{15}$ photons per second over an area covering $\sim 10^{15}$ surface atoms. Classically every atom receives the same trickle of energy and would need many seconds to accumulate $\phi \approx 2$ eV. Experimentally, photoemission begins within a femtosecond and the kinetic energy depends *only* on $\nu$. Einstein's photon hypothesis is the only resolution: each photon either has enough energy individually or nothing happens.

!!! example "Numerical sanity check on the photoelectric equation"
    Caesium has a work function $\phi_{\mathrm{Cs}} \approx 2.1$ eV. Illuminating it with 400 nm violet light ($h\nu = 1240/400 \approx 3.10$ eV) gives $E_{\mathrm{kin}}^{\max} = 3.10 - 2.10 = 1.0$ eV. Doubling the intensity doubles the *count* of photoelectrons but leaves their kinetic energy unchanged. Switching to 700 nm red light ($h\nu \approx 1.77$ eV $< \phi$) produces no photoelectrons at all, no matter how bright the beam. Both predictions are textbook experiments, and both are direct consequences of (4.1.4).

The deeper significance of Einstein's photon is not the photoelectric effect itself — that is one experiment — but the fact that *energy quantisation propagates*. Planck quantised the oscillators in the cavity walls; Einstein quantised the field that travels between them. Once light itself comes in lumps $h\nu$, the symmetry suggesting "matter is also wavelike" is hard to avoid. We will see de Broglie complete that symmetry in §4.1.4.

## 4.1.3 The atom should not exist

The third — and for materials physics most acute — crisis concerns the stability of atoms.

By 1911 Ernest Rutherford's gold-foil experiments had established that an atom consists of a tiny dense positive nucleus surrounded by negative electrons. The natural classical picture is a miniature solar system: electrons orbit the nucleus under the Coulomb attraction, much as planets orbit the sun under gravity. The maths is identical: closed elliptical orbits with energies given by Kepler's laws.

This picture is catastrophically wrong. An electron in circular orbit is an accelerating charge, and an accelerating charge radiates electromagnetic waves — this is exactly how a radio antenna works. The Larmor formula gives the radiated power:

$$P = \frac{e^2 a^2}{6\pi \varepsilon_0 c^3}, \tag{4.1.5}$$

where $a$ is the acceleration. Plug in the numbers for a hydrogen atom: an electron at the Bohr radius $a_0 \approx 0.529$ Å, circling at the velocity required to balance the Coulomb force, has centripetal acceleration $a \sim 9 \times 10^{22}$ m s$^{-2}$. The radiated power is enormous. As it loses energy the electron spirals inward, accelerating still more and radiating still faster. A straightforward integration gives the lifetime of a classical hydrogen atom:

$$\tau \sim 10^{-11}\ \mathrm{s}. \tag{4.1.6}$$

This is not a minor discrepancy. Classical physics predicts that every atom in your body should collapse, releasing a flash of ultraviolet light, in less than a tenth of a nanosecond. Yet the hydrogen atom has been observed, undisturbed, in interstellar clouds for thirteen billion years.

Even worse, the predicted radiation should be a smooth continuum of frequencies — yet what we observe is a sharp line spectrum. Hydrogen emits at very specific wavelengths (the Balmer series in the visible: 656.3 nm, 486.1 nm, 434.0 nm, 410.2 nm, …) and at no others. Empirically, Rydberg had fitted these to

$$\frac{1}{\lambda} = R\left(\frac{1}{n_1^2} - \frac{1}{n_2^2}\right), \quad n_1 < n_2, \tag{4.1.7}$$

with $R = 1.097 \times 10^7$ m$^{-1}$, but no classical mechanism produced anything resembling integer-labelled spectra.

Niels Bohr's 1913 model patched the problem by *fiat*: postulate that the electron is allowed only on certain orbits with quantised angular momentum $L = n\hbar$ (where $\hbar = h/2\pi$), and assert that no radiation is emitted on these orbits, only during jumps between them. The model reproduced the hydrogen spectrum to remarkable accuracy and even predicted the Rydberg constant in terms of fundamental quantities. But it was, transparently, a kludge — it explained nothing about *why* angular momentum should be quantised, and it failed catastrophically for helium.

### The Bohr atom in three lines of algebra

It is worth deriving the Bohr formulae explicitly: even though the model is wrong in detail, the energy scale and length scale it produces are quantitatively right, and they pervade quantum mechanics.

Postulate 1 (mechanical balance): an electron of mass $m_{\mathrm e}$ in a circular orbit of radius $r$ around a proton experiences a Coulomb force balanced by centripetal acceleration,

$$\frac{e^2}{4\pi\varepsilon_0\, r^2} = \frac{m_{\mathrm e} v^2}{r}. \tag{4.1.A}$$

!!! note "Why this step?"
    Both sides are classical: $F = m a$ with the Coulomb force. This is the *only* classical input. Quantisation enters next.

Postulate 2 (Bohr's leap): angular momentum is quantised,

$$L = m_{\mathrm e} v r = n\hbar, \quad n = 1, 2, 3, \ldots \tag{4.1.B}$$

From (4.1.B), $v = n\hbar/(m_{\mathrm e} r)$. Substitute into (4.1.A):

$$\frac{e^2}{4\pi\varepsilon_0\, r^2} = \frac{m_{\mathrm e}}{r}\cdot \frac{n^2\hbar^2}{m_{\mathrm e}^2 r^2} = \frac{n^2 \hbar^2}{m_{\mathrm e} r^3},$$

which rearranges to

$$r_n = \frac{4\pi\varepsilon_0 n^2 \hbar^2}{m_{\mathrm e} e^2} = n^2 a_0, \qquad a_0 \equiv \frac{4\pi\varepsilon_0 \hbar^2}{m_{\mathrm e} e^2} \approx 0.529\ \text{\AA}. \tag{4.1.C}$$

The **Bohr radius** $a_0$ is the natural length scale of the hydrogen atom; it remains the atomic unit of length in modern computational chemistry. The energy at this radius follows from $E = \tfrac12 m_{\mathrm e} v^2 - e^2/(4\pi\varepsilon_0 r)$. Using $\tfrac12 m_{\mathrm e}v^2 = e^2/(8\pi\varepsilon_0 r)$ (half the magnitude of the potential, by the virial theorem), the total energy is

$$E_n = -\frac{e^2}{8\pi\varepsilon_0 r_n} = -\frac{m_{\mathrm e} e^4}{2(4\pi\varepsilon_0)^2 \hbar^2}\cdot \frac{1}{n^2} = -\frac{13.6\ \text{eV}}{n^2}. \tag{4.1.D}$$

!!! example "Why this matches experiment"
    The Rydberg formula (4.1.7) is reproduced by setting $h c/\lambda = E_{n_2} - E_{n_1}$ between two Bohr levels. The Rydberg constant $R$ derived this way is $R = m_{\mathrm e} e^4 / (8\varepsilon_0^2 h^3 c) = 1.0974\times 10^7$ m$^{-1}$, in agreement with the empirical value to four significant figures. Bohr's heuristic is wrong in detail (real hydrogen orbitals are $s, p, d, \ldots$, not classical circles), but it gets the energy levels exactly right. This *cannot* be a coincidence, and de Broglie's standing-wave reinterpretation in §4.1.4 will explain why.

The Bohr model is best regarded as a brilliant interim measure: it grants energy quantisation but stops short of a wave description. Something deeper was needed, and the next idea came from an unlikely source.

## 4.1.4 The de Broglie hypothesis

In 1924 a French aristocrat-turned-physicist, Louis de Broglie, submitted a doctoral thesis containing one of the boldest leaps in twentieth-century physics. Light, Einstein had shown, is both a wave (wavelength $\lambda$) and a particle (momentum $p = h\nu/c = h/\lambda$). What, asked de Broglie, if matter is too?

He posited that any particle of momentum $p$ has an associated wavelength,

$$\lambda = \frac{h}{p}, \tag{4.1.8}$$

the *de Broglie wavelength*. For a tennis ball ($p \sim 1$ kg m s$^{-1}$) this gives $\lambda \sim 10^{-34}$ m, far below any conceivable measurement — which is why we never see tennis balls diffract. For an electron in an atom ($p \sim 10^{-24}$ kg m s$^{-1}$) it gives $\lambda \sim 10^{-10}$ m, comparable to the atom itself.

The Bohr quantisation condition $L = n\hbar$ now acquires a transparent interpretation: a stable orbit is one whose circumference accommodates an integer number of de Broglie wavelengths,

$$2\pi r = n\lambda \quad \Longleftrightarrow \quad pr = n\hbar \quad \Longleftrightarrow \quad L = n\hbar.$$

The "allowed orbits" are *standing waves*. The electron does not orbit at all in the planetary sense — it is a wave wrapped around the nucleus.

Three years later, Davisson and Germer at Bell Labs accidentally scattered electrons off a nickel crystal and recorded the diffraction pattern. Electrons, indubitable particles, were producing the same Bragg peaks that X-rays produce. The de Broglie hypothesis was confirmed. By the 1990s the same experiment had been done with whole molecules of C$_{60}$.

!!! example "Numerical check"
    A 100 eV electron has momentum $p = \sqrt{2m_e E} \approx 5.4 \times 10^{-24}$ kg m s$^{-1}$, giving $\lambda \approx 1.2$ Å. Electron diffraction is therefore a routine technique for probing crystal structure. A 1 eV thermal neutron has $\lambda \approx 0.28$ Å, which is why neutron diffraction works.

### Wave–particle duality and the double slit

The most celebrated experimental embodiment of de Broglie's hypothesis is the **double-slit experiment**. A coherent beam of particles — photons, electrons, neutrons, atoms, or even C$_{60}$ molecules — is directed at a screen pierced by two narrow slits a distance $d$ apart, and the arrivals are recorded on a detector at distance $L \gg d$ behind the slits.

If the particles were classical Newtonian objects, we would expect two blurred images of the slits on the detector — one peak behind each slit. What is observed instead is an interference pattern: alternating bright and dark fringes, with peaks at positions

$$y_n = n\, \frac{\lambda L}{d}, \qquad n = 0, \pm 1, \pm 2, \ldots,$$

exactly the pattern produced by waves of wavelength $\lambda$. For light this is unsurprising and was already understood in Young's day. For *electrons* it was astonishing. And the astonishment deepens once one reduces the beam intensity to the point where only one electron is in the apparatus at any time: the individual arrivals are still detected as point particles (a single dot of phosphorescence on the screen), but the *distribution* of dots accumulated over many runs is the same interference pattern. Each electron seems to know about both slits, even though it is "indivisible".

The Davisson–Germer experiment of 1927 was a controlled accident. Davisson and Germer at Bell Labs were studying low-energy electron scattering from polycrystalline nickel when a vacuum leak forced them to anneal the sample, accidentally producing a nearly single-crystal surface. On resuming the experiment they observed sharp angular peaks in the scattered electron intensity at angles satisfying the Bragg condition $n\lambda = d\sin\theta$, with $\lambda$ given by de Broglie's formula. Electrons were diffracting from the nickel lattice exactly as X-rays would. Independently and almost simultaneously, G.\ P.\ Thomson — son of J.\ J.\ Thomson, who had identified the electron as a particle in 1897 — observed transmission electron diffraction through thin foils. Father and son shared Nobel prizes for proving, respectively, that the electron is a particle and that it is a wave.

By the 1990s the same experiment had been performed on whole C$_{60}$ fullerene molecules in Anton Zeilinger's laboratory: 60 atoms, mass 720 amu, diffracting through a transmission grating with the predicted wavelength. As of the 2020s, molecules of $\sim 25\,000$ amu have been diffracted. There is no known mass scale at which wave behaviour switches off.

### The correspondence principle

If matter is wavelike, why does classical mechanics work so well for tennis balls? The answer is the **correspondence principle**, articulated by Bohr: quantum mechanics must reduce to classical mechanics in the limit where actions are large compared to $\hbar$.

Quantitatively, the relevant small parameter is $\hbar/S$, where $S$ is a characteristic classical action of the system (with dimensions of energy × time, or momentum × length). For a tennis ball of mass $0.06$ kg moving at $20$ m s$^{-1}$ over a court of length $24$ m, the action is $S \sim m v L \sim 30$ J s, while $\hbar \approx 10^{-34}$ J s. The ratio is $10^{-35}$. Wave effects are present in principle but invisible in practice: the de Broglie wavelength is $\lambda = h/(mv) \sim 5\times 10^{-34}$ m, twenty orders of magnitude below the diameter of an atomic nucleus.

For an electron in a hydrogen atom, in contrast, $S \sim p\, a_0 \sim \hbar$ (this is in fact one definition of the atomic scale), and quantum effects are unavoidable. The crossover between the two regimes — the realm where the symbols of classical mechanics begin to lose their meaning and the symbols of quantum mechanics take over — is exactly where modern condensed-matter physics lives.

Two operational forms of the correspondence principle will recur.

- **Large quantum numbers.** For energy levels labelled by an integer $n$, the spacing $E_{n+1} - E_n$ relative to $E_n$ shrinks as $n \to \infty$, and the quantum probability density (averaged over a few wavelengths) approaches the classical position distribution $\rho_{\mathrm{cl}}(x) \propto 1/v_{\mathrm{cl}}(x)$. We will see this explicitly for the harmonic oscillator in §4.4.
- **The $\hbar \to 0$ limit.** Formally taking $\hbar \to 0$ in the Schrödinger equation recovers the Hamilton–Jacobi equation of classical mechanics (the WKB construction). Operationally it is the limit in which the wavefunction becomes a narrow wavepacket whose centroid obeys Newton's laws (Ehrenfest's theorem).

The correspondence principle answers the embarrassing question of why physics is taught in the order it is: because the limit in which classical mechanics is *valid* is the limit in which it was originally *discovered*. Below the atomic scale the limit fails and we must work with the full quantum theory. Above it the quantum theory still applies, but we may legitimately approximate it by classical mechanics — and in doing so we recover the entire edifice of macroscopic physics from a single $\hbar \to 0$ argument.

!!! tip "A unifying picture"
    The four "crises" of §4.1.1–4.1.4 are not independent: they are four faces of the same underlying fact, that action is quantised in units of $\hbar$. Planck's $h\nu$ is the action quantum for an oscillation cycle. Einstein's photon is the action quantum delivered ballistically. Bohr's angular-momentum condition $L = n\hbar$ is the action quantum on a circular orbit. De Broglie's $\lambda = h/p$ is the conversion factor from momentum to wavelength such that $\oint p\, dx = n h$ over a closed orbit. The Schrödinger equation of §4.2 will be the linear PDE that subsumes all four; once we have it, every result above will follow as a worked example.

## 4.1.5 The lesson for materials

Take stock of what these four observations together imply.

- Energy is quantised (Planck).
- Light is both wave and particle (Einstein).
- Electrons in atoms cannot orbit classically (Bohr).
- Matter, like light, has wave character (de Broglie).

For a materials physicist the consequences are stark. Electrons are the glue of chemistry: they form bonds, fill bands, determine optical properties, mediate magnetism. If electrons are waves, then *every* material property is in the end the property of a complicated multi-electron wavefunction defined on the atoms of the solid. There is no escape into classical mechanics.

What we need is an equation. We need an analogue of Newton's $\mathbf F = m\mathbf a$ that governs the evolution of these matter waves. We need a rule that tells us, given a potential $V(\mathbf r)$, what the standing-wave patterns are, what their energies are, and how a wavepacket propagates in time.

That equation was written down by Erwin Schrödinger in late 1925, working at a guesthouse in the Swiss Alps. It is the subject of the next section. We will not derive it — nobody can derive it, in the same way that nobody can derive Newton's laws. It is a postulate, justified by the spectacular agreement of its predictions with experiment over the past century.

But we can prepare ourselves to receive it. Three concepts are essential.

1. **The state of a quantum particle is a complex-valued function** $\psi(\mathbf r, t)$, the *wavefunction*. This is the analogue of the (position, momentum) pair in classical mechanics, but it carries far more information.

2. **The wavefunction obeys a linear partial differential equation.** Linearity is crucial: it means quantum states *superpose*, the way waves on water superpose. It is what makes interference possible.

3. **The connection to experiment is statistical.** $|\psi(\mathbf r, t)|^2$ is the probability density of finding the particle at $\mathbf r$ at time $t$. Quantum mechanics does not predict where a single electron will land on a detector; it predicts the *distribution* of landings over many runs.

With those three points in hand, the equation itself will look almost inevitable. Turn the page.

!!! tip "What to take from this section"
    The historical detour through Planck, Einstein, Bohr and de Broglie is not for nostalgia. Each of the four episodes identified a *specific* failure of classical physics — UV catastrophe, photoelectric effect, atomic stability, electron diffraction — and the *specific* fix needed. In §4.2 we will write down a single equation that contains all four fixes simultaneously. Knowing the historical motivation is not a luxury; it is what makes the equation legible. The same will be true in Chapter 5, where the failure modes of HF (correlation, self-interaction, dispersion) motivate the structure of DFT functionals. Throughout the book: read the failure first, then the fix.
