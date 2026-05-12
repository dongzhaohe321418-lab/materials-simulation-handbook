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

The Bohr model is best regarded as a brilliant interim measure. Something deeper was needed, and the next idea came from an unlikely source.

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
