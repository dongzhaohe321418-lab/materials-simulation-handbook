# 4.1 Why we need quantum mechanics

By the end of the nineteenth century physics looked, to its practitioners, very nearly finished. Newton's laws governed mechanics. Maxwell's equations described the electromagnetic field. Boltzmann had reduced thermodynamics to the statistical motion of atoms. It was widely believed that within a generation every observed phenomenon would be expressible as a corollary of these three pillars. Lord Kelvin, addressing the Royal Institution in 1900, famously remarked that there were only "two small clouds" remaining on the horizon — the failure of the aether-drift experiments, and the puzzle of blackbody radiation. Within twenty-five years those two clouds had grown into the two great storms of modern physics: relativity and quantum mechanics.

This section tells the story of the second storm. We will not develop the theory yet — that is the work of §4.2. Our aim here is more modest, and more important: to make the reader *feel* why a wave description of matter is not optional. Classical mechanics does not just give slightly wrong answers at small scales; it gives qualitatively wrong answers, and in some cases predicts catastrophes that simply do not occur. The arguments below are the same ones that convinced Planck, Einstein, Bohr and de Broglie — they are worth working through carefully even though every physicist alive accepts the conclusion.

!!! info "What problem are we solving?"
    Before any quantum formalism, we need a reason to abandon the physics
    that already works for everything we can see and touch. This page is
    that reason. It collects four nineteenth- and early-twentieth-century
    experiments that classical physics gets *qualitatively* wrong — not by
    a few per cent, but by predicting infinities, instant collapses, or
    effects that simply never happen. The single thread running through all
    four is that **a wave description of matter and light is forced on us**;
    it is not a convenience we adopt for elegance. We are not deriving
    quantum mechanics here (that is §4.2) and we are not proving the wave
    equation — we are establishing that *something* like it must exist.
    Read the failure first; the fix comes in the next section.

!!! note "Plain-language version"
    Classical physics pictures light as a continuous wave and an electron
    as a tiny ball on a track. Four experiments break both pictures at
    once. A hot object would glow with infinite brightness (it does not).
    Dim light should take minutes to kick out an electron (it takes
    femtoseconds). A classical atom should collapse in a fraction of a
    nanosecond (atoms are eternal). And electrons fired at two slits build
    up a wave interference pattern (balls cannot). The only way to fit all
    four is to let energy come in lumps and let matter behave like a wave.

## 4.1.1 The ultraviolet catastrophe

Heat any object — a poker, a star, the tungsten filament of a light bulb — and it glows. The spectrum of the emitted light depends only on the temperature: a "blackbody", an idealised perfect absorber and emitter, radiates with a universal spectral curve $u(\nu, T)$ that gives the energy per unit volume per unit frequency interval.

!!! tip "New vocabulary"
    - **Blackbody** — an idealised object that absorbs all light falling on
      it and re-emits a spectrum that depends only on its temperature.
    - **Spectral energy density** $u(\nu,T)$ — energy stored in the
      radiation field, *per unit volume*, *per unit frequency interval*.
      See the symbol guide below for its units.
    - **Mode** — one standing-wave pattern the cavity can support, like one
      note a guitar string can sound. Each mode is an independent oscillator.
    - **Equipartition** — the classical rule that every independent
      quadratic energy store ("degree of freedom") holds an average energy
      $\tfrac12 k_{\mathrm B}T$ in thermal equilibrium. See the
      [beginner glossary](../undergraduate/glossary-for-beginners.md) for
      *ensemble* and temperature.

Before the formula, here is every symbol that appears in this sub-section, with its SI units.

| Symbol | Meaning | Units |
|---|---|---|
| $u(\nu, T)$ | spectral energy density of the radiation field | J m$^{-3}$ Hz$^{-1}$ |
| $\nu$ | frequency of a mode | Hz $=$ s$^{-1}$ |
| $T$ | absolute temperature | K |
| $c$ | speed of light | m s$^{-1}$ |
| $L$ | side length of the cubical cavity | m |
| $\mathbf k$ | wavevector of a standing-wave mode | m$^{-1}$ |
| $g(\nu)\,d\nu$ | number of modes per unit volume in $[\nu,\nu+d\nu]$ | m$^{-3}$ |
| $h$ | Planck constant, $6.626\times10^{-34}$ | J s |
| $k_{\mathrm B}$ | Boltzmann constant, $1.381\times10^{-23}$ | J K$^{-1}$ |
| $\varepsilon = h\nu$ | energy quantum of a mode | J |

!!! note "Plain-language version"
    The Rayleigh–Jeans argument has two halves. First, *count the modes*:
    how many independent standing waves of each frequency fit inside a box?
    Second, *give each mode its share of thermal energy* using
    equipartition. Multiplying the two gives an energy density. The trouble
    is that the number of modes grows like $\nu^2$ without limit, so the
    total energy diverges. Planck fixes the *second* half — how much energy
    a mode actually carries — not the counting.

Classical electromagnetism makes a very definite prediction for $u(\nu, T)$. The argument, due to Rayleigh and refined by Jeans, runs as follows. Consider a cubical cavity of side $L$ in thermal equilibrium. The electromagnetic field inside can be decomposed into standing-wave modes, each labelled by a wavevector $\mathbf k$. The number of modes per unit volume with frequency between $\nu$ and $\nu + d\nu$ is

$$g(\nu)\, d\nu = \frac{8\pi \nu^2}{c^3}\, d\nu.$$

The factor $8\pi\nu^2/c^3$ is not pulled from the air; it is a careful count of standing waves, worked through line by line below.

??? note "Full derivation: the mode density $g(\nu) = 8\pi\nu^2/c^3$"
    **Step 1 — standing waves in a box.** A field component that vanishes on
    the walls of a cube of side $L$ must be a standing wave, so along each
    axis it is a sine with a whole number of half-wavelengths fitting in the
    box. The allowed wavevector components are therefore

    $$k_x = \frac{n_x \pi}{L},\quad k_y = \frac{n_y \pi}{L},\quad k_z = \frac{n_z \pi}{L},\qquad n_x, n_y, n_z = 1, 2, 3, \ldots \tag{4.1.1a}$$

    Each ordered triple $(n_x, n_y, n_z)$ of *positive* integers labels one
    spatial mode.

    **Step 2 — modes as a grid of points.** Plot the allowed $\mathbf k$ in
    "$k$-space". They form a regular cubic grid with spacing $\pi/L$ along
    each axis, so each mode occupies a small cube of volume $(\pi/L)^3$. The
    number of modes is therefore the available $k$-space volume divided by
    $(\pi/L)^3$.

    **Step 3 — count modes with $|\mathbf k|$ up to $k$.** Because
    $n_x, n_y, n_z$ are all positive, the allowed points fill only the
    *positive octant* — one eighth — of a sphere of radius $k$ in $k$-space.
    The volume of that octant is $\tfrac18 \cdot \tfrac43 \pi k^3$, so the
    number of spatial modes with wavenumber up to $k$ is

    $$N(k) = \frac{\tfrac18 \cdot \tfrac43 \pi k^3}{(\pi/L)^3} = \frac{\tfrac16 \pi k^3 \, L^3}{\pi^3} = \frac{L^3 k^3}{6\pi^2}. \tag{4.1.1b}$$

    **Step 4 — two polarisations.** An electromagnetic wave has two
    independent transverse polarisation states for each spatial mode.
    Multiply by 2:

    $$N(k) = 2 \cdot \frac{L^3 k^3}{6\pi^2} = \frac{L^3 k^3}{3\pi^2}. \tag{4.1.1c}$$

    **Step 5 — change variable from $k$ to $\nu$.** For light in vacuum
    $\nu = ck/2\pi$, i.e. $k = 2\pi\nu/c$, so $k^3 = 8\pi^3\nu^3/c^3$.
    Substituting,

    $$N(\nu) = \frac{L^3}{3\pi^2}\cdot \frac{8\pi^3 \nu^3}{c^3} = \frac{8\pi L^3 \nu^3}{3 c^3}. \tag{4.1.1d}$$

    **Step 6 — differentiate and divide by volume.** The number of modes in
    a thin frequency shell $[\nu, \nu+d\nu]$ is $dN = (dN/d\nu)\,d\nu$, and
    the mode density per unit volume is $g(\nu)\,d\nu = dN/L^3$:

    $$g(\nu) = \frac{1}{L^3}\frac{dN}{d\nu} = \frac{1}{L^3}\cdot \frac{8\pi L^3}{3c^3}\cdot 3\nu^2 = \frac{8\pi \nu^2}{c^3}. \tag{4.1.1e}$$

    The cavity volume $L^3$ has cancelled, as it must: the spectral energy
    density is a property of the radiation, not of the box that holds it.

!!! note "Plain-language version: where the $k_{\mathrm B}T$ per mode comes from"
    Equipartition says each *independent quadratic way of storing energy*
    carries $\tfrac12 k_{\mathrm B}T$ on thermal average. A harmonic
    oscillator stores energy in two such ways at once — kinetic
    ($\propto p^2$) and potential ($\propto x^2$) — so it carries
    $2\times\tfrac12 k_{\mathrm B}T = k_{\mathrm B}T$. Each cavity mode is
    exactly such an oscillator (the field amplitude swings like a mass on a
    spring), hence $k_{\mathrm B}T$ per mode. Crucially, this result does
    *not* depend on the mode's frequency — every mode, however
    high-frequency, gets the same $k_{\mathrm B}T$. That frequency-blindness
    is precisely what breaks at high $\nu$, and precisely what Planck repairs.

By the classical equipartition theorem, each mode is an independent harmonic oscillator and carries an average energy $k_{\mathrm B}T$ in thermal equilibrium. Multiplying gives the Rayleigh–Jeans law,

$$u_{\mathrm{RJ}}(\nu, T) = \frac{8\pi \nu^2}{c^3}\, k_{\mathrm B} T. \tag{4.1.1}$$

At low frequencies — radio waves, microwaves, the red end of the visible spectrum — this formula matches experiment beautifully. But it has a fatal feature: it grows without bound as $\nu \to \infty$. The total radiated energy per unit volume,

$$U = \int_0^\infty u_{\mathrm{RJ}}(\nu, T)\, d\nu = \infty,$$

diverges. This is the "ultraviolet catastrophe": classical physics predicts that *any* warm object should emit an infinite amount of high-frequency radiation. A glowing coal should incinerate the room. It does not.

In December 1900 Max Planck produced a fix. He postulated — at first as a purely mathematical trick — that the energy of each electromagnetic mode of frequency $\nu$ is not continuously variable but comes in discrete lumps of size

$$\varepsilon = h\nu, \tag{4.1.2}$$

where $h \approx 6.626 \times 10^{-34}$ J s. Repeating the equipartition calculation with quantised energies replaces $k_{\mathrm B} T$ by $h\nu /(e^{h\nu/k_{\mathrm B}T} - 1)$. The replacement is not magic; it is a finite sum of a geometric series, derived in full below.

??? note "Full derivation: the Planck mean energy of a mode"
    **Set-up.** Classically a mode could hold *any* energy. Planck allows
    only the discrete values $E_n = n\,h\nu$ for $n = 0, 1, 2, \ldots$ In
    thermal equilibrium the probability of finding the mode in level $n$ is
    given by the Boltzmann factor (normalised so the probabilities sum to
    one):

    $$p_n = \frac{e^{-E_n/k_{\mathrm B}T}}{\displaystyle\sum_{m=0}^{\infty} e^{-E_m/k_{\mathrm B}T}} = \frac{e^{-n h\nu/k_{\mathrm B}T}}{\displaystyle\sum_{m=0}^{\infty} e^{-m h\nu/k_{\mathrm B}T}}. \tag{4.1.3a}$$

    The mean energy is the probability-weighted average $\langle E\rangle = \sum_n E_n p_n$:

    $$\langle E\rangle = \frac{\displaystyle\sum_{n=0}^{\infty} n\,h\nu \, e^{-n h\nu/k_{\mathrm B}T}}{\displaystyle\sum_{n=0}^{\infty} e^{-n h\nu/k_{\mathrm B}T}}. \tag{4.1.3b}$$

    **Shorthand.** Let $x \equiv h\nu/k_{\mathrm B}T$ and $z \equiv e^{-x}$,
    so $0 < z < 1$. Then $\langle E\rangle = h\nu\,\dfrac{\sum_n n z^n}{\sum_n z^n}$.

    **The denominator** is a geometric series. For $|z|<1$,

    $$\sum_{n=0}^{\infty} z^n = \frac{1}{1-z}. \tag{4.1.3c}$$

    **The numerator** is obtained by differentiating the denominator. Apply
    $z\,\dfrac{d}{dz}$ to both sides of (4.1.3c): on the left,
    $z\dfrac{d}{dz}\sum_n z^n = \sum_n n z^n$; on the right,
    $z\dfrac{d}{dz}\dfrac{1}{1-z} = z\cdot\dfrac{1}{(1-z)^2} = \dfrac{z}{(1-z)^2}$. Hence

    $$\sum_{n=0}^{\infty} n\, z^n = \frac{z}{(1-z)^2}. \tag{4.1.3d}$$

    **Take the ratio.** Dividing (4.1.3d) by (4.1.3c),

    $$\frac{\sum_n n z^n}{\sum_n z^n} = \frac{z/(1-z)^2}{1/(1-z)} = \frac{z}{1-z}. \tag{4.1.3e}$$

    **Restore $z = e^{-x}$.** Multiply top and bottom by $e^{x}$:

    $$\frac{z}{1-z} = \frac{e^{-x}}{1-e^{-x}} = \frac{1}{e^{x}-1}. \tag{4.1.3f}$$

    **Result.** Therefore

    $$\langle E\rangle = h\nu\,\frac{z}{1-z} = \frac{h\nu}{e^{h\nu/k_{\mathrm B}T}-1}, \tag{4.1.3g}$$

    which is exactly the replacement for $k_{\mathrm B}T$ quoted in the text.

    **Recover the classical limit.** As $h\nu \to 0$ (equivalently
    $x\to 0$), expand the exponential to first order:
    $e^{x} - 1 = (1 + x + \tfrac12 x^2 + \cdots) - 1 \approx x$. Then

    $$\langle E\rangle = \frac{h\nu}{e^{x}-1} \approx \frac{h\nu}{x} = \frac{h\nu}{h\nu/k_{\mathrm B}T} = k_{\mathrm B}T. \tag{4.1.3h}$$

    So when the quantum step $h\nu$ is tiny compared with the thermal energy
    $k_{\mathrm B}T$, the discreteness is invisible and equipartition's
    $k_{\mathrm B}T$ per mode is recovered. The high-frequency modes, where
    $h\nu \gg k_{\mathrm B}T$, are instead *frozen out*: $\langle E\rangle
    \approx h\nu\, e^{-h\nu/k_{\mathrm B}T} \to 0$. That exponential freeze
    is what tames the ultraviolet catastrophe.

Multiplying the mode density $g(\nu) = 8\pi\nu^2/c^3$ by this mean energy $h\nu/(e^{h\nu/k_{\mathrm B}T}-1)$ gives the Planck law,

$$u(\nu, T) = \frac{8\pi h \nu^3}{c^3}\, \frac{1}{e^{h\nu / k_{\mathrm B} T} - 1}, \tag{4.1.3}$$

reduces to Rayleigh–Jeans at low frequency and cuts off exponentially at high frequency. The total energy is finite and agrees with the Stefan–Boltzmann law. The data fit is perfect.

!!! note "Why this matters for materials"
    Equation (4.1.2) is the first crack in classical physics. Whatever else is true, energy at the atomic scale comes in *discrete amounts*. This is the seed of every band gap, every vibrational quantum, every laser line. When we compute the electronic levels of a solid in Chapter 5 we are computing the descendants of Planck's quanta.

Planck himself was uncomfortable with his own postulate and spent years trying to derive it from classical physics. He failed, because there is no such derivation: nature is not classical.

!!! tip "Reading the Planck spectrum"
    There are two useful limiting forms of (4.1.3). For $h\nu \ll k_{\mathrm B}T$ we may expand the exponential as $e^{h\nu/k_{\mathrm B}T} \approx 1 + h\nu/k_{\mathrm B}T$ and recover Rayleigh–Jeans:
    $$u(\nu, T) \approx \frac{8\pi\nu^2}{c^3}\, k_{\mathrm B}T, \qquad (h\nu \ll k_{\mathrm B}T).$$
    For $h\nu \gg k_{\mathrm B}T$ the denominator $e^{h\nu/k_{\mathrm B}T} - 1 \approx e^{h\nu/k_{\mathrm B}T}$ and the spectrum decays exponentially: this is **Wien's law**, $u(\nu,T) \propto \nu^3 e^{-h\nu/k_{\mathrm B}T}$. The crossover happens at $h\nu \sim k_{\mathrm B}T$, which for $T = 300$ K is $\nu \sim 6\times 10^{12}$ Hz — in the far infrared. This is why room-temperature objects glow in the IR, not the visible: their thermal energy $k_{\mathrm B}T \approx 0.026$ eV is well below visible photon energies $\sim 2$ eV.

The numerical check is instructive. Differentiating $u(\nu, T)$ with respect to $\nu$ and setting the derivative to zero yields **Wien's displacement law**, $h\nu_{\max} \approx 2.82\, k_{\mathrm B}T$. The number $2.82$ is not arbitrary: it is the root of a simple transcendental equation, derived below.

??? note "Full derivation: Wien's constant $2.82$ as a root of $(3-x) = 3e^{-x}$"
    **Reduce to one variable.** Write $u(\nu,T) \propto \nu^3/(e^{x}-1)$
    with $x = h\nu/k_{\mathrm B}T$. At fixed $T$, the peak in $\nu$ is the
    peak of $f(x) = x^3/(e^{x}-1)$, because $\nu \propto x$. We seek the $x$
    where $f'(x) = 0$.

    **Differentiate using the quotient rule.** With numerator $x^3$ and
    denominator $g(x) = e^{x}-1$,

    $$f'(x) = \frac{3x^2(e^{x}-1) - x^3 e^{x}}{(e^{x}-1)^2}. \tag{4.1.W1}$$

    **Set the numerator to zero** (the denominator is never zero for $x>0$):

    $$3x^2(e^{x}-1) - x^3 e^{x} = 0. \tag{4.1.W2}$$

    **Divide by $x^2$** (the $x=0$ root is the trivial minimum, not the peak):

    $$3(e^{x}-1) - x\,e^{x} = 0 \;\Longrightarrow\; 3e^{x} - 3 - x e^{x} = 0 \;\Longrightarrow\; (3 - x)\,e^{x} = 3. \tag{4.1.W3}$$

    **Rearrange** by multiplying through by $e^{-x}$:

    $$3 - x = 3\,e^{-x}. \tag{4.1.W4}$$

    **Solve numerically.** This cannot be solved in closed form, but it is
    a one-line fixed-point iteration $x \mapsto 3(1 - e^{-x})$. Starting from
    $x_0 = 3$ gives $2.852,\ 2.821,\ 2.8214,\ldots$, converging to

    $$x_{\max} = 2.8214\ldots \tag{4.1.W5}$$

    Since $x_{\max} = h\nu_{\max}/k_{\mathrm B}T$, this is exactly
    $h\nu_{\max} \approx 2.82\,k_{\mathrm B}T$.

    *A subtlety worth flagging.* The peak of $u(\nu,T)$ in *frequency* and
    the peak of the corresponding $u(\lambda,T)$ in *wavelength* are **not**
    at the same physical place, because $d\nu$ and $d\lambda$ are related by
    $\nu = c/\lambda$, so $u(\lambda) = u(\nu)\,|d\nu/d\lambda| =
    u(\nu)\,c/\lambda^2$. Redoing the maximisation in $\lambda$ gives the
    *different* root $4.965$ (the familiar $\lambda_{\max}T = 2.898\times
    10^{-3}$ m K). The number you quote depends on whether you bin per unit
    frequency or per unit wavelength; this page uses the frequency form
    throughout.

For the surface of the Sun ($T \approx 5800$ K), $\nu_{\max} \approx 3.4\times 10^{14}$ Hz, i.e.\ $\lambda_{\max} \approx 880$ nm — in the near-infrared, with the peak energy density spilling into the visible. Human vision and the solar spectrum are matched, and the match is set by a single equation that classical physics could not have produced.

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

!!! warning "Common misunderstanding: 'light is *just* a wave'"
    Maxwell's equations, interference and diffraction make light look like a
    pure wave, and for those phenomena it *is* one. The photoelectric effect
    shows the wave picture is incomplete, not wrong: in an absorption or
    emission event, light is delivered in indivisible lumps of energy $h\nu$.
    The correct statement is **both/and**, not either/or — light is a wave
    *and* its energy is quantised. The same will be true of the electron in
    the opposite direction: a particle whose behaviour is governed by a wave.
    Neither "wave" nor "particle" alone is the full story; that is the point
    of *wave–particle duality*.

The deeper significance of Einstein's photon is not the photoelectric effect itself — that is one experiment — but the fact that *energy quantisation propagates*. Planck quantised the oscillators in the cavity walls; Einstein quantised the field that travels between them. Once light itself comes in lumps $h\nu$, the symmetry suggesting "matter is also wavelike" is hard to avoid. We will see de Broglie complete that symmetry in §4.1.4.

## 4.1.3 The atom should not exist

The third — and for materials physics most acute — crisis concerns the stability of atoms.

!!! warning "Common misunderstanding: the electron is *not* a tiny billiard ball"
    The "miniature solar system" picture below is the one this whole
    sub-section is about to demolish. It is a useful target precisely
    because it is wrong: a classical point electron on a fixed circular
    track would radiate and spiral into the nucleus in about $10^{-11}$ s.
    Hold the planetary image loosely — it gets the length and energy scales
    right (next page) but the *mechanism* wrong. The electron is a
    [wavefunction](../undergraduate/glossary-for-beginners.md), not a ball
    on a wire.

The symbols used throughout this sub-section, in SI units, are as follows.

| Symbol | Meaning | Units |
|---|---|---|
| $e$ | elementary charge, $1.602\times10^{-19}$ | C |
| $\varepsilon_0$ | vacuum permittivity, $8.854\times10^{-12}$ | F m$^{-1}$ |
| $m_{\mathrm e}$ | electron mass, $9.109\times10^{-31}$ | kg |
| $r$ | orbit radius | m |
| $v$ | orbital speed | m s$^{-1}$ |
| $a$ | (centripetal) acceleration | m s$^{-2}$ |
| $P$ | radiated power (Larmor) | W |
| $a_0$ | Bohr radius, $0.529\times10^{-10}$ | m |
| $\hbar = h/2\pi$ | reduced Planck constant, $1.055\times10^{-34}$ | J s |
| $L = m_{\mathrm e}vr$ | orbital angular momentum | J s |
| $E_n$ | energy of the $n$-th Bohr level | J (quoted in eV) |
| $R$ | Rydberg constant, $1.097\times10^{7}$ | m$^{-1}$ |

By 1911 Ernest Rutherford's gold-foil experiments had established that an atom consists of a tiny dense positive nucleus surrounded by negative electrons. The natural classical picture is a miniature solar system: electrons orbit the nucleus under the Coulomb attraction, much as planets orbit the sun under gravity. The maths is identical: closed elliptical orbits with energies given by Kepler's laws.

This picture is catastrophically wrong. An electron in circular orbit is an accelerating charge, and an accelerating charge radiates electromagnetic waves — this is exactly how a radio antenna works. The Larmor formula gives the radiated power:

$$P = \frac{e^2 a^2}{6\pi \varepsilon_0 c^3}, \tag{4.1.5}$$

where $a$ is the acceleration. Plug in the numbers for a hydrogen atom: an electron at the Bohr radius $a_0 \approx 0.529$ Å, circling at the velocity required to balance the Coulomb force, has centripetal acceleration $a \sim 9 \times 10^{22}$ m s$^{-2}$. The radiated power is enormous. As it loses energy the electron spirals inward, accelerating still more and radiating still faster. A straightforward integration gives the lifetime of a classical hydrogen atom:

$$\tau \sim 10^{-11}\ \mathrm{s}. \tag{4.1.6}$$

Both the acceleration and the lifetime are worked through, number by number, below — the result $\tau \approx 1.6\times10^{-11}$ s is the one Bohr and his contemporaries found so impossible to live with.

??? note "Full derivation: the acceleration $\sim 9\times10^{22}$ m s$^{-2}$ and the lifetime $\tau \approx 1.6\times10^{-11}$ s"
    **The centripetal acceleration.** From the force balance (the same
    equation as the Bohr Postulate 1 on the next page),
    $\dfrac{e^2}{4\pi\varepsilon_0 r^2} = m_{\mathrm e}a$, so

    $$a = \frac{e^2}{4\pi\varepsilon_0 m_{\mathrm e} r^2}. \tag{4.1.6a}$$

    Putting $r = a_0 = 0.529\times10^{-10}$ m,
    $e = 1.602\times10^{-19}$ C, $\varepsilon_0 = 8.854\times10^{-12}$
    F m$^{-1}$, $m_{\mathrm e} = 9.109\times10^{-31}$ kg:

    $$a = \frac{(1.602\times10^{-19})^2}{4\pi(8.854\times10^{-12})(9.109\times10^{-31})(0.529\times10^{-10})^2} \approx 9.0\times10^{22}\ \text{m s}^{-2}.$$

    **The energy as a function of radius.** A classical circular-orbit
    electron has total energy equal to half the (negative) potential energy
    — this is the direct force-balance result derived on the next page, not
    a separate assumption:

    $$E(r) = -\frac{e^2}{8\pi\varepsilon_0\, r}. \tag{4.1.6b}$$

    Differentiating, a small inward change $dr<0$ *lowers* $E$:

    $$\frac{dE}{dr} = \frac{e^2}{8\pi\varepsilon_0\, r^2}. \tag{4.1.6c}$$

    **Equate energy loss to radiated power.** Energy conservation says
    $dE/dt = -P$, where $P$ is the Larmor power (4.1.5). Using
    $dE/dt = (dE/dr)(dr/dt)$ and substituting $a$ from (4.1.6a) into (4.1.5):

    $$\frac{e^2}{8\pi\varepsilon_0\, r^2}\frac{dr}{dt} = -\frac{e^2 a^2}{6\pi\varepsilon_0 c^3} = -\frac{e^2}{6\pi\varepsilon_0 c^3}\left(\frac{e^2}{4\pi\varepsilon_0 m_{\mathrm e} r^2}\right)^2. \tag{4.1.6d}$$

    **Cancel and tidy.** Cancel the common factor
    $e^2/(\pi\varepsilon_0)$ from both sides, leaving

    $$\frac{1}{8 r^2}\frac{dr}{dt} = -\frac{1}{6 c^3}\cdot \frac{e^4}{(4\pi\varepsilon_0)^2 m_{\mathrm e}^2 r^4}.$$

    Multiply through by $8r^2$:

    $$\frac{dr}{dt} = -\frac{8}{6 c^3}\cdot\frac{e^4}{(4\pi\varepsilon_0)^2 m_{\mathrm e}^2 r^2} = -\frac{4}{3}\,\frac{e^4}{(4\pi\varepsilon_0)^2 m_{\mathrm e}^2 c^3}\,\frac{1}{r^2}. \tag{4.1.6e}$$

    **Separate variables and integrate.** Write $r^2\,dr = -K\,dt$ with the
    constant

    $$K \equiv \frac{4}{3}\,\frac{e^4}{(4\pi\varepsilon_0)^2 m_{\mathrm e}^2 c^3}.$$

    Integrating from the starting radius $a_0$ down to $0$ over the lifetime
    $\tau$:

    $$\int_{a_0}^{0} r^2\,dr = -K\int_0^{\tau} dt \;\Longrightarrow\; -\frac{a_0^3}{3} = -K\tau \;\Longrightarrow\; \tau = \frac{a_0^3}{3K}. \tag{4.1.6f}$$

    **Put in the numbers.** First evaluate
    $K = \tfrac43\,(1.602\times10^{-19})^4 / [(4\pi\cdot8.854\times10^{-12})^2(9.109\times10^{-31})^2(2.998\times10^8)^3]$.
    Working it out, $K \approx 3.17\times10^{-21}$ m$^3$ s$^{-1}$. Then,
    with $a_0^3 = (0.529\times10^{-10})^3 \approx 1.48\times10^{-31}$ m$^3$,

    $$\tau = \frac{1.48\times10^{-31}}{3\times 3.17\times10^{-21}} \approx 1.6\times10^{-11}\ \text{s},$$

    in agreement with the quoted estimate (4.1.6). Classical electrodynamics
    really does predict that hydrogen collapses in about sixteen
    picoseconds.

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

!!! note "Why this step? — the kinetic energy without invoking the virial theorem"
    You do not need to know the virial theorem to get
    $\tfrac12 m_{\mathrm e}v^2 = e^2/(8\pi\varepsilon_0 r)$; it falls straight
    out of the force balance (4.1.A). That equation reads
    $\dfrac{e^2}{4\pi\varepsilon_0 r^2} = \dfrac{m_{\mathrm e}v^2}{r}$.
    Multiply both sides by $r$:

    $$\frac{e^2}{4\pi\varepsilon_0\, r} = m_{\mathrm e} v^2.$$

    Now halve it:

    $$\underbrace{\tfrac12 m_{\mathrm e} v^2}_{\text{kinetic energy}} = \frac{e^2}{8\pi\varepsilon_0\, r}. \tag{4.1.D1}$$

    So the kinetic energy is exactly *half* the magnitude of the potential
    energy $U = -e^2/(4\pi\varepsilon_0 r)$ — which is all the virial theorem
    was asserting here. Adding kinetic and potential then gives the total
    energy used in (4.1.D):

    $$E = \tfrac12 m_{\mathrm e}v^2 + U = \frac{e^2}{8\pi\varepsilon_0 r} - \frac{e^2}{4\pi\varepsilon_0 r} = -\frac{e^2}{8\pi\varepsilon_0 r}. \tag{4.1.D2}$$

    The minus sign matters: a *bound* electron has negative total energy, so
    energy must be *supplied* to ionise it. To reach the final numerical
    value, substitute $r_1 = a_0$ into $E_1 = -e^2/(8\pi\varepsilon_0 a_0)$
    and convert from joules to electronvolts by dividing by
    $e = 1.602\times10^{-19}$ C; the result is $-13.6$ eV.

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

exactly the pattern produced by waves of wavelength $\lambda$. The fringe-spacing formula is short enough to derive in three lines, so we do.

??? note "Full derivation: fringe positions $y_n = n\lambda L/d$ from the path difference"
    **The geometry.** Two slits sit a distance $d$ apart; the detector is a
    distance $L$ behind them, with $L \gg d$. Consider a point on the
    detector a height $y$ above the central axis. The two paths from the
    slits to that point differ in length because one slit is slightly
    farther from the point than the other.

    **The path difference.** For $L \gg d$ the two rays are almost parallel,
    leaving the slits at a common small angle $\theta$ to the axis. Simple
    trigonometry on the little right-angled triangle between the slits gives
    a path difference

    $$\Delta = d\,\sin\theta. \tag{4.1.S1}$$

    **Constructive interference.** Bright fringes occur where the two waves
    arrive *in phase*, i.e. where the path difference is a whole number of
    wavelengths:

    $$d\,\sin\theta = n\lambda, \qquad n = 0, \pm1, \pm2, \ldots \tag{4.1.S2}$$

    **Small-angle approximation.** Because $L \gg d$, the angle $\theta$ is
    tiny, and the height on the detector is $y = L\tan\theta \approx
    L\sin\theta \approx L\theta$. Hence $\sin\theta \approx y/L$. Substitute
    into (4.1.S2):

    $$d\,\frac{y_n}{L} = n\lambda \;\Longrightarrow\; y_n = n\,\frac{\lambda L}{d}, \tag{4.1.S3}$$

    which is the quoted result. Neighbouring bright fringes are therefore
    evenly spaced by $\Delta y = \lambda L/d$: a *direct* readout of the
    wavelength $\lambda$, and hence (via $\lambda = h/p$) of the particle's
    momentum.

For light this is unsurprising and was already understood in Young's day. For *electrons* it was astonishing. And the astonishment deepens once one reduces the beam intensity to the point where only one electron is in the apparatus at any time: the individual arrivals are still detected as point particles (a single dot of phosphorescence on the screen), but the *distribution* of dots accumulated over many runs is the same interference pattern. Each electron seems to know about both slits, even though it is "indivisible".

The Davisson–Germer experiment of 1927 was a controlled accident. Davisson and Germer at Bell Labs were studying low-energy electron scattering from polycrystalline nickel when a vacuum leak forced them to anneal the sample, accidentally producing a nearly single-crystal surface. On resuming the experiment they observed sharp angular peaks in the scattered electron intensity at angles satisfying the Bragg condition $n\lambda = d\sin\theta$, with $\lambda$ given by de Broglie's formula. Electrons were diffracting from the nickel lattice exactly as X-rays would. Independently and almost simultaneously, G.\ P.\ Thomson — son of J.\ J.\ Thomson, who had identified the electron as a particle in 1897 — observed transmission electron diffraction through thin foils. Father and son shared Nobel prizes for proving, respectively, that the electron is a particle and that it is a wave.

By the 1990s the same experiment had been performed on whole C$_{60}$ fullerene molecules in Anton Zeilinger's laboratory: 60 atoms, mass 720 amu, diffracting through a transmission grating with the predicted wavelength. As of the 2020s, molecules of $\sim 25\,000$ amu have been diffracted. There is no known mass scale at which wave behaviour switches off.

### The correspondence principle

If matter is wavelike, why does classical mechanics work so well for tennis balls? The answer is the **correspondence principle**, articulated by Bohr: quantum mechanics must reduce to classical mechanics in the limit where actions are large compared to $\hbar$.

Quantitatively, the relevant small parameter is $\hbar/S$, where $S$ is a characteristic classical action of the system (with dimensions of energy × time, or momentum × length). For a tennis ball of mass $0.06$ kg moving at $20$ m s$^{-1}$ over a court of length $24$ m, the action is $S \sim m v L \sim 30$ J s, while $\hbar \approx 10^{-34}$ J s. The ratio is $10^{-35}$. Wave effects are present in principle but invisible in practice: the de Broglie wavelength is $\lambda = h/(mv) \sim 5\times 10^{-34}$ m, twenty orders of magnitude below the diameter of an atomic nucleus.

For an electron in a hydrogen atom, in contrast, $S \sim p\, a_0 \sim \hbar$ (this is in fact one definition of the atomic scale), and quantum effects are unavoidable. The crossover between the two regimes — the realm where the symbols of classical mechanics begin to lose their meaning and the symbols of quantum mechanics take over — is exactly where modern condensed-matter physics lives.

Two operational forms of the correspondence principle will recur.

- **Large quantum numbers.** For energy levels labelled by an integer $n$, the spacing $E_{n+1} - E_n$ relative to $E_n$ shrinks as $n \to \infty$, and the quantum probability density (averaged over a few wavelengths) approaches the classical position distribution $\rho_{\mathrm{cl}}(x) \propto 1/v_{\mathrm{cl}}(x)$. The $1/v$ here has a simple meaning: a classical particle spends, in each little interval $dx$, a time $dt = dx/v(x)$, so the fraction of its period it is "caught" near $x$ — and hence the probability of finding it there — is proportional to $1/v(x)$. The particle lingers where it moves slowly (near the turning points) and is rarely caught where it is fast (near the centre). We will see this explicitly for the harmonic oscillator in §4.4.
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

!!! tip "Where this appears later"
    The matter-wave idea launched here becomes the working tool of the rest
    of the book. The Schrödinger equation that subsumes these four crises is
    built in [Chapter 4's later sections](index.md) and solved on a
    computer; its many-electron version is the starting point for
    [Chapter 5 (DFT)](../ch05-dft/index.md). The quantised energy levels
    (Planck, Bohr) reappear as the [band structure](../undergraduate/glossary-for-beginners.md)
    of solids, and the wave nature of the electron is what makes
    [reciprocal space](../undergraduate/glossary-for-beginners.md) and
    Brillouin zones the natural language of [Chapter 3b (solid-state)](../ch03b-solid-state/index.md).
    The correspondence-principle limit ($\hbar\to0$) is exactly what lets
    [Chapter 7 (molecular dynamics)](../ch07-md/index.md) treat nuclei as
    classical balls while the electrons remain quantum.

!!! question "Check yourself"
    1. In the Rayleigh–Jeans derivation, two ingredients are multiplied
       together to get $u(\nu,T)$. Name them, and say which one Planck
       changed and which he left untouched.
    2. Show, in one line, that the Planck mean energy
       $\langle E\rangle = h\nu/(e^{h\nu/k_{\mathrm B}T}-1)$ tends to
       $k_{\mathrm B}T$ when $h\nu \ll k_{\mathrm B}T$.
    3. A classical orbiting electron radiates and spirals in. Roughly how
       long does a hydrogen atom survive on this picture, and what one
       experimental fact does this contradict?
    4. The number $2.82$ in Wien's law is the root of which equation? Why is
       the root for the *wavelength* peak a different number?
    5. A 100 eV electron has a de Broglie wavelength of about $1.2$ Å. Why
       does this make electrons useful for studying crystals, while a
       thrown tennis ball never visibly diffracts?
    6. In the double-slit experiment, the fringe spacing is
       $\Delta y = \lambda L/d$. If you halved the slit separation $d$,
       what would happen to the spacing of the bright fringes?

??? success "Answers"
    1. The **mode density** $g(\nu) = 8\pi\nu^2/c^3$ (how many standing
       waves of each frequency fit in the box) and the **mean energy per
       mode**. Planck changed *only the mean energy per mode*, replacing the
       equipartition value $k_{\mathrm B}T$ with
       $h\nu/(e^{h\nu/k_{\mathrm B}T}-1)$. The mode counting is unchanged.
    2. Put $x = h\nu/k_{\mathrm B}T \ll 1$ and use $e^{x}-1 \approx x$:
       $\langle E\rangle = h\nu/(e^{x}-1) \approx h\nu/x = k_{\mathrm B}T$.
    3. About $1.6\times10^{-11}$ s (some tens of picoseconds). This
       contradicts the plain fact that atoms are stable — hydrogen has
       existed undisturbed for billions of years — and that hydrogen emits
       a *sharp line* spectrum, not the smooth continuum a spiralling
       electron would radiate.
    4. The frequency peak satisfies $3 - x = 3e^{-x}$, whose non-trivial
       root is $x = 2.8214\ldots$ The wavelength peak uses a *different*
       quantity, $u(\lambda,T) = u(\nu,T)\,c/\lambda^2$, because converting
       from "per unit frequency" to "per unit wavelength" introduces an
       extra $1/\lambda^2$; maximising that gives the different root
       $4.965\ldots$
    5. The de Broglie wavelength $\lambda = h/p$ is comparable to atomic
       spacings ($\sim 1$ Å) for a 100 eV electron, so it diffracts from the
       lattice and reveals crystal structure. A tennis ball has an enormous
       momentum, so its wavelength is $\sim 10^{-34}$ m — twenty orders of
       magnitude smaller than a nucleus — and no slit or lattice is fine
       enough to reveal the wave (the correspondence principle: $\hbar/S \to
       0$).
    6. The spacing $\Delta y = \lambda L/d$ is *inversely* proportional to
       $d$, so halving $d$ would **double** the fringe spacing.

??? note "Hint"
    For 2, look at the "Recover the classical limit" step of the Planck
    mean-energy derivation. For 4, see the folded Wien derivation and the
    note that follows it. For 5 and 6, the de Broglie relation
    $\lambda = h/p$ and the fringe formula $\Delta y = \lambda L/d$ are all
    you need.

!!! tip "What to take from this section"
    The historical detour through Planck, Einstein, Bohr and de Broglie is not for nostalgia. Each of the four episodes identified a *specific* failure of classical physics — UV catastrophe, photoelectric effect, atomic stability, electron diffraction — and the *specific* fix needed. In §4.2 we will write down a single equation that contains all four fixes simultaneously. Knowing the historical motivation is not a luxury; it is what makes the equation legible. The same will be true in Chapter 5, where the failure modes of HF (correlation, self-interaction, dispersion) motivate the structure of DFT functionals. Throughout the book: read the failure first, then the fix.
