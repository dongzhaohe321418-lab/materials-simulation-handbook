# 2.3 Three Diagrams Every Materials Scientist Reads

Walk into any materials-science seminar and you are likely to see at least one of three diagrams on the projector: a *phase diagram*, a *band structure with density of states*, or a *radial distribution function*. Between them, they encode the thermodynamic, electronic, and structural information that the field cares about most. This section teaches you to read all three, with worked examples and short Python snippets that produce toy versions you can play with.

## Diagram one: the phase diagram

A phase diagram is a map of which phases of a material are stable as a function of thermodynamic variables — usually composition, temperature, and pressure. It is the single most important diagram in classical materials science. Metallurgists draw them on the backs of envelopes.

### What is on the axes

A *binary phase diagram* — two components — typically plots composition (atomic or weight fraction of component B) on the horizontal axis and temperature on the vertical axis, with pressure held fixed at one atmosphere. Each region of the plane is labelled with the stable phase or mixture of phases: a single solid solution, a two-phase mixture, a liquid, or a mixture of solid and liquid.

Three features dominate any binary diagram you will encounter:

- the **liquidus**, the upper boundary above which everything is liquid;
- the **solidus**, the lower boundary below which everything is solid;
- the **two-phase region** between them, in which solid and liquid coexist.

In simple eutectic systems, the liquidus and solidus meet at a single composition called the *eutectic point*, where solid, liquid, and a second solid phase coexist at a single temperature. In other systems, intermetallic compounds appear as vertical lines or narrow regions of stoichiometric stability.

### Worked example: reading a Cu–Zn diagram

The copper–zinc system, which underlies the brasses, is one of the more intricate binary phase diagrams. We will work with a simplified version.

**Figure 2.3.1.** Imagine a binary phase diagram with composition (atomic percent zinc) on the horizontal axis from 0 to 100, and temperature in degrees Celsius on the vertical axis from 0 to 1100. The left edge represents pure copper (melting at 1085 °C); the right edge pure zinc (melting at 420 °C). The liquidus is a single curve descending from 1085 °C on the left to a minimum at around 902 °C near 36 atomic percent zinc, then continuing down to 420 °C on the right. Below the liquidus, several solid regions are labelled: $\alpha$ (FCC Cu solid solution) on the left up to about 32 at.% Zn; $\beta$ (BCC) between 36 and 56 at.% Zn; $\gamma$, $\delta$, $\varepsilon$, $\eta$ at higher Zn content. Two-phase regions sit between them.

Suppose we want to know what phases are present in a brass containing 30 at.% Zn at 600 °C.

**Step 1: locate the point.** Mark the intersection of the vertical line at 30 at.% Zn and the horizontal line at 600 °C.

**Step 2: identify the region.** The point lies inside the $\alpha$ phase field (the FCC solid solution of zinc in copper). So at 600 °C, this composition is a single solid phase — copper-rich FCC with zinc atoms substitutionally dissolved at random.

**Step 3: think about what changes at higher temperature.** Raising the temperature at fixed composition, we eventually cross the $\alpha$ / ($\alpha + L$) boundary, where the alloy partially melts. Above the liquidus, everything is liquid.

**Step 4: the lever rule (for two-phase regions).** If we had instead picked 40 at.% Zn at 900 °C, we would land in the two-phase ($\alpha + L$) region. The fractions of each phase are given by the *lever rule*: draw the horizontal tie-line at 900 °C; the fraction of liquid equals the distance from the alloy composition to the $\alpha$ boundary, divided by the total length of the tie-line. This is the same statement as the conservation of mass.

The lever rule, applied repeatedly, predicts how phase fractions evolve as the alloy is cooled — the basis of equilibrium solidification calculations.

### A toy phase diagram in Python

A real Cu–Zn diagram is computed from a CALPHAD database with dozens of parameters. We can illustrate the structure of a simple eutectic diagram analytically. The idea: each phase has a Gibbs free energy as a function of composition; the equilibrium phase or phase mixture at each $T$ minimises the total free energy by common-tangent construction.

```python
from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt

def gibbs_solid(x: np.ndarray, T: float) -> np.ndarray:
    """Regular-solution free energy of a generic solid phase."""
    R = 8.314  # J / (mol K)
    omega = 25_000.0  # interaction parameter, J / mol
    g_mix = R * T * (x * np.log(x + 1e-12) + (1 - x) * np.log(1 - x + 1e-12))
    return omega * x * (1 - x) + g_mix

def gibbs_liquid(x: np.ndarray, T: float) -> np.ndarray:
    """Ideal-solution free energy of the liquid, shifted by latent heat."""
    R = 8.314
    Tm_A, Tm_B = 1358.0, 693.0  # melting points in K
    L_A, L_B = 13_000.0, 7_300.0  # latent heats in J / mol
    g_pure = (1 - x) * L_A * (1 - T / Tm_A) + x * L_B * (1 - T / Tm_B)
    g_mix = R * T * (x * np.log(x + 1e-12) + (1 - x) * np.log(1 - x + 1e-12))
    return g_pure + g_mix

x = np.linspace(1e-3, 1 - 1e-3, 400)
for T in (500.0, 800.0, 1100.0):
    plt.plot(x, gibbs_solid(x, T), label=f"solid, T={T:.0f} K")
    plt.plot(x, gibbs_liquid(x, T), "--", label=f"liquid, T={T:.0f} K")
plt.xlabel("composition x_B")
plt.ylabel("G (J / mol)")
plt.legend(fontsize=7)
plt.tight_layout()
```

Running this for several temperatures and constructing common tangents between the lowest-lying curves gives the phase boundaries. We will not reproduce a full CALPHAD calculation, but the snippet exposes the underlying ideas: free energy is the master function; phase boundaries follow from minimisation.

!!! note "Beyond binary diagrams"
    Ternary diagrams use the Gibbs triangle to represent three-component compositions; pressure–temperature diagrams omit composition and are used for one-component systems like water or carbon. The same principles apply: equilibrium is set by free-energy minimisation, and phase boundaries follow from common-tangent constructions.

## Diagram two: band structure and density of states

If the phase diagram is the headline diagram of classical materials science, the band structure is its electronic counterpart. It tells you whether a material is a metal, a semiconductor, or an insulator; it suggests where optical transitions will occur; and it is the starting point for understanding electrical and thermal transport.

### What is being plotted

A band structure plots single-electron energies $E_n(\mathbf{k})$ as a function of crystal momentum $\mathbf{k}$, for each band index $n$. The horizontal axis is a one-dimensional path through the *Brillouin zone* — a region of reciprocal space that we will define in Chapter 3 — visiting *high-symmetry points* such as $\Gamma$, X, L, and K. The vertical axis is energy, conventionally with the Fermi level placed at zero.

A *density of states* (DOS) is the integral of the band structure over all $\mathbf{k}$, projected onto energy:

$$
g(E) = \sum_n \int_\mathrm{BZ} \frac{\mathrm{d}^3 \mathbf{k}}{(2\pi)^3} \, \delta(E - E_n(\mathbf{k})).
$$

It records how many states are available at each energy, without keeping the momentum information. Band structures and DOS are usually shown side by side: the band structure on the left, with energy on the vertical axis, and the DOS plotted with its axes swapped to share that vertical scale on the right.

### What to look at

When you see a band structure for the first time, ask these questions in order.

1. **Is there a gap at the Fermi level?** Look at the band that is filled to $E_\mathrm{F}$ (occupied valence band) and the next one above (conduction band). If they touch, the material is a metal. If they are separated by an energy range with no states, the material is a semiconductor or insulator, and the size of the gap is the *band gap*.

2. **If there is a gap, is it direct or indirect?** A *direct* gap has its valence-band maximum and conduction-band minimum at the same $\mathbf{k}$. An *indirect* gap has them at different $\mathbf{k}$ points. Direct-gap materials (GaAs, the perovskites, monolayer MoS$_2$) absorb and emit light efficiently; indirect-gap materials (Si, Ge) need phonons to mediate optical transitions, which makes them poor light emitters.

3. **What is the band curvature near the band edges?** The curvature gives the *effective mass*: $1 / m^* = (1 / \hbar^2) \, \partial^2 E / \partial k^2$. Flat bands mean heavy carriers; steep bands mean light carriers. This determines mobility.

4. **What is the character of the bands?** Modern band-structure plots colour each band by its dominant orbital contribution (s, p, d) or by atomic species. This tells you which atoms and orbitals dominate transport at each energy.

### Worked example: silicon

Silicon is the canonical worked example for a band structure. It is FCC with a two-atom basis (diamond structure), and the relevant high-symmetry path runs $L \to \Gamma \to X \to U \to K \to \Gamma$.

**Figure 2.3.2.** Imagine a band-structure plot with $k$ on the horizontal axis traversing the path L–$\Gamma$–X–U–K–$\Gamma$ from left to right, and energy on the vertical axis from $-12$ eV to $+5$ eV with the Fermi level at $0$. Four valence bands and several conduction bands are drawn. The valence-band maximum is at $\Gamma$ (the point in the middle of the horizontal axis), at $E = 0$. The conduction-band minimum is along the $\Gamma$–X line, about 85% of the way to X. The two are at different $\mathbf{k}$ points: silicon has an *indirect* gap. Adjacent to the band structure on the right, the DOS shows a broad valence-band region from about $-12$ to $0$ eV with structure, a gap, and conduction-band states above $\sim 1.1$ eV.

For silicon, you should be able to read off:

- *Type*: semiconductor with an indirect gap.
- *Gap size*: $1.1$ eV experimentally; standard PBE-DFT gives $0.6$ eV (the usual underestimate); HSE06 hybrid functional gives $\sim 1.2$ eV.
- *Effective masses*: light electrons (steep conduction band), light and heavy holes (two valence bands with different curvatures at $\Gamma$).
- *Orbital character*: valence band is mostly Si 3p; conduction band has mixed 3s and 3p character.

These facts together explain why silicon is the workhorse of microelectronics (a moderate gap is easy to dope and gate) but a poor light emitter (the indirect gap requires phonon assistance).

### A toy band structure in Python

We can construct an instructive band structure with a one-dimensional tight-binding chain. Consider a chain of identical sites with nearest-neighbour hopping $t$:

$$
E(k) = -2t \cos(ka).
$$

```python
import numpy as np
import matplotlib.pyplot as plt

a, t = 1.0, 1.0  # lattice constant, hopping
k = np.linspace(-np.pi / a, np.pi / a, 401)
E = -2 * t * np.cos(k * a)

fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(6, 4))
ax1.plot(k * a / np.pi, E, color="C0")
ax1.set_xlabel("k a / pi")
ax1.set_ylabel("E / t")
ax1.set_title("band structure (1D tight binding)")

# density of states: g(E) = 1 / (pi sqrt((2t)^2 - E^2))
E_grid = np.linspace(-2 * t + 1e-3, 2 * t - 1e-3, 401)
g = 1.0 / (np.pi * np.sqrt((2 * t) ** 2 - E_grid ** 2))
ax2.plot(g, E_grid, color="C1")
ax2.set_xlabel("g(E)")
ax2.set_title("DOS")
plt.tight_layout()
```

The DOS for a one-dimensional band has logarithmic *van Hove singularities* at the band edges $E = \pm 2t$, where the band is flat. Van Hove singularities at the Fermi level are sometimes invoked as drivers of phase transitions; this toy model lets you see one directly.

!!! tip "Reading real band-structure plots"
    When you encounter your first DFT band structure for a real material — silicon, graphene, MoS$_2$, the perovskite of the month — work through the four questions above explicitly. After a dozen materials, the answers will come at a glance. The exercises in this chapter give you a head start.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. What is plotted on each axis of a band structure, and how do you tell at a glance whether a material is a metal or a semiconductor?
    2. What is the difference between a direct and an indirect band gap, and why does silicon's indirect gap make it a poor light emitter?
    3. How is the effective mass of a charge carrier related to the curvature of a band, and what does a flat band imply about mobility?

    If any of these is shaky, re-read the preceding section before continuing.

## Diagram three: the radial distribution function

The third pillar is the *radial distribution function* (RDF), also called the *pair correlation function* and denoted $g(r)$. Unlike the previous two diagrams, the RDF is fundamentally a statistical object: it characterises a *snapshot* (or a time-average of snapshots) of atomic positions.

### Definition

For a homogeneous, isotropic system of $N$ atoms in a volume $V$, the radial distribution function is

$$
g(r) = \frac{V}{N^2} \left\langle \sum_{i \ne j} \delta(r - |\mathbf{r}_i - \mathbf{r}_j|) \right\rangle \cdot \frac{1}{4\pi r^2},
$$

where the angle brackets denote a thermal or time average. In words: $g(r)$ is the probability of finding an atom at distance $r$ from a reference atom, normalised so that $g(r) \to 1$ at large $r$ for a homogeneous system. The quantity $4\pi r^2 g(r) \rho \, \mathrm{d}r$ is the average number of atoms in a shell of thickness $\mathrm{d}r$ at radius $r$.

### What the peaks mean

The shape of $g(r)$ encodes the local environment:

- **Sharp peaks at well-defined distances** indicate a *crystalline* solid. The first peak corresponds to nearest-neighbour distance; subsequent peaks correspond to successive coordination shells.

- **Broadened peaks with structure that persists to several diameters** indicate a *dense liquid* or *amorphous solid*. The first peak is still sharp (atoms cannot interpenetrate), but the second and third peaks are washed out by thermal motion.

- **A single broad peak followed by $g(r) \approx 1$** indicates a *dilute gas* with weak interactions.

- **A peak at very short distances** below the nearest-neighbour shell would indicate covalent bonding or, more commonly, a problem with the simulation (unphysical overlap).

The integral of $g(r)$ up to the first minimum, weighted by $4\pi r^2 \rho$, gives the *coordination number* — the average number of nearest neighbours. For FCC at zero temperature this is 12; for BCC, 8; for diamond, 4.

### Liquid versus crystal: a worked example

**Figure 2.3.3.** Imagine two RDF curves on the same axes. The horizontal axis is $r$ in ångströms from 0 to 10, the vertical axis $g(r)$ from 0 to 4. The crystalline curve has sharp Gaussian-like peaks at $r \approx 2.5$, $3.5$, $4.3$, $5.0$, $\ldots$ Å, each separated by regions of near-zero $g(r)$, corresponding to coordination shells of the FCC lattice. The liquid curve has a first peak at nearly the same position ($\sim 2.5$ Å) but broader, a noticeable trough at $\sim 3$ Å, a second peak at $\sim 5$ Å that is much weaker and broader, and oscillations that decay to $g(r) = 1$ by $r \sim 8$ Å.

The interpretation: the liquid retains short-range order (the nearest-neighbour shell is well defined) but loses long-range order (subsequent shells smear out). This is one of the cleanest experimental signatures of melting — observable in neutron and X-ray scattering, which measure the structure factor $S(q)$, the Fourier transform of $g(r) - 1$.

### A toy RDF in Python

We can compute $g(r)$ from atomic positions in a few lines.

```python
import numpy as np

def compute_rdf(
    positions: np.ndarray,
    box: float,
    n_bins: int = 200,
    r_max: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the radial distribution function in a cubic periodic box."""
    n = positions.shape[0]
    if r_max is None:
        r_max = 0.5 * box
    edges = np.linspace(0.0, r_max, n_bins + 1)
    counts = np.zeros(n_bins, dtype=np.int64)
    for i in range(n - 1):
        d = positions[i + 1 :] - positions[i]
        d -= box * np.round(d / box)  # minimum-image convention
        r = np.linalg.norm(d, axis=1)
        r = r[r < r_max]
        h, _ = np.histogram(r, bins=edges)
        counts += h
    centres = 0.5 * (edges[1:] + edges[:-1])
    shell_volumes = 4.0 * np.pi * centres ** 2 * (edges[1] - edges[0])
    rho = n / box ** 3
    # factor of 2 because each pair was counted once
    g = 2.0 * counts / (n * rho * shell_volumes)
    return centres, g
```

This routine accepts a set of positions in a cubic periodic box and returns the RDF. To test it, generate either (a) atoms on a perfect FCC lattice — you will see sharp shells; or (b) atoms placed at random with a hard-sphere exclusion — you will see a single broad first peak. Real production codes (MDAnalysis, OVITO, ASE) provide more efficient implementations using neighbour lists, but the conceptual content is the snippet above.

!!! note "The structure factor"
    What experiment actually measures is the structure factor $S(q)$, related to $g(r)$ by
    $$
    S(q) = 1 + \rho \int [g(r) - 1] \, e^{-i \mathbf{q} \cdot \mathbf{r}} \, \mathrm{d}^3 \mathbf{r}.
    $$
    For an isotropic system this reduces to a one-dimensional Fourier transform. Both $g(r)$ and $S(q)$ are useful, and either can be computed from the other; conventional practice is to compute $g(r)$ in MD analysis and to compare to $S(q)$ from scattering experiments.

## Bringing it together

These three diagrams partition the kind of information a materials simulation typically produces. The phase diagram is the thermodynamic answer to *what is stable, and where?*. The band structure is the electronic answer to *what does the material do with charge?*. The radial distribution function is the structural answer to *how are the atoms arranged?*. A complete characterisation of a new material typically presents all three.

The next section, on the software ecosystem, points you to the codes that produce each of these outputs in practice. After that, the exercises will check that you can interpret real diagrams as well as toy ones.
