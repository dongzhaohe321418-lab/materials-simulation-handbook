# 3b.4 — The Free Electron Gas and the Sommerfeld Expansion

<figure markdown>
![Free-electron density of states](../assets/figures/ch03b/fig_free_electron_dos.png){ width="600" }
<figcaption>Figure 3b.4.1. The free-electron density of states is \(g(\varepsilon) \propto \sqrt{\varepsilon}\). At \(T = 0\), all states below the Fermi energy \(\varepsilon_F\) are filled (shaded region); for copper this corresponds to \(\varepsilon_F \approx 7\) eV.</figcaption>
</figure>

> *"Take all the electrons in a metal, ignore the ions completely, ignore the interactions completely, and see what you get. The answer is: most of the metal."* — Drude, paraphrased

The simplest band structure is no band structure: ignore the lattice altogether and treat the conduction electrons as a gas of non-interacting fermions in a box. This is the *free electron gas* model, and despite its crudeness it accounts for many properties of simple metals to within tens of percent: electronic specific heat, Pauli paramagnetism, the linear-in-$T$ resistivity coefficient, and the rough magnitude of the Fermi energy. It is also the starting point — quite literally — for the local density approximation in DFT, which uses the exchange–correlation energy of a *uniform* electron gas as its building block.

In this section we develop the free electron gas in three dimensions, derive its density of states, define the Fermi energy and Fermi temperature, and use the Sommerfeld expansion to compute finite-temperature corrections including the famous linear electronic specific heat. We finish with a numerical worked example for copper.

## 3b.4.1 The setup

$N$ non-interacting electrons of mass $m$ live in a cubic box of side $L$, volume $V = L^3$. Impose periodic boundary conditions (the same Born–von Kármán device as in §3b.1). The single-particle eigenstates are plane waves

$$\psi_\mathbf k(\mathbf r) = \frac{1}{\sqrt V}\, e^{i\mathbf k\cdot\mathbf r}, \qquad \mathbf k = \frac{2\pi}{L}(n_x, n_y, n_z), \quad n_i \in \mathbb Z, \tag{3b.4.1}$$

with energy

$$\varepsilon(\mathbf k) = \frac{\hbar^2 |\mathbf k|^2}{2m}. \tag{3b.4.2}$$

Electrons carry spin $1/2$, so each $\mathbf k$ accommodates two electrons. The allowed $\mathbf k$ form a cubic lattice in reciprocal space with spacing $2\pi/L$, hence density

$$\frac{V}{(2\pi)^3}\, d^3 k \tag{3b.4.3}$$

states per unit volume of $\mathbf k$-space (counting spin doubles this).

## 3b.4.2 Filling the Fermi sphere

At zero temperature, the ground state has every state with $\varepsilon < \varepsilon_F$ filled and every state with $\varepsilon > \varepsilon_F$ empty. Because $\varepsilon$ depends only on $|\mathbf k|$, the occupied region is a sphere of radius $k_F$ — the **Fermi sphere**. The Fermi wavevector $k_F$ is fixed by the total number of electrons:

$$N = 2 \cdot \frac{V}{(2\pi)^3} \cdot \frac{4\pi}{3}\, k_F^3, \tag{3b.4.4}$$

where the leading $2$ is the spin degeneracy. Solve for $k_F$:

$$\boxed{\; k_F = (3\pi^2 n)^{1/3}, \quad n := N/V. \;} \tag{3b.4.5}$$

The Fermi energy is

$$\varepsilon_F = \frac{\hbar^2 k_F^2}{2m} = \frac{\hbar^2}{2m}(3\pi^2 n)^{2/3}, \tag{3b.4.6}$$

and the Fermi temperature is

$$T_F := \varepsilon_F / k_B. \tag{3b.4.7}$$

For metals, $T_F$ is typically $10^4 - 10^5$ K — far above room temperature. This is the central observation justifying the *degenerate* electron gas approximation: at any laboratory temperature $T \ll T_F$, the electron gas is essentially in its ground state, and finite-temperature effects show up as small corrections in powers of $T/T_F$.

## 3b.4.3 The density of states

Many thermodynamic quantities depend on the electron distribution only through the *density of states* $g(\varepsilon)$, defined so that $g(\varepsilon) d\varepsilon$ is the number of single-particle states (per unit volume) with energy in $[\varepsilon, \varepsilon + d\varepsilon]$.

Start with the number of states with energy $\le \varepsilon$:

$$N(\varepsilon)/V = 2 \cdot \frac{1}{(2\pi)^3} \cdot \frac{4\pi}{3} k(\varepsilon)^3, \qquad k(\varepsilon) = \sqrt{2m\varepsilon}/\hbar. \tag{3b.4.8}$$

Substituting,

$$\frac{N(\varepsilon)}{V} = \frac{1}{3\pi^2}\left(\frac{2m\varepsilon}{\hbar^2}\right)^{3/2}. \tag{3b.4.9}$$

Differentiating with respect to $\varepsilon$,

$$\boxed{\; g(\varepsilon) = \frac{1}{2\pi^2}\left(\frac{2m}{\hbar^2}\right)^{3/2}\sqrt{\varepsilon}. \;} \tag{3b.4.10}$$

The 3D free-electron DOS grows as $\sqrt{\varepsilon}$. A useful equivalent form, expressing $g$ in terms of $\varepsilon_F$ and $n$:

$$g(\varepsilon_F) = \frac{3n}{2\varepsilon_F}, \qquad g(\varepsilon) = g(\varepsilon_F)\, \sqrt{\varepsilon/\varepsilon_F}. \tag{3b.4.11}$$

A factor that recurs everywhere. In 1D the DOS is $\propto 1/\sqrt\varepsilon$ (diverging at the band bottom); in 2D it is constant; in 3D it is $\propto\sqrt\varepsilon$. The dimensional dependence is a useful sanity check for any DOS you ever compute.

## 3b.4.4 Numerical example: copper

Copper has one $4s$ electron per atom outside the closed $3d^{10}$ shell. The metallic density is $8.96$ g/cm$^3$ with atomic mass $63.55$ g/mol, giving an atomic number density of $8.49\times 10^{28}$ atoms/m$^3$, and hence (one electron per atom) $n = 8.49\times 10^{28}$ m$^{-3}$. Plugging into (3b.4.5) and (3b.4.6):

$$k_F = (3\pi^2 \cdot 8.49\times 10^{28})^{1/3} \approx 1.36\times 10^{10}\text{ m}^{-1}, \tag{3b.4.12}$$

$$\varepsilon_F = \frac{(1.055\times 10^{-34})^2}{2 \cdot 9.109\times 10^{-31}}\cdot (1.36\times 10^{10})^2 \approx 1.13\times 10^{-18}\text{ J} \approx 7.04\text{ eV}, \tag{3b.4.13}$$

$$T_F = \varepsilon_F/k_B \approx 8.16\times 10^4\text{ K}. \tag{3b.4.14}$$

The Fermi temperature of copper is about 82 000 K. At room temperature $T/T_F \approx 0.0037$, vanishingly small. The famous experimental Fermi energy of copper, measured from photoemission, is 7.0 eV — the free electron gas is essentially exact, and this is the reason copper is the textbook simple metal.

## 3b.4.5 Finite-temperature: the Fermi–Dirac distribution

At temperature $T$ the occupation of single-particle state with energy $\varepsilon$ is governed by Fermi–Dirac statistics:

$$f(\varepsilon) = \frac{1}{e^{(\varepsilon - \mu)/k_B T} + 1}, \tag{3b.4.15}$$

with chemical potential $\mu$ determined by particle conservation:

$$n = \int_0^\infty g(\varepsilon)\, f(\varepsilon)\, d\varepsilon. \tag{3b.4.16}$$

At $T = 0$, $f$ is a step function: $f(\varepsilon) = 1$ for $\varepsilon < \mu = \varepsilon_F$ and zero otherwise. At low $T$, the step is *smeared* over an energy window of width $\sim k_B T$ around $\mu$. Only electrons near the Fermi level participate in thermal processes — the deep interior of the Fermi sphere is frozen. This is the central feature of the degenerate electron gas, and it controls every "low-T anomaly" in metals.

## 3b.4.6 The Sommerfeld expansion

To compute the average energy and the specific heat we need integrals of the form

$$I[H] := \int_{-\infty}^\infty H(\varepsilon)\, f(\varepsilon)\, d\varepsilon. \tag{3b.4.17}$$

for slowly varying $H$. The Sommerfeld expansion is a systematic asymptotic series in $k_B T/\mu$. To derive it, integrate by parts using $\int_{-\infty}^\infty H f\, d\varepsilon = -\int_{-\infty}^\infty K(\varepsilon)\, f'(\varepsilon)\, d\varepsilon$ where $K(\varepsilon) = \int_{-\infty}^\varepsilon H(\varepsilon')\, d\varepsilon'$. The derivative $-f'(\varepsilon)$ is sharply peaked at $\mu$ with width $\sim k_B T$. Taylor expand $K(\varepsilon)$ about $\mu$:

$$K(\varepsilon) = K(\mu) + K'(\mu)(\varepsilon - \mu) + \tfrac12 K''(\mu)(\varepsilon - \mu)^2 + \cdots \tag{3b.4.18}$$

Substitute and use $\int(-f') d\varepsilon = 1$, $\int(\varepsilon-\mu)(-f') d\varepsilon = 0$ (the integrand is odd), $\int(\varepsilon-\mu)^2(-f') d\varepsilon = \pi^2(k_B T)^2/3$. The result is

$$\boxed{\; I[H] = \int_{-\infty}^\mu H(\varepsilon)\, d\varepsilon + \frac{\pi^2}{6}(k_B T)^2\, H'(\mu) + O(T^4). \;} \tag{3b.4.19}$$

This is the Sommerfeld expansion. The leading correction is $O(T^2)$, with the coefficient $\pi^2/6$ — a number that appears in every metal property at low temperature.

## 3b.4.7 Specific heat

Apply (3b.4.19) to the energy density:

$$u(T) = \int_0^\infty \varepsilon\, g(\varepsilon)\, f(\varepsilon)\, d\varepsilon. \tag{3b.4.20}$$

Identify $H(\varepsilon) = \varepsilon g(\varepsilon)$. Then $H'(\mu) = g(\mu) + \mu g'(\mu)$. To leading order, $\mu \approx \varepsilon_F$ (corrections are $O(T^2)$ and matter only at next order). Also, the $T=0$ integral $\int_0^{\varepsilon_F} \varepsilon g(\varepsilon)\, d\varepsilon = u(0)$, the zero-temperature energy density. So

$$u(T) = u(0) + \frac{\pi^2}{6}(k_B T)^2\left[g(\varepsilon_F) + \varepsilon_F g'(\varepsilon_F)\right] + O(T^4). \tag{3b.4.21}$$

Differentiating with respect to $T$ to get the specific heat per unit volume,

$$c_v^\text{el}(T) = \frac{\partial u}{\partial T} = \frac{\pi^2}{3} k_B^2 T\left[g(\varepsilon_F) + \varepsilon_F g'(\varepsilon_F)\right]. \tag{3b.4.22}$$

To leading order in $T/T_F$, the second term in the bracket is comparable to the first; but for the free electron gas it is conventional to write the answer as

$$\boxed{\; c_v^\text{el}(T) = \frac{\pi^2}{3} g(\varepsilon_F)\, k_B^2\, T. \;} \tag{3b.4.23}$$

The electronic specific heat is *linear* in $T$, with a coefficient $\gamma := (\pi^2/3) g(\varepsilon_F) k_B^2$ called the **Sommerfeld coefficient**. The lattice (phonon) specific heat is $T^3$ at low $T$ (next section), so at sufficiently low temperature the electronic contribution dominates. Plotting $c_v/T$ vs $T^2$ — a so-called "Sommerfeld plot" — gives a straight line whose intercept is $\gamma$ and slope is the lattice $T^3$ coefficient. This is the experimental standard for measuring $g(\varepsilon_F)$.

The free electron prediction for $\gamma$ in copper, using $g(\varepsilon_F) = 3n/(2\varepsilon_F)$ from (3b.4.11):

$$\gamma_\text{free} = \frac{\pi^2}{3}\cdot \frac{3n}{2\varepsilon_F}\cdot k_B^2 = \frac{\pi^2 n k_B^2}{2\varepsilon_F} \approx 5.0\times 10^{-4}\text{ J K}^{-2}\text{ mol}^{-1}. \tag{3b.4.24}$$

(Per-volume this is $\approx 71$ J K$^{-2}$ m$^{-3}$; multiplication by Cu's molar volume $V_m \approx 7.09\times 10^{-6}$ m$^3$/mol gives the per-mole figure.) The measured value is $7.0\times 10^{-4}$ J K$^{-2}$ mol$^{-1}$. The ratio 1.4 is the *effective mass enhancement* — band structure and electron–phonon coupling combine to make the real Cu electrons slightly heavier than free electrons.

## 3b.4.8 Python: $g(\varepsilon)$, $\varepsilon_F$ for copper

The following script computes $g(\varepsilon)$ from (3b.4.10), inverts (3b.4.16) numerically at $T=0$ to find $\varepsilon_F$ for copper's electron density, and verifies (3b.4.6) directly.

```python
"""Free-electron density of states and Fermi energy for copper."""
from __future__ import annotations
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# Physical constants (SI)
HBAR: float = 1.054_571_817e-34   # J s
M_E: float = 9.109_383_7015e-31   # kg
E_CHARGE: float = 1.602_176_634e-19  # J/eV (for unit conversion)

# Copper number density (one electron per atom)
N_CU: float = 8.49e28  # electrons / m^3

def dos(eps: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """Free-electron DOS per unit volume (states / J / m^3)."""
    pref: float = (1.0 / (2.0 * np.pi**2)) * (2.0 * M_E / HBAR**2) ** 1.5
    return pref * np.sqrt(np.maximum(eps, 0.0))

def electron_density(eps_f: float) -> float:
    """Integrate the DOS from 0 to eps_f (T=0)."""
    eps_grid: npt.NDArray[np.float64] = np.linspace(0.0, eps_f, 10_000)
    return float(np.trapezoid(dos(eps_grid), eps_grid))

def find_fermi_energy(n: float) -> float:
    """Solve electron_density(eps_F) = n for eps_F (Joules)."""
    eps_lo: float = 1e-22   # arbitrary small
    eps_hi: float = 1e-17   # ~ 60 eV; safely above any metal's eps_F
    return brentq(lambda x: electron_density(x) - n, eps_lo, eps_hi, xtol=1e-25)

def main() -> None:
    eps_f_J: float = find_fermi_energy(N_CU)
    eps_f_eV: float = eps_f_J / E_CHARGE
    eps_f_analytic_J: float = (HBAR**2 / (2 * M_E)) * (3 * np.pi**2 * N_CU) ** (2 / 3)
    print(f"eps_F (numerical)  = {eps_f_eV:.4f} eV")
    print(f"eps_F (analytic)   = {eps_f_analytic_J / E_CHARGE:.4f} eV")
    k_f: float = (3 * np.pi**2 * N_CU) ** (1 / 3)
    print(f"k_F                = {k_f:.4e} 1/m")
    print(f"T_F                = {eps_f_J / 1.380649e-23:.2f} K")

    # Plot g(eps)
    eps_plot: npt.NDArray[np.float64] = np.linspace(0.0, 1.5 * eps_f_J, 1000)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(eps_plot / E_CHARGE, dos(eps_plot) * E_CHARGE * 1e-28,
            color="C0", lw=2)
    ax.axvline(eps_f_eV, color="red", linestyle="--",
               label=f"$\\varepsilon_F$ = {eps_f_eV:.2f} eV")
    ax.set_xlabel(r"$\varepsilon$ (eV)")
    ax.set_ylabel(r"$g(\varepsilon)$ ($10^{28}$ states eV$^{-1}$ m$^{-3}$)")
    ax.set_title("Copper: free-electron density of states")
    ax.legend()
    plt.tight_layout()
    plt.savefig("cu_dos.pdf")
    plt.show()

if __name__ == "__main__":
    main()
```

Running this script prints, among other things, $\varepsilon_F \approx 7.04$ eV and $T_F \approx 8.17\times 10^4$ K — matching the by-hand calculation in §3b.4.4 to four significant figures.

## 3b.4.9 The jellium model and the bridge to DFT

If we generalise the free electron gas by adding a uniform positive background of charge density $+en$ to make the system neutral, the resulting system is called **jellium**. The Hamiltonian now contains the Coulomb interaction between electrons. Despite the simplicity of the model, the *exact* exchange–correlation energy per electron is well known as a function of density: it is a smooth, monotonically increasing function $\epsilon_\text{xc}(n)$, computed accurately by Ceperley–Alder Monte Carlo and parametrised by Perdew, Zunger, and others.

The **local density approximation** (LDA) to the exchange–correlation functional in DFT consists of asserting that, *locally*, the exchange–correlation energy density at a point $\mathbf r$ in an *inhomogeneous* electron system is the same as that of a uniform electron gas of the same density $n(\mathbf r)$:

$$E_\text{xc}^\text{LDA}[\rho] = \int \rho(\mathbf r)\, \epsilon_\text{xc}(\rho(\mathbf r))\, d^3 r. \tag{3b.4.25}$$

This is by far the most consequential idea in computational materials science. It is also the reason that every DFT code on Earth uses the free electron gas as its reference, and the reason that LDA does so well for free-electron-like metals (Na, Al, Cu) and so badly for strongly correlated materials.

In Chapter 5 you will see (3b.4.25) written down formally, and in Chapter 6 you will choose a Perdew–Zunger LDA functional in your input file. When that happens, recall: it is *literally* the energy of a Fermi gas, applied pointwise.

!!! warning "When jellium fails"
    LDA inherits the strengths and weaknesses of jellium. It is excellent for slowly varying densities (metals). It is moderate for typical solids. It is poor for systems with strong density variations (van der Waals interactions, atomic densities near nuclei) and disastrous for localised electrons (NiO, MnO, CoO, rare-earth compounds). Each of these is a forward reference to a specific functional improvement in Chapter 5: GGAs (PBE), meta-GGAs (SCAN), hybrids (HSE06), DFT+U.

## Where this is used later

- **Tier 1.** §5.2 (the Kohn–Sham equations specialise to jellium when $V_\text{ext}$ is constant), §5.5 (LDA construction), §6.3 (smearing schemes — Methfessel–Paxton, Fermi–Dirac — all rely on the Sommerfeld picture of states near $\varepsilon_F$), §6.4 (interpreting metallic DOS at $\varepsilon_F$ and computing $\gamma$).
- **Tier 2.** §8.3 (electronic contribution to the Helmholtz free energy at finite temperature, dominated by the Sommerfeld term), §9.6 (transferability of MLIPs to metals: a learnt potential must reproduce the *electronic* equation of state implicit in the free electron gas).
- **Capstone Project 1.** Validate a screened-hybrid DFT calculation against the free-electron-gas prediction for a simple metal — a sanity check that any pre-screening study must pass.

Onwards to phonons (§3b.5), where we leave the electrons behind and quantise the *nuclear* degrees of freedom.
