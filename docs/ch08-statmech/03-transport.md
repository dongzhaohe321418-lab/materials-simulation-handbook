# 8.3 Transport Coefficients

Transport coefficients — diffusion $D$, shear viscosity $\eta$, thermal conductivity $\kappa$, ionic conductivity $\sigma_q$ — describe a system's **response** to a gradient. Diffusion describes mass flux under a concentration gradient ($\mathbf{J} = -D \nabla c$, Fick's first law); viscosity, momentum flux under a velocity gradient; thermal conductivity, heat flux under a temperature gradient.

These are inherently **non-equilibrium** properties. One way to compute them is **non-equilibrium MD** (NEMD): impose a gradient externally, measure the resulting current, take the ratio. We will mention this for completeness but focus on the alternative — extracting transport coefficients from equilibrium MD via the Green-Kubo formalism. Equilibrium fluctuations of currents, it turns out, carry exactly the information about how those currents respond to a small perturbation.

## Green-Kubo for diffusion

For a tagged particle undergoing equilibrium diffusion, the diffusion coefficient is related to the velocity autocorrelation function (VACF) by

$$
D = \frac{1}{3} \int_0^\infty \langle \mathbf{v}(0)\cdot \mathbf{v}(t)\rangle\, dt.
\tag{8.30}
$$

This is the **Green-Kubo** relation for self-diffusion. The factor 3 is the dimensionality (one factor for each Cartesian direction); some references absorb it differently.

**Derivation.** Start from the displacement,

$$
\Delta\mathbf{r}(t) = \int_0^t \mathbf{v}(t')\, dt'.
$$

The mean squared displacement is

$$
\langle |\Delta\mathbf{r}(t)|^2\rangle = \int_0^t \int_0^t \langle \mathbf{v}(t_1)\cdot \mathbf{v}(t_2)\rangle\, dt_1\, dt_2.
\tag{8.31}
$$

Change variables to $\tau = t_2 - t_1$ and $T = (t_1 + t_2)/2$. In equilibrium $\langle \mathbf{v}(t_1)\cdot\mathbf{v}(t_2)\rangle$ depends only on $\tau$, so

$$
\langle |\Delta\mathbf{r}(t)|^2\rangle = 2\int_0^t (t - \tau)\, \langle \mathbf{v}(0)\cdot\mathbf{v}(\tau)\rangle\, d\tau.
\tag{8.32}
$$

For large $t$ (longer than the VACF decay time), the $-\tau$ contribution converges to a constant, and

$$
\langle |\Delta\mathbf{r}(t)|^2\rangle \xrightarrow{t \to \infty} 2t \int_0^\infty \langle \mathbf{v}(0)\cdot\mathbf{v}(\tau)\rangle\, d\tau.
\tag{8.33}
$$

Comparing with the Einstein relation $\langle |\Delta\mathbf{r}(t)|^2\rangle = 6 D t$ (in 3D) gives (8.30). The Green-Kubo and Einstein formulations are equivalent — they are different ways of looking at the same long-time linear growth of the MSD. The Einstein form is integrated over time, the Green-Kubo form is differentiated.

## Fluctuation-dissipation theorem

The Green-Kubo relations are instances of the **fluctuation-dissipation theorem**, which states that the linear response of a system to a small external perturbation is determined by equilibrium fluctuations of the corresponding observable in the unperturbed system.

In its most general form, for an observable $A$ coupled to a perturbing field $f(t)$ via $H = H_0 - f(t) A$, the response

$$
\langle \Delta A(t)\rangle = \int_{-\infty}^t \chi(t - t')\, f(t')\, dt'
$$

is determined by the equilibrium time-correlation $\langle A(0) A(t)\rangle_\mathrm{eq}$ through

$$
\chi(t) = -\beta\, \frac{d}{dt} \langle A(0) A(t)\rangle\, \Theta(t).
\tag{8.34}
$$

For mechanical transport coefficients, the Green-Kubo formulas are explicit applications of (8.34) with the relevant flux as $A$.

The deeper statement: equilibrium fluctuations and non-equilibrium dissipation are two faces of the same microscopic dynamics. The Langevin thermostat is a direct manifestation — its friction and noise satisfy (7.40), the FDT.

## Viscosity

Shear viscosity is the response of momentum flux $\sigma_{xy}$ to a shear gradient $\partial v_x/\partial y$:

$$
\sigma_{xy} = -\eta\, \frac{\partial v_x}{\partial y}.
$$

The Green-Kubo formula:

$$
\eta = \frac{V}{k_B T}\int_0^\infty \langle \sigma_{xy}(0)\, \sigma_{xy}(t)\rangle\, dt,
\tag{8.35}
$$

where $\sigma_{xy}$ is an off-diagonal component of the **stress tensor**:

$$
V \sigma_{xy} = \sum_i m_i v_{ix} v_{iy} + \tfrac{1}{2} \sum_{i \ne j} r_{ij,x}\, F_{ij,y}.
\tag{8.36}
$$

The first term is the kinetic contribution (momentum transport by atoms moving), the second the virial contribution (momentum transport by forces). LAMMPS computes the full stress tensor via `compute stress/atom` and `compute pressure`; the off-diagonal components are what you correlate.

There are three independent off-diagonal stress components ($\sigma_{xy}$, $\sigma_{xz}$, $\sigma_{yz}$), and isotropic averaging combines them for better statistics. The variance of stress autocorrelations is much larger than for VACF — viscosity from Green-Kubo is notoriously noisy. For accurate results, expect to run several nanoseconds and possibly stitch in multiple independent replicas.

## Thermal conductivity

Heat flux $\mathbf{J}_Q$ is defined as the rate of energy transport per unit area. For a single-component system,

$$
\mathbf{J}_Q = \frac{1}{V}\left[\sum_i e_i\, \mathbf{v}_i + \tfrac{1}{2} \sum_{i \ne j} (\mathbf{F}_{ij}\cdot \mathbf{v}_i)\, \mathbf{r}_{ij}\right],
\tag{8.37}
$$

with $e_i$ a per-atom energy (kinetic plus half the pair potential energies of bonds involving $i$). The Green-Kubo formula:

$$
\kappa = \frac{V}{k_B T^2}\int_0^\infty \langle J_{Q,x}(0)\, J_{Q,x}(t)\rangle\, dt.
\tag{8.38}
$$

In LAMMPS, `compute heat/flux` accumulates the relevant quantity per step; you then time-correlate it in post-processing.

Thermal conductivity from Green-Kubo is the most demanding of the transport coefficients: the heat-flux autocorrelation has long tails, especially in solids where phonon transport is correlated over picoseconds to nanoseconds. A 10 ns simulation is a minimum; metals and dielectrics with long-lived phonons may need much longer.

## Convergence diagnostics

A Green-Kubo integral is computed by truncating the upper limit at some finite time $t_\mathrm{max}$:

$$
D(t_\mathrm{max}) = \frac{1}{3} \int_0^{t_\mathrm{max}} \langle \mathbf{v}(0)\cdot\mathbf{v}(t)\rangle\, dt.
\tag{8.39}
$$

Plot $D(t_\mathrm{max})$ as a function of $t_\mathrm{max}$; a converged value shows a plateau. If $D$ continues to drift, the simulation is too short. If $D$ oscillates wildly around a mean, the simulation has insufficient statistics — average over more replicas.

In practice, a common signature of trouble is **long-tail contamination**. The VACF for liquid argon falls off rapidly (sub-picosecond timescale), but contributions to the integral at picosecond scales remain non-negligible due to hydrodynamic backflow — the famous Alder-Wainwright $t^{-3/2}$ tail. Truncating the integral too early misses this tail; integrating too far drowns the signal in noise. Best practice: fit the tail to its expected functional form ($t^{-3/2}$ for fluid self-diffusion; exponential for many viscosity contributions) and extend the integral analytically beyond the noisy region.

```python
import numpy as np

def green_kubo_diffusion(vacf: np.ndarray, dt: float) -> np.ndarray:
    """Running Green-Kubo integral D(t_max) = (1/3) ∫_0^{t_max} <v(0).v(t)> dt.

    Parameters
    ----------
    vacf : (n_steps,) array of <v(0).v(t)> (already summed over species).
    dt : timestep between samples.

    Returns
    -------
    (n_steps,) array of D as function of t_max.
    """
    return np.cumsum(vacf) * dt / 3.0
```

Use this to plot $D(t_\mathrm{max})$. A clean simulation produces a curve that rises rapidly, then plateaus.

## Einstein vs Green-Kubo

For diffusion, the two formulations are equivalent in the long-time limit. Which to use?

| Diagnostic | Einstein (MSD) | Green-Kubo (VACF) |
|---|---|---|
| Unwrapping required | Yes (always use unwrapped) | No |
| Sensitivity to ballistic regime | Discard short times | Automatic |
| Tail contribution explicit | Hidden in MSD slope | Explicit — must be integrated |
| Statistical noise | Moderate | Higher |
| Implementation | One line | FFT autocorrelation |

For self-diffusion in a liquid, Einstein is more robust; use Green-Kubo as a cross-check. For viscosity and thermal conductivity, only Green-Kubo (or NEMD) is practical because there is no natural "Einstein" formulation.

## Non-equilibrium MD as alternative

Direct measurement: impose a gradient, measure the response, take the ratio.

For viscosity, the **Müller-Plathe method** (1997) periodically swaps the momentum of an atom in a high-velocity slab with one in a low-velocity slab, generating a known momentum flux. Measure the resulting velocity profile; viscosity is the slope ratio. Implemented in LAMMPS as `fix viscosity`.

For thermal conductivity, the analogous Müller-Plathe protocol swaps kinetic energy between hot and cold slabs (`fix heat`). The known heat flux divided by the measured temperature gradient gives $\kappa$.

NEMD pros: direct, no autocorrelation integration, often converges faster than Green-Kubo for nasty cases (solids, polymers).

NEMD cons: results may depend on the magnitude of the imposed gradient (need to extrapolate to zero gradient for true linear-response limit); equilibrium ensemble interpretation is less direct.

Modern practice mixes them: run Green-Kubo for a primary estimate, NEMD as cross-validation. If they disagree by more than a factor of 2, both are wrong.

## Ionic conductivity

For electrolytes, the ionic conductivity at zero external field is

$$
\sigma_q = \frac{1}{3 V k_B T} \int_0^\infty \langle \mathbf{J}(0)\cdot \mathbf{J}(t)\rangle\, dt,
\quad
\mathbf{J}(t) = \sum_i q_i\, \mathbf{v}_i(t).
\tag{8.40}
$$

This is the Nernst-Einstein generalisation including all ion-ion correlations. The simple Nernst-Einstein approximation (ignoring cross-correlations) gives

$$
\sigma_q^\mathrm{NE} = \frac{1}{k_B T} \sum_\alpha c_\alpha q_\alpha^2 D_\alpha,
$$

valid when ion-ion correlations are negligible (dilute electrolytes). Concentrated electrolytes need the full (8.40), and the ratio $\sigma_q / \sigma_q^\mathrm{NE}$ — sometimes called the Haven ratio — is a diagnostic of how correlated the ionic motion is.

## What we have

A toolbox for extracting kinetics — diffusion, viscosity, thermal conductivity, ionic conductivity — from equilibrium MD. The same fluctuation-dissipation framework also underlies any linear-response calculation; the dielectric function, magnetic susceptibility, and optical conductivity all derive from analogous correlation functions of the relevant flux.

The chapter closes with [§8.4](04-phase-diagrams.md), which puts the free-energy methods to work for the most demanding application: phase diagrams from atomistic simulation.
