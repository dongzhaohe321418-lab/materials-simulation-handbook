# 8.5 Exercises

Seven exercises with worked solutions. Difficulty markers:

- (★) Routine application.
- (★★) Standard.
- (★★★) Research-flavoured.

---

## Exercise 8.1 (★) — Heat capacity from energy variance

Run an NVT LAMMPS simulation of 4000 LJ argon atoms at $T^* = 1.5, \rho^* = 0.8$ for 200 ps after equilibration. Extract the total energy time series and compute the heat capacity via the fluctuation formula

$$
C_V = \frac{\langle E^2\rangle - \langle E\rangle^2}{k_B T^2}.
$$

Compare to the analytical ideal-gas value $C_V^\mathrm{id} = \tfrac{3}{2} N k_B$ and explain the excess.

**Solution.** Read `step pe ke etotal` from the LAMMPS log. The kinetic energy alone (in NVT with Nosé-Hoover) gives the ideal-gas part: $\langle (KE - \langle KE\rangle)^2\rangle = \tfrac{3}{2}N(k_B T)^2$, so $C_V^\mathrm{kin} = \tfrac{3}{2}N k_B$.

The potential energy contributes the **excess** heat capacity. For LJ at this state point, the potential energy variance gives $C_V^\mathrm{pot}/(N k_B) \approx 1$–2, so the total $C_V/(N k_B) \approx 3$. This matches Verlet's 1968 simulation values and modern reference data for LJ at this density.

```python
import numpy as np
data = np.loadtxt("log.dat", skiprows=2)
step, pe, ke, etot = data.T
n_block = 4000
# Use block averaging to estimate the error
blocks = np.array_split(etot, n_block)
means = np.array([b.mean() for b in blocks])
variances = np.array([b.var() for b in blocks])
Cv = variances.mean() / (kB * T**2)
```

The block-averaged variance estimator (Flyvbjerg-Petersen) gives proper error bars; a single-block variance underestimates the error for autocorrelated data.

---

## Exercise 8.2 (★) — Compressibility from volume variance

In NPT at $T^* = 1.5, P^* = 1.0$, compute the isothermal compressibility $\kappa_T$ from the volume variance. Compare to the experimental compressibility of liquid argon at the corresponding $(T, P)$ in SI units.

**Solution.** $\kappa_T = (\langle V^2\rangle - \langle V\rangle^2)/(k_B T \langle V\rangle)$.

Run NPT with `fix npt iso 1.0 1.0 1.0`, log `vol` every 100 steps for 200 ps. Compute the variance.

In reduced units, $\kappa_T^* = \kappa_T \epsilon/\sigma^3$ converts via $\kappa_T^\mathrm{SI} = \kappa_T^* \sigma^3/\epsilon$. With argon's $\sigma = 3.405$ Å and $\epsilon/k_B = 120$ K, expect $\kappa_T \sim 10^{-9}$ Pa$^{-1}$, agreeing with experimental liquid argon compressibility at the relevant $(T, P)$.

---

## Exercise 8.3 (★★) — FEP for an alchemical swap

Compute the free energy difference for transmuting one argon atom into a "lighter argon" with $\epsilon_2 = 0.8 \epsilon$, $\sigma_2 = \sigma$, in a bath of 999 standard argon atoms at $T^* = 1.5, \rho^* = 0.8$. Use forward and reverse FEP and report both estimates with errors.

**Solution.** Run two simulations:

1. Reference state: 999 standard + 1 standard. At every step compute $\Delta U = U_\mathrm{perturbed} - U_\mathrm{ref}$, the energy change if atom 1's $\epsilon$ were 0.8. Estimate $\Delta A_\mathrm{fwd} = -k_B T \ln \langle e^{-\beta \Delta U}\rangle_\mathrm{ref}$.
2. Perturbed state: 999 standard + 1 light. At every step compute $-\Delta U$ relative to reverting atom 1 back. Estimate $\Delta A_\mathrm{rev}$ analogously; under the Zwanzig identity $\Delta A_\mathrm{fwd} = -\Delta A_\mathrm{rev}$ ideally.

If forward and reverse disagree by more than the combined statistical error, FEP has not converged for this swap. For a single-atom $\epsilon$ change of 20%, $|\beta \Delta U| \sim 0.3$ — well within the FEP-friendly regime, so the two estimates should agree to within $\sim 0.01\, \epsilon$.

The BAR estimator (Bennett 1976) combines forward and reverse data optimally:

```python
import pymbar
delta_A_BAR, err = pymbar.bar(np.array([dU_fwd]), -np.array([dU_rev]))
```

For larger perturbations (e.g., swap argon for krypton with $\epsilon_2 \approx 1.5 \epsilon$), naive FEP fails and either multi-stage windows or BAR are required.

---

## Exercise 8.4 (★★) — Thermodynamic integration for an LJ liquid

Implement the TI for $\mu^\mathrm{ex}$ of LJ at $T^* = 1.5, \rho^* = 0.8$ described in §8.2. Use 11 $\lambda$ values; for each, run 50 ps equilibration and 200 ps production. Plot $\langle U\rangle_\lambda/N$ vs $\lambda$, integrate with Simpson's rule, and compare to the reference value $\mu^\mathrm{ex} \approx -1.60\, \epsilon$.

**Solution.** Follow the LAMMPS script in §8.2 with soft-core LJ for small $\lambda$. The plot of $\langle U\rangle_\lambda/N$ vs $\lambda$ should be a smooth, monotonic curve, rising from $0$ at $\lambda = 0$ to $\langle U\rangle_\mathrm{full}/N \approx -6\, \epsilon$ at $\lambda = 1$.

The integral $\int_0^1 \langle U\rangle_\lambda/N\, d\lambda$ should give $\mu^\mathrm{ex} \approx -1.6 \pm 0.02\, \epsilon$.

Pitfalls:
- Without soft-core, the small-$\lambda$ integrand has a $\lambda^{-3}$ singularity that makes Simpson's rule fail.
- Without enough $\lambda$ points (fewer than 7), Simpson's rule undersamples the curvature.
- Without enough equilibration at each $\lambda$, the integrand is biased.

Compare your value to the Widom particle-insertion estimate at the same state point as an independent cross-check.

---

## Exercise 8.5 (★★) — Green-Kubo diffusion in liquid argon

From an NVE trajectory of 4000 LJ argon atoms at $T^* = 1.5, \rho^* = 0.8$ (use the trajectory from Ex 7.6), compute the velocity autocorrelation function $C_{vv}(t)$ and the running Green-Kubo integral $D(t_\mathrm{max})$. Plot both. Identify the plateau in $D(t_\mathrm{max})$ and compare to the Einstein-relation value from Ex 7.6.

**Solution.** Use the FFT-based VACF computation from §7.6. The VACF shows a positive peak at $t = 0$, drops to negative around $\tau \approx 0.3 \tau$ (the caging time), recovers to small positive at $\tau \approx 1\, \tau$, and decays as $t^{-3/2}$ at long times.

```python
import numpy as np
# velocities: (n_frames, n_atoms, 3) unwrapped from trajectory
vacf = np.mean(np.sum(velocities[0:1] * velocities, axis=-1), axis=-1)
D_t = np.cumsum(vacf) * dt / 3.0
```

The running integral $D(t_\mathrm{max})$ rises rapidly, reaches a near-plateau around $t \sim 5\, \tau$, with a slow drift due to the $t^{-3/2}$ tail. Fitting the tail analytically beyond the noise threshold gives an estimate of the asymptotic $D$.

In LJ reduced units, $D \approx 0.04\, \sigma^2/\tau$ at this state point, consistent with the Einstein value from Ex 7.6 to within a few percent.

---

## Exercise 8.6 (★★★) — Two-phase melting of EAM aluminium

Locate the melting temperature of Mishin-Farkas EAM aluminium ($a_0 = 4.05$ Å) using the two-phase method, with a bracketing accuracy of $\pm 25$ K.

**Solution.** Use the LAMMPS sketch from §8.4 adapted to Al ($a_0 = 4.05$, mass 26.98). Iterate:

1. First run: $T = 900$ K. If solid grows, $T_m > 900$ K.
2. Run at 1000 K. If liquid grows, $T_m < 1000$ K. Bracket: $900 < T_m < 1000$.
3. Run at 950 K. Reduce bracket by half each step.

After 3–4 runs you have $T_m \approx 925 \pm 25$ K for Mishin Al, compared to experimental 933 K — about 1% accurate.

A diagnostic: compute $Q_6$ bond-orientational order parameter for each atom on the fly (LAMMPS `compute orientorder/atom`); plot the fraction of solid-like atoms vs time. At $T_m$ this fraction should be roughly constant.

Subtleties:
- Use `fix npt z 1.0 1.0 1.0` to allow only $z$ to fluctuate. Isotropic NPT distorts the crystal.
- Start with a clean interface: melt the top half by region thermostat to far above $T_m$, then equilibrate everything at the target $T$ before production.
- The interface migrates at typically 1–10 m/s, so a 1 ns simulation moves the interface by 1–10 nm — enough to see clearly in an Al box 4 nm long along $z$.

---

## Exercise 8.7 (★★★) — Viscosity from Green-Kubo

For LJ liquid at $T^* = 1.5, \rho^* = 0.8$, compute the shear viscosity from the Green-Kubo formula (8.35). Run a 5 ns NVT simulation, log the off-diagonal stress tensor every 4 fs, compute the autocorrelation function and its integral. Plot $\eta(t_\mathrm{max})$ and locate the plateau. Compare to the literature LJ viscosity for these conditions ($\eta^* \approx 2.5$ in reduced units).

**Solution.** This is the standard hard exercise. Sample considerations:

```text
compute         stress all pressure NULL virial
fix             sxy all ave/time 4 1 4 c_stress[4] c_stress[5] c_stress[6] file stress.dat
```

In post-processing:

```python
import numpy as np

stress_xy = np.loadtxt("stress.dat", usecols=(1, 2, 3))  # xy, xz, yz
# Use isotropic average over 3 off-diagonal components
acf_xy = autocorr(stress_xy[:, 0])
acf_xz = autocorr(stress_xy[:, 1])
acf_yz = autocorr(stress_xy[:, 2])
acf = (acf_xy + acf_xz + acf_yz) / 3.0
eta_t = np.cumsum(acf) * dt * V / (kB * T)
```

Plot `eta_t` vs `t`. Expect a noisy rise to a plateau around $\eta^* \approx 2.5$ that takes 2–5 $\tau$ to reach. The fluctuations on the plateau are large — viscosity from GK has high variance.

Improvements: (a) average over multiple independent replicas (5–10 short runs are better than one long one); (b) fit the tail to a stretched exponential and extrapolate; (c) cross-check with the Müller-Plathe NEMD method. If the value disagrees by more than 20% with the literature, the simulation is undersampled.

---

## Summary

| Exercise | Skill |
|---|---|
| 8.1 | Heat capacity from variance |
| 8.2 | Compressibility from volume variance |
| 8.3 | Free energy perturbation, forward and reverse |
| 8.4 | Thermodynamic integration over a coupling parameter |
| 8.5 | Green-Kubo diffusion from VACF |
| 8.6 | Two-phase coexistence to locate $T_m$ |
| 8.7 | Viscosity from stress autocorrelation |

These exercises together cover the major numerical workflows of statistical-mechanical analysis from MD. With them mastered, you are ready for [Chapter 9](../ch09-mlip/index.md), where the classical force fields of §7.4 are replaced by machine-learning potentials. The downstream machinery — integrators, thermostats, ensembles, free energy methods, Green-Kubo — applies identically; the only difference is that the forces now come from a neural network rather than from an EAM table.
