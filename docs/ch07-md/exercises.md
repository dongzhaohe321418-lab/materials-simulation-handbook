# 7.7 Exercises

Eight exercises, mixing theory and code. Difficulty markers:

- (★) Routine — straight application of the chapter.
- (★★) Standard — synthesise two ideas.
- (★★★) Research-flavoured — open-ended or requires extension.

Solutions inline.

---

## Exercise 7.1 (★) — Euler vs Verlet for a harmonic oscillator

Implement forward Euler and velocity-Verlet for the 1D harmonic oscillator with $m = k = 1$, $x_0 = 1$, $v_0 = 0$. Run each for $10^5$ steps at $\Delta t = 0.05$. Plot the relative energy drift $(E - E_0)/E_0$ for both integrators on the same axes.

**Solution.**

```python
import numpy as np
import matplotlib.pyplot as plt

def euler(x0, v0, dt, n_steps):
    x = np.zeros(n_steps + 1); v = np.zeros(n_steps + 1)
    x[0], v[0] = x0, v0
    for n in range(n_steps):
        a = -x[n]
        x[n+1] = x[n] + dt * v[n]
        v[n+1] = v[n] + dt * a
    return x, v

def vv(x0, v0, dt, n_steps):
    x = np.zeros(n_steps + 1); v = np.zeros(n_steps + 1)
    x[0], v[0] = x0, v0
    a = -x[0]
    for n in range(n_steps):
        v_half = v[n] + 0.5 * dt * a
        x[n+1] = x[n] + dt * v_half
        a_new = -x[n+1]
        v[n+1] = v_half + 0.5 * dt * a_new
        a = a_new
    return x, v

dt, N = 0.05, 100_000
xE, vE = euler(1.0, 0.0, dt, N)
xV, vV = vv(1.0, 0.0, dt, N)
E_euler = 0.5 * (vE**2 + xE**2); E_vv = 0.5 * (vV**2 + xV**2)

fig, ax = plt.subplots()
ax.semilogy(np.arange(N+1)*dt, np.abs(E_euler - E_euler[0])/E_euler[0],
            label="Euler")
ax.semilogy(np.arange(N+1)*dt, np.abs(E_vv - E_vv[0])/E_vv[0],
            label="Velocity-Verlet")
ax.set(xlabel="t", ylabel="|ΔE|/E_0")
ax.legend(); fig.savefig("ex_7_1.png", dpi=150)
```

Euler grows like $(1 + \omega^2 \Delta t^2)^n$ — exponential blow-up to numerical overflow within a few thousand periods. Velocity-Verlet oscillates with a bounded amplitude $\sim \Delta t^2 / 2$, independent of $n$. This is the symplectic signature in one figure.

---

## Exercise 7.2 (★) — Time-step convergence

Run velocity-Verlet for the same harmonic oscillator at $\Delta t \in \{0.5, 0.2, 0.1, 0.05, 0.02, 0.01\}$. For each, compute the maximum amplitude of energy fluctuations over $10^4$ steps. Plot on a log-log scale and confirm the $\Delta t^2$ scaling.

**Solution.** The maximum amplitude scales as $\Delta t^2$ because Verlet conserves a shadow Hamiltonian $\tilde H = H + O(\Delta t^2)$. The log-log plot should be a straight line of slope 2.

```python
dts = np.array([0.5, 0.2, 0.1, 0.05, 0.02, 0.01])
amps = []
for dt in dts:
    n = int(50 / dt)
    x, v = vv(1.0, 0.0, dt, n)
    E = 0.5 * (v**2 + x**2)
    amps.append(E.max() - E.min())
amps = np.array(amps)
slope, intercept = np.polyfit(np.log(dts), np.log(amps), 1)
print(f"Slope = {slope:.2f}")    # expect ≈ 2.0
```

For very small $\Delta t$ floating-point accumulation noise becomes visible (slope dips below 2). For very large $\Delta t$ the resonance instability $\omega \Delta t > 2$ kicks in.

---

## Exercise 7.3 (★★) — Minimum image convention bug-hunt

The following function is supposed to compute the squared distance between two atoms in a cubic box with the minimum image convention. Find the bug and fix it.

```python
def buggy_distance2(r1, r2, L):
    dr = r1 - r2
    dr = dr % L
    return np.sum(dr ** 2)
```

**Solution.** The modulo operation maps $\Delta x$ into $[0, L)$, not $[-L/2, L/2]$. Two atoms at positions 0.1 and $L - 0.1$ would have $\Delta x \approx L - 0.2$ rather than $-0.2$. The fix:

```python
def distance2(r1, r2, L):
    dr = r1 - r2
    dr -= L * np.round(dr / L)
    return np.sum(dr ** 2)
```

Test: `r1 = [0.1, 0.0, 0.0]`, `r2 = [L - 0.1, 0.0, 0.0]`, $L = 10$. The fixed version returns 0.04, the bug returns 96.04.

---

## Exercise 7.4 (★★) — Berendsen suppresses fluctuations

Run an MD simulation of an ideal Einstein solid (200 independent 3D harmonic oscillators with $\omega = 1$) for $10^5$ steps with $\Delta t = 0.02$, once with a Nosé-Hoover thermostat and once with a Berendsen thermostat, both targeted at $T = 1$ (reduced units; equipartition gives $\langle KE\rangle = 3N/2$). Histogram the instantaneous kinetic temperature. Compare to the analytical canonical prediction $\sigma_T^2/T^2 = 2/(3N) = 1/300$.

**Solution sketch.** Nosé-Hoover will give a histogram whose width matches the canonical prediction (relative $\sigma$ near 0.058). Berendsen will give a histogram of width perhaps 0.01 — much narrower than canonical — because the Berendsen scheme directly damps temperature fluctuations. This is the textbook demonstration of why Berendsen is not used for production sampling. Specific heat from energy variance with Berendsen will be wrong by roughly the same factor squared.

The code skeleton (NH part):

```python
def nh_step(x, v, zeta, dt, Q, gkT):
    a = -x
    v_half = (v + 0.5 * dt * a) / (1 + 0.5 * dt * zeta)
    x_new = x + dt * v_half
    KE_half = 0.5 * np.sum(v_half ** 2)
    zeta_new = zeta + dt * (2 * KE_half - gkT) / Q
    a_new = -x_new
    v_new = (v_half * (1 - 0.5 * dt * zeta_new) + 0.5 * dt * a_new)
    return x_new, v_new, zeta_new
```

A proper Trotter splitting (Martyna-Tuckerman) is preferred for production; this Euler-style update is enough for illustration.

---

## Exercise 7.5 (★★) — Lennard-Jones equilibrium and pressure

Set up a LAMMPS simulation of 4000 LJ argon atoms at the conditions of §7.5 Workflow 1. Compute the equilibrium pressure as a function of density $\rho^* = N\sigma^3/V \in \{0.6, 0.7, 0.8, 0.85, 0.9\}$ at fixed $T^* = 1.5$. Plot $P^*$ vs $\rho^*$ and compare to the LJ equation of state at this temperature.

**Solution.** Run five simulations, one per density, each for $\sim 50$ ps after $\sim 20$ ps equilibration. Extract `press` and `vol` from the LAMMPS log. Convert to reduced units:

$$
P^* = P\sigma^3/\epsilon, \qquad \rho^* = N\sigma^3/V, \qquad T^* = k_B T/\epsilon.
$$

The expected isotherm at $T^* = 1.5$ (well above the critical $T_c^* \approx 1.31$) is monotonic, with $P^* \approx 0.5$ at $\rho^* = 0.6$ rising to $P^* \approx 5$ at $\rho^* = 0.9$. Tabulated values from Johnson, Zollweg, Gubbins (1993) give 1% reference data. Discrepancies larger than 1–2% point to insufficient equilibration, missing tail corrections, or a wrong cutoff.

---

## Exercise 7.6 (★★) — MSD and diffusion in liquid argon

From the trajectory of Exercise 7.5 at $\rho^* = 0.85, T^* = 1.5$, compute the mean squared displacement using the NumPy code in §7.6. Plot $\mathrm{MSD}(t)$ vs $t$; fit a linear region; extract $D$. Compare to the Yeh-Hummer finite-size corrected value at infinite box size.

**Solution.** For LJ argon at the triple-point-like conditions $T^* = 1.5, \rho^* = 0.85$, the literature gives $D^* = D\sqrt{m/(\epsilon \sigma^2)} \approx 0.04$. Convert to physical units: with argon $\sigma = 3.405$ Å, $\epsilon/k_B = 120$ K, $m = 39.95$ amu, the natural time unit is $\tau = \sigma\sqrt{m/\epsilon} \approx 2.16$ ps; hence $D \approx 0.04 \sigma^2 / \tau \approx 0.04 \times 11.6 / 2.16$ Å$^2$/ps $\approx 0.21$ Å$^2$/ps.

After Yeh-Hummer correction with the simulation's viscosity (you can read $\eta$ off Green-Kubo, [Chapter 8](../ch08-statmech/03-transport.md), or pull from the literature) the infinite-box result is about 10% higher than the box-size simulated.

Pitfalls: use unwrapped positions; exclude the first ~1 ps of MSD (ballistic regime) and the last few percent (poor statistics) from the linear fit.

---

## Exercise 7.7 (★★★) — Build a Nosé-Hoover chain thermostat

Implement a Nosé-Hoover **chain** of $M = 3$ thermostats for the Einstein-solid system of Exercise 7.4. Compare the kinetic-energy histogram to (a) a single NH thermostat, (b) the analytical canonical $\chi^2$ distribution with $3N$ degrees of freedom. Show that a single NH gives the wrong distribution for $N = 1$ (a single oscillator) while NH chain recovers it.

**Solution outline.** The chain equations of motion couple thermostat $k$ to thermostat $k+1$:

$$
\dot{\zeta}_1 = \frac{1}{Q_1}\left[2 KE - g k_B T\right] - \zeta_1 \zeta_2,
$$

$$
\dot{\zeta}_k = \frac{1}{Q_k}\left[Q_{k-1} \zeta_{k-1}^2 - k_B T\right] - \zeta_k \zeta_{k+1}, \quad k = 2, \ldots, M-1,
$$

$$
\dot{\zeta}_M = \frac{1}{Q_M}\left[Q_{M-1} \zeta_{M-1}^2 - k_B T\right],
$$

with the physical EOM $\dot{\mathbf{p}}_i = -\nabla_i U - \zeta_1 \mathbf{p}_i$. Use a Trotter splitting: half-step thermostat updates, full-step Verlet, half-step thermostat updates again. For $M = 3$ and the Einstein solid, the resulting kinetic-energy histogram matches the canonical $\chi^2_{3N}$ form to well within the statistical noise of $10^5$ samples. For $N = 1$, a single NH produces a multi-peaked, non-canonical histogram (the Frenkel-Smit textbook has the classic figure); the chain washes this out.

This is also the implementation under the hood of LAMMPS' `fix nvt`, which uses a 3-chain by default.

---

## Exercise 7.8 (★★★) — Cu melting point bracketing

Estimate the EAM melting temperature of copper by running NPT simulations at $T \in \{1100, 1200, 1300, 1400, 1500\}$ K, starting from (a) a perfect crystal and (b) a quenched liquid. Track the equilibrium volume in each case. The midpoint of the hysteresis window is an upper-and-lower bracket for $T_m$. Then refine using a two-phase simulation.

**Solution sketch.**

(a) Crystals are kinetically stable past $T_m$ because the liquid must nucleate; expect the volume to remain crystalline up to perhaps 1450 K before suddenly jumping.

(b) Quenched liquids supercool because the crystal must nucleate; expect them to remain liquid down to perhaps 1100 K before suddenly crystallising — or, with bad luck, all the way to room temperature without crystallising.

The bracket from (a)+(b) is $T_m \in [1200, 1450]$ K. To refine, use the **two-phase method**: take a 16x16x4 supercell, melt half (region $z > L_z/2$) by heating it locally to 2500 K with a region thermostat while keeping the other half cold, equilibrate the two-phase configuration, then run NPT with no temperature ramp. The system spontaneously moves the solid-liquid interface up or down depending on which side $T$ favours; at $T = T_m$ the interface is stationary (within fluctuations). This pinpoints $T_m$ to within a few K for a given potential.

For Mishin Cu EAM, the published two-phase $T_m$ is $\sim 1325$ K versus the experimental 1358 K. Discrepancies of 30–60 K are typical for classical force fields. This is the kind of error that MLIPs ([Chapter 9](../ch09-mlip/index.md)) can shave down to 5–10 K when trained on appropriate DFT reference data.

---

## Summary of what these exercises cover

| Exercise | Skill |
|---|---|
| 7.1 | Implementing and comparing integrators |
| 7.2 | Empirical verification of integrator order |
| 7.3 | Periodic boundary debugging |
| 7.4 | Sampling distributions and thermostat correctness |
| 7.5 | Equation of state from MD |
| 7.6 | Transport coefficient from MSD |
| 7.7 | Implementing a chain thermostat |
| 7.8 | Phase boundaries from MD |

If you can do all eight comfortably, you understand the mechanics of an MD simulation well enough to read research papers critically and to set up your own production calculations. The remaining limitation is the force field itself, which is the subject of Chapter 9.
