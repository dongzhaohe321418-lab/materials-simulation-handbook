# 3b.5 — Phonons

<figure markdown>
![1D monatomic and diatomic chain phonon dispersions](../assets/figures/ch03b/fig_1d_chain_phonons.png){ width="700" }
<figcaption>Figure 3b.5.1. Phonon dispersions of the simplest 1D models. The monatomic chain (left) has a single acoustic branch with a linear (sound-like) regime near \(k=0\). The diatomic chain (right) splits this into an acoustic branch and an optical branch separated by a band gap that opens at the zone boundary when the two masses differ.</figcaption>
</figure>

> *"A solid is a forest of harmonic oscillators."* — every solid state course since 1907

So far we have quantised the electrons and treated the nuclei as fixed in their crystallographic positions. This is the Born–Oppenheimer approximation, and it is excellent for most ground-state electronic structure questions. But the moment we ask anything about *thermal* properties — heat capacity, thermal expansion, sound velocity, lattice stability — we must allow the nuclei to move. The classical lattice vibrations decompose into normal modes; the quanta of those modes are *phonons*. Phonons are the bridge between electronic structure (Tier 1) and finite-temperature statistical mechanics (Ch 8), and they are the most direct benchmark for any interatomic potential, classical or machine-learned.

We build the theory in the standard pedagogical sequence: a 1D monatomic chain, then a 1D diatomic chain (which already exhibits the most important qualitative features — acoustic and optical branches), then generalise to 3D. Working code for the diatomic chain finishes the section.

## 3b.5.1 The 1D monatomic chain

Consider $N$ identical atoms of mass $m$ arranged on a line with equilibrium spacing $a$. Couple each atom to its two nearest neighbours by springs of stiffness $K$. Let $u_n(t)$ be the displacement of atom $n$ from its equilibrium position $na$. Newton's second law for atom $n$ reads

$$m\ddot u_n = K(u_{n+1} - u_n) + K(u_{n-1} - u_n) = K(u_{n+1} + u_{n-1} - 2u_n). \tag{3b.5.1}$$

Impose periodic boundary conditions $u_{n+N} = u_n$ (Born–von Kármán again). The lattice has translational symmetry, so look for solutions of Bloch form — plane waves in *site index*:

$$u_n(t) = A\, e^{i(kna - \omega t)}, \tag{3b.5.2}$$

with $k$ taking $N$ discrete values in the BZ $-\pi/a < k \le \pi/a$. Substituting (3b.5.2) into (3b.5.1) and dividing through by the common exponential,

$$-m\omega^2 A = K\,A\, [e^{ika} + e^{-ika} - 2] = 2KA[\cos(ka) - 1] = -4KA\sin^2(ka/2). \tag{3b.5.3}$$

Solving for $\omega$,

$$\boxed{\; \omega(k) = 2\sqrt{K/m}\, \left|\sin(ka/2)\right|. \;} \tag{3b.5.4}$$

This is the dispersion relation of the 1D monatomic chain. Several features merit attention.

**Long-wavelength limit ($k\to 0$).** $\sin(ka/2)\approx ka/2$ so $\omega \approx \sqrt{K/m}\cdot ka = c_s k$, where

$$c_s = a\sqrt{K/m} \tag{3b.5.5}$$

is the *sound velocity*. The chain supports sound waves with linear dispersion at long wavelengths, as it must — the elastic continuum limit of any 3D crystal also gives linear dispersion (acoustic phonons).

**Zone boundary ($k = \pi/a$).** $\omega_\text{max} = 2\sqrt{K/m}$, a hard maximum frequency. Above this, no propagating modes exist. This is the lattice analogue of a cutoff frequency in a dispersive transmission line.

**Periodicity in $k$.** $\omega(k + 2\pi/a) = \omega(k)$, so all distinct modes lie in the first BZ.

**Number of modes.** Exactly $N$ values of $k$ in the BZ, one for each atom: the number of normal modes equals the number of degrees of freedom, as it must.

**Counting:** in 3D the chain becomes a 3D crystal of $N$ atoms with $3N$ degrees of freedom; the BZ houses $N$ wavevectors, each carrying *three* phonon polarisations (one longitudinal, two transverse for an isotropic medium). The count comes out right.

## 3b.5.2 The 1D diatomic chain — full derivation

Now place two distinct atoms per unit cell: atom A of mass $m_1$ at position $na$ (the cell origin), atom B of mass $m_2$ at position $na + d$ (somewhere inside the cell). Connect each A to each B inside the same cell with spring $K$, and each A to the B of the previous cell also with $K$ (we take the alternating-spring model: same stiffness, alternating bond lengths). Let $u_n$ be the displacement of A in cell $n$, and $v_n$ the displacement of B in cell $n$. The equations of motion are

$$m_1 \ddot u_n = K(v_n + v_{n-1} - 2u_n), \tag{3b.5.6}$$
$$m_2 \ddot v_n = K(u_{n+1} + u_n - 2v_n). \tag{3b.5.7}$$

Use the Bloch ansatz

$$u_n(t) = U\, e^{i(kna - \omega t)}, \qquad v_n(t) = V\, e^{i(kna - \omega t)}. \tag{3b.5.8}$$

Substituting into (3b.5.6),

$$-m_1\omega^2 U = K(V + V e^{-ika} - 2U) = K[V(1 + e^{-ika}) - 2U]. \tag{3b.5.9}$$

Substituting into (3b.5.7),

$$-m_2\omega^2 V = K(U e^{ika} + U - 2V) = K[U(1 + e^{ika}) - 2V]. \tag{3b.5.10}$$

Collect terms into a matrix equation for $(U, V)^T$:

$$\begin{pmatrix} 2K - m_1\omega^2 & -K(1 + e^{-ika}) \\ -K(1 + e^{ika}) & 2K - m_2\omega^2 \end{pmatrix}\begin{pmatrix} U \\ V\end{pmatrix} = 0. \tag{3b.5.11}$$

Nontrivial solutions require the determinant to vanish:

$$(2K - m_1\omega^2)(2K - m_2\omega^2) - K^2|1 + e^{ika}|^2 = 0. \tag{3b.5.12}$$

Use $|1 + e^{ika}|^2 = 2 + 2\cos(ka) = 4\cos^2(ka/2)$. Expand:

$$4K^2 - 2K(m_1 + m_2)\omega^2 + m_1 m_2 \omega^4 = 4K^2 \cos^2(ka/2), \tag{3b.5.13}$$

$$m_1 m_2 \omega^4 - 2K(m_1 + m_2)\omega^2 + 4K^2\sin^2(ka/2) = 0. \tag{3b.5.14}$$

This is a quadratic in $\omega^2$. Solve:

$$\boxed{\; \omega_\pm^2(k) = \frac{K(m_1 + m_2)}{m_1 m_2}\left[1 \pm \sqrt{1 - \frac{4 m_1 m_2}{(m_1 + m_2)^2}\sin^2(ka/2)}\right]. \;} \tag{3b.5.15}$$

The two branches $\omega_-$ and $\omega_+$ are the **acoustic** and **optical** phonons respectively. We inspect each.

**Acoustic branch** ($\omega_-$).
- At $k=0$: $\sin = 0$ so $\omega_- = 0$. The acoustic branch is gapless at the BZ centre.
- Near $k=0$: Taylor expand the square root. The result, after a few lines of algebra, is $\omega_- \approx ka\sqrt{K/[2(m_1 + m_2)]}$ — a linear, gapless dispersion with sound speed $c_s = a\sqrt{K/[2(m_1 + m_2)]}$. The two atoms move *in phase*: $U \approx V$.
- At $k = \pi/a$ (zone boundary): $\sin^2(ka/2) = 1$ so $\omega_-^2 = 2K/\max(m_1, m_2)$ — depends on the *heavier* mass.

**Optical branch** ($\omega_+$).
- At $k=0$: $\omega_+^2 = 2K(m_1 + m_2)/(m_1 m_2) = 2K/\mu$, where $\mu = m_1 m_2/(m_1+m_2)$ is the reduced mass. The optical branch has a *finite* frequency at $\mathbf k = 0$. The atoms move *out of phase*: $m_1 U = -m_2 V$, i.e. the centre of mass is at rest while the two atoms oscillate against each other. This is the physical reason optical phonons in ionic crystals couple to light: they create an oscillating electric dipole.
- At $k = \pi/a$: $\omega_+^2 = 2K/\min(m_1, m_2)$ — depends on the *lighter* mass.

**Gap at the zone boundary.** Between the top of the acoustic branch and the bottom of the optical branch there is a gap whose size is set by the mass difference:

$$\omega_+^2(\pi/a) - \omega_-^2(\pi/a) = 2K\left[\frac{1}{m_1} - \frac{1}{m_2}\right] \quad (m_2 > m_1). \tag{3b.5.16}$$

If $m_1 = m_2$ the gap closes and the two branches join — the diatomic chain *with equal masses* is just a monatomic chain with half the lattice constant, and the optical branch is the *backfold* of the upper half of the original acoustic branch.

## 3b.5.3 3D phonons: the dynamical matrix

In a 3D crystal with $N_\text{at}$ atoms per unit cell, the displacement of atom $a$ in cell $\mathbf R$ along Cartesian direction $\alpha$ is $u_\alpha^a(\mathbf R, t)$. The harmonic potential energy is

$$V = \frac{1}{2}\sum_{\mathbf R \mathbf R'}\sum_{ab}\sum_{\alpha\beta} \Phi_{\alpha\beta}^{ab}(\mathbf R - \mathbf R')\, u_\alpha^a(\mathbf R)\, u_\beta^b(\mathbf R'), \tag{3b.5.17}$$

with $\Phi_{\alpha\beta}^{ab}(\mathbf R)$ the **force constant matrix** — the second derivative of the Born–Oppenheimer energy with respect to atomic displacements, evaluated at equilibrium:

$$\Phi_{\alpha\beta}^{ab}(\mathbf R - \mathbf R') = \frac{\partial^2 E_\text{BO}}{\partial u_\alpha^a(\mathbf R)\, \partial u_\beta^b(\mathbf R')}\bigg|_\text{eq}. \tag{3b.5.18}$$

By translational symmetry, $\Phi$ depends only on $\mathbf R - \mathbf R'$. The equations of motion are

$$m_a \ddot u_\alpha^a(\mathbf R) = -\sum_{\mathbf R' b\beta} \Phi_{\alpha\beta}^{ab}(\mathbf R - \mathbf R') u_\beta^b(\mathbf R'). \tag{3b.5.19}$$

Use the Bloch ansatz $u_\alpha^a(\mathbf R, t) = \epsilon_\alpha^a\, e^{i(\mathbf k\cdot\mathbf R - \omega t)}/\sqrt{m_a}$ (the mass-weighted normalisation is conventional). Substitution gives the eigenvalue problem

$$\omega^2\, \epsilon_\alpha^a = \sum_{b\beta} D_{\alpha\beta}^{ab}(\mathbf k)\, \epsilon_\beta^b, \tag{3b.5.20}$$

with the **dynamical matrix**

$$\boxed{\; D_{\alpha\beta}^{ab}(\mathbf k) = \frac{1}{\sqrt{m_a m_b}} \sum_\mathbf R \Phi_{\alpha\beta}^{ab}(\mathbf R)\, e^{i\mathbf k\cdot\mathbf R}. \;} \tag{3b.5.21}$$

$D(\mathbf k)$ is a $3N_\text{at}\times 3N_\text{at}$ Hermitian matrix. Its eigenvalues are $\omega_\nu^2(\mathbf k)$ — there are $3N_\text{at}$ phonon branches at each $\mathbf k$. Of these, three are acoustic (gapless at $\mathbf k = 0$) and $3N_\text{at} - 3$ are optical. The eigenvectors $\boldsymbol\epsilon_\nu(\mathbf k)$ are the **polarisation vectors**, telling you the relative motion of the atoms in each mode.

The dynamical matrix is a Fourier transform of the force-constant matrix. In practice $\Phi(\mathbf R)$ has short range — it decays rapidly with $|\mathbf R|$ — so only a few neighbour shells of $\Phi$ need be stored. This is why phonon calculations are tractable: a $3N_\text{at}\times 3N_\text{at}$ matrix at each $\mathbf k$.

## 3b.5.4 Force constants from DFT

How do you actually compute $\Phi$? Two ingredients.

**Finite differences.** Displace one atom in a *supercell* by a small amount $\Delta u$, run a self-consistent DFT calculation, read off the forces on every other atom. Each pair of (displacement, force) entries gives one column of the force-constant matrix. Repeat for all symmetry-inequivalent displacements. The python package `phonopy` automates this entire workflow, taking a primitive-cell DFT calculation and producing the full $\Phi(\mathbf R)$, the dynamical matrix at any $\mathbf k$, and the resulting phonon dispersion. This is the workhorse method.

**Density-functional perturbation theory (DFPT).** Compute $\Phi(\mathbf q)$ at any wavevector $\mathbf q$ directly by linear response, without ever building a supercell. This is faster than finite differences and is built into Quantum ESPRESSO (`ph.x`) and ABINIT.

You will use phonopy in Chapter 6 to compute the phonon spectrum of a real material and check that all frequencies are real — imaginary $\omega^2$ indicates a structural instability, and is the first-line tool for detecting that you have not in fact found the ground state structure.

## 3b.5.5 Python: phonon dispersion of a 1D diatomic chain

The following code builds the dynamical matrix (3b.5.11) for the 1D diatomic chain at 200 wavevectors and plots both branches.

```python
"""Phonon dispersion of a 1D diatomic chain.

Builds the dynamical matrix from analytical expressions, diagonalises
at 200 k-points, and plots acoustic + optical branches.
"""
from __future__ import annotations
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

# Parameters (arbitrary units)
M1: float = 1.0      # mass of atom A
M2: float = 3.0      # mass of atom B (heavier)
K_SPRING: float = 1.0  # spring constant
A_LATT: float = 1.0    # lattice constant of the diatomic cell

def dynamical_matrix(k: float) -> npt.NDArray[np.complex128]:
    """Return the 2x2 dynamical matrix at wavevector k.

    Constructed by reading off Eq. (3b.5.11) and dividing by sqrt(m_a m_b).
    """
    D: npt.NDArray[np.complex128] = np.zeros((2, 2), dtype=np.complex128)
    D[0, 0] = (2 * K_SPRING) / M1
    D[1, 1] = (2 * K_SPRING) / M2
    coupling: complex = -K_SPRING * (1.0 + np.exp(-1j * k * A_LATT))
    D[0, 1] = coupling / np.sqrt(M1 * M2)
    D[1, 0] = np.conj(D[0, 1])
    return D

def phonon_frequencies(k: float) -> npt.NDArray[np.float64]:
    """Diagonalise the dynamical matrix; return sqrt of (real) eigenvalues."""
    eigs_omega2: npt.NDArray[np.float64] = np.linalg.eigvalsh(dynamical_matrix(k))
    # Clip tiny negatives from round-off
    eigs_clip: npt.NDArray[np.float64] = np.clip(eigs_omega2, a_min=0.0, a_max=None)
    return np.sqrt(eigs_clip)

def analytical_branches(k: npt.NDArray[np.float64]) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Closed-form Eq. (3b.5.15) for cross-checking."""
    pref: float = K_SPRING * (M1 + M2) / (M1 * M2)
    discr: npt.NDArray[np.float64] = 1.0 - (4 * M1 * M2 / (M1 + M2) ** 2) * np.sin(k * A_LATT / 2) ** 2
    omega_minus: npt.NDArray[np.float64] = np.sqrt(pref * (1.0 - np.sqrt(discr)))
    omega_plus: npt.NDArray[np.float64] = np.sqrt(pref * (1.0 + np.sqrt(discr)))
    return omega_minus, omega_plus

def main() -> None:
    n_k: int = 200
    k_vals: npt.NDArray[np.float64] = np.linspace(-np.pi / A_LATT, np.pi / A_LATT, n_k)
    omega_numeric: npt.NDArray[np.float64] = np.array(
        [phonon_frequencies(float(k)) for k in k_vals]
    )  # shape (n_k, 2)
    omega_minus, omega_plus = analytical_branches(k_vals)

    # Sanity check
    err_acoustic: float = float(np.max(np.abs(omega_numeric[:, 0] - omega_minus)))
    err_optical: float = float(np.max(np.abs(omega_numeric[:, 1] - omega_plus)))
    print(f"Max |numerical - analytical| acoustic = {err_acoustic:.2e}")
    print(f"Max |numerical - analytical| optical  = {err_optical:.2e}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(k_vals * A_LATT / np.pi, omega_numeric[:, 0],
            color="C0", lw=2, label="acoustic (numerical)")
    ax.plot(k_vals * A_LATT / np.pi, omega_numeric[:, 1],
            color="C3", lw=2, label="optical (numerical)")
    ax.plot(k_vals * A_LATT / np.pi, omega_minus,
            color="C0", lw=0, marker=".", markersize=2, label="acoustic (Eq. 3b.5.15)")
    ax.plot(k_vals * A_LATT / np.pi, omega_plus,
            color="C3", lw=0, marker=".", markersize=2, label="optical (Eq. 3b.5.15)")
    ax.axvline(-1.0, color="grey", lw=0.5, linestyle=":")
    ax.axvline(1.0, color="grey", lw=0.5, linestyle=":")
    ax.set_xlabel(r"$k\, a/\pi$")
    ax.set_ylabel(r"$\omega$ (arb. units)")
    ax.set_title(
        f"Diatomic chain phonons (m1={M1}, m2={M2}, K={K_SPRING})"
    )
    ax.legend()
    plt.tight_layout()
    plt.savefig("diatomic_phonons.pdf")
    plt.show()

if __name__ == "__main__":
    main()
```

You should see an acoustic branch vanishing linearly at $k=0$, an optical branch starting at $\omega_+(0) = \sqrt{2K/\mu} \approx \sqrt{2 \cdot 1 / 0.75} \approx 1.633$, and a clear gap at the BZ boundary $k = \pi/a$ of size $\omega_+^2 - \omega_-^2 = 2K(1/m_1 - 1/m_2) = 2(1 - 1/3) = 4/3$. The numerical and analytical results agree to round-off.

## 3b.5.6 Connection to molecular dynamics

In Chapter 7 you will run classical molecular dynamics (MD): integrate Newton's equations for $N$ atoms in a thermostat, generating a trajectory. The atoms vibrate around their equilibrium positions; the time series of atomic positions contains all the lattice-dynamical information of the crystal. Specifically:

1. **Mean-square displacement** of an atom is set by the Bose–Einstein average of $|\boldsymbol\epsilon|^2 / \omega$, summed over modes. In the classical limit ($k_B T \gg \hbar\omega$) this becomes $\langle u^2\rangle = k_B T/(m\omega^2)$ — the equipartition theorem applied to harmonic oscillators.

2. **Velocity autocorrelation function** $C(t) = \langle \mathbf v(0)\cdot\mathbf v(t)\rangle$ is, by the Wiener–Khinchin theorem, the Fourier transform of the *phonon density of states*. So you can extract the phonon DOS from any sufficiently long MD trajectory — a method called the velocity autocorrelation method, or sometimes (in jargon) "thermostat phonons".

3. **Validity of classical MD.** The catch in (1) — that the classical formula assumes $k_B T \gg \hbar\omega$ — is the central restriction on MD. In hydrogen-containing solids, the high-frequency O–H stretches have $\hbar\omega/k_B \sim 4000$ K, far above room temperature. Classical MD therefore *over-populates* these modes and gets the heat capacity too high by a factor of three at low temperatures: equipartition gives $3Nk_B$ instead of the experimental $T^3$ falloff. This is the Einstein/Debye problem of the next section, viewed from the MD side.

In other words: classical MD *is* lattice dynamics, but in the classical limit. Phonons quantise the same physics. They are two views of the same picture.

!!! note "Phonons and MLIPs"
    A machine-learning interatomic potential must reproduce, at minimum, the phonon dispersion of the training material. The reason is that vibrational free energies (Ch 8), thermal expansion, sound velocities, and Raman/infrared spectra all depend on $\omega_\nu(\mathbf k)$. The standard benchmark of an MLIP — used in every Ch 9 paper — is the per-mode phonon error: RMSE in THz between MLIP and DFT, computed at a dense $\mathbf k$-grid in the BZ. State of the art models (MACE, NequIP) routinely achieve $\le 0.2$ THz across all branches.

## Where this is used later

- **Tier 1.** §6.6 (running phonon calculations with phonopy and DFPT), §7.2 (molecular dynamics as classical lattice dynamics in disguise), §7.5 (interatomic potentials must reproduce the lowest few moments of the phonon DOS).
- **Tier 2.** §8.2 (vibrational free energy from phonon DOS, Helmholtz free energy at finite temperature), §8.4 (thermal expansion via Grüneisen parameters), §9.7 (phonon benchmarks for MLIPs), §10.6 (equivariant features and the connection to phonon polarisation vectors).
- **Capstone Project 2.** Computing thermal conductivity in a thermoelectric: the Boltzmann transport equation for phonons takes $\omega_\nu(\mathbf k)$ as input.

Next, §3b.6, where we feed the phonon spectrum into a statistical mechanics machine and recover the specific heat — Einstein for room-temperature intuition, Debye for the low-temperature truth.
