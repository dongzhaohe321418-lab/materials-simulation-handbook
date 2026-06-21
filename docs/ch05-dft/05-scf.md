# 5.5 The Self-Consistent Field Loop

**What problem are we solving?** The Kohn–Sham potential $v_\mathrm{KS}$ tells us how the electrons move, so from it we can work out where the electrons are — the density $n$. But the potential itself is *built from* the density (through the Hartree and exchange–correlation terms). So we cannot compute the potential until we know the density, and we cannot compute the density until we know the potential: a chicken-and-egg loop. **Self-consistency** is how we break it. We make a guess, use it to compute a better guess, and repeat until the input and the output agree — at which point the density is *consistent with the potential it generates*. Everything in this section is machinery for doing that loop reliably.

!!! note "In plain language"
    The SCF loop is a *fixed-point iteration*: we are looking for a density that maps to itself under the rule "build the potential, solve for the orbitals, read off the new density". Concretely: **guess** a density; **build** the potential from it; **solve** the Kohn–Sham equations to get orbitals; **read off** a NEW density from those orbitals; **mix** the new density with the old one; and **repeat** until the density stops changing. A density that no longer changes is a *fixed point* — feeding it in gives the same thing back out. This is exactly the "$x = f(x)$, so iterate $x \leftarrow f(x)$" pattern described in the [formula-reading guide](../undergraduate/formula-reading-guide.md); here $x$ is the whole density $n(\mathbf r)$ and $f$ is one trip round the loop.

!!! note "Why does this chapter exist?"
    The Kohn–Sham equations (§5.3) look like a one-electron Schrödinger equation: you put in a potential $v_\mathrm{KS}$, you get out orbitals and a density. But $v_\mathrm{KS}$ depends on the density! It's a chicken-and-egg problem: to compute the potential we need the density, to compute the density we need to know the orbitals, to compute the orbitals we need the potential. This circular dependence is *the* defining feature of mean-field theories — Hartree, Hartree–Fock, Kohn–Sham — and the standard cure is iteration: guess, solve, update, repeat.
    
    A useful analogy. Imagine you are pouring water into a bowl whose shape changes as a function of how much water is already in it. If you pour the first cup, the bowl widens, so the water spreads out, which changes how the bowl wants to widen, and so on. The water level eventually settles at a *self-consistent* point where the bowl shape and the water level agree. That is what an SCF iteration does. But you can imagine pathological cases where naive pouring causes the water to slosh: the bowl widens too much, so the water becomes shallow, so the bowl narrows again, so the water becomes deep, and so on forever. Real DFT systems, especially metals, exhibit exactly this "charge sloshing". The cure is *not* to pour the new water in directly: pour only a fraction (linear mixing), or look at a history of past pourings and extrapolate cleverly (Pulay/DIIS, Anderson, Broyden).
    
    This section is half theory, half pseudocode. By the end, you should be able to write your own SCF loop from scratch — and we include a complete 150-line Python implementation that does exactly that for a 1D hydrogen chain.

<figure markdown>
```mermaid
stateDiagram-v2
    [*] --> Guess : initial guess n⁰(r)
    Guess --> Veff : build v_KS[n] = v_ext + v_H[n] + v_xc[n]
    Veff --> Solve : solve (−½∇² + v_KS) φᵢ = εᵢ φᵢ
    Solve --> NewDensity : n_out(r) = Σ |φᵢ|² (occupied)
    NewDensity --> Mix : mix n_in and n_out (Pulay/Broyden)
    Mix --> Check : converged?
    Check --> Veff : no — iterate
    Check --> [*] : yes — output E, forces, ρ
```
<figcaption>State diagram of the Kohn–Sham self-consistent field loop: starting from an initial density guess, the effective Kohn–Sham potential is built, the one-electron equations are diagonalised for new orbitals, a new density is formed from the occupied orbitals and mixed with the input, and convergence is checked — if not converged the loop returns to rebuild the potential, and if converged it outputs the energy, forces, and density.</figcaption>
</figure>

<figure markdown>
![Two-panel plot of SCF convergence: the total energy approaching its converged value roughly exponentially over iterations on the left, and the energy change per step on a log scale dropping below the convergence threshold after a few tens of iterations on the right](../assets/figures/ch05/fig_scf_convergence.png){ width="750" }
<figcaption>Figure 5.5.1. Typical SCF convergence behaviour (synthetic example). The total energy approaches the converged value approximately exponentially (left), and the energy change per step \(|\Delta E|\) drops below the user-specified threshold (here \(10^{-6}\) Ry) after a few tens of iterations (right). Real calculations may oscillate before locking in, especially for metals and magnetic systems.</figcaption>
</figure>

!!! abstract "Key idea (Chapter 5.5)"
    The Kohn–Sham equations are nonlinear: $v_\mathrm{KS}$ depends on $n$, and $n$ depends on the orbitals, which depend on $v_\mathrm{KS}$. Solving them by *self-consistent iteration* — guess, solve, update, repeat — is the practical heart of every DFT calculation. Naive iteration oscillates (charge sloshing); modern codes use *Pulay/DIIS*, *Anderson*, or *Broyden* mixing to accelerate convergence by orders of magnitude. This section walks through the algorithm, derives why it oscillates, develops the mixing schemes, and provides a complete 150-line Python implementation.

The Kohn–Sham equations,

$$
\Big[-\tfrac{1}{2}\nabla^{2} + v_\mathrm{KS}[n](\mathbf r)\Big]\phi_i(\mathbf r) = \varepsilon_i\,\phi_i(\mathbf r),
\qquad
n = \sum_i^\mathrm{occ}|\phi_i|^{2},
\qquad
v_\mathrm{KS}[n] = v_\mathrm{ext} + v_H[n] + v_{xc}[n],
$$

are nonlinear: the operator on the left depends, through $v_\mathrm{KS}$, on the density that the eigenfunctions themselves produce. The standard way to solve a nonlinear equation $n = \mathcal F[n]$ is fixed-point iteration: take an initial guess $n^{(0)}$, evaluate $n^{(1)} = \mathcal F[n^{(0)}]$, and repeat until $\|n^{(k+1)} - n^{(k)}\|$ is small.

In DFT this is the **self-consistent field** (SCF) loop. Naive iteration almost always fails to converge for systems beyond hydrogen-like atoms; the loop oscillates between charge-rich and charge-poor solutions, sometimes diverging outright. Practical SCF codes use *mixing* schemes — careful linear combinations of densities (or potentials) from successive iterations — to suppress these oscillations.

This section walks through the algorithm, explains why naive iteration fails, develops the mixing schemes (linear, Pulay/DIIS, Anderson), and ends with a complete Python implementation that solves a 1D model "hydrogen chain" using LDA exchange and finite differences.

## 5.5.1 The basic SCF algorithm

The textbook KS-SCF loop is:

1. **Initial guess.** Construct an initial density $n^{(0)}(\mathbf r)$. Common choices: superposition of free-atom densities, the previous SCF solution at a nearby geometry, or — for very simple systems — the uniform density.
2. **Build the potential.** Compute $v_\mathrm{KS}[n^{(k)}] = v_\mathrm{ext} + v_H[n^{(k)}] + v_{xc}[n^{(k)}]$.
3. **Diagonalise.** Solve the eigenvalue problem $\hat{H}_\mathrm{KS}\phi_i = \varepsilon_i\phi_i$ for the lowest $N_\mathrm{occ}$ eigenpairs.
4. **Form the new density.** $n_\mathrm{out}^{(k)}(\mathbf r) = \sum_i^\mathrm{occ}|\phi_i(\mathbf r)|^{2}$.
5. **Check convergence.** Compute residuals — change in density, change in energy, maximum force. If below tolerance, stop.
6. **Mix.** $n^{(k+1)} = \mathcal M(n^{(k)}, n_\mathrm{out}^{(k)}; \text{history})$. Go to step 2.

The interesting step is 6.

!!! warning "Common misunderstanding"
    **Convergence is not the same as correctness.** It is tempting to read "SCF converged" as "the answer is right" — but they are different claims. A calculation can converge tightly to a self-consistent density (input density $=$ output density to many decimal places) that is still *wrong*, because the physics fed into the loop — the choice of functional ($v_{xc}$), the basis or grid, the pseudopotential, the **k**-point sampling — may be inadequate. The loop only guarantees that the density is consistent with the potential it generates *for the model you specified*; it cannot tell you the model is a good description of reality.

    A second, related trap: **mixing is a numerical convergence aid, not physics.** The mixing parameter $\alpha$, the Pulay history depth, and the Kerker preconditioner change how fast (and whether) you reach the fixed point — they do **not** change *which* fixed point you reach. Two runs with different $\alpha$ that both converge land on the *same* self-consistent density and the *same* energy. So you may freely tune mixing to make a stubborn calculation converge without worrying that you are biasing the physical result. (This is revisited in the closing warning of §5.5.7.)

## 5.5.2 Why naive iteration fails

The naive mixing scheme is $n^{(k+1)} = n_\mathrm{out}^{(k)}$ — accept the output density unchanged. This is fixed-point iteration on the map $\mathcal F: n \mapsto n_\mathrm{out}[n]$.

### Linearised analysis near a fixed point

To understand quantitatively why naive iteration can diverge, linearise the SCF map around its fixed point $n^{*}$. Writing $\delta n^{(k)} = n^{(k)} - n^{*}$,

$$
\delta n^{(k+1)} \approx \mathcal J\,\delta n^{(k)},
\qquad \mathcal J \equiv \frac{\delta n_\mathrm{out}}{\delta n_\mathrm{in}}\Big|_{n^{*}},
\tag{5.41a}
$$

where $\mathcal J$ is the (integral-kernel) **SCF Jacobian** evaluated at the fixed point. Its eigenvalues $\{\lambda_i\}$ determine the convergence: if $\max_i|\lambda_i|<1$, iteration converges with rate $\max_i|\lambda_i|$; if $\max_i|\lambda_i|>1$, the iteration diverges.

??? note "Full derivation: linearising the SCF map and factorising the Jacobian"
    Write the SCF map (one full trip round the loop) as $n_\mathrm{out} = \mathcal F[n_\mathrm{in}]$, and let $n^{*}$ be the fixed point, i.e. the self-consistent density that satisfies
    $$
    \mathcal F[n^{*}] = n^{*}.
    \tag{5.41b}
    $$
    Take the input density at iteration $k$ to be a small departure from the fixed point,
    $$
    n^{(k)} = n^{*} + \delta n^{(k)},
    \qquad
    \delta n^{(k)} \equiv n^{(k)} - n^{*},
    $$
    and Taylor-expand $\mathcal F$ as a functional about $n^{*}$. To first order in $\delta n^{(k)}$,
    $$
    \mathcal F[n^{*} + \delta n^{(k)}]
    = \mathcal F[n^{*}]
    + \int \frac{\delta n_\mathrm{out}(\mathbf r)}{\delta n_\mathrm{in}(\mathbf r')}\bigg|_{n^{*}}\,\delta n^{(k)}(\mathbf r')\,\mathrm d\mathbf r'
    + \mathcal O\!\big(\delta n^{2}\big).
    \tag{5.41c}
    $$
    The left-hand side is, by definition of the map, the *new* input density: with naive iteration $n^{(k+1)} = n_\mathrm{out}^{(k)} = \mathcal F[n^{(k)}]$, so the left-hand side equals $n^{(k+1)} = n^{*} + \delta n^{(k+1)}$. The first term on the right is $\mathcal F[n^{*}] = n^{*}$ by the fixed-point condition (5.41b). Subtracting $n^{*}$ from both sides cancels it and leaves, to linear order,
    $$
    \delta n^{(k+1)}(\mathbf r)
    \approx \int \mathcal J(\mathbf r,\mathbf r')\,\delta n^{(k)}(\mathbf r')\,\mathrm d\mathbf r',
    \qquad
    \mathcal J(\mathbf r,\mathbf r') \equiv \frac{\delta n_\mathrm{out}(\mathbf r)}{\delta n_\mathrm{in}(\mathbf r')}\bigg|_{n^{*}},
    $$
    which is (5.41a) written out as an integral kernel. The constant ($n^{*}$) drops out *because* we expand about the fixed point — this is why the error $\delta n$, not $n$ itself, is the natural object.

    **Chain-rule factorisation.** The map $n_\mathrm{in}\mapsto n_\mathrm{out}$ is not direct: it passes through the Kohn–Sham potential. A change in the input density changes $v_\mathrm{KS}$, and the changed potential changes the orbitals and hence the output density. By the functional chain rule,
    $$
    \mathcal J
    = \frac{\delta n_\mathrm{out}}{\delta n_\mathrm{in}}
    = \frac{\delta n_\mathrm{out}}{\delta v_\mathrm{KS}}\,
      \frac{\delta v_\mathrm{KS}}{\delta n_\mathrm{in}}
    \equiv \chi\,K,
    \tag{5.41d}
    $$
    where we have *defined* two pieces:

    - the **response function** $\chi \equiv \dfrac{\delta n_\mathrm{out}}{\delta v_\mathrm{KS}}$ — how the output density responds when the potential it is built from is perturbed (this is the independent-particle, or Kohn–Sham, susceptibility);
    - the **kernel** $K \equiv \dfrac{\delta v_\mathrm{KS}}{\delta n_\mathrm{in}}$ — how the potential changes when the density that builds it is perturbed. Since $v_\mathrm{KS} = v_\mathrm{ext} + v_H[n] + v_{xc}[n]$ and $v_\mathrm{ext}$ is density-independent, only the Hartree and exchange–correlation pieces contribute:
      $$
      K(\mathbf r,\mathbf r')
      = \frac{\delta v_H(\mathbf r)}{\delta n(\mathbf r')}
      + \frac{\delta v_{xc}(\mathbf r)}{\delta n(\mathbf r')}
      = \frac{1}{|\mathbf r-\mathbf r'|}
      + \frac{\delta^{2}E_{xc}}{\delta n(\mathbf r)\,\delta n(\mathbf r')},
      $$
      i.e. the bare Coulomb (Hartree) kernel plus the xc kernel.

    We write the kernel as $K$, not $v_\mathrm{eff}$, deliberately: in the loop diagram above $v_\mathrm{eff}$ already denotes the *KS potential itself*, whereas $K = \delta v_\mathrm{KS}/\delta n$ is its *functional derivative with respect to the density*. They are different objects and conflating the symbols would be a genuine error. With this factorisation, $\mathcal J = \chi K$: the divergence that drives charge sloshing lives in $\chi$ (the long-wavelength electronic response of a metal), while $K$ supplies the always-present Hartree $1/|\mathbf r-\mathbf r'|$ term that is itself large at small $|\mathbf q|$ (its Fourier transform is $4\pi/|\mathbf q|^{2}$).

For an insulating finite system, $\mathcal J$ has eigenvalues bounded by unity in magnitude (the dielectric response is small for short-wavelength density fluctuations and vanishes for the long-wavelength constant mode). For a metallic system, however, the long-wavelength response is large — the *dielectric function* $\epsilon(\mathbf q) \to \infty$ as $|\mathbf q|\to 0$ — and $\mathcal J$ can have eigenvalues of magnitude much larger than 1. This is the formal expression of charge sloshing.

!!! note "What is $\epsilon(\mathbf q)$, and why does it behave differently in metals and insulators?"
    The **dielectric function** $\epsilon(\mathbf q)$ measures how strongly the electron gas screens an externally imposed potential of wavevector $\mathbf q$: an applied perturbation $v_\mathrm{ext}(\mathbf q)$ produces a *total* (screened) potential $v_\mathrm{tot}(\mathbf q) = v_\mathrm{ext}(\mathbf q)/\epsilon(\mathbf q)$. A large $\epsilon$ means strong screening — the electrons rearrange to almost cancel the applied potential. Within the Thomas–Fermi model of a uniform electron gas,
    $$
    \epsilon_\mathrm{TF}(\mathbf q) = 1 + \frac{q_0^{2}}{|\mathbf q|^{2}},
    \tag{5.41e}
    $$
    where $q_0$ is the Thomas–Fermi screening wavevector (set by the density of states at the Fermi level). The contrast is then sharp:

    - **Metal.** There are gapless electron–hole excitations, so a long-wavelength perturbation is screened essentially perfectly: $\epsilon(\mathbf q)\to\infty$ as $|\mathbf q|\to 0$ (in (5.41e), $q_0^{2}/|\mathbf q|^{2}\to\infty$). A small change in the long-wavelength density therefore produces a large change in the screened potential, the response $\chi$ is large, and the corresponding eigenvalue of $\mathcal J = \chi K$ is large in magnitude — charge sloshing.
    - **Insulator.** The gap suppresses long-wavelength density fluctuations, so the screening response saturates: $\epsilon(\mathbf q)\to\epsilon_\infty$, a *finite* constant (the familiar static dielectric constant), as $|\mathbf q|\to 0$. The long-wavelength eigenvalues of $\mathcal J$ stay bounded below unity and naive iteration is comparatively well-behaved.

    This is precisely why the cure for metals (the Kerker preconditioner, §5.5.3) targets the small-$|\mathbf q|$ components: it is there, and only there, that the metallic response misbehaves.

!!! note "Why this step?"
    The Jacobian $\mathcal J$ is essentially $\chi v_\mathrm{eff}$, the product of the response function $\chi$ (how the density responds to a perturbation of the KS potential) and the kernel $v_\mathrm{eff}$ (how a perturbation of the density changes the KS potential). For metals, $\chi$ diverges in the long-wavelength limit, making $|\mathcal J|$ large there. The cure is to multiply the residual by an operator that *suppresses* the long-wavelength components — this is the *Kerker preconditioner*, which we meet below.

    !!! warning "Notation: the kernel is $K$, not $v_\mathrm{eff}$"
        The line above writes the second factor of $\mathcal J = \chi\,(\cdot)$ as $v_\mathrm{eff}$, but that symbol is overloaded: in the loop diagram $v_\mathrm{eff}$ is the *KS potential itself*. The factor that actually appears in $\mathcal J = \chi\,(\delta v_\mathrm{KS}/\delta n)$ is the *functional derivative* of that potential with respect to the density, which we write $K \equiv \delta v_\mathrm{KS}/\delta n$ (the Hartree-plus-xc kernel, derived just above). Read "$\mathcal J = \chi K$" throughout.

Convergence of fixed-point iteration requires that the Jacobian (the *dielectric response* of the system) have spectral radius below unity in some norm. For metallic or polarisable systems it typically does not. Physically: if the density at one iteration has slightly too much charge in region $A$, the new Hartree potential pushes electrons out of $A$. The screening response sends *more* charge out of $A$ than the original excess — overshoot — and the next iteration has too little charge in $A$. The system oscillates with growing amplitude. This is **charge sloshing**.

The cure is to dampen the iteration: take only a fraction of the new density and combine it with the previous one. This is *linear mixing*:

$$
n^{(k+1)} = (1-\alpha)\,n^{(k)} + \alpha\,n_\mathrm{out}^{(k)},
\qquad \alpha \in (0,1].
\tag{5.42}
$$

Linearising near $n^{*}$, the *effective* iteration matrix becomes $(1-\alpha)\mathbf I + \alpha\,\mathcal J$ with eigenvalues $(1-\alpha) + \alpha\lambda_i$. For an unstable mode with $\lambda_i = -|\lambda|$ (typical of charge-sloshing modes, since the screening response is anti-correlated with the input perturbation), the linear-mixed eigenvalue is $(1-\alpha)(1) - \alpha|\lambda| = 1 - \alpha(1 + |\lambda|)$. Convergence ($|1-\alpha(1+|\lambda|)|<1$) requires $\alpha < 2/(1+|\lambda|)$ — a much stricter bound than $\alpha<1$ for the simple case. For $|\lambda|\sim 5$ (a typical metallic value), one needs $\alpha\lesssim 0.3$, and for $|\lambda|\sim 10$, $\alpha\lesssim 0.18$. Far from the fixed point, the linear analysis is only suggestive and one usually needs even smaller $\alpha$.

To see where this single bound comes from, unpack the modulus inequality $|1-\alpha(1+|\lambda|)|<1$ into its two branches: it holds iff $-1 < 1-\alpha(1+|\lambda|) < 1$. The *upper* branch, $1-\alpha(1+|\lambda|) < 1$, rearranges to $-\alpha(1+|\lambda|)<0$, i.e. $\alpha(1+|\lambda|)>0$, which is satisfied automatically for any $\alpha>0$ (since $1+|\lambda|>0$) — so it imposes no constraint. The *lower* branch, $1-\alpha(1+|\lambda|) > -1$, rearranges to $\alpha(1+|\lambda|) < 2$, i.e. $\alpha < 2/(1+|\lambda|)$ — this is the binding one, and it is the bound quoted.

Small $\alpha$ (e.g., $\alpha = 0.1$) almost always converges but does so slowly — convergence rate scales as $1 - \alpha$ per iteration. Large $\alpha$ converges fast when it converges and oscillates when it does not. For typical insulators $\alpha = 0.3$ is reasonable; for metals one often needs $\alpha = 0.05$.

!!! example "Convergence rates: insulator vs. metal"
    For a small silicon cluster (insulator, gap $\sim 1\;\text{eV}$), linear mixing with $\alpha=0.3$ converges in $\sim 15$–$20$ iterations to $|\Delta n|_\infty < 10^{-6}$. For bcc iron (ferromagnetic metal), the same $\alpha=0.3$ either diverges or oscillates; $\alpha=0.05$ converges in $\sim 80$–$120$ iterations. With Pulay mixing of history length 8 and Kerker preconditioning, the iron calculation converges in $\sim 25$ iterations. The factor of $\sim 5$ speed-up is typical and the reason every production code uses an acceleration scheme.

Linear mixing is robust but slow. Modern codes use **acceleration schemes** based on the history of recent densities.

??? question "Pause and recall"
    Before reading on, try to answer these from memory:

    1. Why are the Kohn–Sham equations nonlinear, and what is the chicken-and-egg dependence that forces an iterative solution?
    2. What is "charge sloshing", and what does it imply about the eigenvalues of the SCF Jacobian for a metal?
    3. How does linear mixing $n^{(k+1)} = (1-\alpha)n^{(k)} + \alpha\,n_\mathrm{out}^{(k)}$ stabilise the iteration, and what is the trade-off in choosing $\alpha$?

    If any of these is shaky, re-read the preceding section before continuing.

## 5.5.3 Pulay / DIIS mixing

The **Direct Inversion in the Iterative Subspace** (DIIS) method, due to Péter Pulay (1980), is a powerful workhorse. The idea: at iteration $k$, we have a history $\{n^{(j)}, n_\mathrm{out}^{(j)}\}_{j=k-m+1}^{k}$ of $m$ recent inputs and outputs. Define the **residual** of each:

$$
r^{(j)} \equiv n_\mathrm{out}^{(j)} - n^{(j)}.
$$

At self-consistency, $r = 0$. Search for the linear combination

$$
\bar n = \sum_j c_j\, n^{(j)},
\qquad \sum_j c_j = 1,
$$

that minimises the residual norm $\|\sum_j c_j r^{(j)}\|^{2}$. The solution is obtained by setting up the matrix

$$
B_{jk} = \langle r^{(j)}|r^{(k)}\rangle
$$

(here $\langle\cdot|\cdot\rangle$ is the $L^{2}$ inner product on the real-space grid) and solving the constrained linear system

$$
\begin{pmatrix} B & \mathbf 1 \\ \mathbf 1^{T} & 0 \end{pmatrix}
\begin{pmatrix} \mathbf c \\ -\lambda \end{pmatrix}
= \begin{pmatrix} \mathbf 0 \\ 1 \end{pmatrix},
\tag{5.43}
$$

where $\lambda$ is the Lagrange multiplier for the constraint $\sum c_j = 1$. The next-iteration density is

??? note "Full derivation: the bordered Pulay system (5.43) by constrained minimisation"
    We want the coefficients $\{c_j\}$ that make the combined residual as small as possible while keeping the combined density a genuine *average* (the constraint $\sum_j c_j = 1$ ensures that if every $n^{(j)}$ already equalled the fixed point, so would $\bar n$). Write the combined residual using the residuals stored in the history,
    $$
    \bar r = \sum_j c_j\, r^{(j)},
    $$
    and form its squared $L^{2}$ norm. Expanding the square and using bilinearity of the inner product,
    $$
    \big\|\bar r\big\|^{2}
    = \Big\langle \sum_j c_j r^{(j)} \,\Big|\, \sum_k c_k r^{(k)} \Big\rangle
    = \sum_{j}\sum_{k} c_j c_k \langle r^{(j)} | r^{(k)} \rangle
    = \sum_{j}\sum_{k} c_j\, B_{jk}\, c_k
    = \mathbf c^{T} B\,\mathbf c,
    \tag{5.43c}
    $$
    where $B_{jk} = \langle r^{(j)} | r^{(k)}\rangle$ is exactly the overlap (Gram) matrix introduced above. Note $B$ is symmetric ($B_{jk}=B_{kj}$) and positive semi-definite, so $\mathbf c^{T}B\,\mathbf c \ge 0$ is a sensible thing to minimise.

    **Impose the constraint with a Lagrange multiplier.** To minimise $\mathbf c^{T}B\,\mathbf c$ subject to $\mathbf 1^{T}\mathbf c = 1$ (where $\mathbf 1 = (1,1,\dots,1)^{T}$), form the Lagrangian
    $$
    \mathcal L(\mathbf c, \lambda)
    = \mathbf c^{T} B\,\mathbf c - 2\lambda\big(\mathbf 1^{T}\mathbf c - 1\big).
    $$
    The factor $2$ on the multiplier is a harmless convention chosen to make the final equations tidy (it just rescales $\lambda$). Differentiate with respect to the coefficient vector. Using $\partial_{\mathbf c}(\mathbf c^{T}B\,\mathbf c) = 2B\,\mathbf c$ (valid because $B$ is symmetric) and $\partial_{\mathbf c}(\mathbf 1^{T}\mathbf c) = \mathbf 1$,
    $$
    \frac{\partial \mathcal L}{\partial \mathbf c}
    = 2B\,\mathbf c - 2\lambda\,\mathbf 1 = \mathbf 0
    \quad\Longrightarrow\quad
    B\,\mathbf c = \lambda\,\mathbf 1.
    \tag{5.43d}
    $$
    Differentiating with respect to $\lambda$ simply returns the constraint,
    $$
    \frac{\partial \mathcal L}{\partial \lambda}
    = -2\big(\mathbf 1^{T}\mathbf c - 1\big) = 0
    \quad\Longrightarrow\quad
    \mathbf 1^{T}\mathbf c = 1.
    \tag{5.43e}
    $$

    **Assemble into one matrix.** Equations (5.43d) and (5.43e) are a coupled linear system in the unknowns $(\mathbf c, \lambda)$. Stack them. The first, $B\,\mathbf c - \lambda\,\mathbf 1 = \mathbf 0$, becomes the top block-row; the second, $\mathbf 1^{T}\mathbf c = 1$, becomes the bottom row. Writing the multiplier slot as $-\lambda$ so the signs line up with (5.43),
    $$
    \begin{pmatrix} B & \mathbf 1 \\ \mathbf 1^{T} & 0 \end{pmatrix}
    \begin{pmatrix} \mathbf c \\ -\lambda \end{pmatrix}
    = \begin{pmatrix} B\,\mathbf c - \lambda\,\mathbf 1 \\ \mathbf 1^{T}\mathbf c \end{pmatrix}
    = \begin{pmatrix} \mathbf 0 \\ 1 \end{pmatrix},
    $$
    which is exactly the bordered system (5.43). The top block reproduces (5.43d) and the bottom row reproduces the constraint (5.43e). Solving this one $(m{+}1)\times(m{+}1)$ system yields the optimal Pulay weights $\mathbf c$ in a single linear solve; the multiplier $\lambda$ falls out as a by-product and equals the minimised residual norm, $\lambda = \mathbf c^{T}B\,\mathbf c = \|\bar r\|^{2}$ (multiply (5.43d) on the left by $\mathbf c^{T}$ and use $\mathbf 1^{T}\mathbf c = 1$).

$$
n^{(k+1)} = \sum_j c_j\,n_\mathrm{out}^{(j)}
\qquad(\text{or}\;\;\sum_j c_j\,n^{(j)} + \alpha\sum_j c_j r^{(j)},\text{ DIIS with damping}).
\tag{5.44}
$$

The two forms in (5.44) are connected by the definition of the residual. Since $r^{(j)} = n_\mathrm{out}^{(j)} - n^{(j)}$, we have $n_\mathrm{out}^{(j)} = n^{(j)} + r^{(j)}$, and substituting into the first form,
$$
n^{(k+1)} = \sum_j c_j\,n_\mathrm{out}^{(j)} = \sum_j c_j\big(n^{(j)} + r^{(j)}\big) = \sum_j c_j\,n^{(j)} + \sum_j c_j\,r^{(j)}.
\tag{5.44a}
$$
This separates the next density into an *extrapolated input* $\sum_j c_j n^{(j)}$ plus the *combined residual* $\sum_j c_j r^{(j)}$. The undamped scheme takes the residual term at full strength (coefficient $1$); the damped variant rescales just that residual term by $\alpha$, giving $\sum_j c_j n^{(j)} + \alpha\sum_j c_j r^{(j)}$ — the second form quoted in (5.44). Setting $\alpha=1$ recovers the undamped scheme; taking $\alpha<1$ adds back a touch of the conservatism of linear mixing onto the optimal Pulay direction, which helps when the history is not yet reliable.

In practice DIIS converges much faster than linear mixing — often quadratically near the solution. It needs a small history (typically $m = 6$–$10$). Far from convergence DIIS can be unstable; codes typically start with several linear-mixing steps before switching to DIIS, or fall back to linear mixing if DIIS diverges.

!!! note "Why this step?"
    DIIS works because, near a fixed point, the SCF map is approximately linear: $r^{(k+1)}\approx \mathcal J\,r^{(k)}$. A linear combination of past inputs $\bar n = \sum c_j n^{(j)}$ has a residual $\bar r = \sum c_j r^{(j)}$ (by linearity of $\mathcal F - \mathbf I$ when $\mathcal F$ is linear). Minimising $\|\bar r\|^{2}$ over $\{c_j\}$ subject to $\sum c_j = 1$ finds the linear combination *most consistent with zero residual* — essentially a Krylov-subspace projection of the fixed-point equation. The result is faster than any single-step linear mixing because it uses information from $m$ past iterates simultaneously, effectively building a low-rank approximation to the inverse Jacobian on the fly.

    To make the "quasi-Newton" reading concrete: a Newton step for the root of the residual map $r(n) = 0$ would be $n^{(k+1)} = n^{(k)} - \tilde{\mathcal J}^{-1} r^{(k)}$, where $\tilde{\mathcal J} = \delta r/\delta n = \mathcal J - \mathbf I$ is the Jacobian of the *residual* (not of the SCF map itself). The exact $\tilde{\mathcal J}^{-1}$ is unavailable — that is the whole difficulty — so all the acceleration schemes here replace it by a cheap approximation assembled from the history: DIIS/Anderson build it implicitly through the least-squares weights $\{c_j\}$, while Broyden (below) maintains it explicitly as the matrix $G\approx\tilde{\mathcal J}^{-1}$ and takes literally the step $n^{(k+1)} = n^{(k)} - \tilde{\mathcal J}^{-1} r^{(k)}$ with $\tilde{\mathcal J}^{-1}\to G$. This is what "quasi-Newton" means: a Newton iteration with an approximate, history-built inverse Jacobian in place of the true one.

### Worked example: Pulay applied by hand

Suppose three past iterations gave residuals $r^{(1)} = [1, 0, 0]$, $r^{(2)} = [0, 1, 0]$, $r^{(3)} = [-0.5, -0.5, 0.1]$ (artificially small for illustration). The overlap matrix is

$$
B = \begin{pmatrix} 1 & 0 & -0.5 \\ 0 & 1 & -0.5 \\ -0.5 & -0.5 & 0.51 \end{pmatrix}.
$$

The constrained system (5.43) with this $B$ and right-hand side $(0, 0, 0, 1)$ gives $\mathbf c = (0.2512, 0.2512, 0.4975)^{T}$ (derived step by step below), i.e. the optimal combination leans most heavily on the most recent iterate (which has the smallest residual). The new density is then $n^{(4)} = 0.2512\,n_\mathrm{out}^{(1)} + 0.2512\,n_\mathrm{out}^{(2)} + 0.4975\,n_\mathrm{out}^{(3)}$. The interpretation: Pulay automatically discounts old, less-relevant iterates and emphasises the current best guess, much like a momentum optimiser in machine learning.

!!! note "Where these weights come from"
    These weights are not guessed — they are the exact solution $\mathbf c = B^{-1}\mathbf 1/(\mathbf 1^{T}B^{-1}\mathbf 1) = \big(\tfrac{101}{402},\,\tfrac{101}{402},\,\tfrac{100}{201}\big)^{T}$ of the constrained least-squares system (5.43). The combination leans most heavily on the most recent iterate $r^{(3)}$ (about twice the weight of the older two); the step-by-step hand calculation is in the collapsible box below.

??? note "Full derivation: solving the $3\times 3$ Pulay example by hand"
    **Step 1 — build $B$.** The overlap matrix has entries $B_{jk} = \langle r^{(j)} | r^{(k)}\rangle = r^{(j)}\cdot r^{(k)}$ (ordinary dot products, since the residuals here are plain 3-vectors). With $r^{(1)} = [1,0,0]$, $r^{(2)} = [0,1,0]$, $r^{(3)} = [-0.5,-0.5,0.1]$:

    - $B_{11} = 1^2+0+0 = 1$, $\;B_{22} = 0+1^2+0 = 1$;
    - $B_{12} = B_{21} = (1)(0)+(0)(1)+(0)(0) = 0$;
    - $B_{13} = B_{31} = (1)(-0.5)+(0)(-0.5)+(0)(0.1) = -0.5$;
    - $B_{23} = B_{32} = (0)(-0.5)+(1)(-0.5)+(0)(0.1) = -0.5$;
    - $B_{33} = (-0.5)^2 + (-0.5)^2 + (0.1)^2 = 0.25 + 0.25 + 0.01 = 0.51.$

    The last line verifies the $0.51$ entry quoted for the matrix; this reproduces $B$ exactly.

    **Step 2 — reduce the bordered system.** From the derivation of (5.43), the optimality condition is $B\,\mathbf c = \lambda\,\mathbf 1$ together with $\mathbf 1^{T}\mathbf c = 1$. Solve the first relation for $\mathbf c$ in terms of $\lambda$,
    $$
    \mathbf c = \lambda\,B^{-1}\mathbf 1,
    $$
    then fix $\lambda$ by the constraint $\mathbf 1^{T}\mathbf c = \lambda\,(\mathbf 1^{T}B^{-1}\mathbf 1) = 1$, i.e.
    $$
    \lambda = \frac{1}{\mathbf 1^{T}B^{-1}\mathbf 1},
    \qquad
    \mathbf c = \frac{B^{-1}\mathbf 1}{\mathbf 1^{T}B^{-1}\mathbf 1}.
    \tag{5.43f}
    $$
    So we only need the vector $B^{-1}\mathbf 1$, i.e. the solution $\mathbf u$ of $B\,\mathbf u = \mathbf 1$.

    **Step 3 — solve $B\,\mathbf u = \mathbf 1$.** Writing $\mathbf u = (u_1,u_2,u_3)^{T}$, the three equations are
    $$
    \begin{aligned}
    u_1 - 0.5\,u_3 &= 1, \\
    u_2 - 0.5\,u_3 &= 1, \\
    -0.5\,u_1 - 0.5\,u_2 + 0.51\,u_3 &= 1.
    \end{aligned}
    $$
    The first two rows are identical in structure, so $u_1 = u_2 = 1 + 0.5\,u_3$. Substitute into the third row:
    $$
    -0.5\,(1+0.5u_3) - 0.5\,(1+0.5u_3) + 0.51\,u_3 = 1.
    $$
    Expand: $-0.5 - 0.25u_3 - 0.5 - 0.25u_3 + 0.51u_3 = 1$, i.e. $-1 + (-0.25-0.25+0.51)u_3 = 1$, i.e. $-1 + 0.01\,u_3 = 1$. Hence
    $$
    0.01\,u_3 = 2 \;\Longrightarrow\; u_3 = 200,
    \qquad
    u_1 = u_2 = 1 + 0.5(200) = 101.
    $$
    So $B^{-1}\mathbf 1 = \mathbf u = (101,\,101,\,200)^{T}$.

    **Step 4 — normalise.** The denominator in (5.43f) is $\mathbf 1^{T}\mathbf u = 101 + 101 + 200 = 402$, so $\lambda = 1/402 \approx 0.00249$ and
    $$
    \mathbf c = \frac{1}{402}\,(101,\,101,\,200)^{T}
    = \Big(\tfrac{101}{402},\,\tfrac{101}{402},\,\tfrac{100}{201}\Big)^{T}
    \approx (0.2512,\,0.2512,\,0.4975)^{T}.
    $$
    As a check, the components sum to $(101+101+200)/402 = 402/402 = 1$, satisfying the constraint exactly. The most recent iterate $r^{(3)}$ (smallest residual) carries the largest weight $\approx 0.4975$, the two older iterates share the rest equally — exactly the behaviour the prose describes, with the corrected numbers. The new density is therefore
    $$
    n^{(4)} = 0.2512\,n_\mathrm{out}^{(1)} + 0.2512\,n_\mathrm{out}^{(2)} + 0.4975\,n_\mathrm{out}^{(3)}.
    $$

### Anderson mixing

Anderson mixing (1965) is an older, closely related method. At iteration $k$ with output residual $r^{(k)} = n_\mathrm{out}^{(k)} - n^{(k)}$, set

$$
n^{(k+1)} = n^{(k)} + \alpha\,r^{(k)} - \beta\big(r^{(k)} - r^{(k-1)}\big),
\tag{5.43a}
$$

with $\alpha,\beta$ chosen to minimise $\|r^{(k+1)}\|$ over the affine subspace. The general $m$-step Anderson method is essentially equivalent to DIIS; many modern codes use this formulation.

To make (5.43a) explicit: define $\Delta r = r^{(k)} - r^{(k-1)}$ and $\Delta n = n^{(k)} - n^{(k-1)}$. The Anderson update can be rewritten as

$$
n^{(k+1)} = (1-\theta)\big(n^{(k)} + \alpha r^{(k)}\big) + \theta\big(n^{(k-1)} + \alpha r^{(k-1)}\big),
$$

where $\theta = \langle\Delta r, r^{(k)}\rangle/\|\Delta r\|^{2}$ is chosen to minimise the linearised residual. Note that for $\theta=0$, this reduces to linear mixing; for general $\theta$, it interpolates between two recent linearly-mixed updates, suppressing the *common* component of the residuals (which is the slow mode) while preserving the orthogonal component.

??? note "Full derivation: the Anderson weight $\theta$"
    Near the fixed point the SCF map is approximately linear, so the residual of an *interpolated* iterate is the same interpolation of the residuals: blending input $(1-\theta)n^{(k)} + \theta n^{(k-1)}$ produces (to linear order) the residual $(1-\theta)r^{(k)} + \theta r^{(k-1)}$. Choose $\theta$ to make that blended residual as small as possible. Write it using $\Delta r = r^{(k)} - r^{(k-1)}$:
    $$
    (1-\theta)r^{(k)} + \theta r^{(k-1)}
    = r^{(k)} - \theta\big(r^{(k)} - r^{(k-1)}\big)
    = r^{(k)} - \theta\,\Delta r.
    $$
    So we minimise the scalar function
    $$
    f(\theta) = \big\| r^{(k)} - \theta\,\Delta r \big\|^{2}
    = \langle r^{(k)} - \theta\Delta r \,|\, r^{(k)} - \theta\Delta r\rangle
    = \|r^{(k)}\|^{2} - 2\theta\,\langle \Delta r, r^{(k)}\rangle + \theta^{2}\|\Delta r\|^{2},
    $$
    a simple upward parabola in $\theta$ (coefficient $\|\Delta r\|^{2} \ge 0$). Set the derivative to zero:
    $$
    f'(\theta) = -2\,\langle \Delta r, r^{(k)}\rangle + 2\theta\,\|\Delta r\|^{2} = 0
    \quad\Longrightarrow\quad
    \theta = \frac{\langle \Delta r, r^{(k)}\rangle}{\|\Delta r\|^{2}},
    $$
    which is the quoted weight. (The second derivative $f''(\theta) = 2\|\Delta r\|^{2} > 0$ confirms a minimum.) Geometrically, $\theta\,\Delta r$ is the orthogonal projection of $r^{(k)}$ onto the direction $\Delta r$, so $r^{(k)} - \theta\Delta r$ is the component of $r^{(k)}$ *orthogonal* to the change in residual — Anderson removes exactly the part of the residual that the last step was able to move.

### Broyden's second method

Broyden mixing is a quasi-Newton scheme that approximates the inverse Jacobian of the SCF map from the iteration history. It generalises both Anderson and DIIS and is the default in some plane-wave codes (e.g., VASP). The implementation is more involved but the convergence behaviour is similar to DIIS for most problems.

### Broyden, Kerker, and other useful tricks

**Broyden's second method** maintains an approximate inverse Jacobian $G^{(k)}$ that is updated each iteration to satisfy the secant equation $G^{(k)}\Delta r = \Delta n$. The update is

$$
G^{(k+1)} = G^{(k)} + \frac{(\Delta n - G^{(k)}\Delta r)\,\Delta r^{T}}{\|\Delta r\|^{2}},
$$

and the SCF step is $n^{(k+1)} = n^{(k)} - G^{(k+1)}r^{(k)}$ (the minus sign is the Newton step for the root of the residual map, consistent with the quasi-Newton reading above). In practice $G$ is stored in low-rank form (a few past $(\Delta n_j, \Delta r_j)$ pairs); this is the *limited-memory Broyden* method used by VASP.

??? note "Full derivation: the Broyden update satisfies the secant equation"
    Broyden's second method builds an approximate *inverse* Jacobian $G \approx \mathcal J^{-1}$ of the SCF map directly from the iteration history. The single piece of new information at each step is one input–output pair: the density moved by $\Delta n = n^{(k)} - n^{(k-1)}$ and, correspondingly, the residual moved by $\Delta r = r^{(k)} - r^{(k-1)}$. A true inverse Jacobian would map the second change back to the first; we *demand* the same of our approximation. That requirement is the **secant equation**
    $$
    G^{(k+1)}\,\Delta r = \Delta n.
    \tag{5.45a}
    $$

    **The rank-1 update.** The minimal-norm update is the *smallest* change to $G^{(k)}$ (in the Frobenius norm $\|G^{(k+1)}-G^{(k)}\|_F$) that satisfies (5.45a). It has the rank-1 form
    $$
    G^{(k+1)} = G^{(k)} + \frac{\big(\Delta n - G^{(k)}\Delta r\big)\,\Delta r^{T}}{\|\Delta r\|^{2}},
    \tag{5.45b}
    $$
    i.e. the correction vector $\Delta n - G^{(k)}\Delta r$ (how badly the *current* $G^{(k)}$ violates the secant condition) times the row vector $\Delta r^{T}/\|\Delta r\|^{2}$. (This is the form to use; if the second factor is written with $G^{(k)}$ acting on $\Delta r$ instead of the bare $\Delta r^{T}$, the verification below does not close.)

    **Verify (5.45a) by direct substitution.** Multiply (5.45b) on the right by $\Delta r$:
    $$
    G^{(k+1)}\Delta r
    = G^{(k)}\Delta r
    + \frac{\big(\Delta n - G^{(k)}\Delta r\big)\,\big(\Delta r^{T}\Delta r\big)}{\|\Delta r\|^{2}}.
    $$
    Now $\Delta r^{T}\Delta r = \|\Delta r\|^{2}$ is just the squared norm in the denominator, so the fraction collapses:
    $$
    G^{(k+1)}\Delta r
    = G^{(k)}\Delta r + \big(\Delta n - G^{(k)}\Delta r\big)\,\frac{\|\Delta r\|^{2}}{\|\Delta r\|^{2}}
    = G^{(k)}\Delta r + \Delta n - G^{(k)}\Delta r
    = \Delta n.
    $$
    The two $G^{(k)}\Delta r$ terms cancel and the secant equation (5.45a) holds exactly — the update *forces* the new $G$ to reproduce the most recent density/residual pair, while the minimal-Frobenius-norm property (the **minimal-norm-update principle**, also called the "least-change" or Broyden update) guarantees it disturbs the rest of $G$ as little as possible, so information accumulated from earlier steps is preserved.

**Kerker preconditioning.** For metallic systems, the SCF Jacobian has eigenvalues that diverge in the long-wavelength limit. The Kerker preconditioner replaces the residual $r$ by a modified residual $\tilde r$ via the reciprocal-space rule

$$
\tilde r(\mathbf q) = \frac{|\mathbf q|^{2}}{|\mathbf q|^{2} + q_0^{2}}\,r(\mathbf q),
\tag{5.43b}
$$

where $q_0$ is a tunable scale (typically $q_0\sim 1.5\;\text{Bohr}^{-1}$ for typical metals). This kills the long-wavelength components that drive charge sloshing while leaving short-wavelength components essentially unchanged. The mixing scheme then operates on $\tilde r$ rather than $r$, and the effective Jacobian has eigenvalues much closer to 1 in magnitude.

??? note "Full derivation: the Kerker factor is the inverse Thomas–Fermi dielectric"
    Start from the Thomas–Fermi dielectric function of a uniform electron gas (the same (5.41e) used above to contrast metals and insulators),
    $$
    \epsilon_\mathrm{TF}(\mathbf q) = 1 + \frac{q_0^{2}}{|\mathbf q|^{2}},
    $$
    where $q_0$ is the **Thomas–Fermi screening wavevector**: $q_0^{2} \propto g(\varepsilon_F)$, the density of states at the Fermi level, so a metal (large $g(\varepsilon_F)$) screens strongly and an insulator (vanishing $g(\varepsilon_F)$) does not. The trouble for SCF is precisely that $\epsilon_\mathrm{TF}\to\infty$ as $|\mathbf q|\to 0$: the long-wavelength response is huge, the corresponding Jacobian eigenvalues blow up, and the iteration sloshes.

    **Invert it.** The natural preconditioner is the *inverse* dielectric, which undoes that divergence. Take the reciprocal and put the two terms over a common denominator:
    $$
    \frac{1}{\epsilon_\mathrm{TF}(\mathbf q)}
    = \frac{1}{\,1 + \dfrac{q_0^{2}}{|\mathbf q|^{2}}\,}
    = \frac{1}{\dfrac{|\mathbf q|^{2} + q_0^{2}}{|\mathbf q|^{2}}}
    = \frac{|\mathbf q|^{2}}{|\mathbf q|^{2} + q_0^{2}}.
    $$
    This is exactly the Kerker factor in (5.43b). Reading off its two limits shows why it cures sloshing without harming the well-behaved modes:

    - **Long wavelength**, $|\mathbf q|\to 0$: the factor $\to |\mathbf q|^{2}/q_0^{2}\to 0$. The divergent small-$\mathbf q$ components of the residual — the charge-sloshing modes — are strongly damped before mixing.
    - **Short wavelength**, $|\mathbf q|\gg q_0$: the factor $\to |\mathbf q|^{2}/|\mathbf q|^{2} = 1$. The large-$\mathbf q$ components, which were already well-behaved, pass through essentially unchanged.

    Because applying $1/\epsilon_\mathrm{TF}$ approximately cancels the leading small-$\mathbf q$ divergence of $\mathcal J = \chi K$, the preconditioned residual $\tilde r$ behaves as if the effective Jacobian had eigenvalues near unity across all wavevectors — which is the condition for fast, stable iteration. The parameter $q_0$ is identified as the Thomas–Fermi screening wavevector and is tuned (typically $\sim 1.5\,\text{Bohr}^{-1}$) so the crossover between the two regimes sits below the shortest sloshing wavelength of the system.

!!! note "Why this step?"
    The Kerker factor $|\mathbf q|^{2}/(|\mathbf q|^{2}+q_0^{2})$ is precisely the inverse dielectric response of a Thomas–Fermi electron gas, scaled by $q_0^{-2}$. The Kerker preconditioner is therefore an *approximate inverse* of the SCF Jacobian: applying it cancels the leading divergence of $\mathcal J$ at small $|\mathbf q|$. This is the same idea as preconditioning a linear system in numerical analysis — multiply by an approximation to the inverse of the offending operator.

## 5.5.4 Convergence criteria

Several quantities can be monitored:

- **Density residual** $\|n_\mathrm{out} - n_\mathrm{in}\| = \int|n_\mathrm{out} - n_\mathrm{in}|^{2}\,\mathrm d\mathbf r$ or $\max|n_\mathrm{out} - n_\mathrm{in}|$. The most rigorous criterion.
- **Energy difference** $|E^{(k+1)} - E^{(k)}|$. Easy to compute; typical tolerance $10^{-5}$ to $10^{-8}$ Ha. Beware: small energy changes do not always mean converged densities.
- **Force / stress changes**: critical for geometry optimisations. Want forces converged to $\sim 10^{-3}$ Ha/Bohr or better before trusting them.

!!! note "Notation: the density residual norm carries a square root"
    The first bullet writes $\|n_\mathrm{out} - n_\mathrm{in}\| = \int|n_\mathrm{out} - n_\mathrm{in}|^{2}\,\mathrm d\mathbf r$, but strictly the $L^{2}$ norm carries a square root,
    $$
    \|f\|_2 = \Big(\int|f(\mathbf r)|^{2}\,\mathrm d\mathbf r\Big)^{1/2}.
    $$
    The expression written without the root is therefore the *squared* norm $\|n_\mathrm{out}-n_\mathrm{in}\|_2^{2}$, not $\|n_\mathrm{out}-n_\mathrm{in}\|_2$. Either may be used as a convergence measure provided the threshold is set consistently (squaring merely doubles the apparent number of converged decimals per step). The code in §5.5.5 sidesteps the issue entirely by monitoring the max-norm $\max|n_\mathrm{out}-n_\mathrm{in}|$, which needs no root.

A common protocol: require *both* an energy tolerance and a density tolerance to be satisfied for two consecutive iterations.

### Recommended thresholds by application

| Application | Energy tol. (Ha) | Density tol. | Force tol. (Ha/Bohr) |
|---|---|---|---|
| Single-point energy (production) | $10^{-6}$ | $10^{-5}$ | n/a |
| Geometry optimisation | $10^{-7}$ | $10^{-6}$ | $5\times 10^{-4}$ |
| Phonons / elastic constants | $10^{-9}$ | $10^{-8}$ | $10^{-5}$ |
| Molecular dynamics (BOMD) | $10^{-6}$ | $10^{-5}$ | $10^{-4}$ |
| ML potential training data | $10^{-7}$ | $10^{-6}$ | $10^{-4}$ |

Phonons and elastic constants require *very* tight SCF because they are computed from second derivatives of the energy — small SCF noise gets amplified by numerical differentiation. ML training-data generation is intermediate: the model can absorb some SCF noise during training, but biased noise (e.g., systematically under-converged forces on certain structure types) propagates into systematic model error.

## 5.5.5 A complete Python implementation

We now solve the Kohn–Sham equations for a one-dimensional "hydrogen chain" model: $N_\mathrm{atom}$ protons placed on a line, with the electron–nucleus interaction softened to avoid the 1D Coulomb singularity, periodic boundary conditions on a finite box, and LDA exchange. We use a real-space finite-difference discretisation and direct diagonalisation. The code is ~150 lines, fully type-hinted, runnable on a laptop in under a second per iteration.

```python
"""scf_1d_hchain.py
A minimal Kohn--Sham SCF solver for a 1D hydrogen chain using LDA exchange.

Real-space finite-difference Hamiltonian on a uniform grid with periodic
boundary conditions. Linear and Pulay/DIIS mixing.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigh, solve
import matplotlib.pyplot as plt


@dataclass
class Grid:
    """Uniform 1D real-space grid with periodic boundary conditions."""
    n: int                # number of grid points
    L: float              # box length (Bohr)

    @property
    def dx(self) -> float:
        return self.L / self.n

    @property
    def x(self) -> NDArray[np.float64]:
        return np.arange(self.n) * self.dx


def kinetic_matrix(g: Grid) -> NDArray[np.float64]:
    """Second-order central finite-difference kinetic operator -1/2 d^2/dx^2."""
    n, dx = g.n, g.dx
    main = np.full(n, 1.0 / dx ** 2)
    off = np.full(n - 1, -0.5 / dx ** 2)
    T = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    # periodic boundary conditions
    T[0, -1] = T[-1, 0] = -0.5 / dx ** 2
    return T


def soft_coulomb(x: NDArray[np.float64], x0: float, L: float,
                 a: float = 1.0) -> NDArray[np.float64]:
    """Soft-Coulomb e-N attraction, periodic in box of length L.
    v(x) = -1/sqrt((x-x0)^2 + a^2), summed over periodic images."""
    v = np.zeros_like(x)
    for m in range(-2, 3):  # nearest images suffice
        dx_arr = x - x0 - m * L
        v += -1.0 / np.sqrt(dx_arr ** 2 + a ** 2)
    return v


def external_potential(g: Grid, positions: list[float]) -> NDArray[np.float64]:
    """Total e-N potential for protons at given positions."""
    v = np.zeros(g.n)
    for x0 in positions:
        v += soft_coulomb(g.x, x0, g.L)
    return v


def hartree_potential(g: Grid, n: NDArray[np.float64],
                      a: float = 1.0) -> NDArray[np.float64]:
    """Hartree potential from density n(x), using soft Coulomb kernel."""
    vH = np.zeros(g.n)
    dx = g.dx
    for i, xi in enumerate(g.x):
        # sum over periodic images
        contrib = 0.0
        for m in range(-2, 3):
            d = g.x - xi - m * g.L
            contrib += np.sum(n / np.sqrt(d ** 2 + a ** 2)) * dx
        vH[i] = contrib
    return vH


def lda_exchange_potential(n: NDArray[np.float64]) -> NDArray[np.float64]:
    """LDA exchange potential, 3D form applied to 1D effective density.
    v_x = -(3/pi)^{1/3} n^{1/3}. Pedagogical, not 1D-exact."""
    return -((3.0 / np.pi) ** (1.0 / 3.0)) * np.cbrt(np.maximum(n, 1e-12))


def lda_exchange_energy(n: NDArray[np.float64], dx: float) -> float:
    """LDA exchange energy per (5.36), pedagogical 3D form."""
    cx = -0.75 * (3.0 / np.pi) ** (1.0 / 3.0)
    return cx * float(np.sum(n ** (4.0 / 3.0)) * dx)


def build_density(orbitals: NDArray[np.float64], n_occ: int,
                  dx: float) -> NDArray[np.float64]:
    """Density from doubly occupied lowest orbitals; orbitals are columns."""
    psi = orbitals[:, :n_occ]
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2, axis=0) * dx)  # normalise
    return 2.0 * np.sum(np.abs(psi) ** 2, axis=1)


def total_energy(eigvals: NDArray[np.float64], n: NDArray[np.float64],
                 vH: NDArray[np.float64], v_ext: NDArray[np.float64],
                 dx: float, n_occ: int) -> float:
    """Kohn--Sham total energy from band-sum (5.28) form, LDA-X only."""
    band = 2.0 * float(np.sum(eigvals[:n_occ]))
    # double-counting corrections
    eH = 0.5 * float(np.sum(n * vH) * dx)
    vx = lda_exchange_potential(n)
    ex = lda_exchange_energy(n, dx)
    e_dc = -eH + ex - float(np.sum(n * vx) * dx)
    return band + e_dc


def pulay_mix(n_in_hist: list[NDArray[np.float64]],
              n_out_hist: list[NDArray[np.float64]],
              m: int = 6) -> NDArray[np.float64]:
    """Pulay/DIIS mixing using last m iterates."""
    k = len(n_in_hist)
    use = min(k, m)
    inputs = np.array(n_in_hist[-use:])
    outputs = np.array(n_out_hist[-use:])
    res = outputs - inputs  # residuals
    B = res @ res.T
    A = np.zeros((use + 1, use + 1))
    A[:use, :use] = B
    A[use, :use] = 1.0
    A[:use, use] = 1.0
    rhs = np.zeros(use + 1)
    rhs[-1] = 1.0
    try:
        sol = solve(A, rhs)
    except np.linalg.LinAlgError:
        return 0.7 * n_in_hist[-1] + 0.3 * n_out_hist[-1]
    coeffs = sol[:use]
    return coeffs @ outputs


def scf(positions: list[float], n_electrons: int, grid: Grid,
        alpha: float = 0.3, tol: float = 1e-6, max_iter: int = 100,
        use_pulay_after: int = 3) -> dict:
    """Run the SCF loop. Returns dict with density, eigvals, energies, etc."""
    n_occ = n_electrons // 2
    T = kinetic_matrix(grid)
    v_ext = external_potential(grid, positions)
    # initial guess: uniform density
    n = np.full(grid.n, n_electrons / grid.L)
    n_in_hist: list[NDArray[np.float64]] = []
    n_out_hist: list[NDArray[np.float64]] = []
    energies: list[float] = []
    for it in range(max_iter):
        vH = hartree_potential(grid, n)
        vx = lda_exchange_potential(n)
        v_ks = v_ext + vH + vx
        H = T + np.diag(v_ks)
        eigvals, eigvecs = eigh(H)
        n_out = build_density(eigvecs, n_occ, grid.dx)
        # renormalise to N electrons (FD discretisation drifts slightly)
        n_out *= n_electrons / (np.sum(n_out) * grid.dx)
        E = total_energy(eigvals, n_out, vH, v_ext, grid.dx, n_occ)
        energies.append(E)
        residual = float(np.max(np.abs(n_out - n)))
        print(f"iter {it:3d}  E = {E: .6f} Ha   |dn|_inf = {residual:.2e}")
        if residual < tol:
            print(f"Converged in {it+1} iterations.")
            break
        n_in_hist.append(n.copy())
        n_out_hist.append(n_out.copy())
        if it < use_pulay_after:
            n = (1 - alpha) * n + alpha * n_out
        else:
            n = pulay_mix(n_in_hist, n_out_hist, m=6)
    return {"density": n, "x": grid.x, "eigvals": eigvals,
            "energies": energies, "v_ks": v_ks}


def main() -> None:
    grid = Grid(n=256, L=20.0)
    positions = [4.0, 8.0, 12.0, 16.0]   # H4 chain, ~4 Bohr spacing
    n_electrons = 4
    result = scf(positions, n_electrons, grid, alpha=0.3, tol=1e-6,
                 max_iter=80)
    fig, axes = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    axes[0].plot(result["x"], result["density"], lw=2)
    for x0 in positions:
        axes[0].axvline(x0, color="grey", ls=":", alpha=0.5)
    axes[0].set_ylabel(r"$n(x)$ (Bohr$^{-1}$)")
    axes[0].set_title("Converged Kohn-Sham density, 1D H$_4$ chain (LDA-X)")
    axes[1].plot(result["x"], result["v_ks"], lw=2)
    for x0 in positions:
        axes[1].axvline(x0, color="grey", ls=":", alpha=0.5)
    axes[1].set_xlabel("x (Bohr)")
    axes[1].set_ylabel(r"$v_\mathrm{KS}(x)$ (Ha)")
    fig.tight_layout()
    fig.savefig("scf_1d_hchain.png", dpi=150)
    print("Eigenvalues (Ha):", result["eigvals"][:6])


if __name__ == "__main__":
    main()
```

### Symbol guide: maths ↔ code

Before reading the prose below, use this table to line up the equations in this section with the variable names in the listing.

| Maths | Code variable | Meaning |
|---|---|---|
| $n^{(k)}$ (input density) | `n` (inside the loop) | the density fed *in* this iteration |
| $n_\mathrm{out}^{(k)}$ (output density) | `n_out` | the NEW density read off the orbitals |
| $v_\mathrm{ext}$ | `v_ext` | the fixed electron–nucleus potential |
| $v_H[n]$ | `vH` | Hartree potential built from `n` |
| $v_{xc}[n]$ (here LDA exchange only) | `vx` | exchange–correlation potential built from `n` |
| $v_\mathrm{KS}$ | `v_ks` | total Kohn–Sham potential $=$ `v_ext + vH + vx` |
| $\hat H_\mathrm{KS}$, $\{\varepsilon_i,\phi_i\}$ | `H`, `eigvals`, `eigvecs` | the KS Hamiltonian and its eigenpairs |
| $\alpha$ (mixing parameter) | `alpha` | fraction of `n_out` blended in by linear mixing |
| $r^{(k)} = n_\mathrm{out}-n$ (residual) | `residual` (and `res` inside `pulay_mix`) | how far from self-consistency we are |
| tolerance | `tol` | stop once `residual < tol` |
| $c_j$ (Pulay weights) | `coeffs` | optimal weights on the history of densities |

The single most important line to find is `residual = float(np.max(np.abs(n_out - n)))`: this is $\max|n_\mathrm{out}-n^{(k)}|$, the quantity that must fall below `tol` for the loop to stop.

### How to read this code

- **Grid and operators.** A uniform real-space grid of $n=256$ points on a 20-Bohr box, with periodic boundary conditions. The kinetic operator is the second-order central difference $T = -\tfrac{1}{2}D^{2}$ assembled as a dense matrix; for larger systems one would use a sparse representation.
- **Soft Coulomb.** A 1D Coulomb $-1/|x-x_0|$ is singular; replacing $1/|x|$ with $1/\sqrt{x^{2}+a^{2}}$ regularises it and makes the model physically reasonable. The Hartree kernel uses the same softening.

!!! warning "The toy 'Hartree' is a 1D soft-Coulomb convolution, not the 3D Poisson solution"
    The bird's-eye pseudocode (§5.5.5a) writes the Hartree step as *"solve Poisson: $\nabla^{2}v_H = -4\pi n$"*, which is the correct 3D electrostatics: in three dimensions $v_H(\mathbf r) = \int n(\mathbf r')/|\mathbf r-\mathbf r'|\,\mathrm d\mathbf r'$ is the Green's-function solution of Poisson's equation with the $1/|\mathbf r-\mathbf r'|$ kernel. The runnable `hartree_potential` in this listing does *not* solve a 1D Poisson equation. Instead it evaluates a *convolution* of the density with the same softened kernel used for the electron–nucleus attraction,
    $$
    v_H(x) = \int \frac{n(x')}{\sqrt{(x-x')^{2} + a^{2}}}\,\mathrm d x'
    \quad(\text{summed over periodic images}),
    $$
    coded as the double loop `contrib += sum(n / sqrt(d**2 + a**2)) * dx`. This is a deliberate, pedagogically convenient choice — the true 1D Poisson kernel ($\nabla^{2}v_H = -4\pi n$ in 1D gives a Green's function $\propto -2\pi|x-x'|$, which grows without bound and behaves quite unlike 3D electrostatics) would not give a sensible bounded potential. So the model borrows the *3D-like* $1/\sqrt{r^{2}+a^{2}}$ interaction throughout, for both $v_\mathrm{ext}$ and $v_H$, rather than literally inverting the 1D Laplacian. Treat the pseudocode's Poisson line as the general 3D recipe and this convolution as its toy-model stand-in.
- **LDA exchange.** We use the 3D LDA exchange potential $v_x \propto -n^{1/3}$ applied to our 1D density. This is *not* the rigorous 1D exchange (which is different functionally), but it is pedagogically standard and produces qualitatively correct behaviour.

!!! note "Where the factor of 2 comes from: the closed-shell assumption"
    Throughout this implementation — in `build_density` (`return 2.0 * ...`), in the band sum `band = 2.0 * sum(eigvals[:n_occ])`, and in the energy derivation above — a factor of $2$ multiplies every sum over the $N_\mathrm{occ}$ spatial orbitals. This is the **closed-shell, spin-unpolarised** assumption: each spatial orbital $\phi_i$ is occupied by *two* electrons, one spin-up and one spin-down, sharing the same spatial wavefunction. With $N$ electrons there are then $N_\mathrm{occ} = N/2$ doubly-occupied orbitals (the code sets `n_occ = n_electrons // 2`), and the density is
    $$
    n(\mathbf r) = \sum_{i,\sigma}^\mathrm{occ}|\phi_{i\sigma}(\mathbf r)|^{2}
    = 2\sum_{i}^{N/2}|\phi_i(\mathbf r)|^{2},
    $$
    since the two spin channels contribute identical $|\phi_i|^{2}$. The same doubling applies to the eigenvalue sum because each orbital energy $\varepsilon_i$ is counted once per occupying electron. This is *only* valid when the system is non-magnetic and has no partially filled shell; spin-polarised or open-shell systems require separate up- and down-spin densities and orbital sets (and the factor $2$ disappears), as flagged in the closing "Pedagogical, not production" note.
- **Mixing.** The first few iterations use linear mixing with $\alpha = 0.3$ to stabilise the history; subsequent iterations switch to Pulay/DIIS using the most recent six densities.
- **Total energy.** Equation (5.28) is computed via the band-sum form: $E = 2\sum_i^\mathrm{occ}\varepsilon_i - U_H + E_x - \int n v_x$, accounting for double-counting between the band-sum $2\sum\varepsilon_i$ and the explicit Hartree/exchange energies.

??? note "Full derivation: the double-counting correction $E = 2\sum_i\varepsilon_i - U_H + E_{xc} - \int n\,v_{xc}$"
    The Kohn–Sham total energy is *not* simply the sum of occupied eigenvalues. The eigenvalues already contain the electron–electron interaction once *inside* each orbital's potential, so naively summing them counts the Hartree and xc interactions twice — hence "double counting". We reconcile the band sum with the true energy functional explicitly. Throughout, the system is closed-shell with doubly-occupied orbitals, so factors of $2$ accompany every sum over the $N_\mathrm{occ}$ spatial orbitals (see the note on this assumption below).

    **The band sum, expanded.** Each occupied orbital satisfies $\big(-\tfrac12\nabla^2 + v_\mathrm{KS}\big)\phi_i = \varepsilon_i\phi_i$. Taking the expectation value in $\phi_i$ (normalised), $\varepsilon_i = \langle\phi_i| -\tfrac12\nabla^2 |\phi_i\rangle + \langle\phi_i| v_\mathrm{KS}|\phi_i\rangle$. Sum over occupied orbitals with the spin factor $2$:
    $$
    2\sum_i^\mathrm{occ}\varepsilon_i
    = \underbrace{2\sum_i^\mathrm{occ}\langle\phi_i| -\tfrac12\nabla^2|\phi_i\rangle}_{\displaystyle T_s}
    + 2\sum_i^\mathrm{occ}\langle\phi_i| v_\mathrm{KS}|\phi_i\rangle.
    $$
    The first group is, by definition, the non-interacting kinetic energy $T_s$. In the second group $v_\mathrm{KS}(\mathbf r)$ is a local multiplicative potential, so $2\sum_i\langle\phi_i|v_\mathrm{KS}|\phi_i\rangle = \int v_\mathrm{KS}(\mathbf r)\,\big[2\sum_i|\phi_i(\mathbf r)|^{2}\big]\,\mathrm d\mathbf r = \int n(\mathbf r)\,v_\mathrm{KS}(\mathbf r)\,\mathrm d\mathbf r$, using $n = 2\sum_i|\phi_i|^{2}$. Substituting $v_\mathrm{KS} = v_\mathrm{ext} + v_H + v_{xc}$,
    $$
    2\sum_i^\mathrm{occ}\varepsilon_i
    = T_s + \int n\,v_\mathrm{ext} + \int n\,v_H + \int n\,v_{xc}.
    \tag{5.46a}
    $$

    **The Hartree term carries a factor 2.** The Hartree potential is $v_H(\mathbf r) = \int \dfrac{n(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r'$, whereas the Hartree *energy* is the half-weighted double integral
    $$
    U_H = \frac{1}{2}\iint \frac{n(\mathbf r)\,n(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r'.
    $$
    The factor $\tfrac12$ corrects for counting each electron pair twice. Comparing, $\int n\,v_H = \iint \dfrac{n(\mathbf r)n(\mathbf r')}{|\mathbf r-\mathbf r'|}\,\mathrm d\mathbf r\,\mathrm d\mathbf r' = 2U_H$. So (5.46a) becomes
    $$
    2\sum_i^\mathrm{occ}\varepsilon_i
    = T_s + \int n\,v_\mathrm{ext} + 2U_H + \int n\,v_{xc}.
    \tag{5.46b}
    $$

    **The true total energy.** The Kohn–Sham energy functional itself is
    $$
    E = T_s + \int n\,v_\mathrm{ext} + U_H + E_{xc},
    \tag{5.46c}
    $$
    with the Hartree energy entering *once* (as $U_H$, not $2U_H$) and the xc *energy* $E_{xc}$ entering (not $\int n\,v_{xc}$, which is a different object — $v_{xc} = \delta E_{xc}/\delta n$).

    **Subtract.** Take (5.46c) minus (5.46b) term by term. The $T_s$ and $\int n\,v_\mathrm{ext}$ pieces are identical and cancel in the difference, leaving
    $$
    E - 2\sum_i\varepsilon_i = \big(U_H - 2U_H\big) + \big(E_{xc} - \textstyle\int n\,v_{xc}\big) = -U_H + E_{xc} - \int n\,v_{xc}.
    $$
    Rearranging gives the working formula used in `total_energy`,
    $$
    \boxed{\,E = 2\sum_i^\mathrm{occ}\varepsilon_i - U_H + E_{xc} - \int n\,v_{xc}\,}.
    \tag{5.46d}
    $$
    Watch the signs: the Hartree term flips from $+2U_H$ in the band sum to a net $-U_H$ (over-counted by $2U_H$, corrected by $+U_H$), and the xc term flips from $+\int n\,v_{xc}$ to a net $E_{xc} - \int n\,v_{xc}$. In the code (LDA exchange only, so $E_{xc}\to E_x$) this is exactly `e_dc = -eH + ex - sum(n*vx)*dx` with `eH` $=U_H$, `ex` $=E_x$, and `sum(n*vx)*dx` $=\int n\,v_x$, added to `band` $=2\sum_i\varepsilon_i$.

### What you should see

Running `python scf_1d_hchain.py` produces:

- Convergence in roughly 10–20 iterations to $|\Delta n|_\infty < 10^{-6}$.
- A density peaked at each proton site, with smooth bonding charge in between — a one-dimensional analogue of the bond-charge build-up in a real H chain.
- Eigenvalues that group into a "band" of four nearly-degenerate states for the four-atom chain — the discrete analogue of the bonding band of an infinite 1D chain.
- A figure `scf_1d_hchain.png` showing the converged density and Kohn–Sham potential.

!!! note "Pedagogical, not production"
    This code is a teaching tool, not a research code. The Hartree integral is $\mathcal O(N^{2})$ in grid size (no FFT), the LDA exchange is the 3D form applied to a 1D problem, and there is no provision for spin polarisation, correlation, or geometry optimisation. Chapter 6 introduces production codes (Quantum ESPRESSO, VASP, ABINIT) that handle all of this rigorously.

## 5.5.5a A bird's-eye view of the algorithm

To consolidate, here is the SCF loop in line-by-line pseudocode form, suitable for translation into any language:

```
input: positions {R_alpha}, charges {Z_alpha}, electron count N
input: tolerance tol, max iterations max_iter, mixing parameter alpha
input: history depth m, Pulay start iteration k_pulay
output: converged density n*, eigenvalues {eps_i}, total energy E

initialise n^(0) = sum of free-atom densities centred at {R_alpha}
compute v_ext = -sum_alpha Z_alpha / |r - R_alpha|
n_hist = [], r_hist = []

for k in 0, 1, ..., max_iter - 1:
    v_H = solve Poisson: nabla^2 v_H = -4 pi n^(k)
    v_xc = delta E_xc / delta n  [evaluated at n^(k)]
    v_KS = v_ext + v_H + v_xc

    H = T_kin + diag(v_KS)
    (eps_i, phi_i) = lowest N_occ eigenpairs of H

    n_out = sum_i f_i |phi_i|^2  [with f_i Fermi-Dirac if metal]
    n_out *= N / int n_out dr  [normalise]

    E = band sum + double-counting corrections [eq. 5.28]

    r = n_out - n^(k)
    if max|r| < tol: break

    n_hist.append(n^(k))
    r_hist.append(r)

    if k < k_pulay:
        n^(k+1) = (1 - alpha) n^(k) + alpha n_out
    else:
        solve eq. 5.43 for c_j given last m residuals
        n^(k+1) = sum_j c_j n_out^(j)
        # optionally apply Kerker preconditioning to r before mixing

return n^(k+1), eps, E
```

The structure is the same in every plane-wave or real-space DFT code. What varies between codes is (i) the basis representation (plane waves, atomic orbitals, real-space grid), (ii) the choice of mixing scheme and acceleration, (iii) the handling of metals and magnetic systems, and (iv) the parallelisation strategy. The conceptual core remains this six-step loop.

## 5.5.6 Practical tips for SCF convergence

When your SCF fails to converge, here are the levers to pull, in approximate order of how often they help.

1. **Smear the occupations.** For metallic systems, the Fermi level sits in a band of states and small changes in the potential cause discontinuous re-occupation between iterations. Replace the integer occupations $f_i \in \{0,2\}$ with smooth Fermi–Dirac (or Gaussian, or Methfessel–Paxton) occupations $f_i = f(\varepsilon_i, \mu, T)$ at an artificial electronic temperature $T$ (typically $kT = 0.05$–$0.2$ eV). This is *always* needed for metals and often helpful for small-gap systems.
2. **Reduce the mixing parameter $\alpha$.** Try $\alpha = 0.1$ or even $0.05$ if oscillating.
3. **Increase the history depth $m$ for DIIS/Anderson.** From 6 to 10 or more.
4. **Mix the potential, not the density** (or vice versa). Some codes mix $v_\mathrm{KS}$; others mix $n$; for hard problems one can flip the choice.
5. **Pre-conditioning.** Long-wavelength density oscillations (charge sloshing in metals) decouple from short-wavelength ones. The **Kerker preconditioner** multiplies the residual in reciprocal space by $|\mathbf q|^{2}/(|\mathbf q|^{2}+q_0^{2})$ to dampen the long-wavelength components that drive oscillation.
6. **Better initial guess.** Restart from the previous geometry's wavefunctions in a relaxation; superposition of atomic densities is good; uniform is the worst.
7. **Tighten the grid / basis.** Sometimes SCF instabilities are artefacts of an under-converged basis or **k**-point sampling.

## 5.5.6a A numerical convergence comparison

To make the convergence behaviour concrete, we compare three mixing schemes on the 1D hydrogen-chain test case (the code above):

| Scheme | Iterations to $|\Delta n|_\infty < 10^{-6}$ | Final energy (Ha) | Wallclock (s) |
|---|---|---|---|
| Linear $\alpha=0.1$ | 78 | $-2.341$ | 11 |
| Linear $\alpha=0.3$ | 28 | $-2.341$ | 4 |
| Linear $\alpha=0.5$ | oscillates | — | — |
| Pulay $m=6$ + linear start | 14 | $-2.341$ | 2 |
| Pulay $m=10$ + linear start | 12 | $-2.341$ | 2 |

The pattern is unsurprising: linear mixing with too small $\alpha$ is robust but slow; with too large $\alpha$ it oscillates; Pulay achieves the same final density in a fraction of the iterations. For production-quality plane-wave calculations on real materials (where each SCF step is $\sim 1\;\text{min}$ rather than $\sim 0.1\;\text{s}$), the speedup matters enormously.

!!! tip "Running it and watching the residual fall"
    Save the listing as `scf_1d_hchain.py` and run `python scf_1d_hchain.py`. Each line the loop prints is one trip round the SCF cycle, for example:

    ```
    iter   3   E = -2.340421    |dn|_inf = 3.50e-03
    ```

    Read it left to right: this is the 4th iteration (`iter 3`, counting from 0); the total energy is currently $-2.340421$ Ha; and `|dn|_inf` $= 3.5\times10^{-3}$ is the residual $\max|n_\mathrm{out}-n^{(k)}|$ — the largest the input and output densities disagree at any grid point. The number to watch is that last column: it should shrink towards zero. When it drops below `tol` $=10^{-6}$, the input and output densities agree everywhere to six decimal places, the loop prints `Converged`, and stops. The energy column settling (e.g. "stopped changing at the sixth decimal", $|\Delta E|<10^{-6}$) is a *convenient secondary signal* that you have reached the fixed point — but it is the density residual that is being tested here. Try changing `alpha=0.3` to `alpha=0.05` and watch the residual fall more slowly but more smoothly; that is linear mixing trading speed for stability.

!!! example "Energy convergence trace"
    Running the SCF code with Pulay mixing and printing $|E^{(k+1)} - E^{(k)}|$ at each iteration produces a sequence like:
    
    ```
    iter   0   E = -2.105674    |dE| = -
    iter   1   E = -2.298471    |dE| = 1.9e-1
    iter   2   E = -2.336891    |dE| = 3.8e-2
    iter   3   E = -2.340421    |dE| = 3.5e-3   <- switch to Pulay
    iter   4   E = -2.340980    |dE| = 5.6e-4
    iter   5   E = -2.341022    |dE| = 4.2e-5
    iter   6   E = -2.341024    |dE| = 2.0e-6
    iter   7   E = -2.341024    |dE| = 8.5e-8
    iter   8   E = -2.341024    |dE| = 3.1e-9
    iter   9   E = -2.341024    |dE| = 9.4e-11  <- converged
    ```
    
    Each Pulay step gains roughly an order of magnitude in $|\Delta E|$ — the characteristic *superlinear* convergence of an acceleration scheme. Once the iterate is close enough to the fixed point for the linearised analysis to apply, convergence is essentially geometric with a very small ratio.

## 5.5.6b Diagnosing a stuck SCF

A field guide to common SCF pathologies and their symptoms:

- *Energy oscillating between two values, density also oscillating.* Charge sloshing in a metallic-like system. Cure: reduce $\alpha$, add Kerker preconditioning, add Fermi smearing.
- *Energy slowly increasing instead of decreasing.* Almost always means the initial guess is in a basin of a different local minimum (magnetic ground state, charge transfer state). Cure: reset, try a different initial density or different magnetic moments.
- *Energy decreasing but density residual stuck.* Numerical noise from under-converged eigensolver or too coarse a grid. Cure: tighten basis or eigensolver tolerance.
- *DIIS divergence after looking converged.* The Pulay history has become near-singular (linearly dependent residuals). Cure: clear history and restart with linear mixing, or use trust-region DIIS.
- *Energy NaN at iteration 1.* Almost always the Hartree integral singular due to a vanishing density at some grid point. Cure: floor the density at $n_\min\sim 10^{-10}$ before evaluating $v_{xc}$.

!!! question "Check yourself"
    1. Why is the Kohn–Sham problem *nonlinear*, even though each individual diagonalisation is a perfectly linear eigenvalue problem?
    2. What does *mixing* actually do, and does it change the final converged density?
    3. Why can naive iteration ($n^{(k+1)} = n_\mathrm{out}^{(k)}$) oscillate or diverge instead of settling down?
    4. If an SCF run converges to a tight tolerance, does that guarantee the answer is physically correct?

    ??? success "Answers"
        1. The diagonalisation is linear *for a fixed potential*, but the potential $v_\mathrm{KS}[n]$ itself depends on the density $n$ that the diagonalisation produces. The map "density $\to$ potential $\to$ orbitals $\to$ new density" is therefore a nonlinear self-referential relation $n = \mathcal F[n]$, which is why it must be solved iteratively rather than in one shot.
        2. Mixing blends the new (output) density with the previous (input) density — e.g. linear mixing keeps only a fraction $\alpha$ of the change — to damp out oscillations and reach the fixed point reliably. It is a *numerical aid only*: it changes the *path* taken to convergence and the *speed*, but not the destination. Any choice of mixing that converges lands on the same self-consistent density.
        3. Because the SCF map can amplify errors: the screening response of the electrons can overshoot, so an excess of charge in one region provokes an over-correction that creates a deficit, then an excess again ("charge sloshing"). Formally, the SCF Jacobian can have eigenvalues with magnitude greater than 1 (especially for metals), so the error grows each step instead of shrinking. Mixing with small enough $\alpha$ brings the effective eigenvalues back inside the unit circle.
        4. No. Convergence means the density is self-consistent *with the chosen functional, basis, pseudopotential and* **k***-points* — not that those choices describe reality. A tightly converged calculation built on an inadequate functional can give a confidently wrong answer. Convergence is necessary but not sufficient; see §5.6.

## 5.5.7 Closing the loop

We now have a complete picture: choose a functional (§5.4), discretise the wavefunctions on a grid (or basis), iterate the SCF loop with sensible mixing, and read out the converged density, total energy, eigenvalues, and (with the Hellmann–Feynman theorem) forces. From these, all the standard observables follow: cohesive energies, lattice constants, elastic moduli, vibrational spectra, band structures.

The next section, §5.6, gives an honest account of where this machine fails and what to do about it. Chapter 6 turns the cogs of a production calculation: plane-wave basis sets, pseudopotentials, **k**-point sampling, and practical convergence testing.

### Summary of §5.5 — what to remember in 3 months

- **The SCF loop**: guess $n^{(0)}$; build $v_\mathrm{KS}$; diagonalise; build new density; mix; repeat until $\|n_\mathrm{out}-n_\mathrm{in}\|<\text{tol}$.
- **Why naive iteration fails**: the SCF Jacobian $\mathcal J = \delta n_\mathrm{out}/\delta n_\mathrm{in}$ can have eigenvalues $|\lambda|>1$, especially for metals (charge sloshing).
- **Linear mixing**: $n^{(k+1)} = (1-\alpha)n^{(k)} + \alpha n_\mathrm{out}^{(k)}$. Robust but slow. Typical $\alpha\sim 0.1$–$0.3$.
- **Pulay/DIIS**: minimise $\|\sum c_j r^{(j)}\|^{2}$ over recent residuals, $\sum c_j = 1$. Converges quadratically near the fixed point.
- **Anderson, Broyden**: closely related acceleration schemes.
- **Kerker preconditioner**: $\tilde r(\mathbf q) = q^{2}/(q^{2}+q_0^{2})\cdot r(\mathbf q)$ suppresses long-wavelength sloshing for metals.
- **Fermi smearing**: replace step occupations with smooth Fermi–Dirac at artificial $T$. Essential for metals (Janak's theorem and Mermin's finite-$T$ DFT justify this).
- **Convergence criteria**: density residual, energy difference, force changes — tighter for phonons and elastic constants than for total energies.
- **Convergence ≠ correctness**: a converged SCF solves the chosen functional self-consistently, not necessarily reality.

!!! note "What Janak's theorem and Mermin's finite-temperature DFT each say"
    Two results underpin the Fermi-smearing trick. *Janak's theorem* states that the Kohn–Sham eigenvalue is the derivative of the total energy with respect to the occupation of that orbital, $\varepsilon_i = \partial E/\partial f_i$, which is what makes *fractional* occupations a well-defined variational quantity rather than an ad hoc device. *Mermin's finite-temperature DFT* generalises the Hohenberg–Kohn theorems to thermal ensembles, showing that the grand potential is a unique functional of the density at temperature $T$; minimising the free energy $F = E - TS$ with smooth Fermi–Dirac occupations is then a rigorous variational principle (the artificial $T$ enters as a real smearing temperature, and zero-temperature observables are recovered by extrapolating $T\to 0$).

!!! warning "Convergence ≠ correctness"
    A converged SCF only means that the equations have been *solved consistently* for the chosen functional, basis, and pseudopotential. It does *not* mean that the answer is physically correct. A PBE calculation can converge beautifully to an unphysical metallic state for NiO, or to a fictitiously short bond for a vdW-bound dimer. Convergence is necessary, not sufficient. The next section addresses the gap between convergence and physical truth.
