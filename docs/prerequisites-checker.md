# Prerequisites Check

This is a twenty-question self-test. The point is not to grade you, but to help you decide where to start reading. The questions are deliberately mixed — some are easy if you have just finished an undergraduate physics module; some require a moment of thought. Each one ends with a folded "Answer" block that gives you the solution and a remediation pointer to a specific section of the book.

Spend no more than five minutes per question; if you genuinely do not know, mark it wrong and move on. The scoring rubric at the bottom tells you which chapter to start at.

---

## Calculus (Questions 1–3)

### Question 1

Compute the derivative

$$ \frac{d}{dx} \left[ x^2 \ln(x) \right] $$

at $x = 1$.

- (a) 0
- (b) 1
- (c) 2
- (d) $e$

??? note "Answer"
    **(b) 1.** By the product rule, $\frac{d}{dx}[x^2 \ln x] = 2x\ln x + x$. At $x=1$, $\ln 1 = 0$, so the result is $1$.

    If you got this wrong, read [Ch 0.3 §2 (the product and chain rules)](ch00-math/03-calculus.md).

### Question 2

Evaluate

$$ \int_0^{\infty} e^{-x^2}\,dx. $$

- (a) $1$
- (b) $\pi$
- (c) $\sqrt{\pi}/2$
- (d) $\sqrt{\pi}$

??? note "Answer"
    **(c) $\sqrt{\pi}/2$.** This is half of the full Gaussian integral $\int_{-\infty}^{\infty} e^{-x^2}\,dx = \sqrt{\pi}$.

    Gaussian integrals appear constantly in statistical mechanics, in basis-set expansions, and inside neural-network activations. If you have not seen this trick, read [Ch 0.3 §4 (Gaussian integrals and gradients)](ch00-math/03-calculus.md).

### Question 3

Short answer. State, in words, what the gradient $\nabla f(\mathbf{x})$ of a scalar function $f : \mathbb{R}^n \to \mathbb{R}$ represents geometrically.

??? note "Answer"
    The vector of partial derivatives, pointing in the direction of steepest ascent of $f$ at the point $\mathbf{x}$; its magnitude is the rate of increase per unit step in that direction. The negative gradient is the direction of steepest descent and is what every optimiser, every force, and every SCF loop is secretly chasing.

    If this felt vague, read [Ch 0.3 §3 (gradients and the variational principle)](ch00-math/03-calculus.md).

---

## Linear algebra (Questions 4–6)

### Question 4

Let $A = \begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$. What are its eigenvalues?

- (a) 1 and 2
- (b) 1 and 3
- (c) 0 and 3
- (d) 2 and 2

??? note "Answer"
    **(b) 1 and 3.** The characteristic equation is $(2-\lambda)^2 - 1 = 0$, giving $\lambda = 1$ or $\lambda = 3$.

    Every electronic-structure calculation is an eigenvalue problem in disguise. If this took you more than a minute, read [Ch 0.2 §4 (eigenvalues and diagonalisation)](ch00-math/02-linear-algebra.md).

### Question 5

A matrix $U \in \mathbb{C}^{n \times n}$ is **unitary** if:

- (a) $U^T = U$
- (b) $U^\dagger U = I$
- (c) $\det(U) = 1$
- (d) all entries of $U$ are real

??? note "Answer"
    **(b) $U^\dagger U = I$.** A unitary matrix has orthonormal columns; its conjugate transpose is its inverse.

    Quantum time evolution, Bloch transformations and irreducible representations are all built on unitary matrices. If you are unsure why this matters, read [Ch 0.2 §6 (inner products, orthogonality, unitarity)](ch00-math/02-linear-algebra.md).

### Question 6

For matrices $A \in \mathbb{R}^{m \times n}$ and $B \in \mathbb{R}^{n \times p}$, which of the following is generally **false**?

- (a) $(AB)^T = B^T A^T$
- (b) $AB$ has shape $m \times p$
- (c) $AB = BA$
- (d) $A(B+C) = AB + AC$ for compatible $C$

??? note "Answer"
    **(c)** is false in general. Matrix multiplication is **not** commutative. (a), (b) and (d) all hold.

    If you ticked (a) or (b), you have a bigger gap than this checker can fix; read [Ch 0.2 §§2–3 (vectors, matrices, multiplication rules)](ch00-math/02-linear-algebra.md) before going further.

---

## Basic quantum mechanics (Questions 7–10)

### Question 7

A particle is in a one-dimensional infinite square well of length $L$. The ground-state energy is:

- (a) $\frac{\hbar^2 \pi^2}{2 m L^2}$
- (b) $\frac{\hbar^2}{2 m L^2}$
- (c) $\frac{\hbar^2 \pi^2}{m L^2}$
- (d) $\frac{\hbar \pi}{2 m L}$

??? note "Answer"
    **(a) $\hbar^2 \pi^2 / (2 m L^2)$.** The allowed energies are $E_n = n^2 \hbar^2 \pi^2 / (2 m L^2)$ with $n = 1, 2, 3, \dots$.

    If you do not recognise this immediately, read [Ch 4.3 (particle in a box, both analytic and numerical)](ch04-quantum/03-particle-in-box.md).

### Question 8

Which of the following is **not** a postulate of quantum mechanics?

- (a) The state of a system is described by a vector in a Hilbert space.
- (b) Observables correspond to Hermitian operators.
- (c) Measurement collapses the state to an eigenstate of the measured observable.
- (d) The wavefunction always factorises into a product of single-particle wavefunctions.

??? note "Answer"
    **(d)** is false. Generic many-particle wavefunctions are entangled and do **not** factorise. This is precisely what makes the many-electron problem hard and motivates DFT.

    Read [Ch 4.5 (many electrons and the exponential wall)](ch04-quantum/05-many-electrons.md) if this surprised you.

### Question 9

Short answer. State the time-independent Schrödinger equation in one line.

??? note "Answer"
    $\hat{H}\psi = E\psi$, where $\hat{H} = -\frac{\hbar^2}{2m}\nabla^2 + V(\mathbf{r})$ for a single particle in a potential. Anything equivalent (Dirac notation $\hat{H}|\psi\rangle = E|\psi\rangle$) is fine.

    If you could not produce this, read [Ch 4.2 (the Schrödinger equation, derived twice)](ch04-quantum/02-schrodinger.md).

### Question 10

The Born–Oppenheimer approximation says, roughly, that:

- (a) Electrons and nuclei have the same kinetic energy.
- (b) Electrons follow the nuclei adiabatically because they are much lighter and faster.
- (c) Nuclei move on classical trajectories while electrons stay frozen.
- (d) The wavefunction is separable in momentum space.

??? note "Answer"
    **(b).** Because $m_e \ll m_N$, electrons respond essentially instantaneously to nuclear motion, so we can solve the electronic problem at fixed nuclear positions and treat the nuclear motion on the resulting potential-energy surface.

    Read [Ch 4.6 (Born–Oppenheimer separation, with a worked example)](ch04-quantum/06-born-oppenheimer.md).

---

## Basic statistical mechanics (Questions 11–13)

### Question 11

In the canonical ensemble, the probability of a microstate $i$ with energy $E_i$ is:

- (a) $e^{-E_i}$
- (b) $e^{-\beta E_i} / Z$, where $\beta = 1/(k_B T)$
- (c) constant for all $i$
- (d) $E_i / Z$

??? note "Answer"
    **(b).** The Boltzmann distribution, with partition function $Z = \sum_i e^{-\beta E_i}$.

    Read [Ch 8.1 (ensembles and partition functions)](ch08-statmech/01-ensembles.md) if you needed to look this up.

### Question 12

The equipartition theorem says that, classically, each quadratic degree of freedom contributes:

- (a) $k_B T$ to the average energy
- (b) $\frac{1}{2} k_B T$ to the average energy
- (c) $\frac{3}{2} k_B T$ to the average energy
- (d) Nothing — only quantum oscillators contribute.

??? note "Answer"
    **(b) $\frac{1}{2} k_B T$.** Hence a classical monatomic ideal gas has $\langle E \rangle = \frac{3}{2} k_B T$ per atom.

    Read [Ch 8.1 §3 (equipartition and its quantum failure mode)](ch08-statmech/01-ensembles.md).

### Question 13

Short answer. What is the second law of thermodynamics, in one sentence?

??? note "Answer"
    The entropy of an isolated system never decreases; or equivalently, no process can convert heat entirely into work without other change. Any reasonable paraphrase is fine.

    If you are rusty, the relevant background is in [Ch 8.1](ch08-statmech/01-ensembles.md) and the appendix [thermodynamics refresher](appendix/A-math-reference.md).

---

## Basic Python and NumPy (Questions 14–16)

### Question 14

What does this snippet print?

```python
import numpy as np
a = np.arange(6).reshape(2, 3)
print(a.sum(axis=0))
```

- (a) `[3 5 7]`
- (b) `[3 12]`
- (c) `[0 1 2 3 4 5]`
- (d) `15`

??? note "Answer"
    **(a) `[3 5 7]`.** `axis=0` sums down columns, giving `[0+3, 1+4, 2+5]`.

    If this confused you, read [Ch 1.2 §3 (NumPy axes and broadcasting)](ch01-python/02-numpy.md).

### Question 15

Which of the following creates a $100 \times 100$ identity matrix?

- (a) `np.ones((100, 100))`
- (b) `np.eye(100)`
- (c) `np.identity(100, 100)`
- (d) `np.diag(100)`

??? note "Answer"
    **(b) `np.eye(100)`.** `np.identity(100)` also works (single argument). `np.diag(100)` is a $1\times 1$ matrix containing $100$.

    Read [Ch 1.2 §2 (NumPy array construction)](ch01-python/02-numpy.md).

### Question 16

Short answer. What is the difference between `a = b` and `a = b.copy()` for a NumPy array `b`?

??? note "Answer"
    `a = b` makes `a` and `b` refer to the **same** array object; mutating `a` mutates `b`. `a = b.copy()` creates an independent array with the same contents. This bites every materials simulation pipeline that updates configurations in-place during an MD loop.

    Read [Ch 1.2 §5 (views, copies and the in-place gotcha)](ch01-python/02-numpy.md) if this caught you out.

---

## Basic crystallography (Questions 17–18)

### Question 17

How many atoms are in the conventional unit cell of an FCC (face-centred cubic) crystal?

- (a) 1
- (b) 2
- (c) 4
- (d) 8

??? note "Answer"
    **(c) 4.** Eight corner atoms each shared by eight cells ($8 \times \frac{1}{8} = 1$) plus six face atoms each shared by two cells ($6 \times \frac{1}{2} = 3$), totalling $4$.

    Read [Ch 3.3 (crystals and lattices, with worked examples)](ch03-atoms/03-crystals.md).

### Question 18

The reciprocal lattice vector $\mathbf{b}_1$ corresponding to direct-lattice vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$ is defined by:

- (a) $\mathbf{b}_1 = \mathbf{a}_1$
- (b) $\mathbf{b}_1 = 2\pi \frac{\mathbf{a}_2 \times \mathbf{a}_3}{\mathbf{a}_1 \cdot (\mathbf{a}_2 \times \mathbf{a}_3)}$
- (c) $\mathbf{b}_1 = \mathbf{a}_2 + \mathbf{a}_3$
- (d) $\mathbf{b}_1 = 2\pi (\mathbf{a}_2 - \mathbf{a}_3)$

??? note "Answer"
    **(b).** The reciprocal lattice vectors are defined so that $\mathbf{a}_i \cdot \mathbf{b}_j = 2\pi \delta_{ij}$.

    If this looked unfamiliar, read [Ch 3.4 (reciprocal space without tears)](ch03-atoms/04-reciprocal.md) and then [Ch 3.5 §1 (Bloch's theorem)](ch03b-solid-state/01-bloch-theorem.md).

---

## Basic machine learning (Questions 19–20)

### Question 19

In supervised learning with mean-squared-error loss, the gradient of the loss with respect to a model parameter $\theta$ tells us:

- (a) The value of $\theta$ at the optimum.
- (b) The direction in which to change $\theta$ to **decrease** the loss.
- (c) The variance of the model's predictions.
- (d) Nothing useful without a learning rate.

??? note "Answer"
    **(b).** Specifically, the **negative** gradient points in the direction of steepest decrease, which is why gradient descent steps in $-\nabla_\theta L$.

    If you got this wrong, read [Ch 1 appendix on optimisation](ch01-python/04-reproducibility.md) and Ch 9 §1.

### Question 20

A **graph neural network** acts on data structured as:

- (a) A regular grid of pixels.
- (b) A sequence of tokens.
- (c) A set of nodes connected by edges, with features on either.
- (d) A single fixed-length vector.

??? note "Answer"
    **(c).** Crystals are natural graphs: nodes are atoms, edges encode neighbour distances. This is the structural reason GNNs work well for materials.

    Read [Ch 10.1 (crystals as graphs)](ch10-gnn/01-graphs.md) for a gentle introduction.

---

## Scoring rubric

Add up your correct answers. Short-answer questions count as correct if you produced the essential content, even if your wording differs from the model answer.

| Score | Recommended starting point |
|---|---|
| 0–8 | Start at **Chapter 0 (Mathematics)** and follow [Path A — Linear](learning-path.md#path-a-linear-1216-weeks-full-coverage). Do not skip the exercises; the prerequisites you are missing are exactly the kind of thing the rest of the book assumes silently. |
| 9–13 | **Skim Chapters 0–3, then start carefully at Chapter 4 (Quantum Mechanics).** Use the cross-reference matrix on the [Learning Path](learning-path.md#cross-reference-matrix) to identify which earlier sections are actually load-bearing for you. This is the most common starting profile. |
| 14–17 | **Light skim of Chapters 0–3, start at Chapter 5 (Density Functional Theory).** You have the maths and the elementary quantum mechanics; the book's value for you starts where it begins teaching electronic-structure machinery. Follow [Path B — Deep core](learning-path.md#path-b-deep-core-1012-weeks-skim-prerequisites). |
| 18–20 | **Jump straight to Chapter 9 (MLIPs)** or pick a Tier-2 project from `projects/` and follow [Path C — Project-driven](learning-path.md#path-c-project-driven-8-weeks). You can treat Chapters 0–8 as reference. If your weak spot was crystallography, Ch 3.5 is the one chapter to dip into before doing periodic-system ML. |

If your score sits awkwardly on a boundary, look at **which categories** you got wrong:

- Wrong in calculus/linear algebra → start lower than your score suggests.
- Wrong in QM/stat mech but fine in maths → start at Ch 4 regardless of score.
- Wrong in Python only → run through Ch 1 over a weekend and start where the rest of your score suggests.
- Wrong in crystallography only → read Ch 3 and Ch 3.5 first, then continue.

Re-take this test after every Tier completes. The remediation pointers will tell you which sections you still owe the book.
