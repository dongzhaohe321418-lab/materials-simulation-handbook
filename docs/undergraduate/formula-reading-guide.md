# How to Read the Formulas

This page does not teach new physics. Its only job is to help you **read** the equations that already appear elsewhere in this handbook. When an equation stops a chapter dead for you, come back here, then return to the chapter.

!!! tip "The core idea"
    An equation is **compressed language**. Someone wrote a sentence — sometimes a whole
    paragraph — and squeezed it into symbols to save space. Reading an equation means doing
    the reverse: unpacking it *slowly*, one symbol at a time, back into ordinary words. You
    are not expected to "see" it all at once. Nobody does. Going slowly is the correct speed.

This is a slower, friendlier companion to the terse reference material; it overlaps on purpose so you have somewhere gentle to land.

## First: what is the equation *for*?

Before reading any symbol, ask one question: **what does this equation compute?** Most equations have a left-hand side and a right-hand side, and usually one side is the thing you *want* and the other side is the *recipe* for getting it.

!!! example "Inputs and outputs"
    In $E = T + V$ (total energy = kinetic + potential), the **output** is $E$ on the left;
    the **inputs** are $T$ and $V$ on the right. Reading direction: "to get $E$, add $T$ and
    $V$." Knowing which symbols are inputs and which is the output is half the battle.

A quick test: ask "if I changed this symbol, would the answer change?" If yes, it is an input. The symbol alone on one side is usually the output.

## Reading subscripts and superscripts

Small marks above and below a symbol are not decoration. They carry meaning, and the same-looking mark can mean different things.

- **Subscript as a label or index.** $n_i$ means "the $n$ belonging to item $i$" — for example the occupation of orbital number $i$. The $i$ picks one out of many.
- **Subscript on a vector.** $\mathbf{r}_i$ is the position vector of atom $i$. Bold says "vector"; the subscript says "which atom".
- **Subscript as a named level.** $\varepsilon_n$ is the energy of level $n$ (the $n$-th eigenvalue). Here $n$ counts states.
- **Superscript as a power.** In $n^{4/3}$ the superscript is an exponent: raise $n$ to the power $4/3$.
- **Superscript as an operation, not a power.** $\psi^*$ means the **complex conjugate** of $\psi$, not "$\psi$ squared". Context tells you which: a number on top is usually a power; a symbol like $*$ or $\dagger$ is an operation.

!!! warning
    $\psi^2$ and $\psi^*$ look similar but are unrelated. When you meet a superscript, pause
    and decide: power, or operation?

## Six kinds of object — and why it matters

Symbols stand for different *kinds* of thing. Telling them apart prevents most confusion. Here is each kind with one example drawn from this book.

- **Scalar** — a single number (with units). Example: the total energy $E$. One value.
- **Vector** — a list of numbers with direction, written in bold. Example: a position $\mathbf{r} = (x, y, z)$. Three numbers describing one point in space.
- **Matrix** — a grid of numbers, indexed by two subscripts. Example: the Hamiltonian matrix $H_{ij}$, where $i$ labels the row and $j$ the column.
- **Function** — a rule that takes an input and returns a number. Example: the electron density $n(\mathbf{r})$ takes a position $\mathbf{r}$ and returns "how much electron is here". Give it a point, get a number.
- **Operator** — a rule that takes a *function* and returns another function. Example: the Hamiltonian operator $\hat{H}$ (the hat is the giveaway) acts on a wavefunction and hands you back a new function.
- **Functional** — a rule that takes a *whole function* and returns a single number. Example: the energy functional $E[n]$ takes the entire density function $n(\mathbf{r})$ and returns one energy value. Note the **square brackets** $[\,\cdot\,]$ — a strong hint that you are looking at a functional.

!!! note "Function vs functional — read this twice"
    A **function** like $n(\mathbf{r})$ eats *one input* (a point) and gives *one number*.
    A **functional** like $E[n]$ eats *an entire function* (the whole density, everywhere at
    once) and gives *one number*. The clue is the brackets: round $(\,)$ for a function,
    square $[\,]$ for a functional. This distinction is the whole reason "density
    **functional** theory" is named the way it is: the energy is a functional *of* the density
    function. Get this and the name stops being mysterious.

## Tracking units

Every physical quantity carries units, and units obey the same algebra as the symbols. This gives you a free error-checker.

The rule: **both sides of an equation must have the same units**, and you can only add quantities with the same units. You cannot add an energy to a length.

!!! example "Units as a safety net"
    Suppose you derive an "energy" and the units come out as Å (a length). You have made a
    mistake *somewhere* — no calculator needed to know that. Conversely, if both sides come
    out in eV, that is a good sign (though not a proof) that the algebra is right. In this book
    energies live in eV or hartree ($1\,\mathrm{Ha} = 27.2114$ eV) and lengths in ångström or
    bohr ($1\,a_0 = 0.529177$ Å); keeping track of which is which prevents silent errors.

## Reading sums and integrals as "add up contributions"

A summation and an integral say the same thing in plain words: **add up many small contributions**.

- $\displaystyle\sum_i x_i$ reads as "add up $x_i$ for every $i$" — first $x_1$, then $x_2$, and so on, totalling them. The $i$ under the $\sum$ is the counter.
- $\displaystyle\int f(\mathbf{r})\,\mathrm{d}\mathbf{r}$ reads as "add up $f$ over all of space". An integral is a sum where the pieces are infinitely small. Think of chopping space into tiny boxes, evaluating $f$ in each box, and totalling — the $\mathrm{d}\mathbf{r}$ is the size of one tiny box.

If you can read $\sum$ and $\int$ both as "total it up", a large fraction of intimidating equations becomes approachable.

## Reading derivatives and gradients as rates of change

Derivatives measure **how fast something changes**.

- $\dfrac{\mathrm{d}f}{\mathrm{d}x}$ is the **slope**: how much $f$ changes when $x$ changes a little. Positive means rising, negative means falling, zero means flat (a peak, valley, or plateau).
- $\nabla f$ (the **gradient**) is the slope generalised to several directions at once. It is a vector pointing "uphill" — the direction in which $f$ increases fastest. A force is often $-\nabla V$: push downhill in energy.
- $\nabla^2 f$ (the **Laplacian**) measures **curvature** — how much $f$ bends. It appears in the kinetic-energy term of quantum mechanics.

A simple mantra: first derivative = slope, second derivative = curvature.

## Reading an eigenvalue equation

The pattern $\hat{A}\psi = a\psi$ appears constantly. In words:

> "When the operator $\hat{A}$ acts on the special function $\psi$, the result is just $\psi$ again, **scaled** by the number $a$."

The function $\psi$ is called an **eigenfunction** and the number $a$ its **eigenvalue**. The key idea is that most functions get *reshaped* by an operator, but eigenfunctions are special: they come back unchanged in *shape*, only stretched by a factor $a$.

## Reading probability expressions

In quantum mechanics, $|\psi(\mathbf{r})|^2$ reads as a **probability density**: "how likely you are to find the particle near point $\mathbf{r}$". It is a density, so to get an actual probability you integrate it over a region — adding up contributions, exactly as above. The conjugate appears because $|\psi|^2 = \psi^*\psi$, which is always real and non-negative, as a probability must be.

## Reading a machine-learning loss

A **loss function** measures **how wrong a model is**, as a single number, so smaller is better. The mean-squared error,

$$
L = \frac{1}{N}\sum_{i=1}^{N}\bigl(y_i - \hat{y}_i\bigr)^2,
$$

reads as: "for each of the $N$ examples, take the gap between the true value $y_i$ and the prediction $\hat{y}_i$, square it (so over- and under-shooting both count as a penalty), then **average** the penalties." Squaring punishes large mistakes more than small ones. Training a model means adjusting it to make this average penalty as small as possible.

---

## Four worked "read-this-equation" examples

For each, read the equation left to right, naming every symbol, then say what the whole thing means.

### 1. A simple energy expression

$$
E = T + V
$$

1. $E$ — a **scalar**, the total energy (output), in eV or hartree.
2. $=$ — "is equal to / is computed as".
3. $T$ — the kinetic energy (energy of motion), an input.
4. $+$ — ordinary addition; allowed because $T$ and $V$ share units.
5. $V$ — the potential energy (energy of position/interaction), an input.

**In words:** "The total energy is the kinetic energy plus the potential energy." A closely related density-functional form is the LDA exchange energy, $E_x \propto \int n(\mathbf{r})^{4/3}\,\mathrm{d}\mathbf{r}$, which reads "add up $n^{4/3}$ over all space" — a single number (a functional of the density), built by integrating a function of $n$ at every point.

### 2. A Schrödinger-style eigenvalue equation

$$
\hat{H}\psi = E\psi
$$

1. $\hat{H}$ — an **operator** (the hat says so); the Hamiltonian, encoding total energy.
2. $\psi$ — a **function**, the wavefunction; the eigenfunction we are looking for.
3. $\hat{H}\psi$ — "let $\hat{H}$ act on $\psi$", producing a new function.
4. $=$ — the two sides are the same function.
5. $E$ — a **scalar** eigenvalue, the energy of this state.
6. $E\psi$ — the original $\psi$ simply scaled by the number $E$.

**In words:** "Acting on $\psi$ with the energy operator gives back the same $\psi$, scaled by the energy $E$." The special functions $\psi$ that satisfy this are the allowed states, and the matching $E$ values are their allowed energies.

### 3. A self-consistent (fixed-point) idea

$$
n_{\text{out}} = F[n_{\text{in}}], \qquad \text{converged when } n_{\text{out}} \approx n_{\text{in}}
$$

1. $n_{\text{in}}$ — a **function**, the density you guess going in.
2. $F[\,\cdot\,]$ — a **functional / procedure** (square brackets): "do the whole calculation on this density". Here it stands for one full pass of the machinery.
3. $n_{\text{out}} = F[n_{\text{in}}]$ — feed the guess in, run the procedure, get an updated density out.
4. $n_{\text{out}} \approx n_{\text{in}}$ — the stopping rule: the density you get out matches the density you put in.

**In words:** "Put a density in, the procedure gives an updated density; repeat, feeding each output back as the next input, until output and input agree." When they agree, the density is **self-consistent** — it reproduces itself. This is exactly the self-consistent field (SCF) loop used in Kohn–Sham DFT (Section 5.5): guess a density, build the potential, solve, get a new density, and iterate to a fixed point. See the chapter overview, [Chapter 5 (DFT)](../ch05-dft/index.md), for the full picture.

### 4. A mean-squared-error ML loss

$$
L = \frac{1}{N}\sum_{i=1}^{N}\bigl(y_i - \hat{y}_i\bigr)^2
$$

1. $L$ — a **scalar**, the loss (output); smaller is better.
2. $\dfrac{1}{N}$ — divide by the number of examples $N$: this makes it an **average**.
3. $\displaystyle\sum_{i=1}^{N}$ — "add up over every example $i$ from $1$ to $N$".
4. $y_i$ — the true value for example $i$ (a label/index subscript).
5. $\hat{y}_i$ — the model's prediction for example $i$ (the hat here means "estimate").
6. $(y_i - \hat{y}_i)$ — the error: how far the prediction is from the truth.
7. $(\,\cdots)^2$ — square it, so the penalty is always positive and big errors count extra.

**In words:** "For each example, square the gap between truth and prediction, then average those squared gaps." It is the **average penalty** for being wrong; training shrinks it.

---

!!! question "When you meet a new equation, ask…"
    Run through this checklist, slowly, every time:

    1. **What does it compute?** Which side is the output, which side is the recipe?
    2. **What kind of object is each symbol** — scalar, vector, matrix, function, operator, or functional? (Watch the hats and the bracket shapes.)
    3. **What do the sub/superscripts mean** — an index, a label, a power, or an operation like $*$?
    4. **Do the units match** on both sides? If not, something is wrong.
    5. **Is there a $\sum$ or $\int$?** If so, read it as "add up contributions".
    6. **Is there a derivative or $\nabla$?** Read it as a slope, a rate of change, or a curvature.
    7. **Is it an eigenvalue equation** $\hat{A}\psi = a\psi$? Then $\psi$ comes back scaled by $a$.
    8. **Can I say the whole thing in one plain English sentence?** If yes, you have read it.
