# 4.5 Many electrons and the exponential wall

So far we have considered a single quantum particle in a one-dimensional potential. Diagonalising a 400×400 tridiagonal matrix on a laptop returns the lowest several eigenstates of the particle-in-a-box or the harmonic oscillator to seven significant figures in a fraction of a second. It is tempting to assume that the same approach, scaled up, will work for a real material — a few hundred electrons in three dimensions on a finer grid. This section explains why that assumption is wrong by an absurd margin.

The "many-body problem" is the central computational difficulty of quantum chemistry and condensed-matter physics. Two ingredients together produce it. First, the wavefunction of $N$ particles is a function not of three coordinates but of $3N$; the configuration space scales linearly with $N$, but the number of grid points needed to *sample* it scales exponentially. Second, the Pauli principle forces the wavefunction to be antisymmetric under particle exchange, which both rules out the simplest product ansatz and dictates a particular algebraic structure (Slater determinants). The interplay is what makes interacting fermion systems hard.

!!! info "What problem are we solving?"
    We can solve the Schrödinger equation for *one* electron on a laptop in
    a fraction of a second (Section 4.3). A real material has hundreds or
    thousands of electrons that all repel one another. The obvious plan —
    "use the same grid method, just with more coordinates" — fails not by a
    little but by an astronomical margin. This section's job is to *measure*
    how badly it fails (the "exponential wall"), and to introduce the two
    structural facts every later method must respect: that the wavefunction
    lives in a $3N$-dimensional space, and that for electrons it must be
    **antisymmetric**. Understanding *why* brute force is impossible is what
    motivates everything that follows — Hartree–Fock (Section 4.7), and
    above all density functional theory in [Chapter 5 (DFT)](../ch05-dft/index.md).

!!! note "Plain-language version"
    A single electron needs three numbers to say where it is: $x$, $y$, $z$.
    Two electrons need six; one hundred electrons need three hundred. The
    wavefunction is one function of *all* of those numbers at once — you
    cannot split it into "electron 1's wavefunction times electron 2's"
    because the electrons push on each other. To store such a function on a
    computer you lay down a grid in every one of those $3N$ directions, and
    the number of grid boxes multiplies. With just ten points per direction
    and ten electrons that is $10^{30}$ boxes — more numbers than could ever
    be written down. That impossibility, not any lack of cleverness, is why
    the rest of this book exists.

!!! tip "New vocabulary"
    - **Wavefunction** $\Psi$ — the object holding everything knowable about
      the system; $|\Psi|^2$ is a probability density. See the
      [beginner glossary](../undergraduate/glossary-for-beginners.md).
    - **Hamiltonian** $\hat H$ — the total-energy operator; solving
      $\hat H\Psi = E\Psi$ is "doing the quantum mechanics". See the
      [beginner glossary](../undergraduate/glossary-for-beginners.md).
    - **Configuration space** — the abstract space of *all* the particle
      coordinates at once; for $N$ particles in 3D it has $3N$ dimensions.
    - **Antisymmetric** — changes sign when you swap any two electrons.
    - **Indistinguishable** — there is no experiment that tells "electron 1"
      from "electron 2"; they carry no individual identity.
    - **Slater determinant** — the antisymmetric wavefunction you build by
      writing orbitals in a determinant; defined in Section 4.5.4.

    | Symbol | Meaning | Units (SI) |
    |---|---|---|
    | $N$, $N_{\mathrm e}$ | number of electrons | — (a count) |
    | $N_{\mathrm n}$ | number of nuclei | — |
    | $\mathbf r_i$ | position of electron $i$ | m |
    | $\mathbf R_I$ | position of nucleus $I$ | m |
    | $\Psi(\mathbf r_1,\dots,\mathbf r_N)$ | many-electron wavefunction | m$^{-3N/2}$ |
    | $\phi_a(\mathbf r)$ | single-electron orbital | m$^{-3/2}$ |
    | $\hat H$ | Hamiltonian (total-energy operator) | J |
    | $m_{\mathrm e}$ | electron mass, $9.109\times10^{-31}$ kg | kg |
    | $M_I$ | mass of nucleus $I$ | kg |
    | $Z_I$ | atomic number of nucleus $I$ | — |
    | $\hbar$ | reduced Planck constant, $1.055\times10^{-34}$ J s | J s |
    | $e$ | elementary charge, $1.602\times10^{-19}$ C | C |
    | $N_g$ | grid points per spatial direction | — |
    | $D$ | dimension of the discretised Hilbert space | — |
    | $\binom{n}{k}$ | binomial coefficient, "$n$ choose $k$" | — |

    A reminder on indices throughout this section: lower-case $i,j$ label
    **electrons**, capital $I,J$ label **nuclei**.

## 4.5.0 Warm-up: non-interacting electrons

Before writing the full Hamiltonian, consider the *limit* in which electron–electron interaction is switched off. The Hamiltonian for $N$ non-interacting electrons in an external potential $v_{\mathrm{ext}}(\mathbf r)$ is

$$\hat H_0 = \sum_{i=1}^N \hat h(\mathbf r_i), \qquad \hat h(\mathbf r) = -\frac{\hbar^2}{2m_{\mathrm e}}\nabla^2 + v_{\mathrm{ext}}(\mathbf r).$$

Because $\hat H_0$ is a *sum* of one-electron operators, each acting on a different coordinate, the eigenvalue equation $\hat H_0 \Psi = E \Psi$ separates. Try the product ansatz

$$\Psi(\mathbf r_1, \ldots, \mathbf r_N) = \phi_{a_1}(\mathbf r_1)\,\phi_{a_2}(\mathbf r_2)\cdots\phi_{a_N}(\mathbf r_N),$$

where each $\phi_{a_i}$ is a single-particle eigenstate of $\hat h$ with eigenvalue $\varepsilon_{a_i}$. Substituting,

$$\hat H_0\,\Psi = \sum_i \hat h(\mathbf r_i) \prod_j \phi_{a_j}(\mathbf r_j) = \sum_i \varepsilon_{a_i}\,\Psi = E\,\Psi,$$

so $E = \sum_i \varepsilon_{a_i}$ is the sum of one-electron eigenvalues. **The non-interacting $N$-electron problem reduces to $N$ independent single-electron problems** — exactly the kind of problem we solved in §4.3 and §4.4.

This is a wonderful state of affairs and it is the structure that underpins band theory (Chapter 3), tight-binding models, and ultimately the *Kohn–Sham* formulation of DFT, which represents an interacting system *as if it were* a non-interacting one in a cleverly chosen effective potential. The trouble, as we are about to see, is that real electrons interact through Coulomb repulsion, and the moment we turn that interaction back on the separability is destroyed.

!!! note "Why this step?"
    The factorisation $\hat H_0 = \sum_i \hat h_i$ works precisely because each $\hat h_i$ acts on coordinate $\mathbf r_i$ only, leaving the others untouched. The product wavefunction is then an eigenfunction by elementary calculus (a partial derivative with respect to $\mathbf r_i$ touches only $\phi_{a_i}$). Adding any term that couples two coordinates — even a single $V(\mathbf r_i, \mathbf r_j)$ for one pair — destroys separability and forces us to deal with the joint $3N$-dimensional problem.

??? note "Full derivation: why the product ansatz separates"
    Take the product trial function and apply the non-interacting Hamiltonian term by term. Write $\Psi = \prod_j \phi_{a_j}(\mathbf r_j)$. The operator $\hat h(\mathbf r_i)$ contains derivatives in $\mathbf r_i$ and a function of $\mathbf r_i$ only, so when it meets the product it passes straight through every factor whose coordinate is *not* $\mathbf r_i$:

    $$\hat h(\mathbf r_i)\,\Psi = \Bigl[\prod_{j\neq i}\phi_{a_j}(\mathbf r_j)\Bigr]\,\hat h(\mathbf r_i)\,\phi_{a_i}(\mathbf r_i).$$

    Because $\phi_{a_i}$ is by assumption an eigenfunction of the *same* one-electron operator $\hat h$ with eigenvalue $\varepsilon_{a_i}$, the bracketed action is just multiplication by a number:

    $$\hat h(\mathbf r_i)\,\phi_{a_i}(\mathbf r_i) = \varepsilon_{a_i}\,\phi_{a_i}(\mathbf r_i),$$

    so

    $$\hat h(\mathbf r_i)\,\Psi = \varepsilon_{a_i}\,\Psi.$$

    Summing over $i$,

    $$\hat H_0\,\Psi = \sum_{i=1}^N \hat h(\mathbf r_i)\,\Psi = \Bigl(\sum_{i=1}^N \varepsilon_{a_i}\Bigr)\Psi,$$

    which is exactly $E\Psi$ with $E=\sum_i\varepsilon_{a_i}$. The crucial property used is that each $\hat h(\mathbf r_i)$ ignores every coordinate except $\mathbf r_i$ — a property the Coulomb repulsion $1/|\mathbf r_i-\mathbf r_j|$ does *not* have, because it depends on two coordinates at once.

!!! warning "Common misunderstanding: do electrons move independently?"
    The clean result above is seductive, but it is the answer to a problem we
    have deliberately *changed*. Real electrons are charged, and the term
    $\hat V_{\mathrm{ee}}=\tfrac12\sum_{i\neq j}e^2/|\mathbf r_i-\mathbf r_j|$
    couples every coordinate to every other. Once it is present:

    - The wavefunction is **not** a product $\phi_{a_1}(\mathbf r_1)\cdots\phi_{a_N}(\mathbf r_N)$, and there is no set of one-electron orbitals whose product is the true ground state.
    - "Electron $i$ is in orbital $a_i$" stops being even approximately true; the electrons are *correlated* — where one goes depends on where the others are.
    - The total energy is **not** a sum of one-electron energies.

    When you later read that Hartree–Fock or Kohn–Sham DFT "puts each electron
    in its own orbital", remember that this is an *approximation chosen for
    tractability*, not a statement that electrons really move independently.
    The error you make by pretending they do is called the **correlation
    energy** (Section 4.5.4) — small in magnitude, but responsible for most
    of chemistry.

## 4.5.1 The full Hamiltonian for a real material

!!! note "Physical picture"
    Everything in a piece of matter that is not a nucleus is an electron, and
    every charge feels every other charge through the Coulomb force. The
    Hamiltonian below is just *the total energy of that swarm*, written one
    contribution at a time: the kinetic energy of the electrons, the kinetic
    energy of the nuclei, and then the three kinds of electrostatic energy —
    electron–electron repulsion, electron–nucleus attraction, and
    nucleus–nucleus repulsion. There is nothing exotic in it; it is the same
    $-\tfrac{\hbar^2}{2m}\nabla^2 + V$ you met for one particle, summed over
    every particle and every pair. What makes it hard is not any single term
    but that the repulsion terms *tie all the coordinates together*.

!!! tip "Gaussian units, briefly"
    The book writes the Coulomb energy of two charges $q_1,q_2$ a distance $r$
    apart as $q_1 q_2/r$. This is the **Gaussian** (CGS) convention, which
    absorbs the SI prefactor $1/4\pi\varepsilon_0$ into the definition of
    charge and keeps the formulae uncluttered. In SI units every Coulomb term
    below would carry an extra factor $1/4\pi\varepsilon_0$; for example
    $\hat V_{\mathrm{ee}} = \tfrac12\sum_{i\neq j} e^2/(4\pi\varepsilon_0|\mathbf r_i-\mathbf r_j|)$.
    The physics is identical — only the bookkeeping of constants differs. We
    follow the source text and use Gaussian units for the rest of this section.

Consider an arbitrary molecule or solid: $N_{\mathrm e}$ electrons (mass $m_{\mathrm e}$, charge $-e$, positions $\mathbf r_i$) and $N_{\mathrm n}$ nuclei (mass $M_I$, charge $+Z_I e$, positions $\mathbf R_I$). The non-relativistic Hamiltonian, in Gaussian units for compactness, is

$$\hat{H} = \hat T_{\mathrm e} + \hat T_{\mathrm n} + \hat V_{\mathrm{ee}} + \hat V_{\mathrm{en}} + \hat V_{\mathrm{nn}}, \tag{4.5.1}$$

with the five terms

$$\hat T_{\mathrm e} = -\frac{\hbar^2}{2m_{\mathrm e}}\sum_{i=1}^{N_{\mathrm e}} \nabla_i^2, \qquad \hat T_{\mathrm n} = -\sum_{I=1}^{N_{\mathrm n}} \frac{\hbar^2}{2M_I} \nabla_I^2, \tag{4.5.2}$$

$$\hat V_{\mathrm{ee}} = \frac{1}{2}\sum_{i\neq j} \frac{e^2}{|\mathbf r_i - \mathbf r_j|}, \qquad \hat V_{\mathrm{en}} = -\sum_{i, I}\frac{Z_I e^2}{|\mathbf r_i - \mathbf R_I|}, \qquad \hat V_{\mathrm{nn}} = \frac{1}{2}\sum_{I\neq J}\frac{Z_I Z_J e^2}{|\mathbf R_I - \mathbf R_J|}. \tag{4.5.3}$$

!!! note "Reading the indices carefully"
    The sums are over different kinds of indices and it is worth pausing to make sure you parse them correctly.
    - $\hat T_{\mathrm e}$ sums over electron labels $i = 1, \ldots, N_{\mathrm e}$. Each term involves the Laplacian $\nabla_i^2$ with respect to the position $\mathbf r_i$ of electron $i$.
    - $\hat T_{\mathrm n}$ sums over nuclear labels $I = 1, \ldots, N_{\mathrm n}$, with $\nabla_I^2$ acting on $\mathbf R_I$. Capital $I$ versus lower-case $i$ is the convention that distinguishes nuclei from electrons throughout this book.
    - $\hat V_{\mathrm{ee}}$ is a double sum over *pairs* of electrons. The factor of $1/2$ corrects double-counting since $(i,j)$ and $(j,i)$ refer to the same pair; $i \neq j$ excludes self-interaction.
    - $\hat V_{\mathrm{en}}$ couples every electron $i$ to every nucleus $I$. There is no factor of $1/2$ because $i$ and $I$ run over different sets.
    - $\hat V_{\mathrm{nn}}$ is a double sum over nuclei pairs; the factor of $1/2$ and $I \neq J$ play the same role as in $\hat V_{\mathrm{ee}}$.
    Memorise this index convention now. The same letters will appear in HF (§4.7) and in every DFT discussion thereafter.

??? note "Full derivation: where the factor of $\tfrac12$ comes from"
    The electron–electron repulsion is a sum over *pairs* of electrons, and each pair must be counted once. There are two equivalent ways to write a sum over pairs, and seeing both fixes the $\tfrac12$ for good.

    **Ordered double sum.** Let $i$ and $j$ each run over all electrons with the single restriction $i\neq j$. The pair $\{1,2\}$ is then hit *twice* — once as $(i,j)=(1,2)$ and once as $(i,j)=(2,1)$ — and both times it contributes the *same* energy, because $|\mathbf r_1-\mathbf r_2| = |\mathbf r_2-\mathbf r_1|$. To count each physical pair once we must therefore divide by $2$:

    $$\hat V_{\mathrm{ee}} = \frac12 \sum_{i\neq j}\frac{e^2}{|\mathbf r_i-\mathbf r_j|} = \frac{e^2}{2}\sum_{i=1}^{N_{\mathrm e}}\sum_{\substack{j=1\\ j\neq i}}^{N_{\mathrm e}}\frac{1}{|\mathbf r_i-\mathbf r_j|}.$$

    **Restricted single sum.** Equivalently, sum only over $i<j$, which visits each pair exactly once and needs *no* prefactor:

    $$\hat V_{\mathrm{ee}} = \sum_{i<j}\frac{e^2}{|\mathbf r_i-\mathbf r_j|}.$$

    The number of distinct pairs is $\binom{N_{\mathrm e}}{2} = \tfrac12 N_{\mathrm e}(N_{\mathrm e}-1)$, so for $N_{\mathrm e}=2$ there is exactly one repulsion term, for $N_{\mathrm e}=3$ there are three, and so on. The same argument fixes the $\tfrac12$ in $\hat V_{\mathrm{nn}}$. There is no $\tfrac12$ in $\hat V_{\mathrm{en}}$ because $i$ (an electron) and $I$ (a nucleus) range over *different* sets, so $(i,I)$ and $(I,i)$ are not "the same pair counted twice" — each electron–nucleus combination occurs once already.

That is it. *Every* property of every material — bond lengths, lattice constants, elastic moduli, band gaps, magnetisation, superconductivity, ferroelectricity, thermal conductivity — is encoded in solving the eigenvalue equation $\hat{H} \Psi = E \Psi$ for this operator. Dirac, having written down a similar Hamiltonian in 1929, declared:

> "The underlying physical laws necessary for the mathematical theory of a large part of physics and the whole of chemistry are thus completely known, and the difficulty is only that the exact application of these laws leads to equations much too complicated to be soluble."

The history of computational materials science is the history of grappling with that "only". In this section we measure the size of the problem and convince ourselves that it cannot be solved by brute force.

## 4.5.2 The wavefunction lives in $3N$-dimensional space

In §4.3 we discretised the wavefunction of a single particle on a 1D grid of, say, $N_g = 400$ points. The wavefunction was represented by 400 complex numbers. The Hamiltonian was a $400\times 400$ matrix — 160 000 entries, of which only $\sim 1200$ were non-zero in the tridiagonal structure. Trivial.

!!! note "Physical picture: what does '$3N$-dimensional space' mean?"
    The space here is **not** the room the electrons sit in — that is always
    just ordinary 3D space. It is *configuration space*: one axis for **every
    coordinate of every electron**. To specify the state of two electrons you
    must give six numbers at once — $(x_1,y_1,z_1,x_2,y_2,z_2)$ — and a single
    point in the 6-dimensional configuration space is one complete "snapshot"
    saying where *both* electrons are simultaneously. The wavefunction assigns
    one complex amplitude to each such snapshot, i.e. to each point of this
    big space.

    For one electron the wavefunction is a cloud you could photograph: a
    function over 3D. For two electrons it is a function over 6D, which you
    cannot draw, because $\Psi(\mathbf r_1,\mathbf r_2)$ tells you how the
    *joint* arrangement is weighted — it knows that if electron 1 is here,
    electron 2 prefers to be over there. That correlation is precisely the
    information a product $\phi(\mathbf r_1)\phi(\mathbf r_2)$ throws away, and
    it is why we cannot collapse the $3N$ axes back down to $3$. The grid we
    are about to lay down has to tile *all* $3N$ axes, and that is the source
    of the catastrophe.

Now consider $N$ electrons in three dimensions. The many-electron wavefunction is a function

$$\Psi(\mathbf r_1, \mathbf r_2, \ldots, \mathbf r_N) \tag{4.5.4}$$

of $3N$ continuous coordinates (we are ignoring spin for the moment). To sample $\Psi$ on a grid with $N_g$ points per spatial direction, we need

$$N_g^{3N} \tag{4.5.5}$$

complex numbers. The Hamiltonian becomes a matrix of size $N_g^{3N} \times N_g^{3N}$. This is the **exponential wall**.

!!! example "Minimal example: two electrons on a 3-point grid"
    Strip the problem to almost nothing. Put space on a grid of just **three
    points** in each of $x,y,z$, so a single electron can sit at any of
    $N_g^3 = 3^3 = 27$ locations. One electron is therefore described by $27$
    amplitudes — a vector of length $27$, completely manageable.

    Now add a *second* electron. A configuration is now a pair "(box for
    electron 1, box for electron 2)", and the number of pairs is

    $$27 \times 27 = 27^2 = 729.$$

    The wavefunction is a $729$-entry array; the Hamiltonian is a
    $729\times729$ matrix. Already we cannot picture it, but a laptop copes.
    The pattern is the key thing: **each electron we add multiplies the count
    by $27$ again.** For $N$ electrons on this toy grid the number of
    amplitudes is

    $$\underbrace{27\times27\times\cdots\times27}_{N\ \text{factors}} = 27^{N} = \bigl(N_g^3\bigr)^N = N_g^{3N},$$

    which is exactly Eq. (4.5.5). Three electrons need $27^3 = 19\,683$;
    five electrons need $27^5 = 14\,348\,907$ (already 14 million on a grid so
    coarse it is useless for chemistry). The multiplication never stops, and
    that *repeated multiplication* — one factor of $N_g^3$ per electron — is
    what the word "exponential" means here. Refine the grid from 3 points to
    a realistic 10 points per direction and the base jumps from $27$ to
    $1000$, giving the table below.

Let us put numbers on it. Take a coarse spatial grid of $N_g = 10$ points per dimension — laughably under-resolved for any real chemistry, but a good lower bound. Then:

| $N$ electrons | grid coords $3N$ | basis states $10^{3N}$ |
|---:|---:|---:|
| 1 | 3 | $10^{3}$ |
| 2 | 6 | $10^{6}$ |
| 3 | 9 | $10^{9}$ |
| 5 | 15 | $10^{15}$ |
| 10 | 30 | $10^{30}$ |
| 20 | 60 | $10^{60}$ |
| 100 | 300 | $10^{300}$ |

The third row, three electrons, already gives a billion-dimensional Hilbert space. The fifth row, ten electrons — roughly a water molecule — gives $10^{30}$ basis functions. A single double-precision complex number occupies 16 bytes, so storing $\Psi$ for ten electrons on a 10×10×10 grid would require $1.6 \times 10^{31}$ bytes. The total digital data created by humanity to date is approximately $10^{23}$ bytes. Just *writing down* the wavefunction of ten electrons on a coarse 1000-point grid would require **a hundred million times the entire world's data storage**.

!!! example "Putting the absurdity in context"
    Take 10 electrons (a water molecule) on a $10\times 10\times 10$ grid: $10^{30}$ basis states. The number of atoms in the Earth is approximately $1.3\times 10^{50}$. So the wavefunction of *water* on a 1000-point coarse grid has fewer entries than there are atoms in the Earth — but only by twenty orders of magnitude, which on a logarithmic scale leaves us essentially nowhere. Refining the grid to chemical resolution (100 points per direction) brings the basis up to $10^{60}$, which is more than the atoms in the *solar system* ($\sim 10^{57}$). The brute-force approach is not just impractical, it is *cosmologically* impossible.

!!! warning "It gets worse"
    Even if we could store $\Psi$, we would still need to diagonalise the Hamiltonian. Dense diagonalisation costs $\mathcal O(D^3)$ operations for a $D\times D$ matrix; sparse Lanczos methods can reach $\mathcal O(D \cdot \text{iterations})$ but still need to apply $\hat{H}$ to a vector of length $D$. For our ten-electron problem, $D = 10^{30}$ — and the fastest supercomputers in 2026 perform of order $10^{18}$ floating-point operations per second. A single matrix–vector product would take $10^{12}$ seconds, roughly thirty thousand years.

This is **not** a problem that will be solved by Moore's law. Doubling our compute power every two years lets us add one electron to the calculation every six years or so. At that rate, going from 10 electrons to 30 (a small organic molecule) would take 120 years of hardware improvement. The Schrödinger equation must be tackled differently.

??? note "Where does 'one electron every six years' come from?"
    Adding one electron multiplies the number of grid amplitudes by $N_g^3$.
    On the $N_g=10$ grid that is a factor of $10^3 = 1000$. Moore's law gives
    a factor-of-2 speed-up roughly every two years, so to gain a factor of
    $1000 = 2^{\log_2 1000} = 2^{\,\approx 10}$ you need about $10$ doublings,
    i.e. about $10\times 2 = 20$ years per electron at fixed grid. (The text's
    "every six years" assumes the coarser toy grid $N_g=3$, where one extra
    electron costs a factor $27\approx 2^{4.75}$, i.e. about $5$ doublings or
    $\sim10$ years; either way the conclusion stands — hardware progress is
    *linear in electrons added* against an *exponential wall*, so it never
    catches up.) The point is not the exact slope but that no realistic rate
    of hardware improvement closes an exponential gap.

!!! question "Check yourself: the exponential wall (pause and recall)"
    Answer these before reading on — they are the load-bearing ideas of the
    whole chapter.

    1. On a grid with $N_g$ points per direction, how many complex amplitudes
       are needed to store the wavefunction of $N$ electrons (ignoring spin)?
       Why is it a *power* of $N_g^3$ and not a multiple of it?
    2. Using a $10$-point grid, how many amplitudes does a single water
       molecule's 10 electrons require, and how does that compare to the
       roughly $10^{23}$ bytes of all digital data ever created?
    3. A colleague says: "We just need a computer ten times bigger to do ten
       electrons instead of nine." What is wrong with this?
    4. Roadmap question: given that brute force is impossible, name the two
       structural facts the rest of the section is about to exploit, and
       which later method ([Chapter 5](../ch05-dft/index.md)) they make
       necessary.

    ??? success "Answer"
        1. $N_g^{3N}$ amplitudes (Eq. 4.5.5). Each electron contributes its
           own factor of $N_g^3$ (one for each of its three coordinate axes,
           each sampled at $N_g$ points), and these factors *multiply*
           because every combination of "where electron 1 is" with "where
           electron 2 is" with … is a distinct configuration. Multiplying a
           fixed base by itself $N$ times is exponential growth in $N$.
        2. $10$ electrons on a $10\times10\times10$ grid is
           $N_g^{3N}=10^{3\times10}=10^{30}$ amplitudes; at $16$ bytes each
           that is $1.6\times10^{31}$ bytes, about $10^{8}$ — a hundred
           million — times all the digital data humanity has ever produced.
        3. The cost is *exponential*, not linear. Going from $9$ to $10$
           electrons multiplies the storage by $N_g^3 = 1000$ on a 10-point
           grid, not by $10$. A computer "ten times bigger" buys you almost
           nothing; you would need one *a thousand times* bigger just for one
           more electron, and another thousandfold for the next.
        4. (i) The wavefunction lives in $3N$-dimensional configuration space
           (Section 4.5.2); (ii) for electrons it must be antisymmetric under
           exchange (Section 4.5.3), which forces the Slater-determinant
           structure (Section 4.5.4). Because storing the exact wavefunction
           is impossible, we are driven to recast the problem in terms of the
           electron *density* — density functional theory,
           [Chapter 5 (DFT)](../ch05-dft/index.md).

    ??? note "Hint"
        For 1, count the axes: each electron carries three coordinates, each
        coordinate is sampled at $N_g$ values, and a configuration fixes them
        all at once. For 3, compare the *factor* you gain (computer size) with
        the *factor* the problem grows by when you add an electron.

## 4.5.3 Pauli antisymmetry

!!! note "Plain-language version: why a sign at all?"
    Electrons are *identical* in the strong quantum sense: no measurement, not
    even in principle, can tell one from another. So if we physically swap two
    electrons, every observable quantity — every probability, every energy —
    must be completely unchanged. Now, all observables come from
    $|\Psi|^2$, never from $\Psi$ itself. Demanding that $|\Psi|^2$ be
    unchanged by a swap does *not* force $\Psi$ to be unchanged: it only
    forces $\Psi$ to pick up a phase factor whose square modulus is one. The
    short argument below shows that for a swap this phase can only be $+1$ or
    $-1$, and that Nature assigns $-1$ to electrons. That single minus sign is
    the Pauli principle, and it controls the entire structure of the periodic
    table.

??? note "Full derivation: indistinguishability forces a sign of $\pm1$"
    Define the **exchange operator** $\hat P_{ij}$ that swaps the coordinates of particles $i$ and $j$:

    $$\hat P_{ij}\,\Psi(\dots,\mathbf r_i,\dots,\mathbf r_j,\dots) = \Psi(\dots,\mathbf r_j,\dots,\mathbf r_i,\dots).$$

    **Step 1 — indistinguishability fixes $|\,$amplitude$\,|$.** Because the two particles are genuinely identical, the swapped configuration is physically the same state, so the probability density cannot change:

    $$\bigl|\hat P_{ij}\Psi\bigr|^2 = |\Psi|^2 .$$

    The most general way two functions can have equal squared modulus everywhere is that they differ by a constant phase $\lambda$ with $|\lambda|=1$:

    $$\hat P_{ij}\,\Psi = \lambda\,\Psi, \qquad |\lambda| = 1. \tag{4.5.6a}$$

    In other words $\Psi$ must be an *eigenfunction* of the exchange operator. (Strictly this is exact only when particles $i,j$ are equivalent throughout, which is the situation we care about; the conclusion is borne out by the full theory.)

    **Step 2 — swapping twice is the identity.** Exchanging the same two particles a second time restores the original labelling:

    $$\hat P_{ij}^{\,2}\,\Psi = \Psi .$$

    But applying Eq. (4.5.6a) twice gives $\hat P_{ij}^{\,2}\Psi = \lambda^2\Psi$. Therefore

    $$\lambda^2 = 1 \quad\Longrightarrow\quad \lambda = \pm 1 .$$

    Only two possibilities survive. (Note we needed $\lambda^2=1$, *not* merely $|\lambda|^2=1$; this is what rules out a generic complex phase and leaves exactly $\pm1$ in three dimensions.)

    **Step 3 — Nature's choice.** Which sign a given particle takes is *not* deducible from non-relativistic quantum mechanics; it is fixed by the **spin–statistics theorem** of relativistic quantum field theory and taken as a postulate here:

    - $\lambda = +1$ for **bosons** (integer spin: photons, $^4$He, …) — *symmetric* wavefunctions;
    - $\lambda = -1$ for **fermions** (half-integer spin: electrons, protons, neutrons, …) — *antisymmetric* wavefunctions.

    Setting $\lambda=-1$ in Eq. (4.5.6a) is exactly Eq. (4.5.6) below. So antisymmetry is *indistinguishability* (which gives $\lambda=\pm1$) plus the *empirical fact* that electrons are fermions (which selects $-1$).

Worse: a generic function $\Psi(\mathbf r_1, \ldots, \mathbf r_N)$ is not even an admissible wavefunction for electrons. The Pauli principle, stated formally, is the postulate that the wavefunction of identical *fermions* (electrons, protons, neutrons, ...) must change sign under exchange of any two particles:

$$\Psi(\ldots, \mathbf r_i, \ldots, \mathbf r_j, \ldots) = -\Psi(\ldots, \mathbf r_j, \ldots, \mathbf r_i, \ldots). \tag{4.5.6}$$

(Strictly the exchange acts on combined space-and-spin coordinates $\mathbf x_i = (\mathbf r_i, \sigma_i)$, but in this section we will keep spin implicit.) Identical *bosons* (photons, $^4$He nuclei, …) take a plus sign instead. This is the **spin-statistics theorem**, derived in relativistic quantum field theory and accepted as a postulate in non-relativistic quantum mechanics.

An immediate corollary: if two electrons are in the same single-particle state ($\mathbf r_i = \mathbf r_j$ and same spin), then (4.5.6) demands $\Psi = -\Psi$, so $\Psi = 0$. This is the exclusion principle: *no two electrons can occupy the same one-electron state*. Stuffing more than two electrons into the lowest level of a particle-in-a-box (one spin-up, one spin-down) is forbidden; the third electron must go in $n = 2$. This is what makes atoms beyond hydrogen possess shell structure, and what makes the periodic table look the way it does.

!!! warning "Common misunderstandings: what antisymmetry does and does not say"
    - **It is a statement about the wavefunction, not about forces.**
      Antisymmetry holds *even with the Coulomb repulsion switched off*. The
      "exchange hole" (below) is a consequence of the minus sign, not of
      electrostatic repulsion — same-spin electrons avoid each other purely
      because $\Psi$ must vanish when their space-and-spin coordinates
      coincide.
    - **The minus sign attaches to a swap of two particles, not to a single
      particle.** $\Psi$ does not "have a sign"; only the *relationship*
      between $\Psi$ before and after exchanging a pair does.
    - **The exclusion principle is the special case of antisymmetry, not an
      extra rule.** "No two electrons in the same state" follows by setting
      the two coordinates equal in Eq. (4.5.6); it is not an independent
      postulate.
    - **It does not forbid two electrons from being at the same point in
      space.** Two electrons of *opposite* spin can sit at the same
      $\mathbf r$ — their antisymmetry is carried by the spin part. Only
      *same-spin* electrons are kept apart in space. (This is why each
      spatial orbital holds two electrons, not one.)
    - **Antisymmetry is not the same as electrons "repelling because they are
      negative".** Coulomb repulsion acts between *all* electrons regardless
      of spin; the exchange effect acts only between *parallel-spin*
      electrons and exists even for neutral, non-interacting fermions.

### Why antisymmetry — a two-electron example

It is worth seeing antisymmetry in action on the smallest possible system. Take two electrons in single-particle orbitals $\phi_a, \phi_b$ (orthonormal, both spatial only — we suppress spin for now). The product ansatz $\Psi_{\mathrm{prod}}(\mathbf r_1, \mathbf r_2) = \phi_a(\mathbf r_1) \phi_b(\mathbf r_2)$ treats the electrons as distinguishable: electron 1 is in $\phi_a$, electron 2 is in $\phi_b$. But electrons are *identical* — there is no fact of the matter about "which one is which" — and the wavefunction must reflect this.

The two possibilities consistent with identicalness are the *symmetric* and *antisymmetric* combinations:

$$\Psi_\pm(\mathbf r_1, \mathbf r_2) = \frac{1}{\sqrt 2}\bigl[\phi_a(\mathbf r_1)\phi_b(\mathbf r_2) \pm \phi_b(\mathbf r_1)\phi_a(\mathbf r_2)\bigr].$$

Under exchange $\mathbf r_1 \leftrightarrow \mathbf r_2$, $\Psi_+$ is unchanged (bosons) and $\Psi_-$ acquires a minus sign (fermions, including electrons).

Now evaluate the joint probability density at $\mathbf r_1 = \mathbf r_2 \equiv \mathbf r$. For the antisymmetric combination,

$$\Psi_-(\mathbf r, \mathbf r) = \frac{1}{\sqrt 2}\bigl[\phi_a(\mathbf r)\phi_b(\mathbf r) - \phi_b(\mathbf r)\phi_a(\mathbf r)\bigr] = 0.$$

The probability of finding both electrons at the same point is exactly zero, regardless of the orbitals $\phi_a, \phi_b$. This is the "exchange hole" — a quantum-mechanical avoidance, present *even without Coulomb repulsion*, between electrons of the same spin. It is the origin of why electrons in different shells of an atom don't collapse together: Pauli, not Coulomb, holds them apart.

For the symmetric (bosonic) combination, by contrast, $\Psi_+(\mathbf r, \mathbf r) = \sqrt 2 \phi_a(\mathbf r)\phi_b(\mathbf r) \neq 0$ — bosons *like* to be together (bunching).

This single calculation, no more elaborate than the algebra above, is the seed from which exchange interactions in magnetism, the Slater determinant in §4.5.4, and the Fermi sea in solid-state physics all grow.

!!! tip "Spin enters the picture"
    With spin included, the full wavefunction is antisymmetric under exchange of the combined space-spin coordinates. This permits two electrons of *opposite* spin to occupy the same spatial orbital (the spatial part is symmetric, the spin part antisymmetric — the *singlet*); but two electrons of *parallel* spin must have an antisymmetric spatial part (and a symmetric *triplet* spin), forcing them into different spatial orbitals. The pattern "two opposite-spin electrons per spatial orbital" is what makes the periodic table look the way it does.

For numerical purposes, antisymmetry restricts the wavefunction to a subspace — the antisymmetric subspace of $L^2(\mathbb R^{3N})$. Its dimension is much smaller than $N_g^{3N}$ (only $\binom{N_g^3}{N}$ ways to occupy $N_g^3$ orbitals with $N$ electrons), but for any non-trivial $N$ this is still astronomical.

!!! example "Minimal example: how much does antisymmetry save?"
    Return to the 3-point grid, where there are $N_g^3 = 27$ single-particle
    boxes, and put two same-spin electrons on it. Counting *all* configurations
    (treating the electrons as distinguishable) gave $27^2 = 729$. But the
    antisymmetric subspace only allows *unordered, distinct* occupations: pick
    $2$ of the $27$ boxes, in no particular order, and never both in the same
    box. The count is

    $$\binom{27}{2} = \frac{27\times 26}{2} = 351.$$

    So antisymmetry removes the $27$ "both electrons in the same box" states
    (forbidden by exclusion) and folds the rest in half (because swapping the
    two labels gives the same physical state), $729 \to (729-27)/2 = 351$. A
    real saving — but only a factor of about two, and it still grows like
    $\binom{N_g^3}{N}$, which for $N_g^3=1000$ and $N=10$ is about $10^{23}$.
    Antisymmetry shrinks the wall; it does not knock it down. The factor of
    $1/N!$ in the Slater determinant of Section 4.5.4 is exactly the
    normalisation bookkeeping behind this "fold in half" (here $1/2! = 1/2$).

### Worked example: helium singlet–triplet splitting

The simplest place in nature where antisymmetry has a *measurable* consequence — beyond the structure of the periodic table itself — is the excited-state spectrum of helium. He has two electrons; in its ground state both occupy the $1s$ orbital, one with spin up and one with spin down. Excite one electron to the $2s$ orbital and you have an electron in $1s$ and an electron in $2s$, with two possible spin configurations:

- **Para-helium (singlet, $S = 0$).** Spatial part symmetric, $\Psi_{\mathrm{sym}}^{\mathrm{space}}(\mathbf r_1, \mathbf r_2) = \tfrac{1}{\sqrt 2}[\phi_{1s}(\mathbf r_1)\phi_{2s}(\mathbf r_2) + \phi_{2s}(\mathbf r_1)\phi_{1s}(\mathbf r_2)]$; spin part antisymmetric, $\tfrac{1}{\sqrt 2}(|\!\uparrow\downarrow\rangle - |\!\downarrow\uparrow\rangle)$.
- **Ortho-helium (triplet, $S = 1$).** Spatial part antisymmetric, $\Psi_{\mathrm{anti}}^{\mathrm{space}} = \tfrac{1}{\sqrt 2}[\phi_{1s}(\mathbf r_1)\phi_{2s}(\mathbf r_2) - \phi_{2s}(\mathbf r_1)\phi_{1s}(\mathbf r_2)]$; spin part one of the three symmetric combinations $|\!\uparrow\uparrow\rangle$, $|\!\downarrow\downarrow\rangle$, $\tfrac{1}{\sqrt 2}(|\!\uparrow\downarrow\rangle + |\!\downarrow\uparrow\rangle)$.

Now compute the Coulomb repulsion expectation value $\langle \Psi | e^2/r_{12} | \Psi\rangle$ for the two cases. After a few lines of algebra,

$$E_\pm = J \pm K, \qquad J = \int |\phi_{1s}(\mathbf r_1)|^2 |\phi_{2s}(\mathbf r_2)|^2\,\frac{e^2}{r_{12}}\,d\mathbf r_1 d\mathbf r_2,$$
$$K = \int \phi_{1s}^*(\mathbf r_1)\phi_{2s}(\mathbf r_1)\,\frac{e^2}{r_{12}}\,\phi_{2s}^*(\mathbf r_2)\phi_{1s}(\mathbf r_2)\,d\mathbf r_1 d\mathbf r_2.$$

The "+" applies to the singlet (symmetric spatial), the "−" to the triplet (antisymmetric spatial). The first integral $J$ is the classical Coulomb energy; the second integral $K$ is the **exchange integral** and is purely quantum-mechanical. Both $J$ and $K$ are positive.

??? note "Full derivation: the 'few lines of algebra' giving $E_\pm = J \pm K$"
    The Coulomb operator $e^2/r_{12}$ acts on the spatial coordinates only, so the spin part of $\Psi$ — being normalised — contributes a factor of $1$ and can be dropped from this expectation value. We work with the normalised spatial parts

    $$\Psi^{\mathrm{space}}_\pm(\mathbf r_1,\mathbf r_2) = \frac{1}{\sqrt2}\bigl[\phi_{1s}(\mathbf r_1)\phi_{2s}(\mathbf r_2) \pm \phi_{2s}(\mathbf r_1)\phi_{1s}(\mathbf r_2)\bigr],$$

    with $+$ for the singlet and $-$ for the triplet. Abbreviate $a \equiv \phi_{1s}$, $b \equiv \phi_{2s}$, and write $1,2$ for $\mathbf r_1,\mathbf r_2$. Then

    $$E_\pm = \big\langle \Psi^{\mathrm{space}}_\pm \big|\, \tfrac{e^2}{r_{12}} \,\big| \Psi^{\mathrm{space}}_\pm \big\rangle = \frac{1}{2}\int \bigl[a(1)b(2) \pm b(1)a(2)\bigr]^*\,\frac{e^2}{r_{12}}\,\bigl[a(1)b(2) \pm b(1)a(2)\bigr]\,d1\,d2.$$

    **Step 1 — expand the bracketed product into four terms.** The $1/2$ multiplies a sum of four integrals:

    $$E_\pm = \frac12\bigl( I_1 + I_2 \pm I_3 \pm I_4 \bigr),$$

    with

    $$I_1 = \int |a(1)|^2\,\frac{e^2}{r_{12}}\,|b(2)|^2\,d1\,d2, \qquad I_2 = \int |b(1)|^2\,\frac{e^2}{r_{12}}\,|a(2)|^2\,d1\,d2,$$
    $$I_3 = \int a^*(1)b^*(2)\,\frac{e^2}{r_{12}}\,b(1)a(2)\,d1\,d2, \qquad I_4 = \int b^*(1)a^*(2)\,\frac{e^2}{r_{12}}\,a(1)b(2)\,d1\,d2.$$

    **Step 2 — identify the direct terms.** $I_1$ and $I_2$ are the *same* number: $I_2$ is just $I_1$ with the integration labels $1\leftrightarrow 2$ relabelled (the integrand and $r_{12}=r_{21}$ are symmetric under it). Each equals the classical Coulomb repulsion between the charge clouds $|a|^2$ and $|b|^2$:

    $$I_1 = I_2 = \int |\phi_{1s}(\mathbf r_1)|^2\,\frac{e^2}{r_{12}}\,|\phi_{2s}(\mathbf r_2)|^2\,d\mathbf r_1 d\mathbf r_2 \equiv J.$$

    **Step 3 — identify the exchange terms.** $I_3$ and $I_4$ are likewise equal to each other (relabel $1\leftrightarrow2$, or note $I_4 = I_3^*$ and the integral is real for real orbitals). Each equals the exchange integral:

    $$I_3 = I_4 = \int \phi_{1s}^*(\mathbf r_1)\phi_{2s}(\mathbf r_1)\,\frac{e^2}{r_{12}}\,\phi_{2s}^*(\mathbf r_2)\phi_{1s}(\mathbf r_2)\,d\mathbf r_1 d\mathbf r_2 \equiv K.$$

    Note carefully: in $J$ the *same* orbital is paired with itself at each point ($|a(1)|^2$, $|b(2)|^2$), whereas in $K$ the two orbitals are *swapped* between the bra and the ket ($a^*(1)b(1)$ at point 1) — this orbital interchange is what makes $K$ a purely quantum, "exchange" object with no classical analogue.

    **Step 4 — collect.** Substituting $I_1=I_2=J$ and $I_3=I_4=K$:

    $$E_\pm = \frac12\bigl( J + J \pm K \pm K \bigr) = \frac12\bigl( 2J \pm 2K \bigr) = J \pm K.$$

    The factor of $\tfrac12$ out front is exactly cancelled by the two equal copies of each integral — which is why the normalisation $1/\sqrt2$ leaves a clean $J\pm K$. The singlet ($+$) sits at $J+K$, the triplet ($-$) at $J-K$, and since $K>0$ the triplet lies lower by $2K$.

The triplet, $E_- = J - K$, lies *below* the singlet, $E_+ = J + K$, by $2K \approx 0.8$ eV in helium — an enormous splitting, measurable spectroscopically (the 1s2s $^1S$ and $^3S$ levels of He) and entirely due to exchange. The triplet wins because the antisymmetric spatial wavefunction has a node at $\mathbf r_1 = \mathbf r_2$ (we saw this in §4.5.3), so the two electrons avoid each other and pay less Coulomb energy.

!!! tip "Hund's rule"
    This is the simplest example of **Hund's rule**: for a given electronic configuration, the term with the *highest spin* lies lowest in energy, because parallel-spin electrons are kept apart by antisymmetry and so pay less Coulomb repulsion. The same mechanism makes magnetism possible in transition metals and rare earths.

## 4.5.4 Slater determinants

!!! info "What problem are we solving?"
    We have just learned that an electronic wavefunction *must* be
    antisymmetric (Section 4.5.3), and that a plain product
    $\phi_a(\mathbf r_1)\phi_b(\mathbf r_2)\cdots$ is *not*. We need a
    recipe that takes any set of one-electron orbitals and assembles them into
    a guaranteed-antisymmetric many-electron wavefunction, with the right
    normalisation, automatically. The **Slater determinant** is that recipe.
    It is the single most important building block in electronic-structure
    theory: Hartree–Fock uses exactly one of them (Section 4.7), and every
    correlated method (FCI, coupled cluster) is a sum of many.

!!! tip "New vocabulary"
    - **Spin-orbital** — a single-electron state that specifies *both* a
      spatial orbital and a spin, e.g. $\chi(\mathbf x) = \phi(\mathbf r)\,\alpha(\sigma)$,
      where $\mathbf x = (\mathbf r,\sigma)$ packages position and spin and
      $\alpha,\beta$ denote spin-up and spin-down. The full antisymmetry of
      Section 4.5.3 is a statement about the combined coordinate $\mathbf x$.
    - **Determinant** — the familiar alternating sum from linear algebra
      (Chapter 0); its defining property, that swapping two rows flips its
      sign, is *precisely* the antisymmetry we need.
    - **Normalisation factor** $1/\sqrt{N!}$ — the constant that makes
      $\int|\Phi|^2 = 1$; derived explicitly below.

The simplest antisymmetric wavefunctions are **Slater determinants**. Given $N$ orthonormal one-electron orbitals $\phi_1, \phi_2, \ldots, \phi_N$ (each a function of one spatial coordinate, and for now a single spin), build

$$\Phi(\mathbf r_1, \ldots, \mathbf r_N) = \frac{1}{\sqrt{N!}}\,\det\!\begin{pmatrix} \phi_1(\mathbf r_1) & \phi_2(\mathbf r_1) & \cdots & \phi_N(\mathbf r_1)\\ \phi_1(\mathbf r_2) & \phi_2(\mathbf r_2) & \cdots & \phi_N(\mathbf r_2)\\ \vdots & \vdots & \ddots & \vdots\\ \phi_1(\mathbf r_N) & \phi_2(\mathbf r_N) & \cdots & \phi_N(\mathbf r_N)\end{pmatrix}. \tag{4.5.7}$$

!!! example "Step-by-step: the Slater determinant for two electrons"
    The general formula (4.5.7) is easiest to absorb at $N=2$. Take two
    orthonormal spin-orbitals $\chi_a$ and $\chi_b$ (each fixes an orbital and
    a spin), and two electrons with combined coordinates $\mathbf x_1,\mathbf x_2$.
    The $N=2$ case of Eq. (4.5.7) is the $2\times2$ determinant

    $$\Phi(\mathbf x_1,\mathbf x_2) = \frac{1}{\sqrt{2!}}\,\det\!\begin{pmatrix} \chi_a(\mathbf x_1) & \chi_b(\mathbf x_1)\\[2pt] \chi_a(\mathbf x_2) & \chi_b(\mathbf x_2)\end{pmatrix}. \tag{4.5.7a}$$

    1. **Expand the $2\times2$ determinant.** Using $\det\begin{pmatrix}p&q\\ r&s\end{pmatrix}=ps-qr$,

       $$\Phi(\mathbf x_1,\mathbf x_2) = \frac{1}{\sqrt 2}\bigl[\chi_a(\mathbf x_1)\,\chi_b(\mathbf x_2) - \chi_b(\mathbf x_1)\,\chi_a(\mathbf x_2)\bigr]. \tag{4.5.7b}$$

       This is exactly the antisymmetric combination $\Psi_-$ we met by hand in Section 4.5.3 — the determinant has *reproduced* it, with no guesswork.

    2. **Swapping the two electrons flips the sign.** Exchange $\mathbf x_1 \leftrightarrow \mathbf x_2$ in (4.5.7b):

       $$\Phi(\mathbf x_2,\mathbf x_1) = \frac{1}{\sqrt 2}\bigl[\chi_a(\mathbf x_2)\chi_b(\mathbf x_1) - \chi_b(\mathbf x_2)\chi_a(\mathbf x_1)\bigr] = -\,\Phi(\mathbf x_1,\mathbf x_2).$$

       In determinant language this is the row swap: exchanging the electrons exchanges *rows* 1 and 2 of the matrix in (4.5.7a), and a determinant changes sign under a row swap (Chapter 0). Antisymmetry, Eq. (4.5.6), is automatic — it is built into the algebra of determinants, not imposed afterwards.

    3. **Equal orbitals give zero (exclusion).** Put $\chi_b = \chi_a$. Then the two *columns* of (4.5.7a) are identical, so $\det = 0$, and directly from (4.5.7b),
       $$\Phi = \frac{1}{\sqrt2}\bigl[\chi_a(\mathbf x_1)\chi_a(\mathbf x_2) - \chi_a(\mathbf x_1)\chi_a(\mathbf x_2)\bigr] = 0.$$
       Two electrons cannot share a spin-orbital — Pauli exclusion, with no extra assumption.

??? note "Full derivation: the $1/\sqrt{N!}$ normalisation"
    Why $1/\sqrt{2}$ for two electrons, and $1/\sqrt{N!}$ in general? Because the un-normalised determinant is a sum of $N!$ orthonormal-product terms, and normalising the *sum* of $N!$ orthonormal pieces brings in $1/\sqrt{N!}$.

    **Two-electron case.** Drop the prefactor and write $\tilde\Phi = \chi_a(\mathbf x_1)\chi_b(\mathbf x_2) - \chi_b(\mathbf x_1)\chi_a(\mathbf x_2)$. Its squared norm is

    $$\langle\tilde\Phi|\tilde\Phi\rangle = \int |\tilde\Phi|^2\,d\mathbf x_1 d\mathbf x_2,$$

    where $\int d\mathbf x$ means integrate over space *and* sum over spin. Expand the modulus-squared into four terms:

    $$|\tilde\Phi|^2 = |\chi_a(\mathbf x_1)|^2|\chi_b(\mathbf x_2)|^2 + |\chi_b(\mathbf x_1)|^2|\chi_a(\mathbf x_2)|^2 - \chi_a^*(\mathbf x_1)\chi_b^*(\mathbf x_2)\chi_b(\mathbf x_1)\chi_a(\mathbf x_2) - (\text{c.c.}).$$

    Integrate term by term, using orthonormality $\int \chi_p^*(\mathbf x)\chi_q(\mathbf x)\,d\mathbf x = \delta_{pq}$:

    - First term: $\bigl(\int|\chi_a|^2\bigr)\bigl(\int|\chi_b|^2\bigr) = 1\cdot 1 = 1.$
    - Second term: likewise $=1.$
    - Third (cross) term: $\bigl(\int \chi_a^*\chi_b\,d\mathbf x_1\bigr)\bigl(\int \chi_b^*... \bigr)$ — each factor is $\int\chi_a^*\chi_b = \delta_{ab} = 0$ since $a\neq b$, so the whole term vanishes. Its complex conjugate vanishes too.

    Hence $\langle\tilde\Phi|\tilde\Phi\rangle = 1 + 1 - 0 - 0 = 2 = 2!$. To make $\langle\Phi|\Phi\rangle = 1$ we divide by $\sqrt{2} = \sqrt{2!}$, which is exactly the factor in (4.5.7a).

    **General $N$.** Expanding the $N\times N$ determinant gives $N!$ terms, one for each permutation $P$ of the orbital labels:

    $$\tilde\Phi = \sum_{P}\mathrm{sgn}(P)\,\prod_{i=1}^N \chi_{P(i)}(\mathbf x_i).$$

    These $N!$ products are mutually orthonormal: two different permutations differ in at least one factor $\int\chi_p^*\chi_q = \delta_{pq}$, which is zero for $p\neq q$, killing every cross term. Each surviving "diagonal" term integrates to $1$, and there are $N!$ of them, so $\langle\tilde\Phi|\tilde\Phi\rangle = N!$. Dividing by $\sqrt{N!}$ normalises it — and this is the same $N!$ that appeared as the "fold by $1/N!$" in the antisymmetric-subspace count of Section 4.5.3.

Two structural facts make Slater determinants the workhorse of quantum chemistry.

1. **They are automatically antisymmetric.** Swapping any two electrons $i, j$ swaps two rows of the determinant, which flips its sign — exactly Pauli's rule (4.5.6).

2. **The exclusion principle is built in.** If two orbitals are equal ($\phi_a = \phi_b$), two columns of the determinant are equal and $\det = 0$. You cannot put two electrons (with the same spin) into the same orbital.

Slater determinants form a complete basis for the antisymmetric Hilbert space: the most general $N$-electron wavefunction can be written as a linear combination of all possible Slater determinants built from any complete one-electron basis $\{\chi_i\}$. The number of such determinants for $N$ electrons in $N_b$ basis orbitals is $\binom{N_b}{N}$ (per spin channel), which is still enormous — but it is the *correct* counting of the antisymmetric subspace.

This expansion is called **full configuration interaction** (FCI). With $N = 20$ electrons and $N_b = 100$ basis orbitals, the number of determinants is $\binom{100}{20} \approx 5\times 10^{20}$. FCI scales factorially with system size and is feasible only for molecules with under ~20 electrons in modest basis sets. For materials, FCI is hopeless.

!!! note "What a Slater determinant *cannot* do"
    A single Slater determinant treats electrons as moving in some average potential created by the others (Hartree–Fock, §4.7). It does not capture *correlation*: the fact that an electron at $\mathbf r$ actively *repels* nearby electrons because of $1/|\mathbf r - \mathbf r'|$, regardless of spin. The energy gap between the best single-determinant solution and the true ground state is called the **correlation energy** and is typically 1% of the total energy — but 100% of the chemistry. Bond-making is correlation. Reaction barriers are correlation. Magnetism is correlation. Reproducing correlation is the central challenge of every electronic-structure method.

## 4.5.5 An order-of-magnitude budget

Let us be concrete. Imagine we want to compute the ground-state wavefunction of a single benzene molecule, $\mathrm C_6 \mathrm H_6$. It has 42 electrons (each carbon contributes 6, each hydrogen 1).

- **Real-space grid.** A modest 3D grid with $N_g = 50$ points per direction (1 Å resolution over a 5 Å box per electron, far too coarse for chemistry) gives configuration space of dimension $50^{3\cdot 42} = 50^{126} \approx 10^{214}$.

- **Atomic-orbital basis.** Standard chemistry basis sets (cc-pVDZ) use about 20 contracted basis functions per heavy atom; benzene has roughly $6\cdot 20 + 6\cdot 5 = 150$ basis functions. The number of Slater determinants for 42 electrons in 150 orbitals (per spin) is $\binom{150}{21}^2 \approx 10^{49}$. FCI in this basis is impossible.

- **CCSD(T)**, the "gold standard" of molecular quantum chemistry, scales as $\mathcal O(N_b^7)$. For benzene with $N_b = 150$, this is $150^7 \approx 1.7 \times 10^{15}$ operations — feasible on a workstation. Result: chemical accuracy (1 kcal/mol) for the benzene total energy.

- **DFT** (Chapter 5) scales as $\mathcal O(N_b^3)$. For benzene with $N_b = 150$ this is $\sim 3 \times 10^6$ operations per SCF step — runs in seconds on a laptop. Result: bond lengths correct to about 0.01 Å.

The contrast is the entire point. Solving the Schrödinger equation directly is impossible. Reducing the problem to a well-chosen mean-field theory (Hartree–Fock) or to a functional of the electron density (DFT) makes routine chemistry tractable. The trade-off is that you no longer compute the full wavefunction; you compute an effective single-particle theory and hope that the part you have thrown away — correlation — is small or systematically correctable.

## 4.5.6 Two strategies for escape

The remaining sections of the chapter introduce two complementary approximations that begin to chip away at the exponential wall.

**The Born–Oppenheimer separation** (§4.6) recognises that nuclei are 1836 times heavier than electrons. We can freeze the nuclei, solve the *electronic* Schrödinger equation at fixed $\{\mathbf R_I\}$, and treat the resulting energy as a potential for the nuclei. The wavefunction factorises:

$$\Psi(\mathbf r, \mathbf R) \approx \chi(\mathbf R)\, \psi(\mathbf r; \mathbf R).$$

This does not reduce the electronic problem — that is still $N_{\mathrm e}$-dimensional — but it removes the nuclear coordinates from the worst-scaling part of the calculation. The PES $E_{\mathrm{BO}}(\mathbf R)$ then becomes the central object of all atomistic simulation, including the classical and machine-learning approaches of Chapters 7–9.

**Hartree–Fock** (§4.7) and its descendants tackle the electronic problem itself by restricting the wavefunction to the simplest antisymmetric form — a single Slater determinant — and minimising the energy variationally. This reduces the many-body problem to a nonlinear set of self-consistent one-electron equations, scaling polynomially rather than exponentially. The missing physics (correlation) must then be reintroduced by post-Hartree–Fock methods or absorbed into an exchange-correlation functional in DFT.

The two strategies are independent and composable: every DFT calculation you will run in Chapter 6 uses Born–Oppenheimer to freeze the nuclei *and* an effective single-determinant ansatz for the electronic ground state. Together they convert an impossible $10^{300}$-dimensional eigenvalue problem into a tractable $10^3$-dimensional one, with errors that are typically a fraction of an electron-volt per atom.

!!! tip "Where this appears later"
    The exponential wall is the reason the next two chapters exist.

    - The Born–Oppenheimer separation is worked out in detail in Section 4.6,
      and its output — the potential energy surface — is the central object
      of [Chapter 7 (molecular dynamics)](../ch07-md/index.md) and the
      machine-learning potentials of [Chapter 9 (MLIPs)](../ch09-mlip/index.md).
    - The single-Slater-determinant idea becomes **Hartree–Fock** in
      Section 4.7.
    - Most important: the impossibility of storing the $3N$-dimensional
      wavefunction is *the* motivation for replacing it with the
      three-dimensional electron *density*. That conceptual leap is
      [Chapter 5 (density functional theory)](../ch05-dft/index.md), the
      workhorse of the whole book. When you read there that "the ground-state
      energy is a functional of the density $n(\mathbf r)$", remember that the
      prize is precisely escaping the $N_g^{3N}$ count established here: the
      density needs only $N_g^3$ numbers, *independent of the number of
      electrons*.

## 4.5.6a The hierarchy of electronic-structure methods

It is useful, before closing this section, to sketch the landscape of methods that have been built to evade the exponential wall. Each method is a different compromise between cost and accuracy, and each is appropriate for a different class of problem.

| Method | Scaling | Captures correlation? | Where it lives |
|---|---|---|---|
| Hartree (mean-field, no exchange) | $\mathcal O(N^3)$ | No | Historical curiosity |
| Hartree–Fock | $\mathcal O(N^4)$ | Exchange only | §4.7; basis for post-HF |
| MP2 | $\mathcal O(N^5)$ | Perturbative | Small molecules |
| CCSD | $\mathcal O(N^6)$ | High accuracy | Benchmark chemistry |
| CCSD(T) | $\mathcal O(N^7)$ | "Gold standard" | $\lesssim 30$ atoms |
| FCI | $\mathcal O(N!)$ | Exact in basis | $\lesssim 15$ electrons |
| DFT (KS) | $\mathcal O(N^3)$ | Approximate (functional) | Chapter 5; workhorse |
| QMC (DMC/VMC) | $\mathcal O(N^{3-4})$ | Exact in principle, stochastic | High accuracy on demand |
| ML potentials | $\mathcal O(N)$ | Inherited from training data | Chapter 9; large MD |

Two observations: (1) every accurate method scales at least as the cube of system size, which is why a method that scales as $N^4$ or worse is described as "expensive"; (2) the only method on the list that scales *exactly* as the exponential wall is FCI, the brute-force diagonalisation we have just argued is impossible. Every other method is some clever rearrangement that captures the dominant correlations in polynomial time.

The trade-off is what we have stressed throughout: solving the Schrödinger equation exactly is forbidden, but solving an *approximate* version of it cheaply is the entire science of electronic structure. The exponential wall is the boundary between physics and computational physics.

## 4.5.7 Take-aways

Three messages from this section should stick.

1. **The Schrödinger equation for any non-trivial material is too large to solve.** The wavefunction lives in an exponentially large space, and storing it is physically impossible even in principle.

2. **Pauli antisymmetry both helps and hurts.** It restricts the wavefunction to a smaller subspace (helpful) but forces a particular algebraic structure — Slater determinants — that complicates calculations (hurtful). Antisymmetry is also the origin of exchange, the most important purely quantum interaction in the many-electron problem.

3. **All electronic-structure methods are approximations.** There is no "exact" method one could in principle run on a sufficiently large computer; the cost of an exact method scales exponentially with system size, full stop. Density functional theory, coupled-cluster theory, quantum Monte Carlo, and the modern neural-network ansätze of variational Monte Carlo are all attempts to capture the most important physics in a polynomially-scaling representation.

The rest of this chapter, and indeed the rest of the book, is the long story of those approximations.

!!! tip "How to think about approximations"
    A useful taxonomy: an *uncontrolled* approximation is one for which we have no a priori estimate of the error and no systematic way to improve it (most density functionals fall here); a *controlled* approximation has a known error bound that shrinks as a parameter is tuned (e.g.\ FCI in a finite basis, with basis-set convergence). The most powerful methods combine the two — for instance, **embedding theories** that solve a small, strongly correlated region with a controlled method (CCSD(T) or FCI) and treat the surrounding environment with a cheaper, uncontrolled one (DFT). The art of modern electronic structure is the *composition* of approximations, each correcting the errors of another, with the system size scaling tamed by careful exploitation of locality, sparsity, and renormalisation.

The remainder of Chapter 4 lays the groundwork. In §4.6 we make the first great simplification — the Born–Oppenheimer separation — which decouples electrons from nuclei. In §4.7 we take the first serious crack at the resulting electronic problem with Hartree–Fock. Then in Chapter 5 we replace the explicit wavefunction with an electron density, the conceptual leap that opened the door to *practical* materials simulation.

!!! question "Check yourself: antisymmetry and Slater determinants"
    1. Starting only from "electrons are indistinguishable", explain why
       exchanging two of them can multiply $\Psi$ by $+1$ or $-1$ but not by,
       say, $\mathrm i$. What extra fact tells us electrons take the $-1$?
    2. Write the two-electron Slater determinant for orthonormal spin-orbitals
       $\chi_a,\chi_b$ and expand it. Verify by hand that swapping the two
       electrons changes its sign.
    3. Set $\chi_b = \chi_a$ in that determinant. What do you get, and which
       physical principle is this?
    4. Where does the normalisation constant $1/\sqrt{2}$ in the two-electron
       determinant come from? (One sentence.)
    5. Can two electrons ever sit at the *same point in space*? Under what
       condition on their spins?

    ??? success "Answer"
        1. Indistinguishability requires $|\Psi|^2$ to be unchanged by the
           swap, so the swap can only attach a phase $\lambda$ with
           $|\lambda|=1$. Applying the swap *twice* returns the original
           state, so $\lambda^2 = 1$, leaving only $\lambda = \pm1$ — a
           generic phase like $\mathrm i$ has $\mathrm i^2=-1\neq1$ and is
           excluded. The extra fact is the spin–statistics theorem: electrons
           have half-integer spin and are therefore fermions, which fixes
           $\lambda=-1$.
        2. $\Phi(\mathbf x_1,\mathbf x_2) = \tfrac{1}{\sqrt2}\bigl[\chi_a(\mathbf x_1)\chi_b(\mathbf x_2) - \chi_b(\mathbf x_1)\chi_a(\mathbf x_2)\bigr]$.
           Swapping $\mathbf x_1\leftrightarrow\mathbf x_2$ turns the bracket
           into $\chi_a(\mathbf x_2)\chi_b(\mathbf x_1)-\chi_b(\mathbf x_2)\chi_a(\mathbf x_1)$,
           which is the negative of the original — sign flipped, as a row swap
           of the determinant requires.
        3. You get $\Phi = 0$ (two equal columns make the determinant vanish).
           This is the Pauli exclusion principle: no two electrons in the same
           spin-orbital.
        4. The un-normalised determinant is a sum of $2! = 2$ orthonormal
           product terms, so its squared norm is $2$; dividing by
           $\sqrt{2}=\sqrt{2!}$ normalises it (general case: $1/\sqrt{N!}$).
        5. Yes — but only if their spins are *opposite*. Antisymmetry then
           lives in the spin part, so the spatial part need not vanish at
           $\mathbf r_1=\mathbf r_2$. Two *same-spin* electrons have a node
           there and cannot coincide.

    ??? note "Hint"
        For 1, recall that $\hat P^2$ is the identity (swapping twice undoes
        the swap). For 4, look at the full-derivation box on the
        $1/\sqrt{N!}$ normalisation: count how many product terms survive when
        you integrate $|\Phi|^2$.
