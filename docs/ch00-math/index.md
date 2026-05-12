# Chapter 0 — Mathematics from Scratch

[Open in Jupyter (browser)](/materials-simulation-handbook/lite/lab/index.html?path=ch00-math.ipynb){ .md-button .md-button--primary }

Computational materials science sits at the intersection of physics, chemistry, computer science, and applied mathematics. Every workflow you will encounter in this handbook — from a simple Lennard-Jones molecular-dynamics run, through density functional theory, to graph neural networks trained on the Materials Project — rests on a small core of mathematical machinery. If that machinery is unfamiliar, the rest of the book will read like a foreign language. This chapter exists so that nobody is shut out for want of a calculus refresher.

We assume only high-school algebra. By the end of the chapter you will be able to manipulate vectors and matrices with confidence, take derivatives and gradients of multivariable functions, recognise where Fourier transforms appear in solid-state physics, and reason probabilistically about thermal ensembles and Bayesian models. None of this is presented as pure mathematics; every concept is grounded in an example you will meet later — the Arrhenius rate law, Hamiltonian diagonalisation, reciprocal-space lattices, Gaussian process regression for materials discovery.

## Why start with mathematics?

There is a temptation in modern computational science to skip straight to the software. Install LAMMPS, run a tutorial, modify the input file, publish a paper. This works until it does not. The moment a simulation gives an unexpected answer — a negative phonon frequency, a non-converged self-consistent cycle, an exploding loss curve during training — the user who treats the underlying mathematics as a black box has no diagnostic tools. The user who can write down the Hamiltonian, expand it in a basis, and reason about eigenvalues has an enormous advantage.

The mathematics in this chapter is not optional rigour. It is the working vocabulary of the field. You will see the same handful of ideas — linear operators, gradients, probability distributions — appear in apparently unrelated contexts. The dispersion relation of a phonon, the band structure of a semiconductor, the principal components of a structural descriptor, and the attention weights of a transformer are all eigenvalue problems in disguise. Recognising this shared structure is what allows a practitioner to move fluidly between sub-fields.

## How to use this chapter

The chapter is divided into five thematic sections plus an exercises file. Each section is self-contained but builds on the previous one. If you are returning to mathematics after a long break, work through them in order. If you are already comfortable with, say, linear algebra, skim the first paragraphs of `02-linear-algebra.md` to check our notation, then move on.

Every section follows the same pattern. Definitions come first, in plain English. A worked derivation or two follows, with each step numbered so that you can refer back. Short Python snippets — using NumPy throughout — show how the same calculations are done numerically. We finish each section with a pointer to where the material is used later in the book.

The exercises file contains eight problems with full worked solutions. Three are conceptual warm-ups, four require a page or so of work, and one is genuinely hard. Difficulty is marked with stars: ★, ★★, ★★★. Solutions are collapsed; resist the urge to expand them before attempting the problem.

## Roadmap

The five core sections are arranged so that each provides a tool used by the next:

1. **Numbers, sets, and functions** (`01-numbers.md`). The bedrock. We establish notation for the number systems and review functions, exponentials, and logarithms, ending with the Arrhenius and Boltzmann factors you will see in Chapters 7 and 8.

2. **Linear algebra** (`02-linear-algebra.md`). Vectors, matrices, determinants, eigenvalues. This is the single most important section. Quantum mechanics (Chapter 4), DFT (Chapter 5), and graph neural networks (Chapter 10) are all built from linear-algebraic primitives.

3. **Calculus and gradients** (`03-calculus.md`). Derivatives, Taylor expansion, partial derivatives, gradients, integrals, and a first encounter with variational thinking. The gradient is the workhorse of geometry optimisation, force evaluation, and machine-learning training.

4. **Complex numbers and Fourier intuition** (`04-complex-fourier.md`). Euler's formula, the FFT, and the position–momentum duality that pervades solid-state physics. Reciprocal lattices, Bloch states, and plane-wave bases all live here.

5. **Probability and statistics** (`05-probability.md`). Random variables, the Gaussian, the Boltzmann distribution, Bayes' rule. Foundations for statistical mechanics (Chapter 8) and Bayesian optimisation (Chapter 11).

Each section ends with a one-line preview of the chapter that consumes it.

!!! note "On notation"
    We use bold lower-case letters $\mathbf{v}$ for vectors, capital letters $A$ for matrices, and a hat $\hat{H}$ for operators. Greek letters carry their conventional physical meanings: $\psi$ for wavefunctions, $\rho$ for density, $\beta = 1/(k_\mathrm{B} T)$ for inverse temperature. Spelling is British English throughout.

!!! warning "What this chapter is not"
    This is not a substitute for a proper course. We skip almost all proofs of existence and uniqueness, and we are cavalier about convergence. If you want rigour, Axler's *Linear Algebra Done Right* and Spivak's *Calculus* are the classics. Our goal is operational fluency, not formal completeness.

By the end of Chapter 0 you will have a working toolkit. By the end of Chapter 1, you will know how to drive that toolkit from Python. From Chapter 2 onwards, we put both to work on actual materials problems.
