# Chapter 3 — Atoms, Bonds, Solids

A solid is a vast collection of atoms held in place by the electromagnetic interactions of their electrons. To simulate a solid, we need to decide what to keep, what to throw away, and what notation to use. That triage is the subject of this chapter.

We begin with single atoms — what they are, what we model explicitly, and what we can hide inside an effective potential. We continue with bonds — the qualitatively different ways atoms stick together, and the consequences for properties and for the choice of simulation method. We end with crystals — how to specify a periodic arrangement of atoms with a few vectors, how to read the conventional notations of the field (Miller indices, Bravais lattices, space groups), and how reciprocal space arises as the natural Fourier-dual description.

This is the most descriptive chapter in the book. The mathematics here is light; the vocabulary is dense. By the end of the chapter you should be able to look at a crystal-structure file, read off its Bravais lattice, name the bonds dominating the cohesion, and write down the reciprocal-lattice vectors and the relevant high-symmetry points in the Brillouin zone. Every chapter from 4 onwards assumes this fluency.

## Why this chapter exists

Computational chemists and materials scientists speak in a shorthand. *We optimised the Si(100)-2×1 surface at the* $\Gamma$ *point with a $4 \times 4 \times 1$ k-mesh and a 50 Ry cutoff.* Every word in that sentence carries technical content; nothing in it is decorative. If even one of the terms is unfamiliar, papers in the field will read as gibberish. This chapter establishes the dictionary.

A reader coming from a physics background may find some of the chapter elementary. A reader coming from chemistry will find the periodic-table material familiar but the reciprocal-space material less so. A reader coming from machine learning may find both halves novel. We have tried to write the chapter so that each section is self-contained: skip what you know, slow down at what you do not.

## What you will know by the end

After this chapter you should be able to:

- explain why core electrons can usually be replaced by pseudopotentials, and what is meant by the *valence* of an element;
- read an electron configuration from the periodic table and predict, qualitatively, what bonding behaviour to expect;
- describe the five canonical bond types (metallic, covalent, ionic, van der Waals, hydrogen) and identify which simulation methods handle each well;
- distinguish a primitive cell from a conventional cell, and identify the 14 Bravais lattices in two dimensions of classification (crystal system and centering);
- read Miller indices for planes and directions and explain what they mean geometrically;
- build silicon, copper, and a rock-salt structure in ASE and inspect their lattice parameters and neighbour distances;
- write down the reciprocal-lattice vectors $\mathbf{b}_i$ from real-space lattice vectors $\mathbf{a}_j$, sketch the Brillouin zone of a simple cubic and an FCC lattice, and locate the high-symmetry points $\Gamma$, X, L, K.

## Roadmap

The chapter has four content sections and an exercise set.

1. **[The atom](01-atom.md)** — What an atom is, in our model. The orbital picture, valence versus core, electron configurations, the periodic table as a pattern of valence electrons. The language used everywhere later.

2. **[Bonding](02-bonding.md)** — The five canonical bond types: metallic, covalent, ionic, van der Waals, hydrogen. Energy scales, length scales, properties predicted, and which simulation methods capture each well or badly.

3. **[Crystals](03-crystals.md)** — Lattices, primitive and conventional cells, the 14 Bravais lattices, Miller indices, common crystal structures (FCC, BCC, HCP, diamond, rock salt, perovskite). Building silicon in ASE.

4. **[Reciprocal space](04-reciprocal.md)** — From periodicity to Fourier transforms to the reciprocal lattice. Brillouin zones, high-symmetry points, why band structures look the way they do. Electron momentum and Bloch's theorem in preview.

5. **[Exercises](exercises.md)** — Seven problems mixing concept and computation.

!!! note "Notation"
    We use $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$ for real-space lattice vectors and $\mathbf{b}_1, \mathbf{b}_2, \mathbf{b}_3$ for reciprocal-lattice vectors. Capital letters $\mathbf{R}$ and $\mathbf{G}$ denote arbitrary direct- and reciprocal-lattice vectors. Atomic positions inside a unit cell are written as fractional coordinates $\boldsymbol{\tau} = (\tau_1, \tau_2, \tau_3)$ with $0 \le \tau_i < 1$. Miller indices use round brackets $(hkl)$ for planes and square brackets $[hkl]$ for directions. The Fermi level is $E_\mathrm{F}$.

!!! tip "Pen and paper"
    The material in Section 3.3, in particular, is much easier to absorb if you sketch the lattices as you go. The exercises ask you to compute neighbour distances and reciprocal vectors by hand at least once. The computer can do them faster, but the hand-calculation builds the geometric intuition that makes later debugging possible.

## The view from Chapter 12

It is worth remembering, as we begin, why this material matters even in a foundation-model future. A graph neural network operating on a crystal structure does not need to know what *FCC* means; it only sees atom types and edges. But the people *training* such networks make decisions — graph cutoffs, message-passing rules, equivariance constraints — that are informed by the same physics this chapter discusses. The most successful MLIP architectures incorporate the same notion of locality, the same separation of bonding scales, and the same symmetry structure that classical materials science codified decades ago. Knowing the classical material is not a historical exercise; it is what lets you decide whether a foundation model has learnt something real or only memorised a benchmark.

With that motivation, we begin with the atom.
