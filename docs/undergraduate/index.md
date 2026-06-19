# Undergraduate Guide

This page is a **guide for reading the existing handbook** — it is not a new course and it does not contain the main technical content. The handbook already covers the physics, the methods and the code. What this undergraduate layer adds, sitting *on top* of the handbook, is a set of reading strategies, vocabulary help, formula help, gentle learning paths, small projects, and self-checks. Think of it as a study companion that walks beside the existing book rather than a second book that competes with it.

If you read nothing else, read this: you are allowed to go slowly, and you are not expected to understand everything on the first pass.

## Who this guide is for

This guide is for a normal undergraduate in materials science, physics, chemistry, or mechanical or chemical engineering. You have seen calculus, vectors and matrices, and some first-year physics, but you may not feel confident with them yet. You may have written a little Python and used a command line once or twice, but you would not call yourself a programmer. You probably have little or no research experience. None of that is a problem. The handbook was written to be learnable, and this layer exists to make the climb less steep.

## What you are expected to know

- **Maths**: A-level or first-year-university level. You should recognise differentiation and integration, vectors, and basic matrix operations, even if you need to look up the details.
- **Physics**: A-level or first-year level — forces, energy, waves, and a first taste of atoms.
- **Python**: enough to read a short script, run it, and change a number. You do not need to write programs from scratch.
- **Command line**: enough to open a terminal, change directory, and run a command someone gives you.

## What you are NOT expected to know

- Any research experience, or how to read a scientific paper.
- Advanced quantum mechanics. A first-year picture of wavefunctions and energy levels is plenty to begin.
- The internal machinery of density functional theory (DFT). You will meet DFT gently in [Chapter 5 (DFT)](../ch05-dft/index.md); you do not need to know how it works on day one.
- Machine learning. The later chapters on machine-learned potentials and graph neural networks build up from the start.

!!! note "It is normal to feel tired"

    Technical material is *dense*. A single paragraph can introduce a new symbol, a new acronym and a new idea all at once, and that is genuinely tiring to absorb. Feeling tired is not a sign that you are not clever enough — it is the expected cost of real learning, and every researcher has felt it. When a section starts to feel like too much, slow down deliberately:

    - **Re-read the plain-language part** at the top of the section before the formal version. The intuition is meant to carry you.
    - **Write the one-sentence summary** of the section in your own words. If you cannot, you have found exactly what to re-read.
    - **Look terms up** in the [beginner glossary](glossary-for-beginners.md), which explains words slowly. (There is also a terser Appendix C glossary in the main book for quick reference once you know a term.)
    - **Do one worked example by hand** before moving on — even a tiny one. Following the numbers yourself fixes an idea far better than reading does.

    Rest, then return. Difficult sections almost always feel easier the second time.

## The three learning layers

The handbook can be read at three depths. You do not have to read all three on your first pass.

- **Layer 1 — Undergraduate Core.** The essential ideas, the intuition, toy examples, basic by-hand calculations, and short pieces of code you can run. This is where you should spend your first read. If you understand Layer 1, you understand what the chapter is *about*.
- **Layer 2 — Advanced Notes.** More formal derivations, deeper theory, and technical caveats. This is genuinely useful, but it is safe to *skim* on a first pass and come back to later when the core idea has settled.
- **Layer 3 — Research Path.** Real workflows, modern methods, capstone-style projects, and the judgement needed to read the literature. This is where the field actually lives, but it is optional on a first read and assumes Layer 1 is comfortable.

!!! tip "Aim for Layer 1 first"

    On a first read of any chapter, your only goal is Layer 1. Treat Layers 2 and 3 as a return visit. Trying to absorb all three at once is the fastest way to feel overwhelmed.

!!! question "Before you move on"

    Before leaving a section, check that you can honestly tick each box. If you cannot tick one, that is your signal for what to revisit — not a reason to give up.

    - [ ] I can state the main idea of this section in one sentence.
    - [ ] I can say what physical system is being discussed (a single atom? a crystal? a box of moving atoms?).
    - [ ] I can define each new piece of vocabulary that appeared.
    - [ ] I can point to the key equation and say what it calculates.
    - [ ] I can explain the simplest example given.
    - [ ] I can run the code, or at least follow what each line does.
    - [ ] I can answer the "Check Yourself" questions in the chapter.

## Where to go next

Once you are comfortable with how to read the handbook, move on to the rest of this undergraduate layer:

- [Learning paths](learning-paths.md) — gentle, time-budgeted routes through the book for self-study.
- [Beginner glossary](glossary-for-beginners.md) — slower, friendlier definitions of the key terms.
- [Formula reading guide](formula-reading-guide.md) — how to read an equation symbol by symbol without panic.
- [Code reading guide](code-reading-guide.md) — how to read and run the handbook's code.
- [Undergraduate projects](undergraduate-projects.md) — small, self-contained projects to practise on.

It is also worth reading two pages from the main Start-Here section:

- [How to use this handbook](../how-to-use.md) — the author's own guidance on navigating the book.
- [Prerequisites checker](../prerequisites-checker.md) — a quick way to see which background topics to brush up on first.
