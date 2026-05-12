# Chapter 14 — Designing Your Own Project

!!! abstract "What this chapter is for"
    Everything before this chapter was about *learning techniques*. From now
    on, the centre of gravity shifts: this chapter is about *doing research*.
    By the end you will have chosen one of five concrete undergraduate-thesis
    projects (in the [`/projects`](https://github.com/) folder of this
    handbook's repository) and you will have a plan to attempt it.

You have read about density functional theory, molecular dynamics,
interatomic potentials, graph neural networks, active learning, foundation
models, and the bridges between them. That is more material than is
covered in most one-year master's programmes. You could comfortably sit a
written exam on any of it.

A written exam is not a thesis. A thesis is a single, coherent piece of
work in which you take a question, decide how to attack it with the tools
you have, and follow the attack through to a defensible answer. The
question may be small. The answer is rarely tidy. The discipline of
seeing a project from start to finish is the central skill of the working
researcher, and it is the skill this chapter is built around.

## What a capstone is, and is not

A capstone is a thesis-style project intended to integrate the technical
material from earlier in the handbook. It is *not* a literature review,
not a homework set, not a methods tutorial. It is your first attempt at
the full loop of research:

1. Pick a question that you actually care about (and that is *answerable*
   with the resources you have).
2. Find out what people have already done.
3. Decide which method, of the dozen you have learnt, is the right one.
4. Set up the calculations carefully.
5. Convince yourself the calculations are converged.
6. Run production and gather data.
7. Analyse what you have.
8. Write it up in a form that someone else could reproduce.

Most of this chapter is about steps 2-8. Step 1 is largely handled for you:
the [five project templates](07-the-five-projects.md) at the end of the
chapter give five different starting points, each calibrated to be tractable
in a typical undergraduate or master's thesis timeframe.

## Why we wrote this chapter

The reason this chapter exists at all is that *learning the methods is not
enough*. We have watched many bright students learn DFT cold, then walk
straight into a project where they:

- Spent six weeks reproducing a result that was already in the literature
  because they did not search for it.
- Reported numbers without converging k-points.
- Quoted MLIP accuracy on the *training set*.
- Spent the last week before the deadline writing the methods section from
  memory, mis-stating their cutoff energy.

None of these failures are mysterious. They are the predictable
consequences of doing research for the first time, without explicit
training in *how research works as a process*. This chapter is the
explicit training.

## Tone of this chapter, and an honest warning

We have tried to write the previous chapters in a balanced way: presenting
methods on their merits, flagging limitations, pointing to alternatives.
This chapter is different. It is mentorship, and it has opinions.

Some of those opinions are well supported across the field. Others are
defensible but contested. We have tried to mark which is which, but if
you are reading this with a supervisor who has thirty years of experience
and disagrees with us on a particular point, *they are very likely right
and we are wrong*. Use what is here as a starting framework, not a
binding rulebook.

## What is in this chapter

[Section 1 — How research works](01-how-research-works.md). The realistic
shape of a six-month thesis project: how time is spent (and mis-spent),
what to expect from group meetings and advisor check-ins, and a handful of
illustrative anecdotes.

[Section 2 — Choosing a problem](02-choosing-a-problem.md). The three
lenses (interest, feasibility, writability) for scoping a tractable
project, with concrete rules of thumb for what fits in a laptop-scale
thesis versus a cluster-scale one.

[Section 3 — Literature search](03-literature-search.md). How to find,
read, and manage the literature relevant to your project. Strategies for
the new-to-the-field reader who has been told "go and learn the literature"
with no further instructions.

[Section 4 — Convergence and validation](04-convergence-and-validation.md).
The non-negotiable starting point of any computational study. Why you must
demonstrate convergence on *your* system, before reporting any results.

[Section 5 — Common pitfalls](05-common-pitfalls.md). The recurring
mistakes that supervisors see in undergraduate and master's projects, with
specific instructions for detecting and avoiding each.

[Section 6 — Writing up](06-writing-up.md). The standard scientific
manuscript structure, and how to write a methods section reproducible
enough that someone else could rerun your work.

[Section 7 — The five project templates](07-the-five-projects.md). Brief
descriptions of the five projects in this handbook's repository, with
prerequisites, expected outcomes, and difficulty markers, plus pointers
to the detailed `projects/0X-name/README.md` for each.

The [exercises](exercises.md) at the end are scoping and critique
exercises: you will be asked to spot the problems in a poorly-written
project description, sketch a convergence study, and draft a methods
section. None of them are about computing a number; all of them are about
the *thinking that precedes the computing*.

## How to use this chapter

A reasonable approach:

- Read [Section 1](01-how-research-works.md) and
  [Section 2](02-choosing-a-problem.md) before you commit to a project.
  They will save you from premature commitment to something unworkable.
- Skim [Section 7](07-the-five-projects.md) and the project READMEs to
  pick a starting point. The five projects span DFT, MD, MLIP training,
  high-throughput screening, and active learning, so there is one in your
  area of comfort.
- Once you have started, return to [Sections 3-6](03-literature-search.md)
  as you reach each phase of the project.
- Use the exercises to check your scoping before you commit to your own
  project plan with your supervisor.

The goal is not for you to finish this chapter and know everything about
research. The goal is for you to *finish your thesis*. That requires
choosing well, planning honestly, working steadily, and writing carefully.
Each of these is harder than learning DFT. None of them is unlearnable.

Begin with [Section 1](01-how-research-works.md).
