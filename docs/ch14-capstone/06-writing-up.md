# Writing Up

Writing the thesis is, for most students, the hardest part of the
project. It is harder than the calculations. It is harder than the
literature search. It is harder than the analysis. And it must be done
in less time than any of them.

We have a strong opinion on this: **writing badly is not a personality
trait, it is a skill deficit**, and the skill is learnable. This section
is the curriculum.

## The standard scientific manuscript structure

The conventional structure of a scientific paper, and by extension a
thesis chapter, is:

- **Abstract**: one paragraph; what you did, what you found, what it
  means.
- **Introduction**: the problem, why it matters, what is known, what
  you will add.
- **Methods**: what you did, in enough detail to reproduce.
- **Results**: what you found, with figures and tables.
- **Discussion**: what the results mean, in context.
- **Conclusion**: brief summary, future directions.

A thesis adds:

- **Literature review**: an extended introduction or a separate
  chapter, depending on length conventions in your department.
- **Appendix**: convergence plots, code listings, supplementary data.
- **Bibliography**: full reference list.

The order matters. Each section serves a specific role and should not
duplicate others. Concretely: the abstract is not a teaser; the
introduction is not a literature review; the discussion is not a
re-statement of results.

We will walk through each section with a focus on the choices that
matter.

## Abstract

The abstract is the most-read part of your thesis. The examiner reads
it before the rest, and many readers will read *only* the abstract. It
must stand alone.

A four-sentence skeleton:

1. The problem and why it matters. ("Predicting defect formation
   energies in alloys is essential for steel design, but few systematic
   studies exist for chromium in austenitic iron.")
2. What you did. ("We computed the formation energy of Cr substitutional
   defects in $\gamma$-Fe over a range of compositions and pressures
   using PBE-DFT with 128-atom supercells.")
3. What you found. ("The formation energy decreases monotonically with
   pressure, with a slope of $-0.04$ eV/GPa, and shows a non-monotonic
   composition dependence reflecting magnetic frustration around
   $x_\mathrm{Cr} = 0.2$.")
4. What it means. ("This identifies the composition regime in which
   defect concentration is suppressed under high pressure, with
   implications for processing of stainless steels for extreme
   environments.")

Keep it to 150-250 words. No citations. No equations unless absolutely
unavoidable. No jargon that an examiner outside your sub-specialty
cannot parse.

## Introduction

The introduction has one job: by the end of it, the reader should
understand what question you are answering and why.

A reliable structure:

- **Paragraph 1**: the broad context. ("Stainless steels are
  used in...")
- **Paragraph 2**: the specific scientific gap. ("However, the role of
  pressure on defect chemistry in these alloys has been studied
  only...")
- **Paragraph 3**: how others have approached related questions.
  ("Recent DFT work on pure iron found...")
- **Paragraph 4**: your specific question and approach. ("In this work
  we extend these calculations to...")
- **Paragraph 5**: a roadmap of the rest of the thesis. ("Section 2
  describes the methods; Section 3 presents the convergence study;
  Section 4 reports the main results; Section 5 discusses
  implications.")

The fifth paragraph is more typical of a thesis than a paper, but it
helps the reader navigate.

What to avoid:

- Vague openings ("Materials are very important in modern technology.").
- Generic reviews ("Many studies have addressed defects in metals.").
- Failing to state the question explicitly.

## Methods

The methods section is *the* reproducibility section. A well-written
methods section allows someone in another group, with no prior contact
with you, to re-run your calculations.

What to include:

- **Software and version**. "Calculations were performed using Quantum
  ESPRESSO v7.2 [@QEref]." Not "we used DFT". Not "Quantum ESPRESSO".
  The *version*.
- **Functional and pseudopotential**. "We used the PBE
  exchange-correlation functional [@PBEref] with PAW pseudopotentials
  from the PSlibrary v1.0.0 [@PSref], specifically the Fe.pbe-spn-kjpaw
  dataset."
- **Numerical parameters**. Cutoff energy, k-grid, smearing scheme and
  width, convergence criteria for SCF and forces.
- **Cell and geometry**. Cell composition, supercell size, geometry
  source (relaxed from experimental input? optimised from scratch?).
- **Protocols and analysis**. Exactly how each reported quantity was
  computed.
- **Hardware**. CPU type, number of cores, typical wall-time. Useful for
  reproducibility and for sizing future work.
- **Code and data availability**. A link to a public repository or
  archive, with a DOI.

A reasonable methods section is 1-3 pages for an undergraduate thesis.
Substantially shorter and you are probably not being reproducible.
Substantially longer and you are probably padding.

### A small methods-section example

> *All calculations were performed using Quantum ESPRESSO v7.2. We used
> the PBE exchange-correlation functional with PAW pseudopotentials
> from the PSlibrary v1.0.0 (`Fe.pbe-spn-kjpaw_psl.0.2.1.UPF` for iron,
> `Cr.pbe-spn-kjpaw_psl.0.3.0.UPF` for chromium). Plane-wave cutoffs
> were 70 Ry for the wavefunctions and 700 Ry for the density. The
> Brillouin zone was sampled with a Monkhorst-Pack $4 \times 4 \times 4$
> k-grid for the 128-atom supercell, with Marzari-Vanderbilt smearing of
> width 0.02 Ry. SCF convergence was set to $10^{-8}$ Ry; force
> convergence for relaxation to $10^{-3}$ Ry/au.*
>
> *Defect formation energies were computed as
> $E_f = E_\mathrm{defect} - E_\mathrm{bulk} - \mu_\mathrm{Cr} + \mu_\mathrm{Fe}$,
> where $\mu_\mathrm{Fe}$ is the energy per atom of bcc Fe at the same
> pressure, and $\mu_\mathrm{Cr}$ is the energy per atom of pure bcc Cr.*
>
> *Calculations were run on the Whatever Cluster at the University of
> Wherever, on Intel Xeon nodes with 32 cores. Typical wall-time per
> relaxation was 6 hours. Input files, output files, and the analysis
> scripts used to produce the figures are available at
> https://doi.org/10.5281/zenodo.1234567 (replace with your own).*

Compare to the bad version: "*We used DFT in VASP. K-grid: 4×4×4.*"
This is uninformative and unreproducible.

## Results

The results section presents what you found. It should be *almost
entirely figures and short paragraphs*, with the figures doing the
heavy lifting.

A few specific recommendations.

### Every figure should answer one question

Each figure has a single, specific question it is answering. If you
cannot state that question in one sentence, the figure is unfocused.

A figure caption should:

- State what is shown.
- State the answer the figure provides.
- Give enough technical detail (axes, units, sample size, error bars)
  that the figure stands alone.

Example: "Figure 3: Defect formation energy of Cr in $\gamma$-Fe as a
function of pressure. The energy decreases monotonically by $\sim 4$
eV over the studied range, with a near-linear slope of $-0.04$ eV/GPa.
Each point is from a single 128-atom supercell DFT calculation;
uncertainties (smaller than the markers) are from k-point convergence."

### Label everything

Axis labels with units. Legend entries with units. Tick marks at
sensible values. No "Figure 3a" that does not appear in the figure.

### Never use the heatmap rainbow palette

The rainbow palette (sometimes called "jet") is *not* perceptually
uniform: identical numerical differences look different in different
parts of the colour scale. It is also inaccessible to colour-blind
readers.

Use a perceptually uniform colormap: viridis, magma, cividis (which is
specifically colour-blind safe), plasma. For diverging data (positive
and negative around zero), use coolwarm, RdBu, or PiYG.

For categorical data (different materials, different conditions), use a
qualitative palette (matplotlib's tab10, ColorBrewer's Set2) and limit
to ~6-8 categories per figure.

### Figures, not screenshots

Produce figures from your analysis scripts and save as PDF (for vector
graphics) or PNG with $\geq 300$ dpi (for raster). Never screenshot a
plot from a Jupyter notebook into a thesis; the result is jagged and
unprofessional.

A reproducibility bonus: include a `figures.py` script in your project
repository that regenerates every figure in the thesis from the raw
data. Future-you, six months later, will need this when the examiner
asks about the y-axis scaling.

## Discussion

The discussion is where you turn results into meaning. This is the
section that most distinguishes a strong thesis from a weak one.

What goes in the discussion:

- **Interpretation of the results**. What physical picture explains
  what you found?
- **Comparison to prior work**. Where your numbers agree, why? Where
  they differ, why?
- **Limitations**. What your results do not tell you. What
  assumptions could break.
- **Implications and future directions**. What would be the next thing
  to do? What design or experimental decisions does your work inform?

A discussion section is not:

- A restatement of the results section in words.
- Vague speculation untethered to the data.
- A grand claim of universal applicability.

A useful exercise: write the discussion section *as if you were the
reviewer*. What are the obvious counter-arguments? What are the
limitations? Address them explicitly. A discussion that openly
acknowledges weaknesses is far more credible than one that pretends
they do not exist.

## Conclusion

The conclusion is short — half a page to a page — and recapitulates
the highlights. It should answer the same question as the abstract,
but with confidence that the reader has now read the supporting
arguments.

Three sentences are often enough:

1. The single most important finding.
2. Why it matters.
3. What comes next.

## References

Use a reference manager. We said this in
[Section 3](03-literature-search.md) and we mean it. Never type a
BibTeX entry by hand. Never let your reference list be a hand-typed
list in the thesis document.

For citation style, follow the convention of your department or of the
journal you are targeting. For materials simulation, common choices:

- **Numeric** (e.g. [1, 2, 3]): more compact, default in many physics
  journals.
- **Author-year** (e.g. Smith et al., 2020): more readable, common in
  reviews and theses.

Be consistent. Use one style throughout.

## Open data and code

Modern computational materials science is increasingly expected to be
*open*: input files, scripts, raw data, all available.

A minimal open-data deposit:

- A Zenodo entry (DOI-stamped) containing:
  - Input files for every calculation reported in the thesis.
  - Selected output files (especially relaxed structures, total
    energies). Full output files for headline results.
  - Analysis scripts that produce the figures.
  - A README explaining the structure of the archive.

Zenodo provides DOIs free of charge and integrates with GitHub. The
process from "I want to share my code" to "I have a DOI to cite in my
thesis" is about 15 minutes.

For the code: a GitHub repository, linked to from the Zenodo deposit,
with a clear README and a sensible licence (MIT or BSD for permissive,
GPL for copyleft).

Citing the handbook: if this textbook has been useful, please cite it.
Citation information is in the handbook's front matter.

## A writing schedule

We promised earlier in this chapter that we would treat writing as a
skill. Skills require practice. Here is a writing schedule for the
last six weeks of a six-month project:

- **Week -6**: outline the thesis. Bullet points for every section.
  Identify which figures are still needed; generate them.
- **Week -5**: draft the methods section. This is the easiest section
  to draft because you already know what you did. It is also a useful
  warm-up.
- **Week -4**: draft the introduction and the literature review. By
  this point your literature notes from earlier in the project should
  be ample.
- **Week -3**: draft the results section. Place every figure with
  caption. Write 2-3 sentences per figure.
- **Week -2**: draft the discussion and conclusion. Now you know what
  the story is, so the discussion will write itself once you commit to
  the story.
- **Week -1**: edit, polish, fix references, generate the final PDF.
  Get someone else to read it.

This schedule allows three to four weeks of *drafting* (which is
slow, generative writing) and one week of *editing* (which is faster,
selective writing). Most students get the proportion the wrong way
round.

!!! tip "Write the methods section first"
    Drafting the methods section first does two things: it gives you a
    quick win (it is the easiest section to write), and it forces you to
    notice if any of your methodological choices are not yet fully
    justified. Better to notice in week -5 than week -1.

## Getting feedback

Get someone else to read your draft. Multiple someones, if possible:

- Your supervisor (will catch scientific errors and gaps).
- A peer in your group (will catch unclear explanations).
- Someone outside materials simulation (will catch unexplained jargon).

The single most useful piece of feedback is: "*what was confusing?*"
Confusion in the reader is your problem to fix, not theirs.

Allow at least a week for feedback turnaround. Supervisors are busy.
Plan accordingly.

## A small note on writing tools

Use a tool that will not break under your thesis.

- **LaTeX**: the default for physics, chemistry, materials science. Steep
  learning curve, but produces beautiful output and is unbreakable for
  long documents with many equations and citations. Pair with a good
  editor (VS Code with LaTeX Workshop, Overleaf for cloud, TeXShop on
  Mac).
- **Markdown + Pandoc** (or Quarto): lighter weight, good for shorter
  documents, can compile to PDF via LaTeX. Increasingly popular for
  reproducible documents that mix code and text.
- **Word / Google Docs**: works for short theses, painful for long
  ones, fragile with equations and citations. Use only if your
  department mandates it.

Whatever you use, put the source under version control. A version-
controlled thesis is one where you can always recover to "yesterday's
working version" when something goes wrong at 2 am the night before
submission. (Things will go wrong at 2 am the night before submission.)

[Section 7](07-the-five-projects.md) presents the five thesis project
templates that this handbook provides.
