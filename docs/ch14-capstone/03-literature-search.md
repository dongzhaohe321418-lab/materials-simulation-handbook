# Finding and Reading the Literature

Most theses are tested at two times. The first is when your supervisor
reads draft 1 and asks "have you seen so-and-so's paper from 2018?" The
second is the viva or examiner's report, when an external researcher
asks "how does your work compare to so-and-so?". Both failures look the
same: a hole in your reading.

A solid literature search has three components: *finding* the right
papers, *reading* them effectively, and *managing* them with a reference
manager. This section covers each.

## Finding the right papers

For a beginning computational materials science researcher, the
following resources are the workhorses.

**Google Scholar** (scholar.google.com) is the default search engine for
academic literature. It indexes essentially all published work in
materials simulation, including theses, preprints, and conference
proceedings. Two features are essential:

- *Cited by*: every paper page shows how many papers cited it, and lets
  you list them. This is your gateway to forward-citation searching:
  given a foundational paper, what subsequent work built on it?
- *Cited references*: every paper's reference list, when the publisher
  makes it available, can be traversed in Scholar. This is
  backward-citation searching.

The combination — forward and backward citation traversal — is the
backbone of literature search. Start from a seed paper; follow citations
in both directions; build out a network of related work.

**arXiv** (arxiv.org) hosts preprints. For materials simulation, the
relevant categories are:

- `cond-mat.mtrl-sci` (condensed matter, materials science)
- `cond-mat.dis-nn` (disordered systems and neural networks — for ML
  papers)
- `physics.chem-ph` (chemical physics)
- `physics.comp-ph` (computational physics)

You can subscribe to the daily arXiv email for any category. For an
active project, a daily skim of `cond-mat.mtrl-sci` is fifteen minutes
that you should not skip.

**Specialised journals**. For materials simulation, the journals you
should know:

- *npj Computational Materials* — open access, broad scope, high-impact
  computational papers across materials.
- *Computational Materials Science* — Elsevier, mid-tier, lots of
  methods papers and applications.
- *Journal of Chemical Physics* and *J. Phys. Chem. C* — for
  electrochemistry, catalysis, surfaces.
- *Physical Review B* and *Physical Review Materials* — APS, broad
  scope, strong on methods.
- *Modelling and Simulation in Materials Science and Engineering* —
  IOP, broad multiscale focus.
- *npj 2D Materials and Applications*, *Acta Materialia*, *Nature
  Computational Science* — depending on your sub-field.

You do not need to read every issue of every journal. But you should
recognise the journal names when they appear in citations, and have a
sense of what each tends to publish.

**Specialised databases**:

- *Materials Project* (materialsproject.org) — DFT data on $10^5$
  compounds, with elastic constants, band structures, formation
  energies.
- *OQMD* and *AFLOW* — similar databases, slightly different
  conventions.
- *NOMAD Repository* — broader, more recent.
- *Open Catalyst Project* — surface adsorption and catalysis.

These are *primary data sources* as well as papers; you can cite the
database entry as well as the paper, and increasingly you should.

## How to find the "best" 10 papers when you start

When you start a project, your supervisor probably hands you three or
four papers as a starting point. From there, you need to grow this seed
set into a working understanding of the field. A practical procedure:

1. **Read the most recent review article** on your sub-topic. Reviews
   are designed exactly for new entrants. They survey the landscape,
   define the language, and cite the key prior work.

2. **From the review, identify 5-10 "key prior" papers**: the original
   methodological papers (the one that introduced DFT+U, the one that
   introduced NEB, etc.) and the most-cited application papers in your
   specific area.

3. **For each key paper, do a forward citation search**: what papers
   cited it? Skim the titles and abstracts. The ones that look
   thematically close to your project go on the reading list.

4. **For each key paper, look at its references**: which earlier work
   does it build on? Add the most relevant of these to your reading
   list.

5. **Cross-reference**: if a paper keeps appearing in the references of
   your reading list, it is almost certainly important. Read it.

6. **Search by methodology**: search Google Scholar for keywords like
   "DFT defect formation energy iron" or "MLIP MgO molecular dynamics".
   This catches papers that the citation network missed.

7. **Skim arXiv**: anything new in your area in the last 12 months may
   not yet be heavily cited but could be highly relevant.

After steps 1-7 you should have roughly 30-50 papers on the reading
list. Of these, 10-15 will be central; the rest are background or
methodological. This is the "best 10 papers" plus context.

The whole process takes one to two weeks, full-time. It is not a
distraction from the project; it *is* the project's first phase.

## Reading a paper

There is a difference between reading a paper *for content* and reading
it *for method*. Both have their place.

### Reading for content

Reading for content is what you do when you want to know what the paper
*found*. Standard order:

1. Title.
2. Abstract.
3. Figures and figure captions (in order — figures often tell the story
   more clearly than the text).
4. Conclusion.
5. Now go back and skim the introduction and discussion, looking for
   the specific claims made and the evidence offered for each.

This is fast. You can read for content at 5-10 papers per day.

### Reading for method

Reading for method is what you do when you want to *reproduce* what the
paper did. Now the relevant sections are:

1. Methods section, line by line.
2. Supplementary information / supplementary methods (often where the
   real details are).
3. Code and data availability statements; if there is a GitHub
   repository linked, clone it and look at the scripts.
4. The parts of results that depend on specific methodological choices.

Reading for method is slow. You may spend an entire afternoon on a
single paper, with the goal of being able to set up a similar
calculation yourself.

### Notes and skim depth

Some papers deserve careful notes; most do not. A heuristic:

- *Annotate / take notes* on the 5-10 papers that are most central to
  your project.
- *Skim and tag* the next 30-40 papers; remember they exist but do not
  invest in deep understanding unless they become relevant later.
- *Track only* the next 100; you should know the title and rough
  content but no more.

The first category goes into your thesis bibliography. The second is
your supporting context. The third is "in case the examiner asks".

## Reference management

Type out your bibliography by hand from a PDF and you will eventually
make a mistake. The mistake will be embarrassing when an examiner finds
it. A reference manager is non-negotiable.

The two free, open-source choices for academic users:

- **Zotero** (zotero.org). The default recommendation. Browser extension
  imports citations from publisher websites and Google Scholar with one
  click. Stores PDFs locally. Exports BibTeX for LaTeX. Has a usable
  desktop and mobile app. Free for unlimited references.

- **JabRef** (jabref.org). Pure BibTeX management; lighter weight but
  fewer integrations.

Commercial options (Mendeley, EndNote, Papers) exist; we have no strong
opinion. Pick one and use it.

### A workflow that works

1. When you find a paper worth saving, *immediately* click the Zotero
   button to import it. Do not "come back later". You will not.

2. In Zotero, tag the paper with one or two keywords relevant to your
   project. ("vacancy", "iron", "DFT", "review", "MLIP-training").

3. Group papers into a thesis-specific *collection*.

4. Export the entire collection as BibTeX from time to time, into a
   `references.bib` file in your thesis repository.

5. In your LaTeX (or Markdown, or Word) document, cite by BibTeX key.
   Most editors will autocomplete from `references.bib`.

This workflow takes ten minutes to set up and saves hours over a
six-month project. It also makes your thesis bibliography
*reproducible*: the same `references.bib` file can be regenerated from
the same Zotero collection.

### What to include in a citation

A BibTeX entry should have at minimum:

- Authors (all of them, not just the first three)
- Year
- Title
- Journal name (full, not abbreviation, unless your style guide says
  otherwise)
- Volume and issue
- Page numbers or article number
- DOI

Auto-imported entries usually have all of this. Check the first few you
import — sometimes the import gets the volume wrong, or capitalises the
title weirdly ("Density Functional Theory" vs. "density functional
theory"). Fix once; it stays fixed.

## Identifying gaps in the literature

The most ambitious goal of a literature search, for a thesis, is to
identify a *gap*: a question that the field has not yet answered and
that your project can answer.

Gaps come in several flavours.

- **Missing data**: nobody has computed the property for *your* specific
  system. (E.g. the formation energy of a particular defect in a
  particular semiconductor.)
- **Methodological gap**: the existing data was computed with a method
  that is now known to be inadequate. (E.g. published values using LDA
  for a system known to need hybrid functionals.)
- **Reconciliation gap**: published results disagree, and the
  disagreement has not been resolved. Your work can settle it.
- **Scaling gap**: properties have been computed for small systems, but
  the relevant phenomenon requires larger systems that are now
  tractable.
- **Combination gap**: properties have been computed separately, but
  not their combination, which would answer a new question.

Reviews often state gaps explicitly: "further work is needed on X", "the
case of Y remains poorly understood". These are starting points.

A more reliable way to find gaps: in the introductions of recent papers,
look for phrases like "however, the role of X has not been
investigated", "to our knowledge, no calculation of Y has been
reported". These sentences exist to *motivate the paper that follows*;
they are signposts to the gap that paper aimed to fill. If the paper
filled it, the gap is closed. If the paper merely identified the gap
and did not fill it, you have a project candidate.

A small warning: a gap that is easy to spot is also easy for others to
spot. If a paper says "this should be done", someone is probably doing
it. Check the most recent arXiv listings; check what the lab that wrote
the paper has done since. You do not want to spend six months only to
find a near-duplicate published last month.

## Tracking the literature over the project

Reading the literature is not a phase that ends in month 2. New papers
will appear throughout your project. Some will affect what you are
doing.

Two habits to develop:

- **Daily arXiv skim**: 10-15 minutes each morning skimming new
  `cond-mat.mtrl-sci` listings. You will spot the rare paper that is
  directly relevant.

- **Google Scholar alert**: set up email alerts for searches related to
  your project ("vacancy formation iron DFT", say). New papers matching
  your search arrive in your inbox.

The cost is small. The benefit, in not missing the one paper that would
have changed your project, is occasionally enormous.

## A common failure mode

A common scenario: a student does an initial literature search in week
1-2, builds a reading list of 30 papers, then *stops reading new
literature for the rest of the project*. By the time of the viva, six
months have passed and significant work has been published that they do
not know about.

The fix is the daily skim and the Scholar alerts. Five minutes a day
prevents the failure mode entirely.

!!! tip "Read a paper a day, every day"
    Once your initial literature phase is done, set yourself the
    discipline of reading *one new paper a day*, however short or
    skimmed. Over a six-month thesis, that is 130 papers. It is more
    than enough to stay current, builds your knowledge of the field
    beyond your specific project, and is a habit that will serve you for
    your career.

[Section 4](04-convergence-and-validation.md) takes us to the moment
when reading stops and calculation begins.
