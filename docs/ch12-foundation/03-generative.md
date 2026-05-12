# 12.3 Generative Models — MatterGen and Friends

The MLIPs of Section 12.2 solve a *forward* problem: given a structure,
predict its energy and forces. The harder and arguably more useful
question runs the other way. Given a target property — a band gap of
exactly $2.3$ eV, a non-magnetic ground state, a formation energy
within $50$ meV/atom of the convex hull — can one *generate* a crystal
structure that satisfies it?

This is the inverse design problem, and it has occupied materials
informatics for two decades. The classical approach — enumerate
candidate compositions, generate plausible structures by prototype
substitution, screen with DFT — has produced real results (the
CALYPSO and USPEX evolutionary searches; the high-throughput
oxide-perovskite scans). Its limitation is that the search space
remains, in any meaningful sense, fixed by the human-designed
prototype library. Truly novel structural motifs are by construction
inaccessible.

The generative-model approach replaces the prototype library with a
learned distribution over structures, parameterised by a neural
network. Once trained, the model samples directly from this
distribution, optionally conditioned on a target property. The dream
is a tool to which one can issue the prompt "give me a thermoelectric
with $ZT > 2$" and receive a ranked list of crystallographic CIFs.
The reality, in 2026, is more nuanced — but real progress has been
made, and this section examines the state of the art.

## The shape of the problem

A crystal is specified by

- a lattice (three vectors $\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3$
  encoded as a $3 \times 3$ matrix $L$),
- $N$ atomic positions $\mathbf{r}_i \in \mathbb{R}^3$, conventionally
  given in fractional coordinates,
- $N$ atomic species $Z_i \in \{1, 2, \ldots, 96\}$.

Three structural features make this a difficult generative target.

First, *the number of atoms is variable*. A generative model must
handle structures with two atoms in the unit cell as gracefully as
those with two hundred. This rules out the simplest fixed-dimensional
architectures (e.g., a vanilla VAE on a flattened coordinate vector).

Second, *the data live on a quotient space*. Two CIFs related by a
lattice symmetry (a rotation of the cell, a permutation of equivalent
atoms, a translation by a lattice vector) represent the same physical
crystal. A generative model that does not respect these symmetries
will spend most of its capacity learning the redundancy, and its
samples will be biased toward the choice of representation in the
training set.

Third, *the support of the distribution is sparse*. The space of
$(L, \mathbf{r}, Z)$ is enormous; the subset corresponding to
thermodynamically stable crystals is a thin manifold within it. A
generative model that puts even a small amount of probability mass on
the wrong region will return overwhelmingly unphysical samples.

The dominant architectural choice that addresses all three
constraints is the *equivariant diffusion model*, and MatterGen
(Microsoft Research, 2024) is its most-cited instantiation.

## Diffusion models, briefly

Diffusion models, introduced by Sohl-Dickstein in 2015 and made
popular by Ho et al. in 2020 (DDPM) and Song et al. (score-based
generation), are now the dominant approach to generative modelling of
continuous data. The idea is to define a *forward* noising process
that progressively destroys structure, and a *reverse* denoising
process — parameterised by a neural network — that reconstructs it.

For a continuous data sample $\mathbf{x}_0$ the forward process is
typically Gaussian:
$$
\mathbf{x}_t = \sqrt{\bar{\alpha}_t}\,\mathbf{x}_0 +
\sqrt{1-\bar{\alpha}_t}\,\boldsymbol{\epsilon},
\qquad
\boldsymbol{\epsilon} \sim \mathcal{N}(0, I),
$$
where the noise schedule $\bar{\alpha}_t$ decreases from $1$ at
$t=0$ to (almost) $0$ at $t=T$. The reverse process is parameterised
by a network $\boldsymbol{\epsilon}_\theta(\mathbf{x}_t, t)$ trained
to predict the noise, equivalently the score
$\nabla_{\mathbf{x}_t} \log p_t(\mathbf{x}_t)$. Sampling proceeds by
starting from pure noise and iteratively denoising.

For crystals the picture must be modified in three ways.

**The lattice $L$ is a $3 \times 3$ matrix** that must remain
positive-definite (no zero or negative volumes) and is naturally
modelled in log-space or via a Cholesky parameterisation. Some
implementations work in the more convenient *symmetric* representation
$LL^\top$.

**The fractional coordinates** $\mathbf{r}_i \in [0, 1)^3$ live on a
torus. Adding Gaussian noise wraps around the unit cell; the
appropriate forward process is a *wrapped* normal or, equivalently, a
Brownian motion on $\mathbb{T}^3$.

**The atomic species $Z_i$** are discrete. The corresponding
generative process is a discrete-state Markov chain — the categorical
diffusion of Austin et al. (2021), or the masked-language-model
analogue used in MatterGen.

The full forward process is a product of these three. The denoising
network must learn to invert all three simultaneously, taking a noisy
$(L_t, \mathbf{r}_t, Z_t)$ and predicting the clean
$(L_0, \mathbf{r}_0, Z_0)$.

## Equivariance, again

The denoising network must respect the symmetries of the data. If
$\mathbf{x}_0$ is a physical crystal and $g$ is a global rotation, the
generative model should be equally likely to produce $\mathbf{x}_0$ or
$g \cdot \mathbf{x}_0$. The clean way to enforce this is to make the
denoising network *equivariant* in the sense of Chapter 9:
transforming the input by $g$ transforms the output by $g$ in the
appropriate representation.

MatterGen's denoiser is a stack of equivariant message-passing layers
of the kind we built in Chapters 9 and 10. The novelty is that it
processes a noisy crystal — atoms with uncertain positions and even
uncertain identities — rather than a clean one. The architecture is
nearly identical to MACE or NequIP; the differences are in the input
features (which now include the time step $t$) and the output
(predictions of the noise rather than the energy).

## MatterGen architecture, sketched

A reader familiar with Chapter 9 can construct a faithful picture from
the following ingredients.

1. **Inputs.** A noisy crystal $(L_t, \mathbf{r}_t, Z_t)$ together with
   the time step $t \in [0, T]$ and an optional conditioning vector
   $\mathbf{c}$ representing the target property.

2. **Element embedding.** Each (possibly noisy) atomic species
   $Z_i$ is mapped to a learnable embedding vector. For the
   *categorical* forward process, the embedding is in fact a soft
   mixture over elements weighted by the current denoising
   distribution.

3. **Time embedding.** The scalar $t$ is passed through a sinusoidal
   embedding and an MLP, producing a feature vector that conditions
   every node and edge representation in the network. This is the
   standard trick from DDPM literature.

4. **Graph construction.** Atoms within a cutoff (typically $5$ Å,
   evaluated with respect to the *current* lattice $L_t$) become
   neighbours. The graph topology changes with the noise level: at
   high $t$, the lattice is so distorted that the cutoff captures
   nearly everything; at low $t$, the graph stabilises to that of the
   clean crystal.

5. **Message passing.** Five or six layers of equivariant message
   passing, in the manner of MACE, update node and edge features.

6. **Heads.** Three readout heads produce, respectively,

   - $\boldsymbol{\epsilon}_L$ — the lattice noise, in the parameterisation
     of choice (a $6$-vector for the symmetric representation),
   - $\boldsymbol{\epsilon}_{\mathbf{r}}$ — the per-atom positional
     noise, an equivariant vector,
   - $\mathbf{p}_Z$ — the predicted clean-species distribution for each
     atom, a categorical over $96$ elements.

7. **Conditioning.** Property targets $\mathbf{c}$ (band gap, formation
   energy, magnetic moment, space group, chemical formula) are encoded
   into a vector that is concatenated to the time embedding, in the
   manner of classifier-free guidance.

Training is standard. For each labelled crystal in the dataset, sample
a time step $t$, apply the forward noising, and minimise a weighted
sum of mean-squared errors on $\boldsymbol{\epsilon}_L$ and
$\boldsymbol{\epsilon}_{\mathbf{r}}$, plus a cross-entropy on
$\mathbf{p}_Z$. MatterGen was trained on Alexandria (about
$6 \times 10^5$ DFT-relaxed structures, after filtering) for several
weeks on a small cluster of A100 GPUs.

## Conditional generation in practice

Sampling from a trained MatterGen, conditioned on a target property,
proceeds as follows.

```python
# mattergen_sample.py — illustrative; the public mattergen package
# exposes a higher-level API that is structurally similar.
from __future__ import annotations
from dataclasses import dataclass

import torch
from ase import Atoms


@dataclass
class GenerationTarget:
    """Property targets supplied as conditioning to the diffusion model."""
    band_gap_eV: float | None = None
    formation_energy_per_atom_eV: float | None = None
    space_group: int | None = None
    chemical_system: tuple[str, ...] | None = None
    n_atoms: int = 8


def sample_one(
    model: torch.nn.Module,
    target: GenerationTarget,
    n_steps: int = 1000,
    guidance_scale: float = 2.0,
    device: str = "cuda",
) -> Atoms:
    """Sample a single crystal from MatterGen, conditioned on `target`.

    Implements classifier-free guidance with `guidance_scale` mixing
    the conditional and unconditional score predictions.
    """
    cond = model.encode_target(target).to(device)
    null = model.encode_target(GenerationTarget(n_atoms=target.n_atoms)).to(device)

    # Initialise from pure noise.
    state = model.sample_prior(n_atoms=target.n_atoms, device=device)

    for t in reversed(range(n_steps)):
        eps_c = model.denoise(state, t, cond)
        eps_u = model.denoise(state, t, null)
        eps = eps_u + guidance_scale * (eps_c - eps_u)
        state = model.reverse_step(state, eps, t)

    return model.to_ase(state)


if __name__ == "__main__":
    from mattergen import load_pretrained  # hypothetical wrapper
    model = load_pretrained("mattergen-v1")
    target = GenerationTarget(
        band_gap_eV=2.5,
        chemical_system=("Sr", "Ti", "O"),
        n_atoms=20,
    )
    samples = [sample_one(model, target) for _ in range(50)]
    for k, atoms in enumerate(samples):
        atoms.write(f"candidate_{k:03d}.cif")
```

Two practical points.

The *guidance scale* trades sample quality for sample diversity. At
$g = 1$ one samples from the conditional distribution directly; at
$g \gg 1$ one biases the samples toward configurations that more
strongly match the target, at the cost of putting more probability
mass on a smaller set of structures. Values between $1.5$ and $3$ are
typical.

The number of denoising steps $T$ governs sample quality. Most
MatterGen variants use $T = 1000$ for inference. Distillation
techniques can reduce this to $T = 50$ at modest quality cost, but
this is not yet standard.

## The validation pipeline

A diffusion model returns *candidate* structures, not verified ones.
A serious workflow will route each candidate through several stages
of validation, exploiting the universal MLIP and DFT machinery of
earlier chapters.

The canonical pipeline is:

1. **Generate.** Sample $N \sim 10^3$–$10^4$ candidates from the
   diffusion model, conditioned on the property target.

2. **Filter for chemical sanity.** Reject candidates with
   unphysical bond lengths (any pair closer than the sum of covalent
   radii minus a tolerance), unphysical coordination numbers, or
   chemistry violating elementary valence-counting rules.

3. **Relax with a universal MLIP.** Use MACE-MP-0 (or CHGNet, or Orb)
   to relax each surviving candidate to a local minimum. Reject those
   that drift far from their initial structure, those that fragment,
   and those whose MLIP-predicted formation energy lies well above the
   convex hull. Typically $\sim 70$% of candidates are filtered out at
   this stage.

4. **Verify with DFT.** Re-relax the most promising candidates (the
   top $\sim 10^2$) with PBE-level DFT. Compute the actual formation
   energy and the target property of interest. This is the only step
   in which one obtains numbers that should be trusted.

5. **Synthesisability screen.** Apply heuristics for synthetic
   accessibility: is there a known precursor route? Is the predicted
   phase stable against decomposition into competing phases that
   appear in the Materials Project convex hull? Is it likely to form
   under reasonable temperature and pressure conditions?

6. **Experimental verification.** The brave step — sending the
   shortlist to a synthesis lab and finding out which compositions
   actually form.

A representative yield, from a 2024 MatterGen paper targeting magnetic
ground states, was roughly: $10^5$ generated candidates, $3 \times 10^4$
survived sanity filtering, $5 \times 10^3$ converged in MLIP relaxation
within $50$ meV/atom of the hull, $\sim 200$ were DFT-verified to lie
on the hull with the requested property, and a handful were attempted
experimentally. Two of those were successfully synthesised.

A 2-in-$10^5$ hit rate sounds modest until one notices that the
classical enumeration baseline, even when restricted to known
prototypes, would have had to evaluate millions of DFT calculations
to find the same hits — if the relevant structures had appeared in
the prototype library at all, which several of them did not.

### The full validation chain, quantified

The yield numbers above conceal a great deal of structure. Each stage
of the pipeline has a characteristic rejection rate, and understanding
where the candidates die clarifies both what the generative model is
good at and what it is not.

A typical breakdown — averaged across several recent MatterGen-class
studies (Zeni et al. 2025; Antunes et al. 2024; Merchant et al. 2023
GNoME) — is approximately:

| Stage | Input | Rejection rate | Survivors | Cumulative survival |
|---|---|---|---|---|
| 1. Generate | — | — | $10^5$ | $100\%$ |
| 2. Chemical sanity | $10^5$ | $\sim 25$% | $7.5 \times 10^4$ | $75\%$ |
| 3. MLIP relaxation converges | $7.5 \times 10^4$ | $\sim 30$% | $5.3 \times 10^4$ | $53\%$ |
| 4. MLIP $E_\mathrm{hull} < 50$ meV/atom | $5.3 \times 10^4$ | $\sim 90$% | $5 \times 10^3$ | $5\%$ |
| 5. DFT relaxation completes | $5 \times 10^3$ | $\sim 5$% | $4.7 \times 10^3$ | $4.7\%$ |
| 6. DFT $E_\mathrm{hull} < 0$ | $4.7 \times 10^3$ | $\sim 95$% | $250$ | $0.25\%$ |
| 7. Property target met | $250$ | $\sim 40$% | $150$ | $0.15\%$ |
| 8. Synthesisability heuristics | $150$ | $\sim 80$% | $30$ | $0.03\%$ |
| 9. Experimentally attempted | $30$ | $\sim 95$% | $1$–$3$ | $\sim 10^{-5}$ |

The headline conclusion is the one already stated: $\sim 95$–$99$% of
generated candidates are rejected before any human looks at them.
This is not a failure mode of the diffusion model; it is its
*intended* behaviour. A model that generated only on-hull structures
would be a model that had memorised the training set. A useful
generative model produces a probability distribution that is *broader*
than the support of stable structures, and the downstream filtering
is what extracts the stable subset.

### Where the rejections come from

The two largest filters are stages 4 (MLIP-hull) and 6 (DFT-hull),
which together remove roughly $99.5$% of candidates. Their failure
modes are partially independent.

**Stage 4 failures** are dominated by candidates with broadly correct
composition and topology but the wrong specific arrangement of atoms.
The MLIP relaxation pulls the structure into a local minimum that
turns out to be $200$–$500$ meV/atom above hull — a familiar enough
crystal but not a competitive one. These are the candidates that
"look right" to a chemist but are thermodynamically dominated by a
known polymorph at the same composition.

**Stage 6 failures** — candidates that pass MLIP screening but fail
DFT — split into two sub-classes. About one-third are MLIP errors:
the foundation MLIP under-predicted the energy and the structure was
not actually below hull. The remaining two-thirds are competing-phase
issues: the MLIP-predicted hull, computed using only Materials
Project compositions, missed a recent addition to the hull or a
mixed-phase decomposition that the candidate is unstable against.

The implication is that the universal MLIPs of §12.2 and the
generative models of this section are deeply coupled. Improvements
in MLIP accuracy translate directly into reduced stage-4 false
positives and stage-6 surprises. Conversely, the most expensive
generative-pipeline failures — candidates that passed all
computational filters but turned out to be wrong in DFT — are
themselves valuable training data for the next-generation MLIP.

### Common failure modes of generated structures

The failure modes at the sanity-filter stage (step 2) are stereotyped
enough to be worth cataloguing.

- **Non-physical bond lengths.** The most common failure: two atoms
  closer together than the sum of their covalent radii minus a
  tolerance. In MatterGen samples, $\sim 15$% of raw outputs have at
  least one bond shorter than $0.7$ Å, an obvious diffusion-process
  artefact when noise at high $t$ collapses two atoms together.

- **Missing atoms.** The categorical diffusion can predict an atom of
  type "padding" or "no atom" with non-zero probability; when an
  atom is dropped from the unit cell the resulting structure is no
  longer charge-balanced and is typically chemically nonsense. The
  fraction varies with the model and is suppressed by careful
  classifier-free guidance.

- **Wrong oxidation states.** A common subtle failure: the model
  generates compositions like NaO$_2$Cl$_3$ that have no
  charge-balanced oxidation-state assignment for any element. The
  Materials Project's BVAnalyzer flags these in seconds; many groups
  filter on this criterion before MLIP relaxation.

- **Implausible coordination geometries.** Octahedral cations placed
  in tetrahedral environments, or vice versa. These are not
  necessarily wrong — high-pressure phases routinely break common
  coordination heuristics — but they are usually a sign that the
  model has wandered far from its training distribution.

- **Hollow regions.** A unit cell with a large region of empty space
  that should contain atoms. Diffusion models sometimes fail to
  populate all required sites, especially at the boundary of the
  cell. Catching these requires a Voronoi-cell or pore-size analysis.

These failures are surprisingly uniform across the published
generative-model families. CDVAE, DiffCSP, MatterGen and the
CIF-language-model variants all produce them, though at different
rates and with different chemical biases.

### Recent benchmarks: MatterGen 2024 numbers

The Zeni et al. *Nature* (2025) paper on MatterGen reports the
following on its largest-scale evaluation. For unconditional
generation, $\sim 78$% of samples pass chemical sanity, $\sim 13$%
relax to within $50$ meV/atom of the convex hull as estimated by
MACE-MP-0, and $\sim 3$% are confirmed on-hull by PBE DFT. For
conditional generation targeting specific space groups, the
on-hull fraction drops to $\sim 1.5$% but the fraction matching the
requested space group reaches $\sim 65$%.

In a head-to-head with DiffCSP (Jiao et al. 2023), MatterGen produces
$1.7\times$ more on-hull structures per generation step but is
$\sim 2\times$ slower per step, yielding roughly comparable on-hull
throughput per GPU-hour. The GNoME effort (Merchant et al. 2023),
using a substitutional rather than generative approach and the
Materials Project DFT pipeline, screened $\sim 2 \times 10^7$
candidate compositions to discover $\sim 4 \times 10^5$ on-hull
structures — a higher absolute yield but at substantially higher
DFT cost per discovery. The two approaches are complementary: GNoME
extends the periodic table coverage of known-prototype chemistry;
MatterGen-class models explore genuinely novel structural motifs at
lower DFT cost per novel motif.

A 2025 follow-up benchmark by the Microsoft team reported that
fine-tuning MatterGen with property-conditioned classifier-free
guidance on $\sim 5 \times 10^4$ DFT-validated samples improved the
on-hull rate on conditional targets by roughly a factor of two over
the base model. This suggests that the generative pipeline is
approaching the regime where the same active-learning ideas of
Chapter 11 will further reduce the rejection rates — a clear
near-term research direction.

## Limitations: what the model does not know

MatterGen and its cousins are powerful but partial tools, and any
serious use requires understanding what they cannot do.

**They do not know about synthesisability.** The training data are
DFT relaxations of structures that have been computed, not synthesised.
The model learns a distribution over *computationally stable* crystals
— a much looser criterion than what an experimentalist means by
"makeable". Candidate structures that are thermodynamically reasonable
but kinetically inaccessible (high-pressure phases, metastable
polymorphs with no known low-temperature route, structures requiring
exotic precursors) appear regularly in the model's output.

**They do not understand chemistry.** The model has learnt
correlations, not principles. It will happily generate structures
with mixed oxidation states that no chemist would propose, or with
unusual coordination geometries that violate established
crystal-chemical heuristics. The valence-counting filter in the
pipeline above catches the worst of these.

**They are trained on a biased corpus.** Alexandria, like
Materials Project, is dominated by oxides and intermetallics; organic
crystals, framework materials and disordered solids are
under-represented. The model's diversity reflects this bias. Asking
for an unusual chemistry pushes the model into the tails of its
training distribution, where sample quality degrades.

**They give no uncertainty.** Each sample is presented as if equally
plausible. Two samples that look qualitatively similar may be very
different in their position relative to the training distribution.
Quantifying this — generative uncertainty, support-aware sampling — is
an active area.

**The conditioning is imperfect.** Asking for a band gap of $2.5$ eV
typically yields a distribution centred near $2.5$ eV but with
significant variance, and the band-gap predictor used during
training was itself a neural network — not DFT, and certainly not
$GW$. The actual measured band gap of a generated structure may
differ by an electron-volt or more.

## Related architectures

MatterGen has cousins, each with its own emphasis.

**CDVAE** (Xie, Fu et al., 2022). The first major equivariant
generative model for crystals. Uses a variational autoencoder rather
than a diffusion model. Simpler but less expressive.

**DiffCSP** (Jiao et al., 2023). A diffusion model targeting
*crystal structure prediction* — given a composition, generate the
ground-state structure. Narrower scope than MatterGen but more
controllable.

**Crystal-LLM** approaches (2024–2025). Tokenise CIF strings and feed
them to a transformer trained as a language model. Surprisingly
competitive on simple chemistries; less principled than equivariant
approaches but trivially conditionable on natural-language prompts.

**EDM, GeoLDM** (for molecules; 2022–2023). The molecular-graph
analogues of MatterGen. Often the architecture is identical; the
difference is in the support (molecules have no lattice).

## Connecting back to the rest of the book

The reader who has worked through Chapters 9, 10 and 11 has all the
machinery in hand. A diffusion model on crystals is

- a denoising network with the equivariant message-passing layers of
  Chapter 9,
- applied to the graph representation of Chapter 10,
- generating candidates that then enter the Bayesian / active-learning
  loop of Chapter 11.

What is genuinely new is the inversion of direction. We have moved
from a regime in which the network maps structures to properties, to
one in which it maps properties (or, more generally, partial
information) to structures. The underlying mathematics — equivariant
neural networks on point clouds — is the same.

The next section steps back from concrete tools and asks what the
remaining open problems look like, and where the field is going in
$2026$ and beyond.
