# 2.2 What Materials Simulation Can and Cannot Do (in 2026)

The field of computational materials science has a public relations problem. Press releases routinely claim that the next breakthrough material has been *discovered in silico*, that artificial intelligence has *predicted* the structure of a new superconductor, or that an algorithm has *designed* a battery. Some of these claims are true. Most are exaggerations. A new student who joins a computational group with expectations set by press releases is going to be disappointed within a month.

This section is the disappointment-prevention service. We will look at what materials simulation can reliably deliver in 2026, what it can deliver only with caveats, and what remains genuinely out of reach. The goal is to leave you with a calibrated sense of when a simulation result is trustworthy and when it is a hopeful guess dressed up in good plots.

## A taxonomy of claims

Before getting into specifics, it helps to distinguish three kinds of statement a simulation might make:

1. **Quantitative predictions for known quantities.** *The lattice constant of silicon is $5.43$ Å.* Simulation can be benchmarked against experiment, errors can be quantified, and a typical error bar can be reported.

2. **Quantitative predictions for unknown quantities.** *The yet-unsynthesised compound X has a band gap of $1.2$ eV.* Here we have no experimental ground truth. The trustworthiness of the prediction depends entirely on benchmarking on similar systems where ground truth is available.

3. **Qualitative or screening claims.** *Compound X is more promising than compound Y for application Z.* These are weaker but often what industrial users actually need. Ranking is easier than absolute prediction.

Different methods have different strengths in each category, and being clear about which kind of claim is being made helps to avoid disappointment.

## What we can do well

### Predict the equilibrium structure of a known material

For an ordered crystalline solid whose composition and approximate structure are known, modern DFT relaxations are robust and accurate. Lattice constants typically agree with experiment to within $1$–$2\%$ using a good GGA functional; meta-GGAs such as SCAN often do better than $1\%$. Phonon spectra computed via DFPT match inelastic neutron scattering data well in regions where anharmonicity is mild. This is the bread and butter of the field, and you should expect a relaxation to *work*.

### Compute mechanical properties

Elastic constants, bulk moduli, and shear moduli are extracted from energy-versus-strain curves and agree with experiment to roughly $5$–$10\%$ for most main-group and transition-metal compounds. Yield strengths, fracture toughness, and fatigue resistance are harder, because they depend on defects whose density is not set by the perfect-crystal calculation; but for ideal-crystal properties the agreement is generally good.

### Screen catalysts at scale

Heterogeneous catalysis on transition-metal surfaces is one of the most mature applications of DFT. The Sabatier principle, scaling relations, and computational descriptors (such as the d-band centre or the OH adsorption energy) allow tens of thousands of candidate surfaces to be ranked. Industrial laboratories routinely use DFT screening to focus experimental effort, and the predictions are good enough that the *ranking* is reliable even when individual absolute energies have errors of $0.1$–$0.3$ eV.

### Predict phase stability from composition (with caveats)

The convex hull of formation enthalpies — the union of stable phases at $T = 0$ — can be constructed from DFT calculations for a fixed composition space. Materials Project, OQMD, AFLOW, and Alexandria have done this for hundreds of thousands of compositions. Stability *at low temperature* is reasonably well predicted. Stability *at high temperature* requires free energies, which depend on phonons and configurational entropy; these are accessible but expensive. The major caveat is that the convex hull is constructed from a finite set of candidate structures: if the true ground-state structure was not in the search, it will not appear on the hull. Crystal-structure prediction algorithms (USPEX, CALYPSO, AIRSS, and recently diffusion models) help, but the search is fundamentally heuristic.

### Compute band gaps to $\sim 0.5$ eV with standard methods

Density functional theory at the GGA level systematically underestimates band gaps, often by a factor of two. Hybrid functionals (HSE06, PBE0) correct most of this and routinely match experimental gaps to within $0.3$–$0.5$ eV across a wide range of semiconductors and insulators. Many-body perturbation theory (GW) does better still — typically within $0.1$–$0.2$ eV — but at $\sim 100\times$ the cost. For screening applications, hybrid DFT is currently the sweet spot.

### Fit machine-learning interatomic potentials that rival DFT

The breakthrough of the past five years. Modern MLIPs — MACE, NequIP, Allegro, SevenNet, the foundation potentials such as MACE-MP-0 and ORB — reach DFT accuracy on energies and forces for many systems while running 100 to 10000 times faster than DFT. This has shifted what is feasible at the atomistic scale: million-atom MD simulations at DFT accuracy are now routine for some chemistries. The remaining challenges are extrapolation to out-of-distribution configurations and incorporation of long-range electrostatics, both areas of active research.

### Generate candidate materials with diffusion models

A development of the past two years. Generative models — equivariant diffusion models such as MatterGen, GNoME, and CrystalDiffusion — can sample from a learned distribution over crystal structures conditioned on desired properties. Several reports of synthesised materials originating from such pipelines have appeared in the literature in 2024–2025. The methodology is real but immature: the *generation* of candidates is now relatively easy; the *validation* (does the structure relax to itself? Is it thermodynamically stable? Can it be synthesised?) remains the bottleneck. Chapter 12 covers this in detail.

## What we can do with substantial caveats

### Predict melting points

Melting is a free-energy crossing between two phases. With a good potential and proper free-energy methods (interface-pinning, the Z method, two-phase coexistence) one can obtain melting points to within $50$ K for many metals. The catch is that *the potential must be good*. Classical EAM potentials for FCC metals get melting points right to $\sim 100$ K; for refractory metals and complex chemistries, errors of several hundred kelvin are common. Foundation MLIPs are starting to change this, but extensive benchmarking remains the user's responsibility.

### Predict defect formation energies and migration barriers

Point defects are a workhorse application of DFT, and formation energies are routinely reported. The caveats are finite-size effects (charged defects in a periodic cell are infamously tricky and require Makov–Payne or Freysoldt corrections), the level of theory (Hubbard $U$ for transition-metal compounds, hybrid functionals for shallow donors in semiconductors), and the choice of chemical potential reference. A defect energy from a careful DFT study is reliable to within $0.1$–$0.3$ eV; one from a sloppy calculation can be off by an electron-volt or more.

### Compute optical spectra

For weakly correlated semiconductors, GW-BSE delivers excellent optical absorption spectra and exciton binding energies. For molecules, TDDFT works well for valence excitations but fails for charge-transfer states and double excitations. For strongly correlated materials, dynamical mean-field theory is needed and is a specialist tool. Routine optical-property prediction for an arbitrary new material is *not* a solved problem.

### Predict reaction pathways

Climbing-image nudged elastic band (CI-NEB) finds saddle points along a known reaction coordinate, and dimer methods explore unknown coordinates. These work well for simple chemistries on metal surfaces. For complex multi-step mechanisms — say, the methane-to-methanol cycle on a non-trivial catalyst — the user must propose plausible intermediates, and the simulation can only test what is proposed. Automated reaction-network discovery is improving but not yet routine.

## What we cannot do

### Replace experiment for genuinely new physics

If a material exhibits behaviour governed by physics absent from your model — a topological phase transition not captured by your Hamiltonian, a quantum spin liquid, an unconventional superconductivity mechanism — your simulation will silently give the wrong answer. Simulation is interpolation within a model; extrapolation beyond it is unreliable. Experiment remains the final arbiter for new physics.

### Predict high-$T_\mathrm{c}$ superconductors

This is the textbook example of the previous point. Despite four decades of effort, we cannot, given a composition and a structure, reliably predict superconducting transition temperatures for unconventional superconductors. Phonon-mediated superconductors (the BCS regime) are tractable: Migdal–Eliashberg theory plus DFT-computed electron–phonon couplings has had real successes, including the prediction of hydrogen-rich high-pressure superconductors that were then synthesised. Cuprates, iron pnictides, and nickelates remain out of reach. Anyone who claims a routine workflow for designing a new room-temperature superconductor is overselling.

### Accurately compute strongly correlated systems

A material is strongly correlated when the kinetic energy of an electron is comparable to or smaller than its Coulomb repulsion with another electron on the same site. DFT, which is a mean-field theory in disguise, fails for such systems. Mott insulators, heavy-fermion compounds, and many transition-metal oxides are in this regime. Hubbard-$U$ corrections help for some properties; dynamical mean-field theory does better but is expensive and requires careful parameter choice. The recent application of quantum embedding methods and neural network wavefunctions is promising but not yet mainstream.

### Predict experimental synthesis yield

Even if you know that a compound is stable, can be synthesised, and has the properties you want, the actual yield of a real laboratory synthesis depends on kinetics, impurities, atmospheric exposure, choice of precursors, and the experimentalist's habits in ways that simulation does not capture. The gap between *predicted to exist* and *successfully made* remains one of the largest in the field. Several recent generative-model papers stumble here: structures are predicted, but only a small fraction can be made, and the simulation provides little guidance for choosing which to attempt.

### Predict failure of real components

Continuum FEM can predict stress and strain in an idealised component. Predicting *when and how* it will fail in service — fatigue cracks, corrosion, stress-corrosion cracking, environmental embrittlement — requires constitutive laws that are themselves uncertain, and microstructure-resolved models that are still research-grade. Aerospace and nuclear engineers rely heavily on experimental certification.

### Predict properties under conditions you have not trained on

This is the elephant in the room for machine-learning approaches. An MLIP trained on near-equilibrium configurations at $300$ K may give catastrophically wrong forces in a melt at $3000$ K. A property predictor trained on perovskites may produce nonsense for spinels. The community is increasingly disciplined about reporting test errors on held-out distributions, but published benchmark numbers usually flatter the model. Healthy scepticism is warranted.

## Where we are right now: a snapshot

The trajectory of the field over the past five years has been dominated by two themes.

**Machine-learning interatomic potentials have become a routine production tool.** Five years ago, fitting an MLIP was a research project in itself. Today, MACE-MP-0 and ORB are foundation potentials that work, out of the box, for most periodic-table chemistries to roughly DFT accuracy, on a million atoms, on a single GPU. This has changed what is conceivable at the atomistic scale.

**Foundation models for materials are emerging.** Beyond MLIPs, foundation models for property prediction (Matformer, MatterSim), generation (MatterGen, GNoME, ADiT), and even direct experiment design are appearing. They are not yet uniformly trustworthy — benchmark performance often does not transfer to new chemistries — but the rate of progress is striking. Chapter 12 traces this in detail.

The frontiers, then, are:

- **Long-range, transferable MLIPs.** Current MLIPs are local. Capturing dispersion, electrostatics, and screening across many chemistries without sacrificing transferability is an open problem.

- **Crystal structure prediction from composition alone.** Sampling the configuration space remains heuristic. Diffusion models offer a probabilistic approach but require validation.

- **Free energies and finite-temperature properties at DFT accuracy.** MLIPs make this feasible computationally; making it accurate (including quantum nuclear effects, configurational entropy, and rare-event sampling) is the next layer of effort.

- **Strongly correlated electrons.** DMFT, neural network wavefunctions, and quantum Monte Carlo are all moving forward, but no single method dominates.

- **Closing the simulation–synthesis loop.** Tools that go from prediction to a successful laboratory recipe — including suggesting precursors, conditions, and characterisation — are in their infancy. Autonomous laboratories and active-learning workflows are the closest current approach.

## How to read papers in 2026

Three small habits will save you a great deal of time:

1. **Look for an out-of-distribution test.** Does the paper benchmark the model on data drawn from a different distribution than the training data? In-distribution test errors of a few meV/atom are common and often mean very little.

2. **Compare to the right baseline.** Is the proposed method compared against the current state of the art (e.g., MACE-MP-0 for MLIPs, $r^2$SCAN for DFT), or against a weak straw man (e.g., LDA, a 2018-vintage GNN)? Many papers improve on outdated baselines.

3. **Look at the absolute error, not the relative.** A *25% improvement* on a benchmark whose absolute error is $50$ meV may matter; the same improvement on a benchmark whose error is $5$ eV is meaningless for any practical purpose.

!!! warning "On hype and humility"
    The temptation in this field, perhaps more than in any other branch of computational science, is to oversell what a method can do. Funding pressure, the dynamics of social media, and a culture of optimism push researchers — including this author — toward extravagant claims. The antidote is not pessimism but discipline: report your error bars, acknowledge where your method fails, and be especially sceptical of any result that claims to solve a problem that has resisted decades of effort. A paper that honestly reports a 70% success rate is a contribution; a paper that claims 99% but is silent on the test distribution should be read with one eyebrow raised.

## Looking ahead

The next section turns from capabilities to communication: the three diagrams every materials scientist reads, and how to interpret them. Once you can read a phase diagram, a band structure, and a radial distribution function, you can interpret the outputs of essentially any simulation in this book.
