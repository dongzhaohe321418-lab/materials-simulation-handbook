# 9.1 Why MLIPs — the accuracy–cost gap

The argument for machine learning interatomic potentials is, at root, an
argument about numbers. To see it cleanly we need to put the two
established alternatives — density functional theory and classical force
fields — side by side and ask what each can and cannot do.

## 9.1.1 The cost of DFT

A Kohn–Sham DFT calculation (Chapter 5) finds the self-consistent
electronic density of a periodic cell by repeatedly diagonalising a
Hamiltonian matrix whose dimension scales with the number of electrons.
For a plane-wave pseudopotential code with $N$ atoms in the cell, the
dominant cost is the orthogonalisation of $M$ Kohn–Sham orbitals at each
of $N_k$ k-points and each of $N_\mathrm{scf}$ self-consistency
iterations:

$$
T_\mathrm{DFT} \;\sim\; N_\mathrm{scf}\, N_k\, M^2 N_\mathrm{pw}
\;\sim\; O(N^3),
$$

where $N_\mathrm{pw}$, the size of the plane-wave basis, itself grows
linearly with $N$ for fixed energy cutoff. The prefactor is large: a
single SCF step on a 100-atom silicon cell with a moderate cutoff costs
roughly a CPU-hour on a modern node, and a converged self-consistent
calculation requires tens of such steps. A single energy-and-force
evaluation for 100 atoms therefore costs an hour or two of wall time on
sixteen cores.

This rules out almost everything one would call molecular dynamics. The
characteristic vibrational period of a stiff bond is around
$20\,\mathrm{fs}$, requiring a time step of $\Delta t \approx
0.5\,\mathrm{fs}$ to resolve it. A picosecond of dynamics is two thousand
steps; a nanosecond is two million. At an hour per step, a nanosecond of
DFT-MD on 100 atoms is two hundred years of computer time. Even the
heroic *ab initio* MD simulations published in the literature rarely
exceed $10\,\mathrm{ps}$ on $\sim\!200$ atoms, and only on the largest
supercomputers.

The scaling makes the picture worse, not better. A 1000-atom cell — the
size needed to study a single dislocation core or the solvation shell of
a polypeptide — is $10^3$ times more expensive than 100 atoms, so a
single force evaluation takes weeks. Linear-scaling DFT methods exist
(ONETEP, BigDFT, CONQUEST) but trade a worse prefactor for the better
asymptotics, and their accuracy on metallic or strongly correlated
systems remains a moving target.

The cost is not waste: it pays for chemical accuracy. A well-converged
GGA calculation reproduces bond energies to roughly
$50\,\mathrm{meV}/\mathrm{atom}$ and forces to a few
$\mathrm{meV}/\text{\AA}$, with no parameters fitted to the system at
hand. Hybrid functionals and beyond-DFT methods such as RPA push the
accuracy further at proportionally higher cost. *Whatever sits at the
end of an MLIP fitting pipeline can only be as good as the DFT data we
feed it.*

## 9.1.2 The transferability of classical force fields

A classical force field writes the potential energy as a sum of hand-coded
terms in interatomic coordinates:

$$
U = \sum_\mathrm{bonds} k_b (r-r_0)^2
  + \sum_\mathrm{angles} k_\theta (\theta-\theta_0)^2
  + \sum_\mathrm{dihedrals} V_n \cos(n\phi - \phi_0)
  + \sum_{i<j} 4\varepsilon_{ij}
    \!\left[\!\left(\tfrac{\sigma_{ij}}{r_{ij}}\right)^{\!12}
            -\!\left(\tfrac{\sigma_{ij}}{r_{ij}}\right)^{\!6}\right]
  + \sum_{i<j} \frac{q_i q_j}{4\pi\varepsilon_0 r_{ij}}.
$$

Each parameter — $k_b, r_0, \varepsilon, \sigma, q$ — is fitted once, to a
small reference dataset, and then frozen. Evaluation costs scale as $O(N)$
or $O(N\log N)$ for the long-range Coulomb part, with a prefactor
measured in nanoseconds per atom per step. A million-atom system runs in
real time on a single GPU.

The price of that speed is transferability. The functional form encodes
strong assumptions: bonds are harmonic, angles are harmonic, Lennard-Jones
captures all non-bonded interactions, charges are fixed. Move outside
the regime in which parameters were fitted — protonate a residue,
dissociate a bond, raise the temperature, change the metal — and the
predictions degrade silently. The harmonic bond term cannot break a
bond at all: stretch it past the inflection point and the force
increases without bound.

ReaxFF, the most ambitious of the reactive force fields, partially
overcomes this by replacing fixed bonds with bond-order variables that
respond to the local environment. The functional form has roughly
fifty parameters per element pair, and bond breaking is built in. In
practice, however, ReaxFF is notoriously brittle. Parameter sets fitted
for hydrocarbons fail on oxidation; parameter sets fitted for the bulk
metal fail on the surface. Cross-parameter coupling means that
re-fitting one element pair degrades others. Years of expert effort
produce parameter files that work in a small chemical neighbourhood and
fail outside it.

The structural reason is simple: any fixed functional form embeds a
prior about which interactions matter, and that prior is wrong as soon
as you visit a chemistry that was not anticipated. To make a force field
*genuinely* transferable across the periodic table you would need a
functional form rich enough to interpolate between all bonding regimes —
ionic, covalent, metallic, van der Waals — without being told in advance
which is operative.

## 9.1.3 The MLIP proposition

Machine learning interatomic potentials adopt exactly this stance. The
energy is still written as a sum of atomic contributions,

$$
U(\{\mathbf{r}_i\}) = \sum_i E_i(\{\mathbf{r}_j\}_{j \in \mathcal{N}(i)}),
$$

with each $E_i$ depending only on the atoms within a cutoff
$r_\mathrm{c}$ of atom $i$. But $E_i$ is no longer a hand-tuned
polynomial: it is a regression model — a neural network, a Gaussian
process, a sparse polynomial in a many-body basis — fitted to reproduce
DFT energies and forces on a reference dataset.

The model has enough flexibility that, given sufficient training data,
it can represent any smooth function of the local environment. The
training data is generated automatically from DFT calculations on small
configurations, and the model learns the chemistry of those
configurations without being told what bonds, angles, or hybridisations
to expect. If the training data covers the configurations the model
will encounter in production, the predictions match DFT to within the
training error — typically $1$–$50\,\mathrm{meV}/\mathrm{atom}$ for
energies, $20$–$100\,\mathrm{meV}/\text{\AA}$ for forces.

At inference time the model is just a function of atomic coordinates,
evaluated by a fixed number of matrix multiplications per atom. The
cost per atom is constant in $N$, and the absolute number is
small. A representative benchmark: MACE on a single NVIDIA A100 GPU
evaluates energy and forces for a $1000$-atom configuration in roughly
$10\,\mathrm{ms}$. That is $10\,\mu\mathrm{s}$ per atom per step,
versus roughly $30\,\mathrm{s}$ per atom per step for DFT — a
$3\times 10^6$ speedup. At $10\,\mathrm{ms}$ per step a researcher
running continuously fits eight million steps a day, or four
nanoseconds of MD on a 1000-atom cell with a $0.5\,\mathrm{fs}$ time
step. That brings nucleation, defect kinetics, electrochemistry,
biomolecular folding — every problem in the picosecond-to-microsecond
regime — within reach at DFT accuracy.

## 9.1.4 A potted history

The idea of fitting interatomic potentials to *ab initio* data is older
than machine learning. Empirical potentials such as Stillinger–Weber
(1985) and Tersoff (1988) were fitted to small DFT or
experimental datasets; embedded-atom-method potentials for metals
(Daw and Baskes, 1984) likewise. What was missing was a generic,
high-dimensional fitting strategy that did not require a chemist to
guess the functional form. That changed in the late 2000s.

**Behler and Parrinello (2007).** Behler and Parrinello proposed the
first modern MLIP. The local environment of each atom is encoded by a
set of *symmetry functions* — sums over neighbours of radial and angular
terms designed to be invariant under rotation and permutation — and the
resulting fixed-length descriptor is fed to a small feed-forward neural
network, one network per chemical element. Training is by stochastic
gradient descent on energy and force errors. The architecture, now
called a Behler–Parrinello neural network (BPNN), demonstrated for the
first time that a neural network could reproduce DFT energies of bulk
silicon over a broad range of densities and phases.

**Gaussian Approximation Potentials (Bartók et al., 2010).** The GAP
framework replaces the neural network with a Gaussian process, regresses
on the SOAP descriptor (introduced in the same paper), and obtains a
fully Bayesian potential with built-in uncertainty estimates. GAP
remains the gold standard for data efficiency on small datasets and
has produced production potentials for tungsten, silicon, amorphous
carbon, and many others.

**SchNet (Schütt et al., 2017).** SchNet was the first widely
adopted *deep* MLIP. Instead of hand-designed symmetry functions, the
network learns a representation directly from pairwise distances via
continuous-filter convolutions. The architecture is still invariant —
all features are scalars — but it generalised the
descriptor-then-regress pattern into a single end-to-end trainable
graph network. SchNet inspired a generation of descendants:
DimeNet, PaiNN, SpookyNet.

**NequIP (Batzner et al., 2022).** NequIP introduced
*equivariant* features to MLIPs. Rather than throwing away
directional information by reducing to scalar descriptors at every
layer, NequIP propagates tensors that transform predictably under
rotation — irreducible representations of $\mathrm{O}(3)$. The
empirical effect is dramatic: on the rMD17 benchmark NequIP reaches the
same accuracy as SchNet with roughly twenty times less training data,
and extrapolates substantially better. Equivariance is now the
dominant paradigm.

**MACE (Batatia et al., 2023).** MACE combines NequIP's equivariant
message passing with the body-order expansion of the Atomic Cluster
Expansion. A single MACE layer captures interactions up to fourth body
order, and two layers suffice for most chemistries. MACE is currently
the most accurate publicly available architecture per parameter and
per training example, and inference is fast enough that million-step
trajectories on thousand-atom cells are routine.

**MACE-MP-0 (Batatia et al., 2023).** MACE-MP-0 is the first
*foundation model* for materials chemistry: a single MACE potential
trained on the Materials Project dataset (roughly $1.5 \times 10^6$
DFT calculations spanning $89$ elements) and released as a checkpoint
that one can apply zero-shot to any chemistry. Out-of-the-box accuracy
is comparable to a well-fitted bespoke potential for many systems, and
the checkpoint can be fine-tuned on a few hundred new configurations
to match the bespoke potential exactly. Foundation models are the
subject of Chapter 12.

## 9.1.5 What MLIPs are not

It is worth stating clearly what MLIPs do *not* solve. They are
interpolators: they reproduce the chemistry contained in their
training set, with smoothly varying confidence in nearby regions.
They extrapolate badly. A potential trained on equilibrium liquid
water will not, in general, describe water under shock compression.
A potential trained on stoichiometric crystals will not, in general,
describe defective ones unless the relevant defects are in the
training data. Training-set design — and its automation via active
learning, Chapter 11 — is the central engineering problem of the
field.

MLIPs also inherit every error of the underlying electronic-structure
theory. A GGA-trained MACE potential reproduces GGA energetics, with
GGA's overbinding, GGA's underestimated band gaps, GGA's poor
description of dispersion. If you need hybrid-DFT accuracy you must
train on hybrid-DFT data, which is roughly thirty times more expensive
per reference calculation. If you need CCSD(T) accuracy you must
train on CCSD(T) data, which is roughly $10^4$ times more expensive.

The MLIP does not absolve you of choosing the right theory; it
amortises the cost of running it.

## 9.1.6 What we will build

Through the rest of this chapter we will build the components needed
to construct an MLIP from scratch and to use a state-of-the-art
implementation in production. The pedagogical sequence is:

1. State the symmetries any potential must satisfy (§9.2).
2. Encode them in a descriptor (§9.3) — Behler–Parrinello, SOAP, ACE.
3. Plug the descriptor into a regression model (§9.4) — BPNN, GAP.
4. Generalise the descriptor and the regressor jointly into an
   equivariant network (§9.5) — NequIP, MACE.
5. Train one (§9.6) and run MD with it.

The mathematical content is mostly linear algebra and a small amount of
group theory; the coding content is plain PyTorch. By the end you will
have trained a usable potential and will be in a position to read the
current literature.
