# Exercises

Difficulty markers: (E) easy / conceptual, (M) medium / requires some
calculation, (H) hard / open-ended or substantial coding.

---

## Exercise 1 (E) — Sequential or concurrent?

For each of the following multiscale problems, decide whether a *sequential*
or *concurrent* coupling is more natural, and justify.

(a) Predicting the temperature-dependent thermal conductivity of a
thermoelectric material from DFT phonons.

(b) Following the propagation of a crack tip in a metal under cyclic loading,
where bond-breaking at the tip is essential.

(c) Computing the diffusion coefficient of a vacancy in copper from
NEB-computed barriers.

(d) Modelling the response of a polymer chain in a flowing solvent.

(e) Studying enzyme catalysis where a single C-H bond is broken in the active
site.

### Solution

(a) **Sequential.** Phonon dispersion is in steady state; once you have the
force constants, you do not need DFT again. The pipeline is DFT → phonon →
Boltzmann transport equation.

(b) **Concurrent.** The crack tip moves; new bonds are broken at each step.
The QM region must travel with the crack tip and the surrounding elastic
solid must respond. This is a classic case for an atomistic-continuum
concurrent method (CADD, bridging-domain).

(c) **Sequential.** All the rates are computed once from DFT; the KMC then
runs without needing DFT again. The pipeline is DFT-NEB → rate catalogue →
KMC.

(d) **Concurrent** if you care about the back-reaction of the solvent flow on
the chain conformation in detail (adaptive resolution scheme), but
**sequential** if you only need an effective drag coefficient from atomistic
simulation to put into a coarse-grained polymer model.

(e) **Concurrent (QM/MM).** The bond-breaking is intrinsically quantum and
localised; the protein environment is classical and polarising. This is the
canonical use case for QM/MM.

---

## Exercise 2 (M) — Error budget

You are predicting the thermal expansion coefficient of an alloy by a
multiscale workflow:

1. DFT to get the harmonic phonons. Quoted accuracy: 5% in mode frequencies.
2. Quasi-harmonic approximation to get $\alpha(T)$.
3. Comparison to experiment at 300 K.

Assume (somewhat simplistically) that a 5% error in phonon frequencies
translates into a 10% error in the thermal expansion coefficient under the
quasi-harmonic approximation.

(a) If the experimental thermal expansion at 300 K is $\alpha_\mathrm{exp} =
12 \times 10^{-6}$ K$^{-1}$, what is the range of values from your DFT
prediction that you should not regard as disagreement with experiment?

(b) The experiment has its own uncertainty, quoted as 2%. How does this
change the range?

(c) Your boss says: "I want to use this workflow to *discover* alloys with
unusual thermal expansion. What is the smallest deviation from a known
material that you can confidently claim?"

### Solution

(a) 10% of $12 \times 10^{-6}$ K$^{-1}$ is $1.2 \times 10^{-6}$ K$^{-1}$.
So values in $[10.8, 13.2] \times 10^{-6}$ K$^{-1}$ are within DFT error.

(b) The experimental uncertainty is $0.24 \times 10^{-6}$ K$^{-1}$. The total
uncertainty (adding in quadrature, treating them as independent Gaussian
errors) is $\sqrt{1.2^2 + 0.24^2} \times 10^{-6} \approx 1.22 \times 10^{-6}$
K$^{-1}$. The experimental error barely moves the needle; the DFT error
dominates.

(c) To claim a discovery you want the predicted deviation to be larger than
the uncertainty. With 10% DFT uncertainty, you need the new alloy to deviate
from the reference by *at least* 10% — roughly $1.2 \times 10^{-6}$ K$^{-1}$
in this example — before the claim is robust. Smaller deviations are within
the noise.

The pedagogical point: *the precision of the workflow sets the smallest
effect you can discover with it*. This is the fundamental limit of any
multiscale workflow, and it is the first thing to estimate before launching a
discovery campaign.

---

## Exercise 3 (M) — A subtractive scheme on paper

A graduate student writes the energy of their two-region system as

$$
E = E^\mathrm{QM}_\mathrm{inner} + E^\mathrm{MM}_\mathrm{outer} + E^\mathrm{MM,QM-MM}_\mathrm{cross}
$$

where the cross term is computed by the MM force field treating QM atoms as
point charges. They also note that their MM force field, applied to the
whole system, would give

$$
E^\mathrm{MM}_\mathrm{full} = E^\mathrm{MM}_\mathrm{inner} + E^\mathrm{MM}_\mathrm{outer} + E^\mathrm{MM,QM-MM}_\mathrm{cross}.
$$

(a) Subtract these two expressions. What does the difference tell you?

(b) Show that the student's additive expression is equivalent to a particular
form of the subtractive ONIOM expression. State the equivalence as an
equation.

(c) When does the equivalence break down? Hint: think about electrostatic
embedding.

### Solution

(a) Subtracting:

$$
E - E^\mathrm{MM}_\mathrm{full} = E^\mathrm{QM}_\mathrm{inner} - E^\mathrm{MM}_\mathrm{inner}
$$

So the student's $E$ equals $E^\mathrm{MM}_\mathrm{full} + (E^\mathrm{QM}_\mathrm{inner} - E^\mathrm{MM}_\mathrm{inner})$.

(b) This is exactly the subtractive ONIOM expression:

$$
E^\mathrm{ONIOM} = E^\mathrm{MM}_\mathrm{full} + E^\mathrm{QM}_\mathrm{inner} - E^\mathrm{MM}_\mathrm{inner}.
$$

So *additive with mechanical embedding* equals *subtractive ONIOM*, exactly.

(c) The equivalence breaks down when the cross-term is no longer purely
classical — i.e. when you use electrostatic embedding. Then
$E^\mathrm{QM}_\mathrm{inner}$ is computed *with* the MM point charges in the
Hamiltonian and its value depends on the MM configuration. The subtractive
expression with vacuum-$E^\mathrm{QM}_\mathrm{inner}$ no longer captures this
polarisation. Electrostatic embedding is intrinsically *additive*.

---

## Exercise 4 (M) — Coarse-graining check

You have run IBI to fit a CG pair potential for water-water interactions in a
mapping of 4 H$_2$O molecules per bead. The fitted potential reproduces the
target $g(r)$ within 1% pointwise.

(a) Name three observables that you should compute as a *test* of the CG
model that are not the pair RDF.

(b) For each, what would success and failure look like?

(c) Suppose the CG model reproduces the RDF and the density correctly but
gives a diffusion coefficient three times too high. What does that tell
you?

### Solution

(a) Three reasonable tests: (i) the angular distribution between
nearest-neighbour bead triplets, (ii) the bulk modulus (pressure response to
small volume change), (iii) the diffusion coefficient measured from MSD.

(b) Success means each is within ~10% of the atomistic value at the same
state point. Failure means a substantial deviation in at least one. The
angular distribution probes three-body correlations, which IBI does not fit
directly. The bulk modulus probes the equation of state. The diffusion
coefficient probes dynamics.

(c) That you have a CG model with correct thermodynamics but wrong kinetics
— exactly as expected, because IBI is fit to a static distribution function
and the projection onto CG coordinates speeds up the dynamics by removing
fast modes. The factor of 3 is in line with typical CG dynamics. You either
accept the time rescaling (and quote results in CG time units), or fit a
memory kernel to recover atomistic kinetics.

---

## Exercise 5 (H) — KMC implementation

Extend the KMC code in [Section 4](04-kmc-phase-field.md) to handle two
*inequivalent* event types: a vacancy in the bulk (barrier $E_a = 0.7$ eV)
and a vacancy adjacent to a fixed solute (barrier $E_a = 0.9$ eV when
hopping *away from* the solute, $0.5$ eV when hopping *towards* it).

(a) Sketch (in code or pseudocode) the modified BKL algorithm. What new data
structure do you need?

(b) Run your modified code at $T = 600$ K with one solute fixed at the
origin on a $50 \times 50$ lattice. Plot the steady-state probability of
finding the vacancy at each lattice site as a function of distance from the
solute.

(c) Predict, before running, the expected behaviour: should the vacancy be
attracted to or repelled by the solute? Justify with a Boltzmann argument.

(d) Why does the *initial condition* matter less here than it would in MD?

### Solution

(a) The modification: at each step, *enumerate all four possible hops from
the current vacancy position*, compute their rates given the current
configuration (which depends on the vacancy's position relative to the
solute), build a cumulative-rate array, draw $\xi_2$, and select the hop by
binary search. The new data structure is a per-step array of rates, or
equivalently a small dictionary mapping local environment to rate.

In pseudocode:

```
while running:
    rates = []
    for direction in [+x, -x, +y, -y]:
        new_pos = pos + direction
        E_a = barrier_for_environment(pos, new_pos, solute_pos)
        rates.append(nu * exp(-E_a / kT))
    R = sum(rates)
    dt = -log(rand()) / R
    j = sample_from(cumulative(rates) / R)
    pos = pos + moves[j]
    t += dt
```

(b) The expected plot: vacancy probability is enhanced near the solute and
decays to a uniform background at large distance. The enhancement factor at
the nearest-neighbour shell is roughly $\exp((E_a^\mathrm{away} -
E_a^\mathrm{toward}) / k_B T) = \exp(0.4 / (k_B \cdot 600)) \approx
\exp(7.7) \approx 2200$, which is enormous: the vacancy is strongly bound
to the solute. In practice, with a fixed solute and good statistics, the
nearest-neighbour shell will be visited far more often than expected by
chance alone.

(c) The Boltzmann argument: in equilibrium, the probability of the vacancy
being at a site of binding energy $E_b$ relative to the bulk is
$\exp(-E_b / k_B T)$. The detailed-balance condition on the rates implies
$E_b = E_a^\mathrm{toward} - E_a^\mathrm{away} = 0.5 - 0.9 = -0.4$ eV
(i.e. the bound state is 0.4 eV more stable than the unbound state). So
strong attraction.

(d) Because KMC samples a stochastic trajectory but the *time-averaged
behaviour* converges to the equilibrium distribution for long runs. The
initial position determines transient behaviour but not the steady state.
In MD, similarly, the initial conditions are washed out by sufficient
equilibration. The KMC advantage is that you can simulate the equilibrium
distribution directly with no equilibration time, because each step *is* an
equilibrium event.

---

## Exercise 6 (H) — Designing a multiscale workflow

You are asked to predict the *fracture toughness* of a polycrystalline
nickel-base superalloy at 800 K. The available computing budget is
substantial but finite (one month of access to a moderate-sized cluster).

Design a multiscale workflow. Your answer should specify:

(a) Which methods will be used at which scales.

(b) What data flows between which methods, and in which direction.

(c) Where the largest uncertainties will arise.

(d) Where you would compare to experiment as sanity checks.

(e) One simplification you would make if the budget were halved.

### Solution

(a) A reasonable workflow:

- **DFT (VASP or QE):** elastic constants of pure Ni and key alloying
  additions (Al, Cr, Co). Surface and grain-boundary energies for low-index
  facets.
- **MLIP (MACE or similar):** trained on DFT data, used to run large-scale
  MD of crack-tip propagation at the atomistic level and to compute
  temperature-dependent elastic constants by stress fluctuations.
- **Phase field:** simulate the polycrystalline microstructure of the
  alloy under representative thermomechanical treatment. Output:
  representative grain size distribution and orientation distribution.
- **Crystal-plasticity FEM:** with single-crystal stiffness from DFT/MLIP,
  slip-system parameters informed by atomistic dislocation simulations,
  and microstructure from the phase-field step. Output: macroscopic
  stress-strain curve and toughness from a notched-specimen calculation.

(b) Data flow (all sequential):

- DFT → MLIP training data
- MLIP-MD → temperature-dependent elastic constants and dislocation
  parameters
- DFT/MLIP → phase-field free-energy and interfacial parameters
- Phase-field → microstructure for CPFEM
- CPFEM → toughness

(c) The largest uncertainties: (i) the MLIP's ability to extrapolate to
crack-tip stress states it was not trained on; (ii) the phase-field free
energy at high temperature (often calibrated against CALPHAD, with its own
uncertainties); (iii) the slip-system hardening parameters in the CPFEM,
which are notoriously hard to extract from atomistic simulations alone.

(d) Sanity checks against experiment: (i) elastic constants of pure Ni at
800 K vs. published data; (ii) recovered lattice parameter and density of
the alloy vs. measurement; (iii) crack-tip stress intensity factor at
fracture under a specific load vs. experimental fracture toughness for a
similar alloy; (iv) yield stress vs. tensile-test data.

(e) Halved-budget simplification: drop the explicit phase-field step and
use a synthetic Voronoi microstructure with grain sizes calibrated from
optical micrographs. Use isotropic continuum plasticity (J2 with hardening)
informed by MLIP-MD on representative polycrystalline samples, rather than
full crystal plasticity. This sacrifices microstructure-sensitivity but
roughly halves the cost and removes the most uncertain element of the
chain (the phase-field calibration).

Marking scheme suggestion: full marks if the answer identifies a
defensible chain of methods, recognises that the MLIP and the slip-system
parameters are the high-uncertainty links, and proposes a sensible
budget-reduction strategy.
