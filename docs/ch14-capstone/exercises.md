# Exercises

These exercises are not computational. They ask you to think about
scoping, critique, and writing — the meta-skills of the chapter. Each
has a worked solution.

---

## Exercise 1 (E) — Critique a poorly-scoped project description

A student writes the following project description:

> *I will use deep learning to discover new battery materials. I will
> train a neural network on the Materials Project database and use it
> to screen for high-capacity cathodes. I will then validate the top
> candidates with DFT. The expected outcome is several novel
> high-performing cathode candidates.*

Identify at least five problems with this scoping. For each, suggest a
specific fix.

### Solution

1. **No specific question.** "Discover new battery materials" is not a
   question. What property? What chemistry (Li, Na, Mg, multivalent)?
   What target value of capacity? Without specifics, success cannot be
   judged.
   **Fix.** "Identify candidate Na-ion cathode materials with a
   theoretical capacity above 250 mAh/g, using PBE-DFT for validation."

2. **No method specifics.** "Deep learning" and "neural network" are not
   methods. Which architecture? GNN? CGCNN? MEGNet? ALIGNN?
   **Fix.** "Use the pre-trained MEGNet model from the Materials Graph
   Library and fine-tune on Na-containing oxides from MP."

3. **No success criterion.** What counts as a "high-performing cathode
   candidate"? How many candidates are enough? What is the failure mode?
   **Fix.** "Success is defined as identifying three candidates with
   DFT-validated capacities exceeding 250 mAh/g and average voltages
   between 2.5 and 4.5 V."

4. **No comparison baseline.** Without comparison to a random or
   simpler-model baseline, the ML pipeline cannot be assessed.
   **Fix.** "We will compare GNN-based screening to a composition-only
   filter (Roost) and random sampling, as baselines."

5. **No feasibility analysis.** How many DFT validations? How many
   GNN-screened candidates? How much compute is needed?
   **Fix.** "We will GNN-screen ~$10^5$ Na compositions, take the top
   100, and DFT-validate the top 10, at an estimated DFT cost of 5
   compounds × 50 CPU-hours = 250 CPU-hours."

6. **(Bonus) No risk register.** What if the GNN's predictions are
   poorly calibrated? What if all top candidates turn out to be already
   known materials?
   **Fix.** Explicit risk register with mitigations.

---

## Exercise 2 (M) — Sketch a convergence study

You are about to compute the formation energy of an interstitial
hydrogen defect in BCC iron. Design a convergence study covering all
relevant parameters. For each parameter, state:

- What values you will test.
- Which property you will monitor for convergence.
- The convergence target (in meV or similar).
- An order-of-magnitude estimate of the compute cost.

### Solution

A reasonable convergence-study plan:

| Parameter | Values to test | Monitor | Target | Cost |
|---|---|---|---|---|
| Plane-wave cutoff $E_\mathrm{cut}$ | 300, 400, 500, 600, 700 eV (or equiv. in Ry) | Defect formation energy in a fixed 54-atom cell | ≤10 meV change | 5 calculations × small cell = ~25 CPU-hours |
| K-grid (Monkhorst-Pack) | $2^3$, $4^3$, $6^3$, $8^3$, $12^3$ for 54-atom cell | Formation energy | ≤10 meV | 5 calculations = ~30 CPU-hours |
| Smearing scheme | Gaussian 0.05 eV, MP-1 0.05 eV, MP-1 0.1 eV | Formation energy | ≤10 meV | 3 calculations = ~15 CPU-hours |
| Supercell size | 16, 54, 128, 250, 432 atoms | Formation energy | ≤20 meV between 250 and 432 | 5 calculations × increasing cost = ~200 CPU-hours |
| Magnetic state | spin-polarised FM ground state vs. NM | Formation energy and total energy | NM should be higher in energy by ~0.5 eV/atom; if not, error | 2 calculations on 128-atom cell = ~40 CPU-hours |

**Order in which to run.** Cutoff and k-grid first (cheap; values
needed for everything else). Then supercell, which is the most
expensive. Smearing and magnetic state can be checked at small supercell
once the others are in hand.

**Total budget**: ~300 CPU-hours for convergence study. Acceptable.

**One subtle pitfall to flag**: the absolute formation energy depends
on the chemical potential of hydrogen, which itself is a separate
DFT calculation (H$_2$ molecule in a large box). That calculation has
its own convergence: the box size must be large enough that the H$_2$
does not interact with its periodic images. Include that as a separate
convergence test.

---

## Exercise 3 (M) — Draft a methods paragraph

You have just finished a project computing the bulk modulus of a binary
intermetallic compound by DFT. Draft a methods paragraph (200-300
words) suitable for a thesis chapter. Include:

- Software and version.
- Functional and pseudopotentials.
- Numerical parameters (k-grid, cutoff, smearing).
- The protocol for computing bulk modulus.
- Hardware.

### Solution

A sample answer:

> *All calculations were performed using Quantum ESPRESSO v7.2
> [@Giannozzi2017] with the PBE generalised-gradient approximation
> exchange-correlation functional [@Perdew1996]. Ultrasoft
> pseudopotentials from the PSlibrary v1.0.0 distribution were used
> for both elements (`Ni.pbe-spn-rrkjus_psl.0.4.1.UPF` and
> `Al.pbe-n-rrkjus_psl.0.1.UPF`).*
>
> *Plane-wave cutoffs of 60 Ry for the wavefunctions and 600 Ry for
> the charge density were employed, with a Marzari-Vanderbilt smearing
> of width 0.02 Ry. The Brillouin zone of the 16-atom conventional
> cubic supercell was sampled with a $12 \times 12 \times 12$
> Monkhorst-Pack k-grid. SCF convergence was set to $10^{-8}$ Ry, and
> force convergence for geometric relaxations to $10^{-3}$ Ry/au.
> Convergence with respect to all numerical parameters is documented
> in Appendix A; the chosen values converge the total energy per atom
> to better than 1 meV.*
>
> *The bulk modulus was computed by fitting the third-order
> Birch-Murnaghan equation of state to the total energy as a function
> of volume, with seven volume points spanning $\pm 5\%$ of the
> equilibrium volume. At each volume, the cell shape and internal
> atomic positions were fully relaxed. The fit was performed using
> the Pymatgen `EOS` module.*
>
> *Calculations were run on the example HPC cluster on Intel Xeon
> nodes with 32 cores per task. Typical wall-time per relaxation was
> 30 minutes. All input files, output files, and analysis scripts are
> archived at example DOI repository.*

Note what this paragraph does:

- Names the software, version, functional, and pseudopotentials
  explicitly.
- Gives every numerical parameter.
- Explains the protocol (Birch-Murnaghan EOS, 7 volume points).
- Mentions where the convergence study lives.
- Names the analysis library used.
- Records hardware and reproducibility metadata.

A reader could in principle reproduce this calculation from the
paragraph alone.

---

## Exercise 4 (M) — Identify the gap

Read the following abstract (composite, not real):

> *We report DFT calculations of the formation energies of intrinsic
> point defects in single-crystal cubic $\alpha$-FeSi$_2$. Vacancy,
> antisite, and interstitial defects on both Fe and Si sublattices are
> considered. PBE calculations in 96-atom supercells with $4\times4
> \times4$ k-grid sampling find that the Si vacancy is the dominant
> defect at silicon-rich conditions. Charged defect calculations are
> not included.*

(a) State two specific gaps that this paper leaves open.

(b) For each gap, sketch a thesis project that would address it.

(c) Estimate the resources each would need.

### Solution

(a) Two clear gaps:

- *Charged defects are not included.* The paper explicitly flags this,
  which is helpful. For a semiconductor or semimetal, charged defects
  often dominate at finite Fermi levels. The omission means the paper's
  conclusions about which defect is dominant may not hold in real
  doped samples.

- *Only PBE is used.* PBE has well-known limitations for semiconductor
  band gaps and for systems with localised d-states. The reported
  formation energies may be functional-dependent. Especially for charge
  transitions, a hybrid functional could reorder defect stabilities.

(b) Two thesis projects:

- *Project A: Charged-defect formation energies in $\alpha$-FeSi$_2$.*
  Re-run the same defect set with explicit charge states, including the
  Freysoldt image-charge correction. Compute formation energies as a
  function of the Fermi level, producing the standard formation-energy
  diagram. Compare to electrical characterisation in the literature.

- *Project B: Functional sensitivity of defect predictions in
  $\alpha$-FeSi$_2$.* Repeat the headline calculations with two
  additional functionals — say PBEsol and HSE06 — and assess how the
  defect ordering changes. Provide a recommended "consensus" picture.

(c) Resources:

- Project A: ~500 CPU-hours (charged defects need larger supercells
  for image-charge corrections). 8-10 weeks.
- Project B: ~1500 CPU-hours (HSE06 is expensive). 10-12 weeks. May
  need to use smaller supercells for the hybrid calculations.

Both are well-scoped thesis projects.

---

## Exercise 5 (H) — Plan a six-month thesis

Choose one of the five projects in [Section 7](07-the-five-projects.md).
Produce a *week-by-week plan* for a six-month thesis attempting it.
Include:

- Weekly goals.
- A specific deliverable (figure, table, computation, draft section)
  for each fortnight.
- Risk register with at least four items and mitigations.
- A "minimum viable thesis" — the simplest version of the result you
  would still write up if everything else fails.

### Solution

A sample answer for **Project 2: Melting Point of Copper via MLIP-Driven MD**.

**Weeks 1-2 — Setup and literature.**

- Read the original two-phase coexistence papers (e.g., Morris et al.
  for solidification methodology).
- Install LAMMPS and ASE; verify a pre-trained MLIP (MACE-MP-0 or
  CHGNet) loads and runs on copper.
- Reproduce the published EAM lattice parameter at 300 K. Deliverable:
  a plot of computed lattice parameter vs. published value.

**Weeks 3-4 — MLIP validation.**

- Compare MLIP and EAM total energies on a set of distorted copper
  structures. Quantify the MLIP's accuracy on the regime that matters.
- Run NPT MD with both potentials at 1000 K and compare lattice
  parameter and density. Deliverable: a methods-section figure with
  MLIP vs. EAM at three temperatures.

**Weeks 5-6 — Two-phase setup.**

- Construct a solid-liquid coexistence cell (one half of the cell
  pre-melted by heating to 2000 K, then re-equilibrated).
- Verify that the interface remains stable in NPT.
- Deliverable: a snapshot of the coexistence cell, plus a plot of cell
  volume vs. time at a fixed coexistence temperature.

**Weeks 7-10 — Production runs (MLIP).**

- Run coexistence MD at 1200, 1300, 1350, 1400 K. Determine whether
  the solid or liquid grows.
- Bisection toward the melting point.
- Deliverable: plot of which phase grew at each temperature; estimate
  of melting point.

**Weeks 11-13 — Production runs (EAM) and comparison.**

- Same protocol with EAM for comparison.
- Compute statistical uncertainty (block averaging over multiple
  trajectories).
- Deliverable: comparison table of MLIP $T_m$, EAM $T_m$, and
  experimental $T_m$.

**Weeks 14-17 — Analysis.**

- Density of states at the interface; identify any pathological
  behaviour.
- Heat-of-fusion calculation as a sanity check.
- Compute the order parameter (Steinhardt $Q_6$, for example) to
  confirm solid vs. liquid identification.
- Deliverable: complete results-section figures and tables.

**Weeks 18-22 — Write up.**

- Draft methods, results, discussion.
- Revisions with supervisor.
- Final figures, bibliography, abstract.

**Weeks 23-24 — Final polish and viva preparation.**

- Submit thesis.
- Practice viva responses.

**Risk register:**

| Risk | P | Impact | Mitigation |
|---|---|---|---|
| Pre-trained MLIP gives unphysical liquid Cu | M | H | Test against EAM in early weeks; if MLIP fails, do thesis with EAM only |
| Coexistence cell does not stabilise | M | M | Use longer equilibration; try larger cells |
| Cluster downtime > 1 month | L | M | Local workstation fallback for EAM; preempt by booking time in advance |
| Insufficient statistics | M | M | Block averaging from the start; running multiple trajectories |
| Melting point off by > 100 K from experiment | M | L | This *is* a possible scientific outcome, not a failure; discuss as such |

**Minimum viable thesis.** If MLIP fails entirely, the EAM-only
two-phase melting-point calculation, with full convergence and a
careful methods section, is itself a complete thesis. It reproduces a
known experimental value and demonstrates the two-phase methodology.
The MLIP comparison becomes a "future work" or a small additional
section.

**Marking note.** Full marks if the plan is realistic, has explicit
deliverables for each fortnight, has a coherent risk register with
mitigations that are sensible, and identifies a minimum viable
thesis. Half marks if any of these are missing.
