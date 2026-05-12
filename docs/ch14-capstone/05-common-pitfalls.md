# Common Pitfalls

We have spent the last four sections on the meta-skills of running a
thesis project. This section is more concrete. It is a catalogue of
specific, recurring mistakes that supervisors see in undergraduate and
master's projects, with diagnostics and fixes for each.

Most of these will eventually be caught — by a careful supervisor, by an
external examiner, by a reviewer. But each catch costs time and
credibility. Pre-empt them.

The list is not exhaustive. It is what we have seen go wrong most often.

## DFT pitfalls

### 1. Not setting `occupations='smearing'` for metals

In a metallic system, the Fermi level cuts through bands. Numerically,
the Brillouin-zone integration of occupation numbers around the Fermi
surface requires *smearing*: assigning fractional occupations to bands
near $E_F$, so that the integration over **k** is numerically stable.

If you forget — i.e. if you use `occupations='fixed'` or an
insulator-style integration scheme on a metal — your calculation either
fails to converge, or it converges to a wrong total energy with
characteristic spurious features (very high SCF iteration counts,
sensitivity to k-grid).

**Detection.** Sudden jumps in total energy as a function of small
perturbations (k-grid, smearing width). Failure to converge SCF in a
reasonable number of iterations. The DOS at the Fermi level being
exactly zero when it should not be.

**Fix.** Set `occupations='smearing'` (VASP/ASE) or
`occupations='smearing'` plus a smearing scheme in QE (Methfessel-Paxton
of order 1 is standard for metals; Gaussian width 0.05-0.1 eV is a
common default; Fermi-Dirac for finite-temperature work).

**Bonus mistake.** Using smearing width *much* larger than $k_B T$ at
your simulation temperature. The smearing width is a numerical
artefact; once converged, smaller is more accurate but more expensive.
A width of 0.05 eV is fine for room-temperature metallic calculations;
0.2 eV is suspicious unless you have a specific reason.

### 2. Using PBE for van der Waals-bonded systems

PBE — and GGAs in general — does not include van der Waals (dispersion)
interactions. For layered materials (graphite, MoS$_2$, h-BN), molecular
crystals, gas-phase adsorption, and any system held together by
dispersion forces, plain PBE gives the wrong interlayer spacing, the
wrong adsorption energy, sometimes even the wrong sign.

**Detection.** Layered systems with too-large interlayer distance.
Adsorption energies that are factor-of-two off literature values.
Molecular crystals that "fly apart" when relaxed.

**Fix.** Use a dispersion-corrected functional:

- **PBE-D3 (Grimme):** the workhorse. Cheap correction, well-validated.
- **PBE-D4:** an improvement on D3.
- **vdW-DF, vdW-DF2:** non-local correlation functionals; more
  principled but heavier.
- **rev-vdW-DF2, optB88-vdW:** widely used variants.

Choose one and report it.

### 3. Under-converged k-grid

The k-grid sampling of the Brillouin zone affects every quantity that
depends on the band structure (essentially every quantity in a periodic
DFT calculation). Under-converged k-grids are the single most common
silent error in undergraduate DFT projects.

**Detection.** Compute the property at several k-grid densities and
plot. If the property is still drifting at your "production" k-grid,
you are under-converged.

**Fix.** Increase the k-grid until the property plateaus. For metals
this may mean ridiculous-feeling grids ($24 \times 24 \times 24$ or
beyond for small primitive cells); for insulators with large lattice
constants, smaller grids suffice. The Monkhorst-Pack rule of thumb is
to scale the grid inversely with the lattice parameter, i.e. the
grid times lattice constant in each direction should be roughly
constant ($\sim$ 30-40 Å for a converged calculation).

### 4. Wrong space group or forgotten symmetry

Most DFT codes can symmetrise the input structure, recognise the space
group, and use the symmetry to reduce the number of independent k-points
and force components. If you skip this — e.g. by setting `ISYM=0` in
VASP or by feeding a slightly distorted structure that breaks
symmetry — you may run perfectly correct calculations that are simply
*more expensive* than necessary, or, worse, calculations that find a
spurious lower-symmetry state.

**Detection.** Run a symmetry analyser (Pymatgen's `SpacegroupAnalyzer`
or VASP's `SYMPREC` output) on your input structure and compare to the
expected space group. If they disagree, your structure is wrong.

**Fix.** Symmetrise the structure before launching the calculation.
Pymatgen's `SpacegroupAnalyzer.get_refined_structure()` is one route.
ASE's `ase.spacegroup.symmetrize` is another.

The exception: for finite-temperature MD or for studying symmetry
breaking, you *want* the symmetry off. Be explicit about which case you
are in.

### 5. Spin polarisation forgotten

Magnetic systems (iron, nickel, manganese-bearing oxides, many defects)
require spin-polarised DFT. Forgetting this gives nonsensical results.

**Detection.** Magnetic system, but calculated magnetic moment is zero.
Or: known ferromagnet computed as a non-magnet because spin polarisation
was off.

**Fix.** Turn on spin polarisation, set initial magnetic moments
sensibly (don't all start at zero — that's a local minimum), and check
that the final state is what you expect.

## MD and MLIP pitfalls

### 6. MD trajectories that haven't equilibrated

You launched a 1 ns simulation, computed your observables from the
whole trajectory, and reported them. But the system spent the first
200 ps drifting away from its initial condition. Your observables
average together the equilibration period and the actual production
period.

**Detection.** Plot the observable as a function of MD time. The
equilibration period is the early part where the observable is still
trending. Discard it.

**Fix.** Define a clear *equilibration* phase (typically 10-20% of the
total trajectory) and discard those frames from the production
averaging. For tricky systems (slow relaxation, hysteresis), longer
equilibration is needed.

### 7. Timestep too large

The MD timestep must be small enough to resolve the fastest motion in
the system. For systems with O-H or N-H bonds, this is roughly 0.5 fs
(or shorter if you do not constrain those bonds). For metals without
hydrogen, 1-2 fs is often fine.

**Detection.** Energy drift in NVE simulations; non-conservation of
total energy. Sudden jumps in temperature. The MD code refuses to
converge or atoms move ridiculous distances per step.

**Fix.** Reduce the timestep. Or constrain the fastest bonds (SHAKE,
RATTLE for H-bonds in classical force fields).

### 8. Train/test split letting polymorphs leak

When training an MLIP, you typically split your DFT data into training
and test sets. The naive split is *by random configuration*. But if
your dataset contains multiple polymorphs of the same composition, or
multiple snapshots from the same MD trajectory, random splitting puts
near-identical configurations into both train and test. The reported
test error is then optimistically biased.

**Detection.** Compute the radial distribution function (or some other
fingerprint) for each configuration; compute the nearest distance from
each test configuration to its nearest training configuration. If this
distance is very small, your test set is too easy.

**Fix.** Split by *polymorph* or by *trajectory*: hold out entire
trajectories, or entire compositions, from the training set. This
gives a much more honest measure of out-of-distribution accuracy.
Report it.

### 9. MLIP MAE reported without baselines

You trained an MLIP and report a test MAE of 5 meV/atom. Good or bad?
Without context, the reader cannot tell. Modern MLIPs on standard
benchmark datasets routinely reach 1-2 meV/atom; 5 meV/atom would be
mediocre. On harder out-of-distribution data, 50 meV/atom would be
respectable.

**Detection.** Re-read your thesis methods section. Does it report an
MAE? Does it compare to *some* baseline (a simpler model, a published
benchmark, the DFT-to-DFT noise on similar data)?

**Fix.** Always report MAE alongside (i) the MAE of a trivial baseline
(e.g. the mean energy), (ii) the MAE of a published comparable model on
the same data, or (iii) the spread of DFT energies themselves on
near-identical configurations. A meaningful MAE is one in context.

## High-throughput and ML pitfalls

### 10. Selecting the test set after fitting

You train an ML model, evaluate it on a test set, observe poor
performance, change the test set, and re-evaluate. This is data
contamination, even if it feels harmless.

**Fix.** Set the test set *before* training. Lock it. If you must
re-split, restart the training from scratch on the new split. Your
final test-set performance must be reported on a split you committed to
before seeing the results.

### 11. Confusing training distribution with reality

A model trained on the Materials Project database is biased by what is
*in* MP: mostly DFT-relaxed inorganic structures. Predictions on
disordered alloys, surfaces, or molecular crystals are out of
distribution and unreliable.

**Detection.** Check whether your prediction target lies inside the
support of the training data. Tools like t-SNE or UMAP plots of the
training-data fingerprints versus the prediction targets give a quick
visualisation.

**Fix.** Either restrict your conclusions to the support of the
training data, or augment the training data, or use uncertainty
estimates (Gaussian processes, MLP ensembles) to flag out-of-domain
predictions.

## General numerical pitfalls

### 12. Mixing units

You read a paper that quotes formation energies in eV. Your DFT code
outputs in Hartree. You convert in your head and miss a factor of 27.
The thesis goes out with the wrong number.

**Detection.** Sanity check every reported number against the order of
magnitude you expect. Defect formation energies in metals are 1-3 eV
typically; if you compute 30 or 0.1, something is wrong with units.

**Fix.** Use a unit-aware library (`ase.units`, `astropy.units`,
`pint`) in your analysis scripts. Never type unit conversions by hand
in code. Always write out the units in your figures and tables.

### 13. Wrong sign convention

DFT codes have different sign conventions for forces (force *on* atom
vs. gradient *of* energy). MLIPs trained on data from one code may
produce wrong-sign forces when interfaced to another. Symptom: MD
that "explodes" instead of equilibrating.

**Fix.** When interfacing two codes, run a small sanity calculation
(say, a 5-atom Lennard-Jones cluster) and verify that both codes give
the same sign of forces.

### 14. Atom indexing off by one

You wrote a script to swap atom 5 with atom 6 to create a defect.
You forgot that Python is 0-indexed and ASE is 0-indexed, but your
input file numbered atoms from 1. Your "defect" is the wrong atom.

**Fix.** Always sanity-check by visualising the structure (in OVITO, VESTA,
or ASE GUI) after constructing it. A 30-second visual inspection
prevents weeks of silently-wrong calculations.

## Reporting pitfalls

### 15. No uncertainty on reported numbers

A formation energy of 1.234 eV — to four significant figures — implies
a precision the calculation does not have. Numerical convergence is
typically 10-20 meV. Pseudopotential and functional choice introduce
0.1-0.5 eV uncertainty. Statistical noise in MD adds more.

**Fix.** Report numbers to a precision consistent with their actual
uncertainty. "1.2 eV", or "1.23 ± 0.05 eV" if you have a defensible
uncertainty estimate. Never quote a result to more digits than you can
defend.

### 16. Cherry-picking results

You ran 20 calculations; 3 gave the result you expected and 17 did not.
You report the 3. This is fraud, even if you did not think of it that
way. The dishonest framing is *the temptation to view results that
support your hypothesis as "the real ones" and inconvenient results as
"errors to debug"*.

**Fix.** Pre-commit to a protocol: which calculations you will run,
what counts as a success, what you will do with outliers. Report the
full distribution of outcomes, not just the favourable subset.

### 17. Confusing precision with accuracy

You ran a calculation with extremely tight numerical convergence (1 μeV
total energy) and reported the result with high confidence. But the
exchange-correlation functional you used has a systematic error of
0.5 eV for this kind of problem. Your precise number is inaccurate.

**Fix.** Distinguish *numerical convergence* (how well your code has
solved the equations you gave it) from *systematic error* (how well
the equations you gave it represent reality). Quote both, when you can,
or at least flag which you are reporting.

## Code and reproducibility pitfalls

### 18. Hard-coded paths and parameters

Your analysis script has lines like
`base_dir = "/home/jdoe/cluster/projects/v3/iron_runs/"`
and `n_atoms = 128`. When the thesis is done and someone else tries to
re-run, the script fails immediately.

**Fix.** Use relative paths, configuration files, or command-line
arguments. Externalise every parameter that might change.

### 19. No version control

You worked on the project for six months in a single directory. There
is no record of when each file was created or modified. You have ten
versions of `analyse.py` (`analyse_v2.py`, `analyse_FINAL.py`,
`analyse_FINAL_v2.py`).

**Fix.** Use git from day one. Commit early, commit often. A thesis
without a git history is a thesis that cannot be safely modified in the
last week before submission.

### 20. No data archiving plan

The thesis is submitted. Six months later, a reviewer asks for the raw
data. You realise that the cluster account has been wiped and you have
no local copy.

**Fix.** Before submitting the thesis, archive the input files,
output files (or at least the essential ones), analysis scripts, and
final figures to a long-term store: an institutional repository, a
Zenodo deposit, a personal hard drive. Cite the archive in the thesis.

## A diagnostic exercise

For each pitfall above, mentally check your current project:

- Have I tested this?
- Where in my thesis does the check appear?
- Could I produce the evidence in 30 seconds if asked?

If you cannot answer yes to all three for a given pitfall, address it
this week. Pitfalls compound; one is annoying, three is a problem, ten
is a thesis that the examiner will reject.

[Section 6](06-writing-up.md) is about turning carefully-converged,
carefully-validated work into a written thesis.
