# Tier-1 code

Importable Python modules extracted from the pedagogical scripts in the
chapter markdown. Each module is testable in isolation; the matching
test files live in `tests/` at the repository root.

The chapter markdown remains the canonical source for prose and
derivations. The code here is for execution and regression testing only.

## Layout

| Path | Source chapter | Contents |
|---|---|---|
| `ch03b/graphene_tb.py` | 3b.3 Tight-binding | Graphene 2x2 TB Hamiltonian, K-Gamma-M-K band path. |
| `ch03b/phonons_diatomic.py` | 3b.5 Phonons | 1D diatomic-chain dynamical matrix, acoustic/optical branches. |
| `ch04/schrodinger_1d.py` | 4.3 Particle in a box | Finite-difference 1D Schrodinger solver. |
| `ch04/harmonic_oscillator.py` | 4.4 Harmonic oscillator | Numerical 1D quantum SHO via FD. |
| `ch05/scf_1d.py` | 5.5 SCF loop | Kohn-Sham SCF in 1D with LDA exchange and Pulay/DIIS mixing. |
| `ch07/velocity_verlet.py` | 7.1 Integration | Velocity-Verlet integrator and a 1D harmonic-oscillator demo. |
| `ch07/trajectory_analysis.py` | 7.6 Trajectory analysis | MSD, RDF, VACF/VDOS and a diffusion-coefficient fit. |
| `ch09/descriptors_g2.py` | 9.3 Descriptors | Behler G^2 radial symmetry function. |
| `ch10/cgcnn.py` | 10.3 CGCNN | CGCNN model in PyTorch Geometric (skipped if torch absent). |
| `ch11/gp.py` | 11.2 Gaussian processes | RBF kernel and GP with marginal-likelihood optimisation. |
| `ch11/acquisition.py` | 11.3 Acquisition | Expected improvement, UCB, and a short BO loop. |

## Running

From the repository root:

```bash
pytest tests/ -v --tb=short
```

Tests that depend on `torch` or `torch_geometric` are skipped
automatically when those libraries are not installed.
