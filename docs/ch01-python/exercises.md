# 1.5 Exercises

Eight problems. Difficulty markers follow [the conventions](../how-to-use.md): <span class="diff-easy">★</span> easy, <span class="diff-med">★★</span> medium, <span class="diff-hard">★★★</span> hard. Worked solutions are hidden in admonitions below each problem — try the problem first, peek when you are stuck.

You will need the `materials-simulation` environment from [Section 1.1](01-install.md) for the later problems.

---

## Exercise 1.1 — Array creation drill <span class="diff-easy">★</span>

Without looking back at Section 1.2, write one line of NumPy each that produces:

(a) a length-100 array of zeros of dtype `float32`;
(b) the 11 numbers $0.0, 0.1, 0.2, \dots, 1.0$;
(c) a $5 \times 5$ identity matrix;
(d) an array of 50 random samples from $\mathcal{N}(0, 1)$, with a seed of 7.

??? success "Solution"

    ```python
    import numpy as np

    a = np.zeros(100, dtype=np.float32)
    b = np.linspace(0.0, 1.0, 11)
    c = np.eye(5)
    d = np.random.default_rng(7).normal(size=50)
    ```

    Notes: `np.arange(0, 1.1, 0.1)` *appears* to work for (b) but is fragile — floating-point error can drop the last point. Use `linspace` whenever you know the number of points.

---

## Exercise 1.2 — Slicing a band-like matrix <span class="diff-easy">★</span>

Construct the $10 \times 10$ tridiagonal matrix with $2$ on the diagonal and $-1$ on the first off-diagonals (this is the discrete Laplacian on a chain). Then, using *one* slicing expression each:

(a) extract its diagonal as a 1D array;
(b) extract its main and super-diagonal as a $2 \times 10$ array (pad with `NaN` where needed);
(c) compute its smallest eigenvalue.

??? success "Solution"

    ```python
    import numpy as np

    n = 10
    A = 2 * np.eye(n) - np.eye(n, k=1) - np.eye(n, k=-1)

    diag = np.diag(A)                                       # (a)

    super_ = np.concatenate([np.diag(A, k=1), [np.nan]])    # (b)
    band = np.vstack([diag, super_])

    smallest = np.linalg.eigvalsh(A).min()                  # (c)
    print(diag, band, smallest)
    ```

    The smallest eigenvalue is approximately $2(1 - \cos(\pi / (n+1))) \approx 0.0810$.

---

## Exercise 1.3 — Plotting the Lennard-Jones potential <span class="diff-easy">★</span>

The Lennard-Jones potential

$$
V(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]
$$

models the interaction between two neutral atoms. Plot $V(r)$ in units of $\epsilon$ versus $r/\sigma$ for $r/\sigma \in [0.9, 3.0]$. Mark the minimum and report its location *analytically*.

??? success "Solution"

    Analytically, $\mathrm{d}V/\mathrm{d}r = 0$ gives $r_{\min} = 2^{1/6}\sigma \approx 1.1225\sigma$, with $V(r_{\min}) = -\epsilon$.

    ```python
    import numpy as np
    import matplotlib.pyplot as plt

    r_over_sigma = np.linspace(0.9, 3.0, 400)
    V_over_eps = 4.0 * (r_over_sigma**-12 - r_over_sigma**-6)

    rmin = 2.0 ** (1.0 / 6.0)

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(r_over_sigma, V_over_eps, color="C0")
    ax.axhline(0, color="0.6", linewidth=0.5)
    ax.scatter([rmin], [-1.0], color="C3", zorder=3)
    ax.annotate(rf"$r_\mathrm{{min}} = 2^{{1/6}}\sigma \approx {rmin:.4f}\sigma$",
                xy=(rmin, -1), xytext=(1.5, -0.5),
                arrowprops=dict(arrowstyle="->", color="0.4"), fontsize=9)

    ax.set_xlabel(r"$r / \sigma$")
    ax.set_ylabel(r"$V / \epsilon$")
    ax.set_ylim(-1.2, 2.0)
    fig.tight_layout()
    fig.savefig("lj.pdf")
    ```

---

## Exercise 1.4 — Distance matrix function <span class="diff-med">★★</span>

Write a function

```python
def distance_matrix(points: np.ndarray) -> np.ndarray:
    ...
```

that takes an $(N, 3)$ array of Cartesian positions and returns the $(N, N)$ matrix of pairwise Euclidean distances. Use broadcasting; no Python loops. Verify on the four corners of a unit tetrahedron that the off-diagonal distances are all equal.

??? success "Solution"

    ```python
    import numpy as np

    def distance_matrix(points: np.ndarray) -> np.ndarray:
        diff = points[:, None, :] - points[None, :, :]   # (N, N, 3)
        return np.linalg.norm(diff, axis=-1)             # (N, N)

    tetra = np.array([
        [ 1,  1,  1],
        [ 1, -1, -1],
        [-1,  1, -1],
        [-1, -1,  1],
    ], dtype=float)
    D = distance_matrix(tetra)
    print(D)
    # Off-diagonal entries should all equal 2*sqrt(2) ≈ 2.8284
    assert np.allclose(D[np.triu_indices(4, k=1)], 2 * np.sqrt(2))
    ```

    For very large $N$, the $(N, N, 3)$ intermediate becomes memory-hungry — `scipy.spatial.distance.cdist` is more memory-efficient. For up to a few thousand atoms, the broadcasting version is faster and simpler.

---

## Exercise 1.5 — Vectorising a naive loop <span class="diff-med">★★</span>

Here is a naive loop that computes the mean nearest-neighbour distance from a set of points to a set of centres:

```python
def naive_mean_nn(points: np.ndarray, centres: np.ndarray) -> float:
    total = 0.0
    for p in points:
        best = float("inf")
        for c in centres:
            d = ((p - c) ** 2).sum() ** 0.5
            if d < best:
                best = d
        total += best
    return total / len(points)
```

Rewrite as a single vectorised function with no Python loops. Time both on $N = 500$ points and centres and confirm at least a 100× speed-up.

??? success "Solution"

    ```python
    import numpy as np
    import time

    def vectorised_mean_nn(points: np.ndarray, centres: np.ndarray) -> float:
        diff = points[:, None, :] - centres[None, :, :]   # (N, M, d)
        d2 = (diff ** 2).sum(axis=-1)                      # (N, M)
        return np.sqrt(d2.min(axis=1)).mean()

    rng = np.random.default_rng(0)
    P = rng.normal(size=(500, 3))
    C = rng.normal(size=(500, 3))

    t0 = time.perf_counter()
    a = naive_mean_nn(P, C)
    t_naive = time.perf_counter() - t0

    t0 = time.perf_counter()
    b = vectorised_mean_nn(P, C)
    t_vec = time.perf_counter() - t0

    assert np.isclose(a, b)
    print(f"naive: {t_naive:.3f} s, vectorised: {t_vec:.4f} s, speed-up: {t_naive / t_vec:.1f}x")
    ```

    Two implementation notes:

    - Squaring then square-rooting at the end is cheaper than computing $\| \mathbf{p} - \mathbf{c} \|$ inside the broadcasting. Taking `sqrt` only after the `min` is also valid because $\sqrt{\cdot}$ is monotonic.
    - For $N \gtrsim 10^4$, the $(N, M)$ intermediate becomes too large. Either chunk the points or use `scipy.spatial.cKDTree`.

---

## Exercise 1.6 — Numerical minimum of a potential <span class="diff-med">★★</span>

Take the Morse potential from Section 1.3 with $D_e = 4.75$ eV, $a = 1.94$ Å$^{-1}$, $r_0 = 0.741$ Å. Use `scipy.optimize.minimize_scalar` to find the minimum *numerically*, starting from $r = 1.2$ Å. Confirm that it matches $r_0$ to six decimal places. Then, *without* changing the optimiser, compute the harmonic-approximation force constant $k = V''(r_0)$ using a central finite difference, and from it the vibrational frequency $\omega = \sqrt{k / \mu}$ for an H$_2$-like reduced mass $\mu = 0.5\,m_\mathrm{H}$.

??? success "Solution"

    ```python
    import numpy as np
    from scipy.optimize import minimize_scalar

    D_e, a, r0 = 4.75, 1.94, 0.741

    def V(r):
        return D_e * (1 - np.exp(-a * (r - r0)))**2 - D_e

    res = minimize_scalar(V, bracket=(0.5, 1.2, 3.0))
    r_star = res.x
    print(f"r* = {r_star:.6f} Å (target {r0})")
    assert abs(r_star - r0) < 1e-6

    # Second derivative by central difference
    h = 1e-4
    Vpp = (V(r_star + h) - 2 * V(r_star) + V(r_star - h)) / h**2
    print(f"V''(r*) = {Vpp:.4f} eV/Å²")

    # Convert k to SI and compute omega.
    # 1 eV/Å² = 16.0218 N/m. m_H = 1.6735e-27 kg, mu = 0.5 m_H.
    k_SI = Vpp * 16.0218
    mu = 0.5 * 1.6735e-27
    omega = np.sqrt(k_SI / mu)          # rad/s
    nu_cm = omega / (2 * np.pi * 3.0e10)
    print(f"omega = {omega:.3e} rad/s, wavenumber ≈ {nu_cm:.0f} cm⁻¹")
    ```

    Analytically, $V''(r_0) = 2 D_e a^2 \approx 35.74$ eV/Å$^2$. The finite difference should reproduce this to four decimal places. The vibrational wavenumber for real H$_2$ is about 4400 cm$^{-1}$; the Morse parameters above are tuned to give that order of magnitude.

---

## Exercise 1.7 — A reproducible mini-project <span class="diff-med">★★</span>

Create a new directory `lj-cluster/` containing:

1. an `environment.yml` with `numpy`, `scipy`, `matplotlib`;
2. a `.gitignore` that excludes `__pycache__`, PDFs, and a `results/` directory;
3. a script `cluster_energy.py` that builds 13 atoms in an icosahedral arrangement at unit spacing, computes the total Lennard-Jones energy ($\epsilon = \sigma = 1$) using the function from Exercise 1.4, and writes the energy together with a `provenance.json` (as in Section 1.4) into `results/run_001/`.

Initialise git, commit, and verify with `git log`.

??? success "Solution"

    ```bash
    mkdir lj-cluster && cd lj-cluster
    git init && git branch -M main
    ```

    `environment.yml`:

    ```yaml
    name: lj-cluster
    channels: [conda-forge]
    dependencies:
      - python=3.11
      - numpy
      - scipy
      - matplotlib
    ```

    `.gitignore`:

    ```gitignore
    __pycache__/
    *.pyc
    results/
    *.pdf
    .DS_Store
    ```

    `cluster_energy.py`:

    ```python
    import json
    import sys
    import platform
    from datetime import datetime, timezone
    from importlib.metadata import version, PackageNotFoundError
    from pathlib import Path
    import numpy as np

    def icosahedron_13(scale: float = 1.0) -> np.ndarray:
        phi = (1 + 5 ** 0.5) / 2
        verts = []
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                verts.append([0, s1, s2 * phi])
                verts.append([s1, s2 * phi, 0])
                verts.append([s2 * phi, 0, s1])
        verts.append([0, 0, 0])
        return scale * np.asarray(verts, dtype=float)

    def distance_matrix(points: np.ndarray) -> np.ndarray:
        diff = points[:, None, :] - points[None, :, :]
        return np.linalg.norm(diff, axis=-1)

    def lj_energy(points: np.ndarray, eps: float = 1.0, sigma: float = 1.0) -> float:
        D = distance_matrix(points)
        i, j = np.triu_indices(len(points), k=1)
        r = D[i, j]
        sr6 = (sigma / r) ** 6
        return float((4 * eps * (sr6 ** 2 - sr6)).sum())

    def provenance() -> dict:
        out = {}
        for p in ("numpy", "scipy", "matplotlib"):
            try:
                out[p] = version(p)
            except PackageNotFoundError:
                out[p] = "missing"
        return {
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "packages": out,
        }

    if __name__ == "__main__":
        points = icosahedron_13(scale=1.0)
        E = lj_energy(points)
        out_dir = Path("results/run_001")
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "energy.txt").write_text(f"{E:.8f}\n")
        (out_dir / "provenance.json").write_text(json.dumps(provenance(), indent=2))
        print(f"Total LJ energy: {E:.6f}")
    ```

    Then:

    ```bash
    git add environment.yml .gitignore cluster_energy.py
    git commit -m "initial 13-atom LJ cluster energy"
    git log --oneline
    ```

    Because `results/` is gitignored, the output files are *not* committed. The recipe to regenerate them is.

---

## Exercise 1.8 — Minimum-image radial distribution function <span class="diff-hard">★★★</span>

Write a function

```python
def rdf(positions: np.ndarray,
        cell: np.ndarray,
        r_max: float,
        n_bins: int) -> tuple[np.ndarray, np.ndarray]:
    ...
```

that computes the **radial distribution function** $g(r)$ of a 3D periodic system from a set of atomic positions and a $3\times 3$ lattice matrix. The minimum-image convention is required; do not assume an orthorhombic cell.

Test on a $4 \times 4 \times 4$ silicon supercell built with `ase.build.bulk("Si", cubic=True).repeat((4,4,4))`. Plot $g(r)$ from 0 to 6 Å with 200 bins. Confirm a first peak near 2.35 Å and a second near 3.84 Å.

State, in two sentences, why a naive Cartesian distance matrix gives the wrong answer here.

??? success "Solution"

    A naive distance matrix only counts displacements within the simulation box. In a periodic crystal, an atom is also a neighbour of every periodic image of every other atom; the *minimum-image convention* picks, for each pair, the periodic copy that lies closest, which is what physical observables like $g(r)$ depend on.

    ```python
    import numpy as np
    import matplotlib.pyplot as plt
    from ase.build import bulk

    def rdf(positions: np.ndarray,
            cell: np.ndarray,
            r_max: float,
            n_bins: int) -> tuple[np.ndarray, np.ndarray]:
        # Convert to fractional coordinates: f = (cell^T)^-1 @ r
        inv_cell_T = np.linalg.inv(cell.T)
        frac = positions @ inv_cell_T.T                     # (N, 3)

        # Pairwise fractional differences, wrapped into [-0.5, 0.5)
        df = frac[:, None, :] - frac[None, :, :]             # (N, N, 3)
        df -= np.round(df)
        # Back to Cartesian
        diff = df @ cell                                     # (N, N, 3)
        D = np.linalg.norm(diff, axis=-1)

        i, j = np.triu_indices(len(positions), k=1)
        d = D[i, j]
        d = d[d < r_max]

        edges = np.linspace(0.0, r_max, n_bins + 1)
        counts, _ = np.histogram(d, bins=edges)
        centres = 0.5 * (edges[:-1] + edges[1:])

        # Normalise: g(r) = counts / (rho * 4 pi r^2 dr * N/2)
        N = len(positions)
        volume = abs(np.linalg.det(cell))
        rho = N / volume
        dr = edges[1] - edges[0]
        shell = 4.0 * np.pi * centres**2 * dr
        ideal_pairs_per_atom = rho * shell
        ideal_total_pairs = ideal_pairs_per_atom * N / 2.0
        g = counts / np.maximum(ideal_total_pairs, 1e-30)
        return centres, g

    si = bulk("Si", cubic=True).repeat((4, 4, 4))
    r, g = rdf(si.get_positions(), si.cell.array, r_max=6.0, n_bins=200)

    fig, ax = plt.subplots(figsize=(4.5, 3))
    ax.plot(r, g, color="C0")
    ax.axhline(1.0, color="0.6", linewidth=0.5, linestyle="--")
    ax.set_xlabel("$r$ (Å)")
    ax.set_ylabel("$g(r)$")
    ax.set_xlim(0, 6)
    fig.tight_layout()
    fig.savefig("si_rdf.pdf")

    first  = r[np.argmax(g)]
    print(f"first peak: {first:.2f} Å")
    ```

    With 4×4×4×8 = 512 atoms the distance tensor uses about 2 MB and is fine on any laptop. Beyond ~5000 atoms the broadcasting version exhausts memory; either chunk the outer index or use a neighbour list such as `ase.neighborlist.NeighborList` or `matscipy.neighbours.neighbour_list`.

    Bonus question to ponder: does your normalisation produce $g(r) \to 1$ at large $r$? Why or why not, for a crystal of 512 atoms?
