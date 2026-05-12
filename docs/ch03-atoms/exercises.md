# Chapter 3 — Exercises

Seven problems, mixing concept and computation. Difficulty: <span class="diff-easy">★</span> easy, <span class="diff-medium">★★</span> medium, <span class="diff-hard">★★★</span> hard. Solutions are collapsed; attempt before peeking.

---

## Exercise 3.1 — Valence electrons and pseudopotentials <span class="diff-easy">★</span>

For each of the following elements, write the ground-state electron configuration, identify which electrons are *valence* in the typical pseudopotential treatment, and state how many valence electrons would appear in a DFT calculation using the standard pseudopotential.

(a) Aluminium (Al, $Z = 13$).
(b) Iron (Fe, $Z = 26$).
(c) Lead (Pb, $Z = 82$).
(d) Oxygen (O, $Z = 8$).

??? success "Solution"
    (a) Al: $1s^2 2s^2 2p^6 3s^2 3p^1 = [\text{Ne}] 3s^2 3p^1$. Valence: $3s^2 3p^1$. Three valence electrons.

    (b) Fe: $[\text{Ar}] 3d^6 4s^2$. The standard valence treatment in modern PAW pseudopotentials is *eight* electrons ($3d^6 4s^2$). For high accuracy in magnetic systems some users freeze fewer states and include the $3p^6$ semi-core, giving fourteen.

    (c) Pb: $[\text{Xe}] 4f^{14} 5d^{10} 6s^2 6p^2$. The valence is usually taken as $6s^2 6p^2$ (four electrons), with the $5d$ in the core. For high precision the $5d$ is sometimes treated as semi-core (fourteen electrons). Spin–orbit coupling is essential.

    (d) O: $1s^2 2s^2 2p^4$. Valence: $2s^2 2p^4$ (six electrons). The $1s^2$ is frozen.

    The pattern: lighter main-group elements have small, unambiguous valence; transition metals require choices about $d$ and (sometimes) $p$ semi-core; heavy elements are dominated by relativistic considerations. The pseudopotential file documents which valence configuration is used.

---

## Exercise 3.2 — Bond-type identification <span class="diff-easy">★</span>

For each of the following materials, identify the *dominant* bond type(s) and predict, qualitatively, (i) whether it is hard or soft, (ii) whether it conducts electricity at room temperature, and (iii) whether standard PBE-DFT without dispersion correction will give an accurate cohesive energy.

(a) Solid argon at 70 K.
(b) Sodium chloride (NaCl).
(c) Solid copper.
(d) Graphite.
(e) Crystalline naphthalene (C$_{10}$H$_8$).
(f) Diamond.

??? success "Solution"
    (a) Argon: pure van der Waals. (i) Very soft. (ii) Insulator (no carriers). (iii) **No.** PBE without vdW correction will give essentially zero binding. Use PBE-D3 or a vdW functional.

    (b) NaCl: ionic. (i) Hard and brittle. (ii) Insulator at room temperature; conducts via ion motion at high temperature or in melt. (iii) **Yes.** Long-range Coulomb is captured naturally in DFT.

    (c) Copper: metallic. (i) Soft to medium; ductile. (ii) Excellent conductor. (iii) **Yes.** Metals are easy for DFT (modulo k-point convergence).

    (d) Graphite: covalent in-plane, vdW out-of-plane. (i) Soft in the layer-stacking direction (cleavage), hard within layers. (ii) Conducts in-plane (semi-metal), poorly out-of-plane. (iii) **No, for the interlayer binding.** Without vdW correction, the interlayer distance and binding will be wrong.

    (e) Naphthalene: covalent within molecule, vdW between molecules. (i) Soft molecular crystal. (ii) Insulator. (iii) **No.** Molecular crystals are dispersion-bound.

    (f) Diamond: pure covalent. (i) The hardest naturally occurring material. (ii) Insulator. (iii) **Yes.** Pure covalent systems with no long range polarisation are well-handled by PBE.

---

## Exercise 3.3 — Building FCC copper and counting neighbours <span class="diff-medium">★★</span>

Using ASE, build a 4 × 4 × 4 conventional supercell of FCC copper. Then:

(a) How many atoms does it contain?

(b) What is the volume of the supercell?

(c) Pick one atom near the centre. How many atoms lie within $3.0$ Å of it (using the minimum image convention)? Compare to the theoretical first-shell coordination number for FCC.

(d) How many atoms lie within $4.0$ Å of it? Identify which shells of the FCC structure are included.

Provide the code and a brief discussion of the output.

??? success "Solution"
    ```python
    import numpy as np
    from ase.build import bulk

    cu = bulk("Cu", crystalstructure="fcc", cubic=True).repeat((4, 4, 4))
    print(f"(a) atoms: {len(cu)}")
    print(f"(b) volume: {cu.get_volume():.2f} Å^3")
    print(f"    lattice parameter: {cu.cell[0, 0] / 4:.3f} Å")

    centre_index = len(cu) // 2
    distances = cu.get_distances(centre_index, range(len(cu)), mic=True)
    distances = distances[distances > 0]

    n_within_3 = int(np.sum(distances < 3.0))
    n_within_4 = int(np.sum(distances < 4.0))
    print(f"(c) atoms within 3.0 Å: {n_within_3}")
    print(f"(d) atoms within 4.0 Å: {n_within_4}")
    ```

    Expected output:

    - (a) $4 \text{ atoms per conventional cell} \times 4^3 = 256$ atoms.
    - (b) Cu lattice parameter is $a = 3.615$ Å (ASE default); volume $= (4 \times 3.615)^3 \approx 3024$ Å$^3$.
    - (c) Theoretical FCC first-shell distance: $a / \sqrt{2} = 3.615 / 1.414 = 2.557$ Å. The next shell (second neighbours) is at $a = 3.615$ Å. So *within 3.0 Å* picks up the first shell only: **12 atoms**.
    - (d) Within 4.0 Å we include shells at 2.557 Å (first, 12 atoms) and 3.615 Å (second, 6 atoms). Total: **18 atoms**.

    The agreement with theoretical coordination numbers (12 for the first FCC shell, 6 for the second) confirms the structure is correct.

---

## Exercise 3.4 — Miller indices <span class="diff-medium">★★</span>

Consider a cubic crystal with lattice parameter $a$.

(a) Draw or describe the (100), (110), and (111) planes. Sketch their orientations relative to a unit cube.

(b) What is the interplanar spacing $d_{hkl}$ for each plane in a cubic crystal? Use the formula $d_{hkl} = a / \sqrt{h^2 + k^2 + \ell^2}$.

(c) For diamond silicon ($a = 5.43$ Å), at what angle would $(111)$ X-ray reflections appear using Cu K$\alpha$ radiation with $\lambda = 1.5406$ Å? Use Bragg's law, $n \lambda = 2 d_{hkl} \sin\theta$, with $n = 1$.

(d) In FCC, the $(100)$ reflection is *forbidden* — it has zero structure factor. Explain why.

??? success "Solution"
    (a) **(100)** is a plane perpendicular to the $x$-axis, passing through the face of the cube. **(110)** is a plane containing the $z$-axis, slicing diagonally through opposite edges of the cube. **(111)** is a plane intersecting all three axes at unit distance — the plane that cuts the cube into a tetrahedron at one corner.

    (b) $d_{100} = a$. $d_{110} = a / \sqrt{2}$. $d_{111} = a / \sqrt{3}$.

    (c) For Si, $d_{111} = 5.43 / \sqrt{3} = 3.135$ Å.
    Bragg: $\sin\theta = \lambda / (2d) = 1.5406 / 6.270 = 0.2457$.
    $\theta = \arcsin(0.2457) = 14.22°$. The diffraction angle is $2\theta = 28.44°$. This is in fact the strongest reflection in the silicon powder pattern.

    (d) The FCC basis has atoms at $(0,0,0)$, $(1/2, 1/2, 0)$, $(1/2, 0, 1/2)$, $(0, 1/2, 1/2)$ in the conventional cell. The structure factor for reflection $(hk\ell)$ is

    $$
    F_{hk\ell} = f \big[1 + e^{i\pi(h+k)} + e^{i\pi(h+\ell)} + e^{i\pi(k+\ell)}\big].
    $$

    For $(100)$ this is $f[1 + e^{i\pi} + 1 + e^{i\pi}] = f[1 - 1 + 1 - 1] = 0$. The reflection is *systematically absent* — a signature of the face centring. In general, FCC reflections are allowed only when $h, k, \ell$ are *all even or all odd*. So $(111)$, $(200)$, $(220)$, $(311)$ are allowed; $(100)$, $(110)$, $(210)$ are forbidden.

---

## Exercise 3.5 — Reciprocal lattice of primitive cubic <span class="diff-medium">★★</span>

Compute the reciprocal-lattice vectors of a simple cubic lattice with conventional vectors $\mathbf{a}_1 = a \hat x$, $\mathbf{a}_2 = a \hat y$, $\mathbf{a}_3 = a \hat z$. Do this by hand using $\mathbf{b}_i \cdot \mathbf{a}_j = 2\pi \delta_{ij}$.

Then, in code, verify your answer using either ASE's `cell.reciprocal()` or direct matrix inversion, and confirm the duality relation numerically.

What is the volume of the reciprocal-space primitive cell?

??? success "Solution"
    Demanding $\mathbf{b}_1 \cdot \mathbf{a}_1 = 2\pi$ and $\mathbf{b}_1 \cdot \mathbf{a}_2 = \mathbf{b}_1 \cdot \mathbf{a}_3 = 0$: $\mathbf{b}_1$ must be parallel to $\hat x$, with magnitude $2\pi/a$.

    By the same logic:
    $$\mathbf{b}_1 = (2\pi/a) \hat x, \quad \mathbf{b}_2 = (2\pi/a) \hat y, \quad \mathbf{b}_3 = (2\pi/a) \hat z.$$

    The reciprocal of simple cubic is simple cubic, with lattice parameter $2\pi / a$. The primitive cell volume in reciprocal space is $(2\pi/a)^3 = (2\pi)^3 / a^3 = (2\pi)^3 / V_\text{direct}$.

    ```python
    import numpy as np

    a = 4.0  # arbitrary
    A = a * np.eye(3)              # rows are a_i
    B = 2 * np.pi * np.linalg.inv(A).T   # rows are b_i
    print("Reciprocal vectors (rows):")
    print(B)
    print(f"|b_i| = {np.linalg.norm(B[0]):.4f}, expected {2*np.pi/a:.4f}")
    print("Duality check (should be identity):")
    print(np.round(A @ B.T / (2 * np.pi), 6))
    ```

    The duality matrix is the identity to machine precision.

---

## Exercise 3.6 — Brillouin zone of FCC silicon <span class="diff-hard">★★★</span>

Use ASE to construct primitive silicon and compute its reciprocal-lattice vectors. Then:

(a) Confirm numerically that the reciprocal lattice of FCC silicon is *body-centred cubic*, with conventional cube edge $4\pi/a$ where $a$ is the FCC conventional lattice parameter.

(b) Compute the volume of the first Brillouin zone of silicon and compare to the simple-cubic Brillouin zone of the same conventional cell. Why is the FCC primitive BZ four times larger than the simple-cubic BZ of the conventional cell?

(c) For a Monkhorst–Pack mesh of $8 \times 8 \times 8$ on the primitive cell, how many k-points are sampled per unit reciprocal-space volume? What density would you need on a $2 \times 2 \times 2$ supercell of the primitive cell to match this sampling density?

??? success "Solution"
    ```python
    import numpy as np
    from ase.build import bulk

    si = bulk("Si")          # primitive FCC, a_lat = 5.43 Å
    a = 5.43

    print("Direct cell:")
    print(si.cell.array)
    print("Reciprocal cell (rows, in 1/Å, with 2pi):")
    B = 2 * np.pi * np.linalg.inv(si.cell.array).T
    print(B)
    print("Lengths of b_i:", np.linalg.norm(B, axis=1))
    print("Expected (BCC primitive, b = sqrt(3) * 2pi / a):",
          np.sqrt(3) * 2 * np.pi / a)
    ```

    (a) The three $\mathbf{b}_i$ should each have length $\sqrt{3} \cdot 2\pi / a$ and should be arranged as the primitive vectors of a BCC lattice with conventional cube edge $4\pi / a$. Numerically (for $a = 5.43$ Å): $|\mathbf{b}_i| \approx 2.004$ Å$^{-1}$; $4\pi/a \approx 2.314$ Å$^{-1}$.

    (b) BZ volume of primitive FCC: $(2\pi)^3 / V_\text{primitive}$. With $V_\text{primitive} = a^3 / 4$, this is $4 (2\pi)^3 / a^3$.

    For the conventional cell (which has volume $a^3$), the would-be simple-cubic BZ has volume $(2\pi)^3 / a^3$ — a factor of 4 smaller. The factor of 4 is precisely the number of primitive cells in one conventional cell. Said differently: when you fold the FCC structure into a conventional cell, you also fold the BZ — replacing the primitive (large) BZ by the conventional (small) one, with bands folded back into it. The two are equivalent representations but the primitive BZ is the natural one for band-structure plots.

    (c) For the primitive cell, $8 \times 8 \times 8 = 512$ k-points cover a BZ volume of $V_\text{BZ}^\text{prim} = 4(2\pi)^3/a^3$. Density: $512 / V_\text{BZ}^\text{prim}$.

    For a $2 \times 2 \times 2$ supercell of the primitive cell, the BZ shrinks by a factor of 8 (cell volume grows 8×). To maintain the same k-point density, we need $512 / 8 = 64$ k-points, i.e. a $4 \times 4 \times 4$ mesh. The general rule: doubling the cell in a direction means *halving* the k-point count along that direction.

---

## Exercise 3.7 — Diagnosing a structure <span class="diff-hard">★★★</span>

You receive a CIF file from a collaborator describing a crystal but with no label. From the file you extract:

- conventional cell: $a = b = c = 4.04$ Å, $\alpha = \beta = \gamma = 90°$;
- two atoms per conventional cell;
- atom positions (fractional): atom 1 at $(0, 0, 0)$, atom 2 at $(1/2, 1/2, 1/2)$;
- both atoms are the same element with $Z = 26$.

(a) Identify the crystal structure and the element.

(b) What is the nearest-neighbour distance? Compare to typical metallic bond lengths.

(c) The collaborator says they obtained this structure by relaxing under high temperature. What might have happened to the well-known low-temperature structure of this element?

(d) If you wanted to simulate the *room-temperature* structure of this element, which crystal structure would you use, and how would the lattice parameter compare?

??? success "Solution"
    (a) Cubic cell, $a = b = c$, two atoms — one at the origin, one at the body centre. This is **body-centred cubic (BCC)**. The element with $Z = 26$ is **iron**.

    At high temperature ($> 1185$ K), iron is in its $\gamma$ phase, which is FCC. At lower temperatures it is $\alpha$-iron (BCC). So this structure is the *low-temperature* $\alpha$-Fe. (The lattice parameter of $\alpha$-Fe is indeed close to 2.87 Å for the primitive BCC cell; the value 4.04 Å here is suspect — it is closer to the *FCC* lattice parameter of $\gamma$-Fe at high temperature, around 3.65 Å, but in a BCC-like description). Let us re-examine.

    Actually, with $a = 4.04$ Å in a 2-atom cubic cell, the nearest-neighbour distance for BCC is $a\sqrt{3}/2 = 3.50$ Å, which is too large for a metallic bond. For real $\alpha$-Fe, $a_\text{BCC} = 2.87$ Å, giving $d_\text{nn} = 2.49$ Å. So either the collaborator has reported an unphysical structure, or the cell is *not* really BCC iron despite looking like one. A more plausible interpretation:

    A 4-atom FCC cell of iron has $a_\text{FCC} \approx 3.65$ Å. Two-atom cubic cells with $a \approx 4$ Å are unusual; one possibility is that the collaborator's relaxation went unstable and produced a phantom structure.

    (b) Reading the puzzle as stated (BCC with $a = 4.04$ Å), the nearest-neighbour distance is $a \sqrt 3 / 2 = 3.50$ Å. Typical Fe–Fe bond lengths are 2.5 Å. The reported structure has bonds far too long — this is unphysical for iron at standard pressure.

    (c) The most likely explanation is that the relaxation has produced an *expanded BCC* structure, perhaps because the calculation was run at very high temperature using a thermostat without proper pressure control, or because the starting structure was a metastable expanded lattice from which the relaxation never recovered. Alternatively, the magnetic configuration may be wrong — non-magnetic Fe in FCC has a different equilibrium volume from ferromagnetic BCC Fe.

    (d) For room-temperature iron, use **BCC** ($\alpha$-Fe) with $a \approx 2.87$ Å and *ferromagnetic* spin configuration (magnetic moment $\sim 2.2 \mu_B$ per atom). DFT with PBE and spin polarisation reproduces this to 1–2%. Forgetting the ferromagnetism gives a totally wrong volume — by a famous historical mistake, early LDA calculations predicted non-magnetic FCC iron to be lower in energy than the experimental ground state. The lesson: for transition metals, always check whether the magnetic state is correct.

    More broadly, this exercise illustrates the most important sanity check on any new structure: *do the nearest-neighbour distances look like real bond lengths?* If not, something is wrong.

---

When you are comfortable with all seven of these, you are ready for Chapter 4, where we will populate this real and reciprocal scaffolding with electrons and begin the descent into quantum mechanics proper.
