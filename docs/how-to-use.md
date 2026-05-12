# How to use this book

A short note on conventions, before you dive in.

## Mathematical notation

- Vectors are bold: $\mathbf{r}$, $\mathbf{k}$. Their magnitudes are italic: $r$, $k$.
- Operators wear a hat: $\hat{H}$, $\hat{T}$, $\hat{V}$.
- $\mathbb{R}$ is the real numbers, $\mathbb{C}$ the complex numbers, $\mathbb{Z}$ the integers.
- We use $\partial_x f$ and $\frac{\partial f}{\partial x}$ interchangeably for partial derivatives.
- $\nabla$ is the gradient; $\nabla^2$ is the Laplacian.
- Inner products are written with bra-ket notation when convenient: $\langle\psi|\phi\rangle = \int \psi^* \phi\, d\tau$.

## Units

Materials simulation uses an unholy mixture of unit systems. We try to be explicit, but the conventions are:

- **Energies** in electronvolts (eV) for chemistry, hartree (Ha) inside DFT codes, joules (J) for thermodynamics. $1\,\mathrm{Ha} = 27.2114\,\mathrm{eV}$.
- **Lengths** in ångströms (Å) for chemistry, bohr ($a_0$) inside DFT codes. $1\,a_0 = 0.529177\,\text{Å}$.
- **Temperatures** in kelvin (K).
- **Times** in femtoseconds (fs) for MD; $1\,\mathrm{fs} = 10^{-15}\,\mathrm{s}$.

A table at the back of [Appendix A](appendix/A-math-reference.md) gives all the conversion factors you will ever need.

## Code conventions

All worked examples are in Python 3.11. We use:

- **NumPy** for numerical arrays
- **Matplotlib** for plots
- **ASE** (Atomic Simulation Environment) for structures and DFT/MD drivers
- **pymatgen** for crystallographic utilities and Materials Project I/O
- **PyTorch** for neural networks
- **PyTorch Geometric** for graph neural networks

When a script can be run as-is, the file is named at the top of the listing and is available in the [`code/`](https://github.com/dongzhaohe321418-lab/materials-simulation-handbook/tree/main/code) directory of the repository.

## Difficulty markers

Exercises are marked:

- <span class="diff-easy">★ easy</span> — warm-up, mostly definitional
- <span class="diff-med">★★ medium</span> — standard problem-set difficulty
- <span class="diff-hard">★★★ hard</span> — research-level; expect to spend a serious afternoon

## Cross-references

A reference like *(Section 5.3)* points to the section with that number in this book. A reference like *(Martin §7.2)* points to a section in an external textbook listed in [Appendix D](appendix/D-reading.md).

## Running the code live in your browser

Every Python code block in this book has a small **▶ Run** button in its top-right corner. Click it and the code executes **inside your browser** — no install, no account, no server. The first run on any page takes ~10 seconds while the WebAssembly Python runtime is fetched and cached; after that, every subsequent run is instant.

What works in-browser: NumPy, SciPy, Matplotlib (all the maths, plotting, and small-system simulation code in Chapters 0–8 and most of 11).

What does NOT work in-browser: PyTorch, MACE, ASE, pymatgen, and other heavy scientific packages. Code blocks that import these (mostly Chapters 6, 9, 10, 12) display a polite note pointing you to [Google Colab](https://colab.research.google.com/), where the full stack is available with one click.

If a code block produces a plot, the figure renders directly under the **▶ Run** button. If it prints text, the text appears there. If it raises an error, the traceback is shown in red — debug it, edit the code, and run again.

!!! tip
    Editing the code in the page directly is not yet supported in the live runner; if you want to tweak a value and re-run, paste the snippet into the Python cell or open the corresponding notebook on Colab.

## Errata and contributions

If you spot a mistake — a typo, an incorrect derivation, a code example that does not run — please open an issue or a pull request on [GitHub](https://github.com/dongzhaohe321418-lab/materials-simulation-handbook). This book gets better when readers push back.
