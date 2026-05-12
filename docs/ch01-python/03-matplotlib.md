# 1.3 Plotting with Matplotlib

A plot is the most efficient way to convince yourself — and a reader — that a calculation did what it should. Matplotlib is the default plotting library in scientific Python, and despite a sometimes idiosyncratic API, it produces publication-quality output and is supported by every journal's figure pipeline.

## The figure / axes mental model

Matplotlib has two layers you must keep separate in your head.

- A **figure** is the whole canvas — the rectangle of pixels (or PDF points) that will eventually be saved or shown.
- An **axes** is a single coordinate system inside the figure — what most non-Matplotlib users call "a plot". A figure can contain one axes or many.

The single most useful function is `plt.subplots`, which returns both:

```python
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(4, 3))
x = np.linspace(0, 2 * np.pi, 200)
ax.plot(x, np.sin(x), label=r"$\sin(x)$")
ax.set_xlabel("$x$")
ax.set_ylabel("$\\sin(x)$")
ax.legend()
fig.tight_layout()
fig.savefig("sin.pdf")
plt.show()
```

The pattern is always: create a figure with one or more axes, draw on the axes, label the axes, save the figure. We almost never use the older `plt.plot(...)` style without first calling `subplots`, because it relies on a hidden "current axes" that becomes confusing the moment you have more than one panel.

!!! tip "Object-oriented or bust"
    Throughout this book we use `fig, ax = plt.subplots(...)` and `ax.plot(...)`, not `plt.plot(...)`. The object-oriented form is explicit and composes cleanly into multi-panel figures.

## Basic line and scatter plots

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-5, 5, 200)
y = np.tanh(x)

fig, ax = plt.subplots(figsize=(4, 3))
ax.plot(x, y, color="C0", linewidth=2, label=r"$\tanh(x)$")
ax.plot(x, np.sign(x), color="C1", linestyle="--", label=r"$\mathrm{sign}(x)$")

ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_xlim(-5, 5)
ax.set_ylim(-1.2, 1.2)
ax.axhline(0, color="0.5", linewidth=0.5)
ax.axvline(0, color="0.5", linewidth=0.5)
ax.legend(frameon=False, loc="lower right")
fig.tight_layout()
```

A scatter plot uses `ax.scatter`, which has a different signature: marker size and colour can be arrays of the same length as the data, encoding extra dimensions.

```python
rng = np.random.default_rng(0)
volumes = rng.uniform(40, 200, size=80)        # Å³
formation_energies = rng.normal(0, 0.5, size=80)  # eV/atom
bandgaps = rng.uniform(0, 4, size=80)          # eV

fig, ax = plt.subplots(figsize=(4.5, 3.5))
sc = ax.scatter(volumes, formation_energies,
                c=bandgaps, cmap="viridis",
                s=30, edgecolor="k", linewidth=0.3)
cbar = fig.colorbar(sc, ax=ax, label="Band gap (eV)")
ax.set_xlabel("Unit-cell volume (Å³)")
ax.set_ylabel("Formation energy (eV/atom)")
fig.tight_layout()
```

The colour bar is itself a tiny axes; `fig.colorbar` adds it next to the parent. Three quantitative variables, encoded clearly.

## Multiple subplots

`plt.subplots(nrows, ncols)` returns a 2D array of axes. The most common form:

```python
fig, axes = plt.subplots(1, 2, figsize=(7, 3), sharey=True)

x = np.linspace(0, 2 * np.pi, 200)
axes[0].plot(x, np.sin(x))
axes[0].set_title("sine")
axes[1].plot(x, np.cos(x))
axes[1].set_title("cosine")

for ax in axes:
    ax.set_xlabel("$x$")
axes[0].set_ylabel("$y$")
fig.tight_layout()
```

`sharey=True` makes both panels use the same $y$-axis — useful for comparison, and it lets Matplotlib hide the duplicate tick labels automatically.

For 2D grids of panels:

```python
fig, axes = plt.subplots(2, 3, figsize=(9, 6))
for i, ax in enumerate(axes.flat):
    ax.plot(np.random.default_rng(i).normal(size=50).cumsum())
    ax.set_title(f"trial {i}")
fig.tight_layout()
```

`axes.flat` flattens the 2D array into a 1D iterator — much cleaner than nested loops.

## Labels, legends, and ticks

Always include axis labels with units. A bare number on an axis is hostile to your reader.

```python
ax.set_xlabel("Wavevector $k$ (Å$^{-1}$)")
ax.set_ylabel("Energy (eV)")
```

LaTeX-style math is enabled by default; wrap it in dollar signs. Use raw strings (`r"$\sin$"`) to avoid escaping issues with backslashes.

`ax.legend()` reads the `label=` keyword from each plotted artist. Recommended options:

```python
ax.legend(frameon=False,        # no box around the legend
          loc="best",            # let Matplotlib pick the position
          fontsize="small")
```

Tick formatting is mostly automatic, but for crystallographic plots you may want explicit ticks (e.g., high-symmetry $k$-points):

```python
ax.set_xticks([0.0, 0.5, 1.0])
ax.set_xticklabels([r"$\Gamma$", "X", "L"])
```

## Logarithmic and twin axes

Switching one or both axes to log scale is a single call:

```python
ax.set_yscale("log")
ax.set_xscale("log")
```

Energy-volume curves and learning curves usually want a log $y$-axis; convergence-with-cutoff plots want both.

A "twin" axes shares an $x$-axis with another and gets its own $y$-axis on the right. Useful for showing two quantities with different units on the same energy or temperature axis:

```python
fig, ax = plt.subplots(figsize=(4.5, 3))
T = np.linspace(100, 1000, 50)
Cv = 3 * 8.314 * (T / (T + 200))           # arbitrary smooth function
S = 50 + 0.05 * (T - 100)

ax.plot(T, Cv, color="C0", label="$C_V$")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("$C_V$ (J K$^{-1}$ mol$^{-1}$)", color="C0")
ax.tick_params(axis="y", labelcolor="C0")

ax2 = ax.twinx()
ax2.plot(T, S, color="C1", label="$S$")
ax2.set_ylabel("$S$ (J K$^{-1}$ mol$^{-1}$)", color="C1")
ax2.tick_params(axis="y", labelcolor="C1")

fig.tight_layout()
```

## Saving publication-quality figures

For a paper or thesis, save as PDF or SVG, not PNG. Vector formats stay crisp at any zoom level and embed cleanly into LaTeX.

```python
fig.savefig("eos.pdf", bbox_inches="tight")
fig.savefig("eos.svg", bbox_inches="tight")

# PNG for slides or web preview — set the DPI high
fig.savefig("eos.png", dpi=300, bbox_inches="tight")
```

Sizing convention: a single-column journal figure is about 3.4 inches wide; a double-column figure is about 7 inches wide. Set `figsize` in inches and font sizes in points and you avoid almost all "the labels are too small" referee comments.

```python
import matplotlib as mpl
mpl.rcParams.update({
    "font.size": 9,
    "axes.labelsize": 9,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "lines.linewidth": 1.5,
})
```

Put this at the top of your plotting script, or save it as a `mplstyle` file and load it with `plt.style.use("paper.mplstyle")`.

## Styles and colours

Matplotlib ships with several built-in styles. Two we recommend trying:

```python
plt.style.use("seaborn-v0_8-darkgrid")   # soft grid, good for slides
plt.style.use("seaborn-v0_8-paper")       # smaller fonts, paper-ready
```

You can list all available styles with `plt.style.available`. The default cycle of line colours, `C0` through `C9`, is colour-blind-friendly (the *tab10* palette).

For custom colours, named CSS colours work (`"firebrick"`, `"steelblue"`); for arbitrary colours pass hex strings (`"#1f77b4"`); for colour maps use `plt.cm.viridis(0.4)` to evaluate the map at a fraction along its range.

!!! warning "Avoid the jet colour map"
    `jet` (and its cousin `rainbow`) introduce false visual structure where there is none and are inaccessible to colour-blind readers. Use **viridis**, **plasma**, or **cividis** for sequential data, and **RdBu_r** or **PiYG** for diverging data.

## A worked example: a 1D potential energy surface

Plot a Morse potential $V(r) = D_e\, [1 - e^{-a(r - r_0)}]^2 - D_e$ for an H–H-like bond, and mark its minimum.

```python
import numpy as np
import matplotlib.pyplot as plt

D_e = 4.75      # eV
a   = 1.94      # 1/Å
r0  = 0.741     # Å

r = np.linspace(0.4, 3.0, 400)
V = D_e * (1 - np.exp(-a * (r - r0)))**2 - D_e

fig, ax = plt.subplots(figsize=(4.5, 3.5))
ax.plot(r, V, color="C0")
ax.axhline(0, color="0.6", linewidth=0.5)
ax.axvline(r0, color="C3", linestyle=":", linewidth=1.0)
ax.scatter([r0], [-D_e], color="C3", zorder=3)
ax.annotate(f"minimum at $r_0 = {r0}$ Å,\n$V = -D_e = {-D_e}$ eV",
            xy=(r0, -D_e), xytext=(1.2, -1.0),
            arrowprops=dict(arrowstyle="->", color="0.4"),
            fontsize=8)

ax.set_xlabel("Bond length $r$ (Å)")
ax.set_ylabel("Potential energy $V$ (eV)")
ax.set_xlim(0.4, 3.0)
ax.set_ylim(-5.2, 5.0)
fig.tight_layout()
fig.savefig("morse.pdf")
```

A reader looking at this figure should understand, without consulting the caption, what is being plotted and where the equilibrium bond length lies.

## A band-structure stub

In Chapter 6 we plot a real DFT band structure. The plotting pattern is identical to what we already know — many lines on one axes — so here is a placeholder you can use to lay out a band-structure figure before you have actual eigenvalues:

```python
import numpy as np
import matplotlib.pyplot as plt

# Fake k-path with three segments: Γ-X-L-Γ
n_per_seg = 40
ks = np.concatenate([np.linspace(0, 1, n_per_seg),
                     np.linspace(1, 2, n_per_seg),
                     np.linspace(2, 3, n_per_seg)])
n_bands = 6
rng = np.random.default_rng(1)
bands = np.sort(rng.normal(size=(len(ks), n_bands))
                * np.linspace(0.2, 1.5, n_bands), axis=1)
bands += np.linspace(-4, 4, n_bands)

fig, ax = plt.subplots(figsize=(4, 4))
for b in range(n_bands):
    ax.plot(ks, bands[:, b], color="C0", linewidth=1.0)

for x in (1, 2):
    ax.axvline(x, color="0.7", linewidth=0.5)
ax.axhline(0, color="0.5", linewidth=0.5, linestyle="--")

ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels([r"$\Gamma$", "X", "L", r"$\Gamma$"])
ax.set_ylabel("Energy $E - E_F$ (eV)")
ax.set_xlim(0, 3)
fig.tight_layout()
```

The vertical lines mark the high-symmetry points and the dashed horizontal line marks the Fermi level. Real band structures use exactly this layout.

## A short style guide

Most of what makes a figure look professional is what you *do not* draw.

- **No chartjunk.** Drop the grid unless it aids reading. Drop the box on the legend. Drop redundant tick marks. Every pixel should carry information.
- **Label every axis with both quantity and unit.** "Energy" is not a label; "Energy (eV)" is.
- **Use colour deliberately.** Two colours is usually enough. Reserve red for "this is important" or "this is wrong"; do not use it for the third line in a series.
- **Match figure size to publication size.** Do not save a giant figure and scale it down in your paper — fonts will be illegible.
- **Use vector formats.** PDF or SVG, not PNG, for anything that will be printed.
- **One idea per figure.** If you find yourself describing four phenomena in one caption, split it into two figures.
- **Caption first, then plot.** Writing the caption forces you to decide what the figure is *for*. The figure should make exactly that point and no other.

The figures in this book obey these rules. Use them as templates. The exercise set at the end of the chapter has you reproduce the Morse potential figure and add a numerical minimum-finding step.
