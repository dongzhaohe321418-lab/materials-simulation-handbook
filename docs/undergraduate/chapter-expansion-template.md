# In-Place Chapter Expansion Template

This page is a template for *improving the chapter pages we already have*. It is not a plan for writing a new book, and it is not a chapter itself. The existing **Materials Simulation Handbook** is the source of truth. Our job in the undergraduate support layer is to add 台阶 (steps, or scaffolding) just *before* the difficult jumps, so that a normal undergraduate can climb where they previously fell off. We expand a page **in place**: we wrap the author's original material in supporting context, and we never replace it.

If you are new to the project, read [How to use this handbook](../how-to-use.md) and the [undergraduate index](index.md) first. For terminology, point readers at the [beginner glossary](glossary-for-beginners.md); for the surrounding plan, see the [expansion roadmap](expansion-roadmap-v1.4.md).

!!! warning "Golden rules — read before editing any chapter"
    - **Do not delete** the original explanation unless it is clearly wrong.
    - **Do not rewrite** existing prose into a new style or voice.
    - **Keep** all existing equations, worked examples, code and figures exactly as they are.
    - **Insert scaffolding *before*** a difficult jump, not after it.
    - **Keep advanced details**, but *label* them (for example with a difficulty tag) so beginners can skip them safely.
    - **Do not over-explain.** Not every sentence needs a gloss.
    - **Focus only** on the places where an undergraduate would actually get stuck.
    - **Keep diffs reviewable.** Small, well-targeted additions are easier to merge than a rewrite.

## How to use the blocks below

What follows is a menu of optional blocks. Add a block **only where it earns its place** — use your judgement. A short, clear page might need two or three of these; a dense derivation-heavy page might need most of them. Forcing every heading onto every page produces bloated, patronising pages, which is exactly what we are trying to avoid. When you do add a block, place it next to the original text it supports, and keep the original text intact above or below it.

---

### Keep the existing section

*Purpose:* the author's original explanation stays in place, unchanged unless it is clearly wrong. This is the anchor everything else hangs from.

```markdown
## The Kohn–Sham equations

[original author text, equations and figures — left exactly as written]
```

---

### What problem are we solving?

*Purpose:* state, in plain language, the question this section answers, before any formalism arrives.

```markdown
!!! info "What problem are we solving?"
    We want the energy of a material without solving the full
    many-electron Schrödinger equation, which is intractable for more
    than a handful of electrons. This section shows the trick that
    makes it possible.
```

---

### Plain-language version

*Purpose:* give the idea once in words, before any equation, so the reader knows what the maths is *for*.

```markdown
!!! note "Plain-language version"
    Density functional theory swaps a hard question ("what is every
    electron doing?") for an easier one ("what is the electron *density*
    doing?"). The density is a single function of position, not a
    function of every electron's coordinates at once.
```

---

### Physical picture

*Purpose:* describe what the atoms, electrons, energy, forces or structures are actually doing, so the symbols later have something concrete to attach to.

```markdown
!!! note "Physical picture"
    Picture the electrons as a smeared-out negative cloud, denser near
    the nuclei. Each electron feels the nuclei *and* the averaged cloud
    of all the others. The forces on the nuclei come from the gradient
    of the total energy as we slide an atom by a small distance.
```

---

### New vocabulary

*Purpose:* define terms *before* they are used heavily. Link to the [beginner glossary](glossary-for-beginners.md) rather than re-defining a term that already lives there.

```markdown
!!! tip "New vocabulary"
    - **Functional** — a rule that takes a *function* in and gives a
      *number* out (an energy, here). More in the
      [beginner glossary](glossary-for-beginners.md).
    - **Self-consistent** — we guess, compute, then re-guess until the
      answer stops changing.
```

---

### Symbol guide

*Purpose:* a small table of every symbol used nearby, with its meaning and units, so no symbol appears undefined. See the [formula reading guide](formula-reading-guide.md) for how to read an equation slowly.

```markdown
| Symbol | Meaning | Units |
|---|---|---|
| $n(\mathbf{r})$ | electron density at position $\mathbf{r}$ | electrons / $a_0^3$ |
| $E_{\text{xc}}[n]$ | exchange–correlation energy functional | Ha |
| $\hat{H}$ | Hamiltonian operator | — |
```

---

### Step-by-step mathematics

*Purpose:* break one big derivation into small, labelled steps. The original compressed derivation stays; this is the slow version beside it.

```markdown
!!! example "Step-by-step: from total energy to forces"
    1. Start from the total energy $E(\{\mathbf{R}_I\})$ as a function
       of the nuclear positions.
    2. The force on atom $I$ is $\mathbf{F}_I = -\,\partial E / \partial \mathbf{R}_I$.
    3. Because the energy is variational, only the *explicit* position
       dependence survives — this is the Hellmann–Feynman result.
```

---

### Minimal example

*Purpose:* a toy model — ideally 1-D or two-particle — that strips the idea to its bones.

```markdown
!!! example "Minimal example: two atoms on a line"
    Take two atoms at positions $x_1$ and $x_2$ on a line, with energy
    $E = \tfrac{1}{2}k\,(x_2 - x_1 - d_0)^2$. The force on atom 2 is
    $-k\,(x_2 - x_1 - d_0)$: a spring pulling the bond towards its
    natural length $d_0$. Every richer model reduces to this near a
    minimum.
```

---

### Worked example

*Purpose:* one small numerical or symbolic example carried through to a number, so the reader sees the machinery turn.

```markdown
!!! example "Worked example"
    With $k = 2\ \mathrm{eV\,\AA^{-2}}$, $d_0 = 1.5\ \mathrm{\AA}$ and
    a bond stretched to $x_2 - x_1 = 1.7\ \mathrm{\AA}$, the force is
    $-2 \times (1.7 - 1.5) = -0.4\ \mathrm{eV\,\AA^{-1}}$ — a restoring
    pull of $0.4\ \mathrm{eV\,\AA^{-1}}$.
```

---

### Code connection

*Purpose:* show where the idea lives in the codebase, and explain the snippet in prose. State the runtime honesty: only NumPy / SciPy / Matplotlib runs live in the browser; ASE, pymatgen, LAMMPS and Quantum ESPRESSO do not.

```markdown
!!! info "Code connection"
    The spring force above is one line of NumPy. This runs live in the
    browser:

    ```python
    import numpy as np
    k, d0 = 2.0, 1.5
    x1, x2 = 0.0, 1.7
    force_on_2 = -k * (x2 - x1 - d0)   # eV / Å
    print(force_on_2)                  # -0.4
    ```

    The full handbook code lives under `code/tier1/`; see the
    [code reading guide](code-reading-guide.md) for how to read it.
```

---

### Common misunderstandings

*Purpose:* list the specific beginner mistakes this section invites, so the reader can recognise and avoid them.

```markdown
!!! warning "Common misunderstandings"
    - The density $n(\mathbf{r})$ is **not** a probability for one
      electron; it integrates to the *total* number of electrons.
    - "Self-consistent" does **not** mean "exact" — the functional is
      still approximate.
```

---

### Check yourself

*Purpose:* three to five short questions that test whether the reader followed the scaffolding. Use the recommended admonition form below.

```markdown
!!! question "Check yourself"
    1. What does the symbol $n(\mathbf{r})$ represent, and what are its units?
    2. Why is the force the *negative* gradient of the energy?
    3. In the two-atom toy model, what happens to the force when the
       bond sits exactly at $d_0$?
```

---

### Short answers or hints

*Purpose:* collapsible answers and hints that support self-study without giving everything away on the page. Keep them folded so the reader chooses when to look.

```markdown
??? success "Answer"
    1. The electron density at position $\mathbf{r}$, in electrons per
       $a_0^3$; it integrates to the total electron count.
    2. The system lowers its energy by moving downhill, so the force
       points down the energy gradient.
    3. The force is zero — the bond is at equilibrium.

??? note "Hint"
    Look again at the step-by-step block: the force is a *gradient*.
```

---

### Where this appears later

*Purpose:* forward links so the reader knows where this idea is used again. Link only to a chapter's overview index, using the exact folder names; write specific numbered sections in prose, never as links.

```markdown
!!! tip "Where this appears later"
    These forces drive every step of molecular dynamics in
    [Chapter 7 (MD)](../ch07-md/index.md), and the same energy
    functional underpins [Chapter 5 (DFT)](../ch05-dft/index.md)
    (see Section 5.3 there).
```

---

## Reviewer checklist

Before approving an in-place expansion, confirm:

- [ ] **No symbol undefined** — every symbol near new text appears in a symbol guide or the [beginner glossary](glossary-for-beginners.md).
- [ ] **No hard equation without an intro** — each difficult equation has a plain-language or physical-picture lead-in.
- [ ] **No code block without explanation** — every snippet is explained in prose, with browser-runnability stated honestly.
- [ ] **At least one self-check per page** — a `!!! question "Check yourself"` block, ideally with folded answers.
- [ ] **Scientific accuracy preserved** — nothing in the original was weakened, contradicted or dumbed down.
- [ ] **British English** throughout, no emoji, no hype words.
- [ ] **Original material intact** — the author's text, equations, figures and code are unchanged (unless clearly wrong), and the diff is small and reviewable.

For the broader plan and what to expand first, see the [expansion roadmap](expansion-roadmap-v1.4.md) and the [undergraduate projects](undergraduate-projects.md) page.
