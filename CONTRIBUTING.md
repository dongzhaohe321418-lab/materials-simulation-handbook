# Contributing

Pull requests and issues are very welcome. A few notes to make things go smoothly.

## Reporting an issue

If you spot a typo, an incorrect derivation, broken code, or a confusing passage, please open an issue describing:

- which chapter / page / section
- what you expected vs. what you found
- a permalink to the line if possible

Small, sharp bug reports are the most useful contribution you can make.

## Proposing a change

For typos and trivial fixes, just open a pull request directly. For larger changes (a rewritten section, a new chapter, a new worked example):

1. Open an issue first to discuss the scope.
2. Fork the repository and create a branch named `chXX-short-description`.
3. Make your changes. Keep prose changes and code changes in separate commits where reasonable.
4. Run `mkdocs serve` locally and verify nothing breaks.
5. Open the pull request, referencing the issue.

## Style

- British English (colour, behaviour, optimise).
- Definitions before derivations; derivations before code; code before exercises.
- No emoji in source files. (The site theme handles iconography.)
- For mathematics, use the existing notation — see [How to use this book](docs/how-to-use.md).
- Code: PEP 8, type hints on public functions, NumPy-style docstrings.

## Code of conduct

Be kind. Argue about the physics, never about the person. See [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
