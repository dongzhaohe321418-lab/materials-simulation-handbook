#!/usr/bin/env python3
"""Verify all key dependencies import cleanly and report versions.

Run this after `conda env create -f environment.yml && conda activate matsim`.
A failure is informative — it tells you which optional dependency is missing,
and the worked examples in the corresponding chapter will still degrade
gracefully (e.g., Materials Project examples fall back to cached data when
`mp-api` is missing or no API key is set).
"""

from __future__ import annotations

import importlib
import os
import sys
from typing import Optional


REQUIRED: list[tuple[str, str]] = [
    ("numpy", "Ch 0, everywhere"),
    ("scipy", "Ch 0, Ch 4"),
    ("matplotlib", "Plotting throughout"),
    ("ase", "Ch 3, Ch 6, Ch 7"),
    ("pymatgen", "Ch 3, Ch 10"),
]

ML_OPTIONAL: list[tuple[str, str]] = [
    ("torch", "Ch 9, Ch 10, Ch 11"),
    ("torch_geometric", "Ch 10 (GNN)"),
    ("mace", "Ch 9 (MLIP training), Ch 12 (foundation model)"),
    ("gpytorch", "Ch 11 (GP)"),
    ("botorch", "Ch 11 (BO)"),
    ("mp_api", "Ch 10 (Materials Project queries) — needs API key"),
]


def try_import(name: str) -> Optional[str]:
    """Return version string if importable, None otherwise."""
    try:
        mod = importlib.import_module(name)
        return getattr(mod, "__version__", "unknown")
    except Exception as exc:
        return f"FAILED: {type(exc).__name__}: {exc}"


def main() -> int:
    print("=" * 60)
    print("Materials Simulation Handbook — environment check")
    print("=" * 60)
    print(f"Python: {sys.version.split()[0]}")
    print()

    print("Required packages")
    print("-" * 60)
    fail_required = 0
    for pkg, where in REQUIRED:
        v = try_import(pkg)
        marker = "OK" if v and not v.startswith("FAILED") else "FAIL"
        print(f"  [{marker}] {pkg:<20s} {v}  ({where})")
        if marker == "FAIL":
            fail_required += 1

    print()
    print("Optional ML / database packages (chapters degrade if missing)")
    print("-" * 60)
    for pkg, where in ML_OPTIONAL:
        v = try_import(pkg)
        marker = "OK" if v and not v.startswith("FAILED") else "SKIP"
        print(f"  [{marker}] {pkg:<20s} {v}  ({where})")

    print()
    print("Materials Project API key")
    print("-" * 60)
    if os.environ.get("MP_API_KEY"):
        print("  [OK] MP_API_KEY is set — live queries enabled.")
    else:
        print("  [SKIP] MP_API_KEY not set — Ch 10 falls back to cached data.")
        print("         Get a free key at https://next-gen.materialsproject.org/api")

    print()
    if fail_required:
        print(f"FAILED: {fail_required} required package(s) missing.")
        return 1
    print("All required packages present.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
