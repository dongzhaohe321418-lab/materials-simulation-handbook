"""Shared pytest fixtures and import-path setup for the tier1 test suite."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# Make the `code/tier1` package importable as the top-level `tier1` package.
_REPO_ROOT = Path(__file__).resolve().parent.parent
_CODE_DIR = _REPO_ROOT / "code"
if str(_CODE_DIR) not in sys.path:
    sys.path.insert(0, str(_CODE_DIR))


@pytest.fixture(scope="session")
def rng() -> np.random.Generator:
    return np.random.default_rng(seed=12345)


@pytest.fixture(scope="session")
def small_grid_size() -> int:
    return 64


@pytest.fixture(scope="session")
def box_length_si() -> float:
    return 1.0e-9
