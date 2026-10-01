"""Shared test configuration.

``@pytest.mark.needs_oceanmesh`` marks a test that exercises the laboratory's
oceanmesh fork (GPL-3.0; not installable from PyPI or conda-forge, see
docs/USER_GUIDE.md §2). Where it is not installed -- GitHub's CI -- those
tests are skipped with that reason instead of failing on the import.
"""

from __future__ import annotations

import importlib

import pytest


def _importable(name: str) -> bool:
    try:
        importlib.import_module(name)
    except ImportError:
        return False
    return True


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "needs_oceanmesh: needs the laboratory's oceanmesh fork installed")


def pytest_collection_modifyitems(config, items):
    # the optional backend is imported only when a selected test needs it:
    # its import has side effects (a matplotlib cache) that a pure selection
    # should not depend on (review of the extend tools, round 14 F8)
    marked = [item for item in items if "needs_oceanmesh" in item.keywords]
    if not marked or _importable("oceanmesh"):
        return
    skip = pytest.mark.skip(reason="needs the laboratory's oceanmesh fork, not installed")
    for item in marked:
        item.add_marker(skip)
