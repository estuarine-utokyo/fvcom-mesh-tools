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
    config.addinivalue_line(
        "markers", "needs_ocsmesh: needs OCSMesh, an optional private-use backend "
                   "outside the default environment (THIRD_PARTY_NOTICES.md)")


@pytest.hookimpl(hookwrapper=True, trylast=True)
def pytest_collection_modifyitems(config, items):
    # The optional backend is imported only when a test that survived
    # selection (-k, -m) needs it: its import has side effects (a matplotlib
    # cache) that a pure selection should not depend on (review of the
    # extend tools, rounds 14 F8, 15 F5). As a wrapper, this runs after the
    # selection hooks have removed the deselected items.
    yield
    for marker, module, reason in (
            ("needs_oceanmesh", "oceanmesh", "needs the laboratory's oceanmesh fork, not "
                                             "installed"),
            ("needs_ocsmesh", "ocsmesh", "needs OCSMesh (optional, private use), not "
                                         "installed")):
        marked = [item for item in items if marker in item.keywords]
        if not marked or _importable(module):
            continue
        skip = pytest.mark.skip(reason=reason)
        for item in marked:
            item.add_marker(skip)
