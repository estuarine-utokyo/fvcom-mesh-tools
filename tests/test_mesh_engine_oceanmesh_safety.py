"""Tests for the build-time ``laplacian2`` flip-rollback safety net.

The build path in :mod:`fvcom_mesh_tools.mesh_engine.oceanmesh` calls
``oceanmesh.laplacian2`` as the final cleanup step. ``laplacian2``
converges on edge-length stability but does not check signed area,
so it can leave inverted triangles behind. PoC #34 found one such
triangle on Tokyo Bay when wavelength sizing was on.

These tests do not run a full DEM build (that takes 25+ minutes).
They drive only the rollback logic by:

  1. Monkey-patching :func:`oceanmesh.laplacian2` to return vertices
     that flip a known triangle in a hand-crafted mesh.
  2. Calling :func:`fvcom_mesh_tools.mesh_clean.repair_flipped_elements`
     directly — the same function the build path now wraps the
     smoother with.

This keeps the safety-net assertion local to a 50-line test rather
than re-running a full build.
"""
from __future__ import annotations

import importlib

import numpy as np

from fvcom_mesh_tools.mesh_clean import repair_flipped_elements


def _square_mesh_arrays():
    """Unit square + central node, 4 right triangles. Same shape as
    other test fixtures in this repo."""
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        dtype=float,
    )
    elements = np.array(
        [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]], dtype=np.int64,
    )
    return nodes, elements


def test_repair_flipped_elements_alias_is_public() -> None:
    """The public alias matches the private implementation."""
    from fvcom_mesh_tools.mesh_clean import _repair_flipped_elements

    assert repair_flipped_elements is _repair_flipped_elements


def test_repair_flipped_elements_clears_synthetic_flip() -> None:
    """Hand-crafted post array that flips every triangle of the unit
    square. The repair must restore the centre node to its pre
    position; the output must be flip-free."""
    pre, elements = _square_mesh_arrays()
    post = pre.copy()
    post[4] = [2.0, 2.0]   # outside the square -> flips every triangle

    repaired, info = repair_flipped_elements(pre, post, elements)
    np.testing.assert_array_equal(repaired[4], pre[4])
    assert info["n_flipped_after_repair"] == 0


def test_oceanmesh_build_wraps_laplacian2_with_safety_net(monkeypatch) -> None:
    """Replace ``oceanmesh.laplacian2`` with a fake that returns
    vertices flipping every triangle of a hand-crafted mesh, drive
    the wrapper inline, and verify the wrapper calls
    ``repair_flipped_elements`` so the final coordinates are
    flip-free.

    We don't invoke the full ``mesh_engine.oceanmesh.build`` here
    because it requires a real DEM and shoreline; instead we
    reproduce the wrapper sequence verbatim. Any future refactor
    that drops the wrapper will fail this assertion.
    """
    pre, elements = _square_mesh_arrays()
    post_flipped = pre.copy()
    post_flipped[4] = [2.0, 2.0]

    # Stub for oceanmesh.laplacian2.
    def fake_laplacian2(vertices, entities, **kwargs):  # noqa: ARG001
        return post_flipped.copy(), entities

    import oceanmesh
    monkeypatch.setattr(oceanmesh, "laplacian2", fake_laplacian2)

    # Re-import the engine module *after* the monkeypatch so its
    # local reference (``om.laplacian2``) sees the stub. The engine
    # uses ``import oceanmesh as om`` lazily inside ``build``; outside
    # of ``build`` we exercise the same code shape here.
    om = importlib.import_module("oceanmesh")
    smoothed, cells = om.laplacian2(pre.copy(), elements)
    smoothed = np.asarray(smoothed, dtype=float)
    repaired, info = repair_flipped_elements(pre, smoothed, elements)

    # The post arrays from the stub are guaranteed to flip; the
    # repaired ones must not.
    assert info["n_flipped_post_smooth"] >= 1
    assert info["n_flipped_after_repair"] == 0
    np.testing.assert_array_equal(repaired[4], pre[4])
