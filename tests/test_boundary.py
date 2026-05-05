"""Tests for fvcom_mesh_tools.algorithms.boundary."""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.algorithms import (
    boundary_edges_from_tris,
    chain_edges_to_loops,
    classify_boundaries_by_bbox,
    classify_outer_loop_by_bbox,
    outer_loop,
)
from fvcom_mesh_tools.io import Fort14Mesh


def _square_with_hole_mesh() -> Fort14Mesh:
    """Annulus mesh: 4 outer corners, 4 inner corners, 8 triangles.

    Outer ring: nodes 0..3 at +/-1 box. Inner ring (hole): nodes 4..7 at
    a smaller box centred on the origin. Eight triangles bridge inner
    and outer rings.
    """
    nodes = np.array(
        [
            [-1.0, -1.0],  # 0
            [1.0, -1.0],   # 1
            [1.0, 1.0],    # 2
            [-1.0, 1.0],   # 3
            [-0.3, -0.3],  # 4
            [0.3, -0.3],   # 5
            [0.3, 0.3],    # 6
            [-0.3, 0.3],   # 7
        ],
        dtype=np.float64,
    )
    elements = np.array(
        [
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
        dtype=np.int64,
    )
    return Fort14Mesh(
        title="annulus",
        nodes=nodes,
        depths=np.zeros(8),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[],
    )


def test_boundary_edges_two_loops() -> None:
    mesh = _square_with_hole_mesh()
    edges = boundary_edges_from_tris(mesh.elements)
    # Two rings of 4 edges each.
    assert edges.shape == (8, 2)


def test_chain_edges_to_loops_returns_two_closed_loops() -> None:
    mesh = _square_with_hole_mesh()
    edges = boundary_edges_from_tris(mesh.elements)
    loops = chain_edges_to_loops(edges)
    assert len(loops) == 2
    for loop in loops:
        assert loop[0] == loop[-1]
        # 4 distinct nodes + closing duplicate.
        assert loop.size == 5


def test_outer_loop_picks_largest_ring() -> None:
    mesh = _square_with_hole_mesh()
    loops = chain_edges_to_loops(boundary_edges_from_tris(mesh.elements))
    outer = outer_loop(loops, mesh.nodes)
    # Outer ring contains nodes 0..3 (in some rotation) plus closing dup.
    assert set(outer.tolist()) == {0, 1, 2, 3}


def test_classify_outer_loop_by_bbox_splits_open_and_land() -> None:
    """Loop along the +/-1 square: clip the bbox to the +y side so only the
    top edge is "open" and bottom + the two sides are "land"."""
    nodes = np.array(
        [
            [-1.0, -1.0],
            [1.0, -1.0],
            [1.0, 1.0],
            [-1.0, 1.0],
        ],
        dtype=np.float64,
    )
    outer = np.array([0, 1, 2, 3, 0], dtype=np.int64)
    # bbox covers only the top edge (y=1) ; tol=1e-6 keeps it strict.
    bbox = (-1.0, 0.99, 1.0, 1.01)
    open_segs, land_segs = classify_outer_loop_by_bbox(
        outer, nodes, bbox=bbox, tol=1e-6,
    )
    assert len(open_segs) == 1
    assert len(land_segs) == 1
    # Open segment should contain the two top nodes (2, 3).
    assert {2, 3}.issubset(set(open_segs[0].tolist()))
    # Adjacent segments share an endpoint.
    assert open_segs[0][0] in land_segs[0].tolist()
    assert open_segs[0][-1] in land_segs[0].tolist()


def test_classify_outer_loop_handles_all_open() -> None:
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]],
        dtype=np.float64,
    )
    outer = np.array([0, 1, 2, 3, 0], dtype=np.int64)
    open_segs, land_segs = classify_outer_loop_by_bbox(
        outer, nodes, bbox=(0.0, 0.0, 1.0, 1.0), tol=1e-6,
    )
    assert len(open_segs) == 1
    assert open_segs[0].size == 4
    assert land_segs == []


def test_classify_boundaries_by_bbox_flags_island_as_land() -> None:
    """Annulus mesh with a generous bbox covering only one side of the
    outer ring: the outer ring is half open / half land, and the hole
    ring is reported as a closed land segment."""
    mesh = _square_with_hole_mesh()
    # bbox covers the +x side of the outer ring.
    bbox = (0.99, -1.0, 1.01, 1.0)
    open_segs, land_bnds = classify_boundaries_by_bbox(
        mesh, bbox=bbox, tol=1e-6, land_ibtype=20,
    )
    assert len(open_segs) >= 1
    assert any(20 == ib for ib, _ in land_bnds)
    # At least one of the land segments must consist purely of inner-ring
    # nodes (the island).
    inner = {4, 5, 6, 7}
    assert any(set(seg.tolist()).issubset(inner) for _, seg in land_bnds)


def test_classify_boundaries_by_bbox_no_loops_returns_empty() -> None:
    """Degenerate single-triangle mesh: every edge is a boundary edge,
    forming exactly one loop. The bbox covers the whole triangle, so we
    expect a single open segment and no land segments."""
    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    elements = np.array([[0, 1, 2]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="tri", nodes=nodes, depths=np.zeros(3), elements=elements,
        open_boundaries=[], land_boundaries=[],
    )
    open_segs, land_bnds = classify_boundaries_by_bbox(
        mesh, bbox=(0.0, 0.0, 1.0, 1.0), tol=1e-6,
    )
    assert len(open_segs) == 1
    assert land_bnds == []
