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


def test_classify_outer_loop_merges_short_coast_gap() -> None:
    """An open-open-land-open-open ring with a 1-node land run between two
    open runs should merge into a single open segment when
    open_merge_coast_gap >= 1."""
    # 6-node ring: 5 nodes near the bbox top edge, 1 node intruding.
    nodes = np.array(
        [
            [-1.0, 1.0],   # 0 open
            [-0.5, 1.0],   # 1 open
            [0.0, 0.5],    # 2 land (intrusion)
            [0.5, 1.0],    # 3 open
            [1.0, 1.0],    # 4 open
            [0.0, -1.0],   # 5 land (south)
        ],
        dtype=np.float64,
    )
    outer = np.array([0, 1, 2, 3, 4, 5, 0], dtype=np.int64)
    bbox = (-1.5, 1.0, 1.5, 1.0)

    # Without merging: 2 open segments.
    open_no_merge, _ = classify_outer_loop_by_bbox(
        outer, nodes, bbox=bbox, tol=0.05,
    )
    assert len(open_no_merge) == 2

    # With merging gap >= 1: a single open segment containing all five
    # bbox-touching nodes plus the bridged intrusion.
    open_merged, _ = classify_outer_loop_by_bbox(
        outer, nodes, bbox=bbox, tol=0.05, open_merge_coast_gap=1,
    )
    assert len(open_merged) == 1
    bridged_set = set(open_merged[0].tolist())
    assert {0, 1, 2, 3, 4}.issubset(bridged_set)


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


# ---------------------------------------------------------------------------
# boundary_snap (exact conformity)
# ---------------------------------------------------------------------------


def test_snap_boundary_to_polylines_snaps_smooth_and_caps_far():
    from shapely.geometry import LineString

    from fvcom_mesh_tools.algorithms.boundary_snap import (
        snap_boundary_to_polylines,
    )

    # 3x3-node grid, 1000 m cells; the "true coastline" runs 80 m
    # north of the bottom edge (smooth -> snappable) and one far line
    # 5 km away that must be ignored by the cap.
    n = 3
    nodes = np.array([[i * 1000.0, j * 1000.0]
                      for j in range(n) for i in range(n)])
    elements = []
    for j in range(n - 1):
        for i in range(n - 1):
            a, b = j * n + i, j * n + i + 1
            c, d = (j + 1) * n + i + 1, (j + 1) * n + i
            elements.append([a, b, c])
            elements.append([a, c, d])
    mesh = Fort14Mesh(
        title="snap",
        nodes=nodes,
        depths=np.full(n * n, 5.0),
        elements=np.asarray(elements),
        open_boundaries=[np.array([5, 8])],  # right side: excluded
        land_boundaries=[(20, np.array([8, 7, 6, 3, 0, 1, 2, 5]))],
    )
    coast = LineString([(-500.0, 80.0), (2500.0, 80.0)])
    far = LineString([(-500.0, -5000.0), (2500.0, -5000.0)])
    out, info = snap_boundary_to_polylines(mesh, [coast, far])

    # Bottom-row land nodes 0,1,2 sit 80 m from the coast (< 0.6*h)
    # -> snapped exactly onto y=80.
    assert np.allclose(out.nodes[[0, 1, 2], 1], 80.0)
    # OBC nodes untouched.
    assert np.allclose(out.nodes[[5, 8]], mesh.nodes[[5, 8]])
    assert info["n_snapped"] >= 3
    # The toy grid's top/side nodes have no nearby reference line, so
    # global percentiles only need to not regress.
    assert info["dist_after_p50_m"] <= info["dist_before_p50_m"] + 1e-9
    # Input untouched.
    assert np.allclose(mesh.nodes[0], [0.0, 0.0])


def test_snap_nodes_to_segment_projects_and_clamps():
    from fvcom_mesh_tools.algorithms.boundary_snap import snap_nodes_to_segment

    nodes = np.array([
        [0.0, 30.0], [1000.0, -40.0], [2000.0, 25.0],
        [0.0, 1000.0], [1000.0, 1000.0], [2000.0, 1000.0],
    ])
    elements = np.array([[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]])
    mesh = Fort14Mesh(
        title="seg",
        nodes=nodes,
        depths=np.full(6, 5.0),
        elements=elements,
        open_boundaries=[np.array([0, 1, 2])],
        land_boundaries=[(20, np.array([2, 5, 4, 3, 0]))],
    )
    out, info = snap_nodes_to_segment(
        mesh, [0, 1, 2], (0.0, 0.0), (2000.0, 0.0),
    )
    assert info["n_snapped"] == 3
    assert np.allclose(out.nodes[[0, 1, 2], 1], 0.0)
    assert np.allclose(out.nodes[[0, 1, 2], 0], [0.0, 1000.0, 2000.0])
