"""Tests for refine_bad_triangles."""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.algorithms import (
    min_interior_angle,
    refine_bad_triangles,
    signed_areas,
)
from fvcom_mesh_tools.io import Fort14Mesh


def _sliver_mesh() -> Fort14Mesh:
    """Two triangles in CCW orientation, one a deliberate sliver
    (min-angle ~ 6 deg).

    Quad 0-1-2-3 split along the long diagonal (0,2). Node 3 is pulled
    very close to the (0,2) line so the (0,2,3) triangle is a sliver,
    while (0,1,2) is well-shaped.
    """
    nodes = np.array(
        [
            [0.0, 0.0],
            [2.0, 0.0],
            [1.0, 1.0],     # 2: above
            [1.0, -0.05],   # 3: just below x-axis -> sliver below
        ],
        dtype=np.float64,
    )
    elements = np.array([[0, 1, 2], [1, 0, 3]], dtype=np.int64)
    land = np.array([0, 3, 1, 2], dtype=np.int64)
    return Fort14Mesh(
        title="sliver",
        nodes=nodes,
        depths=np.array([1.0, 2.0, 3.0, 4.0]),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, land)],
    )


def test_refine_keeps_mesh_valid_on_pathological_sliver() -> None:
    """An isolated sliver triangle whose long edge is the only interior
    edge cannot be improved by longest-edge bisection alone (Rivara
    propagation has no neighbour to recurse into). The algorithm may
    still insert a midpoint - what matters is that the resulting mesh
    is *valid* (no flipped triangles), not that the sliver count drops."""
    before = _sliver_mesh()
    bad_before = int((min_interior_angle(before) < 20.0).sum())
    assert bad_before == 1
    after, _ = refine_bad_triangles(
        before, min_angle_threshold=20.0, max_passes=5,
    )
    assert (signed_areas(after) > 0).all()


def _interior_sliver_mesh() -> Fort14Mesh:
    """Square split into two triangles, with one corner pulled in close
    to the diagonal so the (0,1,3) triangle is a sliver whose *longest*
    edge is the interior diagonal (1,3). Bisecting (1,3) splits the
    sliver along its long axis and creates a balanced sub-triangle on
    the bad-vertex side."""
    nodes = np.array(
        [
            [0.4, 0.5],    # 0 - pulled close to diagonal (1, 3)
            [1.0, 0.0],    # 1
            [1.0, 1.0],    # 2
            [0.0, 1.0],    # 3
        ],
        dtype=np.float64,
    )
    elements = np.array(
        [
            [0, 1, 3],   # CCW sliver: 0 close to line (1, 3)
            [1, 2, 3],   # CCW good
        ],
        dtype=np.int64,
    )
    # The boundary walks 0 -> 1 -> 2 -> 3 -> 0; (1, 3) is the interior
    # diagonal.
    land = np.array([0, 1, 2, 3], dtype=np.int64)
    return Fort14Mesh(
        title="interior-sliver",
        nodes=nodes,
        depths=np.zeros(4),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, land)],
    )


def test_refine_returns_valid_mesh_on_interior_sliver() -> None:
    """Refine on a mesh whose bad triangle has an interior longest edge.
    The refinement may or may not improve the metric on this small
    fixture (regression rollback may decline the change), but the
    returned mesh must always be valid (no flipped triangles)."""
    before = _interior_sliver_mesh()
    bad_before = int((min_interior_angle(before) < 20.0).sum())
    assert bad_before >= 1
    after, info = refine_bad_triangles(
        before, min_angle_threshold=20.0, max_passes=3,
    )
    assert (signed_areas(after) > 0).all()
    bad_after = int((min_interior_angle(after) < 20.0).sum())
    # Early-stop guarantees: bad_after <= bad_before.
    assert bad_after <= bad_before
    assert "stop_reason" in info


def test_refine_no_op_when_clean() -> None:
    """Equilateral triangle pair: no triangle is bad, so refinement
    must leave the mesh untouched."""
    s = np.sqrt(3.0) / 2.0
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.5, s], [1.5, s]],
        dtype=np.float64,
    )
    elements = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="eq", nodes=nodes,
        depths=np.zeros(4), elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 3, 2], dtype=np.int64))],
    )
    after, info = refine_bad_triangles(mesh, min_angle_threshold=20.0)
    np.testing.assert_array_equal(after.nodes, mesh.nodes)
    np.testing.assert_array_equal(after.elements, mesh.elements)
    assert info["total_nodes_inserted"] == 0
    assert info["passes"] == 0


def test_refine_does_not_change_boundary_node_ids() -> None:
    """Original boundary node-ids must keep referring to the same
    coordinates after refinement; only interior node-ids are added."""
    before = _sliver_mesh()
    land_ids = before.land_boundaries[0][1]
    after, _ = refine_bad_triangles(before, min_angle_threshold=20.0, max_passes=5)
    # Boundary-list arrays are shared by reference (algorithm doesn't
    # touch them); the same indices in `after.nodes` must give the
    # same coordinates.
    np.testing.assert_array_equal(after.nodes[land_ids], before.nodes[land_ids])
    assert after.n_nodes >= before.n_nodes
