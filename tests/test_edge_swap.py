"""Tests for fvcom_mesh_tools.algorithms.edge_swap."""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.algorithms import (
    alpha_quality,
    min_interior_angle,
    signed_areas,
    swap_edges_for_quality,
)
from fvcom_mesh_tools.io import Fort14Mesh


def _bad_diagonal_quad() -> Fort14Mesh:
    """Convex quadrilateral split along the wrong diagonal.

    The four nodes form a thin "kite" (1, 0)-(2, 0.5)-(1, 1)-(0, 0.5).
    Splitting along the long diagonal (0, 2) produces two reasonable
    triangles. The mesh is constructed instead with the short diagonal
    (1, 3), producing two slivers - exactly the case the swap is meant
    to fix.
    """
    nodes = np.array(
        [
            [1.0, 0.0],   # 0
            [2.0, 0.5],   # 1
            [1.0, 1.0],   # 2
            [0.0, 0.5],   # 3
        ],
        dtype=np.float64,
    )
    # Bad diagonal: (1, 3). Two triangles share edge (1, 3).
    elements = np.array(
        [
            [0, 1, 3],
            [1, 2, 3],
        ],
        dtype=np.int64,
    )
    # Outer ring as one land segment; no open boundary.
    land = np.array([0, 1, 2, 3], dtype=np.int64)
    return Fort14Mesh(
        title="bad-diagonal",
        nodes=nodes,
        depths=np.zeros(4),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, land)],
    )


def test_swap_fixes_bad_diagonal() -> None:
    before = _bad_diagonal_quad()
    after, info = swap_edges_for_quality(before, max_iters=5)
    # Exactly one swap should resolve this.
    assert info["total_swaps"] == 1
    # Resulting triangles must share the long diagonal (0, 2).
    pairs = {tuple(sorted(t)) for t in after.elements.tolist()}
    assert pairs == {(0, 1, 2), (0, 2, 3)}


def test_swap_improves_min_angle_distribution() -> None:
    before = _bad_diagonal_quad()
    a_before = min_interior_angle(before).min()
    q_before = alpha_quality(before).min()
    after, _ = swap_edges_for_quality(before, max_iters=5)
    a_after = min_interior_angle(after).min()
    q_after = alpha_quality(after).min()
    assert a_after > a_before
    assert q_after > q_before


def test_swap_no_change_when_already_optimal() -> None:
    """Equilateral pair: any swap would create one obtuse triangle, so
    the algorithm must leave it alone."""
    s = np.sqrt(3.0) / 2.0
    nodes = np.array(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, s],
            [1.5, s],
        ],
        dtype=np.float64,
    )
    elements = np.array([[0, 1, 2], [1, 3, 2]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="eq",
        nodes=nodes,
        depths=np.zeros(4),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, np.array([0, 1, 3, 2], dtype=np.int64))],
    )
    after, info = swap_edges_for_quality(mesh, max_iters=5)
    assert info["total_swaps"] == 0
    np.testing.assert_array_equal(after.elements, mesh.elements)


def test_swap_does_not_flip_triangles() -> None:
    before = _bad_diagonal_quad()
    after, _ = swap_edges_for_quality(before, max_iters=5)
    assert (signed_areas(after) > 0).all()


def test_swap_preserves_node_set() -> None:
    """Swapping must not introduce or remove nodes from elements."""
    before = _bad_diagonal_quad()
    nodes_before = set(before.elements.ravel().tolist())
    after, _ = swap_edges_for_quality(before, max_iters=5)
    nodes_after = set(after.elements.ravel().tolist())
    assert nodes_before == nodes_after


def test_swap_mapping_with_many_triangles() -> None:
    """Regression: the edge -> triangle mapping uses ``k % NE``, not
    ``k // 3``. With more than 3 triangles the two formulas diverge,
    so this case catches a bug that the small test fixtures cannot."""
    # 6-triangle mesh: a 2x2 quadrilateral grid split into two triangles
    # per cell, with one cell deliberately mis-split so a swap exists.
    nodes = np.array(
        [
            [0.0, 0.0], [1.0, 0.0], [2.0, 0.0],
            [0.0, 1.0], [1.0, 1.0], [2.0, 1.0],
            [0.0, 2.0], [1.0, 2.0], [2.0, 2.0],
        ],
        dtype=np.float64,
    )
    elements = np.array(
        [
            # Lower-left cell: standard split (0, 1, 4) (0, 4, 3).
            [0, 1, 4],
            [0, 4, 3],
            # Lower-right cell: also standard (1, 2, 5) (1, 5, 4).
            [1, 2, 5],
            [1, 5, 4],
            # Upper-left cell: standard (3, 4, 7) (3, 7, 6).
            [3, 4, 7],
            [3, 7, 6],
            # Upper-right cell: deliberately mis-split, so a swap exists
            # that would convert the diagonal (4, 8) to (5, 7).
            [4, 5, 8],
            [4, 8, 7],
        ],
        dtype=np.int64,
    )
    land = np.array([0, 1, 2, 5, 8, 7, 6, 3], dtype=np.int64)
    mesh = Fort14Mesh(
        title="grid",
        nodes=nodes,
        depths=np.zeros(9),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, land)],
    )
    after, info = swap_edges_for_quality(mesh, max_iters=5)
    # Bug regression: this used to crash with IndexError on the third
    # block of edges because of a wrong row->triangle mapping.
    assert info["total_swaps"] >= 0
    assert (signed_areas(after) > 0).all()
    # Node set unchanged.
    assert set(after.elements.ravel().tolist()) == set(elements.ravel().tolist())
