"""Tests for the Laplacian smoother and quality metrics."""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.algorithms import (
    alpha_quality,
    laplacian_smooth,
    min_interior_angle,
    signed_areas,
)
from fvcom_mesh_tools.io import Fort14Mesh


def _sliver_mesh() -> Fort14Mesh:
    """3x3 grid of points with a deliberately offset interior node.

    Outer 8 nodes form a unit square boundary; node 4 is the central
    interior node, pulled toward one corner so the surrounding triangles
    are slivers. Boundaries cover the entire outer ring as land. After
    smoothing, node 4 should drift back toward the centroid.
    """
    nodes = np.array(
        [
            [0.0, 0.0],   # 0  (corner)
            [0.5, 0.0],   # 1
            [1.0, 0.0],   # 2  (corner)
            [0.0, 0.5],   # 3
            [0.10, 0.10], # 4  interior, intentionally off-centre
            [1.0, 0.5],   # 5
            [0.0, 1.0],   # 6  (corner)
            [0.5, 1.0],   # 7
            [1.0, 1.0],   # 8  (corner)
        ],
        dtype=np.float64,
    )
    elements = np.array(
        [
            [0, 1, 4],
            [1, 2, 4],
            [2, 5, 4],
            [5, 8, 4],
            [8, 7, 4],
            [7, 6, 4],
            [6, 3, 4],
            [3, 0, 4],
        ],
        dtype=np.int64,
    )
    # Whole outer ring as a single land segment.
    land = np.array([0, 1, 2, 5, 8, 7, 6, 3], dtype=np.int64)
    return Fort14Mesh(
        title="sliver",
        nodes=nodes,
        depths=np.zeros(9),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[(0, land)],
    )


def test_laplacian_smooth_improves_quality_on_sliver() -> None:
    before = _sliver_mesh()
    q_before = alpha_quality(before).mean()
    after, info = laplacian_smooth(before, n_iters=20, alpha=0.5)
    q_after = alpha_quality(after).mean()
    assert q_after > q_before
    assert info["degenerate_remaining"] == 0


def test_laplacian_smooth_keeps_boundary_fixed() -> None:
    before = _sliver_mesh()
    land = before.land_boundaries[0][1]
    after, _ = laplacian_smooth(before, n_iters=5, alpha=0.5)
    np.testing.assert_array_equal(after.nodes[land], before.nodes[land])


def test_laplacian_smooth_no_op_when_alpha_tiny() -> None:
    before = _sliver_mesh()
    after, info = laplacian_smooth(before, n_iters=1, alpha=1e-12)
    np.testing.assert_allclose(after.nodes, before.nodes, atol=1e-10)
    # Still records 1 iteration of work.
    assert info["n_iters"] == 1


def test_laplacian_smooth_does_not_flip_triangles() -> None:
    before = _sliver_mesh()
    after, info = laplacian_smooth(before, n_iters=20, alpha=1.0, prevent_flips=True)
    assert (signed_areas(after) > 0).all()
    assert info["degenerate_remaining"] == 0


def test_laplacian_smooth_central_node_drifts_to_centroid() -> None:
    """alpha=1, many iterations: node 4 (the only movable node) should
    converge to the mean of its 8 neighbours (the outer ring centroid)."""
    before = _sliver_mesh()
    # Outer ring is a unit square, centroid (0.5, 0.5).
    after, _ = laplacian_smooth(before, n_iters=50, alpha=1.0, prevent_flips=False)
    np.testing.assert_allclose(after.nodes[4], (0.5, 0.5), atol=1e-9)


def test_min_interior_angle_equilateral_is_60_deg() -> None:
    nodes = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3.0) / 2.0]],
        dtype=np.float64,
    )
    elements = np.array([[0, 1, 2]], dtype=np.int64)
    mesh = Fort14Mesh(
        title="eq",
        nodes=nodes,
        depths=np.zeros(3),
        elements=elements,
        open_boundaries=[],
        land_boundaries=[],
    )
    np.testing.assert_allclose(min_interior_angle(mesh)[0], 60.0, atol=1e-9)
    np.testing.assert_allclose(alpha_quality(mesh)[0], 1.0, atol=1e-9)
