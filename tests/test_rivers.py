"""Tests for the river-inflow loader and segment-splitting algorithm."""

from __future__ import annotations

import numpy as np
import pytest

from fvcom_mesh_tools.algorithms import add_river_inflow_segments
from fvcom_mesh_tools.io import Fort14Mesh, load_river_points


def _coast_mesh() -> Fort14Mesh:
    """Five-node land segment along a 1-degree coast at y=0; one open
    segment on the south. The land segment runs left-to-right and is a
    convenient target for a single river-mouth point."""
    nodes = np.array(
        [
            [139.0, 35.0],   # 0 land
            [139.2, 35.0],   # 1 land
            [139.4, 35.0],   # 2 land
            [139.6, 35.0],   # 3 land
            [139.8, 35.0],   # 4 land
            [139.0, 34.5],   # 5 open
            [139.8, 34.5],   # 6 open
        ],
        dtype=np.float64,
    )
    elements = np.array(
        [
            [0, 5, 1],
            [1, 5, 6],
            [1, 6, 2],
            [2, 6, 3],
            [3, 6, 4],
        ],
        dtype=np.int64,
    )
    return Fort14Mesh(
        title="coast",
        nodes=nodes,
        depths=np.zeros(7),
        elements=elements,
        open_boundaries=[np.array([5, 6], dtype=np.int64)],
        land_boundaries=[(20, np.array([0, 1, 2, 3, 4], dtype=np.int64))],
    )


def test_add_river_inflow_splits_parent_segment() -> None:
    mesh = _coast_mesh()
    point = np.array([[139.4, 35.05]])  # closest land node is index 2
    out, info = add_river_inflow_segments(
        mesh, point, n_nodes_per_river=3, river_ibtype=21,
    )
    assert len(info["rivers"]) == 1
    assert info["rivers"][0]["snapped_node"] == 2
    # Original 1 segment becomes 3: prefix(land), river, suffix(land)
    assert len(out.land_boundaries) == 3
    ibtypes = [ib for ib, _ in out.land_boundaries]
    assert 21 in ibtypes
    river_seg = next(ids for ib, ids in out.land_boundaries if ib == 21)
    assert river_seg.size == 3
    # Prefix + river + suffix should sum to the original size.
    total = sum(ids.size for _, ids in out.land_boundaries)
    assert total == mesh.land_boundaries[0][1].size


def test_add_river_inflow_skips_when_outside_tolerance() -> None:
    mesh = _coast_mesh()
    far_pt = np.array([[180.0, 0.0]])  # equator, far from Tokyo
    out, info = add_river_inflow_segments(
        mesh, far_pt, n_nodes_per_river=3, river_ibtype=21, snap_tol_m=10000.0,
    )
    assert len(info["skipped"]) == 1
    assert info["rivers"] == []
    # Boundary structure unchanged.
    assert len(out.land_boundaries) == 1


def test_add_river_inflow_clamps_at_segment_endpoint() -> None:
    """River centred at the very last node of a 5-node segment with
    n=3 should produce: prefix(2 nodes) + river(3 nodes), no suffix."""
    mesh = _coast_mesh()
    end_pt = np.array([[139.8, 35.05]])  # closest node is index 4 (last)
    out, info = add_river_inflow_segments(
        mesh, end_pt, n_nodes_per_river=3, river_ibtype=21,
    )
    assert info["rivers"][0]["snapped_node"] == 4
    # Two segments: prefix and river (no suffix).
    assert len(out.land_boundaries) == 2
    river_seg = next(ids for ib, ids in out.land_boundaries if ib == 21)
    assert river_seg.tolist() == [2, 3, 4]


def test_add_river_inflow_multiple_points_split_independently() -> None:
    mesh = _coast_mesh()
    pts = np.array([[139.2, 35.05], [139.6, 35.05]])
    out, info = add_river_inflow_segments(
        mesh, pts, n_nodes_per_river=1, river_ibtype=21,
    )
    assert len(info["rivers"]) == 2
    # Original segment of 5 nodes split twice -> at most 5 segments.
    assert len(out.land_boundaries) >= 3
    river_segs = [ids for ib, ids in out.land_boundaries if ib == 21]
    assert len(river_segs) == 2


def test_load_river_points_csv(tmp_path) -> None:
    p = tmp_path / "rivers.csv"
    p.write_text("name,lon,lat\nSumida,139.79,35.65\nTama,139.78,35.55\n")
    pts = load_river_points([p])
    assert pts.shape == (2, 2)
    np.testing.assert_allclose(pts[0], [139.79, 35.65])


def test_load_river_points_missing_file_raises(tmp_path) -> None:
    missing = tmp_path / "nope.csv"
    with pytest.raises(FileNotFoundError):
        load_river_points([missing])
