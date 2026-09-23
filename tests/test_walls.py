"""Walls: a line of edges made into boundary, and nothing else changed."""

from __future__ import annotations

import numpy as np
import pytest

from fvcom_mesh_tools.walls import split_along_walls, wall_edges_from_path


def grid(n=9, h=100.0):
    g = np.arange(n) * h
    gx, gy = np.meshgrid(g, g)
    nodes = np.column_stack([gx.ravel(), gy.ravel()])
    tri = []
    for j in range(n - 1):
        for i in range(n - 1):
            a = j * n + i
            tri += [[a, a + 1, a + n + 1], [a, a + n + 1, a + n]]
    return nodes, np.asarray(tri, dtype=np.int64), n


def node(i, j, n):
    return j * n + i


def boundary_edges(tri):
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    return u[c == 1]


def areas(xy, tri):
    a = xy[tri[:, 1]] - xy[tri[:, 0]]
    b = xy[tri[:, 2]] - xy[tri[:, 0]]
    return 0.5 * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])


def test_a_pier_from_the_coast_is_split_everywhere_but_its_tip():
    xy, tri, n = grid()
    path = [node(4, j, n) for j in range(0, 5)]          # coast at j=0, tip at j=4
    out, t2, copy_of, rep = split_along_walls(xy, tri, wall_edges_from_path(path))
    assert rep["n_free_tips"] == 1 and rep["free_tips"] == [node(4, 4, n)]
    assert rep["n_copies"] == 4, "root and three interior nodes; not the tip"
    assert np.array_equal(out[copy_of], out), "a copy sits exactly on its original"
    assert np.allclose(areas(out, t2), areas(xy, tri)), "no element changed shape"
    assert rep["n_components"] == 1, "a pier does not cut the domain"
    # every wall edge is now two boundary edges
    b = {tuple(r) for r in np.sort(copy_of[boundary_edges(t2)], axis=1).tolist()}
    for a, c in wall_edges_from_path(path).tolist():
        assert (min(a, c), max(a, c)) in b


def test_a_detached_breakwater_has_two_tips_and_makes_an_island():
    xy, tri, n = grid()
    path = [node(2, j, n) for j in range(2, 7)]
    out, t2, copy_of, rep = split_along_walls(xy, tri, wall_edges_from_path(path))
    assert rep["n_free_tips"] == 2
    assert rep["n_copies"] == 3
    assert len(boundary_edges(t2)) == len(boundary_edges(tri)) + 2 * 4, (
        "the island ring is the four wall edges, once per side")


def test_a_junction_of_k_walls_gets_k_minus_one_copies():
    xy, tri, n = grid()
    c = node(4, 4, n)
    walls = np.vstack([
        wall_edges_from_path([node(4, j, n) for j in range(2, 7)]),
        wall_edges_from_path([node(i, 4, n) for i in range(2, 7)]),
    ])
    _, _, copy_of, rep = split_along_walls(xy, tri, walls)
    assert int((copy_of == c).sum()) - 1 == 3
    assert rep["max_sectors_at_a_node"] == 4
    assert rep["n_free_tips"] == 4


def test_a_wall_from_coast_to_coast_is_reported_as_cutting_the_domain():
    xy, tri, n = grid()
    path = [node(4, j, n) for j in range(0, n)]
    _, _, _, rep = split_along_walls(xy, tri, wall_edges_from_path(path))
    assert rep["n_components"] == 2
    assert rep["n_free_tips"] == 0


def test_element_order_is_kept_so_element_data_needs_no_map():
    xy, tri, n = grid()
    path = [node(4, j, n) for j in range(0, 5)]
    _, t2, copy_of, _ = split_along_walls(xy, tri, wall_edges_from_path(path))
    assert np.array_equal(copy_of[t2], tri), (
        "mapping every reference back through copy_of must give the input")


def test_what_cannot_be_a_wall_is_refused():
    xy, tri, n = grid()
    with pytest.raises(ValueError, match="not an edge of the mesh"):
        split_along_walls(xy, tri, [[node(0, 0, n), node(3, 3, n)]])
    with pytest.raises(ValueError, match="only an interior edge"):
        split_along_walls(xy, tri, [[node(0, 0, n), node(1, 0, n)]])
    with pytest.raises(ValueError, match="joins a node to itself"):
        split_along_walls(xy, tri, [[5, 5]])
    with pytest.raises(ValueError, match="does not have"):
        split_along_walls(xy, tri, [[0, 10_000]])


def test_no_wall_is_no_change():
    xy, tri, _ = grid()
    out, t2, copy_of, rep = split_along_walls(xy, tri, np.zeros((0, 2)))
    assert np.array_equal(out, xy) and np.array_equal(t2, tri)
    assert rep["n_copies"] == 0 and rep["n_components"] == 1


# --- extraction --------------------------------------------------------------

def _extract(land, h0):
    from fvcom_mesh_tools.patch import filter_shoreline
    from fvcom_mesh_tools.walls import extract_walls

    area, _ = filter_shoreline(land, h0, elements_per_feature=2)
    return area, *extract_walls(land, area, h0)


def test_a_pier_becomes_one_wall_rooted_in_its_quay_and_reaching_its_tip():
    import shapely

    quay = shapely.box(0.0, 0.0, 400.0, 300.0)
    pier = shapely.box(190.0, 300.0, 202.0, 800.0)      # 12 m x 500 m
    area, walls, rep = _extract(shapely.union_all([quay, pier]), 30.0)
    assert len(walls) == 1, rep
    w = walls[0]
    ends = np.asarray(w.coords)[[0, -1]]
    root, tip = sorted(ends, key=lambda p: p[1])
    assert abs(root[1] - 300.0) < 1.0, "the root must be ON the quay"
    assert abs(tip[1] - 800.0) < 8.0, "the wall must reach the real tip"
    assert np.allclose(np.asarray(w.coords)[:, 0], 196.0, atol=2.0), "centreline"
    assert shapely.equals(area.buffer(0), quay.buffer(0)) or \
        abs(area.area - quay.area) < 0.01 * quay.area, "the quay stays land"
    assert rep["footprint_given_to_water_m2"] == pytest.approx(pier.area, rel=0.05)


def test_an_l_shaped_breakwater_is_one_wall_with_its_corner():
    import shapely

    arm1 = shapely.box(0.0, 0.0, 300.0, 10.0)
    arm2 = shapely.box(290.0, 0.0, 300.0, 200.0)
    _, walls, rep = _extract(shapely.union_all([arm1, arm2]), 30.0)
    assert len(walls) == 1, rep
    assert walls[0].length == pytest.approx(300 + 200 - 10, rel=0.08)


def test_a_wide_structure_stays_land_and_a_stub_is_dropped():
    import shapely

    wide = shapely.box(0.0, 0.0, 70.0, 400.0)            # 70 m >= 2 x 30
    stub = shapely.box(500.0, 0.0, 508.0, 20.0)          # 20 m < L_min = 30
    area, walls, rep = _extract(shapely.union_all([wide, stub]), 30.0)
    assert not walls, rep
    assert area.area == pytest.approx(wide.area, rel=0.01)
    assert rep["n_dropped"] >= 1


def test_mitre_joins_leave_a_square_quay_whole():
    """Round joins shave each convex corner into a crescent."""
    import shapely

    from fvcom_mesh_tools.patch import filter_shoreline

    quay = shapely.box(0.0, 0.0, 400.0, 300.0)
    area, rep = filter_shoreline(quay, 30.0, elements_per_feature=2)
    assert rep["land_lost_m2"] < 1.0


def test_extraction_holds_at_utm_coordinates():
    """In UTM the input is 3.9e6 m with points 0.5 m apart on straight lines.

    qhull stopped on the first real harbour with "a wide merge error"; the
    joggle that fixes that leaves a small loop half-way along a straight
    axis, which cut one pier into two walls until the skeleton was made a
    tree.  Both only happen far from the origin.
    """
    import shapely
    from shapely.affinity import translate

    quay = shapely.box(0.0, 0.0, 400.0, 300.0)
    pier = shapely.box(190.0, 300.0, 202.0, 800.0)
    land = translate(shapely.union_all([quay, pier]), 393000.0, 3909000.0)
    _, walls, rep = _extract(land, 30.0)
    assert len(walls) == 1, rep
    assert walls[0].length == pytest.approx(500.0, abs=8.0)


def test_walls_that_cross_or_butt_are_noded_at_a_shared_vertex():
    import shapely

    from fvcom_mesh_tools.walls import node_walls

    a = shapely.LineString([(0, 0), (100, 0)])
    b = shapely.LineString([(50, -40), (50, 40)])          # crosses a
    c = shapely.LineString([(80, 3), (80, 60)])            # stops 3 m short of a
    out = node_walls([a, b, c], snap_m=5.0)
    ends = [tuple(np.round(np.asarray(g.coords)[e], 6)) for g in out for e in (0, -1)]
    assert ends.count((50.0, 0.0)) == 4, "the crossing is a 4-way vertex"
    assert ends.count((80.0, 0.0)) == 3, "the T is a 3-way vertex"
