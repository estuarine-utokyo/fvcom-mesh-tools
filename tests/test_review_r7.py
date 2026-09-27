"""Coastline rules review, round 7."""

import numpy as np
import shapely
from shapely.geometry import LineString, Point, Polygon, box

from fvcom_mesh_tools.patch import (
    _source_substring,
    hole_polygon,
    island_rings,
    rim_constraints,
    rim_repair,
    select_patch,
    unresolvable_water,
)
from fvcom_mesh_tools.walls import close_wall_pockets
from tests.test_patch import grid_mesh


def _const(v):
    return lambda q: np.full(len(np.atleast_2d(q)), float(v))


def _island_mesh():
    nodes, elements = grid_mesh(15, 15)
    cen = nodes[elements].mean(axis=1)
    land = ((cen[:, 0] > 600) & (cen[:, 0] < 800) & (cen[:, 1] > 600) & (cen[:, 1] < 800))
    elements = elements[~land]
    return nodes, elements, select_patch(nodes, elements, shapely.box(300, 300, 1100, 1100))


def test_resolve_leaves_a_free_island_to_the_source():
    nodes, elements, sel = _island_mesh()
    source = box(650, 650, 750, 750)
    rc = rim_constraints(nodes, sel, size=_const(30.0), coastline="resolve",
                         shoreline=[source.exterior])
    assert rc["n_islands_left_to_source"] == 1
    hole = hole_polygon(rc["pfix"], rc["egfix"])
    assert hole.contains(Point(620, 620))          # base island land is no longer kept
    rings, rep = island_rings(source, hole, _const(30.0))
    assert rep["n_islands_added"] == 1
    assert shapely.Polygon(rings[0]).area == np.float64(10000.0) or \
        abs(shapely.Polygon(rings[0]).area - 10000.0) < 50.0


def test_preserve_still_keeps_a_free_island():
    nodes, elements, sel = _island_mesh()
    rc = rim_constraints(nodes, sel, size=_const(30.0), coastline="preserve")
    assert rc["n_islands_left_to_source"] == 0
    assert not hole_polygon(rc["pfix"], rc["egfix"]).contains(Point(620, 620))


def _stepped_size():
    def size(q):
        x = np.atleast_2d(q)[:, 0]
        return np.where(x < 70, 100.0, np.where(x < 120, 40.0, 200.0))
    size.size_bounds = lambda g, step: (40.0, 100.0 if g.bounds[2] < 120 else 200.0)
    return size


def test_continuous_width_refuses_a_bounded_field_without_a_maximum():
    added, rep = unresolvable_water(box(-1000, -4000, 0, 4000), _stepped_size(),
                                    box(0, -1000, 5, 1000),
                                    radius_factor=0.75, min_h=60.0, spacing=1.0)
    assert added.is_empty and rep["unbounded_pad"]


def test_continuous_width_checks_annuli_out_to_the_maximum():
    size = _stepped_size()
    size.size_max = 200.0
    added, rep = unresolvable_water(box(-1000, -4000, 0, 4000), size,
                                    box(0, -1000, 5, 1000),
                                    radius_factor=0.75, min_h=60.0, spacing=1.0)
    assert added.is_empty and not rep["unbounded_pad"]


def _cone(c, cap):
    c = np.asarray(c, dtype=float)
    return lambda q: np.minimum(cap, 30.0 + np.linalg.norm(np.atleast_2d(q)[:, :2] - c, axis=1))


def test_pocket_without_a_floor_is_not_closed_over_a_fine_spot():
    wall = LineString([(0, 0), (0, 500), (150, 500), (150, 200), (40, 200)])
    _, _, rep = close_wall_pockets([wall], box(-2000, -2000, 2000, 0), _cone((75, 350), 190.0),
                                   min_h=60.0)
    assert rep["n_pockets_closed"] == 0


def test_slit_without_a_floor_keeps_a_fine_basin():
    pts = np.array([(0, 0), (500, 0), (500, 500), (260, 500), (260, 540), (325, 540),
                    (325, 840), (175, 840), (175, 540), (240, 540), (240, 500), (0, 500)],
                   dtype=float)
    e = np.array([(k, (k + 1) % len(pts)) for k in range(len(pts))])
    p, eg, _b, _m, rep = rim_repair(pts, e, np.full(len(pts), -1), Polygon(pts),
                                    _cone((250, 690), 190.0), operations=("slits",),
                                    retreat_tips=False)
    assert hole_polygon(p, eg).contains(Point(250, 690))


def test_short_edge_removal_keeps_a_valid_ring():
    p = np.array([(0, 0), (10, 0), (20, 0), (10, 1)], dtype=float)
    e = np.array([(0, 1), (1, 2), (2, 3), (3, 0)])
    out, eg, _b, _m, rep = rim_repair(p, e, np.array([0, 1, 2, -1]), Polygon(p), _const(30.0))
    hole = hole_polygon(out, eg)
    assert hole.is_valid and hole.area > 0


def test_short_edge_removal_still_removes_a_point():
    pfix = np.array([[0.0, 0.0], [297.0, 0.0], [300.0, 0.0], [300.0, 300.0], [0.0, 300.0]])
    egfix = np.array([[k, (k + 1) % 5] for k in range(5)])
    out, eg, _b, _m, rep = rim_repair(pfix, egfix, np.array([0, -1, 2, 3, 4]),
                                      Polygon(pfix), _const(30.0))
    assert len(out) == 4


def test_arc_choice_sees_the_whole_stretch_and_not_the_ring_start():
    pts = np.vstack([np.column_stack([np.full(160, -20.0), np.linspace(100, 101, 160)]),
                     [(0, 1000), (1000, 1000), (1000, 0), (-20, 0)]])
    for ring in (LineString([(0, 1000), (0, 0), (1000, 0), (1000, 1000), (0, 1000)]),
                 LineString([(0, 0), (1000, 0), (1000, 1000), (0, 1000), (0, 0)])):
        got = _source_substring(pts, ring)
        assert got is not None
        assert abs(LineString(got).length - 3900.0) < 1.0
