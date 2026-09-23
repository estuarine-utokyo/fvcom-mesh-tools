"""The hires branch: its fork, its weight, and the things it must refuse."""

from __future__ import annotations

import numpy as np
import pytest
import shapely

from fvcom_mesh_tools.bathy_patch import (
    blend_weights,
    edge_slopes,
    feasibility_margin,
    interface_lines,
    patch_depths,
)


def square_patch():
    """A 1000 m square cut out of a 3000 m mesh, so there is a real interface."""
    n = 7
    g = np.linspace(0.0, 3000.0, n)
    gx, gy = np.meshgrid(g, g)
    nodes = np.column_stack([gx.ravel(), gy.ravel()])
    tri = []
    for j in range(n - 1):
        for i in range(n - 1):
            a = j * n + i
            tri += [[a, a + 1, a + n + 1], [a, a + n + 1, a + n]]
    tri = np.asarray(tri, dtype=np.int64)
    hole = shapely.box(1000.0, 1000.0, 2000.0, 2000.0)
    cen = nodes[tri].mean(axis=1)
    inside = shapely.contains(hole, shapely.points(cen[:, 0], cen[:, 1]))
    return nodes, tri, np.flatnonzero(~inside), hole


def test_the_weight_is_one_on_the_core_and_zero_at_the_interface():
    nodes, tri, retained, hole = square_patch()
    iface = interface_lines(nodes, tri, retained)
    assert not iface.is_empty, "a cut in the middle of a mesh has an interface"
    core = shapely.box(1350.0, 1350.0, 1650.0, 1650.0)

    on_core = np.array([[1500.0, 1500.0], [1350.0, 1500.0], [1650.0, 1650.0]])
    assert np.allclose(blend_weights(on_core, [core], iface), 1.0)

    on_iface = np.array(shapely.get_coordinates(iface))
    assert np.allclose(blend_weights(on_iface, [core], iface), 0.0)

    between = np.array([[1200.0, 1500.0], [1100.0, 1500.0], [1050.0, 1500.0]])
    w = blend_weights(between, [core], iface)
    assert np.all(np.diff(w) < 0), "the weight has to fall towards the interface"
    assert np.all((w > 0) & (w < 1))


def test_the_interface_is_the_actual_cut_not_the_analytic_buffer():
    """select_patch takes whole faces and grows; the weight must follow THAT.

    Here the cut is grown by one more ring of faces than the declared square,
    so a weight keyed to the square's buffer would not vanish where the mesh
    actually meets the frozen zone.
    """
    nodes, tri, retained, hole = square_patch()
    grown = np.array([k for k in retained
                      if not shapely.intersects(
                          hole.buffer(500.0),
                          shapely.Point(nodes[tri[k]].mean(axis=0)))])
    iface = interface_lines(nodes, tri, grown)
    core = shapely.box(1350.0, 1350.0, 1650.0, 1650.0)
    w = blend_weights(np.array(shapely.get_coordinates(iface)), [core], iface)
    assert np.allclose(w, 0.0), (
        "w must vanish on the interface the cut HAS, not the one it was asked for")


def test_the_coastline_is_not_forced_to_zero():
    """A free coastline inside the cut is hole boundary and must take the source.

    Weighting it to zero would give the newly resolved shore its BASE depths,
    which is the opposite of what the branch is for.
    """
    nodes, tri, retained, hole = square_patch()
    rim_nodes = np.unique(tri[retained])
    coast = rim_nodes[nodes[rim_nodes][:, 1] <= 1000.0]
    full = interface_lines(nodes, tri, retained)
    minus_coast = interface_lines(nodes, tri, retained, coast_nodes=coast)
    assert minus_coast.length < full.length, "the coast must have been dropped"


def test_a_patch_with_no_interface_takes_the_source_everywhere():
    nodes, tri, retained, hole = square_patch()
    empty = interface_lines(nodes, tri, np.zeros(0, dtype=np.int64))
    core = shapely.box(1350.0, 1350.0, 1650.0, 1650.0)
    w = blend_weights(np.array([[1500.0, 1500.0], [2900.0, 2900.0]]), [core], empty)
    assert np.allclose(w, 1.0)


def test_overlapping_regions_need_no_separate_rule():
    nodes, tri, retained, hole = square_patch()
    iface = interface_lines(nodes, tri, retained)
    a = shapely.box(1300.0, 1300.0, 1600.0, 1600.0)
    b = shapely.box(1500.0, 1500.0, 1800.0, 1800.0)
    shared = np.array([[1550.0, 1550.0]])
    assert np.allclose(blend_weights(shared, [a, b], iface), 1.0)


def test_the_feasibility_margin_reproduces_the_counterexample():
    """The review's counterexample, which is why the continuity claim went.

    Base 3 m, survey 300 m, ramp over 330 m: an edge 30 m inside the interface
    carries r = 27/33 = 0.818, so a 30 m target is not feasible at rmax = 0.2.
    """
    m = feasibility_margin(target_h_m=30.0, transition_m=330.0,
                           depth_m=16.5, difference_m=297.0, rmax=0.2)
    assert not m["ok"]
    assert m["allowed_edge_m"] == pytest.approx(2.0 * 16.5 * 330.0 * 0.2 / 297.0)
    assert m["expected_r_at_target"] > 0.2
    easy = feasibility_margin(target_h_m=30.0, transition_m=2000.0,
                              depth_m=4.0, difference_m=0.5, rmax=0.2)
    assert easy["ok"]


def test_a_flat_difference_is_always_feasible():
    m = feasibility_margin(30.0, 1000.0, 4.0, 0.0)
    assert m["ok"] and np.isinf(m["allowed_edge_m"])


def test_the_slope_survives_an_intertidal_node():
    """|dh|/L is defined where |dh|/(hi+hj) is not, and both are reported."""
    nodes = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0], [100.0, 100.0]])
    tri = np.array([[0, 1, 3], [0, 3, 2]], dtype=np.int64)
    depths = np.array([5.0, -6.0, 1.0, 2.0])       # node 1 is 6 m above the datum
    rep = edge_slopes(nodes, tri, depths, [False, True, True, True])
    assert rep["n_r_undefined"] >= 1, "an intertidal pair must be flagged, not hidden"
    assert rep["all_changed"]["slope_max"] > 0
    assert rep["seam"]["n"] > 0, "node 0 is retained, so its edges are the seam"


def test_scope_core_with_a_ramp_is_refused():
    with pytest.raises(ValueError, match="not a ramp"):
        patch_depths(np.zeros((1, 2)), np.zeros((1, 2)), np.ones(1),
                     regions=[], interface=None, scope="core", blend="ramp")


def test_the_blend_meets_the_base_exactly_where_the_weight_is_zero():
    """No real product is read: the ladder is replaced by a known field.

    What is under test is the arithmetic of the ramp, and the property that
    matters is that a w = 0 node receives the base value EXACTLY -- that is
    what lets the frozen depths be met rather than approached.
    """
    import fvcom_mesh_tools.dem.tokyo_bay as tb

    nodes, tri, retained, hole = square_patch()
    iface = interface_lines(nodes, tri, retained)
    core = shapely.box(1350.0, 1350.0, 1650.0, 1650.0)
    xy = np.vstack([np.array(shapely.get_coordinates(iface))[:3],
                    np.array([[1500.0, 1500.0]])])
    base = np.array([7.0, 8.0, 9.0, 10.0])
    real = tb.sample
    tb.sample = lambda lon, lat, extrapolate=True: (
        np.full(len(np.atleast_1d(lon)), 99.0),
        np.zeros(len(np.atleast_1d(lon)), dtype=np.int8),
        np.zeros(len(np.atleast_1d(lon))))
    try:
        depth, rep = patch_depths(np.zeros((4, 2)), xy, base,
                                  regions=[core], interface=iface)
    finally:
        tb.sample = real
    assert np.allclose(depth[:3], base[:3]), "w = 0 must give the base exactly"
    assert depth[3] == pytest.approx(99.0), "w = 1 must give the source exactly"
    assert rep["n_at_w0"] == 3 and rep["n_at_w1"] == 1
