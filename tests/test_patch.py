"""Tests for the patch generator.

The claims worth testing here are contract claims, not numerical ones: that
the cut is a surface with a boundary, that the rim comes back as a closed walk
with every segment constrained, that stitching preserves the retained mesh
exactly, and that the repair pass cannot touch what it does not own.  A mesh
that looks plausible and quietly moved a frozen node is the failure this
module exists to make impossible, so most of these tests assert on identity
rather than on tolerance.
"""

from __future__ import annotations

import numpy as np
import pytest
import shapely

from fvcom_mesh_tools.patch import (
    ambient_size_field,
    boundary_rings,
    coastline_points,
    effective_gradation,
    hole_polygon,
    improve_patch,
    patch_sizing,
    rim_constraints,
    select_patch,
    stitch_patch,
    verify_patch,
)


def grid_mesh(nx: int = 9, ny: int = 9, h: float = 100.0):
    """A right-triangulated rectangle: nodes (N,2) in metres, elements (M,3) CCW."""
    gx, gy = np.meshgrid(np.arange(nx) * h, np.arange(ny) * h)
    nodes = np.column_stack([gx.ravel(), gy.ravel()])
    tri = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            b, c, d = a + 1, a + nx, a + nx + 1
            tri += [[a, b, d], [a, d, c]]
    return nodes, np.asarray(tri, dtype=np.int64)


def ccw(nodes, elements):
    a = nodes[elements[:, 1]] - nodes[elements[:, 0]]
    b = nodes[elements[:, 2]] - nodes[elements[:, 0]]
    return a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]


# ---------------------------------------------------------------- selection


def test_grid_fixture_is_counter_clockwise():
    nodes, elements = grid_mesh()
    assert (ccw(nodes, elements) > 0).all()


def test_select_takes_whole_triangles_and_reports_the_overreach():
    nodes, elements = grid_mesh()
    foot = shapely.Point(400.0, 400.0).buffer(150.0)
    sel = select_patch(nodes, elements, foot)
    assert sel.n_removed > 0
    # the requested reach is 150 m; whole triangles push past it
    assert sel.report["selection_reach_m"] > 150.0
    # every removed element's centroid really was inside
    cen = nodes[elements[sel.removed]].mean(axis=1)
    assert shapely.contains(foot, shapely.points(cen[:, 0], cen[:, 1])).all()


def test_rim_is_a_closed_walk_with_degree_two_everywhere():
    nodes, elements = grid_mesh()
    sel = select_patch(nodes, elements, shapely.Point(400, 400).buffer(150.0))
    deg = np.zeros(len(nodes), dtype=int)
    np.add.at(deg, sel.rim_edges.ravel(), 1)
    assert set(np.unique(deg[deg > 0]).tolist()) == {2}
    assert len(sel.rings) == 1
    assert len(sel.rings[0]) == len(sel.rim_edges)


def test_interface_rim_nodes_are_all_frozen():
    """The structural fact the coastline logic relies on."""
    nodes, elements = grid_mesh()
    sel = select_patch(nodes, elements, shapely.Point(400, 400).buffer(150.0))
    frozen = set(sel.frozen_nodes.tolist())
    interface = sel.rim_edges[~sel.physical_rim]
    assert interface.size
    assert set(interface.ravel().tolist()) <= frozen


def test_selection_reaching_the_edge_marks_physical_rim():
    nodes, elements = grid_mesh()
    sel = select_patch(nodes, elements, shapely.Point(0, 400).buffer(250.0))
    assert sel.physical_rim.any()
    assert sel.report["n_physical_rim_edges"] > 0


def test_empty_selection_is_refused():
    nodes, elements = grid_mesh()
    with pytest.raises(ValueError, match="selects no elements"):
        select_patch(nodes, elements, shapely.Point(-9999, -9999).buffer(10.0))


def test_taking_everything_is_refused():
    nodes, elements = grid_mesh()
    with pytest.raises(ValueError, match="entire mesh"):
        select_patch(nodes, elements, shapely.box(-1e4, -1e4, 1e4, 1e4))


def test_cut_reaching_the_open_boundary_is_refused():
    nodes, elements = grid_mesh()
    obc = np.flatnonzero(nodes[:, 1] == 0.0)
    with pytest.raises(ValueError, match="open boundary"):
        select_patch(nodes, elements, shapely.Point(400, 50).buffer(120.0),
                     open_boundary_nodes=obc)


def test_open_boundary_guard_band_refuses_a_near_miss():
    nodes, elements = grid_mesh()
    obc = np.flatnonzero(nodes[:, 1] == 0.0)
    foot = shapely.Point(400, 400).buffer(120.0)
    select_patch(nodes, elements, foot, open_boundary_nodes=obc)  # fine bare
    with pytest.raises(ValueError, match="guard band"):
        select_patch(nodes, elements, foot, open_boundary_nodes=obc,
                     obc_guard_m=1000.0)


def test_cut_that_severs_the_mesh_is_refused():
    nodes, elements = grid_mesh()
    band = shapely.box(-10.0, 320.0, 1e4, 400.0)   # all of element row 3
    with pytest.raises(ValueError, match="splits the retained mesh"):
        select_patch(nodes, elements, band)


def test_pinch_is_repaired_by_growing_the_cut():
    """Two diagonal squares touching at a corner: the shared vertex has four
    rim edges, so no closed walk exists until the cut grows."""
    nodes, elements = grid_mesh()
    foot = shapely.union(shapely.Point(250, 250).buffer(110.0),
                         shapely.Point(450, 450).buffer(110.0))
    sel = select_patch(nodes, elements, foot)
    deg = np.zeros(len(nodes), dtype=int)
    np.add.at(deg, sel.rim_edges.ravel(), 1)
    assert set(np.unique(deg[deg > 0]).tolist()) == {2}


# -------------------------------------------------------------------- rings


def test_nested_rings_alternate_hole_parity():
    outer = np.array([[0, 0], [10, 0], [10, 10], [0, 10]], dtype=float)
    inner = np.array([[3, 3], [7, 3], [7, 7], [3, 7]], dtype=float)
    xy = np.vstack([outer, inner])
    edges = np.array([[0, 1], [1, 2], [2, 3], [0, 3],
                      [4, 5], [5, 6], [6, 7], [4, 7]])
    rings, is_hole = boundary_rings(xy, edges)
    assert len(rings) == 2
    assert sorted(is_hole) == [False, True]


def test_rim_with_a_bad_degree_is_rejected():
    xy = np.array([[0, 0], [1, 0], [2, 0]], dtype=float)
    with pytest.raises(ValueError, match="exactly two rim edges"):
        boundary_rings(xy, np.array([[0, 1], [1, 2]]))


# ---------------------------------------------------------------- coastline


def test_preserve_keeps_every_point_on_the_base_polyline():
    pts = np.array([[0, 0], [100, 30], [250, 10], [400, 60]], dtype=float)
    out = coastline_points(pts, 20.0, mode="preserve")
    line = shapely.LineString(pts)
    assert float(shapely.distance(shapely.points(out), line).max()) < 1e-6


def test_every_mode_keeps_the_frozen_endpoints():
    pts = np.array([[0, 0], [100, 30], [250, 10], [400, 60]], dtype=float)
    for mode in ("preserve", "spline"):
        out = coastline_points(pts, 25.0, mode=mode)
        assert np.allclose(out[0], pts[0])
        assert np.allclose(out[-1], pts[-1])


def test_spacing_follows_a_size_field():
    pts = np.array([[0.0, 0.0], [1000.0, 0.0]])

    def size(q):
        return 10.0 + 0.2 * np.asarray(q)[:, 0]

    out = coastline_points(pts, size, mode="preserve")
    step = np.linalg.norm(np.diff(out, axis=0), axis=1)
    assert step[0] < step[-1] / 3  # fine at the start, coarse at the end
    assert np.allclose(out[:, 1], 0.0)


def test_constant_size_gives_even_spacing():
    pts = np.array([[0.0, 0.0], [300.0, 0.0]])
    out = coastline_points(pts, 30.0, mode="preserve")
    step = np.linalg.norm(np.diff(out, axis=0), axis=1)
    assert np.allclose(step, step[0], rtol=0.25)
    assert abs(step.mean() - 30.0) < 10.0


def test_resample_without_a_shoreline_is_refused():
    pts = np.array([[0.0, 0.0], [100.0, 0.0]])
    with pytest.raises(ValueError, match="needs the source shoreline"):
        coastline_points(pts, 10.0, mode="resample")


def test_unknown_coastline_mode_is_refused():
    pts = np.array([[0.0, 0.0], [100.0, 0.0]])
    with pytest.raises(ValueError, match="unknown coastline mode"):
        coastline_points(pts, 10.0, mode="recut")


def test_resample_beyond_the_tolerance_is_refused():
    pts = np.array([[0.0, 0.0], [400.0, 0.0]])
    detour = shapely.LineString([[0, 0], [200, 300], [400, 0]])
    with pytest.raises(ValueError, match="beyond the"):
        coastline_points(pts, 20.0, mode="resample", shoreline=detour,
                         tolerance_m=50.0)


def test_resample_within_the_tolerance_follows_the_source():
    pts = np.array([[0.0, 0.0], [400.0, 0.0]])
    detour = shapely.LineString([[0, 0], [200, 20], [400, 0]])
    out = coastline_points(pts, 20.0, mode="resample", shoreline=detour,
                           tolerance_m=50.0)
    assert out[:, 1].max() > 5.0
    assert np.allclose(out[0], pts[0]) and np.allclose(out[-1], pts[-1])


def test_a_nonpositive_size_is_refused():
    pts = np.array([[0.0, 0.0], [100.0, 0.0]])
    with pytest.raises(ValueError, match="finite and positive"):
        coastline_points(pts, lambda q: np.zeros(len(q)), mode="preserve")


# --------------------------------------------------------------------- rim


def rim_of(nodes, elements, foot, **kw):
    sel = select_patch(nodes, elements, foot)
    return sel, rim_constraints(nodes, sel, size=60.0, coastline="preserve", **kw)


def test_rim_constraints_emit_one_segment_per_rim_point():
    nodes, elements = grid_mesh()
    sel, rc = rim_of(nodes, elements, shapely.Point(400, 400).buffer(150.0))
    assert rc["n_egfix"] == rc["n_pfix"]
    assert rc["egfix"].max() < rc["n_pfix"]
    # every fixed point is used by exactly two segments: a closed ring
    deg = np.zeros(rc["n_pfix"], dtype=int)
    np.add.at(deg, rc["egfix"].ravel(), 1)
    assert set(np.unique(deg).tolist()) == {2}


def test_frozen_rim_nodes_keep_their_coordinates_and_identity():
    nodes, elements = grid_mesh()
    sel, rc = rim_of(nodes, elements, shapely.Point(400, 400).buffer(150.0))
    kept = rc["pfix_base"] >= 0
    assert np.allclose(rc["pfix"][kept], nodes[rc["pfix_base"][kept], :2])
    assert set(sel.frozen_nodes.tolist()) >= set(rc["pfix_base"][kept].tolist())


def test_an_interior_cut_replaces_no_coastline():
    nodes, elements = grid_mesh()
    _, rc = rim_of(nodes, elements, shapely.Point(400, 400).buffer(150.0))
    assert rc["n_coastline_nodes_replaced"] == 0
    assert (rc["pfix_base"] >= 0).all()


def test_a_cut_reaching_the_edge_replaces_the_free_stretch():
    nodes, elements = grid_mesh()
    sel = select_patch(nodes, elements, shapely.Point(0, 400).buffer(250.0))
    rc = rim_constraints(nodes, sel, size=25.0, coastline="preserve")
    assert rc["n_coastline_nodes_replaced"] > 0
    assert rc["n_coastline_nodes_new"] > rc["n_coastline_nodes_replaced"]


def test_hole_polygon_closes_and_matches_the_cut_area():
    nodes, elements = grid_mesh()
    sel, rc = rim_of(nodes, elements, shapely.Point(400, 400).buffer(150.0))
    poly = hole_polygon(rc["pfix"], rc["egfix"])
    assert poly.is_valid
    a = nodes[elements[sel.removed]]
    area = 0.5 * np.abs((a[:, 1, 0] - a[:, 0, 0]) * (a[:, 2, 1] - a[:, 0, 1])
                        - (a[:, 1, 1] - a[:, 0, 1]) * (a[:, 2, 0] - a[:, 0, 0]))
    assert poly.area == pytest.approx(area.sum(), rel=1e-6)


def test_unclosed_rim_is_refused():
    with pytest.raises(ValueError, match="exactly two rim edges|do not close"):
        hole_polygon(np.array([[0.0, 0.0], [1.0, 0.0]]), np.array([[0, 1]]))


# ------------------------------------------------------------------ sizing


def test_ambient_field_recovers_a_uniform_edge_length():
    nodes, elements = grid_mesh(h=100.0)
    amb = ambient_size_field(nodes, elements)
    # right triangles: legs 100, hypotenuse 141
    assert 100.0 < float(np.median(amb)) < 142.0


def test_smoothing_flattens_the_ambient_field():
    nodes, elements = grid_mesh(11, 11)
    nodes = nodes.copy()
    nodes[:, 0] *= 1.0 + 0.6 * (nodes[:, 1] / nodes[:, 1].max())
    e = np.unique(np.sort(np.vstack([elements[:, [0, 1]], elements[:, [1, 2]],
                                     elements[:, [2, 0]]]), axis=1), axis=0)
    ln = np.linalg.norm(nodes[e[:, 0]] - nodes[e[:, 1]], axis=1)
    raw = ambient_size_field(nodes, elements, smooth_passes=0)
    smooth = ambient_size_field(nodes, elements, smooth_passes=20)
    slope = lambda h: np.abs(h[e[:, 0]] - h[e[:, 1]]) / ln  # noqa: E731
    assert slope(smooth).max() < slope(raw).max()
    assert np.median(smooth) == pytest.approx(np.median(raw), rel=0.15)


def test_sizing_is_the_target_inside_and_the_base_outside():
    nodes, elements = grid_mesh(15, 15)
    core = shapely.Point(700, 700).buffer(100.0)
    fh = patch_sizing(nodes, elements, [(core, 20.0, 400.0)], distmesh_scale=1.0)
    assert fh(np.array([[700.0, 700.0]]))[0] == pytest.approx(20.0)
    far = fh(np.array([[100.0, 100.0]]))[0]
    assert far > 90.0  # the base size out there, not the target


def test_sizing_never_coarsens_the_base_mesh():
    nodes, elements = grid_mesh(15, 15)
    core = shapely.Point(700, 700).buffer(100.0)
    fh = patch_sizing(nodes, elements, [(core, 5000.0, 400.0)], distmesh_scale=1.0)
    amb = ambient_size_field(nodes, elements)
    q = np.column_stack([np.linspace(200, 1200, 50), np.full(50, 700.0)])
    assert (fh(q) <= amb.max() + 1e-9).all()


def test_distmesh_scale_divides_the_field():
    nodes, elements = grid_mesh(15, 15)
    core = shapely.Point(700, 700).buffer(100.0)
    a = patch_sizing(nodes, elements, [(core, 20.0, 400.0)], distmesh_scale=1.0)
    b = patch_sizing(nodes, elements, [(core, 20.0, 400.0)], distmesh_scale=1.2)
    q = np.array([[700.0, 700.0], [900.0, 700.0]])
    assert np.allclose(a(q) / 1.2, b(q))


def test_effective_gradation_flags_a_steep_ramp():
    nodes, elements = grid_mesh(15, 15, h=100.0)
    core = shapely.Point(700, 700).buffer(100.0)
    easy = effective_gradation(nodes, elements, [(core, 20.0, 2000.0)])
    hard = effective_gradation(nodes, elements, [(core, 20.0, 150.0)])
    assert easy["ramp_below_reference"]
    assert not hard["ramp_below_reference"]
    assert hard["max_effective_gradation"] > easy["max_effective_gradation"]


# ---------------------------------------------------------------- stitching


def replay_patch(nodes, elements, foot):
    """Cut, then fill the hole with exactly what was removed.

    The identity fill: the stitched mesh must come back equal to the input.
    Anything the stitch gets wrong shows up here with nothing else in the way.
    """
    sel = select_patch(nodes, elements, foot)
    rc = rim_constraints(nodes, sel, size=1e9, coastline="preserve")
    taken = np.unique(elements[sel.removed])
    local = {int(g): i for i, g in enumerate(taken)}
    ptri = np.array([[local[int(v)] for v in row]
                     for row in elements[sel.removed]], dtype=np.int64)
    return sel, rc, nodes[taken, :2], ptri


def test_identity_fill_reproduces_the_base_mesh():
    nodes, elements = grid_mesh()
    depths = 5.0 + nodes[:, 0] / 1000.0
    foot = shapely.Point(400, 400).buffer(150.0)
    sel, rc, pn, pt = replay_patch(nodes, elements, foot)
    out, oute, outd, nm, rep = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                            rc["pfix"], rc["pfix_base"])
    assert rep["n_nodes"] == len(nodes)
    assert rep["n_elements"] == len(elements)
    # The node numbering changes -- deleted rows are refilled at the end -- so
    # the meshes are compared by geometry, which is what "the same mesh" means.
    cell = lambda xy, t: {tuple(sorted(  # noqa: E731
        (round(x, 6), round(y, 6)) for x, y in xy[r])) for r in t.tolist()}
    assert cell(out, oute) == cell(nodes, elements)
    keep = sel.frozen_nodes
    assert np.allclose(out[nm[keep]], nodes[keep, :2])
    assert np.allclose(outd[nm[keep]], depths[keep])


def test_verify_passes_on_the_identity_fill():
    from fvcom_mesh_tools.patch import boundary_after_patch

    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    out, oute, outd, nm, st = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                           rc["pfix"], rc["pfix_base"])
    ver = verify_patch(nodes, depths, elements, sel, out, oute, outd, nm,
                       expected_boundary=boundary_after_patch(
                           elements, sel, rc, nm, st["pfix_new"]))
    assert ver["ok"]
    assert ver["n_frozen_moved"] == 0
    assert ver["n_retained_faces_missing"] == 0


def test_a_lost_fixed_point_is_reported_not_snapped():
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    pn = pn.copy()
    pn[0] += 50.0  # the fill moved a constrained point
    with pytest.raises(ValueError, match="constrained points are absent"):
        stitch_patch(nodes, elements, depths, sel, pn, pt,
                     rc["pfix"], rc["pfix_base"])


def test_verify_catches_a_moved_frozen_node():
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    out, oute, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                          rc["pfix"], rc["pfix_base"])
    out = out.copy()
    out[nm[sel.frozen_nodes[0]]] += 0.01
    ver = verify_patch(nodes, depths, elements, sel, out, oute, outd, nm)
    assert not ver["ok"]
    assert ver["n_frozen_moved"] == 1


def test_verify_catches_a_missing_retained_face():
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    out, oute, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                          rc["pfix"], rc["pfix_base"])
    ver = verify_patch(nodes, depths, elements, sel, out, oute[1:], outd, nm)
    assert not ver["ok"]
    assert ver["n_retained_faces_missing"] >= 1


def test_verify_catches_a_changed_open_boundary():
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    obc = np.flatnonzero(nodes[:, 1] == nodes[:, 1].max())
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 200).buffer(150.0))
    out, oute, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                          rc["pfix"], rc["pfix_base"])
    good = verify_patch(nodes, depths, elements, sel, out, oute, outd, nm,
                        open_boundaries=[obc])
    assert good["open_boundary_unchanged"]
    out = out.copy()
    out[nm[obc[0]]] += 1.0
    bad = verify_patch(nodes, depths, elements, sel, out, oute, outd, nm,
                       open_boundaries=[obc])
    assert not bad["open_boundary_unchanged"]
    assert not bad["ok"]


def test_verify_refuses_a_non_finite_coordinate():
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    out, oute, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                          rc["pfix"], rc["pfix_base"])
    out = out.copy()
    out[nm[sel.frozen_nodes[0]]] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        verify_patch(nodes, depths, elements, sel, out, oute, outd, nm)


def test_depths_of_new_nodes_come_from_the_base_field():
    nodes, elements = grid_mesh()
    depths = 5.0 + nodes[:, 0] / 100.0
    sel = select_patch(nodes, elements, shapely.Point(400, 400).buffer(150.0))
    rc = rim_constraints(nodes, sel, size=1e9, coastline="preserve")
    taken = np.unique(elements[sel.removed])
    local = {int(g): i for i, g in enumerate(taken)}
    ptri = np.array([[local[int(v)] for v in row]
                     for row in elements[sel.removed]], dtype=np.int64)
    pn = nodes[taken, :2].copy()
    # add one genuinely new vertex by splitting nothing: reuse an interior point
    _, _, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, ptri,
                                     rc["pfix"], rc["pfix_base"])
    assert np.allclose(outd[nm[sel.frozen_nodes]], depths[sel.frozen_nodes])


# ------------------------------------------------------------------ repair


def perturbed_patch():
    """A grid whose middle nodes have been jittered into poor triangles."""
    nodes, elements = grid_mesh(11, 11)
    nodes = nodes.copy()
    rng = np.random.default_rng(3)
    inner = ((nodes[:, 0] > 150) & (nodes[:, 0] < 850)
             & (nodes[:, 1] > 150) & (nodes[:, 1] < 850))
    nodes[inner] += rng.uniform(-38.0, 38.0, (int(inner.sum()), 2))
    return nodes, elements, inner


def test_repair_never_moves_a_node_it_does_not_own():
    nodes, elements, inner = perturbed_patch()
    movable = inner.copy()
    movable[np.flatnonzero(inner)[:5]] = False  # deliberately withheld
    out, _, _ = improve_patch(nodes, elements, movable,
                              np.ones(len(elements), dtype=bool))
    assert np.allclose(out[~movable], nodes[~movable])


def test_repair_never_rewrites_an_immutable_face():
    nodes, elements, inner = perturbed_patch()
    mutable = np.zeros(len(elements), dtype=bool)
    mutable[: len(elements) // 2] = True
    _, out, _ = improve_patch(nodes, elements, inner, mutable)
    assert np.array_equal(np.sort(out[~mutable], axis=1),
                          np.sort(elements[~mutable], axis=1))


def test_repair_improves_the_worst_angle_and_never_inverts():
    nodes, elements, inner = perturbed_patch()
    from fvcom_mesh_tools.patch import _angles_deg

    before = _angles_deg(nodes, elements).min()
    out_xy, out_t, info = improve_patch(nodes, elements, inner,
                                        np.ones(len(elements), dtype=bool))
    assert info["min_angle_deg"] >= before
    assert (ccw(out_xy, out_t) > 0).all()


def test_repair_is_a_no_op_when_nothing_may_move():
    nodes, elements, _ = perturbed_patch()
    out_xy, out_t, info = improve_patch(
        nodes, elements, np.zeros(len(nodes), dtype=bool),
        np.zeros(len(elements), dtype=bool))
    assert info["n_flips"] == 0 and info["n_moves"] == 0
    assert np.allclose(out_xy, nodes) and np.array_equal(out_t, elements)


def test_sliding_stays_on_the_boundary_chain():
    nodes, elements, _ = perturbed_patch()
    _u, _c = np.unique(np.sort(np.vstack([elements[:, [0, 1]], elements[:, [1, 2]],
                                          elements[:, [2, 0]]]), axis=1),
                       axis=0, return_counts=True)
    on_b = np.zeros(len(nodes), dtype=bool)
    on_b[np.unique(_u[_c == 1])] = True
    ring = shapely.LineString(nodes[np.append(_u[_c == 1][:, 0], _u[_c == 1][0, 0])])
    out, _, _ = improve_patch(nodes, elements, np.zeros(len(nodes), dtype=bool),
                              np.ones(len(elements), dtype=bool), slidable=on_b)
    moved = np.linalg.norm(out - nodes, axis=1) > 1e-9
    assert not moved[~on_b].any()
    if moved.any():
        assert float(shapely.distance(shapely.points(out[moved]), ring).max()) < 1e-6


# ------------------------------------------ what the adversarial review found
#
# gpt-6-astra reviewed the first implementation (2026-09-22,
# docs/local_refine_implementation_review.md) and supplied nine failing tests.
# All nine reproduced. These are the regressions for the fixes; each one
# failed before the fix and the review's own wording is kept where it is
# sharper than mine.


def test_a_duplicated_retained_face_is_rejected():
    """A set comparison is happy to see the same face twice."""
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    out, oute, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                          rc["pfix"], rc["pfix_base"])
    doubled = np.vstack([oute, oute[0]])
    ver = verify_patch(nodes, depths, elements, sel, out, doubled, outd, nm)
    assert not ver["ok"]
    assert ver["n_extra_faces"] == 1
    assert ver["n_nonmanifold_edges"] > 0


def test_frozen_means_exact_not_within_a_micrometre():
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    out, oute, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                          rc["pfix"], rc["pfix_base"])
    v = nm[sel.frozen_nodes[0]]
    out = out.copy()
    outd = outd.copy()
    out[v, 0] += 5e-7
    outd[v] += 5e-7
    ver = verify_patch(nodes, depths, elements, sel, out, oute, outd, nm)
    assert not ver["ok"]
    assert not ver["frozen_exact"]


def test_nested_rings_are_not_filled_in():
    """Four nested squares plus a disjoint one: area 57, not 101."""
    xy = np.vstack([np.array([[a, a], [b, a], [b, b], [a, b]], float)
                    for a, b in [(0, 10), (1, 9), (2, 8), (3, 7), (20, 21)]])
    edges = np.array([[i, i // 4 * 4 + (i + 1) % 4] for i in range(len(xy))])
    assert boundary_rings(xy, edges)[1] == [False, True, False, True, False]
    assert hole_polygon(xy, edges).area == pytest.approx(57.0)


def test_a_physical_island_inside_the_cut_stitches():
    """Its rim nodes are free, and calling them frozen made the stitch raise."""
    nodes, elements = grid_mesh()
    cen = nodes[elements].mean(axis=1)
    land = ((cen[:, 0] > 300) & (cen[:, 0] < 400)
            & (cen[:, 1] > 300) & (cen[:, 1] < 400))
    elements = elements[~land]
    sel = select_patch(nodes, elements, shapely.Point(350, 350).buffer(220.0))
    assert sum(sel.ring_is_hole) == 1
    rc = rim_constraints(nodes, sel, size=1e9, coastline="preserve")
    taken = np.unique(elements[sel.removed])
    local = np.full(len(nodes), -1, dtype=np.int64)
    local[taken] = np.arange(len(taken))
    out, oute, outd, nm, st = stitch_patch(
        nodes, elements, np.ones(len(nodes)), sel, nodes[taken, :2],
        local[elements[sel.removed]], rc["pfix"], rc["pfix_base"])
    from fvcom_mesh_tools.patch import boundary_after_patch

    assert verify_patch(nodes, np.ones(len(nodes)), elements, sel,
                        out, oute, outd, nm,
                        expected_boundary=boundary_after_patch(
                            elements, sel, rc, nm, st["pfix_new"]))["ok"]


def test_a_zero_width_transition_stays_finite():
    """Ambient already at the target is legitimate, and 0/0 is not."""
    nodes, elements = grid_mesh()
    fh = patch_sizing(nodes, elements,
                      [(shapely.Point(400, 400).buffer(100.0), 200.0, 0.0)],
                      distmesh_scale=1.0)
    q = np.array([[400.0, 400.0], [700.0, 700.0]])
    assert np.isfinite(fh(q)).all()


def test_the_repaired_mesh_stays_inside_the_valence_gate():
    """The guard counted incident faces over a subset, not over the mesh, and
    read [5,5,7,7] where the mesh had [5,5,7,9].  Valence is now tracked
    globally; it is a cost on a flip rather than a veto, because forbidding
    every intermediate excess also blocks the sequences that end below the
    limit, so what is asserted is the finished mesh."""
    nodes, elements, inner = perturbed_patch()
    for rounds in (1, 30):
        _, out_t, _ = improve_patch(nodes, elements, inner,
                                    np.ones(len(elements), dtype=bool),
                                    rounds=rounds, max_valence=8)
        assert np.bincount(out_t.ravel(), minlength=len(nodes)).max() <= 8


def test_a_vertex_pinch_in_the_retained_mesh_is_not_accepted():
    """Four faces round one vertex, middle two cut: the rest meet at a point."""
    xy = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [-1, 1], [-1, 0]], float)
    tri = np.array([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]])
    with pytest.raises(ValueError, match="splits the retained mesh|entire mesh"):
        select_patch(xy, tri, shapely.Point(0, 2 / 3).buffer(0.4))


def test_an_isolated_retained_face_is_absorbed_into_the_cut():
    """A spike joined only at its vertices; the cut takes it rather than refuse."""
    nodes, elements = grid_mesh()
    sel = select_patch(nodes, elements, shapely.Point(400, 200).buffer(150.0))
    # every retained face now shares an edge with another retained face
    sub = elements[~sel.removed]
    e = np.sort(np.vstack([sub[:, [0, 1]], sub[:, [1, 2]], sub[:, [2, 0]]]), axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    owner = np.tile(np.arange(len(sub)), 3)
    order = np.lexsort((e[:, 1], e[:, 0]))
    e2, owner = e[order], owner[order]
    k = np.flatnonzero(np.all(e2[:-1] == e2[1:], axis=1))
    has = np.zeros(len(sub), dtype=bool)
    has[owner[k]] = True
    has[owner[k + 1]] = True
    assert has.all()
    assert sel.report["n_grown_by_repair"] > 0


def test_a_declared_bbox_stays_a_bbox():
    """A square's corners are all equidistant from its centre."""
    from fvcom_mesh_tools.refine import RefineRegion

    box = RefineRegion({"name": "square", "target_h_m": 30,
                        "geometry": {"bbox": [139.78, 35.32, 139.79, 35.32816]}})
    circle = RefineRegion({"name": "disc", "target_h_m": 30,
                           "geometry": {"circle": {"center": [139.78, 35.32],
                                                   "radius_m": 300}}})
    assert box.kind == "bbox" and box.circle is None
    assert circle.kind == "circle"
    assert circle.circle == (139.78, 35.32, 300.0)


def test_depths_follow_a_node_the_repair_moved():
    from fvcom_mesh_tools.patch import refresh_depths

    nodes, elements = grid_mesh(5, 5)
    depths = 5.0 + nodes[:, 0] / 100.0
    moved = np.zeros(len(nodes), dtype=bool)
    moved[12] = True
    after = nodes.copy()
    after[12] = [250.0, 250.0]
    out, n_outside = refresh_depths(nodes, elements, depths, after, moved, depths)
    assert n_outside == 0
    assert out[12] == pytest.approx(7.5)
    assert np.array_equal(out[~moved], depths[~moved])


def test_sliding_really_slides_and_stays_on_the_given_curve():
    """The first version of this test never passed slide_on, so it proved
    nothing: improve_patch disables sliding when no curve is supplied.

    The curve is densified first, because a node sitting ON a curve vertex is
    pinned -- see the corner test below -- and a curve made of exactly the
    mesh's own boundary nodes would pin every one of them."""
    nodes, elements = grid_mesh(11, 11)
    nodes = nodes.copy()
    _u, _c = np.unique(np.sort(np.vstack([elements[:, [0, 1]], elements[:, [1, 2]],
                                          elements[:, [2, 0]]]), axis=1),
                       axis=0, return_counts=True)
    b = _u[_c == 1]
    # Bunch the boundary nodes up ALONG the rectangle, leaving the outline
    # exactly where it was: that is a spacing a slide can improve and a move
    # off the curve cannot.
    rng = np.random.default_rng(11)
    lo, hi = 0.0, nodes.max()
    for v in np.unique(b):
        for ax in (0, 1):
            if lo < nodes[v, ax] < hi:
                nodes[v, ax] += rng.uniform(-35.0, 35.0)
    # The curve is the rectangle itself, whose only vertices are its four
    # corners: the nodes sit between them, as a resampled node sits between
    # two vertices of a source shoreline.
    side = nodes.max()
    curves = [shapely.LineString([(0, 0), (side, 0), (side, side),
                                  (0, side), (0, 0)])]
    on_b = np.zeros(len(nodes), dtype=bool)
    on_b[np.unique(b)] = True
    out, _, _ = improve_patch(nodes, elements, np.zeros(len(nodes), dtype=bool),
                              np.ones(len(elements), dtype=bool),
                              slidable=on_b, slide_on=curves)
    moved = np.linalg.norm(out - nodes, axis=1) > 1e-9
    assert not moved[~on_b].any()
    assert moved.any(), "no node slid, so this test would prove nothing"
    every = shapely.MultiLineString([np.asarray(c.coords) for c in curves])
    assert float(shapely.distance(shapely.points(out[moved]), every).max()) < 1e-6
    # and the four corners, being vertices of the curve, are pinned
    for corner in ((0, 0), (side, 0), (side, side), (0, side)):
        v = int(np.argmin(np.linalg.norm(nodes - corner, axis=1)))
        assert np.array_equal(out[v], nodes[v])


def test_a_slide_cannot_cut_a_corner():
    """Staying on the curve is not the same as leaving the curve where it was.

    A node slid past a vertex is still exactly on the curve, and the polyline
    has lost the corner: measured 0 m off the curve and 52.5 m of Hausdorff
    movement in the boundary itself. Curve vertices are pinned and every
    other node is confined to its own span, so the polyline is invariant.
    """
    from fvcom_mesh_tools.patch import _edge_table

    nodes, elements = grid_mesh(5, 5)
    nodes = nodes.copy()
    nodes[22, 1] = 500.0  # a real corner on the boundary
    u, c = _edge_table(elements)
    rings, _ = boundary_rings(nodes, u[c == 1])
    curve = shapely.LineString(nodes[np.append(rings[0], rings[0][0])])
    slidable = np.zeros(len(nodes), dtype=bool)
    slidable[22] = True
    out, out_t, _ = improve_patch(
        nodes, elements, np.zeros(len(nodes), dtype=bool),
        np.zeros(len(elements), dtype=bool),   # no flips: isolate the slide
        slidable=slidable, slide_on=[curve], rounds=1)
    assert np.array_equal(out[22], nodes[22])
    u2, c2 = _edge_table(out_t)
    rings2, _ = boundary_rings(out, u2[c2 == 1])
    after = shapely.LineString(out[np.append(rings2[0], rings2[0][0])])
    assert curve.hausdorff_distance(after) < 1e-9


# ------------------------------------- what the SECOND adversarial review found
#
# gpt-6-astra reviewed commit 8caca27 (docs/local_refine_implementation_review_2.md)
# and supplied nine more failing assertions. All nine reproduced.


def test_a_missing_patch_face_is_rejected():
    """Every other invariant survives a hole in the water intact.

    Deleting one interior patch triangle left the frozen zone exact, no
    retained face missing, no interface split, no extra face, no non-manifold
    edge, no inversion and no orphan -- and 5,000 m2 of water gone. The
    boundary is what tells, and only if the caller says what it should be.
    """
    from fvcom_mesh_tools.patch import boundary_after_patch

    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    foot = shapely.Point(400, 400).buffer(150.0)
    sel, rc, pn, pt = replay_patch(nodes, elements, foot)
    out, oute, outd, nm, st = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                           rc["pfix"], rc["pfix_base"])
    want = boundary_after_patch(elements, sel, rc, nm, st["pfix_new"])
    good = verify_patch(nodes, depths, elements, sel, out, oute, outd, nm,
                        expected_boundary=want)
    assert good["ok"] and good["boundary_checked"]

    # an interior patch face, all of whose vertices are used elsewhere
    n_ret = len(sel.retained)
    interior = next(k for k in range(n_ret, len(oute))
                    if ((oute == oute[k][0]).sum() > 1
                        and (oute == oute[k][1]).sum() > 1
                        and (oute == oute[k][2]).sum() > 1))
    holed = np.delete(oute, interior, axis=0)
    bad = verify_patch(nodes, depths, elements, sel, out, holed, outd, nm,
                       expected_boundary=want)
    assert not bad["ok"]
    assert bad["n_unexpected_boundary_edges"] > 0


def test_verify_is_not_ok_when_it_did_not_check_coverage():
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    out, oute, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                          rc["pfix"], rc["pfix_base"])
    blind = verify_patch(nodes, depths, elements, sel, out, oute, outd, nm)
    assert not blind["boundary_checked"]
    assert not blind["ok"], "an unchecked claim is not a satisfied one"


def test_an_island_taken_whole_carries_its_own_curve():
    """Otherwise it is measured against somebody else's coastline.

    A 15x15 grid with one cell removed as land and a cut that takes the
    island and part of the mainland: the island had no curve at all, so the
    driver measured its unmoved vertices 700 m from a mainland stretch and
    rejected identity geometry.
    """
    nodes, elements = grid_mesh(15, 15)
    cen = nodes[elements].mean(axis=1)
    land = ((cen[:, 0] > 600) & (cen[:, 0] < 700)
            & (cen[:, 1] > 600) & (cen[:, 1] < 700))
    elements = elements[~land]
    sel = select_patch(nodes, elements,
                       shapely.box(-1.0, 400.0, 850.0, 900.0))
    rc = rim_constraints(nodes, sel, size=1e9, coastline="preserve")
    assert sum(sel.ring_is_hole) == 1
    assert len(rc["curves"]) >= 2, "the island ring registered no curve"
    curves = shapely.MultiLineString([np.asarray(c) for c in rc["curves"]])
    new = rc["pfix"][rc["pfix_base"] < 0]
    assert float(shapely.distance(shapely.points(new), curves).max()) < 1e-9


def test_the_source_is_chosen_by_the_whole_stretch():
    """A line that touches one endpoint and leaves wins a closest-pair test."""
    pts = np.array([[0.0, 0.0], [5.0, 0.0], [10.0, 0.0]])
    alongside = shapely.LineString([(0, 0.1), (5, 0.1), (10, 0.1)])
    diagonal = shapely.LineString([(0, 0), (10, -8)])
    out = coastline_points(pts, 1.0, mode="resample",
                           shoreline=[diagonal, alongside], tolerance_m=1.0)
    assert float(np.abs(out[:, 1]).max()) < 0.5


def test_a_sole_retained_face_is_not_a_spike():
    """It has no edge-neighbour because there is nothing left to be one."""
    nodes = np.array([[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]])
    elements = np.array([[0, 1, 2], [0, 2, 3]])
    sel = select_patch(nodes, elements, shapely.Point(75.0, 25.0).buffer(30.0))
    assert sel.n_removed == 1
    assert len(sel.retained) == 1


def test_rings_that_touch_are_refused_rather_than_unioned():
    """A valid Polygon is not evidence that it is the domain the rim asked for."""
    xy = np.array([[0, 0], [2, 0], [2, 2], [0, 2],
                   [2, 0], [4, 0], [4, 2], [2, 2]], dtype=float)
    edges = np.array([[0, 1], [1, 2], [2, 3], [3, 0],
                      [4, 5], [5, 6], [6, 7], [7, 4]])
    with pytest.raises(ValueError, match="rings touch or overlap"):
        hole_polygon(xy, edges)


def test_a_new_edge_is_not_bounded_by_the_base_r_factor():
    """The same-triangle bound does not survive a change of connectivity.

    Every base edge here is within r = 0.2 and a new edge between two
    interpolated depths is at 0.349. The earlier docstring promised the
    bound globally; this pins the counterexample so it cannot come back.
    """
    from fvcom_mesh_tools.refine import depths_from_base

    nodes, elements = grid_mesh(3, 2, 1.0)
    depths = 2.0 * 1.5 ** nodes[:, 0]
    u, _ = _edge_table_for_test(elements)
    base_r = (np.abs(depths[u[:, 0]] - depths[u[:, 1]])
              / (depths[u[:, 0]] + depths[u[:, 1]])).max()
    new, outside = depths_from_base(nodes, elements, depths,
                                    np.array([[0.1, 0.5], [1.9, 0.5]]))
    assert outside == 0
    assert base_r <= 0.2
    assert abs(new[0] - new[1]) / new.sum() > 0.34


def _edge_table_for_test(elements):
    from fvcom_mesh_tools.patch import _edge_table

    return _edge_table(elements)


def test_a_patch_is_not_blamed_for_the_base_mesh_s_own_failures():
    """The goto2023 production mesh fails C1 at one element 18 km from Futtsu.

    The contract freezes that element, so an absolute gate made every seed
    report "QA 20/21" for a defect the patch is forbidden to touch. What the
    patch is answerable for is a violation involving something it made.
    """
    from types import SimpleNamespace

    from fvcom_mesh_tools.patch import introduced_violations

    n_ret = 100
    elements = np.arange(3 * 150).reshape(150, 3) % 90
    checks = [
        SimpleNamespace(check_id="c1_min_angle", status="fail",
                        requirement=">= 30 deg", observed="min = 28.99",
                        offenders=[{"kind": "element", "id": 7},      # retained
                                   {"kind": "element", "id": 120}]),  # patch
        SimpleNamespace(check_id="c4_area_change", status="fail",
                        requirement="<= 0.5", observed="max = 0.6",
                        offenders=[{"kind": "edge", "elements": [3, 9]},
                                   {"kind": "edge", "elements": [9, 130]}]),
        SimpleNamespace(check_id="c2_max_angle", status="pass",
                        requirement="", observed="", offenders=[]),
    ]
    new = introduced_violations(checks, n_ret, elements)
    assert [v["check"] for v in new] == ["c1_min_angle", "c4_area_change"]
    assert new[0]["id"] == 120
    assert new[1]["elements"] == [9, 130]


def test_a_node_offender_needs_the_connectivity_to_be_attributed():
    from types import SimpleNamespace

    from fvcom_mesh_tools.patch import introduced_violations

    elements = np.array([[0, 1, 2], [1, 2, 3], [2, 3, 4]])
    checks = [SimpleNamespace(check_id="c5_valence", status="fail",
                              requirement="<= 8", observed="max = 9",
                              offenders=[{"kind": "node", "id": 0}])]
    # node 0 is only in element 0, which is retained -> inherited
    assert introduced_violations(checks, 1, elements) == []
    # without the connectivity it cannot be attributed, and is the patch's
    assert len(introduced_violations(checks, 1)) == 1


# ------------------------------------------------- several regions at once
#
# Overlapping fisheries are an ordinary input: two rights over the same
# water. In a refinement they do not conflict, because a target is a ceiling.


def test_a_target_is_a_ceiling_where_regions_overlap():
    """The finest wins, so every region gets at least what it declared."""
    nodes, elements = grid_mesh(15, 15)
    fine = shapely.Point(600, 700).buffer(150.0)
    coarse = shapely.Point(800, 700).buffer(150.0)
    fh = patch_sizing(nodes, elements,
                      [(fine, 20.0, 400.0), (coarse, 60.0, 400.0)],
                      distmesh_scale=1.0)
    shared = np.array([[700.0, 700.0]])          # inside both
    assert fine.contains(shapely.Point(shared[0]))
    assert coarse.contains(shapely.Point(shared[0]))
    assert fh(shared)[0] == pytest.approx(20.0)
    # In the coarse region but outside the fine one: no coarser than its own
    # 60 m ceiling, and finer than that because the fine region's transition
    # reaches it. A ceiling is satisfied by anything below it.
    outside_fine = fh(np.array([[850.0, 700.0]]))[0]
    assert 20.0 < outside_fine <= 60.0


def test_priority_does_not_coarsen_a_patch():
    """It did, and measuring the field is what showed that to be wrong.

    A core imposing a size the surrounding ramp disagrees with is a well or
    a step: with a coarse core winning by priority the measured field slope
    ran 0.50 to 1.93 against a C4 reference of 0.414, and the fill could not
    mesh it.
    """
    nodes, elements = grid_mesh(15, 15)
    fine = shapely.Point(600, 700).buffer(150.0)
    coarse = shapely.Point(800, 700).buffer(150.0)
    shared = np.array([[700.0, 700.0]])
    with_priority = patch_sizing(
        nodes, elements,
        [(fine, 20.0, 400.0, 0.0), (coarse, 60.0, 400.0, 99.0)],
        distmesh_scale=1.0)
    assert with_priority(shared)[0] == pytest.approx(20.0)


def test_region_conflicts_reports_the_shared_water():
    from fvcom_mesh_tools.patch import region_conflicts

    fine = shapely.Point(600, 700).buffer(150.0)
    coarse = shapely.Point(800, 700).buffer(150.0)
    apart = shapely.Point(2000, 2000).buffer(100.0)
    rep = region_conflicts([(fine, 20.0, 400.0, 0.0),
                            (coarse, 60.0, 400.0, 1.0),
                            (apart, 30.0, 400.0, 0.0)],
                           ["fine", "coarse", "apart"])
    assert rep["any_overlap"]
    assert len(rep["overlapping_pairs"]) == 1
    pair = rep["overlapping_pairs"][0]
    assert sorted(pair["regions"]) == ["coarse", "fine"]
    assert pair["effective_target_h_m"] == 20.0
    assert rep["finer_than_declared"]["coarse"]["gets_h_m"] == 20.0
    assert "fine" not in rep["finer_than_declared"]
    assert rep["priority_ignored"], "differing priorities must be called out"


def test_regions_that_do_not_touch_do_not_conflict():
    from fvcom_mesh_tools.patch import region_conflicts

    a = shapely.Point(0, 0).buffer(100.0)
    b = shapely.Point(1000, 0).buffer(100.0)
    rep = region_conflicts([(a, 20.0, 400.0), (b, 30.0, 400.0)], ["a", "b"])
    assert not rep["any_overlap"]
    assert rep["finer_than_declared"] == {}
    assert not rep["priority_ignored"]


def test_the_field_slope_is_measured_not_derived():
    """effective_gradation is a per-region formula; this is the field."""
    from fvcom_mesh_tools.patch import field_gradation

    nodes, elements = grid_mesh(15, 15)
    core = shapely.Point(700, 700).buffer(120.0)
    gentle = patch_sizing(nodes, elements, [(core, 30.0, 900.0)],
                          distmesh_scale=1.0)
    steep = patch_sizing(nodes, elements, [(core, 30.0, 60.0)],
                         distmesh_scale=1.0)
    foot = core.buffer(900.0)
    g = field_gradation(gentle, foot, spacing=20.0)
    t = field_gradation(steep, foot, spacing=20.0)
    assert g["max_slope"] < g["c4_reference_gradation"]
    assert t["max_slope"] > g["max_slope"]
    assert t["fraction_above_reference"] > 0


# --------------------------------------------- what the FOURTH review found
#
# gpt-6-astra's fourth pass ran out of model capacity before writing a
# report, but not before leaving twelve failing probes. All twelve
# reproduced. These are the regressions.


def _check(check_id="c1_min_angle", n=1, offenders=None):
    from types import SimpleNamespace

    return SimpleNamespace(check_id=check_id, status="fail", requirement=">= 30",
                           observed=f"{n} failures", n_violations=n,
                           offenders=offenders if offenders is not None else [])


def test_a_failing_check_that_names_nobody_is_the_patch_s():
    """The most dangerous shape of this bug: QA fails, nothing is attributed,
    and "0 introduced" accepts the mesh.

    `node_index_valid` fails with `offenders=[]`, and the attribution loop
    then yields nothing at all.
    """
    from fvcom_mesh_tools.patch import introduced_violations

    out = introduced_violations([_check("node_index_valid", n=1)], 10)
    assert len(out) == 1
    assert out[0]["kind"] == "unattributed"


def test_violations_a_check_counted_but_did_not_name_are_the_patch_s():
    """run_qa truncates its offender list; the unlisted ones are unproven."""
    from fvcom_mesh_tools.patch import introduced_violations

    listed = [{"kind": "element", "id": i} for i in range(3)]
    out = introduced_violations([_check(n=10, offenders=listed)], 10)
    assert [v["kind"] for v in out] == ["unattributed"]
    assert out[0]["n_unattributed"] == 7


def test_an_offender_this_cannot_place_is_not_waved_through():
    """An id that is a pair used to raise TypeError inside the attribution."""
    from fvcom_mesh_tools.patch import introduced_violations

    tri = np.array([[0, 1, 2], [0, 2, 3]])
    out = introduced_violations(
        [_check("no_duplicate_nodes", n=1,
                offenders=[{"kind": "node_pair", "id": [0, 1]}])], 2, tri)
    assert len(out) == 1


def test_limit_rfactor_accepts_a_flat_bottom():
    """`rfactor_limit: base` passes the base's own worst r, which is 0 for a
    constant-depth base, and that used to raise."""
    from fvcom_mesh_tools.refine import limit_rfactor

    elements = np.array([[0, 1, 2], [0, 2, 3]])
    depths = np.full(4, 8.0)
    out, info = limit_rfactor(elements, depths, np.array([0, 1, 1, 0], dtype=bool), 0.0)
    assert info["converged"]
    assert np.array_equal(out, depths)


def test_an_unfixable_r_factor_edge_is_not_called_converged():
    """A new edge between two frozen nodes cannot be fixed here, and saying
    so is the caller's only chance to notice it."""
    from fvcom_mesh_tools.refine import limit_rfactor

    depths = np.array([3.0, 4.5, 6.75, 4.5])
    after = np.array([[0, 1, 2], [0, 2, 3]])   # new connectivity: 0-2 is new
    _, info = limit_rfactor(after, depths, np.zeros(4, dtype=bool), 0.2,
                            depth_min=3.0, depth_max=6.75)
    assert not info["converged"]
    assert info["n_over_rmax_frozen_pair"] == 1
    assert info["n_over_rmax_movable"] == 0


def test_field_gradation_measures_the_gradient_not_the_axes():
    """A ramp rising equally in x and y has slope sqrt(2) * 0.35, and axis
    differences report 0.35 -- a factor sqrt(2) gentler than it is."""
    from fvcom_mesh_tools.patch import field_gradation

    box = shapely.box(0, 0, 100, 100)
    rep = field_gradation(lambda p: 100 + 0.35 * np.asarray(p)[:, 0]
                          + 0.35 * np.asarray(p)[:, 1], box, spacing=10.0)
    assert rep["max_slope"] == pytest.approx(np.hypot(0.35, 0.35), rel=1e-6)


def test_a_region_swallowed_by_another_s_transition_is_reported():
    """Only core-to-core overlap was looked at, so a 5 m core 200 m away with
    a 1 km transition swallowed a 90 m core and nothing was said.

    A ramp runs from its target to the BASE mesh's size, so what the
    neighbour's transition is worth here is a fact about the base mesh: with
    it the claim is the field's, and without it there is no claim to make.
    """
    from fvcom_mesh_tools.patch import region_conflicts

    fine = shapely.Point(500, 700).buffer(50.0)
    coarse = shapely.Point(700, 700).buffer(30.0)
    regions = [(fine, 5.0, 1000.0), (coarse, 90.0, 100.0)]
    rep = region_conflicts(regions, ["fine", "coarse"],
                           base_size=lambda p: np.full(len(p), 300.0))
    assert not rep["any_overlap"], "their cores do not touch"
    assert rep["transitions_evaluated"]
    assert "coarse" in rep["finer_than_declared"]
    assert rep["finer_than_declared"]["coarse"]["gets_h_m"] < 90.0
    assert rep["finer_than_declared"]["coarse"]["by_core_overlap"] is False
    assert rep["finer_than_declared"]["coarse"]["because_of"] == ["fine"]

    blind = region_conflicts(regions, ["fine", "coarse"])
    assert not blind["transitions_evaluated"]
    assert blind["finer_than_declared"] == {}, (
        "without the base mesh a ramp's value is unknown, and an unknown is "
        "not a claim")


def test_the_conflict_report_is_the_field_not_a_formula():
    """A 12 m region 480 m from a 5 m one keeps its 12 m in the field the
    mesher is handed; a report with its own ramp formula said 9.6 m."""
    from fvcom_mesh_tools.patch import (
        base_size_field,
        patch_sizing,
        region_conflicts,
    )

    nodes, elements = grid_mesh(15, 15)
    regions = [(shapely.box(490, 690, 510, 710), 5.0, 1000.0),
               (shapely.box(990, 690, 1010, 710), 12.0, 100.0)]
    here = np.array([[1000.0, 700.0]])
    alone = patch_sizing(nodes, elements, [regions[1]], distmesh_scale=1.0)(here)
    joint = patch_sizing(nodes, elements, regions, distmesh_scale=1.0)(here)
    assert alone[0] == pytest.approx(12.0) and joint[0] == pytest.approx(12.0)
    rep = region_conflicts(regions, ["fine", "coarse"],
                           base_size=base_size_field(nodes, elements))
    assert "coarse" not in rep["finer_than_declared"]


# ------------------------------------------------- attribution, fourth review


def test_attribution_is_not_capped_with_the_display_list():
    """A cap on what is PRINTED is not a cap on what is attributed.

    An unchanged 73x73 mesh of 30-30-120 triangles fails C1 10,368 times.
    Reading the display list, capped at 10,000, turned the remaining 368 into
    "unattributed" violations and rejected a mesh with no patch in it at all.
    Raising the cap only moves the failure, so the two lists are separate.
    """
    from fvcom_mesh_tools.io.fort14 import Fort14Mesh
    from fvcom_mesh_tools.patch import introduced_violations
    from fvcom_mesh_tools.qa import run_qa

    nodes, elements = grid_mesh(73, 73)
    nodes[:, 1] *= 0.2                       # thin triangles: C1 fails everywhere
    mesh = Fort14Mesh(title="flat", nodes=nodes, depths=np.full(len(nodes), 8.0),
                      elements=elements, open_boundaries=[], land_boundaries=[])
    qa = run_qa(mesh, max_offenders=10_000)
    c1 = next(c for c in qa.checks if c.check_id == "c1_min_angle")
    assert c1.n_violations > len(c1.offenders), "this needs a truncated list"
    assert len(c1.offender_ids) == c1.n_violations
    assert introduced_violations([c1], len(elements), elements) == [], (
        "every element is retained, so the patch introduced nothing")


def test_attribution_still_refuses_a_check_that_named_nobody():
    """The cap fix must not restore the empty-list acceptance hole."""
    from types import SimpleNamespace

    from fvcom_mesh_tools.patch import introduced_violations

    check = SimpleNamespace(check_id="node_index_valid", status="fail",
                            requirement="in range", observed="3 bad",
                            n_violations=3, offenders=[], offender_ids=[])
    out = introduced_violations([check], 100, np.zeros((100, 3), dtype=int))
    assert [v["kind"] for v in out] == ["unattributed"]
    assert out[0]["n_unattributed"] == 3


@pytest.mark.parametrize("ident", ["0", 0.5, -1, True, 10_000])
def test_an_id_that_cannot_be_looked_up_does_not_exonerate(ident):
    """`0.5 < n_retained` is True and meant nothing; `-1` indexes from the end."""
    from types import SimpleNamespace

    from fvcom_mesh_tools.patch import introduced_violations

    nodes, elements = grid_mesh(5, 5)
    check = SimpleNamespace(check_id="c1_min_angle", status="fail",
                            requirement=">= 30", observed="bad", n_violations=1,
                            offenders=[{"kind": "element", "id": ident}],
                            offender_ids=[{"kind": "element", "id": ident}])
    assert introduced_violations([check], 2, elements), (
        f"id {ident!r} is not an element of this mesh and cannot clear the patch")


def test_a_field_thinner_than_the_lattice_is_unknown_not_flat():
    """Masking before differencing left one row with NaN neighbours.

    np.gradient across it is NaN, the finite filter dropped every sample, and
    a field rising 1 m per metre reported a slope of zero.
    """
    from fvcom_mesh_tools.patch import field_gradation

    rep = field_gradation(lambda p: 100.0 + np.asarray(p)[:, 0],
                          shapely.box(0, 0, 100, 20), spacing=10.0)
    assert rep["max_slope"] == pytest.approx(1.0)
    # One row of samples has no neighbour across the box, so the lattice is
    # refined until it has one rather than reporting a lower bound as the
    # gradient (fifth review). At the requested spacing there were nine.
    assert rep["spacing_m"] < 10.0 and rep["n_samples"] > 9
    coarse = field_gradation(lambda p: 100.0 + np.asarray(p)[:, 0],
                             shapely.box(0, 0, 100, 20), spacing=10.0,
                             refine=0)
    assert coarse["n_samples"] == 9
    assert coarse["max_slope"] == pytest.approx(1.0), (
        "even without refinement the one-sided bound must not read as zero")

    # A footprint thinner than the lattice is now refined until it fits
    # (fifth review); what must never happen is a confident zero.
    thin = field_gradation(lambda p: 100.0 + np.asarray(p)[:, 0],
                           shapely.box(0, 0, 100, 1), spacing=10.0)
    assert thin["measured"] and thin["spacing_m"] < 10.0
    assert thin["max_slope"] == pytest.approx(1.0)

    blind = field_gradation(lambda p: 100.0 + np.asarray(p)[:, 0],
                            shapely.box(0, 0, 100, 1), spacing=10.0, refine=0)
    assert blind["measured"] is False
    assert blind["max_slope"] is None, "an unmeasured slope is not a gentle one"


# ------------------------------------------------- resolution, fifth review


def _wedge_mesh(fine_h=30.0, coarse_h=300.0, n=40):
    """Two square patches side by side: the left fine, the right coarse."""
    left, _ = grid_mesh(n, n, fine_h)
    right, _ = grid_mesh(5, 5, coarse_h)
    right = right + np.array([fine_h * (n - 1) + coarse_h, 0.0])
    nodes = np.vstack([left, right])
    tri = []
    for block, start, nx, ny in ((left, 0, n, n),
                                 (right, len(left), 5, 5)):
        for j in range(ny - 1):
            for i in range(nx - 1):
                a = start + j * nx + i
                tri += [[a, a + 1, a + nx + 1], [a, a + nx + 1, a + nx]]
    return nodes, np.asarray(tri, dtype=np.int64)


def test_resolution_is_measured_over_the_area_not_over_the_edges():
    """A request half refined and half untouched passed the edge median.

    The fine half owns almost every edge inside the region, so the median
    edge length is the fine half's. On the delivered Futtsu mesh a request
    with 46 % of its water at the base size read as delivered (fifth
    review). Asking the water instead of the edges sees it.
    """
    from fvcom_mesh_tools.patch import region_resolution

    nodes, elements = _wedge_mesh()
    fine = shapely.box(100.0, 100.0, 1000.0, 1000.0)
    both = shapely.box(100.0, 100.0, 2000.0, 1000.0)

    good = region_resolution(nodes, elements, fine, 30.0)
    # A right-triangulated square of side h has the area of half a square,
    # so its equivalent edge is 1.07 h: this grid is not an equilateral fill.
    assert good["median_ratio"] == pytest.approx(1.07, abs=0.02)
    assert good["covered_fraction"] == pytest.approx(1.0)

    half = region_resolution(nodes, elements, both, 30.0)
    assert half["covered_fraction"] < 0.95, (
        "half of this request is at the base size and must not read as covered")

    # The edge statistic is the one that was fooled, and it still is: it is
    # reported beside the area measure rather than gated on.
    e = np.unique(np.sort(np.vstack([elements[:, [0, 1]], elements[:, [1, 2]],
                                     elements[:, [2, 0]]]), axis=1), axis=0)
    mid = 0.5 * (nodes[e[:, 0]] + nodes[e[:, 1]])
    length = np.linalg.norm(nodes[e[:, 0]] - nodes[e[:, 1]], axis=1)
    inside = np.asarray(shapely.contains(both, shapely.points(mid[:, 0],
                                                              mid[:, 1])))
    assert np.median(length[inside]) < 1.05 * 30.0


def test_a_region_outside_the_mesh_is_reported_and_is_a_miss():
    from fvcom_mesh_tools.patch import region_resolution

    nodes, elements = grid_mesh(9, 9, 30.0)
    away = shapely.box(10_000.0, 10_000.0, 10_500.0, 10_500.0)
    rep = region_resolution(nodes, elements, away, 30.0)
    assert rep["outside_mesh_fraction"] == pytest.approx(1.0)
    assert rep["covered_fraction"] == 0.0
    assert rep["median_ratio"] is None, "water that does not exist has no size"

    half_out = shapely.box(100.0, 100.0, 400.0, 10_000.0)
    rep = region_resolution(nodes, elements, half_out, 30.0)
    assert 0.0 < rep["outside_mesh_fraction"] < 1.0
    assert rep["median_ratio"] is not None


def test_the_conflict_area_of_a_thin_region_is_not_a_single_point():
    """One representative point was weighted as the whole polygon.

    A 100 x 10 m rectangle inside a 10 km x 10 m one is 1 % of it; the
    area-derived lattice missed the shape entirely and the answer came back
    as 0 % or 100 % depending on where the thin part sat (fifth review).
    """
    from fvcom_mesh_tools.patch import region_conflicts

    long_thin = shapely.box(0.0, 0.0, 10_000.0, 10.0)
    for x in (0.0, 4950.0):
        fine = shapely.box(x, 0.0, x + 100.0, 10.0)
        rep = region_conflicts([(long_thin, 30.0, 0.0), (fine, 5.0, 0.0)],
                               ["coarse", "fine"],
                               base_size=lambda p: np.full(len(p), 300.0))
        got = rep["finer_than_declared"]["coarse"]
        assert got["fraction"] == pytest.approx(0.01, abs=0.005), (
            f"the fine region covers 1 % of the coarse one, not "
            f"{got['fraction']}")
        assert got["area_is_estimated"] and got["area_sampled"]


def test_a_transverse_ramp_in_a_narrow_channel_is_not_flat():
    """Neither axis had support across a 40 m channel at 25 m spacing, so a
    field rising 1 m per metre across it measured as zero (fifth review)."""
    from fvcom_mesh_tools.patch import field_gradation

    channel = shapely.box(0.0, 0.0, 1000.0, 40.0)
    rep = field_gradation(lambda p: 100.0 + np.asarray(p)[:, 1], channel)
    assert rep["measured"]
    assert rep["max_slope"] == pytest.approx(1.0, rel=1e-6)
    assert rep["spacing_m"] < 25.0, "the lattice has to be refined to see it"


def test_an_offender_that_cannot_be_hashed_is_unplaceable_not_a_crash():
    """A dict id raised TypeError before it could be refused."""
    from types import SimpleNamespace

    from fvcom_mesh_tools.patch import introduced_violations

    _, elements = grid_mesh(5, 5)
    off = {"kind": "element", "id": {"index": 0}}
    check = SimpleNamespace(check_id="c1_min_angle", status="fail",
                            requirement=">= 30", observed="bad",
                            n_violations=1, offenders=[off],
                            offender_ids=[off])
    assert introduced_violations([check], 2, elements)


def test_a_fallback_identity_keeps_what_places_the_offender():
    """`kind` and `id` alone cannot place a C4 edge; its elements can.

    The fallback dropped them and turned a wholly inherited edge into an
    introduced violation, with no truncation involved (fifth review).
    """
    from fvcom_mesh_tools.patch import introduced_violations
    from fvcom_mesh_tools.qa import QACheck

    _, elements = grid_mesh(5, 5)
    check = QACheck("c4_area_change", "quality", True, False, "<= 0.5",
                    "max = 0.6", 1,
                    offenders=[{"kind": "edge", "id": [0, 2],
                                "elements": [0, 1]}])
    assert check.offender_ids[0]["elements"] == [0, 1]
    assert introduced_violations([check], len(elements), elements) == []


def test_two_obc_segments_with_the_same_pair_do_not_collide():
    """The decoration key omitted the segment, so both records named one."""
    from fvcom_mesh_tools.patch import _offender_key

    a = {"kind": "obc_pair", "id": [1, 3], "segment": 0}
    b = {"kind": "obc_pair", "id": [1, 3], "segment": 1}
    assert _offender_key(a) != _offender_key(b)


# --- coastline: resolve, the hires fork -------------------------------------

def _wiggly_source(n=400):
    """A source shoreline with detail far finer than any element here."""
    import shapely

    t = np.linspace(0.0, 1000.0, n)
    return shapely.LineString(np.column_stack([t, 20.0 * np.sin(t / 25.0)]))


def test_resolve_follows_the_source_past_the_departure_veto():
    """`resample` refuses a departure over the tolerance; `resolve` is the branch
    where the coastline is SUPPOSED to move, so it has no veto at all."""
    from fvcom_mesh_tools.patch import coastline_points

    src = _wiggly_source()
    base = np.array([[0.0, 0.0], [500.0, 0.0], [1000.0, 0.0]])
    with pytest.raises(ValueError, match="beyond the"):
        coastline_points(base, 25.0, mode="resample", shoreline=src, tolerance_m=5.0)
    out = coastline_points(base, 25.0, mode="resolve", shoreline=src, tolerance_m=5.0)
    assert len(out) > len(base)
    assert np.abs(out[1:-1, 1]).max() > 5.0, "it must have left the base polyline"
    assert np.allclose(out[0], base[0]) and np.allclose(out[-1], base[-1]), (
        "the frozen anchors are still frozen")


def test_resolve_refuses_a_stretch_the_source_does_not_cover():
    """Silently subdividing the base is how a fidelity check certifies itself."""
    import shapely

    from fvcom_mesh_tools.patch import coastline_curve, coastline_points

    # Both endpoints of the stretch project to the SAME station on the source,
    # so there is no substring between them and _source_substring returns None.
    source = shapely.LineString([[0.0, 0.0], [1000.0, 0.0]])
    base = np.array([[500.0, 100.0], [500.0, 200.0]])
    with pytest.raises(ValueError, match="no source component"):
        coastline_points(base, 25.0, mode="resolve", shoreline=source)
    with pytest.raises(ValueError, match="no source component"):
        coastline_curve(base, "resolve", source)
    # the same stretch on the `resample` branch falls back to the base, which
    # is the behaviour this branch refuses rather than a bug in that one
    out = coastline_points(base, 25.0, mode="resample", shoreline=source,
                           tolerance_m=1e9)
    assert np.allclose(out[0], base[0]) and np.allclose(out[-1], base[-1])


def test_resolve_keeps_core_detail_that_the_median_rule_destroys():
    """The simplifier's tolerance is a quarter of the MEDIAN size over the whole
    substring, so a stretch fine in the core and coarse in the transition loses
    its core detail at the transition's scale (review P2-11)."""
    from fvcom_mesh_tools.patch import coastline_points

    src = _wiggly_source()
    base = np.array([[0.0, 0.0], [1000.0, 0.0]])

    def size(q):
        # 10 m over the first fifth, 400 m after it
        return np.where(np.asarray(q)[:, 0] < 200.0, 10.0, 400.0)

    coarse = coastline_points(base, size, mode="resample", shoreline=src,
                              tolerance_m=1e9)
    fine = coastline_points(base, size, mode="resolve", shoreline=src,
                            tolerance_m=1e9)
    in_core = lambda p: p[(p[:, 0] > 0) & (p[:, 0] < 200.0)]      # noqa: E731
    assert len(in_core(fine)) > len(in_core(coarse)), (
        "the pointwise tolerance has to keep what the median rule threw away")


def test_resolve_follows_the_source_only_where_the_mesh_is_fine_enough():
    """Measured twice, in opposite directions, on the real patches.

    A 300 m fishery 2 km offshore put its whole coastline in a 400-1700 m
    field and `resolve` replaced 17 base nodes with 13 -- COARSER than the
    base.  A region on the shore resolved its core cleanly and produced two
    6-degree elements 6 km away at Kimitsu port.  Both are the same thing:
    walking a detailed shoreline at the coarse end of a transition.
    """
    from fvcom_mesh_tools.patch import coastline_points

    src = _wiggly_source()
    # base vertices 100 m apart, so 10 m of mesh is finer than the base and
    # 400 m is coarser
    base = np.column_stack([np.arange(0.0, 1001.0, 100.0), np.zeros(11)])

    def size(q):
        return np.where(np.asarray(q)[:, 0] < 500.0, 10.0, 400.0)

    out = coastline_points(base, size, mode="resolve", shoreline=src,
                           tolerance_m=1e9)
    west, east = out[out[:, 0] < 500.0], out[out[:, 0] > 500.0]
    assert np.abs(west[:, 1]).max() > 5.0, (
        "where the mesh is finer than the base, the source must be followed")
    assert np.allclose(east[:, 1], 0.0, atol=1e-9), (
        "where it is coarser, the base polyline must be kept exactly")
    assert len(out) > len(base)


def test_resolve_never_returns_fewer_points_than_the_base_stretch_had():
    """The coarsening that started this: 17 base nodes replaced by 13."""
    from fvcom_mesh_tools.patch import coastline_points

    src = _wiggly_source()
    base = np.column_stack([np.arange(0.0, 1001.0, 50.0), np.zeros(21)])
    out = coastline_points(base, 900.0, mode="resolve", shoreline=src,
                           tolerance_m=1e9)
    assert len(out) >= len(base), (
        f"a 900 m walk returned {len(out)} points for a 21-point base stretch; "
        "subdividing keeps every original vertex and walking does not")


def test_the_size_field_outside_the_mesh_is_a_cliff_or_a_continuation():
    """`resolve` moves the boundary off the base mesh, and what is out there
    decides whether the mesher is handed a step.

    Measured on the real patches: with the field maximum outside, the hires
    runs report a size-field slope of 35-42 against a C4 reference of 0.414,
    where the same patch on the default branch reports 0.375.

    The mesh here has to be GRADED for the question to exist: on a uniform
    grid the field's maximum IS its local value, so both fills agree and the
    first version of this test proved nothing.
    """
    from fvcom_mesh_tools.patch import base_size_field

    # columns 20 m apart on the left, 400 m apart on the right
    xs = np.concatenate([[0.0], np.cumsum(np.linspace(20.0, 400.0, 14))])
    ys = np.linspace(0.0, 400.0, 5)
    gx, gy = np.meshgrid(xs, ys)
    nodes = np.column_stack([gx.ravel(), gy.ravel()])
    nx, ny = len(xs), len(ys)
    tri = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            tri += [[a, a + 1, a + nx + 1], [a, a + nx + 1, a + nx]]
    tri = np.asarray(tri, dtype=np.int64)

    fine = np.array([[10.0, 200.0]])                  # in the 20 m columns
    beyond = np.array([[10.0, 430.0]])                # 30 m past the top edge
    f_max = base_size_field(nodes, tri)
    f_near = base_size_field(nodes, tri, outside="nearest")
    assert np.allclose(f_max(fine), f_near(fine)), (
        "inside the mesh the two must be the same field")
    local = float(f_max(fine)[0])
    step = abs(float(f_max(beyond)[0]) - local)
    smooth = abs(float(f_near(beyond)[0]) - local)
    assert step > 10 * max(smooth, 1e-9), (
        f"the maximum fill should be a cliff here: {step:.1f} m against "
        f"{smooth:.1f} m, on a field whose local value is {local:.1f} m")
    with pytest.raises(ValueError, match="outside must be"):
        base_size_field(nodes, tri, outside="zero")


def test_a_replacement_that_would_cross_the_frozen_boundary_is_refused():
    """A moved boundary can cross the one it does not own.

    Measured: a hires run with a clean size field produced exactly one
    crossing pair -- a retained boundary edge against a resolved one -- and
    matplotlib's TriFinder refused the mesh with "Triangulation is invalid",
    while verify_patch had passed it. Subdividing cannot cross anything,
    because it stays on the base polyline.
    """
    from fvcom_mesh_tools.patch import _crosses_boundary

    xy = np.array([[0.0, 0.0], [100.0, 0.0], [200.0, 0.0],   # the stretch
                   [50.0, -50.0], [50.0, 50.0]])             # an edge across it
    idx = np.array([0, 1, 2])
    other = np.array([[3, 4]])
    own = np.array([[0, 1], [1, 2]])
    straight = xy[idx]
    assert not _crosses_boundary(straight, xy, idx, own), (
        "the stretch's own edges are shared endpoints, not crossings")
    assert _crosses_boundary(straight, xy, idx, np.vstack([own, other]))
    away = np.array([[0.0, 0.0], [100.0, 200.0], [200.0, 0.0]])
    assert not _crosses_boundary(away, xy, idx, np.vstack([own, other]))
    assert not _crosses_boundary(straight, xy, idx, None)
