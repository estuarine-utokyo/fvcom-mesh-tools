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
    with pytest.raises(ValueError, match="do not close"):
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


def test_effective_gradation_flags_a_ramp_c4_cannot_carry():
    nodes, elements = grid_mesh(15, 15, h=100.0)
    core = shapely.Point(700, 700).buffer(100.0)
    easy = effective_gradation(nodes, elements, [(core, 20.0, 2000.0)])
    hard = effective_gradation(nodes, elements, [(core, 20.0, 150.0)])
    assert easy["within_c4"]
    assert not hard["within_c4"]
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
    nodes, elements = grid_mesh()
    depths = np.full(len(nodes), 8.0)
    sel, rc, pn, pt = replay_patch(nodes, elements,
                                   shapely.Point(400, 400).buffer(150.0))
    out, oute, outd, nm, _ = stitch_patch(nodes, elements, depths, sel, pn, pt,
                                          rc["pfix"], rc["pfix_base"])
    ver = verify_patch(nodes, depths, elements, sel, out, oute, outd, nm)
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
