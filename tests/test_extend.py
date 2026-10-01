"""Extension of a base mesh: sizing composition, merge, frozen contract."""

import numpy as np
import pytest

from fvcom_mesh_tools.extend import (
    band_field,
    compose_sizing,
    graded_up,
    merge_outer,
    verify_frozen_base,
)
from fvcom_mesh_tools.io.fort14 import Fort14Mesh


def _lattice(n=41, h=1000.0):
    x, y = np.meshgrid(np.arange(n) * h, np.arange(n) * h)
    return x, y


def _grade_ok(v, x, y, grade):
    for a, b in ((v[:, 1:], v[:, :-1]), (v[1:, :], v[:-1, :])):
        if np.max(np.abs(a - b)) > grade * 1000.0 + 1e-6:
            return False
    return True


def test_graded_up_is_the_smallest_feasible_field_above():
    x, y = _lattice()
    v = np.full(x.shape, 100.0)
    v[20, 20] = 5000.0
    up = graded_up(v, x, y, 0.2)
    assert up[20, 20] == 5000.0 and np.all(up >= v) and _grade_ok(up, x, y, 0.2)
    assert up[20, 25] == pytest.approx(5000.0 - 0.2 * 5000.0)


def test_compose_respects_floor_band_and_gradation():
    x, y = _lattice()
    ambient = 300.0 + 0.3 * y                      # grows away from y = 0
    floor = np.where(x > 30_000, 4000.0, 0.0)      # deep water on the east
    band = band_field(x, y, [[0, 20_000], [40_000, 20_000]], [1500.0, 1500.0], 800.0)
    h, rep = compose_sizing(ambient, x, y, grade=0.2, floor=floor, bands=[band])
    on = np.isfinite(band)
    assert np.allclose(h[on], 1500.0)             # the band is set, up or down
    assert _grade_ok(h, x, y, 0.2 * 1.0000001)
    off_band = ~on & (np.abs(y - 20_000) > 15_000)
    assert np.all(h[off_band] >= floor[off_band] - 1e-6)
    assert rep["band_0_cells"] == int(on.sum())


def test_band_field_interpolates_along_the_line_and_refuses_bad_input():
    x, y = _lattice()
    b = band_field(x, y, [[0, 0], [40_000, 0]], [1000.0, 3000.0], 10.0)
    assert b[0, 0] == pytest.approx(1000.0) and b[0, 40] == pytest.approx(3000.0)
    assert b[0, 20] == pytest.approx(2000.0) and np.isnan(b[5, 5])
    with pytest.raises(ValueError, match="one target per point"):
        band_field(x, y, [[0, 0], [1, 1]], [1.0], 5.0)


def _base():
    """Two triangles over the unit square; the east side (1,2) is the open boundary."""
    nodes = np.array([[0, 0], [1000, 0], [1000, 1000], [0, 1000]], float)
    elems = np.array([[0, 1, 2], [0, 2, 3]])
    return Fort14Mesh("base", nodes, np.array([5.0, 6.0, 7.0, 8.0]), elems,
                      [np.array([1, 2])], [(0, np.array([2, 3, 0, 1]))])


def _outer():
    """Two triangles over the square to the east, sharing x = 1000."""
    nodes = np.array([[1000, 0], [2000, 0], [2000, 1000], [1000, 1000.2]], float)
    elems = np.array([[0, 2, 1], [0, 3, 2]])      # opposite orientation on purpose
    return nodes, elems


def test_merge_keeps_the_base_first_and_turns_outer_elements():
    base = _base()
    on, oe = _outer()
    m = merge_outer(base, on, oe, [0, 3], [1, 2], [1, 2])
    assert m.n_nodes == 6 and m.n_elements == 4
    assert np.array_equal(m.nodes[:4], base.nodes) and np.array_equal(m.elements[:2], base.elements)
    a = m.nodes[m.elements]
    area = 0.5 * ((a[:, 1, 0] - a[:, 0, 0]) * (a[:, 2, 1] - a[:, 0, 1])
                  - (a[:, 2, 0] - a[:, 0, 0]) * (a[:, 1, 1] - a[:, 0, 1]))
    assert np.all(np.sign(area) == np.sign(area[0]))
    assert m.open_boundaries[0].tolist() == [4, 5] and np.isnan(m.depths[4:]).all()
    rep = verify_frozen_base(m, base, [1, 2])
    assert rep["n_interface_edges"] == 1


def test_merge_refuses_interface_nodes_that_do_not_coincide():
    on, oe = _outer()
    on[3, 1] = 1003.0
    with pytest.raises(ValueError, match="differ from the base"):
        merge_outer(_base(), on, oe, [0, 3], [1, 2], [1, 2])


def test_verify_catches_a_moved_base_node_and_an_unshared_interface():
    base = _base()
    on, oe = _outer()
    m = merge_outer(base, on, oe, [0, 3], [1, 2], [1, 2])
    moved = m.nodes.copy()
    moved[0, 0] += 1e-9
    with pytest.raises(ValueError, match="node coordinates"):
        verify_frozen_base(Fort14Mesh("m", moved, m.depths, m.elements, [], []), base, [1, 2])
    with pytest.raises(ValueError, match="not shared"):
        verify_frozen_base(m, base, [0, 3])


def test_land_segments_split_loops_at_the_open_boundary():
    from fvcom_mesh_tools.extend import land_segments

    # a 3x3-node square (8 boundary nodes) with a hole-free interior
    t = []
    for j in range(2):
        for i in range(2):
            a = j * 3 + i
            t += [[a, a + 1, a + 4], [a, a + 4, a + 3]]
    t = np.array(t)
    runs = land_segments(t, [np.array([2, 5, 8])])          # the east side is open
    assert len(runs) == 1 and runs[0][0] == 0
    run = runs[0][1].tolist()
    assert run[0] in (2, 8) and run[-1] in (2, 8) and 5 not in run and len(run) == 7
    islands = land_segments(t, [])
    assert len(islands) == 1 and islands[0][0] == 1 and len(islands[0][1]) == 9


def test_rfactor_smooth_free_leaves_fixed_nodes_alone():
    from fvcom_mesh_tools.extend import rfactor_smooth_free

    # a chain 0-1-2-3: node 0 fixed at 10 m, the rest free and deep
    h0 = np.array([10.0, 200.0, 400.0, 800.0])
    ei, ej = np.array([0, 1, 2]), np.array([1, 2, 3])
    free = np.array([False, True, True, True])
    h, it, r = rfactor_smooth_free(h0, ei, ej, free, rmax=0.2, hmin=3.0)
    assert h[0] == 10.0 and r <= 0.2 + 1e-6
    assert np.all(np.abs(h[ei] - h[ej]) / (h[ei] + h[ej]) <= 0.2 + 1e-6)
    # an edge between two fixed nodes is not touched and not counted
    h, it, r = rfactor_smooth_free(np.array([10.0, 100.0]), np.array([0]), np.array([1]),
                                   np.array([False, False]), rmax=0.2, hmin=3.0)
    assert h.tolist() == [10.0, 100.0] and it == 0


def test_trim_lone_corners_drops_a_cape_tip_but_not_a_kept_or_frozen_one():
    from fvcom_mesh_tools.extend import trim_lone_corners

    # a hexagonal fan around node 0 (ring 1..6) and a spike on edge 1-2
    fan = [[0, i, i % 6 + 1] for i in range(1, 7)]
    t = np.array(fan + [[1, 7, 2]])
    out, mut, rep = trim_lone_corners(t, np.ones(len(t), bool))
    assert rep["n_elements_dropped"] == 1 and not (out == 7).any() and len(out) == 6
    assert rep["lone_nodes_left"] == []
    out, mut, rep = trim_lone_corners(t, np.r_[np.ones(6, bool), False])
    assert len(out) == 7 and rep["lone_nodes_left"] == [7]
    out, mut, rep = trim_lone_corners(t, np.ones(len(t), bool), keep_nodes=[7])
    assert len(out) == 7


def test_trim_lone_corners_leaves_a_lone_corner_that_is_not_a_spike():
    from fvcom_mesh_tools.extend import trim_lone_corners

    # a lone triangle: its sides are all boundary, dropping it would orphan more
    t = np.array([[0, 1, 2]])
    out, mut, rep = trim_lone_corners(t, np.ones(1, bool))
    assert len(out) == 1 and sorted(rep["lone_nodes_left"]) == [0, 1, 2]


def test_rfactor_smooth_free_meets_the_cap_and_refuses_the_infeasible():
    """Review F3: a cap after smoothing broke r; an infeasible limit returned silently."""
    from fvcom_mesh_tools.extend import rfactor_smooth_free

    ei, ej = np.array([0]), np.array([1])
    h, it, r = rfactor_smooth_free(np.array([100.0, 10.0]), ei, ej, np.array([False, True]),
                                   rmax=0.2, hmin=3.0, hmax=80.0)
    assert h[1] <= 80.0 and r <= 0.2 + 1e-6
    with pytest.raises(ValueError, match="not reached"):
        # 100 m fixed beside a free node capped at 10 m: r = 0.82 at best
        rfactor_smooth_free(np.array([100.0, 10.0]), ei, ej, np.array([False, True]),
                            rmax=0.2, hmin=3.0, hmax=10.0, max_iter=50)
    with pytest.raises(ValueError, match="not reached"):
        rfactor_smooth_free(np.array([1.0, 3.0]), ei, ej, np.array([False, True]),
                            rmax=0.2, hmin=3.0, max_iter=50)


def test_compose_reports_a_band_set_below_the_floor():
    """Review F4: a band below the time-step floor must be visible to the caller."""
    x, y = _lattice()
    floor = np.full(x.shape, 4000.0)
    band = band_field(x, y, [[0, 20_000], [40_000, 20_000]], [1500.0, 1500.0], 800.0)
    h, rep = compose_sizing(np.full(x.shape, 5000.0), x, y, grade=0.2, floor=floor,
                            bands=[band])
    assert rep["band_0_below_floor_cells"] == int(np.isfinite(band).sum())


def test_verify_refuses_an_outer_element_overlapping_the_base():
    """Review F10: two owners on the same side of an interface edge is an overlap."""
    base = _base()
    on = np.array([[1000, 0], [1000, 1000], [500, 500]], float)   # third vertex inside the base
    oe = np.array([[0, 1, 2]])
    m = merge_outer(base, on, oe, [0, 1], [1, 2], [0, 1])
    with pytest.raises(ValueError, match="same side"):
        verify_frozen_base(m, base, [1, 2])


def test_land_segments_keeps_a_single_land_edge():
    """Review F11: a land run of one edge was dropped."""
    from fvcom_mesh_tools.extend import land_segments

    runs = land_segments(np.array([[0, 1, 2]]), [np.array([0, 1, 2])])
    assert len(runs) == 1 and sorted(runs[0][1].tolist()) == [0, 2]


def test_rfactor_limiter_round2_edges():
    """Round 2 F14: convergence on the last pass; F15: NaN depths and bad controls."""
    from fvcom_mesh_tools.extend import rfactor_smooth_free

    h, it, r = rfactor_smooth_free(np.array([10.0, 100.0]), np.array([0]), np.array([1]),
                                   np.array([True, True]), rmax=0.2, hmin=3.0, max_iter=1)
    assert r == pytest.approx(0.2) and h.tolist() == pytest.approx([44.0, 66.0])
    with pytest.raises(ValueError, match="finite and positive"):
        rfactor_smooth_free(np.array([10.0, np.nan]), np.array([0]), np.array([1]),
                            np.array([False, True]), rmax=0.2, hmin=3.0)
    with pytest.raises(ValueError, match="bad controls"):
        rfactor_smooth_free(np.array([10.0, 20.0]), np.array([0]), np.array([1]),
                            np.array([False, True]), rmax=1.5, hmin=3.0)


def test_compose_refuses_bands_that_cannot_both_hold():
    """Round 2 F17: a later band silently overrode an earlier one."""
    x, y = np.meshgrid(np.arange(5) * 1000.0, np.arange(5) * 1000.0)
    b1 = np.where(y == 0, 1000.0, np.nan)
    b2 = np.where(y == 1000, 5000.0, np.nan)
    with pytest.raises(ValueError, match="cannot hold"):
        compose_sizing(np.full(x.shape, 3000.0), x, y, grade=0.2, bands=[b1, b2])


def test_check_no_overlap_finds_overlap_away_from_the_interface():
    """Round 2 F9: an outer element overlapping the base elsewhere passed."""
    from fvcom_mesh_tools.extend import check_no_overlap

    base = _base()
    on, oe = _outer()
    m = merge_outer(base, on, oe, [0, 3], [1, 2], [1, 2])
    assert check_no_overlap(m, base.n_elements)["n_outer_touching_base"] >= 1
    nodes = np.vstack([m.nodes, [[500.0, 500.0]]])              # a node inside the base
    bad = Fort14Mesh("m", nodes, np.r_[m.depths, 1.0],
                     np.vstack([m.elements, [[4, 5, 6]]]), m.open_boundaries, [])
    with pytest.raises(ValueError, match="overlap the base"):
        check_no_overlap(bad, base.n_elements)


def test_band_floor_count_is_taken_on_the_final_field():
    """Round 3 F1: a band smoothed below its floor counted 0; F3: iterators."""
    x, y = np.meshgrid(np.arange(5) * 100.0, np.arange(5) * 100.0)
    b0 = np.where(y == 0, 1000.0, np.nan)
    b1 = np.where(y == 400, np.array([1000, 1040, 1040, 1040, 1040.0])[None, :]
                  * np.ones((5, 1)), np.nan)
    floor = np.where(y == 400, b1, 0.0)
    h, rep = compose_sizing(np.full(x.shape, 1000.0), x, y, grade=0.2, floor=floor,
                            bands=iter([b0, b1]))
    assert rep["band_1_below_floor_cells"] >= 1
    assert "band_1_max_rel_deviation" in rep


def test_frozen_base_is_checked_bit_for_bit():
    """Review round 8 F12: -0.0 for +0.0 passed."""
    from dataclasses import replace

    import numpy as np
    import pytest

    from fvcom_mesh_tools.extend import verify_frozen_base
    from fvcom_mesh_tools.io.fort14 import Fort14Mesh

    nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    base = Fort14Mesh("b", nodes, np.array([0.0, 5.0, 5.0]), np.array([[0, 1, 2]]), [], [])
    with pytest.raises(ValueError, match="depths changed"):
        verify_frozen_base(replace(base, depths=np.array([-0.0, 5.0, 5.0])), base, [])
    with pytest.raises(ValueError, match="coordinates changed"):
        verify_frozen_base(replace(base, nodes=nodes.astype(np.float32)), base, [])


def test_trim_lone_corners_does_not_drop_mutually_supporting_elements():
    """Review round 9 F9: both triangles of a square were dropped, then the
    next round reduced an empty array."""
    import numpy as np

    from fvcom_mesh_tools.extend import trim_lone_corners

    out, mut, rep = trim_lone_corners(np.array([[0, 1, 2], [0, 2, 3]]), [True, True])
    assert len(out) == 1 and rep["n_elements_dropped"] == 1
    out, _, _ = trim_lone_corners(np.empty((0, 3), np.int64), np.empty(0, bool))
    assert out.shape == (0, 3)


def test_free_depth_bounds_hold_without_any_live_edge():
    """Review round 9 F10."""
    import numpy as np

    from fvcom_mesh_tools.extend import rfactor_smooth_free

    h, it, r = rfactor_smooth_free([5, 5, 1], np.array([0]), np.array([1]),
                                   [False, False, True], rmax=0.2, hmin=3, hmax=4)
    assert list(h) == [5, 5, 3] and it == 0


def test_merge_refuses_fractional_indices_and_a_float32_base():
    """Review round 9 F4 and F6."""
    import numpy as np
    import pytest

    from fvcom_mesh_tools.extend import merge_outer
    from fvcom_mesh_tools.io.fort14 import Fort14Mesh

    base = Fort14Mesh("b", np.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
                      np.full(4, 5.0), np.array([[0, 1, 2], [0, 2, 3]]), [np.array([1, 2])], [])
    outer = np.array([[1.0, 0.0], [1.0, 1.0], [2.0, 0.0], [2.0, 1.0]])
    with pytest.raises(ValueError, match="whole number"):
        merge_outer(base, outer, np.array([[0.0, 2.0, 1.9], [2.0, 3.0, 1.0]]),
                    [0, 1], [1, 2], [2, 3])
    b32 = Fort14Mesh("b", base.nodes.astype(np.float32), base.depths.astype(np.float32),
                     base.elements, base.open_boundaries, [])
    with pytest.raises(ValueError, match="float64"):
        merge_outer(b32, outer, np.array([[0, 2, 1], [2, 3, 1]]), [0, 1], [1, 2], [2, 3])


def test_frozen_base_elements_are_checked_by_dtype_and_bytes():
    """Review round 10 F10."""
    from dataclasses import replace

    import numpy as np
    import pytest

    from fvcom_mesh_tools.extend import verify_frozen_base
    from fvcom_mesh_tools.io.fort14 import Fort14Mesh

    base = Fort14Mesh("b", np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]), np.full(3, 5.0),
                      np.array([[0, 1, 2]], np.int64), [], [])
    with pytest.raises(ValueError, match="elements changed"):
        verify_frozen_base(replace(base, elements=base.elements.astype(np.int32)), base, [])


def test_trim_report_describes_the_returned_mesh():
    """Review round 10 F11."""
    import numpy as np

    from fvcom_mesh_tools.extend import trim_lone_corners

    out, _, rep = trim_lone_corners(np.array([[0, 1, 2], [0, 2, 3]]), [True, True], max_rounds=1)
    assert out.tolist() == [[0, 2, 3]]
    assert rep["lone_nodes_left"] == [0, 2, 3] and rep["round_limit_reached"]


def test_rounding_keeps_depths_inside_their_bounds():
    """Review round 10 F12."""
    import numpy as np

    from fvcom_mesh_tools.extend import round_depths_inside

    assert round_depths_inside([3.0000004], 3.0000004)[0] >= 3.0000004
    assert round_depths_inside([3.0000006], 3.0, 3.0000006)[0] <= 3.0000006
    assert list(round_depths_inside(np.array([3.1234567, 7.0]), 3.0)) == [3.123457, 7.0]
