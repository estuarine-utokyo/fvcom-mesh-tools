"""The step after a hires refinement: floor, cap, smooth -- and touch nothing else."""

from __future__ import annotations

import numpy as np
import pytest

from fvcom_mesh_tools.cli.refine_depths import finish_depths


def strip(n=6):
    """A strip of triangles, so every node has neighbours on an edge."""
    xs = np.arange(n, dtype=float) * 100.0
    nodes = np.column_stack([np.repeat(xs, 2), np.tile([0.0, 100.0], n)])
    tri = []
    for i in range(n - 1):
        a = 2 * i
        tri += [[a, a + 1, a + 3], [a, a + 3, a + 2]]
    return nodes, np.asarray(tri, dtype=np.int64)


def test_the_floor_and_the_cap_apply_only_to_the_patch():
    _, tri = strip()
    h = np.array([1.0, 2.0, 0.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 400.0, 11.0, 12.0])
    movable = np.zeros(12, dtype=bool)
    movable[2:10] = True
    out, rep = finish_depths(h, tri, movable, hmin=3.0, hmax=300.0, rfactor=0.2)
    assert out[0] == 1.0 and out[1] == 2.0, (
        "a frozen depth below the floor is the BASE's and must not move")
    assert rep["n_frozen_below_hmin"] == 2
    assert out[2] >= 3.0 and out[9] <= 300.0
    assert rep["max_frozen_change_m"] == 0.0


def test_a_frozen_node_above_the_datum_stops_the_step_with_a_reason():
    """The limiter refuses depths it is not defined on, including frozen ones.

    goto2023 is floored at 3 m so this cannot happen there, but a base that
    carries its own tidal flats would need flooring itself -- a change to the
    model, not to the patch -- and the message has to say that rather than
    leave the limiter to complain about a mesh the caller did not choose.
    """
    _, tri = strip()
    h = np.array([1.0, -2.0] + [5.0] * 10)
    movable = np.zeros(12, dtype=bool)
    movable[2:] = True
    with pytest.raises(ValueError, match="FROZEN node"):
        finish_depths(h, tri, movable, hmin=3.0, hmax=300.0, rfactor=0.2)


def test_the_smoothing_is_what_the_floor_makes_possible():
    """`|dh|/(hi+hj)` is a sign error dressed as a gradient above the datum."""
    _, tri = strip()
    h = np.array([4.0, 4.0, -1.0, -0.5, 8.0, 8.0, 20.0, 20.0,
                  20.0, 20.0, 20.0, 20.0])
    movable = np.zeros(12, dtype=bool)
    movable[2:10] = True
    out, rep = finish_depths(h, tri, movable, hmin=3.0, hmax=300.0, rfactor=0.2)
    assert (out[movable] > 0).all()
    assert rep["rfactor_report"]["converged"], (
        "with a floor the limiter has a feasible set on every movable edge")
    e = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]],
                                     tri[:, [2, 0]]]), axis=1), axis=0)
    r = np.abs(out[e[:, 0]] - out[e[:, 1]]) / (out[e[:, 0]] + out[e[:, 1]])
    touched = movable[e].any(axis=1)
    assert r[touched].max() <= 0.2 + 1e-9


def test_a_run_of_all_frozen_depths_is_reported_not_fixed():
    _, tri = strip()
    h = np.array([1.0, 10.0] + [5.0] * 10)      # a frozen pair at r = 0.82
    movable = np.zeros(12, dtype=bool)
    movable[4:] = True
    out, rep = finish_depths(h, tri, movable, hmin=3.0, hmax=300.0, rfactor=0.2)
    assert rep["rfactor_report"]["n_over_rmax_frozen_pair"] >= 1
    assert not rep["rfactor_report"]["converged"]
    assert out[0] == 1.0 and out[1] == 10.0


@pytest.mark.parametrize("kw,msg", [
    ({"hmin": 0.0}, "hmin"),
    ({"hmin": -1.0}, "hmin"),
    ({"hmax": 1.0}, "hmax"),
])
def test_invalid_bounds_are_refused(kw, msg):
    _, tri = strip()
    h = np.full(12, 5.0)
    args = {"hmin": 3.0, "hmax": 300.0, "rfactor": 0.2, **kw}
    with pytest.raises(ValueError, match=msg):
        finish_depths(h, tri, np.ones(12, dtype=bool), **args)


def test_non_finite_depths_are_refused():
    _, tri = strip()
    h = np.full(12, 5.0)
    h[3] = np.nan
    with pytest.raises(ValueError, match="finite"):
        finish_depths(h, tri, np.ones(12, dtype=bool), hmin=3.0, hmax=300.0,
                      rfactor=0.2)


def test_convergence_is_judged_at_the_depth_file_s_own_precision():
    """The base is written to six decimals and sits ON its own bound.

    `TokyoBay_dep_m7001tp_rfac0p2_cap300.dat` has a max r of exactly 0.2000
    and 497 of its 8,858 edges sit on it; read back from six decimals they
    come out at 0.2 + 3e-8, and limit_rfactor's 1e-9 test reported the patch
    NOT CONVERGED for 374 edges it is forbidden to touch.
    """
    _, tri = strip()
    # two frozen nodes exactly on the bound, as the file would store them
    h = np.array([round(4.0, 6), round(4.0 * 1.5, 6)] + [6.0] * 10)
    movable = np.zeros(12, dtype=bool)
    movable[2:] = True
    out, rep = finish_depths(h, tri, movable, hmin=3.0, hmax=300.0, rfactor=0.2)
    assert rep["max_r_frozen_pair"] == pytest.approx(0.2, abs=1e-9)
    assert rep["converged_at_write_precision"], (
        "an edge sitting on the bound is not an edge over it")
    assert rep["n_over_movable_at_tolerance"] == 0
    assert rep["write_precision_tolerance"] > 0
