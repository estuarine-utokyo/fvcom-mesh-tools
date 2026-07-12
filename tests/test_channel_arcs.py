"""Tests for channel_arcs: arc ordering and the barrier-safe
corridor carve (metric coordinates; metric_scale=(1, 1))."""

import numpy as np
import pytest
import shapely
from shapely.geometry import box

from fvcom_mesh_tools.channel_arcs import (
    arc_from_points,
    bank_chains,
    carve_channel_corridor,
    snap_arc_to_channel,
)

SCALE = (1.0, 1.0)
DOM = box(0.0, 0.0, 10_000.0, 10_000.0)


def _land_band_with_slot(slot_x0=4_900.0, slot_x1=5_100.0):
    """Land band y in [4000, 6000] with a narrow through-slot."""
    band = box(0.0, 4_000.0, 10_000.0, 6_000.0)
    slot = box(slot_x0, 3_900.0, slot_x1, 6_100.0)
    return band.difference(slot)


ARC = np.array([[5_000.0, 3_000.0], [5_000.0, 5_000.0],
                [5_000.0, 7_000.0]])


def test_arc_from_points_orders_scattered_curve():
    t = np.linspace(0.0, np.pi, 40)
    pts = np.column_stack([t * 3_000.0,
                           1_000.0 * np.sin(t)])
    rng = np.random.default_rng(0)
    arc = arc_from_points(pts[rng.permutation(40)],
                          smooth_passes=0)
    ends = {tuple(np.round(arc[0])), tuple(np.round(arc[-1]))}
    assert tuple(np.round(pts[0])) in ends
    assert tuple(np.round(pts[-1])) in ends
    # ordered path is not much longer than the true curve
    d = np.linalg.norm(np.diff(arc, axis=0), axis=1).sum()
    d_true = np.linalg.norm(np.diff(pts, axis=0), axis=1).sum()
    assert d < 1.35 * d_true


def test_carve_widens_through_slot():
    land = _land_band_with_slot()
    new_land, info = carve_channel_corridor(
        land, ARC, 700.0, min_gap_m=150.0, metric_scale=SCALE,
        domain_poly=DOM)
    # slot widened: a point 250 m off-axis is now water
    assert not new_land.covers(shapely.Point(4_750.0, 5_000.0))
    assert not new_land.covers(shapely.Point(5_250.0, 5_000.0))
    # far land untouched
    assert new_land.covers(shapely.Point(2_000.0, 5_000.0))
    assert info["arc_on_land_m"] < 1.0


def test_carve_keeps_gap_to_parallel_basin():
    # isolated basin inside the band, 200 m east of the corridor's
    # would-be edge: the carve must stop min_gap short of it
    land = _land_band_with_slot().difference(
        box(5_300.0, 4_500.0, 5_800.0, 5_500.0))
    new_land, _ = carve_channel_corridor(
        land, ARC, 700.0, min_gap_m=150.0, metric_scale=SCALE,
        domain_poly=DOM)
    # protective land strip between corridor and basin survives
    assert new_land.covers(shapely.Point(5_225.0, 5_000.0))
    # western side still fully widened
    assert not new_land.covers(shapely.Point(4_750.0, 5_000.0))


def test_carve_refuses_arc_across_solid_barrier():
    land = box(0.0, 4_000.0, 10_000.0, 6_000.0)  # no slot
    with pytest.raises(RuntimeError, match="crosses land"):
        carve_channel_corridor(
            land, ARC, 700.0, min_gap_m=150.0, metric_scale=SCALE,
            domain_poly=DOM)


def test_detected_arc_corner_clip_is_preserved_not_pierced():
    # arc clips the band/bank corner within tolerance: with
    # carve_crossings=False the corner must stay LAND (no
    # fabricated passage), while the slot still widens and the
    # ends stay connected through the real channel
    land = _land_band_with_slot()
    arc = np.array([[5_150.0, 3_800.0], [5_000.0, 5_000.0],
                    [5_000.0, 6_200.0]])
    new_land, info = carve_channel_corridor(
        land, arc, 700.0, min_gap_m=150.0, metric_scale=SCALE,
        domain_poly=DOM, arc_on_land_tol_m=400.0,
        carve_crossings=False)
    assert info["arc_on_land_m"] > 20.0
    # the clipped corner region is preserved as land
    assert new_land.covers(shapely.Point(5_115.0, 4_120.0))
    # the slot is still widened elsewhere
    assert not new_land.covers(shapely.Point(4_750.0, 5_500.0))


def test_manual_edit_still_carves_crossing():
    # explicit opt-in (pier drawn as land): the same corner clip
    # IS carved with carve_crossings=True
    land = _land_band_with_slot()
    arc = np.array([[5_150.0, 3_800.0], [5_000.0, 5_000.0],
                    [5_000.0, 6_200.0]])
    new_land, _ = carve_channel_corridor(
        land, arc, 700.0, min_gap_m=150.0, metric_scale=SCALE,
        domain_poly=DOM, arc_on_land_tol_m=400.0,
        carve_crossings=True)
    assert not new_land.covers(shapely.Point(5_115.0, 4_120.0))


def test_thin_barrier_inside_corridor_width_is_protected():
    # a ditch runs parallel 50 m beyond the slot bank, WITHIN the
    # corridor half-width: the wall between slot and ditch must
    # survive (fabricated-passage regression, breach at
    # 139.8991/35.3703)
    land = _land_band_with_slot().difference(
        box(5_150.0, 4_200.0, 5_250.0, 5_800.0))
    new_land, _ = carve_channel_corridor(
        land, ARC, 700.0, min_gap_m=150.0, metric_scale=SCALE,
        domain_poly=DOM, carve_crossings=False)
    # the 50 m wall stays land
    assert new_land.covers(shapely.Point(5_125.0, 5_000.0))
    # west side still widened
    assert not new_land.covers(shapely.Point(4_750.0, 5_000.0))


def test_carve_refuses_anisotropic_scale():
    with pytest.raises(ValueError, match="anisotropic"):
        carve_channel_corridor(
            _land_band_with_slot(), ARC, 700.0, min_gap_m=150.0,
            metric_scale=(1.0, 2.0), domain_poly=DOM)


def test_arc_needs_two_points():
    with pytest.raises(ValueError):
        arc_from_points(np.array([[0.0, 0.0]]))


def _land_kinked_slot():
    """Slot whose axis is OFFSET from the guide arc (center x=5200,
    width 300) with a full-width pier/bridge blockage midway (the
    Haneda D-runway pattern: real passage drawn as land)."""
    band = box(0.0, 4_000.0, 10_000.0, 6_000.0)
    slot = box(5_050.0, 3_900.0, 5_350.0, 6_100.0)
    bridge = box(5_050.0, 4_950.0, 5_350.0, 5_050.0)
    return band.difference(slot).union(bridge)


def test_snap_recovers_offset_center_width_and_pinch():
    land = _land_kinked_slot()
    guide = np.array([[5_000.0, 3_200.0], [5_000.0, 5_000.0],
                      [5_000.0, 6_800.0]])
    snap = snap_arc_to_channel(land, guide, metric_scale=SCALE,
                               step_m=150.0, smooth_passes=0)
    mid = np.abs(snap["arc"][:, 1] - 5_500.0) < 400.0
    assert np.allclose(snap["arc"][mid, 0], 5_200.0, atol=30.0)
    assert np.allclose(snap["width_m"][mid], 300.0, atol=30.0)
    # the 100 m-long full-width blockage is marked un-snapped and
    # its width interpolated back to ~300 from the clean neighbours
    pinch = np.abs(snap["arc"][:, 1] - 5_000.0) < 45.0
    assert pinch.any()
    assert (~snap["snapped"][pinch]).all()
    assert (snap["width_m"][pinch] > 250.0).all()


def test_bank_chains_sit_inside_water():
    land = _land_kinked_slot()
    guide = np.array([[5_000.0, 3_200.0], [5_000.0, 5_000.0],
                      [5_000.0, 6_800.0]])
    snap = snap_arc_to_channel(land, guide, metric_scale=SCALE,
                               step_m=150.0, smooth_passes=0)
    pfix, egfix = bank_chains(snap, spacing_m=350.0,
                              metric_scale=SCALE)
    assert len(pfix) >= 4 and len(egfix) >= 2
    # chains: every edge joins two points on the SAME bank
    for a, b in egfix:
        assert abs(pfix[a][0] - pfix[b][0]) < 40.0
    # clean-station bank nodes sit strictly inside the slot water
    for x, y in pfix:
        if abs(y - 5_000.0) > 120.0:      # outside the pinch
            assert 5_050.0 < x < 5_350.0


def test_carve_variable_width_removes_only_pinch():
    land = _land_kinked_slot()
    guide = np.array([[5_000.0, 3_200.0], [5_000.0, 5_000.0],
                      [5_000.0, 6_800.0]])
    snap = snap_arc_to_channel(land, guide, metric_scale=SCALE,
                               step_m=150.0)
    new_land, info = carve_channel_corridor(
        land, snap["arc"], snap["width_carve_m"], min_gap_m=150.0,
        metric_scale=SCALE, domain_poly=DOM)
    # bridge pinch removed
    assert not new_land.covers(shapely.Point(5_120.0, 5_000.0))
    # real banks preserved outside the pinch (faithful coastline)
    assert new_land.covers(shapely.Point(5_030.0, 5_600.0))
    assert new_land.covers(shapely.Point(5_370.0, 5_600.0))
    # land actually removed is small (just the bridge notch)
    assert info["land_removed_m2"] < 50_000.0


def test_split_choke_edges_doubles_the_section():
    from fvcom_mesh_tools.algorithms.obc_finish import (
        split_choke_edges,
    )
    from fvcom_mesh_tools.io import Fort14Mesh

    # 4x2-node channel strip: interior edge A2-B1 has both ends on
    # opposite banks with a long boundary detour -> choke
    A = [(0.0, 0.0), (600.0, 0.0), (1200.0, 0.0), (1800.0, 0.0)]
    B = [(0.0, 500.0), (600.0, 500.0), (1200.0, 500.0),
         (1800.0, 500.0)]
    nodes = np.asarray(A + B)
    els = np.asarray([
        [0, 1, 4], [1, 5, 4], [1, 2, 5],
        [2, 6, 5], [2, 3, 6], [3, 7, 6]])
    mesh = Fort14Mesh(
        title="t", nodes=nodes, depths=np.full(8, 5.0),
        elements=els, open_boundaries=[], land_boundaries=[])
    out, info = split_choke_edges(mesh)
    assert info["split"] == 1
    assert info["edges"] == [(2, 5)]
    assert out.n_nodes == 9
    assert out.n_elements == 8
    # midpoint inserted on the choke edge
    assert np.allclose(out.nodes[8], [900.0, 250.0])
    # all elements keep positive (CCW) area
    p = out.nodes
    t = out.elements
    ar = 0.5 * ((p[t[:, 1], 0] - p[t[:, 0], 0])
                * (p[t[:, 2], 1] - p[t[:, 0], 1])
                - (p[t[:, 2], 0] - p[t[:, 0], 0])
                * (p[t[:, 1], 1] - p[t[:, 0], 1]))
    assert (ar > 0).all()
