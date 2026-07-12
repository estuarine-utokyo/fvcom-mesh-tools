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
