"""Tests for channel_arcs: arc ordering and the barrier-safe
corridor carve (metric coordinates; metric_scale=(1, 1))."""

import numpy as np
import pytest
import shapely
from shapely.geometry import box

from fvcom_mesh_tools.channel_arcs import (
    arc_from_points,
    carve_channel_corridor,
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
