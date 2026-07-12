"""Tests for the geometry-stage narrow-channel policy."""

import numpy as np
import pytest
from shapely.geometry import Point, box
from shapely.ops import unary_union

from fvcom_mesh_tools.prep.channel_policy_geom import (
    apply_channel_policy_to_land,
)

H = 350.0  # mesh-size floor (m)
DOMAIN = box(0, 0, 20000, 10000)
OBC = (1000, 5000)  # in the left/main water


def _water(land):
    return DOMAIN.difference(land)


class TestThroughChannel:
    def test_through_slot_gets_widened_not_closed(self):
        # island pair with a 300 m slot; water wraps around, so the
        # slot connects main to main -> widen by pushing banks
        land = unary_union([
            box(8000, 2000, 12000, 4850),
            box(8000, 5150, 12000, 8000),
        ])
        new_land, info = apply_channel_policy_to_land(
            land, DOMAIN, h_mesh_m=H, obc_point=OBC)
        assert len(info["widened"]) == 1
        assert len(info["closed"]) == 0
        w = _water(new_land)
        # banks pushed: >=600 m of water across the old 300 m slot
        assert w.contains(Point(10000, 4700))
        assert w.contains(Point(10000, 5300))

    def test_big_basin_channel_widened(self):
        # basin of ~12 cells behind a 300 m slot -> keep + widen
        land = unary_union([
            box(8000, 7000, 12000, 10000),
        ]).difference(unary_union([
            box(9000, 8000, 9800, 8800),      # 0.64 km2 ~ 12 cells
            box(9300, 7000, 9600, 8000),      # 300 m slot
        ]))
        new_land, info = apply_channel_policy_to_land(
            land, DOMAIN, h_mesh_m=H, obc_point=OBC)
        assert len(info["widened"]) == 1
        assert _water(new_land).contains(Point(9450, 8400))


class TestCloseCases:
    def test_dead_end_inlet_closed(self):
        land = box(15000, 8000, 20000, 10000).difference(
            box(17000, 8000, 17300, 9500))    # blind 300 m inlet
        new_land, info = apply_channel_policy_to_land(
            land, DOMAIN, h_mesh_m=H, obc_point=OBC)
        assert len(info["closed"]) == 1
        assert not _water(new_land).contains(Point(17150, 9000))

    def test_small_basin_closed_with_channel(self):
        land = box(2000, 7000, 6000, 10000).difference(
            unary_union([
                box(3000, 8000, 3500, 8500),  # 0.25 km2 ~ 4.7 cells
                box(3150, 7000, 3450, 8000),  # 300 m slot
            ]))
        new_land, info = apply_channel_policy_to_land(
            land, DOMAIN, h_mesh_m=H, obc_point=OBC)
        assert len(info["closed"]) == 1
        w = _water(new_land)
        assert not w.contains(Point(3250, 8250))   # basin filled
        assert not w.contains(Point(3300, 7500))   # slot filled


class TestRiverMouth:
    def test_river_touching_one_body_is_closed(self):
        # a long 300 m river entering from land: touches the main
        # body exactly once -> close, no matter how big main is
        land = box(8000, 0, 20000, 10000).difference(
            box(8000, 4850, 20000, 5150))    # river to the edge
        new_land, info = apply_channel_policy_to_land(
            land, DOMAIN, h_mesh_m=H, obc_point=OBC)
        assert len(info["widened"]) == 0
        assert len(info["closed"]) >= 1
        assert not _water(new_land).contains(Point(16000, 5000))

    def test_phantom_water_is_dropped(self):
        # land coverage gap far from the sea: an isolated "water"
        # pocket must not participate at all
        land = box(8000, 0, 20000, 10000).difference(
            box(15000, 2000, 18000, 4000))   # inland no-data hole
        new_land, info = apply_channel_policy_to_land(
            land, DOMAIN, h_mesh_m=H, obc_point=OBC)
        assert info["n_phantom_water_dropped"] >= 1
        # the hole is untouched land-side (not analysed, not edited)
        assert info["n_narrow"] == 0 or all(
            r["center"][0] < 15000 for r in
            info["widened"] + info["closed"])


class TestNetworkCases:
    def test_chained_corridor_stays_open(self):
        # Haneda case: slot-pocket-slot chain through an island
        # group; each link alone touches wide water once, but the
        # NETWORK connects main to main -> widen the whole chain
        land = unary_union([
            box(8000, 2000, 12000, 8000),
        ]).difference(unary_union([
            box(8000, 4850, 9500, 5150),
            box(9500, 4750, 10000, 5250),   # pocket ~4.7 cells (<6)
            box(10000, 4850, 12000, 5150),
        ]))
        new_land, info = apply_channel_policy_to_land(
            land, DOMAIN, h_mesh_m=H, obc_point=OBC)
        # ONE network spanning slot-pocket-slot, widened as a whole
        assert len(info["widened"]) == 1
        assert info["widened"][0]["n_members"] >= 2
        assert len(info["closed"]) == 0
        w = _water(new_land)
        assert w.contains(Point(8700, 5000))
        assert w.contains(Point(9750, 5000))   # pocket kept
        assert w.contains(Point(11000, 5000))

    def test_breakwater_gap_not_widened(self):
        # detached breakwater off a straight beach: the gap behind
        # it is 'through' topologically but saves no distance and
        # must NOT erode the beach
        land = unary_union([
            box(0, 8000, 20000, 10000),          # straight coast
            box(9000, 7700, 9800, 7900),         # breakwater
        ])
        new_land, info = apply_channel_policy_to_land(
            land, DOMAIN, h_mesh_m=H, obc_point=OBC)
        assert len(info["widened"]) == 0
        # gap filled, beach untouched
        assert not _water(new_land).contains(Point(9400, 7950))
        assert new_land.contains(Point(9400, 8100))


class TestValidation:
    def test_anisotropic_scale_raises(self):
        with pytest.raises(ValueError, match="anisotropic"):
            apply_channel_policy_to_land(
                box(0, 0, 1, 1), DOMAIN, h_mesh_m=H,
                obc_point=OBC, metric_scale=(1.0, 2.0))

    def test_overlarge_h_raises(self):
        with pytest.raises(RuntimeError, match="ALL water"):
            apply_channel_policy_to_land(
                box(0, 0, 1, 1), box(0, 0, 4000, 4000),
                h_mesh_m=1e5, obc_point=(2000, 2000))
