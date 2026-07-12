"""Arc-based mathematical waterway detection (no reference mesh:
OSM-only input). Mirrors the owner's keep/close rules previously
tested on the polygon-morphology policy, now with barrier-safe
arc carving."""

import pytest
from shapely.geometry import Point, box
from shapely.ops import unary_union

from fvcom_mesh_tools.waterways import (
    apply_waterway_policy,
    detect_waterways,
)

H = 350.0
DOMAIN = box(0, 0, 20000, 10000)
OBC = (1000, 5000)
SCALE = (1.0, 1.0)


def _run(land, **kw):
    recs = detect_waterways(land, DOMAIN, h_mesh_m=H,
                            obc_point=OBC, metric_scale=SCALE,
                            **kw)
    new_land, info = apply_waterway_policy(
        land, DOMAIN, recs, h_mesh_m=H, metric_scale=SCALE)
    return recs, new_land, info


def _water(land):
    return DOMAIN.difference(land)


class TestKeep:
    def test_through_slot_widened_to_two_rows(self):
        land = unary_union([
            box(8000, 2000, 12000, 4850),
            box(8000, 5150, 12000, 8000),
        ])
        recs, new_land, info = _run(land)
        assert info["kept"] == 1 and info["closed"] == 0
        assert not info["blocked"]
        w = _water(new_land)
        # 300 m slot -> ~600 m: banks pushed on both sides
        assert w.contains(Point(10000, 4720))
        assert w.contains(Point(10000, 5280))
        # far banks untouched
        assert new_land.contains(Point(10000, 4300))
        assert new_land.contains(Point(10000, 5700))

    def test_big_basin_channel_kept(self):
        land = unary_union([
            box(8000, 7000, 12000, 10000),
        ]).difference(unary_union([
            box(9000, 8000, 9800, 8800),
            box(9300, 7000, 9600, 8000),
        ]))
        recs, new_land, info = _run(land)
        assert info["kept"] == 1
        assert _water(new_land).contains(Point(9450, 8400))
        assert _water(new_land).contains(Point(9450, 7500))


class TestClose:
    def test_dead_end_inlet_closed(self):
        land = box(15000, 8000, 20000, 10000).difference(
            box(17000, 8000, 17300, 9500))
        recs, new_land, info = _run(land)
        assert info["closed"] == 1 and info["kept"] == 0
        assert not _water(new_land).contains(Point(17150, 9000))

    def test_small_basin_closed_with_channel(self):
        land = box(2000, 7000, 6000, 10000).difference(
            unary_union([
                box(3000, 8000, 3500, 8500),
                box(3150, 7000, 3450, 8000),
            ]))
        recs, new_land, info = _run(land)
        assert info["closed"] == 1
        w = _water(new_land)
        assert not w.contains(Point(3250, 8250))
        assert not w.contains(Point(3300, 7500))

    def test_river_touching_one_body_is_closed(self):
        land = box(8000, 0, 20000, 10000).difference(
            box(8000, 4850, 20000, 5150))
        recs, new_land, info = _run(land)
        assert info["kept"] == 0 and info["closed"] >= 1
        assert not _water(new_land).contains(Point(16000, 5000))

    def test_breakwater_gap_closed_beach_untouched(self):
        land = unary_union([
            box(0, 8000, 20000, 10000),
            box(9000, 7700, 9800, 7900),
        ])
        recs, new_land, info = _run(land)
        assert info["kept"] == 0
        assert not _water(new_land).contains(Point(9400, 7950))
        assert new_land.contains(Point(9400, 8100))


class TestNetworksAndSafety:
    def test_chained_corridor_kept_whole(self):
        land = unary_union([
            box(8000, 2000, 12000, 8000),
        ]).difference(unary_union([
            box(8000, 4850, 9500, 5150),
            box(9500, 4750, 10000, 5250),
            box(10000, 4850, 12000, 5150),
        ]))
        recs, new_land, info = _run(land)
        assert info["kept"] == 1 and info["closed"] == 0
        w = _water(new_land)
        for x in (8700, 9750, 11000):
            assert w.contains(Point(x, 5000))
        # widened to ~2 rows mid-chain
        assert w.contains(Point(11000, 4750))

    def test_widening_respects_parallel_basin_barrier(self):
        # a protected basin runs parallel 200 m south of the slot:
        # widening must keep >= min_gap of land, not merge them
        land = unary_union([
            box(8000, 2000, 12000, 4850),
            box(8000, 5150, 12000, 8000),
        ]).difference(box(8500, 4450, 11500, 4650))
        recs, new_land, info = _run(land)
        w = _water(new_land)
        # slot widened northwards at least
        assert w.contains(Point(10000, 5280))
        # the land sliver between slot and basin survives
        assert new_land.intersection(
            box(8600, 4650, 11400, 4850)).area > 0

    def test_phantom_water_ignored(self):
        land = box(8000, 0, 20000, 10000).difference(
            box(15000, 2000, 18000, 4000))
        recs, new_land, info = _run(land)
        # the inland no-data hole must not be analysed or edited
        assert _water(new_land).contains(Point(16500, 3000)) or \
            new_land.contains(Point(16500, 3000))
        for r in recs:
            assert not r["geometry"].contains(Point(16500, 3000))

    def test_anisotropic_scale_raises(self):
        with pytest.raises(ValueError, match="anisotropic"):
            detect_waterways(box(0, 0, 1, 1), DOMAIN, h_mesh_m=H,
                             obc_point=OBC,
                             metric_scale=(1.0, 2.0))
