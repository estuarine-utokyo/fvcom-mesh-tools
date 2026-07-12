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

    def test_unwidenable_channel_is_closed_not_left_narrow(self):
        # a protected basin runs parallel 200 m south of the slot:
        # widening is capped by the barrier at ~1.3 rows. Owner
        # rule: a channel that cannot reach ~two rows is NOT
        # meshed at all -- the slot is filled, the basin and the
        # barrier survive untouched.
        land = unary_union([
            box(8000, 2000, 12000, 4850),
            box(8000, 5150, 12000, 8000),
        ]).difference(box(8500, 4450, 11500, 4650))
        recs, new_land, info = _run(land)
        assert len(info["blocked"]) == 1
        assert info["blocked"][0].get("closed") is True
        assert "two rows are not attainable" in \
            info["blocked"][0]["reason"]
        w = _water(new_land)
        # slot filled, not left one cell wide
        assert not w.contains(Point(10000, 5000))
        # basin and the protecting barrier survive
        assert w.contains(Point(10000, 4550))
        assert new_land.intersection(
            box(8600, 4650, 11400, 4850)).area > 0

    def test_big_pocket_chain_links_are_kept(self):
        # Keihin-canal pattern: main - slot - BIG pocket - slot -
        # BIG pocket. The second slot touches no main water, only
        # two big anchors -- it must be KEPT (its closure severed
        # the real Keihin canal, run 6185186).
        land = box(8000, 2000, 20000, 8000).difference(unary_union([
            box(8000, 4850, 10000, 5150),      # slot from main
            box(10000, 4300, 11200, 5700),     # big pocket A
            box(11200, 4850, 13000, 5150),     # link slot
            box(13000, 4300, 14200, 5700),     # big pocket B
        ]))
        recs, new_land, info = _run(land)
        w = _water(new_land)
        # the whole chain stays open end to end
        for x in (9000, 10600, 12000, 13600):
            assert w.contains(Point(x, 5000)), x
        # link slots widened, not closed
        closed_centers = [r["geometry"].representative_point()
                          for r in recs if r["action"] == "close"]
        assert not any(11200 < c.x < 13000 and 4700 < c.y < 5300
                       for c in closed_centers)

    def test_too_thin_through_channel_not_resolved(self):
        # a 120 m ditch (0.34 h) connects main to main: the sample
        # never resolves such channels -- keeping it forced a
        # 600 m corridor across land with one-wide remnants
        # (owner 2026-07-12). It must be closed, not widened.
        land = unary_union([
            box(8000, 2000, 12000, 4940),
            box(8000, 5060, 12000, 8000),
        ])
        recs, new_land, info = _run(land)
        assert info["kept"] == 0
        # ditch filled, banks untouched at width scale
        assert not _water(new_land).contains(Point(10000, 5000))
        assert new_land.contains(Point(10000, 4700))
        assert new_land.contains(Point(10000, 5300))

    def test_bridge_chained_canal_kept_and_opened(self):
        # Daishi-canal pattern: a through canal chopped by a
        # 120 m road strip drawn as land. The two halves must
        # chain into ONE through network and the bridge strip
        # must be carved open.
        land = unary_union([
            box(8000, 2000, 12000, 4850),
            box(8000, 5150, 12000, 8000),
        ]).union(box(9900, 4850, 10020, 5150))   # bridge strip
        recs, new_land, info = _run(land)
        assert info["kept"] == 1 and info["closed"] == 0
        assert info["bridges_opened"] >= 1
        w = _water(new_land)
        # continuous end to end THROUGH the ex-bridge
        for x in (8700, 9960, 11000):
            assert w.contains(Point(x, 5000)), x

    def test_parallel_canals_do_not_chain_across_land(self):
        # two PARALLEL dead-end slots 200 m apart share a long
        # frontage: they must NOT chain into one fake through
        # network (frontage gate), and no bridge may be opened
        # between them
        land = box(8000, 2000, 20000, 8000).difference(unary_union([
            box(8000, 4850, 11000, 5150),      # slot A (dead end)
            box(8000, 5350, 11000, 5650),      # slot B parallel
        ]))
        recs, new_land, info = _run(land)
        assert info["bridges_opened"] == 0
        # A and B stay separate records, never one network
        boxA, boxB = box(8500, 4850, 10500, 5150), \
            box(8500, 5350, 10500, 5650)
        for r in recs:
            g = r["geometry"]
            assert not (g.intersects(boxA) and g.intersects(boxB))

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
