"""normalize_unresolved_water (owner rule 2026-07-15, OW05
Urayasu): water we decided not to resolve is LAND for later
geometry decisions. Fills must cover (a) fully disconnected
components and (b) basins reachable only through sub-floor
passages -- and must NEVER touch the main body, water behind a
resolvable passage, or water joined through a carved corridor
(keep_tubes)."""

import numpy as np
import pytest
from shapely.geometry import Point, box
from shapely.ops import unary_union

from fvcom_mesh_tools.waterways import (
    _junction_bridge_pairs,
    normalize_unresolved_water,
)

H = 350.0
DOMAIN = box(0, 0, 20000, 10000)
OBC = (1000, 5000)
SCALE = (1.0, 1.0)


def _basin_land(slot_w: float):
    """Sealed 7.4 x 1.7 km basin along the top edge, sole
    entrance = a ``slot_w``-wide slot in the south wall."""
    return unary_union([
        box(6000, 8000, 10000 - slot_w / 2, 8300),   # south wall W
        box(10000 + slot_w / 2, 8000, 14000, 8300),  # south wall E
        box(6000, 8300, 6300, 10000),                # west wall
        box(13700, 8300, 14000, 10000),              # east wall
    ])


def _run(land, **kw):
    return normalize_unresolved_water(
        land, DOMAIN, h_mesh_m=H, obc_point=OBC,
        metric_scale=SCALE, **kw)


def _covers(fills, x, y):
    p = Point(x, y)
    return any(f.covers(p) or f.distance(p) < 1.0 for f in fills)


class TestBasins:
    def test_subfloor_entrance_basin_is_filled(self):
        fills, info = _run(_basin_land(slot_w=200.0))
        assert info["basin_parts_filled"] >= 1
        assert _covers(fills, 10000, 9000)          # basin core
        assert not _covers(fills, 10000, 3000)      # open bay
        assert info["area_filled_ha"] > 100.0

    def test_wide_entrance_basin_is_kept(self):
        # 1200 m entrance (> h): resolvable, nothing fills
        fills, info = _run(_basin_land(slot_w=1200.0))
        assert not _covers(fills, 10000, 9000)
        assert info["basin_parts_filled"] == 0

    def test_through_passage_untouched(self):
        # a 200 m top strip AND a 2 km bottom passage around one
        # obstacle: the narrow strip joins two main-side regions,
        # so nothing may fill
        land = box(8000, 2000, 12000, 9800)
        fills, info = _run(land)
        assert info["basin_parts_filled"] == 0
        assert info["components_filled"] == 0


class TestComponents:
    def test_disconnected_pond_is_filled(self):
        # land strips wall off the NE corner completely
        land = unary_union([
            box(15000, 6000, 15400, 10000),
            box(15000, 6000, 20000, 6400),
        ])
        fills, info = _run(land)
        assert info["components_filled"] == 1
        assert _covers(fills, 18000, 8000)
        assert not _covers(fills, 5000, 5000)


class TestKeepTubes:
    def test_tube_preserves_marginal_corridor_basin(self):
        # sub-floor entrance, but a carved-corridor tube runs
        # through it (marginal-kept branch): basin must survive
        arc = np.array([[10000.0, 5000.0], [10000.0, 9000.0]])
        w = np.array([200.0, 200.0])
        fills, info = _run(_basin_land(slot_w=200.0),
                           keep_tubes=[(arc, w)])
        assert not _covers(fills, 10000, 9000)
        assert info["basin_parts_filled"] == 0

    def test_tube_over_wall_does_not_bridge(self):
        # tube crossing SOLID land must not fabricate a link: the
        # basin has NO water entrance at all (separate component,
        # filled regardless of the tube)
        land = unary_union([
            box(6000, 8000, 14000, 8300),
            box(6000, 8300, 6300, 10000),
            box(13700, 8300, 14000, 10000),
        ])
        arc = np.array([[10000.0, 5000.0], [10000.0, 9000.0]])
        w = np.array([400.0, 400.0])
        fills, info = _run(land, keep_tubes=[(arc, w)])
        assert info["components_filled"] == 1
        assert _covers(fills, 10000, 9000)


class TestJunctionBridges:
    """F9-c4 lesson (run 6210307): per-branch arcs stop at
    junctions; the junction hole must be bridged at corridor
    width or the through path is left to realization roulette."""

    def test_junction_gap_bridged_at_narrower_width(self):
        a1 = np.array([[0.0, 0.0], [1000.0, 0.0]])
        a2 = np.array([[1400.0, 0.0], [2400.0, 0.0]])
        done = [(a1, np.array([700.0, 700.0])),
                (a2, np.array([600.0, 600.0]))]
        prs = _junction_bridge_pairs(done, 350.0, (1.0, 1.0))
        assert len(prs) == 1     # gap 400 in (0.35h, 2h]
        pa, pb, wj = prs[0]
        assert wj == 600.0
        assert {tuple(pa), tuple(pb)} == {(1000.0, 0.0),
                                          (1400.0, 0.0)}

    def test_touching_and_far_endpoints_skipped(self):
        a1 = np.array([[0.0, 0.0], [1000.0, 0.0]])
        a2 = np.array([[1000.0, 0.0], [2000.0, 0.0]])  # gap 0
        a3 = np.array([[3000.0, 0.0], [4000.0, 0.0]])  # gap 1000
        done = [(a, np.array([700.0, 700.0]))
                for a in (a1, a2, a3)]
        assert _junction_bridge_pairs(done, 350.0,
                                      (1.0, 1.0)) == []

    def test_network_geom_guard_rejects_cross_land_pair(self):
        # parallel branches across a land spit: endpoints 400 m
        # apart, but the connecting segment leaves the network
        # water -> no bridge (run 6218497: fabricated H5-d3
        # passage)
        a1 = np.array([[0.0, 0.0], [1000.0, 0.0]])
        a2 = np.array([[1400.0, 0.0], [2400.0, 0.0]])
        done = [(a1, np.array([700.0, 700.0])),
                (a2, np.array([600.0, 600.0]))]
        water = unary_union([
            box(-100, -200, 1100, 200),
            box(1300, -200, 2500, 200)])   # gap NOT in water
        assert _junction_bridge_pairs(
            done, 350.0, (1.0, 1.0), network_geom=water) == []
        joined = box(-100, -200, 2500, 200)  # junction in water
        assert len(_junction_bridge_pairs(
            done, 350.0, (1.0, 1.0), network_geom=joined)) == 1


class TestLoudFailures:
    def test_no_water_raises(self):
        with pytest.raises(RuntimeError):
            _run(DOMAIN.buffer(1.0))

    def test_anisotropic_scale_raises(self):
        with pytest.raises(ValueError):
            normalize_unresolved_water(
                box(1, 1, 2, 2), DOMAIN, h_mesh_m=H,
                obc_point=OBC, metric_scale=(1.0, 2.0))
