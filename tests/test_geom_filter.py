"""Tests for filter_multipolygon_by_area."""

from __future__ import annotations

import pytest

from fvcom_mesh_tools.io import filter_multipolygon_by_area

shapely = pytest.importorskip("shapely")  # noqa: F841


def _multipoly():
    """Return a multipolygon with one ~1 km^2 polygon (with a tiny hole and a
    1-km-side hole) and one ~25 m x 25 m blip - in lon/lat, so the area
    helper has to project to evaluate."""
    from shapely.geometry import MultiPolygon, Polygon
    # A degree at lat=35 is ~111 km lat / 91 km lon. Use a 0.01 deg square
    # which is roughly 910 m x 1110 m ~ 1 km^2.
    big_outer = [(139.0, 35.0), (139.01, 35.0), (139.01, 35.01), (139.0, 35.01)]
    # Tiny hole: 0.0001 deg x 0.0001 deg ~ 9.1 m x 11.1 m ~ 100 m^2
    tiny_hole = [
        (139.001, 35.001), (139.0011, 35.001),
        (139.0011, 35.0011), (139.001, 35.0011),
    ]
    # Big hole: 0.005 deg x 0.005 deg ~ 455 m x 555 m ~ 250000 m^2
    big_hole = [
        (139.002, 35.002), (139.007, 35.002),
        (139.007, 35.007), (139.002, 35.007),
    ]
    big = Polygon(big_outer, holes=[tiny_hole, big_hole])
    blip = Polygon([(140.0, 35.0), (140.00025, 35.0),
                    (140.00025, 35.00025), (140.0, 35.00025)])
    return MultiPolygon([big, blip])


def test_no_filter_returns_input_unchanged() -> None:
    mp = _multipoly()
    out = filter_multipolygon_by_area(mp)
    assert len(out.geoms) == len(mp.geoms)
    assert sum(len(list(p.interiors)) for p in out.geoms) == sum(
        len(list(p.interiors)) for p in mp.geoms
    )


def test_polygon_threshold_drops_blip() -> None:
    mp = _multipoly()
    # Threshold 1000 m^2 -> drops the ~750 m^2 blip; keeps the 1 km^2 big.
    out = filter_multipolygon_by_area(mp, min_polygon_area_m2=1000.0)
    assert len(out.geoms) == 1


def test_island_threshold_keeps_big_hole_drops_tiny() -> None:
    mp = _multipoly()
    out = filter_multipolygon_by_area(mp, min_island_area_m2=1000.0)
    # Same number of outer polygons.
    assert len(out.geoms) == len(mp.geoms)
    # Big polygon should keep only the big hole; tiny hole filtered.
    big = out.geoms[0]
    assert len(list(big.interiors)) == 1


def test_combined_filter_drops_blip_and_tiny_hole() -> None:
    mp = _multipoly()
    out = filter_multipolygon_by_area(
        mp, min_polygon_area_m2=1000.0, min_island_area_m2=1000.0,
    )
    assert len(out.geoms) == 1
    assert len(list(out.geoms[0].interiors)) == 1


def test_negative_thresholds_rejected() -> None:
    mp = _multipoly()
    with pytest.raises(ValueError):
        filter_multipolygon_by_area(mp, min_polygon_area_m2=-1.0)
    with pytest.raises(ValueError):
        filter_multipolygon_by_area(mp, min_island_area_m2=-1.0)


def test_polygon_input_promoted() -> None:
    """A bare Polygon should be promoted to a length-1 MultiPolygon."""
    from shapely.geometry import Polygon
    poly = Polygon([(139.0, 35.0), (139.01, 35.0),
                    (139.01, 35.01), (139.0, 35.01)])
    out = filter_multipolygon_by_area(poly, min_polygon_area_m2=100.0)
    assert len(out.geoms) == 1
