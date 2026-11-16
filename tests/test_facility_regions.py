"""fmesh-facility-regions: facility near fields -> refinement regions -> recipe."""

import json
import math

import pytest
from shapely.geometry import Polygon, box

from fvcom_mesh_tools import facility_regions as fr
from fvcom_mesh_tools.cli import facility_regions as cli

# Near Futtsu: 1e-3 deg is about 91 m in lon and 111 m in lat.
LON0, LAT0 = 139.80, 35.33


def rect_lonlat(x0, y0, w_m, h_m):
    """A w x h metre rectangle at (x0, y0) metres from (LON0, LAT0)."""
    dlon = 1.0 / (111320.0 * math.cos(math.radians(LAT0)))
    dlat = 1.0 / 110574.0
    return box(LON0 + x0 * dlon, LAT0 + y0 * dlat, LON0 + (x0 + w_m) * dlon,
               LAT0 + (y0 + h_m) * dlat)


def test_utm_zone():
    assert fr.utm_crs_for(139.8, 35.3) == "EPSG:32654"
    assert fr.utm_crs_for(-70.0, -33.0) == "EPSG:32719"
    with pytest.raises(ValueError):
        fr.utm_crs_for(0.0, 85.0)


def test_length_scale_rules():
    p = box(0, 0, 100, 20)
    assert fr.length_scale(p, "long_side") == pytest.approx(100.0)
    assert fr.length_scale(p, "sqrt_area") == pytest.approx(math.sqrt(2000.0))
    assert fr.length_scale(p, "diameter") == pytest.approx(math.hypot(100, 20))
    # rotation does not change the long side
    rot = Polygon([(0, 0), (70.71, 70.71), (56.57, 84.85), (-14.14, 14.14)])
    assert fr.length_scale(rot) == pytest.approx(100.0, rel=1e-3)
    with pytest.raises(ValueError):
        fr.length_scale(p, "perimeter")
    with pytest.raises(ValueError):
        fr.length_scale(Polygon(), "long_side")


def test_close_facilities_merge_far_ones_do_not():
    # A and B 60 m apart, each 50 m long: near fields of 50 m touch -> one region.
    # C is 1 km away -> its own region.
    facs = [("A", rect_lonlat(0, 0, 50, 20)), ("B", rect_lonlat(110, 0, 50, 20)),
            ("C", rect_lonlat(1100, 0, 30, 30))]
    regions, ff, crs = fr.facility_regions(facs, prefix="t")
    assert crs == "EPSG:32654"
    assert [r.name for r in regions] == ["t_c01", "t_c02"]
    assert sorted(regions[0].members) == ["A", "B"] and regions[1].members == ["C"]
    assert [round(f.buffer_m) for f in ff] == [50, 50, 30]
    # the region contains each member's full near field
    for f in ff:
        reg = regions[0] if f.fid in "AB" else regions[1]
        assert reg.polygon_m.buffer(1e-6).contains(f.polygon_m.buffer(f.buffer_m * 0.99))


def test_factor_and_min_buffer():
    facs = [("A", rect_lonlat(0, 0, 50, 20)), ("B", rect_lonlat(200, 0, 50, 20))]
    regions, _, _ = fr.facility_regions(facs)            # 50 m fields, 150 m gap
    assert len(regions) == 2
    regions, ff, _ = fr.facility_regions(facs, factor=2.0)   # 100 m fields touch
    assert len(regions) == 1
    regions, ff, _ = fr.facility_regions(facs, factor=0.5, min_buffer_m=80.0)
    assert [round(f.buffer_m) for f in ff] == [80, 80] and len(regions) == 1
    with pytest.raises(ValueError):
        fr.facility_regions(facs, factor=0.0)
    with pytest.raises(ValueError):
        fr.facility_regions(facs, factor=float("nan"))
    with pytest.raises(ValueError):
        fr.facility_regions([("A", facs[0][1]), ("A", facs[1][1])])


def test_enclosed_water_is_filled_unless_kept():
    # four bars around a 200 m square: their near fields enclose water
    facs = [("S", rect_lonlat(0, 0, 400, 20)), ("N", rect_lonlat(0, 380, 400, 20)),
            ("W", rect_lonlat(0, 0, 20, 400)), ("E", rect_lonlat(380, 0, 20, 400))]
    filled, _, _ = fr.facility_regions(facs, factor=0.2)
    kept, _, _ = fr.facility_regions(facs, factor=0.2, fill_holes=False)
    assert len(filled) == len(kept) == 1
    assert len(filled[0].polygon_m.interiors) == 0
    assert len(kept[0].polygon_m.interiors) == 1
    assert filled[0].area_m2 > kept[0].area_m2


def test_summary_estimates():
    facs = [("A", rect_lonlat(0, 0, 50, 20))]
    regions, ff, _ = fr.facility_regions(facs)
    rep = fr.summary(regions, ff, target_h_m=5.0, ambient_h_m=50.0, gradation=0.165,
                     within_m=regions[0].polygon_m.buffer(10.0))
    assert rep["transition_width_m"] == pytest.approx(45.0 / 0.165)
    assert rep["core_outside_within_frac"] == pytest.approx(0.0, abs=1e-9)
    assert 0.0 < rep["hole_outside_within_frac"] < 1.0
    expect = regions[0].area_m2 / (math.sqrt(3) / 4 * 25.0)
    assert rep["core_elements_estimate"] == round(expect)


def _geojson(path, polys, props=None):
    feats = [{"type": "Feature", "properties": {"id": k + 1, **(props or {})},
              "geometry": json.loads(json.dumps(p.__geo_interface__))}
             for k, p in enumerate(polys)]
    path.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))


def test_read_facilities_filters_and_multipolygon(tmp_path):
    from shapely.geometry import MultiPolygon

    a, b = rect_lonlat(0, 0, 50, 20), rect_lonlat(500, 0, 50, 20)
    f = tmp_path / "f.geojson"
    feats = [{"type": "Feature", "properties": {"id": 7, "ground": "x"},
              "geometry": MultiPolygon([a, b]).__geo_interface__},
             {"type": "Feature", "properties": {"id": 8, "ground": "y"},
              "geometry": a.__geo_interface__}]
    f.write_text(json.dumps({"type": "FeatureCollection", "features": json.loads(
        json.dumps(feats))}))
    got = fr.read_facilities(f, {"ground": "x"})
    assert [g[0] for g in got] == ["7.0", "7.1"]
    with pytest.raises(ValueError, match="no facility"):
        fr.read_facilities(f, {"ground": "z"})


def test_cli_recipe_loads_in_refine(tmp_path):
    """The recipe the CLI writes is accepted by refine.load_refine, one region
    per cluster, each resolving to exactly one polygon."""
    from fvcom_mesh_tools import refine
    from fvcom_mesh_tools.sizing import _geometry_from_file

    nets = tmp_path / "nets.geojson"
    _geojson(nets, [rect_lonlat(0, 0, 50, 20), rect_lonlat(110, 0, 50, 20),
                    rect_lonlat(1100, 0, 30, 30)])
    lic = tmp_path / "lic.geojson"
    _geojson(lic, [rect_lonlat(-500, -500, 2200, 1000)], {"ground": "g"})
    base = tmp_path / "case_grd.dat"
    dep = tmp_path / "case_dep.dat"
    base.write_text("Node Number = 0\nCell Number = 0\n")
    dep.write_text("Node Number = 0\n")
    out = tmp_path / "recipes"
    rc = cli.main(["--facilities", str(nets), "--target-h-m", "5", "--ambient-h-m", "50",
                   "--within", str(lic), "--within-where", "ground=g", "--prefix", "nets",
                   "--base-mesh", str(base), "--base-depth", str(dep),
                   "--out-geojson", str(out / "nets_regions.geojson"),
                   "--out-recipe", str(out / "nets_5m.yaml")])
    assert rc == 0
    cfg = refine.load_refine(out / "nets_5m.yaml")
    names = [r.name for r in cfg["refine"]]
    assert names == ["nets_c01", "nets_c02"]
    assert all(r.target_h_m == 5.0 for r in cfg["refine"])
    assert cfg["coastline"] == "preserve"
    for name in names:
        geom, _ = _geometry_from_file({"file": str(out / "nets_regions.geojson"),
                                       "where": {"region": name}})
        assert geom.is_valid and geom.area > 0
    rep = json.loads((out / "nets_5m.report.json").read_text())
    assert rep["facilities"] == 3 and rep["regions"] == 2
    assert rep["core_outside_within_frac"] == pytest.approx(0.0, abs=1e-9)
