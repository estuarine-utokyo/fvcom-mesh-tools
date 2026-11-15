"""New-region parameters of the fmesh-pipeline path (2026-11-16): recipe CRS,
region water files, polygon engine and size field, channel protection,
lone nodes / land boundaries, depth products, recipe rules."""

from __future__ import annotations

import json

import numpy as np
import pytest

from fvcom_mesh_tools.cli import pipeline as pl
from fvcom_mesh_tools.io import Fort14Mesh, read_fort14, write_fort14

EPSG = 32631          # UTM 31N (Volta mouth)
X0, Y0 = 300_000.0, 640_000.0


def _lattice(nx, ny, dx=100.0, x0=X0, y0=Y0):
    xs, ys = np.meshgrid(x0 + dx * np.arange(nx), y0 + dx * np.arange(ny))
    P = np.column_stack([xs.ravel(), ys.ravel()])
    T = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            T.append([a, a + 1, a + nx])
            T.append([a + 1, a + nx + 1, a + nx])
    return P, np.asarray(T, np.int64)


def _mesh(P, T, opens=()):
    return Fort14Mesh(title="t", nodes=np.asarray(P, float), depths=np.full(len(P), 5.0),
                      elements=np.asarray(T, np.int64),
                      open_boundaries=[np.asarray(o, np.int64) for o in opens],
                      land_boundaries=[])


# ------------------------------------------------------------------ recipe CRS
def test_recipe_epsg_sources():
    assert pl._recipe_epsg({"crs": "EPSG:32631"}) == 32631
    assert pl._recipe_epsg({"crs": 32631}) == 32631
    assert pl._recipe_epsg({"finish": {"utm_epsg": 32654}}) == 32654       # legacy key
    assert pl._recipe_epsg({"prep": {"bbox": [139.4, 34.9, 140.3, 35.9]}}) == 32654
    assert pl._recipe_epsg({"prep": {"bbox": [-0.15, 5.0, 1.45, 6.5]}}) == 32631
    assert pl._recipe_epsg({"prep": {"bbox": [151.0, -34.0, 151.5, -33.5]}}) == 32756


def test_recipe_epsg_errors():
    with pytest.raises(SystemExit):
        pl._recipe_epsg({})                                  # no silent Tokyo default
    with pytest.raises(SystemExit):
        pl._recipe_epsg({"crs": 32631, "finish": {"utm_epsg": 32654}})
    with pytest.raises(SystemExit):
        pl._recipe_epsg({"crs": "UTM31"})


def test_library_functions_require_crs(tmp_path):
    from fvcom_mesh_tools.autofinish import apply_directives
    from fvcom_mesh_tools.finishing import finish_constrained_mesh
    from fvcom_mesh_tools.obc_tools import assign_west_south_obc
    from fvcom_mesh_tools.qa import run_qa

    P, T = _lattice(4, 4)
    m = _mesh(P, T)
    with pytest.raises(ValueError, match="utm_epsg"):
        finish_constrained_mesh(m, tmp_path / "x.shp", tmp_path, log=lambda *a: None)
    with pytest.raises(ValueError, match="utm_epsg"):
        assign_west_south_obc(m, log=lambda *a: None)
    with pytest.raises(ValueError, match="utm_epsg"):
        apply_directives(m, [], log=lambda *a: None)
    shp = tmp_path / "land.gpkg"
    import geopandas as gpd
    from shapely.geometry import box

    gpd.GeoDataFrame(geometry=[box(0, 0, 1, 1)], crs=4326).to_file(shp)
    with pytest.raises(ValueError, match="utm_epsg"):
        run_qa(m, land_solid_shp=shp, channel_check=False)


# ------------------------------------------------------------ recipe structure
def test_section_and_recipe_rules(tmp_path):
    import yaml

    f = tmp_path / "r.yaml"
    f.write_text(yaml.safe_dump({"version": "v3", "obc": {"isobath_m": 125},
                                 "pipeline": {"name": "x", "out_dir": str(tmp_path)}}))
    assert pl._load_recipe(f, section="pipeline")["name"] == "x"
    with pytest.raises(SystemExit):
        pl._load_recipe(f, section="missing")
    with pytest.raises(SystemExit, match="obc"):
        pl._check_recipe({"domain": {"obc_line": [[0, 0], [1, 0]]}, "finish": {}})
    pl._check_recipe({"domain": {"obc_line": [[0, 0], [1, 0]]}, "finish": {}, "obcfinal": {}})
    pl._check_recipe({"finish": {}})                         # no line: nothing to protect


def test_path_expands_env(monkeypatch):
    monkeypatch.setenv("FMESH_TEST_ROOT", "/data/x")
    assert str(pl._path("${FMESH_TEST_ROOT}/a.tif")) == "/data/x/a.tif"
    monkeypatch.delenv("FMESH_TEST_ROOT")
    with pytest.raises(SystemExit):
        pl._path("${FMESH_TEST_ROOT}/a.tif")


# ----------------------------------------------------------- region water files
def test_merge_water_files(tmp_path):
    import geopandas as gpd
    from shapely.geometry import box

    from fvcom_mesh_tools.prep.shoreline import merge_water_files

    a, b = tmp_path / "gh.gpkg", tmp_path / "tg.gpkg"
    gpd.GeoDataFrame({"fclass": ["water"], "osm_id": ["1"]}, geometry=[box(0.5, 5.8, 0.6, 5.9)],
                     crs=4326).to_file(a)
    gpd.GeoDataFrame({"fclass": ["riverbank"], "osm_id": ["2"], "extra": [1]},
                     geometry=[box(1.2, 6.0, 1.3, 6.1)], crs=4326).to_file(b)
    out = merge_water_files([a, b], (0.0, 5.0, 2.0, 7.0), tmp_path / "m.gpkg")
    g = gpd.read_file(out)
    assert sorted(g.fclass) == ["riverbank", "water"] and "extra" not in g.columns
    with pytest.raises(FileNotFoundError):
        merge_water_files([tmp_path / "nope.shp"], (0, 5, 2, 7), tmp_path / "m2.gpkg")


# ------------------------------------------------------------------ size field
def test_size_field_background_regions_gradation():
    from pyproj import Transformer
    from shapely.geometry import box

    from fvcom_mesh_tools.mesh_engine.polygon import build_size_field

    dom = box(X0, Y0, X0 + 10_000, Y0 + 5_000)
    inv = Transformer.from_crs(f"EPSG:{EPSG}", "EPSG:4326", always_xy=True)
    lo0, la0 = inv.transform(X0 + 4000, Y0 + 2000)
    lo1, la1 = inv.transform(X0 + 5000, Y0 + 3000)
    reg = [{"name": "channel", "geometry": {"bbox": [lo0, la0, lo1, la1]}, "target_h_m": 100.0}]
    sz = build_size_field(dom, crs_epsg=EPSG, background_m=1000.0, raster_m=100.0,
                          regions=reg, gradation=0.2, log=lambda *a: None)
    h, dx = sz["h"], sz["dx"]
    assert h.min() == pytest.approx(100.0)
    assert h.max() == pytest.approx(1000.0)
    assert np.abs(np.diff(h, axis=1)).max() <= 0.2 * dx + 1e-9   # graded
    assert sz["report"]["regions"][0]["cells"] > 0
    with pytest.raises(ValueError):
        build_size_field(dom, crs_epsg=EPSG, log=lambda *a: None)  # neither raster nor background


def test_size_field_raster(tmp_path):
    import rasterio
    from rasterio.transform import from_origin
    from shapely.geometry import box

    from fvcom_mesh_tools.mesh_engine.polygon import build_size_field

    tif = tmp_path / "size.tif"
    a = np.full((60, 120), 500.0, np.float32)
    a[:, :10] = np.nan
    with rasterio.open(tif, "w", driver="GTiff", height=60, width=120, count=1, dtype="float32",
                       crs=f"EPSG:{EPSG}", transform=from_origin(X0 - 1000, Y0 + 5000, 100, 100),
                       nodata=np.nan) as dst:
        dst.write(a, 1)
    dom = box(X0, Y0, X0 + 8000, Y0 + 4000)
    sz = build_size_field(dom, crs_epsg=EPSG, raster=tif, max_m=2000.0, log=lambda *a: None)
    assert np.all(np.diff(sz["ys"]) > 0) and sz["h"].max() == pytest.approx(2000.0)
    with pytest.raises(ValueError, match="CRS"):
        build_size_field(dom, crs_epsg=32654, raster=tif, max_m=2000.0, log=lambda *a: None)
    with pytest.raises(ValueError, match="cover"):
        build_size_field(box(X0, Y0, X0 + 50_000, Y0 + 4000), crs_epsg=EPSG, raster=tif,
                         max_m=2000.0, log=lambda *a: None)


# ------------------------------------------------------------- polygon engine
def _domain():
    from shapely.geometry import LineString, Polygon

    # sea box 6 x 3 km whose south edge is the open boundary, plus a 300 m wide
    # channel 3 km long going north
    sea = Polygon([(X0, Y0), (X0 + 6000, Y0), (X0 + 6000, Y0 + 3000), (X0 + 3150, Y0 + 3000),
                   (X0 + 3150, Y0 + 6000), (X0 + 2850, Y0 + 6000), (X0 + 2850, Y0 + 3000),
                   (X0, Y0 + 3000)])
    obc = LineString([(X0, Y0), (X0 + 6000, Y0)])
    chan = Polygon([(X0 + 2850, Y0 + 3000), (X0 + 3150, Y0 + 3000), (X0 + 3150, Y0 + 6000),
                    (X0 + 2850, Y0 + 6000)])
    return sea, obc, chan


def test_polygon_engine_open_string_and_components():
    from fvcom_mesh_tools.mesh_engine.polygon import (
        build_size_field,
        generate_polygon_mesh,
    )

    water, obc, chan = _domain()
    sz = build_size_field(water, crs_epsg=EPSG, background_m=400.0, raster_m=50.0,
                          log=lambda *a: None)
    res = generate_polygon_mesh(water, sz, obc_line=obc, ladder_k=1.0, max_iter=30,
                                keep_polygons=chan, log=lambda *a: None)
    p, t, o = res["points"], res["cells"], res["open_string"]
    assert len(o) >= 10
    assert np.abs(p[o, 1] - Y0).max() < 1e-6                  # the open string is on the line
    cnt = np.bincount(t.ravel(), minlength=len(p))
    lone = cnt == 1
    lone[o] = False
    assert not lone.any()                                      # no single-element nodes
    post = [r for r in res["components"] if r["stage"] == "post-clean"]
    assert any(r["kept"] for r in post)
    assert all(r["kept"] or (r["area_km2"] < 0.5 and not r["in_keep"]) for r in post)


# --------------------------------------------------------- channel protection
def _strip_mesh(n=10):
    """A one-element-wide strip (a channel two nodes across): every interior
    element is a bridge with all three nodes on the boundary, the cap class
    the legacy finish deletes, which splits the channel."""
    P = np.array([(100.0 * i, 0.0) for i in range(n)] + [(100.0 * i, 100.0) for i in range(n)])
    T = []
    for i in range(n - 1):
        T.append([i, i + 1, n + i])
        T.append([i + 1, n + i + 1, n + i])
    return P, np.asarray(T, np.int64)


def test_cap_deletion_guard_and_component_policy():
    import shapely
    from shapely.geometry import box

    from fvcom_mesh_tools.finishing import (
        _cap_deletion_allowed,
        keep_components_protected,
        remove_elements_keep_boundaries,
    )

    P, T = _strip_mesh()
    m = _mesh(P, T, opens=[[0, 10]])                 # the west end is the open boundary
    mid = len(T) // 2
    pm = shapely.contains_xy(box(-10, -10, 2000, 200), P[:, 0], P[:, 1])
    ok = _cap_deletion_allowed(m, pm)
    assert not ok[mid]                               # protected bridge: kept
    assert ok[-1]                                    # the east-end ear may go
    assert not ok[0]                                 # two open-boundary nodes: kept
    assert _cap_deletion_allowed(m, None)[mid]       # unprotected: legacy behaviour
    m2, tag = remove_elements_keep_boundaries(m, np.arange(len(T)) != mid, pm)
    assert [list(o) for o in m2.open_boundaries] == [[0, 10]]
    with pytest.raises(RuntimeError, match="protected"):
        keep_components_protected(m2, tag, policy="fail", log=lambda *a: None)
    m3, info, _ = keep_components_protected(m2, tag, policy="report", log=lambda *a: None)
    assert len(info["dropped"]) == 1 and info["dropped"][0]["protected"]
    assert m3.n_elements == mid
    assert m3.nodes[m3.open_boundaries[0]].tolist() == [[0.0, 0.0], [0.0, 100.0]]


def test_remove_lone_nodes_and_land_rebuild():
    from fvcom_mesh_tools.finishing import rebuild_land_boundaries, remove_lone_nodes

    P, T = _lattice(4, 4)
    # node 3 (corner) is used by one element; node 0 too, but it is an OBC node
    m = _mesh(P, T, opens=[[0, 1, 2]])
    m2, info = remove_lone_nodes(m, log=lambda *a: None)
    cnt = np.bincount(m2.elements.ravel(), minlength=m2.n_nodes)
    o = set(m2.open_boundaries[0].tolist())
    assert all(cnt[v] > 1 for v in range(m2.n_nodes) if v not in o)
    assert info["elements_deleted"] >= 1
    m3 = rebuild_land_boundaries(m2)
    ob = m3.open_boundaries[0]
    (ib, land), = m3.land_boundaries
    assert ib == 20 and land[0] == ob[-1] and land[-1] == ob[0]


def test_finish_keeps_open_boundary_and_protected_channel(tmp_path):
    import geopandas as gpd
    from pyproj import Transformer
    from shapely.geometry import LineString

    from fvcom_mesh_tools.finishing import finish_constrained_mesh
    from fvcom_mesh_tools.mesh_engine.polygon import (
        build_size_field,
        generate_polygon_mesh,
    )

    water, obc, chan = _domain()
    sz = build_size_field(water, crs_epsg=EPSG, background_m=150.0, raster_m=50.0,
                          log=lambda *a: None)
    res = generate_polygon_mesh(water, sz, obc_line=obc, max_iter=30, log=lambda *a: None)
    m = _mesh(res["points"], res["cells"], opens=[res["open_string"]])
    import shapely

    cen0 = m.nodes[m.elements].mean(axis=1)
    n_chan0 = int(shapely.contains_xy(chan, cen0[:, 0], cen0[:, 1]).sum())
    assert n_chan0 > 20
    rings = [LineString(water.exterior.coords)]
    shp = tmp_path / "shore.shp"
    gpd.GeoDataFrame(geometry=rings, crs=EPSG).to_crs(4326).to_file(shp)
    out, info = finish_constrained_mesh(m, shp, tmp_path / "w", utm_epsg=EPSG,
                                        optimize_budget_s=20.0, protect_polygons=chan,
                                        max_valence=8, remove_lone=True, freeze_obc_band=True,
                                        log=lambda *a: None)
    assert set(info["obc_band"]) == {"restored", "refused"}
    assert np.isfinite(out.depths).all() and np.allclose(out.depths, 5.0)  # tags undone
    assert len(out.open_boundaries) == 1 and len(out.open_boundaries[0]) >= 3
    assert np.bincount(out.elements.ravel()).max() <= 8
    cen = out.nodes[out.elements].mean(axis=1)
    n_chan = int(shapely.contains_xy(chan, cen[:, 0], cen[:, 1]).sum())
    assert n_chan >= 0.8 * n_chan0                                     # channel kept
    assert out.land_boundaries and out.land_boundaries[0][0] == 20
    del Transformer


# --------------------------------------------------------------- depth stage
def test_depths_stage_products(tmp_path):
    P, T = _lattice(3, 3)
    m = _mesh(P, T, opens=[[0, 1, 2]])
    f14 = tmp_path / "final.14"
    write_fort14(m, f14)
    dep = tmp_path / "x_dep.dat"
    dep.write_text(f"Node Number = {len(P)}\n"
                   + "".join(f"{x:.4f} {y:.4f} {7.5 + i}\n" for i, (x, y) in enumerate(P)))
    rec = {"name": "x", "depths": {"product": "external_dep_file", "path": str(dep)}}
    out = pl._stage_depths(rec, tmp_path, {"final_mesh": str(f14)}, lambda *a: None)
    got = read_fort14(out["final_mesh"]).depths
    assert np.allclose(got, 7.5 + np.arange(len(P)))
    bad = tmp_path / "bad_dep.dat"
    bad.write_text(f"Node Number = {len(P)}\n"
                   + "".join(f"{x + 5:.4f} {y:.4f} 1\n" for x, y in P))
    with pytest.raises(SystemExit, match="another mesh"):
        pl._stage_depths({"name": "x", "depths": {"product": "external_dep_file",
                                                  "path": str(bad)}},
                         tmp_path, {"final_mesh": str(f14)}, lambda *a: None)
    assert pl._stage_depths({"name": "x", "depths": {"product": "none"}}, tmp_path,
                            {"final_mesh": str(f14)}, lambda *a: None) == {"depths_product": "none"}
    with pytest.raises(SystemExit):
        pl._stage_depths({"name": "x", "depths": {"product": "gebco"}}, tmp_path,
                         {"final_mesh": str(f14)}, lambda *a: None)


# --------------------------------------------- pipeline: domain + polygon build
def test_pipeline_domain_build_qa(tmp_path):
    import geopandas as gpd
    import yaml

    water, obc, chan = _domain()
    gp = tmp_path / "dom.gpkg"
    gpd.GeoDataFrame(geometry=[water], crs=EPSG).to_file(gp, layer="water")
    gpd.GeoDataFrame(geometry=[obc], crs=EPSG).to_file(gp, layer="obc")
    gpd.GeoDataFrame(geometry=[chan], crs=EPSG).to_file(gp, layer="channel")
    rec = {"name": "t", "out_dir": str(tmp_path / "out"), "crs": f"EPSG:{EPSG}",
           "domain": {"water": {"path": str(gp), "layer": "water"},
                      "obc_line": {"path": str(gp), "layer": "obc"}},
           "protect_polygons": [{"path": str(gp), "layer": "channel"}],
           "build": {"engine": "oceanmesh-polygon", "max_iter": 30,
                     "size": {"background_m": 400.0, "raster_m": 50.0}},
           "depths": {"product": "none"},
           "qa": {"lang": "en"}}
    f = tmp_path / "r.yaml"
    f.write_text(yaml.safe_dump({"pipeline": rec}))
    assert pl.main([str(f), "--section", "pipeline"]) == 0
    prov = json.loads((tmp_path / "out" / "pipeline_provenance.json").read_text())
    assert prov["artifacts"]["raw_mesh_crs"] == f"EPSG:{EPSG}"
    m = read_fort14(prov["artifacts"]["raw_mesh"])
    assert len(m.open_boundaries) == 1 and m.nodes[:, 0].min() > 1e5      # metres
    line = json.loads((tmp_path / "out" / "domain" / "obc_line_lonlat.json").read_text())
    assert len(line) == 2 and 0.0 < line[0][0] < 3.0                       # lon/lat


def test_polygon_engine_lost_component_and_expected_drop():
    from shapely.geometry import box

    from fvcom_mesh_tools.mesh_engine.polygon import (
        build_size_field,
        generate_polygon_mesh,
    )

    water, obc, _chan = _domain()
    pond = box(X0 + 8000, Y0, X0 + 9500, Y0 + 1500)          # 2.25 km2, not connected
    w2 = water.union(pond)
    sz = build_size_field(w2, crs_epsg=EPSG, background_m=400.0, raster_m=50.0,
                          log=lambda *a: None)
    with pytest.raises(RuntimeError, match="dropped"):
        generate_polygon_mesh(w2, sz, obc_line=obc, max_iter=20, log=lambda *a: None)
    res = generate_polygon_mesh(w2, sz, obc_line=obc, max_iter=20,
                                expected_drop=pond.buffer(100.0), log=lambda *a: None)
    assert any(r["expected_drop"] and not r["kept"] for r in res["components"])
    assert res["points"][:, 0].max() < X0 + 7000                  # the pond is not meshed
