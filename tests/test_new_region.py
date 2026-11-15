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
    ok = _cap_deletion_allowed(m, pm[T].any(axis=1))   # r2 #1: an element mask
    assert not ok[mid]                               # protected bridge: kept
    assert not ok[-1]                                # review r1 #2: no ear in protection either
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
    assert info["elements_deleted"] + info["flipped"] >= 1   # r3 #2: a flip comes first
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
                            {"final_mesh": str(f14)}, lambda *a: None) == {
        "depths_product": "none", "final_mesh": str(f14)}
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
           "obc": {"trim": 1, "snap": False},
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


def test_obc_perpendicularity_and_obc_lock_key(tmp_path):
    from fvcom_mesh_tools.finishing import obc_perpendicularity

    P, T = _lattice(5, 4)
    m = _mesh(P, T, opens=[[0, 1, 2, 3, 4]])
    dev = obc_perpendicularity(m)
    assert dev.shape == (5,) and np.allclose(dev, 0.0)          # vertical edges exist
    m.nodes = m.nodes.copy()
    m.nodes[5:10, 0] += 60.0                                     # shear the first ring
    assert obc_perpendicularity(m).max() > 20.0
    f14 = tmp_path / "f.14"
    write_fort14(_mesh(P, T, opens=[[0, 1, 2, 3, 4]]), f14)
    with pytest.raises(SystemExit, match="obc_lock"):
        pl._stage_finishing({"name": "x", "finishing": {"auto": {"obc_lock": "free"}}},
                            tmp_path, {"final_mesh": str(f14)}, lambda *a: None)
    out = pl._stage_finishing({"name": "x", "finishing": {"auto": {"obc_lock": "none"}}},
                              tmp_path, {"final_mesh": str(f14)}, lambda *a: None)
    assert np.allclose(read_fort14(out["final_mesh"]).nodes, P)  # clean lattice: no patch


def test_site_fallback_not_needed_on_clean_mesh():
    P, T = _lattice(5, 4)
    m = _mesh(P, T, opens=[[0, 1, 2, 3, 4]])
    out, rec = pl._finishing_site_fallback({"crs": EPSG}, m, {}, lambda *a: None)
    assert out is m and rec["outcome"] == "not needed"
    score, bad = pl._gate_violations(m)
    assert not any(k in bad for k in pl._SITE_FALLBACK_GATES)


# ------------------------------------------- review newregion r1 (regressions)
def _thin_strip(n=19, width=30.0):
    """One-element-wide strip 100 m per cell, ``width`` m across: every node is
    ON the outline and every element a sliver cap (min angle < 20 deg)."""
    P = np.array([(100.0 * i, 0.0) for i in range(n)] + [(100.0 * i, width) for i in range(n)])
    T = []
    for i in range(n - 1):
        T.append([i, i + 1, n + i])
        T.append([i + 1, n + i + 1, n + i])
    return P, np.asarray(T, np.int64)


def test_r1_1_protection_includes_outline(tmp_path, monkeypatch):
    """#1: a channel meshed bank to bank is protected although no node is
    strictly inside the polygon; finish deletes none of its elements."""
    from shapely.geometry import box

    import fvcom_mesh_tools.finishing as fin

    P, T = _thin_strip()
    water = box(0, 0, 1800, 30)
    assert fin.protected_node_mask(P, water).all()
    monkeypatch.setattr(fin, "phase_h_optimize", lambda m, **k: (m, {}))
    monkeypatch.setattr(fin, "build_coastline_projector", lambda *a, **k: None)
    import geopandas as gpd

    shp = tmp_path / "s.shp"
    from shapely.geometry import LineString

    gpd.GeoDataFrame(geometry=[LineString(water.exterior.coords)], crs=EPSG
                     ).to_crs(4326).to_file(shp)
    m = _mesh(P, T, opens=[[0, 19]])
    out, info = fin.finish_constrained_mesh(m, shp, tmp_path / "w", utm_epsg=EPSG,
                                            protect_polygons=water, log=lambda *a: None)
    assert out.n_elements == len(T) and info["caps_removed"] == 0


def test_r1_2_lone_nodes_respect_protection_and_obc():
    """#2: a fully protected strip loses nothing to lone-node removal (it was
    eaten from the far end); the leftover lone node fails or is reported."""
    from shapely.geometry import box

    from fvcom_mesh_tools.finishing import remove_lone_nodes

    P, T = _strip_mesh(10)                           # 900 x 100 m, OBC at the west end
    m = _mesh(P, T, opens=[[0, 10]])
    prot = box(0, 0, 900, 100)
    out, info = remove_lone_nodes(m, protect_polygons=prot, policy="report",
                                  log=lambda *a: None)
    # r3 #2: the far-end lone node is now fixed by a split, never by deleting water
    assert info["elements_deleted"] == 0 and info["lone_kept"] == 0
    assert out.n_elements >= m.n_elements
    one = _mesh(np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]]), np.array([[0, 1, 2]]))
    with pytest.raises(RuntimeError, match="single-element"):
        remove_lone_nodes(one, protect_polygons=box(0, 0, 100, 100), policy="fail",
                          log=lambda *a: None)
    # unprotected: elements holding open-boundary nodes are never deleted
    P3, T3 = _strip_mesh(4)                          # 300 x 100 m, OBC on the west column
    m3 = _mesh(P3, T3, opens=[[0, 4]])
    out3, _ = remove_lone_nodes(m3, policy="report", log=lambda *a: None)
    assert out3.n_elements >= 1 and len(out3.open_boundaries[0]) == 2


def test_r1_3_polygon_build_requires_obc_restoration():
    """#3: a polygon build with an open-boundary line may trim the open string;
    a recipe without obc/obcfinal is refused."""
    rec = {"domain": {"obc_line": [[0.1, 5.7], [1.2, 6.0]]},
           "build": {"engine": "oceanmesh-polygon"}, "qa": {}}
    with pytest.raises(SystemExit, match="trim"):
        pl._check_recipe(rec)
    pl._check_recipe(dict(rec, obc={}))


def test_r1_3_polygon_engine_reports_open_cover():
    from fvcom_mesh_tools.mesh_engine.polygon import build_size_field, generate_polygon_mesh

    water, obc, _chan = _domain()
    sz = build_size_field(water, crs_epsg=EPSG, background_m=400.0, raster_m=50.0,
                          log=lambda *a: None)
    res = generate_polygon_mesh(water, sz, obc_line=obc, max_iter=20, log=lambda *a: None)
    cov = res["open_cover"]
    assert cov["line_m"] == pytest.approx(6000.0)
    assert 0.0 <= cov["uncovered_m"] <= 6000.0 and "protected_repair_area_m2" in res


def test_r1_4_band_restore_interpolates_new_node_depths():
    """#4: a node created by a split keeps the interpolated depth (51 m between
    1 m and 101 m), not the nearest original node's depth."""
    from fvcom_mesh_tools.finishing import _restore_obc_band, _tag_obc_band

    P = np.array([[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]])
    T = np.array([[0, 1, 2], [0, 2, 3]])
    m = _mesh(P, T, opens=[[0, 1]])
    m.depths = np.array([1.0, 101.0, 101.0, 1.0])
    state = _tag_obc_band(m)
    tags = m.depths.copy()
    # an optimizer split of edge 0-1: new node 4 at the midpoint, tag = mean of the tags
    m2 = Fort14Mesh(title="t", nodes=np.vstack([P, [[50.0, 0.0]]]),
                    depths=np.append(tags, 0.5 * (tags[0] + tags[1])),
                    elements=np.array([[0, 4, 2], [4, 1, 2], [0, 2, 3]]),
                    open_boundaries=[np.array([0, 4, 1])], land_boundaries=[])
    out, _ = _restore_obc_band(m2, state, log=lambda *a: None)
    assert out.depths[4] == pytest.approx(51.0)
    assert np.allclose(out.depths[:4], [1.0, 101.0, 101.0, 1.0])


def test_r1_5_expected_drop_is_an_area_fraction():
    """#5: a 2.25 km2 pond finely meshed in a small corner is NOT an expected
    drop when expected_drop covers only that corner (4 % of its area)."""
    from pyproj import Transformer
    from shapely.geometry import box

    from fvcom_mesh_tools.mesh_engine.polygon import build_size_field, generate_polygon_mesh

    water, obc, _chan = _domain()
    pond = box(X0 + 8000, Y0, X0 + 9500, Y0 + 1500)
    corner = box(X0 + 8000, Y0, X0 + 8300, Y0 + 300)
    w2 = water.union(pond)
    inv = Transformer.from_crs(f"EPSG:{EPSG}", "EPSG:4326", always_xy=True)
    lo0, la0 = inv.transform(*corner.bounds[:2])
    lo1, la1 = inv.transform(*corner.bounds[2:])
    reg = [{"name": "corner", "geometry": {"bbox": [lo0, la0, lo1, la1]}, "target_h_m": 40.0}]
    sz = build_size_field(w2, crs_epsg=EPSG, background_m=700.0, raster_m=20.0, regions=reg,
                          log=lambda *a: None)
    with pytest.raises(RuntimeError, match="dropped"):
        generate_polygon_mesh(w2, sz, obc_line=obc, max_iter=15, expected_drop=corner,
                              log=lambda *a: None)


def test_r1_6_raster_edge_half_cell_keeps_edge_values():
    """#6: inside the raster footprint but outside the outermost centres the
    size is the edge value, not the global maximum."""
    from fvcom_mesh_tools.mesh_engine.polygon import size_interpolator

    xs = ys = np.arange(50.0, 1000.0, 100.0)
    h = np.full((len(ys), len(xs)), 400.0)
    h[:, 0] = 100.0
    hat = size_interpolator({"xs": xs, "ys": ys, "h": h})
    assert hat([[50.0, 500.0]])[0] == pytest.approx(100.0)
    assert hat([[0.0, 500.0]])[0] == pytest.approx(100.0)
    assert hat([[1000.0, 1000.0]])[0] == pytest.approx(400.0)


def test_r1_7_site_fallback_refuses_a_new_failing_gate(monkeypatch):
    """#7: fewer C1 offenders never pay for a newly failing other gate."""
    import fvcom_mesh_tools.algorithms.boundary_snap as bs
    import fvcom_mesh_tools.site_session as ss

    P, T = _lattice(5, 4)
    m = _mesh(P, T, opens=[[0, 1, 2, 3, 4]])
    seq = iter([(1002, {"c1_min_angle": 1002}),
                (1001, {"c1_min_angle": 1000, "obc_perpendicularity": 1})])
    monkeypatch.setattr(pl, "_gate_violations", lambda mesh: next(seq))
    monkeypatch.setattr(ss, "apply_site_operators", lambda mesh, lines, **k: (mesh, []))
    monkeypatch.setattr(bs, "load_polylines", lambda *a, **k: [])
    monkeypatch.setattr(pl, "_shoreline", lambda *a: "x")
    out, rec = pl._finishing_site_fallback({"crs": EPSG}, m, {}, lambda *a: None)
    assert rec["outcome"] == "refused" and out is m


def test_r1_8_depths_use_the_latest_mesh(tmp_path):
    """#8: the depth stage takes the mesh the last mesh stage wrote (a polished
    mesh after obc), and product none propagates it."""
    P, T = _lattice(3, 3)
    old = _mesh(P[:6], np.array([[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]]), opens=[[0, 1, 2]])
    new = _mesh(P, T, opens=[[0, 1, 2]])
    fo, fn = tmp_path / "obc.14", tmp_path / "polished.14"
    write_fort14(old, fo)
    write_fort14(new, fn)
    dep = tmp_path / "d_dep.dat"
    dep.write_text(f"Node Number = {len(P)}\n" + "".join(f"{x} {y} 3.0\n" for x, y in P))
    art = {"obc_mesh": str(fo), "polish_mesh": str(fn), "latest_mesh": str(fn)}
    out = pl._stage_depths({"name": "x", "depths": {"product": "external_dep_file",
                                                    "path": str(dep)}}, tmp_path, art,
                           lambda *a: None)
    assert read_fort14(out["final_mesh"]).n_nodes == len(P)
    out2 = pl._stage_depths({"name": "x", "depths": {"product": "none"}}, tmp_path, art,
                            lambda *a: None)
    assert out2["final_mesh"] == str(fn)


def test_r1_9_merge_water_files_projected_source(tmp_path):
    """#9: a projected (UTM) water file is filtered with the bbox in ITS CRS."""
    import geopandas as gpd
    from shapely.geometry import box

    from fvcom_mesh_tools.prep.shoreline import merge_water_files

    a = tmp_path / "utm.gpkg"
    gpd.GeoDataFrame({"fclass": ["water"]}, geometry=[box(250_000, 640_000, 251_000, 641_000)],
                     crs=EPSG).to_file(a)
    out = merge_water_files([a], (0.0, 5.0, 2.0, 7.0), tmp_path / "m.gpkg")
    assert len(gpd.read_file(out)) == 1
    out2 = merge_water_files([a], (3.0, 5.0, 4.0, 7.0), tmp_path / "m2.gpkg")
    assert len(gpd.read_file(out2)) == 0


def test_r1_10_prep_receives_the_recipe_crs(tmp_path, monkeypatch):
    """#10: an explicit recipe crs reaches fmesh-prep-shoreline; none, none."""
    from fvcom_mesh_tools.cli import prepshoreline

    seen = []
    monkeypatch.setattr(prepshoreline, "main", lambda argv: seen.append(list(argv)) or 0)
    base = {"prep": {"bbox": [-1.5, 5.0, 0.5, 6.5]}}
    pl._stage_prep(dict(base, crs=32630), tmp_path, lambda *a: None)
    pl._stage_prep(base, tmp_path, lambda *a: None)
    i = seen[0].index("--utm-epsg")
    assert seen[0][i + 1] == "32630" and "--utm-epsg" not in seen[1]


def test_r1_11_build_records_its_crs(tmp_path, monkeypatch):
    """#11: every build records the CRS of the raw mesh, so a stale projected
    marker from an earlier polygon build cannot skip the projection."""
    import fvcom_mesh_tools.finishing as fin
    from fvcom_mesh_tools.cli import buildmesh

    ll = _mesh(np.array([[0.5, 5.8], [0.51, 5.8], [0.5, 5.81]]), np.array([[0, 1, 2]]))

    def fake_build(argv):
        write_fort14(ll, argv[1])
        return 0

    monkeypatch.setattr(buildmesh, "main", fake_build)
    rec = {"name": "x", "crs": EPSG, "build": {"engine": "oceanmesh", "dem": str(tmp_path),
                                                "coastline": "c.shp", "hmin_m": 100}}
    out = pl._stage_build(rec, tmp_path, {"raw_mesh_crs": f"EPSG:{EPSG}"}, lambda *a: None)
    assert out["raw_mesh_crs"] == "EPSG:4326"
    got = {}

    def fake_finish(mesh, *a, **k):
        got["x"] = mesh.nodes[:, 0].copy()
        return mesh, {}

    monkeypatch.setattr(fin, "finish_constrained_mesh", fake_finish)
    pl._stage_finish(rec, tmp_path, dict(out), lambda *a: None)
    assert got["x"].min() > 1e5                                    # projected to metres


# ------------------------------------------- review newregion r2 (regressions)
def _rect_lattice(nx, ny, dx, dy):
    xs, ys = np.meshgrid(dx * np.arange(nx), dy * np.arange(ny))
    P = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    T = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            T.append([a, a + 1, a + nx])
            T.append([a + 1, a + nx + 1, a + nx])
    return P, np.asarray(T, np.int64)


def test_r2_1_protection_by_element_overlap(tmp_path, monkeypatch):
    """#1: a protected strip narrower than the mesh spacing has no vertex in it;
    its elements are still protected and finish deletes none of them."""
    import geopandas as gpd
    from shapely.geometry import LineString, box

    import fvcom_mesh_tools.finishing as fin

    P, T = _rect_lattice(4, 4, 300.0, 100.0)
    strip = box(0, 20, 900, 80)
    assert not fin.protected_node_mask(P, strip).any()
    assert fin.protected_element_mask(P, T, strip).sum() == 6
    monkeypatch.setattr(fin, "phase_h_optimize", lambda m, **k: (m, {}))
    monkeypatch.setattr(fin, "build_coastline_projector", lambda *a, **k: None)
    shp = tmp_path / "s.shp"
    gpd.GeoDataFrame(geometry=[LineString(box(0, 0, 900, 300).exterior.coords)], crs=EPSG
                     ).to_crs(4326).to_file(shp)
    m = _mesh(P, T, opens=[[12, 13, 14, 15]])
    out, info = fin.finish_constrained_mesh(m, shp, tmp_path / "w", utm_epsg=EPSG,
                                            protect_polygons=strip, log=lambda *a: None)
    assert fin.protected_overlap_m2(out.nodes, out.elements, strip) == pytest.approx(strip.area)


def test_r2_2_siteops_cannot_erode_protected_water():
    """#2: the site operators (collapse at a protected bank) keep the protected
    coverage; the reviewed code lost 5,000 m2 here."""
    from shapely.geometry import box

    from fvcom_mesh_tools.finishing import protected_overlap_m2
    from fvcom_mesh_tools.site_session import apply_site_operators

    P, T = _rect_lattice(4, 4, 100.0, 100.0)
    P[0] = (80.0, 30.0)
    water = box(0, 0, 300, 300)
    m = _mesh(P, T, opens=[[12, 13, 14, 15]])
    before = protected_overlap_m2(m.nodes, m.elements, water)
    from shapely.geometry import LineString

    lines = [LineString(water.exterior.coords)]
    out, _log = apply_site_operators(m, lines, passes=1, protect_polygons=water,
                                     log=lambda *a: None)
    assert protected_overlap_m2(out.nodes, out.elements, water) >= before - 1.0


def test_r2_3_small_protected_component_is_not_dropped_silently():
    """#3: a protected pond below max_dropped_component_km2 fails unless it is
    an expected drop."""
    from shapely.geometry import box

    from fvcom_mesh_tools.mesh_engine.polygon import build_size_field, generate_polygon_mesh

    water, obc, _chan = _domain()
    pond = box(X0 + 8000, Y0, X0 + 8600, Y0 + 600)                # 0.36 km2 < 0.5
    w2 = water.union(pond)
    sz = build_size_field(w2, crs_epsg=EPSG, background_m=200.0, raster_m=50.0,
                          log=lambda *a: None)
    with pytest.raises(RuntimeError, match="protected water"):
        generate_polygon_mesh(w2, sz, obc_line=obc, max_iter=15, keep_polygons=pond,
                              log=lambda *a: None)
    res = generate_polygon_mesh(w2, sz, obc_line=obc, max_iter=15, keep_polygons=pond,
                                expected_drop=pond.buffer(50.0), log=lambda *a: None)
    assert res["points"][:, 0].max() < X0 + 7000
    # the relative tolerance: all of its protected water may go only with frac >= 1
    res = generate_polygon_mesh(w2, sz, obc_line=obc, max_iter=15, keep_polygons=pond,
                                protected_drop_tol_frac=1.0, log=lambda *a: None)
    assert res["points"][:, 0].max() < X0 + 7000                  # dropped, with consent
    assert any(r["in_keep"] for r in res["components"] if r["stage"] == "pre-clean")


def _mock_generator(monkeypatch, P, T, water_bounds=(-150.0, 0.0)):
    """The generator returns (P, T) in the engine's local frame (shifted by the
    origin it derives from the water bounds)."""
    import oceanmesh
    import oceanmesh.clean

    o = np.floor(np.asarray(water_bounds) / 1000.0) * 1000.0 - 5000.0
    monkeypatch.setattr(oceanmesh, "generate_mesh", lambda *a, **k: (P - o, T.copy()))
    monkeypatch.setattr(oceanmesh.clean, "om2d_default_clean",
                        lambda p, t, **k: (p, t))


def _ear_lattice():
    """5 x 5 lattice at 300 m plus one ear element at the west side whose
    centroid lies outside a thin protected strip it overlaps."""
    P, T = _rect_lattice(5, 5, 300.0, 300.0)
    P = np.vstack([P, [[-150.0, 150.0]]])
    T = np.vstack([T, [[0, 25, 5]]])
    return P, T


def test_r2_5_repair_loss_by_intersection_area(monkeypatch):
    """#5: the ear's overlap with the protected strip counts although its
    centroid is outside; with a zero threshold the build fails."""
    from shapely.geometry import box

    from fvcom_mesh_tools.mesh_engine.polygon import build_size_field, generate_polygon_mesh

    P, T = _ear_lattice()
    _mock_generator(monkeypatch, P, T)
    water = box(-150, 0, 1200, 1200)
    strip = box(-150, 0, 150, 60)
    sz = build_size_field(water, crs_epsg=EPSG, background_m=300.0, raster_m=100.0,
                          log=lambda *a: None)
    res = generate_polygon_mesh(water, sz, keep_polygons=strip, log=lambda *a: None)
    assert res["protected_repair_area_m2"] > 1000.0
    with pytest.raises(RuntimeError, match="max_protected_ear_area_km2"):
        generate_polygon_mesh(water, sz, keep_polygons=strip, max_protected_ear_area_km2=0.0,
                              log=lambda *a: None)


def test_r2_4_pipeline_forwards_build_guards(tmp_path, monkeypatch):
    """#4: build.max_protected_ear_area_km2 (and the drop policy) reach the
    engine; an excessive repair loss fails the pipeline build."""
    import geopandas as gpd
    from shapely.geometry import box

    P, T = _ear_lattice()
    _mock_generator(monkeypatch, P, T)
    gp = tmp_path / "d.gpkg"
    gpd.GeoDataFrame(geometry=[box(-150, 0, 1200, 1200)], crs=EPSG).to_file(gp, layer="water")
    gpd.GeoDataFrame(geometry=[box(-150, 0, 150, 60)], crs=EPSG).to_file(gp, layer="strip")
    rec = {"name": "x", "crs": EPSG, "domain": {"water": {"path": str(gp), "layer": "water"}},
           "protect_polygons": [{"path": str(gp), "layer": "strip"}],
           "build": {"engine": "oceanmesh-polygon", "max_protected_ear_area_km2": 0.0,
                     "size": {"background_m": 300.0, "raster_m": 100.0}}}
    art = pl._stage_domain(rec, tmp_path, {}, lambda *a: None)
    with pytest.raises(RuntimeError, match="max_protected_ear_area_km2"):
        pl._stage_build_polygon(rec, tmp_path, art, lambda *a: None)
    rec["build"]["max_protected_ear_area_km2"] = 1.0
    out = pl._stage_build_polygon(rec, tmp_path, art, lambda *a: None)
    assert out["raw_protected_repair_m2"] > 1000.0


def test_r2_6_band_restore_checks_every_c4_pair():
    """#6: restoring a band node may not turn passing pairs into C4 failures
    hidden behind an existing worse pair."""
    from fvcom_mesh_tools.finishing import _restore_obc_band, _tag_obc_band
    from fvcom_mesh_tools.qa import run_qa

    P, T = _rect_lattice(4, 4, 100.0, 100.0)
    P0 = P.copy()
    P0[[5, 6]] += (-40.0, -40.0)
    m0 = _mesh(P0, T, opens=[[0, 1, 2, 3]])
    state = _tag_obc_band(m0)
    P1 = P0.copy()
    P1[5] = P[5]                                     # the optimizer moved node 5 back
    m1 = _mesh(P1, T, opens=[[0, 1, 2, 3]])
    m1.depths = m0.depths.copy()

    def c4(m):
        r = run_qa(m, coords="metric", channel_check=False)
        return [c.n_violations for c in r.checks if c.check_id == "c4_area_change"][0]

    n_before = c4(m1)
    out, info = _restore_obc_band(m1, state, log=lambda *a: None)
    assert c4(out) <= n_before


def test_r2_7_figures_do_not_reproject_a_lonlat_raw_mesh(tmp_path, monkeypatch):
    """#7: a lon/lat raw mesh (recorded EPSG:4326) is plotted where it is."""
    import fvcom_mesh_tools.plotting as plotting

    ll = _mesh(np.array([[139.75, 35.45], [139.76, 35.45], [139.75, 35.46]]),
               np.array([[0, 1, 2]]))
    f = tmp_path / "raw.14"
    write_fort14(ll, f)
    seen = []
    monkeypatch.setattr(plotting, "plot_mesh_overview",
                        lambda mesh, out, **k: seen.append(mesh.nodes.copy()))
    pl._stage_figures({"crs": 32654, "figures": {}}, tmp_path,
                      {"raw_mesh": str(f), "raw_mesh_crs": "EPSG:4326"}, lambda *a: None)
    assert seen and np.allclose(seen[0][0], [139.75, 35.45])


def test_r2_8_prep_receives_a_legacy_crs_override(tmp_path, monkeypatch):
    """#8: finish.utm_epsg different from prep's automatic zone reaches prep;
    when they agree (Tokyo) the prep arguments are unchanged."""
    from fvcom_mesh_tools.cli import prepshoreline

    seen = []
    monkeypatch.setattr(prepshoreline, "main", lambda argv: seen.append(list(argv)) or 0)
    pl._stage_prep({"prep": {"bbox": [-0.15, 5.0, 1.45, 6.5]}, "finish": {"utm_epsg": 32630}},
                   tmp_path, lambda *a: None)
    pl._stage_prep({"prep": {"bbox": [139.4, 34.9, 140.3, 35.9]},
                    "finish": {"utm_epsg": 32654}}, tmp_path, lambda *a: None)
    i = seen[0].index("--utm-epsg")
    assert seen[0][i + 1] == "32630" and "--utm-epsg" not in seen[1]


# ------------------------------------------- review newregion r3 (regressions)
def test_r3_1_finely_meshed_protected_pond_is_measured():
    """#1: a separate protected pond holding MOST of the elements is still a
    dropped component; it fails instead of vanishing (no OBC reaches it)."""
    from pyproj import Transformer
    from shapely.geometry import LineString, box

    from fvcom_mesh_tools.mesh_engine.polygon import build_size_field, generate_polygon_mesh

    sea = box(X0, Y0, X0 + 3000, Y0 + 1333.33)
    obc = LineString([(X0, Y0), (X0 + 3000, Y0)])
    pond = box(X0 + 4000, Y0, X0 + 4600, Y0 + 600)
    w = sea.union(pond)
    inv = Transformer.from_crs(f"EPSG:{EPSG}", "EPSG:4326", always_xy=True)
    lo0, la0 = inv.transform(*pond.buffer(100).bounds[:2])
    lo1, la1 = inv.transform(*pond.buffer(100).bounds[2:])
    reg = [{"name": "pond", "geometry": {"bbox": [lo0, la0, lo1, la1]}, "target_h_m": 60.0}]
    sz = build_size_field(w, crs_epsg=EPSG, background_m=333.0, raster_m=20.0, regions=reg,
                          log=lambda *a: None)
    with pytest.raises(RuntimeError, match="protected water"):
        generate_polygon_mesh(w, sz, obc_line=obc, max_iter=30, ladder_k=1.0,
                              clean_dj_cutoff=0.0005, keep_polygons=pond,
                              log=lambda *a: None)


def _fan_mesh():
    P = np.array([[0.0, 0.0], [1.0, -1.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0],
                  [50.0, 50.0]])
    T = np.array([[0, 1, 5], [1, 2, 5], [2, 3, 5], [3, 4, 5], [4, 0, 5]])
    return P, T


def test_r3_2_weld_keeps_open_boundary_and_protected_water(tmp_path, monkeypatch):
    """#2: the 2 m weld may not merge an open-boundary node or take protected
    water; the reviewed code welded (1,-1) into (0,0)."""
    import geopandas as gpd
    from shapely.geometry import LineString, Polygon

    import fvcom_mesh_tools.finishing as fin

    P, T = _fan_mesh()
    poly = Polygon(P[:5])
    monkeypatch.setattr(fin, "build_coastline_projector", lambda *a, **k: None)
    shp = tmp_path / "s.shp"
    gpd.GeoDataFrame(geometry=[LineString(list(poly.exterior.coords))], crs=EPSG
                     ).to_crs(4326).to_file(shp)
    m = _mesh(P, T, opens=[[1, 2]])
    out, info = fin.finish_constrained_mesh(m, shp, tmp_path / "w", utm_epsg=EPSG,
                                            optimize_budget_s=5.0, protect_polygons=poly,
                                            freeze_obc_band=True, log=lambda *a: None)
    assert np.allclose(out.nodes[out.open_boundaries[0]], [[1.0, -1.0], [100.0, 0.0]])
    assert fin.protected_overlap_m2(out.nodes, out.elements, poly) == pytest.approx(poly.area)


def test_r3_3_patches_cannot_lose_protected_water(monkeypatch):
    """#3: a patch that moves a protected bank node inward is reverted, even
    when it gains water elsewhere (geometric loss, not netted)."""
    from shapely.geometry import box

    import fvcom_mesh_tools.autofinish as af

    P, T = _rect_lattice(4, 4, 100.0, 100.0)
    m = _mesh(P, T, opens=[[12, 13, 14, 15]])
    water = box(0, 0, 300, 300)

    def fake_execute(nodes, cells, patches, obc_nodes=None, log=print, **k):
        n2 = np.asarray(nodes, float).copy()
        n2[1] = (100.0, 20.0)      # bank node pulled in: protected water lost ...
        n2[7] = (320.0, 100.0)     # ... while the east bank gains water
        return n2, [{"id": 0, "outcome": "fixed [test]"}]

    monkeypatch.setattr(af, "execute_patches", fake_execute)
    out, led = pl._execute_patches_guarded(m, [{"id": 0}], [12, 13, 14, 15], 20.0,
                                           lambda *a: None, protect=water, tol_m2=0.0)
    assert np.allclose(out.nodes, P) and led[0]["outcome"].startswith("reverted")


def _r3_lattice():
    P, T = _rect_lattice(5, 5, 100.0, 100.0)
    T = [list(t) for t in T]
    T[0], T[1] = [0, 1, 6], [0, 6, 5]
    T[-2], T[-1] = [18, 19, 24], [18, 24, 23]
    return P + [X0, Y0], np.asarray(T, np.int64)


def _fallback_setup(tmp_path, protect):
    import geopandas as gpd
    from shapely.geometry import LineString, box

    rect = box(X0, Y0, X0 + 400, Y0 + 400)
    shp = tmp_path / "shore.shp"
    gpd.GeoDataFrame(geometry=[LineString(list(rect.exterior.coords))], crs=EPSG
                     ).to_crs(4326).to_file(shp)
    rec = {"crs": EPSG, "finishing": {"auto": {"obc_perp_max_deg": 10.0}}}
    if protect:
        gp = tmp_path / "p.gpkg"
        gpd.GeoDataFrame(geometry=[rect], crs=EPSG).to_file(gp)
        rec["protect_polygons"] = [{"path": str(gp)}]
    return rec, {"domain_shoreline": str(shp)}


def test_r3_4_fallback_keeps_the_open_boundary(tmp_path):
    """#4: without protection the fallback may not return a mesh without its
    open boundary (the reviewed code returned open_boundaries == [])."""
    P, T = _r3_lattice()
    P[8] = (X0 + 343.126008, Y0 + 41.208301)
    P[7] = (X0 + 214.748498, Y0 + 22.858384)
    m = _mesh(P, T, opens=[[0, 1, 2, 3, 4]])
    rec, art = _fallback_setup(tmp_path, protect=False)
    out, info = pl._finishing_site_fallback(rec, m, art, lambda *a: None)
    assert len(out.open_boundaries) == 1
    assert np.allclose(out.nodes[out.open_boundaries[0]], P[[0, 1, 2, 3, 4]])


def test_r3_5_fallback_respects_the_configured_perpendicularity(tmp_path):
    """#5: the whole fallback is held to obc_perp_max_deg (10 deg here), per
    node against the original deviations (all 0 deg)."""
    from fvcom_mesh_tools.finishing import obc_perpendicularity

    P, T = _r3_lattice()
    P[17] = (X0 + 232.394916693, Y0 + 262.252995411)
    P[11] = (X0 + 79.252582251, Y0 + 161.262146908)
    m = _mesh(P, T, opens=[[0, 1, 2, 3, 4]])
    assert np.allclose(obc_perpendicularity(m), 0.0)
    rec, art = _fallback_setup(tmp_path, protect=True)
    out, info = pl._finishing_site_fallback(rec, m, art, lambda *a: None)
    assert obc_perpendicularity(out).max() <= 10.0 + 1e-9


def test_r3_2_weld_merges_coincident_open_boundary_nodes(tmp_path, monkeypatch):
    """#2 (consequence): the OBC assignment can leave two COINCIDENT open-
    boundary nodes (VR v3 job 7591832); the guarded weld merges such a pair
    (a repair, the geometry is unchanged) but never a distinct OBC node."""
    import geopandas as gpd
    from scipy.spatial import cKDTree
    from shapely.geometry import LineString, box

    import fvcom_mesh_tools.finishing as fin

    P, T = _rect_lattice(3, 3, 100.0, 100.0)
    P = np.vstack([P, P[[1]]])                          # node 9 coincides with node 1
    T = T.copy()
    T[0] = [0, 9, 3]                                    # the first element uses the duplicate
    m = _mesh(P, T, opens=[[0, 9, 1, 2]])
    monkeypatch.setattr(fin, "phase_h_optimize", lambda m, **k: (m, {}))
    monkeypatch.setattr(fin, "build_coastline_projector", lambda *a, **k: None)
    shp = tmp_path / "s.shp"
    gpd.GeoDataFrame(geometry=[LineString(box(0, 0, 200, 200).exterior.coords)], crs=EPSG
                     ).to_crs(4326).to_file(shp)
    out, info = fin.finish_constrained_mesh(m, shp, tmp_path / "w", utm_epsg=EPSG,
                                            protect_polygons=box(0, 0, 200, 200),
                                            freeze_obc_band=True, log=lambda *a: None)
    assert not cKDTree(out.nodes).query_pairs(1e-3)
    assert np.allclose(out.nodes[out.open_boundaries[0]], [[0, 0], [100, 0], [200, 0]])


def test_r3_2_lone_node_fixed_by_flip_without_losing_water():
    """#2 (consequence): a single-element node in protected water is fixed by
    flipping the opposite edge; coverage and element count are unchanged."""
    from shapely.geometry import box

    from fvcom_mesh_tools.finishing import protected_overlap_m2, remove_lone_nodes

    P, T = _rect_lattice(3, 3, 100.0, 100.0)         # nodes 0 and 8 are single-element nodes
    m = _mesh(P, T)
    water = box(0, 0, 200, 200)
    out, info = remove_lone_nodes(m, protect_polygons=water, policy="fail", log=lambda *a: None)
    assert info["flipped"] == 2 and info["elements_deleted"] == 0 and info["lone_kept"] == 0
    assert out.n_elements == len(T)
    assert protected_overlap_m2(out.nodes, out.elements, water) == pytest.approx(4e4)


def test_r3_2_lone_node_fixed_by_split_when_flip_would_make_another():
    """A flip that would leave a or b in one element is replaced by a split of
    the opposite edge (VR job 7592080: the lone node at (279934, 658992))."""
    from shapely.geometry import box

    from fvcom_mesh_tools.finishing import protected_overlap_m2, remove_lone_nodes

    P = np.array([[0.0, 0.0], [100.0, 0.0], [100.0, 100.0], [0.0, 100.0]])
    T = np.array([[0, 1, 3], [1, 2, 3]])
    water = box(0, 0, 100, 100)
    out, info = remove_lone_nodes(_mesh(P, T), protect_polygons=water, policy="fail",
                                  log=lambda *a: None)
    cnt = np.bincount(out.elements.ravel(), minlength=out.n_nodes)
    assert (cnt >= 2).all() and info["elements_deleted"] == 0
    assert protected_overlap_m2(out.nodes, out.elements, water) == pytest.approx(1e4)


def _region_line_case(tmp_path):
    """A region recipe (polygon build + domain line) whose water is a 400 x 200 m
    rectangle with the open line along its south edge, and a matching mesh."""
    import geopandas as gpd
    from pyproj import Transformer
    from shapely.geometry import box

    P, T = _rect_lattice(5, 3, 100.0, 100.0)
    P = P + [X0, Y0]
    inv = Transformer.from_crs(f"EPSG:{EPSG}", "EPSG:4326", always_xy=True)
    line = [list(inv.transform(X0, Y0)), list(inv.transform(X0 + 400, Y0))]
    wp = tmp_path / "water.gpkg"
    gpd.GeoDataFrame(geometry=[box(X0, Y0, X0 + 400, Y0 + 200)], crs=EPSG).to_file(wp)
    rec = {"crs": EPSG, "domain": {"obc_line": line},
           "build": {"engine": "oceanmesh-polygon"}}
    return rec, {"domain_water": str(wp)}, P, T


def test_r3_2_obc_stages_check_the_open_line_coverage(tmp_path):
    """#2 (r3) / r4 #3: the obc stages check the open string against the wet
    extent of the supplied line (finish/polish may move open nodes)."""
    rec, art, P, T = _region_line_case(tmp_path)
    full = _mesh(P, T, opens=[[0, 1, 2, 3, 4]])
    assert pl._check_open_cover(rec, tmp_path, art, full, {}, lambda *a: None)["uncovered_m"] < 1
    short = _mesh(P, T, opens=[[0, 1]])                          # 300 m short > 250 m default
    with pytest.raises(SystemExit, match="no longer matches"):
        pl._check_open_cover(rec, tmp_path, art, short, {}, lambda *a: None)
    legacy = {"crs": EPSG, "prep": {"obc_line": rec["domain"]["obc_line"]}}
    assert pl._check_open_cover(legacy, tmp_path, {}, short, {}, lambda *a: None) is None



# ------------------------------------------- review newregion r4 (regressions)
def test_r4_1_obc_aftercare_may_not_delete_protected_water():
    """#1: the reviewer's 10 x 3 lattice: aftercare deleted 4 protected
    elements (20,000 m2) and returned; now the repair is refused."""
    from pyproj import Transformer
    from shapely.geometry import box

    from fvcom_mesh_tools.obc_tools import assign_west_south_obc

    P, T = _rect_lattice(10, 3, 100.0, 100.0)
    P = P + [300000.0, 640000.0]
    m = _mesh(P, T)
    inv = Transformer.from_crs("EPSG:32631", "EPSG:4326", always_xy=True)
    line = [inv.transform(300000.0, 640000.0), inv.transform(300900.0, 640000.0)]
    water = box(300000, 640000, 300900, 640200)
    with pytest.raises(RuntimeError, match="protected water"):
        assign_west_south_obc(m, utm_epsg=32631, snap=False, max_move_m=60.0,
                              obc_line_lonlat=line, protect_polygons=water,
                              log=lambda *a: None)


def test_r4_2_coverage_guard_needs_no_runtime_metadata(tmp_path):
    """#2: a resumed obc/obcfinal run has only persisted artifacts: the guard
    uses the domain (a path artifact), fails loudly without it, and the
    coverage record is persisted (dict artifacts survive provenance)."""
    rec, art, P, T = _region_line_case(tmp_path)
    full = _mesh(P, T, opens=[[0, 1, 2, 3, 4]])
    with pytest.raises(SystemExit, match="domain water"):
        pl._check_open_cover(rec, tmp_path, {}, full, {}, lambda *a: None)
    assert "open_cover" in pl._STRUCTURED_ARTIFACTS and "raw_open_cover" in pl._STRUCTURED_ARTIFACTS


def test_r4_2_open_cover_survives_a_resumed_run(tmp_path):
    """#2: run domain+build+obc, then --only obc: the coverage record is in the
    persisted artifacts both times."""
    import geopandas as gpd
    import yaml

    water, obc, chan = _domain()
    gp = tmp_path / "dom.gpkg"
    gpd.GeoDataFrame(geometry=[water], crs=EPSG).to_file(gp, layer="water")
    gpd.GeoDataFrame(geometry=[obc], crs=EPSG).to_file(gp, layer="obc")
    rec = {"name": "t", "out_dir": str(tmp_path / "out"), "crs": f"EPSG:{EPSG}",
           "domain": {"water": {"path": str(gp), "layer": "water"},
                      "obc_line": {"path": str(gp), "layer": "obc"}},
           "build": {"engine": "oceanmesh-polygon", "max_iter": 30,
                     "size": {"background_m": 400.0, "raster_m": 50.0}},
           "obc": {"trim": 1, "snap": False}}
    f = tmp_path / "r.yaml"
    f.write_text(yaml.safe_dump({"pipeline": rec}))
    assert pl.main([str(f), "--section", "pipeline"]) == 0
    prov = json.loads((tmp_path / "out" / "pipeline_provenance.json").read_text())
    assert prov["artifacts"]["open_cover"]["stage"] == "obc"
    assert pl.main([str(f), "--section", "pipeline", "--only", "obc"]) == 0
    prov = json.loads((tmp_path / "out" / "pipeline_provenance.json").read_text())
    assert prov["artifacts"]["open_cover"]["stage"] == "obc"


def test_r4_3_a_shortened_open_string_fails_and_coverage_may_not_retreat(tmp_path):
    """#3: the allowance is the wet extent with an explicit end tolerance, not
    the trimmed raw mesh; a later stage may not retreat from earlier coverage."""
    rec, art, P, T = _region_line_case(tmp_path)
    short = _mesh(P, T, opens=[[0, 1, 2, 3]])                    # 100 m short at the east
    with pytest.raises(SystemExit, match="end_tol_m"):
        pl._check_open_cover(rec, tmp_path, art, short, {"end_tol_m": 50.0}, lambda *a: None)
    ok = pl._check_open_cover(rec, tmp_path, art, short, {"end_tol_m": 150.0}, lambda *a: None)
    assert ok["uncovered_m"] == pytest.approx(100.0)
    art2 = dict(art, open_cover={"stage": "obc", "uncovered_m": 0.0})
    with pytest.raises(SystemExit, match="retreated"):
        pl._check_open_cover(rec, tmp_path, art2, short, {"end_tol_m": 150.0}, lambda *a: None)


def test_r4_4_region_siteops_keeps_the_open_boundary(tmp_path):
    """#4: standalone siteops of a region recipe without protection keeps (and
    validates) the open boundary; the reviewed code cleared it."""
    P, T = _r3_lattice()
    m = _mesh(P, T, opens=[[0, 1, 2, 3, 4]])
    f = tmp_path / "obc.14"
    write_fort14(m, f)
    rec, art = _fallback_setup(tmp_path, protect=False)
    rec.update(name="x", build={"engine": "oceanmesh-polygon"},
               domain={"obc_line": [[0.1, 5.7], [0.2, 5.7]]}, siteops={"passes": 1})
    out = pl._stage_siteops(rec, tmp_path, dict(art, obc_mesh=str(f)), lambda *a: None)
    m2 = read_fort14(out["siteops_mesh"])
    assert len(m2.open_boundaries) == 1
    assert np.allclose(m2.nodes[m2.open_boundaries[0]], P[[0, 1, 2, 3, 4]])


def test_r4_4_validation_accepts_a_coincident_open_boundary_pair():
    """The OBC assignment can leave two coincident open-boundary nodes; the
    string validation treats them as one node (VR job 7592472)."""
    from fvcom_mesh_tools.finishing import validate_open_boundary

    P, T = _rect_lattice(3, 3, 100.0, 100.0)
    P = np.vstack([P, P[[1]]])                          # node 9 coincides with node 1
    T = T.copy()
    T[0] = [0, 9, 3]
    m = _mesh(P, T, opens=[[0, 9, 1, 2]])
    assert validate_open_boundary(m, m, same_positions=True) == []
