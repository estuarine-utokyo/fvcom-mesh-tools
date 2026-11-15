"""fmesh-pipeline: recipe-driven end-to-end mesh construction.

Replaces the PoC-era sed-cloned notebook scripts with a single
declarative entry point:

    fmesh-pipeline recipes/tokyo_bay_v5.yaml [--only STAGE[,STAGE]]

Stages run in fixed order, each at most ONCE (no convergence loops),
each optional by presence of its key in the recipe:

    prep     -> fmesh-prep-shoreline (land opening + skeleton seeds)
    domain   -> a region-supplied water polygon + open-boundary line
    build    -> fmesh-buildmesh (oceanmesh engine, CDT boundary), or the
                projected polygon engine (``engine: oceanmesh-polygon``)
    finish   -> finishing.finish_constrained_mesh (UTM)
    obc / siteops / polish / obcfinal / finishing -> quality chain
    depths   -> depth product (m7001_production | external_dep_file | none)
    qa       -> fvcom_mesh_tools.qa.run_qa report (ja/en)
    figures  -> per-stage overview/zoom PNGs

The projected CRS is the top-level ``crs`` (EPSG code or "EPSG:nnnnn");
without it the legacy ``finish.utm_epsg``, else the UTM zone of
``prep.bbox`` or of the domain polygon is used. ``protect_polygons``
(top level) names water that no stage may cut off (channels).

A pipeline_provenance.json in the output directory records the
recipe, git state and per-stage results.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path
from typing import Any


# jsasaki 2026-11-16: optional section (a block of a region repository's recipe)
def _load_recipe(path: Path, section: str | None = None) -> dict[str, Any]:
    import yaml

    with path.open() as f:
        recipe = yaml.safe_load(f)
    if section is not None:
        if not isinstance(recipe, dict) or not isinstance(recipe.get(section), dict):
            raise SystemExit(f"recipe {path} has no mapping {section!r}")
        recipe = recipe[section]
    if not isinstance(recipe, dict):
        raise SystemExit(f"recipe {path} is not a mapping")
    for key in ("name", "out_dir"):
        if key not in recipe:
            raise SystemExit(f"recipe missing required key: {key!r}")
    return recipe


#---> jsasaki 2026-11-16: new-region support (crs, domain polygon, protection, paths)
def _parse_epsg(value) -> int:
    s = str(value).strip().upper()
    if s.startswith("EPSG:"):
        s = s[5:]
    try:
        code = int(s)
    except ValueError:
        raise SystemExit(f"crs must be an EPSG code (e.g. 32631 or 'EPSG:32631'), not {value!r}")
    return code


def _recipe_epsg(recipe: dict[str, Any]) -> int:
    """The projected CRS of the recipe: ``crs``; else the legacy
    ``finish.utm_epsg``; else the UTM zone of ``prep.bbox``; else that of
    the domain polygon. There is no Tokyo fallback."""
    if recipe.get("crs") is not None:
        code = _parse_epsg(recipe["crs"])
        legacy = (recipe.get("finish") or {}).get("utm_epsg")
        if legacy is not None and int(legacy) != code:
            raise SystemExit(f"crs {code} and finish.utm_epsg {legacy} disagree")
        return code
    legacy = (recipe.get("finish") or {}).get("utm_epsg")
    if legacy is not None:
        return int(legacy)
    from fvcom_mesh_tools.prep.shoreline import auto_utm_epsg

    bbox = (recipe.get("prep") or {}).get("bbox")
    if bbox:
        return auto_utm_epsg(0.5 * (bbox[0] + bbox[2]), 0.5 * (bbox[1] + bbox[3]))
    dom = (recipe.get("domain") or {}).get("water")
    if dom:
        g = _read_geoms(dom).to_crs(4326)
        c = g.union_all().centroid
        return auto_utm_epsg(c.x, c.y)
    raise SystemExit("recipe names no crs, finish.utm_epsg, prep.bbox or domain.water: "
                     "the projected CRS is unknown")


def _path(value) -> Path:
    """A recipe path with environment variables expanded; an unresolved
    variable fails loudly (GENKAI rule)."""
    import os

    s = os.path.expandvars(str(value))
    if "$" in s:
        raise SystemExit(f"unresolved environment variable in path: {s}")
    return Path(s)


def _read_geoms(spec):
    """GeoDataFrame from ``{path, layer?}`` (or a bare path string)."""
    import geopandas as gpd

    if isinstance(spec, (str, Path)):
        spec = {"path": spec}
    p = _path(spec["path"])
    if not p.exists():
        raise SystemExit(f"geometry file {p} does not exist")
    g = gpd.read_file(p, layer=spec.get("layer")) if spec.get("layer") else gpd.read_file(p)
    if g.crs is None:
        raise SystemExit(f"geometry file {p} has no CRS")
    return g


def _protect_geometry(recipe, epsg):
    """Union of ``protect_polygons`` in the recipe CRS (None if absent).
    Items: ``{path, layer?, buffer_m?}`` or a lon/lat coordinate list."""
    items = recipe.get("protect_polygons")
    if not items:
        return None
    import geopandas as gpd
    from shapely.geometry import Polygon
    from shapely.ops import unary_union

    parts = []
    for it in items:
        if isinstance(it, dict):
            g = _read_geoms(it).to_crs(epsg)
            geom = unary_union(list(g.geometry))
            if it.get("buffer_m"):
                geom = geom.buffer(float(it["buffer_m"]))
        else:
            geom = gpd.GeoSeries([Polygon([tuple(q) for q in it])], crs=4326
                                 ).to_crs(epsg).iloc[0]
        parts.append(geom)
    return unary_union(parts)


def _shoreline(recipe, artifacts):
    """Engineered shoreline for projectors (lon/lat file)."""
    return Path(artifacts.get("projector_lines")
                or artifacts.get("land_opened")
                or artifacts.get("domain_shoreline")
                or _path(recipe["build"]["coastline"]))


def _check_recipe(recipe):
    """Recipe rules the stages rely on."""
    has_line = bool((recipe.get("prep") or {}).get("obc_line")
                    or (recipe.get("domain") or {}).get("obc_line"))
    edits = [k for k in ("finish", "siteops", "polish") if k in recipe]
    # jsasaki 2026-11-16: review newregion r1 #3: the polygon build may trim the open string at
    # R4 corners (its ends become land): it needs a stage that re-derives the open boundary too
    if has_line and (recipe.get("build") or {}).get("engine") == "oceanmesh-polygon":
        edits = ["build (oceanmesh-polygon)"] + edits
    if has_line and edits and not ("obc" in recipe or "obcfinal" in recipe):
        raise SystemExit(f"the recipe has an open-boundary line and {edits}, which may "
                         "drop, move or trim open-boundary nodes: add the obc and/or obcfinal "
                         "stage")
#<--- jsasaki 2026-11-16


def _git_describe(repo: Path) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(repo), "describe", "--always", "--dirty"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:  # noqa: BLE001 - provenance is best-effort
        return "unknown"


# jsasaki 2026-11-16: review newregion r2 #8: CRS prep must use when the recipe overrides its zone
def _prep_crs_override(recipe):
    """The recipe CRS when it is explicit (``crs``) or when the legacy
    ``finish.utm_epsg`` differs from the UTM zone prep derives from its bbox;
    else None (prep keeps its automatic zone, as before)."""
    from fvcom_mesh_tools.prep.shoreline import auto_utm_epsg

    if recipe.get("crs") is not None:
        return _recipe_epsg(recipe)
    legacy = (recipe.get("finish") or {}).get("utm_epsg")
    bbox = (recipe.get("prep") or {}).get("bbox")
    if legacy is not None and bbox:
        auto = auto_utm_epsg(0.5 * (bbox[0] + bbox[2]), 0.5 * (bbox[1] + bbox[3]))
        if int(legacy) != auto:
            return int(legacy)
    return None


def _stage_prep(recipe, out_dir, log):
    from fvcom_mesh_tools.cli import prepshoreline

    cfg = recipe["prep"]
    prep_dir = out_dir / "prep"
    argv = [
        "--bbox", *[str(v) for v in cfg["bbox"]],
        "--out-dir", str(prep_dir),
        "--r-open", str(cfg.get("r_open_m", 150.0)),
        "--min-island", str(cfg.get("min_island_area_m2", 3.6e5)),
        "--px", str(cfg.get("skeleton_px_m", 50.0)),
        "--half-width",
        str(cfg.get("skeleton_half_width_m", [140.0, 460.0])[0]),
        str(cfg.get("skeleton_half_width_m", [140.0, 460.0])[1]),
    ]
    if cfg.get("cache_dir"):
        argv += ["--cache-dir", str(cfg["cache_dir"])]
    if cfg.get("min_water_area") is not None:
        argv += ["--min-water-area", str(cfg["min_water_area"])]
    v6_outer = bool(cfg.get("outer_smooth_r_m"))
    if cfg.get("obc_line") and not v6_outer:
        argv += ["--obc-line"] + [
            str(c) for pt in cfg["obc_line"] for c in pt
        ]
    ir = (recipe.get("build", {}) or {}).get("interest_region")
    if ir:
        argv += ["--interest-region"] + [
            str(c) for pt in ir for c in pt
        ]
        if cfg.get("outside_simplify_m") is not None:
            argv += ["--outside-simplify",
                     str(cfg["outside_simplify_m"])]
    if not cfg.get("skeleton", True):
        argv += ["--no-skeleton"]
    # jsasaki 2026-11-16: region water polygons (list) instead of the Kanto default
    for w in cfg.get("water_shp") or []:
        argv += ["--water-shp", str(_path(w))]
    # jsasaki 2026-11-16: review newregion r1 #10: an explicit recipe crs reaches prep too (a
    # recipe without one keeps prep's own UTM zone, the measured Tokyo behaviour)
    # jsasaki 2026-11-16: review newregion r2 #8: also the legacy finish.utm_epsg when it differs
    # from prep's automatic zone (equal for Tokyo: argv unchanged)
    if _prep_crs_override(recipe) is not None:
        argv += ["--utm-epsg", str(_prep_crs_override(recipe))]
    rc = prepshoreline.main(argv)
    if rc:
        raise SystemExit(f"prep stage failed (rc={rc})")
    out = {
        "land_opened": str(prep_dir / "land_opened.shp"),
        "skeleton_seeds": str(prep_dir / "skeleton_seeds.shp"),
    }
    r_sm = cfg.get("outer_smooth_r_m")
    # jsasaki 2026-11-16: the recipe CRS replaces the literal 32654 below
    epsg = _recipe_epsg(recipe) if r_sm else None
    if r_sm:
        # v6 ordering: smooth FIRST, extend the OBC ends against the
        # SMOOTHED coast (the mouth coasts in the final mesh follow
        # it; extending against the detailed coast left the junction
        # corner ~1 element off the effective line), then cut the
        # wall into BOTH coastlines with the SAME effective line.
        import json as _json

        import geopandas as gpd

        from fvcom_mesh_tools.prep.shoreline import (
            cut_domain_at_obc_line,
            extend_obc_ends_perpendicular,
            projector_lines,
            smooth_land_per_polygon,
        )

        inner_poly = (recipe.get("build", {}) or {}).get(
            "nests", {}).get("inner", {}).get("polygon")
        land = gpd.read_file(out["land_opened"])  # no wall yet
        r_in = cfg.get("inner_smooth_r_m")
        if r_in:
            # sample-like profile: round coves/piers on the DETAILED
            # coastline too (goto2023's SMS-simplified shoreline)
            # jsasaki 2026-11-16: explicit recipe CRS (was the function's 32654 default)
            land = smooth_land_per_polygon(
                land, r_m=float(r_in),
                min_area_m2=float(cfg.get("min_island_area_m2",
                                          3.6e5)),
                utm_epsg=epsg,
            )
            land.to_file(out["land_opened"])
            log(f"[prep] inner coastline smoothed r={r_in:g} m "
                f"({len(land)} polygons)")
        # jsasaki 2026-11-16: explicit recipe CRS (was the function's 32654 default)
        outer_gdf = smooth_land_per_polygon(land, r_m=float(r_sm), utm_epsg=epsg)
        outer_pre_cut = outer_gdf.copy()
        eff = None
        if cfg.get("obc_line"):
            # jsasaki 2026-11-16: review newregion r1 #10: explicit recipe crs (None = UTM zone of
            # the line, as before)
            eff = extend_obc_ends_perpendicular(
                [tuple(q) for q in cfg["obc_line"]], outer_gdf,
                # jsasaki 2026-11-16: review newregion r2 #8: legacy override too
                utm_epsg=_prep_crs_override(recipe),
            )
            bbox = tuple(cfg["bbox"])
            land = cut_domain_at_obc_line(land, eff, bbox)
            outer_gdf = cut_domain_at_obc_line(outer_gdf, eff, bbox)
            land.to_file(out["land_opened"])
            prov_path = prep_dir / "provenance.json"
            prov = (_json.loads(prov_path.read_text())
                    if prov_path.exists() else {})
            prov["obc_line_effective"] = [list(q) for q in eff]
            from fvcom_mesh_tools.prep.shoreline import (
                junction_tail_constraints,
            )

            # jsasaki 2026-11-16: explicit recipe CRS (was the function's 32654 default)
            jc = junction_tail_constraints(
                eff, outer_pre_cut,
                tail_len_m=float(cfg.get("junction_tail_m", 1500.0)),
                utm_epsg=epsg,
            )
            if jc:
                prov["junction_constraints"] = jc
                log(f"[prep] junction tail constraints: {len(jc)}")
            prov_path.write_text(_json.dumps(prov, indent=2))
            log(f"[prep] wall cut into both coastlines "
                f"({len(eff)} effective pts)")
        outer_path = prep_dir / "land_outer.shp"
        outer_gdf.to_file(outer_path)
        out["land_outer"] = str(outer_path)
        log(f"[prep] land_outer.shp: {len(outer_gdf)} polygons "
            f"(r={r_sm:g} m)")
        solid_p = prep_dir / "land_solid.shp"
        if solid_p.exists():
            # tripwire: land converted to water by prep, deep
            # interior only, > 0.25 km2 -> warn with grid cells
            from shapely import unary_union

            from fvcom_mesh_tools.gridref import GridRef

            # jsasaki 2026-11-16: recipe CRS (was a literal 32654)
            solid = gpd.read_file(solid_p).to_crs(epsg)
            opened_u = unary_union(list(
                land.to_crs(epsg).geometry))
            lost = unary_union(list(solid.geometry)).difference(
                opened_u).buffer(-200.0)
            grid = GridRef(*cfg["bbox"])
            import geopandas as _gpd2
            for geom in getattr(lost, "geoms", [lost]):
                if geom.is_empty or geom.area < 2.5e5:
                    continue
                # jsasaki 2026-11-16: recipe CRS (was a literal 32654)
                cen = _gpd2.GeoSeries([geom.centroid], crs=epsg
                                      ).to_crs(4326).iloc[0]
                try:
                    cell = grid.point_to_cell(cen.x, cen.y)
                except ValueError:
                    cell = "?"
                log(f"[prep] WARNING: {geom.area/1e6:.2f} km2 of "
                    f"solid land became water near {cell} "
                    f"({cen.x:.3f},{cen.y:.3f})")
        if inner_poly:
            proj_gdf = projector_lines(land, outer_gdf, inner_poly)
            proj_path = prep_dir / "projector_lines.shp"
            proj_gdf.to_file(proj_path)
            out["projector_lines"] = str(proj_path)
            log(f"[prep] projector_lines.shp: {len(proj_gdf)} lines")
    return out


#---> jsasaki 2026-11-16: region-supplied domain and the projected polygon engine
def _domain_obc_lonlat(recipe):
    """``domain.obc_line``: lon/lat pairs, or ``{path, layer?}`` holding one
    line (any CRS), or a csv with lon,lat columns -> [(lon, lat), ...]."""
    spec = (recipe.get("domain") or {}).get("obc_line")
    if not spec:
        return None
    if isinstance(spec, list):
        return [tuple(map(float, q)) for q in spec]
    p = _path(spec["path"] if isinstance(spec, dict) else spec)
    if p.suffix.lower() == ".csv":
        import pandas as pd

        from fvcom_mesh_tools.base_recipe import read_open_boundary

        df = pd.read_csv(p)
        if {"lon", "lat"} <= set(df.columns):
            return [(float(a), float(b)) for a, b in zip(df.lon, df.lat)]
        return read_open_boundary(p)
    from shapely.ops import linemerge, unary_union

    g = _read_geoms(spec if isinstance(spec, dict) else {"path": p}).to_crs(4326)
    line = unary_union(list(g.geometry))
    if line.geom_type == "MultiLineString":
        line = linemerge(line)
    if line.geom_type != "LineString":
        raise SystemExit(f"domain.obc_line {p}: expected one line, got {line.geom_type}")
    return [tuple(c[:2]) for c in line.coords]


def _stage_domain(recipe, out_dir, artifacts, log):
    """Region-supplied domain: the water polygon to mesh (``domain.water``,
    any CRS) and the open-boundary line (``domain.obc_line``). Writes the
    polygon in the recipe CRS and its boundary as lon/lat lines (the
    shoreline the finishing projectors slide along)."""
    import geopandas as gpd
    from shapely.geometry import LineString
    from shapely.ops import unary_union

    cfg = recipe["domain"]
    epsg = _recipe_epsg(recipe)
    ddir = out_dir / "domain"
    ddir.mkdir(parents=True, exist_ok=True)
    water = unary_union(list(_read_geoms(cfg["water"]).to_crs(epsg).geometry))
    if water.is_empty or water.geom_type not in ("Polygon", "MultiPolygon"):
        raise SystemExit("domain.water is not a polygon")
    wpath = ddir / "water_crs.gpkg"
    gpd.GeoDataFrame(geometry=[water], crs=epsg).to_file(wpath, driver="GPKG")
    rings = [r for q in getattr(water, "geoms", [water])
             for r in [q.exterior] + list(q.interiors)]
    spath = ddir / "shoreline_ll.shp"
    gpd.GeoDataFrame(geometry=[LineString(r.coords) for r in rings], crs=epsg
                     ).to_crs(4326).to_file(spath)
    out = {"domain_water": str(wpath), "domain_shoreline": str(spath)}
    line = _domain_obc_lonlat(recipe)
    if line:
        lpath = ddir / "obc_line_lonlat.json"
        lpath.write_text(json.dumps([list(q) for q in line]))
        out["domain_obc_line"] = str(lpath)
    log(f"[domain] water {water.area / 1e6:.1f} km2 in EPSG:{epsg}, {len(rings)} rings; "
        f"open-boundary line: {len(line) if line else 0} vertices")
    return out


# jsasaki 2026-11-16: review newregion r2 #3: expected_drop = listed geometries + (optionally) the
# water the mesh cannot resolve
def _expected_drop_geometry(recipe, cfg, water, obc, epsg, log):
    """``build.expected_drop`` items plus, with ``build.expected_drop_unresolved_m: w``,
    the water not reached by the open-boundary component of the water opened
    at width ``w`` (channels, necks and ponds narrower than the mesh can carry)."""
    from shapely.ops import unary_union

    geoms = []
    listed = _protect_geometry({"protect_polygons": cfg.get("expected_drop")}, epsg)
    if listed is not None:
        geoms.append(listed)
    w = cfg.get("expected_drop_unresolved_m")
    if w:
        w = float(w)
        opened = water.buffer(-w / 2).buffer(w / 2)
        parts = [g for g in getattr(opened, "geoms", [opened]) if not g.is_empty]
        if obc is not None:
            main = [g for g in parts if g.distance(obc) <= w]
        else:
            main = [max(parts, key=lambda g: g.area)] if parts else []
        unresolved = water.difference(unary_union(main)) if main else water
        log(f"[polygon] unresolved water at {w:g} m: {unresolved.area / 1e6:.2f} km2 "
            "(expected drop)")
        geoms.append(unresolved)
    return unary_union(geoms) if geoms else None


def _stage_build_polygon(recipe, out_dir, artifacts, log):
    """``build.engine: oceanmesh-polygon`` (mesh_engine.polygon)."""
    import numpy as np
    from pyproj import Transformer
    from shapely.geometry import LineString
    from shapely.ops import unary_union

    from fvcom_mesh_tools.finishing import rebuild_land_boundaries
    from fvcom_mesh_tools.io import Fort14Mesh, write_fort14
    from fvcom_mesh_tools.mesh_engine.polygon import (
        build_size_field,
        generate_polygon_mesh,
    )

    cfg = recipe["build"]
    epsg = _recipe_epsg(recipe)
    if "domain_water" not in artifacts:
        raise SystemExit("engine oceanmesh-polygon needs the domain stage (domain.water)")
    import geopandas as gpd

    water = unary_union(list(gpd.read_file(artifacts["domain_water"]).geometry))
    sz = cfg.get("size") or {}
    regions = sz.get("regions")
    if sz.get("regions_file"):
        from fvcom_mesh_tools.sizing import load_sizing

        regions = (regions or []) + load_sizing(_path(sz["regions_file"]))["regions"]
    size = build_size_field(
        water, crs_epsg=epsg,
        raster=_path(sz["raster"]) if sz.get("raster") else None,
        background_m=sz.get("background_m"), raster_m=float(sz.get("raster_m", 50.0)),
        regions=regions, hmin_m=sz.get("hmin_m"), max_m=sz.get("max_m"),
        gradation=sz.get("gradation"), log=log)
    line = _obc_line_from(recipe, out_dir)
    obc = None
    if line:
        fwd = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
        ax, ay = fwd.transform([q[0] for q in line], [q[1] for q in line])
        obc = LineString(list(zip(ax, ay)))
    res = generate_polygon_mesh(
        water, size, obc_line=obc,
        ladder_k=float(cfg.get("obc_ladder_k", 0.0)),
        dm_scale=float(sz.get("dm_scale", 1.0)),
        max_iter=int(cfg.get("max_iter", 60)), seed=int(cfg.get("seed", 0)),
        clean_dj_cutoff=float(cfg.get("clean_dj_cutoff", 0.25)),
        max_dropped_component_km2=float(cfg.get("max_dropped_component_km2", 0.5)),
        keep_polygons=_protect_geometry(recipe, epsg),
        expected_drop=_expected_drop_geometry(recipe, cfg, water, obc, epsg, log),
        # jsasaki 2026-11-16: review newregion r2 #4: the threshold and drop policy reach the engine
        max_protected_ear_area_km2=cfg.get("max_protected_ear_area_km2"),
        protected_drop_policy=str(cfg.get("protected_drop_policy", "fail")),
        protected_drop_tol_m2=float(cfg.get("protected_drop_tol_m2", 1.0)),
        protected_drop_tol_frac=float(cfg.get("protected_drop_tol_frac", 0.01)),
        log=log)
    pts, cls = res["points"], res["cells"]
    mesh = Fort14Mesh(
        title=f"{recipe['name']} polygon raw (EPSG:{epsg})",
        nodes=np.asarray(pts, float),
        depths=np.full(len(pts), float(cfg.get("placeholder_depth_m", 10.0))),
        elements=np.asarray(cls, np.int64),
        open_boundaries=[np.asarray(res["open_string"], np.int64)]
        if len(res["open_string"]) else [],
        land_boundaries=[])
    mesh = rebuild_land_boundaries(mesh)
    raw14 = out_dir / f"{recipe['name']}_raw.14"
    write_fort14(mesh, raw14)
    (out_dir / "build_components.json").write_text(json.dumps(res["components"], indent=1))
    (out_dir / "build_size.json").write_text(json.dumps(size["report"], indent=1))
    log(f"[pipeline] polygon raw -> {raw14} (EPSG:{epsg} metres)")
    # jsasaki 2026-11-16: review newregion r1 #3: record the open-line coverage of the raw mesh
    return {"raw_mesh": str(raw14), "raw_mesh_crs": f"EPSG:{epsg}",
            "raw_open_cover": res["open_cover"],
            "raw_protected_repair_m2": res["protected_repair_area_m2"],
            # jsasaki 2026-11-16: review newregion r2 #5: repair loss outside expected_drop
            "raw_protected_repair_unexpected_m2": res["protected_repair_unexpected_m2"]}
#<--- jsasaki 2026-11-16


def _stage_build(recipe, out_dir, artifacts, log):
    import os

    from fvcom_mesh_tools.cli import buildmesh

    cfg = recipe["build"]
    # jsasaki 2026-11-16: projected polygon engine for region-supplied domains
    if cfg.get("engine") == "oceanmesh-polygon":
        return _stage_build_polygon(recipe, out_dir, artifacts, log)
    raw14 = out_dir / f"{recipe['name']}_raw.14"
    coast = artifacts.get("land_opened") or cfg.get("coastline")
    dem = os.path.expandvars(str(cfg["dem"]))
    if "$" in dem:
        raise SystemExit(
            f"unresolved environment variable in dem path: {dem} "
            "(is DATA_DIR set? GENKAI rule: fail loudly)"
        )
    if cfg.get("engine") == "oceanmesh-multiscale":
        import numpy as np

        from fvcom_mesh_tools.io import Fort14Mesh, write_fort14
        from fvcom_mesh_tools.mesh_engine.multiscale import (
            build_multiscale,
        )

        nests = cfg.get("nests", {})
        jc = None
        prov_path = out_dir / "prep" / "provenance.json"
        if prov_path.exists():
            import json as _json

            jc = _json.loads(prov_path.read_text()).get(
                "junction_constraints")
        pts, cls, dep = build_multiscale(
            land_detailed_shp=Path(artifacts["land_opened"]),
            land_outer_shp=Path(artifacts["land_outer"]),
            dem_path=Path(dem),
            bbox=tuple(recipe["prep"]["bbox"]),
            inner_polygon=nests["inner"]["polygon"],
            outer=nests.get("outer"),
            inner={k: v for k, v in nests.get("inner", {}).items()
                   if k != "polygon"},
            courant=cfg.get("courant"),
            seed=int(cfg.get("seed", 0)),
            junction_constraints=jc,
            log=log,
        )
        mesh = Fort14Mesh(
            title=f"{recipe['name']} multiscale raw",
            nodes=np.asarray(pts, dtype=float),
            depths=np.asarray(dep, dtype=float),
            elements=np.asarray(cls, dtype=np.int64),
            open_boundaries=[],
            land_boundaries=[],
        )
        write_fort14(mesh, raw14)
        log(f"[pipeline] multiscale raw -> {raw14}")
        # jsasaki 2026-11-16: review newregion r1 #11: every build records its CRS
        return {"raw_mesh": str(raw14), "raw_mesh_crs": "EPSG:4326"}
    argv = [
        dem, str(raw14),
        "--engine", cfg.get("engine", "oceanmesh"),
    ]
    prep_bbox = (recipe.get("prep", {}) or {}).get("bbox")
    if prep_bbox:
        argv += ["--region", *[str(v) for v in prep_bbox]]
    obc_line = (recipe.get("prep", {}) or {}).get("obc_line")
    if obc_line and cfg.get("obc_coarsen", True):
        argv += ["--om-obc-coarsen-line"] + [
            str(c) for pt in obc_line for c in pt
        ]
        argv += ["--om-obc-coarsen-size-m",
                 str(cfg.get("obc_coarsen_size_m", 1600.0)),
                 "--om-obc-coarsen-radius-m",
                 str(cfg.get("obc_coarsen_radius_m", 10000.0))]
    if cfg.get("courant_sizing"):
        argv += ["--om-courant-sizing",
                 "--om-courant-target",
                 str(cfg.get("courant_target", 0.5)),
                 "--om-courant-timestep",
                 str(cfg.get("courant_timestep_s", 16.0))]
    if cfg.get("interest_region"):
        argv += ["--om-interest-region"] + [
            str(c) for pt in cfg["interest_region"] for c in pt
        ]
        argv += ["--om-outside-min-m",
                 str(cfg.get("outside_min_m", 1000.0)),
                 "--om-outside-blend-m",
                 str(cfg.get("outside_blend_m", 5000.0))]
    argv += [
        "--coastline", str(coast),
        "--hmin", str(cfg["hmin_m"]),
        "--hmax", str(cfg.get("hmax_m", 1500.0)),
        "--zmax", str(cfg.get("zmax", 0.0)),
        "--om-seed", str(cfg.get("seed", 0)),
    ]
    if cfg.get("shoreline_h0_m"):
        argv += ["--om-shoreline-h0-m", str(cfg["shoreline_h0_m"])]
    if cfg.get("minimum_area_mult"):
        argv += ["--om-minimum-area-mult", str(cfg["minimum_area_mult"])]
    if cfg.get("enforce_hmin_floor", True):
        argv += ["--om-enforce-hmin-floor"]
    if cfg.get("constrain_boundary", True):
        argv += ["--om-constrain-boundary"]
    seeds = artifacts.get("skeleton_seeds")
    if seeds and Path(seeds).exists():
        argv += ["--om-high-fidelity-lines", str(seeds)]
    for extra in cfg.get("extra_args", []):
        argv.append(str(extra))
    rc = buildmesh.main(argv)
    if rc:
        raise SystemExit(f"build stage failed (rc={rc})")
    # jsasaki 2026-11-16: review newregion r1 #11: every build records its CRS
    return {"raw_mesh": str(raw14), "raw_mesh_crs": "EPSG:4326"}


def _stage_finish(recipe, out_dir, artifacts, log):
    import numpy as np
    from pyproj import Transformer

    from fvcom_mesh_tools.finishing import finish_constrained_mesh
    from fvcom_mesh_tools.io import read_fort14, write_fort14

    cfg = recipe.get("finish", {}) or {}
    # jsasaki 2026-11-16: recipe CRS (crs / finish.utm_epsg / derived), no 32654 fallback
    utm = _recipe_epsg(recipe)
    mesh = read_fort14(Path(artifacts["raw_mesh"]))
    # jsasaki 2026-11-16: a raw mesh already in the recipe CRS (polygon engine) is not reprojected
    # jsasaki 2026-11-16: review newregion r1 #11: transform FROM the CRS the build recorded
    # (legacy provenance without it: lon/lat)
    src_crs = artifacts.get("raw_mesh_crs") or "EPSG:4326"
    if src_crs != f"EPSG:{utm}":
        tr = Transformer.from_crs(src_crs, f"EPSG:{utm}", always_xy=True)
        x, y = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        mesh.nodes = np.column_stack([x, y])
    # jsasaki 2026-11-16: shoreline lookup includes the region-supplied domain
    shoreline = _shoreline(recipe, artifacts)
    mesh, finfo = finish_constrained_mesh(
        mesh, shoreline, out_dir / "work",
        optimize_budget_s=float(cfg.get("optimize_budget_s", 600.0)),
        weld_tol_m=float(cfg.get("weld_tol_m", 2.0)),
        utm_epsg=utm,
        #---> jsasaki 2026-11-16: channel protection, valence, lone nodes (all off by default)
        protect_polygons=_protect_geometry(recipe, utm),
        on_dropped_protected=str(cfg.get("on_dropped_protected", "fail")),
        max_valence=cfg.get("max_valence"),
        valence_min_angle_deg=float(cfg.get("valence_min_angle_deg", 0.0)),
        remove_lone=bool(cfg.get("remove_lone_nodes", False)),
        freeze_obc_band=bool(cfg.get("freeze_obc_band", False)),
        # jsasaki 2026-11-16: review newregion r3 #2: protected-loss tolerance (default 0 m2)
        protected_loss_tol_m2=float(cfg.get("protected_loss_tol_m2", 0.0)),
        #<--- jsasaki 2026-11-16
        log=log,
    )
    out14 = out_dir / f"{recipe['name']}_finished.14"
    write_fort14(mesh, out14)
    return {"finished_mesh": str(out14), "finish_info": finfo}


def _obc_line_from(recipe, out_dir):
    """Effective OBC line: prep may extend the ends perpendicular to
    the coast; prefer its recorded result."""
    prov = out_dir / "prep" / "provenance.json"
    if prov.exists():
        import json as _json

        d = _json.loads(prov.read_text())
        eff = d.get("obc_line_effective")
        if eff:
            return [tuple(q) for q in eff]
    # jsasaki 2026-11-16: a region-supplied domain names its own open-boundary line
    line = (recipe.get("prep", {}) or {}).get("obc_line")
    return line if line else _domain_obc_lonlat(recipe)


#---> jsasaki 2026-11-16: review newregion r3 #2: the obc stages must give back the open line
def _check_open_cover(recipe, out_dir, artifacts, mesh, cfg, log):
    """For a region-built mesh (``raw_open_cover`` recorded by the polygon
    build): the open string must cover the supplied line at least as far as
    the raw mesh did (+1 m) and every open node lie within ``max_move_m`` of
    it; else the stage fails (QA does not see a shortened open boundary)."""
    import shapely
    from pyproj import Transformer

    raw = artifacts.get("raw_open_cover")
    line_ll = _obc_line_from(recipe, out_dir)
    if not raw or not line_ll or not mesh.open_boundaries:
        return None
    fwd = Transformer.from_crs("EPSG:4326", f"EPSG:{_recipe_epsg(recipe)}", always_xy=True)
    x, y = fwd.transform([q[0] for q in line_ll], [q[1] for q in line_ll])
    line = shapely.LineString(list(zip(x, y)))
    o = mesh.open_boundaries[0]
    pts = shapely.points(mesh.nodes[o])
    pos = sorted(line.project(q) for q in (pts[0], pts[-1]))
    unc = float(line.length - (pos[1] - pos[0]))
    off = float(shapely.distance(pts, line).max())
    rec = {"uncovered_m": unc, "max_offset_m": off, "raw_uncovered_m": raw["uncovered_m"]}
    log(f"[obc] open line covered {pos[0]:.0f}-{pos[1]:.0f} of {line.length:.0f} m, "
        f"max offset {off:.1f} m")
    lim = float(cfg.get("max_move_m", 600.0))
    if unc > raw["uncovered_m"] + 1.0 or off > lim:
        raise SystemExit(f"open boundary no longer matches the line: {rec} (max_move_m {lim:g})")
    return rec
#<--- jsasaki 2026-11-16


def _stage_obc(recipe, out_dir, artifacts, log):
    from fvcom_mesh_tools.io import read_fort14, write_fort14
    from fvcom_mesh_tools.obc_tools import assign_west_south_obc

    cfg = recipe.get("obc", {}) or {}
    src = Path(artifacts.get("finished_mesh") or artifacts["raw_mesh"])
    mesh = read_fort14(src)
    # jsasaki 2026-11-16: region-supplied shoreline as a further fallback
    shoreline = artifacts.get("land_opened") \
        or artifacts.get("domain_shoreline") \
        or recipe["build"].get("coastline")
    # jsasaki 2026-11-16: recipe CRS, no 32654 fallback
    utm = int(cfg["utm_epsg"]) if cfg.get("utm_epsg") is not None else _recipe_epsg(recipe)
    mesh, info = assign_west_south_obc(
        mesh,
        utm_epsg=utm,
        protect_polygons=_protect_geometry(recipe, utm),
        band_deg=float(cfg.get("band_deg", 0.012)),
        shoreline_shp=shoreline,
        coast_tol_m=float(cfg.get("coast_tol_m", 500.0)),
        trim=int(cfg.get("trim", 1)),
        max_move_m=float(cfg.get("max_move_m", 600.0)),
        min_depth_m=cfg.get("min_depth_m"),
        obc_line_lonlat=_obc_line_from(recipe, out_dir),
        # jsasaki 2026-11-16: snap off when the build put the nodes on the line already
        snap=bool(cfg.get("snap", True)),
        log=log,
    )
    # jsasaki 2026-11-16: review newregion r3 #2: the open boundary must still cover the line
    _check_open_cover(recipe, out_dir, artifacts, mesh, cfg, log)
    out14 = out_dir / f"{recipe['name']}_obc.14"
    write_fort14(mesh, out14)
    return {"obc_mesh": str(out14), "obc_info": info}


def _stage_siteops(recipe, out_dir, artifacts, log):
    import json as _json

    from fvcom_mesh_tools.algorithms.boundary_snap import load_polylines
    from fvcom_mesh_tools.io import read_fort14, write_fort14
    from fvcom_mesh_tools.site_session import apply_site_operators

    cfg = recipe.get("siteops", {}) or {}
    src = Path(artifacts.get("obc_mesh")
               or artifacts.get("finished_mesh")
               or artifacts["raw_mesh"])
    mesh = read_fort14(src)
    # jsasaki 2026-11-16: recipe CRS and shoreline lookup (no 32654 fallback)
    utm = _recipe_epsg(recipe)
    shoreline = _shoreline(recipe, artifacts)
    lines = load_polylines(shoreline, to_crs=utm)
    min_edge = cfg.get("min_edge_m")
    if min_edge:
        from fvcom_mesh_tools.site_session import collapse_short_edges

        # jsasaki 2026-11-16: review newregion r2 #2: protected water guarded in the sweep
        mesh, n_short = collapse_short_edges(
            mesh, lines, min_len=float(min_edge),
            protect_polygons=_protect_geometry(recipe, utm), log=log,
        )
    # jsasaki 2026-11-16: protected water for the final component step
    mesh, edit_log = apply_site_operators(
        mesh, lines, passes=int(cfg.get("passes", 2)),
        protect_polygons=_protect_geometry(recipe, utm), log=log,
    )
    from fvcom_mesh_tools.mesh_clean import compact_nodes

    mesh, _ = compact_nodes(mesh)
    out14 = out_dir / f"{recipe['name']}_siteops.14"
    write_fort14(mesh, out14)
    (out_dir / "siteops_edit_log.json").write_text(
        _json.dumps(edit_log, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    ops = {}
    for r in edit_log:
        ops[r["op"]] = ops.get(r["op"], 0) + 1
    log(f"[pipeline] siteops: {ops}")
    return {"siteops_mesh": str(out14), "siteops_ops": ops}


def _stage_export(recipe, out_dir, artifacts, log):
    import numpy as np
    from pyproj import Transformer

    from fvcom_mesh_tools.io import read_fort14
    from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case

    cfg = recipe.get("export", {}) or {}
    src = Path(artifacts.get("final_mesh")
               or artifacts.get("obc_mesh")
               or artifacts["raw_mesh"])
    mesh = read_fort14(src)
    if not mesh.open_boundaries or not len(mesh.open_boundaries[0]):
        log("[pipeline] export SKIPPED: mesh has no open boundary")
        return {}
    # jsasaki 2026-11-16: recipe CRS (no 32654 fallback)
    utm = _recipe_epsg(recipe)
    tr = Transformer.from_crs(f"EPSG:{utm}", "EPSG:4326",
                              always_xy=True)
    _lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    obc = mesh.open_boundaries[0]
    sponge = [
        (int(v), float(cfg.get("sponge_radius_m", 3000.0)),
         float(cfg.get("sponge_coeff", 0.001)))
        for v in obc
    ]
    case_dir = Path(cfg.get("case_dir", out_dir / "fvcom_inputs"))
    written = export_fvcom_case(
        mesh, case_dir, cfg.get("casename", recipe["name"]),
        obc_type=int(cfg.get("obc_type", 1)),
        cor=lat, sponge=sponge,
    )
    for k, pth in written.items():
        log(f"[pipeline] export {k}: {pth}")
    del np
    return {"case_dir": str(case_dir)}


def _stage_polish(recipe, out_dir, artifacts, log):
    from fvcom_mesh_tools.finishing import finish_constrained_mesh
    from fvcom_mesh_tools.io import read_fort14, write_fort14

    cfg = recipe.get("polish", {}) or {}
    src = Path(artifacts.get("siteops_mesh")
               or artifacts["finished_mesh"])
    mesh = read_fort14(src)
    # jsasaki 2026-11-16: recipe CRS and shoreline lookup (no 32654 fallback)
    shoreline = _shoreline(recipe, artifacts)
    utm = _recipe_epsg(recipe)
    mesh, finfo = finish_constrained_mesh(
        mesh, shoreline, out_dir / "work",
        optimize_budget_s=float(cfg.get("optimize_budget_s", 300.0)),
        weld_tol_m=float(cfg.get("weld_tol_m", 2.0)),
        utm_epsg=utm,
        #---> jsasaki 2026-11-16: channel protection, valence, lone nodes (all off by default)
        protect_polygons=_protect_geometry(recipe, utm),
        on_dropped_protected=str(cfg.get("on_dropped_protected", "fail")),
        max_valence=cfg.get("max_valence"),
        valence_min_angle_deg=float(cfg.get("valence_min_angle_deg", 0.0)),
        remove_lone=bool(cfg.get("remove_lone_nodes", False)),
        freeze_obc_band=bool(cfg.get("freeze_obc_band", False)),
        # jsasaki 2026-11-16: review newregion r3 #2: protected-loss tolerance (default 0 m2)
        protected_loss_tol_m2=float(cfg.get("protected_loss_tol_m2", 0.0)),
        #<--- jsasaki 2026-11-16
        log=log,
    )
    out14 = out_dir / f"{recipe['name']}_polished.14"
    write_fort14(mesh, out14)
    return {"polish_mesh": str(out14), "polish_info": finfo}


def _stage_finishing(recipe, out_dir, artifacts, log):
    """Stage-2 automated finishing (docs/STAGE2_DESIGN.md):
    detect -> plan -> execute, single pass, per-patch atomic."""
    from fvcom_mesh_tools.autofinish import (
        detect_violations,
        execute_patches,
        plan_patches,
        write_ledger,
    )
    from fvcom_mesh_tools.io import read_fort14, write_fort14

    cfg = (recipe.get("finishing") or {}).get("auto") or {}
    src = Path(artifacts.get("final_mesh")
               or artifacts.get("obc_mesh")
               or artifacts["raw_mesh"])
    mesh = read_fort14(src)
    dcfg = (recipe.get("finishing") or {}).get("directives")
    dledger = []
    if dcfg:
        from fvcom_mesh_tools.autofinish import apply_directives
        from fvcom_mesh_tools.gridref import GridRef

        bbox = recipe["prep"]["bbox"]
        grid = GridRef(bbox[0], bbox[1], bbox[2], bbox[3])
        for d in dcfg:
            if "cell" in d and "polygon" not in d:
                d["polygon"] = grid.cell_polygon(d["cell"])
                log(f"[finishing] directive cell {d['cell']} -> "
                    "polygon")

        # jsasaki 2026-11-16: recipe CRS (no 32654 fallback)
        utm = _recipe_epsg(recipe)
        mesh, dledger = apply_directives(mesh, dcfg,
                                         utm_epsg=utm, log=log)
    obc = [int(v) for b in mesh.open_boundaries for v in b]
    det = detect_violations(mesh.nodes, mesh.elements)
    # jsasaki 2026-11-16: review newregion r3 #3: the mesh before the stage, for its water balance
    from fvcom_mesh_tools.finishing import _snapshot
    stage_in = _snapshot(mesh)
    #---> jsasaki 2026-11-16: obc_lock "none": patches near the open boundary are
    # worked too (open-boundary nodes never move); a patch that pushes an open-boundary
    # node over obc_perp_max_deg (or a violating one further) is reverted
    obc_lock = str(cfg.get("obc_lock", "kring"))
    if obc_lock not in ("kring", "none"):
        raise SystemExit(f"finishing.auto.obc_lock must be 'kring' or 'none', not {obc_lock!r}")
    patches = plan_patches(
        mesh.nodes, mesh.elements, det,
        obc_nodes=obc if obc_lock == "kring" else None,
        kring=int(cfg.get("patch_kring", 2)),
        max_patches=int(cfg.get("max_patches", 50)),
    )
    log(f"[finishing] {len(patches)} patches planned")
    ledger = []
    # jsasaki 2026-11-16: review newregion r3 #3: with protected water every patch is guarded
    protect = _protect_geometry(recipe, _recipe_epsg(recipe)) \
        if recipe.get("protect_polygons") else None
    if patches and (obc_lock == "none" or protect is not None):
        # jsasaki 2026-11-16: review newregion r1 #7: shared with the site fallback
        mesh, ledger = _execute_patches_guarded(
            mesh, patches, obc, float(cfg.get("obc_perp_max_deg", 20.0)), log,
            protect=protect, tol_m2=float(cfg.get("protected_loss_tol_m2", 0.0)))
    elif patches:
    #<--- jsasaki 2026-11-16
        mesh.nodes, ledger = execute_patches(
            mesh.nodes, mesh.elements, patches, obc_nodes=obc,
            log=log,
        )
    #---> jsasaki 2026-11-16: one bounded site-operator fallback for what the patches left
    if cfg.get("site_fallback"):
        mesh, fb = _finishing_site_fallback(recipe, mesh, artifacts, log)
        ledger = list(ledger) + [fb]
    #<--- jsasaki 2026-11-16
    #---> jsasaki 2026-11-16: review newregion r3 #3: protected water lost / gained by the stage
    if protect is not None:
        from fvcom_mesh_tools.finishing import protected_loss_m2

        lost = protected_loss_m2(stage_in, mesh, protect)
        gained = protected_loss_m2(mesh, stage_in, protect)
        log(f"[finishing] protected water: {lost:.1f} m2 lost, {gained:.1f} m2 gained")
        ledger = list(ledger) + [{"op": "protected_water", "outcome": "measured",
                                  "lost_m2": lost, "gained_m2": gained}]
    #<--- jsasaki 2026-11-16
    write_ledger(list(dledger) + list(ledger or patches),
                 out_dir / "finishing_ledger.json")
    out14 = Path(artifacts.get("final_mesh")
                 or out_dir / f"{recipe['name']}_finishing.14")
    write_fort14(mesh, out14)
    # jsasaki 2026-11-16: the site-fallback record is not a patch
    unfixed = [q for q in ledger
               if "outcome" in q and q.get("op") not in ("site_fallback", "protected_water")
               and not q["outcome"].startswith("fixed")]
    if unfixed:
        log(f"[finishing] {len(unfixed)} patches NOT fixed "
            "(see ledger)")
    return {"finishing_mesh": str(out14),
            "final_mesh": str(out14)}


#---> jsasaki 2026-11-16: site-operator fallback of the finishing stage
# jsasaki 2026-11-16: review newregion r1 #7: guarded patch execution, shared with the fallback
# jsasaki 2026-11-16: review newregion r3 #3: protected water guarded per patch as well
def _execute_patches_guarded(mesh, patches, obc, pmax, log, protect=None, tol_m2=0.0):
    """Execute autofinish patches one at a time; a patch that pushes an
    open-boundary node over ``pmax`` deg (or a violating one further), or
    that loses protected water (``protect``; geometric loss over ``tol_m2``,
    not netted against gains), is reverted. Returns ``(mesh, ledger)``."""
    import numpy as np

    from fvcom_mesh_tools.autofinish import execute_patches
    from fvcom_mesh_tools.finishing import _snapshot, obc_perpendicularity, protected_loss_m2

    ledger = []
    dev0 = obc_perpendicularity(mesh)
    for q in patches:
        nodes_try, led = execute_patches(mesh.nodes, mesh.elements, [q], obc_nodes=obc,
                                         log=log)
        before = _snapshot(mesh) if protect is not None else None
        keep_nodes = mesh.nodes
        mesh.nodes = nodes_try
        dev1 = obc_perpendicularity(mesh)
        why = None
        if (dev1 > np.maximum(dev0, pmax) + 1e-9).any():
            why = f"open-boundary perpendicularity over {pmax:g} deg"
        elif protect is not None:
            lost = protected_loss_m2(before, mesh, protect)
            if lost > tol_m2:
                why = f"{lost:.1f} m2 of protected water lost"
        if why:
            mesh.nodes = keep_nodes
            for r in led:
                if "outcome" in r and r["outcome"].startswith("fixed"):
                    r["outcome"] = f"reverted ({why})"
            log(f"[finishing] patch {q.get('id')} reverted: {why}")
        else:
            dev0 = dev1
        ledger.extend(led)
    return mesh, ledger



def _gate_violations(mesh) -> tuple[int, dict[str, int]]:
    """Total offender count over the failed QA gates (metric coordinates,
    no channel check) and the per-check counts."""
    from fvcom_mesh_tools.qa import run_qa

    rep = run_qa(mesh, coords="metric", channel_check=False, max_offenders=1)
    bad = {c.check_id: int(c.n_violations) for c in rep.checks
           if c.gate and not c.skipped and not c.passed}
    # jsasaki 2026-11-16: review newregion r1 #7: no finite penalty; new failures refuse instead
    return sum(bad.values()), bad


# gates the fallback may trade against each other; any other gate newly failing (or failing
# with more offenders) refuses the fallback (review newregion r1 #7)
_SITE_FALLBACK_GATES = ("c1_min_angle", "c2_max_angle", "c4_area_change")


def _finishing_site_fallback(recipe, mesh, artifacts, log):
    """``finishing.auto.site_fallback``: when C1/C2/C4 offenders remain after
    the patches, run ONE pass of the site operators (equalize / split /
    collapse; open-boundary nodes protected), then one more detect -> plan ->
    execute round, and keep the result only if the failed-gate offender count
    drops and no other gate starts failing. Single shot, no loop."""
    import copy

    from fvcom_mesh_tools.algorithms.boundary_snap import load_polylines
    from fvcom_mesh_tools.autofinish import (
        detect_violations,
        plan_patches,
    )
    from fvcom_mesh_tools.mesh_clean import compact_nodes
    from fvcom_mesh_tools.site_session import apply_site_operators

    s0, bad0 = _gate_violations(mesh)
    if not any(k in bad0 for k in _SITE_FALLBACK_GATES):
        return mesh, {"op": "site_fallback", "outcome": "not needed", "before": bad0}
    utm = _recipe_epsg(recipe)
    lines = load_polylines(_shoreline(recipe, artifacts), to_crs=utm)
    #---> jsasaki 2026-11-16: review newregion r3 #3/#4/#5: references taken BEFORE any edit
    from fvcom_mesh_tools.finishing import (
        _snapshot,
        obc_perpendicularity,
        protected_loss_m2,
        validate_open_boundary,
    )

    acfg = (recipe.get("finishing") or {}).get("auto", {}) or {}
    pmax = float(acfg.get("obc_perp_max_deg", 20.0))
    tol = float(acfg.get("protected_loss_tol_m2", 0.0))
    protect = _protect_geometry(recipe, utm)
    ref = _snapshot(mesh)
    dev_ref = obc_perpendicularity(mesh)
    #<--- jsasaki 2026-11-16
    trial = copy.deepcopy(mesh)
    # jsasaki 2026-11-16: review newregion r3 #4: the open boundary survives without protection too
    trial, elog = apply_site_operators(trial, lines, passes=1,
                                       protect_polygons=protect, keep_open_boundary=True,
                                       log=log)
    trial, _ = compact_nodes(trial)
    obc = [int(v) for b in trial.open_boundaries for v in b]
    det = detect_violations(trial.nodes, trial.elements)
    pts = plan_patches(trial.nodes, trial.elements, det, obc_nodes=None)
    # jsasaki 2026-11-16: review newregion r1 #7: fallback patches get the perpendicularity guard
    if pts:
        # jsasaki 2026-11-16: review newregion r3 #3: protected water guarded per patch
        trial, _led = _execute_patches_guarded(trial, pts, obc, pmax, lambda *a: None,
                                               protect=protect, tol_m2=tol)
    s1, bad1 = _gate_violations(trial)
    # jsasaki 2026-11-16: review newregion r1 #7: no gate outside C1/C2/C4 may start failing or get
    # worse; only then the offender totals are compared
    other_worse = [k for k in bad1 if k not in _SITE_FALLBACK_GATES
                   and bad1[k] > bad0.get(k, 0)]
    #---> jsasaki 2026-11-16: review newregion r3 #4/#5/#3: the WHOLE trial against the original:
    # open boundary intact (same nodes, order, boundary edges), every node's perpendicularity
    # within max(original, obc_perp_max_deg), no protected water lost
    import numpy as np

    obc_bad = validate_open_boundary(ref, trial, same_positions=True)
    if obc_bad:
        other_worse.append(f"open boundary {obc_bad}")
    elif len(dev_ref):
        dev1 = obc_perpendicularity(trial)
        if (dev1 > np.maximum(dev_ref, pmax) + 1e-9).any():
            other_worse.append(f"open-boundary perpendicularity over {pmax:g} deg")
    if protect is not None:
        lost = protected_loss_m2(ref, trial, protect)
        if lost > tol:
            other_worse.append(f"{lost:.1f} m2 of protected water lost")
    #<--- jsasaki 2026-11-16
    ok = not other_worse and s1 < s0
    # jsasaki 2026-11-16: review newregion r1 #7: name the gates that got worse
    log(f"[finishing] site fallback: {bad0} -> {bad1} "
        f"({'kept' if ok else 'refused'}{'; worse: ' + str(other_worse) if other_worse else ''})")
    rec = {"op": "site_fallback", "outcome": "kept" if ok else "refused",
           "before": bad0, "after": bad1,
           "ops": sorted({r.get("op") for r in elog if isinstance(r, dict)})}
    return (trial if ok else mesh), rec
#<--- jsasaki 2026-11-16


def _stage_obcfinal(recipe, out_dir, artifacts, log):
    from fvcom_mesh_tools.io import read_fort14, write_fort14
    from fvcom_mesh_tools.obc_tools import assign_west_south_obc

    cfg = recipe.get("obcfinal", recipe.get("obc", {})) or {}
    src = Path(artifacts.get("polish_mesh")
               or artifacts.get("siteops_mesh")
               or artifacts["obc_mesh"])
    mesh = read_fort14(src)
    # jsasaki 2026-11-16: region-supplied shoreline; recipe CRS (no 32654 fallback)
    shoreline = artifacts.get("land_opened") \
        or artifacts.get("domain_shoreline") \
        or recipe["build"].get("coastline")
    mesh, info = assign_west_south_obc(
        mesh,
        utm_epsg=_recipe_epsg(recipe),
        protect_polygons=_protect_geometry(recipe, _recipe_epsg(recipe)),
        shoreline_shp=shoreline,
        coast_tol_m=float(cfg.get("coast_tol_m", 500.0)),
        trim=int(cfg.get("trim", 1)),
        max_move_m=float(cfg.get("max_move_m", 600.0)),
        min_depth_m=cfg.get("min_depth_m"),
        perp_seed=9900,
        snap=False,
        obc_line_lonlat=_obc_line_from(recipe, out_dir),
        log=log,
    )
    eff = _obc_line_from(recipe, out_dir)
    proj_lines = artifacts.get("projector_lines")
    if eff and proj_lines and cfg.get("junction_realize", False):
        from pyproj import Transformer
        from shapely.geometry import LineString

        from fvcom_mesh_tools.algorithms.boundary_snap import (
            load_polylines,
        )
        from fvcom_mesh_tools.obc_tools import (
            realize_perpendicular_junctions,
        )

        # jsasaki 2026-11-16: recipe CRS (no 32654 fallback)
        utm = _recipe_epsg(recipe)
        trj = Transformer.from_crs("EPSG:4326", f"EPSG:{utm}",
                                   always_xy=True)
        ax, ay = trj.transform([q[0] for q in eff],
                               [q[1] for q in eff])
        arc = LineString(list(zip(ax, ay)))
        lines = load_polylines(Path(proj_lines), to_crs=utm)
        mesh, jinfo = realize_perpendicular_junctions(
            mesh, arc, lines,
            max_move_m=float(cfg.get(
                "junction_realize_max_move_m", 1800.0)),
            log=log,
        )
        info["junction_realized"] = jinfo
    # jsasaki 2026-11-16: review newregion r3 #2: the open boundary must still cover the line
    _check_open_cover(recipe, out_dir, artifacts, mesh, cfg, log)
    out14 = out_dir / f"{recipe['name']}_final.14"
    write_fort14(mesh, out14)
    return {"final_mesh": str(out14), "obcfinal_info": info}


#---> jsasaki 2026-11-16: pluggable depth products (new regions bring their own)
DEPTH_PRODUCTS = ("m7001_production", "external_dep_file", "none")
# jsasaki 2026-11-16: review newregion r1 #8: artifact keys of the meshes a stage writes
_MESH_KEYS = ("depths_mesh", "finishing_mesh", "final_mesh", "polish_mesh", "siteops_mesh",
              "obc_mesh", "finished_mesh", "raw_mesh")


def _read_dep_file(path):
    """Node depths from an FVCOM ``_dep.dat`` (``Node Number = N`` then
    ``x y h`` rows) or a fort.14: returns (xy[N, 2], h[N])."""
    import numpy as np

    lines = Path(path).read_text().splitlines()
    head = lines[0].lower()
    if "node" in head and "=" in head:
        n = int(lines[0].split("=")[1])
        a = np.array([ln.split()[:3] for ln in lines[1:1 + n]], float)
        return a[:, :2], a[:, 2]
    from fvcom_mesh_tools.io import read_fort14

    m = read_fort14(Path(path))
    return m.nodes, m.depths


def _stage_depths(recipe, out_dir, artifacts, log):
    """``depths.product``: ``none`` keeps the depths the mesh carries;
    ``external_dep_file`` reads ``depths.path`` (FVCOM _dep.dat or fort.14)
    whose nodes must be this mesh's nodes in order (within
    ``depths.match_tol_m``, default 0.01 m); ``m7001_production`` runs the
    Tokyo Bay M7001 recipe (dem.m7001.production_depths)."""
    import numpy as np
    from pyproj import Transformer

    from fvcom_mesh_tools.io import read_fort14, write_fort14

    cfg = recipe.get("depths") or {}
    product = str(cfg.get("product", "none"))
    if product not in DEPTH_PRODUCTS:
        raise SystemExit(f"depths.product must be one of {DEPTH_PRODUCTS}, not {product!r}")
    # jsasaki 2026-11-16: review newregion r1 #8: the mesh the LAST mesh stage wrote, and it is
    # propagated as final_mesh (also for product none)
    src = Path(artifacts.get("latest_mesh") or artifacts.get("final_mesh")
               or artifacts.get("obc_mesh") or artifacts["raw_mesh"])
    if product == "none":
        log(f"[depths] product none: {src.name} keeps its depths")
        return {"depths_product": "none", "final_mesh": str(src)}
    mesh = read_fort14(src)
    if product == "external_dep_file":
        xy, h = _read_dep_file(_path(cfg["path"]))
        tol = float(cfg.get("match_tol_m", 0.01))
        if len(xy) != mesh.n_nodes:
            raise SystemExit(f"depths file has {len(xy)} nodes, the mesh {mesh.n_nodes}: "
                             "it belongs to another mesh")
        dev = float(np.abs(xy - mesh.nodes).max())
        if dev > tol:
            raise SystemExit(f"depths file node coordinates differ from the mesh by {dev:.3g} m "
                             f"(> {tol:g}): it belongs to another mesh")
        mesh.depths = np.asarray(h, float)
    else:
        from fvcom_mesh_tools.dem.m7001 import production_depths

        tr = Transformer.from_crs(f"EPSG:{_recipe_epsg(recipe)}", "EPSG:4326",
                                  always_xy=True)
        lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        mesh.depths, rep_ = production_depths(np.asarray(lon), np.asarray(lat), mesh.elements)
        (out_dir / "depths_report.json").write_text(json.dumps(rep_, indent=1, default=str))
    out14 = out_dir / f"{recipe['name']}_depths.14"
    write_fort14(mesh, out14)
    log(f"[depths] {product}: {mesh.depths.min():.2f}..{mesh.depths.max():.2f} m -> {out14}")
    return {"depths_product": product, "depths_mesh": str(out14), "final_mesh": str(out14)}
#<--- jsasaki 2026-11-16


def _stage_qa(recipe, out_dir, artifacts, log):
    from fvcom_mesh_tools.io import read_fort14
    from fvcom_mesh_tools.qa import format_report, run_qa

    cfg = recipe.get("qa", {}) or {}
    target = Path(artifacts.get("final_mesh")
                  or artifacts.get("obc_mesh")
                  or artifacts["raw_mesh"])
    mesh = read_fort14(target)
    land_solid = out_dir / "prep" / "land_solid.shp"
    report = run_qa(mesh, name=target.name, path=target,
                    max_offenders=int(cfg.get("max_offenders", 1000)),
                    land_solid_shp=(land_solid if land_solid.exists()
                                    else None),
                    # jsasaki 2026-11-16: recipe CRS (no 32654 fallback)
                    utm_epsg=_recipe_epsg(recipe),
                    land_interior_m=float(
                        cfg.get("land_interior_m", 200.0)),
                    land_strict_polygon=(recipe.get("build", {})
                                         .get("nests", {})
                                         .get("inner", {})
                                         .get("polygon")))
    txt = format_report(report, lang=cfg.get("lang", "ja"))
    (out_dir / "qa_report.txt").write_text(txt, encoding="utf-8")
    (out_dir / "qa_report.json").write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    skip = set(cfg.get("skip_checks") or [])
    fails = [c.check_id for c in report.checks
             if c.gate and not c.skipped and not c.passed]
    counted = [c for c in fails if c not in skip]
    ignored = [c for c in fails if c in skip]
    log(f"[pipeline] QA: {len(counted)} gate failures"
        + (f" (+{len(ignored)} out-of-scope: {ignored})"
           if ignored else "")
        + f" (report -> {out_dir / 'qa_report.txt'})")
    return {"qa_failures": len(counted),
            "qa_out_of_scope": ignored}


def _stage_figures(recipe, out_dir, artifacts, log):
    import numpy as np
    from pyproj import Transformer

    from fvcom_mesh_tools.io import read_fort14
    from fvcom_mesh_tools.plotting import plot_mesh_overview

    cfg = recipe.get("figures", {}) or {}
    views = cfg.get("views", {"full": None})
    coast = cfg.get("coast_bbox")
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for key in ("raw_mesh", "finished_mesh", "obc_mesh",
                "siteops_mesh", "polish_mesh", "final_mesh"):
        path = artifacts.get(key)
        if not path or not Path(path).exists():
            continue
        mesh = read_fort14(Path(path))
        # jsasaki 2026-11-16: review newregion r2 #7: the raw mesh is transformed FROM its recorded
        # CRS and not at all when that is EPSG:4326; later meshes are in the recipe CRS when metric
        src_crs = (artifacts.get("raw_mesh_crs") if key == "raw_mesh" else None) \
            or (f"EPSG:{_recipe_epsg(recipe)}" if np.abs(mesh.nodes[:, 0]).max() > 1000.0
                else "EPSG:4326")
        if src_crs != "EPSG:4326":
            # jsasaki 2026-11-16: recipe CRS (was a literal EPSG:32654)
            tr = Transformer.from_crs(src_crs, "EPSG:4326",
                                      always_xy=True)
            lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
            mesh.nodes = np.column_stack([lon, lat])
        for view, zb in views.items():
            out = fig_dir / f"{Path(path).stem}_{view}.png"
            plot_mesh_overview(
                mesh, out, crs="EPSG:4326", cell_m=None,
                coast=tuple(coast) if coast else None,
                zoom=tuple(zb) if zb else None,
                dpi=int(cfg.get("dpi", 220)),
                title=f"{Path(path).stem} - {view}",
            )
            written.append(str(out))
    log(f"[pipeline] figures: {len(written)} -> {fig_dir}")
    return {"figures": written}


STAGES = [
    ("prep", _stage_prep),
    # jsasaki 2026-11-16: region-supplied domain (water polygon + open-boundary line)
    ("domain", _stage_domain),
    ("build", _stage_build),
    ("finish", _stage_finish),
    ("obc", _stage_obc),
    ("siteops", _stage_siteops),
    ("polish", _stage_polish),
    ("obcfinal", _stage_obcfinal),
    ("finishing", _stage_finishing),
    # jsasaki 2026-11-16: pluggable depth products
    ("depths", _stage_depths),
    ("qa", _stage_qa),
    ("export", _stage_export),
    ("figures", _stage_figures),
]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="fmesh-pipeline",
                                description=__doc__)
    p.add_argument("recipe", type=Path, help="Recipe YAML file.")
    p.add_argument(
        "--only", type=str, default=None,
        help="Comma-separated subset of stages to run "
             "(artifacts of earlier stages are read from out_dir).",
    )
    # jsasaki 2026-11-16: a region repository may keep the pipeline recipe as one block of its own
    p.add_argument(
        "--section", type=str, default=None,
        help="Read the pipeline recipe from this top-level key of the YAML "
             "(e.g. 'pipeline' in a region repository's mesh recipe).",
    )
    args = p.parse_args(argv)

    # jsasaki 2026-11-16: --section and the recipe rule check
    recipe = _load_recipe(args.recipe.resolve(), section=args.section)
    _check_recipe(recipe)
    out_dir = Path(recipe["out_dir"]).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    only = set(args.only.split(",")) if args.only else None

    prov_path = out_dir / "pipeline_provenance.json"
    artifacts: dict[str, Any] = {}
    if prov_path.exists():
        artifacts = json.loads(
            prov_path.read_text()
        ).get("artifacts", {})

    def log(msg: str) -> None:
        print(msg, flush=True)

    results: dict[str, Any] = {}
    t0 = time.perf_counter()
    for name, fn in STAGES:
        if name not in recipe:
            continue
        if only is not None and name not in only:
            continue
        log(f"[pipeline] === stage {name} ===")
        t1 = time.perf_counter()
        out = fn(recipe, out_dir, artifacts, log) if name != "prep" \
            else fn(recipe, out_dir, log)
        artifacts.update(out or {})
        # jsasaki 2026-11-16: review newregion r1 #8: remember the mesh the last mesh stage wrote
        for key in _MESH_KEYS:
            if (out or {}).get(key):
                artifacts["latest_mesh"] = out[key]
                break
        results[name] = {"wall_s": round(time.perf_counter() - t1, 1),
                         **{k: v for k, v in (out or {}).items()
                            if not isinstance(v, list)}}

    prov_path.write_text(json.dumps({
        "recipe_file": str(args.recipe.resolve()),
        "recipe": recipe,
        "git": _git_describe(Path(__file__).resolve().parents[3]),
        "artifacts": {k: v for k, v in artifacts.items()
                      if isinstance(v, str)},
        "results": results,
        "total_wall_s": round(time.perf_counter() - t0, 1),
    }, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    log(f"[pipeline] provenance -> {prov_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
