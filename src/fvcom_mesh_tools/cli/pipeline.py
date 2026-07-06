"""fmesh-pipeline: recipe-driven end-to-end mesh construction.

Replaces the PoC-era sed-cloned notebook scripts with a single
declarative entry point:

    fmesh-pipeline recipes/tokyo_bay_v5.yaml [--only STAGE[,STAGE]]

Stages run in fixed order, each at most ONCE (no convergence loops),
each optional by presence of its key in the recipe:

    prep     -> fmesh-prep-shoreline (land opening + skeleton seeds)
    build    -> fmesh-buildmesh (oceanmesh engine, CDT boundary)
    finish   -> finishing.finish_constrained_mesh (UTM)
    qa       -> fvcom_mesh_tools.qa.run_qa report (ja/en)
    figures  -> per-stage overview/zoom PNGs

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


def _load_recipe(path: Path) -> dict[str, Any]:
    import yaml

    with path.open() as f:
        recipe = yaml.safe_load(f)
    if not isinstance(recipe, dict):
        raise SystemExit(f"recipe {path} is not a mapping")
    for key in ("name", "out_dir"):
        if key not in recipe:
            raise SystemExit(f"recipe missing required key: {key!r}")
    return recipe


def _git_describe(repo: Path) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(repo), "describe", "--always", "--dirty"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:  # noqa: BLE001 - provenance is best-effort
        return "unknown"


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
    rc = prepshoreline.main(argv)
    if rc:
        raise SystemExit(f"prep stage failed (rc={rc})")
    out = {
        "land_opened": str(prep_dir / "land_opened.shp"),
        "skeleton_seeds": str(prep_dir / "skeleton_seeds.shp"),
    }
    r_sm = cfg.get("outer_smooth_r_m")
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
        outer_gdf = smooth_land_per_polygon(land, r_m=float(r_sm))
        outer_pre_cut = outer_gdf.copy()
        eff = None
        if cfg.get("obc_line"):
            eff = extend_obc_ends_perpendicular(
                [tuple(q) for q in cfg["obc_line"]], outer_gdf,
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

            jc = junction_tail_constraints(
                eff, outer_pre_cut,
                tail_len_m=float(cfg.get("junction_tail_m", 1500.0)),
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

            solid = gpd.read_file(solid_p).to_crs(32654)
            opened_u = unary_union(list(
                land.to_crs(32654).geometry))
            lost = unary_union(list(solid.geometry)).difference(
                opened_u).buffer(-200.0)
            grid = GridRef(*cfg["bbox"])
            import geopandas as _gpd2
            for geom in getattr(lost, "geoms", [lost]):
                if geom.is_empty or geom.area < 2.5e5:
                    continue
                cen = _gpd2.GeoSeries([geom.centroid], crs=32654
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


def _stage_build(recipe, out_dir, artifacts, log):
    import os

    from fvcom_mesh_tools.cli import buildmesh

    cfg = recipe["build"]
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
        return {"raw_mesh": str(raw14)}
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
    return {"raw_mesh": str(raw14)}


def _stage_finish(recipe, out_dir, artifacts, log):
    import numpy as np
    from pyproj import Transformer

    from fvcom_mesh_tools.finishing import finish_constrained_mesh
    from fvcom_mesh_tools.io import read_fort14, write_fort14

    cfg = recipe.get("finish", {}) or {}
    utm = int(cfg.get("utm_epsg", 32654))
    mesh = read_fort14(Path(artifacts["raw_mesh"]))
    tr = Transformer.from_crs("EPSG:4326", f"EPSG:{utm}", always_xy=True)
    x, y = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
    mesh.nodes = np.column_stack([x, y])
    shoreline = Path(artifacts.get("projector_lines")
                     or artifacts.get("land_opened")
                     or recipe["build"]["coastline"])
    mesh, finfo = finish_constrained_mesh(
        mesh, shoreline, out_dir / "work",
        optimize_budget_s=float(cfg.get("optimize_budget_s", 600.0)),
        weld_tol_m=float(cfg.get("weld_tol_m", 2.0)),
        utm_epsg=utm,
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
    return (recipe.get("prep", {}) or {}).get("obc_line")


def _stage_obc(recipe, out_dir, artifacts, log):
    from fvcom_mesh_tools.io import read_fort14, write_fort14
    from fvcom_mesh_tools.obc_tools import assign_west_south_obc

    cfg = recipe.get("obc", {}) or {}
    src = Path(artifacts.get("finished_mesh") or artifacts["raw_mesh"])
    mesh = read_fort14(src)
    shoreline = artifacts.get("land_opened") \
        or recipe["build"].get("coastline")
    mesh, info = assign_west_south_obc(
        mesh,
        utm_epsg=int(cfg.get("utm_epsg",
                             (recipe.get("finish") or {})
                             .get("utm_epsg", 32654))),
        band_deg=float(cfg.get("band_deg", 0.012)),
        shoreline_shp=shoreline,
        coast_tol_m=float(cfg.get("coast_tol_m", 500.0)),
        trim=int(cfg.get("trim", 1)),
        max_move_m=float(cfg.get("max_move_m", 600.0)),
        min_depth_m=cfg.get("min_depth_m"),
        obc_line_lonlat=_obc_line_from(recipe, out_dir),
        log=log,
    )
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
    utm = int((recipe.get("finish") or {}).get("utm_epsg", 32654))
    shoreline = Path(artifacts.get("projector_lines")
                     or artifacts.get("land_opened")
                     or recipe["build"]["coastline"])
    lines = load_polylines(shoreline, to_crs=utm)
    min_edge = cfg.get("min_edge_m")
    if min_edge:
        from fvcom_mesh_tools.site_session import collapse_short_edges

        mesh, n_short = collapse_short_edges(
            mesh, lines, min_len=float(min_edge), log=log,
        )
    mesh, edit_log = apply_site_operators(
        mesh, lines, passes=int(cfg.get("passes", 2)), log=log,
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
    utm = int((recipe.get("finish") or {}).get("utm_epsg", 32654))
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
    shoreline = Path(artifacts.get("projector_lines")
                     or artifacts.get("land_opened")
                     or recipe["build"]["coastline"])
    utm = int((recipe.get("finish") or {}).get("utm_epsg", 32654))
    mesh, finfo = finish_constrained_mesh(
        mesh, shoreline, out_dir / "work",
        optimize_budget_s=float(cfg.get("optimize_budget_s", 300.0)),
        weld_tol_m=float(cfg.get("weld_tol_m", 2.0)),
        utm_epsg=utm,
        log=log,
    )
    out14 = out_dir / f"{recipe['name']}_polished.14"
    write_fort14(mesh, out14)
    return {"polish_mesh": str(out14), "polish_info": finfo}


def _stage_finishing(recipe, out_dir, artifacts, log):
    """Stage-2 automated finishing (docs/STAGE2_DESIGN.md):
    detect -> plan -> execute, single pass, per-patch atomic."""
    from fvcom_mesh_tools.finishing import (
        detect_violations,
        execute_patches,
        plan_patches,
        write_ledger,
    )
    from fvcom_mesh_tools.io import read_fort14, write_fort14

    cfg = (recipe.get("finishing") or {}).get("auto") or {}
    src = Path(artifacts.get("obc_mesh")
               or artifacts.get("finished_mesh")
               or artifacts["raw_mesh"])
    mesh = read_fort14(src)
    dcfg = (recipe.get("finishing") or {}).get("directives")
    dledger = []
    if dcfg:
        from fvcom_mesh_tools.finishing import apply_directives
        from fvcom_mesh_tools.gridref import GridRef

        bbox = recipe["prep"]["bbox"]
        grid = GridRef(bbox[0], bbox[1], bbox[2], bbox[3])
        for d in dcfg:
            if "cell" in d and "polygon" not in d:
                d["polygon"] = grid.cell_polygon(d["cell"])
                log(f"[finishing] directive cell {d['cell']} -> "
                    "polygon")

        utm = int((recipe.get("finish") or {})
                  .get("utm_epsg", 32654))
        mesh, dledger = apply_directives(mesh, dcfg,
                                         utm_epsg=utm, log=log)
    obc = [int(v) for b in mesh.open_boundaries for v in b]
    det = detect_violations(mesh.nodes, mesh.elements)
    patches = plan_patches(
        mesh.nodes, mesh.elements, det, obc_nodes=obc,
        kring=int(cfg.get("patch_kring", 2)),
        max_patches=int(cfg.get("max_patches", 50)),
    )
    log(f"[finishing] {len(patches)} patches planned")
    ledger = []
    if patches:
        mesh.nodes, ledger = execute_patches(
            mesh.nodes, mesh.elements, patches, obc_nodes=obc,
            log=log,
        )
    write_ledger(list(dledger) + list(ledger or patches),
                 out_dir / "finishing_ledger.json")
    out14 = out_dir / f"{recipe['name']}_finishing.14"
    write_fort14(mesh, out14)
    unfixed = [q for q in ledger
               if "outcome" in q
               and not q["outcome"].startswith("fixed")]
    if unfixed:
        log(f"[finishing] {len(unfixed)} patches NOT fixed "
            "(see ledger)")
    return {"finishing_mesh": str(out14)}


def _stage_obcfinal(recipe, out_dir, artifacts, log):
    from fvcom_mesh_tools.io import read_fort14, write_fort14
    from fvcom_mesh_tools.obc_tools import assign_west_south_obc

    cfg = recipe.get("obcfinal", recipe.get("obc", {})) or {}
    src = Path(artifacts.get("finishing_mesh")
               or artifacts.get("polish_mesh")
               or artifacts.get("siteops_mesh")
               or artifacts["obc_mesh"])
    mesh = read_fort14(src)
    shoreline = artifacts.get("land_opened") \
        or recipe["build"].get("coastline")
    mesh, info = assign_west_south_obc(
        mesh,
        utm_epsg=int((recipe.get("finish") or {})
                     .get("utm_epsg", 32654)),
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

        utm = int((recipe.get("finish") or {})
                  .get("utm_epsg", 32654))
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
    out14 = out_dir / f"{recipe['name']}_final.14"
    write_fort14(mesh, out14)
    return {"final_mesh": str(out14), "obcfinal_info": info}


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
                    utm_epsg=int((recipe.get("finish") or {})
                                 .get("utm_epsg", 32654)))
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
        if np.abs(mesh.nodes[:, 0]).max() > 1000.0:
            tr = Transformer.from_crs("EPSG:32654", "EPSG:4326",
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
    ("build", _stage_build),
    ("finish", _stage_finish),
    ("obc", _stage_obc),
    ("siteops", _stage_siteops),
    ("polish", _stage_polish),
    ("obcfinal", _stage_obcfinal),
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
    args = p.parse_args(argv)

    recipe = _load_recipe(args.recipe.resolve())
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
