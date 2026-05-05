"""``fmesh-buildmesh`` CLI: DEM -> OCSMesh -> classified, perp-fixed fort.14.

Pipeline (single command, no MATLAB intermediate):

    DEM (NetCDF/GeoTIFF) -> ocsmesh.Geom(zmax)
                          -> ocsmesh.Hfun(hmin, hmax)
                          -> MeshDriver(engine="gmsh").run()
                          -> mesh.interpolate(raster) for depth
                          -> classify boundaries by DEM bbox proximity
                          -> (optional) edge-swap + Laplacian quality pass
                          -> (optional) align_open_boundary_first_ring
                          -> fort.14 (FVCOM/ADCIRC convention, depth +down)

The intent is "give us a fort.14 the FVCOM harness can actually load".
Quality post-processing is best-effort: it monotonically improves mean
alpha-quality but plateaus at a level that depends on the initial size
function (PoC #10). Driving the bad-element fraction to zero requires
adaptive sizing, which ``docs/python_pipeline_gap_analysis.md`` tracks.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms import (
    align_open_boundary_first_ring,
    alpha_quality,
    classify_boundaries_by_bbox,
    laplacian_smooth,
    min_interior_angle,
    signed_areas,
    swap_edges_for_quality,
)
from fvcom_mesh_tools.io import (
    Fort14Mesh,
    filter_multipolygon_by_area,
    load_coastline_as_lines,
    write_fort14,
)

EARTH_R_M = 6_371_000.0


def _deg_per_metre(lat_deg: float) -> float:
    """Conservative degrees-per-metre at ``lat_deg`` (longitude direction)."""
    return 1.0 / (EARTH_R_M * np.cos(np.deg2rad(lat_deg)) * np.pi / 180.0)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-buildmesh",
        description=(
            "Generate a FVCOM-ready fort.14 from a single DEM. "
            "Mesh is produced by OCSMesh + gmsh; depths are interpolated "
            "from the same DEM; the outer ring is split into open and "
            "land segments by proximity to the DEM bounding box; an "
            "optional first-ring perpendicularity fix is applied."
        ),
    )
    p.add_argument("dem", type=Path, help="DEM raster (GeoTIFF / NetCDF / etc).")
    p.add_argument("output", type=Path, help="Output fort.14 path.")
    p.add_argument(
        "--hmin", type=float, default=200.0,
        help="Minimum element size in metres (default: 200).",
    )
    p.add_argument(
        "--hmax", type=float, default=5000.0,
        help="Maximum element size in metres (default: 5000).",
    )
    p.add_argument(
        "--zmax", type=float, default=0.0,
        help="Geom: water domain is the region where DEM <= zmax (default: 0.0).",
    )
    p.add_argument(
        "--engine", choices=["gmsh"], default="gmsh",
        help="OCSMesh engine (default: gmsh).",
    )
    p.add_argument(
        "--interp-method", choices=["spline", "linear", "nearest"], default="linear",
        help="Depth-interpolation method (default: linear).",
    )
    p.add_argument(
        "--bbox-tol-m", type=float, default=None,
        help=(
            "Distance tolerance for 'on the DEM bbox' open-boundary "
            "classification, in metres. Default: 0.75 * hmin."
        ),
    )
    p.add_argument(
        "--land-ibtype", type=int, default=20,
        help=(
            "ibtype to write for every land segment in fort.14 "
            "(default: 20, matching the Tokyo Bay reference)."
        ),
    )
    p.add_argument(
        "--open-merge-coast-gap", type=int, default=0, metavar="NODES",
        help=(
            "Bridge short coast intrusions between two open arcs into "
            "one open segment when the intrusion is shorter than NODES "
            "(default 0 = no merging). Useful when the DEM bbox is "
            "rectangular but the coastline pokes into the bbox in a "
            "few places, splitting one geometric open arc into many."
        ),
    )
    p.add_argument(
        "--min-polygon-area-m2", type=float, default=0.0, metavar="M2",
        help=(
            "Drop wet-domain polygons whose metric area is below this "
            "threshold before meshing. Useful for stripping isolated "
            "single-pixel water bodies. 0 = keep everything (default)."
        ),
    )
    p.add_argument(
        "--min-island-area-m2", type=float, default=0.0, metavar="M2",
        help=(
            "Drop holes (islands) inside the wet-domain polygons whose "
            "metric area is below this threshold. 0 = keep everything "
            "(default). Typical: 10000-100000 m^2 to keep islands "
            "smaller than a few mesh cells."
        ),
    )
    p.add_argument(
        "--coastline", type=Path, action="append", default=[], metavar="PATH",
        help=(
            "Shapefile / GeoJSON of coastline polylines (or polygons; "
            "rings will be flattened) to drive Hfun.add_feature for "
            "coastline-aware sizing. Pass multiple times to combine "
            "sources. Inputs are reprojected to EPSG:4326 and clipped "
            "to the DEM bbox before use."
        ),
    )
    p.add_argument(
        "--coast-target-size", type=float, default=None, metavar="METRES",
        help=(
            "Target element size on the coastline polylines (default: "
            "--hmin). Sizing expands away from the coast at "
            "--coast-expansion-rate."
        ),
    )
    p.add_argument(
        "--coast-expansion-rate", type=float, default=0.005, metavar="RATE",
        help=(
            "Expansion rate passed to OCSMesh's Hfun.add_feature. "
            "Larger -> sizing relaxes faster away from the coast "
            "(default: 0.005)."
        ),
    )
    p.add_argument(
        "--quality-pass", type=int, default=0, metavar="ROUNDS",
        help=(
            "Run ROUNDS alternating edge-swap + Laplacian-smooth passes "
            "before perpfix. 0 disables (default). 6 is a good upper "
            "bound (PoC #10); returns plateau quickly."
        ),
    )
    p.add_argument(
        "--smooth-iters", type=int, default=5,
        help="Smoothing iterations per --quality-pass round (default: 5).",
    )
    p.add_argument(
        "--smooth-alpha", type=float, default=0.5,
        help="Damping factor for the per-round smoothing pass (default: 0.5).",
    )
    p.add_argument(
        "--no-perpfix", action="store_true",
        help="Skip the open-boundary first-ring perpendicularity correction.",
    )
    p.add_argument(
        "--perpfix-iters", type=int, default=1,
        help="Iterations of the perpendicularity fix (default: 1).",
    )
    p.add_argument(
        "--title", type=str, default=None,
        help="Title to write on the first line of fort.14 (default: derived from DEM).",
    )
    p.add_argument(
        "--quiet", action="store_true",
        help="Suppress the progress / summary output.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.dem.exists():
        print(f"DEM not found: {args.dem}", file=sys.stderr)
        return 2
    if args.hmin <= 0 or args.hmax < args.hmin:
        print("--hmin must be > 0 and --hmax must be >= --hmin.", file=sys.stderr)
        return 2

    args.output.parent.mkdir(parents=True, exist_ok=True)

    log = (lambda *a, **k: None) if args.quiet else print  # noqa: E731

    # OCSMesh / pyproj imports are deferred so the rest of the package
    # stays importable in environments without OCSMesh / gmsh installed.
    from ocsmesh import Geom, Hfun, MeshDriver, Raster
    from pyproj import CRS, Transformer

    log(f"[buildmesh] DEM: {args.dem}")
    log(f"[buildmesh] hmin={args.hmin:g} m  hmax={args.hmax:g} m  zmax={args.zmax:g}")

    raster = Raster(str(args.dem))
    xmin, ymin, xmax, ymax = raster.bbox.bounds
    log(f"[buildmesh] DEM bbox: x[{xmin:.6f}, {xmax:.6f}]  y[{ymin:.6f}, {ymax:.6f}]")

    t0 = time.perf_counter()
    geom = Geom(raster, zmax=args.zmax)

    if args.min_polygon_area_m2 > 0 or args.min_island_area_m2 > 0:
        log(
            f"[buildmesh] geom filter: min_polygon={args.min_polygon_area_m2:g} m^2, "
            f"min_island={args.min_island_area_m2:g} m^2"
        )
        t_filter = time.perf_counter()
        mp = geom.get_multipolygon()
        n_polys_before = len(mp.geoms)
        n_holes_before = sum(len(list(p.interiors)) for p in mp.geoms)
        mp_filt = filter_multipolygon_by_area(
            mp,
            src_crs=raster.crs,
            min_polygon_area_m2=args.min_polygon_area_m2,
            min_island_area_m2=args.min_island_area_m2,
        )
        n_polys_after = len(mp_filt.geoms)
        n_holes_after = sum(len(list(p.interiors)) for p in mp_filt.geoms)
        log(
            f"[buildmesh] geom filter: polygons {n_polys_before} -> {n_polys_after}, "
            f"holes {n_holes_before} -> {n_holes_after} "
            f"({time.perf_counter() - t_filter:.2f} s)"
        )
        if n_polys_after == 0:
            print(
                "Geom filter dropped every polygon; "
                "loosen --min-polygon-area-m2.",
                file=sys.stderr,
            )
            return 4
        # Replace the raster geom with a polygon-backed one driven by
        # the filtered multipolygon. Mesher uses this for the wet
        # domain; Hfun keeps the raster for sizing. OCSMesh requires
        # an explicit CRS for the polygon variant.
        geom = Geom(mp_filt, crs=raster.crs)

    hfun = Hfun(raster, hmin=args.hmin, hmax=args.hmax)

    if args.coastline:
        coast = load_coastline_as_lines(
            args.coastline, bbox=(xmin, ymin, xmax, ymax),
        )
        n_lines = len(coast.geoms)
        if n_lines == 0:
            log(
                "[buildmesh] coastline: no features inside DEM bbox; "
                "skipping add_feature."
            )
        else:
            target = args.coast_target_size or args.hmin
            log(
                f"[buildmesh] coastline: {n_lines} line strings; "
                f"target_size={target:g} m  expansion_rate={args.coast_expansion_rate:g}"
            )
            t_feat = time.perf_counter()
            hfun.add_feature(
                feature=coast,
                expansion_rate=args.coast_expansion_rate,
                target_size=target,
            )
            log(f"[buildmesh] coastline: add_feature took {time.perf_counter() - t_feat:.2f} s")

    driver = MeshDriver(geom, hfun=hfun, engine_name=args.engine)
    mesh = driver.run()
    t_gen = time.perf_counter() - t0
    log(f"[buildmesh] mesh generation ({args.engine}): {t_gen:.2f} s")

    log(f"[buildmesh] interpolating depths (method={args.interp_method}) ...")
    t1 = time.perf_counter()
    mesh.interpolate(raster, method=args.interp_method)
    log(f"[buildmesh] depth interpolation: {time.perf_counter() - t1:.2f} s")

    coords = np.asarray(mesh.coord, dtype=np.float64)
    elements = np.asarray(mesh.triangles, dtype=np.int64)
    if elements.size == 0:
        print("OCSMesh produced zero triangles; aborting.", file=sys.stderr)
        return 3
    raw_values = np.asarray(mesh.value, dtype=np.float64).ravel()
    if raw_values.size != coords.shape[0]:
        # Defensive guard: degenerate interpolations have surfaced empty
        # value arrays in past OCSMesh releases.
        depths = np.zeros(coords.shape[0], dtype=np.float64)
    else:
        # OCSMesh's GRD writer writes -value, which corresponds to ADCIRC's
        # "positive depth = below MSL" convention when DEM elevation is
        # +up. Replicate that here so our writer sees the same sign.
        depths = -raw_values

    # OCSMesh internally meshes in a metric CRS (UTM by default) and
    # exposes coords there, even when the input DEM is geographic.
    # Project back to EPSG:4326 so our fort.14 matches the legacy
    # reference's lon/lat convention.
    src_crs = mesh.crs
    if src_crs is not None and not CRS(src_crs).equals(CRS.from_epsg(4326)):
        log(f"[buildmesh] projecting coords {src_crs} -> EPSG:4326")
        transformer = Transformer.from_crs(
            src_crs, CRS.from_epsg(4326), always_xy=True,
        )
        lon, lat = transformer.transform(coords[:, 0], coords[:, 1])
        coords = np.column_stack([lon, lat])

    title = args.title or f"fmesh-buildmesh {args.dem.name}"
    f14 = Fort14Mesh(
        title=title,
        nodes=coords[:, :2].copy(),
        depths=depths,
        elements=elements.copy(),
        open_boundaries=[],
        land_boundaries=[],
    )

    # OCSMesh / gmsh emit triangles in clockwise order; the ADCIRC
    # fort.14 convention is counter-clockwise. Detect and flip.
    sa = signed_areas(f14)
    if (sa < 0).mean() > 0.5:
        log("[buildmesh] flipping triangle winding (CW -> CCW)")
        f14.elements = f14.elements[:, [0, 2, 1]].copy()

    bbox_tol_m = args.bbox_tol_m if args.bbox_tol_m is not None else 0.75 * args.hmin
    lat_mid = 0.5 * (ymin + ymax)
    tol_deg = bbox_tol_m * _deg_per_metre(lat_mid)
    log(f"[buildmesh] bbox classify tol: {bbox_tol_m:g} m ({tol_deg:.2e} deg)")

    open_segs, land_bnds = classify_boundaries_by_bbox(
        f14,
        bbox=(xmin, ymin, xmax, ymax),
        tol=tol_deg,
        land_ibtype=args.land_ibtype,
        open_merge_coast_gap=args.open_merge_coast_gap,
    )
    f14.open_boundaries = open_segs
    f14.land_boundaries = land_bnds
    log(
        f"[buildmesh] classified: open={len(open_segs)} segments "
        f"({sum(s.size for s in open_segs)} nodes), "
        f"land={len(land_bnds)} segments "
        f"({sum(s.size for _, s in land_bnds)} nodes)"
    )

    if args.quality_pass > 0:
        log(
            f"[buildmesh] quality pass: rounds={args.quality_pass}, "
            f"smooth_iters={args.smooth_iters}, smooth_alpha={args.smooth_alpha}"
        )
        q_before = float(alpha_quality(f14).mean())
        bad_before = float((min_interior_angle(f14) < 20).mean()) * 100
        for r in range(args.quality_pass):
            f14, swap_info = swap_edges_for_quality(f14, max_iters=10)
            f14, smooth_info = laplacian_smooth(
                f14,
                n_iters=args.smooth_iters,
                alpha=args.smooth_alpha,
                prevent_flips=True,
            )
            log(
                f"[buildmesh]   round {r + 1}: swaps={swap_info['total_swaps']:,}  "
                f"smooth_reverts={int(sum(smooth_info['reverts'])):,}"
            )
        q_after = float(alpha_quality(f14).mean())
        bad_after = float((min_interior_angle(f14) < 20).mean()) * 100
        log(
            f"[buildmesh] quality pass: alpha {q_before:.4f} -> {q_after:.4f}, "
            f"frac<20deg {bad_before:.2f}% -> {bad_after:.2f}%"
        )

    if not args.no_perpfix and open_segs:
        log(f"[buildmesh] perpfix: aligning first-ring (iters={args.perpfix_iters})")
        nodes_before = f14.nodes.copy()
        f14, info = align_open_boundary_first_ring(
            f14,
            alpha=1.0,
            n_iters=args.perpfix_iters,
            smooth_iters=0,
            segment_index=0,
        )
        log(f"[buildmesh] perpfix moved {info['moved']:,} interior nodes")

        # The length-preserving perpendicular projection may flip
        # narrow triangles. Iterate the revert: each round undoes any
        # moved node that participates in a flipped triangle, then re-
        # checks signed areas. After the quality pass the surrounding
        # geometry is tighter, so a single round can be insufficient.
        n_revert_total = 0
        for _ in range(5):
            bad_tris = signed_areas(f14) <= 0
            if not bool(bad_tris.any()):
                break
            moved_mask = np.any(f14.nodes != nodes_before, axis=1)
            bad_nodes = np.zeros(f14.n_nodes, dtype=bool)
            bad_nodes[np.unique(f14.elements[bad_tris].ravel())] = True
            revert = bad_nodes & moved_mask
            if not bool(revert.any()):
                break
            f14.nodes[revert] = nodes_before[revert]
            n_revert_total += int(revert.sum())
        if n_revert_total:
            log(f"[buildmesh] perpfix: reverted {n_revert_total} moves to avoid flips")
        still_bad = int((signed_areas(f14) <= 0).sum())
        if still_bad:
            # Safety net: if iterative revert cannot clear the flips,
            # roll back perpfix entirely. The mesh is then exactly as
            # the quality pass left it - which is at least valid.
            log(
                f"[buildmesh] WARN: {still_bad} flipped triangles remain; "
                f"undoing perpfix to keep the mesh valid."
            )
            f14.nodes = nodes_before

    write_fort14(f14, args.output)
    log(f"[buildmesh] wrote {args.output}")
    log(
        f"[buildmesh] NP={f14.n_nodes:,}  NE={f14.n_elements:,}  "
        f"depth p50={float(np.nanpercentile(f14.depths, 50)):.3f} m  "
        f"depth max={float(np.nanmax(f14.depths)):.3f} m"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
