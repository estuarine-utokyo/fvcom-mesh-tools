"""``fmesh-buildmesh`` CLI: DEM -> mesh engine -> classified, perp-fixed fort.14.

Pipeline (single command, no MATLAB intermediate):

    DEM (NetCDF/GeoTIFF)
        -> dem.bbox.read              (lon/lat extent for boundary classify)
        -> mesh_engine.build(engine)  (oceanmesh DistMesh; ocsmesh+gmsh deprecated)
        -> dem.interp.at_points       (per-node depth, mesher-agnostic)
        -> classify boundaries by DEM bbox proximity
        -> (optional) edge-swap + Laplacian quality pass
        -> (optional) longest-edge bisection refine
        -> align_open_boundary_first_ring (perpfix)
        -> fort.14 (FVCOM/ADCIRC convention, depth +down)

Both mesh engines emit ``(points, cells)`` in EPSG:4326 lon/lat; everything
downstream is mesher-agnostic. The intent is "give us a fort.14 the FVCOM
harness can actually load". Quality post-processing is best-effort: it
monotonically improves mean alpha-quality but plateaus at a level that
depends on the initial size function (PoC #10). Driving the bad-element
fraction to zero requires adaptive sizing, which
``docs/python_pipeline_gap_analysis.md`` tracks.

.. note::
   ``--engine ocsmesh`` is **deprecated** and emits a
   :class:`DeprecationWarning` plus an stderr notice when selected.
   It is retained for one release to give downstream callers time
   to migrate. Production meshes should use ``--engine oceanmesh``
   (the default). See ``docs/engine_complementarity.md`` for the
   rationale (PoC #30 quality gap and Triangle-backend limitation).
   Library-level use of ocsmesh (``ops.combine_mesh``, ``utils``,
   ``Raster``) is unaffected.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms import (
    add_river_inflow_segments,
    align_open_boundary_first_ring,
    alpha_quality,
    classify_boundaries_by_bbox,
    laplacian_smooth,
    min_interior_angle,
    refine_bad_triangles,
    signed_areas,
    swap_edges_for_quality,
)
from fvcom_mesh_tools.io import (
    Fort14Mesh,
    load_river_points,
    write_fort14,
)

EARTH_R_M = 6_371_000.0

# Threshold below which a positive signed area is treated as a flipped
# triangle. Belt-and-suspenders with the writer's ``:.15f`` precision:
# any triangle with ``sa <= EPSILON_SA`` is at risk of being pancaked
# or flipped by FP rounding through any downstream writer/reader, so
# we revert it before write. 1e-12 deg^2 at mid-latitudes is ~0.01 m^2
# - a sliver, not a real element.
EPSILON_SA = 1.0e-12


def _deg_per_metre(lat_deg: float) -> float:
    """Conservative degrees-per-metre at ``lat_deg`` (longitude direction)."""
    return 1.0 / (EARTH_R_M * np.cos(np.deg2rad(lat_deg)) * np.pi / 180.0)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="fmesh-buildmesh",
        description=(
            "Generate a FVCOM-ready fort.14 from a single DEM. "
            "Mesh is produced by the selected --engine "
            "(oceanmesh DistMesh by default, ocsmesh+gmsh for fast "
            "drafts); depths are interpolated from the same DEM; the "
            "outer ring is split into open and land segments by "
            "proximity to the DEM bounding box; an optional first-ring "
            "perpendicularity fix is applied."
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
        "--engine", choices=["oceanmesh", "ocsmesh"], default="oceanmesh",
        help=(
            "Mesh-generation backend (default: oceanmesh). 'oceanmesh' "
            "is the OceanMesh2D Python port - higher quality "
            "(alpha~0.96) but slower (~25 min on Tokyo Bay). 'ocsmesh' "
            "is DEPRECATED and slated for removal in a future release "
            "(see docs/engine_complementarity.md): it uses OCSMesh+gmsh "
            "- alpha~0.85, max valence 26, and PoC #30 confirmed the "
            "Triangle escape hatch is unusable. Library use of ocsmesh "
            "(ops.combine_mesh, utils, Raster) is unaffected."
        ),
    )
    p.add_argument(
        "--interp-method", choices=["linear", "nearest"], default="linear",
        help="Depth-interpolation method (default: linear).",
    )
    p.add_argument(
        "--om-slope-parameter", type=float, default=20.0,
        help=(
            "[oceanmesh] slope_parameter for "
            "bathymetric_gradient_sizing_function (default: 20)."
        ),
    )
    p.add_argument(
        "--om-gradation", type=float, default=0.15,
        help=(
            "[oceanmesh] gradation for enforce_mesh_gradation "
            "(default: 0.15, OceanMesh2D default)."
        ),
    )
    p.add_argument(
        "--om-max-iter", type=int, default=50,
        help="[oceanmesh] DistMesh max iterations (default: 50).",
    )
    p.add_argument(
        "--om-seed", type=int, default=0,
        help="[oceanmesh] DistMesh PRNG seed (default: 0, deterministic).",
    )
    p.add_argument(
        "--om-no-bathymetric-gradient", action="store_true",
        help=(
            "[oceanmesh] Skip the bathymetric-gradient sizing function. "
            "Useful when the DEM is too coarse to drive slope-based "
            "sizing meaningfully."
        ),
    )
    p.add_argument(
        "--om-minimum-area-mult", type=float, default=4.0,
        help=(
            "[oceanmesh] minimum_area_mult forwarded to om.Shoreline "
            "(default: 4.0). Inner-shoreline features smaller than "
            "minimum_area_mult * h0^2 are dropped; raise to coalesce "
            "small islets."
        ),
    )
    p.add_argument(
        "--om-wavelength-sizing", action="store_true",
        help=(
            "[oceanmesh] Add wavelength_sizing_function (CFL / shallow-"
            "water celerity) to the size composition: "
            "dx ∝ T·sqrt(g·h)/wl. Off by default. Useful when the "
            "gradient-based sizing under-resolves shoaling regions where "
            "∇h is small but h is too — typical in inner bays / harbours "
            "where the FVCOM CFL condition would otherwise force a "
            "small dt."
        ),
    )
    p.add_argument(
        "--om-wavelength-period", type=float, default=44712.0,
        help=(
            "[oceanmesh] Reference period in seconds for "
            "--om-wavelength-sizing. Default 44712.0 ≈ M2 (12.42 h). "
            "Halve to require a finer mesh / shorter dt."
        ),
    )
    p.add_argument(
        "--om-wavelength-grid-spacing", type=int, default=100,
        help=(
            "[oceanmesh] wl parameter (cells per wavelength) for "
            "--om-wavelength-sizing. Default 100; corresponds to "
            "dt = T/wl ≈ 7.5 min for M2 — a comfortable FVCOM time "
            "step. Raise for finer meshes / shorter dt; lower for "
            "draft work."
        ),
    )
    p.add_argument(
        "--om-courant-sizing", action="store_true",
        help=(
            "[oceanmesh] Add a per-cell Courant-bound sizing "
            "contribution. Sets the upper sizing envelope so the mesh "
            "respects approximate Courant number "
            "C = c_char * dt / dx <= --om-courant-target at the "
            "requested --om-courant-timestep. Where wavelength sizing "
            "ties dx to a wavelength (a property of the dynamics), "
            "Courant sizing ties dx to an explicit dt (a property of "
            "the solver) — the two compose by min."
        ),
    )
    p.add_argument(
        "--om-courant-target", type=float, default=0.7,
        help=(
            "[oceanmesh] Upper bound on the approximate Courant number "
            "for --om-courant-sizing. Default 0.7 (FVCOM-friendly "
            "preset). Lower for more conservative meshes."
        ),
    )
    p.add_argument(
        "--om-courant-timestep", type=float, default=5.0,
        help=(
            "[oceanmesh] Reference time step in seconds for "
            "--om-courant-sizing. Default 5.0 — a comfortable production "
            "FVCOM coastal step. Lower if the FVCOM run actually uses a "
            "smaller dt; raising it relaxes the upper sizing envelope."
        ),
    )
    p.add_argument(
        "--om-courant-wave-amplitude", type=float, default=2.0,
        help=(
            "[oceanmesh] Wave amplitude (metres) for the linear-wave-"
            "theory regime threshold used by --om-courant-sizing. "
            "Default 2.0 m matches the OceanMesh2D / ocsmesh reference."
        ),
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
        "--river-inflow-points", type=Path, action="append", default=[],
        metavar="PATH",
        help=(
            "Vector / CSV file of river-mouth points (lon/lat). Each "
            "point is snapped to the nearest land-boundary node and "
            "the surrounding land segment is split so the river "
            "occupies its own segment with --river-ibtype. May be "
            "passed multiple times to combine sources."
        ),
    )
    p.add_argument(
        "--river-segment-nodes", type=int, default=5, metavar="N",
        help="Land-boundary nodes per river segment, centred on the snap (default 5).",
    )
    p.add_argument(
        "--river-ibtype", type=int, default=21,
        help="ibtype for river segments in fort.14 (default 21 = FVCOM discharge).",
    )
    p.add_argument(
        "--river-snap-tol-m", type=float, default=None, metavar="METRES",
        help="Skip river points whose nearest land node is farther than this.",
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
        "--refine-min-angle", type=float, default=0.0, metavar="DEG",
        help=(
            "After --quality-pass, run longest-edge bisection on every "
            "triangle whose minimum interior angle is below this "
            "threshold. 0 disables (default). Typical: 15-20 deg."
        ),
    )
    p.add_argument(
        "--refine-max-passes", type=int, default=5,
        help="Cap on refine passes (default 5; early-stops on regression).",
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


_OCSMESH_DEPRECATION_MESSAGE = (
    "--engine ocsmesh is DEPRECATED and slated for removal in a future "
    "release. Use --engine oceanmesh (the default) for production "
    "meshes. See docs/engine_complementarity.md for the rationale: "
    "ocsmesh+gmsh produces alpha~0.85 / max valence 26 (vs 0.96 / 9 "
    "for oceanmesh on the same input), and ocsmesh's Triangle backend "
    "rejects raster-driven varying sizing (PoC #30), so gmsh cannot "
    "be cheaply replaced. Library-level use of ocsmesh "
    "(ops.combine_mesh, utils, Raster) is unaffected."
)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.engine == "ocsmesh":
        # Visible on stderr regardless of warning filters; also emit a
        # DeprecationWarning so any harness that promotes warnings to
        # errors picks it up.
        import warnings

        print(
            f"[buildmesh] WARNING: {_OCSMESH_DEPRECATION_MESSAGE}",
            file=sys.stderr,
        )
        warnings.warn(
            _OCSMESH_DEPRECATION_MESSAGE,
            DeprecationWarning,
            stacklevel=2,
        )

    if not args.dem.exists():
        print(f"DEM not found: {args.dem}", file=sys.stderr)
        return 2
    if args.hmin <= 0 or args.hmax < args.hmin:
        print("--hmin must be > 0 and --hmax must be >= --hmin.", file=sys.stderr)
        return 2
    if args.om_wavelength_sizing:
        if args.om_wavelength_period <= 0:
            print("--om-wavelength-period must be > 0.", file=sys.stderr)
            return 2
        if args.om_wavelength_grid_spacing < 1:
            print("--om-wavelength-grid-spacing must be >= 1.", file=sys.stderr)
            return 2
    if args.om_courant_sizing:
        if not 0.0 < args.om_courant_target <= 1.5:
            print(
                "--om-courant-target must be in (0, 1.5] (typical "
                "FVCOM-friendly value 0.7).",
                file=sys.stderr,
            )
            return 2
        if args.om_courant_timestep <= 0:
            print("--om-courant-timestep must be > 0.", file=sys.stderr)
            return 2
        if args.om_courant_wave_amplitude <= 0:
            print(
                "--om-courant-wave-amplitude must be > 0.", file=sys.stderr,
            )
            return 2

    args.output.parent.mkdir(parents=True, exist_ok=True)

    log = (lambda *a, **k: None) if args.quiet else print  # noqa: E731

    log(f"[buildmesh] DEM: {args.dem}")
    log(
        f"[buildmesh] engine={args.engine}  "
        f"hmin={args.hmin:g} m  hmax={args.hmax:g} m  zmax={args.zmax:g}"
    )

    # DEM bbox drives boundary classification regardless of engine
    # choice. Reading it through the dem subpackage keeps rasterio
    # confined to one entry point.
    from fvcom_mesh_tools.dem.bbox import read as read_dem_bbox

    xmin, ymin, xmax, ymax = read_dem_bbox(args.dem)
    log(f"[buildmesh] DEM bbox: x[{xmin:.6f}, {xmax:.6f}]  y[{ymin:.6f}, {ymax:.6f}]")

    # Engine-specific kwargs. The dispatcher forwards these as
    # **engine_kwargs to the per-engine adapter.
    if args.engine == "oceanmesh":
        if not args.coastline:
            print(
                "oceanmesh engine requires at least one --coastline.",
                file=sys.stderr,
            )
            return 5
        if args.min_polygon_area_m2 > 0 or args.min_island_area_m2 > 0:
            log(
                "[buildmesh] note: --min-polygon-area-m2 / "
                "--min-island-area-m2 are ocsmesh-only and ignored "
                "by the oceanmesh engine."
            )
        engine_kwargs = dict(
            zmax=args.zmax,
            slope_parameter=args.om_slope_parameter,
            gradation=args.om_gradation,
            max_iter=args.om_max_iter,
            seed=args.om_seed,
            use_bathymetric_gradient=not args.om_no_bathymetric_gradient,
            minimum_area_mult=args.om_minimum_area_mult,
            use_wavelength_sizing=args.om_wavelength_sizing,
            wavelength_period_s=args.om_wavelength_period,
            wavelength_grid_spacing=args.om_wavelength_grid_spacing,
            use_courant_sizing=args.om_courant_sizing,
            courant_target=args.om_courant_target,
            courant_timestep_s=args.om_courant_timestep,
            courant_wave_amplitude_m=args.om_courant_wave_amplitude,
        )
    else:  # ocsmesh; argparse choices guard the value space
        engine_kwargs = dict(
            zmax=args.zmax,
            min_polygon_area_m2=args.min_polygon_area_m2,
            min_island_area_m2=args.min_island_area_m2,
            coast_target_size=args.coast_target_size,
            coast_expansion_rate=args.coast_expansion_rate,
        )

    from fvcom_mesh_tools.dem.interp import at_points as interp_dem_at_points
    from fvcom_mesh_tools.mesh_engine import build as build_engine

    t0 = time.perf_counter()
    try:
        points, cells = build_engine(
            args.engine,
            dem_path=args.dem,
            coastline_paths=list(args.coastline),
            bbox=(xmin, ymin, xmax, ymax),
            hmin_m=args.hmin,
            hmax_m=args.hmax,
            log=log,
            **engine_kwargs,
        )
    except ValueError as e:
        print(str(e), file=sys.stderr)
        return 4
    except RuntimeError as e:
        print(str(e), file=sys.stderr)
        return 3
    log(
        f"[buildmesh] mesh generation ({args.engine}): "
        f"{time.perf_counter() - t0:.2f} s"
    )

    log(f"[buildmesh] interpolating depths (method={args.interp_method}) ...")
    t1 = time.perf_counter()
    depths = interp_dem_at_points(args.dem, points, method=args.interp_method)
    log(f"[buildmesh] depth interpolation: {time.perf_counter() - t1:.2f} s")

    coords = points
    elements = cells

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

    if args.river_inflow_points:
        river_pts = load_river_points(args.river_inflow_points)
        log(
            f"[buildmesh] river inflow: {river_pts.shape[0]} points; "
            f"n_per_river={args.river_segment_nodes}, "
            f"ibtype={args.river_ibtype}"
        )
        f14, river_info = add_river_inflow_segments(
            f14,
            river_pts,
            n_nodes_per_river=args.river_segment_nodes,
            river_ibtype=args.river_ibtype,
            snap_tol_m=args.river_snap_tol_m,
        )
        for r in river_info["rivers"]:
            log(
                f"[buildmesh]   river @ ({r['point'][0]:.4f}, {r['point'][1]:.4f}) "
                f"-> node {r['snapped_node']}, dist={r['dist_m']:.0f} m, "
                f"river_n={r['river_n_nodes']}"
            )
        for s in river_info["skipped"]:
            log(
                f"[buildmesh]   river @ {s['point']} SKIPPED "
                f"(dist={s['dist_m']:.0f} m > snap_tol)"
            )
        log(
            f"[buildmesh] river inflow: now "
            f"{len(f14.land_boundaries)} land segments "
            f"({sum(1 for ib, _ in f14.land_boundaries if ib == args.river_ibtype)} "
            f"with ibtype={args.river_ibtype})"
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

    if args.refine_min_angle > 0:
        log(
            f"[buildmesh] refine: min-angle threshold {args.refine_min_angle:g} deg, "
            f"max_passes={args.refine_max_passes}"
        )
        bad_pre = float((min_interior_angle(f14) < args.refine_min_angle).mean()) * 100
        np_pre = f14.n_nodes
        f14, refine_info = refine_bad_triangles(
            f14,
            min_angle_threshold=args.refine_min_angle,
            max_passes=args.refine_max_passes,
        )
        bad_post = float((min_interior_angle(f14) < args.refine_min_angle).mean()) * 100
        log(
            f"[buildmesh] refine: passes={refine_info['passes']} "
            f"({refine_info['stop_reason']}); "
            f"frac<{args.refine_min_angle:g}deg {bad_pre:.2f}% -> {bad_post:.2f}%; "
            f"NP {np_pre:,} -> {f14.n_nodes:,} "
            f"(+{refine_info['total_nodes_inserted']:,} nodes)"
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
        # narrow triangles or pancake near-degenerate ones. Iterate
        # the revert: each round undoes any moved node that
        # participates in a flipped or near-degenerate triangle, then
        # re-checks signed areas. After the quality pass the
        # surrounding geometry is tighter, so a single round can be
        # insufficient. We use ``EPSILON_SA`` (not 0) so we also catch
        # triangles whose positive ``sa`` is small enough to round to
        # zero through the fort.14 writer.
        n_revert_total = 0
        for _ in range(5):
            bad_tris = signed_areas(f14) <= EPSILON_SA
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
        still_bad = int((signed_areas(f14) <= EPSILON_SA).sum())
        if still_bad:
            # Safety net: if iterative revert cannot clear the flips,
            # roll back perpfix entirely. The mesh is then exactly as
            # the quality pass left it - which is at least valid.
            log(
                f"[buildmesh] WARN: {still_bad} flipped or near-degenerate "
                f"triangles remain; undoing perpfix to keep the mesh valid."
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
