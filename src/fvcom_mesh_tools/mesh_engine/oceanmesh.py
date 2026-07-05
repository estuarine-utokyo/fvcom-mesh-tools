"""``oceanmesh`` engine adapter for ``fmesh-buildmesh``.

The OceanMesh2D Python port. Combines up to three sizing functions
via :func:`oceanmesh.compute_minimum`:

  * :func:`oceanmesh.feature_sizing_function` — coastline-aware sizing
    (always on).
  * :func:`oceanmesh.bathymetric_gradient_sizing_function` —
    depth-gradient sizing (on by default).
  * :func:`oceanmesh.wavelength_sizing_function` — CFL/celerity
    sizing ``dx ∝ T·√(g·h)/wl`` (off by default).

The composed sizing is then smoothed by
:func:`oceanmesh.enforce_mesh_gradation` and fed to DistMesh via
:func:`oceanmesh.generate_mesh`. Output is cleaned up with the standard
``make_mesh_boundaries_traversable`` / ``delete_boundary_faces`` /
``laplacian2`` post-processing chain. The build-time ``laplacian2``
call is wrapped with
:func:`fvcom_mesh_tools.mesh_clean.repair_flipped_elements` so any
inverted triangle that the unsupervised smoother might leave behind
(PoC #34 found 1 such triangle on Tokyo Bay when wavelength sizing
was on) is rolled back to its pre-smoothing state.

Used as the *primary* engine because PoC #18 showed alpha mean 0.96 and
``frac<20deg`` 0.03 % on Tokyo Bay vs. 0.85 / 1.13 % for OCSMesh+gmsh.
The trade is wall-clock: ~26 min vs. ~40 s on the same problem.
"""

from __future__ import annotations

import shutil
import tempfile
from collections.abc import Callable
from pathlib import Path

import numpy as np

WGS84_PRJ_WKT = (
    'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",'
    '6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",'
    "0.017453292519943295]]"
)


def _stage_shapefile(src: Path, dst_dir: Path) -> Path:
    """Copy a shapefile into ``dst_dir`` and ensure a WGS84 .prj sidecar.

    oceanmesh.Shoreline requires a CRS-tagged shapefile; some sources
    (notably MLIT C23 main file) ship without a ``.prj``. We never
    mutate the user's data tree.
    """
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src.name
    for ext in (".shp", ".shx", ".dbf", ".cpg"):
        s = src.with_suffix(ext)
        if s.exists():
            shutil.copyfile(s, dst.with_suffix(ext))
    prj = src.with_suffix(".prj")
    if prj.exists():
        shutil.copyfile(prj, dst.with_suffix(".prj"))
    else:
        dst.with_suffix(".prj").write_text(WGS84_PRJ_WKT, encoding="utf-8")
    return dst


def _m_to_deg_lat(m: float) -> float:
    return m / 110_574.0


def _m_to_deg_lon(m: float, lat: float) -> float:
    return m / (111_320.0 * float(np.cos(np.deg2rad(lat))))


def _meters_per_degree_at_lat(lat_deg: float) -> float:
    """Metres-per-degree at ``lat_deg``, **matching oceanmesh's
    internal convention**: the upstream sizing functions feed the
    latitude expressed in *degrees* directly into ``np.cos``, so we
    do the same to keep the metres-to-degrees scaling consistent
    across all sizing components when they are merged via
    ``compute_minimum``. The numerical value is therefore not the
    physically correct WGS84 metres-per-degree; it is the same scale
    factor the wavelength / bathymetric-gradient functions apply.
    """
    return (
        111132.92
        - 559.82 * np.cos(2 * lat_deg)
        + 1.175 * np.cos(4 * lat_deg)
        - 0.0023 * np.cos(6 * lat_deg)
    )


def courant_sizing_function(
    dem,
    *,
    target_courant: float = 0.7,
    timestep_s: float = 5.0,
    wave_amplitude_m: float = 2.0,
    min_edgelength: float | None = None,
    max_edge_length: float | None = None,
    gravity: float = 9.81,
    crs: str | int = "EPSG:4326",
):
    """Per-cell sizing such that the approximate Courant number is
    capped at ``target_courant`` for the given ``timestep_s``.

    The characteristic celerity is approximated from linear long-wave
    theory, exactly as in OceanMesh2D's MATLAB Courant constraint and
    the ``ocsmesh.add_courant_num_constraint`` reference. For a depth
    ``h`` (positive metres):

    * ocean (``h > nu``):
      ``c = nu * sqrt(g / h) + sqrt(g * h)`` (particle velocity from
      linear wave theory + long-wave celerity).
    * overland (``h <= nu``):
      ``c = 2 * sqrt(g * nu)`` (the linear approximation breaks down
      when ``h ~ nu`` so we use the standard overland surrogate).

    Maximum element size such that ``C = c * dt / dx <= target_C``:

        ``dx_max = c * dt / target_C``.

    Composes via :func:`oceanmesh.compute_minimum` alongside feature /
    gradient / wavelength sizing functions; the final mesh respects
    ``C <= target_C`` everywhere.

    The algorithm is a few-line analytical formula based on the
    documented OceanMesh2D recipe and is implemented here from first
    principles — no code is borrowed from ``ocsmesh`` (CC0) or
    ``oceanmesh`` (GPL-3.0).

    Parameters
    ----------
    dem
        :class:`oceanmesh.DEM` instance.
    target_courant
        Upper bound on the approximate Courant number. Default 0.7
        matches the FVCOM-friendly preset documented in
        ``docs/architecture.md`` § 2.
    timestep_s
        Reference time step (seconds). The size function is scaled so
        that ``C <= target_courant`` at this step. Default 5.0 s — a
        comfortable production FVCOM step in coastal applications.
    wave_amplitude_m
        Wave amplitude / surface-elevation amplitude in metres. The
        linear-wave-theory regime threshold (``h > wave_amplitude_m``)
        and the particle-velocity coefficient. Default 2.0 m matches
        ``ocsmesh.add_courant_num_constraint``.
    min_edgelength, max_edge_length
        Optional clamp on the output sizing. Units must match the
        ``crs`` (degrees for ``EPSG:4326``).
    gravity
        Gravitational acceleration in m / s². Default 9.81.
    crs
        Coordinate reference of the returned :class:`oceanmesh.Grid`.
        ``EPSG:4326`` (default) returns sizes in degrees, matching the
        rest of the oceanmesh sizing chain.

    Returns
    -------
    :class:`oceanmesh.Grid`
        Sizing grid with values in degrees (geographic) or metres
        (projected). Its ``hmin`` attribute is set if
        ``min_edgelength`` was supplied.
    """
    import oceanmesh as om

    if target_courant <= 0:
        raise ValueError(f"target_courant must be > 0, got {target_courant}")
    if timestep_s <= 0:
        raise ValueError(f"timestep_s must be > 0, got {timestep_s}")
    if wave_amplitude_m <= 0:
        raise ValueError(
            f"wave_amplitude_m must be > 0, got {wave_amplitude_m}"
        )

    lon, lat = dem.create_grid()
    tmpz = dem.eval((lon, lat))
    abs_h = np.abs(np.asarray(tmpz, dtype=float))
    abs_h = np.where(abs_h < 1.0, 1.0, abs_h)  # safety floor (matches om)

    nu = float(wave_amplitude_m)
    deep_mask = abs_h > nu
    sqrt_gh = np.sqrt(gravity * abs_h)
    u_mag_deep = nu * np.sqrt(gravity / abs_h)
    u_mag_shallow = np.sqrt(gravity * nu)
    char_vel = np.where(
        deep_mask, u_mag_deep + sqrt_gh, 2.0 * u_mag_shallow,
    )

    dx_max_m = char_vel * float(timestep_s) / float(target_courant)

    grid = om.Grid(
        bbox=dem.bbox, dx=dem.dx, dy=dem.dy,
        extrapolate=True, values=0.0, crs=crs,
    )
    if crs in ("EPSG:4326", 4326):
        mean_latitude = float(np.mean(dem.bbox[2:]))
        meters_per_degree = _meters_per_degree_at_lat(mean_latitude)
        grid.values = dx_max_m / meters_per_degree
        grid.dx = dem.dx
        grid.dy = dem.dy
    else:
        grid.values = dx_max_m

    if min_edgelength is not None:
        grid.values = np.where(
            grid.values < min_edgelength, min_edgelength, grid.values,
        )
        grid.hmin = min_edgelength
    if max_edge_length is not None:
        grid.values = np.where(
            grid.values > max_edge_length, max_edge_length, grid.values,
        )

    grid.build_interpolant()
    return grid


def build(
    *,
    dem_path: Path,
    coastline_paths: list[Path],
    bbox: tuple[float, float, float, float],
    hmin_m: float,
    hmax_m: float,
    zmax: float = 0.0,
    slope_parameter: float = 20.0,
    filter_quotient: int = 50,
    gradation: float = 0.15,
    max_iter: int = 50,
    seed: int = 0,
    min_qual: float = 0.15,
    use_bathymetric_gradient: bool = True,
    minimum_area_mult: float = 4.0,
    use_wavelength_sizing: bool = False,
    wavelength_period_s: float = 44712.0,   # M2 ≈ 12.42 h
    wavelength_grid_spacing: int = 100,
    use_courant_sizing: bool = False,
    courant_target: float = 0.7,
    courant_timestep_s: float = 5.0,
    courant_wave_amplitude_m: float = 2.0,
    high_fidelity: bool = False,
    high_fidelity_lines: Path | None = None,
    shoreline_h0_m: float | None = None,
    enforce_hmin_floor: bool = False,
    constrain_boundary: bool = False,
    obc_coarsen_line: list | None = None,
    obc_coarsen_size_m: float = 1600.0,
    obc_coarsen_radius_m: float = 10000.0,
    interest_region: list | None = None,
    outside_min_m: float = 1000.0,
    outside_blend_m: float = 5000.0,
    log: Callable[[str], None] = print,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a mesh with oceanmesh + DistMesh.

    Parameters
    ----------
    dem_path
        DEM raster path (must be CRS-tagged; see ``fmesh-subset-dem`` to
        prepare global DEMs).
    coastline_paths
        List of vector shapefile paths (lines or polygons). The first
        one is fed to ``om.Shoreline``; if you need to combine multiple
        sources, pre-merge them with :mod:`geopandas` first.
    bbox
        ``(minlon, minlat, maxlon, maxlat)`` in degrees; defines the
        meshing domain (independent of the DEM raster bounds).
    hmin_m, hmax_m
        Target minimum / maximum element size in metres (converted to
        degrees internally for the EPSG:4326 grid).
    slope_parameter, filter_quotient
        Tunables for ``bathymetric_gradient_sizing_function`` (see
        oceanmesh docs). Defaults follow OceanMesh2D recipes.
    gradation
        Maximum ratio of element-size change between neighbours
        (``enforce_mesh_gradation``). 0.15 is the OceanMesh2D default.
    max_iter, seed
        DistMesh iterations and PRNG seed (default seed=0 makes the
        output deterministic between runs).
    min_qual
        Threshold for ``delete_boundary_faces``: triangles with
        normalised quality below this get peeled off the boundary.
    use_bathymetric_gradient
        Skip the bathymetric-gradient sizing function and rely solely
        on coastline feature sizing. Useful when the DEM is too coarse
        for meaningful slope information.
    use_wavelength_sizing
        Add a CFL/celerity-based sizing function via
        ``oceanmesh.wavelength_sizing_function``: ``dx ∝ T·√(g·h) / wl``.
        Off by default. Useful when the gradient-based sizing under-
        resolves shoaling regions where ``∇h`` is small but ``h`` is
        too — typical in inner bays / harbours where the FVCOM CFL
        condition would otherwise force a small dt.
    wavelength_period_s
        Reference period in seconds for ``wavelength_sizing_function``
        (only used when ``use_wavelength_sizing`` is True). Default
        44712.0 s ≈ M2 (12.42 h).
    wavelength_grid_spacing
        ``wl`` parameter (number of cells per wavelength) for
        ``wavelength_sizing_function``. Default 100 — corresponds to
        ``dt = T/wl`` ≈ 7.5 min for M2, a comfortable FVCOM time step.
    use_courant_sizing
        Add a per-cell Courant-bound sizing contribution via
        :func:`courant_sizing_function`. Off by default. Sets the
        upper sizing envelope so the mesh respects ``C <= courant_target``
        at the requested ``courant_timestep_s``. Where
        ``wavelength_sizing_function`` ties dx to a *wavelength*
        (a property of the dynamics), Courant ties dx to an *explicit
        time step* (a property of the solver) — the two compose by
        ``compute_minimum``. PoC #39 quantifies the trade.
    courant_target
        Upper bound on the approximate Courant number used by
        Phase E sizing. Default 0.7 matches the FVCOM-friendly preset
        in ``docs/architecture.md`` § 2.
    courant_timestep_s
        Target time step in seconds. Default 5.0 s — a comfortable
        production FVCOM coastal step.
    courant_wave_amplitude_m
        Wave amplitude / surface-elevation amplitude (metres) used by
        the linear-wave-theory regime threshold. Default 2.0 m
        matches the OceanMesh2D / ``ocsmesh`` reference value.
    minimum_area_mult
        Forwarded to ``om.Shoreline``. Inner-shoreline features
        smaller than ``minimum_area_mult * h0**2`` (with ``h0`` being
        the per-CRS minimum edge length) are dropped. Default 4.0
        matches oceanmesh; raise it to filter out more islets when the
        coastline shapefile is over-detailed.
    log
        Logging hook; defaults to ``print``.
    """
    import oceanmesh as om

    if not coastline_paths:
        raise ValueError("oceanmesh engine requires at least one --coastline.")
    minlon, minlat, maxlon, maxlat = bbox
    om_bbox = (minlon, maxlon, minlat, maxlat)
    lat_mid = 0.5 * (minlat + maxlat)
    hmin_deg = float(min(_m_to_deg_lat(hmin_m), _m_to_deg_lon(hmin_m, lat_mid)))
    hmax_deg = float(_m_to_deg_lat(hmax_m))
    log(
        f"[oceanmesh] bbox={om_bbox}  hmin={hmin_m:g} m -> "
        f"{hmin_deg:.6f} deg  hmax={hmax_m:g} m -> {hmax_deg:.6f} deg"
    )

    region = om.Region(extent=om_bbox, crs=4326)

    with tempfile.TemporaryDirectory(prefix="fmesh_om_") as td:
        td_path = Path(td)
        if len(coastline_paths) > 1:
            log(
                f"[oceanmesh] WARN: passing only the first of "
                f"{len(coastline_paths)} coastline files to om.Shoreline; "
                "merge upstream if multiple sources are needed."
            )
        coast_staged = _stage_shapefile(Path(coastline_paths[0]), td_path)
        # Decouple the DOMAIN's shoreline resolution from the element
        # size: om.Shoreline simplifies the coast to its h0, and the
        # signed-distance domain inherits that simplification — the
        # root cause of the ~h0/6 boundary-conformity floor (PoC
        # #60-62). A smaller shoreline_h0_m keeps the domain faithful
        # while sizing stays at hmin (DistMesh's boundary projection
        # then lands coarse boundary nodes ON the detailed line).
        if shoreline_h0_m is None:
            shore_h0_deg = hmin_deg
        else:
            shore_h0_deg = float(min(
                _m_to_deg_lat(shoreline_h0_m),
                _m_to_deg_lon(shoreline_h0_m, lat_mid),
            ))
        log(
            f"[oceanmesh] reading shoreline {coast_staged.name}  "
            f"minimum_area_mult={minimum_area_mult:g}  "
            f"shoreline_h0={shore_h0_deg:.6f} deg ..."
        )
        shore = om.Shoreline(
            str(coast_staged), region, shore_h0_deg,
            minimum_area_mult=minimum_area_mult,
        )
        sdf = om.signed_distance_function(shore)

        log(f"[oceanmesh] reading DEM {Path(dem_path).name} ...")
        dem = om.DEM(str(dem_path), bbox=region, crs=4326)

        log("[oceanmesh] feature_sizing_function ...")
        edge_feat = om.feature_sizing_function(
            shore, sdf, max_edge_length=hmax_deg, crs=4326,
        )
        if enforce_hmin_floor:
            # feature_sizing has no lower bound (its min_edge_length
            # argument is NOT a floor — passing hmin there collapsed
            # the mesh in PoC #65), so nearshore sizes fall to
            # ~hmin/4 and boundary vertices pack at that scale.
            # Clamping the grid VALUES afterwards is the safe floor:
            # with a detailed domain (shoreline_h0_m) the boundary
            # vertices then space at hmin ON the detailed line.
            n_below = int((edge_feat.values < hmin_deg).sum())
            edge_feat.values = np.maximum(edge_feat.values, hmin_deg)
            log(
                f"[oceanmesh] hmin floor: raised {n_below:,} "
                f"feature-sizing cells to {hmin_deg:.6f} deg"
            )
        edge_components = [edge_feat]
        if use_bathymetric_gradient:
            log("[oceanmesh] bathymetric_gradient_sizing_function ...")
            edge_grad = om.bathymetric_gradient_sizing_function(
                dem,
                slope_parameter=slope_parameter,
                filter_quotient=filter_quotient,
                min_edge_length=hmin_deg,
                max_edge_length=hmax_deg,
                crs=4326,
            )
            edge_components.append(edge_grad)
        if use_wavelength_sizing:
            log(
                f"[oceanmesh] wavelength_sizing_function "
                f"(period={wavelength_period_s:g} s, "
                f"wl={wavelength_grid_spacing}) ..."
            )
            edge_wave = om.wavelength_sizing_function(
                dem,
                wl=int(wavelength_grid_spacing),
                period=float(wavelength_period_s),
                min_edgelength=hmin_deg,
                max_edge_length=hmax_deg,
                crs=4326,
            )
            edge_components.append(edge_wave)
        if use_courant_sizing:
            log(
                f"[oceanmesh] courant_sizing_function "
                f"(C_target={courant_target:g}, "
                f"dt={courant_timestep_s:g} s, "
                f"nu={courant_wave_amplitude_m:g} m) ..."
            )
            edge_courant = courant_sizing_function(
                dem,
                target_courant=float(courant_target),
                timestep_s=float(courant_timestep_s),
                wave_amplitude_m=float(courant_wave_amplitude_m),
                min_edgelength=hmin_deg,
                max_edge_length=hmax_deg,
                crs=4326,
            )
            edge_components.append(edge_courant)
        if len(edge_components) == 1:
            combined = edge_feat
        else:
            combined = om.compute_minimum(edge_components)
        if obc_coarsen_line:
            # Reference-mesh policy (goto2023, user 2026-07-05): the
            # open boundary carries COARSE simple elements (~1.6 km
            # spacing at the arc) because the boundary zone is
            # numerically delicate and outside the region of
            # interest. Impose a size floor ramping from
            # obc_coarsen_size_m at the arc down to nothing at
            # obc_coarsen_radius_m inland; bathymetric-gradient
            # refinement over the Uraga canyon otherwise pins the
            # mouth at hmin.
            import shapely
            from shapely.geometry import LineString

            arcline = LineString(
                [(float(q[0]), float(q[1])) for q in obc_coarsen_line]
            )
            xg, yg = combined.create_grid()
            pts_g = shapely.points(xg.ravel(), yg.ravel())
            d_deg = shapely.distance(pts_g, arcline).reshape(xg.shape)
            deg_per_m = float(_m_to_deg_lat(1.0))
            d_m = d_deg / deg_per_m
            size_deg = obc_coarsen_size_m * deg_per_m
            ramp = size_deg * np.clip(
                1.0 - d_m / obc_coarsen_radius_m, 0.0, 1.0,
            )
            n_raised = int((combined.values < ramp).sum())
            combined.values = np.maximum(combined.values, ramp)
            log(
                f"[oceanmesh] OBC coarsening: raised {n_raised:,} "
                f"cells (size {obc_coarsen_size_m:g} m at the arc, "
                f"radius {obc_coarsen_radius_m:g} m)"
            )
        if interest_region:
            # Reference-mesh architecture (goto2023 MATLAB scripts:
            # 3-level nest, fine ONLY inside the interest polygon;
            # outside is 1-2 km+). Outside the polygon the sizing is
            # floored at outside_min_m, blended over outside_blend_m.
            import shapely
            from shapely.geometry import Polygon

            poly = Polygon(
                [(float(q[0]), float(q[1])) for q in interest_region]
            )
            xg2, yg2 = combined.create_grid()
            pts2 = shapely.points(xg2.ravel(), yg2.ravel())
            d_out_deg = shapely.distance(pts2, poly).reshape(xg2.shape)
            deg_per_m2 = float(_m_to_deg_lat(1.0))
            d_out_m = d_out_deg / deg_per_m2
            floor2 = (outside_min_m * deg_per_m2) * np.clip(
                d_out_m / max(outside_blend_m, 1.0), 0.0, 1.0,
            )
            n_up = int((combined.values < floor2).sum())
            combined.values = np.maximum(combined.values, floor2)
            log(
                f"[oceanmesh] interest-region floor: raised {n_up:,} "
                f"cells outside the polygon "
                f"(min {outside_min_m:g} m, blend {outside_blend_m:g} m)"
            )
        edge = om.enforce_mesh_gradation(combined, gradation=gradation)

        pfix = None
        egfix = None
        if constrain_boundary:
            # OM2D-parity: the shoreline's own polylines (h0-detail)
            # resampled at local h become pfix + CONSTRAINED EDGE
            # chains (egfix) — the boundary of the triangulation IS
            # the shoreline, no cleanup retreat, no post-hoc snap.
            pfix, egfix = om.shoreline_to_fixed_points(
                shore, edge, return_edges=True,
            )
            log(
                f"[oceanmesh] constrained boundary: {len(pfix):,} "
                f"fixed points, {len(egfix):,} fixed edges"
            )
            if high_fidelity_lines is not None:
                # Extra interior seed lines (e.g. water skeleton):
                # points only, no edge constraints.
                import geopandas as gpd

                gdf = gpd.read_file(high_fidelity_lines)
                raw_lines = []
                for geom in gdf.geometry:
                    if geom is None or geom.is_empty:
                        continue
                    boundary = (
                        geom.boundary
                        if geom.geom_type.endswith("Polygon") else geom
                    )
                    geoms = (
                        boundary.geoms
                        if hasattr(boundary, "geoms") else [boundary]
                    )
                    for g in geoms:
                        raw_lines.append(
                            np.asarray(g.coords, dtype=float)
                        )
                seed_pts = om.polylines_to_fixed_points(raw_lines, edge)
                if len(seed_pts):
                    pfix = np.vstack([pfix, seed_pts])
                log(
                    f"[oceanmesh] + {len(seed_pts):,} interior seed "
                    "points (no edge constraints)"
                )
            if len(pfix) == 0:
                pfix = None
                egfix = None
        elif high_fidelity_lines is not None:
            # Constrain to the RAW (unsimplified) vectors: Shoreline's
            # mainland/inner are h0-simplified, so shore-derived pfix
            # cannot beat the simplification floor (PoC #61: ~50 m at
            # 300 m). Equivalent to OceanMesh2D high_fidelity=2 with
            # local-h resampling.
            import geopandas as gpd

            gdf = gpd.read_file(high_fidelity_lines)
            raw_lines = []
            for geom in gdf.geometry:
                if geom is None or geom.is_empty:
                    continue
                boundary = (
                    geom.boundary if geom.geom_type.endswith("Polygon") else geom
                )
                geoms = (
                    boundary.geoms if hasattr(boundary, "geoms") else [boundary]
                )
                for g in geoms:
                    raw_lines.append(np.asarray(g.coords, dtype=float))
            pfix = om.polylines_to_fixed_points(raw_lines, edge)
            log(
                f"[oceanmesh] high-fidelity (raw lines): {len(pfix):,} "
                f"fixed points from {len(raw_lines):,} polylines"
            )
            if len(pfix) == 0:
                pfix = None
        elif high_fidelity:
            # OceanMesh2D V6.0 #264 port: fix resampled shoreline points
            # into the DistMesh iteration so the boundary lies exactly
            # on the (locally-resampled) shoreline.
            pfix = om.shoreline_to_fixed_points(shore, edge)
            log(
                f"[oceanmesh] high-fidelity: {len(pfix):,} fixed "
                "shoreline points"
            )
            if len(pfix) == 0:
                pfix = None

        log(f"[oceanmesh] generate_mesh (max_iter={max_iter}, seed={seed}) ...")
        points, cells = om.generate_mesh(
            sdf, edge, max_iter=max_iter, seed=seed, pfix=pfix,
            egfix=egfix,
        )
        log(f"[oceanmesh] raw output: NP={points.shape[0]:,} NE={cells.shape[0]:,}")

        log("[oceanmesh] cleanup pipeline ...")
        points, cells = om.make_mesh_boundaries_traversable(points, cells)
        if not constrain_boundary:
            points, cells = om.delete_faces_connected_to_one_face(points, cells)
            points, cells = om.delete_boundary_faces(points, cells, min_qual=min_qual)
        else:
            # The quality-based boundary deleters are exactly the
            # "cleanup retreat" (PoC #77/#92: boundary pulled ~0.7 h
            # inside). With CDT-constrained chains the boundary
            # triangles are mesh1d-spaced by construction; keep only
            # the manifold guarantee.
            log(
                "[oceanmesh] constrained boundary: skipping "
                "delete_boundary_faces / one-face deletion"
            )

        # ``om.laplacian2`` converges on edge-length stability but does
        # not check signed area, so it can leave a few inverted
        # triangles behind (PoC #34 surfaced 1 such triangle when
        # wavelength sizing was on). Wrap with the same flip-rollback
        # used by Phase G in mesh_clean.
        from fvcom_mesh_tools.mesh_clean import repair_flipped_elements

        pre = np.asarray(points, dtype=float).copy()
        pfix_idx = None
        if pfix is not None:
            # Cleanup may have renumbered/dropped vertices; re-locate
            # the surviving fixed points by exact coordinate match so
            # the smoother cannot move them off the shoreline.
            from scipy.spatial import cKDTree

            d, idx = cKDTree(pre).query(pfix)
            pfix_idx = np.unique(idx[d < 1e-9])
            log(
                f"[oceanmesh] high-fidelity: locking {len(pfix_idx):,}"
                f"/{len(pfix):,} surviving fixed points in laplacian2"
            )
        smoothed, cells = om.laplacian2(points, cells, pfix=pfix_idx)
        smoothed = np.asarray(smoothed, dtype=float)
        cells_arr = np.asarray(cells, dtype=int)
        repaired, repair_info = repair_flipped_elements(pre, smoothed, cells_arr)
        if repair_info["n_flipped_post_smooth"] > 0:
            log(
                f"[oceanmesh] laplacian2 repair: rolled back "
                f"{repair_info['n_nodes_rolled_back']} node(s) to clear "
                f"{repair_info['n_flipped_post_smooth']} flipped triangle(s)"
                f"{' (full rollback)' if repair_info['full_rollback'] else ''}"
            )
        points = repaired

    return np.asarray(points, dtype=float), np.asarray(cells, dtype=np.int64)
