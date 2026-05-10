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
        log(
            f"[oceanmesh] reading shoreline {coast_staged.name}  "
            f"minimum_area_mult={minimum_area_mult:g} ..."
        )
        shore = om.Shoreline(
            str(coast_staged), region, hmin_deg,
            minimum_area_mult=minimum_area_mult,
        )
        sdf = om.signed_distance_function(shore)

        log(f"[oceanmesh] reading DEM {Path(dem_path).name} ...")
        dem = om.DEM(str(dem_path), bbox=region, crs=4326)

        log("[oceanmesh] feature_sizing_function ...")
        edge_feat = om.feature_sizing_function(
            shore, sdf, max_edge_length=hmax_deg, crs=4326,
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
        if len(edge_components) == 1:
            edge = om.enforce_mesh_gradation(edge_feat, gradation=gradation)
        else:
            edge = om.enforce_mesh_gradation(
                om.compute_minimum(edge_components), gradation=gradation,
            )

        log(f"[oceanmesh] generate_mesh (max_iter={max_iter}, seed={seed}) ...")
        points, cells = om.generate_mesh(sdf, edge, max_iter=max_iter, seed=seed)
        log(f"[oceanmesh] raw output: NP={points.shape[0]:,} NE={cells.shape[0]:,}")

        log("[oceanmesh] cleanup pipeline ...")
        points, cells = om.make_mesh_boundaries_traversable(points, cells)
        points, cells = om.delete_faces_connected_to_one_face(points, cells)
        points, cells = om.delete_boundary_faces(points, cells, min_qual=min_qual)

        # ``om.laplacian2`` converges on edge-length stability but does
        # not check signed area, so it can leave a few inverted
        # triangles behind (PoC #34 surfaced 1 such triangle when
        # wavelength sizing was on). Wrap with the same flip-rollback
        # used by Phase G in mesh_clean.
        from fvcom_mesh_tools.mesh_clean import repair_flipped_elements

        pre = np.asarray(points, dtype=float).copy()
        smoothed, cells = om.laplacian2(points, cells)
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
