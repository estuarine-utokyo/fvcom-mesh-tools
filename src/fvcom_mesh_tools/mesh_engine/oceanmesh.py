"""``oceanmesh`` engine adapter for ``fmesh-buildmesh``.

The OceanMesh2D Python port. Combines :func:`oceanmesh.feature_sizing_function`
(coastline-aware sizing) with :func:`oceanmesh.bathymetric_gradient_sizing_function`
(depth-gradient sizing), ties them together with
:func:`oceanmesh.enforce_mesh_gradation`, then runs DistMesh via
:func:`oceanmesh.generate_mesh`. Output is cleaned up with the standard
``make_mesh_boundaries_traversable`` / ``delete_boundary_faces`` /
``laplacian2`` post-processing chain.

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
        log(f"[oceanmesh] reading shoreline {coast_staged.name} ...")
        shore = om.Shoreline(str(coast_staged), region, hmin_deg)
        sdf = om.signed_distance_function(shore)

        log(f"[oceanmesh] reading DEM {Path(dem_path).name} ...")
        dem = om.DEM(str(dem_path), bbox=region, crs=4326)

        log("[oceanmesh] feature_sizing_function ...")
        edge_feat = om.feature_sizing_function(
            shore, sdf, max_edge_length=hmax_deg, crs=4326,
        )
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
            edge = om.enforce_mesh_gradation(
                om.compute_minimum([edge_feat, edge_grad]), gradation=gradation,
            )
        else:
            edge = om.enforce_mesh_gradation(edge_feat, gradation=gradation)

        log(f"[oceanmesh] generate_mesh (max_iter={max_iter}, seed={seed}) ...")
        points, cells = om.generate_mesh(sdf, edge, max_iter=max_iter, seed=seed)
        log(f"[oceanmesh] raw output: NP={points.shape[0]:,} NE={cells.shape[0]:,}")

        log("[oceanmesh] cleanup pipeline ...")
        points, cells = om.make_mesh_boundaries_traversable(points, cells)
        points, cells = om.delete_faces_connected_to_one_face(points, cells)
        points, cells = om.delete_boundary_faces(points, cells, min_qual=min_qual)
        points, cells = om.laplacian2(points, cells)

    return np.asarray(points, dtype=float), np.asarray(cells, dtype=np.int64)
