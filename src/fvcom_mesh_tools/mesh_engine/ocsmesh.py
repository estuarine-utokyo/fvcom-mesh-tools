"""``ocsmesh`` engine adapter for ``fmesh-buildmesh``.

Wraps OCSMesh + gmsh (CC0 + GPL-2.0+ runtime). Faster (~40 s on
Tokyo Bay vs ~25 min for the oceanmesh engine) but lower quality
(alpha mean ~0.85 vs 0.96; PoC #18). Useful as the *draft* engine
during parameter iteration and as the backend for
``fmesh-mesh-combine --strategy {overlap,neighbor}``.

Like the oceanmesh adapter, this module emits ``(points, cells)`` in
EPSG:4326 lon/lat. Depth interpolation, boundary classification,
quality passes, perpfix, and ``fort.14`` writing are the caller's
responsibility (see ``cli/buildmesh.py``).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.io import (
    filter_multipolygon_by_area,
    load_coastline_as_lines,
)


def build(
    *,
    dem_path: Path,
    coastline_paths: list[Path],
    bbox: tuple[float, float, float, float],
    hmin_m: float,
    hmax_m: float,
    zmax: float = 0.0,
    min_polygon_area_m2: float = 0.0,
    min_island_area_m2: float = 0.0,
    coast_target_size: float | None = None,
    coast_expansion_rate: float = 0.005,
    log: Callable[[str], None] = print,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a mesh with OCSMesh + gmsh.

    Parameters
    ----------
    dem_path
        DEM raster path. Must be CRS-tagged (use ``fmesh-subset-dem``
        to prepare a CF-tagged GeoTIFF from a global DEM).
    coastline_paths
        Optional coastline shapefiles. If non-empty, the lines are
        clipped to ``bbox``, projected to EPSG:4326, and fed to
        ``Hfun.add_feature`` for coastline-aware sizing.
    bbox
        ``(minlon, minlat, maxlon, maxlat)`` in degrees. Used for
        clipping coastline polylines; OCSMesh derives the wet domain
        from the raster geometry directly.
    hmin_m, hmax_m
        Target minimum / maximum element size in metres.
    zmax
        Geom: water domain is the region where DEM <= zmax (default 0).
    min_polygon_area_m2, min_island_area_m2
        Pre-mesh polygon / hole filtering applied to the wet domain
        before meshing. ``0`` disables (default).
    coast_target_size
        Target element size on the coastline polylines (default
        ``hmin_m``).
    coast_expansion_rate
        Expansion rate passed to ``Hfun.add_feature``. Larger values
        relax sizing faster away from the coast.
    log
        Logging hook; defaults to ``print``.

    Returns
    -------
    points, cells
        ``points`` is ``(NP, 2)`` lon/lat in EPSG:4326. ``cells`` is
        ``(NE, 3)`` 0-based vertex indices.

    Raises
    ------
    ValueError
        If the geom filter drops every polygon.
    RuntimeError
        If OCSMesh produces zero triangles.
    """
    from ocsmesh import Geom, Hfun, MeshDriver, Raster
    from pyproj import CRS, Transformer

    raster = Raster(str(dem_path))
    geom = Geom(raster, zmax=zmax)

    if min_polygon_area_m2 > 0 or min_island_area_m2 > 0:
        log(
            f"[ocsmesh] geom filter: min_polygon={min_polygon_area_m2:g} m^2, "
            f"min_island={min_island_area_m2:g} m^2"
        )
        mp = geom.get_multipolygon()
        n_polys_before = len(mp.geoms)
        n_holes_before = sum(len(list(p.interiors)) for p in mp.geoms)
        mp_filt = filter_multipolygon_by_area(
            mp,
            src_crs=raster.crs,
            min_polygon_area_m2=min_polygon_area_m2,
            min_island_area_m2=min_island_area_m2,
        )
        n_polys_after = len(mp_filt.geoms)
        n_holes_after = sum(len(list(p.interiors)) for p in mp_filt.geoms)
        log(
            f"[ocsmesh] geom filter: polygons "
            f"{n_polys_before} -> {n_polys_after}, "
            f"holes {n_holes_before} -> {n_holes_after}"
        )
        if n_polys_after == 0:
            raise ValueError(
                "Geom filter dropped every polygon; "
                "loosen --min-polygon-area-m2."
            )
        # Replace the raster geom with a polygon-backed one driven by
        # the filtered multipolygon. Mesher uses this for the wet
        # domain; Hfun keeps the raster for sizing. OCSMesh requires
        # an explicit CRS for the polygon variant.
        geom = Geom(mp_filt, crs=raster.crs)

    hfun = Hfun(raster, hmin=hmin_m, hmax=hmax_m)

    if coastline_paths:
        coast = load_coastline_as_lines(coastline_paths, bbox=bbox)
        n_lines = len(coast.geoms)
        if n_lines == 0:
            log(
                "[ocsmesh] coastline: no features inside DEM bbox; "
                "skipping add_feature."
            )
        else:
            target = coast_target_size if coast_target_size is not None else hmin_m
            log(
                f"[ocsmesh] coastline: {n_lines} line strings; "
                f"target_size={target:g} m  expansion_rate={coast_expansion_rate:g}"
            )
            hfun.add_feature(
                feature=coast,
                expansion_rate=coast_expansion_rate,
                target_size=target,
            )

    log("[ocsmesh] generate_mesh (gmsh) ...")
    driver = MeshDriver(geom, hfun=hfun, engine_name="gmsh")
    mesh = driver.run()

    coords = np.asarray(mesh.coord, dtype=np.float64)
    cells = np.asarray(mesh.triangles, dtype=np.int64)
    if cells.size == 0:
        raise RuntimeError("OCSMesh produced zero triangles.")

    # OCSMesh internally meshes in a metric CRS (UTM by default) and
    # exposes coords there even when the input DEM is geographic.
    # Project back to EPSG:4326 to honour the engine contract.
    src_crs = mesh.crs
    if src_crs is not None and not CRS(src_crs).equals(CRS.from_epsg(4326)):
        log(f"[ocsmesh] projecting coords {src_crs} -> EPSG:4326")
        transformer = Transformer.from_crs(
            src_crs, CRS.from_epsg(4326), always_xy=True,
        )
        lon, lat = transformer.transform(coords[:, 0], coords[:, 1])
        coords = np.column_stack([lon, lat])

    return coords, cells
