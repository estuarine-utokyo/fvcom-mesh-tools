"""OM2D-parity multiscale build engine (v6 stage 1).

Two nests from per-nest coastline sources (the validated PoC #124
configuration): an outer nest over the whole walled domain built
from the per-polygon-smoothed coastline, and an inner nest on the
interest polygon built from the detailed engineered coastline.
Sizing runs through oceanmesh.finalize_sizing (banded bounds +
two-sided CFL); generation through generate_multiscale_mesh.

Returns the raw mesh in EPSG:4326 like the single-scale engine, so
all downstream stages (finish/obc/siteops/polish/qa/export) apply
unchanged.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

DEG = 1.0 / 111194.9266

__all__ = ["build_multiscale"]


def build_multiscale(
    *,
    land_detailed_shp: Path,
    land_outer_shp: Path,
    dem_path: Path,
    bbox: tuple[float, float, float, float],
    inner_polygon: list,
    outer: dict[str, Any] | None = None,
    inner: dict[str, Any] | None = None,
    courant: dict[str, Any] | None = None,
    seed: int = 0,
    junction_constraints: list | None = None,
    log=print,
) -> tuple[np.ndarray, np.ndarray]:
    """Build the two-nest multiscale mesh; returns (points, cells)
    in EPSG:4326.

    ``outer``/``inner`` sizing dicts (defaults = PoC #124):
      outer: h0_m 1000, hmin_m 1000, hmax_m 2000, wl 30, grade 0.10
      inner: h0_m 100, hmin_m 400, hmax_m 500, r 3, grade 0.10
    ``courant``: {"timestep": 16.0, "max": 0.5} (finalize CFL).
    """
    import oceanmesh as om
    from oceanmesh import Region

    o = {"h0_m": 1000.0, "hmin_m": 1000.0, "hmax_m": 2000.0,
         "wl": 30, "r": 3, "grade": 0.10}
    o.update(outer or {})
    i = {"h0_m": 100.0, "hmin_m": 400.0, "hmax_m": 500.0,
         "r": 3, "grade": 0.10}
    i.update(inner or {})
    cfl = {"timestep": 16.0, "max": 0.5}
    cfl.update(courant or {})

    region_outer = Region((bbox[0], bbox[2], bbox[1], bbox[3]), 4326)
    poly_inner = np.asarray(
        list(inner_polygon) + [inner_polygon[0]], dtype=float
    )

    log(f"[ms] outer Shoreline h0={o['h0_m']:g} m "
        f"({Path(land_outer_shp).name})")
    shore_o = om.Shoreline(str(land_outer_shp), region_outer.bbox,
                           o["h0_m"] * DEG)
    sdf_o = om.signed_distance_function(shore_o)
    dem = om.DEM(str(dem_path), bbox=region_outer)

    comps_o = [om.feature_sizing_function(
        shore_o, sdf_o, r=o["r"], max_edge_length=o["hmax_m"] * DEG)]
    if o.get("wl"):
        comps_o.append(om.wavelength_sizing_function(
            dem, wl=o["wl"], min_edgelength=o["hmin_m"] * DEG,
            max_edge_length=o["hmax_m"] * DEG))
    edge_o, dt_o = om.finalize_sizing(
        comps_o, dem=dem,
        hmin=o["hmin_m"], max_edge_length=o["hmax_m"],
        gradation=o["grade"], courant=cfl,
    )
    log(f"[ms] outer sizing done (dt={dt_o})")

    log(f"[ms] inner Shoreline h0={i['h0_m']:g} m, polygon boubox "
        f"({Path(land_detailed_shp).name})")
    shore_i = om.Shoreline(str(land_detailed_shp), poly_inner,
                           i["h0_m"] * DEG)
    sdf_i = om.signed_distance_function(shore_i)
    edge_i, dt_i = om.finalize_sizing(
        [om.feature_sizing_function(
            shore_i, sdf_i, r=i["r"],
            max_edge_length=i["hmax_m"] * DEG)],
        dem=dem,
        hmin=i["hmin_m"], max_edge_length=i["hmax_m"],
        gradation=i["grade"], courant=cfl,
    )
    log(f"[ms] inner sizing done (dt={dt_i})")

    ms_kw = {}
    if junction_constraints:
        pf = np.asarray(
            [q for pair in junction_constraints for q in pair],
            dtype=float,
        )
        eg = np.asarray(
            [[2 * k, 2 * k + 1]
             for k in range(len(junction_constraints))], dtype=int,
        )
        ms_kw = {"pfix": pf, "egfix": eg}
        log(f"[ms] junction constraints: {len(pf)} pfix, "
            f"{len(eg)} egfix")
    log("[ms] generate_multiscale_mesh ...")
    points, cells = om.generate_multiscale_mesh(
        [sdf_o, sdf_i], [edge_o, edge_i], seed=seed, **ms_kw,
    )
    log(f"[ms] raw NP={len(points):,} NE={len(cells):,}")

    points, cells = om.make_mesh_boundaries_traversable(points, cells)
    points, cells = om.delete_faces_connected_to_one_face(points, cells)
    points, cells = om.laplacian2(points, cells)

    # keep only the largest connected component (wall-band remnants)
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    e = np.vstack([cells[:, [0, 1]], cells[:, [1, 2]],
                   cells[:, [2, 0]]])
    A = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])),
                   shape=(len(points), len(points)))
    ncc, lab = connected_components(A, directed=False)
    if ncc > 1:
        keep_lab = np.argmax(np.bincount(lab))
        cells = cells[(lab[cells] == keep_lab).all(axis=1)]
        from oceanmesh.fix_mesh import fix_mesh

        points, cells, _ = fix_mesh(points, cells, delete_unused=True)
        log(f"[ms] kept largest of {ncc} components "
            f"-> NP={len(points):,}")

    from oceanmesh.mesh_improve import area_length_quality

    q = area_length_quality(points, cells)
    log(f"[ms] clean NP={len(points):,} NE={len(cells):,} "
        f"AL-qual min/mean = {q.min():.3f}/{q.mean():.3f}")

    # provisional positive-down depths (depth design itself is out
    # of scope; downstream stages expect non-empty depths)
    depths = om.interp_bathymetry(points, cells, dem,
                                  method="cell-averaging",
                                  min_depth=1.0)
    return points, cells, depths
