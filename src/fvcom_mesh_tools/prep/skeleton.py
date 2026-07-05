"""Water medial-axis seed lines for narrow-but-essential water.

At a true hmin floor the DistMesh seed lattice equals the floor, so
water narrower than ~2 floors receives no initial points and drops
out silently (mechanism M3 in docs/DESIGN_HISTORY.md). This module
rasterizes the water mask, extracts the medial axis where the local
half-width sits in a configurable band, and vectorizes the ridge into
polylines. Fed through ``--om-high-fidelity-lines`` they become
interior fixed seed points at local-h spacing — the automated version
of SMS hand-meshing 300-600 m harbor basins with 1-2 element rows
(independently validated by ADMESH+, Kang & Kubatko 2024, GMD).
"""

from __future__ import annotations

__all__ = ["water_skeleton_lines"]


def water_skeleton_lines(
    land_gdf,
    *,
    px_m: float = 50.0,
    half_width_range_m: tuple[float, float] = (140.0, 460.0),
    min_chain_px: int = 5,
    simplify_deg: float = 3e-4,
    utm_epsg: int | None = None,
):
    """Medial-axis polylines of the water complement of ``land_gdf``
    (EPSG:4326 polygons), kept where the local half-width is within
    ``half_width_range_m``. Returns LineStrings in EPSG:4326.

    ``px_m`` is the raster cell size (keep it well below the lower
    half-width bound); chains shorter than ``min_chain_px`` pixels
    are dropped.
    """
    import geopandas as gpd
    import numpy as np
    from shapely import contains_xy, make_valid, unary_union
    from shapely.geometry import LineString, box
    from skimage.morphology import medial_axis

    from fvcom_mesh_tools.prep.shoreline import auto_utm_epsg

    gdf = land_gdf
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)
    land_ll = unary_union([
        make_valid(g) for g in gdf.to_crs(4326).geometry
        if g is not None and not g.is_empty
    ])
    if utm_epsg is None:
        c = land_ll.centroid
        utm_epsg = auto_utm_epsg(c.x, c.y)
    land = gpd.GeoSeries([land_ll], crs=4326).to_crs(utm_epsg).iloc[0]

    minx, miny, maxx, maxy = land.bounds
    water = box(minx, miny, maxx, maxy).difference(land)

    nx = int((maxx - minx) / px_m) + 1
    ny = int((maxy - miny) / px_m) + 1
    xs = minx + px_m * (np.arange(nx) + 0.5)
    ys = miny + px_m * (np.arange(ny) + 0.5)
    xx, yy = np.meshgrid(xs, ys)
    mask = contains_xy(water, xx.ravel(), yy.ravel()).reshape(ny, nx)

    skel, dist = medial_axis(mask, return_distance=True)
    half_w = dist * px_m
    lo, hi = half_width_range_m
    keep = skel & (half_w >= lo) & (half_w <= hi)

    pts = {(int(i), int(j)) for i, j in zip(*np.where(keep))}

    def _nbrs(p):
        i, j = p
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if (di or dj) and (i + di, j + dj) in pts:
                    yield (i + di, j + dj)

    lines = []
    visited: set = set()
    for p in sorted(pts):
        if p in visited:
            continue
        chain = [p]
        visited.add(p)
        ends = [q for q in _nbrs(p) if q not in visited]
        for direction in ends[:2]:
            cur, prev = direction, p
            seg = []
            while cur and cur not in visited:
                visited.add(cur)
                seg.append(cur)
                nxt = [q for q in _nbrs(cur)
                       if q != prev and q not in visited]
                prev, cur = cur, (nxt[0] if nxt else None)
            if direction is ends[0]:
                chain = list(reversed(seg)) + chain
            else:
                chain = chain + seg
        if len(chain) >= min_chain_px:
            coords = [(minx + px_m * (j + 0.5), miny + px_m * (i + 0.5))
                      for i, j in chain]
            lines.append(LineString(coords))

    out = gpd.GeoDataFrame(geometry=lines, crs=utm_epsg).to_crs(4326)
    if simplify_deg > 0 and len(out):
        out.geometry = out.geometry.simplify(simplify_deg)
    return out.reset_index(drop=True)
