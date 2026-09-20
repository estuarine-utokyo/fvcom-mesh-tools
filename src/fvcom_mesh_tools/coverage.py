"""Find original water omitted by a mesh, in a caller-declared metric domain.

Optional dependencies: the existing io-vector and prep extras. This is a
geometric screening diagnostic, not a guarantee that a quality mesh can be built.
"""

from __future__ import annotations

import numpy as np

from fvcom_mesh_tools.io import Fort14Mesh


def _polygons(geometry):
    if geometry.geom_type == "Polygon":
        if not geometry.is_empty:
            yield geometry
    elif hasattr(geometry, "geoms"):
        for part in geometry.geoms:
            yield from _polygons(part)


def _summary(poly):
    from shapely.geometry import mapping

    point = poly.representative_point()
    return {
        "area_m2": float(poly.area),
        "bounds_m": list(poly.bounds),
        "xy_m": [point.x, point.y],
        "geometry": mapping(poly),
    }


def coverage_gaps(
    mesh: Fort14Mesh,
    original_land,
    domain,
    *,
    policy_land=None,
    pixel_size_m: float | None = None,
    min_area_m2: float = 100.0,
    max_raster_cells: int = 2_000_000,
    corridor_fraction: float = 0.8,
    h_target_m: float | None = None,
):
    """Return JSON-ready ranked gaps and separate policy-filled polygons.

    All inputs MUST share a projected CRS in metres. ``domain`` must already
    exclude the seaward side of the open boundary; never infer it from the mesh
    hull (which would hide missing bays). Policy fill is original water covered
    by ``policy_land`` but not by elements. No fill footprints are inferred from
    approximate centres in normalization logs.

    ``h_target_m`` caps the local size by the policy's coastal target: the
    nearest EXISTING triangle can be a coarse open-water one far from the gap,
    which would make an unmeshed inlet look unresolvable although the sizing
    field would have refined it (owner 2026-09-20).

    h is the median edge length of the nearest triangle (polygon distance),
    sampled along a raster medial axis, then median-aggregated per region.
    Width is twice the exact distance to the gap boundary along that axis.
    The median w/h suppresses tapered tips and tiny wide pockets; p10 and p90
    and the full sampled profile remain available. Rows use w/h conservatively
    (no equilateral-triangle packing bonus). Clear defects have median w/h >= 2
    and area >= 4*h**2; one-row regions are marginal. Accessibility is reported
    independently: a wide basin behind a narrow neck still contains missing water.

    Corridor tests erode ORIGINAL water minus policy fill by rows*h/2 plus a
    half-pixel diagonal, using the region median h. Four-connected raster paths
    must lead from actual covered water to ``corridor_fraction`` of the maximum
    straight-line distance from the gap/mesh interface. Thus a wide chamber
    behind a narrow neck cannot pass just because its median width is large.
    Isolated regions instead require one eroded component to contain that
    fraction of medial-axis samples. Length is total medial-axis network length;
    inland reach is straight-line distance to the mouth, not river chainage.

    Default pixels are min(10 m, nearest h/12). Oversized rasters and subpixel
    regions are explicitly unmeasured, never silently coarsened or cleared.
    """
    import shapely
    from scipy.ndimage import label
    from shapely.geometry import Polygon, box
    from shapely.ops import unary_union
    from shapely.strtree import STRtree
    from skimage.morphology import medial_axis

    if mesh.n_elements == 0:
        raise ValueError("coverage requires a nonempty mesh")
    if not np.isfinite(mesh.nodes).all():
        raise ValueError("mesh coordinates must be finite")
    if pixel_size_m is not None and (not np.isfinite(pixel_size_m) or pixel_size_m <= 0):
        raise ValueError("pixel_size_m must be finite and positive")
    if not np.isfinite(min_area_m2) or min_area_m2 < 0:
        raise ValueError("min_area_m2 must be finite and nonnegative")
    if max_raster_cells < 1 or not 0 < corridor_fraction <= 1:
        raise ValueError("invalid raster limit or corridor fraction")
    for geom in (original_land, domain, policy_land):
        if geom is not None and (not geom.is_valid or geom.geom_type not in
                                 ("Polygon", "MultiPolygon")):
            raise ValueError("land and domain must be valid polygonal geometries")
    if domain.is_empty:
        raise ValueError("domain must not be empty")
    triangles = np.array([Polygon(p) for p in mesh.nodes[mesh.elements]], dtype=object)
    if any(p.area <= 0 or not p.is_valid for p in triangles):
        raise ValueError("mesh contains degenerate triangles")
    tree = STRtree(triangles)
    vertices = mesh.nodes[mesh.elements]
    sizes = np.median(np.linalg.norm(vertices - np.roll(vertices, 1, axis=1), axis=2), axis=1)
    covered = unary_union(triangles).intersection(domain)
    original_water = domain.difference(original_land)
    omitted = original_water.difference(covered)
    filled = omitted.intersection(policy_land) if policy_land is not None else Polygon()
    gaps = omitted.difference(filled)
    water = original_water.difference(policy_land) if policy_land is not None else original_water
    records = []
    small_area = 0.0
    small_count = 0
    for poly in _polygons(gaps):
        if poly.area < min_area_m2:
            small_area += poly.area
            small_count += 1
            continue
        rec = _summary(poly)
        mouth = poly.boundary.intersection(covered.boundary)
        connected = mouth.length > 1e-6
        rec.update(touches_mesh=connected, mouth_length_m=float(mouth.length))
        h0 = float(sizes[tree.nearest(poly.representative_point())])
        ds = pixel_size_m if pixel_size_m is not None else min(10.0, h0 / 12)
        # Padding provides meshed-water seeds for the erosion test at the mouth.
        x0, y0, x1, y1 = poly.buffer(3 * h0).bounds
        nx, ny = int(np.ceil((x1 - x0) / ds)), int(np.ceil((y1 - y0) / ds))
        rec.update(pixel_size_m=float(ds), classification="unmeasured")
        if nx * ny > max_raster_cells:
            rec["measurement_note"] = f"raster needs {nx * ny} cells; limit {max_raster_cells}"
            records.append(rec)
            continue
        xx, yy = np.meshgrid(x0 + (np.arange(nx) + 0.5) * ds,
                             y0 + (np.arange(ny) + 0.5) * ds)
        mask = shapely.contains_xy(poly, xx, yy)
        axis = medial_axis(mask, rng=0)
        ay, ax = np.nonzero(axis)
        if not len(ax):
            rec["measurement_note"] = "gap is below raster resolution"
            records.append(rec)
            continue
        points = shapely.points(xx[axis], yy[axis])
        widths = 2 * shapely.distance(points, poly.boundary)
        nearest = tree.nearest(points)
        local_h = sizes[nearest]
        if h_target_m is not None:
            local_h = np.minimum(local_h, float(h_target_m))
        ratios = widths / local_h
        h = float(np.median(local_h))
        ratio = float(np.median(ratios))
        length = 0.0
        for dy, dx in ((0, 1), (1, -1), (1, 0), (1, 1)):
            by, bx = ay + dy, ax + dx
            valid = (by < ny) & (bx >= 0) & (bx < nx)
            length += np.count_nonzero(axis[by[valid], bx[valid]]) * ds * np.hypot(dy, dx)
        gap_points = shapely.points(xx[mask], yy[mask])
        distances = shapely.distance(gap_points, mouth) if connected else None
        reach = float(np.max(distances)) if connected else None
        # Clearance to real water banks, not the artificial gap/mesh mouth.
        wet = shapely.contains_xy(water, xx, yy)
        wet_points = shapely.points(xx[wet], yy[wet])
        clearance = np.zeros_like(xx)
        # Clip banks beyond the largest tested radius to keep distance queries local.
        # Artificial clip edges cannot affect either erosion threshold.
        margin = 2 * h + 2 * ds
        local_water = water.intersection(box(x0 - margin, y0 - margin,
                                             x0 + nx * ds + margin, y0 + ny * ds + margin))
        clearance[wet] = shapely.distance(wet_points, local_water.boundary)
        seeds = shapely.contains_xy(covered.intersection(water), xx, yy)
        corridors = {}
        for rows in (1, 2):
            labels, _ = label(wet & (clearance >= rows * h / 2 + ds / np.sqrt(2)))
            if connected:
                seed_labels = np.unique(labels[seeds])
                accessible = np.isin(labels[mask], seed_labels[seed_labels > 0])
                achieved = float(np.max(distances[accessible])) if accessible.any() else 0.0
                fraction = achieved / reach if reach else 0.0
            else:
                counts = np.bincount(labels[axis])
                fraction = float(counts[1:].max() / len(ax)) if len(counts) > 1 else 0.0
                achieved = None
            corridors[str(rows)] = {
                "feasible": bool(fraction >= corridor_fraction),
                "reach_fraction": float(fraction), "reach_m": achieved,
            }
        # A long inlet that is wide at its mouth and tapers inland has a small
        # MEDIAN w/h while still omitting water that the target size resolves.
        # Measure how much of the medial axis is wide enough for one and two
        # rows and classify on that, not on the median alone.
        step = ds  # medial-axis samples are one pixel apart
        wide_1h_m = float(np.count_nonzero(ratios >= 1.0) * step)
        wide_2h_m = float(np.count_nonzero(ratios >= 2.0) * step)
        classification = "subtarget"
        if ratio >= 1 or wide_1h_m >= 2 * h:
            classification = "marginal"
        if (ratio >= 2 and poly.area >= 4 * h**2) or wide_2h_m >= 2 * h:
            classification = "clear-defect"
        rec.update(
            classification=classification, h_m=h,
            h_range_m=[float(local_h.min()), float(local_h.max())],
            width_p10_m=float(np.quantile(widths, 0.1)),
            width_p50_m=float(np.median(widths)), width_p90_m=float(np.quantile(widths, 0.9)),
            w_h_p10=float(np.quantile(ratios, 0.1)), w_h_p50=ratio,
            w_h_p90=float(np.quantile(ratios, 0.9)), rows_p50=int(np.floor(ratio)),
            wide_1h_m=wide_1h_m, wide_2h_m=wide_2h_m,
            length_m=float(length), inland_reach_m=reach,
            distance_to_mesh_m=float(poly.distance(covered)), corridors=corridors,
            width_profile=[{"xy_m": [float(x), float(y)], "width_m": float(w),
                            "h_m": float(hh), "element_id": int(i) + 1}
                           for x, y, w, hh, i in zip(xx[axis], yy[axis], widths, local_h, nearest)],
        )
        records.append(rec)
    priority = {"clear-defect": 0, "marginal": 1, "unmeasured": 2, "subtarget": 3}
    records.sort(key=lambda r: (priority[r["classification"]], not r["touches_mesh"],
                               -r["area_m2"] * r.get("w_h_p50", 0)))
    for i, rec in enumerate(records, 1):
        rec["id"] = i
    policy = sorted((_summary(p) for p in _polygons(filled)), key=lambda p: -p["area_m2"])
    for i, rec in enumerate(policy, 1):
        rec.update(id=f"P{i}", classification="policy-filled")
    return {
        "method": {
            "units": "metres; area in square metres", "width_quantile": 0.5,
            "width_note": "median medial-axis width suppresses tapered ends; p10/p90 also reported",
            "h_note": "nearest triangle median edge, median over gap medial axis",
            "length_note": "total raster medial-axis network length, including branches",
            "inland_note": "max sampled straight-line distance to mesh interface; isolated=null",
            "corridor_note": "4-connected conservative erosion, constant region median h",
            "corridor_required_fraction": corridor_fraction,
            "clear_defect_rule": "median w/h >= 2 and area >= 4*h^2; accessibility separate",
            "min_area_m2": min_area_m2,
        },
        "regions": records, "policy_filled": policy,
        "clear_defects": sum(r["classification"] == "clear-defect" for r in records),
        "unmeasured": sum(r["classification"] == "unmeasured" for r in records),
        "below_area_threshold": {"count": small_count, "area_m2": float(small_area)},
        "original_water_area_m2": float(original_water.area),
        "omitted_water_area_m2": float(omitted.area),
        "policy_filled_area_m2": float(filled.area),
    }
