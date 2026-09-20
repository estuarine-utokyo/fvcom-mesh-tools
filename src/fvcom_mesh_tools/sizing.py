"""Declarative regional sizing, in metres on a geographic lon/lat lattice.

YAML schema (unknown keys are errors)::

    coastal_target_m: 290
    max_edge_length_m: 1400
    gradation: 0.165
    cfl: {dt_s: 15, cr: 0.45}
    hmin_m: 10                 # optional absolute lower bound, NOT coastal target
    regions:
      - name: port
        geometry: {bbox: [139.7, 35.4, 139.8, 35.5]}
        target_h_m: 30
        transition_m: 500     # optional linear blend outside geometry
        priority: 0           # higher wins; ties choose smaller target

Geometry also accepts a GeoJSON Polygon, including holes. Coordinates are
EPSG:4326; sizes are never converted from degrees here. Distances/areas use a
local equirectangular projection (111000 m/degree, cosine of mean latitude).
Targets replace the ambient field within their footprint, including coarsening.
The final eight-neighbour shortest-path lower envelope enforces
|h_i-h_j| <= gradation * distance(i,j). It can lower a regional target further.
A scalar hmin is enforced; CFL is diagnostic only. dt = cr*h/sqrt(9.81*depth)
(and dt at Cr=1 is also reported). Dry-only/unsampled statistics are null.
Requires the existing pipeline YAML and vector Shapely extras; no GPL imports.
"""
from __future__ import annotations

import heapq
from pathlib import Path

import numpy as np


def _positive(value, name, *, zero=False):
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise ValueError(f"{name} must be numeric")
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if not np.isfinite(value) or (value < 0 if zero else value <= 0):
        raise ValueError(f"invalid {name}: {value}")
    return value


def _keys(obj, required, optional=()):
    if not isinstance(obj, dict) or not set(required) <= obj.keys():
        raise ValueError(f"expected mapping with {required}")
    if obj.keys() - set(required) - set(optional):
        raise ValueError(f"unknown keys: {obj.keys() - set(required) - set(optional)}")


def _geometry(spec):
    from shapely.geometry import box, shape

    if not isinstance(spec, dict):
        raise ValueError("geometry must be a mapping")
    if "bbox" in spec:
        _keys(spec, ["bbox"])
        b = np.asarray(spec["bbox"], dtype=float)
        if b.shape != (4,) or not np.isfinite(b).all() or np.any(b[:2] >= b[2:]):
            raise ValueError("bbox must be [west, south, east, north]")
        geom = box(*b)
    else:
        _keys(spec, ["type", "coordinates"])
        if spec["type"] != "Polygon":
            raise ValueError("only GeoJSON Polygon is supported")
        try:
            rings = spec["coordinates"]
            if not rings or any(len(ring) < 4 or list(ring[0]) != list(ring[-1])
                                for ring in rings):
                raise ValueError("polygon rings must be closed with at least four positions")
            if any(np.asarray(ring).ndim != 2 or np.asarray(ring).shape[1] != 2
                   or not np.isfinite(np.asarray(ring, float)).all() for ring in rings):
                raise ValueError("polygon coordinates must be finite lon/lat pairs")
            geom = shape(spec)
        except Exception as exc:
            raise ValueError("invalid polygon") from exc
    if geom.is_empty or not geom.is_valid or geom.area <= 0:
        raise ValueError("geometry must be a valid nonempty polygon")
    w, s, e, n = geom.bounds
    if w < -180 or e > 180 or s <= -90 or n >= 90 or e-w > 180:
        raise ValueError("geometry must be local geographic coordinates, without dateline crossing")
    return geom


def _regions(regions):
    if not isinstance(regions, list):
        raise ValueError("regions must be a list")
    names = set()
    for r in regions:
        _keys(r, ["name", "geometry", "target_h_m"], ["transition_m", "priority"])
        if not isinstance(r["name"], str) or not r["name"].strip() or r["name"] in names:
            raise ValueError("region names must be nonempty and unique")
        names.add(r["name"])
        _geometry(r["geometry"])
        _positive(r["target_h_m"], "target_h_m")
        _positive(r.get("transition_m", 0), "transition_m", zero=True)
        if isinstance(r.get("priority", 0), bool) or not isinstance(
            r.get("priority", 0), (float, int)
        ) or not np.isfinite(
            r.get("priority", 0)
        ):
            raise ValueError("priority must be finite numeric")


def load_sizing(path):
    """Read and strictly validate a YAML sizing recipe; paths resolve absolutely."""
    import yaml

    class UniqueLoader(yaml.SafeLoader):
        pass

    def mapping(loader, node):
        result = {}
        for key_node, value_node in node.value:
            key = loader.construct_object(key_node, deep=True)
            if not isinstance(key, str) or key in result:
                raise ValueError("recipe keys must be unique strings")
            result[key] = loader.construct_object(value_node, deep=True)
        return result

    UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, mapping)
    with Path(path).resolve().open() as stream:
        cfg = yaml.load(stream, Loader=UniqueLoader)
    _keys(cfg, ["coastal_target_m", "max_edge_length_m", "gradation", "cfl", "regions"],
          ["hmin_m", "one_wide"])
    from fvcom_mesh_tools.one_wide import parse_one_wide
    if "one_wide" in cfg:
        cfg["one_wide"] = parse_one_wide(cfg["one_wide"])
    for key in ("coastal_target_m", "max_edge_length_m", "gradation", "hmin_m"):
        if key in cfg:
            cfg[key] = _positive(cfg[key], key)
    if cfg["coastal_target_m"] > cfg["max_edge_length_m"]:
        raise ValueError("coastal target exceeds max edge length")
    if cfg.get("hmin_m", 0) > cfg["max_edge_length_m"]:
        raise ValueError("hmin exceeds max edge length")
    _keys(cfg["cfl"], ["dt_s", "cr"])
    for key in ("dt_s", "cr"):
        cfg["cfl"][key] = _positive(cfg["cfl"][key], key)
    _regions(cfg["regions"])
    return cfg


def _limit(values, x, y, grade):
    """Multi-source Dijkstra lower envelope on the eight-neighbour lattice."""
    out = values.copy()
    queue = [(float(v), i, j) for (i, j), v in np.ndenumerate(out)]
    heapq.heapify(queue)
    while queue:
        value, i, j = heapq.heappop(queue)
        if value > out[i, j]:
            continue
        for di, dj in ((-1, -1), (-1, 0), (-1, 1), (0, -1),
                       (0, 1), (1, -1), (1, 0), (1, 1)):
            a, b = i+di, j+dj
            if 0 <= a < out.shape[0] and 0 <= b < out.shape[1]:
                bound = value + grade * np.hypot(x[a, b]-x[i, j], y[a, b]-y[i, j])
                if bound < out[a, b]:
                    out[a, b] = bound
                    heapq.heappush(queue, (bound, a, b))
    return out


def apply_sizing_regions(values_m, lon, lat, regions, *, gradation, hmin_m=None, cfl=None):
    """Return a new field and JSON-compatible diagnostics; never mutate inputs.

    lon/lat must be same-shape 2-D lattice arrays. Empty regions still enforce
    gradation; callers requiring an exact no-op must skip this function.
    """
    import shapely
    from shapely.affinity import scale, translate

    values = np.asarray(values_m, dtype=float)
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    if (values.ndim != 2 or not values.size or lon.shape != values.shape
            or lat.shape != values.shape or not np.isfinite([values, lon, lat]).all()
            or np.any(values <= 0) or np.any(abs(lat) >= 90) or np.any(abs(lon) > 180)):
        raise ValueError("expected finite positive 2-D sizes and matching geographic arrays")
    grade = _positive(gradation, "gradation", zero=True)
    lower = 0 if hmin_m is None else _positive(hmin_m, "hmin_m")
    _regions(regions)
    speed = None
    if cfl is not None:
        _keys(cfl, ["dt_s", "cr", "depth_m"])
        dt, cr = _positive(cfl["dt_s"], "dt_s"), _positive(cfl["cr"], "cr")
        depth = np.asarray(cfl["depth_m"], float)
        if depth.shape != values.shape or not np.isfinite(depth).all() or np.any(depth < 0):
            raise ValueError("depth_m must be finite, nonnegative and match values")
        speed = np.sqrt(9.81*depth)
    sx = 111000*np.cos(np.deg2rad(lat.mean()))
    origin_x, origin_y = float(lon.mean()), float(lat.mean())
    x, y = (lon-origin_x)*sx, (lat-origin_y)*111000
    points = shapely.points(x, y)
    out = values.copy()
    masks = {}
    areas = {}
    # Weakest first: priority wins, then smallest target. Name ensures stable ties.
    ordered = sorted(regions, key=lambda r: (r.get("priority", 0), -r["target_h_m"], r["name"]))
    for r in ordered:
        poly = scale(translate(_geometry(r["geometry"]), -origin_x, -origin_y),
                     xfact=sx, yfact=111000, origin=(0, 0))
        dist = shapely.distance(points, poly)
        inside = shapely.covers(poly, points)
        masks[r["name"]], areas[r["name"]] = inside, float(poly.area)
        transition = r.get("transition_m", 0)
        weight = np.maximum(0, 1-dist/transition) if transition else inside.astype(float)
        candidate = values*(1-weight) + r["target_h_m"]*weight
        out = np.where(weight > 0, candidate, out)
    out = _limit(np.maximum(out, lower), x, y, grade)

    def minimum(a):
        return float(a.min()) if a.size else None

    def timestep(mask, cr_value):
        wet = mask & (speed > 0)
        return minimum(cr_value*out[wet]/speed[wet])

    report = {"regions": [], "implied_dt_s": None, "implied_dt_cr1_s": None}
    for r in regions:
        mask = masks[r["name"]]
        achieved = out[mask]
        item = {"name": r["name"], "area_m2": areas[r["name"]],
                "sample_count": int(mask.sum()), "requested_target_m": r["target_h_m"],
                "achieved_min_m": minimum(achieved),
                "achieved_median_m": float(np.median(achieved)) if achieved.size else None,
                "hmin_yield_fraction": float(r["target_h_m"] < lower) if achieved.size else None,
                "cfl_floor_min_m": None, "cfl_floor_max_m": None,
                "target_below_cfl_fraction": None, "implied_dt_s": None,
                "implied_dt_cr1_s": None, "requested_dt_s": None}
        if speed is not None:
            floor = speed[mask]*dt/cr
            item.update(cfl_floor_min_m=minimum(floor),
                        cfl_floor_max_m=float(floor.max()) if floor.size else None,
                        target_below_cfl_fraction=float(np.mean(r["target_h_m"] < floor))
                        if floor.size else None,
                        implied_dt_s=timestep(mask, cr), implied_dt_cr1_s=timestep(mask, 1),
                        requested_dt_s=minimum(cr*r["target_h_m"]/speed[mask & (speed > 0)]))
        report["regions"].append(item)
    if speed is not None:
        report["implied_dt_s"] = timestep(np.ones(values.shape, bool), cr)
        report["implied_dt_cr1_s"] = timestep(np.ones(values.shape, bool), 1)
    return out, report
