"""Refinement regions around small installations: the facility-scale policy.

The policy (owner, 2026-11-16, first for the Futtsu nori nets): a mesh that
already resolves a fishery at an intermediate size (e.g. 50 m) is refined to a
fine size (e.g. 5 m) only at each facility and its near field, and stays at
the intermediate size everywhere else, including the rest of the licensed
area. The near field of a facility is **as wide as the facility is long**: a
facility of length scale ``L`` is buffered by ``factor * L`` (``factor`` = 1
by default). Facilities whose near fields touch are merged, because two cores
closer than their transitions cannot be meshed apart anyway, and each merged
cluster becomes one region of a local-refinement recipe
(:mod:`fvcom_mesh_tools.refine`, ``docs/local_refine.md``).

Nothing here is specific to nori: any set of facility polygons (aquaculture
rafts, fish cages, outfall diffusers, piers) goes in, and the result is a
GeoJSON with one polygon per cluster plus the recipe that names them.

Length scale rules (``rule``):

``long_side``
    the long side of the minimum rotated rectangle (default; for an elongated
    net block, its length);
``sqrt_area``
    the square root of the area (an isotropic scale);
``diameter``
    the largest distance between two vertices.

Lengths are measured in a metric CRS: the UTM zone of the facilities' centroid
unless one is given.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

RULES = ("long_side", "sqrt_area", "diameter")


@dataclass
class Facility:
    fid: Any
    polygon_m: Any  # shapely Polygon in the metric CRS
    length_m: float
    buffer_m: float


@dataclass
class FacilityRegion:
    name: str
    polygon_m: Any
    members: list = field(default_factory=list)

    @property
    def area_m2(self) -> float:
        return float(self.polygon_m.area)


def utm_crs_for(lon: float, lat: float) -> str:
    """The EPSG code of the UTM zone that contains (lon, lat)."""
    if not (-180.0 <= lon <= 180.0 and -80.0 <= lat <= 84.0):
        raise ValueError(f"no UTM zone for lon={lon}, lat={lat}")
    zone = min(int((lon + 180.0) // 6.0) + 1, 60)
    return f"EPSG:{32600 + zone if lat >= 0 else 32700 + zone}"


def length_scale(polygon_m, rule: str = "long_side") -> float:
    """The length scale of one facility polygon (metric coordinates)."""
    if rule not in RULES:
        raise ValueError(f"rule must be one of {RULES}, got {rule!r}")
    if polygon_m.is_empty or polygon_m.area <= 0.0:
        raise ValueError("facility polygon is empty or has no area")
    if rule == "sqrt_area":
        return math.sqrt(polygon_m.area)
    if rule == "diameter":
        xy = np.asarray(polygon_m.exterior.coords)
        d = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
        return float(d.max())
    rect = polygon_m.minimum_rotated_rectangle
    xy = np.asarray(rect.exterior.coords)
    return float(np.hypot(*np.diff(xy, axis=0).T).max())


def read_facilities(path, where: dict | None = None, id_field: str | None = "id"):
    """``[(id, lon/lat Polygon), ...]`` from a GeoJSON file.

    MultiPolygon features contribute one facility per part (ids get a
    ``.<part>`` suffix). ``where`` keeps features whose properties equal every
    given value (compared as strings, so a CLI ``k=v`` filter works).
    """
    from shapely.geometry import shape

    data = json.loads(Path(path).read_text())
    crs = (data.get("crs") or {}).get("properties", {}).get("name", "")
    if crs and "4326" not in crs and "CRS84" not in crs:
        raise ValueError(f"{path}: GeoJSON must be lon/lat (WGS84), found crs {crs!r}")
    out = []
    for k, feat in enumerate(data.get("features", [])):
        props = feat.get("properties") or {}
        if where and any(str(props.get(key)) != str(val) for key, val in where.items()):
            continue
        geom = shape(feat["geometry"])
        fid = props.get(id_field, k) if id_field else k
        parts = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
        if any(p.geom_type != "Polygon" for p in parts):
            raise ValueError(f"{path}: feature {fid} is a {geom.geom_type}, not a polygon")
        for j, p in enumerate(parts):
            out.append((fid if len(parts) == 1 else f"{fid}.{j}", p))
    if not out:
        raise ValueError(f"{path}: no facility matches {where or '(no filter)'}")
    return out


def facility_regions(
    facilities: Sequence[tuple[Any, Any]],
    *,
    factor: float = 1.0,
    rule: str = "long_side",
    min_buffer_m: float = 0.0,
    crs: str | None = None,
    fill_holes: bool = True,
    prefix: str = "facilities",
    quad_segs: int = 8,
) -> tuple[list[FacilityRegion], list[Facility], str]:
    """Buffer each facility by ``max(factor * L, min_buffer_m)`` and merge.

    Returns the regions (largest first, named ``<prefix>_c01`` ...), the
    facilities with their length and buffer, and the metric CRS used.
    ``fill_holes`` closes water enclosed by a cluster's near fields: it is
    surrounded by fine cells on every side and a coarse island there would
    only add transitions.
    """
    import pyproj
    from shapely.geometry import Polygon
    from shapely.ops import transform, unary_union

    if not facilities:
        raise ValueError("no facilities")
    if not (factor >= 0.0 and math.isfinite(factor)):
        raise ValueError(f"factor must be a finite number >= 0, got {factor}")
    if not (min_buffer_m >= 0.0 and math.isfinite(min_buffer_m)):
        raise ValueError(f"min_buffer_m must be a finite number >= 0, got {min_buffer_m}")
    if factor == 0.0 and min_buffer_m == 0.0:
        raise ValueError("factor and min_buffer_m are both 0: the regions would be the "
                         "facilities themselves, with no near field")
    ids = [fid for fid, _ in facilities]
    if len(set(map(str, ids))) != len(ids):
        raise ValueError("facility ids must be unique")
    if crs is None:
        c = unary_union([p for _, p in facilities]).centroid
        crs = utm_crs_for(c.x, c.y)
    fwd = pyproj.Transformer.from_crs("EPSG:4326", crs, always_xy=True).transform

    facs = []
    for fid, poly in facilities:
        pm = transform(fwd, poly)
        if not pm.is_valid:
            raise ValueError(f"facility {fid}: invalid polygon")
        L = length_scale(pm, rule)
        facs.append(Facility(fid, pm, L, max(factor * L, min_buffer_m)))
    merged = unary_union([f.polygon_m.buffer(f.buffer_m, quad_segs=quad_segs) for f in facs])
    parts = list(merged.geoms) if merged.geom_type == "MultiPolygon" else [merged]
    if fill_holes:
        parts = [Polygon(p.exterior) for p in parts]
    parts.sort(key=lambda p: -p.area)
    width = max(2, len(str(len(parts))))
    regions = []
    for k, p in enumerate(parts, 1):
        members = [f.fid for f in facs if p.intersects(f.polygon_m)]
        regions.append(FacilityRegion(f"{prefix}_c{k:0{width}d}", p, members))
    return regions, facs, crs


def to_lonlat(geom_m, crs: str):
    import pyproj
    from shapely.ops import transform

    inv = pyproj.Transformer.from_crs(crs, "EPSG:4326", always_xy=True).transform
    return transform(inv, geom_m)


def summary(
    regions: Sequence[FacilityRegion],
    facs: Sequence[Facility],
    *,
    target_h_m: float,
    ambient_h_m: float | None = None,
    gradation: float | None = None,
    within_m=None,
) -> dict[str, Any]:
    """What the regions will cost, before any meshing.

    Core elements are estimated for equilateral cells (``area / (sqrt(3)/4 h^2)``).
    With ``ambient_h_m`` and ``gradation`` the transition width
    ``(ambient - target) / gradation`` is reported, and with ``within_m`` (a
    metric polygon, e.g. the licence) the share of core and of core +
    transition that falls outside it.
    """
    from shapely.ops import unary_union

    core = unary_union([r.polygon_m for r in regions])
    out: dict[str, Any] = {
        "facilities": len(facs),
        "length_m": {"min": min(f.length_m for f in facs),
                     "median": float(np.median([f.length_m for f in facs])),
                     "max": max(f.length_m for f in facs)},
        "regions": len(regions),
        "core_area_km2": core.area / 1e6,
        "core_elements_estimate": int(round(core.area / (math.sqrt(3) / 4 * target_h_m**2))),
        "per_region": [{"name": r.name, "area_km2": r.area_m2 / 1e6,
                        "facilities": [str(m) for m in r.members]} for r in regions],
    }
    if ambient_h_m is not None and gradation is not None:
        width = max(ambient_h_m - target_h_m, 0.0) / gradation
        hole = core.buffer(width)
        out["transition_width_m"] = width
        out["hole_area_km2"] = hole.area / 1e6
        if within_m is not None:
            out["core_outside_within_frac"] = core.difference(within_m).area / core.area
            out["hole_outside_within_frac"] = hole.difference(within_m).area / hole.area
    return out


def write_regions_geojson(path, regions: Sequence[FacilityRegion], crs: str,
                          meta: dict | None = None) -> None:
    from shapely.geometry import mapping

    feats = []
    for r in regions:
        feats.append({"type": "Feature",
                      # a plain string: list-valued properties do not survive
                      # every GeoJSON reader (refine reads this with geopandas)
                      "properties": {"region": r.name, "area_m2": round(r.area_m2, 1),
                                     "facilities": ",".join(str(m) for m in r.members)},
                      "geometry": mapping(to_lonlat(r.polygon_m, crs))})
    doc = {"type": "FeatureCollection", "metric_crs": crs, "features": feats}
    if meta:
        doc["generated_by"] = meta
    Path(path).write_text(json.dumps(doc, indent=1) + "\n")


def _rel(target, start_dir) -> str:
    return os.path.relpath(Path(target).resolve(), Path(start_dir).resolve())


def recipe_text(
    recipe_path,
    regions_geojson,
    regions: Iterable[FacilityRegion],
    *,
    base_mesh,
    target_h_m: float,
    gradation: float,
    dt_expected_s: float,
    base_depth=None,
    base_obc=None,
    header: str = "",
) -> str:
    """A ``recipes/refine`` YAML whose regions are the clusters, by file.

    Paths are written relative to the recipe, as :func:`refine.load_refine`
    resolves them. The coastline is kept (``preserve``): the base already
    carries the coastline its own refinement chose.
    """
    here = Path(recipe_path).resolve().parent
    lines = [f"# {ln}".rstrip() for ln in header.strip().splitlines()] if header else []
    lines += [f"base_mesh: {_rel(base_mesh, here)}"]
    if base_depth is not None:
        lines.append(f"base_depth: {_rel(base_depth, here)}")
    if base_obc is not None:
        lines.append(f"base_obc: {_rel(base_obc, here)}")
    lines += ["", f"dt_expected_s: {dt_expected_s}        # advisory: an alert, not a veto",
              f"gradation: {gradation}", "coastline: preserve", "", "refine:"]
    geo = _rel(regions_geojson, here)
    for r in regions:
        lines += [f"  - name: {r.name}",
                  "    geometry:",
                  f"      file: {geo}",
                  f"      where: {{region: {r.name}}}",
                  f"    target_h_m: {target_h_m:g}",
                  "    priority: 0"]
    return "\n".join(lines) + "\n"
