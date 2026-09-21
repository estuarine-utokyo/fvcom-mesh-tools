"""Local refinement of an existing mesh: specification, and the checks that
must pass before any meshing is attempted.

The sizing recipe (``sizing.py``) declares regions that change the SIZING
FIELD, which means the mesh is rebuilt everywhere -- DistMesh is a global
relaxation, so nothing outside the region survives unchanged.  A refinement
recipe says something different: **keep this mesh, and re-cut only this
region and the ground it needs to blend over**.

Three zones, and the contract for each:

``core``
    inside the declared geometry; achieved edge length ``target_h_m``.
``transition``
    the annulus the gradation needs in order to reach the ambient size;
    remeshed, so it changes.
``frozen``
    everything else; node coordinates and connectivity **identical** to the
    base mesh.  :func:`frozen_changes` checks it rather than asserting it.

The transition width is not a free choice.  With a linear size gradation ``g``
the edge length grows as ``h(r) = h_in + g * (r - r_core)``, so reaching the
ambient size takes

    W = (h_ambient - h_target) / g

-- 1,939 m for a 30 m target in a 350 m ambient field at g = 0.165, which is
six times the radius of a 300 m fishery.  Declare the target and let the width
follow; an explicit ``transition_m`` is checked against this and rejected when
it is too small to be reached at the recipe's gradation.

The cost that bites is the time step.  FVCOM integrates with one global
external step, so the smallest element anywhere sets it for the whole run:
``dt = h / sqrt(g * H)``.  A 30 m target over 4 m of water allows 4.8 s where
the present mesh allows 11.9 s.  ``dt_floor_s`` is therefore mandatory, and
:func:`preflight` refuses a region that would break it -- before spending the
compute, not after.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from fvcom_mesh_tools.sizing import _geometry, _keys, _positive

__all__ = [
    "GRAVITY_M_S2",
    "RefineRegion",
    "frozen_changes",
    "load_refine",
    "preflight",
    "transition_width_m",
]

GRAVITY_M_S2 = 9.81


def transition_width_m(target_h_m: float, ambient_h_m: float, gradation: float) -> float:
    """Distance a linear gradation needs to climb from ``target`` to ``ambient``."""
    if ambient_h_m <= target_h_m:
        return 0.0
    return (float(ambient_h_m) - float(target_h_m)) / float(gradation)


class RefineRegion:
    """One declared refinement: geometry in lon/lat, plus its target and rules."""

    def __init__(self, spec: dict[str, Any]):
        _keys(spec, ["name", "geometry", "target_h_m"],
              ["transition_m", "priority", "touch_coast"])
        name = spec["name"]
        if not isinstance(name, str) or not name.strip():
            raise ValueError("region name must be a nonempty string")
        self.name = name
        self.geometry = _geometry(spec["geometry"])
        self.target_h_m = _positive(spec["target_h_m"], "target_h_m")
        self.transition_m = (
            _positive(spec["transition_m"], "transition_m") if "transition_m" in spec else None
        )
        priority = spec.get("priority", 0)
        if isinstance(priority, bool) or not isinstance(priority, (int, float)) \
                or not np.isfinite(priority):
            raise ValueError("priority must be finite numeric")
        self.priority = float(priority)
        touch = spec.get("touch_coast", False)
        if not isinstance(touch, bool):
            raise ValueError("touch_coast must be true or false")
        self.touch_coast = touch

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"<RefineRegion {self.name} target={self.target_h_m:g} m>"


def load_refine(path) -> dict[str, Any]:
    """Read and strictly validate a YAML refinement recipe.

    Schema (unknown keys are errors)::

        base_mesh: outputs/.../sample_repro_final.14
        dt_floor_s: 4.5
        gradation: 0.165
        refine:
          - name: futtsu_nori
            geometry: {circle: {center: [139.7881, 35.3228], radius_m: 300}}
            target_h_m: 30
            touch_coast: false
            priority: 0
    """
    import yaml

    class UniqueLoader(yaml.SafeLoader):
        pass

    def mapping(loader, node):
        result: dict[str, Any] = {}
        for key_node, value_node in node.value:
            key = loader.construct_object(key_node, deep=True)
            if not isinstance(key, str) or key in result:
                raise ValueError("recipe keys must be unique strings")
            result[key] = loader.construct_object(value_node, deep=True)
        return result

    UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, mapping)
    with Path(path).resolve().open() as stream:
        cfg = yaml.load(stream, Loader=UniqueLoader)
    _keys(cfg, ["base_mesh", "dt_floor_s", "gradation", "refine"])
    base = Path(cfg["base_mesh"])
    if not base.is_absolute():
        base = (Path(path).resolve().parent / base).resolve()
    if not base.exists():
        raise ValueError(f"base_mesh not found: {base}")
    cfg["base_mesh"] = base
    cfg["dt_floor_s"] = _positive(cfg["dt_floor_s"], "dt_floor_s")
    cfg["gradation"] = _positive(cfg["gradation"], "gradation")
    if not isinstance(cfg["refine"], list) or not cfg["refine"]:
        raise ValueError("refine must be a nonempty list")
    regions = [RefineRegion(r) for r in cfg["refine"]]
    if len({r.name for r in regions}) != len(regions):
        raise ValueError("region names must be unique")
    cfg["refine"] = regions
    return cfg


def _to_metres(geom, lat0: float):
    """Project a lon/lat geometry to local metres about the equator meridian."""
    from shapely.affinity import affine_transform

    cos = float(np.cos(np.radians(lat0)))
    return affine_transform(geom, [111000.0 * cos, 0.0, 0.0, 111000.0, 0.0, 0.0])


def preflight(
    region: RefineRegion,
    *,
    gradation: float,
    dt_floor_s: float,
    ambient_h_m: float,
    depth_of,
    land=None,
    samples: int = 4000,
) -> dict[str, Any]:
    """Decide whether a region can be met, before any meshing happens.

    ``depth_of(lon, lat)`` returns positive-down depths in metres (the same
    bathymetry the mesh will carry, floors already applied). ``land`` is the
    land polygon in lon/lat; it is required when ``touch_coast`` is false.

    Returns a report. Raises ``ValueError`` when the region cannot be built as
    declared -- a failure here costs a second, a failure after meshing costs
    the run.
    """
    import shapely

    geom = region.geometry
    lat0 = float(geom.centroid.y)
    needed = transition_width_m(region.target_h_m, ambient_h_m, gradation)
    if region.transition_m is not None and region.transition_m + 1e-9 < needed:
        raise ValueError(
            f"{region.name}: transition_m={region.transition_m:g} m is shorter than the "
            f"{needed:.0f} m that gradation {gradation:g} needs to climb from "
            f"{region.target_h_m:g} m to {ambient_h_m:g} m; raise it, raise the target, "
            "or raise the gradation")
    width = region.transition_m if region.transition_m is not None else needed

    # Sample the core on a lattice dense enough that a 300 m disc gets
    # thousands of points; the deepest sample sets the time step.
    w, s, e, n = geom.bounds
    k = int(np.ceil(np.sqrt(samples)))
    gx, gy = np.meshgrid(np.linspace(w, e, k), np.linspace(s, n, k))
    inside = shapely.contains(geom, shapely.points(gx.ravel(), gy.ravel()))
    if not inside.any():
        raise ValueError(f"{region.name}: geometry contains no sample points")
    plon, plat = gx.ravel()[inside], gy.ravel()[inside]

    on_land = np.zeros(plon.size, dtype=bool)
    if land is not None:
        on_land = shapely.contains(land, shapely.points(plon, plat))
        if not region.touch_coast and on_land.any():
            raise ValueError(
                f"{region.name}: {int(on_land.sum())} of {plon.size} core samples are on "
                "land and touch_coast is false; move the region, shrink it, or set "
                "touch_coast: true to remesh the coastline inside it")
    wet = ~on_land
    if not wet.any():
        raise ValueError(f"{region.name}: the core is entirely on land")

    depth = np.asarray(depth_of(plon[wet], plat[wet]), dtype=float)
    if not np.isfinite(depth).all() or (depth <= 0).any():
        raise ValueError(f"{region.name}: depths must be finite and positive-down")
    dt = region.target_h_m / np.sqrt(GRAVITY_M_S2 * depth.max())
    if dt < dt_floor_s:
        allowed = dt_floor_s * np.sqrt(GRAVITY_M_S2 * depth.max())
        raise ValueError(
            f"{region.name}: target {region.target_h_m:g} m over {depth.max():.2f} m of "
            f"water allows dt = {dt:.2f} s, below the floor {dt_floor_s:g} s; the "
            f"smallest target this water permits is {allowed:.0f} m")

    area = _to_metres(geom, lat0).area
    outer = _to_metres(geom.buffer(width / 111000.0), lat0).area
    per_elem = np.sqrt(3.0) / 4.0
    n_core = area / (per_elem * region.target_h_m ** 2)
    # The transition's element count, integrated over the annulus with a
    # linearly growing size, is well approximated by its area over the mean
    # of the size at both ends squared.
    h_mid = 0.5 * (region.target_h_m + ambient_h_m)
    n_trans = (outer - area) / (per_elem * h_mid ** 2)
    n_was = outer / (per_elem * ambient_h_m ** 2)
    return {
        "name": region.name,
        "target_h_m": region.target_h_m,
        "ambient_h_m": float(ambient_h_m),
        "gradation": float(gradation),
        "transition_m": float(width),
        "transition_required_m": float(needed),
        "core_area_m2": float(area),
        "affected_area_m2": float(outer),
        "core_samples": int(plon.size),
        "core_on_land": int(on_land.sum()),
        "core_depth_min_m": float(depth.min()),
        "core_depth_max_m": float(depth.max()),
        "dt_s": float(dt),
        "dt_floor_s": float(dt_floor_s),
        "elements_core": float(n_core),
        "elements_transition": float(n_trans),
        "elements_replaced": float(n_was),
        "elements_added": float(n_core + n_trans - n_was),
    }


def frozen_changes(base_nodes, new_nodes, affected_mask, tol_m: float = 1e-6) -> dict[str, Any]:
    """How far the mesh moved where it promised not to.

    ``affected_mask`` marks the nodes inside core + transition. Every other
    node must keep its coordinates; the count of those that did not is the
    number this contract is judged by.
    """
    a = np.asarray(base_nodes, dtype=float)[:, :2]
    b = np.asarray(new_nodes, dtype=float)[:, :2]
    mask = np.asarray(affected_mask, dtype=bool)
    if a.shape != b.shape or mask.shape[0] != a.shape[0]:
        raise ValueError("frozen check needs matching node arrays and a mask over them")
    moved = np.linalg.norm(b - a, axis=1)
    frozen = ~mask
    bad = frozen & (moved > tol_m)
    return {
        "n_frozen": int(frozen.sum()),
        "n_moved_in_frozen": int(bad.sum()),
        "max_move_in_frozen_m": float(moved[frozen].max()) if frozen.any() else 0.0,
        "ok": bool(not bad.any()),
    }
