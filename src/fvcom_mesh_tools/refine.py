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
``dt = L / sqrt(g * H)`` with ``L`` the minimum altitude.  A 30 m target over
4 m of water allows 4.1 s where the present mesh allows 11.9 s.

**The time step does not veto a region** (owner 2026-09-22).  A fishery is
given -- its position and its required resolution are inputs, not preferences
-- so a recipe that cannot meet a time step is still the recipe.  What the
caller needs is to be *told*, loudly, before the run: ``dt_expected_s`` is the
step the caller was counting on, and :func:`preflight` raises an alert in its
report when the region will not deliver it.  Refusal is reserved for things
that make the operation impossible, not expensive.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from fvcom_mesh_tools.sizing import _geometry, _keys, _positive

__all__ = [
    "GRAVITY_M_S2",
    "ALTITUDE_OVER_EDGE",
    "depths_from_base",
    "RefineRegion",
    "frozen_changes",
    "hole_clearance",
    "load_refine",
    "preflight",
    "transition_width_m",
]

GRAVITY_M_S2 = 9.81

#: The reported time step uses the MINIMUM ALTITUDE of a triangle, not an edge
#: (notebook 392, and coast_fit guards both).  For an equilateral triangle the
#: altitude is sqrt(3)/2 of the edge, so a target edge length buys only 0.866
#: of the dt a naive edge estimate suggests -- which is how a 30 m target over
#: 4.15 m of water came to be advertised as 4.70 s when the honest figure is
#: 4.07 s, below its own 4.5 s floor.
ALTITUDE_OVER_EDGE = float(np.sqrt(3.0) / 2.0)


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
        dt_expected_s: 4.5      # advisory: an alert, not a veto
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
    _keys(cfg, ["base_mesh", "dt_expected_s", "gradation", "refine"])
    base = Path(cfg["base_mesh"])
    if not base.is_absolute():
        base = (Path(path).resolve().parent / base).resolve()
    if not base.exists():
        raise ValueError(f"base_mesh not found: {base}")
    cfg["base_mesh"] = base
    cfg["dt_expected_s"] = _positive(cfg["dt_expected_s"], "dt_expected_s")
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
    dt_expected_s: float,
    ambient_h_m: float,
    depth_of,
    land=None,
    samples: int = 4000,
) -> dict[str, Any]:
    """Decide whether a region can be met, before any meshing happens.

    ``depth_of(lon, lat)`` returns positive-down depths in metres (the same
    bathymetry the mesh will carry, floors already applied). ``land`` is the
    land polygon in lon/lat; it is **required** when ``touch_coast`` is false,
    and omitting it is an error rather than a silent pass.

    This checks the CORE only. The hole the generator actually cuts is the
    core plus its transition, which is several times larger, so a core clear
    of land says nothing about the hole -- the first Futtsu recipe passed here
    while its 2,239 m hole reached a coastline 1,014 m from the centre.
    :func:`hole_clearance` is the test for that, and needs the mesh.

    Returns a report. ``ValueError`` is reserved for a region that cannot be
    built at all -- an empty or wholly dry geometry, a transition too short for
    the gradation, a missing land polygon. A time step below ``dt_expected_s``
    sets ``dt_alert`` in the report instead: the region is still built, and the
    caller is told what it will cost.
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

    if land is None and not region.touch_coast:
        raise ValueError(
            f"{region.name}: touch_coast is false, so a land polygon is required "
            "to check it")
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
    c = np.sqrt(GRAVITY_M_S2 * depth.max())
    dt_edge = region.target_h_m / c
    dt = ALTITUDE_OVER_EDGE * dt_edge          # the measure that is reported
    alert = None
    if dt < dt_expected_s:
        would_need = dt_expected_s * c / ALTITUDE_OVER_EDGE
        alert = (
            f"{region.name}: target {region.target_h_m:g} m over {depth.max():.2f} m of "
            f"water allows dt = {dt:.2f} s by minimum altitude ({dt_edge:.2f} s by "
            f"shortest edge), against the {dt_expected_s:g} s expected -- the run will "
            f"cost {dt_expected_s / dt:.1f}x the external steps. Keeping {dt_expected_s:g} s "
            f"would need a {would_need:.0f} m target. This is an equilateral upper bound: "
            "a legal 30-30-120 cell halves the altitude and halves the step again.")

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
        "dt_by_shortest_edge_s": float(dt_edge),
        "dt_measure": "minimum altitude of an equilateral cell / sqrt(g*Hmax)",
        "dt_expected_s": float(dt_expected_s),
        "dt_alert": alert,
        "dt_step_cost_factor": float(dt_expected_s / dt) if dt > 0 else float("inf"),
        "elements_core": float(n_core),
        "elements_transition": float(n_trans),
        "elements_replaced": float(n_was),
        "elements_added": float(n_core + n_trans - n_was),
    }


def depths_from_base(base_nodes, base_elements, base_depths, new_nodes):
    """Depths for the patched mesh, taken from the base mesh's own field.

    The bathymetry is **not** refined (owner 2026-09-22). The base mesh is the
    topography actually being simulated; a refined region has to sit on the
    same seabed, so new nodes get the base field evaluated at their position
    and every retained node keeps its exact value. Refining the bathymetry
    too is a separate question for later.

    This also settles what looked like the hardest open problem. Re-sampling a
    survey inside the patch would have produced a mixed-source depth field and
    an r-factor constraint problem across the seam -- a mutable node beside a
    frozen depth ``H`` must satisfy ``H/1.5 <= h <= 1.5H``, and two frozen
    neighbours can make that infeasible. Interpolation cannot: a value inside
    a base element lies between that element's own vertex depths, so for two
    points in the same element

        |h_a - h_b| / (h_a + h_b)  <=  (h_max - h_min) / (h_max + h_min)

    which is the r of the base edge joining that element's deepest and
    shallowest vertices. **Refinement can only leave the r-factor where it was
    or improve it**, never worsen it.

    Linear interpolation over the base triangulation; points outside it fall
    back to the nearest base node, and the count is reported.
    """
    from scipy.spatial import cKDTree

    xy = np.asarray(base_nodes, dtype=float)[:, :2]
    tri = np.asarray(base_elements, dtype=np.int64)
    dep = np.asarray(base_depths, dtype=float)
    if dep.shape[0] != xy.shape[0]:
        raise ValueError("base depths must be one per base node")
    if not np.isfinite(dep).all():
        raise ValueError("base depths must be finite")
    new = np.asarray(new_nodes, dtype=float)[:, :2]

    from matplotlib.tri import LinearTriInterpolator, Triangulation

    mtri = Triangulation(xy[:, 0], xy[:, 1], tri)
    out = np.asarray(LinearTriInterpolator(mtri, dep)(new[:, 0], new[:, 1]))
    outside = ~np.isfinite(out)
    if outside.any():
        out = np.array(out, dtype=float)
        out[outside] = dep[cKDTree(xy).query(new[outside])[1]]
    return np.asarray(out, dtype=float), int(outside.sum())


def hole_clearance(nodes, elements, region, *, transition_m, open_boundaries=()):
    """What the hole would actually cut, on the real mesh.

    :func:`preflight` judges the core; the generator cuts core + transition,
    which is several times larger. This selects the elements whose centroid
    falls inside that footprint and reports what they touch -- the physical
    boundary, the open boundary, and how far the selection reaches beyond the
    analytic envelope, which it always does because whole triangles are taken.

    ``nodes`` are in the mesh CRS (metres) and ``region`` is the geometry in
    the same CRS. Returns a report; the caller decides what is acceptable.
    """
    import shapely

    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    footprint = region.buffer(float(transition_m))
    centroid = xy[tri].mean(axis=1)
    inside = shapely.contains(footprint, shapely.points(centroid[:, 0], centroid[:, 1]))
    sel = tri[inside]
    if not len(sel):
        return {"n_selected": 0, "reaches_boundary": False, "reaches_open_boundary": False}

    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    mesh_boundary = set(map(tuple, u[c == 1].tolist()))
    es = np.sort(np.vstack([sel[:, [0, 1]], sel[:, [1, 2]], sel[:, [2, 0]]]), axis=1)
    us, cs = np.unique(es, axis=0, return_counts=True)
    rim = us[cs == 1]
    physical = [tuple(x) for x in rim.tolist() if tuple(x) in mesh_boundary]
    interface = [tuple(x) for x in rim.tolist() if tuple(x) not in mesh_boundary]
    obc = set()
    for seg in open_boundaries:
        obc.update(np.asarray(seg, dtype=np.int64).ravel().tolist())
    touched = set(np.unique(sel).tolist())
    lengths = np.linalg.norm(xy[rim[:, 0]] - xy[rim[:, 1]], axis=1) if len(rim) else np.zeros(0)
    reach = shapely.distance(
        shapely.points(xy[np.unique(sel), 0], xy[np.unique(sel), 1]), region.centroid)
    return {
        "n_selected": int(inside.sum()),
        "n_rim_edges": int(len(rim)),
        "n_interface_edges": int(len(interface)),
        "n_physical_boundary_edges": int(len(physical)),
        "reaches_boundary": bool(physical),
        "reaches_open_boundary": bool(touched & obc),
        "n_open_boundary_nodes": int(len(touched & obc)),
        "rim_edge_min_m": float(lengths.min()) if len(lengths) else 0.0,
        "rim_edge_median_m": float(np.median(lengths)) if len(lengths) else 0.0,
        "rim_edge_max_m": float(lengths.max()) if len(lengths) else 0.0,
        "selection_reach_m": float(reach.max()) if len(reach) else 0.0,
        "requested_reach_m": float(shapely.distance(region.centroid, region.boundary)
                                   + transition_m),
    }


def frozen_changes(base_nodes, new_nodes, affected_mask, tol_m: float = 1e-6) -> dict[str, Any]:
    """How far the mesh moved where it promised not to.

    ``affected_mask`` marks the nodes inside core + transition. Every other
    node must keep its coordinates; the count of those that did not is the
    number this contract is judged by.

    This is a COORDINATE check on arrays in row correspondence, which only
    holds while the node numbering is unchanged. A patch inserts and deletes
    nodes and renumbers them, so the generator must supply an old-to-new node
    map and this check must be applied through it -- and it is only one of the
    invariants the contract needs: retained connectivity and orientation,
    depths, boundary membership and order must be compared too.
    """
    a = np.asarray(base_nodes, dtype=float)[:, :2]
    b = np.asarray(new_nodes, dtype=float)[:, :2]
    mask = np.asarray(affected_mask, dtype=bool)
    if a.shape != b.shape or mask.shape[0] != a.shape[0]:
        raise ValueError("frozen check needs matching node arrays and a mask over them")
    if not np.isfinite(b).all():
        raise ValueError("the new coordinates contain non-finite values")
    moved = np.linalg.norm(b - a, axis=1)
    frozen = ~mask
    # NaN fails every comparison, so `moved > tol` would pass a node whose
    # coordinate was destroyed; the finiteness check above is what catches it.
    bad = frozen & (moved > tol_m)
    return {
        "n_frozen": int(frozen.sum()),
        "n_moved_in_frozen": int(bad.sum()),
        "max_move_in_frozen_m": float(moved[frozen].max()) if frozen.any() else 0.0,
        "ok": bool(not bad.any()),
    }
