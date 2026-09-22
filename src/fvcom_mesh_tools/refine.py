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
    "COASTLINE_MODES",
    "GRAVITY_M_S2",
    "ALTITUDE_OVER_EDGE",
    "depths_from_base",
    "limit_rfactor",
    "RefineRegion",
    "frozen_changes",
    "hole_clearance",
    "load_refine",
    "preflight",
    "transition_width_m",
]

#: How the coastline inside the hole is treated.  The hole reaches the coast
#: whenever the transition does, which is normal and not a reason to refuse:
#: the coastline is simply meshed at the target size like any other boundary.
#: What differs is where the new boundary nodes are placed.
#:
#: ``preserve``  on the existing segments only.  The polyline is geometrically
#:               identical -- subdividing a segment does not move it.
#: ``resample``  along the SOURCE shoreline at the target size.  This is the
#:               default because it is the only one that improves fidelity:
#:               the base polyline is 300-600 m between nodes and sits up to
#:               68.5 m from the real coast at its segment midpoints, and a
#:               30 m resample recovers that for free where we are refining
#:               anyway.
#: ``spline``    a smooth resample of the base polyline, for when there is no
#:               usable source shoreline.  It eases a sharp corner but adds no
#:               information.
#:
#: Every mode is bounded by ``coastline_tolerance_m`` against the base
#: polyline, and none of them touches the coastline outside the hole.
COASTLINE_MODES = ("preserve", "resample", "spline")

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
              ["transition_m", "priority"])
        name = spec["name"]
        if not isinstance(name, str) or not name.strip():
            raise ValueError("region name must be a nonempty string")
        self.name = name
        self.geometry = _geometry(spec["geometry"])
        # What the recipe DECLARED, kept alongside the polygon.  A consumer
        # that has to project the region needs this: a bbox's four corners
        # are all about equidistant from its centre, so guessing "circle"
        # from the spread of vertex radii turns every near-square box into a
        # disc (review finding 14, 2026-09-22).
        geom_spec = spec["geometry"]
        self.kind = ("circle" if "circle" in geom_spec
                     else "bbox" if "bbox" in geom_spec
                     else "file" if "file" in geom_spec else "polygon")
        self.circle = (
            (float(geom_spec["circle"]["center"][0]),
             float(geom_spec["circle"]["center"][1]),
             float(geom_spec["circle"]["radius_m"]))
            if self.kind == "circle" else None)
        # Where a polygon came from, so the report can say it. A fishery
        # boundary read from a file is the case this exists for: "region
        # futtsu_nori" is not enough to reproduce a run, and the file, the
        # filter and the row are.
        self.source = None
        if self.kind == "file":
            from fvcom_mesh_tools.sizing import _geometry_from_file

            _, self.source = _geometry_from_file(geom_spec)
        self.target_h_m = _positive(spec["target_h_m"], "target_h_m")
        self.transition_m = (
            _positive(spec["transition_m"], "transition_m") if "transition_m" in spec else None
        )
        priority = spec.get("priority", 0)
        if isinstance(priority, bool) or not isinstance(priority, (int, float)) \
                or not np.isfinite(priority):
            raise ValueError("priority must be finite numeric")
        self.priority = float(priority)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"<RefineRegion {self.name} target={self.target_h_m:g} m>"


def load_refine(path) -> dict[str, Any]:
    """Read and strictly validate a YAML refinement recipe.

    Schema (unknown keys are errors)::

        base_mesh: .../TokyoBay_grd.dat    # or a fort.14
        base_depth: .../TokyoBay_dep_m7001tp_rfac0p2_cap300.dat
        base_obc: .../TokyoBay_obc.dat
        dt_expected_s: 4.5            # advisory: an alert, not a veto
        gradation: 0.165
        coastline: preserve           # preserve | resample | spline
        rfactor_limit: base           # base | off | a number in (0, 1)
        coastline_tolerance_m: 100    # max departure from the base polyline
        refine:
          - name: futtsu_nori
            geometry: {circle: {center: [139.7881, 35.3228], radius_m: 300}}
            target_h_m: 30
            priority: 0

    ``geometry`` is also a GeoJSON ``Polygon``, a ``bbox``, or a polygon read
    from a file -- ``{file: fishery.geojson, where: {...}, index: 0,
    buffer_m: 25}`` -- which is how a real fishery boundary arrives.
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
    _keys(cfg, ["base_mesh", "dt_expected_s", "gradation", "refine"],
          ["base_depth", "base_obc", "coastline", "coastline_tolerance_m",
           "rfactor_limit"])

    def _resolve(key, required=True):
        if key not in cfg or cfg[key] is None:
            if required:
                raise ValueError(f"{key} is required")
            return None
        # `~` is how a recipe names a file in the user's checkout of the model
        # repository, which is where a production base lives.
        q = Path(cfg[key]).expanduser()
        if not q.is_absolute():
            q = (Path(path).resolve().parent / q).resolve()
        if not q.exists():
            raise ValueError(f"{key} not found: {q}")
        return q

    base = _resolve("base_mesh")
    cfg["base_mesh"] = base
    # A finished FVCOM case is three files, and the depth file is the one that
    # matters: `_grd.dat` carries a depth column from whenever it was made,
    # and the baseline names a different one by tag.  On goto2023 node 1 that
    # is 4.31 m in the grd against 7.16 m in the b12 baseline's
    # `TokyoBay_dep_m7001tp_rfac0p2_cap300.dat`.  So a .dat base must say
    # which depths it means; a fort.14 carries its own and must not.
    cfg["base_depth"] = _resolve("base_depth", required=base.suffix == ".dat")
    cfg["base_obc"] = _resolve("base_obc", required=False)
    if base.suffix != ".dat" and cfg["base_depth"] is not None:
        raise ValueError("base_depth applies to an FVCOM _grd.dat base; a "
                         "fort.14 carries its own depths")
    cfg["dt_expected_s"] = _positive(cfg["dt_expected_s"], "dt_expected_s")
    cfg["gradation"] = _positive(cfg["gradation"], "gradation")
    # `preserve` is the default: the operation refines an existing mesh and
    # leaves its coastline where it is (owner, 2026-09-22).  `resample` is
    # for the case where the base polyline is known to be a poor rendering of
    # a source shoreline that is available -- it buys fidelity and costs the
    # guarantee that the coastline did not move.
    mode = cfg.setdefault("coastline", "preserve")
    if mode not in COASTLINE_MODES:
        raise ValueError(f"coastline must be one of {COASTLINE_MODES}, got {mode!r}")
    cfg["coastline_tolerance_m"] = _positive(
        cfg.get("coastline_tolerance_m", 100.0), "coastline_tolerance_m")
    # The base is a product with a property -- m7001tp_rfac0p2_cap300 means
    # r <= 0.2 on every edge -- and interpolation does not inherit it. `base`
    # asks for whatever the base itself achieves, a number asks for that, and
    # `off` accepts the new edges as interpolation leaves them.
    rl = cfg.setdefault("rfactor_limit", "base")
    if rl not in ("base", "off"):
        cfg["rfactor_limit"] = _positive(rl, "rfactor_limit")
        if not 0.0 < cfg["rfactor_limit"] < 1.0:
            raise ValueError("rfactor_limit must be in (0, 1), 'base' or 'off'")
    if not isinstance(cfg["refine"], list) or not cfg["refine"]:
        raise ValueError("refine must be a nonempty list")
    # A relative geometry file resolves against the RECIPE, like base_mesh,
    # not against whatever directory the run happens to start in. A recipe
    # that only works from the repository root is not a recipe.
    here = Path(path).resolve().parent
    for r in cfg["refine"]:
        g = r.get("geometry")
        if isinstance(g, dict) and "file" in g:
            q = Path(str(g["file"])).expanduser()
            g["file"] = str(q if q.is_absolute() else (here / q).resolve())
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
    bathymetry the mesh will carry). ``land`` is the land polygon in lon/lat;
    when given, the report says how much of the core is dry.

    This looks at the CORE only. The hole the generator cuts is the core plus
    its transition, several times larger, and it commonly reaches the coast --
    the Futtsu core clears land by 714 m while its 2,239 m hole overlaps the
    coastline by 1,225 m. That is normal: the coastline is meshed at the target
    size like any other boundary, under the recipe's ``coastline`` mode.
    :func:`hole_clearance` measures what the hole actually touches.

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

    on_land = np.zeros(plon.size, dtype=bool)
    land_alert = None
    if land is not None:
        on_land = shapely.contains(land, shapely.points(plon, plat))
        if on_land.any():
            # The region is given, so a dry patch of it is news, not grounds
            # for refusal -- unless there is no water at all to mesh.
            land_alert = (
                f"{region.name}: {int(on_land.sum())} of {plon.size} core samples "
                f"({100 * on_land.mean():.0f} %) are on land; the core will be meshed "
                "only where there is water")
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
            f"would need a {would_need:.0f} m target. This is an equilateral upper "
            "bound at Cr = 1 with no velocity allowance, and it is sampled over the "
            "CORE only: a legal 30-30-120 cell of the same side has 1/sqrt(3) of the "
            "altitude, and the transition can be deeper than anything sampled here. "
            "Measure the achieved step on the finished mesh.")

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
        "land_alert": land_alert,
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
    shallowest vertices.

    That bound holds **only between two points of the same base element**.
    Refinement changes the connectivity, so a new edge can join points from
    different base elements and is not bounded by it: on a 3x2 grid with
    column depths 2, 3 and 4.5 every base edge has r <= 0.2, and the new edge
    between the interpolated 2.1 and 4.35 has r = 0.349 (second review,
    finding 8).  The earlier revisions of this docstring claimed the bound
    globally; they were wrong.  If the r-factor is a requirement, measure it
    on the finished mesh -- the QA battery does not gate it.

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
    # Interpolation cannot leave the convex hull of its inputs, so anything
    # outside the base range is rounding, and rounding below the base minimum
    # is not harmless: the base mesh sits exactly on the 2 m depth floor, so
    # 1.9999999999999998 is a QA failure that used to be hidden by an
    # 11-digit write format.
    return np.clip(np.asarray(out, dtype=float), dep.min(), dep.max()), \
        int(outside.sum())


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


def limit_rfactor(elements, depths, movable, rmax: float, *,
                  depth_min: float | None = None, depth_max: float | None = None,
                  rounds: int = 200):
    """Bring new nodes' depths inside the base mesh's own r-factor limit.

    The base is a product with a property: ``m7001tp_rfac0p2_cap300`` means
    every edge of the goto2023 mesh satisfies ``|hi-hj|/(hi+hj) <= 0.2``, and
    inheriting its bathymetry ought to inherit that too. Interpolation does
    not deliver it. The bound proved in :func:`depths_from_base` holds between
    two points of the SAME base element, and refinement changes the
    connectivity: measured on the Futtsu patch, ten new edges exceed 0.2 and
    the worst reaches 0.3075, while every wholly frozen edge stays at 0.2.

    Only ``movable`` depths change, so **no base depth moves**: a frozen node
    keeps the value the model runs with. What moves is the value this code
    chose for a node the base never had, and the report says by how much --
    that number is the price of the property, and the caller should look at
    it rather than trust it.

    ``r <= rmax`` is exactly ``max/min <= (1+rmax)/(1-rmax)``. Each violating
    edge is pulled to that ratio: the one movable end onto the bound, or both
    ends symmetrically about their geometric mean. That is a Gauss-Seidel
    sweep, not a solve, and it can fail -- a movable node between two frozen
    depths more than ``R`` apart has an empty feasible set, and the report
    says so instead of pretending otherwise.
    """
    tri = np.asarray(elements, dtype=np.int64)
    h = np.array(depths, dtype=float)
    h0 = h.copy()
    free = np.asarray(movable, dtype=bool)
    if not 0.0 <= rmax < 1.0:
        raise ValueError("rmax must be in [0, 1)")
    # rmax = 0 is what a constant-depth base asks for -- "no jump at all" --
    # and `rfactor_limit: base` passes the base's own worst r, so a flat
    # bottom reached this with 0 and raised (fourth review).
    ratio = (1.0 + rmax) / (1.0 - rmax)
    lo = float(depth_min) if depth_min is not None else -np.inf
    hi_cap = float(depth_max) if depth_max is not None else np.inf

    e = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]],
                                     tri[:, [2, 0]]]), axis=1), axis=0)
    n_rounds = 0
    for n_rounds in range(1, rounds + 1):
        r = np.abs(h[e[:, 0]] - h[e[:, 1]]) / (h[e[:, 0]] + h[e[:, 1]])
        bad = np.flatnonzero(r > rmax + 1e-12)
        bad = bad[free[e[bad]].any(axis=1)]
        if not bad.size:
            break
        for k in bad:
            i, j = int(e[k, 0]), int(e[k, 1])
            if h[i] < h[j]:
                i, j = j, i                      # i is the deeper end
            if free[i] and free[j]:
                g = float(np.sqrt(h[i] * h[j]))
                h[i] = np.clip(g * np.sqrt(ratio), lo, hi_cap)
                h[j] = np.clip(g / np.sqrt(ratio), lo, hi_cap)
            elif free[i]:
                h[i] = np.clip(h[j] * ratio, lo, hi_cap)
            elif free[j]:
                h[j] = np.clip(h[i] / ratio, lo, hi_cap)
    r = np.abs(h[e[:, 0]] - h[e[:, 1]]) / (h[e[:, 0]] + h[e[:, 1]])
    touched = free[e].any(axis=1)
    over = r > rmax + 1e-9
    moved = np.abs(h - h0)
    # `converged` is about EVERY edge, not only the ones this could move. An
    # edge between two frozen nodes is unfixable here, and a patch can create
    # one that the base never had -- new connectivity joining two retained
    # nodes that were not neighbours. Reporting that as converged because
    # nothing could be done about it is how it would reach a mesh unnoticed
    # (fourth review). The caller decides whether a frozen-pair violation is
    # inherited; it has the base to compare against and this does not.
    return h, {
        "rmax": float(rmax),
        "rounds": int(n_rounds),
        "converged": bool(not over.any()),
        "n_edges_over_rmax": int(over.sum()),
        "n_over_rmax_movable": int((over & touched).sum()),
        "n_over_rmax_frozen_pair": int((over & ~touched).sum()),
        "frozen_pair_edges_over_rmax": e[over & ~touched].tolist(),
        "max_r_on_new_edges": float(r[touched].max()) if touched.any() else 0.0,
        "max_r_frozen": float(r[~touched].max()) if (~touched).any() else 0.0,
        "n_depths_changed": int((moved > 1e-9).sum()),
        "max_depth_change_m": float(moved.max()),
        "max_frozen_depth_change_m": float(moved[~free].max())
        if (~free).any() else 0.0,
    }
