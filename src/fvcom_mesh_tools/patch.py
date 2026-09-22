"""The patch generator: cut a hole in an existing mesh and fill it finer.

``refine.py`` declares a refinement and refuses the impossible ones before any
meshing happens.  This module does the mechanical work that follows, and it is
deliberately free of any meshing library: it *selects* the faces to remove,
*extracts* the rim the filler must honour, *stitches* the filled patch back in,
and *verifies* that everything outside the patch is bit-for-bit what it was.
The fill itself -- DistMesh with the rim as constrained edges -- lives in a
notebook, because the only implementation available is GPL (see
``docs/local_refine.md`` and the licence policy in ``CLAUDE.md``).

The contract that makes this worth doing is the frozen zone.  A sizing-region
rebuild produces *a* mesh with a fine fishery in it; this produces *the same
mesh* with a fine fishery in it, so a run against the base mesh and a run
against the patched mesh differ by the patch and nothing else.  Every function
here exists to keep that claim checkable rather than asserted:
:func:`select_patch` records exactly which faces it takes,
:func:`stitch_patch` carries an old-to-new node map, and
:func:`verify_patch` compares coordinates, connectivity, depths and boundary
membership *through* that map.

Two node classes run through the whole module and are worth naming once:

``frozen``
    used by at least one retained element.  Its coordinate, its depth and the
    retained faces around it are part of the contract, so it cannot move and
    cannot be deleted.  Every interface rim node is frozen by definition.
``free``
    inside the hole, or on the hole's stretch of the physical boundary with no
    retained element left to hold it.  These may be moved, replaced or dropped
    -- which is what lets the coastline inside the hole be re-meshed at the
    target size while the coastline outside it is untouched.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

__all__ = [
    "PatchSelection",
    "introduced_violations",
    "boundary_after_patch",
    "refresh_depths",
    "ambient_size_field",
    "base_size_field",
    "boundary_rings",
    "coastline_points",
    "hole_polygon",
    "improve_patch",
    "effective_gradation",
    "field_gradation",
    "patch_sizing",
    "region_conflicts",
    "region_resolution",
    "rim_constraints",
    "select_patch",
    "stitch_patch",
    "verify_patch",
]


# --------------------------------------------------------------------------
# selection
# --------------------------------------------------------------------------


@dataclass
class PatchSelection:
    """Which faces the patch replaces, and the rim it must honour.

    All node ids are indices into the BASE mesh.  ``rings`` are closed walks
    of base node ids around each connected piece of the hole (first id is not
    repeated at the end); ``ring_is_hole`` marks the walks that bound a piece
    of retained mesh left standing inside the cut.
    """

    removed: np.ndarray
    retained: np.ndarray
    frozen_nodes: np.ndarray
    free_nodes: np.ndarray
    rim_edges: np.ndarray
    physical_rim: np.ndarray
    rings: list[np.ndarray]
    ring_is_hole: list[bool]
    report: dict[str, Any] = field(default_factory=dict)

    @property
    def n_removed(self) -> int:
        return int(self.removed.sum())


def _edge_table(tri: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unique undirected edges of a triangulation and how many faces use each."""
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    return np.unique(e, axis=0, return_counts=True)


def _isolated_faces(tri: np.ndarray, removed: np.ndarray) -> np.ndarray:
    """Retained faces with no retained edge-neighbour: spikes on the rim."""
    keep = np.flatnonzero(~removed)
    if keep.size < 2:
        # One face has no edge-neighbour because there is nothing left to be
        # a neighbour, not because it is a spike.  Absorbing it turns a valid
        # if tiny result into "the footprint removes the entire mesh"
        # (second review, finding 7).
        return np.zeros(len(tri), dtype=bool)
    sub = tri[keep]
    e = np.sort(np.vstack([sub[:, [0, 1]], sub[:, [1, 2]], sub[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(sub)), 3)
    order = np.lexsort((e[:, 1], e[:, 0]))
    e, owner = e[order], owner[order]
    k = np.flatnonzero(np.all(e[:-1] == e[1:], axis=1))
    has = np.zeros(len(sub), dtype=bool)
    has[owner[k]] = True
    has[owner[k + 1]] = True
    out = np.zeros(len(tri), dtype=bool)
    out[keep[~has]] = True
    return out


def _rim_degree(rim: np.ndarray, n_nodes: int) -> np.ndarray:
    deg = np.zeros(n_nodes, dtype=np.int64)
    if len(rim):
        np.add.at(deg, rim.ravel(), 1)
    return deg


def select_patch(
    nodes,
    elements,
    footprint,
    *,
    open_boundary_nodes=(),
    obc_guard_m: float = 0.0,
    max_repair_rounds: int = 50,
) -> PatchSelection:
    """Choose the faces to remove, and make the resulting hole meshable.

    ``footprint`` is the core-plus-transition geometry in the mesh CRS.  An
    element is taken when its centroid lies inside it, which means **the cut
    is not the footprint**: whole triangles go, so the hole reaches past the
    envelope -- 2,615 m against a requested 2,239 m on the Futtsu case.  The
    report carries what was actually taken; the circle in the recipe is not
    the thing the contract talks about.

    Centroid selection alone can leave a hole that is not a surface with a
    boundary: a vertex where the removed faces form two separate fans has
    four rim edges, and no closed walk exists through it.  The repair is to
    keep taking every face at such a vertex until each rim node has exactly
    two rim edges.  This only ever *grows* the hole, so it terminates, and
    the growth is reported rather than hidden -- it can push the cut towards
    the open boundary, which is checked afterwards, not before.

    Raises ``ValueError`` when the cut reaches the open boundary (with its
    guard band), empties the mesh, or disconnects the retained mesh.
    """
    import shapely

    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    n_nodes = xy.shape[0]
    centroid = xy[tri].mean(axis=1)
    removed = np.asarray(
        shapely.contains(footprint, shapely.points(centroid[:, 0], centroid[:, 1])),
        dtype=bool,
    )
    if not removed.any():
        raise ValueError("the footprint selects no elements; nothing to refine")

    grown = 0
    resolved = False
    for _ in range(max_repair_rounds):
        sel = tri[removed]
        u, c = _edge_table(sel)
        rim = u[c == 1]
        deg = _rim_degree(rim, n_nodes)
        bad = np.where(deg > 2)[0]
        # A retained face whose three neighbours are all taken hangs off the
        # mesh by its vertices alone.  That is the same defect as a pinch seen
        # from the other side, and the same repair fixes it: take the spike
        # too.  A cut at (400, 200) on the 9x9 test grid leaves exactly one
        # (review finding 7, 2026-09-22).
        spike = _isolated_faces(tri, removed)
        if not bad.size and not spike.any():
            resolved = True
            break
        add = (np.isin(tri, bad).any(axis=1) & ~removed) | spike
        if not add.any():
            raise ValueError(
                f"the cut pinches at {bad.size} vertices and taking their faces "
                "does not resolve it; the footprint straddles a one-element-wide "
                "neck")
        removed |= add
        grown += int(add.sum())
    if not resolved:
        # One more look before giving up: the last round may have fixed it.
        sel = tri[removed]
        u, c = _edge_table(sel)
        if (_rim_degree(u[c == 1], n_nodes) > 2).any() \
                or _isolated_faces(tri, removed).any():
            raise ValueError("the cut could not be made manifold within "
                             f"{max_repair_rounds} rounds")

    if removed.all():
        raise ValueError("the footprint removes the entire mesh")

    retained = tri[~removed]
    frozen = np.unique(retained)
    sel = tri[removed]
    u, c = _edge_table(sel)
    rim = u[c == 1]
    all_u, all_c = _edge_table(tri)
    mesh_boundary = {tuple(x) for x in all_u[all_c == 1].tolist()}
    is_physical = np.array([tuple(x) in mesh_boundary for x in rim.tolist()], dtype=bool)

    obc = np.unique(np.asarray(list(open_boundary_nodes), dtype=np.int64)) \
        if len(open_boundary_nodes) else np.empty(0, dtype=np.int64)
    if obc.size:
        taken = np.unique(sel)
        hit = np.intersect1d(taken, obc)
        if hit.size:
            raise ValueError(
                f"the cut reaches the open boundary ({hit.size} OBC nodes removed); "
                "the open boundary is an input and must not be re-meshed")
        if obc_guard_m > 0:
            from scipy.spatial import cKDTree

            d = cKDTree(xy[taken]).query(xy[obc])[0]
            if float(d.min()) < obc_guard_m:
                raise ValueError(
                    f"the cut comes within {d.min():.0f} m of the open boundary, "
                    f"inside the {obc_guard_m:g} m guard band; an element merely "
                    "incident to an OBC node can change its orthogonality even "
                    "though no OBC coordinate moves")

    _require_connected(retained, frozen)

    rings, is_hole = boundary_rings(xy, rim)
    taken_nodes = np.unique(sel)
    free = np.setdiff1d(taken_nodes, frozen)
    reach = np.linalg.norm(xy[taken_nodes] - np.asarray(footprint.centroid.coords)[0],
                           axis=1)
    lengths = np.linalg.norm(xy[rim[:, 0]] - xy[rim[:, 1]], axis=1)
    report = {
        "n_elements_removed": int(removed.sum()),
        "n_elements_retained": int(len(retained)),
        "n_grown_by_repair": grown,
        "n_rim_edges": int(len(rim)),
        "n_physical_rim_edges": int(is_physical.sum()),
        "n_interface_rim_edges": int((~is_physical).sum()),
        "n_frozen_rim_nodes": int(np.intersect1d(np.unique(rim), frozen).size),
        "n_free_nodes": int(free.size),
        "n_rings": len(rings),
        "n_interior_islands": int(sum(is_hole)),
        "rim_edge_min_m": float(lengths.min()) if len(lengths) else 0.0,
        "rim_edge_median_m": float(np.median(lengths)) if len(lengths) else 0.0,
        "rim_edge_max_m": float(lengths.max()) if len(lengths) else 0.0,
        "selection_reach_m": float(reach.max()) if reach.size else 0.0,
    }
    return PatchSelection(removed=removed, retained=retained, frozen_nodes=frozen,
                          free_nodes=free, rim_edges=rim, physical_rim=is_physical,
                          rings=rings, ring_is_hole=is_hole, report=report)


def _require_connected(retained: np.ndarray, frozen: np.ndarray) -> None:
    """Refuse a cut that severs a piece of the mesh from the rest.

    The graph is over FACES joined by shared EDGES, not over nodes joined by
    edges.  A node graph calls two triangles that meet at a single point
    connected, and they are not: that is a non-manifold pinch, FVCOM will not
    run on it, and the promise made to the caller is edge connectivity.  A
    four-triangle fan with its two middle faces removed is the smallest case
    that a node graph waves through (review finding 7, 2026-09-22).
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    if not len(retained):
        raise ValueError("the cut leaves no elements")
    e = np.sort(np.vstack([retained[:, [0, 1]], retained[:, [1, 2]],
                           retained[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(retained)), 3)
    order = np.lexsort((e[:, 1], e[:, 0]))
    e, owner = e[order], owner[order]
    k = np.flatnonzero(np.all(e[:-1] == e[1:], axis=1))
    n = len(retained)
    g = coo_matrix((np.ones(len(k)), (owner[k], owner[k + 1])), shape=(n, n))
    ncomp, _ = connected_components(g, directed=False)
    if ncomp != 1:
        raise ValueError(
            f"the cut splits the retained mesh into {ncomp} pieces; "
            "a refinement may not disconnect the domain")
    del frozen


def boundary_rings(xy, rim_edges) -> tuple[list[np.ndarray], list[bool]]:
    """Order a set of rim edges into closed walks, and say which are islands.

    Every rim node has exactly two rim edges (:func:`select_patch` enforces
    it), so the walk is unambiguous.  A walk is flagged as an island when it
    winds the opposite way to the walk that contains it -- that is a piece of
    the base mesh left standing inside the cut, and the filler must treat it
    as a hole in its domain rather than as more water to mesh.
    """
    import shapely

    rim = np.asarray(rim_edges, dtype=np.int64)
    if not len(rim):
        return [], []
    nbr: dict[int, list[int]] = {}
    for a, b in rim.tolist():
        nbr.setdefault(a, []).append(b)
        nbr.setdefault(b, []).append(a)
    bad = [v for v, w in nbr.items() if len(w) != 2]
    if bad:
        raise ValueError(f"{len(bad)} rim nodes do not have exactly two rim edges")

    seen: set[int] = set()
    rings: list[np.ndarray] = []
    for start in nbr:
        if start in seen:
            continue
        walk = [start]
        seen.add(start)
        prev, cur = None, start
        while True:
            a, b = nbr[cur]
            nxt = a if a != prev else b
            if nxt == start:
                break
            walk.append(nxt)
            seen.add(nxt)
            prev, cur = cur, nxt
        rings.append(np.asarray(walk, dtype=np.int64))

    pts = np.asarray(xy, dtype=float)[:, :2]
    polys = [shapely.Polygon(pts[r]) for r in rings]
    areas = np.array([p.area for p in polys])
    # Parents first (descending area) so a ring's parent is already decided,
    # but the parent is looked up smallest-container-first (ascending): with
    # three levels of nesting the parity only comes out right against the
    # IMMEDIATE parent.
    order = np.argsort(areas)
    is_hole = [False] * len(rings)
    for k in order[::-1]:
        for j in order:
            if j == k or areas[j] <= areas[k]:
                continue
            if shapely.contains(polys[j], polys[k]):
                is_hole[k] = not is_hole[j]
                break
    return rings, is_hole


# --------------------------------------------------------------------------
# the rim the filler must honour
# --------------------------------------------------------------------------


def coastline_points(
    pts: np.ndarray,
    size,
    *,
    mode: str = "resample",
    shoreline=None,
    tolerance_m: float = 100.0,
) -> np.ndarray:
    """New boundary points along one free stretch of the coastline.

    ``pts`` is the base polyline for the stretch, endpoints first and last and
    **both frozen** -- they are where the stretch meets retained mesh, so they
    are returned unchanged whatever the mode.  Only the interior is replaced.

    ``size`` is the ACHIEVED edge length wanted, and it is a function of
    position, not the region's target.  Spacing the whole stretch at the
    target is the mistake that looks harmless and is not: the Futtsu hole
    reaches 2.4 km from a 300 m core, so a 30 m coastline would sit in a
    sizing field asking for 400 m triangles, and DistMesh resolves that
    contradiction with slivers -- measured, a minimum angle of 0.8 degrees and
    an element quality of 0.002 against a 30-degree gate.

    ``preserve`` subdivides the base segments, which cannot move the polyline
    at all.  ``resample`` walks the SOURCE shoreline between the projections of
    the two endpoints, which is the only mode that adds information: the base
    polyline is 300-600 m between nodes and sits up to 68.5 m off the real
    coast at its segment midpoints.  ``spline`` fits a smooth curve through the
    base nodes, for when no source shoreline is available; it eases a sharp
    corner but invents the easing.

    Every mode is checked against ``tolerance_m``: the result may not depart
    from the base polyline by more than that, because the recipe's promise is
    a finer coastline in the same place, not a different coastline.
    """
    import shapely

    pts = np.asarray(pts, dtype=float)[:, :2]
    if len(pts) < 2:
        raise ValueError("a coastline stretch needs at least two points")
    if mode not in ("preserve", "resample", "spline"):
        raise ValueError(f"unknown coastline mode {mode!r}")
    base = shapely.LineString(pts)

    if mode == "resample":
        if shoreline is None:
            raise ValueError("coastline: resample needs the source shoreline")
        out = _resample_on_source(pts, shoreline, size)
    elif mode == "spline":
        out = _spline_resample(pts, size)
    else:
        out = _subdivide(pts, size)

    off = shapely.distance(shapely.points(out[1:-1]), base) if len(out) > 2 \
        else np.zeros(0)
    if off.size and float(off.max()) > tolerance_m:
        raise ValueError(
            f"coastline: {mode} departs {off.max():.0f} m from the base polyline, "
            f"beyond the {tolerance_m:g} m tolerance; raise coastline_tolerance_m "
            "or use coastline: preserve")
    out[0], out[-1] = pts[0], pts[-1]
    return out


def _size_at(size, q: np.ndarray) -> np.ndarray:
    """Evaluate a size that may be a constant or a field."""
    q = np.atleast_2d(np.asarray(q, dtype=float))[:, :2]
    if callable(size):
        h = np.asarray(size(q), dtype=float).ravel()
    else:
        h = np.full(len(q), float(size))
    if not np.isfinite(h).all() or (h <= 0).any():
        raise ValueError("the coastline size field must be finite and positive")
    return h


def _walk(path: np.ndarray, size) -> np.ndarray:
    """Place points along a polyline at the LOCAL size, both ends kept.

    The step is re-read at every point, so a stretch that runs from beside the
    core out to the ambient field is fine at one end and coarse at the other
    -- which is the only spacing DistMesh can then honour without slivers.
    The last interval is stretched or dropped rather than left short: a stub
    next to a 400 m neighbour is the same defect in miniature.
    """
    path = np.asarray(path, dtype=float)[:, :2]
    seg = np.linalg.norm(np.diff(path, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(s[-1])
    if total <= 0:
        return path[[0, -1]].copy()

    def at(u):
        return np.column_stack([np.interp(u, s, path[:, 0]),
                                np.interp(u, s, path[:, 1])])

    stations = [0.0]
    u = 0.0
    while True:
        h = float(_size_at(size, at(np.array([u])))[0])
        nxt = u + h
        if nxt >= total - 0.5 * h:
            break
        stations.append(nxt)
        u = nxt
    stations.append(total)
    return at(np.asarray(stations))


def _subdivide(pts: np.ndarray, size) -> np.ndarray:
    """Insert points on the existing segments at about the local size.

    Every original vertex is KEPT and only interior points are added, so the
    polyline is geometrically identical -- subdividing a segment is not a
    displacement.  Re-walking the stretch by arc length instead would drop
    vertices wherever the local size is coarser than the original spacing,
    and on the Futtsu stretch that quietly straightened 18 base nodes into 13
    and moved the coastline 74.5 m.  That is the opposite of what
    ``preserve`` promises.
    """
    pts = np.asarray(pts, dtype=float)[:, :2]
    out = [pts[0]]
    for a, b in zip(pts[:-1], pts[1:]):
        if float(np.linalg.norm(b - a)) <= 0:
            continue
        # walk each ORIGINAL segment separately: the step still follows the
        # size field, but a segment's two ends are always emitted
        out.extend(_walk(np.vstack([a, b]), size)[1:])
    return np.asarray(out, dtype=float)


def _spline_resample(pts: np.ndarray, size) -> np.ndarray:
    """A cubic parametric spline through the base nodes, sampled at the local size."""
    from scipy.interpolate import CubicSpline

    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    if s[-1] <= 0:
        return pts.copy()
    cs = CubicSpline(s, pts, axis=0)
    dense = np.asarray(cs(np.linspace(0.0, s[-1], max(200, 20 * len(pts)))), dtype=float)
    return _walk(dense, size)


def coastline_curve(pts: np.ndarray, mode: str, shoreline=None) -> np.ndarray:
    """The curve one stretch of coastline was cut from.

    Whatever the mode, the new boundary nodes lie on this curve, so it is also
    the only curve they may later be slid along.  The base mesh's own
    ``land_boundaries`` list is NOT that curve: on ``sample_repro_final.14``
    only 820 of its 939 consecutive pairs are boundary edges and it jumps up
    to 2,867 m, so a line built from it crosses open water and a node
    projected onto it lands in the sea.
    """
    pts = np.asarray(pts, dtype=float)[:, :2]
    if mode == "resample":
        piece = _source_substring(pts, shoreline)
        if piece is not None:
            return piece

    if mode == "spline":
        from scipy.interpolate import CubicSpline

        t = np.concatenate([[0.0], np.cumsum(
            np.linalg.norm(np.diff(pts, axis=0), axis=1))])
        if t[-1] > 0:
            cs = CubicSpline(t, pts, axis=0)
            return np.asarray(cs(np.linspace(0.0, t[-1], max(200, 20 * len(pts)))))
    return pts


def _source_substring(pts: np.ndarray, shoreline):
    """The piece of the source shoreline between a stretch's two endpoints.

    ``shoreline`` is one LineString, or several to choose from.  When it is
    several, the one nearest THIS stretch is used -- not the one nearest the
    free rim as a whole, which is a closest-pair distance and can hand every
    stretch the ring that only one of them lies on (review finding 11,
    2026-09-22).  ``None`` when there is no usable piece.
    """
    import shapely
    from shapely.ops import substring

    if shoreline is None:
        return None
    lines = _slide_lines(shoreline)
    if not lines:
        return None
    here = shapely.LineString(pts)
    # By the fit along the WHOLE stretch, not by the closest pair.  A
    # candidate that merely touches one endpoint and then departs wins a
    # min-distance contest: a line running diagonally away from (0,0) beat
    # one that stays 0.1 from the entire base stretch (second review,
    # finding 3).
    probe = shapely.points(np.column_stack([
        np.interp(np.linspace(0, 1, 32), np.linspace(0, 1, len(pts)), pts[:, 0]),
        np.interp(np.linspace(0, 1, 32), np.linspace(0, 1, len(pts)), pts[:, 1])]))
    line = min(lines, key=lambda ln: float(shapely.distance(probe, ln).max()))
    del here
    s0 = line.project(shapely.Point(pts[0]))
    s1 = line.project(shapely.Point(pts[-1]))
    if s0 == s1:
        return None
    piece = substring(line, s0, s1)
    coords = np.asarray(piece.coords, dtype=float)[:, :2]
    if len(coords) < 2:
        return None
    # Walking the wrong way round a closed ring gives a piece far longer than
    # the stretch it replaces; fall back rather than swap in half a coastline.
    if piece.length > 3.0 * shapely.LineString(pts).length:
        return None
    return coords


def _resample_on_source(pts: np.ndarray, shoreline, size,
                        simplify_frac: float = 0.25) -> np.ndarray:
    """Re-space one stretch along the source shoreline between its endpoints.

    The source is simplified to a fraction of the LOCAL element size first.
    A shoreline digitised at metres cannot be represented by 400 m elements,
    and walking it at 400 m produces chords that turn sharply against each
    other: measured on the Futtsu patch, every surviving QA failure was a
    coastal element 1.8-2.3 km out, where the transition is coarse and the
    coast is not.  Simplifying is not throwing detail away -- there is no
    room for it at that size -- and what remains is still bounded by
    ``coastline_tolerance_m`` against the base polyline.
    """
    import shapely

    coords = _source_substring(pts, shoreline)
    if coords is None:
        return _subdivide(pts, size)
    h = float(np.median(_size_at(size, coords)))
    simple = shapely.simplify(shapely.LineString(coords), simplify_frac * h)
    out = np.asarray(simple.coords, dtype=float)[:, :2]
    return _walk(out if len(out) >= 2 else coords, size)


def rim_constraints(
    nodes,
    selection: PatchSelection,
    *,
    size,
    coastline: str = "resample",
    shoreline=None,
    tolerance_m: float = 100.0,
) -> dict[str, Any]:
    """Build the fixed points and constrained segments the filler needs.

    Returns ``pfix`` (coordinates), ``egfix`` (index pairs into ``pfix``, one
    per rim segment), ``pfix_base`` (the base node id for each pfix row, or -1
    for a new coastline point), and ``curves`` -- the coastline each new point
    was cut from, with ``curve_of_pfix`` mapping the row to its curve.  A
    later repair pass may slide those points, and this is the only curve they
    may be slid along.  ``pfix`` alone is not enough: it fixes
    positions, and it is ``egfix`` that forces the segments between them into
    the triangulation through the CDT (``mesh_generator.py:1077``).  Without
    the segments an unconstrained Delaunay can bridge a concave rim, and the
    819 m interface edges measured on the Futtsu rim are exactly the kind of
    chord it bridges.
    """
    xy = np.asarray(nodes, dtype=float)[:, :2]
    frozen = set(selection.frozen_nodes.tolist())

    pts: list[np.ndarray] = []
    base_id: list[int] = []
    segs: list[tuple[int, int]] = []
    curves: list[np.ndarray] = []
    curve_of: dict[int, int] = {}
    n_resampled = 0
    n_new = 0

    # A free rim node's two rim edges are both on the physical boundary, and
    # this is structural rather than lucky: an interface edge is shared with a
    # retained face, so BOTH its endpoints are frozen by definition.  Hence
    # every maximal run of free rim nodes is a stretch of coastline bounded by
    # two frozen nodes, and it is exactly the set of stretches that may be
    # re-cut.  Nothing else on the rim can move.
    for ring in selection.rings:
        ring = np.asarray(ring, dtype=np.int64)
        m = len(ring)
        is_free = np.array([int(v) not in frozen for v in ring], dtype=bool)
        ring_rows: list[int] = []

        if is_free.all():
            # An island taken whole: no frozen node anchors it, so there is no
            # stretch to anchor a resample between, and its shape is kept as
            # it stands.  Its nodes are emitted as NEW points even though the
            # coordinates are the base ones -- they are free, no retained face
            # holds them, and calling them frozen makes stitch_patch look for
            # them in a node map that only carries survivors (review finding
            # 6, 2026-09-22).
            ring_rows = [_push(pts, base_id, xy[v], -1) for v in ring]
            # Its own closed curve, and its own nodes bound to it.  Without
            # this the island has no provenance at all: the driver measures
            # every new boundary node against the curves it was given, so an
            # island that did not move measured 700 m from a mainland stretch
            # and the whole run was rejected (second review, finding 1).
            curves.append(np.asarray(xy[np.append(ring, ring[0])], dtype=float))
            for row in ring_rows:
                curve_of[row] = len(curves) - 1
            n_new += len(ring)
        else:
            start = int(np.argmax(~is_free))
            order = [(start + k) % m for k in range(m)]
            i = 0
            while i < m:
                v = order[i]
                ring_rows.append(_push(pts, base_id, xy[ring[v]], int(ring[v])))
                j = i + 1
                run: list[int] = []
                while j < m and is_free[order[j]]:
                    run.append(order[j])
                    j += 1
                if run:
                    idx = np.concatenate([[ring[v]], ring[run], [ring[order[j % m]]]])
                    new = coastline_points(xy[idx], size, mode=coastline,
                                           shoreline=shoreline,
                                           tolerance_m=tolerance_m)
                    curves.append(coastline_curve(xy[idx], coastline, shoreline))
                    for q in new[1:-1]:
                        curve_of[_push(pts, base_id, q, -1)] = len(curves) - 1
                        ring_rows.append(len(pts) - 1)
                        n_new += 1
                    n_resampled += len(run)
                i = j

        for k in range(len(ring_rows)):
            segs.append((ring_rows[k], ring_rows[(k + 1) % len(ring_rows)]))

    pfix = np.asarray(pts, dtype=float)
    egfix = np.asarray(segs, dtype=np.int64)
    return {
        "pfix": pfix,
        "egfix": egfix,
        "pfix_base": np.asarray(base_id, dtype=np.int64),
        "curves": curves,
        "curve_of_pfix": curve_of,
        "n_pfix": int(len(pfix)),
        "n_egfix": int(len(egfix)),
        "n_coastline_nodes_replaced": n_resampled,
        "n_coastline_nodes_new": n_new,
        "coastline_mode": coastline,
    }


def _push(pts: list, base_id: list, q, bid: int) -> int:
    pts.append(np.asarray(q, dtype=float))
    base_id.append(int(bid))
    return len(pts) - 1


def hole_polygon(pfix: np.ndarray, egfix: np.ndarray):
    """The domain the filler meshes, built from the FINAL rim.

    Not the union of the removed triangles: a resampled coastline leaves that
    polygon by up to the tolerance, and meshing the old outline would put the
    new boundary nodes outside the domain.

    Nor the union of every face ``polygonize`` returns.  That fills the
    exclusions: four nested squares plus a disjoint one come back as area 101
    where the domain is 57, and an island inside the cut would be handed to
    the filler as water with its coastline reduced to an interior line
    (review finding 1, 2026-09-22).  Rings are assembled by nesting parity,
    each shell carrying the holes whose immediate parent it is.
    """
    import shapely

    rings, is_hole = boundary_rings(pfix, egfix)
    if not rings:
        raise ValueError("the rim segments do not close a polygon")
    polys = [shapely.Polygon(np.asarray(pfix, dtype=float)[r]) for r in rings]
    areas = np.array([p.area for p in polys])
    order = np.argsort(areas)
    shells = []
    for k in range(len(rings)):
        if is_hole[k]:
            continue
        holes = []
        for j in order:
            if not is_hole[j] or areas[j] >= areas[k]:
                continue
            # the hole belongs to its IMMEDIATE parent, the smallest shell
            # that contains it
            parent = next((m for m in order
                           if not is_hole[m] and areas[m] > areas[j]
                           and shapely.contains(polys[m], polys[j])), None)
            if parent == k:
                holes.append(np.asarray(pfix, dtype=float)[rings[j]])
        shells.append(shapely.Polygon(
            np.asarray(pfix, dtype=float)[rings[k]], holes))
    out = shapely.union_all(shells)
    if out.is_empty:
        raise ValueError("the rim segments do not close a polygon")
    # A valid Polygon is not evidence that it is the domain the rim asked
    # for.  Two rings touching along an edge with distinct node ids union
    # into a rectangle and two units of constraint simply vanish from its
    # boundary (second review, finding 6).  The lengths must agree.
    want = float(sum(shapely.LineString(np.asarray(pfix, dtype=float)[
        np.append(r, r[0])]).length for r in rings))
    got = float(shapely.length(shapely.boundary(out)))
    if abs(got - want) > 1e-6 * max(1.0, want):
        raise ValueError(
            f"the assembled hole boundary is {got:.6g} m against {want:.6g} m "
            "of rim constraints; the rings touch or overlap, which this does "
            "not represent")
    return out


# --------------------------------------------------------------------------
# the sizing the filler uses
# --------------------------------------------------------------------------


def ambient_size_field(nodes, elements, *, smooth_passes: int = 20) -> np.ndarray:
    """The base mesh's own edge length at each node, smoothed.

    This is what the patch must grow into.  A scalar "ambient" is a fiction at
    the rim: measured on the Futtsu selection the rim edges run 180 / 467 /
    819 m against a nominal 350 m field, so a field that ignores the local
    mesh and stops at 350 leaves DistMesh fighting the rim it is pinned to.

    The raw per-node mean is too rough to mesh against.  On the base mesh its
    slope runs to 0.686 -- above the 0.414 that C4 allows between neighbours,
    and 26 % of edges exceed the recipe's own 0.165 -- because a node-mean
    over a fan is not a smooth function.  Twenty Jacobi passes bring the p90
    slope from 0.243 to 0.094 and the maximum to 0.206 while leaving the
    median size where it was (421 -> 439 m), which is the difference between
    a field a mesh generator can follow and one it cannot.
    """
    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    e = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]],
                                     tri[:, [2, 0]]]), axis=1), axis=0)
    ln = np.linalg.norm(xy[e[:, 0]] - xy[e[:, 1]], axis=1)
    tot = np.zeros(len(xy))
    cnt = np.zeros(len(xy))
    np.add.at(tot, e.ravel(), np.repeat(ln, 2))
    np.add.at(cnt, e.ravel(), 1.0)
    out = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
    if np.isnan(out).any():
        out[np.isnan(out)] = np.nanmedian(out)
    for _ in range(int(smooth_passes)):
        tot[:] = 0.0
        cnt[:] = 0.0
        np.add.at(tot, e.ravel(), out[e[:, ::-1]].ravel())
        np.add.at(cnt, e.ravel(), 1.0)
        out = 0.5 * out + 0.5 * (tot / np.maximum(cnt, 1.0))
    return out


def patch_sizing(
    nodes,
    elements,
    regions,
    *,
    distmesh_scale: float = 1.2,
):
    """A callable ``h(points) -> edge length`` for the hole, in mesh CRS metres.

    ``regions`` is a list of ``(geometry_in_mesh_crs, target_h_m, width_m)``
    or ``(geometry, target, width, priority)``.  Inside a region the size is
    its target; over the next ``width_m`` it ramps linearly to **the base
    mesh's own size at that point**, and beyond it is the base mesh's size.
    It never exceeds the base size, so a refinement can only refine.

    **Overlaps, and why a target is a ceiling.**  In a refinement a target
    says "no coarser than this here", not "exactly this here", so wherever
    regions meet the size is the SMALLEST any of them asks for.  Everyone
    gets at least what they declared and nobody is surprised.

    The other rule -- a core imposing its own target, with ``priority``
    deciding overlaps -- was implemented first and is wrong, which measuring
    the field showed.  A 60 m core 600 m from a 30 m core sits where the 30 m
    region's ramp wants 129 m, and imposing 60 dips the field 129 -> 60 ->
    129 over 200 m; imposing 120 makes it jump the other way, and the
    measured slope went 0.50, 0.92, 1.93 for targets of 60, 90 and 120 m
    against a C4 reference of 0.414.  Coarsening on purpose is a thing a
    SIZING recipe does, where the whole mesh is rebuilt; a patch is clipped
    by the base size anyway and can only refine.  ``priority`` therefore has
    no effect here, and ``region_conflicts`` says so rather than ignoring it.

    The obvious rule -- ``min(target + gradation * distance, ambient)`` -- is
    the one that fails, and it fails at the seam.  Around Futtsu the declared
    0.165 gradation reaches 422 m at the rim while the base mesh there runs to
    811 m, so the patch arrives beside a retained element four times its area:
    C4 measured 0.649 against a 0.5 gate, and six edges over the limit.
    Widening the hole until the ramp catches up is not available either -- at
    3.7 km the cut crosses the Futtsu spit and severs the mesh.  Landing on
    the ambient field by construction costs nothing and removes the failure
    mode; what varies instead is the LOCAL slope, which the caller should
    check with :func:`effective_gradation` against what C4 allows.
    """
    base_size = base_size_field(nodes, elements)
    geoms = [_region_geometry(r) for r in regions]

    def h(points):
        p = np.atleast_2d(np.asarray(points, dtype=float))[:, :2]
        base = base_size(p)
        out = base.copy()
        for g in geoms:
            out = np.minimum(out, _region_contribution(g, p, base))
        return np.minimum(out, base) / distmesh_scale

    return h


def base_size_field(nodes, elements):
    """``f(points) -> the base mesh's own edge length there``, in metres.

    The ambient field of :func:`ambient_size_field`, interpolated, with the
    field's maximum outside the mesh.  Shared by :func:`patch_sizing` and
    :func:`region_conflicts` so that what the report evaluates is the field
    the mesher is given, not a second formula that resembles it.
    """
    from matplotlib.tri import LinearTriInterpolator, Triangulation

    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    amb = ambient_size_field(xy, tri)
    interp = LinearTriInterpolator(Triangulation(xy[:, 0], xy[:, 1], tri), amb)
    amb_max = float(np.nanmax(amb))

    def f(points):
        p = np.atleast_2d(np.asarray(points, dtype=float))[:, :2]
        out = np.asarray(interp(p[:, 0], p[:, 1]), dtype=float)
        return np.where(np.isfinite(out), out, amb_max)

    return f


def _region_geometry(region):
    """``(boundary, polygon, target, width, priority)`` for one region."""
    import shapely

    g, t, w, pr = _region4(region)
    edge = shapely.boundary(g) if g.geom_type in ("Polygon", "MultiPolygon") else g
    return edge, g, t, w, pr


def _region_contribution(geom, points, base):
    """What one region asks for at ``points``, given the base size there.

    Its target in the core, ramping linearly to the base mesh's own size
    over ``width``.  This is the one expression of the rule; everything that
    needs to know what a region asks for calls it, so a report cannot drift
    away from the field (fourth review: the report used ``target/(1-u)`` and
    named a region finer than the field ever made it).
    """
    import shapely

    edge, poly, target, width, _pr = geom
    p = np.atleast_2d(np.asarray(points, dtype=float))[:, :2]
    pt = shapely.points(p[:, 0], p[:, 1])
    d = np.where(np.asarray(shapely.contains(poly, pt)), 0.0,
                 shapely.distance(pt, edge))
    # A zero width is legitimate -- it means the base mesh is already at or
    # below the target, so there is nothing to ramp -- and d/0 makes the
    # whole field NaN, which DistMesh accepts and then produces nothing from.
    u = np.clip(d / width, 0.0, 1.0) if width > 0 else (d > 0).astype(float)
    return target + (np.asarray(base, dtype=float) - target) * u


def _sample_points(geom, n_target: int = 400,
                   rounds: int = 12) -> tuple[np.ndarray, bool]:
    """A lattice of points inside ``geom``, and whether it is an area sample.

    The spacing follows the area, which is wrong for a shape that is thin:
    a 100 x 10 m rectangle got a 15.8 m lattice, no row of it landed inside,
    and the single representative point was then weighted as the whole
    polygon -- reporting 1 % of the water as 100 % of it (fifth review).
    So the lattice is halved until it actually resolves the shape, and when
    even that fails the caller is told the answer is one point, not an area.
    """
    import shapely

    xmin, ymin, xmax, ymax = geom.bounds
    area = float(geom.area)
    spacing = max(np.sqrt(area / max(n_target, 1)), 1e-12) if area > 0 else 1.0
    for _ in range(rounds):
        gx, gy = np.meshgrid(np.arange(xmin, xmax + spacing, spacing),
                             np.arange(ymin, ymax + spacing, spacing))
        p = np.column_stack([gx.ravel(), gy.ravel()])
        keep = np.asarray(shapely.contains(geom,
                                           shapely.points(p[:, 0], p[:, 1])))
        if keep.sum() >= n_target // 4:
            return p[keep], True
        spacing *= 0.5
    rep = shapely.get_coordinates(geom.representative_point())
    return (np.vstack([p[keep], rep]) if keep.any() else rep), False


def _region4(region):
    """``(geom, target, width[, priority])`` -> a four-tuple."""
    if len(region) == 4:
        g, t, w, pr = region
    elif len(region) == 3:
        g, t, w = region
        pr = 0.0
    else:
        raise ValueError("a region is (geometry, target_h_m, width_m[, priority])")
    return g, float(t), float(w), float(pr)


def field_gradation(h, footprint, spacing: float | None = None,
                    max_area_change: float = 0.5, *,
                    refine: int = 6, max_samples: int = 2_000_000,
                    ) -> dict[str, Any]:
    """The slope the sizing field ACTUALLY has, measured, not derived.

    :func:`effective_gradation` reports ``(ambient - target) / width`` for
    each region separately.  That is a formula about one region, and the
    field is not one region: it has an ambient term the formula omits, and
    where several regions meet it has whatever their minimum has.  A 60 m
    core 600 m from a 30 m core sits in a place the 30 m region's ramp wants
    at 129 m, so the field dips 129 -> 60 -> 129 over 200 m -- a slope no
    per-region number reports, and one that makes the fill hard enough that
    ten seeds all left a gate violation.

    So this samples ``h`` on a lattice over the footprint and differences it.
    The comparison number is the same C4 reference: two similar neighbouring
    triangles one size apart sit exactly on the 0.5 gate at ``g = 0.414``.
    It is still a diagnostic -- C4 is an area ratio across an edge, and it is
    gated on the finished mesh -- but it is a diagnostic about the field that
    exists rather than one that was intended.
    """
    import shapely

    xmin, ymin, xmax, ymax = footprint.bounds
    if spacing is None:
        spacing = max(25.0, min(xmax - xmin, ymax - ymin) / 200.0)
    limit = 1.0 / np.sqrt(1.0 - max_area_change) - 1.0
    gx, gy = np.meshgrid(np.arange(xmin, xmax + spacing, spacing),
                         np.arange(ymin, ymax + spacing, spacing))
    inside = np.asarray(shapely.contains(footprint,
                                         shapely.points(gx.ravel(), gy.ravel()))
                        ).reshape(gx.shape)
    unknown = {
        "spacing_m": float(spacing),
        "n_samples": int(inside.sum()),
        "n_partial_samples": 0,
        "max_slope": None, "p99_slope": None,
        "c4_reference_gradation": float(limit),
        "fraction_above_reference": None,
        "measured": False,
        "note": "no lattice sample has a neighbour inside the footprint; the "
                "slope is unknown, not zero -- pass a smaller spacing",
    }
    if not inside.any():
        # Nothing of the footprint is on the lattice: refine before giving
        # up, for the same reason a one-sided sample is refined below.
        if refine > 0:
            return field_gradation(h, footprint, spacing / 2.0,
                                   max_area_change, refine=refine - 1,
                                   max_samples=max_samples)
        return unknown
    size = np.asarray(h(np.column_stack([gx.ravel(), gy.ravel()])),
                      dtype=float).reshape(gx.shape)

    def derivative(field, ok, axis):
        """One derivative per sample, from neighbours INSIDE the footprint.

        Both ends of a difference have to be in the footprint. Sampling
        outside it means sampling outside the MESH, where the ambient
        interpolator returns its maximum: differencing across that edge
        reported a slope of 24 on a field whose p99 is 0.37. And masking
        first is no better -- np.gradient then returns NaN for the whole
        stencil, every sample is dropped and a 100 + x field over a 20 m
        tall box read as flat (fourth review). So the stencil is trimmed
        rather than the field: central where both neighbours are in,
        one-sided where one is, unknown where neither.
        """
        f = np.moveaxis(field, axis, 0)
        m = np.moveaxis(ok, axis, 0)
        d = (f[1:] - f[:-1]) / spacing          # between adjacent samples
        good = m[1:] & m[:-1]
        tot = np.zeros_like(f)
        cnt = np.zeros_like(f)
        for sl_a, sl_b in ((slice(1, None), slice(None)), (slice(0, -1), slice(None))):
            tot[sl_a] += np.where(good, d, 0.0)[sl_b]
            cnt[sl_a] += good[sl_b]
        with np.errstate(invalid="ignore", divide="ignore"):
            out = np.where(cnt > 0, tot / np.maximum(cnt, 1), np.nan)
        return np.moveaxis(out, 0, axis), np.moveaxis(cnt > 0, 0, axis)

    dy, have_y = derivative(size, inside, 0)
    dx, have_x = derivative(size, inside, 1)
    partial = int((inside & ~(have_x & have_y)).sum())
    have = (have_x | have_y) & inside
    if not have.any():
        if refine > 0:
            return field_gradation(h, footprint, spacing / 2.0,
                                   max_area_change, refine=refine - 1,
                                   max_samples=max_samples)
        return unknown
    partial_fraction = float(partial) / max(int(inside.sum()), 1)
    # Refine for a shape the lattice cannot see across, not for the handful
    # of slivers every real footprint has along its boundary: the Futtsu
    # hole has 7 partial samples out of 27,035, and refining for those takes
    # the lattice from 27 thousand points to 27 million.
    if partial_fraction > 0.01 and refine > 0 and 4 * gx.size <= max_samples:
        # A sample with support on one axis only reports the magnitude of
        # what it can see, which is a LOWER bound: a 40 m wide channel at
        # the default 25 m spacing has no transverse neighbour anywhere, and
        # a field rising 1 m per metre across it measured as flat (fifth
        # review). Halving the lattice is what makes the missing direction
        # visible, so it is halved rather than qualified in a footnote.
        finer = field_gradation(h, footprint, spacing / 2.0,
                                max_area_change, refine=refine - 1,
                                max_samples=max_samples)
        if finer.get("measured") and not finer.get("n_partial_samples"):
            return finer
    # An axis with no support contributes nothing rather than a NaN: the
    # result is then a LOWER bound on the slope, and the count of such
    # samples is reported rather than buried.
    slopes = np.hypot(np.where(have_x, dx, 0.0), np.where(have_y, dy, 0.0))[have]
    slopes = slopes[np.isfinite(slopes)]
    if not slopes.size:
        return unknown
    return {
        "spacing_m": float(spacing),
        "n_samples": int(inside.sum()),
        "n_partial_samples": partial,
        "max_slope": float(slopes.max()),
        "p99_slope": float(np.percentile(slopes, 99)),
        "c4_reference_gradation": float(limit),
        "fraction_above_reference": float((slopes > limit).mean()),
        "measured": True,
        "partial_support": bool(partial),
    }


def region_resolution(nodes, elements, geometry, target_h_m: float, *,
                      spacing: float | None = None,
                      coverage_tolerance: float = 1.25,
                      max_samples: int = 200_000) -> dict[str, Any]:
    """What a declared region actually got, measured over its AREA.

    Edge statistics inside a region are counted per edge, so the finely
    meshed water contributes most of them: a request whose northern half is
    at 30 m and whose southern half was never refined at all still shows a
    median of 28.5 m, because the fine half owns 1,100 of the edges and the
    coarse half owns eight. A patch that leaves 46 % of the requested water
    at the base size passed the edge-median gate on the delivered Futtsu mesh
    (fifth review). So the question is asked of the WATER, not of the edges:
    sample the region on a lattice, find the element covering each sample,
    and ask how big that element is.

    The size of an element is its **equivalent edge**, the edge of the
    equilateral triangle with the same area, which is what ``target_h_m``
    means as a resolution. The longest edge would be a different and
    stricter contract: it runs about 1.2 times the equivalent edge on a
    DistMesh fill, so gating on it would reject meshes that deliver the
    requested resolution.

    Samples that fall outside the mesh are reported, not counted: a fishery
    boundary may legitimately run onto land, and water that does not exist
    cannot be refined.
    """
    import shapely
    from matplotlib.tri import Triangulation

    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    target = float(target_h_m)
    xmin, ymin, xmax, ymax = geometry.bounds
    if spacing is None:
        spacing = max(target / 3.0, 1e-6)
    # A big region at a fine target would otherwise ask for a lattice nobody
    # can afford; the coarser spacing is reported with the result.
    span = max(xmax - xmin, 1e-9) * max(ymax - ymin, 1e-9)
    if span / (spacing * spacing) > max_samples:
        spacing = float(np.sqrt(span / max_samples))
    gx, gy = np.meshgrid(np.arange(xmin, xmax + spacing, spacing),
                         np.arange(ymin, ymax + spacing, spacing))
    p = np.column_stack([gx.ravel(), gy.ravel()])
    keep = np.asarray(shapely.contains(geometry,
                                       shapely.points(p[:, 0], p[:, 1])))
    p = p[keep]
    if not p.shape[0]:
        # A region thinner than the lattice: its representative point is the
        # one place it certainly covers, and one sample is better than a
        # claim of nothing.
        p = shapely.get_coordinates(geometry.representative_point())
    finder = Triangulation(xy[:, 0], xy[:, 1], tri).get_trifinder()
    found = finder(p[:, 0], p[:, 1])
    inside = found >= 0
    u = xy[tri[:, 1]] - xy[tri[:, 0]]
    v = xy[tri[:, 2]] - xy[tri[:, 0]]
    area = 0.5 * np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
    h_eq = np.sqrt(4.0 * area / np.sqrt(3.0))
    out: dict[str, Any] = {
        "target_h_m": target,
        "spacing_m": float(spacing),
        "n_samples": int(p.shape[0]),
        "n_outside_mesh": int((~inside).sum()),
        "outside_mesh_fraction": float((~inside).mean()),
        "coverage_tolerance": float(coverage_tolerance),
    }
    if not inside.any():
        out.update(median_m=None, p90_m=None, max_m=None, median_ratio=None,
                   p90_ratio=None, max_ratio=None, covered_fraction=0.0)
        return out
    h = h_eq[found[inside]]
    r = h / target
    out.update(
        median_m=float(np.median(h)), p90_m=float(np.percentile(h, 90)),
        max_m=float(h.max()), median_ratio=float(np.median(r)),
        p90_ratio=float(np.percentile(r, 90)), max_ratio=float(r.max()),
        covered_fraction=float((r <= coverage_tolerance).mean()),
    )
    return out


def region_conflicts(regions, names=None, *, base_size=None) -> dict[str, Any]:
    """Which declared cores overlap, and what the shared water gets.

    Overlapping fisheries are an ordinary input -- two rights over the same
    water -- and in a refinement they do not conflict: a target is a ceiling,
    so the shared water takes the smallest of them and every region gets at
    least what it asked for.  What is worth reporting is that one region's
    water will be FINER than it declared, because that costs time steps and
    elements somebody has to pay for.

    ``base_size`` is the base mesh's own size as a function of position, from
    :func:`base_size_field`.  A region's ramp runs from its target to THAT,
    so without it a transition's value is not a number this function knows:
    it then reports only what a core overlap makes certain and says as much
    in ``transitions_evaluated``.  Pass it whenever the base mesh is at hand
    -- it is what makes the report the field's rather than a formula's.

    ``priority_ignored`` is true when the recipe sets different priorities on
    overlapping regions: in a refinement priority cannot change the outcome,
    and saying so beats ignoring the key.
    """
    import shapely

    r4 = [_region4(x) for x in regions]
    names = list(names) if names else [f"region_{i}" for i in range(len(r4))]
    pairs, finer = [], {}
    priority_ignored = False
    for a in range(len(r4)):
        for b in range(a + 1, len(r4)):
            inter = shapely.intersection(r4[a][0], r4[b][0])
            if inter.is_empty or inter.area <= 0:
                continue
            if r4[a][3] != r4[b][3]:
                priority_ignored = True
            pairs.append({
                "regions": [names[a], names[b]],
                "overlap_m2": float(inter.area),
                "fraction_of": {names[a]: float(inter.area / r4[a][0].area),
                                names[b]: float(inter.area / r4[b][0].area)},
                "targets_h_m": {names[a]: r4[a][1], names[b]: r4[b][1]},
                "effective_target_h_m": min(r4[a][1], r4[b][1]),
            })
    # A region comes out finer than it declared for two reasons, and only one
    # of them is an overlap of CORES: a neighbour's TRANSITION can reach it
    # and be finer there than its own target.  A 5 m core 200 m away with a
    # 1 km transition swallows a 90 m core whole, and looking only at core
    # intersections reported nothing at all (third review).  So the test is
    # the field: what this region would get alone, against what it gets --
    # sampled from the same expression the mesher is handed, because a
    # second expression that merely resembles it named a 12 m region as
    # getting 9.6 m in a field that gave it 12 (fourth review).
    geoms = [_region_geometry(x) for x in regions]
    for i, name in enumerate(names):
        cores = [j for j in range(len(r4)) if j != i
                 and not shapely.intersection(r4[i][0], r4[j][0]).is_empty]
        if base_size is None:
            # Without the base mesh a ramp's value is unknown, so only a core
            # overlap by a finer target is certain.
            smaller = [j for j in cores if r4[j][1] < r4[i][1] - 1e-9]
            if not smaller:
                continue
            taken = shapely.intersection(
                r4[i][0], shapely.union_all([r4[j][0] for j in smaller]))
            finer[name] = {
                "core_overlap_area_m2": float(taken.area),
                "area_m2": float(taken.area),
                "fraction": float(taken.area / r4[i][0].area),
                "own_target_h_m": r4[i][1],
                "gets_h_m": float(min(r4[j][1] for j in smaller)),
                "because_of": sorted(names[j] for j in smaller),
                "by_core_overlap": True,
                "area_is_estimated": False,
                "area_sampled": True,
            }
            continue
        p, is_area_sample = _sample_points(r4[i][0])
        base = np.asarray(base_size(p), dtype=float)
        contrib = np.array([_region_contribution(g, p, base) for g in geoms])
        alone = np.minimum(contrib[i], base)
        joint = np.minimum(contrib.min(axis=0), base)
        hit = joint < alone - 1e-9
        if not hit.any():
            continue
        others = [names[j] for j in range(len(r4)) if j != i
                  and (contrib[j][hit] < alone[hit] - 1e-9).any()]
        # The exact core overlap is geometry and is reported as such; the
        # part a neighbour's transition reaches is an estimate from the
        # samples, and says so rather than dressing one point as an area.
        exact = shapely.intersection(
            r4[i][0], shapely.union_all([r4[j][0] for j in cores])) \
            if cores else shapely.Polygon()
        fraction = float(hit.mean()) if is_area_sample else None
        finer[name] = {
            "core_overlap_area_m2": float(exact.area),
            "area_m2": (float(r4[i][0].area) * fraction
                        if fraction is not None else float(exact.area)),
            "fraction": fraction,
            "own_target_h_m": r4[i][1],
            "gets_h_m": float(joint[hit].min()),
            "because_of": sorted(others),
            "by_core_overlap": bool(cores),
            "n_samples": int(p.shape[0]),
            "area_is_estimated": True,
            "area_sampled": bool(is_area_sample),
        }
    return {
        "overlapping_pairs": pairs,
        "finer_than_declared": finer,
        "any_overlap": bool(pairs),
        "priority_ignored": priority_ignored,
        "transitions_evaluated": base_size is not None,
    }


def effective_gradation(nodes, elements, regions) -> dict[str, Any]:
    """The steepest size ramp :func:`patch_sizing` will actually ask for.

    The declared gradation sets the WIDTH; the slope that results is
    ``(local base size - target) / width``, and where the base mesh is coarse
    that is steeper than declared.  It has a hard limit rather than a taste:
    neighbouring elements one size apart differ in area by
    ``1 - 1/(1+g)^2``, so the C4 gate of 0.5 is reached at ``g = 0.414``.
    Past that the patch cannot pass QA however well it is meshed.
    """
    xy = np.asarray(nodes, dtype=float)[:, :2]
    amb = ambient_size_field(xy, np.asarray(elements, dtype=np.int64))
    worst = 0.0
    per_region = {}
    for geom, target, width, _pr in map(_region4, regions):
        import shapely

        d = shapely.distance(shapely.points(xy[:, 0], xy[:, 1]), geom)
        inside = d < width
        g = float(((amb[inside] - target) / width).max()) if inside.any() else 0.0
        per_region[str(target)] = g
        worst = max(worst, g)
    return {
        "max_effective_gradation": worst,
        # A reference number, not a verdict.  The implemented field is
        # target + (B(x) - target) * d(x) / W, whose gradient also has an
        # ambient term this does not measure, and the 1 - 1/(1+g)^2 relation
        # assumes similar neighbouring triangles where C4 is an area ratio
        # across a shared edge.  C4 itself is gated on the finished mesh.
        "c4_reference_gradation": float(1.0 / np.sqrt(1.0 - 0.5) - 1.0),
        "ramp_below_reference": bool(worst <= 1.0 / np.sqrt(1.0 - 0.5) - 1.0),
        "per_region": per_region,
    }


# --------------------------------------------------------------------------
# stitching
# --------------------------------------------------------------------------


def stitch_patch(
    base_nodes,
    base_elements,
    base_depths,
    selection: PatchSelection,
    patch_nodes,
    patch_elements,
    pfix,
    pfix_base,
    *,
    tol_m: float = 1e-6,
):
    """Put the filled patch back into the base mesh.

    Returns ``(nodes, elements, depths, node_map, report)``.  ``node_map`` is
    the old-to-new node index for every base node, ``-1`` where the base node
    was inside the hole and did not survive.  Every check in
    :func:`verify_patch` runs through this map, because the patch renumbers
    everything and a coordinate comparison on raw row order would be
    meaningless.

    Retained nodes come first and in their base order, so the map is the
    identity for everything before the first deleted node -- which keeps
    diffs of the written fort.14 readable.
    """
    from scipy.spatial import cKDTree

    xy = np.asarray(base_nodes, dtype=float)[:, :2]
    dep = np.asarray(base_depths, dtype=float)
    pnodes = np.asarray(patch_nodes, dtype=float)[:, :2]
    ptri = np.asarray(patch_elements, dtype=np.int64)
    pfix = np.asarray(pfix, dtype=float)[:, :2]
    pfix_base = np.asarray(pfix_base, dtype=np.int64)

    keep = selection.frozen_nodes
    node_map = np.full(len(xy), -1, dtype=np.int64)
    node_map[keep] = np.arange(len(keep), dtype=np.int64)

    # Which patch vertex is which fixed point?  Exact, not nearest: a fixed
    # point that DistMesh lost shows up here as a missing row rather than as a
    # silent snap onto its neighbour.
    d, idx = cKDTree(pnodes).query(pfix)
    lost = np.where(d > tol_m)[0]
    if lost.size:
        raise ValueError(
            f"{lost.size} of {len(pfix)} constrained points are absent from the "
            f"filled patch (worst offset {d.max():.3g} m); the fill dropped the "
            "rim it was pinned to")
    if np.unique(idx).size != len(idx):
        raise ValueError("two constrained points map to the same patch vertex; "
                         "the fill merged the rim")

    patch_map = np.full(len(pnodes), -1, dtype=np.int64)
    frozen_rows = pfix_base >= 0
    patch_map[idx[frozen_rows]] = node_map[pfix_base[frozen_rows]]
    if (patch_map[idx[frozen_rows]] < 0).any():
        raise ValueError("a constrained rim point refers to a base node that the "
                         "cut removed; the selection and the rim disagree")

    fresh = np.where(patch_map < 0)[0]
    patch_map[fresh] = len(keep) + np.arange(len(fresh), dtype=np.int64)

    nodes = np.vstack([xy[keep], pnodes[fresh]])
    elements = np.vstack([node_map[selection.retained], patch_map[ptri]])
    if (elements < 0).any():
        raise ValueError("an element references a node that was not carried over")

    # Depths come from the FULL base triangulation, not the retained part of
    # it: a new node inside the hole is, by construction, over base elements
    # that were removed, so interpolating on the retained mesh alone would
    # send every one of them to the nearest-node fallback.
    from fvcom_mesh_tools.refine import depths_from_base

    if len(fresh):
        new_dep, n_outside = depths_from_base(
            xy, np.asarray(base_elements, dtype=np.int64), dep, pnodes[fresh])
    else:
        new_dep, n_outside = np.zeros(0), 0
    depths = np.concatenate([dep[keep], new_dep])

    report = {
        # Where every constrained rim point ended up.  boundary_after_patch
        # used to recover this by looking the coordinates up in the finished
        # mesh, which made it callable only BEFORE the repair slid anything.
        # A precondition nobody can see is a defect waiting for a refactor
        # (third review, finding 1).
        "pfix_new": patch_map[idx],
        "n_nodes": int(len(nodes)),
        "n_elements": int(len(elements)),
        "n_nodes_retained": int(len(keep)),
        "n_nodes_new": int(len(fresh)),
        "n_nodes_deleted": int(len(xy) - len(keep)),
        "n_new_depths_outside_base": int(n_outside),
    }
    return nodes, elements, depths, node_map, report


def refresh_depths(base_nodes, base_elements, base_depths, nodes, moved, depths):
    """Re-read the base bathymetry wherever the repair moved a node.

    ``stitch_patch`` evaluates the base field at the positions the fill
    produced; :func:`improve_patch` then moves some of them, and a depth left
    behind at the old position is not the base field at the delivered
    coordinate.  On a ``5 + x`` field a node moved from (0.2, 0.3) to (1, 1)
    keeps 5.2 where 6.0 is right (review finding 9, 2026-09-22).  Frozen
    nodes never move, so their depths are untouched either way.
    """
    from fvcom_mesh_tools.refine import depths_from_base

    depths = np.array(depths, dtype=float)
    moved = np.asarray(moved, dtype=bool)
    if not moved.any():
        return depths, 0
    fresh, n_outside = depths_from_base(base_nodes, base_elements, base_depths,
                                        np.asarray(nodes, dtype=float)[moved])
    depths[moved] = fresh
    return depths, int(n_outside)


def boundary_after_patch(base_elements, selection, rc, node_map, pfix_new):
    """The boundary edge set the patched mesh must have, in final node ids.

    Two parts and nothing else: the base mesh's boundary edges that the cut
    did not take, and the rim segments that are supposed to BE boundary --
    the coastline chains, not the interface ones.  ``pfix_new`` is
    ``stitch_patch``'s report entry of the same name: where each constrained
    rim point landed, by node id, which the repair does not change.  Comparing the finished
    mesh against this is what notices a face that is simply gone: every other
    invariant survives a missing patch triangle intact.
    """
    tri = np.asarray(base_elements, dtype=np.int64)
    nm = np.asarray(node_map, dtype=np.int64)
    u, c = _edge_table(tri)
    base_boundary = u[c == 1]
    taken = _edge_table(tri[selection.removed])[0]
    taken_set = {tuple(x) for x in taken.tolist()}
    kept = [tuple(sorted(nm[list(e)].tolist())) for e in base_boundary.tolist()
            if tuple(e) not in taken_set]
    if any(v < 0 for e in kept for v in e):
        raise ValueError("a retained boundary edge lost a node in the patch")

    base_id = np.asarray(rc["pfix_base"], dtype=np.int64)
    row = np.asarray(pfix_new, dtype=np.int64)
    if row.shape[0] != base_id.shape[0]:
        raise ValueError("pfix_new must give one output node per pfix row")
    interface = {tuple(sorted(e)) for e, phys
                 in zip(selection.rim_edges.tolist(), selection.physical_rim)
                 if not phys}
    coast = []
    for i, j in np.asarray(rc["egfix"], dtype=np.int64).tolist():
        if base_id[i] >= 0 and base_id[j] >= 0 and \
                tuple(sorted((int(base_id[i]), int(base_id[j])))) in interface:
            continue
        coast.append(tuple(sorted((int(row[i]), int(row[j])))))
    return set(kept) | set(coast)


def introduced_violations(qa_checks, n_retained_elements: int,
                          elements=None) -> list[dict]:
    """The QA failures this patch is answerable for.

    A refinement may not be held to a standard its base does not meet. The
    goto2023 production mesh -- the `current` hydro baseline's own grid --
    fails C1 at element 2101 with a 28.99 degree angle, 18 km from Futtsu.
    The contract says that element is frozen, so the patch keeps the failure
    and would be blamed for it by an absolute gate: every seed came back
    "QA 20/21" for a defect it is forbidden to touch.

    An offender is INHERITED when every element it involves is a retained
    one, because a retained element is bit-identical to the base. Everything
    else -- a patch element, a seam edge between a patch and a retained
    element, a frozen node whose fan now contains a patch face -- is the
    patch's, and is what this returns.

    ``qa_checks`` are ``QAReport.checks``.  Attribution reads
    ``check.offender_ids`` -- every offender, whatever ``max_offenders`` was
    -- and not ``check.offenders``, which is the display list and is capped.
    Reading the capped list made a check that named 10,000 of its 10,368
    inherited failures produce 368 "unattributed" entries, and so rejected an
    UNCHANGED mesh (fourth review).  A raised cap only moves that failure;
    the display list and the gate are two different lists.

    An offender this cannot place is counted as the patch's, and so is every
    violation a check counted but did not identify: attribution is a claim,
    and the absence of one is not a claim of innocence.  ``elements`` is the
    patched connectivity, needed to attribute a node offender (C5 valence):
    a frozen node is inherited only while every face around it is retained.
    Without it a node offender is always counted as the patch's, which errs
    towards blaming the patch.
    """
    tri = None if elements is None else np.asarray(elements, dtype=np.int64)
    n_elements = None if tri is None else int(tri.shape[0])
    out: list[dict] = []
    for check in qa_checks:
        if getattr(check, "status", "") != "fail":
            continue
        # Display first, identities second: the entries reported are the
        # readable ones, the gate is decided on all of them.
        shown = list(getattr(check, "offenders", []) or [])
        identified = list(getattr(check, "offender_ids", None) or shown)
        detail = {}
        for off in shown:
            detail.setdefault(_offender_key(off), off)
        for off in identified:
            kind = off.get("kind")
            ident = off.get("id")
            if kind == "element":
                elems = [_element_id(ident, n_elements)]
            elif kind == "edge" and "elements" in off:
                elems = [_element_id(e, n_elements) for e in off["elements"]]
            elif kind == "node" and tri is not None:
                node = _element_id(ident, int(tri.max()) + 1 if tri.size else 0)
                elems = None if node is None else \
                    np.flatnonzero((tri == node).any(axis=1)).tolist()
            else:
                # An offender this cannot place -- a kind it does not know, an
                # id that is a pair or not an index at all, or no connectivity
                # to look it up in -- is the patch's. Attribution is a claim,
                # and the absence of one is not a claim of innocence.
                elems = None
            if elems is not None and elems and all(
                    e is not None and e < n_retained_elements for e in elems):
                continue
            record = detail.get(_offender_key(off), off)
            out.append({"check": check.check_id, "requirement": check.requirement,
                        "observed": check.observed, **record})
        # A failing check that identified nobody, or fewer than it counted,
        # cannot be shown to be the base's -- node_index_valid fails with no
        # offenders at all, and an empty list read as "0 introduced" once
        # accepted a mesh with a failing gate (fourth review).
        counted = int(getattr(check, "n_violations", 0) or 0)
        unlisted = max(0, counted - len(identified)) or (0 if identified else 1)
        if unlisted:
            out.append({"check": check.check_id, "requirement": check.requirement,
                        "observed": check.observed, "kind": "unattributed",
                        "n_unattributed": unlisted,
                        "note": "the check identified fewer offenders than it "
                                "counted, so they cannot be shown to be the "
                                "base's"})
    return out


def _offender_key(off: dict):
    """What identifies an offender across the two lists.

    The id alone is not enough: a C4 edge offender is identified by the pair
    of elements it sits between and need not carry an id at all, and keying
    on the id alone handed every such edge the first one's detail.
    """
    return (off.get("kind"), _ident_key(off.get("id")),
            _ident_key(off.get("elements")), _ident_key(off.get("segment")))


def _ident_key(ident):
    """A hashable form of an offender id (an index, or a pair of them)."""
    if ident is None:
        return None
    if isinstance(ident, (list, tuple, np.ndarray)):
        return tuple(np.asarray(ident).ravel().tolist())
    try:
        hash(ident)
    except TypeError:
        # A malformed id this cannot hash is still an offender to report;
        # a dict one raised TypeError before it could be refused (fifth
        # review). Unplaceable is the answer, not an exception.
        return ("unhashable", repr(ident))
    return ident


def _element_id(ident, n: int | None) -> int | None:
    """``ident`` as an element index, or None when it is not one.

    A non-integer id, a bool, a negative index or one past the end of the
    mesh is not an element this can place, and a patch is not exonerated by
    an id nobody can look up: ``0.5`` and ``-1`` both compared happily
    against ``n_retained_elements`` and cleared the patch (fourth review).
    """
    if isinstance(ident, (bool, np.bool_)) or not isinstance(
            ident, (int, np.integer)):
        return None
    i = int(ident)
    if i < 0 or (n is not None and i >= n):
        return None
    return i


# --------------------------------------------------------------------------
# verification
# --------------------------------------------------------------------------


def verify_patch(
    base_nodes,
    base_depths,
    base_elements,
    selection: PatchSelection,
    nodes,
    elements,
    depths,
    node_map,
    *,
    open_boundaries=(),
    expected_boundary=None,
) -> dict[str, Any]:
    """Check the frozen zone really is frozen, through the node map.

    Several separate claims, because "the mesh outside is unchanged" is
    several claims wearing one coat: the retained nodes kept their
    coordinates, they kept their depths, the retained faces are all still
    there with the same vertices, their orientation did not flip, the
    interface segments are still shared rather than split, and the open
    boundary is the same list in the same order, and the boundary is the one
    the patch was built to have.  A patch can satisfy all but one of those
    and still break the model.  ``area_change_fraction`` is
    reported rather than gated: it moves legitimately when the coastline is
    resampled (-0.008 % on the Futtsu patch) and is the one number that
    notices a face quietly dropped.
    """
    xy = np.asarray(base_nodes, dtype=float)[:, :2]
    dep0 = np.asarray(base_depths, dtype=float)
    new_xy = np.asarray(nodes, dtype=float)[:, :2]
    new_dep = np.asarray(depths, dtype=float)
    tri = np.asarray(elements, dtype=np.int64)
    nm = np.asarray(node_map, dtype=np.int64)

    if not np.isfinite(new_xy).all() or not np.isfinite(new_dep).all():
        raise ValueError("the patched mesh contains non-finite coordinates or depths")

    keep = selection.frozen_nodes
    # EXACT, not within a tolerance.  These values are copied, not computed:
    # a frozen coordinate that differs by 5e-7 m has not been copied, it has
    # been recomputed, and a contract that says bit-for-bit has to be checked
    # bit-for-bit (review finding 3, 2026-09-22).
    moved = np.linalg.norm(new_xy[nm[keep]] - xy[keep], axis=1)
    ddep = np.abs(new_dep[nm[keep]] - dep0[keep])
    frozen_exact = bool(np.array_equal(new_xy[nm[keep]], xy[keep])
                        and np.array_equal(new_dep[nm[keep]], dep0[keep]))

    # MULTISET, not set: a duplicated retained face keeps every set-based
    # check happy while laying a second element on top of the first.
    from collections import Counter

    want = Counter(tuple(sorted(r)) for r in nm[selection.retained].tolist())
    have = Counter(tuple(sorted(r)) for r in tri.tolist())
    missing = {f for f, n_ in want.items() if have[f] < n_}
    # Only over the retained faces: the patch's own faces are new and are
    # supposed to be here.  What this counts is a retained face appearing
    # more often than it did -- a second element laid on the first.
    extra = sum(have[f] - n_ for f, n_ in want.items() if have[f] > n_)

    a = new_xy[tri[:, 1]] - new_xy[tri[:, 0]]
    b = new_xy[tri[:, 2]] - new_xy[tri[:, 0]]
    area2 = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]

    obc_ok = True
    for seg in open_boundaries:
        s = np.asarray(seg, dtype=np.int64)
        if (nm[s] < 0).any() or not np.array_equal(
                new_xy[nm[s]], xy[s]):
            obc_ok = False

    dup = 0
    if len(new_xy):
        from scipy.spatial import cKDTree
        pairs = cKDTree(new_xy).query_pairs(1e-3, output_type="ndarray")
        dup = int(len(pairs))
    orphan = int(len(new_xy) - np.unique(tri).size)

    # Interface segments must survive as INTERIOR edges.  A fill can insert a
    # vertex on a constrained edge; on the coastline that is harmless -- the
    # chord is where it was -- but on an interface segment it leaves a hanging
    # node, the retained face on the far side keeps the unsplit edge, and the
    # mesh is no longer conforming although every face and every coordinate
    # still checks out.
    u, c = np.unique(np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]],
                                        tri[:, [2, 0]]]), axis=1),
                     axis=0, return_counts=True)
    interior = {tuple(x) for x in u[c == 2].tolist()}
    iface = selection.rim_edges[~selection.physical_rim]
    mapped = nm[iface] if len(iface) else np.empty((0, 2), dtype=np.int64)
    split = [e for e in np.sort(mapped, axis=1).tolist() if tuple(e) not in interior]
    # No edge may carry three faces.  A duplicated or overlapping element
    # shows up here and nowhere else.
    nonmanifold = int((c > 2).sum())

    # Coverage.  Every gate above is about the faces that are there; none of
    # them notices a face that is not.  Deleting one interior patch triangle
    # left the frozen zone exact, no retained face missing, no interface
    # split, no extra face, no non-manifold edge, no inversion, no orphan --
    # and a 5,000 m2 hole in the water (second review, finding 2).  The
    # boundary is what tells: a mesh that covers what it promised has exactly
    # the boundary edges it was built to have.
    boundary = {tuple(x) for x in u[c == 1].tolist()}
    unexpected: set = set()
    absent: set = set()
    if expected_boundary is not None:
        want_b = {tuple(sorted(e)) for e in expected_boundary}
        unexpected = boundary - want_b
        absent = want_b - boundary

    # Area is the blunt instrument that catches a hole nothing else notices:
    # a dropped face keeps every other invariant intact.
    def _area(xy_, t_):
        u_ = xy_[t_[:, 1]] - xy_[t_[:, 0]]
        v_ = xy_[t_[:, 2]] - xy_[t_[:, 0]]
        return float(np.abs(u_[:, 0] * v_[:, 1] - u_[:, 1] * v_[:, 0]).sum() / 2)

    base_area = _area(xy, np.asarray(base_elements, dtype=np.int64))
    new_area = _area(new_xy, tri)

    return {
        "n_frozen_nodes": int(len(keep)),
        "max_frozen_move_m": float(moved.max()) if len(moved) else 0.0,
        "n_frozen_moved": int((moved > 0).sum()),
        "max_frozen_depth_change_m": float(ddep.max()) if len(ddep) else 0.0,
        "n_retained_faces_missing": int(len(missing)),
        "n_inverted_elements": int((area2 <= 0).sum()),
        "n_duplicate_nodes": dup,
        "n_orphan_nodes": orphan,
        "n_interface_segments_split": int(len(split)),
        "boundary_checked": bool(expected_boundary is not None),
        "n_unexpected_boundary_edges": int(len(unexpected)),
        "n_missing_boundary_edges": int(len(absent)),
        "n_extra_faces": int(extra),
        "n_nonmanifold_edges": nonmanifold,
        "frozen_exact": frozen_exact,
        "open_boundary_unchanged": bool(obc_ok),
        "area_change_fraction": float((new_area - base_area) / base_area)
        if base_area > 0 else 0.0,
        # An unchecked claim is not a satisfied one.  Without
        # expected_boundary the coverage question was never asked, and a
        # Boolean success that means "everything I was asked to look at was
        # fine" is too easy to read as a certificate (third review).
        "ok": bool(expected_boundary is not None
                   and frozen_exact and not missing and (area2 > 0).all()
                   and extra == 0 and nonmanifold == 0 and not unexpected
                   and not absent
                   and dup == 0 and orphan == 0 and obc_ok and not split),
    }


# --------------------------------------------------------------------------
# seam repair
# --------------------------------------------------------------------------


def _angles_deg(xy: np.ndarray, tri: np.ndarray) -> np.ndarray:
    """Interior angles, one row per element, degrees."""
    p = xy[tri]
    out = np.empty((len(tri), 3))
    for k, (a, b, c) in enumerate(((0, 1, 2), (1, 2, 0), (2, 0, 1))):
        u, v = p[:, b] - p[:, a], p[:, c] - p[:, a]
        nu = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)
        out[:, k] = np.degrees(np.arccos(np.clip(
            (u * v).sum(axis=1) / np.where(nu > 0, nu, 1.0), -1.0, 1.0)))
    return out


def _areas(xy: np.ndarray, tri: np.ndarray) -> np.ndarray:
    a = xy[tri[:, 1]] - xy[tri[:, 0]]
    b = xy[tri[:, 2]] - xy[tri[:, 0]]
    return 0.5 * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])


def _health(xy, tri, faces, adj, min_angle, max_angle, max_area_change) -> float:
    """The worst gate margin alone; see :func:`_scores`."""
    return _scores(xy, tri, faces, adj, min_angle, max_angle, max_area_change)[0]


def _scores(xy, tri, faces, adj, min_angle, max_angle, max_area_change):
    """How much room the worst gate in a neighbourhood has left, as a fraction.

    One number for angles and area jump together, normalised so that 1.0 is
    exactly on the gate: maximising it is the same as pushing every gate in
    the neighbourhood away from its limit, and a move that buys a C1 gain with
    a C4 loss scores no better than the loss.

    ``adj`` is the face-to-face table, and it is not an optimisation.  C4 is
    measured across an edge, and the face on the far side of that edge need
    not touch the node being moved: the two surviving C4 failures on the
    Futtsu patch were each a patch element beside a RETAINED element twice
    its area, invisible to a score that looked only at the moved node's own
    fan.  A neighbourhood therefore means the faces plus everything they
    share an edge with.
    """
    if not len(faces):
        return np.inf, np.inf
    t = tri[faces]
    ar = _areas(xy, t)
    if (ar <= 0).any():
        return -np.inf, -np.inf
    ang = _angles_deg(xy, t)
    parts = [(ang / min_angle).ravel(),
             ((180.0 - ang) / (180.0 - max_angle)).ravel()]
    nb = adj[faces]
    have = nb >= 0
    if have.any():
        mine = np.repeat(np.abs(ar), 3).reshape(-1, 3)[have]
        theirs = np.abs(_areas(xy, tri[nb[have]]))
        change = np.abs(mine - theirs) / np.maximum(mine, theirs)
        parts.append((1.0 - change) / (1.0 - max_area_change))
    all_scores = np.concatenate(parts)
    # The hard score is the worst margin; the soft one saturates, so raising
    # a bad quantity counts and polishing an already-good one does not.  A
    # move is taken when the hard score improves, or when it holds and the
    # soft one improves: the hard min alone freezes out every node whose fan
    # minimum is set by an element it cannot help, which is how the greedy
    # pass stalls at 22.96 deg on one run and 30.01 deg on the next.
    return float(all_scores.min()), float(np.minimum(all_scores, 1.5).sum())


def _face_adjacency(tri) -> np.ndarray:
    """(nfaces, 3) neighbouring face across each edge, -1 on the boundary."""
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(tri)), 3)
    slot = np.repeat(np.arange(3)[None, :], len(tri), axis=0).T.ravel()
    order = np.lexsort((e[:, 1], e[:, 0]))
    e, owner, slot = e[order], owner[order], slot[order]
    adj = np.full((len(tri), 3), -1, dtype=np.int64)
    k = np.flatnonzero(np.all(e[:-1] == e[1:], axis=1))
    adj[owner[k], slot[k]] = owner[k + 1]
    adj[owner[k + 1], slot[k + 1]] = owner[k]
    return adj


def _better(after, before, soft: bool) -> bool:
    """Is this candidate an improvement?

    Strictly better on the worst margin; or, with ``soft``, no worse on it
    and better on the saturated sum.  Either way never worse on the worst
    margin, so the pass stays monotone in the quantity the gates are written
    on.

    Soft is not free and is not the default.  Accepting a lateral move
    changes the trajectory, and a greedy pass that wanders can finish in a
    worse basin than one that does not: on this patch the strict pass
    reached 30.01 deg and the soft pass 27.46, while on a neighbouring
    configuration the strict pass stalled at 22.96 and soft carried it to
    27.46.  The driver runs strict first and reaches for soft only when the
    gates are still unmet.
    """
    if after[0] > before[0] + 1e-9:
        return True
    return bool(soft and after[0] >= before[0] - 1e-12
                and after[1] > before[1] + 1e-9)


def _incidence(tri, n_nodes):
    """node -> incident face ids, as one sorted array plus per-node slices."""
    faces = np.repeat(np.arange(len(tri)), 3)
    who = tri.ravel()
    order = np.argsort(who, kind="stable")
    who, faces = who[order], faces[order]
    idx = np.arange(n_nodes)
    return faces, np.searchsorted(who, idx), np.searchsorted(who, idx, side="right")


def _interior_edges(tri):
    """Every edge with two faces, as (a, b, face0, face1)."""
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(tri)), 3)
    order = np.lexsort((e[:, 1], e[:, 0]))
    e, owner = e[order], owner[order]
    k = np.flatnonzero(np.all(e[:-1] == e[1:], axis=1))
    return e[k, 0], e[k, 1], owner[k], owner[k + 1]


def improve_patch(
    nodes,
    elements,
    movable,
    mutable_faces,
    *,
    slidable=None,
    slide_on=None,
    min_angle_deg: float = 30.0,
    max_angle_deg: float = 130.0,
    max_area_change: float = 0.5,
    max_valence: int = 8,
    only_below: float = 1.15,
    soft: bool = False,
    rounds: int = 30,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Repair the seam without touching anything the patch does not own.

    The seam is where a fill pinned to an immutable rim meets a mesh that was
    already finished against the same gates: the base mesh sits at min angle
    30.01 deg and area change 0.500 exactly, so there is no margin to absorb a
    new neighbour, and the first patch came out at 26.1 deg with seven area
    jumps over the limit.  Those elements are all patch elements -- the
    retained mesh is untouched -- so they can be repaired in place.

    Two operations, both restricted:

    * an **edge flip**, allowed only when both incident faces are mutable.  A
      flip rewrites connectivity, and rewriting a retained face would break
      the contract even if every coordinate stayed put.
    * a **node move**, allowed only for nodes in ``movable``, and for nodes in
      ``slidable`` only along ``slide_on``, the coastline curve those nodes
      were cut from.  This matters because the worst elements in the first
      patch had two of their three vertices on the coast and no interior
      vertex left to move.  Sliding along the mesh's own polyline instead is
      the version that looks equivalent and is not: each node stays on the
      old line, but the chords between them cut the corners, and on this
      patch that moved the coastline 59.9 m.  Projecting every candidate back
      onto the source curve leaves the coastline exactly as faithful as
      ``coastline_tolerance_m`` already made it.

    The candidate positions are not the Laplacian centroid.  The fill's own
    smoother has already driven every interior node there, so on this mesh a
    Laplacian sweep proposed a zero-length move for every node and the pass
    did nothing at all; the offending node was 0.000 m from its centroid and a
    15 m step off it raised the neighbourhood health from 0.664 to 0.705.
    What is searched instead is a ring of directions at a few fractions of the
    local edge length.

    Every candidate is scored by :func:`_health` and accepted only when it
    strictly improves the worst gate in its neighbourhood, so the pass is
    monotone and cannot trade one failure for another.  Neighbourhoods already
    clear of every gate by ``only_below`` are skipped: most of a patch is
    fine, DistMesh made it that way on purpose, and touching it would be both
    slow and gratuitous.  It is not guaranteed
    to reach the gates -- when it cannot, that is a fact about the cut, and
    the caller should widen the transition rather than lower the gate.
    """
    xy = np.array(nodes, dtype=float)[:, :2]
    tri = np.array(elements, dtype=np.int64)
    mutable = np.asarray(mutable_faces, dtype=bool)
    lines = _slide_lines(slide_on)
    can_slide = (np.asarray(slidable, dtype=bool) if slidable is not None
                 and lines else np.zeros(len(xy), dtype=bool))
    slide = _boundary_neighbours(xy, tri, len(xy), can_slide) \
        if can_slide.any() else {}
    # Each sliding node is bound to ONE curve, here and for good.  Choosing it
    # per candidate would let a node hop between curves, and locating along a
    # MultiLineString measures the concatenation, so a projection can land on
    # another island with nothing raised.  A node is also bound to the SPAN
    # between the two curve vertices that bracket it, and a node sitting on a
    # curve vertex does not slide at all: staying on the curve is not the same
    # as leaving it where it was.  A node slid past a corner is still exactly
    # on the curve and the polyline has lost the corner -- measured on a
    # single perturbed grid node, 0 m off the curve and 52.5 m of Hausdorff
    # movement in the boundary itself (review finding 10, reconfirmed
    # 2026-09-22).  Since `preserve` keeps every original vertex as a node,
    # pinning the vertices and confining the rest to their own span makes the
    # polyline exactly invariant.
    curve_of = {}
    for v in slide:
        ln = min(lines, key=lambda c, q=xy[v]: c.distance(shapely_point(q)))
        span = _vertex_span(ln, xy[v])
        if span is not None:
            curve_of[v] = (ln, *span)
    slide = {v: w for v, w in slide.items() if v in curve_of}
    # A pinned node is not merely un-slidable, it is immovable: leaving it in
    # `movable` sends it down the free-ring branch of _candidates and off the
    # coastline altogether, which is worse than the corner-cutting this pins
    # it to prevent.
    allowed = np.zeros(len(xy), dtype=bool)
    allowed[list(slide)] = True
    movable = np.flatnonzero((np.asarray(movable, dtype=bool) & ~can_slide)
                             | allowed)
    n_flip = n_move = 0

    for _ in range(rounds):
        changed = False

        # --- flips -------------------------------------------------------
        ea, eb, f0s, f1s = _interior_edges(tri)
        both = mutable[f0s] & mutable[f1s]
        inc, lo, hi = _incidence(tri, len(xy))
        adj = _face_adjacency(tri)
        # Valence is the incident-FACE count over the whole mesh, matching
        # qa.py's C5 gate.  Counting it over the candidate's own fan gave
        # [5,5,7,7] where the mesh had [5,5,7,9] and let a flip through
        # against max_valence=8 (review finding 4, 2026-09-22).
        val = np.bincount(tri.ravel(), minlength=len(xy))
        for a, b, f0, f1 in zip(ea[both], eb[both], f0s[both], f1s[both]):
            a, b, f0, f1 = int(a), int(b), int(f0), int(f1)
            if not ({a, b} <= set(tri[f0].tolist())
                    and {a, b} <= set(tri[f1].tolist())):
                # An earlier accepted flip in this sweep rewrote one of these
                # faces, so the edge list no longer describes the mesh.
                continue
            c = int(np.setdiff1d(tri[f0], [a, b])[0])
            d = int(np.setdiff1d(tri[f1], [a, b])[0])
            if not _convex_quad(xy, a, b, c, d):
                # A concave quad's flip overlaps itself, and both halves can
                # still come out positively oriented -- an area test alone
                # would wave it through.
                continue
            # Valence is tracked over the whole mesh, not over the
            # candidate's own fan -- counting the fan gave [5,5,7,7] where
            # the mesh had [5,5,7,9] and let a flip through against a limit
            # of 8 (review finding 4, 2026-09-22).  It is a cost here rather
            # It is a veto, so the gate cannot be breached at any point.
            # Allowing an intermediate excess does buy reach -- the greedy
            # pass can route through it -- but it leaves nodes at 9 that no
            # later flip can bring down, and on this patch `preserve` failed
            # C5 at every seed that way.  A guarantee beats a heuristic.
            if val[c] + 1 > max_valence or val[d] + 1 > max_valence:
                continue
            faces = np.unique(np.concatenate(
                [inc[lo[v]:hi[v]] for v in (a, b, c, d)]))
            before = _scores(xy, tri, faces, adj, min_angle_deg, max_angle_deg,
                             max_area_change)
            if before[0] >= only_below:
                continue
            keep0, keep1 = tri[f0].copy(), tri[f1].copy()
            tri[f0] = _ccw(xy, np.array([c, d, b]))
            tri[f1] = _ccw(xy, np.array([d, c, a]))
            after = _scores(xy, tri, faces, _face_adjacency(tri),
                            min_angle_deg, max_angle_deg, max_area_change)
            if _better(after, before, soft):
                n_flip += 1
                changed = True
                # Accepted: every precomputed table is now stale for the
                # vertices this touched, so they are rebuilt.  Flips are few
                # (6-71 on the Futtsu patch) and a stale fan is how the
                # "strictly improving" claim stopped being true.
                val[a] -= 1
                val[b] -= 1
                val[c] += 1
                val[d] += 1
                inc, lo, hi = _incidence(tri, len(xy))
                adj = _face_adjacency(tri)
            else:
                tri[f0], tri[f1] = keep0, keep1

        # --- moves -------------------------------------------------------
        # Connectivity is unchanged by a move, so one incidence table serves
        # the whole sweep.
        inc, lo, hi = _incidence(tri, len(xy))
        adj = _face_adjacency(tri)
        for v in movable:
            faces = inc[lo[v]:hi[v]]
            if not len(faces):
                continue
            before = _scores(xy, tri, faces, adj, min_angle_deg, max_angle_deg,
                             max_area_change)
            if before[0] >= only_below:
                continue
            keep = xy[v].copy()
            scale = float(np.linalg.norm(
                xy[tri[faces]].reshape(-1, 2) - keep, axis=1).mean())
            best, best_p = before, None
            for q in _candidates(xy, tri, faces, v, keep, scale,
                                 slide.get(int(v)), curve_of.get(int(v))):
                xy[v] = q
                sc = _scores(xy, tri, faces, adj, min_angle_deg, max_angle_deg,
                             max_area_change)
                if _better(sc, best, soft):
                    best, best_p = sc, q
            xy[v] = keep if best_p is None else best_p
            if best_p is not None:
                n_move += 1
                changed = True

        if not changed:
            break

    # Valence cleanup.  Intermediate excess is allowed above because
    # forbidding it blocks the sequences that end below the limit -- but
    # ending above it is a C5 failure, so any node still over the gate gets
    # one more round of flips whose only job is to bring it down, taken
    # whenever they do not cost anything the other gates measure.
    n_valence_fixed = 0
    for _ in range(rounds):
        val = np.bincount(tri.ravel(), minlength=len(xy))
        over = np.flatnonzero(val > max_valence)
        if not over.size:
            break
        ea, eb, f0s, f1s = _interior_edges(tri)
        adj = _face_adjacency(tri)
        inc, lo, hi = _incidence(tri, len(xy))
        fixed_any = False
        for a, b, f0, f1 in zip(ea, eb, f0s, f1s):
            a, b, f0, f1 = int(a), int(b), int(f0), int(f1)
            if not (mutable[f0] and mutable[f1]):
                continue
            if val[a] <= max_valence and val[b] <= max_valence:
                continue
            if not ({a, b} <= set(tri[f0].tolist())
                    and {a, b} <= set(tri[f1].tolist())):
                continue
            c = int(np.setdiff1d(tri[f0], [a, b])[0])
            d = int(np.setdiff1d(tri[f1], [a, b])[0])
            if not _convex_quad(xy, a, b, c, d):
                continue
            if val[c] + 1 > max_valence or val[d] + 1 > max_valence:
                continue
            faces = np.unique(np.concatenate(
                [inc[lo[v]:hi[v]] for v in (a, b, c, d)]))
            before = _scores(xy, tri, faces, adj, min_angle_deg, max_angle_deg,
                             max_area_change)
            keep0, keep1 = tri[f0].copy(), tri[f1].copy()
            tri[f0] = _ccw(xy, np.array([c, d, b]))
            tri[f1] = _ccw(xy, np.array([d, c, a]))
            after = _scores(xy, tri, faces, _face_adjacency(tri), min_angle_deg,
                            max_angle_deg, max_area_change)[0]
            # Spending margin is allowed here, failing a gate is not. The main
            # pass may not lose ground anywhere; this one has a gate to
            # satisfy -- a node above max_valence -- and a flip that takes a
            # comfortable neighbourhood from 1.5 to 1.2 still passes
            # everything. Refusing it only to keep the margin is how a C5
            # violation survived every seed.
            if after >= min(before[0], 1.0) - 1e-12:
                n_valence_fixed += 1
                n_flip += 1
                fixed_any = True
                val[a] -= 1
                val[b] -= 1
                val[c] += 1
                val[d] += 1
                inc, lo, hi = _incidence(tri, len(xy))
                adj = _face_adjacency(tri)
            else:
                tri[f0], tri[f1] = keep0, keep1
        if not fixed_any:
            break

    ang = _angles_deg(xy, tri)
    return xy, tri, {
        "n_flips": n_flip,
        "n_moves": n_move,
        "n_slidable_used": len(slide),
        "n_valence_flips": n_valence_fixed,
        "max_valence": int(np.bincount(tri.ravel(), minlength=len(xy)).max()),
        "min_angle_deg": float(ang.min()),
        "max_angle_deg": float(ang.max()),
    }


def _vertex_span(line, q, eps: float = 1e-6):
    """The arc-length window a point may slide in without crossing a vertex.

    ``None`` when the point IS a vertex of the curve, which is the case that
    must not move: every corner of a `preserve` coastline is a mesh node, so
    pinning them is what keeps the polyline identical rather than merely
    keeping each node on it.
    """
    import shapely

    coords = np.asarray(line.coords, dtype=float)[:, :2]
    stations = np.concatenate(
        [[0.0], np.cumsum(np.linalg.norm(np.diff(coords, axis=0), axis=1))])
    s0 = float(shapely.line_locate_point(line, shapely_point(q)))
    tol = max(eps, 1e-9 * stations[-1])
    if np.abs(stations - s0).min() <= tol:
        return None
    lo = float(stations[stations < s0].max())
    hi = float(stations[stations > s0].min())
    return lo, hi


def shapely_point(q):
    import shapely

    return shapely.Point(float(q[0]), float(q[1]))


def _slide_lines(slide_on) -> list:
    """Normalise the slide geometry to a list of LineStrings."""
    if slide_on is None:
        return []
    import shapely

    items = slide_on if isinstance(slide_on, (list, tuple)) else [slide_on]
    out = []
    for it in items:
        if it is None:
            continue
        g = it if hasattr(it, "geom_type") else shapely.LineString(
            np.asarray(it, dtype=float)[:, :2])
        out.extend(list(g.geoms) if g.geom_type == "MultiLineString" else [g])
    return out


def _candidates(xy, tri, faces, v, keep, scale, slide_pair, curve):
    """Where a node might go: along its coastline, or a ring of small steps."""
    if slide_pair is not None and curve is not None:
        import shapely

        line, lo, hi = curve
        a, b = slide_pair
        for other in (a, b):
            for f in (0.05, 0.12, 0.25, 0.4):
                q = shapely_point(keep + f * (xy[other] - keep))
                station = float(shapely.line_locate_point(line, q))
                # Clamped inside its own span: a slide re-spaces the
                # coastline, it does not reshape it.
                station = min(max(station, lo), hi)
                yield np.asarray(shapely.line_interpolate_point(
                    line, station).coords[0])
        return
    ring = tri[faces].ravel()
    yield xy[ring[ring != v]].mean(axis=0)
    for ang in np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False):
        step = np.array([np.cos(ang), np.sin(ang)])
        for f in (0.03, 0.08, 0.16, 0.3):
            yield keep + f * scale * step


def _boundary_neighbours(xy, tri, n_nodes, slidable) -> dict[int, tuple[int, int]]:
    """The two boundary neighbours of each slidable node.

    They set the DIRECTION and length scale of a trial slide; where the node
    actually lands is decided by projecting that trial back onto the source
    coastline, so these two only have to bracket the node along the chain.
    """
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    u, c = np.unique(e, axis=0, return_counts=True)
    b = u[c == 1]
    nbr: dict[int, list[int]] = {}
    for x, y in b.tolist():
        nbr.setdefault(x, []).append(y)
        nbr.setdefault(y, []).append(x)
    out: dict[int, tuple[int, int]] = {}
    for v, w in nbr.items():
        if len(w) != 2 or v >= n_nodes or not slidable[v]:
            continue
        out[v] = (w[0], w[1])
    return out


def _convex_quad(xy, a, b, c, d) -> bool:
    """Is a-c-b-d convex?  Only then does flipping a-b to c-d tile the quad."""
    e = xy[d] - xy[c]
    sa = e[0] * (xy[a][1] - xy[c][1]) - e[1] * (xy[a][0] - xy[c][0])
    sb = e[0] * (xy[b][1] - xy[c][1]) - e[1] * (xy[b][0] - xy[c][0])
    return bool(sa * sb < 0)


def _ccw(xy, t: np.ndarray) -> np.ndarray:
    """Vertex order with positive area; fort.14 wants counter-clockwise."""
    u, v = xy[t[1]] - xy[t[0]], xy[t[2]] - xy[t[0]]
    return t if (u[0] * v[1] - u[1] * v[0]) > 0 else t[[0, 2, 1]]


