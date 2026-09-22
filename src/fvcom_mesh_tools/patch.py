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
    "ambient_size_field",
    "boundary_rings",
    "coastline_points",
    "hole_polygon",
    "improve_patch",
    "effective_gradation",
    "patch_sizing",
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
    for _ in range(max_repair_rounds):
        sel = tri[removed]
        u, c = _edge_table(sel)
        rim = u[c == 1]
        deg = _rim_degree(rim, n_nodes)
        bad = np.where(deg > 2)[0]
        if not bad.size:
            break
        touch = np.isin(tri, bad).any(axis=1)
        add = touch & ~removed
        if not add.any():
            raise ValueError(
                f"the cut pinches at {bad.size} vertices and taking their faces "
                "does not resolve it; the footprint straddles a one-element-wide "
                "neck")
        removed |= add
        grown += int(add.sum())
    else:
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

    Edge connectivity, not vertex: two triangles meeting at a single point are
    a non-manifold mesh, not a connected one, and FVCOM will not run on it.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    if not len(retained):
        raise ValueError("the cut leaves no elements")
    e = np.vstack([retained[:, [0, 1]], retained[:, [1, 2]], retained[:, [2, 0]]])
    n = int(frozen.max()) + 1
    g = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n, n))
    # Nodes no retained element uses form their own singleton components, so
    # only the labels of the used nodes are counted.
    lab = connected_components(g, directed=False)[1][frozen]
    if np.unique(lab).size != 1:
        raise ValueError(
            f"the cut splits the retained mesh into {np.unique(lab).size} pieces; "
            "a refinement may not disconnect the domain")


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

    ``shoreline`` must be a single LineString -- the caller picks the ring the
    stretch belongs to, because picking it here by nearest distance would
    silently jump rings at a strait.  ``None`` when there is no usable piece.
    """
    import shapely
    from shapely.ops import substring

    if shoreline is None:
        return None
    line = shapely.LineString(np.asarray(shoreline.coords, dtype=float)[:, :2]) \
        if hasattr(shoreline, "coords") else shoreline
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


def _resample_on_source(pts: np.ndarray, shoreline, size) -> np.ndarray:
    """Re-space one stretch along the source shoreline between its endpoints."""
    coords = _source_substring(pts, shoreline)
    if coords is None:
        return _subdivide(pts, size)
    return _walk(coords, size)


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
            # stretch to anchor a resample between.  Kept as it stands.
            ring_rows = [_push(pts, base_id, xy[v], int(v)) for v in ring]
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
    """
    import shapely
    from shapely.ops import polygonize, unary_union

    lines = [shapely.LineString([pfix[a], pfix[b]]) for a, b in egfix.tolist()]
    polys = list(polygonize(unary_union(lines)))
    if not polys:
        raise ValueError("the rim segments do not close a polygon")
    return unary_union(polys)


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

    ``regions`` is a list of ``(geometry_in_mesh_crs, target_h_m, width_m)``.
    Inside a region the size is its target; over the next ``width_m`` it ramps
    linearly to **the base mesh's own size at that point**, and beyond it is
    the base mesh's size.  It never exceeds the base size, so a refinement can
    only refine.

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
    import shapely
    from matplotlib.tri import LinearTriInterpolator, Triangulation

    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    amb = ambient_size_field(xy, tri)
    interp = LinearTriInterpolator(Triangulation(xy[:, 0], xy[:, 1], tri), amb)
    amb_max = float(np.nanmax(amb))
    geoms = [(shapely.boundary(g) if g.geom_type in ("Polygon", "MultiPolygon") else g,
              g, float(h), float(w)) for g, h, w in regions]

    def h(points):
        p = np.atleast_2d(np.asarray(points, dtype=float))[:, :2]
        pt = shapely.points(p[:, 0], p[:, 1])
        base = np.asarray(interp(p[:, 0], p[:, 1]), dtype=float)
        base = np.where(np.isfinite(base), base, amb_max)
        out = base.copy()
        for edge, poly, target, width in geoms:
            d = shapely.distance(pt, edge)
            d = np.where(shapely.contains(poly, pt), 0.0, d)
            u = np.clip(d / width, 0.0, 1.0)
            out = np.minimum(out, target + (base - target) * u)
        return np.minimum(out, base) / distmesh_scale

    return h


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
    for geom, target, width in regions:
        import shapely

        d = shapely.distance(shapely.points(xy[:, 0], xy[:, 1]), geom)
        inside = d < width
        g = float(((amb[inside] - target) / width).max()) if inside.any() else 0.0
        per_region[str(target)] = g
        worst = max(worst, g)
    return {
        "max_effective_gradation": worst,
        "c4_limit_gradation": float(1.0 / np.sqrt(1.0 - 0.5) - 1.0),
        "within_c4": bool(worst <= 1.0 / np.sqrt(1.0 - 0.5) - 1.0),
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
        "n_nodes": int(len(nodes)),
        "n_elements": int(len(elements)),
        "n_nodes_retained": int(len(keep)),
        "n_nodes_new": int(len(fresh)),
        "n_nodes_deleted": int(len(xy) - len(keep)),
        "n_new_depths_outside_base": int(n_outside),
    }
    return nodes, elements, depths, node_map, report


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
    tol_m: float = 1e-6,
    open_boundaries=(),
) -> dict[str, Any]:
    """Check the frozen zone really is frozen, through the node map.

    Several separate claims, because "the mesh outside is unchanged" is
    several claims wearing one coat: the retained nodes kept their
    coordinates, they kept their depths, the retained faces are all still
    there with the same vertices, their orientation did not flip, the
    interface segments are still shared rather than split, and the open
    boundary is the same list in the same order.  A patch can satisfy all but
    one of those and still break the model.  ``area_change_fraction`` is
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
    moved = np.linalg.norm(new_xy[nm[keep]] - xy[keep], axis=1)
    ddep = np.abs(new_dep[nm[keep]] - dep0[keep])

    want = {tuple(sorted(r)) for r in nm[selection.retained].tolist()}
    have = {tuple(sorted(r)) for r in tri.tolist()}
    missing = want - have

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
        "n_frozen_moved": int((moved > tol_m).sum()),
        "max_frozen_depth_change_m": float(ddep.max()) if len(ddep) else 0.0,
        "n_retained_faces_missing": int(len(missing)),
        "n_inverted_elements": int((area2 <= 0).sum()),
        "n_duplicate_nodes": dup,
        "n_orphan_nodes": orphan,
        "n_interface_segments_split": int(len(split)),
        "open_boundary_unchanged": bool(obc_ok),
        "area_change_fraction": float((new_area - base_area) / base_area)
        if base_area > 0 else 0.0,
        "ok": bool((moved <= tol_m).all() and (ddep <= tol_m).all()
                   and not missing and (area2 > 0).all()
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
        return np.inf
    t = tri[faces]
    ar = _areas(xy, t)
    if (ar <= 0).any():
        return -np.inf
    ang = _angles_deg(xy, t)
    score = min(float(ang.min() / min_angle),
                float((180.0 - ang.max()) / (180.0 - max_angle)))
    nb = adj[faces]
    have = nb >= 0
    if have.any():
        mine = np.repeat(np.abs(ar), 3).reshape(-1, 3)[have]
        theirs = np.abs(_areas(xy, tri[nb[have]]))
        change = np.abs(mine - theirs) / np.maximum(mine, theirs)
        score = min(score, float(((1.0 - change) / (1.0 - max_area_change)).min()))
    return score


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
    movable = np.flatnonzero(np.asarray(movable, dtype=bool) | can_slide)
    # Each sliding node is bound to ONE curve, here and for good.  Choosing it
    # per candidate would let a node hop between curves, and locating along a
    # MultiLineString measures the concatenation, so a projection can land on
    # another island with nothing raised.
    curve_of = {v: min(lines, key=lambda ln, q=xy[v]: ln.distance(shapely_point(q)))
                for v in slide}
    n_flip = n_move = 0

    for _ in range(rounds):
        changed = False

        # --- flips -------------------------------------------------------
        ea, eb, f0s, f1s = _interior_edges(tri)
        both = mutable[f0s] & mutable[f1s]
        inc, lo, hi = _incidence(tri, len(xy))
        adj = _face_adjacency(tri)
        touched: set[int] = set()
        for a, b, f0, f1 in zip(ea[both], eb[both], f0s[both], f1s[both]):
            a, b, f0, f1 = int(a), int(b), int(f0), int(f1)
            if f0 in touched or f1 in touched:
                continue  # the precomputed tables would be stale
            c = int(np.setdiff1d(tri[f0], [a, b])[0])
            d = int(np.setdiff1d(tri[f1], [a, b])[0])
            if not _convex_quad(xy, a, b, c, d):
                # A concave quad's flip overlaps itself, and both halves can
                # still come out positively oriented -- an area test alone
                # would wave it through.
                continue
            faces = np.unique(np.concatenate(
                [inc[lo[v]:hi[v]] for v in (a, b, c, d)]))
            before = _health(xy, tri, faces, adj, min_angle_deg, max_angle_deg,
                             max_area_change)
            if before >= only_below:
                continue
            keep0, keep1 = tri[f0].copy(), tri[f1].copy()
            tri[f0] = _ccw(xy, np.array([c, d, b]))
            tri[f1] = _ccw(xy, np.array([d, c, a]))
            if (_health(xy, tri, faces, _face_adjacency(tri), min_angle_deg,
                        max_angle_deg, max_area_change) > before + 1e-9
                    and _valence_ok(tri, faces, (a, b, c, d), max_valence)):
                n_flip += 1
                changed = True
                touched.add(f0)
                touched.add(f1)
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
            before = _health(xy, tri, faces, adj, min_angle_deg, max_angle_deg,
                             max_area_change)
            if before >= only_below:
                continue
            keep = xy[v].copy()
            scale = float(np.linalg.norm(
                xy[tri[faces]].reshape(-1, 2) - keep, axis=1).mean())
            best, best_p = before, None
            for q in _candidates(xy, tri, faces, v, keep, scale,
                                 slide.get(int(v)), curve_of.get(int(v))):
                xy[v] = q
                sc = _health(xy, tri, faces, adj, min_angle_deg, max_angle_deg,
                             max_area_change)
                if sc > best + 1e-9:
                    best, best_p = sc, q
            xy[v] = keep if best_p is None else best_p
            if best_p is not None:
                n_move += 1
                changed = True

        if not changed:
            break

    ang = _angles_deg(xy, tri)
    return xy, tri, {
        "n_flips": n_flip,
        "n_moves": n_move,
        "n_slidable_used": len(slide),
        "min_angle_deg": float(ang.min()),
        "max_angle_deg": float(ang.max()),
    }


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

        a, b = slide_pair
        for other in (a, b):
            for f in (0.05, 0.12, 0.25, 0.4):
                q = shapely_point(keep + f * (xy[other] - keep))
                yield np.asarray(shapely.line_interpolate_point(
                    curve, shapely.line_locate_point(curve, q)).coords[0])
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


def _valence_ok(tri, faces, verts, max_valence: int) -> bool:
    """Valence is the incident-face count, matching the C5 gate in qa.py."""
    sub = tri[faces]
    for v in verts:
        if int((sub == v).sum()) > max_valence:
            return False
    return True
