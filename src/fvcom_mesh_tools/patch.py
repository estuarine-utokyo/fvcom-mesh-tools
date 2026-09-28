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

    A cut can also leave the retained mesh in pieces: a circle across the
    canals of a base that draws them (Odaiba on this project's own base)
    parts the water beyond it from the rest.  Every piece borders the hole,
    so the fill joins them again, and they are kept as they are -- the
    transition serves efficiency, not detail (owner, 2026-09-28).  Taking
    them into the hole instead re-cut a canal 2.7 km from the region at the
    transition's size, closed it, and left a land spit the mesh could not
    carry.  ``n_retained_pieces`` reports them; the refined mesh is checked
    to be one piece after stitching.

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

    n_pieces = _retained_pieces(retained, rim)

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
        "n_retained_pieces": n_pieces,
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


def _face_components(faces: np.ndarray) -> np.ndarray:
    """Component label of each face, faces joined by shared edges."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    e = np.sort(np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(len(faces)), 3)
    order = np.lexsort((e[:, 1], e[:, 0]))
    e, owner = e[order], owner[order]
    k = np.flatnonzero(np.all(e[:-1] == e[1:], axis=1))
    n = len(faces)
    g = coo_matrix((np.ones(len(k)), (owner[k], owner[k + 1])), shape=(n, n))
    return connected_components(g, directed=False)[1]


def _retained_pieces(retained: np.ndarray, rim: np.ndarray) -> int:
    """How many edge-connected pieces the retained mesh is in; each must
    border the hole, which the fill then joins them through.

    A piece touching the hole only at a vertex is a pinch, which the repair
    has already resolved; a piece with no rim edge at all would stay apart.
    """
    if not len(retained):
        raise ValueError("the cut leaves no elements")
    lab = _face_components(retained)
    n = int(lab.max()) + 1
    if n == 1:
        return 1
    rim_set = {tuple(x) for x in np.sort(np.asarray(rim), axis=1).tolist()}
    e = np.sort(np.vstack([retained[:, [0, 1]], retained[:, [1, 2]],
                           retained[:, [2, 0]]]), axis=1)
    on_rim = np.array([tuple(x) in rim_set for x in e.tolist()]).reshape(3, -1).any(axis=0)
    touching = set(lab[on_rim].tolist())
    if len(touching) != n:
        raise ValueError(
            f"the cut splits the retained mesh into {n} pieces, "
            f"{n - len(touching)} of them not on the hole; a refinement may not "
            "disconnect the domain")
    return n


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
    fine_h: float | None = None,
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
    if mode not in ("preserve", "resample", "spline", "resolve"):
        raise ValueError(f"unknown coastline mode {mode!r}")
    base = shapely.LineString(pts)

    if mode == "resolve":
        if shoreline is None:
            raise ValueError("coastline: resolve needs the source shoreline")
        out = _resolve_stretch(pts, shoreline, size, fine_h=fine_h)
    elif mode == "resample":
        if shoreline is None:
            raise ValueError("coastline: resample needs the source shoreline")
        out = _resample_on_source(pts, shoreline, size)
    elif mode == "spline":
        out = _spline_resample(pts, size)
    else:
        out = _subdivide(pts, size)

    # `resolve` is the hires branch, and there the coastline is SUPPOSED to
    # move: the base polyline runs 300-600 m between nodes and sits up to
    # 68.5 m from OSM at its segment midpoints, and recovering that is the
    # request.  The departure is measured and reported by the driver instead
    # (owner, 2026-09-23).  Every other mode keeps the veto unchanged.
    if mode == "resolve":
        out[0], out[-1] = pts[0], pts[-1]
        return out

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
    if mode in ("resample", "resolve"):
        piece = _source_substring(pts, shoreline)
        if piece is not None:
            return piece
        if mode == "resolve":
            # Silently returning the base here is how a fidelity check
            # certifies "zero departure from the source" without ever having
            # looked at the source (review P2-12, 2026-09-23).  On the hires
            # branch the instruction is to follow OSM, so a stretch OSM does
            # not cover is a refusal, not a fallback.
            raise ValueError(
                "coastline: resolve found no source component for a stretch; "
                "the shoreline does not cover it, or the match was rejected "
                "as running the wrong way round a ring")

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
    # Every vertex, and points evenly spaced by ARC LENGTH: spaced by vertex
    # index, 32 probes all fell in a dense 1 m cluster of a 3.9 km stretch
    # and the fit never saw the rest of it (review, round 7).
    probe = shapely.points(np.vstack([
        pts[:, :2],
        shapely.get_coordinates(shapely.line_interpolate_point(
            here, np.linspace(0.0, here.length, 65)))]))
    line = min(lines, key=lambda ln: float(shapely.distance(probe, ln).max()))
    del here
    s0 = line.project(shapely.Point(pts[0]))
    s1 = line.project(shapely.Point(pts[-1]))
    if s0 == s1:
        return None
    piece = substring(line, s0, s1)
    # A closed ring has two arcs between the endpoints, and the stretch may
    # be the long one.  At Yokohama a base island (the Daikoku pier) kept
    # one coast edge in the frozen zone; the rest of its coast was one
    # stretch from one end of that edge round to the other, and the direct
    # arc between them was that edge itself -- so the island vanished from
    # the rim and was meshed as water.  The arc that FITS the stretch is the
    # one it replaces; the longer one is taken only when it fits clearly
    # better, so a stretch that the shorter arc follows is unchanged.  Which
    # arc is "direct" depends on where the ring happens to start; shorter and
    # longer do not (review, round 7).
    if line.is_closed and line.length > 0:
        L = line.length
        if s0 < s1:
            other = [substring(line, s0, 0.0), substring(line, L, s1)]
        else:
            other = [substring(line, s0, L), substring(line, 0.0, s1)]
        oc = [np.asarray(g.coords, dtype=float)[:, :2] for g in other if g.length > 0]
        if oc:
            wrap = shapely.LineString(np.vstack([oc[0], *[c[1:] for c in oc[1:]]]))
            fit = {id(g): float(shapely.distance(probe, g).max()) for g in (piece, wrap)}
            if abs(piece.length - wrap.length) <= 1e-9 * L:
                # equal arcs: the one that fits, and on a tie the one whose
                # midpoint comes first -- never "whichever did not cross the
                # ring's start" (review, round 8)
                def mid(g):
                    return tuple(np.round(g.interpolate(0.5, normalized=True).coords[0], 6))
                piece = min((piece, wrap), key=lambda g: (round(fit[id(g)], 9), mid(g)))
            else:
                short, long_ = (piece, wrap) if piece.length < wrap.length else (wrap, piece)
                piece = long_ if fit[id(long_)] < 0.5 * fit[id(short)] else short
    coords = np.asarray(piece.coords, dtype=float)[:, :2]
    if len(coords) < 2:
        return None
    # Walking the wrong way round a closed ring gives a piece far longer than
    # the stretch it replaces; fall back rather than swap in half a coastline.
    if piece.length > 3.0 * shapely.LineString(pts).length:
        return None
    return coords


def _unusable_replacement(new: np.ndarray, xy: np.ndarray, idx,
                          boundary_edges, placed, where: list | None = None) -> str | None:
    """Why this replacement cannot be used, or ``None`` when it can.

    Three ways a moved boundary breaks the mesh, and all three were met on
    real patches:

    ``self``      the new polyline touches or crosses ITSELF.  A pier 20 m
                  wide resolved at 30 m does this: the walk goes up one side
                  and back down the other and the two sides interleave.  The
                  ring is then not simple, and `hole_polygon` fails inside
                  GEOS with "side location conflict", which is a sentence
                  about topology and not about the port.
    ``other``     it crosses a stretch already placed on this rim.
    ``frozen``    it crosses the base mesh's boundary somewhere it does not
                  own -- a shared endpoint is a touch, not a crossing, so the
                  edges this stretch REPLACES are excluded.

    Subdividing cannot do any of the three, because it stays on a polyline
    that was already part of a valid mesh.
    """
    if len(new) < 2:
        return None
    import shapely

    line = shapely.LineString(np.asarray(new, dtype=float)[:, :2])
    if not line.is_simple:
        return "self"
    if placed:
        tree = shapely.STRtree(placed)
        for q in tree.query(line):
            if shapely.crosses(line, placed[q]) or shapely.overlaps(line, placed[q]):
                return "other"
    if boundary_edges is None:
        return None
    e = np.asarray(boundary_edges, dtype=np.int64).reshape(-1, 2)
    mine = set(np.asarray(idx, dtype=np.int64).ravel().tolist())
    others = e[[not (i in mine and j in mine) for i, j in e.tolist()]]
    if not len(others):
        return None
    segs = [shapely.LineString(xy[[i, j], :2]) for i, j in others.tolist()]
    tree = shapely.STRtree(segs)
    hit = [q for q in tree.query(line) if shapely.crosses(line, segs[q])]
    if hit:
        if where is not None:
            # where, so a caller can take the ground there into the hole
            where.extend(np.asarray(shapely.intersection(line, segs[q])
                                    .representative_point().coords)[0].tolist()
                         for q in hit)
        return "frozen"
    return None


def filter_shoreline(land, h0: float, *, elements_per_feature: float = 3.0,
                     land_width_factor: float | None = None):
    """Remove what a mesh of size ``h0`` cannot resolve, and say what went.

    This is the judgement the declared grid size implies, and it has to be
    made explicitly because nothing else makes it.  `oceanmesh.Shoreline`
    culls by AREA -- islands under ``minimum_area_mult * h0**2`` -- which does
    not touch a long thin one: the Kimitsu pier is about 20 m across and
    700 m long, 14,000 m2 against a 3,600 m2 threshold, and a 30 m mesh
    cannot carry it whatever its area.

    Width is what matters, and a morphological opening and closing is what
    measures it.  ``buffer(-r).buffer(r)`` deletes land narrower than ``2r``;
    the same pair the other way round deletes water narrower than ``2r``.

    ``elements_per_feature`` is how many elements a surviving feature must
    have room for, so ``r = elements_per_feature * h0 / 2``.  It is not 1:
    a channel exactly ``h0`` wide has a node on each bank and no element
    between them, and the first run of this filter at ``r = h0/2`` collapsed
    43 of 262 fixed points onto 21 vertices, in pairs 30 m apart at a 30 m
    target.  Three is what OceanMesh2D's feature sizing asks for across a
    channel, and it is the smallest number that leaves an element with
    neighbours on both sides.

    ``land_width_factor``, when given, sets the LAND threshold on its own:
    land narrower than ``land_width_factor * h0`` goes, while water keeps
    ``elements_per_feature``.  The mesh lives in the water, so a narrow pier
    can still be meshed round as land: its outline is a coastline, and the
    only thing its width limits is the edge across its end, which C1 needs
    at half an element or more (owner, 2026-09-24).  ``None`` keeps one
    threshold for both, as before.

    Returns ``(filtered, report)``.  The report is the point: a coastline
    that quietly lost its piers is worse than one that says it did.
    """
    import shapely

    if not (np.isfinite(h0) and h0 > 0):
        raise ValueError("h0 must be finite and positive")
    if not (np.isfinite(elements_per_feature) and elements_per_feature > 0):
        raise ValueError("elements_per_feature must be finite and positive")
    if land_width_factor is not None and not (
            np.isfinite(land_width_factor) and land_width_factor > 0):
        raise ValueError("land_width_factor must be finite and positive")
    r = 0.5 * float(elements_per_feature) * float(h0)
    r_land = r if land_width_factor is None else 0.5 * float(land_width_factor) * float(h0)
    before = shapely.union_all(
        [land] if hasattr(land, "geom_type") else list(land))
    # MITRE joins: a port is rectilinear, and round joins shave every convex
    # corner of a quay into a crescent -- land lost that was never narrow,
    # and a sliver that the wall extraction would read as a structure.
    opened = before.buffer(-r_land, join_style="mitre").buffer(r_land, join_style="mitre")
    closed = opened.buffer(r, join_style="mitre").buffer(-r, join_style="mitre")
    out = shapely.make_valid(closed)

    def _rings(g):
        return sum(1 + len(q.interiors) for q in getattr(g, "geoms", [g])
                   if not q.is_empty and q.geom_type == "Polygon")

    report = {
        "h0_m": float(h0),
        "elements_per_feature": float(elements_per_feature),
        "removes_features_narrower_than_m": float(2.0 * r),
        "removes_land_narrower_than_m": float(2.0 * r_land),
        "fills_water_narrower_than_m": float(2.0 * r),
        "area_before_m2": float(before.area),
        "area_after_m2": float(out.area),
        "area_removed_m2": float(before.area - out.area),
        "area_removed_fraction": float((before.area - out.area) / before.area)
        if before.area > 0 else 0.0,
        "rings_before": _rings(before),
        "rings_after": _rings(out),
        "land_lost_m2": float(shapely.difference(before, out).area),
        "water_lost_m2": float(shapely.difference(out, before).area),
        "perimeter_before_m": float(before.length),
        "perimeter_after_m": float(out.length),
    }
    return out, report


def unresolvable_water(land, size_field, footprint, *, radius_factor=1.0,
                       min_h=0.0, spacing=10.0, levels_ratio=1.1):
    """The water no disc of the LOCAL element size reaches, as polygons.

    Water is kept where it is covered by a disc of radius
    ``radius_factor * h(c)`` whose centre ``c`` is at least that far from
    land; the rest is returned.  h is the local size itself, not the lower
    bound of an octave band, so a channel is judged at the size of the
    elements it will get and nothing is cut where two bands meet: the octave
    bands leave water two lower bounds wide -- ONE element where the elements
    are twice the bound -- and a straight cut where the next band closes it
    (Funabashi, Yokohama; USER_GUIDE section 11, band seams).

    Only water that touches ONE body of land is returned: water between two
    is a strait, and closing it would join them.

    Computed on a raster of cell ``spacing`` over ``footprint``: the distance
    to land, the centres, and the discs re-grown from them in radius levels
    ``levels_ratio`` apart (each level grown at its UPPER radius, so water is
    kept rather than lost to the quantisation).  Only where ``h > min_h``.
    The pieces are smoothed by one cell and simplified, so the coast they
    make is not a staircase; pieces smaller than a disc of ``min_h`` go.
    Returns ``(polygon, report)``.
    """
    import shapely
    from rasterio import features
    from rasterio.transform import from_origin
    from scipy import ndimage

    s = float(spacing)
    if not (np.isfinite(s) and s > 0):
        raise ValueError("spacing must be finite and positive")
    if not (np.isfinite(radius_factor) and radius_factor > 0):
        raise ValueError("radius_factor must be finite and positive")
    if not (np.isfinite(levels_ratio) and levels_ratio > 1.0):
        raise ValueError("levels_ratio must be finite and above 1")
    empty = shapely.Polygon()
    rep = {"spacing_m": s, "n_levels": 0, "lost_cells": 0, "area_m2": 0.0,
           "n_pieces": 0, "n_straits_left_open": 0, "n_entrances_left_open": 0,
           "n_detached_left_open": 0, "unbounded_pad": False}
    land_u = shapely.union_all([land] if hasattr(land, "geom_type") else list(land))
    # no footprint, or no coast to measure from: nothing can be judged
    if footprint is None or footprint.is_empty or land_u.is_empty:
        return empty, rep

    def sizes(x, y):
        v = np.asarray(size_field(np.column_stack([x.ravel(), y.ravel()])),
                       dtype=float).reshape(x.shape)
        if not np.isfinite(v).all() or (v <= 0).any():
            raise ValueError("the size field is not finite and positive over the footprint")
        return v

    # The raster reaches past the footprint by the largest disc, so a disc
    # centred outside it still counts (review, round 1): the footprint is
    # where water may be CLOSED, not where it may be measured from.  The
    # largest disc is the one over the PADDED box, not over the footprint --
    # sizes grow outward, and a disc centred 86 m out reached water the
    # footprint's own sizes said nothing could (review, round 6) -- so the
    # pad grows until it holds the largest disc centred inside it (a proved
    # upper bound, size_upper_bound).  Beyond it a disc is smaller than its
    # distance to the footprint while radius_factor times the field's slope
    # stays under 1 (the ramps: 0.75 x 0.165).  A pad that does not settle
    # proves nothing, and then nothing is closed.
    fx0, fy0, fx1, fy1 = footprint.bounds
    h_step = max(4.0 * s, float(min_h))

    def disc_over(p_):
        box_ = shapely.box(fx0 - p_, fy0 - p_, fx1 + p_, fy1 + p_)
        try:
            return radius_factor * size_upper_bound(box_, size_field, h_step) + 2.0 * s
        except ValueError as e:
            raise ValueError("the size field is not finite and positive over the "
                             "footprint") from e

    # A field that states its global maximum (size_max: patch_sizing,
    # base_size_field) is checked outward in annuli doubling to that
    # maximum's disc: a centre beyond the pad P is P from the footprint, so
    # the discs up to 2P must be under P, and so on.  A settled pad alone
    # proved nothing for a field that steps up outside it (review, round 7:
    # 100 -> 40 -> 200 m).  A field that brings its own bounds (size_bounds)
    # makes no slope promise, so without size_max nothing is closed; a plain
    # callable has the slope premise (radius_factor x slope < 1).
    r_max = getattr(size_field, "size_max", None)
    r_max = None if r_max is None else radius_factor * float(r_max)
    if r_max is None and hasattr(size_field, "size_bounds"):
        rep["unbounded_pad"] = True
        return empty, rep
    pad, last, grew = 2.0 * s, 0.0, np.inf
    settled = False
    for _ in range(16):
        # a pad that grows as fast as it did the step before diverges; stop
        # before the box (and the raster) outgrows memory
        if (fx1 - fx0 + 2 * pad) * (fy1 - fy0 + 2 * pad) / (s * s) > 5e7:
            break
        need = disc_over(pad)
        if need <= pad:
            ring = pad
            while r_max is not None and ring < r_max and disc_over(2.0 * ring) <= ring:
                ring *= 2.0
            if r_max is None or ring >= r_max:
                settled = True
                break
            need = disc_over(2.0 * ring)        # an annulus out there reaches in
            grew, last = np.inf, 0.0
        elif need - last >= grew:
            break
        else:
            grew, last = need - last, need
        pad = 1.25 * need     # overshoot: a growing field never meets need exactly
    if not settled:
        rep["unbounded_pad"] = True
        return empty, rep
    # The raster stays where the footprint's own sizes put it, and a larger
    # pad adds whole cells around it: moving the origin re-quantises the land,
    # which moved every recipe's coast by 50-95 m for a pad that only needed
    # to ADD disc centres.
    sx, sy = np.meshgrid(np.linspace(fx0, fx1, 24), np.linspace(fy0, fy1, 24))
    pad0 = radius_factor * float(sizes(sx, sy).max()) + 2.0 * s
    pad = pad0 + np.ceil(max(0.0, pad - pad0) / s) * s
    x0, y0, x1, y1 = fx0 - pad, fy0 - pad, fx1 + pad, fy1 + pad
    nx, ny = int(np.ceil((x1 - x0) / s)) + 1, int(np.ceil((y1 - y0) / s)) + 1
    transform = from_origin(x0, y1, s, s)
    is_land = features.rasterize([(land_u, 1)], out_shape=(ny, nx), transform=transform,
                                 fill=0, all_touched=False, dtype="uint8").astype(bool)
    if not is_land.any():
        return empty, rep
    inside = features.rasterize([(footprint, 1)], out_shape=(ny, nx), transform=transform,
                                fill=0, dtype="uint8").astype(bool)
    gx = x0 + (np.arange(nx) + 0.5) * s
    gy = y1 - (np.arange(ny) + 0.5) * s
    mx, my = np.meshgrid(gx, gy)
    h = sizes(mx, my)
    r = radius_factor * h
    sea = ~is_land                       # water anywhere on the raster
    water = inside & sea                 # water that may be closed
    if not water.any():                  # all land: nothing to close
        return empty, rep
    dist = ndimage.distance_transform_edt(sea) * s
    centres = sea & (dist >= r)
    covered = np.zeros_like(water)
    rmin, rmax = float(r[sea].min()), float(r[sea].max())
    lo = max(rmin, s)
    while lo <= rmax * levels_ratio:
        hi = lo * levels_ratio
        c = centres & (r >= lo) & (r < hi)
        if c.any():
            covered |= ndimage.distance_transform_edt(~c) * s <= hi
            rep["n_levels"] += 1
        lo = hi
    # centres whose radius is under one cell cover themselves
    covered |= centres & (r < max(rmin, s))
    lost = water & ~covered & (h > min_h)
    rep["lost_cells"] = int(lost.sum())
    if not lost.any():
        return empty, rep
    polys = [shapely.geometry.shape(g) for g, v in
             features.shapes(lost.astype("uint8"), mask=lost, transform=transform) if v == 1]
    g = shapely.union_all(polys)
    g = shapely.simplify(g.buffer(s, join_style="round").buffer(-s, join_style="round"), s)
    g = shapely.intersection(g, footprint)
    min_area = np.pi * max(min_h, s) ** 2 / 4.0
    pieces = [q for q in getattr(g, "geoms", [g]) if q.geom_type == "Polygon"
              and q.area >= min_area]
    # Only water that ends at ONE body of land: a dead end, a strip along a
    # bank.  Water between two bodies is a strait, and closing it joins them
    # -- the rim is resolved from the land's own components, and at Odaiba a
    # join left a stretch with no source component.
    comps = [q for q in getattr(land_u, "geoms", [land_u]) if not q.is_empty]
    tree = shapely.STRtree(comps)
    # ...and only water whose closing cuts nothing off: an entrance whose
    # basin behind is wide enough for elements touches one body of land too,
    # and closing it left the basin a lake (review, round 1)
    sea_v = shapely.difference(shapely.box(x0, y0, x1, y1), land_u)
    kept = []
    for q in pieces:
        ring = q.buffer(2.0 * s)
        touching = [k for k in tree.query(ring) if shapely.intersects(comps[k], ring)]
        if len(touching) > 1:
            rep["n_straits_left_open"] += 1
            continue
        if not touching:
            # water no land is next to is not a dead end or a strip -- only a
            # size field that jumps can leave it, and closing it would make
            # an island out of open water (review, round 2)
            rep["n_detached_left_open"] += 1
            continue
        # judged on the piece grown by two cells: smoothing leaves it a hair
        # short of the banks, and the hairline of water joined the two sides
        qd = q.buffer(2.0 * s)
        rest = shapely.difference(shapely.intersection(sea_v, qd.buffer(4.0 * s)), qd)
        sides = [w for w in getattr(rest, "geoms", [rest])
                 if w.geom_type == "Polygon" and w.area > s * s]
        if len(sides) > 1 and _splits_water(sea_v, qd, sides):
            rep["n_entrances_left_open"] += 1
            continue
        kept.append(q)
    out = shapely.union_all(kept) if kept else empty
    rep["area_m2"] = float(out.area)
    rep["n_pieces"] = len(kept)
    return out, rep


def _splits_water(sea, piece, sides) -> bool:
    """True when taking ``piece`` out of ``sea`` separates two of ``sides``."""
    import shapely

    left = shapely.difference(sea, piece)
    parts = [w for w in getattr(left, "geoms", [left]) if w.geom_type == "Polygon"]
    owner = set()
    for sd in sides:
        p = sd.representative_point()
        for k, w in enumerate(parts):
            if w.buffer(1e-6).contains(p):
                owner.add(k)
                break
    return len(owner) > 1


def filter_shoreline_local(land, size_field, h0: float, footprint, *,
                           elements_per_feature: float = 2.0,
                           land_width_factor: float | None = None,
                           land_width_max_band: int | None = None,
                           spacing: float | None = None,
                           keep_land=None,
                           keep_water=None,
                           continuous_width: bool = False):
    """:func:`filter_shoreline` at the LOCAL element size, not one size.

    One size for the whole hole let 60-90 m features and walls survive where
    the transition's elements are 150-400 m, and they came back as the worst
    elements of the Kimitsu port mesh (4.7 deg on a coastline 2.3 km from the
    region).  The owner's rule (2026-09-23): the size decides, wherever it is
    -- which inside the region is h0 itself, so the region is unchanged.

    The size field is sampled over ``footprint`` and cut into octave bands,
    ``[h0 * 2**k, h0 * 2**(k+1))``; each band takes the land filtered at its
    LOWER bound, so no band removes a feature the finest element in it could
    carry.  The pieces are joined and filtered once more at ``h0``, which
    removes the slivers a seam between two bands leaves and nothing else --
    every band already removed more than that.  (Quarter-octave bands were
    tried for Funabashi and moved the Kimitsu transition coast by 211 m.)

    ``land_width_factor`` applies in bands up to ``land_width_max_band``
    (all bands when None); coarser bands keep one threshold for land and
    water.  Out in the coarse
    transition, land kept down to half a 240 m element put a 200 m spike
    of OSM land against the frozen interface's 440-600 m edges, and an
    element of 17.8 deg between them.

    Returns ``(filtered, report)`` with a line per band.
    """
    import shapely
    from rasterio import features
    from rasterio.transform import from_origin

    # an iterator of polygons would be spent by the first band (review, round 3)
    if not hasattr(land, "geom_type"):
        land = list(land)
    if not (np.isfinite(h0) and h0 > 0):
        raise ValueError("h0 must be finite and positive")
    if footprint is None or footprint.is_empty:
        raise ValueError("the footprint is empty: there is no local size to filter at")
    s = float(spacing if spacing is not None else h0)
    if not (np.isfinite(s) and s > 0):
        raise ValueError("spacing must be finite and positive")
    x0, y0, x1, y1 = footprint.bounds
    nx, ny = int(np.ceil((x1 - x0) / s)) + 1, int(np.ceil((y1 - y0) / s)) + 1
    gx = x0 + (np.arange(nx) + 0.5) * s
    gy = y1 - (np.arange(ny) + 0.5) * s
    mx, my = np.meshgrid(gx, gy)
    h = np.asarray(size_field(np.column_stack([mx.ravel(), my.ravel()])), dtype=float)
    # checked before the h0 clip, which made a zero or negative size look
    # like the finest band (review, round 5)
    if h.shape not in ((nx * ny,), (nx * ny, 1)) or not np.isfinite(h).all() or (h <= 0).any():
        raise ValueError("the size field is not finite and positive over the footprint")
    h = h.reshape(ny, nx)
    band = np.floor(np.log2(np.maximum(h, h0) / h0)).astype(np.int32)
    transform = from_origin(x0, y1, s, s)
    zones: dict[int, list] = {}
    for geom, value in features.shapes(band, transform=transform):
        zones.setdefault(int(value), []).append(shapely.geometry.shape(geom))
    pieces, rows = [], []
    land_union = shapely.union_all([land] if hasattr(land, "geom_type") else list(land)) \
        if keep_land is not None or keep_water is not None else None

    def mostly(q, keep):
        return q.geom_type == "Polygon" and q.area > 0 \
            and shapely.intersection(q, keep).area >= 0.5 * q.area
    for k in sorted(zones):
        hk = h0 * 2.0 ** k
        zone = shapely.intersection(shapely.union_all(zones[k]), footprint)
        if zone.is_empty:
            continue
        lw = land_width_factor if (land_width_max_band is None
                                   or k <= land_width_max_band) else None
        fk, rk = filter_shoreline(land, hk, elements_per_feature=elements_per_feature,
                                  land_width_factor=lw)
        if keep_land is not None and lw is None:
            # In the coarse bands land under two elements goes like water
            # does, and where both are unresolvable the order decides: at
            # Yokohama two 150 m quay blocks the base mesh has as land went
            # first and left 0.19 km2 of new water behind.  A piece of land
            # the band removes stays land if most of it is ``keep_land``
            # (land in the source AND in the base) -- the WHOLE piece: kept
            # only where the base has it, a Funabashi peninsula was cut in
            # two along the base's coast, and the half that went left a wall
            # and a 27 deg pocket against the half that stayed.  The narrow
            # water round a kept piece is closed by the h0 pass below.
            removed = shapely.difference(land_union, fk)
            back = [q for q in getattr(removed, "geoms", [removed]) if mostly(q, keep_land)]
            if back:
                fk = shapely.union_all([fk, *back])
        if keep_water is not None and lw is None:
            # The same for water: a piece of water the band fills stays water
            # if most of it is ``keep_water`` (water in the base too).  The
            # transition serves efficiency, not detail (owner, 2026-09-28):
            # on this project's base, which draws the Odaiba channels, a
            # 240 m band closed one inside the hole while the frozen mesh
            # beyond kept it, and the seam got a spit the fill could not
            # carry (8 violations).
            filled = shapely.difference(fk, land_union)
            back = [q for q in getattr(filled, "geoms", [filled]) if mostly(q, keep_water)]
            if back:
                fk = shapely.difference(fk, shapely.union_all(back))
        part = shapely.intersection(fk, zone)
        pieces.append(part)
        rows.append({"band": k, "h_m": hk, "zone_km2": float(zone.area / 1e6),
                     "removes_narrower_than_m": rk["removes_features_narrower_than_m"],
                     "removes_land_narrower_than_m": rk["removes_land_narrower_than_m"]})
    # land outside the sampled footprint is kept exactly as the finest band has it
    f0, _ = filter_shoreline(land, h0, elements_per_feature=elements_per_feature,
                             land_width_factor=land_width_factor)
    pieces.append(shapely.difference(f0, footprint))
    joined = shapely.union_all([q for q in pieces if not q.is_empty])
    cw_rep = None
    if continuous_width:
        # Water judged at the local size itself, in the coarse zone only
        # (h > 2 h0): no octave lower bound and no band cut there.  The bands
        # close water under 2 hk, which is between h and 2 h of the elements
        # actually there -- 1.5 h on average; judged at two whole elements
        # it closed 1.26 km2 of the Odaiba port and 0.53 km2 of Yokohama's.
        # So the continuous rule keeps the bands' average: water narrower than
        # 0.75 of ``elements_per_feature`` elements goes.
        cw, cw_rep = unresolvable_water(joined, size_field, footprint,
                                        radius_factor=0.375 * float(elements_per_feature),
                                        min_h=2.0 * h0, spacing=h0 / 3.0)
        if keep_water is not None and not cw.is_empty:
            kept = [q for q in getattr(cw, "geoms", [cw]) if not mostly(q, keep_water)]
            cw_rep["n_kept_as_base_water"] = len(getattr(cw, "geoms", [cw])) - len(kept)
            cw = shapely.union_all(kept) if kept else shapely.Polygon()
        if not cw.is_empty:
            joined = shapely.union_all([joined, cw])
    out, rep = filter_shoreline(joined, h0, elements_per_feature=elements_per_feature,
                                land_width_factor=land_width_factor)
    before = shapely.union_all([land] if hasattr(land, "geom_type") else list(land))
    rep = {**rep, "bands": rows,
           "continuous_width": cw_rep,
           "land_lost_m2": float(shapely.difference(before, out).area),
           "water_lost_m2": float(shapely.difference(out, before).area),
           "area_before_m2": float(before.area), "area_after_m2": float(out.area)}
    return out, rep


def _resolve_stretch(pts: np.ndarray, shoreline, size, *, fine_h=None) -> np.ndarray:
    """Follow the source, over the whole stretch.

    An earlier version kept the base polyline wherever the local element size
    was coarser than the base's own spacing, because resolving a coastline in
    a 400-1700 m field replaced 17 base nodes with 13.  The owner settled it
    (2026-09-23): the transition is treated like the region, the source is
    OSM throughout, and the base is not consulted for the SHAPE anywhere
    inside the hole.  A coarse part of the stretch is therefore a coarse
    sampling of OSM, not a copy of the base -- fewer nodes, but one
    provenance.

    What decides whether a feature survives is :func:`filter_shoreline`,
    applied to the source before any of this, at the declared grid size.
    """
    return _resample_on_source(np.asarray(pts, dtype=float)[:, :2], shoreline,
                               size, pointwise=True, require=True, fine_h=fine_h)


def _resample_on_source(pts: np.ndarray, shoreline, size,
                        simplify_frac: float = 0.25, *,
                        pointwise: bool = False,
                        require: bool = False,
                        fine_h: float | None = None) -> np.ndarray:
    """Re-space one stretch along the source shoreline between its endpoints.

    The source is simplified to a fraction of the LOCAL element size first.
    A shoreline digitised at metres cannot be represented by 400 m elements,
    and walking it at 400 m produces chords that turn sharply against each
    other: measured on the Futtsu patch, every surviving QA failure was a
    coastal element 1.8-2.3 km out, where the transition is coarse and the
    coast is not.  Simplifying is not throwing detail away -- there is no
    room for it at that size.

    ``pointwise`` is what the hires branch needs.  The default tolerance is a
    quarter of the MEDIAN size over the whole substring, and a stretch that is
    fine in the core and coarse in the transition therefore has its core
    detail simplified away at the transition's scale (review P2-11).  The
    pointwise form splits the substring into octave bands of the local size
    and simplifies each with a quarter of that band's SMALLEST size, so detail
    survives wherever there are elements small enough to carry it.
    """
    import shapely

    coords = _source_substring(pts, shoreline)
    if coords is None:
        if require:
            raise ValueError(
                "coastline: resolve found no source component for a stretch")
        return _subdivide(pts, size)
    if not pointwise:
        h = float(np.median(_size_at(size, coords)))
        simple = shapely.simplify(shapely.LineString(coords), simplify_frac * h)
        out = np.asarray(simple.coords, dtype=float)[:, :2]
        return _walk(out if len(out) >= 2 else coords, size)

    h = _size_at(size, coords)
    # One band per octave of the local size.  Runs are simplified separately
    # and rejoined; the shared vertex at each join is kept exactly once, so
    # the curve stays continuous and no band can move another band's ends.
    band = np.floor(np.log2(h / h.min())).astype(np.int64)
    cuts = np.flatnonzero(np.diff(band)) + 1
    out = [coords[0]]
    for a, b in zip(np.concatenate([[0], cuts]),
                    np.concatenate([cuts, [len(coords)]])):
        run = coords[a:min(b + 1, len(coords))]
        if len(run) < 2:
            continue
        tol = simplify_frac * float(h[a:b].min())
        simple = np.asarray(
            shapely.simplify(shapely.LineString(run), tol).coords,
            dtype=float)[:, :2]
        out.extend((simple if len(simple) >= 2 else run)[1:])
    out = np.asarray(out, dtype=float)
    # corners kept: an arc-length walk cut across every pier narrower than
    # two elements, and piers down to half an element are land now
    return _corner_walk(out if len(out) >= 2 else coords, size, fine_h=fine_h)


def rim_constraints(
    nodes,
    selection: PatchSelection,
    *,
    size,
    coastline: str = "resample",
    boundary_edges=None,
    shoreline=None,
    tolerance_m: float = 100.0,
    fine_h: float | None = None,
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

    A ring with no frozen node -- an island inside the hole -- has nothing to
    anchor a stretch.  It is copied as it stands, except under ``resolve``
    with a ``shoreline``: then it is left out (``n_islands_left_to_source``)
    and the caller adds the source's land inside the hole with
    :func:`island_rings`.
    """
    xy = np.asarray(nodes, dtype=float)[:, :2]
    frozen = set(selection.frozen_nodes.tolist())

    pts: list[np.ndarray] = []
    base_id: list[int] = []
    segs: list[tuple[int, int]] = []
    curves: list[np.ndarray] = []
    curve_of: dict[int, int] = {}
    n_resampled = 0
    n_uncrossed = 0
    placed: list = []
    kept_because: dict[str, int] = {}
    kept_at: list = []
    n_new = 0
    n_islands_left = 0
    # A replaced stretch may cross neither the base's own boundary nor the
    # interface between the hole and the retained mesh: at Funabashi on this
    # project's base a resolved stretch cut across an interface edge, and the
    # rim came out crossing itself (the boundary check never saw it).
    iface = np.asarray(selection.rim_edges, dtype=np.int64).reshape(-1, 2)[
        ~np.asarray(selection.physical_rim, dtype=bool)]
    check_edges = iface if boundary_edges is None else np.vstack(
        [np.asarray(boundary_edges, dtype=np.int64).reshape(-1, 2), iface])
    # ...but not the coast of a base island the source replaces: it is not
    # on the rim at all (below), and a pier meeting it was kept on the base
    # line with its whole harbour (Kimitsu, on this project's base)
    if coastline == "resolve" and shoreline is not None:
        left = {int(v) for r in selection.rings
                if all(int(v) not in frozen for v in r) for v in r}
        if left:
            check_edges = check_edges[~np.isin(check_edges, list(left)).all(axis=1)]

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

        if is_free.all() and coastline == "resolve" and shoreline is not None:
            # Under `resolve` the source decides, and a copied base island
            # kept a 200 m square where the source had a 100 m one (review,
            # round 7).  Left out here: the caller adds the source's own land
            # inside the hole as islands at the local size (island_rings,
            # which the driver runs next), and the land check covers what it
            # refuses.  A base island the source does not have is water.
            n_islands_left += 1
            continue
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
                    new = coastline_points(xy[idx], size, mode=coastline, fine_h=fine_h,
                                           shoreline=shoreline,
                                           tolerance_m=tolerance_m)
                    why = _unusable_replacement(new, xy, idx, check_edges,
                                                placed, where=kept_at)
                    if why is not None:
                        new = _subdivide(xy[idx], size)
                        n_uncrossed += 1
                        kept_because[why] = kept_because.get(why, 0) + 1
                        curves.append(np.asarray(xy[idx], dtype=float))
                    elif coastline == "resolve":
                        # The DELIVERED polyline, not the source substring:
                        # `resolve` follows the source over the whole stretch
                        # (_resolve_stretch), simplified and cut at the local
                        # size, and the curve is what the seam repair may
                        # SLIDE these nodes along.  Sliding along the raw
                        # substring would pull a node off the delivered
                        # coastline into detail the local size removed.
                        curves.append(np.asarray(new, dtype=float))
                    else:
                        curves.append(coastline_curve(xy[idx], coastline, shoreline))
                    if len(new) > 1:
                        import shapely as _sh

                        placed.append(_sh.LineString(
                            np.asarray(new, dtype=float)[:, :2]))
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
        "n_stretches_kept_to_avoid_a_crossing": n_uncrossed,
        "kept_because": kept_because,
        "kept_at": [[round(float(x), 1), round(float(y), 1)] for x, y in kept_at],
        "n_coastline_nodes_replaced": n_resampled,
        "n_coastline_nodes_new": n_new,
        "n_islands_left_to_source": n_islands_left,
        "coastline_mode": coastline,
    }


def _push(pts: list, base_id: list, q, bid: int) -> int:
    pts.append(np.asarray(q, dtype=float))
    base_id.append(int(bid))
    return len(pts) - 1


def _corner_walk(pts: np.ndarray, size, *, closed: bool = False,
                 fine_h: float | None = None) -> np.ndarray:
    """Resample a polyline at the local size, keeping its corners.

    ``_walk`` places stations by arc length and so cuts every corner it
    passes, which on a pier narrower than two elements cuts across the pier:
    its stations go out along one side and come back along the other, and
    the chord between them runs over the tip.  Here every vertex of the
    (already simplified) line is a candidate; an edge under half an element
    loses the end that turns less, never a stretch's own two ends; then

    * a STEP of under 0.75 of an element -- an edge whose ends turn by more
      than 60 deg in OPPOSITE directions -- also becomes one point;
    * a TIP -- an edge under 0.75 of an element whose two ends both turn by
      more than 60 deg the same way, the end of a pier 15-22 m wide at a
      30 m size -- becomes one point at its midpoint (owner, 2026-09-24):
      an end that short puts an element under 30 deg beside it;
    * a SPIKE -- a vertex whose two edges are both under 0.75 of an element
      and which turns by more than 60 deg -- goes;
    * long edges are subdivided, every vertex kept.

    ``closed`` treats ``pts`` as a ring (not repeated at the end) and
    returns it the same way. ``fine_h`` limits the step (0.5-0.75 h) and
    spike rules to where the local size is at most ``fine_h``.
    """
    pts = np.asarray(pts, dtype=float)[:, :2].copy()
    # EXACT repetition only: allclose's relative tolerance at a UTM northing
    # of 3.9e6 counts two corners 30 m apart as one, and a 100 x 30 m quay
    # block lost half its area to it (review F2)
    if closed and len(pts) > 1 and np.array_equal(pts[0], pts[-1]):
        pts = pts[:-1]

    def fixed(k):
        return not closed and k in (0, len(pts) - 1)

    def turn(k):
        """Signed turn at vertex k, in degrees (nan at an open end)."""
        if fixed(k):
            return np.nan
        a = pts[k] - pts[k - 1]
        b = pts[(k + 1) % len(pts)] - pts[k]
        cr = a[0] * b[1] - a[1] * b[0]
        return float(np.degrees(np.arctan2(cr, a @ b)))

    def edges():
        n = len(pts)
        return [(k, (k + 1) % n) for k in range(n if closed else n - 1)]

    def fine(h):
        # The step and spike rules serve element quality at the TARGET size.
        # Applied in the coarse transition they collapsed 100 m jogs among
        # 200 m elements that already met QA, and moved the Futtsu coast up
        # to 47 m off OSM for nothing. ``fine_h`` bounds where they act.
        return fine_h is None or h <= fine_h

    made: set = set()
    changed = True
    while changed and len(pts) > (3 if closed else 2):
        changed = False
        # short edges: drop the end that turns less
        cand = []
        for i, j in edges():
            L = float(np.linalg.norm(pts[j] - pts[i]))
            h = float(_size_at(size, 0.5 * (pts[i] + pts[j])[None])[0])
            if L < 0.5 * h and not (fixed(i) and fixed(j)):
                cand.append((L / h, i, j))
        if cand:
            _, i, j = min(cand)
            if fixed(i):
                drop = j
            elif fixed(j):
                drop = i
            elif min(abs(turn(i)), abs(turn(j))) < 20.0:
                # one end is on a straight run: it goes, and the shape with it
                drop = i if abs(turn(i)) <= abs(turn(j)) else j
            else:
                # a STEP, both ends corners: dropping either cut a 14 m step
                # in a quay into a 165 m chord across the water beside it.
                # The two become one point between them instead.
                pts[i] = 0.5 * (pts[i] + pts[j])
                drop = j
            pts = np.delete(pts, drop, axis=0)
            changed = True
            continue
        # tips: a short end between two turns the same way
        for i, j in edges():
            if fixed(i) or fixed(j):
                continue
            L = float(np.linalg.norm(pts[j] - pts[i]))
            h = float(_size_at(size, 0.5 * (pts[i] + pts[j])[None])[0])
            ti, tj = turn(i), turn(j)
            if L < 0.75 * h and abs(ti) > 60.0 and abs(tj) > 60.0 \
                    and np.sign(ti) != np.sign(tj) and fine(h):
                # a STEP of up to 0.75 of an element: a 19.6 m jog in an
                # Odaiba quay left an element of 29.7 deg at 30 m. The two
                # corners become one point between them.
                pts[i] = 0.5 * (pts[i] + pts[j])
                pts = np.delete(pts, j, axis=0)
                changed = True
                break
            if L < 0.75 * h and abs(ti) > 60.0 and abs(tj) > 60.0 \
                    and np.sign(ti) == np.sign(tj):
                # The sides stay parallel up to half an element short of the
                # end and close from there onto the end's midpoint.  Joining
                # the root corners straight to the midpoint turned a 20 x 70 m
                # pier into a triangle.
                n = len(pts)
                p_, q_ = (i - 1) % n, (j + 1) % n
                back = []
                for c_, o_ in ((i, p_), (j, q_)):
                    side = pts[o_] - pts[c_]
                    ls = float(np.linalg.norm(side))
                    d = min(0.5 * h, 0.5 * ls)
                    back.append(pts[c_] + side / max(ls, 1e-12) * d)
                mid = 0.5 * (pts[i] + pts[j])
                new = [back[0], mid, back[1]]
                made.update(tuple(q) for q in new)   # a gable is not a spike
                if j == 0:                      # the end wraps round a ring
                    pts = np.vstack([pts[1:i], new])
                else:
                    pts = np.vstack([pts[:i], new, pts[j + 1:]])
                changed = True
                break
        if changed:
            continue
        # spikes: one vertex with two short edges, turning sharply -- a small
        # triangular bump of 16 and 20 m on the Odaiba coast made an element
        # of 29.7 deg at 30 m. The vertex goes, the chord replaces it.
        n = len(pts)
        for k in range(n):
            if fixed(k) or (not closed and k in (0, n - 1)) or tuple(pts[k]) in made:
                continue
            a_ = float(np.linalg.norm(pts[k] - pts[k - 1]))
            b_ = float(np.linalg.norm(pts[(k + 1) % n] - pts[k]))
            h = float(_size_at(size, pts[k][None])[0])
            if a_ < 0.75 * h and b_ < 0.75 * h and abs(turn(k)) > 60.0 and fine(h):
                pts = np.delete(pts, k, axis=0)
                changed = True
                break
    if closed:
        return _subdivide(np.vstack([pts, pts[:1]]), size)[:-1]
    return _subdivide(pts, size)


def _ring_at_size(ring: np.ndarray, size, *, fine_h=None) -> np.ndarray:
    """A closed ring resampled at the local size, its corners kept.

    Douglas-Peucker at a tenth of the local element drops the wiggles, then
    :func:`_corner_walk` merges short edges, points narrow tips and
    subdivides.  Returns the points once each, not closed.
    """
    import shapely

    ring = np.asarray(ring, dtype=float)[:, :2]
    h_min = float(np.min(_size_at(size, ring)))
    g = shapely.simplify(shapely.LinearRing(ring), 0.1 * h_min)
    return _corner_walk(np.asarray(g.coords, dtype=float)[:-1], size, closed=True,
                        fine_h=fine_h)


def rim_repair(pfix, egfix, pfix_base, water, size, *, protect=(),
               min_edge_factor=0.5, gap_factor=1.0, min_angle_deg=60.0,
               operations=("short_edges", "slits", "angles"), rounds=2,
               focus=None, focus_factor=1.5, retreat_tips=True, size_floor=None):
    """Check the finished rim against the local size and repair what fails.

    Each rule upstream (the filter, the corner walk, blunting, rooting) is
    right on the geometry it was given, but a later one can leave what an
    earlier one would have refused.  This pass looks at the rim as the fill
    will get it, and repairs it by a fixed set of operations in a fixed
    order, ``rounds`` times (review 6, 2026-09-25):

    ``short_edges``  an edge shorter than ``min_edge_factor`` of the local
        element loses one of its ends -- a free point (not frozen, two
        constrained edges, not in ``protect``) whose neighbours are joined
        instead -- if the new edge crosses no other rim edge, the triangle
        handed between land and water is no wider than the short edge, and
        no water angle is sharpened below ``min_angle_deg``.  (Yokohama:
        blunting stopped at 0.9 of a quay edge and left 3.0 m beside its
        cap among 30 m elements; 13.1 and 19.5 deg.)
    ``slits``  a point closer than ``gap_factor`` of an element to a rim
        edge it is not next to, across water, is a throat (one element by
        default: a slot a blunting cap leaves 14-33 m wide among 30 m
        elements was missed at half an element, Yokohama): the throat is
        cut, and the dead end beyond it becomes land if no disc of one
        element fits in it and nothing on its outline is frozen or a wall
        root.  No angle rule sees this -- the pier corner at a 3.3 m slit
        is a 311 deg water angle (Yokohama, review 6).
    ``angles``  :func:`blunt_acute_corners` once more on the finished rim,
        wall roots protected: rooting and the operations above can make a
        water corner under ``min_angle_deg``.

    ``retreat_tips`` False leaves a tip where it is (the slit is reported).

    ``size_floor`` is the smallest size the field has anywhere; the slit and
    retreat decisions use a guaranteed lower bound of the size over the
    water or land they change hands (:func:`size_lower_bound`), not the size
    at a few vertices (review, round 4).  Default: the field's own
    ``size_min``; for a field without one, half the smallest size at the rim
    points sets the grid step only (:func:`resolve_size_floor`).

    ``focus`` ((m, 2) points, e.g. the QA offenders of a failed seed)
    limits the short-edge and slit operations to within ``focus_factor``
    local elements of one of them, so a retry with looser thresholds touches
    only where the mesh failed.

    Nothing is forced: what no operation can repair is reported.  ``water``
    is the hole the rim bounds; ``size`` maps (n, 2) points to the local
    element size.  Returns ``(pfix, egfix, pfix_base, remap, report)``;
    ``remap[old] = new`` (-1 for a removed point) for the input's points, so
    that references into the rim (wall roots) can follow it.
    """
    import shapely
    from shapely.ops import nearest_points

    pfix = np.asarray(pfix, dtype=float)[:, :2].copy()
    egfix = np.asarray(egfix, dtype=np.int64).copy()
    pfix_base = np.asarray(pfix_base, dtype=np.int64).copy()
    n_in = len(pfix)
    floor, floor_ok = resolve_size_floor(size, size_floor, pfix) if len(pfix) else (1.0, False)
    ident = np.arange(n_in)             # current index of each input point
    protect = set(int(k) for k in protect)
    removed, refused, slits, slits_left, retreated = [], [], [], [], []
    angles_rep: list = []

    def key(xy):
        return (round(float(xy[0]), 6), round(float(xy[1]), 6))

    protected_xy = {key(pfix[k]) for k in protect}
    focus_pts = None if focus is None else np.asarray(focus, dtype=float).reshape(-1, 2)

    def in_focus(xy):
        if focus_pts is None:
            return True
        if not len(focus_pts):
            return False
        h = float(np.asarray(size(np.asarray([xy])), dtype=float)[0])
        return bool(np.min(np.linalg.norm(focus_pts - xy, axis=1)) <= focus_factor * h)

    def is_protected(k):
        return key(pfix[k]) in protected_xy

    def nbrs_of():
        nb: dict = {}
        for a, b in egfix.tolist():
            nb.setdefault(a, []).append(b)
            nb.setdefault(b, []).append(a)
        return nb

    def free(k, nb):
        return pfix_base[k] < 0 and len(nb.get(k, [])) == 2 and not is_protected(k)

    def compact(drop):
        """Delete points ``drop``; keep ``ident`` in step."""
        nonlocal pfix, egfix, pfix_base, ident
        keep = np.ones(len(pfix), dtype=bool)
        keep[list(drop)] = False
        m = np.full(len(pfix), -1, dtype=np.int64)
        m[keep] = np.arange(int(keep.sum()))
        egfix = m[egfix]
        pfix, pfix_base = pfix[keep], pfix_base[keep]
        ident = np.where(ident >= 0, m[np.maximum(ident, 0)], -1)

    def angle_at(xy, p_xy, n_xy, poly):
        up, un = p_xy - xy, n_xy - xy
        lp, ln = float(np.linalg.norm(up)), float(np.linalg.norm(un))
        if lp <= 0 or ln <= 0:
            return None
        ang = float(np.degrees(np.arccos(np.clip(up @ un / (lp * ln), -1, 1))))
        bis = up / lp + un / ln
        if np.linalg.norm(bis) < 1e-9:
            return 180.0
        probe = xy + bis / np.linalg.norm(bis) * 0.05 * min(lp, ln)
        return ang if poly.contains(shapely.Point(probe)) else 360.0 - ang

    # ------------------------------------------------------------ short edges
    def short_edges():
        nonlocal egfix
        # refused in THIS pass only: an edge refused before may be removable
        # once blunting or a slit changed its neighbours (review, round 12)
        done_refused: set = set()
        for _ in range(len(pfix)):
            nb = nbrs_of()
            L = np.linalg.norm(pfix[egfix[:, 0]] - pfix[egfix[:, 1]], axis=1)
            mid = 0.5 * (pfix[egfix[:, 0]] + pfix[egfix[:, 1]])
            h = np.asarray(size(mid), dtype=float)
            order = [k for k in np.argsort(L / h) if L[k] < min_edge_factor * h[k]
                     and in_focus(mid[k])]
            order = [k for k in order
                     if tuple(sorted((key(pfix[egfix[k, 0]]), key(pfix[egfix[k, 1]]))))
                     not in done_refused]
            if not order:
                return
            k = order[0]
            a, b = int(egfix[k, 0]), int(egfix[k, 1])
            options = []
            for x in (a, b):
                if not free(x, nb):
                    continue
                u, w = nb[x]
                # u and w already joined: the ring is a triangle, and taking
                # x would leave two points and a doubled edge (review, round 1)
                if u == w or w in nb.get(u, []):
                    continue
                du, dw = pfix[u] - pfix[x], pfix[w] - pfix[x]
                tri = abs(float(du[0] * dw[1] - du[1] * dw[0])) / 2.0
                if tri > 0.5 * float(L[k]) * float(h[k]):
                    continue
                new = shapely.LineString([pfix[u], pfix[w]])
                others = [shapely.LineString(pfix[[c, d]]) for c, d in egfix.tolist()
                          if not ({c, d} & {x, u, w})]
                if others and shapely.MultiLineString(others).intersects(new):
                    continue
                # the ring it leaves must still be a polygon: three points
                # in a line passed every test above and returned a ring of
                # zero area (review, round 7)
                eg_try = np.vstack([egfix[~np.isin(egfix, [x]).any(axis=1)], [[u, w]]])
                if not _rim_edit_ok(pfix, egfix, pfix, eg_try):
                    continue
                ok = True
                for y, far in ((u, w), (w, u)):
                    other = [q for q in nb.get(y, []) if q != x]
                    if len(other) != 1:
                        continue
                    before = angle_at(pfix[y], pfix[other[0]], pfix[x], water)
                    after = angle_at(pfix[y], pfix[other[0]], pfix[far], water)
                    if after is None or (after < min_angle_deg and
                                         (before is None or after < before - 1e-6)):
                        ok = False
                if ok:
                    options.append((tri, x, u, w))
            if not options:
                edge_xy = (key(pfix[a]), key(pfix[b]))
                done_refused.add(tuple(sorted(edge_xy)))
                if all(tuple(sorted(r["edge_xy"])) != tuple(sorted(edge_xy)) for r in refused):
                    refused.append({"edge_xy": edge_xy, "length_m": round(float(L[k]), 2),
                                    "h_m": round(float(h[k]), 1),
                                    "at": [round(float(v), 1) for v in mid[k]]})
                continue
            tri, x, u, w = min(options)
            removed.append({"at": [round(float(v), 1) for v in pfix[x]],
                            "edge_m": round(float(L[k]), 2), "h_m": round(float(h[k]), 1),
                            "area_m2": round(tri, 1)})
            egfix = np.vstack([egfix[~np.isin(egfix, [x]).any(axis=1)], [[u, w]]])
            compact([x])

    def retreat(v, q, h, nb, poly):
        """Move point v straight away from q until the gap is no throat.

        To just past the throat test (``gap_factor`` elements): placed at one
        element with the test at one, rounding left it a hair inside and the
        same tip was stepped back 2,180 times by nothing; with the retry's
        looser test it could never get out (Yokohama, after review round 2).
        A step that gets nowhere is refused.
        """
        nonlocal pfix
        if not free(v, nb):
            return False
        d = pfix[v] - q
        dist = float(np.linalg.norm(d))
        target = 1.001 * max(gap_factor, 1.0) * h
        if dist <= 0 or target - dist < 1e-3:
            return False
        new = q + d / dist * target
        h = target
        u, w = nb[v]
        # it must go INTO the land (the pier), not into the water
        if poly.contains(shapely.Point(new)):
            return False
        # ...and hand over only land no element fits in: stepping back by
        # 1.5 elements took 5,000 m2 off a 120 m pier (review, round 3)
        handed = shapely.make_valid(shapely.Polygon([pfix[u], pfix[v], pfix[w], new]))
        if land_an_element_fits(handed, size, floor, certified=floor_ok):
            return False
        edges = shapely.MultiLineString([[pfix[u], new], [new, pfix[w]]])
        others = [shapely.LineString(pfix[[c, e2]]) for c, e2 in egfix.tolist()
                  if not ({c, e2} & {v, u, w})]
        if others and shapely.MultiLineString(others).intersects(edges):
            return False
        for y, far in ((u, w), (w, u)):
            other = [x for x in nb.get(y, []) if x != v]
            if len(other) != 1:
                continue
            a0 = angle_at(pfix[y], pfix[other[0]], pfix[v], poly)
            a1 = angle_at(pfix[y], pfix[other[0]], new, poly)
            if a1 is None or (a1 < min_angle_deg and (a0 is None or a1 < a0 - 1e-6)):
                return False
        a0 = angle_at(pfix[v], pfix[u], pfix[w], poly)
        a1 = angle_at(new, pfix[u], pfix[w], poly)
        if a1 is None or (a1 < min_angle_deg and (a0 is None or a1 < a0 - 1e-6)):
            return False
        # ...and every ring keeps its role: a tip stepped back over a small
        # lake turned its water into an island (review, round 10)
        p_try = pfix.copy()
        p_try[v] = new
        if not _rim_edit_ok(pfix, egfix, p_try, egfix):
            return False
        retreated.append({"at": [round(float(c), 1) for c in pfix[v]],
                          "by_m": round(h - dist, 1), "gap_was_m": round(dist, 2)})
        pfix = p_try
        return True

    # ------------------------------------------------------------------ slits
    def slit_once(poly):
        """Close ONE slit; True if something changed."""
        nonlocal pfix, egfix, pfix_base, ident
        nb = nbrs_of()
        segs = shapely.linestrings(np.stack([pfix[egfix[:, 0]], pfix[egfix[:, 1]]], axis=1))
        tree = shapely.STRtree(segs)
        h_all = np.asarray(size(pfix), dtype=float)
        cands = []
        for v in range(len(pfix)):
            if len(nb.get(v, [])) != 2 or not in_focus(pfix[v]):
                continue
            near = {v, *nb[v]}
            for u in list(nb[v]):
                near.update(nb.get(u, []))
            pt = shapely.Point(pfix[v])
            for k in tree.query(pt.buffer(gap_factor * h_all[v])):
                a, b = int(egfix[k, 0]), int(egfix[k, 1])
                if a in near or b in near:
                    continue
                d = float(segs[k].distance(pt))
                if d < gap_factor * h_all[v]:
                    cands.append((d / h_all[v], v, k))
        for _, v, k in sorted(cands):
            a, b = int(egfix[k, 0]), int(egfix[k, 1])
            q = np.asarray(nearest_points(segs[k], shapely.Point(pfix[v]))[0].coords[0])
            throat = shapely.LineString([pfix[v], q])
            if throat.length <= 0 or not poly.buffer(1e-6).contains(throat):
                continue
            # the two ways round from v to the edge (a, b): each, closed by
            # the throat, bounds one side
            sides = []
            for first in nb[v]:
                path, prev, cur = [v], v, first
                while cur not in (a, b) and len(path) <= len(pfix):
                    path.append(cur)
                    nxt = [x for x in nb.get(cur, []) if x != prev]
                    if len(nxt) != 1:
                        path = None
                        break
                    prev, cur = cur, nxt[0]
                if path is None or cur not in (a, b):
                    continue
                path.append(cur)
                ring = np.vstack([pfix[path], q[None]])
                g = shapely.make_valid(shapely.Polygon(ring))
                sides.append((float(g.area), path, g))
            if len(sides) != 2:
                continue                    # not one loop: an island, or a junction
            area, path, g = min(sides, key=lambda t: t[0])
            inner = path[1:-1]              # the points that go
            end = path[-1]                  # the edge end on the dead-end side
            if not inner or shapely.difference(g, poly.buffer(1e-6)).area > 1e-6 * max(area, 1.0):
                continue                    # the side is not water in the hole
            # no element may fit ANYWHERE in the dead end, judged on a
            # guaranteed lower bound of the size over it, not at its
            # vertices: a finer basin behind a coarse mouth was closed
            # (review, round 4)
            h = size_lower_bound(g, size, floor, certified=floor_ok) if not g.is_empty \
                else float(np.min(h_all[path]))
            if not g.buffer(-h).is_empty:
                # Closing would cut off water an element fits in -- a pier
                # whose tip nearly touches the quay across (Yokohama, 3.3 m).
                # The tip steps back instead, straight away from the quay,
                # until the gap is one element.
                if retreat_tips and retreat(v, q, float(h_all[v]), nb, poly):
                    return True
                slits_left.append({"at": [round(float(c), 1) for c in pfix[v]],
                                   "why": "wide enough for an element, and the tip "
                                          "could not step back"})
                continue
            if any(pfix_base[x] >= 0 or is_protected(x) for x in inner + [end]):
                slits_left.append({"at": [round(float(c), 1) for c in pfix[v]],
                                   "why": "frozen or wall root inside"})
                continue
            other_end = b if end == a else a
            drop = set(inner) | {end}
            keep_e = ~np.isin(egfix, list(drop)).any(axis=1)
            # q snaps to the far end when it is that close already -- but the
            # snapped chord is not the throat that was checked, and one
            # crossed the edge beside a protected corner (review, round 8):
            # each candidate must leave the ring a simple polygon
            cands = [(len(pfix), [q])]
            if float(np.linalg.norm(q - pfix[other_end])) <= 0.25 * h:
                cands.insert(0, (other_end, []))
            chosen = None
            for qi, new_pts in cands:
                new_e = [[v, qi]] + ([[qi, other_end]] if qi != other_end else [])
                eg_try = np.vstack([egfix[keep_e], np.asarray(new_e, dtype=np.int64)])
                p_try = np.vstack([pfix, np.asarray(new_pts)]) if new_pts else pfix
                if not _rim_edit_ok(pfix, egfix, p_try, eg_try):
                    continue
                try:                            # nor cross another ring
                    hole_polygon(p_try, eg_try)
                except ValueError:
                    continue
                chosen = (eg_try, p_try, bool(new_pts))
                break
            if chosen is None:
                slits_left.append({"at": [round(float(c), 1) for c in pfix[v]],
                                   "why": "closing it would not leave a simple ring"})
                continue
            egfix, p_try, added = chosen
            if added:
                pfix = p_try
                pfix_base = np.concatenate([pfix_base, [-1]])
            slits.append({"at": [round(float(c), 1) for c in pfix[v]],
                          "throat_m": round(float(throat.length), 2),
                          "area_m2": round(area, 1), "h_m": round(h, 1)})
            compact(drop)
            return True
        return False

    # ----------------------------------------------------------------- angles
    def angles(poly):
        nonlocal pfix, egfix, pfix_base, ident
        prot = [k for k in range(len(pfix)) if is_protected(k)]
        before = {key(pfix[k]): k for k in range(len(pfix))}
        p2, e2, b2, rep = blunt_acute_corners(pfix, egfix, pfix_base, poly, size,
                                              min_angle_deg, protect=prot)
        if rep["n_corners_blunted"] == 0:
            return
        angles_rep.append(rep)
        m = np.full(len(pfix), -1, dtype=np.int64)
        for j, xy in enumerate(p2):
            k = before.get(key(xy))
            if k is not None:
                m[k] = j
        ident = np.where(ident >= 0, m[np.maximum(ident, 0)], -1)
        pfix, egfix, pfix_base = np.asarray(p2, dtype=float)[:, :2], np.asarray(e2), \
            np.asarray(b2)

    def n_done():
        return (len(removed), len(slits), len(retreated),
                sum(r["n_corners_blunted"] for r in angles_rep))

    for _ in range(max(1, int(rounds))):
        n_before = n_done()
        if "short_edges" in operations:
            short_edges()
        if "slits" in operations:
            for _i in range(len(pfix)):
                if not slit_once(hole_polygon(pfix, egfix)):
                    break
        if "angles" in operations and focus_pts is None:
            angles(hole_polygon(pfix, egfix))
        if n_done() == n_before:
            break
    if protect and (ident[sorted(protect)] < 0).any():
        raise RuntimeError("rim_repair removed a protected point")
    # what is left is read off the rim returned, not the refusals on the
    # way: a later blunting took an edge the report still named (round 12)
    left = []
    if len(egfix):
        L = np.linalg.norm(pfix[egfix[:, 0]] - pfix[egfix[:, 1]], axis=1)
        mid = 0.5 * (pfix[egfix[:, 0]] + pfix[egfix[:, 1]])
        h = np.asarray(size(mid), dtype=float)
        left = [{"length_m": round(float(L[k]), 2), "h_m": round(float(h[k]), 1),
                 "at": [round(float(v), 1) for v in mid[k]]}
                for k in np.argsort(L / h) if L[k] < min_edge_factor * h[k]
                and in_focus(mid[k])]
    report = {"n_points_removed": len(removed), "removed": removed[:50],
              "n_tips_stepped_back": len(retreated), "stepped_back": retreated[:50],
              "n_short_edges_left": len(left), "short_edges_left": left[:50],
              "n_short_edges_refused": len(refused),
              "n_slits_closed": len(slits), "slits": slits[:50],
              "slits_left": slits_left[:50],
              "n_corners_blunted": sum(r["n_corners_blunted"] for r in angles_rep),
              "corners": [a for r in angles_rep for a in r["at"]][:50]}
    return pfix, egfix, pfix_base, ident, report


def _covering_grid(geom, step: float):
    """Grid points of spacing ``step`` within ``c`` of ``geom``, and ``c``.

    Every point of ``geom`` lies within ``c = step / sqrt(2)`` of one of
    them.  Keeping only the points INSIDE ``geom`` loses that along its
    edges, where the nearest grid point is outside (review, round 5: a
    100 x 1 m box kept none, and its vertices read 47 m from a 30 m spot).
    """
    import shapely

    c = step / np.sqrt(2.0)
    x0, y0, x1, y1 = geom.bounds
    gx, gy = np.meshgrid(np.arange(x0, x1 + step, step), np.arange(y0, y1 + step, step))
    grid = np.column_stack([gx.ravel(), gy.ravel()])
    shapely.prepare(geom)
    keep = shapely.dwithin(geom, shapely.points(grid[:, 0], grid[:, 1]), c * (1 + 1e-9))
    return grid[np.asarray(keep, dtype=bool)], c


def _size_bounds(geom, size, h_floor: float, slope: float) -> tuple[float, float]:
    """``(lo, hi)`` of ``size`` over ``geom``; see :func:`size_lower_bound`."""
    import shapely

    if not (np.isfinite(h_floor) and h_floor > 0):
        raise ValueError("h_floor must be finite and positive")
    if not (np.isfinite(slope) and slope >= 0):
        raise ValueError("slope must be finite and non-negative")
    step = 0.25 * h_floor
    own = getattr(size, "size_bounds", None)
    if own is not None:
        lo, hi = own(geom, step)
        if not (np.isfinite(lo) and np.isfinite(hi) and 0 < lo <= hi):
            raise ValueError("the size field is not finite and positive over the geometry")
        return lo, hi
    grid, c = _covering_grid(geom, step)
    pts = np.vstack([q for q in (shapely.get_coordinates(geom)[:, :2], grid) if len(q)])
    h = np.asarray(size(pts), dtype=float).reshape(-1)
    if h.shape != (len(pts),) or not np.isfinite(h).all() or (h <= 0).any():
        raise ValueError("the size field is not finite and positive over the geometry")
    return float(h.min() - slope * c), float(h.max() + slope * c)


def size_lower_bound(geom, size, h_floor: float, *, slope: float = 1.0,
                     certified: bool = True) -> float:
    """A size no point of ``geom`` is below -- guaranteed, not sampled.

    A field that knows its own structure says so through a ``size_bounds``
    attribute (:func:`patch_sizing` and :func:`base_size_field` do: their
    bounds come from element vertices and candidate nearest nodes, and hold
    without any slope premise -- the ``nearest`` extension is discontinuous,
    review, round 5).  Any other field is read on a grid of ``h_floor / 4``
    that COVERS ``geom`` -- every point of it within half a diagonal of a
    reading -- and on its vertices, and the smallest reading is lowered by
    ``slope`` times that half diagonal: a field that changes by at most
    ``slope`` per metre cannot hide anything smaller between readings.
    ``h_floor`` -- the smallest size the field has anywhere -- is the answer
    when the readings allow no better (reviews, rounds 2-5).

    ``certified=False`` says ``h_floor`` is only an estimate (it then sets
    the grid step and nothing else): the bound is not raised to it, and may
    come out at or below zero, which every caller reads as "cannot prove
    absence" (review, round 7: half the smallest size at a pocket's vertices
    raised a bound over a 30 m spot to 95 m).
    """
    lo = _size_bounds(geom, size, h_floor, slope)[0]
    return float(max(h_floor, lo) if certified else lo)


def resolve_size_floor(size, size_floor=None, points=None) -> tuple[float, bool]:
    """``(floor, certified)``: the smallest size ``size`` has anywhere.

    ``size_floor`` when the caller knows it; else the field's own
    ``size_min`` (:func:`patch_sizing`, :func:`base_size_field`); else only
    an estimate, half the smallest size at ``points``, which the bounds use
    for their grid step and never as a floor (review, round 7).
    """
    if size_floor is not None:
        f = float(size_floor)
        if not (np.isfinite(f) and f > 0):
            raise ValueError("size_floor must be finite and positive")
        return f, True
    m = getattr(size, "size_min", None)
    if m is not None:
        return float(m), True
    if points is not None and len(points):
        v = np.asarray(size(np.atleast_2d(np.asarray(points, dtype=float))[:, :2]),
                       dtype=float)
        if len(v) and np.isfinite(v).all() and (v > 0).all():
            return 0.5 * float(v.min()), False
    raise ValueError("no size floor: pass size_floor, or a field with size_min")


def size_upper_bound(geom, size, h_floor: float, *, slope: float = 1.0) -> float:
    """A size no point of ``geom`` is above; the mirror of
    :func:`size_lower_bound`, for a test that must hold at the COARSEST
    element on ``geom`` (an island at least one element in area, clear of
    the coast by half of one) -- a lower bound passes those too easily
    (review, round 5).
    """
    return float(_size_bounds(geom, size, h_floor, slope)[1])


def _polygons(geom):
    """Every polygon in ``geom``, however deeply collections nest them."""
    if geom is None or geom.is_empty:
        return
    if geom.geom_type == "Polygon":
        yield geom
    for g in getattr(geom, "geoms", ()):
        yield from _polygons(g)


def land_an_element_fits(land, size, h_min: float, *, factor: float = 0.5,
                         slope: float = 1.0, certified: bool = True) -> list:
    """Where ``land`` may hold a disc of ``factor`` local elements -- empty
    only when it certainly holds none.

    Absence is proved, not sampled: for each polygon a guaranteed lower bound
    of the size over it (:func:`size_lower_bound`) is found, and if eroding
    the polygon by ``factor`` times that bound leaves nothing, no disc of
    ``factor * size(c)`` fits anywhere.  Otherwise the polygon is reported --
    it may hold one, and an unresolved case must not pass (review, round 4:
    a sampled search missed an in-centre with 31 m of clearance against
    30 m, and one point's size missed a 25 m disc in round 3).  ``h_min`` is
    the smallest size the field has anywhere.  Returns one point per
    reported polygon.
    """
    if not (np.isfinite(factor) and factor > 0):
        raise ValueError("factor must be finite and positive")
    out = []
    # make_valid nests polygons in collections; every one is land (round 5)
    for g in _polygons(land):
        if certified and g.buffer(-factor * h_min).is_empty:
            continue                                   # nothing fits even at the floor
        lo = size_lower_bound(g, size, h_min, slope=slope, certified=certified)
        core = g.buffer(-factor * lo)
        if core.is_empty:
            continue
        out.append(np.round(np.asarray(core.representative_point().coords)[0], 1).tolist())
    return out


def _clear_of(ring, other, size, factor):
    """``(ok, gap, h, possible)``: whether every edge of the closed ``ring``
    is ``factor`` local elements away from ``other``; the tightest edge's gap
    and size; and whether every edge is at least ``factor`` of the FINEST of
    those sizes away -- a gap under that carries no element at all.

    Judged edge by edge at the size where the edge meets the gap -- the
    point of the edge nearest ``other``, bounded over the piece -- not at
    the finest size anywhere
    on the ring: a 30 m corner let a lake stand 50 m from a coast where the
    elements are 170 m (review, round 13).  Not at the middle of the gap or
    along it: from an edge on the far side the shortest line crosses the
    island itself, and on land the field takes a coarse node's value --
    both refused a Funabashi island 159 m clear of the rim that meshes
    cleanly.
    """
    import shapely

    ring = np.asarray(ring, dtype=float)[:, :2]
    # Pieces of a quarter of the finest size on the ring: one foot per long
    # edge lands anywhere along a parallel gap, and one fell at a 30 m end
    # of a gap that needs 200 m elements in its middle (review, round 14).
    step = 0.25 * float(np.min(np.asarray(size(ring), dtype=float)))
    dense = shapely.get_coordinates(shapely.segmentize(shapely.LinearRing(ring), step))[:-1]
    ring = dense if len(dense) >= 3 else ring
    segs = shapely.linestrings(np.stack([ring, np.roll(ring, -1, axis=0)], axis=1))
    gaps = shapely.shortest_line(segs, other)
    d = np.asarray(shapely.length(gaps), dtype=float)
    foot = shapely.get_coordinates(shapely.get_point(gaps, 0))
    # the largest size anywhere on the piece: the most read at its foot and
    # its ends, plus a slope of 1 m/m over its length -- samples alone left
    # a 31 m element between two 30 m readings (review, round 15).  A field
    # that jumps (steeper than that) is outside this bound.
    ends = np.vstack([ring, np.roll(ring, -1, axis=0)])
    h_ends = np.asarray(size(ends), dtype=float).reshape(2, -1)
    length = np.asarray(shapely.length(segs), dtype=float)
    h_foot = np.asarray(size(foot), dtype=float)
    h = np.maximum(h_foot, h_ends.max(axis=0)) + length
    k = int(np.argmin(d - factor * h))
    # impossible only on what was read: the finest element seen anywhere
    return (bool((d >= factor * h).all()), float(d[k]), float(h[k]),
            bool((d >= factor * float(min(h_foot.min(), h_ends.min()))).all()))


def island_rings(land, water, size, clearance_factor=0.5, *, fine_h=None):
    """The land wholly inside the water the patch meshes, as rings to add.

    The rim is cut from the BASE mesh's coastline, re-drawn along the source
    where the source runs; land the base never had -- a quay block standing
    in the harbour, left detached when the width filter turned the narrow
    pier that joined it to the shore into a wall -- has no stretch to be
    re-drawn from, and was meshed as water.  Every polygon of ``land`` that
    lies inside ``water`` with at least ``clearance_factor`` of a local
    element to spare (:func:`_clear_of`) is returned as a closed ring at the
    local size, and so is each lake in it, as water.  One closer than that
    is left out (reported in ``skipped``) when its gap could carry no
    element at all (under ``clearance_factor`` of the finest size along it)
    or no element fits on it; otherwise leaving it out would mesh it as
    water -- or fill a lake an element fits in as land -- so it is kept and
    named in ``tight``, and the QA gate judges the gap.

    Each polygon's outline -- resampled or the source's -- is chosen with
    its lakes, but that choice can make a NEIGHBOUR impossible: a sibling
    island, or an island standing in its lake (review, round 22).  So the
    whole pass is repeated, preferring the source outline for the placed
    rings next to anything refused, while that loses less (bounded by the
    number of polygons), and the pass losing least area is returned.

    Returns ``(rings, report)``; each ring is an (n, 2) array, not closed.
    """
    import shapely

    polys = sorted(_polygons(land), key=lambda q: -shapely.Polygon(q.exterior).area)
    src_keys = [frozenset((round(float(x), 6), round(float(y), 6)) for x, y in r.coords)
                for g in polys for r in (g.exterior, *g.interiors)]

    def owner(ring):
        # the polygon one of whose source rings this ring was cut from --
        # the nearest by Hausdorff distance; a resampled ring may stand a
        # little outside its source
        lr = shapely.LinearRing(ring)
        return min(range(len(polys)), key=lambda i: min(
            shapely.hausdorff_distance(lr, r)
            for r in (polys[i].exterior, *polys[i].interiors)))

    prefer: set = set()
    best = None
    for _pass in range(len(polys) + 1):
        rings, rep, refused = _island_pass(polys, water, size, clearance_factor,
                                           fine_h, prefer)
        loss = float(sum(q.area for q in refused))
        if best is None or loss < best[0] - 1e-6:
            best = (loss, rings, rep)
        if loss <= 0:
            break
        new = set()
        for q in refused:
            reach = 2.0 * float(np.max(_size_at(size, np.asarray(q.exterior.coords)[:, :2])))
            for r in rings:
                key = frozenset((round(float(x), 6), round(float(y), 6)) for x, y in r)
                if key in src_keys:
                    continue                      # already the source outline
                if shapely.LinearRing(r).distance(q) < reach:
                    i = owner(r)
                    if i not in prefer:
                        new.add(i)
        if not new:
            break
        prefer |= new
    _loss, rings, rep = best
    rep["n_passes"] = _pass + 1
    return rings, rep


def _preference_search(n, first, run, *, exhaustive_limit=8):
    """The best ``(cost, outcome)`` over per-lake preference tuples, and the
    last reason a group was invalid.

    ``run(prefs)`` returns ``((cost, outcome, ...), None)`` or ``(None,
    why)``.  Up to ``exhaustive_limit`` lakes every tuple is tried, ``first``
    ones first; beyond, a breadth-first search over single flips from
    ``first``: deduplicated when queued, charged per evaluation (4 n^2),
    expanding invalid groups and sideways moves alike -- a single-flip
    search stopped where two flips at once keep a lake, duplicates spent the
    budget, and an earlier shell's score stopped this one's search (reviews,
    rounds 21-22).  Stops at the first cost of 0.
    """
    import itertools
    from collections import deque

    best, why = None, None
    if n <= exhaustive_limit:
        order = list(first) + [o for o in itertools.product((True, False), repeat=n)
                               if o not in first]
        for prefs in order:
            got, why_g = run(prefs)
            if got is None:
                why = why_g
                continue
            if best is None or got[0] < best[0]:
                best = (got[0], got[1])
            if best[0] == 0:
                break
        return best, why
    queue, queued = deque(first), set(first)
    for _eval in range(4 * n * n):
        if not queue:
            break
        prefs = queue.popleft()
        got, why_g = run(prefs)
        if got is None:
            why = why_g
            expand = True
        else:
            if best is None or got[0] < best[0]:
                best = (got[0], got[1])
            if got[0] == 0:
                break
            expand = got[0] <= best[0]
        if expand:
            for i in range(n):
                nxt = prefs[:i] + (not prefs[i],) + prefs[i + 1:]
                if nxt not in queued:
                    queued.add(nxt)
                    queue.append(nxt)
    return best, why


def _island_pass(polys, water, size, clearance_factor, fine_h, prefer):
    """One pass of :func:`island_rings`: ``(rings, report, refused)``, where
    ``refused`` are the source islands and lakes it could not deliver;
    ``prefer`` holds the indices of ``polys`` to try source outline first.
    """
    import shapely

    # a constant, or a field of any output shape, as before (review, round
    # 14); the bounds a field states travel with it
    field_ = size

    def size(q):
        return _size_at(field_, q)

    for attr in ("size_min", "size_max", "size_bounds"):
        if hasattr(field_, attr):
            setattr(size, attr, getattr(field_, attr))
    if not callable(field_):
        size.size_min = size.size_max = float(field_)

    rings, skipped = [], []
    n_islands = n_lakes = 0
    edge = shapely.boundary(water)
    placed: list = []                   # (resampled polygon, source polygon, report point)
    filled: list = []                   # source lakes refused, meshed as land
    tight: dict = {}                    # placed index -> why it is kept although close
    refused: list = []                  # source islands and lakes not delivered
    floor, floor_ok = resolve_size_floor(
        size, None, np.vstack([shapely.get_coordinates(g)[:, :2] for g in polys])) \
        if polys else (1.0, False)

    def at(geom):
        c = np.asarray(geom.representative_point().coords)[0]
        return [round(float(c[0]), 1), round(float(c[1]), 1)]

    def apart(poly, src):
        """A ring may not touch one already placed (two lakes resampled
        apart came out overlapping, review round 13), and it must sit inside
        exactly the placed rings its source sits inside -- and they in it:
        two sibling lakes came out one inside the other (round 14)."""
        for q, q_src, _at in placed:
            if poly.exterior.intersects(q.exterior):
                return False
            if poly.within(q) != src.within(q_src) or q.within(poly) != q_src.within(src):
                return False
        return True

    def choose(cands, src, inside, fits, what):
        """The first candidate ring that may be delivered, as
        ``(ring, why_tight_or_None)``, or ``(None, why)``.

        Each candidate -- the resampled ring, then the source outline -- is
        judged as delivered (the source's 50 m gap came back 8 m after
        resampling, review round 16): a valid polygon ``inside`` its
        container, apart from every placed ring, and clear of the rim and of
        every placed ring by ``_clear_of`` -- and every placed ring clear of
        IT, since sizes differ on the two sides of a gap (round 17).  Too
        close for the local element but an element fits on it and the gap
        could carry one: kept, and named -- leaving out land an element fits
        on meshes it as water, which the driver refuses (a Funabashi island
        159 m from the rim where the element is 367 m meshes cleanly), and
        filling a lake an element fits in makes water land; the QA gate
        judges the gap.  A placed ring the new one makes tight is named too.
        """
        why = "no valid ring at the local size"
        others = shapely.union_all([edge, *[q.exterior for q, _s, _a in placed]])
        for cand in cands:
            if len(cand) < 3:
                continue
            poly = shapely.Polygon(cand)
            if not (poly.is_valid and inside(poly)
                    and not poly.exterior.intersects(others) and apart(poly, src)):
                continue
            ok, gap, h_gap, possible = _clear_of(cand, others, size, clearance_factor)
            if not ok:
                why = f"{what} {gap:.1f} m from the rim or another ring where the " \
                      f"element is {h_gap:.0f} m"
                if not possible or not fits():
                    continue
            back, fine = {}, True
            for i, (q, _s, _a) in enumerate(placed):
                q_ok, q_gap, q_h, q_possible = _clear_of(
                    np.asarray(q.exterior.coords)[:-1], poly.exterior, size, clearance_factor)
                if not q_possible:
                    fine = False
                    break
                if not q_ok:
                    back[i] = (q_gap, q_h)
            if not fine:
                why = f"{what}: a ring already placed would be too close to it"
                continue
            for i, (q_gap, q_h) in back.items():
                tight.setdefault(i, {"at": placed[i][2], "why": (
                    f"{q_gap:.1f} m from a ring added later where the element is "
                    f"{q_h:.0f} m")})
            return cand, (None if ok else why)
        return None, why

    # every polygon, however make_valid nested it; lines are no land (round
    # 6).  Outside in, so an island is met after the lake it stands in.
    # By the SHELL's area: net area leaves out the lakes, and an island in a
    # large lake came before the island around it (review, round 14).
    for gi, g in enumerate(polys):
        if not shapely.intersects(water, g):
            continue
        if any(f.contains(g) for f in filled):
            # its lake is land now: its own coast would turn it into water
            # by the nesting parity (review, round 13)
            skipped.append({"at": at(g), "why": "inside a lake meshed as land"})
            continue
        if not shapely.within(g, water):
            inside = float(shapely.area(shapely.intersection(water, g)))
            if inside > 0:
                # this land is meshed as WATER where it lies inside the hole;
                # the area says whether that is a sliver or a pier block
                skipped.append({"at": at(g), "why": "crosses the rim",
                                "area_inside_m2": round(inside, 1)})
            continue
        ext = np.asarray(g.exterior.coords, dtype=float)[:-1, :2]
        src = shapely.Polygon(g.exterior)
        lakes_src = [shapely.Polygon(h_) for h_ in g.interiors]
        lakes_u = shapely.union_all(lakes_src) if lakes_src else shapely.Polygon()

        def in_water(poly):
            # its outline in the water, and the island in the water outside
            # its own lakes: a lake may hold a rim island (review, round 18)
            return poly.exterior.within(water) and poly.difference(lakes_u).within(water)

        # The shell and its lakes are chosen TOGETHER: a resampled shell
        # committed first filled a lake the source shell keeps (review,
        # round 18), and one lake's outline made a sibling fail that another
        # choice keeps (rounds 19-21).  Each shell candidate is tried with
        # per-lake outline preferences, searched as described below, and the
        # group filling fewest lakes is kept.  A group
        # whose land (the shell minus the lakes it keeps) holds any of the
        # rim -- a rim island inside a filled lake turned to water, whatever
        # its size (rounds 19-20) -- is no group at all.
        state = (len(rings), len(placed), len(filled), len(skipped), dict(tight))
        why = "no valid ring at the local size"

        def reset():
            del rings[state[0]:], placed[state[1]:], filled[state[2]:], skipped[state[3]:]
            tight.clear()
            tight.update(state[4])

        def run_group(shell_cand, prefs):
            """The group for one shell and per-lake preferences, as
            ``(n_filled, outcome, filled_idx)``, or ``(None, why)``."""
            reset()
            r, why_r = choose((shell_cand,), src, in_water,
                              lambda: bool(land_an_element_fits(g, size, floor,
                                                                certified=floor_ok)),
                              "island")
            if r is None:
                return None, why_r
            if why_r:
                tight[len(placed)] = {"at": at(g), "why": why_r}
            rings.append(r)
            placed.append((shapely.Polygon(r), src, at(g)))
            shell = shapely.Polygon(r)
            kept, filled_idx = [], []
            for j, lake_p in enumerate(lakes_src):
                lake = np.asarray(lake_p.exterior.coords, dtype=float)[:-1, :2]
                resampled = _ring_at_size(lake, size, fine_h=fine_h)
                lr, why_l = choose(
                    (resampled, lake) if prefs[j] else (lake, resampled), lake_p,
                    lambda poly: shell.contains(poly),
                    lambda: bool(land_an_element_fits(lake_p, size, floor,
                                                      certified=floor_ok)),
                    "lake")
                if lr is None:
                    filled.append(lake_p)
                    skipped.append({"at": at(lake_p), "why": f"{why_l}; meshed as land",
                                    "lake_area_m2": round(float(lake_p.area), 1)})
                    filled_idx.append(j)
                    continue
                if why_l:
                    tight[len(placed)] = {"at": at(lake_p), "why": why_l}
                rings.append(lr)
                placed.append((shapely.Polygon(lr), lake_p, at(lake_p)))
                kept.append(shapely.Polygon(lr))
            land_now = shell.difference(shapely.union_all(kept)) if kept else shell
            if land_now.intersects(edge) or land_now.difference(water).area > 1e-6:
                return None, "island: a lake it cannot keep holds part of the rim"
            return (len(filled_idx), (list(rings[state[0]:]), list(placed[state[1]:]),
                                      list(filled[state[2]:]), list(skipped[state[3]:]),
                                      dict(tight), len(kept)), filled_idx), None

        n_l = len(lakes_src)
        first = [tuple([True] * n_l), tuple([False] * n_l)]
        shell_cands = (_ring_at_size(ext, size, fine_h=fine_h), ext)
        if gi in prefer:                 # the outer pass asks for the source
            first.reverse()
            shell_cands = shell_cands[::-1]
        best = None
        for shell_cand in shell_cands:
            got, why_s = _preference_search(
                n_l, first, lambda prefs, sc=shell_cand: run_group(sc, list(prefs)))
            if why_s:
                why = why_s
            if got is not None and (best is None or got[0] < best[0]):
                best = got
            if best is not None and best[0] == 0:
                break
        del rings[state[0]:], placed[state[1]:], filled[state[2]:], skipped[state[3]:]
        tight.clear()
        tight.update(state[4])
        if best is None:
            skipped.append({"at": at(g), "why": why})
            refused.append(src)
            continue
        b_rings, b_placed, b_filled, b_skipped, b_tight, b_kept = best[1]
        refused.extend(b_filled)
        rings.extend(b_rings)
        placed.extend(b_placed)
        filled.extend(b_filled)
        skipped.extend(b_skipped)
        tight.clear()
        tight.update(b_tight)
        n_islands += 1
        n_lakes += b_kept
    tight = [tight[i] for i in sorted(tight)]
    return rings, {"n_islands_added": n_islands, "n_lakes_added": n_lakes,
                   "skipped": skipped, "tight": tight,
                   "n_island_points": int(sum(len(r) for r in rings))}, refused


def _rim_layout(pfix, egfix):
    """``[(vertex keys, polygon, first vertex)]`` per rim ring, or None when
    the edges do not form closed rings of three points or more."""
    import shapely

    try:
        rings, _ = boundary_rings(pfix, egfix)
    except ValueError:
        return None
    xy = np.asarray(pfix, dtype=float)[:, :2]
    out = []
    for r in rings:
        if len(r) < 3:
            return None
        keys = {(round(float(x), 6), round(float(y), 6)) for x, y in xy[r]}
        out.append((keys, shapely.Polygon(xy[r]), xy[r[0]]))
    return out


def _rim_edit_ok(p_old, e_old, p_new, e_new) -> bool:
    """Whether an edit of the rim leaves it a rim with the same meaning.

    Every ring a simple polygon of some area (three points in a line passed
    every edge test once, review round 7); no two rings touching or crossing
    (a fold ran the coast through an island, round 11); and every ring --
    matched to its old self by the vertices they share -- inside exactly the
    same rings as before, so no land becomes water or water land (an island
    left outside its shell, round 9; a lake and an island that swapped roles
    under an unchanged depth count, round 11).
    """
    import shapely

    old, new = _rim_layout(p_old, e_old), _rim_layout(p_new, e_new)
    if old is None or new is None or len(old) != len(new):
        return False
    for _k, poly, _v in new:
        if not (poly.is_valid and poly.area > 0):
            return False
    lines = [poly.exterior for _k, poly, _v in new]
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            if lines[i].intersects(lines[j]):
                return False
    match = []
    for keys, _p, _v in new:
        shared = [len(keys & k0) for k0, _p0, _v0 in old]
        j = int(np.argmax(shared))
        if shared[j] == 0:
            return False
        match.append(j)
    if len(set(match)) != len(match):
        return False

    def containers(layout):
        return [{j for j, (_k, q, _v) in enumerate(layout)
                 if j != i and q.contains(shapely.Point(v))}
                for i, (_k, _p, v) in enumerate(layout)]

    c_old, c_new = containers(old), containers(new)
    return all({match[j] for j in c_new[i]} == c_old[match[i]] for i in range(len(new)))


def blunt_acute_corners(pfix, egfix, pfix_base, water, size, min_angle_deg=60.0, *,
                        protect=()):
    """Cut off every coastline corner sharper than ``min_angle_deg`` on the water side.

    The water between two boundary lines that meet at under 60 degrees holds
    ONE element if every angle is to stay at 30 or more, so the corner node is
    in a single element -- and FVCOM never updates such a node: on the walled
    Kimitsu harbour a 44 deg corner of the resolved coastline kept an M2
    amplitude of exactly 0, and bisecting its element left angles of 22 deg.

    From the corner the coastline is walked one local element ``size`` along
    each side, and the two points reached are joined by a chord; the corner
    and any point passed on the way are removed.  Walking a fixed fraction of
    the first edge instead left a 7.4 m chord where that edge was short, and
    the external step fell from 2.46 to 1.37 s.  The walk stops early at a
    point it may not remove -- a frozen one (``pfix_base >= 0``) or one with
    other than two constrained edges -- and, so that the chord stays inside
    the water, at 0.9 of a side it would otherwise run past.

    ``water`` is the hole the rim bounds; ``size`` maps (n, 2) points to the
    local element size.  Returns ``(pfix, egfix, pfix_base, report)``.
    """
    import shapely

    pfix = np.asarray(pfix, dtype=float).copy()
    egfix = np.asarray(egfix, dtype=np.int64).copy()
    pfix_base = np.asarray(pfix_base, dtype=np.int64).copy()
    # ``protect``: indices never removed (wall roots); held by position,
    # because the indices shift as points go
    protected = {tuple(np.round(pfix[int(k)], 6)) for k in protect}

    def cut_one():
        """Find ONE acute water corner in the current rim and cut it.

        One at a time, against the rim as it now stands: planning every cut
        on the original adjacency let a later cut delete an end point an
        earlier cut's new edges used, and two adjacent 45 deg corners came
        back with edges to node -1 (review F1).
        """
        nonlocal pfix, egfix, pfix_base
        nbrs: dict = {}
        for a, b in egfix.tolist():
            nbrs.setdefault(a, []).append(b)
            nbrs.setdefault(b, []).append(a)

        def free(k):
            return (pfix_base[k] < 0 and len(nbrs.get(k, [])) == 2
                    and tuple(np.round(pfix[k], 6)) not in protected)

        def bend(k, prev):
            nxt = [q for q in nbrs[k] if q != prev][0]
            a, b = pfix[k] - pfix[prev], pfix[nxt] - pfix[k]
            return float(np.degrees(np.arccos(np.clip(
                a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12), -1, 1))))

        def walk(v, first, d):
            """Walk from v through `first` for arc length d.

            Returns (point, stop, passed, at_vertex). The walk stops AT a
            corner (a bend over 30 deg) it reaches: passing one put the chord
            over land on the far side of it, and the Odaiba 57 deg corner was
            left uncut because of it.
            """
            prev, cur, left, passed = v, first, d, []
            while True:
                seg = float(np.linalg.norm(pfix[cur] - pfix[prev]))
                if seg >= left or not free(cur):
                    t = min(left, 0.9 * seg) / seg
                    return pfix[prev] + t * (pfix[cur] - pfix[prev]), cur, passed, False
                if bend(cur, prev) > 30.0:
                    return pfix[cur].copy(), cur, passed, True
                left -= seg
                passed.append(cur)
                nxt = [k for k in nbrs[cur] if k != prev][0]
                prev, cur = cur, nxt

        def water_angle(v, p, n):
            up, un = pfix[p] - pfix[v], pfix[n] - pfix[v]
            lp, ln = float(np.linalg.norm(up)), float(np.linalg.norm(un))
            if lp <= 0 or ln <= 0:
                return None
            ang = float(np.degrees(np.arccos(np.clip(up @ un / (lp * ln), -1, 1))))
            bis = up / lp + un / ln
            if np.linalg.norm(bis) < 1e-9:
                return None
            probe = pfix[v] + bis / np.linalg.norm(bis) * 0.05 * min(lp, ln)
            return ang if water.contains(shapely.Point(probe)) else 360.0 - ang

        def open_frozen(v):
            nonlocal pfix, egfix, pfix_base
            p, n = nbrs[v]
            ang = water_angle(v, p, n)
            if ang is None or ang >= min_angle_deg:
                return None
            for drop, keep in ((p, n), (n, p)):
                if not free(drop):
                    continue
                nxt = [k for k in nbrs[drop] if k != v][0]
                if nxt == keep:
                    continue
                chord = shapely.LineString([pfix[v], pfix[nxt]])
                # Opening the angle gives the triangle v-drop-nxt to the
                # water, so the chord runs over land by construction; what it
                # may not do is cross another stretch of the rim.
                others = [shapely.LineString(pfix[[a, b]]) for a, b in egfix.tolist()
                          if not ({a, b} & {v, drop, nxt})]
                if others and shapely.MultiLineString(others).intersects(chord):
                    continue
                opened = water_angle(v, nxt, keep)
                if opened is None or opened < min_angle_deg:
                    continue            # a point on a straight run opens nothing
                keep_e = ~np.isin(egfix, [drop]).any(axis=1)
                eg_try = np.vstack([egfix[keep_e], [[v, nxt]]])
                if not _rim_edit_ok(pfix, egfix, pfix, eg_try):
                    continue            # three points in a line (review, round 7)
                egfix = eg_try
                live = np.setdiff1d(np.arange(len(pfix)), [drop])
                remap = np.full(len(pfix), -1, dtype=np.int64)
                remap[live] = np.arange(len(live))
                pfix, pfix_base, egfix = pfix[live], pfix_base[live], remap[egfix]
                return ang, pfix[remap[v]].copy(), float(chord.length), 1
            return None

        for v in range(len(pfix)):
            if pfix_base[v] >= 0 and len(nbrs.get(v, [])) == 2:
                # A FROZEN corner cannot be cut, but its angle can be opened:
                # where resolved coast meets frozen coast at under 60 deg, the
                # free point next to the junction goes and the junction joins
                # the point after it (Odaiba: 57 deg -> 115 deg; left alone
                # the corner held one element, split into two of 28.5 deg).
                got = open_frozen(v)
                if got is not None:
                    return got
                continue
            if not free(v):
                continue
            p, n = nbrs[v]
            up, un = pfix[p] - pfix[v], pfix[n] - pfix[v]
            lp, ln = float(np.linalg.norm(up)), float(np.linalg.norm(un))
            if lp <= 0 or ln <= 0:
                continue
            ang = float(np.degrees(np.arccos(np.clip(up @ un / (lp * ln), -1, 1))))
            # the angle between the edges is the WATER's only if the bisector
            # points into the water; otherwise the water has 360 minus it
            bis = up / lp + un / ln
            if ang >= min_angle_deg or np.linalg.norm(bis) < 1e-9:
                continue
            probe = pfix[v] + bis / np.linalg.norm(bis) * 0.05 * min(lp, ln)
            if not water.contains(shapely.Point(probe)):
                continue
            d = float(size(pfix[v][None])[0])
            a_xy, a_stop, a_pass, a_at = walk(v, p, d)
            b_xy, b_stop, b_pass, b_at = walk(v, n, d)
            removed = {v, *a_pass, *b_pass}
            if a_stop in removed or b_stop in removed or a_stop == b_stop:
                continue
            # the ring must keep three points: blunting a small triangle
            # took it down to two, and the next hole_polygon failed
            # (review, round 1)
            ring, todo = {v}, [v]
            while todo:
                for q in nbrs.get(todo.pop(), []):
                    if q not in ring:
                        ring.add(q)
                        todo.append(q)
            n_new = sum(1 for at_ in (a_at, b_at) if not at_)
            if len(ring) - len(removed) + n_new < 3:
                continue
            chord = shapely.LineString([a_xy, b_xy])
            if not water.buffer(1e-6).contains(chord):
                continue
            # a chord end AT a corner is that corner, not a new point beside it
            ia = len(pfix)
            new_xy, ends = [], []
            for xy_, stop_, at_ in ((a_xy, a_stop, a_at), (b_xy, b_stop, b_at)):
                if at_:
                    ends.append((stop_, None))
                else:
                    ends.append((stop_, ia + len(new_xy)))
                    new_xy.append(xy_)
            path = []
            for stop_, idx in ends:
                path.append(idx if idx is not None else stop_)
            new_e = []
            (sa, ia_), (sb, ib_) = ends
            if ia_ is not None:
                new_e.append([sa, ia_])
            new_e.append([path[0], path[1]])
            if ib_ is not None:
                new_e.append([ib_, sb])
            keep_e = ~np.isin(egfix, list(removed)).any(axis=1)
            eg_try = np.vstack([egfix[keep_e], np.asarray(new_e, dtype=np.int64)])
            p_try = np.vstack([pfix, np.asarray(new_xy)]) if new_xy else pfix
            if not _rim_edit_ok(pfix, egfix, p_try, eg_try):
                continue
            egfix = eg_try
            if new_xy:
                pfix = p_try
                pfix_base = np.concatenate([pfix_base, np.full(len(new_xy), -1)])
            live = np.setdiff1d(np.arange(len(pfix)), list(removed))
            remap = np.full(len(pfix), -1, dtype=np.int64)
            remap[live] = np.arange(len(live))
            pfix, pfix_base, egfix = pfix[live], pfix_base[live], remap[egfix]
            return ang, a_xy, float(chord.length), len(removed)
        return None

    corners = []
    n_removed = 0
    for _ in range(max(1, len(pfix))):
        got = cut_one()
        if got is None:
            break
        corners.append(got[:3])
        n_removed += got[3]
    if (egfix < 0).any():               # the invariant F1 broke; never hand it on
        raise RuntimeError("blunt_acute_corners left an edge to a deleted point")
    report = {"n_corners_blunted": len(corners),
              "angles_deg": [round(c[0], 1) for c in corners],
              "chord_m": [round(c[2], 1) for c in corners],
              "n_points_removed": n_removed,
              "at": [[round(float(x), 1) for x in c[1]] for c in corners]}
    return pfix, egfix, pfix_base, report


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
    try:
        out = shapely.union_all(shells)
    except Exception as exc:
        # GEOS says "side location conflict at <x> <y>", which is a sentence
        # about topology and not about the mesh.  A rim ring that is not
        # simple is what produces it, and a shoreline feature narrower than
        # the local element size is what produces that.
        bad = [k for k, p in enumerate(polys) if not p.is_valid
               or not p.exterior.is_simple]
        raise ValueError(
            f"the rim does not bound a valid polygon ({exc}); "
            f"{len(bad)} of {len(polys)} ring(s) are not simple. A coastline "
            "feature narrower than the local element size is the usual cause "
            "-- notebooks/423_site_survey.py measures that by eroding the "
            "water polygon") from exc
    if out.is_empty:
        raise ValueError("the rim segments do not close a polygon")
    if not out.is_valid:
        # union_all can hand back an invalid polygon rather than raise, and
        # a repair that made one was accepted (review, round 8)
        raise ValueError("the rim does not bound a valid polygon ("
                         f"{shapely.is_valid_reason(out)})")
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
    outside: str = "max",
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
    import shapely

    base_size = base_size_field(nodes, elements, outside=outside)
    geoms = [_region_geometry(r) for r in regions]

    def h(points):
        p = np.atleast_2d(np.asarray(points, dtype=float))[:, :2]
        base = base_size(p)
        out = base.copy()
        for g in geoms:
            out = np.minimum(out, _region_contribution(g, p, base))
        return np.minimum(out, base) / distmesh_scale

    def size_bounds(geom, step):
        # A region contributes target + (base - target) * u, u = clip(d / width)
        # rising with the distance d from it: over ``geom`` u is at least its
        # value at the nearest point, and the base at least its own lower
        # bound, so where that bound is over the target the contribution is
        # at least target + (lo - target) * u_min, and elsewhere at least the
        # base bound itself.  Contributions only lower the base, whose upper
        # bound is the field's.  No slope premise (review, round 5).
        lo_b, hi = base_size.size_bounds(geom, step)
        lo = lo_b
        for _edge, poly, target, width, _pr in geoms:
            d = float(shapely.distance(poly, geom))
            u = min(d / width, 1.0) if width > 0 else float(d > 0)
            if lo_b > target:
                lo = min(lo, float(target) + (lo_b - float(target)) * u)
        return lo / distmesh_scale, hi / distmesh_scale

    h.size_bounds = size_bounds
    # a contribution lies between its target and the base size
    h.size_min = min([base_size.size_min] + [float(g[2]) for g in geoms]) / distmesh_scale
    h.size_max = base_size.size_max / distmesh_scale
    return h


def base_size_field(nodes, elements, *, outside: str = "max"):
    """``f(points) -> the base mesh's own edge length there``, in metres.

    The ambient field of :func:`ambient_size_field`, interpolated.  Shared by
    :func:`patch_sizing` and :func:`region_conflicts` so that what the report
    evaluates is the field the mesher is given, not a second formula that
    resembles it.

    ``outside`` decides what happens beyond the base mesh's footprint.

    ``max``      the field's maximum, which is the historical behaviour and
                 is safe while the hole stays inside the mesh -- and it does,
                 because the rim is built from base edges.
    ``nearest``  the nearest in-mesh value.  ``coastline: resolve`` moves the
                 boundary onto the source shoreline, so parts of the hole end
                 up OUTSIDE the base mesh, and there ``max`` is a cliff: the
                 hires runs measured a field slope of 35-42 against a C4
                 reference of 0.414, where the same patch on the default
                 branch measures 0.375.  A mesher handed a size field with a
                 step resolves it with slivers, and that is where the
                 7.30-degree elements at Kimitsu port came from.
    """
    from matplotlib.tri import LinearTriInterpolator, Triangulation

    if outside not in ("max", "nearest"):
        raise ValueError(f"outside must be 'max' or 'nearest', got {outside!r}")
    xy = np.asarray(nodes, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    amb = ambient_size_field(xy, tri)
    interp = LinearTriInterpolator(Triangulation(xy[:, 0], xy[:, 1], tri), amb)
    amb_max = float(np.nanmax(amb))
    tree = None
    if outside == "nearest":
        from scipy.spatial import cKDTree

        tree = cKDTree(xy)

    def f(points):
        p = np.atleast_2d(np.asarray(points, dtype=float))[:, :2]
        out = np.asarray(interp(p[:, 0], p[:, 1]), dtype=float)
        bad = ~np.isfinite(out)
        if not bad.any():
            return out
        if tree is None:
            return np.where(bad, amb_max, out)
        out = np.array(out, dtype=float)
        out[bad] = amb[tree.query(p[bad])[1]]
        return out

    index = {}

    def _index():
        import shapely
        from scipy.spatial import cKDTree

        if not index:
            index["tri"] = shapely.STRtree(shapely.polygons(xy[tri]))
            e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
            ue, cnt = np.unique(e, axis=0, return_counts=True)
            index["rim"] = shapely.STRtree(shapely.linestrings(xy[ue[cnt == 1]]))
            index["nodes"] = tree if tree is not None else cKDTree(xy)
        return index

    def size_bounds(geom, step):
        """``(lo, hi)``: bounds of this field over ``geom``, proved.

        Inside the mesh the field is linear on each element, so it lies
        between the values at the vertices of the elements ``geom`` meets.
        Outside, ``nearest`` takes one node's value and is discontinuous
        across the nodes' Voronoi edges -- no slope bounds it (review, round
        5) -- so every node that can be nearest to a point of ``geom`` is
        found: each such point lies within ``c`` of a grid point ``g`` (grid
        of ``step``, half-diagonal ``c``) that is either outside the mesh or
        within ``c`` of its boundary, and its nearest node lies within
        ``d(g) + 2c`` of ``g``.
        """
        import shapely

        ix = _index()
        vals = []
        hit = ix["tri"].query(geom, predicate="intersects")
        if len(hit):
            vals.append(amb[tri[hit].ravel()])
        g, c = _covering_grid(geom, step)
        if len(g):
            gp = shapely.points(g[:, 0], g[:, 1])
            inside = np.isfinite(np.asarray(interp(g[:, 0], g[:, 1]), dtype=float))
            near_rim = np.zeros(len(g), dtype=bool)
            near_rim[np.unique(ix["rim"].query(gp, predicate="dwithin", distance=c)[0])] = True
            out = ~inside | near_rim
            if out.any():
                if tree is None:
                    vals.append(np.array([amb_max]))
                else:
                    go = g[out]
                    d = ix["nodes"].query(go)[0]
                    hits = ix["nodes"].query_ball_point(go, d + 2.0 * c)
                    vals.append(amb[np.unique(np.concatenate(
                        [np.asarray(k, dtype=np.int64) for k in hits]))])
        v = np.concatenate(vals) if vals else np.zeros(0)
        v = v[np.isfinite(v)]
        if not len(v):
            raise ValueError("the size field has no value over the geometry")
        return float(v.min()), float(v.max())

    f.size_bounds = size_bounds
    # certified global bounds: every value is a node's (inside, a convex
    # combination of three; outside, one node's, or the maximum)
    f.size_min, f.size_max = float(np.nanmin(amb)), amb_max
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
    # Which triangle each sample falls in, asked of the triangles directly.
    # matplotlib's trapezoid map was used first and refuses two things a
    # mesh with walls legitimately has: coincident nodes (a split puts two
    # nodes at one point) and, once the repair slides a wall node on one
    # side only, two boundary polylines on one line whose vertices no longer
    # match.  A spatial index over the triangles minds neither.
    polys = shapely.polygons(xy[tri])
    tree = shapely.STRtree(polys)
    pi, ti = tree.query(shapely.points(p[:, 0], p[:, 1]), predicate="intersects")
    found = np.full(p.shape[0], -1, dtype=np.int64)
    # first hit per sample; a sample on a shared edge belongs to either side
    order = np.argsort(pi, kind="stable")
    pi, ti = pi[order], ti[order]
    first = np.r_[True, pi[1:] != pi[:-1]] if len(pi) else np.zeros(0, bool)
    found[pi[first]] = ti[first]
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
    copy_of=None,
) -> dict[str, Any]:
    """Check the frozen zone really is frozen, through the node map.

    ``copy_of`` declares the nodes a wall split duplicated on purpose
    (:func:`fvcom_mesh_tools.walls.split_along_walls`): a coincident pair is
    a defect unless both are copies of the same node.

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
        if copy_of is not None and len(pairs):
            co = np.asarray(copy_of, dtype=np.int64)
            pairs = pairs[co[pairs[:, 0]] != co[pairs[:, 1]]]
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
    # One piece: a cut may leave the retained mesh in pieces the fill joins
    # (select_patch), and this is where "joins" is checked.
    n_components = int(_face_components(tri).max()) + 1 if len(tri) else 0

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
        # where, so a failed contract can be looked at: the counts alone
        # sent Yokohama's one missing edge per seed to a guessing game
        "unexpected_boundary_edges_at": [
            new_xy[list(e)].mean(axis=0).round(1).tolist() for e in sorted(unexpected)[:10]],
        "missing_boundary_edges_at": [
            new_xy[list(e)].mean(axis=0).round(1).tolist() for e in sorted(absent)[:10]],
        "n_extra_faces": int(extra),
        "n_nonmanifold_edges": nonmanifold,
        "n_components": n_components,
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
                   and dup == 0 and orphan == 0 and obc_ok and not split
                   and n_components == 1),
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
            if val[a] <= 2 or val[b] <= 2:
                # a node left in ONE element is one FVCOM never updates
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
            if val[a] <= 2 or val[b] <= 2:
                # a node left in ONE element is one FVCOM never updates
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
