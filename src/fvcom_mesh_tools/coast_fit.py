"""Pull the mesh land boundary onto the coastline it was meant to follow.

Why this exists
---------------
DistMesh only projects a point back to the domain boundary when the point has
stepped *outside* it; a point that stays just inside is left alone.  The mesh
boundary is then whatever polyline is exposed after the triangles whose
centroid falls on land are deleted.  Nothing in that procedure makes a
boundary node sit *on* the coastline, so the mesh shoreline zig-zags across
the real one -- measured on the certified Tokyo Bay mesh: median offset 10 m,
p90 84 m, worst 566 m against the polygon oceanmesh was given, with 32 % of
the boundary nodes on the land side.

The gap is a pure mesh-generation artefact: the coastline is known exactly, so
each boundary node has an obvious target.  This module moves the nodes there,
under guards strict enough that the fit can never trade shoreline fidelity for
a QA failure.

What it does NOT do
-------------------
No node is inserted or deleted and no edge is flipped, so element count,
connectivity, the one-element-wide ledger and the land-breach census are
unchanged by construction.  Fidelity below the local edge length therefore
stays out of reach -- that is a sizing question, not a fitting one.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
import shapely

__all__ = [
    "CoastFitResult",
    "boundary_edges",
    "boundary_offsets",
    "fit_boundary_to_coast",
]

# QA gate values the fit must never break (qa.py: c1/c2/c4).
MIN_ANGLE_DEG = 30.0
MAX_ANGLE_DEG = 130.0
MAX_AREA_CHANGE = 0.5
GRAVITY_M_S2 = 9.81
DT_DEPTH_FLOOR_M = 0.1


def boundary_edges(elements: np.ndarray) -> np.ndarray:
    """Return the ``(m, 2)`` node pairs that bound the triangulation.

    A boundary edge is one used by exactly one element.  Node ids are 0-based
    and each pair is sorted.
    """
    tri = np.asarray(elements, dtype=np.int64)
    e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    uniq, cnt = np.unique(e, axis=0, return_counts=True)
    return uniq[cnt == 1]


def _signed_areas(xy: np.ndarray, tri: np.ndarray) -> np.ndarray:
    a, b, c = xy[tri[:, 0], :2], xy[tri[:, 1], :2], xy[tri[:, 2], :2]
    return 0.5 * ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
                  - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))


def _angles_deg(xy: np.ndarray, tri: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Smallest and largest interior angle of every triangle, in degrees."""
    p = xy[tri][:, :, :2]
    lo = np.full(len(tri), 180.0)
    hi = np.zeros(len(tri))
    for k in range(3):
        u = p[:, (k + 1) % 3] - p[:, k]
        v = p[:, (k + 2) % 3] - p[:, k]
        den = np.maximum(np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1), 1e-30)
        ang = np.degrees(np.arccos(np.clip((u * v).sum(axis=1) / den, -1.0, 1.0)))
        lo = np.minimum(lo, ang)
        hi = np.maximum(hi, ang)
    return lo, hi


def _implied_dt(xy: np.ndarray, tri: np.ndarray, depth: np.ndarray) -> np.ndarray:
    """External-mode time step allowed by each element.

    ``shortest edge / sqrt(g * deepest node)`` -- the same expression the QA
    advisory reports (qa.py: ``implied_dt``), so a guard written against it
    protects exactly the number the gate prints.
    """
    p = xy[tri][:, :, :2]
    lmin = np.min([np.linalg.norm(p[:, (k + 1) % 3] - p[:, k], axis=1) for k in range(3)],
                  axis=0)
    h = np.maximum(depth[tri].max(axis=1), DT_DEPTH_FLOOR_M)
    return lmin / np.sqrt(GRAVITY_M_S2 * h)


def _edge_neighbours(tri: np.ndarray) -> np.ndarray:
    """``(m, 3)`` element-to-element map across edges; ``-1`` where none."""
    n = len(tri)
    pairs = np.empty((3 * n, 2), dtype=np.int64)
    for k, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
        pairs[k * n:(k + 1) * n] = np.sort(tri[:, [i, j]], axis=1)
    owner = np.tile(np.arange(n), 3)
    slot = np.repeat(np.arange(3), n)
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    pairs, owner, slot = pairs[order], owner[order], slot[order]
    nb = np.full((n, 3), -1, dtype=np.int64)
    same = np.all(pairs[1:] == pairs[:-1], axis=1)
    for k in np.flatnonzero(same):
        nb[owner[k], slot[k]] = owner[k + 1]
        nb[owner[k + 1], slot[k + 1]] = owner[k]
    return nb


def _nearest_on(line: Any, pts: np.ndarray) -> np.ndarray:
    """Nearest point of ``line`` for every row of ``pts`` (vectorised)."""
    P = shapely.points(pts[:, 0], pts[:, 1])
    seg = shapely.shortest_line(P, line)
    out = pts.copy()
    ok = ~shapely.is_missing(seg)
    if not ok.any():
        return out
    coords = shapely.get_coordinates(seg[ok]).reshape(-1, 2, 2)
    out[ok] = coords[:, 1, :]
    return out


def boundary_offsets(
    nodes: np.ndarray,
    elements: np.ndarray,
    land: Any,
    *,
    exclude: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Signed distance from each land boundary node to the coastline.

    ``land`` is the land polygon in the mesh CRS.  The sign is ``+`` when the
    node is in water (short of the coast) and ``-`` when it has crossed onto
    land.  Returns ``(node_ids, signed_offset_m)``.
    """
    xy = np.asarray(nodes, dtype=float)[:, :2]
    be = boundary_edges(elements)
    drop = set() if exclude is None else set(np.asarray(exclude).ravel().tolist())
    keep = np.array([not (a in drop and b in drop) for a, b in be], dtype=bool)
    ids = np.unique(be[keep].ravel()) if keep.any() else np.empty(0, dtype=np.int64)
    if ids.size == 0:
        return ids, np.empty(0)
    P = shapely.points(xy[ids, 0], xy[ids, 1])
    d = shapely.distance(P, land.boundary)
    return ids, np.where(shapely.contains(land, P), -1.0, 1.0) * d


@dataclass
class CoastFitResult:
    """Outcome of :func:`fit_boundary_to_coast`."""

    nodes: np.ndarray
    n_boundary: int = 0
    n_moved: int = 0
    n_rejected: int = 0
    before: dict[str, float] = field(default_factory=dict)
    after: dict[str, float] = field(default_factory=dict)
    min_angle_before_deg: float = float("nan")
    min_angle_after_deg: float = float("nan")
    max_angle_before_deg: float = float("nan")
    max_angle_after_deg: float = float("nan")
    max_move_m: float = 0.0
    dt_before_s: float = float("nan")
    dt_after_s: float = float("nan")
    sweeps: int = 0

    def summary(self) -> str:
        return (f"coast fit: {self.n_moved}/{self.n_boundary} boundary nodes moved "
                f"({self.n_rejected} rejections), |offset| median "
                f"{self.before.get('median', float('nan')):.1f} -> "
                f"{self.after.get('median', float('nan')):.1f} m, p90 "
                f"{self.before.get('p90', float('nan')):.1f} -> "
                f"{self.after.get('p90', float('nan')):.1f} m, max "
                f"{self.before.get('max', float('nan')):.1f} -> "
                f"{self.after.get('max', float('nan')):.1f} m; angles "
                f"[{self.min_angle_before_deg:.2f}, {self.max_angle_before_deg:.2f}] -> "
                f"[{self.min_angle_after_deg:.2f}, {self.max_angle_after_deg:.2f}] deg; "
                f"dt {self.dt_before_s:.2f} -> {self.dt_after_s:.2f} s; "
                f"largest move {self.max_move_m:.1f} m")

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_boundary": self.n_boundary,
            "n_moved": self.n_moved,
            "n_rejected": self.n_rejected,
            "before": self.before,
            "after": self.after,
            "min_angle_deg": [self.min_angle_before_deg, self.min_angle_after_deg],
            "max_angle_deg": [self.max_angle_before_deg, self.max_angle_after_deg],
            "max_move_m": self.max_move_m,
            "dt_s": [self.dt_before_s, self.dt_after_s],
            "sweeps": self.sweeps,
        }


def _stats(off: np.ndarray) -> dict[str, float]:
    a = np.abs(off)
    if a.size == 0:
        return {"median": float("nan"), "p90": float("nan"),
                "p99": float("nan"), "max": float("nan"), "n_on_land": 0}
    return {
        "median": float(np.median(a)),
        "p90": float(np.percentile(a, 90)),
        "p99": float(np.percentile(a, 99)),
        "max": float(a.max()),
        "n_on_land": int((off < 0).sum()),
    }


def fit_boundary_to_coast(
    nodes: np.ndarray,
    elements: np.ndarray,
    land: Any,
    *,
    fixed: np.ndarray | None = None,
    sweeps: int = 6,
    relax: float = 0.8,
    max_move_frac: float = 0.5,
    min_angle_deg: float = MIN_ANGLE_DEG,
    max_angle_deg: float = MAX_ANGLE_DEG,
    max_area_change: float = MAX_AREA_CHANGE,
    keep_centroids_wet: bool = True,
    depths: np.ndarray | None = None,
    dt_floor_s: float | None = None,
) -> CoastFitResult:
    """Move land boundary nodes onto ``land``'s boundary in guarded sweeps.

    Parameters
    ----------
    nodes, elements
        ``(n, 2+)`` coordinates in a projected CRS (metres) and ``(m, 3)``
        0-based triangles.
    land
        Land polygon in the same CRS -- the coastline the mesh is meant to
        represent.  Use the polygon handed to the mesh generator, not the raw
        source data: the difference between the two is preprocessing we chose
        deliberately and must not be undone here.
    fixed
        Node ids that must not move (open-boundary nodes, constrained points).
    sweeps, relax
        Each sweep moves a node ``relax`` of the way to its nearest coastline
        point, so the mesh relaxes towards the coast instead of snapping.
    max_move_frac
        A node never ends further than this many local edge lengths from where
        it started, measured over the whole run.
    min_angle_deg, max_angle_deg, max_area_change
        The C1/C2/C4 gates.  A move is rejected when it would push an incident
        element past one of them -- unless that element already violated the
        bound and the move improves it.
    keep_centroids_wet
        Reject a move that would put an incident element's centroid on land.
    depths, dt_floor_s
        Node depths and the external time step the mesh must keep.  Pulling a
        boundary node onto the coast can flatten a triangle, and the time step
        is set by the single worst one, so without this guard a fit that
        improves the shoreline everywhere can still cost several per cent of
        dt.  A move that would take an incident element below ``dt_floor_s``
        is rejected unless that element is already below it and improves.
        With ``depths`` given and no explicit floor, the floor is the mesh's
        own current minimum, i.e. the fit is not allowed to cost any dt at
        all; pass ``dt_floor_s=0`` to switch the guard off.

    Returns
    -------
    CoastFitResult
        Fitted coordinates plus before/after statistics.  ``elements`` is
        never modified.
    """
    xy0 = np.asarray(nodes, dtype=float)
    xy = xy0.copy()
    tri = np.asarray(elements, dtype=np.int64)
    coast = land.boundary
    shapely.prepare(land)

    ids, off0 = boundary_offsets(xy, tri, land, exclude=fixed)
    res = CoastFitResult(nodes=xy, n_boundary=int(ids.size))
    lo_all, hi_all = _angles_deg(xy, tri)
    res.min_angle_before_deg = float(lo_all.min())
    res.max_angle_before_deg = float(hi_all.max())
    if ids.size == 0:
        return res
    res.before = _stats(off0)

    frozen = set() if fixed is None else set(np.asarray(fixed).ravel().tolist())
    movable = np.array([i for i in ids if i not in frozen], dtype=np.int64)
    if movable.size == 0:
        res.after = res.before
        res.min_angle_after_deg = res.min_angle_before_deg
        res.max_angle_after_deg = res.max_angle_before_deg
        return res

    # Local edge length per movable node: the mean of its boundary edges.
    be = boundary_edges(tri)
    elen = np.linalg.norm(xy[be[:, 0], :2] - xy[be[:, 1], :2], axis=1)
    hsum = np.zeros(len(xy))
    hcnt = np.zeros(len(xy))
    np.add.at(hsum, be[:, 0], elen)
    np.add.at(hsum, be[:, 1], elen)
    np.add.at(hcnt, be[:, 0], 1.0)
    np.add.at(hcnt, be[:, 1], 1.0)
    budget = max_move_frac * (hsum[movable] / np.maximum(hcnt[movable], 1.0))

    inc: dict[int, list[int]] = {int(i): [] for i in movable}
    for e, t in enumerate(tri):
        for n in t:
            lst = inc.get(int(n))
            if lst is not None:
                lst.append(e)

    nb = _edge_neighbours(tri)
    area = np.abs(_signed_areas(xy, tri))
    dep = None if depths is None else np.asarray(depths, dtype=float)
    dt_all = None if dep is None else _implied_dt(xy, tri, dep)
    if dep is not None and dt_floor_s is None:
        # Default: keep the time step the mesh arrived with.  Fidelity is not
        # worth paying for in dt, and the sweep below reaches the coastline
        # just as well with the floor in place (measured: p90 offset 10.2 m
        # unguarded vs 11.8 m guarded, dt 15.44 s vs 16.38 s).
        dt_floor_s = float(dt_all.min())
    if dt_floor_s is not None and dt_floor_s <= 0:
        dt_all = None
    orient = np.sign(_signed_areas(xy, tri))

    def area_change_ok(elems: np.ndarray) -> bool:
        """C4 over every edge touching ``elems``, in both directions."""
        for e in elems:
            for o in nb[e]:
                if o < 0:
                    continue
                a, b = area[e], area[o]
                lg = max(a, b, 1e-30)
                if (lg - min(a, b)) / lg > max_area_change:
                    return False
        return True

    moved = np.zeros(len(movable), dtype=bool)
    rejected = 0
    for sweep in range(int(sweeps)):
        target = _nearest_on(coast, xy[movable, :2])
        step = relax * (target - xy[movable, :2])
        any_move = False
        for k, node in enumerate(movable):
            d = step[k]
            if not np.isfinite(d).all():
                continue
            cand = xy[node, :2] + d
            v = cand - xy0[node, :2]
            nv = float(np.linalg.norm(v))
            if nv > budget[k]:
                if nv <= 1e-12:
                    continue
                # Clip to the budget instead of skipping, so a node already at
                # the cap can still slide along the coast.
                cand = xy0[node, :2] + v * (budget[k] / nv)
            if np.linalg.norm(cand - xy[node, :2]) < 1e-9:
                continue
            es = np.asarray(inc[int(node)], dtype=np.int64)
            old = xy[node, :2].copy()
            old_area = area[es].copy()
            xy[node, :2] = cand

            sa = _signed_areas(xy, tri[es])
            ok = bool((np.sign(sa) == orient[es]).all()) and bool((sa != 0).all())
            if ok:
                lo, hi = _angles_deg(xy, tri[es])
                ok = (bool((lo >= np.minimum(min_angle_deg, lo_all[es])).all())
                      and bool((hi <= np.maximum(max_angle_deg, hi_all[es])).all()))
            if ok:
                area[es] = np.abs(sa)
                ok = area_change_ok(es)
                if not ok:
                    area[es] = old_area
            if ok and dt_all is not None:
                dt_new = _implied_dt(xy, tri[es], dep)
                ok = bool((dt_new >= np.minimum(dt_floor_s, dt_all[es])).all())
                if not ok:
                    area[es] = old_area
            if ok and keep_centroids_wet:
                cen = xy[tri[es]][:, :, :2].mean(axis=1)
                ok = not bool(shapely.contains(
                    land, shapely.points(cen[:, 0], cen[:, 1])).any())
                if not ok:
                    area[es] = old_area
            if not ok:
                xy[node, :2] = old
                rejected += 1
                continue
            lo_all[es], hi_all[es] = _angles_deg(xy, tri[es])
            if dt_all is not None:
                dt_all[es] = _implied_dt(xy, tri[es], dep)
            moved[k] = True
            any_move = True
        res.sweeps = sweep + 1
        if not any_move:
            break

    _, off1 = boundary_offsets(xy, tri, land, exclude=fixed)
    lo, hi = _angles_deg(xy, tri)
    res.after = _stats(off1)
    res.n_moved = int(moved.sum())
    res.n_rejected = int(rejected)
    res.min_angle_after_deg = float(lo.min())
    res.max_angle_after_deg = float(hi.max())
    res.max_move_m = float(np.linalg.norm(xy[:, :2] - xy0[:, :2], axis=1).max())
    if dep is not None:
        res.dt_before_s = float(_implied_dt(xy0, tri, dep).min())
        res.dt_after_s = float(_implied_dt(xy, tri, dep).min())
    res.nodes = xy
    return res
