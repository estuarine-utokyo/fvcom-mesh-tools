"""Extend a finished base mesh outward: sizing, merge and the frozen contract.

The owner's decision (2026-09-29): the Tokyo Bay base mesh is adopted **as
it is** and a wider mesh is made by meshing only the sea outside it.  The
base's open boundary becomes an interior line: its nodes are fixed in the
outer generation (pfix/egfix), the outer mesh is merged onto them, and the
base's nodes and elements come through byte for byte.

Sizing (all in metres on a lon/lat lattice with metric ``x``, ``y``):

* the ambient field (coast distance, max edge) is limited to the gradation;
* the time-step floor ``dt * sqrt(g * depth) / Cr`` is raised into a field
  that is itself gradation-feasible (a graded dilation), and the ambient
  field may not go below it -- the extension must not be what limits dt;
* along each constrained line (the base interface, the new open boundary)
  the size is **set** to the line's own spacing on a band, and graded away
  from it both upward and downward.

A maximum or minimum of gradation-feasible fields is gradation-feasible, so
the result needs no further limiting.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from fvcom_mesh_tools.io.fort14 import Fort14Mesh
from fvcom_mesh_tools.sizing import _limit

__all__ = [
    "band_field",
    "check_no_overlap",
    "land_segments",
    "compose_sizing",
    "graded_up",
    "merge_outer",
    "rfactor_smooth_free",
    "trim_lone_corners",
    "verify_frozen_base",
]


#: Largest relative departure of a band from its target sizes that is accepted.
BAND_TOLERANCE = 0.05


def graded_up(values, x, y, grade):
    """The smallest gradation-feasible field that is >= ``values``."""
    return -_limit(-np.asarray(values, float), x, y, grade)


def band_field(x, y, line_xy, targets, half_width_m):
    """``targets`` along a polyline, spread over a band; NaN elsewhere.

    ``line_xy`` is the line in the lattice's metric coordinates and
    ``targets`` one size per vertex; a lattice point within ``half_width_m``
    of the line takes the target interpolated at its projection.
    """
    import shapely

    line_xy = np.asarray(line_xy, float)
    targets = np.asarray(targets, float)
    if len(line_xy) < 2 or len(targets) != len(line_xy):
        raise ValueError("a band needs a line of two points or more and one target per point")
    line = shapely.LineString(line_xy)
    pts = shapely.points(np.ravel(x), np.ravel(y))
    d = shapely.distance(pts, line)
    near = d <= half_width_m
    s_vertex = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(line_xy, axis=0), axis=1))]
    out = np.full(np.size(x), np.nan)
    s = shapely.line_locate_point(line, pts[near])
    out[near] = np.interp(s, s_vertex, targets)
    return out.reshape(np.shape(x))


def compose_sizing(ambient, x, y, *, grade, floor=None, bands=()):
    """The final size field and a report; see the module docstring."""
    h = _limit(np.asarray(ambient, float), x, y, grade)
    report = {}
    bands = [np.asarray(b, float) for b in bands]     # traversed twice (review r3 F3)
    if floor is not None:
        up = graded_up(floor, x, y, grade)
        report["raised_by_floor"] = int((up > h).sum())
        h = np.maximum(h, up)
    for k, band in enumerate(bands):
        on = np.isfinite(band)
        if not on.any():
            raise ValueError(f"band {k} covers no lattice point")
        lo = _limit(np.where(on, band, np.inf), x, y, grade)
        hi = graded_up(np.where(on, band, -np.inf), x, y, grade)
        h = np.minimum(np.maximum(h, hi), lo)
        report[f"band_{k}_cells"] = int(on.sum())
    # every band must come out at its own target: two bands closer than their
    # sizes allow under the gradation cannot both hold, and the later one
    # would silently win (review round 2 F17). A band whose own sizes change
    # faster along the line than the gradation allows is smoothed a little
    # (0.9 % on the Tokyo Bay interface, 2026-10-01); up to BAND_TOLERANCE
    # that is accepted and reported, beyond it refused.
    for k, band in enumerate(bands):
        on = np.isfinite(band)
        if floor is not None:
            # where a band leaves the final field below the time-step floor
            # (review F4; counted on the field that comes out, round 3 F1).
            # The caller decides what to do: meshes are made from the real
            # depths, and the time step is the depth stage's business (owner,
            # 2026-10-01), so 446 reports it
            report[f"band_{k}_below_floor_cells"] = int(
                (on & (h < np.asarray(floor) - 1e-6)).sum())
        dev = float(np.max(np.abs(h[on] - band[on]) / band[on]))
        report[f"band_{k}_max_rel_deviation"] = dev
        if dev > BAND_TOLERANCE:
            raise ValueError(f"band {k} cannot hold its sizes (off by up to {dev:.1%}): "
                             "its own sizes vary faster than the gradation allows, or "
                             "another band is too close")
    if floor is not None:
        report["below_floor_fraction"] = float(np.mean(h < np.asarray(floor) - 1e-6))
    return h, report


def _signed_areas(nodes, elements):
    a = nodes[elements]
    return 0.5 * ((a[:, 1, 0] - a[:, 0, 0]) * (a[:, 2, 1] - a[:, 0, 1])
                  - (a[:, 2, 0] - a[:, 0, 0]) * (a[:, 1, 1] - a[:, 0, 1]))


def merge_outer(base: Fort14Mesh, outer_nodes, outer_elements, interface_outer,
                interface_base, outer_open, *, tol_m: float = 0.5) -> Fort14Mesh:
    """Append an outer mesh to ``base`` through their shared interface nodes.

    ``interface_outer[i]`` (an outer node) and ``interface_base[i]`` (a base
    node) are the same point, at most ``tol_m`` apart; the base coordinate
    is kept.  The base's nodes and elements stay first and unchanged; outer
    elements are turned to the base's orientation.  ``outer_open`` is the new
    open boundary in outer indices; it becomes the only open boundary.
    Depths of the outer nodes are NaN until a depth stage fills them.

    The base must hold float64 coordinates and depths (what the readers
    give): the merged arrays are float64, and the base is carried bit for
    bit (review round 9 F6). Every index is checked to be a whole number in
    range before it is used (round 9 F4).
    """
    from fvcom_mesh_tools.io.fvcom_native import _indices

    if np.asarray(base.nodes).dtype != np.float64 or np.asarray(base.depths).dtype != np.float64:
        raise ValueError(f"the base must be float64 (nodes {np.asarray(base.nodes).dtype}, "
                         f"depths {np.asarray(base.depths).dtype})")
    outer_nodes = np.asarray(outer_nodes, float)
    if outer_nodes.ndim != 2 or outer_nodes.shape[1] < 2 or not np.isfinite(outer_nodes).all():
        raise ValueError("outer nodes must be finite (N, 2) coordinates")
    n_out = len(outer_nodes)
    outer_elements = _indices(outer_elements, n_out, "outer elements", ndim=2)
    if outer_elements.shape[1:] != (3,):
        raise ValueError(f"outer elements must be (NE, 3), not {outer_elements.shape}")
    io_ = _indices(interface_outer, n_out, "interface_outer")
    ib = _indices(interface_base, base.n_nodes, "interface_base")
    outer_open = _indices(outer_open, n_out, "outer_open")
    if io_.shape != ib.shape or len(io_) < 2:
        raise ValueError("the interface needs two nodes or more, paired one to one")
    gap = np.hypot(*(outer_nodes[io_] - base.nodes[ib, :2]).T)
    if gap.max() > tol_m:
        raise ValueError(f"interface nodes differ from the base by up to {gap.max():.3f} m")
    nb = base.n_nodes
    new = np.setdiff1d(np.arange(len(outer_nodes)), io_)
    index = np.full(len(outer_nodes), -1, np.int64)
    index[io_] = ib
    index[new] = nb + np.arange(len(new))
    elems = index[outer_elements]
    sign_base = np.sign(np.median(_signed_areas(base.nodes[:, :2], base.elements)))
    nodes = np.vstack([base.nodes[:, :2], outer_nodes[new]])
    flip = np.sign(_signed_areas(nodes, elems)) != sign_base
    elems[flip] = elems[flip][:, [0, 2, 1]]
    depths = np.r_[np.asarray(base.depths, float), np.full(len(new), np.nan)]
    merged = replace(base, nodes=nodes, depths=depths,
                     elements=np.vstack([base.elements, elems]),
                     open_boundaries=[index[np.asarray(outer_open, np.int64)]],
                     land_boundaries=[])
    return merged


def verify_frozen_base(merged: Fort14Mesh, base: Fort14Mesh, interface_base) -> dict:
    """Is the base inside ``merged`` exactly as it was?  Raises if not.

    Checks that the first ``NP`` nodes and ``NE`` elements are the base's,
    bit for bit, and that every base open-boundary edge is now interior
    (one base element and one outer element on it).
    """
    nb, eb = base.n_nodes, base.n_elements

    def same_bits(a, b):
        # dtype and bit pattern: numerical equality lets -0.0 pass for +0.0
        # (review round 8 F12)
        a, b = np.ascontiguousarray(a), np.ascontiguousarray(b)
        return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()

    if not same_bits(merged.nodes[:nb, :2], base.nodes[:, :2]):
        raise ValueError("base node coordinates changed")
    me, be = np.asarray(merged.elements[:eb]), np.asarray(base.elements)
    if not (me.dtype.kind in "iu" and be.dtype.kind in "iu" and np.array_equal(me, be)):
        raise ValueError("base elements changed")
    if not same_bits(merged.depths[:nb], base.depths):
        raise ValueError("base depths changed")
    ib = np.asarray(interface_base, np.int64)
    e = np.sort(np.vstack([merged.elements[:, [0, 1]], merged.elements[:, [1, 2]],
                           merged.elements[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(merged.n_elements), 3)
    keys = e[:, 0] * merged.n_nodes + e[:, 1]
    xy = merged.nodes[:, :2]
    for a, b in zip(ib[:-1], ib[1:]):
        k = min(a, b) * merged.n_nodes + max(a, b)
        who = owner[keys == k]
        if len(who) != 2 or (who < eb).sum() != 1:
            raise ValueError(f"interface edge {a}-{b} is not shared by one base and one "
                             f"outer element (elements {who.tolist()})")
        # the two elements must lie on opposite sides of the edge, or the outer
        # one overlaps the base (review F10)
        side = []
        for e_ in who:
            c = [v for v in merged.elements[e_] if v not in (a, b)][0]
            side.append(np.sign((xy[b, 0] - xy[a, 0]) * (xy[c, 1] - xy[a, 1])
                                - (xy[b, 1] - xy[a, 1]) * (xy[c, 0] - xy[a, 0])))
        if side[0] * side[1] >= 0:
            raise ValueError(f"interface edge {a}-{b}: the base and the outer element are "
                             "on the same side (they overlap)")
    return {"n_base_nodes": nb, "n_base_elements": eb,
            "n_interface_edges": int(len(ib) - 1),
            "n_nodes": merged.n_nodes, "n_elements": merged.n_elements}


def check_no_overlap(merged: Fort14Mesh, n_base_elements: int, rel_tol: float = 1e-9) -> dict:
    """No new element may cover any part of the base (review round 2 F9).

    The shared interface is a line, so an outer element may touch the base
    footprint only along it: an intersection with positive area (beyond
    ``rel_tol`` of the element's area, for round-off) is an overlap, wherever
    it is -- not only at the interface edges ``verify_frozen_base`` checks.
    Raises on the first overlap; returns counts.
    """
    import shapely

    xy = merged.nodes[:, :2]
    base = shapely.union_all(shapely.polygons(xy[merged.elements[:n_base_elements]]))
    outer = shapely.polygons(xy[merged.elements[n_base_elements:]])
    tree = shapely.STRtree(outer)
    cand = np.unique(tree.query(base, predicate="intersects"))
    if len(cand):
        area = shapely.area(shapely.intersection(outer[cand], base))
        own = shapely.area(outer[cand])
        bad = cand[area > rel_tol * own]
        if len(bad):
            k = int(bad[0])
            raise ValueError(f"{len(bad)} new element(s) overlap the base, e.g. element "
                             f"{n_base_elements + k} ({float(area[cand == k][0]):.3g} m2)")
    return {"n_outer_touching_base": int(len(cand))}


def land_segments(elements, open_chains) -> list[tuple[int, np.ndarray]]:
    """The land boundary runs: every boundary loop minus its open-boundary edges.

    Returns ``(ibtype, nodes)`` runs as ``Fort14Mesh.land_boundaries`` wants
    them: 0 for a run on a loop that carries an open boundary (the mainland),
    1 for a closed loop without one (an island).  An open chain's end nodes
    appear in the land runs too, as in ADCIRC.
    """
    from fvcom_mesh_tools.io.fvcom_native import boundary_loops

    open_edges = set()
    for c in open_chains:
        c = np.asarray(c, np.int64)
        open_edges |= {frozenset((int(a), int(b))) for a, b in zip(c[:-1], c[1:])}
    out = []
    for loop in boundary_loops(np.asarray(elements)):
        loop = [int(v) for v in loop]
        n = len(loop)
        is_open = [frozenset((loop[i], loop[(i + 1) % n])) in open_edges for i in range(n)]
        if not any(is_open):
            out.append((1, np.array(loop + [loop[0]], np.int64)))
            continue
        start = is_open.index(True)
        run: list[int] = []
        for k in range(1, n + 1):
            i = (start + k) % n
            if is_open[i]:
                if run:      # one land edge is a run too (review F11)
                    out.append((0, np.array(run + [loop[i]], np.int64)))
                run = []
            else:
                run.append(loop[i])
        if run:
            out.append((0, np.array(run + [loop[(start + n) % n]], np.int64)))
    return out


def rfactor_smooth_free(h0, ei, ej, free, *, rmax, hmin, hmax=None, max_iter=5000):
    """r-factor limiter that moves only the ``free`` nodes.

    As ``dem.m7001.rfactor_smooth`` (Beckmann-Haidvogel), but a node that is
    not free keeps its depth: on an edge with one fixed end the free end
    takes the whole correction, and an edge with two fixed ends is left as
    it is (its r is the base's own).  ``hmin`` and ``hmax`` bound the free
    depths *during* the smoothing, so the result satisfies both (a cap
    applied afterwards could break the r-factor; review F3). Returns
    ``(depth, iterations, max r over edges with a free end)``; raises when the
    limit is not reached -- infeasible (e.g. a fixed 1 m node beside a free
    node held at 3 m) or not converged within ``max_iter``.
    """
    h = np.asarray(h0, float).copy()
    free = np.asarray(free, bool)
    ei, ej = np.asarray(ei), np.asarray(ej)
    # finite positive depths and sane controls, or NaN slips through the
    # r > limit test (review round 2 F15)
    if not (0 < rmax < 1 and np.isfinite(hmin) and hmin > 0
            and (hmax is None or (np.isfinite(hmax) and hmax >= hmin))
            and int(max_iter) >= 1):
        raise ValueError(f"bad controls: rmax {rmax}, hmin {hmin}, hmax {hmax}, "
                         f"max_iter {max_iter}")
    live = free[ei] | free[ej]
    ei, ej = ei[live], ej[live]
    fi, fj = free[ei], free[ej]
    used = np.unique(np.r_[ei, ej, np.flatnonzero(free)])
    if not (np.isfinite(h[used]).all() and (h[used] > 0).all()):
        raise ValueError("depths on the limited edges and free nodes must be finite and "
                         "positive")
    # the bounds hold for every free node, on an edge or not (review round 9
    # F10)
    if hmax is not None:
        h = np.where(free, np.minimum(h, hmax), h)
    h = np.where(free, np.maximum(h, hmin), h)
    if not len(ei):
        return h, 0, 0.0
    for it in range(int(max_iter)):
        hi, hj = h[ei], h[ej]
        r = np.abs(hi - hj) / (hi + hj)
        bad = r > rmax + 1e-9
        if not bad.any():
            return h, it, float(r.max())
        excess = np.where(bad, np.abs(hi - hj) - rmax * (hi + hj), 0.0)
        sgn = np.sign(hi - hj)                       # +1: i is the deeper end
        both = fi & fj
        # a free end moves by half the excess when the other end moves too,
        # by all of it (scaled for the r-denominator) when the other is fixed
        di = np.where(both, excess / 2, excess / (1 + rmax))
        dj = np.where(both, excess / 2, excess / (1 + rmax))
        add = np.zeros_like(h)
        cnt = np.zeros_like(h)
        np.add.at(add, ei[fi & bad], -(sgn * di)[fi & bad])
        np.add.at(cnt, ei[fi & bad], 1.0)
        np.add.at(add, ej[fj & bad], (sgn * dj)[fj & bad])
        np.add.at(cnt, ej[fj & bad], 1.0)
        step = add / np.where(cnt > 0, cnt, 1.0)
        h = np.where(free, np.clip(h + step, hmin, np.inf if hmax is None else hmax), h)
    hi, hj = h[ei], h[ej]
    r = np.abs(hi - hj) / (hi + hj)
    if np.isfinite(r).all() and r.max() <= rmax + 1e-9:
        return h, int(max_iter), float(r.max())     # met on the last pass (round 2 F14)
    k = int(np.argmax(r))
    raise ValueError(f"r-factor limit {rmax} not reached in {max_iter} iterations: "
                     f"r = {r[k]:.4f} on edge {int(ei[k])}-{int(ej[k])} "
                     f"(depths {hi[k]:.3f}, {hj[k]:.3f} m); the depth bounds may make "
                     "it infeasible")


def trim_lone_corners(elements, mutable, keep_nodes=(), max_rounds=20):
    """Drop the element under a node that no other element shares.

    A cape one element wide ends in a node that sits in a single element;
    FVCOM never updates such a node. Bisecting the element (as
    ``walls.open_lone_corners`` does for wall bends) leaves two slivers at a
    cape tip, so here the tip is not resolved: the element goes, as the
    resolution principle says for what the element size cannot carry. Only
    ``mutable`` elements are dropped, never one holding a node in
    ``keep_nodes`` (the open boundary). Repeats while new lone nodes appear.

    Returns ``(elements, mutable, report)``.
    """
    t = np.asarray(elements, np.int64)
    mut = np.asarray(mutable, bool)
    keep = set(int(v) for v in keep_nodes)
    dropped, left = 0, []
    for _ in range(max_rounds):
        if len(t) == 0:
            break
        count = np.bincount(t.ravel(), minlength=int(t.max()) + 1)
        edge_count: dict[tuple[int, int], int] = {}
        for a, b in np.sort(np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]]), axis=1).tolist():
            edge_count[(a, b)] = edge_count.get((a, b), 0) + 1
        drop = np.zeros(len(t), bool)
        left = []
        for v in np.flatnonzero(count == 1).tolist():
            k = int(np.flatnonzero((t == v).any(axis=1))[0])
            a, b = sorted(int(x) for x in t[k] if x != v)
            # only a spike: the side facing the lone node is shared, so the
            # element's removal leaves no new lone node behind
            # Removals are judged against what survives: each one takes its
            # edges out of the count, so two elements that each lean on the
            # other's shared side are not both dropped (review round 9 F9).
            if (not drop[k] and mut[k] and not (set(t[k].tolist()) & keep)
                    and edge_count[(a, b)] == 2):
                drop[k] = True
                for e in ((t[k, 0], t[k, 1]), (t[k, 1], t[k, 2]), (t[k, 2], t[k, 0])):
                    edge_count[tuple(sorted(int(x) for x in e))] -= 1
            elif not drop[k]:
                left.append(v)
        if not drop.any():
            break
        t, mut = t[~drop], mut[~drop]
        dropped += int(drop.sum())
    return t, mut, {"n_elements_dropped": dropped, "lone_nodes_left": left}
