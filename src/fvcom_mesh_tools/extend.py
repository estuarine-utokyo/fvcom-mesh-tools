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
    "land_segments",
    "compose_sizing",
    "graded_up",
    "merge_outer",
    "rfactor_smooth_free",
    "verify_frozen_base",
]


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
    if floor is not None:
        up = graded_up(floor, x, y, grade)
        report["raised_by_floor"] = int((up > h).sum())
        h = np.maximum(h, up)
    for k, band in enumerate(bands):
        band = np.asarray(band, float)
        on = np.isfinite(band)
        if not on.any():
            raise ValueError(f"band {k} covers no lattice point")
        lo = _limit(np.where(on, band, np.inf), x, y, grade)
        hi = graded_up(np.where(on, band, -np.inf), x, y, grade)
        h = np.minimum(np.maximum(h, hi), lo)
        report[f"band_{k}_cells"] = int(on.sum())
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
    """
    outer_nodes = np.asarray(outer_nodes, float)
    outer_elements = np.asarray(outer_elements, np.int64)
    io_ = np.asarray(interface_outer, np.int64)
    ib = np.asarray(interface_base, np.int64)
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
    if not np.array_equal(merged.nodes[:nb, :2], base.nodes[:, :2]):
        raise ValueError("base node coordinates changed")
    if not np.array_equal(merged.elements[:eb], base.elements):
        raise ValueError("base elements changed")
    if not np.array_equal(merged.depths[:nb], base.depths):
        raise ValueError("base depths changed")
    ib = np.asarray(interface_base, np.int64)
    e = np.sort(np.vstack([merged.elements[:, [0, 1]], merged.elements[:, [1, 2]],
                           merged.elements[:, [2, 0]]]), axis=1)
    owner = np.tile(np.arange(merged.n_elements), 3)
    keys = e[:, 0] * merged.n_nodes + e[:, 1]
    for a, b in zip(ib[:-1], ib[1:]):
        k = min(a, b) * merged.n_nodes + max(a, b)
        who = owner[keys == k]
        if len(who) != 2 or (who < eb).sum() != 1:
            raise ValueError(f"interface edge {a}-{b} is not shared by one base and one "
                             f"outer element (elements {who.tolist()})")
    return {"n_base_nodes": nb, "n_base_elements": eb,
            "n_interface_edges": int(len(ib) - 1),
            "n_nodes": merged.n_nodes, "n_elements": merged.n_elements}


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
                if len(run) > 1:
                    out.append((0, np.array(run + [loop[i]], np.int64)))
                run = []
            else:
                run.append(loop[i])
        if len(run) > 1:
            out.append((0, np.array(run + [loop[(start + n) % n]], np.int64)))
    return out


def rfactor_smooth_free(h0, ei, ej, free, *, rmax, hmin, max_iter=5000):
    """r-factor limiter that moves only the ``free`` nodes.

    As ``dem.m7001.rfactor_smooth`` (Beckmann-Haidvogel), but a node that is
    not free keeps its depth: on an edge with one fixed end the free end
    takes the whole correction, and an edge with two fixed ends is left as
    it is (its r is the base's own).  Returns ``(depth, iterations, max r
    over edges with a free end)``.
    """
    h = np.asarray(h0, float).copy()
    free = np.asarray(free, bool)
    ei, ej = np.asarray(ei), np.asarray(ej)
    live = free[ei] | free[ej]
    ei, ej = ei[live], ej[live]
    fi, fj = free[ei], free[ej]
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
        h = np.where(free, np.maximum(h + step, hmin), h)
    hi, hj = h[ei], h[ej]
    return h, int(max_iter), float((np.abs(hi - hj) / (hi + hj)).max())
