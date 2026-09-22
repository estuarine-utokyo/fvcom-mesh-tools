"""Renumber a mesh's nodes and elements for memory locality.

**Why this exists, measured rather than assumed.** FVCOM decomposes the
domain itself at run time -- ``setup_domain.F`` calls ``DOMDEC``, which calls
METIS with ``ncommon = 2`` -- so nothing this package writes changes how the
mesh is partitioned. What it does change is the order *inside* each rank:
``genmap.F`` builds every rank's local numbering by walking the GLOBAL ids
from 1 to NGL and appending the ones that belong to it::

      DO I=1,NGL
         IF(EL_PID(I) == MYID) THEN
            N = N + 1
            NTEMP(N) = I

so each rank's arrays are the file's order, filtered. The file's numbering is
therefore the memory layout, and a mesh whose neighbours are far apart in
index space is a mesh whose neighbours are far apart in cache.

**And this project makes that worse.** A local refinement keeps the retained
elements first, in base order, and appends the patch -- which is what makes
the frozen zone verifiable, and is worth keeping. But every edge between a
new node and a retained one then spans thousands of indices. Measured on the
Banzu case: the goto2023 base has a bandwidth of 80 over 3,210 nodes, and the
refined mesh has **3,596** over 5,409. Renumbering is not a polish step here;
it puts back what the patch took out.

**The base is already well numbered**, which is the other half of the point:
RCM on goto2023 makes its bandwidth *worse*, 80 -> 128. So every function
here reports what it would do, and :func:`renumber_mesh` refuses to make a
mesh worse unless told to.
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = [
    "METHODS",
    "element_permutation",
    "hilbert_permutation",
    "locality_stats",
    "morton_permutation",
    "open_boundary_seeds",
    "rcm_permutation",
    "renumber_mesh",
]

METHODS = ("rcm", "morton", "hilbert")


# --------------------------------------------------------------------------
# what we are trying to improve
# --------------------------------------------------------------------------


def locality_stats(n_nodes: int, elements, xy=None) -> dict[str, Any]:
    """How far apart, in index space, the mesh's neighbours are.

    ``bandwidth`` is the classical one -- the largest ``|i - j|`` over the
    edges -- and it is the number a single bad pair ruins, so the mean and
    the 99th percentile are reported beside it: they are what a cache
    actually experiences. ``element_span`` is the same question asked of an
    element's three nodes, which is what FVCOM gathers on every cell.
    """
    tri = np.asarray(elements, dtype=np.int64)
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    e = np.unique(np.sort(e, axis=1), axis=0)
    d = np.abs(e[:, 0] - e[:, 1])
    span = tri.max(axis=1) - tri.min(axis=1)
    out = {
        "n_nodes": int(n_nodes),
        "n_elements": int(tri.shape[0]),
        "bandwidth": int(d.max()) if d.size else 0,
        "mean_edge_span": float(d.mean()) if d.size else 0.0,
        "p99_edge_span": float(np.percentile(d, 99)) if d.size else 0.0,
        # The profile of the adjacency matrix: summed over rows, the
        # distance to that row's leftmost entry -- what a skyline solver
        # stores, and a proxy for how much of each row's neighbourhood sits
        # in one cache line. Filled in below, where the edges exist.
        "profile": 0,
        "mean_element_span": float(span.mean()) if span.size else 0.0,
        "max_element_span": int(span.max()) if span.size else 0,
    }
    if d.size:
        lo = np.arange(n_nodes, dtype=np.int64)
        np.minimum.at(lo, e[:, 1], e[:, 0])
        np.minimum.at(lo, e[:, 0], e[:, 1])
        out["profile"] = int((np.arange(n_nodes) - lo).sum())
    return out


# --------------------------------------------------------------------------
# the orderings
# --------------------------------------------------------------------------


def _adjacency(n_nodes: int, elements):
    """CSR-ish node adjacency: ``(neighbours, start, stop)`` per node."""
    tri = np.asarray(elements, dtype=np.int64)
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    e = np.unique(np.vstack([e, e[:, ::-1]]), axis=0)
    e = e[np.argsort(e[:, 0], kind="stable")]
    ids = np.arange(n_nodes)
    return (e[:, 1], np.searchsorted(e[:, 0], ids),
            np.searchsorted(e[:, 0], ids, side="right"))


def rcm_permutation(n_nodes: int, elements, seeds=None) -> np.ndarray:
    """Reverse Cuthill-McKee over the node adjacency graph.

    ``seeds`` is where the sweep starts. This is not a detail: it is the
    difference between this and what SMS does, and it is worth more than the
    algorithm. SMS's workflow -- select the open-boundary nodestring, then
    renumber -- starts the sweep at the open boundary, and that is how the
    production base was numbered: measured on goto2023, seeding at the OBC
    gives a bandwidth of **78** against the mesh's own **80**, while letting
    scipy pick its own pseudo-peripheral node gives **128**. On the refined
    Banzu mesh the two agree (194 against 195), so the seeded sweep costs
    nothing where it does not help.

    Seeding at the WHOLE boundary -- coastline as well as open boundary --
    is a different thing and a bad one: a front 781 nodes wide produces a
    bandwidth of 781, ten times the mesh's own, whatever order the front is
    in. A selection marks where to start, not what to number first.

    ``seeds=None`` hands the choice to scipy. Returns ``perm`` such that
    ``perm[k]`` is the OLD id of the node that becomes new id ``k``.

    The reversal is the "R" of RCM and leaves every metric here unchanged
    -- ``|i - j|`` does not care which end it is measured from -- but it
    decides which end the seed lands at. Reversed, an OBC-seeded sweep puts
    the open boundary LAST, which is where the production base has it.
    """
    from collections import deque

    tri = np.asarray(elements, dtype=np.int64)
    if seeds is None:
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import reverse_cuthill_mckee

        e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
        e = np.vstack([e, e[:, ::-1]])
        a = coo_matrix((np.ones(e.shape[0], dtype=np.int8), (e[:, 0], e[:, 1])),
                       shape=(n_nodes, n_nodes)).tocsr()
        return np.asarray(reverse_cuthill_mckee(a, symmetric_mode=True),
                          dtype=np.int64)

    nbr, lo, hi = _adjacency(n_nodes, tri)
    deg = hi - lo
    seen = np.zeros(n_nodes, dtype=bool)
    order: list[int] = []
    seeds = [int(v) for v in np.unique(np.asarray(seeds, dtype=np.int64))]
    if any(v < 0 or v >= n_nodes for v in seeds):
        raise ValueError("a seed is not a node of this mesh")
    q = deque(sorted(seeds, key=lambda v: deg[v]))
    for v in q:
        seen[v] = True
    while True:
        while q:
            v = q.popleft()
            order.append(v)
            nb = [int(w) for w in nbr[lo[v]:hi[v]] if not seen[w]]
            for w in sorted(nb, key=lambda w: deg[w]):
                seen[w] = True
                q.append(w)
        rest = np.flatnonzero(~seen)
        if not rest.size:
            break
        # A mesh in several pieces, or a seed set that cannot reach all of
        # it: carry on from the lowest-degree node left rather than stop.
        nxt = int(rest[np.argmin(deg[rest])])
        seen[nxt] = True
        q.append(nxt)
    return np.asarray(order, dtype=np.int64)[::-1].copy()


def _boundary_nodes(elements) -> np.ndarray:
    """Every node on a boundary edge -- coastline and open boundary alike."""
    tri = np.asarray(elements, dtype=np.int64)
    e = np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    u, c = np.unique(np.sort(e, axis=1), axis=0, return_counts=True)
    return np.unique(u[c == 1])


def _quantise(xy, bits: int) -> np.ndarray:
    """Coordinates onto a ``2**bits`` integer lattice, per axis."""
    p = np.asarray(xy, dtype=float)[:, :2]
    lo = p.min(axis=0)
    span = np.maximum(p.max(axis=0) - lo, 1e-30)
    q = np.floor((p - lo) / span * (2 ** bits - 1e-9)).astype(np.int64)
    return np.clip(q, 0, 2 ** bits - 1)


def _interleave(q: np.ndarray, bits: int) -> np.ndarray:
    """Morton (Z-order) key from a 2-column integer lattice."""
    key = np.zeros(q.shape[0], dtype=np.int64)
    for b in range(bits):
        key |= ((q[:, 0] >> b) & 1) << (2 * b)
        key |= ((q[:, 1] >> b) & 1) << (2 * b + 1)
    return key


def morton_permutation(xy, bits: int = 16) -> np.ndarray:
    """Sort by Morton (Z) order of the coordinates.

    Cheap, and it has one well-known flaw: the Z curve jumps the full width
    of the domain every time it crosses a power-of-two boundary, so a few
    pairs of neighbours land far apart however fine the lattice is. The
    textbook remedy is :func:`hilbert_permutation` -- which on the meshes
    here does not help: on a 33x33 grid Morton's mean edge span is 22.4 and
    Hilbert's 25.3. Both are offered because both are asked for; neither
    beat RCM on anything measured.
    """
    return np.argsort(_interleave(_quantise(xy, bits), bits), kind="stable")


def _hilbert_key(q: np.ndarray, bits: int) -> np.ndarray:
    """Hilbert index of a 2-column integer lattice (Skilling's transform)."""
    x = q[:, 0].copy()
    y = q[:, 1].copy()
    rx = np.zeros_like(x)
    ry = np.zeros_like(x)
    d = np.zeros_like(x)
    s = np.int64(2 ** (bits - 1))
    while s > 0:
        rx = ((x & s) > 0).astype(np.int64)
        ry = ((y & s) > 0).astype(np.int64)
        d += s * s * ((3 * rx) ^ ry)
        # rotate the quadrant so the curve stays continuous across it
        swap = ry == 0
        flip = swap & (rx == 1)
        x[flip] = s - 1 - x[flip]
        y[flip] = s - 1 - y[flip]
        xs = x[swap].copy()
        x[swap] = y[swap]
        y[swap] = xs
        s //= 2
    return d


def hilbert_permutation(xy, bits: int = 16) -> np.ndarray:
    """Sort by Hilbert-curve order of the coordinates.

    The Hilbert curve visits every cell of the lattice with unit steps, so
    two points close on the curve are close in space and -- unlike Morton --
    the converse is meant to hold far more often. Measured on structured
    grids of 17, 33 and 65 a side it did not beat Morton's mean edge span,
    and neither came close to RCM's. Kept because it is the other thing
    asked for, and because a curve that loses on this mesh may win on the
    next one; the command measures rather than assumes.
    """
    return np.argsort(_hilbert_key(_quantise(xy, bits), bits), kind="stable")


def element_permutation(elements, node_perm=None, xy=None,
                        method: str = "node") -> np.ndarray:
    """An element order to match a node order.

    Elements are renumbered too because ``EL_PID`` is indexed by the global
    element id and ``genmap.F`` walks it in the same way: an element order
    that disagrees with the node order costs the gather on every cell.
    ``method="node"`` sorts by the element's smallest NEW node id, which
    follows whatever the nodes did; ``"centroid"`` sorts by the Hilbert
    order of the centroid instead.
    """
    tri = np.asarray(elements, dtype=np.int64)
    if method == "centroid":
        if xy is None:
            raise ValueError("centroid ordering needs coordinates")
        cent = np.asarray(xy, dtype=float)[:, :2][tri].mean(axis=1)
        return hilbert_permutation(cent)
    if method != "node":
        raise ValueError(f"unknown element ordering {method!r}")
    if node_perm is None:
        key = tri.min(axis=1)
    else:
        inv = np.empty(node_perm.shape[0], dtype=np.int64)
        inv[np.asarray(node_perm, dtype=np.int64)] = np.arange(
            node_perm.shape[0])
        key = inv[tri].min(axis=1)
    return np.argsort(key, kind="stable")


# --------------------------------------------------------------------------
# applying it
# --------------------------------------------------------------------------


def open_boundary_seeds(mesh) -> np.ndarray | None:
    """The mesh's open-boundary nodes, or None when it has none."""
    segs = [np.asarray(s, dtype=np.int64) for s in mesh.open_boundaries]
    segs = [s for s in segs if s.size]
    return np.unique(np.concatenate(segs)) if segs else None


def renumber_mesh(mesh, method: str = "rcm", *, elements: str = "node",
                  seeds: Any = "obc", only_if_better: bool = True):
    """Return ``(mesh, report)`` with nodes and elements renumbered.

    Everything indexed by a node or an element moves with it: the
    coordinates, the depths, the connectivity, the open-boundary lists and
    the land-boundary lists. ``report["node_perm"]`` is the OLD id of each
    new node, which is what a caller needs to carry anything else -- a
    sponge file, a node map, an observation index -- across.

    ``only_if_better`` keeps the mesh as it is when the ordering would not
    reduce the mean edge span. The production base is already well numbered
    and RCM takes its bandwidth from 80 to 128; a tool that always "improves"
    is a tool that cannot be trusted to have measured anything.
    """
    import dataclasses

    xy = np.asarray(mesh.nodes, dtype=float)
    tri = np.asarray(mesh.elements, dtype=np.int64)
    before = locality_stats(xy.shape[0], tri)
    if method == "rcm":
        if isinstance(seeds, str):
            if seeds == "obc":
                seed_ids = open_boundary_seeds(mesh)
            elif seeds == "auto":
                seed_ids = None
            elif seeds == "boundary":
                seed_ids = _boundary_nodes(tri)
            else:
                raise ValueError(f"unknown seed set {seeds!r}")
        else:
            seed_ids = None if seeds is None else np.asarray(seeds, dtype=np.int64)
        node_perm = rcm_permutation(xy.shape[0], tri, seed_ids)
    elif method == "morton":
        node_perm = morton_permutation(xy)
    elif method == "hilbert":
        node_perm = hilbert_permutation(xy)
    else:
        raise ValueError(f"unknown method {method!r}; expected one of {METHODS}")
    inv = np.empty(node_perm.shape[0], dtype=np.int64)
    inv[node_perm] = np.arange(node_perm.shape[0], dtype=np.int64)
    new_tri = inv[tri]
    elem_perm = element_permutation(tri, node_perm, xy, elements)
    new_tri = new_tri[elem_perm]
    after = locality_stats(xy.shape[0], new_tri)

    report = {
        "method": method,
        "element_order": elements,
        "seeds": (seeds if isinstance(seeds, str) else "explicit")
        if method == "rcm" else None,
        "before": before,
        "after": after,
        "applied": True,
        "bandwidth_ratio": (after["bandwidth"] / before["bandwidth"]
                            if before["bandwidth"] else 1.0),
        "mean_edge_span_ratio": (after["mean_edge_span"] / before["mean_edge_span"]
                                 if before["mean_edge_span"] else 1.0),
        "p99_edge_span_ratio": (after["p99_edge_span"] / before["p99_edge_span"]
                                if before["p99_edge_span"] else 1.0),
    }
    # BOTH the mean and the tail have to improve. The mean alone is not
    # enough: Morton on the goto2023 base takes the mean from 32.6 to 27.8
    # and the 99th percentile from 74 to 610, because the Z curve jumps the
    # width of the domain at every power-of-two boundary. A mesh whose
    # typical neighbour is closer and whose worst neighbours are eight times
    # further away has not been improved.
    if only_if_better and (report["mean_edge_span_ratio"] >= 1.0
                           or report["p99_edge_span_ratio"] >= 1.0):
        report["applied"] = False
        report["note"] = (
            f"{method} would not improve this mesh (mean edge span "
            f"{before['mean_edge_span']:.1f} -> {after['mean_edge_span']:.1f}, "
            f"p99 {before['p99_edge_span']:.1f} -> {after['p99_edge_span']:.1f})"
            "; it is already well numbered")
        report["node_perm"] = np.arange(xy.shape[0], dtype=np.int64)
        report["element_perm"] = np.arange(tri.shape[0], dtype=np.int64)
        return mesh, report

    out = dataclasses.replace(
        mesh,
        nodes=xy[node_perm].copy(),
        depths=np.asarray(mesh.depths, dtype=float)[node_perm].copy(),
        elements=new_tri,
        open_boundaries=[inv[np.asarray(s, dtype=np.int64)]
                         for s in mesh.open_boundaries],
        land_boundaries=[(int(ib), inv[np.asarray(s, dtype=np.int64)])
                         for ib, s in mesh.land_boundaries],
    )
    report["node_perm"] = node_perm
    report["element_perm"] = elem_perm
    return out, report
