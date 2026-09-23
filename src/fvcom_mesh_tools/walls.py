"""Walls: turn a line of mesh edges into domain boundary.

An interior edge does not stop water in FVCOM.  Scalars and elevation live on
nodes in median-dual control volumes (``tge.F``: each face runs from an
element centroid to an edge midpoint), and the control volume of a node ON a
line contains parts of the elements on both sides of it -- water crossing the
line there never crosses a face.  Momentum lives on elements, and the shared
edge is the face between them.  So a breakwater laid as a line of edges is
not a breakwater.

What stops water is boundary.  :func:`split_along_walls` duplicates the nodes
along each wall so that the elements on one side use one copy and those on the
other side use the other; the two sides then share no edge, and FVCOM finds
the wall topologically -- ``NBE = 0`` on both sides, ``ISBCE = 1`` -- exactly
as it finds a coastline.  That is also how FVCOM's own ``THIN_DAM`` represents
a wall; what ``THIN_DAM`` adds is the re-stitching above a crest, i.e.
overtopping, which a breakwater treated as impermeable does not need
(owner, 2026-09-23).  A split mesh without ``THIN_DAM`` needs no FVCOM change.

Which nodes get copies follows from the fan of elements around each node,
cut by the wall edges at that node into sectors:

* a node inside a wall: its closed fan is cut twice -> two sectors, one copy;
* a node where a wall meets the coast: its open fan is cut once -> two
  sectors, one copy, and the coastline continues on both sides;
* a free tip: its closed fan is cut once -> still ONE sector, no copy -- the
  boundary turns through 360 degrees there;
* a junction of k walls: k sectors, k - 1 copies.

The copies share coordinates exactly, which is the point and also the thing a
duplicate-node check will object to, so every copy is returned with the node
it was copied from.
"""

from __future__ import annotations

from typing import Any

import numpy as np

__all__ = ["split_along_walls", "wall_edges_from_path"]


def _edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def wall_edges_from_path(path) -> np.ndarray:
    """Consecutive node pairs of a polyline given as mesh node ids."""
    p = np.asarray(path, dtype=np.int64).ravel()
    if p.size < 2:
        return np.zeros((0, 2), dtype=np.int64)
    return np.column_stack([p[:-1], p[1:]])


def split_along_walls(nodes, elements, wall_edges) -> tuple[np.ndarray, np.ndarray,
                                                             np.ndarray, dict[str, Any]]:
    """Cut the mesh along ``wall_edges`` so they become boundary on both sides.

    ``wall_edges`` are ``(M, 2)`` node-id pairs, each of which must be an
    INTERIOR edge of the mesh -- shared by exactly two elements.  An edge on
    the boundary already is refused: it has only one side to separate.

    Returns ``(nodes, elements, copy_of, report)``.  ``nodes`` has the copies
    appended; ``copy_of[k]`` is the original node of row ``k`` (itself for an
    original), so depths and any other per-node field are carried across as
    ``field[copy_of]``.  Element ORDER is unchanged -- only the ids they
    reference move -- so element-keyed data needs no mapping at all.
    """
    xy = np.asarray(nodes, dtype=float)
    orig = np.array(elements, dtype=np.int64, copy=True)
    # Sectors are found on the ORIGINAL connectivity and the relabelling is
    # written to a copy.  Doing both on one array made each node see the ids
    # its neighbours had already been given -- edges the table did not have.
    tri = orig.copy()
    n0 = len(xy)
    walls = np.asarray(wall_edges, dtype=np.int64).reshape(-1, 2)
    if walls.size and (walls.min() < 0 or walls.max() >= n0):
        raise ValueError("a wall edge references a node the mesh does not have")
    if (walls[:, 0] == walls[:, 1]).any():
        raise ValueError("a wall edge joins a node to itself")

    # every edge and the elements that use it
    edge_elems: dict[tuple[int, int], list[int]] = {}
    for k, (a, b, c) in enumerate(orig.tolist()):
        for u, v in ((a, b), (b, c), (c, a)):
            edge_elems.setdefault(_edge_key(u, v), []).append(k)
    wall_set = {_edge_key(int(a), int(b)) for a, b in walls.tolist()}
    for e in wall_set:
        users = edge_elems.get(e)
        if users is None:
            raise ValueError(f"wall edge {e} is not an edge of the mesh; walls "
                             "must be meshed as constrained edges first")
        if len(users) != 2:
            raise ValueError(f"wall edge {e} is used by {len(users)} element(s); "
                             "only an interior edge has two sides to separate")

    touched = sorted({n for e in wall_set for n in e})
    fan: dict[int, list[int]] = {}
    for k, row in enumerate(orig.tolist()):
        for n in row:
            fan.setdefault(n, []).append(k)

    extra_xy: list[np.ndarray] = []
    copy_src: list[int] = []
    sectors_at: dict[int, int] = {}
    for v in touched:
        elems = fan.get(v, [])
        # elements around v are joined across every edge (v, w) that is NOT
        # a wall; the connected pieces are the sectors
        parent = {k: k for k in elems}

        def find(k):
            while parent[k] != k:
                parent[k] = parent[parent[k]]
                k = parent[k]
            return k

        for k in elems:
            for w in orig[k].tolist():
                if w == v:
                    continue
                e = _edge_key(v, w)
                if e in wall_set:
                    continue
                for j in edge_elems[e]:
                    if j != k and j in parent:
                        ra, rb = find(k), find(j)
                        if ra != rb:
                            parent[ra] = rb
        groups: dict[int, list[int]] = {}
        for k in elems:
            groups.setdefault(find(k), []).append(k)
        sectors = sorted(groups.values(), key=min)
        sectors_at[v] = len(sectors)
        # the first sector keeps the original id; each further one gets a copy
        for sector in sectors[1:]:
            new_id = n0 + len(copy_src)
            copy_src.append(v)
            extra_xy.append(xy[v])
            for k in sector:
                tri[k][orig[k] == v] = new_id

    copy_of = np.concatenate([np.arange(n0), np.asarray(copy_src, dtype=np.int64)])
    out_xy = np.vstack([xy, np.asarray(extra_xy).reshape(-1, xy.shape[1])]) \
        if extra_xy else xy.copy()

    # every wall edge must now be boundary on BOTH sides: two distinct
    # boundary edges whose ends are copies of the same pair
    e2 = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    u, c = np.unique(e2, axis=0, return_counts=True)
    bnd = u[c == 1]
    orig_pairs = np.sort(copy_of[bnd], axis=1)
    per_wall = {e: 0 for e in wall_set}
    for a, b in orig_pairs.tolist():
        if (a, b) in per_wall:
            per_wall[(a, b)] += 1
    not_split = [e for e, n in per_wall.items() if n != 2]
    if not_split:
        raise ValueError(f"{len(not_split)} wall edge(s) did not become boundary on "
                         f"both sides, first {not_split[:3]}")
    if (c > 2).any():
        raise ValueError("the split produced a non-manifold edge")

    # How many pieces the mesh is now in.  A wall from coast to coast cuts
    # the domain in two, which may be what was meant (a closed dock) and may
    # not (a breakwater whose OSM gap was lost); either way it is said.
    ek: dict[tuple[int, int], list[int]] = {}
    for k, (a, b, cc) in enumerate(tri.tolist()):
        for x, y in ((a, b), (b, cc), (cc, a)):
            ek.setdefault(_edge_key(x, y), []).append(k)
    root = list(range(len(tri)))

    def _find(k):
        while root[k] != k:
            root[k] = root[root[k]]
            k = root[k]
        return k

    for users in ek.values():
        if len(users) == 2:
            ra, rb = _find(users[0]), _find(users[1])
            if ra != rb:
                root[ra] = rb
    n_components = len({_find(k) for k in range(len(tri))})

    tips = [v for v, s in sectors_at.items() if s == 1]
    return out_xy, tri, copy_of, {
        "n_wall_edges": len(wall_set),
        "n_wall_nodes": len(touched),
        "n_copies": len(copy_src),
        "n_free_tips": len(tips),
        "free_tips": tips,
        "max_sectors_at_a_node": max(sectors_at.values()) if sectors_at else 0,
        "n_components": n_components,
        "pairs": [[int(s), int(n0 + k)] for k, s in enumerate(copy_src)],
    }
