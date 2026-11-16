"""The seam repair's local adjacency update equals a full rebuild."""

import numpy as np
from scipy.spatial import Delaunay

from fvcom_mesh_tools.patch import (
    _ccw,
    _convex_quad,
    _face_adjacency,
    _flip_in_place,
    _interior_edges,
    _unflip,
    improve_patch,
)


def _mesh(n=400, seed=0):
    rng = np.random.default_rng(seed)
    xy = rng.random((n, 2)) * 100.0
    return xy, Delaunay(xy).simplices.astype(np.int64)


def test_flip_in_place_matches_rebuild_and_unflip_restores():
    xy, tri = _mesh()
    adj = _face_adjacency(tri)
    rng = np.random.default_rng(1)
    n_done = 0
    for _ in range(300):
        ea, eb, f0s, f1s = _interior_edges(tri)
        k = int(rng.integers(len(ea)))
        a, b, f0, f1 = int(ea[k]), int(eb[k]), int(f0s[k]), int(f1s[k])
        c = int(np.setdiff1d(tri[f0], [a, b])[0])
        d = int(np.setdiff1d(tri[f1], [a, b])[0])
        if not _convex_quad(xy, a, b, c, d):
            continue
        tri0, adj0 = tri.copy(), adj.copy()
        undo = _flip_in_place(tri, adj, f0, f1, _ccw(xy, np.array([c, d, b])),
                              _ccw(xy, np.array([d, c, a])))
        np.testing.assert_array_equal(adj, _face_adjacency(tri))
        if rng.random() < 0.5:
            _unflip(tri, adj, undo)
            np.testing.assert_array_equal(tri, tri0)
            np.testing.assert_array_equal(adj, adj0)
        n_done += 1
    assert n_done > 100


def test_improve_patch_on_a_perturbed_mesh_still_improves():
    xy, tri = _mesh(600, seed=3)
    rng = np.random.default_rng(4)
    hull = np.zeros(len(xy), dtype=bool)
    hull[np.unique(Delaunay(xy).convex_hull)] = True
    movable = ~hull
    xy2 = xy.copy()
    xy2[movable] += rng.normal(0, 0.3, (movable.sum(), 2))
    mutable = np.ones(len(tri), dtype=bool)
    out_xy, out_tri, info = improve_patch(xy2, tri, movable, mutable)
    p = out_xy[out_tri]
    area2 = ((p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
             - (p[:, 1, 1] - p[:, 0, 1]) * (p[:, 2, 0] - p[:, 0, 0]))
    assert (area2 > 0).all()
    assert info["n_flips"] + info["n_moves"] > 0
