"""Linear interpolation over a mesh that may have coincident nodes.

A finished mesh can carry two nodes at the same coordinates: where the
coastline pinches (two land pieces touching at a point) the finishing splits
the non-manifold vertex so each fan has its own node, and FVCOM is fine with
that. Matplotlib's ``TrapezoidMapTriFinder`` is not -- it refuses the whole
triangulation as invalid -- so every ``LinearTriInterpolator`` over such a
mesh failed (the Futtsu 50 m mesh: 7 split pinch nodes by the port; the 5 m
refinement on it stopped in its first minute, job 137981).

For interpolation, coincident nodes are merged onto one point (the first in
node order) and the elements are renumbered onto it; a split vertex's fans
then share the point again, which is a valid triangulation. The value at a
merged point is its first node's value: split nodes are copies, and where a
later smoothing moved them apart (by centimetres at Futtsu) the choice is
confined to the pinch itself.
"""

from __future__ import annotations

import numpy as np


def merged_triangulation(xy, elements):
    """``(Triangulation, first)``: coincident nodes merged; ``first[k]`` is the
    node whose value point ``k`` takes."""
    from matplotlib.tri import Triangulation

    xy = np.asarray(xy, dtype=float)[:, :2]
    tri = np.asarray(elements, dtype=np.int64)
    uniq, first, inv = np.unique(xy, axis=0, return_index=True, return_inverse=True)
    inv = np.asarray(inv).ravel()
    if len(uniq) == len(xy):
        return Triangulation(xy[:, 0], xy[:, 1], tri), np.arange(len(xy))
    t = inv[tri]
    keep = (t[:, 0] != t[:, 1]) & (t[:, 1] != t[:, 2]) & (t[:, 0] != t[:, 2])
    return Triangulation(uniq[:, 0], uniq[:, 1], t[keep]), first


def linear_interpolator(xy, elements, values):
    """A ``LinearTriInterpolator`` of per-node ``values`` that tolerates
    coincident nodes (see the module docstring)."""
    from matplotlib.tri import LinearTriInterpolator

    tri, first = merged_triangulation(xy, elements)
    return LinearTriInterpolator(tri, np.asarray(values, dtype=float)[first])
