"""Interpolation over a mesh with a split pinch vertex (coincident nodes)."""

import numpy as np

from fvcom_mesh_tools.refine import depths_from_base
from fvcom_mesh_tools.tri_interp import linear_interpolator, merged_triangulation

# Two fans that touch at (1, 1): the lower-left square and the upper-right
# square, each of two triangles. The pinch vertex is split into nodes 2 and 7.
XY = np.array([[0, 0], [1, 0], [1, 1], [0, 1],      # lower-left square, pinch = 2
               [2, 1], [2, 2], [1, 2], [1, 1]],     # upper-right square, pinch = 7
              dtype=float)
TRI = np.array([[0, 1, 2], [0, 2, 3], [7, 4, 5], [7, 5, 6]])
VAL = np.array([0.0, 1.0, 2.0, 1.0, 3.0, 4.0, 3.0, 2.0])


def test_merged_interpolation_matches_each_fan():
    f = linear_interpolator(XY, TRI, VAL)
    q = np.array([[0.75, 0.25], [0.25, 0.75], [1.5, 1.5], [1.0, 1.0]])
    got = np.asarray(f(q[:, 0], q[:, 1]))
    # value = x + y in the lower square, x + y in the upper one too
    np.testing.assert_allclose(got, q.sum(1), atol=1e-12)
    assert np.ma.is_masked(f(1.75, 0.25))   # outside both fans


def test_no_duplicates_is_unchanged():
    tri, first = merged_triangulation(XY[:7], TRI[:2])
    assert len(tri.x) == 7 and (first == np.arange(7)).all()


def test_depths_from_base_with_split_pinch():
    dep, outside = depths_from_base(XY, TRI, VAL + 3.0, np.array([[0.5, 0.25], [1.5, 1.75]]))
    np.testing.assert_allclose(dep, [3.75, 6.25])
    assert outside == 0
