# Why does matplotlib's TriFinder refuse a mesh that verify_patch accepts?
#
# "Triangulation is invalid" comes from the trapezoid map, which needs a
# simple, non-overlapping triangulation.  verify_patch checks a different and
# weaker list.  This says which of the usual causes is present.
#
#   python notebooks/422_mesh_validity.py <mesh.14>
import sys
from pathlib import Path

import numpy as np
import shapely

from fvcom_mesh_tools.io.fort14 import read_fort14

m = read_fort14(Path(sys.argv[1]).resolve())
xy, tri = m.nodes[:, :2], m.elements
print(f"NP={m.n_nodes:,} NE={m.n_elements:,}")

u, inv, cnt = np.unique(np.round(xy, 9), axis=0, return_inverse=True,
                        return_counts=True)
print(f"exactly duplicate nodes: {int((cnt > 1).sum())} position(s), "
      f"{int(cnt[cnt > 1].sum() - (cnt > 1).sum())} extra node(s)")
if (cnt > 1).any():
    for k in np.flatnonzero(cnt > 1)[:5]:
        ids = np.flatnonzero(inv == k)
        print(f"  {u[k]} <- nodes {ids.tolist()}")

a = xy[tri[:, 1]] - xy[tri[:, 0]]
b = xy[tri[:, 2]] - xy[tri[:, 0]]
area2 = a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]
print(f"zero-area elements: {int((np.abs(area2) < 1e-9).sum())}; "
      f"clockwise: {int((area2 < 0).sum())}; area min {0.5 * np.abs(area2).min():.3e} m2")
deg = np.array([np.degrees(np.arccos(np.clip(
    np.sum((xy[tri[:, (i + 1) % 3]] - xy[tri[:, i]])
           * (xy[tri[:, (i + 2) % 3]] - xy[tri[:, i]]), axis=1)
    / (np.linalg.norm(xy[tri[:, (i + 1) % 3]] - xy[tri[:, i]], axis=1)
       * np.linalg.norm(xy[tri[:, (i + 2) % 3]] - xy[tri[:, i]], axis=1)), -1, 1)))
    for i in range(3)])
_worst = int(np.unravel_index(deg.argmin(), deg.shape)[1])
print(f"min angle {deg.min():.3f} deg in element {_worst}")

e = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
ue, ce = np.unique(e, axis=0, return_counts=True)
print(f"edges used by >2 faces (non-manifold): {int((ce > 2).sum())}")
bnd = ue[ce == 1]
ring = shapely.MultiLineString([shapely.LineString(xy[[i, j]]) for i, j in bnd.tolist()])
print(f"boundary edges {len(bnd):,}; boundary is simple: {ring.is_simple}")
# The trapezoid map's real complaint is usually overlapping faces, which show
# up as boundary segments that cross each other.
from shapely.ops import unary_union  # noqa: E402

merged = unary_union(ring)
print(f"boundary merges to {merged.geom_type} with "
      f"{len(getattr(merged, 'geoms', [merged]))} part(s)")
crossings = 0
idx = shapely.STRtree([shapely.LineString(xy[[i, j]]) for i, j in bnd.tolist()])
segs = [shapely.LineString(xy[[i, j]]) for i, j in bnd.tolist()]
for k, s in enumerate(segs):
    for j in idx.query(s):
        if j <= k:
            continue
        if shapely.crosses(s, segs[j]):
            crossings += 1
            if crossings <= 5:
                print(f"  boundary segments {bnd[k].tolist()} and "
                      f"{bnd[j].tolist()} CROSS")
print(f"crossing boundary segment pairs: {crossings}")
