# Find the exact NaN-part whose edge set loses a crossing.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path as P
from oceanmesh import DEM, Region, Shoreline
from oceanmesh import edges as om_edges

DSET = P(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
dem_probe = DEM(str(DSET / "PostSandyNCEI.nc"))
reg = Region(dem_probe.bbox, 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  30.0 / 111e3)
q = np.array([-74.100, 40.777])

def parts(arr):
    a = np.asarray(arr, dtype=float)
    idx = np.where(np.isnan(a[:, 0]))[0]
    start = 0
    for stop in list(idx) + [len(a)]:
        seg = a[start:stop]; start = stop + 1
        if len(seg) >= 3:
            yield seg

def cross_edges(seg, q):
    # crossings east of q using OUR edge builder on this part alone
    part = np.vstack((seg, [[np.nan, np.nan]]))
    e = om_edges.get_poly_edges(part)
    an = np.nan_to_num(part)
    p1 = an[e[:, 0]]; p2 = an[e[:, 1]]
    ca = p1[:, 1] > q[1]; cb = p2[:, 1] > q[1]
    cr = ca != cb
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = p1[:, 0] + (q[1] - p1[:, 1]) * (p2[:, 0] - p1[:, 0]) \
            / (p2[:, 1] - p1[:, 1])
    return int((cr & (xi > q[0])).sum())

def cross_true(seg, q):
    # ground truth: explicit closed ring
    ring = np.vstack((seg, seg[0]))
    p1 = ring[:-1]; p2 = ring[1:]
    ca = p1[:, 1] > q[1]; cb = p2[:, 1] > q[1]
    cr = ca != cb
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = p1[:, 0] + (q[1] - p1[:, 1]) * (p2[:, 0] - p1[:, 0]) \
            / (p2[:, 1] - p1[:, 1])
    return int((cr & (xi > q[0])).sum())

for label, arr in (("mainland", shore.mainland),
                   ("inner", shore.inner)):
    for k, seg in enumerate(parts(arr)):
        a_, b_ = cross_edges(seg, q), cross_true(seg, q)
        if a_ != b_:
            closed = np.allclose(seg[0], seg[-1])
            d = np.hypot(*(seg[0] - seg[-1]))
            print(f"[hunt] {label} part {k}: ours={a_} true={b_} "
                  f"len={len(seg)} closed={closed} "
                  f"first-last-gap={d:.3e}", flush=True)
            print(f"       first={seg[0]} last={seg[-1]}",
                  flush=True)
print("[hunt] done", flush=True)
