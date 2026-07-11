# Dissect one parity-broken row: dump every edge crossing the ray
# for a red-row point, compare with shapely ground truth.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from oceanmesh import Shoreline
from oceanmesh import edges as om_edges
from oceanmesh.geometry import inpoly2

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox_poly = np.array([
    [-71.6, 42.7], [-64, 30], [-80, 24], [-85, 38], [-71.6, 42.7]])
shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox_poly, 1e3*DEG)
outer = np.asarray(shore.boubox)
inner = np.asarray(shore.inner)
poly = np.vstack((outer, inner))
e = om_edges.get_poly_edges(poly)
pn = np.nan_to_num(poly)

# find a red row on the earlier lattice
x = np.linspace(-80.4, -78.0, 1200)
y = np.linspace(24.2, 26.2, 1000)
X, Y = np.meshgrid(x, y)
q = np.column_stack([X.ravel(), Y.ravel()])
ins, _ = inpoly2(q, pn, e)
G = ins.reshape(X.shape)
row_frac = G.mean(axis=1)
# rows fully red inside the mostly-blue band 24.9..26.1
cand = [(i, f) for i, f in enumerate(row_frac)
        if y[i] > 24.95 and f < 0.30]
print("[row] suspicious rows:", [(i, round(y[i], 5), round(f, 2))
                                 for i, f in cand[:8]], flush=True)
if cand:
    i = cand[0][0]
    qy = y[i]
    qx = -79.0
    import shapely.geometry as sg
    # ground truth via shapely on the outer ring only
    print(f"[row] probing point ({qx}, {qy})", flush=True)
    x1 = pn[e[:, 0], 0]; y1 = pn[e[:, 0], 1]
    x2 = pn[e[:, 1], 0]; y2 = pn[e[:, 1], 1]
    ca = y1 > qy; cb = y2 > qy
    cross = ca != cb
    dy = y2 - y1
    with np.errstate(divide="ignore", invalid="ignore"):
        xint = x1 + (qy - y1) * (x2 - x1) / dy
    left = cross & (xint > -85.0) & (xint < qx)
    print(f"[row] crossings left of point: {int(left.sum())}",
          flush=True)
    for j in np.where(left)[0]:
        print(f"  edge {j}: ({x1[j]:.8f},{y1[j]:.10f}) -> "
              f"({x2[j]:.8f},{y2[j]:.10f}) xint={xint[j]:.8f} "
              f"dy={dy[j]:.3e}", flush=True)
    # vertex ties near qy
    vt = np.where(np.abs(pn[:, 1] - qy) < 1e-9)[0]
    print(f"[row] vertices with y within 1e-9 of row: {len(vt)}",
          flush=True)
    for v in vt[:10]:
        print(f"  vert {v}: ({pn[v,0]:.10f},{pn[v,1]:.12f}) "
              f"dy_exact={pn[v,1]-qy:.3e}", flush=True)
    # row neighbours for context
    for di in (-1, 1):
        insd, _ = inpoly2(np.array([[qx, y[i+di]]]), pn, e)
        print(f"[row] neighbour row y={y[i+di]:.6f}: inside={bool(insd[0])}",
              flush=True)
    insd, _ = inpoly2(np.array([[qx, qy]]), pn, e)
    print(f"[row] the point itself: inside={bool(insd[0])}", flush=True)
