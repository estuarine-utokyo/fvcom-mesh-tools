# REAL-DATA equivalence test: old Cython kernel vs dengwirda inpoly
# on the exact (node, edge) arrays the SDF uses for nest 3.
import os, time
import numpy as np
os.environ.setdefault("MPLBACKEND", "Agg")
import oceanmesh as om
from oceanmesh import Shoreline
from oceanmesh.geometry.point_in_polygon import inpoly2 as kernel_inpoly
from oceanmesh import edges as om_edges
from inpoly import inpoly2 as ext_inpoly

DEG = 1.0/111e3
import numpy as np
x3 = [139.142156,139.142156,139.144724,139.152351,139.164815,139.18178,139.202839,139.227552,139.25548,139.286203,139.31934,139.354542,139.391505,139.429958,139.469663,139.510408,139.552007,139.594291,139.637105,139.680308,139.723765,139.767348,139.810932,139.854389,139.897592,139.940406,139.98269,140.024289,140.065034,140.104739,140.143192,140.180155,140.215357,140.248493,140.279217,140.307145,140.331858,140.352917,140.369882,140.382346,140.389973,140.392541,140.392541,139.142156]
y3 = [35.859592,35.225989,35.182486,35.139579,35.097818,35.057672,35.019512,34.983607,34.95014,34.919216,34.890891,34.865177,34.842062,34.82152,34.803513,34.788004,34.774954,34.764328,34.756096,34.750234,34.746724,34.745556,34.746724,34.750234,34.756096,34.764328,34.774954,34.788004,34.803513,34.82152,34.842062,34.865177,34.890891,34.919216,34.95014,34.983607,35.019512,35.057672,35.097818,35.139579,35.182486,35.225989,35.859592,35.859592]
bbox3 = np.column_stack([x3, y3])
shore = Shoreline("outputs/tb_varres_3r/land_osm_wide.shp", bbox3, 100.0*DEG)
poly = np.vstack((shore.inner, shore.boubox))
node = np.nan_to_num(poly)
edge = om_edges.get_poly_edges(poly)
print(f"[equiv] nodes {len(node):,} edges {len(edge):,}", flush=True)
rng = np.random.default_rng(7)
q = np.column_stack([rng.uniform(139.15, 140.39, 300000),
                     rng.uniform(34.75, 35.86, 300000)])
os.environ["OCEANMESH_INPOLY_METHOD"] = "cython"
t0=time.time(); s_old, b_old = kernel_inpoly(q, node, edge); t1=time.time()
print(f"[equiv] cython kernel: {t1-t0:.2f} s", flush=True)
os.environ["OCEANMESH_INPOLY_METHOD"] = ""

t0=time.time(); s_ext, b_ext = ext_inpoly(q, node, edge); t1=time.time()
print(f"[equiv] ext:    {t1-t0:.2f} s", flush=True)

# cleaned-input variant: drop unreferenced nodes (incl. the
# nan_to_num (0,0) dummies), remap edges, drop zero-length edges
used = np.unique(edge)
remap = np.full(len(node), -1, int)
remap[used] = np.arange(len(used))
node_c = node[used]
edge_c = remap[edge]
seglen = np.linalg.norm(node_c[edge_c[:,0]] - node_c[edge_c[:,1]], axis=1)
edge_c = edge_c[seglen > 0]
print(f"[equiv] cleaned: nodes {len(node_c):,} edges {len(edge_c):,} "
      f"(dropped {len(edge)-len(edge_c)} zero-len)", flush=True)
t0=time.time(); s_ext2, _ = ext_inpoly(q, node_c, edge_c); t1=time.time()
print(f"[equiv] ext(clean): {t1-t0:.2f} s", flush=True)
mism2 = s_old != s_ext2
print(f"[equiv] CLEANED mismatches: {int(mism2.sum()):,} / {len(q):,}",
      flush=True)
mism = s_old != s_ext
print(f"[equiv] STAT mismatches: {int(mism.sum()):,} / {len(q):,} "
      f"({100*mism.mean():.3f}%)", flush=True)
# brute-force crossing-number referee on mismatches
def brute(qpts, node, edge):
    x1 = node[edge[:,0],0]; y1 = node[edge[:,0],1]
    x2 = node[edge[:,1],0]; y2 = node[edge[:,1],1]
    out = np.zeros(len(qpts), bool)
    for i,(qx,qy) in enumerate(qpts):
        cond = (y1 > qy) != (y2 > qy)
        with np.errstate(divide='ignore', invalid='ignore'):
            xint = x1 + (qy - y1)*(x2 - x1)/(y2 - y1)
        out[i] = (np.count_nonzero(cond & (qx < xint)) % 2) == 1
    return out

import importlib
import oceanmesh.geometry.point_in_polygon as pip
importlib.reload(pip)
t0=time.time(); s_nb, b_nb = pip._inpoly_numba(q, node, edge, 5e-14); t1=time.time()
print(f"[equiv] numba(cold+jit): {t1-t0:.2f} s", flush=True)
t0=time.time(); s_nb, b_nb = pip._inpoly_numba(q, node, edge, 5e-14); t1=time.time()
print(f"[equiv] numba(warm): {t1-t0:.2f} s", flush=True)
mn = s_old != s_nb
print(f"[equiv] NUMBA vs kernel mismatches: {int(mn.sum()):,} / {len(q):,}",
      flush=True)

sel = np.random.default_rng(1).choice(len(q), 300, replace=False)
ref = brute(q[sel], node, edge)
print(f"[equiv] REFEREE(300 random): cython {100*(ref==s_old[sel]).mean():.1f}% "
      f"| numba {100*(ref==s_nb[sel]).mean():.1f}%", flush=True)

if mism.any():
    pts = q[mism][:8]
    print("[equiv] sample mismatch pts:", np.round(pts, 4).tolist(),
          flush=True)
    # distance of mismatch pts to nearest node (boundary proximity?)
    from scipy.spatial import cKDTree
    d, _ = cKDTree(poly[~np.isnan(poly[:,0])]).query(q[mism])
    print(f"[equiv] mismatch dist-to-shoreline deg: p50={np.percentile(d,50):.5f} "
          f"p90={np.percentile(d,90):.5f} max={d.max():.4f}", flush=True)
