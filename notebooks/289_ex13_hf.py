# Ladder step 13: Example_13 NZ high-fidelity (mesh1d pfix/egfix
# chains; cleanup and egfix-healing OFF per meshgen.m:537,865).
import os, sys, logging, time
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(900, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import Region, Shoreline
from oceanmesh.mesh1d import mesh1d
from oceanmesh import edges as om_edges
from oceanmesh.geometry import inpoly2

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
NZ2 = DS / "New_Zealand_Marine_Environment_Classification_Feature_Layers-shp"
OUT = Path("outputs/om2d_examples/ex13hf")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

# nest 1: NZ coarse
reg1 = Region((166.0, 176.0, -48.0, -40.0), 4326)
sh1 = Shoreline(str(DS/"GSHHS_shp/f/GSHHS_f_L1.shp"), reg1.bbox,
                1e3*DEG)
sdf1 = om.signed_distance_function(sh1)
f1 = om.feature_sizing_function(sh1, sdf1, r=3,
                                max_edge_length=100e3*DEG)
g1, _ = om.finalize_sizing([f1], shoreline=sh1, hmin=1e3,
                           max_edge_length=100e3,
                           max_edge_length_nearshore=5e3,
                           gradation=0.15)
print(f"[ex13] nest1 +{time.time()-t0:.0f}s", flush=True)

# nest 2: Lyttelton high-fidelity
bbox2 = (172.61348324, 172.85932924, -43.67102284, -43.51383230)
reg2 = Region(bbox2, 4326)
sh2 = Shoreline(str(NZ2/"Coastline_epsg4326.shp"), reg2.bbox,
                50.0*DEG)
sdf2 = om.signed_distance_function(sh2)
f2 = om.feature_sizing_function(sh2, sdf2, r=5,
                                max_edge_length=500*DEG)
g2, _ = om.finalize_sizing([f2], shoreline=sh2, hmin=50.0,
                           max_edge_length=500.0,
                           max_edge_length_nearshore=250.0,
                           gradation=0.05)
print(f"[ex13] nest2 +{time.time()-t0:.0f}s", flush=True)

# composite fh (deepest wins) for mesh1d
ring2 = np.array([[bbox2[0], bbox2[2]], [bbox2[1], bbox2[2]],
                  [bbox2[1], bbox2[3]], [bbox2[0], bbox2[3]],
                  [bbox2[0], bbox2[2]]])

def fh_comp(q):
    v = np.asarray(g1.eval(q), dtype=float)
    m = ((q[:,0]>=bbox2[0])&(q[:,0]<=bbox2[1])
         &(q[:,1]>=bbox2[2])&(q[:,1]<=bbox2[3]))
    if m.any():
        v[m] = np.asarray(g2.eval(q[m]), dtype=float)
    return v

# HF chains: in-box runs of sh2 parts -> mesh1d
np.random.seed(0)
parts = []
for arr in (sh2.mainland, sh2.inner):
    if arr is None or not len(arr):
        continue
    a = np.asarray(arr)
    idx = np.where(~np.isfinite(a[:, 0]))[0]
    s = 0
    for e in list(idx) + [len(a)]:
        seg = a[s:e]; s = e + 1
        if len(seg) < 2:
            continue
        inb = ((seg[:,0]>=bbox2[0])&(seg[:,0]<=bbox2[1])
               &(seg[:,1]>=bbox2[2])&(seg[:,1]<=bbox2[3]))
        j = 0
        while j < len(seg):
            if not inb[j]:
                j += 1; continue
            k = j
            while k < len(seg) and inb[k]:
                k += 1
            if k - j >= 2:
                parts.append(seg[j:k])
            j = k
pfix_l, egfix_l, off = [], [], 0
for part in parts:
    pf, ef = mesh1d(part, fh_comp, 50.0*DEG)
    if pf is None:
        continue
    pfix_l.append(pf); egfix_l.append(ef + off); off += len(pf)
pfix = np.vstack(pfix_l) if pfix_l else None
egfix = np.vstack(egfix_l) if egfix_l else None
if pfix is not None:
    # fixgeo2 analogue: weld coincident chain endpoints, remap
    # edges, drop degenerate/duplicate edges
    key = np.round(pfix / 1e-9).astype(np.int64)
    _, uidx, inv = np.unique(key, axis=0, return_index=True,
                             return_inverse=True)
    pfix = pfix[uidx]
    egfix = inv[egfix]
    egfix = egfix[egfix[:, 0] != egfix[:, 1]]
    egfix = np.unique(np.sort(egfix, axis=1), axis=0)
    print(f"[ex13] welded pfix {len(pfix):,} egfix {len(egfix):,}",
          flush=True)
print(f"[ex13] HF chains: {0 if pfix is None else len(pfix):,} pfix, "
      f"{len(parts)} parts +{time.time()-t0:.0f}s", flush=True)

p, t = om.generate_multiscale_mesh(
    [sdf1, sdf2], [g1, g2], max_iter=100, seed=0,
    gradation=[0.15, 0.05], pfix=pfix, egfix=egfix,
    cleanup="none", heal_fixed_edges_every=0)
print(f"[ex13] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p.npy", p); np.save(OUT/"t.npy", t)
np.save(OUT/"pfix.npy", pfix if pfix is not None else np.zeros((0,2)))
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=171 +lat_0=-44 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[ex13] qual min={q.min():.4f} mean={q.mean():.4f}", flush=True)
