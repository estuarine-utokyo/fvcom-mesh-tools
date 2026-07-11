import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import Shoreline
from oceanmesh.fix_mesh import simp_qual

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox = (166.0, 176.0, -48.0, -40.0)
shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox, 1e3*DEG)
sdf = om.signed_distance_function(shore)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=100e3*DEG)
grid, _ = om.finalize_sizing([feat], shoreline=shore, hmin=1e3,
                             max_edge_length=100e3,
                             max_edge_length_nearshore=5e3,
                             gradation=0.35)
p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=0,
                        cleanup="none")
CPP = np.cos(np.deg2rad(-44.0))
pp = p.copy(); pp[:, 0] *= CPP
q = simp_qual(pp, t)
bad = np.where(q < 0.01)[0]
print(f"[degen] NT={len(t):,} q<0.01: {len(bad)}", flush=True)
for e in bad[:12]:
    tri = pp[t[e]]
    L = np.linalg.norm(np.roll(tri, -1, 0) - tri, axis=1) / DEG
    c = p[t[e]].mean(0)
    on_frame = (abs(c[0]-166) < .01 or abs(c[0]-176) < .01
                or abs(c[1]+48) < .01 or abs(c[1]+40) < .01)
    d = sdf.eval(p[t[e]])
    print(f"[degen] q={q[e]:.5f} at ({c[0]:.3f},{c[1]:.3f}) "
          f"frame={on_frame} bars_km={np.round(L/1e3,2)} "
          f"sdf={np.round(d/DEG/1e3,2)}km", flush=True)

np.save("outputs/om2d_examples/ex1_nz/p_raw.npy", p)
np.save("outputs/om2d_examples/ex1_nz/t_raw.npy", t)
om.write_fort14("outputs/om2d_examples/ex1_nz/ours_raw.14", p, t,
                depth=np.ones(len(p)))
print("[degen] raw exported", flush=True)
