# (A) offshore transect through each finalize stage: which stage
#     floors the sizing at ~1.2 km?
# (B) why does collapse_thin_triangles not collapse the q=0.019
#     interior sliver?
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox = (166.0, 176.0, -48.0, -40.0)
shore = Shoreline(str(OM2D / "datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox, 1e3 * DEG)
sdf = om.signed_distance_function(shore)
feat = om.feature_sizing_function(shore, sdf, r=3,
                                  max_edge_length=100e3 * DEG)
# transect: due west from Hokitika coast (~170.9,-42.7) into Tasman
pts = np.column_stack([np.linspace(170.8, 167.0, 12),
                       np.full(12, -42.7)])

def probe(grid, tag):
    v = grid.eval(pts) / DEG
    print(f"[stage] {tag:14s}: " +
          " ".join(f"{x:8.0f}" for x in v), flush=True)

probe(feat, "feature")
from oceanmesh.finalize import (enforce_nearshore_max_edge,
                                compute_minimum)
from oceanmesh import enforce_mesh_gradation
g = compute_minimum([feat]); probe(g, "min-combine")
g = enforce_nearshore_max_edge(g, shore, 5e3 * DEG); probe(g, "ns-cap")
g.values[g.values > 100e3 * DEG] = 100e3 * DEG; probe(g, "max_el")
g.values[g.values < 1e3 * DEG] = 1e3 * DEG
g.build_interpolant(); probe(g, "hmin floor")
g = enforce_mesh_gradation(g, gradation=0.35); probe(g, "gradation")

# (B) collapse probe
from oceanmesh.fix_mesh import simp_qual
from oceanmesh.mesh_improve import collapse_thin_triangles
OUT = Path("outputs/om2d_examples/ex1_nz")
p = np.load(OUT / "p.npy"); t = np.load(OUT / "t.npy")
q = simp_qual(p, t)
w = int(np.argmin(q))
print(f"[collapse] worst q={q[w]:.4f} tri verts=\n{p[t[w]]}", flush=True)
p2, t2 = collapse_thin_triangles(p.copy(), t.copy(), min_qual=0.25)
q2 = simp_qual(p2, t2)
print(f"[collapse] after: NT {len(t)}->{len(t2)} min {q2.min():.4f} "
      f"n(q<0.25): {int((q<0.25).sum())}->{int((q2<0.25).sum())}",
      flush=True)
