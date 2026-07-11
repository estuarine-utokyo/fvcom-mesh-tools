# Ladder step 12a: Example_11 BASE mesh (NZ, h0=1km, max_el_ns).
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

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex11rm")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()
reg = Region((166.0, 176.0, -48.0, -40.0), 4326)
sh = Shoreline(str(DS/"GSHHS_shp/f/GSHHS_f_L1.shp"), reg.bbox,
               1e3*DEG)
sdf = om.signed_distance_function(sh)
f = om.feature_sizing_function(sh, sdf, r=3,
                               max_edge_length=100e3*DEG)
g, _ = om.finalize_sizing([f], shoreline=sh, hmin=1e3,
                          max_edge_length=100e3,
                          max_edge_length_nearshore=5e3,
                          gradation=0.35)
print(f"[ex11] sizing +{time.time()-t0:.0f}s", flush=True)
p, t = om.generate_mesh(sdf, g, max_iter=100, seed=0)
print(f"[ex11] base NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p_base.npy", p); np.save(OUT/"t_base.npy", t)
import pickle
with open(OUT/"fh.pkl", "wb") as fh_:
    pickle.dump({"bbox": g.bbox, "dx": g.dx, "dy": g.dy,
                 "values": np.asarray(g.values)}, fh_)
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=171 +lat_0=-44 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[ex11] qual min={q.min():.4f} mean={q.mean():.4f}", flush=True)
