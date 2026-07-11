# Ladder completeness: Example_5_JBAY (no weirs; dt=2 explicit).
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
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DSET = OM2D / "datasets/PostSandyNCEI"
OUT = Path("outputs/om2d_examples/ex5jbay")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()
reg = Region((-73.97, -73.75, 40.5, 40.68), 4326)
sh = Shoreline(str(DSET/"PostSandyNCEI.shp"), reg.bbox, 15.0*DEG)
sdf = om.signed_distance_function(sh)
dem = DEM(str(DSET/"PostSandyNCEI.nc"), bbox=reg)
f = om.feature_sizing_function(
    sh, sdf, r=3, max_edge_length=1e3*DEG,
    lattice_anchor=(dem.bbox[0], dem.bbox[2]))
g, _ = om.finalize_sizing([f], dem=dem, shoreline=sh, hmin=15.0,
                          max_edge_length=1e3, gradation=0.15,
                          courant={"timestep": 2.0})
print(f"[ex5] sizing +{time.time()-t0:.0f}s", flush=True)
p, t = om.generate_mesh(sdf, g, max_iter=100, seed=0)
print(f"[ex5] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p.npy", p); np.save(OUT/"t.npy", t)
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=-73.86 +lat_0=40.59 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[ex5] qual min={q.min():.4f} mean={q.mean():.4f}", flush=True)
