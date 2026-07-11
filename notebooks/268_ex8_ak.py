# Ladder step 9: Example_8_AK (polygon boubox from file, pure
# distance sizing 'dis', no DEM in edgefx).
import os, sys, logging, time
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(900, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex8ak")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

mm = loadmat(str(DS/"ak_outerpoly.mat"), squeeze_me=False)
key = [k for k in mm if not k.startswith("__")][0]
poly = np.asarray(mm[key], dtype=float)
print(f"[ex8] boubox var '{key}' shape {poly.shape} "
      f"lon [{np.nanmin(poly[:,0]):.2f},{np.nanmax(poly[:,0]):.2f}] "
      f"lat [{np.nanmin(poly[:,1]):.2f},{np.nanmax(poly[:,1]):.2f}]",
      flush=True)

sh = Shoreline(str(DS/"GSHHS_shp/f/GSHHS_f_L1.shp"), poly, 5e3*DEG)
sdf = om.signed_distance_function(sh)
f = om.distance_sizing_function(sh, rate=0.25,
                                max_edge_length=50e3*DEG)
g, _ = om.finalize_sizing([f], shoreline=sh, hmin=5e3,
                          max_edge_length=50e3, gradation=0.25)
print(f"[ex8] sizing +{time.time()-t0:.0f}s", flush=True)
p, t = om.generate_mesh(sdf, g, max_iter=100, seed=0)
print(f"[ex8] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p.npy", p); np.save(OUT/"t.npy", t)
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=-155 +lat_0=60 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[ex8] qual min={q.min():.4f} mean={q.mean():.4f}", flush=True)
