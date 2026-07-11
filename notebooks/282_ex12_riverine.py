# Ladder step 11: Example_12_Riverine (Pearl River delta;
# explicit dt=5 CFL floor; polygon boubox).
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
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex12riv")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

bbox = np.array([[112, 20], [112, 24], [116, 24], [116, 20],
                 [112, 20]], dtype=float)
reg = Region((112.0, 116.0, 20.0, 24.0), 4326)
sh = Shoreline(str(DS/"GSHHS_shp/f/GSHHS_f_L1.shp"), bbox, 100.0*DEG)
sdf = om.signed_distance_function(sh)
dem = DEM(str(DS/"SRTM15+.nc"), bbox=reg)
f = om.feature_sizing_function(
    sh, sdf, r=3, max_edge_length=10e3*DEG,
    lattice_anchor=(dem.bbox[0], dem.bbox[2]))
g, _ = om.finalize_sizing([f], dem=dem, shoreline=sh, hmin=100.0,
                          max_edge_length=10e3, gradation=0.3,
                          courant={"timestep": 5.0})
print(f"[ex12] sizing +{time.time()-t0:.0f}s", flush=True)
p, t = om.generate_mesh(sdf, g, max_iter=100, seed=0)
print(f"[ex12] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p.npy", p); np.save(OUT/"t.npy", t)
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=114 +lat_0=22 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[ex12] qual min={q.min():.4f} mean={q.mean():.4f}", flush=True)
