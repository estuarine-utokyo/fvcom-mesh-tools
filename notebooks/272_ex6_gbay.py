# Ladder step 10: Example_6_GBAY (US Medium Shoreline + measured
# DEM + channels ch=0.1).
import os, sys, logging, time
os.environ["OM_TRACE_IMPROVE"] = "1"
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
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex6gbay")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

mm = loadmat(str(DS/"ECGC_Thalwegs.mat"), squeeze_me=False)
channels = [np.asarray(a, dtype=float)
            for a in mm["pts2"].ravel() if np.size(a) >= 4]

bbox = (-95.40, -94.4, 29.14, 30.09)
reg = Region(bbox, 4326)
sh = Shoreline(
    str(DS/"US_Medium_Shoreline/us_medium_shoreline_polygon.shp"),
    reg.bbox, 60.0*DEG)
sh.detect_inpoly_flip(
    str(DS/"GSHHS_shp/l/GSHHS_l_L1.shp"))
sdf = om.signed_distance_function(sh)
dem = DEM(str(DS/"galveston_13_mhw_2007.nc"), bbox=reg)
f = om.feature_sizing_function(sh, sdf, r=3,
                               max_edge_length=1e3*DEG,
                               lattice_anchor=(dem.bbox[0],
                                               dem.bbox[2]))
c = om.channel_sizing_function(
    dem, channels, ch=0.1,
    min_edge_length_channel=100.0,
    angle_of_reslope=60.0,
    min_edge_length=60.0,
    max_edge_length=1e3,
    dx=60.0*DEG)
_use_ch = os.environ.get("EX6_CHANNELS", "1") == "1"
g, _ = om.finalize_sizing([f, c] if _use_ch else [f],
                          dem=dem, shoreline=sh,
                          hmin=60.0, max_edge_length=1e3,
                          gradation=0.25)
print(f"[ex6] sizing +{time.time()-t0:.0f}s", flush=True)
p, t = om.generate_mesh(sdf, g, max_iter=100, seed=0)
print(f"[ex6] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p.npy", p); np.save(OUT/"t.npy", t)
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=-94.9 +lat_0=29.6 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[ex6] qual min={q.min():.4f} mean={q.mean():.4f}", flush=True)
