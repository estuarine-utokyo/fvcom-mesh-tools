# Ladder step 8: Example_9_TBAY (channels via thalwegs + the
# enforceMin=0/1 pair).
import os, sys, logging, time
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(900, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex9tbay")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

mm = loadmat(str(DS/"ECGC_Thalwegs.mat"), squeeze_me=False)
channels = [np.asarray(a, dtype=float)
            for a in mm["pts2"].ravel() if np.size(a) >= 4]
print(f"[ex9] {len(channels)} thalweg polylines", flush=True)

def build(bbox, with_channels):
    reg = Region(bbox, 4326)
    sh = Shoreline(str(DS/"GSHHS_shp/f/GSHHS_f_L1.shp"),
                   reg.bbox, 100.0*DEG)
    sdf = om.signed_distance_function(sh)
    dem = DEM(str(DS/"SRTM15+.nc"), bbox=reg)
    f = om.feature_sizing_function(sh, sdf, r=3,
                                   max_edge_length=1e3*DEG)
    fns = [f]
    if with_channels:
        c = om.channel_sizing_function(
            dem, channels, ch=1.0,
            min_edge_length_channel=100.0,
            angle_of_reslope=60.0,
            min_edge_length=100.0,
            max_edge_length=1e3,
            dx=100.0*DEG)
        fns.append(c)
    g, _ = om.finalize_sizing(fns, dem=dem, shoreline=sh,
                              hmin=100.0, max_edge_length=1e3,
                              gradation=0.35)
    return sdf, g

d1, g1 = build((-83.0, -82.0, 27.0, 28.5), True)
print(f"[ex9] nest1 (channels) +{time.time()-t0:.0f}s", flush=True)
d2, g2 = build((-82.8, -82.4, 27.25, 28.25), False)
print(f"[ex9] nest2 +{time.time()-t0:.0f}s", flush=True)

for tag, emin in (("nomin", False), ("min", True)):
    p, t = om.generate_multiscale_mesh(
        [d1, d2], [g1, g2], max_iter=100, seed=0,
        gradation=0.35, enforce_min=emin)
    np.save(OUT/f"p_{tag}.npy", p); np.save(OUT/f"t_{tag}.npy", t)
    from oceanmesh.fix_mesh import simp_qual
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:4326",
        "+proj=tmerc +lon_0=-82.5 +lat_0=27.75 +ellps=WGS84 +units=m",
        always_xy=True)
    xx, yy = tr.transform(p[:, 0], p[:, 1])
    q = simp_qual(np.column_stack([xx, yy]), t)
    print(f"[ex9] {tag}: NP={len(p):,} NT={len(t):,} "
          f"qual min={q.min():.4f} mean={q.mean():.4f} "
          f"+{time.time()-t0:.0f}s", flush=True)
    # rebuild sizing grids fresh for the second pass (enforce_min
    # mutates the finer grids in place)
    if tag == "nomin":
        d1, g1 = build((-83.0, -82.0, 27.0, 28.5), True)
        d2, g2 = build((-82.8, -82.4, 27.25, 28.25), False)
print("[ex9] done", flush=True)
