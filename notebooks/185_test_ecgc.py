# Ladder step 2c: port of OM2D Tests/TestECGC.m — US East Coast
# polygon domain, feature+wavelength sizing with two-sided Courant
# bounds (dt=50 s, Cr in [0.5, 6.0]), then bound_courant_number
# (dt=50, max_cr=4, min_cr=0.25) and assert CalcCFL in bounds.
import os, sys, logging, time
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
bbox_poly = np.array([
    [-71.6, 42.7], [-64, 30], [-80, 24], [-85, 38], [-71.6, 42.7]])
min_el, max_el, wl, grade, R = 1e3, 50e3, 30, 0.15, 3
dt, cr_min, cr_max = 50.0, 0.5, 6.0

t0 = time.time()
reg = Region((float(bbox_poly[:,0].min()), float(bbox_poly[:,0].max()),
              float(bbox_poly[:,1].min()), float(bbox_poly[:,1].max())),
             4326)
shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                  bbox_poly, min_el*DEG)
sdf = om.signed_distance_function(shore)
dem = DEM(str(OM2D/"datasets/SRTM15+.nc"), bbox=reg)
print(f"[ecgc] geodata+dem +{time.time()-t0:.0f}s", flush=True)
feat = om.feature_sizing_function(shore, sdf, r=R,
                                  max_edge_length=max_el*DEG)
wlg = om.wavelength_sizing_function(dem, wl=wl,
                                    min_edgelength=min_el*DEG,
                                    max_edge_length=max_el*DEG)
grid, dt_used = om.finalize_sizing(
    [feat, wlg], dem=dem, shoreline=shore, hmin=min_el,
    max_edge_length=max_el, gradation=grade,
    courant={"timestep": dt, "max": cr_max, "min": cr_min},
)
print(f"[ecgc] sizing +{time.time()-t0:.0f}s", flush=True)
p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=0)
print(f"[ecgc] mesh NP={len(p):,} +{time.time()-t0:.0f}s", flush=True)
b = om.interp_bathymetry(p, t, dem, method="cell-averaging")
p2, t2, b2 = om.bound_courant_number(p, t, b, 50.0, cr_max=4.0,
                                     cr_min=0.25, max_iter=20)
cr = om.calc_cfl(p2, t2, b2, dt=50.0)
hi, lo = int((cr > 4.0).sum()), int((cr < 0.25).sum())
print(f"[ecgc] NP {len(p):,}->{len(p2):,} Cr max={cr.max():.2f} "
      f"min={cr.min():.2f} viol>4: {hi} viol<0.25: {lo}", flush=True)
print(f"[ecgc] {'PASSED' if hi == 0 and lo == 0 else 'FAILED'}",
      flush=True)

OUTD = Path("outputs/om2d_examples/ecgc")
OUTD.mkdir(parents=True, exist_ok=True)
np.save(OUTD / "p.npy", p2)
np.save(OUTD / "t.npy", t2)
np.save(OUTD / "b.npy", b2)
print("[ecgc] saved", flush=True)
