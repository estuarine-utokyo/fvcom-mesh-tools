import os, sys, logging, time
import numpy as np
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
DEG = 1.0 / 111e3
dem_p = DS/"PR_1arcsec/pr_1s.nc"
dp = DEM(str(dem_p)); bbox = dp.bbox
print(f"[pr] bbox={bbox}", flush=True)
reg = Region(tuple(bbox), 4326)
sh = Shoreline(str(DS/"PR_1arcsec/pr_1s_0m_contour.shp"), reg.bbox,
               30.0*DEG)
sdf = om.signed_distance_function(sh)
dem = DEM(str(dem_p), bbox=reg)
f = om.feature_sizing_function(sh, sdf, r=-5,
                               max_edge_length=1000*DEG)
for nm, g_ in (("feat", f),):
    v = np.asarray(g_.values, dtype=float) / DEG
    print(f"[pr] {nm}: p5/50/95 = {np.nanpercentile(v,5):.0f}/"
          f"{np.nanpercentile(v,50):.0f}/"
          f"{np.nanpercentile(v,95):.0f} m  min={np.nanmin(v):.0f}",
          flush=True)
g, dta = om.finalize_sizing([f], dem=dem, shoreline=sh,
                            hmin=30.0, max_edge_length=1000.0,
                            gradation=0.25,
                            courant={"timestep": 0.0})
v = np.asarray(g.values, dtype=float) / DEG
print(f"[pr] FINAL: p5/50/95 = {np.nanpercentile(v,5):.0f}/"
      f"{np.nanpercentile(v,50):.0f}/{np.nanpercentile(v,95):.0f} m "
      f"dt={dta:.2f}", flush=True)
acc = np.minimum(1.0, (30*DEG/np.maximum(np.asarray(g.values,float),1e-12))**2)
print(f"[pr] mean acceptance = {np.nanmean(acc)*100:.2f}%", flush=True)
