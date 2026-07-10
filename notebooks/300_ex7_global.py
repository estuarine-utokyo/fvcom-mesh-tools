# Ladder final: Example_7_Global (4-20 km, stereo generation).
import os, sys, logging, time
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(900, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from oceanmesh.region import to_stereo

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex7glob")
OUT.mkdir(parents=True, exist_ok=True)
SCR = Path(os.environ.get("EX7_SCR", "/tmp"))
DEG = 1.0 / 111e3
t0 = time.time()

# --- merge GSHHS L1 + L6 (lonlat) once
mll = SCR / "gshhs_l1l6.shp"
mst = SCR / "gshhs_l1l6_stereo.shp"
if not mll.exists():
    g1 = gpd.read_file(DS/"GSHHS_shp/f/GSHHS_f_L1.shp")
    g6 = gpd.read_file(DS/"GSHHS_shp/f/GSHHS_f_L6.shp")
    gm = gpd.GeoDataFrame(pd.concat([g1, g6], ignore_index=True),
                          crs=g1.crs)
    gm.to_file(mll)
    print(f"[ex7] merged lonlat shp ({len(gm)}) +{time.time()-t0:.0f}s",
          flush=True)
    # stereo projection of every exterior ring
    recs = []
    for geom in gm.geometry:
        if geom is None or geom.is_empty:
            continue
        polys = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
        for pg in polys:
            xy = np.asarray(pg.exterior.coords)
            u, v = to_stereo(xy[:, 0], xy[:, 1])
            if np.isfinite(u).all() and np.isfinite(v).all() and len(u) > 3:
                recs.append(Polygon(np.column_stack([u, v])))
    gs = gpd.GeoDataFrame(geometry=recs)
    gs.to_file(mst)
    print(f"[ex7] stereo shp ({len(recs)}) +{time.time()-t0:.0f}s",
          flush=True)

bbox = (-180.0, 180.0, -89.0, 90.0)
reg = Region(bbox, 4326)
sh = Shoreline(str(mll), reg.bbox, 4e3*DEG)
sdf = om.signed_distance_function(sh)
dem = DEM(str(DS/"SRTM15+.nc"), bbox=reg)
print(f"[ex7] shoreline+dem +{time.time()-t0:.0f}s", flush=True)
f = om.feature_sizing_function(sh, sdf, r=3,
                               max_edge_length=20e3*DEG)
w = om.wavelength_sizing_function(dem, wl=30,
                                  min_edgelength=4e3*DEG,
                                  max_edge_length=20e3*DEG,
                                  grid_dx=4e3*DEG)
# Example_7 passes 'slp' WITHOUT 'fl'; edgefx.m defval=0 => the
# slope filter is OFF (gradient of the raw decimated bathymetry).
# The barotropic/50 filter we first used smoothed the 10-19 km
# transition band ~12% coarser than the golden field.
s = om.bathymetric_gradient_sizing_function(
    dem, slope_parameter=10,
    min_edge_length=4e3*DEG, max_edge_length=20e3*DEG,
    type_of_filter="none", grid_dx=4e3*DEG)
g, dta = om.finalize_sizing([f, w, s], dem=dem, shoreline=sh,
                            hmin=4e3, max_edge_length=20e3,
                            gradation=0.25,
                            courant={"timestep": 0.0})
print(f"[ex7] sizing done dt={dta:.2f} +{time.time()-t0:.0f}s",
      flush=True)
# NOTE: no extra stereo re-gradation — edgefx.m grades ONCE on the
# lonlat lattice (with the cyclical x-padding now inside finalize's
# limgradStruct); the fork's "reinject stereographic gradient" pass
# is not part of OM2D and added ~15% nodes (diag 6180327).

sh2 = Shoreline(str(mst), reg.bbox, 4e3*DEG, stereo=True)
dom = om.signed_distance_function(sh2)
p, t = om.generate_mesh(dom, g, stereo=True, max_iter=100, seed=0)
print(f"[ex7] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p_stereo.npy", p); np.save(OUT/"t.npy", t)
from oceanmesh.region import to_lat_lon
lon, lat = to_lat_lon(p[:, 0], p[:, 1])
np.save(OUT/"p.npy", np.column_stack([lon, lat]))
print("[ex7] saved", flush=True)
