# Ladder step 2b: port of OM2D Tests/TestInterp.m — mesh the same
# Gulf coast tile against 3 CUDEM lattices (equal/unequal/varying
# spacing), CA-interp bathymetry, and require the 3 mesh volumes
# to agree within 0.5% relative.
import os, sys, logging
import numpy as np
logging.basicConfig(level=logging.WARNING)
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DEG = 1.0 / 111e3
min_el, max_el, grade = 40.0, 500.0, 0.25
dis = grade
DEMS = ["CUDEM_equal.nc", "CUDEM_unequal.nc", "CUDEM_varying.nc"]

volumes = []
for name in DEMS:
    f = OM2D / "datasets/CUDEMS" / name
    import xarray as xr
    ds = xr.open_dataset(f)
    lonv = [v for v in ds.coords if "lon" in v.lower() or v == "x"][0]
    latv = [v for v in ds.coords if "lat" in v.lower() or v == "y"][0]
    bbox = (float(ds[lonv].min()), float(ds[lonv].max()),
            float(ds[latv].min()), float(ds[latv].max()))
    reg = Region(bbox, 4326)
    shore = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                      bbox, min_el*DEG)
    sdf = om.signed_distance_function(shore)
    grid0 = om.distance_sizing_function(shore, rate=dis,
                                        max_edge_length=max_el*DEG)
    grid, _ = om.finalize_sizing([grid0], shoreline=shore,
                                 hmin=min_el,
                                 max_edge_length=max_el,
                                 gradation=grade)
    p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=0)
    dem = DEM(str(f), bbox=reg)
    b = om.interp_bathymetry(p, t, dem, method="cell-averaging")
    a2 = p[t]
    areas = 0.5*np.abs(
        (a2[:,1,0]-a2[:,0,0])*(a2[:,2,1]-a2[:,0,1])
        - (a2[:,2,0]-a2[:,0,0])*(a2[:,1,1]-a2[:,0,1]))
    bc = b[t].mean(axis=1)
    vol = float((areas*bc).sum())
    volumes.append(vol)
    print(f"[interp] {name}: NP={len(p):,} volume={vol:.6e}",
          flush=True)
volumes = np.array(volumes)
vd = (volumes.max()-volumes.min())/volumes.mean()
ok = abs(vd) <= 0.005
print(f"[interp] max relative volume difference: {vd:.5f} "
      f"(<=0.005) {'PASSED' if ok else 'FAILED'}", flush=True)
