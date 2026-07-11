# Ladder step 3: 1:1 translation of Examples/Example_2_NY.m
# (PostSandyNCEI shoreline + DEM, h0=30 m, max_el_ns=240 m,
# max_el=1 km, dt=2 s CFL bound, grade 0.20, R=3, utm).
import os, sys, logging, time
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
DEG = 1.0 / 111e3
min_el, max_el, max_el_ns = 30.0, 1e3, 240.0
dt, grade, R = 2.0, 0.20, 3

t0 = time.time()
# geodata without bbox: extents come from the DEM (geodata.m)
dem_probe = DEM(str(DSET / "PostSandyNCEI.nc"))
x0, x1, y0, y1 = dem_probe.bbox
reg = Region((x0, x1, y0, y1), 4326)
print(f"[ex2] dem bbox {dem_probe.bbox} +{time.time()-t0:.0f}s",
      flush=True)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"),
                  reg.bbox, min_el * DEG)
sdf = om.signed_distance_function(shore)
dem = DEM(str(DSET / "PostSandyNCEI.nc"), bbox=reg)
feat = om.feature_sizing_function(shore, sdf, r=R,
                                  max_edge_length=max_el * DEG)
grid, dt_used = om.finalize_sizing(
    [feat], dem=dem, shoreline=shore, hmin=min_el,
    max_edge_length=max_el,
    max_edge_length_nearshore=max_el_ns,
    gradation=grade,
    courant={"timestep": dt},
)
print(f"[ex2] sizing +{time.time()-t0:.0f}s", flush=True)
p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=0)
print(f"[ex2] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
b = om.interp_bathymetry(p, t, dem, method="cell-averaging",
                         min_depth=1.0)
OUT = Path("outputs/om2d_examples/ex2_ny")
OUT.mkdir(parents=True, exist_ok=True)
np.save(OUT / "p.npy", p)
np.save(OUT / "t.npy", t)
np.save(OUT / "b.npy", b)
from oceanmesh.fix_mesh import simp_qual
q = simp_qual(p, t)
print(f"[ex2] qual min={q.min():.4f} mean={q.mean():.4f} saved "
      f"+{time.time()-t0:.0f}s", flush=True)
