# Ladder step 6: FULL Example_3_ECGC (2 nests, real-data
# single-loop multiscale).
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
DSET = OM2D / "datasets/PostSandyNCEI"
OUT = Path("outputs/om2d_examples/ex3full")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

# ---- nest 1: ECGC quad, GSHHS + SRTM15+
bbox1 = np.array([[-71.6, 42.7], [-64, 30], [-80, 24],
                  [-85, 38], [-71.6, 42.7]])
reg1 = Region((float(bbox1[:,0].min()), float(bbox1[:,0].max()),
               float(bbox1[:,1].min()), float(bbox1[:,1].max())),
              4326)
sh1 = Shoreline(str(OM2D/"datasets/GSHHS_shp/f/GSHHS_f_L1.shp"),
                bbox1, 1e3*DEG)
sdf1 = om.signed_distance_function(sh1)
dem1 = DEM(str(OM2D/"datasets/SRTM15+.nc"), bbox=reg1)
f1 = om.feature_sizing_function(sh1, sdf1, r=3,
                                max_edge_length=50e3*DEG)
w1 = om.wavelength_sizing_function(dem1, wl=30,
                                   min_edgelength=1e3*DEG,
                                   max_edge_length=50e3*DEG)
g1, dt1 = om.finalize_sizing(
    [f1, w1], dem=dem1, shoreline=sh1, hmin=1e3,
    max_edge_length=50e3, gradation=0.35,
    courant={"timestep": 0.0},
)
print(f"[ex3f] nest1 sizing dt_auto={dt1:.2f}s "
      f"+{time.time()-t0:.0f}s", flush=True)

# ---- nest 2: NY polygon, PostSandy
bbox2 = np.array([[-74.25, 40.5], [-73.75, 40.55],
                  [-73.75, 41], [-74, 41], [-74.25, 40.5]])
reg2 = Region((-74.25, -73.75, 40.5, 41.0), 4326)
sh2 = Shoreline(str(DSET/"PostSandyNCEI.shp"), bbox2, 30.0*DEG)
sdf2 = om.signed_distance_function(sh2)
dem2 = DEM(str(DSET/"PostSandyNCEI.nc"), bbox=reg2)
f2 = om.feature_sizing_function(sh2, sdf2, r=3,
                                max_edge_length=1e3*DEG)
w2 = om.wavelength_sizing_function(dem2, wl=30,
                                   min_edgelength=30*DEG,
                                   max_edge_length=1e3*DEG)
g2, dt2 = om.finalize_sizing(
    [f2, w2], dem=dem2, shoreline=sh2, hmin=30.0,
    max_edge_length=1e3, max_edge_length_nearshore=240.0,
    gradation=0.35, courant={"timestep": 0.0},
)
print(f"[ex3f] nest2 sizing dt_auto={dt2:.2f}s "
      f"+{time.time()-t0:.0f}s", flush=True)

p, t = om.generate_multiscale_mesh(
    [sdf1, sdf2], [g1, g2], max_iter=100, seed=0,
    gradation=0.35,
)
print(f"[ex3f] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p.npy", p); np.save(OUT/"t.npy", t)

# combined interp: SRTM everywhere, PostSandy overwrites inside
# nest 2 (OM2D interp(m,{gdat1 gdat2}) sequential semantics)
b = om.interp_bathymetry(p, t, dem1, method="cell-averaging",
                         min_depth=1.0, nan_fill=True)
from oceanmesh import edges as om_edges
from oceanmesh.geometry import inpoly2
ring2 = np.vstack([bbox2, [[np.nan, np.nan]]])
e2 = om_edges.get_poly_edges(ring2)
m2, _ = inpoly2(p, np.nan_to_num(ring2), e2)
if m2.any():
    b2 = om.interp_bathymetry(p[m2], t[:0], dem2, method="nearest",
                              min_depth=1.0, nan_fill=False)
print(f"[ex3f] interp done +{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"b.npy", b)
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=-74.5 +lat_0=33 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[ex3f] qual min={q.min():.4f} mean={q.mean():.4f} saved",
      flush=True)
