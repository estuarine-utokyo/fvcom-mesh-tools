# Ladder step 7: Example_4_PRVI (4 nests, fs=-5 auto, slp, Rossby
# filter fl=-50, dt=0, itmax=50).
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
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
OUT = Path("outputs/om2d_examples/ex4prvi")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
WL, GRADE, R, SLP, FLQ = 30, 0.25, -5, 15, 50
t0 = time.time()

def build_nest(shp, dem_path, h0_m, bbox=None):
    if bbox is None:
        dp = DEM(str(dem_path))
        bbox = dp.bbox
    reg = Region(tuple(bbox) if not isinstance(bbox, np.ndarray)
                 else (float(bbox[:,0].min()), float(bbox[:,0].max()),
                       float(bbox[:,1].min()), float(bbox[:,1].max())),
                 4326)
    bb_for_shore = bbox if isinstance(bbox, np.ndarray) else reg.bbox
    sh = Shoreline(str(shp), bb_for_shore, h0_m * DEG)
    sdf = om.signed_distance_function(sh)
    dem = DEM(str(dem_path), bbox=reg)
    f = om.feature_sizing_function(sh, sdf, r=R,
                                   max_edge_length=1000 * DEG)
    w = om.wavelength_sizing_function(dem, wl=WL,
                                      min_edgelength=h0_m * DEG,
                                      max_edge_length=1000 * DEG)
    s = om.bathymetric_gradient_sizing_function(
        dem, slope_parameter=SLP, filter_quotient=FLQ,
        min_edge_length=h0_m * DEG,
        max_edge_length=1000 * DEG,
        type_of_filter="barotropic")
    g, dt_a = om.finalize_sizing(
        [f, w, s], dem=dem, shoreline=sh, hmin=h0_m,
        max_edge_length=1000.0, gradation=GRADE,
        courant={"timestep": 0.0})
    print(f"[ex4] nest {shp} dt_auto={dt_a:.2f}s "
          f"+{time.time()-t0:.0f}s", flush=True)
    return sdf, g

doms, grids = [], []
d1, g1 = build_nest(DS/"GSHHS_shp/f/GSHHS_f_L1.shp",
                    DS/"SRTM15+.nc", 1000.0,
                    bbox=(-100.0, -53.0, 5.0, 52.5))
doms.append(d1); grids.append(g1)
for shp, dem_p, h0 in (
    (DS/"PR_1arcsec/pr_1s_0m_contour.shp",
     DS/"PR_1arcsec/pr_1s.nc", 30.0),
    (DS/"USVI_1arcsec/usvi_0m_contour.shp",
     DS/"USVI_1arcsec/usvi_1_mhw_2014.nc", 30.0),
    (DS/"SanJuan_1-9arcsec/sj_0contour_closed.shp",
     DS/"SanJuan_1-9arcsec/san_juan_19_prvd02_2015.nc", 10.0),
):
    d_, g_ = build_nest(shp, dem_p, h0)
    doms.append(d_); grids.append(g_)

p, t = om.generate_multiscale_mesh(doms, grids, max_iter=50,
                                   seed=0, gradation=GRADE)
print(f"[ex4] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p.npy", p); np.save(OUT/"t.npy", t)
from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=-76.5 +lat_0=28 +ellps=WGS84 +units=m",
    always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[ex4] qual min={q.min():.4f} mean={q.mean():.4f} saved",
      flush=True)
