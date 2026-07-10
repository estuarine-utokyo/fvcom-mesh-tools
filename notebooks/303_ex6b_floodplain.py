# Ladder step 14: Example_6b_GBAY_w_floodplain (two-stage build:
# underwater mesh -> makens auto -> extractFixedConstraints ->
# overland mesh against the 10m-LMSL contour with pfix/egfix).
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
OUT = Path("outputs/om2d_examples/ex6bfp")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

mm = loadmat(str(DS/"ECGC_Thalwegs.mat"), squeeze_me=False)
channels = [np.asarray(a, dtype=float)
            for a in mm["pts2"].ravel() if np.size(a) >= 4]

# --- STAGE 1: underwater mesh (Example_6 config, banded max_el/grade) ---
bbox = (-95.40, -94.4, 29.14, 30.09)
reg = Region(bbox, 4326)
sh = Shoreline(
    str(DS/"US_Medium_Shoreline/us_medium_shoreline_polygon.shp"),
    reg.bbox, 60.0*DEG)
sh.detect_inpoly_flip(str(DS/"GSHHS_shp/l/GSHHS_l_L1.shp"))
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
max_el = np.array([[1e3, -np.inf, 0.0], [500.0, 0.0, np.inf]])
grade = np.array([[0.25, -np.inf, 0.0], [0.05, 0.0, np.inf]])
g, _ = om.finalize_sizing([f, c], dem=dem, shoreline=sh,
                          hmin=60.0, max_edge_length=max_el,
                          gradation=grade)
print(f"[ex6b] sizing +{time.time()-t0:.0f}s", flush=True)
puw, tuw = om.generate_mesh(sdf, g, max_iter=100, seed=0)
print(f"[ex6b] stage1 NP={len(puw):,} NT={len(tuw):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"puw.npy", puw); np.save(OUT/"tuw.npy", tuw)

# --- makens auto + extractFixedConstraints ---
# depth for classification comes from gdat's DEM interpolant
# (msh.m: eb_depth = gdat.Fb(eb_mid)); clamp queries into the DEM
# bbox to emulate griddedInterpolant's nearest extrapolation
eps = 1e-9
qx = np.clip(puw[:, 0], dem.bbox[0] + eps, dem.bbox[1] - eps)
qy = np.clip(puw[:, 1], dem.bbox[2] + eps, dem.bbox[3] - eps)
dep = -np.asarray(dem.eval((qx, qy)), dtype=float)  # positive-down
bc = om.make_bc_auto(puw, tuw, depth=dep, shoreline=sh,
                     classifier="both")
pfix, egfix = om.extract_fixed_constraints(puw, tuw, bc["open"])
print(f"[ex6b] constraints pfix={len(pfix):,} egfix={len(egfix):,} "
      f"open_strings={len(bc['open'])} +{time.time()-t0:.0f}s",
      flush=True)
np.save(OUT/"pfix.npy", pfix); np.save(OUT/"egfix.npy", egfix)

# --- STAGE 2: overland mesh (10m-LMSL contour), same fh ---
sh2 = Shoreline(
    str(DS/"CoastalReliefModel_10m_LMSL/us_coastalreliefmodel_10mLMSL.shp"),
    reg.bbox, 60.0*DEG)
# MATLAB flips this dataset (golden log: "Shapefile inpoly is
# inconsistent with GHSSS test file"); our raw-segment parity vote
# lands at 43/100 (borderline, threshold 50) because the .m votes
# with the BFS-closed obj.outer which this port does not build.
# Golden stage-2 mesh confirms the domain: bay+gulf+floodplain IN,
# >10m inland OUT.
sh2.inpoly_flip = True
sdf2 = om.signed_distance_function(sh2)
print(f"[ex6b] stage2 shoreline +{time.time()-t0:.0f}s", flush=True)
p2, t2 = om.generate_mesh(sdf2, g, max_iter=100, seed=0,
                          pfix=pfix, egfix=egfix)
print(f"[ex6b] stage2 NP={len(p2):,} NT={len(t2):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT/"p.npy", p2); np.save(OUT/"t.npy", t2)

from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326",
    "+proj=tmerc +lon_0=-94.9 +lat_0=29.6 +ellps=WGS84 +units=m",
    always_xy=True)
for tag, pp, tt in [("stage1", puw, tuw), ("stage2", p2, t2)]:
    xx, yy = tr.transform(pp[:, 0], pp[:, 1])
    q = simp_qual(np.column_stack([xx, yy]), tt)
    print(f"[ex6b] {tag} qual min={q.min():.4f} mean={q.mean():.4f}",
          flush=True)
