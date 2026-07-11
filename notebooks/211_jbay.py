# Ladder step 4: 1:1 translation of Example_5b_JBAY_w_weirs.m
# (2 weirs as faux islands + pfix/egfix constraints, h0=15 m).
import os, sys, logging, time
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
from scipy.io import loadmat
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from oceanmesh.weirs import build_weir_geometry, match_weir_nodes

DSET = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/PostSandyNCEI"))
WST = Path(os.path.expanduser(
    "~/Github/OceanMesh2D/datasets/weirs_struct.mat"))
OUT = Path("outputs/om2d_examples/jbay")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
bbox = (-73.97, -73.75, 40.5, 40.68)
min_el, max_el, dt, grade, R = 15.0, 1e3, 2.0, 0.15, 3

t0 = time.time()
w = loadmat(str(WST), squeeze_me=True)["weirs"]
pfix_all, egfix_all, pairs_all, rings = [], [], [], []
off = 0
for wi in np.atleast_1d(w):
    cl = np.column_stack([np.atleast_1d(wi["X"]),
                          np.atleast_1d(wi["Y"])])
    pf, ef, ib = build_weir_geometry(cl, float(wi["width"]),
                                     float(wi["min_ele"]))
    pfix_all.append(pf)
    egfix_all.append(ef + off)
    off += len(pf)
    pairs_all.append(ib)
    rings.append(np.vstack([pf, pf[0], [[np.nan, np.nan]]]))
pfix = np.vstack(pfix_all)
egfix = np.vstack(egfix_all)
print(f"[jbay] weirs: pfix {len(pfix)} egfix {len(egfix)}",
      flush=True)

reg = Region(bbox, 4326)
shore = Shoreline(str(DSET / "PostSandyNCEI.shp"), reg.bbox,
                  min_el * DEG)
# geodata.m:293-303: append weir faux islands to inner AFTER
# shoreline processing (no densify/smooth on the weir rings)
shore.inner = np.vstack(
    [np.asarray(shore.inner).reshape(-1, 2),
     np.array([[np.nan, np.nan]])] + rings
)
sdf = om.signed_distance_function(shore)
dem = DEM(str(DSET / "PostSandyNCEI.nc"), bbox=reg)
feat = om.feature_sizing_function(shore, sdf, r=R,
                                  max_edge_length=max_el * DEG)
weir_cfg = [{"crestline": np.column_stack(
                 [np.atleast_1d(wi["X"]), np.atleast_1d(wi["Y"])]),
             "min_ele_m": float(wi["min_ele"]),
             "width_m": float(wi["width"])}
            for wi in np.atleast_1d(w)]
grid, _ = om.finalize_sizing(
    [feat], dem=dem, shoreline=shore, hmin=min_el,
    max_edge_length=max_el, gradation=grade,
    courant={"timestep": dt}, weirs=weir_cfg,
)
print(f"[jbay] sizing +{time.time()-t0:.0f}s", flush=True)
import sys as _sys
_seed = int(os.environ.get("JBAY_SEED", "0"))
_raw_only = os.environ.get("JBAY_RAW_ONLY") == "1"
if _raw_only:
    p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=_seed,
                            pfix=pfix, egfix=egfix,
                            cleanup="none")
    OUT_ = Path("outputs/om2d_examples/jbay")
    OUT_.mkdir(parents=True, exist_ok=True)
    np.save(OUT_ / f"p_raw_fresh_s{_seed}.npy", p)
    np.save(OUT_ / f"t_raw_fresh_s{_seed}.npy", t)
    with open(OUT_ / f"raw_fresh_s{_seed}.14", "w") as f:
        f.write(f"raw fresh seed {_seed}\n{len(t)} {len(p)}\n")
        for i, (x, y) in enumerate(p, 1):
            f.write(f"{i} {x:.10f} {y:.10f} 0.0\n")
        for i, (a_, b_, c_) in enumerate(t + 1, 1):
            f.write(f"{i} 3 {a_} {b_} {c_}\n")
        f.write("0\n0\n0\n0\n")
    print(f"[jbay] raw saved NP={len(p):,}", flush=True)
    raise SystemExit(0)
p, t = om.generate_mesh(sdf, grid, max_iter=100, seed=_seed,
                        pfix=pfix, egfix=egfix)
print(f"[jbay] seed={_seed}", flush=True)
print(f"[jbay] mesh NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
b = om.interp_bathymetry(p, t, dem, method="cell-averaging",
                         min_depth=1.0, nan_fill=True)
np.save(OUT / "p.npy", p); np.save(OUT / "t.npy", t)
np.save(OUT / "b.npy", b)

# TestJBAY criteria
Re2 = 111.0 ** 2
X, Y = p[:, 0], p[:, 1]
x1, y1 = X[t[:, 0]], Y[t[:, 0]]
x2, y2 = X[t[:, 1]], Y[t[:, 1]]
x3, y3 = X[t[:, 2]], Y[t[:, 2]]
pa = 0.5 * np.abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))
cosf = np.cos(np.deg2rad((y1 + y2 + y3) / 3))
area = float((pa * cosf).sum() * Re2)
bc = b[t].mean(axis=1) / 1e3
vol = float((pa * cosf * bc).sum() * Re2)
from pyproj import Transformer
_tr = Transformer.from_crs(
    "EPSG:4326", "+proj=tmerc +lon_0=-73.86 +lat_0=40.59 "
    "+ellps=WGS84 +units=m", always_xy=True)
xx, yy = _tr.transform(p[:, 0], p[:, 1])
from oceanmesh.fix_mesh import simp_qual
q = simp_qual(np.column_stack([xx, yy]), t)
# weir pairing (bd.nbou==2 analog)
n_ok = 0
for ib in pairs_all:
    try:
        match_weir_nodes(p, ib)
        n_ok += 1
    except Exception as e:
        print(f"[jbay] weir pairing failed: {e}", flush=True)
print(f"[jbay] NP={len(p):,} (target 25,190 +-5%) "
      f"NT={len(t):,} (45,270 +-5%)", flush=True)
print(f"[jbay] area={area:.1f} km2 (211 +-1%) "
      f"vol={vol:.3f} km3 (2.075 +-1%)", flush=True)
print(f"[jbay] qual min={q.min():.4f} (>=0.25) "
      f"weirs paired={n_ok}/2", flush=True)
okNP = abs(len(p) - 25190) / 25190 <= 0.05
okNT = abs(len(t) - 45270) / 45270 <= 0.05
okA = abs(area - 211) / 211 <= 0.01
okV = abs(vol - 2.075) / 2.075 <= 0.01
okQ = q.min() >= 0.25
print(f"[jbay] {'PASSED' if all([okNP, okNT, okA, okV, okQ, n_ok==2]) else 'FAILED'} "
      f"(NP:{okNP} NT:{okNT} A:{okA} V:{okV} Q:{okQ})", flush=True)
