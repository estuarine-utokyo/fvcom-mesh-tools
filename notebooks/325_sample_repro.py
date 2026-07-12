# Reproduce the goto2023 SAMPLE (TokyoBay_grd.dat) with the
# faithful OM2D stack: bay-only domain, the sample's OBC arc as a
# CONSTRAINED input line (pfix+egfix), uniform-Courant sizing
# (h = sqrt(gH)*dt/Cr clipped to [hmin, maxel]), OSM coastline.
import os, sys, logging, time
import faulthandler
faulthandler.enable()
faulthandler.dump_traceback_later(600, repeat=True)
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
OUT = Path("outputs/sample_repro")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3
t0 = time.time()

# knobs (single-pass targeted tuning against the sample bands)
H0 = float(os.environ.get("SR_H0", 290.0))          # m
MAXEL = float(os.environ.get("SR_MAXEL", 1400.0))    # m
GRADE = float(os.environ.get("SR_GRADE", 0.165))
DT = float(os.environ.get("SR_DT", 15.0))            # s
CRMIN = float(os.environ.get("SR_CRMIN", 0.45))
FS = float(os.environ.get("SR_FS", 3.0))

# The OBC is an INPUT: the goto2023 sample's 13-node smooth arc
# (TokyoBay_obc.dat + TokyoBay_grd.dat, EPSG:32654 -> 4326),
# NW (Miura/Otsu coast) -> SE (mid-channel corner). Mesh nodes are
# CONSTRAINED onto this line via pfix+egfix.
OBC_ARC = np.array([
    [139.6713, 35.1396], [139.6737, 35.1288], [139.6772, 35.1168],
    [139.6816, 35.1031], [139.6871, 35.0877], [139.6946, 35.0705],
    [139.7000, 35.0576], [139.7069, 35.0445], [139.7134, 35.0327],
    [139.7216, 35.0184], [139.7289, 35.0047], [139.7373, 34.9916],
    [139.7497, 34.9750]])
OBC_SEG = np.column_stack([np.arange(12), np.arange(1, 13)])

# bay-only domain POLYGON matching the goto2023 sample geometry:
# the SW crossing follows the OBC arc, extended into Miura land at
# the NW end; the southern closure runs east at lat~34.973 to the
# Boso coast (the sample carries this as an artificial land line).
poly = np.array(
    [[139.83, 34.973], [140.12, 34.973], [140.12, 35.75],
     [139.60, 35.75], [139.60, 35.20], [139.6642, 35.1546]]
    + OBC_ARC.tolist()
    + [[139.83, 34.973]])
bbox = (139.60, 140.12, 34.96, 35.75)
reg = Region(bbox, 4326)

# geometry-stage narrow-channel policy (owner 2026-07-12): decide
# channel fates on the SHORELINE, before meshing -- deterministic
# and minimum-mesh-size preserving. Through / big-basin channels
# get their banks pushed into land (width -> 1.8x h_mesh) so two
# STANDARD-size rows fit; dead-ends and small basins are closed.
import geopandas as gpd
from shapely.geometry import Polygon as _Poly
from shapely.ops import unary_union as _uu
from fvcom_mesh_tools.prep.channel_policy_geom import (
    apply_channel_policy_to_land,
)

CH_SHP = OUT / "land_channel_adj.shp"
_land_g = gpd.read_file("outputs/tb_varres_3r/land_osm_wide.shp")
_dom = _Poly(poly)
_cosw = float(np.cos(np.deg2rad(35.35)))
_new_land, chinfo = apply_channel_policy_to_land(
    _uu(list(_land_g.geometry)), _dom,
    h_mesh_m=H0 * 1.2, obc_point=tuple(OBC_ARC[6]),
    metric_scale=(111e3 * _cosw, 111e3))
print(f"[sr] channel policy (geometry stage): "
      f"{len(chinfo['widened'])} widened, "
      f"{len(chinfo['closed'])} closed "
      f"(of {chinfo['n_narrow']} narrow corridors)", flush=True)
for r in chinfo["widened"] + chinfo["closed"]:
    print(f"[sr]   {r['action']}: ({r['center'][0]:.3f}, "
          f"{r['center'][1]:.3f}) area={r['area_cells']:.1f} "
          f"extent={r['extent_cells']:.1f} cells "
          f"members={r['n_members']}", flush=True)
_geoms = list(_new_land.geoms) if hasattr(_new_land, "geoms")     else [_new_land]
gpd.GeoDataFrame(geometry=_geoms, crs=_land_g.crs).to_file(CH_SHP)

sh = Shoreline(str(CH_SHP), poly, H0 * DEG)
sdf = om.signed_distance_function(sh)
dem = DEM(str(OM2D / "datasets/TokyoBay/dem/SRTM15_kanto_15s.nc"),
          bbox=reg, nc_reader="coords")
print(f"[sr] inputs +{time.time()-t0:.0f}s", flush=True)

f = om.feature_sizing_function(
    sh, sdf, r=FS, max_edge_length=MAXEL * DEG,
    lattice_anchor=(dem.bbox[0], dem.bbox[2]))
# SR_MODE=courant: the sample's measured recipe (anatomy 2026-07-11,
# rev.3 after the E/W coast split):
#   The sample has a TILTED, UNIFORM-BASE mouth zone -- coastal
#   base ~930 (mesh) everywhere south of the line through
#   (139.695, 35.215)-(139.795, 35.155)  [W coast coarse to 35.215,
#   E coast only to ~35.155; same latitude, different sides ->
#   not a latitude function], and ~447 (mesh) in the bay. Inside
#   each zone:  h = min(base + g*d_coast, cap),
#   cap = 1030 (bay) / MAXEL (mouth);  limgrad(g) then smooths the
#   zone seam over ~3 km exactly like the sample's 437-665
#   transition bands. Field = mesh/1.2 (DistMesh calibration).
ZBASE = float(os.environ.get("SR_ZBASE", 630.0))   # mouth-zone base, m
ZW = float(os.environ.get("SR_ZW_LAT", 35.215))    # zone lat at 139.695
ZE = float(os.environ.get("SR_ZE_LAT", 35.17))     # zone lat at 139.795
DMAX = float(os.environ.get("SR_DMAX", 1030.0))    # bay far-field cap, m
if os.environ.get("SR_MODE", "courant") == "courant":
    fdst = om.distance_sizing_function(sh, rate=GRADE,
                                       max_edge_length=None)
    vals = np.ma.filled(np.ma.asarray(fdst.values), MAXEL * DEG)
    d_m = (vals / DEG - H0) / GRADE            # metric coast distance
    lon_g, lat_g = fdst.create_grid()
    lat_line = ZW + (ZE - ZW) / 0.10 * (lon_g - 139.695)
    in_zone = lat_g < lat_line
    base_m = np.where(in_zone, ZBASE, H0)
    cap_m = np.where(in_zone, MAXEL, DMAX)
    h_m = np.minimum(base_m + GRADE * d_m, cap_m)
    fdst.values = h_m * DEG
    fdst.build_interpolant()
    f = fdst
    print(f"[sr] mode=courant (two-zone distance field: bay {H0:.0f}"
          f"/{DMAX:.0f}, mouth {ZBASE:.0f}/{MAXEL:.0f}, seam "
          f"({139.695},{ZW})-({139.795},{ZE}), g={GRADE})", flush=True)
g, dt_eff = om.finalize_sizing(
    [f], dem=dem, shoreline=sh, hmin=H0,
    max_edge_length=MAXEL, gradation=GRADE, courant=None)
print(f"[sr] sizing done +{time.time()-t0:.0f}s", flush=True)

# OBC boundary-band (fvcom_mesh_tools.obc_band, default ON):
# smooth inner guide line at k*local-target-size offsets + a
# boundary-priority size corridor along the whole southern
# crossing, applied AFTER limgrad. General rule (sample-calibrated
# K=1.25); no sample-specific numbers remain.
if os.environ.get("SR_OBC_LADDER", "on") == "on":
    from fvcom_mesh_tools.obc_band import (
        apply_corridor, build_obc_band, corridor_targets)
    h_arc_m = np.asarray(g.eval(OBC_ARC)).ravel() / DEG * 1.2
    # end sizes: a field eval AT an arc end that abuts a coast or
    # an artificial closure is contaminated by the coastal halo
    # (TB SE corner: ~820 m where the boundary-adjacent water is
    # 1680 m class). No silent guessing -- override explicitly.
    if os.environ.get("SR_OBC_H0"):
        h_arc_m[0] = float(os.environ["SR_OBC_H0"])
    if os.environ.get("SR_OBC_H1"):
        h_arc_m[-1] = float(os.environ["SR_OBC_H1"])
    band = build_obc_band(OBC_ARC, h_arc_m, k_offset=1.25,
                          skip_ends=1)
    PFIX, SEGS = band["pfix"], band["egfix"]
    print(f"[sr] OBC band: offsets "
          f"{band['offsets_m'].min():.0f}-"
          f"{band['offsets_m'].max():.0f} m (K=1.25 x local size)",
          flush=True)
    closure = np.array([[139.7497, 34.9750], [139.83, 34.973]])
    h_coast_m = float(np.asarray(
        g.eval(np.array([[139.82, 34.974]]))).ravel()[0]) / DEG * 1.2
    # corridor target = the BAND size (K x local), not the ambient
    # local size -- the boundary band must stay one class coarser
    pts_m, tgt_m = corridor_targets(
        OBC_ARC, band["offsets_m"], closure_ll=closure,
        h_closure_end_m=h_coast_m)
    lon_g2, lat_g2 = g.create_grid()
    g.values, n_up = apply_corridor(
        lon_g2, lat_g2, np.asarray(g.values, dtype=float),
        pts_m, tgt_m, grade=GRADE,
        arc_mean_lat=float(OBC_ARC[:, 1].mean()))
    g.build_interpolant()
    print(f"[sr] boundary corridor: raised {n_up} lattice cells "
          f"(post-limgrad, boundary-priority)", flush=True)
else:
    PFIX, SEGS = OBC_ARC, OBC_SEG

# the built-in msh.clean('default') runs with pfix nodes pinned and
# (fork feature) egfix-carrying faces excluded from the boundary
# deletion loop, so the constrained OBC line survives the clean
from scipy.spatial import cKDTree
import shapely
from fvcom_mesh_tools.algorithms.obc_finish import (
    prune_one_wide_protected,
)

arc_line = shapely.LineString(OBC_ARC)
TOL = 1e-6  # deg (~0.1 m): constrained nodes are exact


def _loop_area(pts, lp):
    x, y = pts[lp, 0], pts[lp, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1))
                       - np.dot(y, np.roll(x, -1)))


def build_mesh(g):
    """One full generate->clean->bc pass on the current sizing."""
    p, t = om.generate_mesh(sdf, g, max_iter=60, seed=0,
                            pfix=PFIX, egfix=SEGS)
    ne0 = len(t)
    p, t = prune_one_wide_protected(p, t, PFIX)
    p, t = om.make_mesh_boundaries_traversable(p, t)
    print(f"[sr] 1-wide pruning (pfix-protected): NE {ne0:,} -> "
          f"{len(t):,}", flush=True)
    print(f"[sr] mesh NP={len(p):,} NE={len(t):,} "
          f"+{time.time()-t0:.0f}s", flush=True)
    d_arc, arc_idx = cKDTree(p).query(OBC_ARC)
    if (d_arc > 1e-8).any():
        bad = np.where(d_arc > 1e-8)[0]
        raise RuntimeError(
            f"OBC arc nodes lost during cleanup: arc rows "
            f"{bad.tolist()} (offsets {d_arc[bad]} deg).")
    b = om.interp_bathymetry(p, t, dem, method="cell-averaging",
                             min_depth=2.0)
    loops = om.boundary_loops(t)
    outer = max(loops, key=lambda lp: abs(_loop_area(p, lp)))
    arc_set = set(arc_idx.tolist())
    if not arc_set.issubset(set(outer.tolist())):
        raise RuntimeError("OBC arc nodes are not all on the outer "
                           "boundary loop -- inspect the mesh.")
    start = int(np.where(outer == arc_idx[0])[0][0])
    ring = np.roll(outer, -start)
    if shapely.distance(shapely.Point(p[ring[1]]), arc_line) > TOL:
        ring = np.roll(ring[::-1], 1)
    stop = int(np.where(ring == arc_idx[-1])[0][0])
    open_str = ring[:stop + 1]
    off = np.array([shapely.distance(shapely.Point(p[v]), arc_line)
                    for v in open_str])
    if (off > TOL).any() or not arc_set.issubset(
            set(open_str.tolist())):
        raise RuntimeError(
            "outer-loop run between the OBC arc ends leaves the "
            f"constrained line (max offset {off.max():.2e} deg).")
    land_str = np.append(ring[stop:], ring[0])
    bc = {"open": [open_str], "land": [land_str],
          "island": [lp for lp in loops if lp is not outer]}
    print(f"[sr] bc: open={len(bc['open'])} land={len(bc['land'])} "
          f"island={len(bc['island'])}", flush=True)
    for kk, s2 in enumerate(bc["open"]):
        print(f"[sr]  open[{kk}]: {len(s2)} nodes at "
              f"({p[s2, 0].min():.3f}-{p[s2, 0].max():.3f}, "
              f"{p[s2, 1].min():.3f}-{p[s2, 1].max():.3f})",
              flush=True)
    return p, t, b, bc


p, t, b, bc = build_mesh(g)

om.write_fort14(str(OUT / "sample_repro.14"), p, t, depth=b,
                boundaries=bc)
np.save(OUT / "p.npy", p); np.save(OUT / "t.npy", t)
# kickoff: FVCOM production mesh is cartesian UTM54N (EPSG:32654)
from pyproj import Transformer as _T
_tr = _T.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
_xu, _yu = _tr.transform(p[:, 0], p[:, 1])
om.write_fort14(str(OUT / "sample_repro_utm.14"),
                np.column_stack([_xu, _yu]), t, depth=b,
                boundaries=bc)
print("[sr] wrote UTM54N fort.14", flush=True)

# per-node Courant at DT (msh.CalcCFL port) -- the design target is
# Cr <= 0.5 everywhere with dt = DT
cr = om.calc_cfl(p, t, b, dt=DT)
print(f"[sr] Cr(dt={DT}) p50/p90/p99/max = "
      f"{np.percentile(cr, [50, 90, 99]).round(3)} {cr.max():.3f} "
      f"(n>0.5: {(cr > 0.5).sum()}, n>0.6: {(cr > 0.6).sum()})",
      flush=True)

from oceanmesh.fix_mesh import simp_qual
from pyproj import Transformer
tr = Transformer.from_crs("EPSG:4326", "EPSG:32654", always_xy=True)
xx, yy = tr.transform(p[:, 0], p[:, 1])
q = simp_qual(np.column_stack([xx, yy]), t)
print(f"[sr] qual min={q.min():.3f} mean={q.mean():.3f}", flush=True)

# band table vs the sample
e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
e = np.unique(np.sort(e, axis=1), axis=0)
lat = 0.5 * (p[e[:, 0], 1] + p[e[:, 1], 1])
L = np.hypot((p[e[:, 0], 0] - p[e[:, 1], 0])
             * np.cos(np.deg2rad(lat)),
             p[e[:, 0], 1] - p[e[:, 1], 1]) / DEG
print("[sr] bands (sample: 1130/1146/541/547/534/518):", flush=True)
for lo, hi in [(34.95, 35.10), (35.10, 35.20), (35.20, 35.30),
               (35.30, 35.40), (35.40, 35.55), (35.55, 35.70)]:
    s = (lat >= lo) & (lat < hi)
    if s.sum():
        print(f"[sr]  {lo:.2f}-{hi:.2f}: p10/50/90 = "
              f"{np.percentile(L[s], [10, 50, 90]).round(0)}",
              flush=True)
