# STAGE-1 ACCEPTANCE: 1:1 Python translation of OceanMesh2D
# Tokyo_Bay/scripts/mesh_wide_varres_3r.m using the oceanmesh port.
# Same inputs (GSHHS_f_L1 / coastline_2 / SRTM15), same nest
# polygons, same edgefx parameters, same pipeline order
# (geodata -> edgefx -> meshgen.build -> interp -> make_bc).
# NO walls, NO cuts, NO morphological preprocessing.
import faulthandler
import logging
import os
import sys
from pathlib import Path

import numpy as np

faulthandler.enable()
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")

os.environ.setdefault("MPLBACKEND", "Agg")

import oceanmesh as om  # noqa: E402
from oceanmesh import DEM, Region, Shoreline  # noqa: E402

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
GSHHS = OM2D / "datasets/GSHHS_shp/f/GSHHS_f_L1.shp"
# user policy: the engineered OSM coastline for the bay/Kanto
# nests (coastline_2 was the .m's hand-made composite)
COAST2 = Path("outputs/tb_varres_3r/land_osm_wide.shp")
# sliced DEMs (the .m used 'SRTM15_V2.3_sliced.nc'; full-global
# SRTM15+ made nest-1 sizing a 170M-cell computation)
SRTM_PACIFIC = OM2D / "datasets/TokyoBay/dem/SRTM15_pacific_4min.nc"
# user policy (this stage): the same DEM family the .m used —
# SRTM15 (sliced); M7001 comes later with the production pipeline
SRTM_KANTO = OM2D / "datasets/TokyoBay/dem/SRTM15_kanto_15s.nc"
OUT = Path("outputs/tb_varres_3r")
OUT.mkdir(parents=True, exist_ok=True)
DEG = 1.0 / 111e3  # OM2D uses /111e3 throughout
MESH_ITMAX = 50

# ---- region polygons (verbatim from the .m file) --------------------
x1 = [157.4634361, 158.84262217, 160.046858, 161.07987092, 161.94757837,
      162.65742188, 163.21770287, 163.63704393, 163.92386584, 164.08625707,
      164.1323057, 164.06882521, 163.9009549, 163.63627547, 163.31281843,
      162.89825466, 162.38857843, 161.80134206, 161.13903484, 160.4052709,
      159.60336662, 158.7360147, 157.80598566, 156.81590455, 155.7687459,
      154.66746366, 153.5151905, 152.3154099, 151.0718658, 149.78866673,
      148.4702309, 147.12135316, 145.74715292, 144.35301788, 142.94462003,
      141.52778842, 140.10846886, 138.69263968, 137.2862313, 135.89504337,
      134.52470652, 133.18055169, 131.86766327, 130.59072034, 129.35412832,
      128.16190388, 127.01767503, 125.92498131, 124.88679107, 123.90615375,
      122.98581257, 122.12855282, 121.60387397, 115.64257669, 158.59463098,
      167.27851682, 157.4634361]
y1 = [51.36354136, 50.35540912, 49.27890879, 48.14482943, 46.96329793,
      45.74371019, 44.4946652, 43.2239559, 41.93859546, 40.64492778,
      39.34870289, 38.0551143, 36.76910533, 35.56167205, 34.36768078,
      33.12965026, 31.91772707, 30.7320848, 29.57667543, 28.45502315,
      27.37062929, 26.32729387, 25.32881846, 24.37921688, 23.48234609,
      22.64223232, 21.86295783, 21.14846084, 20.50262818, 19.92912188,
      19.4314301, 19.01266254, 18.67552551, 18.42237308, 18.25486888,
      18.17421142, 18.18096748, 18.27508101, 18.45587251, 18.7220402,
      19.07182143, 19.5027674, 20.01215433, 20.59669911, 21.25297733,
      21.97724182, 22.76552456, 23.61401919, 24.51851667, 25.47512814,
      26.47984282, 27.52881366, 28.22438211, 39.405273, 75.05163401,
      62.25276725, 51.36354136]
x2 = [138.196907945361, 139.784782123012, 141.294666567004,
      139.718998032166, 138.285893837899, 138.196907945361]
y2 = [34.6003740729489, 33.5207179317026, 35.2698037998702,
      35.9357551145740, 35.2841759386230, 34.6003740729489]
x3 = [139.142156, 139.142156, 139.144724, 139.152351, 139.164815,
      139.18178, 139.202839, 139.227552, 139.25548, 139.286203,
      139.31934, 139.354542, 139.391505, 139.429958, 139.469663,
      139.510408, 139.552007, 139.594291, 139.637105, 139.680308,
      139.723765, 139.767348, 139.810932, 139.854389, 139.897592,
      139.940406, 139.98269, 140.024289, 140.065034, 140.104739,
      140.143192, 140.180155, 140.215357, 140.248493, 140.279217,
      140.307145, 140.331858, 140.352917, 140.369882, 140.382346,
      140.389973, 140.392541, 140.392541, 139.142156]
y3 = [35.859592, 35.225989, 35.182486, 35.139579, 35.097818,
      35.057672, 35.019512, 34.983607, 34.95014, 34.919216,
      34.890891, 34.865177, 34.842062, 34.82152, 34.803513,
      34.788004, 34.774954, 34.764328, 34.756096, 34.750234,
      34.746724, 34.745556, 34.746724, 34.750234, 34.756096,
      34.764328, 34.774954, 34.788004, 34.803513, 34.82152,
      34.842062, 34.865177, 34.890891, 34.919216, 34.95014,
      34.983607, 35.019512, 35.057672, 35.097818, 35.139579,
      35.182486, 35.225989, 35.859592, 35.859592]
bbox_01 = np.column_stack([x1, y1])
bbox_02 = np.column_stack([x2, y2])
bbox_03 = np.column_stack([x3, y3])

NESTS = [
    dict(shp=GSHHS, bbox=bbox_01, h0=10e3, min_el=10e3, max_el=50e3,
         wl=30, dt=0.0, grade=0.3, R=3, slp=50, fl=-50),
    dict(shp=COAST2, bbox=bbox_02, h0=1e3, min_el=1e3, max_el=2e3,
         wl=30, dt=0.0, grade=0.1, R=3, slp=50, fl=-50),
    dict(shp=COAST2, bbox=bbox_03, h0=1e2, min_el=1e2, max_el=5e2,
         wl=30, dt=0.0, grade=0.1, R=3, slp=50, fl=-50),
]

sdfs, edges, gdats = [], [], []
import time as _time
_t0 = _time.time()

def _tick(msg):
    print(f"[time +{_time.time()-_t0:6.0f}s] {msg}", flush=True)

for k, nz in enumerate(NESTS, 1):
    _tick(f"nest{k} geodata start")
    print(f"[nest{k}] geodata: {Path(nz['shp']).name} "
          f"h0={nz['h0']:g} m", flush=True)
    reg = Region((float(nz['bbox'][:, 0].min()),
                  float(nz['bbox'][:, 0].max()),
                  float(nz['bbox'][:, 1].min()),
                  float(nz['bbox'][:, 1].max())), 4326)
    shore = Shoreline(str(nz['shp']), nz['bbox'], nz['h0'] * DEG)
    sdf = om.signed_distance_function(shore)
    dem = DEM(str(SRTM_PACIFIC if k == 1 else SRTM_KANTO),
              bbox=reg)
    print(f"[nest{k}] edgefx: fs={nz['R']} wl={nz['wl']} "
          f"slp={nz['slp']} fl={nz['fl']} max_el={nz['max_el']:g} "
          f"dt={nz['dt']:g} g={nz['grade']}", flush=True)
    comps = [
        om.feature_sizing_function(
            shore, sdf, r=nz['R'],
            max_edge_length=nz['max_el'] * DEG),
        om.wavelength_sizing_function(
            dem, wl=nz['wl'], min_edgelength=nz['min_el'] * DEG,
            max_edge_length=nz['max_el'] * DEG),
        om.bathymetric_gradient_sizing_function(
            dem, slope_parameter=nz['slp'],
            filter_quotient=abs(nz['fl']),
            type_of_filter="barotropic",
            min_edge_length=nz['min_el'] * DEG,
            max_edge_length=nz['max_el'] * DEG),
    ]
    grid, dt_used = om.finalize_sizing(
        comps, dem=dem,
        hmin=nz['min_el'], max_edge_length=nz['max_el'],
        gradation=nz['grade'],
        courant={"timestep": nz['dt'], "max": 0.5},
    )
    print(f"[nest{k}] finalize done (auto dt = {dt_used})",
          flush=True)
    _tick(f"nest{k} sizing done")
    sdfs.append(sdf)
    edges.append(grid)
    gdats.append((shore, dem))

_tick("meshgen start")
print("[meshgen] build (itmax=%d) ..." % MESH_ITMAX, flush=True)
p, t = om.generate_multiscale_mesh(sdfs, edges, max_iter=MESH_ITMAX,
                                   seed=0)
print(f"[meshgen] NP={len(p):,} NE={len(t):,}", flush=True)
_tick("meshgen done")
p, t = om.make_mesh_boundaries_traversable(p, t)
p, t = om.delete_faces_connected_to_one_face(p, t)
from oceanmesh.mesh_improve import area_length_quality  # noqa: E402

q = area_length_quality(p, t)
print(f"[clean] NP={len(p):,} NE={len(t):,} AL min/mean = "
      f"{q.min():.3f}/{q.mean():.3f}", flush=True)

# m = interp(m, gdats, mindepth 0.01)
dem3 = gdats[2][1]
b = om.interp_bathymetry(p, t, dem3, method="cell-averaging",
                         min_depth=0.01)
# m = make_bc(m,'auto',gdat_01,'both')
bc = om.make_bc_auto(p, t, depth=b, classifier="both",
                     depth_lim=10.0, cut_lim=10)
om.write_fort14(str(OUT / "tb_varres_3regions.14"), p, t, depth=b,
                boundaries=bc)
print("[out] tb_varres_3regions.14 written", flush=True)

# quick bay-window band stats vs sample for the log
e = np.vstack([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
e = np.unique(np.sort(e, axis=1), axis=0)
L = np.linalg.norm(p[e[:, 0]] - p[e[:, 1]], axis=1) / DEG
lat = 0.5 * (p[e[:, 0], 1] + p[e[:, 1], 1])
lonm = 0.5 * (p[e[:, 0], 0] + p[e[:, 1], 0])
inbay = (lonm > 139.55) & (lonm < 140.2)
for lo, hi in [(34.95, 35.10), (35.10, 35.20), (35.20, 35.30),
               (35.30, 35.40), (35.40, 35.55), (35.55, 35.70)]:
    s = inbay & (lat >= lo) & (lat < hi)
    if s.sum():
        print(f"[band] {lo:.2f}-{hi:.2f}: p10/p50/p90 = "
              f"{np.percentile(L[s], [10, 50, 90]).round(0)}",
              flush=True)
