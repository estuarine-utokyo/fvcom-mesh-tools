# TB nest1/nest2: our FINAL sizing fields vs the MATLAB dumps
# (ml_fh1.mat / ml_fh2.mat), banded by ML value and by depth.
import os, sys, re, logging
import numpy as np
import h5py
logging.basicConfig(level=logging.WARNING, stream=sys.stdout)
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline
from oceanmesh.finalize import _dem_on_grid

DEG = 1.0 / 111e3
OM2D = os.path.expanduser("~/Github/OceanMesh2D")
OUT = Path("outputs/tb_varres_3r")
src = open("notebooks/150_tb_varres_3r_translation.py").read()


def grab(name):
    m = re.search(rf"^{name} = \[(.*?)\]", src, re.S | re.M)
    return [float(v) for v in re.findall(r"-?\d+\.?\d*", m.group(1))]


x2 = [138.196907945361, 139.784782123012, 141.294666567004,
      139.718998032166, 138.285893837899, 138.196907945361]
y2 = [34.6003740729489, 33.5207179317026, 35.2698037998702,
      35.9357551145740, 35.2841759386230, 34.6003740729489]
CFG = [
    ("nest1", np.column_stack([grab("x1"), grab("y1")]),
     OM2D + "/datasets/GSHHS_shp/f/GSHHS_f_L1.shp",
     OM2D + "/datasets/TokyoBay/dem/SRTM15_pacific_4min.nc",
     10e3, 50e3, 0.3, "ml_fh1.mat",
     [(0, 12e3), (12e3, 20e3), (20e3, 35e3), (35e3, 50.1e3)]),
    ("nest2", np.column_stack([x2, y2]),
     "outputs/tb_varres_3r/land_osm_wide.shp",
     OM2D + "/datasets/TokyoBay/dem/SRTM15_kanto_15s.nc",
     1e3, 2e3, 0.1, "ml_fh2.mat",
     [(0, 1.2e3), (1.2e3, 1.6e3), (1.6e3, 2.01e3)]),
]

for tag, poly, shp, demf, h0, mx, gr, mlf, bands in CFG:
    reg = Region((poly[:, 0].min(), poly[:, 0].max(),
                  poly[:, 1].min(), poly[:, 1].max()), 4326)
    sh = Shoreline(shp, poly, h0 * DEG)
    sdf = om.signed_distance_function(sh)
    dem = DEM(demf, bbox=reg, nc_reader="coords")
    f = om.feature_sizing_function(
        sh, sdf, r=3, max_edge_length=mx * DEG,
        lattice_anchor=(dem.bbox[0], dem.bbox[2]))
    w = om.wavelength_sizing_function(
        dem, wl=30, min_edgelength=h0 * DEG,
        max_edge_length=mx * DEG, grid_dx=h0 * DEG)
    s = om.bathymetric_gradient_sizing_function(
        dem, slope_parameter=50, filter_quotient=50,
        type_of_filter="barotropic", min_edge_length=h0 * DEG,
        max_edge_length=mx * DEG, grid_dx=h0 * DEG)
    g, dta = om.finalize_sizing(
        [f, w, s], dem=dem, shoreline=sh, hmin=h0,
        max_edge_length=mx, gradation=gr,
        courant={"timestep": 0.0})
    print(f"[{tag}] dt = {dta:.4f}", flush=True)

    with h5py.File(OUT / mlf, "r") as fml:
        xg = np.asarray(fml["xg"]).ravel()
        yg = np.asarray(fml["yg"]).ravel()
        hh = np.asarray(fml["hh"])
    if hh.shape == (len(yg), len(xg)):
        hh = hh.T
    step = max(1, len(xg) // 900)
    xs = xg[::step]; ys = yg[::step]
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    q = np.column_stack([X.ravel(), Y.ravel()])
    mlv = hh[::step, ::step]
    ours = np.asarray(g.eval(q)).reshape(X.shape) / DEG
    z = _dem_on_grid(f, dem)
    from scipy.interpolate import RegularGridInterpolator
    xf, yf = f.create_grid()
    Fz = RegularGridInterpolator((xf[:, 0], yf[0, :]), z,
                                 bounds_error=False)
    zz = Fz(q).reshape(X.shape)
    r = ours / mlv
    fin = np.isfinite(r) & (mlv > 0) & np.isfinite(zz) & (zz < 0)
    print(f"[{tag}] ocean ratio p10/50/90 = "
          f"{np.percentile(r[fin], [10, 50, 90]).round(3)}", flush=True)
    for lo, hi in bands:
        m = fin & (mlv >= lo) & (mlv < hi)
        if m.sum() > 100:
            print(f"[{tag}]  ml [{lo/1e3:g}k,{hi/1e3:g}k): "
                  f"p50={np.percentile(r[m], 50):.3f} "
                  f"n={int(m.sum()):,}", flush=True)
    for lo, hi, zt in [(-11000, -2000, "deep"), (-2000, -200, "slope"),
                       (-200, 0, "shelf")]:
        m = fin & (zz >= lo) & (zz < hi)
        if m.sum() > 100:
            print(f"[{tag}]  z {zt:6s}: p50={np.percentile(r[m], 50):.3f} "
                  f"ml_p50={np.percentile(mlv[m], 50):,.0f} "
                  f"ours_p50={np.percentile(ours[m], 50):,.0f}", flush=True)
print("[fh12] done", flush=True)
