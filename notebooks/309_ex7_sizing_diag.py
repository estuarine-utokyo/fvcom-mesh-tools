# Ex7 node-explosion diagnosis: rebuild the global sizing, dump the
# field BEFORE and AFTER the stereo gradation pass, and compute the
# expected-node integral N = (2/sqrt(3)) * sum(dA/fh^2) over wet
# cells for each stage. Golden final NP is 1,372,623; our run seeded
# 2.74M and split to 6.0M — locate which stage doubles N.
import os, sys, logging, time
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh import DEM, Region, Shoreline

OM2D = Path(os.path.expanduser("~/Github/OceanMesh2D"))
DS = OM2D / "datasets"
SCR = Path(os.environ.get("EX7_SCR",
    os.path.expanduser("~/Github/fvcom-mesh-tools/outputs/om2d_examples/ex7glob")))
DEG = 1.0 / 111e3
t0 = time.time()

mll = SCR / "gshhs_l1l6.shp"
bbox = (-180.0, 180.0, -89.0, 90.0)
reg = Region(bbox, 4326)
sh = Shoreline(str(mll), reg.bbox, 4e3*DEG)
sdf = om.signed_distance_function(sh)
dem = DEM(str(DS/"SRTM15+.nc"), bbox=reg)
print(f"[diag7] shoreline+dem +{time.time()-t0:.0f}s", flush=True)
f = om.feature_sizing_function(sh, sdf, r=3,
                               max_edge_length=20e3*DEG)
w = om.wavelength_sizing_function(dem, wl=30,
                                  min_edgelength=4e3*DEG,
                                  max_edge_length=20e3*DEG,
                                  grid_dx=4e3*DEG)
s = om.bathymetric_gradient_sizing_function(
    dem, slope_parameter=10, filter_quotient=50,
    min_edge_length=4e3*DEG, max_edge_length=20e3*DEG,
    type_of_filter="barotropic", grid_dx=4e3*DEG)
g, dta = om.finalize_sizing([f, w, s], dem=dem, shoreline=sh,
                            hmin=4e3, max_edge_length=20e3,
                            gradation=0.25,
                            courant={"timestep": 0.0})
print(f"[diag7] sizing dt={dta:.2f} +{time.time()-t0:.0f}s", flush=True)


def integral(grid, tag):
    xs, ys = grid.create_vectors()
    # subsample x for the sdf pass to bound cost; dy full
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    q = np.column_stack([X.ravel(), Y.ravel()])
    d = np.asarray(sdf.eval(q))
    wet = (d < 0).reshape(X.shape)
    h = np.asarray(grid.values, dtype=float)
    m = wet & np.isfinite(h) & (h > 0)
    # physical cell areas: dx*cos(lat)*dy
    dxm = float(grid.dx) * 111e3
    dym = float(grid.dy) * 111e3
    dA = dxm * dym * np.cos(np.deg2rad(np.clip(Y, -89.9, 89.9)))
    hm = h * 111e3  # degrees -> metres
    N = (2/np.sqrt(3)) * np.sum(dA[m] / hm[m]**2)
    print(f"[diag7] {tag}: expected N = {N:,.0f} "
          f"(wet cells {int(m.sum()):,}, h p10/50/90 = "
          f"{np.percentile(hm[m],10):,.0f}/{np.percentile(hm[m],50):,.0f}/"
          f"{np.percentile(hm[m],90):,.0f} m)", flush=True)
    return N


np.save(SCR/"g_before_stereo.npy", np.asarray(g.values, dtype=float))
N1 = integral(g, "after finalize (pre-stereo-gradation)")
g2 = om.enforce_mesh_gradation(g, gradation=0.25, stereo=True)
np.save(SCR/"g_after_stereo.npy", np.asarray(g2.values, dtype=float))
N2 = integral(g2, "after stereo gradation")
print(f"[diag7] ratio N2/N1 = {N2/N1:.3f}", flush=True)
print(f"[diag7] done +{time.time()-t0:.0f}s", flush=True)
