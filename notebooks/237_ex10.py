# Ladder step 5: Example_10_Multiscale_Smoother translation.
import os, sys, logging, time
import numpy as np
logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                    format="%(levelname)s %(name)s: %(message)s")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.path.insert(0, os.path.expanduser("~/Github/oceanmesh"))
from pathlib import Path
import oceanmesh as om
from oceanmesh.grid import Grid
from oceanmesh.signed_distance_function import create_bbox

DEG = 1.0 / 111e3
OUT = Path("outputs/om2d_examples/ex10")
OUT.mkdir(parents=True, exist_ok=True)

h1 = 1e3 * DEG          # inner unit box, 1 km
h2 = 1e4 * DEG          # outer box, 10 km
bbox1 = (0.0, 1.0, 0.0, 1.0)
bbox2 = (-1.0, 2.0, -1.0, 2.0)

def uniform_grid(bbox, h):
    g = Grid(bbox=bbox, dx=h, extrapolate=True, values=float(h),
             crs=4326)
    g.hmin = h
    # OM2D edgefx applies gradation 0.2; a uniform field is
    # already grade-0 so limiting is a no-op
    g.build_interpolant()
    return g

dom2 = create_bbox(bbox2)
dom1 = create_bbox(bbox1)
dom2.crs = "EPSG:4326"
dom1.crs = "EPSG:4326"
# box nests: the covering (nest footprint) is the box itself
dom2.covering = dom2.domain
dom1.covering = dom1.domain
g2 = uniform_grid(bbox2, h2)
g1 = uniform_grid(bbox1, h1)

t0 = time.time()
p, t = om.generate_multiscale_mesh(
    [dom2, dom1], [g2, g1], max_iter=100, seed=0,
    qual_tol=0.0025, cleanup="none", gradation=0.2,
)
print(f"[ex10] raw NP={len(p):,} NT={len(t):,} "
      f"+{time.time()-t0:.0f}s", flush=True)
np.save(OUT / "p_raw.npy", p); np.save(OUT / "t_raw.npy", t)
from oceanmesh.clean import om2d_default_clean
pc, tc = om2d_default_clean(p, t)
from oceanmesh.fix_mesh import simp_qual
q = simp_qual(pc, tc)
print(f"[ex10] clean NP={len(pc):,} NT={len(tc):,} "
      f"qual min={q.min():.4f} mean={q.mean():.4f}", flush=True)
np.save(OUT / "p.npy", pc); np.save(OUT / "t.npy", tc)
print("[ex10] saved", flush=True)
