# The delivered mesh on its own: the whole mesh, the patch, and any close-ups.
#
# plotting.plot_mesh_views draws it in the project's fixed colours: thin grey
# interior edges, every solid boundary black -- coastline, quay, and a wall
# represented as a line alike -- and the open boundary red.  The same figures
# come from `fmesh-plot-views`; this adds the patch's own extent.
#
#   FMESH_VIEWS=port:391900:394100:3908350:3910600 \
#       python notebooks/434_final_mesh.py <refinement output dir>
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.io.fort14 import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import plot_mesh_views  # noqa: E402

out = Path(sys.argv[1]).resolve()
m = read_fort14(next(out.glob("*.14")))
fc = np.load(out / "fill_constraints.npz")
lo, hi = fc["pfix"].min(axis=0), fc["pfix"].max(axis=0)
pad = 0.05 * float(max(hi - lo))
views = {"patch": (lo[0] - pad, hi[0] + pad, lo[1] - pad, hi[1] + pad)}
# close-ups, "name:x0:x1:y0:y1" plus-separated (qsub -v splits on commas)
for item in filter(None, os.environ.get("FMESH_VIEWS", "").split("+")):
    name, *box = item.split(":")
    views[name] = tuple(float(v) for v in box)
for png in plot_mesh_views(m, out, views):
    print(f"wrote {png}")
