# The delivered mesh on its own: the whole patch, and any close-ups asked for.
#
# Drawn with plotting.draw_mesh, so the colours are the project's fixed ones:
# thin grey interior edges, every solid boundary black -- coastline, quay,
# and a wall represented as a line alike -- and the open boundary red.
# One PNG per view at MESH_PNG_DPI.
#
#   FMESH_VIEWS=port:391900:394100:3908350:3910600 \
#       python notebooks/434_final_mesh.py <refinement output dir>
import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from fvcom_mesh_tools.io.fort14 import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import (  # noqa: E402
    MESH_PNG_DPI,
    boundary_legend,
    draw_mesh,
)

out = Path(sys.argv[1]).resolve()
m = read_fort14(next(out.glob("*.14")))
fc = np.load(out / "fill_constraints.npz")
lo, hi = fc["pfix"].min(axis=0), fc["pfix"].max(axis=0)
pad = 0.05 * float(max(hi - lo))
# The whole patch always; close-ups from FMESH_VIEWS, "name:x0:x1:y0:y1+..."
# (plus-separated, since qsub -v splits variables on commas)
views = {"patch": (lo[0] - pad, hi[0] + pad, lo[1] - pad, hi[1] + pad)}
for item in filter(None, os.environ.get("FMESH_VIEWS", "").split("+")):
    name, *box = item.split(":")
    views[name] = tuple(float(v) for v in box)
for name, (x0, x1, y0, y1) in views.items():
    fig, ax = plt.subplots(figsize=(9, 9 * (y1 - y0) / (x1 - x0)))
    lw = 0.9 if name == "patch" else 1.4
    n = draw_mesh(ax, m.nodes, m.elements, m.open_boundaries,
                  mesh_lw=0.12 if name == "patch" else 0.25, boundary_lw=lw)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m, EPSG:32654)")
    ax.set_ylabel("y (m, EPSG:32654)")
    ax.set_title(f"{m.title}  NP={m.n_nodes:,} NE={m.n_elements:,}  [{name}]", fontsize=9)
    boundary_legend(ax, open_boundary=n["n_open"] > 0, loc="upper right", fontsize=8)
    fig.tight_layout()
    png = out / f"final_mesh_{name}.png"
    fig.savefig(png, dpi=MESH_PNG_DPI)
    plt.close(fig)
    print(f"wrote {png}")
