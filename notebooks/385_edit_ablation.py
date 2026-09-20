"""Edit ablation: what the four human-judgment edits actually do.

Three columns per site -- goto2023 sample | certified chain (edits ON)
| same chain with SR_EDITS_EXCLUDE set to all four edits (job 115210,
outputs/sample_repro_noedits/). One row per edit site plus the two
sites the ablation changed most. The QA gate, the connectivity
comparator and the one-wide ledger all stay green without the edits,
so the difference is only visible as geometry.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import sys  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from pyproj import Transformer  # noqa: E402
from shapely.ops import unary_union  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fvcom_mesh_tools.plotting import add_atlas_grid, use_readable_style  # noqa: E402

use_readable_style()
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs/sample_repro"
NOED = ROOT / "outputs/sample_repro_noedits"
COSW = float(np.cos(np.deg2rad(35.35)))
tr = Transformer.from_crs(32654, 4326, always_xy=True)


def load14_ll(p):
    ls = Path(p).read_text().splitlines()
    ne, nn = map(int, ls[1].split()[:2])
    nod = np.array([ls[2 + i].split()[1:3] for i in range(nn)], float)
    tri = np.array([ls[2 + nn + i].split()[2:5] for i in range(ne)], int) - 1
    lon, lat = tr.transform(nod[:, 0], nod[:, 1])
    return np.column_stack([lon, lat]), tri


MESHES = [
    ("goto2023 sample", *load14_ll(OUT / "sample_original.14")),
    ("edits ON (certified)", *load14_ll(OUT / "sample_repro_final.14")),
    ("edits OFF (ablation)", *load14_ll(NOED / "sample_repro_final.14")),
]
land = unary_union(list(gpd.read_file(ROOT / "outputs/tb_varres_3r/land_osm_wide.shp").geometry))
gl = gpd.GeoSeries([land], crs="EPSG:4326")

SITES = [
    ("edit_001 Haneda D-runway", 139.7995, 35.5262, 0.020,
     "OSM draws the pier passage as land; the edit restores it as water"),
    ("edit_003 west domain edge", 139.6080, 35.6200, 0.048,
     "Tama river cut by the artificial edge: sliver forced to land"),
    ("edit_004/005 Urayasu", 139.8285, 35.6330, 0.024,
     "harbour -> land + 700 m axis-centred channel"),
]

fig, axes = plt.subplots(len(SITES), 3, figsize=(6.4 * 3, 6.9 * len(SITES)))
for r, (sid, cx, cy, half, how) in enumerate(SITES):
    for c, (title, po, tri) in enumerate(MESHES):
        ax = axes[r, c]
        gl.plot(ax=ax, color="0.88", edgecolor="0.6", linewidth=0.5)
        ax.triplot(po[:, 0], po[:, 1], tri, color="steelblue", linewidth=0.55, alpha=0.9)
        ax.plot(cx, cy, marker="o", ms=26, mfc="none", mec="darkorange", mew=2.0, zorder=7)
        ax.set_xlim(cx - half / COSW, cx + half / COSW)
        ax.set_ylim(cy - half, cy + half)
        ax.set_aspect(1.0 / COSW)
        add_atlas_grid(ax, crs="EPSG:4326")
        ttl = f"{sid}  {title}"
        if c == 1:
            ttl += f"\n{how}"
        ax.set_title(ttl, fontsize=13)
fig.tight_layout()
out = ROOT / "outputs/figures/385_edit_ablation.png"
fig.savefig(out, dpi=120, bbox_inches="tight")
print(f"[385] saved {out}", flush=True)
