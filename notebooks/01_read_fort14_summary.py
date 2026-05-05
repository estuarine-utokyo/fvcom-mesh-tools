"""Read tb_futtsu20220311.14 and emit a summary + thumbnail plot.

This is the first PoC for fvcom-mesh-tools' Python migration. It loads the
reference Tokyo Bay mesh through the existing oceanmesh-tools parser
(installed as a sibling editable package) and prints node/element/boundary
statistics. A native loader will replace the sibling import in a later step.

Inputs:
    data/mesh/reference/tokyo_bay/tb_futtsu20220311.14   (symlink into $DATA_DIR)

Outputs:
    outputs/01_tb_futtsu20220311_summary.txt
    outputs/01_tb_futtsu20220311_mesh.png
"""

from __future__ import annotations

from pathlib import Path
from statistics import mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from oceanmesh_tools.io.fort14 import parse_fort14


REPO_ROOT = Path(__file__).resolve().parent.parent
MESH_PATH = REPO_ROOT / "data" / "mesh" / "reference" / "tokyo_bay" / "tb_futtsu20220311.14"
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "01_tb_futtsu20220311_summary.txt"
MESH_PNG = OUT_DIR / "01_tb_futtsu20220311_mesh.png"


def summarize(mesh) -> str:
    xs = [x for (_, x, _, _) in mesh.nodes]
    ys = [y for (_, _, y, _) in mesh.nodes]
    ds = [d for (_, _, _, d) in mesh.nodes]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    dmin, dmax, dmean = min(ds), max(ds), mean(ds)

    lines = [
        f"file:               {MESH_PATH}",
        f"title:              {mesh.title!r}",
        f"nodes (NP):         {mesh.n_nodes:,}",
        f"elements (NE):      {mesh.n_elements:,}",
        f"bbox (lon, lat):    [{xmin:.6f}, {ymin:.6f}] -> [{xmax:.6f}, {ymax:.6f}]",
        f"bbox extent:        dlon={xmax - xmin:.4f} deg   dlat={ymax - ymin:.4f} deg",
        f"depth (file value): min={dmin:.4f}   max={dmax:.4f}   mean={dmean:.4f}",
        f"open boundaries:    {len(mesh.open_boundaries)} segments,"
        f"  total nodes={sum(len(b) for b in mesh.open_boundaries):,}",
        f"land boundaries:    {len(mesh.land_boundaries)} segments,"
        f"  total nodes={sum(len(b) for b in mesh.land_boundaries):,}",
    ]
    if mesh.open_boundaries:
        sizes = [len(b) for b in mesh.open_boundaries]
        lines.append(f"  open seg sizes:   min={min(sizes)}  max={max(sizes)}  mean={mean(sizes):.1f}")
    if mesh.land_boundaries:
        sizes = [len(b) for b in mesh.land_boundaries]
        lines.append(f"  land seg sizes:   min={min(sizes)}  max={max(sizes)}  mean={mean(sizes):.1f}")
    return "\n".join(lines)


def thumbnail(mesh, png_path: Path) -> None:
    """Quick coarse mesh sketch (no triangulation; just node scatter + boundaries)."""
    xs = np.fromiter((x for (_, x, _, _) in mesh.nodes), dtype=np.float64, count=mesh.n_nodes)
    ys = np.fromiter((y for (_, _, y, _) in mesh.nodes), dtype=np.float64, count=mesh.n_nodes)
    node_xy_by_id = {nid: (x, y) for (nid, x, y, _d) in mesh.nodes}

    fig, ax = plt.subplots(figsize=(8, 7), dpi=120)
    ax.scatter(xs, ys, s=0.2, color="0.6", linewidths=0, label=f"nodes (n={mesh.n_nodes})")

    for seg in mesh.open_boundaries:
        coords = np.array([node_xy_by_id[i] for i in seg if i in node_xy_by_id])
        if len(coords):
            ax.plot(coords[:, 0], coords[:, 1], "-", color="tab:blue", lw=1.5)
    for seg in mesh.land_boundaries:
        coords = np.array([node_xy_by_id[i] for i in seg if i in node_xy_by_id])
        if len(coords):
            ax.plot(coords[:, 0], coords[:, 1], "-", color="tab:red", lw=0.6)

    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title(f"{MESH_PATH.name}  —  NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}")
    ax.grid(True, lw=0.3, color="0.85")
    proxy_open = plt.Line2D([], [], color="tab:blue", lw=1.5, label="open boundary")
    proxy_land = plt.Line2D([], [], color="tab:red", lw=0.8, label="land boundary")
    ax.legend(handles=[proxy_open, proxy_land], loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not MESH_PATH.exists():
        raise SystemExit(f"mesh not found: {MESH_PATH}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[01] reading {MESH_PATH}")
    mesh = parse_fort14(MESH_PATH)

    text = summarize(mesh)
    print(text)
    SUMMARY_TXT.write_text(text + "\n", encoding="utf-8")
    print(f"[01] wrote {SUMMARY_TXT}")

    thumbnail(mesh, MESH_PNG)
    print(f"[01] wrote {MESH_PNG}")


if __name__ == "__main__":
    main()
