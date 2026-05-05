"""Read tb_futtsu20220311.14 and emit a summary + thumbnail plot.

PoC #1 for the OceanMesh2D -> Python migration. Loads the reference Tokyo
Bay mesh through fvcom-mesh-tools' native fort.14 reader and prints node /
element / boundary statistics.

Inputs:
    data/mesh/reference/tokyo_bay/tb_futtsu20220311.14   (symlink into $DATA_DIR)

Outputs:
    outputs/01_tb_futtsu20220311_summary.txt
    outputs/01_tb_futtsu20220311_mesh.png
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from fvcom_mesh_tools.io import Fort14Mesh, read_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
MESH_PATH = REPO_ROOT / "data" / "mesh" / "reference" / "tokyo_bay" / "tb_futtsu20220311.14"
OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "01_tb_futtsu20220311_summary.txt"
MESH_PNG = OUT_DIR / "01_tb_futtsu20220311_mesh.png"


def summarize(mesh: Fort14Mesh) -> str:
    xmin, ymin, xmax, ymax = mesh.bbox
    d = mesh.depths

    n_open_nodes = sum(len(b) for b in mesh.open_boundaries)
    n_land_nodes = sum(len(ids) for _, ids in mesh.land_boundaries)
    ibtype_counter = Counter(ib for ib, _ in mesh.land_boundaries)
    ibtype_str = ", ".join(f"{k}: {v}" for k, v in sorted(ibtype_counter.items()))

    lines = [
        f"file:               {MESH_PATH}",
        f"title:              {mesh.title.strip()!r}",
        f"nodes (NP):         {mesh.n_nodes:,}",
        f"elements (NE):      {mesh.n_elements:,}",
        f"bbox (lon, lat):    [{xmin:.6f}, {ymin:.6f}] -> [{xmax:.6f}, {ymax:.6f}]",
        f"bbox extent:        dlon={xmax - xmin:.4f} deg   dlat={ymax - ymin:.4f} deg",
        f"depth (file value): min={d.min():.4f}   max={d.max():.4f}   mean={d.mean():.4f}",
        f"open boundaries:    {len(mesh.open_boundaries)} segments / {n_open_nodes:,} nodes",
        f"land boundaries:    {len(mesh.land_boundaries)} segments / {n_land_nodes:,} nodes",
        f"land ibtype dist:   {{{ibtype_str}}}",
    ]
    if mesh.land_boundaries:
        sizes = sorted((len(ids) for _, ids in mesh.land_boundaries), reverse=True)
        lines.append(f"land top-5 sizes:   {sizes[:5]}")
        lines.append(f"land bottom-5:      {sizes[-5:]}")
    return "\n".join(lines)


def thumbnail(mesh: Fort14Mesh, png_path: Path) -> None:
    """Quick mesh sketch: node scatter + per-segment boundary polylines."""
    fig, ax = plt.subplots(figsize=(8, 7), dpi=120)
    ax.scatter(
        mesh.nodes[:, 0], mesh.nodes[:, 1],
        s=0.2, color="0.6", linewidths=0,
    )

    for seg in mesh.open_boundaries:
        coords = mesh.nodes[seg]
        ax.plot(coords[:, 0], coords[:, 1], "-", color="tab:blue", lw=1.5)
    for _ibtype, seg in mesh.land_boundaries:
        coords = mesh.nodes[seg]
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
    mesh = read_fort14(MESH_PATH)

    text = summarize(mesh)
    print(text)
    SUMMARY_TXT.write_text(text + "\n", encoding="utf-8")
    print(f"[01] wrote {SUMMARY_TXT}")

    thumbnail(mesh, MESH_PNG)
    print(f"[01] wrote {MESH_PNG}")


if __name__ == "__main__":
    main()
