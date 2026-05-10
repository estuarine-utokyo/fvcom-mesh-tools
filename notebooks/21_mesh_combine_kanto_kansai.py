"""PoC #21: stitch Tokyo Bay and Osaka Bay into one fort.14.

Demonstrates ``fmesh-mesh-combine`` with the ``disjoint`` strategy:
the two oceanmesh-generated regional meshes (PoCs #19 and #20) live
on opposite ends of Honshu and never share a node, so we just
concatenate them with index renumbering and let the boundaries carry
forward verbatim.

Outputs:
    outputs/21_kanto_kansai.14
    outputs/21_kanto_kansai_summary.txt
    outputs/21_kanto_kansai_mesh.png
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from fvcom_mesh_tools.algorithms import (  # noqa: E402
    alpha_quality,
    min_interior_angle,
    signed_areas,
)
from fvcom_mesh_tools.cli.meshcombine import main as meshcombine_main  # noqa: E402
from fvcom_mesh_tools.io import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
TOKYO = REPO_ROOT / "outputs" / "19_tokyo_bay_oceanmesh.14"
OSAKA = REPO_ROOT / "outputs" / "20_osaka_bay_oceanmesh.14"
OUT_DIR = REPO_ROOT / "outputs"
F14_OUT = OUT_DIR / "21_kanto_kansai.14"
SUMMARY_TXT = OUT_DIR / "21_kanto_kansai_summary.txt"
MESH_PNG = OUT_DIR / "21_kanto_kansai_mesh.png"


def plot_mesh(mesh, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(11, 6), dpi=120)
    ax.triplot(
        mesh.nodes[:, 0], mesh.nodes[:, 1], mesh.elements,
        color="0.7", lw=0.1,
    )
    plotted_open = False
    for seg in mesh.open_boundaries:
        ax.plot(
            mesh.nodes[seg, 0], mesh.nodes[seg, 1], "-",
            color="tab:red", lw=1.0,
            label="open" if not plotted_open else None,
        )
        plotted_open = True
    plotted_land = False
    plotted_river = False
    for ib, seg in mesh.land_boundaries:
        if ib == 21:
            ax.plot(
                mesh.nodes[seg, 0], mesh.nodes[seg, 1], "-",
                color="tab:blue", lw=2.0,
                label="river (ibtype=21)" if not plotted_river else None,
            )
            plotted_river = True
        else:
            ax.plot(
                mesh.nodes[seg, 0], mesh.nodes[seg, 1], "-",
                color="0.4", lw=0.4,
                label=f"land (ibtype={ib})" if not plotted_land else None,
            )
            plotted_land = True
    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title("PoC #21 fmesh-mesh-combine disjoint  Kanto + Kansai")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def main() -> None:
    for p in (TOKYO, OSAKA):
        if not p.exists():
            raise SystemExit(
                f"required input missing: {p}. "
                "Run notebooks/19_oceanmesh_full_pipeline.py and "
                "notebooks/20_osaka_bay_oceanmesh.py first."
            )
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tokyo = read_fort14(TOKYO)
    osaka = read_fort14(OSAKA)

    print("[21] running fmesh-mesh-combine --strategy disjoint ...")
    t0 = time.perf_counter()
    rc = meshcombine_main([
        str(TOKYO), str(OSAKA), str(F14_OUT),
        "--strategy", "disjoint",
        "--title", "Kanto + Kansai PoC #21 (disjoint combine)",
    ])
    wall = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"fmesh-mesh-combine exited {rc}")
    print(f"[21] wall: {wall:.2f} s")

    combined = read_fort14(F14_OUT)
    expected_np = tokyo.n_nodes + osaka.n_nodes
    expected_ne = tokyo.n_elements + osaka.n_elements

    n_open = len(combined.open_boundaries)
    n_land = len(combined.land_boundaries)
    n_river_segs = sum(1 for ib, _ in combined.land_boundaries if ib == 21)
    flipped = int((signed_areas(combined) <= 0).sum())
    q = alpha_quality(combined)
    a = min_interior_angle(combined)

    summary_lines = [
        f"input A: {TOKYO}  (NP={tokyo.n_nodes:,} NE={tokyo.n_elements:,})",
        f"input B: {OSAKA}  (NP={osaka.n_nodes:,} NE={osaka.n_elements:,})",
        "strategy: disjoint",
        f"wall:    {wall:.2f} s",
        "",
        f"NP={combined.n_nodes:,}  NE={combined.n_elements:,}  flipped={flipped}",
        f"NP expected={expected_np:,}  NE expected={expected_ne:,}",
        "",
        "[Boundary structure preserved from inputs]",
        f"  open boundaries     : {n_open} (sum of inputs: "
        f"{len(tokyo.open_boundaries) + len(osaka.open_boundaries)})",
        f"  land segments total : {n_land} (sum of inputs: "
        f"{len(tokyo.land_boundaries) + len(osaka.land_boundaries)})",
        f"    ibtype=21 (river) : {n_river_segs}",
        "",
        "[Quality (combined)]",
        f"  alpha mean          : {q.mean():.4f}",
        f"  min-angle p50       : {np.percentile(a, 50):.2f}",
        f"  min-angle < 20 deg  : {(a < 20).mean() * 100:.2f} %",
    ]
    summary = "\n".join(summary_lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[21] wrote {SUMMARY_TXT}")

    plot_mesh(combined, MESH_PNG)
    print(f"[21] wrote {MESH_PNG}")

    assert combined.n_nodes == expected_np, (
        f"NP mismatch: {combined.n_nodes} != {expected_np}"
    )
    assert combined.n_elements == expected_ne, (
        f"NE mismatch: {combined.n_elements} != {expected_ne}"
    )
    assert flipped == 0, f"{flipped} flipped triangles after combine"
    assert n_open == (
        len(tokyo.open_boundaries) + len(osaka.open_boundaries)
    ), "open boundary count mismatch"
    assert n_land == (
        len(tokyo.land_boundaries) + len(osaka.land_boundaries)
    ), "land boundary count mismatch"
    assert n_river_segs == sum(
        1 for ib, _ in tokyo.land_boundaries if ib == 21
    ) + sum(
        1 for ib, _ in osaka.land_boundaries if ib == 21
    ), "river ibtype=21 count mismatch"


if __name__ == "__main__":
    main()
