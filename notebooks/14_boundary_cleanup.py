"""PoC #14: combined island filter + open-segment merge.

PoC #13 demonstrated --min-{polygon,island}-area-m2 cuts land
segments from 166 to 23 with no quality regression. The remaining
mismatch with the reference fort.14 is the open-boundary count: our
classifier emits 3 segments where the reference has 1 contiguous arc.
PoC #14 sweeps --open-merge-coast-gap and confirms that combining
both options produces a fort.14 whose *shape* (segment counts) is
within striking distance of the reference, while keeping the same
quality numbers as PoC #12 / #13.

Outputs:
    outputs/14_boundary_cleanup_islands_only.14
    outputs/14_boundary_cleanup_islands_open5.14
    outputs/14_boundary_cleanup_islands_open50.14
    outputs/14_boundary_cleanup_summary.txt
    outputs/14_boundary_cleanup_meshes.png
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

from fvcom_mesh_tools.algorithms import (  # noqa: E402
    alpha_quality,
    min_interior_angle,
    signed_areas,
)
from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main  # noqa: E402
from fvcom_mesh_tools.io import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO_ROOT / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO_ROOT / "outputs"
F14_ISLANDS = OUT_DIR / "14_boundary_cleanup_islands_only.14"
F14_OPEN5 = OUT_DIR / "14_boundary_cleanup_islands_open5.14"
F14_OPEN50 = OUT_DIR / "14_boundary_cleanup_islands_open50.14"
SUMMARY_TXT = OUT_DIR / "14_boundary_cleanup_summary.txt"
MESHES_PNG = OUT_DIR / "14_boundary_cleanup_meshes.png"


def _run(label: str, out: Path, *, open_merge: int) -> float:
    args = [
        str(DEM), str(out),
        "--hmin", "200", "--hmax", "5000", "--zmax", "0.0",
        "--interp-method", "linear", "--land-ibtype", "20",
        "--coastline", str(COASTLINE),
        "--coast-target-size", "200", "--coast-expansion-rate", "0.005",
        "--min-polygon-area-m2", "1000000",
        "--min-island-area-m2", "100000",
        "--open-merge-coast-gap", str(open_merge),
        "--quality-pass", "6", "--smooth-iters", "5", "--smooth-alpha", "0.5",
        "--perpfix-iters", "1", "--quiet",
    ]
    print(f"[14] {label}: running ...")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    t = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh ({label}) exited {rc}")
    print(f"[14] {label}: {t:.2f} s")
    return t


def _summary(label: str, mesh, t: float) -> str:
    n_open = len(mesh.open_boundaries)
    n_open_nodes = sum(int(b.size) for b in mesh.open_boundaries)
    n_land = len(mesh.land_boundaries)
    n_land_nodes = sum(int(b.size) for _, b in mesh.land_boundaries)
    flipped = int((signed_areas(mesh) <= 0).sum())
    q = alpha_quality(mesh)
    a = min_interior_angle(mesh)
    return (
        f"=== {label} (wall {t:.2f} s) ===\n"
        f"  NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  flipped={flipped}\n"
        f"  open : {n_open} segments / {n_open_nodes:,} nodes\n"
        f"  land : {n_land} segments / {n_land_nodes:,} nodes\n"
        f"  alpha mean={q.mean():.4f}  alpha<0.3={(q < 0.3).mean() * 100:.2f} %  "
        f"min-angle<20={(a < 20).mean() * 100:.2f} %"
    )


def plot_meshes(meshes: dict, png: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), dpi=120, sharex=True, sharey=True)
    titles = list(meshes.keys())
    for ax, label in zip(axes, titles):
        m = meshes[label]
        ax.triplot(m.nodes[:, 0], m.nodes[:, 1], m.elements, color="0.4", lw=0.15)
        for _, seg in m.land_boundaries:
            ax.plot(m.nodes[seg, 0], m.nodes[seg, 1], "-", color="tab:gray", lw=0.6)
        for seg in m.open_boundaries:
            ax.plot(m.nodes[seg, 0], m.nodes[seg, 1], "-", color="tab:red", lw=1.4)
        n_open = len(m.open_boundaries)
        n_land = len(m.land_boundaries)
        ax.set_aspect("equal")
        ax.set_title(f"{label}\nopen={n_open}  land={n_land}  NP={m.n_nodes:,}")
        ax.set_xlabel("lon (deg)")
        ax.grid(True, lw=0.3, color="0.9")
    axes[0].set_ylabel("lat (deg)")
    fig.suptitle("Boundary-cleanup sweep", y=1.02)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def main() -> None:
    if not DEM.exists() or not COASTLINE.exists():
        raise SystemExit("DEM or coastline not found")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    runs = {
        "islands_only":   (F14_ISLANDS, 0),
        "islands+open5":  (F14_OPEN5,   5),
        "islands+open50": (F14_OPEN50, 50),
    }
    walls = {label: _run(label, out, open_merge=om) for label, (out, om) in runs.items()}
    meshes = {label: read_fort14(out) for label, (out, _om) in runs.items()}

    text = (
        "DEM:        " + str(DEM) + "\n"
        "coastline:  " + str(COASTLINE) + "\n"
        "min_polygon=1e6 m^2  min_island=1e5 m^2\n"
        "Reference: NP=95,551 NE=182,603 open 1/193 land 54/8,464\n\n"
        + "\n\n".join(_summary(label, mesh, walls[label]) for label, mesh in meshes.items())
    )
    print(text)
    SUMMARY_TXT.write_text(text + "\n", encoding="utf-8")
    print(f"[14] wrote {SUMMARY_TXT}")

    plot_meshes(meshes, MESHES_PNG)
    print(f"[14] wrote {MESHES_PNG}")

    for label, mesh in meshes.items():
        assert (signed_areas(mesh) > 0).all(), f"{label}: flipped"


if __name__ == "__main__":
    main()
