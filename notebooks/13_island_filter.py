"""PoC #13: filter tiny outer polygons and tiny island holes.

PoC #11 / #12 left ~165 land segments (mostly islands inside the
main bay polygon) plus several detached single-pixel water bodies in
the wet-domain multipolygon. The reference mesh has 54 land
segments. This PoC sweeps two filter thresholds:

* ``--min-polygon-area-m2`` drops detached water polygons (DEM
  artefacts: e.g. small puddles outside the bay).
* ``--min-island-area-m2`` drops small holes inside surviving
  polygons (islands not worth resolving at the chosen ``hmin``).

Each configuration also runs the coastline-aware Hfun + 6-round
quality pass to keep the comparison apples-to-apples with PoC #12.

Outputs:
    outputs/13_filter_none.14
    outputs/13_filter_islands_1ha.14
    outputs/13_filter_islands_10ha.14
    outputs/13_island_filter_summary.txt
    outputs/13_island_filter_meshes.png
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
    REPO_ROOT
    / "data"
    / "coastline"
    / "tokyo_bay"
    / "MLIT_C23"
    / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO_ROOT / "outputs"
F14_NONE = OUT_DIR / "13_filter_none.14"
F14_1HA = OUT_DIR / "13_filter_islands_1ha.14"
F14_10HA = OUT_DIR / "13_filter_islands_10ha.14"
SUMMARY_TXT = OUT_DIR / "13_island_filter_summary.txt"
MESHES_PNG = OUT_DIR / "13_island_filter_meshes.png"

HMIN_M = 200.0
HMAX_M = 5000.0
QP_ROUNDS = 6


def _run(label: str, out: Path, *, min_poly: float, min_island: float) -> float:
    args = [
        str(DEM), str(out),
        "--hmin", str(HMIN_M),
        "--hmax", str(HMAX_M),
        "--zmax", "0.0",
        "--interp-method", "linear",
        "--land-ibtype", "20",
        "--coastline", str(COASTLINE),
        "--coast-target-size", "200",
        "--coast-expansion-rate", "0.005",
        "--quality-pass", str(QP_ROUNDS),
        "--smooth-iters", "5",
        "--smooth-alpha", "0.5",
        "--perpfix-iters", "1",
        "--quiet",
    ]
    if min_poly > 0:
        args += ["--min-polygon-area-m2", str(min_poly)]
    if min_island > 0:
        args += ["--min-island-area-m2", str(min_island)]
    print(f"[13] {label}: running fmesh-buildmesh ...")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    t = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh ({label}) exited {rc}")
    print(f"[13] {label}: {t:.2f} s")
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
        f"  alpha mean={q.mean():.4f}  alpha<0.3 frac={(q < 0.3).mean() * 100:.2f} %\n"
        f"  min-angle p50={float(a.mean()):.2f}  "
        f"min-angle<20 frac={(a < 20).mean() * 100:.2f} %"
    )


def plot_meshes(meshes: dict, png: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), dpi=120, sharex=True, sharey=True)
    titles = list(meshes.keys())
    for ax, label in zip(axes, titles):
        m = meshes[label]
        ax.triplot(m.nodes[:, 0], m.nodes[:, 1], m.elements, color="0.4", lw=0.15)
        for _, seg in m.land_boundaries:
            ax.plot(m.nodes[seg, 0], m.nodes[seg, 1], "-", color="tab:gray", lw=0.4)
        for seg in m.open_boundaries:
            ax.plot(m.nodes[seg, 0], m.nodes[seg, 1], "-", color="tab:red", lw=0.8)
        ax.set_aspect("equal")
        ax.set_title(f"{label} (NP={m.n_nodes:,})")
        ax.set_xlabel("lon (deg)")
        ax.grid(True, lw=0.3, color="0.9")
    axes[0].set_ylabel("lat (deg)")
    fig.suptitle("Island-area filter sweep", y=1.02)
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def main() -> None:
    if not DEM.exists() or not COASTLINE.exists():
        raise SystemExit("DEM or coastline not found")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    runs = {
        "none":          (F14_NONE,  0.0,           0.0),
        "islands_1ha":   (F14_1HA,   1_000_000.0,   1.0e4),    # 1 km^2 / 1 ha
        "islands_10ha":  (F14_10HA,  1_000_000.0,   1.0e5),    # 1 km^2 / 10 ha
    }
    walls: dict[str, float] = {}
    for label, (out, mp, mi) in runs.items():
        walls[label] = _run(label, out, min_poly=mp, min_island=mi)

    meshes = {label: read_fort14(out) for label, (out, _mp, _mi) in runs.items()}
    summaries = [
        _summary(label, mesh, walls[label]) for label, mesh in meshes.items()
    ]
    text = "\n\n".join(summaries) + "\n\n"
    text += "[Reference for comparison] NP=95,551 NE=182,603 open 1/193 land 54/8,464\n"
    print(text)
    SUMMARY_TXT.write_text(text, encoding="utf-8")
    print(f"[13] wrote {SUMMARY_TXT}")

    plot_meshes(meshes, MESHES_PNG)
    print(f"[13] wrote {MESHES_PNG}")

    for label, mesh in meshes.items():
        assert (signed_areas(mesh) > 0).all(), f"{label}: flipped triangles"


if __name__ == "__main__":
    main()
