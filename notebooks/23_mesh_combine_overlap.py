"""PoC #23: nested-resolution combine via fmesh-mesh-combine --strategy overlap.

Validates the OCSMesh-backed ``--strategy overlap`` path by stitching a
coarse Tokyo Bay outer (hmin=1000 m) with a fine northern-bay inner
(hmin=200 m, bbox 139.78-140.0 x 35.55-35.75) into a single fort.14.
The merge wraps ``ocsmesh.ops.merge_overlapping_meshes``: the first
input is the *background* (coarse), subsequent inputs are
*foregrounds* (fine) that get carved into the background and the
seam is re-meshed.

PoC #21 covered ``--strategy disjoint``. This is the missing
complementary PoC for the OCSMesh-backed ``overlap`` strategy that
the CLI exposed but had no real-data validation.

Notes:
- The overlap strategy uses OCSMesh ``MeshData`` internally, which
  carries no boundary structure; the combined output therefore has
  0 open and 0 land boundaries. Re-classify with the
  ``fmesh-buildmesh``-style bbox classifier downstream if FVCOM
  ingestion needs them.
- Both inputs are built without the open-boundary perpfix (the
  combine throws boundaries away anyway), so the per-input
  generation times here are *not* directly comparable with PoC #19.

Outputs:
    outputs/23_tokyo_inner_dem.tif
    outputs/23_tokyo_outer.14
    outputs/23_tokyo_inner.14
    outputs/23_overlap_combined.14
    outputs/23_overlap_summary.txt
    outputs/23_overlap_mesh.png
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
from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main  # noqa: E402
from fvcom_mesh_tools.cli.meshcombine import main as meshcombine_main  # noqa: E402
from fvcom_mesh_tools.dem.subset import to_geotiff  # noqa: E402
from fvcom_mesh_tools.io import read_fort14  # noqa: E402
from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO_ROOT / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)

OUT_DIR = REPO_ROOT / "outputs"
INNER_DEM = OUT_DIR / "23_tokyo_inner_dem.tif"
OUTER_F14 = OUT_DIR / "23_tokyo_outer.14"
INNER_F14 = OUT_DIR / "23_tokyo_inner.14"
COMBINED_F14 = OUT_DIR / "23_overlap_combined.14"
SUMMARY_TXT = OUT_DIR / "23_overlap_summary.txt"
MESH_PNG = OUT_DIR / "23_overlap_mesh.png"

# Inner sub-region: northern Tokyo Bay, where Sumida / Arakawa /
# Edogawa enter. ~22 km x 22 km.
INNER_BBOX = (139.78, 35.55, 140.0, 35.75)

# Resolution levels. Outer:inner ratio is 5x; the seam re-mesh has
# to grade between them.
HMIN_OUTER = 1000.0
HMAX_OUTER = 10000.0
HMIN_INNER = 200.0
HMAX_INNER = 1000.0

EARTH_R_M = 6_371_000.0


def edge_lengths_m(mesh) -> np.ndarray:
    """Approximate edge length in metres on a sphere."""
    nodes = mesh.nodes
    elems = mesh.elements
    pairs = np.vstack([
        elems[:, [0, 1]],
        elems[:, [1, 2]],
        elems[:, [2, 0]],
    ])
    pairs.sort(axis=1)
    pairs = np.unique(pairs, axis=0)
    p1 = nodes[pairs[:, 0]]
    p2 = nodes[pairs[:, 1]]
    lat_mid_rad = np.deg2rad(0.5 * (p1[:, 1] + p2[:, 1]))
    dx = (
        (p2[:, 0] - p1[:, 0])
        * np.cos(lat_mid_rad)
        * np.deg2rad(1.0)
        * EARTH_R_M
    )
    dy = (p2[:, 1] - p1[:, 1]) * np.deg2rad(1.0) * EARTH_R_M
    return np.hypot(dx, dy)


def plot(combined, png: Path) -> None:
    fig, axes = plt.subplots(
        1, 2, figsize=(13, 6), dpi=120,
        gridspec_kw={"width_ratios": [1.5, 1]},
    )
    ax = axes[0]
    ax.triplot(
        combined.nodes[:, 0], combined.nodes[:, 1], combined.elements,
        color="0.6", lw=0.15,
    )
    minlon, minlat, maxlon, maxlat = INNER_BBOX
    ax.plot(
        [minlon, maxlon, maxlon, minlon, minlon],
        [minlat, minlat, maxlat, maxlat, minlat],
        "-", color="tab:red", lw=1.4, label="inner bbox",
    )
    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title(
        f"PoC #23 overlap combine  "
        f"NP={combined.n_nodes:,} NE={combined.n_elements:,}"
    )
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    el = edge_lengths_m(combined)
    ax.hist(np.log10(el), bins=60, color="0.5")
    ax.axvline(
        np.log10(HMIN_INNER), color="tab:blue", ls="--", lw=1,
        label=f"inner hmin={HMIN_INNER:.0f} m",
    )
    ax.axvline(
        np.log10(HMIN_OUTER), color="tab:orange", ls="--", lw=1,
        label=f"outer hmin={HMIN_OUTER:.0f} m",
    )
    ax.axvline(
        np.log10(HMAX_OUTER), color="tab:green", ls="--", lw=1,
        label=f"outer hmax={HMAX_OUTER:.0f} m",
    )
    ax.set_xlabel("log10 edge length (m)")
    ax.set_ylabel("edge count")
    ax.set_title("edge-length distribution (log10 m)")
    ax.legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)


def build(
    dem: Path, hmin: float, hmax: float, out_f14: Path, *, label: str,
) -> float:
    """Build a mesh via ``fmesh-buildmesh --engine oceanmesh``; return wall."""
    args = [
        str(dem), str(out_f14),
        "--engine", "oceanmesh",
        "--hmin", str(hmin), "--hmax", str(hmax), "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE),
        "--om-seed", "0",
        "--om-max-iter", "50",
        "--quality-pass", "0",
        "--refine-min-angle", "0",
        "--no-perpfix",
        "--quiet",
    ]
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    dt = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh ({label}) exited {rc}")
    return dt


def main() -> None:
    for p in (DEM, COASTLINE):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[23] subset DEM -> {INNER_DEM.name}  bbox={INNER_BBOX}")
    info = to_geotiff(DEM, INNER_DEM, INNER_BBOX)
    print(f"[23]   shape={info['shape']}  crs={info['crs']}")

    print(
        f"[23] build outer (background)  "
        f"hmin={HMIN_OUTER:.0f}  hmax={HMAX_OUTER:.0f} ..."
    )
    t_outer = build(DEM, HMIN_OUTER, HMAX_OUTER, OUTER_F14, label="outer")
    outer = read_fort14(OUTER_F14)
    print(
        f"[23]   outer NP={outer.n_nodes:,} NE={outer.n_elements:,} "
        f"({t_outer:.2f} s)"
    )

    print(
        f"[23] build inner (foreground)  "
        f"hmin={HMIN_INNER:.0f}  hmax={HMAX_INNER:.0f} ..."
    )
    t_inner = build(INNER_DEM, HMIN_INNER, HMAX_INNER, INNER_F14, label="inner")
    inner = read_fort14(INNER_F14)
    print(
        f"[23]   inner NP={inner.n_nodes:,} NE={inner.n_elements:,} "
        f"({t_inner:.2f} s)"
    )

    print("[23] fmesh-mesh-combine --strategy overlap (defaults) ...")
    t0 = time.perf_counter()
    rc = meshcombine_main([
        str(OUTER_F14), str(INNER_F14), str(COMBINED_F14),
        "--strategy", "overlap",
        "--title", "PoC #23 Tokyo Bay coarse + northern fine (overlap)",
    ])
    t_combine = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"fmesh-mesh-combine exited {rc}")
    combined = read_fort14(COMBINED_F14)

    flipped = int((signed_areas(combined) <= 0).sum())
    q = alpha_quality(combined)
    a = min_interior_angle(combined)
    el = edge_lengths_m(combined)
    el_outer = edge_lengths_m(outer)
    el_inner = edge_lengths_m(inner)

    summary_lines = [
        (
            f"input outer (background): {OUTER_F14.name}  "
            f"NP={outer.n_nodes:,} NE={outer.n_elements:,}  "
            f"({t_outer:.2f} s)"
        ),
        (
            f"input inner (foreground): {INNER_F14.name}  "
            f"NP={inner.n_nodes:,} NE={inner.n_elements:,}  "
            f"({t_inner:.2f} s)"
        ),
        (
            f"strategy: overlap (defaults: buffer_size=0.0075, "
            f"buffer_domain=0.002, min_int_ang=30, adjacent_layers=0)  "
            f"combine wall: {t_combine:.2f} s"
        ),
        "",
        (
            f"NP={combined.n_nodes:,}  NE={combined.n_elements:,}  "
            f"flipped={flipped}"
        ),
        (
            "  open boundaries (overlap drops these): "
            f"{len(combined.open_boundaries)}"
        ),
        (
            "  land boundaries (overlap drops these): "
            f"{len(combined.land_boundaries)}"
        ),
        "",
        "[Quality (combined)]",
        f"  alpha mean      : {q.mean():.4f}",
        f"  alpha p10       : {np.percentile(q, 10):.4f}",
        f"  min-angle p50   : {np.percentile(a, 50):.2f} deg",
        f"  frac<20deg      : {(a < 20).mean() * 100:.2f} %",
        "",
        "[Edge length, m, p5/p25/p50/p75/p95]",
        (
            "  outer    : "
            f"{np.percentile(el_outer, 5):.0f} / "
            f"{np.percentile(el_outer, 25):.0f} / "
            f"{np.percentile(el_outer, 50):.0f} / "
            f"{np.percentile(el_outer, 75):.0f} / "
            f"{np.percentile(el_outer, 95):.0f}"
        ),
        (
            "  inner    : "
            f"{np.percentile(el_inner, 5):.0f} / "
            f"{np.percentile(el_inner, 25):.0f} / "
            f"{np.percentile(el_inner, 50):.0f} / "
            f"{np.percentile(el_inner, 75):.0f} / "
            f"{np.percentile(el_inner, 95):.0f}"
        ),
        (
            "  combined : "
            f"{np.percentile(el, 5):.0f} / "
            f"{np.percentile(el, 25):.0f} / "
            f"{np.percentile(el, 50):.0f} / "
            f"{np.percentile(el, 75):.0f} / "
            f"{np.percentile(el, 95):.0f}"
        ),
    ]
    summary = "\n".join(summary_lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[23] wrote {SUMMARY_TXT}")

    plot(combined, MESH_PNG)
    print(f"[23] wrote {MESH_PNG}")

    if flipped > 0:
        raise SystemExit(f"{flipped} flipped triangles in combined mesh")
    if combined.n_elements == 0:
        raise SystemExit("combined mesh is empty")


if __name__ == "__main__":
    main()
