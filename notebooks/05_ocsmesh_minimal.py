"""PoC #5: minimal OCSMesh pipeline -> fort.14.

End-to-end smoke test of an OCSMesh-driven mesh-generation pipeline that
mirrors what OceanMesh2D does in MATLAB at the coarsest possible level:

    DEM raster -> Geom (zmax<=0 land/water mask)
                -> Hfun (uniform hmin/hmax in metres)
                -> MeshDriver(engine="gmsh")
                -> EuclideanMesh2D
                -> fort.14 (via OCSMesh GRD writer)

The output is then re-read with our own :func:`read_fort14` to verify the
file is parseable end-to-end and to record NP / NE / depth statistics. No
boundary classification is performed here: that is a known feature gap to
be quantified in PoC #6 / Phase 4.

Outputs:
    outputs/05_tokyo_bay_minimal.14
    outputs/05_ocsmesh_minimal_summary.txt
    outputs/05_ocsmesh_minimal_mesh.png
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from fvcom_mesh_tools.io import read_fort14  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM_PATH = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
OUT_DIR = REPO_ROOT / "outputs"
FORT14_OUT = OUT_DIR / "05_tokyo_bay_minimal.14"
SUMMARY_TXT = OUT_DIR / "05_ocsmesh_minimal_summary.txt"
MESH_PNG = OUT_DIR / "05_ocsmesh_minimal_mesh.png"

# Coarse uniform sizing (metres). The DEM CRS is geographic (EPSG:4326), but
# OCSMesh treats hmin/hmax in the same units as Hfun's output mesh; the
# defaults inside OCSMesh handle the lon/lat -> metric conversion.
HMIN_M = 200.0
HMAX_M = 5000.0
ZMAX = 0.0  # everything at or below msl is meshed (water domain)


def _depth_stats(name: str, x: np.ndarray) -> str:
    return (
        f"{name}: n={x.size:,}  min={np.nanmin(x):.4f}  "
        f"p50={np.nanpercentile(x, 50):.4f}  p95={np.nanpercentile(x, 95):.4f}  "
        f"max={np.nanmax(x):.4f}  mean={np.nanmean(x):.4f}"
    )


def plot_mesh(nodes: np.ndarray, elements: np.ndarray, depths: np.ndarray, png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 7), dpi=120)
    tri = ax.tripcolor(
        nodes[:, 0], nodes[:, 1], elements,
        facecolors=depths[elements].mean(axis=1),
        cmap="viridis_r", edgecolors="0.85", linewidth=0.05,
    )
    fig.colorbar(tri, ax=ax, label="depth (m, +down)")
    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title("OCSMesh minimal pipeline: Tokyo Bay")
    fig.tight_layout()
    fig.savefig(png, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not DEM_PATH.exists():
        raise SystemExit(f"DEM not found: {DEM_PATH}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # OCSMesh imports are deferred so the rest of the package stays
    # importable in environments without OCSMesh installed.
    from ocsmesh import Geom, Hfun, MeshDriver, Raster

    print(f"[05] DEM:  {DEM_PATH}")
    print(f"[05] hmin={HMIN_M:g} m  hmax={HMAX_M:g} m  zmax={ZMAX:g} m")

    t0 = time.perf_counter()
    raster = Raster(str(DEM_PATH))
    geom = Geom(raster, zmax=ZMAX)
    hfun = Hfun(raster, hmin=HMIN_M, hmax=HMAX_M)
    driver = MeshDriver(geom, hfun=hfun, engine_name="gmsh")
    mesh = driver.run()
    t_gen = time.perf_counter() - t0
    print(f"[05] mesh generation: {t_gen:.2f} s")

    # Write GRD (= fort.14) directly with OCSMesh.
    if FORT14_OUT.exists():
        FORT14_OUT.unlink()
    mesh.write(str(FORT14_OUT), format="grd", overwrite=True)
    print(f"[05] wrote {FORT14_OUT}")

    # Round-trip through our own parser to verify parseability.
    f14 = read_fort14(FORT14_OUT)
    print(f"[05] re-read NP={f14.n_nodes:,}  NE={f14.n_elements:,}")

    n_open_segs = len(f14.open_boundaries)
    n_open_nodes = sum(len(b) for b in f14.open_boundaries)
    n_land_segs = len(f14.land_boundaries)
    n_land_nodes = sum(len(ids) for _, ids in f14.land_boundaries)

    lines = [
        f"DEM:  {DEM_PATH}",
        f"hmin={HMIN_M:g} m  hmax={HMAX_M:g} m  zmax={ZMAX:g} m",
        "engine: gmsh",
        f"mesh generation wall time: {t_gen:.2f} s",
        "",
        "[Mesh size]",
        f"  NP={f14.n_nodes:,}  NE={f14.n_elements:,}",
        "",
        "[Bounding box (lon/lat deg)]",
        f"  xmin={f14.bbox[0]:.6f}  ymin={f14.bbox[1]:.6f}  "
        f"xmax={f14.bbox[2]:.6f}  ymax={f14.bbox[3]:.6f}",
        "",
        "[Depth column (sign as written by OCSMesh)]",
        _depth_stats("  depth", f14.depths),
        "",
        "[Boundaries written by OCSMesh GRD writer]",
        f"  open  segments={n_open_segs}  total nodes={n_open_nodes}",
        f"  land  segments={n_land_segs}  total nodes={n_land_nodes}",
        "",
        "[NOTE]",
        "  Boundary classification is not performed by the minimal Geom/Hfun "
        "pipeline.",
        "  PoC #6 quantifies the parity gap vs the reference fort.14 and "
        "Phase 4 lists",
        "  the missing capabilities (boundary tagging, channel widening, "
        "river inflow nodes, etc.).",
    ]
    summary = "\n".join(lines)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[05] wrote {SUMMARY_TXT}")

    plot_mesh(f14.nodes, f14.elements, f14.depths, MESH_PNG)
    print(f"[05] wrote {MESH_PNG}")


if __name__ == "__main__":
    main()
