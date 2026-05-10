"""PoC #18: head-to-head benchmark of OCSMesh+gmsh vs oceanmesh.

Same Tokyo Bay inputs as PoC #16 (DEM, MLIT C23 coastline) and the
same target sizes (hmin=200 m, hmax=5000 m). The OCSMesh number is
read from the existing ``outputs/16_river_inflow_summary.txt``; this
script builds a comparable mesh with oceanmesh + DistMesh, applying
the OceanMesh2D-style recipe with feature_sizing_function plus
bathymetric_gradient_sizing_function.

Comparison metrics (computed identically on both, in lon/lat space):

* NP, NE
* alpha_quality (mean and frac < 0.3)
* min interior angle (p50 and frac < 20 deg)
* flipped triangles (signed area <= 0)
* wall-clock time

Outputs:
    outputs/18_oceanmesh_tokyo_bay.png
    outputs/18_oceanmesh_metrics.json
    outputs/18_comparison_summary.txt
"""

from __future__ import annotations

import json
import shutil
import tempfile
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import oceanmesh as om  # noqa: E402
import rasterio  # noqa: E402

from fvcom_mesh_tools.plotting import MESH_PNG_DPI  # noqa: E402

WGS84_PRJ = (
    'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",'
    '6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",'
    "0.017453292519943295]]"
)


def stage_shapefile_with_prj(src_shp: Path, dst_dir: Path) -> Path:
    """Copy a shapefile into ``dst_dir`` and ensure a WGS84 .prj sidecar.

    MLIT C23 ships ``C23-06_TOKYOBAY.shp`` without a ``.prj`` (only the
    OUTER/INNER split files carry one), and oceanmesh.Shoreline requires
    a CRS-tagged input. We materialise a CRS-tagged copy in a writable
    directory rather than mutating the user's data tree.
    """
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst_shp = dst_dir / src_shp.name
    for ext in (".shp", ".shx", ".dbf", ".cpg"):
        s = src_shp.with_suffix(ext)
        if s.exists():
            shutil.copyfile(s, dst_shp.with_suffix(ext))
    prj = src_shp.with_suffix(".prj")
    dst_prj = dst_shp.with_suffix(".prj")
    if prj.exists():
        shutil.copyfile(prj, dst_prj)
    else:
        dst_prj.write_text(WGS84_PRJ, encoding="utf-8")
    return dst_shp

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO_ROOT / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
OUT_DIR = REPO_ROOT / "outputs"
PNG = OUT_DIR / "18_oceanmesh_tokyo_bay.png"
METRICS_JSON = OUT_DIR / "18_oceanmesh_metrics.json"
SUMMARY_TXT = OUT_DIR / "18_comparison_summary.txt"
OCSMESH_SUMMARY = OUT_DIR / "16_river_inflow_summary.txt"

HMIN_M = 200.0
HMAX_M = 5000.0


def m_to_deg_lat(m: float) -> float:
    return m / 110_574.0


def m_to_deg_lon(m: float, lat: float) -> float:
    return m / (111_320.0 * np.cos(np.deg2rad(lat)))


def compute_metrics(points: np.ndarray, cells: np.ndarray) -> dict:
    """alpha-quality, min-interior-angle, flipped count - lon/lat planar."""
    p0 = points[cells[:, 0]]
    p1 = points[cells[:, 1]]
    p2 = points[cells[:, 2]]
    e1 = np.linalg.norm(p1 - p0, axis=1)
    e2 = np.linalg.norm(p2 - p1, axis=1)
    e3 = np.linalg.norm(p0 - p2, axis=1)
    sa = 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    )
    A = np.abs(sa)
    eps = 1e-30
    alpha = 4.0 * np.sqrt(3.0) * A / (e1 ** 2 + e2 ** 2 + e3 ** 2 + eps)
    a, b, c = e1, e2, e3
    cos_a = np.clip((b ** 2 + c ** 2 - a ** 2) / (2.0 * b * c + eps), -1.0, 1.0)
    cos_b = np.clip((a ** 2 + c ** 2 - b ** 2) / (2.0 * a * c + eps), -1.0, 1.0)
    cos_c = np.clip((a ** 2 + b ** 2 - c ** 2) / (2.0 * a * b + eps), -1.0, 1.0)
    cos_max = np.maximum(np.maximum(cos_a, cos_b), cos_c)
    min_angle = np.degrees(np.arccos(cos_max))
    return {
        "NP": int(points.shape[0]),
        "NE": int(cells.shape[0]),
        "alpha_mean": float(alpha.mean()),
        "alpha_lt_03_pct": float((alpha < 0.3).mean() * 100.0),
        "min_angle_p50_deg": float(np.percentile(min_angle, 50)),
        "frac_lt_20deg_pct": float((min_angle < 20.0).mean() * 100.0),
        "flipped": int((sa <= 0).sum()),
    }


def parse_ocsmesh_summary(path: Path) -> dict:
    """Pull NP/NE/alpha/frac20/flipped/wall from the PoC #16 summary."""
    text = path.read_text(encoding="utf-8")
    out: dict = {}
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("NP="):
            parts = s.replace(",", "").split()
            for tok in parts:
                if "=" in tok:
                    k, v = tok.split("=")
                    if k == "NP":
                        out["NP"] = int(v)
                    elif k == "NE":
                        out["NE"] = int(v)
                    elif k == "flipped":
                        out["flipped"] = int(v)
        elif "alpha mean" in s:
            out["alpha_mean"] = float(s.split(":")[-1])
        elif "alpha < 0.3" in s:
            out["alpha_lt_03_pct"] = float(s.split(":")[-1].rstrip(" %"))
        elif "min-angle p50" in s:
            out["min_angle_p50_deg"] = float(s.split(":")[-1])
        elif "min-angle < 20 deg" in s:
            out["frac_lt_20deg_pct"] = float(s.split(":")[-1].rstrip(" %"))
        elif s.startswith("wall:"):
            out["wall_s"] = float(s.split(":")[-1].rstrip(" s"))
    return out


def main() -> None:
    for p in (DEM, COASTLINE, OCSMESH_SUMMARY):
        if not p.exists():
            raise SystemExit(f"missing input: {p}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # DEM bbox = our domain
    with rasterio.open(DEM) as r:
        bb = r.bounds
        bbox = (float(bb.left), float(bb.right), float(bb.bottom), float(bb.top))
    print(f"[18] Tokyo Bay bbox (lon, lat): {bbox}")

    lat_mid = 0.5 * (bbox[2] + bbox[3])
    hmin_lat = m_to_deg_lat(HMIN_M)
    hmin_lon = m_to_deg_lon(HMIN_M, lat_mid)
    hmax_lat = m_to_deg_lat(HMAX_M)
    hmin_deg = float(min(hmin_lat, hmin_lon))
    hmax_deg = float(hmax_lat)
    print(
        f"[18] hmin={HMIN_M:.0f} m -> {hmin_deg:.6f} deg; "
        f"hmax={HMAX_M:.0f} m -> {hmax_deg:.6f} deg"
    )

    t_total0 = time.perf_counter()

    region = om.Region(extent=bbox, crs=4326)
    print("[18] reading shoreline ...")
    t0 = time.perf_counter()
    with tempfile.TemporaryDirectory(prefix="poc18_coast_") as td:
        coast_path = stage_shapefile_with_prj(COASTLINE, Path(td))
        shore = om.Shoreline(str(coast_path), region, hmin_deg)
    print(f"[18]   shoreline: {time.perf_counter() - t0:.2f} s")

    sdf = om.signed_distance_function(shore)

    print("[18] reading DEM ...")
    t0 = time.perf_counter()
    dem = om.DEM(str(DEM), bbox=region, crs=4326)
    print(f"[18]   DEM: {time.perf_counter() - t0:.2f} s")

    print("[18] building feature_sizing_function ...")
    t0 = time.perf_counter()
    edge_feat = om.feature_sizing_function(
        shore, sdf, max_edge_length=hmax_deg, crs=4326,
    )
    print(f"[18]   feature: {time.perf_counter() - t0:.2f} s")

    print("[18] building bathymetric_gradient_sizing_function ...")
    t0 = time.perf_counter()
    edge_grad = om.bathymetric_gradient_sizing_function(
        dem,
        slope_parameter=20.0,
        filter_quotient=50,
        min_edge_length=hmin_deg,
        max_edge_length=hmax_deg,
        crs=4326,
    )
    print(f"[18]   gradient: {time.perf_counter() - t0:.2f} s")

    edge = om.enforce_mesh_gradation(
        om.compute_minimum([edge_feat, edge_grad]), gradation=0.15,
    )

    print("[18] running generate_mesh (DistMesh) ...")
    t0 = time.perf_counter()
    points, cells = om.generate_mesh(sdf, edge, max_iter=50, seed=0)
    print(f"[18]   generate_mesh: {time.perf_counter() - t0:.2f} s")

    print("[18] mesh cleanup ...")
    points, cells = om.make_mesh_boundaries_traversable(points, cells)
    points, cells = om.delete_faces_connected_to_one_face(points, cells)
    points, cells = om.delete_boundary_faces(points, cells, min_qual=0.15)
    points, cells = om.laplacian2(points, cells)

    wall = time.perf_counter() - t_total0
    print(f"[18] total wall: {wall:.2f} s")

    # Ensure CCW so signed_area is consistently positive (lon=x, lat=y).
    p0 = points[cells[:, 0]]
    p1 = points[cells[:, 1]]
    p2 = points[cells[:, 2]]
    sa = 0.5 * (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p2[:, 0] - p0[:, 0]) * (p1[:, 1] - p0[:, 1])
    )
    if (sa < 0).mean() > 0.5:
        cells = cells[:, [0, 2, 1]].copy()

    om_metrics = compute_metrics(points, cells)
    om_metrics["wall_s"] = wall
    METRICS_JSON.write_text(json.dumps(om_metrics, indent=2), encoding="utf-8")

    ocs_metrics = parse_ocsmesh_summary(OCSMESH_SUMMARY)

    def fmt(d: dict, key: str, unit: str = "") -> str:
        v = d.get(key)
        return "n/a" if v is None else (
            f"{v:.4f}{unit}" if isinstance(v, float) else f"{v:,}{unit}"
        )

    rows = [
        ("NP",                  "NP",                ""),
        ("NE",                  "NE",                ""),
        ("alpha mean",          "alpha_mean",        ""),
        ("alpha < 0.3",         "alpha_lt_03_pct",   " %"),
        ("min-angle p50",       "min_angle_p50_deg", " deg"),
        ("min-angle < 20 deg",  "frac_lt_20deg_pct", " %"),
        ("flipped",             "flipped",           ""),
        ("wall",                "wall_s",            " s"),
    ]
    lines = []
    lines.append("Tokyo Bay benchmark: OCSMesh+gmsh (PoC #16) vs oceanmesh+DistMesh")
    lines.append(f"DEM:       {DEM.name}")
    lines.append(f"coastline: {COASTLINE.name}")
    lines.append(f"hmin/hmax: {HMIN_M:.0f} m / {HMAX_M:.0f} m")
    lines.append("")
    header = f"{'metric':<22}{'OCSMesh':>16}{'oceanmesh':>16}"
    lines.append(header)
    lines.append("-" * len(header))
    for label, key, unit in rows:
        lines.append(
            f"{label:<22}{fmt(ocs_metrics, key, unit):>16}{fmt(om_metrics, key, unit):>16}"
        )
    summary = "\n".join(lines)
    print()
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")

    fig, ax = plt.subplots(figsize=(9, 7), dpi=120)
    ax.triplot(points[:, 0], points[:, 1], cells, color="0.5", lw=0.15)
    ax.set_aspect("equal")
    ax.set_xlabel("lon (deg)")
    ax.set_ylabel("lat (deg)")
    ax.set_title("oceanmesh + DistMesh: Tokyo Bay")
    fig.tight_layout()
    fig.savefig(PNG, bbox_inches="tight", dpi=MESH_PNG_DPI)
    plt.close(fig)
    print(f"[18] wrote {PNG}")


if __name__ == "__main__":
    main()
