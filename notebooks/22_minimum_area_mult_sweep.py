"""PoC #22: --om-minimum-area-mult sweep on the Tokyo Bay shoreline.

Counts how many inner-shoreline polygons survive ``om.Shoreline`` with
different ``minimum_area_mult`` values, without running DistMesh.
This is the cheap way to demonstrate that the new
``--om-minimum-area-mult`` flag actually controls the islet density
that ends up in the eventual mesh.

Outputs:
    outputs/22_min_area_mult_sweep.txt
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import oceanmesh as om

from fvcom_mesh_tools.mesh_engine.oceanmesh import _stage_shapefile

REPO_ROOT = Path(__file__).resolve().parent.parent
COASTLINE = (
    REPO_ROOT / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
SUMMARY_TXT = REPO_ROOT / "outputs" / "22_min_area_mult_sweep.txt"
BBOX = (139.5650, 140.1716, 35.1015, 35.8561)
HMIN_M = 200.0
LAT_MID = 0.5 * (BBOX[2] + BBOX[3])
HMIN_DEG = HMIN_M / max(110_574.0, 111_320.0 * float(np.cos(np.deg2rad(LAT_MID))))


def main() -> None:
    if not COASTLINE.exists():
        raise SystemExit(f"missing coastline: {COASTLINE}")
    SUMMARY_TXT.parent.mkdir(parents=True, exist_ok=True)

    region = om.Region(extent=BBOX, crs=4326)

    import tempfile

    rows = []
    rows.append(
        "minimum_area_mult sweep on Tokyo Bay (MLIT C23, hmin=200 m)"
    )
    rows.append(
        f"  hmin = {HMIN_M:g} m -> h0 = {HMIN_DEG:.6f} deg "
        f"-> threshold = mult * h0^2 (deg^2)"
    )
    rows.append("")
    rows.append(f"{'mult':>6}  {'#inner':>8}  {'#mainland':>10}  threshold (deg^2)")
    rows.append("-" * 48)

    for mult in (1.0, 4.0, 25.0, 100.0, 500.0, 2000.0):
        with tempfile.TemporaryDirectory(prefix="poc22_") as td:
            staged = _stage_shapefile(COASTLINE, Path(td))
            shore = om.Shoreline(
                str(staged), region, HMIN_DEG, minimum_area_mult=mult,
            )
            n_inner = sum(
                1
                for poly in np.split(
                    shore.inner,
                    np.where(np.isnan(shore.inner[:, 0]))[0] + 1,
                )
                if poly.size > 0
            ) if shore.inner.size > 0 else 0
            n_main = sum(
                1
                for poly in np.split(
                    shore.mainland,
                    np.where(np.isnan(shore.mainland[:, 0]))[0] + 1,
                )
                if poly.size > 0
            ) if shore.mainland.size > 0 else 0
            threshold = mult * (HMIN_DEG ** 2)
            rows.append(f"{mult:>6.1f}  {n_inner:>8d}  {n_main:>10d}  {threshold:.3e}")

    summary = "\n".join(rows)
    print(summary)
    SUMMARY_TXT.write_text(summary + "\n", encoding="utf-8")
    print(f"[22] wrote {SUMMARY_TXT}")


if __name__ == "__main__":
    main()
