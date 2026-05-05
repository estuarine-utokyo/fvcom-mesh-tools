"""PoC #11: ``fmesh-buildmesh --quality-pass 6`` end-to-end smoke test.

Verifies that the CLI integration of edge-swap + smoothing produces
the same mesh quality the standalone PoC #10 driver produced when
applied to the PoC #7 output. Both should land at frac<20deg ~ 17 %
and mean alpha ~ 0.75.

Outputs:
    outputs/11_tokyo_bay_quality_cli.14
    outputs/11_buildmesh_quality_summary.txt
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from fvcom_mesh_tools.algorithms import (
    alpha_quality,
    min_interior_angle,
    open_bdy_perpendicularity,
    signed_areas,
)
from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main
from fvcom_mesh_tools.io import read_fort14

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
OUT_DIR = REPO_ROOT / "outputs"
FORT14 = OUT_DIR / "11_tokyo_bay_quality_cli.14"
SUMMARY = OUT_DIR / "11_buildmesh_quality_summary.txt"

HMIN_M = 200.0
HMAX_M = 5000.0
QUALITY_ROUNDS = 6


def _stats(name: str, x: np.ndarray, fmt: str = ".4f") -> str:
    if x.size == 0:
        return f"{name}: empty"
    return (
        f"{name}: n={x.size:,}  "
        f"min={np.nanmin(x):{fmt}}  "
        f"p05={np.nanpercentile(x, 5):{fmt}}  "
        f"p50={np.nanpercentile(x, 50):{fmt}}  "
        f"p95={np.nanpercentile(x, 95):{fmt}}  "
        f"max={np.nanmax(x):{fmt}}  "
        f"mean={np.nanmean(x):{fmt}}"
    )


def main() -> None:
    if not DEM.exists():
        raise SystemExit(f"DEM not found: {DEM}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[11] running fmesh-buildmesh --quality-pass {QUALITY_ROUNDS}...")
    t0 = time.perf_counter()
    rc = buildmesh_main([
        str(DEM), str(FORT14),
        "--hmin", str(HMIN_M),
        "--hmax", str(HMAX_M),
        "--zmax", "0.0",
        "--interp-method", "linear",
        "--land-ibtype", "20",
        "--quality-pass", str(QUALITY_ROUNDS),
        "--smooth-iters", "5",
        "--smooth-alpha", "0.5",
        "--perpfix-iters", "1",
    ])
    t_total = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"fmesh-buildmesh exited {rc}")
    print(f"[11] wall time: {t_total:.2f} s")

    f14 = read_fort14(FORT14)
    q = alpha_quality(f14)
    a = min_interior_angle(f14)
    perp = open_bdy_perpendicularity(f14, segment_index=0)
    flipped = int((signed_areas(f14) <= 0).sum())

    lines = [
        f"DEM:    {DEM}",
        f"output: {FORT14}",
        f"quality_pass={QUALITY_ROUNDS}  smooth_iters=5  smooth_alpha=0.5",
        f"wall time: {t_total:.2f} s",
        "",
        f"NP={f14.n_nodes:,}  NE={f14.n_elements:,}",
        f"flipped triangles: {flipped:,}",
        "",
        "[Triangle alpha-quality (1 = equilateral)]",
        _stats("  alpha", q),
        f"  frac alpha < 0.3: {(q < 0.3).mean() * 100:.2f} %",
        "",
        "[Triangle minimum interior angle (deg)]",
        _stats("  min-angle", a, ".2f"),
        f"  frac min-angle < 20 deg: {(a < 20).mean() * 100:.2f} %",
        "",
        "[Open-boundary perpendicularity (deg from 90)]",
        _stats("  perp", perp, ".4f"),
    ]
    summary = "\n".join(lines)
    print(summary)
    SUMMARY.write_text(summary + "\n", encoding="utf-8")
    print(f"[11] wrote {SUMMARY}")

    assert flipped == 0, f"{flipped} flipped triangles after CLI quality pass"
    print("[11] sanity checks PASSED")


if __name__ == "__main__":
    main()
