"""PoC #36: ``--om-max-iter`` sweep on Tokyo Bay.

Baseline (``--om-max-iter 50``) is PoC #19's setup; this PoC repeats
the build at progressively smaller iteration caps to characterise the
speed / quality trade-off. The headline question: at what iteration
count does oceanmesh's wall-clock match the deprecated ocsmesh+gmsh
draft turnaround (~40 s on the same input), and what does that cost
in mesh quality?

The result feeds the decision on whether ``--engine ocsmesh`` can be
removed in favour of "oceanmesh with reduced iterations for draft
work" (the rationale already in ``docs/engine_complementarity.md``).

All settings other than ``--om-max-iter`` are held fixed and match
PoC #19 (river-inflow points included, ``--om-seed 0`` for
deterministic runs, ``--quality-pass 0`` so the post-DistMesh
quality pass does not mask iteration-count effects, no perpfix).

Outputs:
    outputs/36_om_iter_<N>.14
    outputs/36_om_max_iter_summary.txt
    outputs/36_om_max_iter_summary.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.quality import compute_metrics

REPO = Path(__file__).resolve().parent.parent
DEM = REPO / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
RIVERS = REPO / "data" / "rivers" / "tokyo_bay" / "tokyo_bay_rivers.csv"
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "36_om_max_iter_summary.txt"
SUMMARY_JSON = OUT_DIR / "36_om_max_iter_summary.json"

HMIN_M = 200.0
HMAX_M = 5000.0

ITER_VALUES = [50, 25, 10, 5]


def _build(out_path: Path, iters: int) -> tuple[float, int]:
    """Run fmesh-buildmesh with --om-max-iter ``iters``. Returns
    ``(wall_seconds, exit_code)``."""
    args = [
        str(DEM), str(out_path),
        "--engine", "oceanmesh",
        "--hmin", str(HMIN_M), "--hmax", str(HMAX_M), "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE),
        "--om-slope-parameter", "20",
        "--om-gradation", "0.15",
        "--om-max-iter", str(iters),
        "--om-seed", "0",
        "--open-merge-coast-gap", "50",
        "--river-inflow-points", str(RIVERS),
        "--river-segment-nodes", "5",
        "--river-ibtype", "21",
        "--river-snap-tol-m", "5000",
        "--quality-pass", "0",
        "--refine-min-angle", "0",
        "--land-ibtype", "20",
        "--perpfix-iters", "1",
    ]
    print(f"[36] iters={iters}  → {out_path.name}")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    wall = time.perf_counter() - t0
    print(f"[36] iters={iters}  wall={wall:.1f} s  rc={rc}")
    return wall, rc


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not (DEM.exists() and COASTLINE.exists() and RIVERS.exists()):
        raise SystemExit("required inputs missing")

    runs: list[dict] = []
    for iters in ITER_VALUES:
        out_path = OUT_DIR / f"36_om_iter_{iters:02d}.14"
        wall, rc = _build(out_path, iters)
        if rc != 0:
            print(f"[36] iters={iters} build failed; skipping metrics")
            runs.append({"iters": iters, "wall": wall, "rc": rc})
            continue
        mesh = read_fort14(out_path)
        metrics = compute_metrics(mesh)
        runs.append({
            "iters": iters,
            "wall": wall,
            "rc": rc,
            "path": str(out_path.resolve()),
            "metrics": metrics,
        })
        print(
            f"[36] iters={iters}  NP={metrics['n_nodes']:,}  "
            f"NE={metrics['n_elements']:,}  "
            f"alpha={metrics['alpha_mean']:.4f}  "
            f"frac<20°={metrics['frac_lt_20deg']:.4%}  "
            f"max_v={metrics['max_valence']}"
        )

    # Side-by-side table for the human reader.
    headers = ["iters", "wall_s", "NP", "NE", "alpha", "alpha_p05",
               "min_p05°", "frac<20°", "max_v", "n_overconn", "n_flipped"]
    lines = [
        "PoC #36 --om-max-iter sweep on Tokyo Bay (PoC #19 settings)",
        "",
        "  " + "  ".join(h.rjust(11) for h in headers),
        "  " + "  ".join("-" * 11 for _ in headers),
    ]
    for r in runs:
        if r["rc"] != 0:
            lines.append(
                f"  {r['iters']:>11d}  {r['wall']:>11.1f}  "
                + "  ".join("FAILED".rjust(11) for _ in headers[2:])
            )
            continue
        m = r["metrics"]
        lines.append("  " + "  ".join([
            f"{r['iters']:>11d}",
            f"{r['wall']:>11.1f}",
            f"{m['n_nodes']:>11,}",
            f"{m['n_elements']:>11,}",
            f"{m['alpha_mean']:>11.4f}",
            f"{m['alpha_p05']:>11.4f}",
            f"{m['min_angle_p05_deg']:>11.2f}",
            f"{m['frac_lt_20deg'] * 100:>10.4f}%",
            f"{m['max_valence']:>11d}",
            f"{m['n_overconnected']:>11d}",
            f"{m['n_flipped']:>11d}",
        ]))

    SUMMARY_JSON.write_text(json.dumps({"runs": runs}, indent=2),
                            encoding="utf-8")
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print(f"\nwrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
