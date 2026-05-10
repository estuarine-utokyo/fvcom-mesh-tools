"""PoC #34: ``--om-wavelength-sizing`` A/B comparison.

Generates two Tokyo-Bay meshes with identical settings except for
``--om-wavelength-sizing``: PoC #19's gradient-only baseline vs.
the new wavelength + gradient combination. Reports

    * unified quality metrics (alpha, frac<20°, max_valence, ...)
      via ``fvcom_mesh_tools.quality.compute_metrics``;
    * NP / NE diff;
    * shallow-cell density: count of triangles whose centroid lies
      in <= 5 m of water, normalised by area — quantifies whether
      the wavelength-driven sizing actually refines the shoaling
      regions;
    * minimum CFL-feasible dt at C=0.7 for the worst-case (smallest
      dx / largest celerity) triangle.

The two presets:

    * baseline: feature + bathymetric_gradient (PoC #19 settings).
    * wavelength: + wavelength_sizing_function with
      ``period=44712 s (M2)``, ``wl=100`` (so the implied dt is
      ~7.5 min — a comfortable FVCOM time step).

Both are identical in everything else (hmin / hmax / coastline /
gradation / seed) so the only confounder is the new sizing.

Outputs:
    outputs/34_tokyo_bay_baseline.14
    outputs/34_tokyo_bay_wavelength.14
    outputs/34_wavelength_sizing_summary.txt
    outputs/34_wavelength_sizing_summary.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

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
F14_BASELINE = OUT_DIR / "34_tokyo_bay_baseline.14"
F14_WAVELENGTH = OUT_DIR / "34_tokyo_bay_wavelength.14"
SUMMARY_TXT = OUT_DIR / "34_wavelength_sizing_summary.txt"
SUMMARY_JSON = OUT_DIR / "34_wavelength_sizing_summary.json"

HMIN_M = 200.0
HMAX_M = 5000.0


def _common_args(out_path: Path) -> list[str]:
    """The PoC #19-equivalent argv (without --om-wavelength-sizing)."""
    return [
        str(DEM), str(out_path),
        "--engine", "oceanmesh",
        "--hmin", str(HMIN_M), "--hmax", str(HMAX_M), "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE),
        "--river-inflow-points", str(RIVERS),
        "--river-segment-nodes", "5",
        "--river-ibtype", "21",
        "--river-snap-tol-m", "5000",
        "--open-merge-coast-gap", "50",
        "--om-seed", "0",
        "--om-max-iter", "50",
        "--quality-pass", "0",  # leave quality pass off so we measure the
                                 # sizing-driven mesh only
        "--perpfix-iters", "1",
        "--land-ibtype", "20",
    ]


def _run_buildmesh(out_path: Path, *, wavelength: bool) -> float:
    args = _common_args(out_path)
    if wavelength:
        args += [
            "--om-wavelength-sizing",
            "--om-wavelength-period", "44712.0",     # M2
            "--om-wavelength-grid-spacing", "100",   # implied dt ~ 7.5 min
        ]
    print(f"[34] running fmesh-buildmesh -> {out_path.name} "
          f"(wavelength={wavelength})")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    wall = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh exit {rc}")
    print(f"[34] {out_path.name} wall: {wall:.1f} s")
    return wall


def _shallow_cell_metrics(mesh) -> dict[str, float]:
    """Count triangles whose centroid sits in <=5 m of water."""
    if mesh.n_elements == 0:
        return {"n_shallow_le_5m": 0, "frac_shallow_le_5m": 0.0,
                "min_cfl_dt_s": 0.0}
    # Centroid depth = mean of element-vertex depths (positive down).
    centroid_depths = mesh.depths[mesh.elements].mean(axis=1)
    shallow_mask = centroid_depths <= 5.0
    n_shallow = int(shallow_mask.sum())
    frac_shallow = float(n_shallow / mesh.n_elements)

    # Crude min-CFL estimate at C=0.7. dx in metres needs lat conversion;
    # use median-per-element edge length and centre-of-mass latitude.
    # We compute |edges| in metres at the local latitude of each centroid.
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    centroids = (p0 + p1 + p2) / 3.0
    lat0 = centroids[:, 1]
    deg2m_lat = 110_574.0
    deg2m_lon = 111_320.0 * np.cos(np.deg2rad(lat0))
    def _len_m(a, b):
        dlon = (b[:, 0] - a[:, 0]) * deg2m_lon
        dlat = (b[:, 1] - a[:, 1]) * deg2m_lat
        return np.sqrt(dlon**2 + dlat**2)
    e01 = _len_m(p0, p1)
    e12 = _len_m(p1, p2)
    e20 = _len_m(p2, p0)
    dx_min_m = np.minimum(np.minimum(e01, e12), e20)
    h_eff = np.maximum(np.abs(centroid_depths), 0.5)  # avoid 0
    celerity = np.sqrt(9.81 * h_eff)
    cfl_dt = 0.7 * dx_min_m / celerity
    return {
        "n_shallow_le_5m": n_shallow,
        "frac_shallow_le_5m": frac_shallow,
        "min_cfl_dt_s": float(np.percentile(cfl_dt, 1)),  # 1 %tile worst
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not (DEM.exists() and COASTLINE.exists() and RIVERS.exists()):
        raise SystemExit("required inputs missing")

    wall_b = _run_buildmesh(F14_BASELINE, wavelength=False)
    wall_w = _run_buildmesh(F14_WAVELENGTH, wavelength=True)

    mb = read_fort14(F14_BASELINE)
    mw = read_fort14(F14_WAVELENGTH)
    qb = compute_metrics(mb)
    qw = compute_metrics(mw)
    sb = _shallow_cell_metrics(mb)
    sw = _shallow_cell_metrics(mw)

    print()
    print("=== quality (baseline vs wavelength) ===")
    keys_q = [
        "n_nodes", "n_elements",
        "alpha_mean", "alpha_p05", "alpha_p50",
        "min_angle_p05_deg", "min_angle_p50_deg",
        "frac_lt_20deg", "max_valence", "n_overconnected",
        "n_flipped", "n_components", "n_disjoint_elems",
    ]
    keys_s = ["n_shallow_le_5m", "frac_shallow_le_5m", "min_cfl_dt_s"]
    for k in keys_q:
        print(f"  {k:<22}  {qb[k]!s:>16}  {qw[k]!s:>16}")
    for k in keys_s:
        print(f"  {k:<22}  {sb[k]!s:>16}  {sw[k]!s:>16}")
    print(f"  {'wall_seconds':<22}  {wall_b:>16.1f}  {wall_w:>16.1f}")

    payload = {
        "baseline": {"path": str(F14_BASELINE.resolve()),
                     "wall_seconds": wall_b,
                     "metrics": qb, "shallow": sb},
        "wavelength": {"path": str(F14_WAVELENGTH.resolve()),
                       "wall_seconds": wall_w,
                       "metrics": qw, "shallow": sw},
        "delta": {
            **{k: float(qw[k] - qb[k]) for k in keys_q
               if isinstance(qb[k], (int, float))
               and not (isinstance(qb[k], float) and np.isnan(qb[k]))},
            **{k: float(sw[k] - sb[k]) for k in keys_s},
            "wall_seconds": float(wall_w - wall_b),
        },
        "config": {
            "hmin_m": HMIN_M, "hmax_m": HMAX_M,
            "wavelength_period_s": 44712.0,
            "wavelength_grid_spacing": 100,
        },
    }
    SUMMARY_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    summary_lines = [
        f"PoC #34 wavelength_sizing A/B on Tokyo Bay (hmin={HMIN_M:g} m).",
        "",
        f"{'metric':<22}  {'baseline':>16}  {'wavelength':>16}",
        "  " + "-" * 56,
    ]
    for k in keys_q + keys_s:
        summary_lines.append(f"  {k:<22}  {qb.get(k, 0)!s:>16}  "
                             f"{qw.get(k, 0)!s:>16}")
    summary_lines.append(
        f"  {'wall_seconds':<22}  {wall_b:>16.1f}  {wall_w:>16.1f}"
    )
    SUMMARY_TXT.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(f"\nwrote {SUMMARY_TXT}\nwrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
