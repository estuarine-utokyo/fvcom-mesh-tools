"""PoC #39: ``--om-courant-sizing`` A/B comparison on Tokyo Bay.

Validates the new Courant-bounded sizing primitive on the same input
PoC #19 / PoC #34 used. Builds two meshes:

    * baseline: PoC #19 settings (feature + bathymetric_gradient
      sizing only).
    * courant:  + ``--om-courant-sizing --om-courant-target 0.7
                  --om-courant-timestep 10.0``.

The interesting metric is the *post-hoc* minimum CFL-feasible time
step at C=0.7. The Courant primitive sets an upper sizing envelope
``dx_max(h, dt) = c_char(h, nu) * dt / C`` so that the resulting mesh
satisfies ``C <= 0.7`` at the requested ``dt``. We expect:

  * NP / NE: similar to or slightly larger than baseline (the
    constraint can only refine, never coarsen, when ``min_edgelength``
    floors the sizing — and only refines in mid-depth water where
    ``c_char * dt / C < hmax``).
  * alpha / frac<20°: comparable to baseline (sizing is per-cell so
    no skewing artefact from the constraint itself).
  * 1 %-ile CFL-feasible dt at C=0.7: in the baseline this depends
    on whatever the gradient sizing happened to allocate. In the
    courant build it should be **>= 10 s** by construction (modulo
    DistMesh imperfections in the 1 %-ile tail).

Tokyo Bay is shallow (max ~30 m in the bay proper), so the Courant
upper envelope at ``dt=10 s`` only bites in mid-depth regions where
the baseline sizing was relatively coarse. ``dt=10 s`` is chosen
above the realistic FVCOM coastal step (~5 s) deliberately so the
constraint actively binds in this otherwise tame basin.

Outputs:
    outputs/39_tokyo_bay_courant_dt10.14
    outputs/39_courant_sizing_summary.txt
    outputs/39_courant_sizing_summary.json
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
F14_BASELINE = OUT_DIR / "19_tokyo_bay_oceanmesh.14"     # reuse PoC #19
F14_COURANT = OUT_DIR / "39_tokyo_bay_courant_dt10.14"
SUMMARY_TXT = OUT_DIR / "39_courant_sizing_summary.txt"
SUMMARY_JSON = OUT_DIR / "39_courant_sizing_summary.json"

HMIN_M = 200.0
HMAX_M = 5000.0
COURANT_TARGET = 0.7
COURANT_DT_S = 10.0
COURANT_NU_M = 2.0


def _common_args(out_path: Path) -> list[str]:
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
        "--quality-pass", "0",
        "--perpfix-iters", "1",
        "--land-ibtype", "20",
    ]


def _run_courant_build(out_path: Path) -> float:
    args = _common_args(out_path)
    args += [
        "--om-courant-sizing",
        "--om-courant-target", str(COURANT_TARGET),
        "--om-courant-timestep", str(COURANT_DT_S),
        "--om-courant-wave-amplitude", str(COURANT_NU_M),
    ]
    print(f"[39] running fmesh-buildmesh -> {out_path.name}  "
          f"(courant C={COURANT_TARGET} dt={COURANT_DT_S} s)")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    wall = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh exit {rc}")
    print(f"[39] {out_path.name} wall: {wall:.1f} s")
    return wall


def _post_hoc_courant_metrics(mesh) -> dict[str, float]:
    """Compute, per element: minimum edge length in metres, char_vel
    from depth, and the implied minimum CFL-feasible dt at C=0.7.
    Returns the 1 %-ile worst-case dt and depth-stratified ratios.
    """
    if mesh.n_elements == 0:
        return {"min_cfl_dt_s_p01": 0.0, "min_cfl_dt_s_p50": 0.0,
                "n_below_target_dt": 0}
    p0 = mesh.nodes[mesh.elements[:, 0]]
    p1 = mesh.nodes[mesh.elements[:, 1]]
    p2 = mesh.nodes[mesh.elements[:, 2]]
    centroids = (p0 + p1 + p2) / 3.0
    centroid_depths = mesh.depths[mesh.elements].mean(axis=1)

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

    # Apply the same characteristic-velocity formula the sizing fn uses.
    abs_h = np.maximum(np.abs(centroid_depths), 1.0)
    nu = COURANT_NU_M
    deep = abs_h > nu
    sqrt_gh = np.sqrt(9.81 * abs_h)
    u_mag = np.where(deep, nu * np.sqrt(9.81 / abs_h),
                     np.sqrt(9.81 * nu))
    char_vel = np.where(deep, u_mag + sqrt_gh, 2.0 * u_mag)

    cfl_dt = COURANT_TARGET * dx_min_m / char_vel
    return {
        "min_cfl_dt_s_p01": float(np.percentile(cfl_dt, 1)),
        "min_cfl_dt_s_p50": float(np.percentile(cfl_dt, 50)),
        "n_below_target_dt": int((cfl_dt < COURANT_DT_S).sum()),
    }


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not (DEM.exists() and COASTLINE.exists() and RIVERS.exists()):
        raise SystemExit("required inputs missing")
    if not F14_BASELINE.exists():
        raise SystemExit(
            f"baseline missing: {F14_BASELINE} — run PoC #19 first"
        )

    print(f"[39] reusing baseline {F14_BASELINE.name}")
    wall_c = _run_courant_build(F14_COURANT)

    mb = read_fort14(F14_BASELINE)
    mc = read_fort14(F14_COURANT)
    qb = compute_metrics(mb)
    qc = compute_metrics(mc)
    cb = _post_hoc_courant_metrics(mb)
    cc = _post_hoc_courant_metrics(mc)

    keys_q = [
        "n_nodes", "n_elements", "alpha_mean", "alpha_p05",
        "min_angle_p05_deg", "frac_lt_20deg",
        "max_valence", "n_overconnected", "n_flipped",
    ]
    keys_c = ["min_cfl_dt_s_p01", "min_cfl_dt_s_p50", "n_below_target_dt"]

    payload = {
        "config": {
            "hmin_m": HMIN_M, "hmax_m": HMAX_M,
            "courant_target": COURANT_TARGET,
            "courant_dt_s": COURANT_DT_S,
            "courant_wave_amplitude_m": COURANT_NU_M,
        },
        "baseline": {"path": str(F14_BASELINE.resolve()),
                     "metrics": qb, "courant": cb},
        "courant": {"path": str(F14_COURANT.resolve()),
                    "wall_seconds": wall_c,
                    "metrics": qc, "courant": cc},
    }
    SUMMARY_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = [
        "PoC #39 --om-courant-sizing A/B on Tokyo Bay",
        f"hmin={HMIN_M:g} m  hmax={HMAX_M:g} m  "
        f"C_target={COURANT_TARGET}  dt={COURANT_DT_S} s  "
        f"nu={COURANT_NU_M} m",
        "",
        f"  {'metric':<22}  {'baseline':>16}  {'courant':>16}  {'delta':>10}",
        "  " + "-" * 70,
    ]
    for k in keys_q:
        delta = qc[k] - qb[k] if isinstance(qb[k], (int, float)) else "-"
        lines.append(
            f"  {k:<22}  {qb[k]!s:>16}  {qc[k]!s:>16}  {delta!s:>10}"
        )
    for k in keys_c:
        delta = cc[k] - cb[k]
        lines.append(
            f"  {k:<22}  {cb[k]!s:>16}  {cc[k]!s:>16}  {delta!s:>10}"
        )
    lines.append(
        f"  {'wall_seconds_courant':<22}  "
        f"{'-':>16}  {wall_c:>16.1f}  {'-':>10}"
    )

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print(f"\nwrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
