"""PoC #26: ablate OCSMesh+gmsh post-processing to find the cause of
over-connected nodes in PoC #16.

Background: PoC #25 traced 440 over-connected nodes (max valence 26) in
``16_tokyo_bay_with_rivers.14`` to the OCSMesh+gmsh engine path; PoC
#19 with ``--engine oceanmesh`` on the same inputs produced only 3
over-connected nodes (max valence 9). This PoC ablates the two main
post-processing stages — longest-edge bisection refinement
(``--refine-min-angle``) and the swap+smooth quality pass
(``--quality-pass``) — to find which one is responsible.

Variants (engine fixed to ocsmesh, all other flags identical to PoC #16):

    1. baseline   : --quality-pass 6 --refine-min-angle 20
    2. no_refine  : --quality-pass 6 --refine-min-angle 0
    3. no_qpass   : --quality-pass 0 --refine-min-angle 20
    4. both_off   : --quality-pass 0 --refine-min-angle 0

For each variant we report:

    * NP, NE, mesh-build wall time
    * over-connected node count + max valence (default cap 8)
    * disjoint, dead-end, thin, thin-chain, unreachable counts (so we
      do not regress other detectors when shedding post-processing)

The variant with the lowest over-connected count and an otherwise
acceptable mesh identifies the upstream cause and gives a concrete
recommendation for OCSMesh-engine users until repair lands.

Outputs:
    outputs/26_tokyo_<variant>.14
    outputs/26_overconn_ablation_summary.txt
"""

from __future__ import annotations

import time
from pathlib import Path

from fvcom_mesh_tools.cli.buildmesh import main as buildmesh_main
from fvcom_mesh_tools.diagnostics import run_diagnostics
from fvcom_mesh_tools.io import read_fort14

REPO_ROOT = Path(__file__).resolve().parent.parent
DEM = REPO_ROOT / "data" / "bathymetry" / "tokyo_bay" / "dem_00_01_change.nc"
COASTLINE = (
    REPO_ROOT / "data" / "coastline" / "tokyo_bay"
    / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
)
RIVERS = REPO_ROOT / "data" / "rivers" / "tokyo_bay" / "tokyo_bay_rivers.csv"

OUT_DIR = REPO_ROOT / "outputs"
SUMMARY_TXT = OUT_DIR / "26_overconn_ablation_summary.txt"

VARIANTS: dict[str, dict[str, str]] = {
    "baseline":  {"quality_pass": "6", "refine_min_angle": "20"},
    "no_refine": {"quality_pass": "6", "refine_min_angle": "0"},
    "no_qpass":  {"quality_pass": "0", "refine_min_angle": "20"},
    "both_off":  {"quality_pass": "0", "refine_min_angle": "0"},
}


def _common_args(out_f14: Path) -> list[str]:
    """Engine + sizing + coastline + river flags identical to PoC #16."""
    return [
        str(DEM), str(out_f14),
        "--engine", "ocsmesh",
        "--hmin", "200", "--hmax", "5000", "--zmax", "0.0",
        "--interp-method", "linear",
        "--coastline", str(COASTLINE),
        "--coast-target-size", "200",
        "--coast-expansion-rate", "0.005",
        "--min-polygon-area-m2", "1000000",
        "--min-island-area-m2", "100000",
        "--open-merge-coast-gap", "50",
        "--river-inflow-points", str(RIVERS),
        "--river-segment-nodes", "5",
        "--river-ibtype", "21",
        "--river-snap-tol-m", "5000",
        "--land-ibtype", "20",
        "--perpfix-iters", "1",
        "--quiet",
    ]


def run_variant(name: str, opts: dict[str, str]) -> dict[str, object]:
    out_f14 = OUT_DIR / f"26_tokyo_{name}.14"
    args = _common_args(out_f14) + [
        "--quality-pass", opts["quality_pass"],
        "--refine-min-angle", opts["refine_min_angle"],
    ]
    print(f"\n[26] === variant={name} qp={opts['quality_pass']} "
          f"refine={opts['refine_min_angle']} ===")
    t0 = time.perf_counter()
    rc = buildmesh_main(args)
    wall = time.perf_counter() - t0
    if rc != 0:
        raise SystemExit(f"buildmesh ({name}) exited {rc}")

    mesh = read_fort14(out_f14)
    report = run_diagnostics(mesh, name=name, path=out_f14)

    return {
        "name": name,
        "wall_s": wall,
        "n_nodes": mesh.n_nodes,
        "n_elements": mesh.n_elements,
        "max_valence": int(report.valence.max()),
        "n_overconn_default8": int(report.overconnected_flag.sum()),
        "n_disjoint": int(report.disjoint_flag.sum()),
        "n_dead_end": int(report.dead_end_flag.sum()),
        "n_thin": int(report.thin_flag.sum()),
        "n_thin_chain3": int(report.thin_chain_flag.sum()),
        "n_unreachable": int(report.unreachable_flag.sum()),
    }


def main() -> None:
    for p in (DEM, COASTLINE, RIVERS):
        if not p.exists():
            raise SystemExit(f"required input missing: {p}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for name, opts in VARIANTS.items():
        rows.append(run_variant(name, opts))

    header = (
        f"{'variant':<10s}{'wall(s)':>8s}{'NP':>10s}{'NE':>10s}"
        f"{'max_v':>8s}{'oc>8':>7s}{'disjt':>7s}{'deadE':>7s}"
        f"{'thin':>7s}{'thch3':>7s}{'unrch':>7s}"
    )
    sep = "-" * len(header)

    table_lines = [header, sep]
    for r in rows:
        table_lines.append(
            f"{r['name']:<10s}{r['wall_s']:>8.1f}{r['n_nodes']:>10,d}"
            f"{r['n_elements']:>10,d}{r['max_valence']:>8d}"
            f"{r['n_overconn_default8']:>7d}{r['n_disjoint']:>7d}"
            f"{r['n_dead_end']:>7d}{r['n_thin']:>7d}"
            f"{r['n_thin_chain3']:>7d}{r['n_unreachable']:>7d}"
        )

    text = "\n".join([
        "PoC #26: OCSMesh+gmsh post-processing ablation on Tokyo Bay (rivers).",
        f"DEM: {DEM}",
        f"coastline: {COASTLINE}",
        f"rivers: {RIVERS}",
        "engine: ocsmesh (fixed); other flags identical to PoC #16.",
        "All counts use the default fmesh-mesh-check thresholds "
        "(max_nbr_elem=8, min_thin_chain=3).",
        "",
        *table_lines,
        "",
        "Reference (PoC #16, ocsmesh, qp=6, refine=20): max_v=26, oc>8=440",
        "Reference (PoC #19, oceanmesh same inputs):    max_v=9,  oc>8=3",
    ])
    print()
    print(text)
    SUMMARY_TXT.write_text(text + "\n", encoding="utf-8")
    print(f"\n[26] wrote {SUMMARY_TXT}")


if __name__ == "__main__":
    main()
