"""PoC #59k-5 — crack the v3 residual floor (C1=37/C2=2/C4=21).

#59k-4 cleared every junction/manifold issue but ``phase_h_finish``
stalled: three seeds accepted zero moves. The residual sits on
100 m coastline detail where BOUNDARY nodes are frozen — exactly the
failure mode the two unused escalation levers address:

* ``phase_h_optimize`` Pass C (2-step lookahead) — cracks fail
  elements single operators cannot;
* ``phase_h_finish`` with a **coastline projector**: boundary nodes
  become movable ALONG the C23 coastline (tangent-clamped Gaussian +
  snap-back), so coastline fidelity is preserved by construction
  while the boundary-locked triangles gain a degree of freedom.

Projector-in-UTM trick: the C23 shapefile is reprojected to
EPSG:32654 and saved WITHOUT CRS metadata (so the loader keeps the
metric coordinates), and the snap tolerance is passed pre-scaled for
``mean_latitude_deg=0`` (deg-per-metre = 1/111194.9), making the
internal degree conversion an identity on metres.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from pyproj import Transformer

from fvcom_mesh_tools.algorithms import align_open_boundary_local
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.io.fvcom_native import export_fvcom_case
from fvcom_mesh_tools.mesh_clean import compact_nodes
from fvcom_mesh_tools.mesh_clean_phase_h import (
    build_coastline_projector,
    phase_h_finish,
    phase_h_optimize,
)
from fvcom_mesh_tools.qa import format_report, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59k_v3_gate_passed.14"   # stalled 59k-4 state
OUT = REPO / "outputs" / "59k_v3_gate_passed.14"
QA_JSON = REPO / "outputs" / "59k_v3_gate_passed_qa.json"
LOG_JSON = REPO / "outputs" / "59k5_cycle_log.json"
EXPORT_DIR = REPO / "outputs" / "fvcom_inputs_v3"
CASENAME = "tokyo_bay_v3"
C23 = REPO / "data" / "coastline" / "tokyo_bay" / "MLIT_C23" / "C23-06_TOKYOBAY.shp"
C23_UTM_DIR = REPO / "outputs" / "c23_utm54_nocrs"

MAX_CYCLES = 8
QUALITY_IDS = {"c1_min_angle", "c2_max_angle", "c4_area_change", "c5_valence"}
FATAL_OR_QUALITY = {
    "node_index_valid", "ccw_all_elements", "no_isolated_elements",
    "r4_mixed_boundary", "manifold_boundary", "no_duplicate_nodes",
    "no_orphan_nodes", "no_tiny_area", "isbce2_authentic",
    "obc_on_boundary", "obc_chain_adjacency", "obc_interior_neighbor",
    "obc_ordering", "c1_min_angle", "c2_max_angle", "c4_area_change",
    "c5_valence", "single_component", "obc_reachable", "min_depth_clip",
}

# Identity-on-metres projector scaling (see module docstring).
_DEG_PER_M_AT_EQ = 1.0 / 111_194.92664455873
SNAP_M = 500.0


def _utm_projector():
    import geopandas as gpd

    C23_UTM_DIR.mkdir(parents=True, exist_ok=True)
    shp = C23_UTM_DIR / "c23_utm54.shp"
    if not shp.exists():
        gdf = gpd.read_file(C23).set_crs(4326, allow_override=True).to_crs(32654)
        gdf = gdf.set_crs(None, allow_override=True)
        gdf.to_file(shp)
    return build_coastline_projector(
        [shp],
        max_snap_distance_m=SNAP_M / _DEG_PER_M_AT_EQ,
        mean_latitude_deg=0.0,
    )


def _failed(report):
    return {c.check_id for c in report.checks
            if c.gate and not c.skipped and not c.passed}


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    print(f"[59k5] input: NP={mesh.n_nodes:,} NE={mesh.n_elements:,} "
          f"obc={[len(s) for s in mesh.open_boundaries]}", flush=True)
    projector = _utm_projector()
    print(f"[59k5] coastline projector ready (snap {SNAP_M:.0f} m)", flush=True)

    log = []
    report = run_qa(mesh, name=OUT.name, path=OUT)
    prev_quality_total = None
    for cycle in range(MAX_CYCLES):
        failed = _failed(report)
        counts = {c.check_id: int(c.n_violations) for c in report.checks
                  if c.gate and not c.skipped and c.n_violations}
        print(f"[59k5] cycle {cycle}: failed = {sorted(failed)} {counts}",
              flush=True)
        log.append({"cycle": cycle, "failed": sorted(failed), "counts": counts})
        if not failed:
            break
        if failed & QUALITY_IDS:
            qtot = sum(counts.get(k, 0) for k in QUALITY_IDS)
            t1 = time.perf_counter()
            mesh, finfo = phase_h_finish(
                mesh, seed=200 + cycle, coastline_projector=projector,
            )
            mesh, cinfo = compact_nodes(mesh)
            print(f"[59k5] cycle {cycle}: projector-finish "
                  f"{finfo.get('before')} -> {finfo.get('after')} "
                  f"({time.perf_counter() - t1:.0f} s, compacted "
                  f"{cinfo['n_orphans_removed']})", flush=True)
            report = run_qa(mesh, name=OUT.name, path=OUT)
            new_counts = {c.check_id: int(c.n_violations)
                          for c in report.checks
                          if c.gate and not c.skipped and c.n_violations}
            new_qtot = sum(new_counts.get(k, 0) for k in QUALITY_IDS)
            if new_qtot >= qtot and prev_quality_total == qtot:
                # Projector finish stalled too — one lookahead
                # optimize round, then loop.
                t1 = time.perf_counter()
                mesh, oinfo = phase_h_optimize(
                    mesh,
                    min_angle_target=30.0,
                    max_angle_target=130.0,
                    lookahead_enabled=True,
                    pass_f_enabled=True,
                    pass_g_enabled=True,
                    pass_g_min_angle_target=30.0,
                    max_outer_rounds=3,
                    coastline_projector=projector,
                )
                mesh, cinfo = compact_nodes(mesh)
                print(f"[59k5] cycle {cycle}: lookahead optimize "
                      f"({time.perf_counter() - t1:.0f} s, compacted "
                      f"{cinfo['n_orphans_removed']})", flush=True)
                report = run_qa(mesh, name=OUT.name, path=OUT)
            prev_quality_total = qtot
        elif "obc_perpendicularity" in failed:
            mesh, pinfo = align_open_boundary_local(mesh, seed=9200 + cycle)
            print(f"[59k5] cycle {cycle}: perp fixer accepted="
                  f"{pinfo['accepted_total']} remaining="
                  f"{len(pinfo['remaining'])}", flush=True)
            report = run_qa(mesh, name=OUT.name, path=OUT)
            if pinfo["accepted_total"] == 0:
                break
        else:
            print(f"[59k5] cycle {cycle}: non-quality residual "
                  f"{sorted(failed)} — stopping", flush=True)
            break

    mesh.title = "PoC 59k UTM54N Tokyo Bay v3 (100 m coastline-fidelity)"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(json.dumps(report.to_dict(), indent=2, ensure_ascii=False),
                       encoding="utf-8")
    LOG_JSON.write_text(json.dumps(log, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    print(format_report(report, lang="ja"), flush=True)

    failed = _failed(report)
    if not (failed & FATAL_OR_QUALITY):
        tr = Transformer.from_crs("EPSG:32654", "EPSG:4326", always_xy=True)
        _lon, lat = tr.transform(mesh.nodes[:, 0], mesh.nodes[:, 1])
        obc = mesh.open_boundaries[0]
        sponge = [(int(v), 3000.0, 0.001) for v in obc]
        written = export_fvcom_case(
            mesh, EXPORT_DIR, CASENAME, obc_type=1, cor=lat, sponge=sponge,
        )
        for k, p in written.items():
            print(f"[59k5] export {k}: {p}", flush=True)
    else:
        print(f"[59k5] EXPORT SKIPPED: {sorted(failed & FATAL_OR_QUALITY)}",
              flush=True)
    print(f"[59k5] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
