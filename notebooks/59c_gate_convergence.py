"""PoC #59c — alternate perpfix / phase_h_finish until the QA gate passes.

PoC #59b reached C1 = C2 = C4 = C5 = 0 in metric space but left two
gates failing:

* ``no_orphan_nodes`` (4) — ``phase_h_finish``'s ``vertex_remove``
  keeps removed vertices in the node array; fixed here with the new
  ``compact_nodes`` utility.
* ``obc_perpendicularity`` (22 nodes, worst 41.0°) — the finish stage
  moves first-ring nodes. #59b's cycle logic exited early because the
  orphan failure masked the perp-only branch.

This PoC runs the proper convergence loop: QA -> (quality failing ->
finish + compact | perp failing -> damped perpfix) -> QA, bounded.

Output: outputs/59c_gate_passed.14 (+ _qa.json, cycle log).
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from fvcom_mesh_tools.algorithms import align_open_boundary_first_ring
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean import compact_nodes
from fvcom_mesh_tools.mesh_clean_phase_h import phase_h_finish
from fvcom_mesh_tools.qa import format_report, run_qa

REPO = Path(__file__).resolve().parents[1]
SRC = REPO / "outputs" / "59b_utm54n_finished.14"
OUT = REPO / "outputs" / "59c_gate_passed.14"
QA_JSON = REPO / "outputs" / "59c_gate_passed_qa.json"
CYCLES_JSON = REPO / "outputs" / "59c_cycle_log.json"

MAX_CYCLES = 6
QUALITY_IDS = {"c1_min_angle", "c2_max_angle", "c4_area_change", "c5_valence"}


def _failed_gate_ids(report) -> set[str]:
    return {
        c.check_id for c in report.checks
        if c.gate and not c.skipped and not c.passed
    }


def _perpfix_all(mesh, *, alpha: float, n_iters: int):
    for k in range(len(mesh.open_boundaries)):
        mesh, _ = align_open_boundary_first_ring(
            mesh, alpha=alpha, n_iters=n_iters,
            smooth_iters=2, smooth_alpha=0.3, segment_index=k,
        )
    return mesh


def main() -> int:
    t0 = time.perf_counter()
    mesh = read_fort14(SRC)
    mesh, cinfo = compact_nodes(mesh)
    print(f"[59c] compact_nodes: {cinfo}", flush=True)

    log: list[dict] = []
    report = run_qa(mesh, name=OUT.name, path=OUT)
    for cycle in range(MAX_CYCLES):
        failed = _failed_gate_ids(report)
        counts = {
            c.check_id: int(c.n_violations)
            for c in report.checks if c.gate and not c.skipped and c.n_violations
        }
        print(f"[59c] cycle {cycle}: failed = {sorted(failed)} {counts}", flush=True)
        log.append({"cycle": cycle, "failed_gates": sorted(failed), "counts": counts})
        if not failed:
            break
        if failed & QUALITY_IDS:
            mesh, finfo = phase_h_finish(mesh, seed=50 + cycle)
            mesh, cinfo = compact_nodes(mesh)
            print(
                f"[59c] cycle {cycle}: finish {finfo.get('before')} -> "
                f"{finfo.get('after')}; compacted {cinfo['n_orphans_removed']}",
                flush=True,
            )
        elif "obc_perpendicularity" in failed:
            mesh = _perpfix_all(mesh, alpha=0.6, n_iters=2)
            print(f"[59c] cycle {cycle}: damped perpfix applied", flush=True)
        else:
            print(f"[59c] cycle {cycle}: unhandled residual — stopping", flush=True)
            break
        report = run_qa(mesh, name=OUT.name, path=OUT)

    mesh.title = "PoC 59c UTM54N gate-converged Tokyo Bay mesh"
    write_fort14(mesh, OUT)
    QA_JSON.write_text(
        json.dumps(report.to_dict(), indent=2, ensure_ascii=False), encoding="utf-8",
    )
    CYCLES_JSON.write_text(
        json.dumps(log, indent=2, ensure_ascii=False), encoding="utf-8",
    )
    print(format_report(report, lang="ja"), flush=True)
    print(f"[59c] wrote {OUT}", flush=True)
    print(f"[59c] wall: {time.perf_counter() - t0:.1f} s", flush=True)
    return 0 if report.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
