"""PoC #47b: Pass D parameter sweep on v3 residual.

PoC #47a (default v1 settings) reported 0 accepts because the
default gate is strict and the boundary rejection is conservative.
This sweep quantifies what relaxations would yield, informing v2.

Sweeps two axes:

  alpha_target           ∈ {0.95, 0.85, 0.75}
  reject_boundary        ∈ {True, False}

Each combination produces an accept/reject count and an
addressable-fail share. Read-only and fast (~1 s per config) — no
mesh written.

Outputs:
    outputs/47b_phase_h_pass_d_sweep_summary.txt
    outputs/47b_phase_h_pass_d_sweep_summary.json
"""
from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path

from fvcom_mesh_tools.algorithms.quality import (
    alpha_quality,
    min_interior_angle,
)
from fvcom_mesh_tools.io import read_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    DEFAULT_PATCH_MAX_CLUSTER_SIZE,
    DEFAULT_PATCH_MIN_CLUSTER_SIZE,
    _attempt_patch_recdt,
    _find_fail_clusters,
)

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "43_phase_h_v3_optimized.14"
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "47b_phase_h_pass_d_sweep_summary.txt"
SUMMARY_JSON = OUT_DIR / "47b_phase_h_pass_d_sweep_summary.json"

MIN_ANGLE_TARGET = 20.0
MIN_CLUSTER_SIZE = DEFAULT_PATCH_MIN_CLUSTER_SIZE
MAX_CLUSTER_SIZE = DEFAULT_PATCH_MAX_CLUSTER_SIZE

ALPHA_SWEEP = (0.95, 0.85, 0.75)
BOUNDARY_SWEEP = (True, False)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mesh = read_fort14(INPUT)
    alpha = alpha_quality(mesh)
    min_ang = min_interior_angle(mesh)
    fail = (alpha < 0.95) | (min_ang < MIN_ANGLE_TARGET)
    n_fail = int(fail.sum())
    print(
        f"[47b] input: NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  "
        f"fail={n_fail:,}",
        flush=True,
    )

    # Cluster list is fixed (depends only on the fail mask, not on
    # the per-element gate used inside Pass D's accept logic).
    clusters = _find_fail_clusters(
        mesh,
        alpha_target=0.95, min_angle_target=MIN_ANGLE_TARGET,
        min_cluster_size=MIN_CLUSTER_SIZE,
        max_cluster_size=MAX_CLUSTER_SIZE,
    )
    n_clusters = len(clusters)
    print(
        f"[47b] {n_clusters:,} clusters in "
        f"[{MIN_CLUSTER_SIZE}, {MAX_CLUSTER_SIZE}]",
        flush=True,
    )

    sweep_results: list[dict] = []
    for alpha_t in ALPHA_SWEEP:
        for reject_bnd in BOUNDARY_SWEEP:
            t0 = time.perf_counter()
            n_accepted = 0
            n_fail_addressable = 0
            reason_hist: Counter = Counter()
            for cluster in clusters:
                new_mesh, _info, reason = _attempt_patch_recdt(
                    mesh, cluster,
                    alpha_target=alpha_t,
                    min_angle_target=MIN_ANGLE_TARGET,
                    reject_boundary_clusters=reject_bnd,
                )
                if new_mesh is not None:
                    n_accepted += 1
                    n_fail_addressable += int(cluster.size)
                else:
                    reason_hist[reason] += 1
            wall = time.perf_counter() - t0
            row = {
                "alpha_target": alpha_t,
                "reject_boundary_clusters": reject_bnd,
                "n_accepted": n_accepted,
                "n_fail_addressable": n_fail_addressable,
                "share_of_clusters": (
                    n_accepted / max(n_clusters, 1)
                ),
                "share_of_residual_fails": (
                    n_fail_addressable / max(n_fail, 1)
                ),
                "reject_reasons": dict(reason_hist),
                "wall_seconds": float(wall),
            }
            sweep_results.append(row)
            print(
                f"  α≥{alpha_t}  reject_bnd={reject_bnd}  "
                f"accepted={n_accepted}/{n_clusters} "
                f"({row['share_of_clusters']:.1%})  "
                f"addressable_fails={n_fail_addressable}/{n_fail} "
                f"({row['share_of_residual_fails']:.1%})  "
                f"wall={wall:.2f}s",
                flush=True,
            )

    SUMMARY_JSON.write_text(json.dumps({
        "input": str(INPUT),
        "n_nodes": int(mesh.n_nodes),
        "n_elements": int(mesh.n_elements),
        "n_fail": n_fail,
        "n_clusters_in_bounds": n_clusters,
        "alpha_sweep": list(ALPHA_SWEEP),
        "boundary_sweep": list(BOUNDARY_SWEEP),
        "min_angle_target": MIN_ANGLE_TARGET,
        "min_cluster_size": MIN_CLUSTER_SIZE,
        "max_cluster_size": MAX_CLUSTER_SIZE,
        "sweep_results": sweep_results,
    }, indent=2))

    lines = [
        f"PoC #47b — Pass D parameter sweep on {INPUT.name}",
        f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  fail={n_fail:,}",
        f"{n_clusters:,} clusters in "
        f"[{MIN_CLUSTER_SIZE}, {MAX_CLUSTER_SIZE}]",
        f"min_angle_target={MIN_ANGLE_TARGET}°",
        "",
        f"  {'α_target':>8}  {'reject_bnd':>10}  "
        f"{'accepted':>10}  {'cluster_%':>10}  "
        f"{'fail_%':>8}  {'top reject reason':<28}",
        "  " + "-" * 90,
    ]
    for row in sweep_results:
        top_reason = next(iter(
            sorted(
                row["reject_reasons"].items(),
                key=lambda kv: -kv[1],
            )
        ), ("none", 0))
        lines.append(
            f"  {row['alpha_target']:>8.2f}  "
            f"{str(row['reject_boundary_clusters']):>10}  "
            f"{row['n_accepted']:>10,}  "
            f"{row['share_of_clusters']:>10.1%}  "
            f"{row['share_of_residual_fails']:>8.1%}  "
            f"{top_reason[0]}={top_reason[1]:,}"
        )
    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
