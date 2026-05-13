"""PoC #47a: Pass D (cluster patch re-CDT) dry-run on v3 residual.

The Pass D design (``docs/patch_re_cdt_design.md``) and the
cluster-structure analysis (PoC #47a-prep) together project a
51.4 % addressable share of the v3 residual. This PoC measures the
actual operator-level yield read-only — for every fail cluster in
``outputs/43_phase_h_v3_optimized.14`` we run
``_apply_patch_recdt`` once and record whether it would accept.

We report:

* number of clusters tried, accepted, rejected (with size-bucket
  histograms);
* total fail count *addressable* by Pass D — i.e. the sum of
  ``cluster_size`` over accepted clusters;
* the cumulative fail-reduction projection.

This is read-only — no mesh is written. The bare ``_apply_patch_recdt``
helper is called directly with ``reject_boundary_clusters=True`` to
match the default v1 policy.

Outputs:
    outputs/47a_phase_h_pass_d_dry_run_summary.txt
    outputs/47a_phase_h_pass_d_dry_run_summary.json
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
SUMMARY_TXT = OUT_DIR / "47a_phase_h_pass_d_dry_run_summary.txt"
SUMMARY_JSON = OUT_DIR / "47a_phase_h_pass_d_dry_run_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 20.0
MIN_CLUSTER_SIZE = DEFAULT_PATCH_MIN_CLUSTER_SIZE
MAX_CLUSTER_SIZE = DEFAULT_PATCH_MAX_CLUSTER_SIZE
REJECT_BOUNDARY = True


def _size_bucket(size: int) -> str:
    if size <= 9:
        return f"size_{size}"
    decade = (size // 10) * 10
    return f"size_{decade}_{decade + 9}"


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")

    t_load = time.perf_counter()
    mesh = read_fort14(INPUT)
    alpha = alpha_quality(mesh)
    min_ang = min_interior_angle(mesh)
    fail = (alpha < ALPHA_TARGET) | (min_ang < MIN_ANGLE_TARGET)
    n_fail = int(fail.sum())
    t_load = time.perf_counter() - t_load

    print(
        f"[47a] input: NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}",
        flush=True,
    )
    print(
        f"[47a] fail elements: {n_fail:,} "
        f"({n_fail / max(mesh.n_elements, 1):.4%})  "
        f"load wall: {t_load:.2f} s",
        flush=True,
    )

    t_clust = time.perf_counter()
    clusters = _find_fail_clusters(
        mesh,
        alpha_target=ALPHA_TARGET, min_angle_target=MIN_ANGLE_TARGET,
        min_cluster_size=MIN_CLUSTER_SIZE,
        max_cluster_size=MAX_CLUSTER_SIZE,
    )
    t_clust = time.perf_counter() - t_clust
    n_clusters = len(clusters)
    n_fail_in_clusters = sum(c.size for c in clusters)
    print(
        f"[47a] {n_clusters:,} clusters in "
        f"[{MIN_CLUSTER_SIZE}, {MAX_CLUSTER_SIZE}] "
        f"covering {n_fail_in_clusters:,} fails "
        f"({n_fail_in_clusters / max(n_fail, 1):.1%})  "
        f"cluster build: {t_clust:.2f} s",
        flush=True,
    )

    n_accepted = 0
    n_rejected = 0
    n_fail_addressed = 0
    accept_size_hist: Counter = Counter()
    reject_size_hist: Counter = Counter()
    reject_reason_hist: Counter = Counter()
    sample_accepts: list[dict] = []
    sample_rejects: list[dict] = []

    t_dry = time.perf_counter()
    for i, cluster in enumerate(clusters):
        size = int(cluster.size)
        new_mesh, info, reason = _attempt_patch_recdt(
            mesh, cluster,
            alpha_target=ALPHA_TARGET,
            min_angle_target=MIN_ANGLE_TARGET,
            reject_boundary_clusters=REJECT_BOUNDARY,
        )
        if new_mesh is not None and info is not None:
            n_accepted += 1
            n_fail_addressed += size
            accept_size_hist[_size_bucket(size)] += 1
            if len(sample_accepts) < 20:
                sample_accepts.append({
                    "cluster_index": i,
                    "cluster_size": size,
                    "rim_size": int(info["rim_size"]),
                    "n_new_elements": int(info["n_new_elements"]),
                    "n_interior_orphaned": int(info["n_interior_orphaned"]),
                    "external_rim_fail_before": int(
                        info["external_rim_fail_before"]
                    ),
                    "external_rim_fail_after": int(
                        info["external_rim_fail_after"]
                    ),
                })
        else:
            n_rejected += 1
            reject_size_hist[_size_bucket(size)] += 1
            reject_reason_hist[reason] += 1
            if len(sample_rejects) < 20:
                sample_rejects.append({
                    "cluster_index": i,
                    "cluster_size": size,
                    "reason": reason,
                    "cluster_eids_head": [int(e) for e in cluster[:6]],
                })
        if (i + 1) % 500 == 0:
            elapsed = time.perf_counter() - t_dry
            rate = (i + 1) / elapsed
            eta_min = (n_clusters - i - 1) / rate / 60.0
            print(
                f"  [{i+1:>5d}/{n_clusters}] "
                f"accepted={n_accepted}  rejected={n_rejected}  "
                f"addressable_fails={n_fail_addressed}  "
                f"rate={rate:.0f}/s  eta={eta_min:.1f}min",
                flush=True,
            )
    t_dry = time.perf_counter() - t_dry

    addressable_share = n_fail_addressed / max(n_fail, 1)
    cluster_accept_share = n_accepted / max(n_clusters, 1)

    print(
        f"[47a] done in {t_dry:.0f}s "
        f"({n_clusters / max(t_dry, 1e-9):.0f} clusters/s)",
        flush=True,
    )
    print(
        f"  accepted        : {n_accepted:>6,} / {n_clusters:>6,} "
        f"({cluster_accept_share:.1%} of clusters)",
        flush=True,
    )
    print(
        f"  addressable fail: {n_fail_addressed:>6,} / {n_fail:>6,} "
        f"({addressable_share:.1%} of all v3 residual fails)",
        flush=True,
    )

    pair_lines = []
    pair_lines.append(
        f"  {'bucket':<14}  {'accepted':>10}  {'rejected':>10}"
    )
    pair_lines.append("  " + "-" * 38)
    all_buckets = sorted(
        set(accept_size_hist.keys()) | set(reject_size_hist.keys()),
        key=lambda s: int(s.split("_")[1]),
    )
    for b in all_buckets:
        pair_lines.append(
            f"  {b:<14}  "
            f"{accept_size_hist.get(b, 0):>10,}  "
            f"{reject_size_hist.get(b, 0):>10,}"
        )

    reason_lines = [
        f"  {'reason':<28}  {'count':>8}  {'share':>7}",
        "  " + "-" * 50,
    ]
    for reason, count in reject_reason_hist.most_common():
        reason_lines.append(
            f"  {reason:<28}  {count:>8,}  "
            f"{count / max(n_rejected, 1):>7.1%}"
        )

    summary_lines = [
        f"PoC #47a — Pass D dry-run on {INPUT.name}",
        f"NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}  "
        f"fail={n_fail:,}",
        f"cluster bounds: [{MIN_CLUSTER_SIZE}, {MAX_CLUSTER_SIZE}]  "
        f"reject_boundary_clusters={REJECT_BOUNDARY}",
        f"alpha_target={ALPHA_TARGET}  min_angle_target={MIN_ANGLE_TARGET}°",
        "",
        f"clusters tried:   {n_clusters:,}",
        f"clusters accepted:{n_accepted:,} ({cluster_accept_share:.1%})",
        f"clusters rejected:{n_rejected:,} ({1 - cluster_accept_share:.1%})",
        f"fails covered by accepted clusters: "
        f"{n_fail_addressed:,} / {n_fail:,} ({addressable_share:.1%})",
        f"dry-run wall: {t_dry:.0f} s",
        "",
        "per-size-bucket breakdown:",
        *pair_lines,
        "",
        "reject reason histogram:",
        *reason_lines,
    ]
    SUMMARY_TXT.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    SUMMARY_JSON.write_text(json.dumps({
        "input": str(INPUT),
        "n_nodes": int(mesh.n_nodes),
        "n_elements": int(mesh.n_elements),
        "n_fail": n_fail,
        "alpha_target": ALPHA_TARGET,
        "min_angle_target": MIN_ANGLE_TARGET,
        "min_cluster_size": MIN_CLUSTER_SIZE,
        "max_cluster_size": MAX_CLUSTER_SIZE,
        "reject_boundary_clusters": REJECT_BOUNDARY,
        "n_clusters_in_bounds": n_clusters,
        "n_fail_in_clusters": n_fail_in_clusters,
        "n_accepted": n_accepted,
        "n_rejected": n_rejected,
        "n_fail_addressable": n_fail_addressed,
        "fail_addressable_share": float(addressable_share),
        "cluster_accept_share": float(cluster_accept_share),
        "accept_size_histogram": dict(accept_size_hist),
        "reject_size_histogram": dict(reject_size_hist),
        "reject_reason_histogram": dict(reject_reason_hist),
        "sample_accepts": sample_accepts,
        "sample_rejects": sample_rejects,
        "wall_seconds": {
            "load": float(t_load),
            "cluster_build": float(t_clust),
            "dry_run": float(t_dry),
        },
    }, indent=2))
    print()
    print("\n".join(summary_lines))
    print()
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
