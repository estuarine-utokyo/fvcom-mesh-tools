"""PoC #47a-prep: cluster-structure analysis of the v3 residual.

The Pass D design (``docs/patch_re_cdt_design.md``) assumes that the
v3 residual is dominated by **cluster-scale** defects — connected
components of fail elements under face-face adjacency that a single
1-ring local edit cannot lift. This script *measures* that
assumption read-only before any Pass D code is written.

For each connected component ``C`` of fail elements in
``outputs/43_phase_h_v3_optimized.14``:

  size = |C|             — element count
  share = size / total_fail
  cum_share = sum of shares for clusters of size ≤ N

We report the histogram (size 1, 2, 3, 4-9, 10-29, 30+), the
cumulative coverage of total fail count by each bucket, and the
distribution of penalties within clusters.

If most fails live in clusters of size 1 (lone) — Pass D buys
little; the 1-ring lookahead path (v4.1) is the right tool.
If clusters of size ≥ 3 carry the majority — Pass D is plausibly
the right next operator class.

Outputs (read-only):
    outputs/47a_cluster_structure_summary.txt
    outputs/47a_cluster_structure_summary.json
"""
from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components

from fvcom_mesh_tools.algorithms.quality import (
    alpha_quality,
    min_interior_angle,
)
from fvcom_mesh_tools.diagnostics import face_face_adjacency
from fvcom_mesh_tools.io import read_fort14

REPO = Path(__file__).resolve().parent.parent
INPUT = REPO / "outputs" / "43_phase_h_v3_optimized.14"
OUT_DIR = REPO / "outputs"
SUMMARY_TXT = OUT_DIR / "47a_cluster_structure_summary.txt"
SUMMARY_JSON = OUT_DIR / "47a_cluster_structure_summary.json"

ALPHA_TARGET = 0.95
MIN_ANGLE_TARGET = 20.0

# Bucket boundaries: [1], [2], [3], [4-9], [10-29], [30+]
BUCKETS: tuple[tuple[str, int, int | None], ...] = (
    ("size = 1 (lone fails)", 1, 1),
    ("size = 2 (adjacent pair)", 2, 2),
    ("size = 3", 3, 3),
    ("size 4-9", 4, 9),
    ("size 10-29", 10, 29),
    ("size >= 30", 30, None),
)


def _bucket_label(size: int) -> str:
    for label, lo, hi in BUCKETS:
        if hi is None and size >= lo:
            return label
        if hi is not None and lo <= size <= hi:
            return label
    return "size = 0 (unreachable)"


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists():
        raise SystemExit(f"input missing: {INPUT}")

    t0 = time.perf_counter()
    mesh = read_fort14(INPUT)
    alpha = alpha_quality(mesh)
    min_ang = min_interior_angle(mesh)
    fail = (alpha < ALPHA_TARGET) | (min_ang < MIN_ANGLE_TARGET)
    fail_eids = np.where(fail)[0]
    n_fail = int(fail_eids.size)
    print(f"[47a] input: NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}",
          flush=True)
    print(f"[47a] fail elements: {n_fail:,} "
          f"({n_fail / max(mesh.n_elements, 1):.4%})",
          flush=True)
    if n_fail == 0:
        raise SystemExit("[47a] no fail elements — nothing to analyse")

    print("[47a] building face-face adjacency ...", flush=True)
    adj = face_face_adjacency(mesh.elements)

    # Restrict adjacency to the fail subgraph by zeroing rows / cols
    # of non-fail elements. We do this via a node-set mask on the row
    # / column indices of the CSR matrix.
    fail_mask = fail.astype(bool)
    fail_sub = adj[fail_mask][:, fail_mask]  # (n_fail, n_fail) CSR
    print(
        f"[47a] fail-subgraph: nnz={fail_sub.nnz:,} "
        f"(symmetric → ~{fail_sub.nnz // 2:,} unique edges)",
        flush=True,
    )

    n_components, labels = connected_components(
        fail_sub, directed=False, return_labels=True,
    )
    # labels[i] is the cluster id of the i-th fail element, in order
    # of np.where(fail) (i.e. ascending element id). We never need
    # the global-id → cluster-id reverse map for the summary; only
    # cluster sizes and top-cluster stats.
    cluster_sizes = np.bincount(labels, minlength=n_components)
    print(
        f"[47a] {n_components:,} clusters; "
        f"size stats min={cluster_sizes.min()}  "
        f"max={cluster_sizes.max()}  "
        f"mean={cluster_sizes.mean():.2f}  "
        f"median={int(np.median(cluster_sizes))}",
        flush=True,
    )

    # Histogram by bucket, plus coverage of total fail count.
    bucket_n_clusters: Counter = Counter()
    bucket_n_fails: Counter = Counter()
    for size in cluster_sizes:
        label = _bucket_label(int(size))
        bucket_n_clusters[label] += 1
        bucket_n_fails[label] += int(size)

    # Cluster penalty / α stats per cluster for the top-N biggest.
    a_pen = np.maximum(0.0, ALPHA_TARGET - alpha) ** 2
    g_pen = np.maximum(0.0, MIN_ANGLE_TARGET - min_ang) ** 2 / 100.0
    pen = a_pen + g_pen
    top_clusters: list[dict] = []
    order = np.argsort(-cluster_sizes)[:20]
    for cid in order:
        members = fail_eids[labels == cid]
        if members.size == 0:
            continue
        top_clusters.append({
            "cluster_id": int(cid),
            "size": int(members.size),
            "alpha_min": float(alpha[members].min()),
            "alpha_mean": float(alpha[members].mean()),
            "min_angle_min_deg": float(min_ang[members].min()),
            "min_angle_mean_deg": float(min_ang[members].mean()),
            "penalty_sum": float(pen[members].sum()),
            "penalty_max": float(pen[members].max()),
        })

    elapsed = time.perf_counter() - t0
    print(f"[47a] analysis wall: {elapsed:.1f} s", flush=True)

    # Pretty TXT summary.
    lines = [
        f"PoC #47a-prep — cluster-structure analysis of {INPUT.name}",
        f"input: NP={mesh.n_nodes:,}  NE={mesh.n_elements:,}",
        (
            f"fail elements: {n_fail:,} "
            f"({n_fail / max(mesh.n_elements, 1):.4%})"
        ),
        f"clusters under face-face adjacency: {n_components:,}",
        (
            f"cluster size:  min={cluster_sizes.min()}  "
            f"max={cluster_sizes.max()}  "
            f"mean={cluster_sizes.mean():.2f}  "
            f"median={int(np.median(cluster_sizes))}"
        ),
        f"analysis wall: {elapsed:.2f} s",
        "",
        f"  {'bucket':<26}  {'#clusters':>10}  {'#fails':>10}  "
        f"{'share':>8}  {'cumul':>8}",
        "  " + "-" * 70,
    ]
    cumul = 0
    for label, _lo, _hi in BUCKETS:
        nc = bucket_n_clusters.get(label, 0)
        nf = bucket_n_fails.get(label, 0)
        cumul += nf
        share = nf / max(n_fail, 1)
        cumul_share = cumul / max(n_fail, 1)
        lines.append(
            f"  {label:<26}  {nc:>10,}  {nf:>10,}  "
            f"{share:>7.1%}  {cumul_share:>7.1%}"
        )
    lines.append("")
    lines.append("Top-20 largest clusters (size, α_min, min_angle_min, Σ penalty):")
    lines.append(
        f"  {'cluster':>8}  {'size':>5}  {'α_min':>7}  "
        f"{'ang_min':>7}  {'Σpen':>8}"
    )
    for c in top_clusters:
        lines.append(
            f"  {c['cluster_id']:>8}  {c['size']:>5}  "
            f"{c['alpha_min']:>7.4f}  {c['min_angle_min_deg']:>7.2f}  "
            f"{c['penalty_sum']:>8.4f}"
        )

    SUMMARY_TXT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    SUMMARY_JSON.write_text(json.dumps({
        "input": str(INPUT),
        "n_nodes": int(mesh.n_nodes),
        "n_elements": int(mesh.n_elements),
        "alpha_target": ALPHA_TARGET,
        "min_angle_target": MIN_ANGLE_TARGET,
        "n_fail": n_fail,
        "n_clusters": int(n_components),
        "cluster_size_min": int(cluster_sizes.min()),
        "cluster_size_max": int(cluster_sizes.max()),
        "cluster_size_mean": float(cluster_sizes.mean()),
        "cluster_size_median": float(np.median(cluster_sizes)),
        "bucket_n_clusters": dict(bucket_n_clusters),
        "bucket_n_fails": dict(bucket_n_fails),
        "top_clusters": top_clusters,
        "wall_seconds": float(elapsed),
    }, indent=2))
    print()
    print("\n".join(lines))
    print()
    print(f"wrote {SUMMARY_TXT}")
    print(f"wrote {SUMMARY_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
