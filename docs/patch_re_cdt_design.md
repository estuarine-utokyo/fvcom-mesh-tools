# Phase H Pass D — cluster-scale patch re-CDT (design sketch)

**Status:** v1 implemented; PoC #47a dry-run confirms 0 accepts at
``alpha_target=0.95`` with the default boundary-reject policy.
PoC #47b parameter sweep (below, § 1.3) attributes the failure to
two structural barriers — Delaunay quality and boundary handling —
and motivates a v2 with adaptive spine points + boundary segment
book-keeping.

## 1. Motivation

After Phase H v3 on the Tokyo-Bay reference mesh, **39 % of the
sampled v3 residual is unfixable by 2-step lookahead** (PoC #44,
n=1 000, seed=42). The PoC #44 caveat hypothesised that this
population is dominated by **cluster-scale defects**: multiple
adjacent fail elements whose joint quality cannot be lifted by any
sequence of single-vertex / single-edge edits without dragging a
non-fail neighbour into fail status. Even v4.1's strict
``target_exits_fail`` gate cannot crack them because op1 ∈
``{smooth_node}`` modifies only one vertex — insufficient to
re-position a cluster.

PoC #46 (the v4.1 strict-gate end-to-end run) is expected to
reduce the residual relative to v3 but not below the
1-ring-fundamental floor measured by PoC #44b. Whatever fraction of
the v3 residual remains is the candidate pool for **patch re-CDT**
— retriangulate a contiguous fail-element patch in one shot using
the rim polygon, accept only if every new patch element passes the
per-element gate.

This is what SMS users do by hand: identify a bad neighbourhood,
delete the local mesh inside a chosen boundary, re-mesh.

### 1.1 v4.1 theoretical ceiling (PoC #44b)

The companion read-only dry-run
``notebooks/44b_phase_h_target_exits_dry_run_poc.py`` measured the
v4.1 ``target_exits_fail`` gate on the same n=1,000 random sample
of the v3 residual that PoC #44 used. The result is striking:

  partition                                   count   share
  ------------------------------------------- ------- -------
  1-step strict fixable (sanity ✓)                 0    0.0 %
  op1-only (smooth_node alone exits E from fail)  20    2.0 %
  2-step (smooth_node + smooth_node)              26    2.6 %
  unfixable                                      954   95.4 %

The v4.1 driver's theoretical ceiling on this benchmark is
**4.6 %** of the v3 residual — extrapolated to 11,015 fails,
≈ 510 elements. PoC #45's 61 % "fixable" under the
``union_penalty`` gate dropped to 4.6 % under the strict
``target_exits_fail`` gate; the gap shows how many accepts the
v4 gate let through that did not actually fix the target.

**Implication for Pass D**: Pass C v4.1 closes at most ~510 of the
11,015 v3 residual fails. The remaining ~10,500 (95.4 %) are the
real market for Pass D — and the cluster-structure analysis below
shows 51.4 % of them sit in clusters of size ≥ 3 where no
sequence of 1-ring edits can help.

### 1.2 Empirical cluster structure (v3 residual)

The cluster-scale assumption was tested directly. Running the
read-only analysis ``notebooks/47a_cluster_structure_analysis.py``
on ``outputs/43_phase_h_v3_optimized.14`` (11,015 fail elements in
46,665 total) finds **5,298 connected components** of the fail
subgraph under face-face adjacency, with size 1 to 20:

    bucket                  #clusters   #fails   share   cumul
    -----------------------------------------------------------
    size = 1 (lone)             2,812    2,812   25.5 %  25.5 %
    size = 2 (adj. pair)        1,269    2,538   23.0 %  48.6 %
    size = 3                      499    1,497   13.6 %  62.2 %
    size 4-9                      655    3,405   30.9 %  93.1 %
    size 10-29                     63      763    6.9 %  100.0 %
    size >= 30                      0        0    0.0 %  100.0 %

The largest cluster has 20 elements (alpha_min 0.72, min_angle_min
27.5°). No mega-clusters threaten the design's ``max_cluster_size``
cap.

* **25.5 %** lone fails — Pass D buys nothing here; 1-ring
  lookahead (v4.1 Pass C) is the right tool.
* **23.0 %** adjacent pairs — Pass D *might* help (pair re-CDT on
  the 4-vertex rim is the smallest non-trivial case); needs
  validation in the dry-run PoC.
* **51.4 %** clusters of size ≥ 3 — Pass D's strong target.
  These cannot be fixed by any sequence of single-vertex / single-
  edge edits and are the addressable market for cluster re-CDT.

The findings confirm the design hypothesis: a non-trivial majority
of the residual sits in cluster components, the cluster size
distribution is well-bounded (max 20, mean 2.08, median 1) and a
``max_cluster_size`` of 100 leaves no cluster unhandled. Pass D's
upper bound is 74.5 % of residual fails (everything in size ≥ 2);
the realistic target is the 62.2 % cumulative coverage at size ≥ 3.

### 1.3 PoC #47a dry-run + #47b sweep (v1 result)

PoC #47a (``notebooks/47a_phase_h_pass_d_dry_run_poc.py``) tries
``_attempt_patch_recdt`` on every fail cluster in
``[3, 100]`` size range. With the v1 defaults
(``alpha_target=0.95``, ``reject_boundary_clusters=True``):

* **0 / 1 217 clusters accepted (0.0 %)**
* Reject reason histogram:
  * ``rim_on_boundary`` — 685 (56.3 %) — clusters whose rim sits
    on a Tokyo-Bay coastline segment
  * ``gate1_alpha`` — 523 (43.0 %) — pure-Delaunay re-mesh of the
    rim polygon produces at least one triangle below
    ``alpha_target``
  * ``rim_walk_failed`` — 9 (0.7 %) — non-simply-connected rim
  * ``gate2_rim_regression`` — 0 (the strict gate never fires;
    whenever Gate 1 passes, Gate 2 passes too)

PoC #47b (``notebooks/47b_phase_h_pass_d_sweep.py``) sweeps
``alpha_target ∈ {0.95, 0.85, 0.75}`` × ``reject_boundary ∈
{True, False}`` to quantify the two barriers:

    α_target  reject_bnd   accepted   cluster_%   fail_%   top reject
    -------   ----------   --------   ---------   ------   -----------
    0.95      True              0       0.0 %     0.0 %   rim_on_boundary=685
    0.95      False             0       0.0 %     0.0 %   gate1_alpha=1208
    0.85      True            468      38.5 %    18.9 %   rim_on_boundary=685
    0.85      False           919      75.5 %    35.5 %   gate1_alpha=289
    0.75      True            500      41.1 %    21.1 %   rim_on_boundary=685
    0.75      False         1 125      92.4 %    45.8 %   gate1_alpha=83

**Conclusion: Pass D v1 as designed is non-viable on Tokyo Bay.**
The pure-Delaunay rim re-mesh cannot lift cluster quality above
``alpha=0.95`` (43 % of clusters fail Gate 1), and the conservative
boundary-reject policy blocks another 56 % (rim on coast). The
intersection blocks the rest.

To unlock the design's projected 51 % addressable share, v2 needs
both:

* **Adaptive spine points** (Steiner / density-driven interior
  nodes), borrowed from PoC #37 Stage 2 medial-axis: lifts Gate 1
  pass rate by densifying the re-mesh away from the rim.
* **Boundary segment book-keeping**: walk the cluster's rim, track
  which rim edges coincide with an open / land segment, preserve
  the segment ordering through the re-mesh so the boundary lists
  stay valid. Lifts ``rim_on_boundary`` rejections.

Without these v2 enhancements, ``--phase-h-patch-recdt`` is left
off by default and documented as a known-low-yield baseline.

## 2. Hypothesis

Cluster fails persist after Pass A / B / C because the optimal
local placement of the cluster's interior vertices is
*non-unique* — any single move toward the optimum drags a neighbour
below threshold. A simultaneous re-triangulation of the cluster
(deleting all interior vertices and recomputing the connectivity)
can land in a global minimum the local solver cannot reach.

Stage 2 medial-axis re-meshing (Phase E, PoC #37) already
demonstrates the principle in a different context (under-resolved
channels). PoC #37 used a Constrained Delaunay Triangulation on a
patch ``(rim ∪ spine)`` and achieved 3-5× the local element count
while strictly improving channel resolution. Pass D borrows the
patch-extraction + re-CDT machinery but targets *quality* clusters
(low-α + low-min-angle) rather than width clusters (low w/h).

## 3. Operator design

### 3.1 Cluster detector

Input: a fort.14 mesh + the per-element fail mask.

```text
fail_eids = where(alpha < α_target ∨ min_angle < min_angle_target)
adj       = face_face_adjacency(mesh)
G_fail    = subgraph of adj induced by fail_eids
clusters  = connected_components(G_fail)
```

Output: a list of cluster element-id sets. Each cluster is a
connected component of the fail-subgraph under face-face adjacency.

Filter parameters (with sensible defaults from PoC #37 sweep):

| Parameter | Default | Purpose |
|---|---:|---|
| ``min_cluster_size`` | 3 | Skip lone fails (already handled by Pass C lookahead) |
| ``max_cluster_size`` | 100 | Skip huge clusters where re-CDT might re-introduce as many fails as it heals |
| ``alpha_target`` | 0.95 | Inherits the global gate |
| ``min_angle_target`` | 20.0 | Inherits the global gate |

PoC #37 found 10 as the sweet spot for *under-resolved* channel
patches; *quality* clusters are likely smaller. Start with 3
(per-element-greedy could not crack them) and tune via PoC #47.

### 3.2 Patch extraction

For each cluster's element set ``C``:

1. Collect every node referenced by any element in ``C``.
2. Classify each node as **interior** (every incident element is in
   ``C``) or **rim** (at least one incident element is outside
   ``C``).
3. Walk the rim by following non-``C`` ↔ ``C`` boundary edges → a
   single ordered rim polygon (assumed simply connected; reject and
   skip the cluster otherwise).
4. Snapshot the original block: rim node coordinates, interior
   node coordinates (will be re-positioned or deleted), depths,
   boundary-segment membership for any rim node that lies on an
   open / land boundary.

### 3.3 Re-triangulation

Reuse ``_retriangulate_patch(rim_xy, spine_xy, n_rim)`` from
``mesh_clean.py`` exactly as Stage 2 medial-axis does. Pass D's
spine is empty by default (let the patch be a pure rim Delaunay)
because we are not enforcing a w/h target. As a v2 enhancement,
allow the caller to supply **density-based spine points** sampled
from a sizing field for finer resolution; for v1 keep the spine
empty.

Output: a new triangle list indexed into ``[rim_xy; spine_xy]``,
already validated for CCW orientation and rim-edge preservation
(re-uses Stage 2's checks).

### 3.4 Acceptance gate (Pass D)

Each candidate patch produces a candidate mesh ``m_after`` (cluster
elements replaced by new triangles, interior nodes orphaned). The
gate is **stricter than v4.1**:

```text
accept iff
    all new patch elements pass per-element gate
        (α(e) >= α_target ∧ min_angle(e) >= min_angle_target
         for every e in new patch)
  AND
    no new fail element introduced *outside* the patch
        (fail count over rim 1-ring elements ≤ before)
```

The first half is the SMS standard ("the patch is repaired"). The
second half stops the patch from spraying new fails across the rim
boundary — which is exactly the 2-ring drift PoC #45 surfaced for
v4. Together they guarantee Pass D is monotonic on fail count.

If the patch fails the gate, reject and try the next cluster
unchanged. No partial accept (no "accept the better triangles and
keep the bad ones" — that re-introduces local complexity Pass D is
designed to escape).

### 3.5 Boundary handling

Rim nodes on an open / land segment must keep their boundary
membership. The re-CDT does not insert / move rim nodes (the rim
is fixed by extraction), so the only book-keeping is preserving
the segment arrays around the cluster. For v1 Pass D **rejects any
cluster whose rim crosses a boundary segment** to dodge the
book-keeping — boundary-overlapping clusters are typically already
addressed by Phase E's medial-axis path, so the residual we are
targeting is interior-dominated.

## 4. Driver integration

Pass D fits after Pass C in ``phase_h_optimize``:

```text
For round in range(max_outer_rounds):
    Pass A: batch Gauss-Seidel smooth   (existing)
    Pass B: per-element greedy 1-step  (existing)
    Pass C: 2-step lookahead           (v4.1, opt-in)
    Pass D: cluster patch re-CDT       (NEW, opt-in)
    If no pass accepted anything in this round: break
```

Configuration knobs on ``phase_h_optimize`` and the CLIs:

| Kwarg / flag | Default | Purpose |
|---|---|---|
| ``cluster_re_cdt_enabled`` / ``--phase-h-patch-recdt`` | ``False`` | Master switch |
| ``min_cluster_size`` / ``--phase-h-min-cluster-size`` | ``3`` | Filter out lone fails |
| ``max_cluster_size`` / ``--phase-h-max-cluster-size`` | ``100`` | Skip mega-clusters |
| ``max_patches_per_round`` / ``--phase-h-max-patches-per-round`` | ``1000`` | Per-round accept cap |

Pass D writes its own histogram in ``info``:

```python
info["patch_recdt_accepted"]: int           # number of patches retriangulated
info["patch_recdt_rejected_by_gate"]: int
info["patch_recdt_skipped_boundary"]: int
info["patch_recdt_skipped_size"]: int
info["patch_recdt_rim_size_hist"]: dict[int, int]  # for cost tracking
```

## 5. Reuse plan (concrete file map)

| New / changed | Path | Origin |
|---|---|---|
| `_face_face_adjacency` | already exists in `diagnostics.py` | reuse |
| `_retriangulate_patch` | already exists in `mesh_clean.py` | reuse |
| `_extract_cluster_rim_polygon` | new in `mesh_clean_phase_h.py` | analogous to PoC #37's ring extraction |
| `_apply_patch_recdt(mesh, cluster_eids, ...)` | new — wraps extract + retriangulate + apply | analogous to `_apply_vertex_remove` × N |
| `_patch_recdt_round(mesh, ...)` | new | analogous to `_topology_round` / `_lookahead_round` |
| `phase_h_optimize` | extended with Pass D | this PR |
| CLI flags | `meshclean.py`, `meshpipeline.py` | mirror existing Phase H flags |
| Tests | extend `test_mesh_clean_phase_h.py` | direct re-CDT calls + integration |

Net new LOC estimate: ~400 (driver + helpers + tests). Comparable to
the v4.1 Pass C addition.

## 6. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Re-CDT introduces new fails outside the rim ring | Gate (3.4) second half — reject patches that grow neighbouring fail count |
| Mega-cluster (e.g. an entire over-resolved peninsula) is hit and Delaunay struggles | ``max_cluster_size`` cap, plus `_retriangulate_patch` already returns ``None`` on degenerate / non-convex rims |
| Cluster rim crosses an open / land boundary | v1 rejects such clusters (5.5); v2 can lift this with explicit segment book-keeping mirroring `_apply_edge_split_boundary` |
| Boundary-tangent smoothness lost when adjacent boundary smooth-fails recouple | Pass A in the next outer round will re-smooth tangents; if instability persists, gate the Pass D accept by `n_flipped`-equality (no flips introduced) |
| Orphan interior nodes leak into the saved fort.14 | Already handled by `Fort14Mesh.write` — orphans are kept in the node array but referenced by zero elements (matches Phase H v1 `vertex_remove`'s pattern) |
| Pass D wall-time blow-up | `_retriangulate_patch` is O(rim²) at worst; cap rim size via `max_cluster_size`. Each patch is faster than a `vertex_remove` of the same neighbourhood (one Delaunay vs many) |

## 7. Validation plan

PoC #47 (after PoC #46 has fully reported):

* Input: ``outputs/46_phase_h_v4_1_optimized.14`` (v4.1's residual)
  if PoC #46 succeeds, else fall back to
  ``outputs/43_phase_h_v3_optimized.14``.
* Run ``phase_h_optimize(lookahead_enabled=True,
  cluster_re_cdt_enabled=True)`` end-to-end.
* Report: NP, NE delta, alpha p05, min_angle p05, frac<20°, fail
  count, wall, patch accept/reject breakdown.
* A/B/C/D table vs v3 / v4.1 / v4.

Pre-implementation safety check (analogous to PoC #44 / #44b):

* Dry-run PoC #47a: extract all fail clusters in the v3 / v4.1
  residual under face-face adjacency; for each, try the patch
  re-CDT and count "would Pass D accept?". Read-only.
* This bounds Pass D's potential before any code is wired into the
  driver. If the dry-run yield is <10 % of residual elements,
  consider the operator class a dead end and look elsewhere
  (simulated annealing on the cluster's interior nodes, sizing-
  function-aware re-mesh, etc.).

## 8. Out of scope (deferred to Pass E or later)

* Cluster splitting (merging across small no-fail elements that
  separate two large fail clusters) — handle as a v2 extension.
* Density-adaptive spine generation (PoC #37-style) — punt to v2.
* Multi-pass within a single Pass D round (re-cluster after each
  accept) — start with a static cluster list per round.
* Boundary-crossing clusters — punt to v2.
* Acceptance gates other than the strict per-element + no-regression
  combination — see PoC #45 for evidence that loose gates thrash.

## 9. Open questions

1. Should Pass D run *before* Pass C in the outer loop? — Pass C
   targets isolated 1-ring fails; Pass D handles the cluster
   residual. Running Pass D first leaves Pass C to clean up the
   smaller scraps Pass D leaves behind. Probably yes; will A/B in
   PoC #47.
2. Should we limit ``max_patches_per_round`` to keep wall-time
   bounded? — Yes; default 1000 mirrors Pass B/C caps. Each patch
   is ~10-50 ms (Delaunay on 10-30 points), so 1000 accepts ≈
   10-50 s of pure re-CDT compute per round, modest.
3. Is there benefit to allowing per-element gate **relaxation** for
   patch interiors (e.g. alpha >= 0.90 for new patch elements)? —
   Probably not; the design rationale is to escape into the
   global optimum, which should comfortably clear 0.95. If the
   gate is consistently unreachable, the cluster is genuinely
   under-resolved and belongs in Phase E, not Pass D.
