# Phase H Finishing Chain — design notes

**Status:** PoC-level, validated on Tokyo-Bay (PoC #54c → #58k,
2026-05-14). 77 → 1 violation on 87 k elements with C1 / C2 / C5
= 0 and 1 residual C4 = 0.536. Promotion to a production module
is item #1 in [project open candidates](../docs/) /
`MEMORY.md`.

## 1. Motivation

After Phase H deterministic passes A → F (smooth + topology +
lookahead + cluster reCDT + C4-aware smoothing + C1-aware
smoothing) the Tokyo-Bay reference mesh stalled at **50 violations
/ 87 k elements** (PoC #57 Stage 1 — 2 C1 + 48 C4). Diagnostics on
the residual (PoC #58a / #58c) confirmed the deterministic
operators had **converged to a fixed point**: each Pass A-G uses
a fixed proposal mechanism (Laplacian centroid or a specific
operator inventory) under a strict count-gate, and the remaining
fails sat in geometries where every proposal landed on a
neighbouring fail.

## 2. The chain

```
                                      mesh
                                       │
       PoC #54c   build-time:          │
                  oceanmesh g = 0.10,  ▼
                  thin_chain = none  ┌────────────────────┐
                                     │  Phase G (build)   │
                                     └─────────┬──────────┘
                                               │
       PoC #56b   ──────────────────►          │
                  keep_components              │
                  (min_elements = 5)           │
                                               │
                                               ▼
                                     ┌────────────────────┐
                                     │   Phase H A..F     │  ← PoC #57 Stage 1
                                     │  (deterministic)   │     50 violations
                                     └─────────┬──────────┘
                                               │
       PoC #58d   ──────────────────►          │
                  stochastic local fixer       │
                  (seed = 42, move-only)       │
                                               │
                                               ▼
                                     ┌────────────────────┐
                                     │  4 violations      │  ← seed = 42
                                     │  4.4 ± 0.49 (5     │     mean across
                                     │   seeds, PoC #58h) │     seeds 42..46
                                     └─────────┬──────────┘
                                               │
       PoC #58j   ──────────────────►          │
                  vertex_remove(force) on      │
                  C1 interior obtuse vertex    │
                  + stochastic cleanup         │
                                               │
                                               ▼
                                     ┌────────────────────┐
                                     │  2 violations      │  ← C1 = 0
                                     │  (C1 = 0, C4 = 2)  │
                                     └─────────┬──────────┘
                                               │
       PoC #58k   ──────────────────►          │
                  vertex_remove(force) on      │
                  C4 fail-pair interior        │
                  vertex + stochastic cleanup  │
                                               │
                                               ▼
                                     ┌────────────────────┐
                                     │  1 violation       │  ← coastline floor
                                     │  (C4 = 1,          │
                                     │   ratio = 0.536)   │
                                     └────────────────────┘
```

## 3. The stochastic local fixer (PoC #58d / #58e)

The single most important component. For each current fail element:

```python
rng = numpy.random.default_rng(seed=42)
for try_i in range(MAX_TRIES_PER_FAIL):  # 500
    op = "move" if rng.random() < 0.7 else "insert"
    if op == "move":
        v = rng.choice(fail_element.vertices)
        delta = rng.standard_normal(2) * 0.3 * mean_local_edge_length
        if boundary_node_mask[v]:
            delta = project_onto_boundary_tangent(v, delta)
        new_pos = mesh.nodes[v] + delta
        # check signed area + C1/C2/C4/C5 on n2e[v] only
        if local_quality_ok and target_eid_no_longer_fails:
            accept
    elif op == "insert":
        bary = sample_simplex(min=0.05)
        # split into 3 sub-triangles
        # check local quality (currently rejects 100% of the time
        # because sub-tri angles geometrically near 30°)
```

**Critical design constants:**

* Local patch = **the 1-ring of the moved vertex** (`n2e[v]`), not
  "fail element + edge buddies". This is the smallest correct
  patch (exactly the elements whose geometry changes when v moves).
  An earlier draft used the looser patch and silently introduced
  C1/C2 regressions (50 → 9 vs corrected 50 → 5).
* σ = 0.3 × mean local edge length. Smaller → can't escape; larger
  → triangles flip.
* MAX_TRIES_PER_FAIL = 500. Many fixable fails clear at 50-100
  tries, but the long tail genuinely needs the budget.
* Seed-fixed `numpy.random.default_rng(seed=42)` makes every run
  bit-identical. PoC #58h verified the result is robust across
  seeds 42-46 (mean 4.4 / median 4.0 / std 0.49 violations).

**Why move dominates insert:** PoC #58e/#58f wired the insert
operator at 30 % dispatch weight; PoC #58f loosened the gate to
`BARYCENTRIC_MIN = 0.05` + count-comparison perimeter C4 check.
Still **0 inserts accepted in ~14 k attempts**. Root cause is
the sub-triangle min-angle gate: random barycentric insertion
into a fail element forces at least one of the 3 sub-triangles
to have a vertex angle near 30°. For C1 fails specifically, the
thin angle is preserved in every sub-triangle (geometric
impossibility). A useful insert operator would need a smarter
proposal (Ruppert circumcenter or edge-swap-equivalent).

## 4. Chaining with `vertex_remove`

After the stochastic fixer plateaus, the residual is structurally
beyond local random perturbation. The next layer composes
`_apply_vertex_remove(force=True)` with another stochastic
cleanup:

```python
for fail_eid in remaining_fails:
    for vertex in interior_vertices_of(fail_eid):
        snapshot = clone(mesh)
        new_mesh, info = _apply_vertex_remove(
            snapshot, vertex, n2e[vertex],
            force=True,
        )
        if new_mesh is None:
            continue
        # Delaunay retriangulation often introduces 1-2 new
        # marginal fails near the removed rim.  The stochastic
        # fixer absorbs them.
        new_mesh, _ = run_stochastic_fixer(new_mesh, seed=42)
        if total_violations(new_mesh) < total_violations(mesh):
            mesh = new_mesh
            break
```

**The compositional principle:** topology changes
(`vertex_remove`) free geometry that the stochastic fixer cannot
escape on its own, and the stochastic cleanup absorbs the local
Delaunay artefacts that the topology change introduces. Neither
half clears the residual alone — `vertex_remove` with a strict
count-gate rejects the +2 regression from the bare Delaunay rim,
and the stochastic fixer plateaus at ~4. The chain mirrors what
SMS users do manually ("delete the bad node, let the smoother fix
the new gap"), expressed as deterministic code.

## 5. Structural floor

After 1 vertex_remove on C1 (PoC #58j) and 1 on C4 (PoC #58k), the
Tokyo-Bay residual is **1 C4 fail at ratio 0.536**.  The element
pair shares an internal edge but **all 4 vertices are on the
coastline** and the quadrilateral they form is **concave** — so
`_apply_edge_swap` is rejected on signed area (the alternate
diagonal gives CW triangles, PoC #58l).  Closing this last fail
would require **editing the coastline** (collapsing one of the 4
boundary nodes onto a neighbour), which changes the geometry
input and is out of scope for the automatic mesh repair pipeline.

## 6. Open questions for promotion

* **Single seed vs. seed sweep at production time.** PoC #58h
  showed seeds 42-46 produce residuals in [4, 5] on the same
  input. Should `phase_h_finish` accept a single seed, or run
  multiple seeds and pick the lowest residual?
* **Vertex_remove ordering policy.** PoC #58k iterates C4 fail
  pairs and tries the highest-valence interior vertex first.
  Greedy on the first improving candidate. Should it instead
  exhaustively try every candidate and pick the global minimum?
  (cost ≈ 8 × current).
* **Insert operator design.** The 0-accept rate of random
  barycentric is a known geometric fact. Adding a Ruppert
  circumcenter proposal might unlock the all-boundary concave
  case (insert near the concave corner instead of inside the
  element).
* **Coastline editing.** Out of scope today, but the
  structural-floor concept suggests an optional Stage 6
  "boundary-node collapse on concave pairs" pass for users who
  prefer mesh quality over coastline fidelity.
