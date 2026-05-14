# Phase H finishing chain — user guide

`phase_h_finish` is the production entry point for the FVCOM
mesh-quality finishing chain that drives an already-meshed
`fort.14` to **the structural floor of the auto-pipeline**:
every FVCOM manual criterion (C1 min angle ≥ 30°, C2 max angle
≤ 130°, C4 adjacent-element area ratio ≤ 0.5, C5 valence ≤ 8)
satisfied except for a handful of irreducible coastline-geometry
artefacts.

On the Tokyo-Bay reference mesh, this chain drove **77 violations
to 1** (87 049 elements ≈ 0.00115 %) — see [Results](#results)
below.

This document is the operational reference; for the design
rationale and the empirical journey that established the recipe,
see [`phase_h_finishing_chain.md`](phase_h_finishing_chain.md).

## When to use it

Run `phase_h_finish` **after** the deterministic Phase H pipeline
(`phase_h_optimize` or `fmesh-mesh-pipeline --max-iters 4
--phase-h`) has converged. Its job is to clean up the residual
that the deterministic passes cannot escape:

* deterministic passes are fixed local-edit functions of the
  current geometry, so they hit fixed points;
* `phase_h_finish` breaks those fixed points with a seed-fixed
  stochastic local-search fixer that explores positions the
  deterministic proposers never propose, followed by targeted
  `vertex_remove` + stochastic cleanup over the surviving fail
  pairs.

If your residual after `phase_h_optimize` is already zero on the
FVCOM gates you care about, you do not need `phase_h_finish`.

## Quick start

```python
from fvcom_mesh_tools.io import read_fort14, write_fort14
from fvcom_mesh_tools.mesh_clean_phase_h import (
    build_coastline_projector,
    phase_h_finish,
)

mesh = read_fort14("after_phase_h_optimize.14")
projector = build_coastline_projector(
    ["data/coastline/tokyo_bay/MLIT_C23/C23-06_TOKYOBAY.shp"],
    max_snap_distance_m=500.0,
    mean_latitude_deg=float(mesh.nodes[:, 1].mean()),
)
out, info = phase_h_finish(
    mesh,
    seed=42,                         # bit-reproducible
    coastline_projector=projector,   # tangent + projection for boundary nodes
)
write_fort14(out, "after_phase_h_finish.14")
print(info["before"], "→", info["after"], f"(Δtotal {info['delta_total']:+d})")
```

Wall time is ~2-3 minutes for a 50-violation, 87 k-element mesh
on a single login-node core. Memory use is modest (peak ~500 MB
on Tokyo-Bay).

## Results

Cumulative improvement on the Tokyo-Bay reference mesh
(87 049 elements, FVCOM strict gates):

| Stage | Input | C1 | C2 | C4 | C5 | total |
|-------|-------|---:|---:|---:|---:|------:|
| PoC #54c (`g=0.10`, `thin_chain=none`) | build-time | 8 | 1 | 68 | 0 | **77** |
| PoC #56b (`keep_components(min=5)`) | post-process | 2 | 0 | 68 | 0 | **70** |
| PoC #57 Stage 1 (`phase_h_optimize` w/ Pass F) | A-F | 2 | 0 | 48 | 0 | **50** |
| PoC #58e (stochastic fixer, seed=42) | finish stage 1 | 1 | 0 | 3 | 0 | **4** |
| PoC #58j (vertex_remove chain on C1) | finish stage 2 | 0 | 0 | 2 | 0 | **2** |
| PoC #58k (vertex_remove chain on C4) | finish stage 3 | 0 | 0 | 1 | 0 | **1** |

PoC #58l confirmed the final fail (one C4 edge at ratio 0.536) is
**geometrically irreducible** without coastline editing: every
vertex of the failing element pair lies on the coastline, and the
quadrilateral they form is concave, so the alternate diagonal
gives CW (negative-area) triangles — `edge_swap` is rejected on
validity.

Running `phase_h_finish` end-to-end on the PoC #57 Stage 1 input
(50 violations) reproduces ~2 violations in a single call (the
chained PoC path reaches 1 with extra hand-tuned vertex_remove
ordering — see [Known limitations](#known-limitations)).

## Algorithm summary

`phase_h_finish` runs three stages in sequence. The same `seed`
flows through all three so the run is bit-reproducible.

1. **Stochastic local fixer.** For each currently-failing element,
   try up to `max_tries_per_fail` (default 500) random Gaussian
   perturbations of a randomly-chosen vertex of the triangle.
   Accept iff the perturbed mesh satisfies every FVCOM criterion
   in the **1-ring of the moved vertex** (the exact set of
   elements whose geometry changes). Boundary nodes get a
   1-D Gaussian along the boundary tangent, clamped to
   `[0.05, 0.95]` of the segment and projected onto the coastline
   if a projector is supplied. Outer-pass loop terminates on
   no-progress.

2. **C1 `vertex_remove` chain.** For each surviving C1 fail
   element, try `_apply_vertex_remove(force=True)` on each of its
   interior (non-boundary) vertices in descending-valence order,
   followed by the stochastic fixer on the resulting mesh. Accept
   the first candidate that strictly decreases the global total.
   The compositional principle is: the `vertex_remove` frees
   geometry the fixer cannot escape, and the fixer absorbs the
   local Delaunay artefacts the topology change introduces.
   Iterate until no candidate improves.

3. **C4 `vertex_remove` chain.** Same as (2) but targeting C4
   fail-pair interior vertices.

## Parameters

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `seed` | `42` | RNG seed flowing through every random draw; the entire run is bit-reproducible at fixed `seed`. |
| `min_angle_target` | `30.0` | FVCOM C1 strict threshold (degrees). |
| `max_angle_target` | `130.0` | FVCOM C2 strict threshold (degrees). |
| `area_ratio_target` | `0.5` | FVCOM C4 threshold on `\|A_i - A_j\| / max(A_i, A_j)`. |
| `max_valence` | `8` | FVCOM C5 legacy cap. |
| `alpha_target` | `0.95` | Alpha-quality target inside `_apply_vertex_remove` (bypassed via `force=True` in this driver). |
| `max_tries_per_fail` | `500` | Per-element random-move budget in Stage 1 + cleanup. |
| `perturbation_sigma` | `0.30` | Gaussian σ as a fraction of the mean local edge length. |
| `max_outer_passes` | `5` | Stochastic-fixer outer-pass cap (each walks every current fail once). |
| `coastline_projector` | `None` | Optional callable mapping a proposed 2-D position to the nearest coastline polyline. Returned by `build_coastline_projector`. |

The defaults match the PoC #58d / #58j / #58k chain that produced
the Tokyo-Bay 77 → 1 result. Tuning is rarely needed — see
"Tuning notes" below.

## Diagnostics

`phase_h_finish` returns `(updated_mesh, info)` where `info` is a
dict with:

```python
info = {
    "seed": 42,
    "min_angle_target": 30.0,
    "max_angle_target": 130.0,
    # ... every parameter echoed back ...
    "before":         {"NP": 49_036, "NE": 86_965, "C1": 2, "C2": 0, "C4": 48, "C5": 0},
    "stage1_stochastic": {"n_outer_passes": 2, "n_fixed": 86, "n_stuck": 18},
    "after_stage1":   {"NP": 49_036, "NE": 86_965, "C1": 1, "C2": 0, "C4": 4,  "C5": 0},
    "stage2_c1_vertex_remove": [
        # one record per attempted candidate (accepted or rejected):
        {"target_kind": "C1", "vertex": 46_234, "result": "applied",
         "before_total": 5, "after_total": 5, "after_counts": {...},
         "info": {"vertex_removed": 46_234, "rim_size": 7, "n_new_elements": 5,
                  "cleanup_stats": {...}}},
        # ...
    ],
    "after_stage2":   {"NP": 49_036, "NE": 86_965, "C1": 1, "C2": 0, "C4": 4,  "C5": 0},
    "stage3_c4_vertex_remove": [...],
    "after":          {"NP": 49_036, "NE": 86_961, "C1": 0, "C2": 0, "C4": 2,  "C5": 0},
    "delta_total":    -48,
}
```

The per-stage `after_*` counts let you locate which stage made
progress on which criterion. The `stage{2,3}_*` records are
useful for post-mortem: each candidate's `before_total` and
`after_total` make the greedy decision auditable.

## Known limitations

These are documented limitations of the current chain that you
should be aware of when interpreting results or designing
follow-up work.

### 1. Structural coastline floor (irreducible)

The final residual on the Tokyo-Bay benchmark — 1 C4 fail at
ratio 0.536 — is **geometrically irreducible without coastline
editing**:

* every vertex of the failing element pair lies on the
  coastline;
* the quadrilateral they form is concave, so the alternate
  diagonal that `edge_swap` would propose gives a CW
  (negative-area) triangle and is rejected on validity;
* `vertex_remove` cannot touch boundary vertices (would change
  the coastline geometry).

Closing this last fail requires moving a boundary node to a
different position on (or off) the coastline — i.e., a
coastline-data edit, which is out of scope for automatic mesh
repair. We treat 0.00115 % residual as "effectively zero" for
the SMS-phase-out goal.

### 2. Greedy vertex_remove ordering may stop short of the optimum

`phase_h_finish` accepts the first improving candidate per
iteration. On the Tokyo-Bay benchmark the chain reaches 2
violations, whereas the hand-tuned PoC pipeline (PoC #58j →
#58k) reaches 1 via a specific candidate ordering. The
difference is one violation (0.00115 % of the mesh) and falls
inside the stochastic-fixer seed-variance band ([4, 5] across
seeds 42–46 in PoC #58h).

A future revision could replace the greedy "first improving wins"
with an exhaustive "try every candidate and pick the global
best per iteration" search. The cost would be ~8× more cleanup
runs per iteration but should close the gap.

### 3. The insert operator never accepts under the current gate

PoC #58e / #58f wired a "random barycentric vertex insertion"
operator into the stochastic fixer at 30 % dispatch weight. It
accepted **zero** of ~14 000 attempts because the strict
30°-min-angle gate forces at least one of the three sub-triangles
to land near the threshold for any random barycentric position.
We removed the insert dispatch from the production driver — the
chain is move-only at Stage 1.

A useful insert operator would need a smarter proposal (Ruppert
circumcenter, edge-swap-equivalent) — see future work.

### 4. Chain assumes a reasonable starting mesh

`phase_h_finish` is a *finishing* step, not a meshing step. It
assumes the input mesh:

* has correct boundary topology (open / land segments listed in
  the `fort.14` correspond to actual boundary edges);
* has at most a few percent of elements failing the gates (the
  500-try budget per fail and 5-outer-pass cap presume "the
  hard work is done");
* has been through `phase_h_optimize` or an equivalent
  deterministic finishing pipeline.

If the input has 30+ % fails (e.g., raw OCSMesh + gmsh output),
run `fmesh-mesh-pipeline` first to bring it into the
finishing-chain regime.

### 5. RNG path depends on seed-and-mesh-state

The chain is bit-reproducible at a fixed seed *for a fixed input
mesh*. Two different inputs with the same nominal residual count
will follow different random streams (because the seed sequence
hits different elements at different times) and converge to
different residuals. The PoC #58h seed sweep confirmed the
variance is small ([4, 5] on the Tokyo-Bay input) — but on a
new mesh, plan to either:

* trust seed=42 as the canonical run; or
* sweep seeds {42, 43, 44, 45, 46} and pick the lowest residual.

## Tuning notes

The defaults reproduce the Tokyo-Bay 50 → 2 result; tuning is
rarely needed. Knobs that *can* matter:

| Knob | Effect of raising | When to consider |
|------|-------------------|------------------|
| `max_tries_per_fail` | larger search budget per fail; ~linear wall-time cost | residual contains "stubborn" pairs that the default 500-try budget misses |
| `max_outer_passes` | more chances for Stage 1 to re-visit newly-touched fails | residual has a chain of fails where fixing one exposes another |
| `perturbation_sigma` | larger random jumps; may flip triangles more often | residual is dominated by boundary-tangent corner cases where small moves don't escape |
| `seed` | different random stream; same expected residual distribution | the default seed gives an unlucky outlier on your input |

Lowering `min_angle_target` below 30° or raising
`area_ratio_target` above 0.5 is **not recommended** — those are
the FVCOM strict gates the chain is designed to satisfy.

## Future work (open candidates)

This is the prioritised follow-up list for the finishing chain.
See [`MEMORY.md`](../MEMORY.md) entries for the project's current
open-work menu and [`phase_h_finishing_chain.md`](phase_h_finishing_chain.md)
for the design-level open questions.

1. **FVCOM run-time tolerance test on the 1-violation mesh.**
   Does FVCOM accept a 0.536 area ratio (vs the 0.5 manual gate)
   as a soft warning, or does it crash / refuse to run? If clean,
   the auto-pipeline result is "production-ready" for SMS
   replacement. If FVCOM rejects, coastline editing becomes the
   next required step.
2. **Generality validation on Osaka / Ariake / Sendai / Mikawa.**
   77 → 1 is a single-mesh result. Three other basins with
   different geometry character (Ariake wide-shallow, Sendai
   open-Pacific, Mikawa multi-channel mouth) would confirm the
   chain works as a general method, not Tokyo-Bay coincidence.
3. **Exhaustive vs greedy vertex_remove ordering.** Replace the
   first-improving-wins policy with a "try every candidate per
   iteration and accept the global best". Cost ~8×; expected
   payoff is closing the 1-vs-2 gap on Tokyo-Bay-class meshes.
4. **Ruppert-style insert operator.** Replace the (currently
   0-accept) random-barycentric insert with a circumcenter-based
   proposal that places the new node on the perpendicular
   bisector of the fail-element's longest edge. May make the
   insert operator productive on cluster residuals that the
   move-only fixer cannot escape.
5. **Coastline editing pass (Stage 4).** For users who prefer
   mesh quality over strict coastline fidelity, an optional
   stage that collapses one of the 4 vertices of an
   all-boundary concave C4 pair onto a neighbour would close
   the structural floor. The collapse changes the coastline
   slightly (typically a few cm to a few m in geographic
   coordinates) — defensible for hydrodynamic simulation but
   not for cartography.
6. **Integrate `phase_h_finish` into `fmesh-mesh-pipeline`.**
   Today the pipeline rungs stop at `phase_h_optimize`; adding a
   "rung 4 = `phase_h_finish`" gates with the existing
   `fmesh-mesh-quality` thresholds would make the finishing
   chain accessible from the CLI without writing Python.
7. **Public CLI: `fmesh-phase-h-finish`.** A thin click-based
   wrapper around `phase_h_finish` so users can run the chain on
   an existing `fort.14` without touching Python (cf.
   `fmesh-mesh-clean`).

## References

* [`docs/phase_h_finishing_chain.md`](phase_h_finishing_chain.md)
  — design rationale, ASCII flow diagram, critical constants,
  compositional principle, structural-floor argument.
* `MEMORY.md` entry "stochastic-local-fixer-breakthrough" —
  empirical journey notes (this session's record of how the
  chain was discovered and refined).
* Notebooks `notebooks/58d_...py`, `58e_...py`, `58j_...py`,
  `58k_...py`, `58l_...py` — original PoC scripts. Useful for
  reproducing the journey or running variants of the chain.
* Implementation: `src/fvcom_mesh_tools/mesh_clean_phase_h.py`,
  function `phase_h_finish`.
* Tests: `tests/test_mesh_clean_phase_h.py`, the
  `phase_h_finish_*` test group.
