# Stage 2: Automated Finishing (SMS replacement) — Design

Status: DRAFT for approval (2026-07-06). Owner: fvcom-mesh-tools.
Prerequisite: OM2D port P1-P4 complete in the oceanmesh fork
(remesh_patch, merge_meshes, mesh_improve suite, boundary tools).

## 1. Goal

Replace the manual SMS finishing pass that produced goto2023 with
a **recipe-driven, QA-driven, per-patch local rebuild layer** that
is domain-agnostic and fully reproducible. This capability exists
in neither SMS-automation form nor OceanMesh2D; it is new work.

Non-goals: global regeneration loops (policy: no
iterate-until-gates-pass), depth-field design, interactive UI.

## 2. What SMS hand-finishing actually did -> operator mapping

| SMS manual action            | Stage-2 operator                              |
|------------------------------|-----------------------------------------------|
| redraw a bad cluster         | `remesh_patch` (reconstructed + target sizing) |
| merge/split nodes, swap edges| `collapse_thin_triangles`, `flip_edge`, `bound_connectivity` (P2) |
| smooth a neighbourhood       | `direct_smoother_lur` on the patch (P2)       |
| nudge boundary nodes         | boundary-slide along coastline / OBC arc (whitelist) |
| local resolution override    | recipe *directive* -> `remesh_patch(target_h)` |
| delete junk (pinch, spikes)  | structural ops (pinch removal, dj components)  |

Already automated upstream of stage 2 (v6): OBC arc snap +
perpendicular junction (egfix), projector snap, CFL sizing.
Stage 2 therefore focuses on QUALITY-CLUSTER repair + USER
DIRECTIVES.

## 3. Architecture (detector -> planner -> executor -> reporter)

```
QA report (json)          recipe finishing: directives
      |                          |
  [Detectors]                    |
  C1/C2/C4/C5 sites, pinch,      |
  disjoint, valence              |
      \________________________ /
                |
            [Planner]
  cluster violations -> patches (k-ring dilation, k=2;
  overlapping patches merged); classify each patch:
    micro  (1-2 elems, interior)  -> P2 micro ops
    patch  (cluster >= 3 elems)   -> remesh_patch
    bound  (touches coast/OBC)    -> whitelist ops only
                |
            [Executor]  (UTM, per-patch atomic)
  apply operator -> LOCAL acceptance: patch+halo gates
  (C1>=30deg, C2<=130deg, C4<=2.0 incl. cross-halo pairs,
  C5<=8, no flips, boundary nodes on their lines) ->
  accept or REVERT THIS PATCH (others unaffected)
                |
            [Reporter]
  per-patch ledger json + before/after zoom figures;
  unfixed patches surfaced for user judgment
```

Design rules (binding, from project history):
- **Single pass**: each detected patch is attempted once per
  pipeline run. No re-detection loop after execution; the QA stage
  reports what remains.
- **Cross-halo acceptance** (lesson of PoC #129 v1): acceptance
  evaluates pairs spanning the patch/halo boundary, not only
  inside.
- **Boundary whitelist**: coastline/OBC nodes move only along
  their line or are deleted; junction/OBC geometry is stage-1.5
  property and stage 2 must not touch OBC members.
- **Determinism**: fixed seed per patch id; same input -> same
  output.
- **Inset priority in stitching**: `merge_meshes` keeps the
  rebuilt patch verbatim and smooths only the seam band.

## 4. Recipe interface

```yaml
finishing:
  auto:
    checks: [c1, c2, c4, c5, pinch, valence]   # detector set
    patch_kring: 2
    max_patches: 50            # hard cap, no silent loops
  directives:                  # SMS-style local edits, declarative
    - polygon: [[139.76, 35.54], ...]
      target_h_m: 150          # local refine (port channel)
    - polygon: [[139.7, 35.0], ...]
      op: coarsen
      target_h_m: 1200
```

Directives execute BEFORE auto-repair (they change sizing;
auto-repair then heals their seams if needed).

## 5. Pipeline placement

`build -> finish -> obc -> [stage2 finishing] -> obcfinal -> qa ->
export`. Replaces today's `siteops` + `polish` + ad-hoc notebooks
(#129 class). Those stages remain available behind recipe flags
during transition, then deprecate.

## 6. Acceptance criteria

1. **Tokyo Bay parity**: recipes/tokyo_bay_v6.yaml with
   `finishing.auto` reproduces all-gates-pass with the ad-hoc
   sitefix notebook deleted; mesh differences confined to patch
   neighbourhoods.
2. **Generality**: an unseen second domain (candidate: Ise Bay,
   same OSM/DEM tooling) runs prep->build->stage2->qa with NO
   domain-specific code and reaches 0 quality-gate failures or a
   reported residual list.
3. **Directive demo**: one refine directive (Tokyo port channel
   150 m) meshes, passes local gates, and survives M2.

## 7. Implementation plan

- M1 detector+planner (analysis only, ledger + figures; no mesh
  writes) — validates clustering on v6.
- M2 executor: micro ops + remesh_patch + acceptance/revert.
- M3 recipe directives.
- M4 generality run (Ise Bay) + A/B note.

Each milestone: pytest (happy + error), English docs, single
review gate with the user before the next milestone.
