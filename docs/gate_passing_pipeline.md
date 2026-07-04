# The gate-passing pipeline — from Phase-H endpoint to FVCOM-ready inputs

> **Superseded (2026-07-05):** this documents the v1-era
> pipeline (post-hoc OBC arc + local perp fixing on the
> fine-coastal mesh). The current architecture — CDT-constrained
> boundary, land-opening preprocess, skeleton seeding, session
> cleanup — is described in `DESIGN_HISTORY.md`.


Status snapshot 2026-07-04. This documents the PoC #59 chain that took
the Phase-H quality endpoint (`outputs/58l_chained.14`, the "77 → 1"
mesh) to the first mesh that passes **every** `fmesh-mesh-qa` gate,
plus the export to FVCOM native inputs. Companion docs:
`fvcom_source_constraints.md` (the rules), `fvcom_mesh_kickoff.md`
(the acceptance spec), `phase_h_finishing_chain.md` (how 77 → 1 was
reached in lon/lat space).

## Why the "finished" 58l mesh failed 10 / 21 gates

The first `fmesh-mesh-qa` run exposed three classes of problems the
quality-focused Phase A–H work never gated on:

1. **Coordinate space.** Phase H optimized angles in raw lon/lat; at
   35.5°N the cos(lat) anisotropy shifts angles by up to ~5°. FVCOM
   production builds are CARTESIAN, so the QA gate evaluates a local
   metric projection — where 58l had 55 C1 violations (invisible in
   lon/lat) — and, conversely, the lon/lat-space "irreducible" C4
   floor of 1 does not exist in metric space.
2. **Open-boundary structure.** The lineage's OBC sat on the eastern
   140.10E data-clip line near the Chiba coast: 15 R4 elements
   (`tge.F` PSTOP), 11 fake-ISBCE=2 elements, necked OBC nodes,
   perpendicularity to 77.6°. FVCOM would refuse to start.
3. **Never-applied housekeeping.** 44 components / 4,878 disjoint
   elements, 2 orphan nodes, 9,688 nodes above the 2 m depth clip.

## The chain (all seeds fixed, wall times on 1-2 GENKAI cores)

| PoC | Action | Result |
|---|---|---|
| 59a | keep largest component → iterative R4/fake-ISBCE=2 deletion + bbox boundary rebuild → 2 m clip (lon/lat) | all structural gates pass, 8.7 s |
| 59b | EPSG:32654 projection → perpfix → `phase_h_finish` (metric) | {C1 79, C2 4, C4 41} → **0/0/0/0** in 12 s |
| 59c | convergence loop (finish ⇄ damped global perpfix) + `compact_nodes` | 20/21 gates; global perpfix ping-pongs with finish at the eastern clip |
| 59d | **OBC relocated to the Uraga Channel transect** (southern DEM edge 35.10N, 50 nodes, depths to 338 m); eastern clip closed as land | 20/21; perp residual down to 8 nodes / 29.0° |
| 59e | **quality-gated local perp fixer** (move only violating nodes' first-ring partners; accept iff perp fixed AND 1-ring C1/C2/C4 hold AND neighbours stay legal) | **21/21 PASS**, 8/8 nodes in one pass, 5.8 s |

Final mesh: `outputs/59e_gate_passed.14` — 45,436 nodes / 81,603
elements, UTM 54N, min angle exactly ≥ 30°, max angle 116.5°,
C4 ≤ 0.499, valence ≤ 8, OBC best-edge deviation ≤ 19.9°
(mean 10.6°), implied external Δt 1.05 s.

Deliverables: `outputs/fvcom_inputs/tokyo_bay_v1_{grd,dep,obc,cor,spg}.dat`
+ `tokyo_bay_v1.2dm` + `tokyo_bay_v1_provenance.md`
(via `fmesh-export-fvcom --cor crs --crs EPSG:32654`); figures under
`outputs/figures/` (via `fmesh-plot-mesh`).

## Lessons now encoded in the toolkit

* **Gate in the space FVCOM runs in.** `fmesh-mesh-qa` evaluates
  metric geometry; any future finishing must run after projection
  (or the build should work in UTM throughout — recipe-layer item).
* **Global perpfix and quality finishing fight each other.** The
  fix is locality: `phase_h_finish` for quality, then the #59e-style
  targeted mover for the few perp violators. Promote #59e into a
  toolkit function (`align_open_boundary_local`?) — Phase-2 item.
* **`phase_h_finish` leaves orphan nodes** (vertex_remove keeps the
  node array); `mesh_clean.compact_nodes` repairs, and the native
  writers refuse to emit orphans. Consider compacting inside
  `phase_h_finish` itself (backward-compatible: NP shrinks).
* **OBC placement is a first-class decision, not a bbox
  side-effect.** The bbox-proximity classifier silently produced an
  unusable OBC; domain forensics (ring walk + lat banding) found the
  real transect. The recipe layer must carry an explicit
  open-boundary specification.

## FVCOM run-time acceptance — PASSED (PoC #59g)

Cold-start, zero-forcing smoke run of the exported `tokyo_bay_v1`
inputs on the TB-FVCOM hydro binary (FVCOM 5.1 lineage, 4 MPI ranks,
EXTSTEP 0.5 s ≤ the QA-implied 1.05 s, 2 simulated minutes,
`OBC_ELEVATION_FORCING_ON=F`, `BOTTOM_ROUGHNESS_KIND='constant'`,
2-line uniform `sigma.dat`, nml cloned from
`hydro/tuning/new_bc/tb_pge_run.nml`):

* `ISONB SETTING / NISBCE/LISBCE/EPOR: COMPLETE` — every `tge.F` and
  `SETUP_OBC` PSTOP that `fmesh-mesh-qa` gates on was accepted;
* all 24 internal steps ran to `END_DATE`, normal completion
  (`TADA!`), no NaN;
* at-rest flow stats sane (max internal velocity ~2 cm/s = the
  expected sigma-PGE truncation signal, free surface ≤ 1 cm).

The mesh-QA gate's prediction of FVCOM startup acceptance is now
verified against the real binary; the auto-pipeline mesh runs with
**no manual SMS step anywhere in the chain**. Case directory:
`outputs/fvcom_smoke/` (job script `notebooks/59g_fvcom_smoke.pjsub`).
Note the production namelist confirms `PROJECTION_REFERENCE =
'proj=utm +zone=54 +datum=WGS84'` — the EPSG:32654 datum choice
matches the existing TB-FVCOM convention.

## Open items toward the kickoff Definition of Done

* **User decisions to confirm** (also in the provenance file):
  eastern 140.10E clip closed as land wall (vs secondary clamped
  OBC); whether the OBC should span the entire southern 35.10N cut
  (incl. the shallow Sagami-side shelf) instead of only the deep
  channel run; OBC type 1 uniform.
* recipe.yaml layer (§3), depth-field regeneration from the
  DATA_INVENTORY precedence stack (the current depths predate it),
  narrow-channel inclusion criterion (§7.4 — w/h currently
  advisory), boundary conformity vs coastline polygons (§9-2),
  OBC-line straightening (one ~200 m notch remains at x≈390.5 km —
  gate-passing but §7.2 prefers a smooth arc).
