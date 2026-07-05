# Design History: Boundary Conformity, Feature Inclusion, and Quality

A post-mortem ledger of every approach tried in the PoC #59-#93
series (2026-07-04/05) for the three coupled problems of coarse
(~300 m) FVCOM meshing of Tokyo Bay:

- **P1 Conformity** — mesh boundary must lie ON the coastline where
  the coastline is representable at the mesh scale;
- **P2 Feature inclusion** — Tokyo Port, river mouths and other
  essential-but-narrow water must be meshed, "forcibly if necessary"
  (shrinking or erasing artificial structures), sacrificing fidelity;
- **P3 Quality gates** — the FVCOM manual's SMS criteria (min angle
  >= 30 deg, max <= 130 deg, area change <= 0.5, valence <= 8, OBC
  triangle edge normal to the boundary), which the hand-finished
  goto2023 reference mesh satisfies almost exactly.

Failed approaches are retained deliberately: each failure isolated a
mechanism, and the final architecture is the sum of those mechanisms.

## 1. The four mechanisms that explain everything

These took ~30 PoCs to isolate; every failure below traces to one of
them.

**M1 — feature sizing has no floor (oceanmesh).**
`feature_sizing_function` refines to ~h0/4 near coastal features, so
every "300 m" build was actually a 75-150 m coastal mesh (~60k
nodes). The TRUE hmin=300 Tokyo Bay mesh is ~5-6k nodes. Fix: clamp
the sizing grid VALUES after construction
(`--om-enforce-hmin-floor`). The function's own `min_edge_length`
argument is NOT a floor (passing hmin there collapses the mesh —
PoC #65 was misread as a bug for half a day).

**M2 — the cleanup retreat.**
The persistent boundary-to-coastline offset is ~0.7 x the local
element size at EVERY resolution (60c: 51 m @ 75-150 m; #63: 19.8 m
@ 50-100 m; #77: 292 m @ 300 m). It is NOT shoreline simplification:
post-generation cleaning (`delete_boundary_faces`,
one-face deletion) removes the sawtooth boundary slivers and thereby
pulls the boundary ~1 element inside the domain's zero line.
Consequence: generation without edge constraints can never conform;
and any conformity mechanism must also disarm these cleaners.

**M3 — the DistMesh seed lattice equals the sizing minimum.**
With a true 300 m floor, initial points are laid ~300 m apart, so
water narrower than ~2 floors (600 m) receives no seed points and
silently drops out — Tokyo Port basins, river channels. Narrow
essential water therefore needs either explicit interior seed
points or artificial widening.

**M4 — the Python port lacks constrained Delaunay (egfix).**
MATLAB OceanMesh2D's high-fidelity mode pins resampled shoreline
points AND EDGES into the triangulation (MATLAB
`delaunayTriangulation` supports constraints natively). The Python
port had only unconstrained CGAL DT, so fixed points drift topologically
even when their positions are locked — chains do not survive as
edges. This is the root reason all point-based (pfix-only) and all
post-hoc approaches hit ceilings.

## 2. Approach ledger

### P1 Conformity — failures first

| # | Approach | Result | Failure mechanism |
|---|---|---|---|
| #61 | pfix from `Shoreline.mainland` (simplified) | p50 46 m (baseline 51 m) | mainland/inner are h0-simplified; constraints cannot beat the simplification floor |
| #62 | pfix from RAW OSM polylines | p50 56 m (WORSE) | fixed points off the SDF domain get culled with their triangles; domain geometry, not constraint mechanics, governs (M2/M4) |
| #63 | domain detail h0=100 m, sizing 300 m | p50 19.8 m but 201k elements | works because coastal elements shrink (M1 unfloored); boundary vertices pack at h0 scale — cost, not conformity, fails |
| #69 | snap ALL boundary nodes (frac 1.2), repair after | p50 8e-11 m; mesh wrecked (1,559 violations; optimize deleted 29% of nodes) | snap-then-globally-grind ignores that conformity edits must be locally validated (SMS discipline) |
| #70 | per-node quality-gated snap | only 37% snapped, p50 65 m | gate judged one node while neighbours had not moved: collective motion misread as shear |
| #71 | chain-collective snap (40-node chains) | 9% of chains accepted | all-or-nothing: one bad quad rolled back 39 good nodes |
| #72 | + bisection narrowing | 49% nodes snapped, p50 57 m | correct granularity, but ~half the nodes cannot conform by MOTION alone at C1 >= 30 (need topology change) |
| #73 | + on-line slide of snapped neighbours | +3% only; 2 CCW flips (bug: slide gate must cover the slid node's FULL incident ring, not just the chain patch) | marginal: the motion ceiling again |
| #79 | boundary strip extrusion (fill the retreat gap with a new triangle row ON the line) | p50 1.7e-10 m, fatal gates PASS — but C1/C4 debris at strip seams; extrusion v1 trimmed whole runs (58 lost) until segment-level trimming (v2: 62 strips) | workable but is a patch OVER M2 rather than a fix OF it; quality debt (~230 violations) handled downstream |
| **#92** | **CDT egfix: shoreline chains constrained into the triangulation** | **p50 7 mm at generation** | — (this is the fix of M4; requires disarming M2's cleaners, see #92 v2 below) |

Verdict: post-hoc conformity (snap/slide/extrude) is a dead end the
user correctly called: it produced jagged boundaries and unbounded
repair debt. Generation-time edge constraints are the answer;
everything else in the snap family is retained only as small QA
repair operators.

Sub-failure worth remembering — **#92 v2**: with egfix in place the
boundary STILL sat 229 m inside. The CDT held the chains through
every retriangulation; the post-generation cleaners then deleted the
constrained boundary row (M2 again). In constrained mode only
`make_mesh_boundaries_traversable` may run.

### P2 Feature inclusion

| # | Approach | Result | Failure mechanism |
|---|---|---|---|
| #64-#66, #76-#77 | morphological opening of the WATER (erode/dilate water r=150 m) to delete sub-grid channels | Tokyo Port and all river channels VANISHED; user: "実用にならない" | opening removes thin WATER — exactly the essential features; also shifts 300-600 m-scale water boundaries so even smooth-coast fidelity degraded (#67 metric) |
| #75/#76 | raw or opened domain + hmin floor | ~5.7k nodes but nearshore/port empty | M3: 300 m seed lattice cannot seed <600 m water regardless of domain preprocessing |
| #90 | morphological opening of the LAND (erode/dilate land r=150 m) | port basins + Arakawa survive; thin piers/breakwaters/islets erased | — (correct operand; matches the goto2023 hand policy "erase thin artificial structures, keep the water") |
| #91 | + water-skeleton seeding (medial axis where half-width 140-460 m -> interior pfix lines) | 300-600 m basins/channels meshed | — (independently mirrors ADMESH+ 2024; see literature survey) |

Verdict: the operand of the morphology was the whole story — open
the LAND, never the water; then seed the remaining narrow-but-wide-
enough water explicitly (M3).

### P3 Quality repair operators (on residual defects)

| Operator | Target | Result |
|---|---|---|
| patch re-CDT (#81) | any quality cluster | 1/10 accepted; also taught the stale-index lesson: an accepted edit renumbers elements — later sites must re-resolve targets by NODE ids |
| `collapse_edge` (#83) | boundary needle-ears (2 boundary edges + short edge) | C1 18 -> 6, C2 cleared: the right tool for needles |
| `split_edge_pair` / `grade_region` (#84-#85) | size-cliff C4 pairs | ~0 accepted: every split halves angles -> new C1; splitting is the WRONG tool at a size cliff |
| `equalize_pair` (#86) | size-cliff C4 pairs | 9 accepted: move the two apex nodes (big toward the shared edge, small away) — the actual SMS "drag the interior node" move |
| cap-triangle deletion (#93) | CDT boundary slivers (all 3 nodes on the chain) | 538 removed with conformity INTACT (deletion exposes the constrained edges) |

### Process failures (equally load-bearing)

- **Iterate-until-gates-pass does not converge** (59k5 hung 6.8 h in
  lookahead; #68 ground >1 h silently; 60d hit its elapse limit).
  SMS practice never expects convergence either. Replaced by the
  session model: every stage runs ONCE, bulk optimization gets a
  10-minute budget + stagnation stop (accepts < 2% of first sweep),
  heartbeat logging makes alive-vs-hung decidable from the job log,
  residual sites are reported with figures and edited individually.
- **Aggregate fidelity percentiles mislead** (#67): p50 over all
  boundary nodes mixes representable and sub-grid sections; judge
  conformity on smooth sections, and judge engineered meshes against
  the ENGINEERED shoreline, not raw OSM.
- **Two changes per experiment cost a day** (#65 conflated the
  min_edge_length argument with the island threshold; the argument
  was blamed for what was actually the true-coarse node count).
- **Stage PNGs are mandatory** (user directive): every pipeline
  stage saves a mesh figure so failures are locatable by inspection
  (see `outputs/figures/87_gallery`, `93_s*`).

## 3. Current architecture (v5) and reference targets

Pipeline: land-opening (r=150 m) -> water-skeleton seed lines ->
oceanmesh build with h0=100 m domain detail + hmin floor +
`--om-constrain-boundary` (CDT egfix chains) -> constrained-mode
cleanup (manifold fix only) -> cap-sliver deletion + weld + one
budgeted projector optimize -> per-site operator suite -> QA gates
-> export -> M2 tidal test.

State at 2026-07-05 morning: conformity p50 0.27 m, C1 244 / C2 9 /
flips 0 at 10,012 nodes. Reference (goto2023, OM2D+SMS,
hand-finished): 3,210 nodes, edges 327/561/1777 m (p1/p50/p99), min
angle 29.0 abs, area change <= 0.50, valence <= 8, one smooth OBC
arc at the Uraga narrows (south limit 34.97 N).

Open items: residual narrow-channel slivers (collapse/equalize,
constraint-aware); weld/optimize keeping boundary nodes exactly
on-line; OBC arc via domain-polygon cut; node budget (10k vs 3.2k:
interior growth + skeleton density); rivers beyond the Arakawa; full
QA + M2 on the v5 mesh.

Related documents: `LITERATURE_SURVEY_AUTOMESH.md` (no published
FVCOM automesher; RiverMapper pseudo-channels and ADMESH+ width
classification as adoptable prior art), `gate_passing_pipeline.md`
(v1-era pipeline; superseded by the architecture above),
`fvcom_source_constraints.md` (the QA gate derivations).

## 4. The Yokosuka seam campaign (2026-07-05/06): a cascade failure

**User verdict (2026-07-06): the mesh BEFORE this campaign was much
better; the entire fix stack was withdrawn and `src/` reverted to
the review17 state (commit `2e4d89d`, CFL floor). The 14 campaign
commits remain in history as the failure ledger.**

### The original defect (a preprocessing bug)

`simplify_outside_region` (introduced at `a187715`, review16) cut
each land polygon at the interest-region boundary and Douglas-
Peucker-simplified the OUTSIDE piece — including the CUT-EDGE
vertices. The two halves then no longer shared a seam, leaving a
sliver strip (hundreds of metres wide, km long) that belongs to no
land polygon. In the engineered-shoreline product that strip is
WATER, so the mesher meshed it faithfully: "mesh on OSM land" at
Yokosuka/Kurihama, exactly where the interest boundary crosses
land. Root cause class: PREPROCESSING — the mesher and the OBC
machinery were never at fault.

### The cascade: each fix made the open boundary worse

1. Seam-exact rewrite + morphological smoothing (r=400 m) +
   set_precision (`cfaf866`..`d0f91f5`, reviews 18-23): fixed the
   land sliver, but MOVED THE COAST near the OBC junctions. The
   junction geometry shift destabilised the (latent, fragile)
   OBC list machinery — coordinate-keyed membership went stale
   after siteops/polish motion and the OBC stopped reaching the
   coasts (review23; unmarked boundary = artificial wall in FVCOM;
   caught by the user).
2. Repair surgeries on the first row then cascaded (reviews 24-33):
   - geometric membership derivation + ring-path ordering: held
     (coast-to-coast restored) but exposed off-line stragglers;
   - R4 centroid splits: children re-flag -> 140 splits, C1 317;
   - fake-open splits: same cascade class;
   - on-line edge collapse: gates too tight, 1-2 firings only;
   - corridor purge + strip re-extrusion: extrusion laddered only
     2 strips into a purged 22 km corridor -> d p50 = 1 km.
3. Net effect vs review17: OBC quality strictly worse at every
   step; reverted wholesale.

### Lessons (binding)

- **Fix preprocessing bugs IN PREPROCESSING.** A broken input
  (fake-water sliver) must be repaired where it was created; every
  downstream compensation layer added a new failure mode.
- **Boundary-touching surgery cascades.** Deleting elements at the
  OBC retreats the boundary; splitting them re-flags children;
  "repair" operators that run inside a convergence-style loop are
  the anti-pattern the no-convergence-loop doctrine forbids —
  meta-level included (iterating OPERATOR DESIGNS against one
  residual is the same trap).
- **Geometry changes near junctions invalidate downstream
  assumptions.** Any prep change that moves the coast near an
  OBC junction must be followed by a full-pipeline visual check of
  the junctions before anything else is tuned.
- **Revert early.** The correct response after the second failed
  surgery was a wholesale revert; instead nine rebuilds were spent.

### The planned correct fix (not yet implemented)

Seam-free formulation: simplify the coastline as RING ARCS, not
area pieces — split each polygon ring into arcs inside/outside the
interest region, DP-simplify only the outside arcs with their
endpoints (the region-boundary crossing points) pinned, and
reassemble the ring. No area cut -> no seam -> no fake water, and
the smoothing/precision/validity patches become unnecessary.
