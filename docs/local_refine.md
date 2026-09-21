# Local refinement of an existing mesh — design

**Status: specification implemented and the first worked example buildable;
the generator is designed but not written.** This document is
revised as the design moves, and records the decisions and their reasons so
that a later change is made knowingly rather than by accident.

Revision history is the git history of this file.
`docs/local_refine_review.md` is an adversarial review of revision 1
(gpt-6-astra, 2026-09-22); every blocker it raised was reproduced
independently before this revision was written. Revision 2 corrects four
factual errors and records what the review showed is still missing.

### What revision 1 got wrong

| claim in revision 1 | what is true |
|---|---|
| the Futtsu recipe is eligible | it is — but not for the reason given. Only the **core** was checked against land. The hole is 2,239 m and overlaps the coastline by 1,225 m; the `touch_coast: false` rule wrongly made that a refusal, when meshing a coastline at the target size is ordinary work (revision 3 replaces the flag with a `coastline` mode) |
| the region allows dt = 4.70 s | that is the **shortest-edge** figure. The reported measure is minimum altitude: **4.07 s**, below its own 4.5 s floor |
| the base mesh carries production depths | it does not: min 2.00 m (not the 3 m floor), max 735.3 m (not the 300 m cap), 2,258 edges above r = 0.2 |
| `pfix` constrains the rim | `pfix` fixes positions only. Segments need `egfix`, which is what drives the CDT (`mesh_generator.py:1077`) |
| the rim "already equals the ambient size" | measured rim edges are **180 / 467 / 819 m**, not 350 m. The claim was tautological and the compatibility it asserted is false |

## 1. What problem this solves

A user wants one part of the domain at a much finer resolution than the rest —
a fishery, a port basin, a discharge site — and wants **the rest of the mesh
left alone**.

The existing `regions:` block in a sizing recipe (`src/fvcom_mesh_tools/
sizing.py`, `docs/sizing_regions.md`) does not do this. It edits the **sizing
field**, after which the whole mesh is regenerated. DistMesh is a global
relaxation: perturbing the field anywhere moves nodes everywhere. The result
is a different mesh, not the old mesh with a refined patch.

So there are two operations, sharing one vocabulary for geometry and targets:

| | outside the region | cost | status |
|---|---|---|---|
| **A. sizing region** | changes | one full rebuild | implemented (`sizing.py`) |
| **B. local refinement** | **unchanged, and checked** | one patch | this document |

B is what the rest of this document specifies.

## 2. The three zones and the contract

| zone | definition | guarantee |
|---|---|---|
| **core** | inside the declared geometry | achieved edge length ≈ `target_h_m` |
| **transition** | the annulus the gradation needs to reach the ambient size | remeshed; changes |
| **frozen** | everything else | node coordinates and connectivity **identical** to the base mesh |

**The bathymetry is not refined.** The base mesh is the topography actually
being simulated, so a refined region sits on the same seabed: retained nodes
keep their exact depths and new nodes take the base field interpolated at
their position (`refine.depths_from_base`). Refining the bathymetry as well is
a separate question, deliberately left for later — mixing a finer seabed into
this change would make the effect of resolving the fishery impossible to
separate from the effect of changing the depths.

The frozen guarantee is the point of the exercise, so it is **checked, not
asserted**. `refine.frozen_changes()` is a start and not yet sufficient: it
compares coordinates in row correspondence, which only holds while the node
numbering is unchanged, and a patch inserts, deletes and renumbers nodes. The
contract needs an explicit old-to-new node map and comparison of retained
connectivity and orientation, depths, and boundary membership and order —
identity under a map, not equality of a serialized table. (A destroyed
coordinate used to pass, because `NaN > tol` is false; that is fixed and
tested.)

The footprint in the contract is also **not** the circle in the recipe. Whole
triangles are removed, so the actual cut reaches further: on the Futtsu
selection, 2,615 m against a requested 2,239 m. The generator must report the
removed-face set and its departure from the requested envelope before filling,
and any enlargement must stay inside a declared maximum.

## 3. Two facts the caller cannot see, and must not have to guess

### 3.1 The transition width is not a free choice

With a linear size gradation `g`, the edge length grows as
`h(r) = h_target + g · (r − r_core)`, so reaching the ambient size takes

```
W = (h_ambient − h_target) / g
```

For the certified Tokyo Bay chain (`g = 0.165`, ambient ≈ 350 m):

| target | elements across | transition W | W / core radius (300 m) |
|---|---|---|---|
| 15 m | 20.2 | 1,909 m | 6.4 |
| **30 m** | **15.7** | **1,939 m** | **6.5** |
| 50 m | 12.4 | 1,697 m | 5.7 |
| 100 m | 7.8 | 1,394 m | 4.6 |

**The transition is several times the size of the region itself.** The user
declares the target; the width follows. An explicit `transition_m` is checked
against `W` and rejected when the gradation cannot cover it, with the three
ways out named in the message (widen it, coarsen the target, steepen the
gradation).

Steepening the gradation shortens the transition but spends the C4 margin.
C4 is `(A_large − A_small)/A_large ≤ 0.5`, i.e. an **area** ratio ≤ 2; across a
shared edge that is an **altitude** ratio ≤ 2, not a universal edge ratio
≤ √2 (which would assume similar triangles). A size-gradient bound alone does
not certify C4. `g = 0.33` halves the width; beyond that, expect QA failures.

### 3.2 The time step is the cost that bites

FVCOM integrates with **one global external step**, so the smallest element
anywhere sets it for the whole run. **Which length goes in the numerator
matters**: the 11.9 s quoted for the certified mesh is the **minimum altitude**
of a triangle (notebook 392; `coast_fit` guards altitude and shortest edge
both), while the shortest edge gives a larger, flattering number.

```
dt = L / sqrt(g · H),  L = minimum altitude
```

For an equilateral cell the altitude is `sqrt(3)/2 = 0.866` of the edge, so a
target edge length buys 13 % less dt than the naive estimate. That is exactly
the margin revision 1 spent without noticing: 30 m over 4.15 m of water gives
4.70 s by edge and **4.07 s by altitude**, against a 4.5 s floor.

Worse, the quality gates do not bound the altitude. A 30°–30°–120° triangle
with two 30 m sides passes C1 and C2 and has a minimum altitude of 15 m — half
the equilateral value, so half the time step. **`preflight` assumes equilateral
cells and therefore reports an upper bound, not a guarantee**; the achieved dt
must be measured on the finished mesh.

Minimum altitude of an equilateral cell, in seconds:

| depth | h = 30 m | h = 50 m | h = 100 m | h = 350 m |
|---|---|---|---|---|
| 2 m | 5.86 | 9.78 | 19.55 | 68.4 |
| 3 m | 4.79 | 7.98 | 15.96 | 55.9 |
| 4.15 m | **4.07** | 6.79 | 13.57 | 47.5 |
| 5 m | 3.71 | 6.18 | 12.37 | 43.3 |
| 10 m | 2.62 | 4.37 | 8.75 | 30.6 |

The certified mesh allows 11.9 s. A 30 m region costs a factor of 2–4,
**depending entirely on the depth of that region**: shallow water is cheap to
refine, deep water is not. The step count is not the whole cost — added
elements, shape and solver settings matter too.

**The time step does not veto a region** (owner 2026-09-22). A fishery is
*given*: its position and the resolution it needs are inputs, not preferences,
so a recipe that cannot hold a time step is still the recipe. `dt_expected_s`
is therefore advisory — `refine.preflight()` raises an **alert** in its report,
naming the step the region will actually allow, the factor by which the run
gets longer, and the target that would have kept the expected step. Refusal is
reserved for what makes the operation impossible, not for what makes it
expensive.

### 3.3 Node count, for completeness

Core elements ≈ `area / (√3/4 · h²)`. A 300 m disc at 30 m is ~725 elements;
its transition adds ~830. Against a base mesh of 8,252 elements that is
**+16 %** for one small fishery, and the transition is the larger half.
Enlarging the core is comparatively cheap — the transition width does not grow
with it.

## 4. Specification

A refinement recipe is YAML, version-controlled, and lives in
`recipes/refine/`. Environment variables (`SR_*`) remain for sweeps; a durable
specification belongs in a file.

```yaml
base_mesh: ../../outputs/verify_409.115302/fit/sample_repro_final.14
dt_expected_s: 4.5            # advisory: an alert, not a veto
gradation: 0.165
coastline: resample           # preserve | resample | spline
coastline_tolerance_m: 100    # max departure from the base polyline

refine:
  - name: futtsu_nori
    geometry:
      circle: {center: [139.7881, 35.3228], radius_m: 300}
    target_h_m: 30            # achieved edge length, as elsewhere in the chain
    priority: 0
```

### Keys

| key | meaning |
|---|---|
| `base_mesh` | the mesh to patch; relative paths resolve against the recipe |
| `dt_expected_s` | mandatory to state, advisory in effect: an alert when the region will not deliver it |
| `gradation` | size growth per metre through the transition |
| `refine[].name` | unique, non-empty |
| `refine[].geometry` | `circle`, `bbox` or GeoJSON `Polygon`, all lon/lat |
| `refine[].target_h_m` | **achieved** edge length, matching `SR_H_TARGET=achieved` |
| `refine[].transition_m` | optional; derived when absent, checked when present |
| `coastline` | `preserve`, `resample` (default) or `spline`; see below |
| `coastline_tolerance_m` | max departure from the base polyline, default 100 m |
| `refine[].priority` | higher wins an overlap; ties take the smaller target |

Unknown keys, non-finite values and invalid geometry are errors.

### Geometry forms

`circle: {center: [lon, lat], radius_m: N}` is the most direct way to say
"this fishery, roughly here" and is the recommended form for exploration.
`bbox` and GeoJSON `Polygon` cover the rest. For production work a real
boundary — a fishery-right polygon read from GeoJSON or a shapefile — is more
reproducible than typed coordinates; **reading geometry from a file is not yet
implemented** and is the first extension to make.

### The coastline inside the hole

The hole reaches the coast whenever the **transition** does, and that is
normal — the core need not go anywhere near land. Measured on Futtsu: the core
clears the coastline by 714 m while the 2,239 m hole overlaps it by 1,225 m,
picking up 12 coastline edges (4,236 m, lengths 180 / 321 / 597 m). This is
not a reason to refuse anything; the coastline there is simply meshed at the
target size, as any boundary is.

What has to be decided is **where the new boundary nodes go**. Keeping the
existing boundary node list is not an option: a 30 m cell against an untouched
597 m boundary edge would need its apex 172 m away, so it cannot be 30 m, and
the coarse boundary row would break C4 against the fine cells behind it.

| mode | where new nodes go | the sharp corner | fidelity to the real coast |
|---|---|---|---|
| `preserve` | on the existing segments only | kept as it is | unchanged |
| **`resample`** (default) | along the **source** shoreline at the target size | follows the data | **improves** |
| `spline` | on a smooth curve through the base polyline | eased | unchanged |

`preserve` keeps the polyline geometrically identical — subdividing a segment
does not move it — and is the strictest option. `resample` is the default
because it is the only one that *adds* information: the base polyline runs
300–600 m between nodes and sits up to **68.5 m** from the source shoreline at
its segment midpoints (median 5.1 m), and a 30 m resample recovers that for
free exactly where the mesh is being refined anyway. `spline` exists for a
patch with no usable source shoreline; it eases a corner but invents the
easing.

Measured on the Futtsu chain, 10 of its 11 interior vertices turn by 164–180°,
so straight subdivision would be unproblematic along almost all of it; a
single vertex turns by **42.4°**, and that one corner is the whole argument
for `resample` over `preserve` — the source data says whether the coast is
really that sharp or whether 300 m sampling made it so.

Every mode is bounded by `coastline_tolerance_m` (default 100 m) against the
base polyline and verified against it. **No mode touches the coastline outside
the hole**, and the two endpoints where the hole's coastline chain meets the
frozen coastline are fixed, so the polyline joins exactly.

Re-cutting the coastline from the original OSM data — as opposed to resampling
the shoreline polygon the generator was already given — is **not** part of
this operation.

## 5. The patch generator (designed, not yet written)

1. **Select.** Build `core ⊕ transition` in the mesh CRS (EPSG:32654). Mark
   every element whose centroid falls inside it.
2. **Cut.** Delete those elements. What remains is the base mesh with a hole;
   the hole's rim is a closed chain of **existing** nodes.
3. **Fill.** Generate a triangulation inside the hole with
   `h(r) = target + g · (r − r_core)`, clipped above by the ambient field, with
   the rim nodes as fixed points (`pfix`) **and every rim segment as a
   constrained edge (`egfix`)**. `pfix` fixes positions only; the segments
   between them are forced into the triangulation by the CDT, and only `egfix`
   builds those constraints (`mesh_generator.py:1077`). Without them an
   unconstrained Delaunay can bridge a concave rim, and centroid rejection
   does not prove the survivors lie inside the hole.

   The rim does **not** match the ambient size. Measured on the Futtsu
   selection, rim edges run 180 / 467 / 819 m against a 350 m ambient target.
   An 819 m immutable edge cannot carry an approximately equilateral 350 m
   cell: for all angles ≥ 30° its opposite vertex must be ≥ 236 m away and
   another side ≥ 473 m. The sizing field is soft enough to allow this, but
   nothing about the rim is "matched by construction".
4. **Stitch.** Concatenate, deduplicate the rim nodes, renumber.
5. **Verify.** `fmesh-mesh-qa` (21 gates), `frozen_changes()` = 0, the 342
   connectivity comparator, the 364 width gate, and the implied dt against
   `dt_floor_s`.

### Constraints the generator must respect (revision 2)

- The hole must not reach the **open boundary**: OBC nodes are an input.
  Refuse if `core ⊕ transition` intersects the OBC arc.
- With `touch_coast: false` the hole must not reach the coastline either, so
  the land boundary is untouched.
- Every rim node must survive as a distinct vertex. The fork's named error
  (`generate_mesh`, commit `529a462`) fires when the nearest-vertex map is not
  injective after pruning; that is a symptom, and proximity to the target size
  is not the only cause (lost points and topology changes do it too). Requiring
  a small correspondence distance as well is what proves the surviving vertex
  is the right one.
- **The selection is not the analytic footprint.** Whole triangles are taken,
  so the cut reaches past the envelope — measured on Futtsu, 2,615 m against a
  requested 2,239 m, with 22 selected vertices outside the disc. The footprint
  the contract talks about must be the actual set of removed faces, reported
  before filling, not the circle in the recipe.
- **The cut need not be one simply connected hole.** It can have several
  components, pinch at a vertex, or enclose retained islands; a coast-reaching
  cut has an interface chain joined to a physical boundary chain, not a closed
  water-only loop. Each case needs its own handling or an explicit refusal.
- **Seam quality is a whole-mesh property.** A seam vertex's valence counts
  retained neighbours too, so a patch-only valence ≤ 8 is not enough; a corner
  of total angle θ admits k new triangles only if `30k ≤ θ ≤ 130k`. An element
  merely *incident* to an OBC node can change its orthogonality even when every
  OBC coordinate is fixed, so the OBC needs a guard band, not just
  non-intersection.

### Where the code goes

The package is Apache-2.0 and **must not import oceanmesh (GPL)**. The split:

- `src/fvcom_mesh_tools/refine.py` — specification, pre-flight, selection,
  stitching, verification. License-clean.
- `notebooks/4xx_local_refine.py` — the DistMesh call, as notebook 325 already
  does. Notebooks are not part of the distributed package.

## 6. Worked case: Futtsu nori area — REJECTED as declared

`recipes/refine/futtsu_nori.yaml`. The centre came from the data, not by eye:
the Futtsu tidal flat is the largest connected patch of water shallower than
1.5 m (T.P.) off Futtsu — 4.02 km², lon 139.768–139.821, lat 35.304–35.331 —
and the declared centre lies 654 m north of its edge.

| quantity | value | verdict |
|---|---|---|
| centre / radius | (139.7881, 35.3228) / 300 m | |
| target | 30 m | |
| core depth | 3.00 – 4.15 m | |
| dt by shortest edge | 4.70 s | flattering |
| **dt by minimum altitude** | **4.07 s** | **alert: 1.1x the expected steps** |
| target that would hold 4.5 s | 33 m | |
| transition | 1,939 m | |
| core clears the coastline by | 714 m | the fishery itself is offshore |
| transition overlaps the coastline by | 1,225 m | normal; meshed at the target size |
| coastline edges on the rim | 12 (4,236 m) | `coastline: resample` |
| selected elements / rim edges | 173 / 47 | |
| rim edge min / median / max | 180 / 467 / 819 m | not 350 m |
| selection reach vs requested | 2,615 m vs 2,239 m | |

Of 1,200 candidate cells 650–750 m north of the flat with the **core** clear of
land, 112 met a 4.5 s floor **computed from the shortest edge**. Both filters
were wrong: the core is the wrong body to test for land, and the edge is the
wrong length for dt. Under the corrected rules **this recipe is buildable as declared**. The dt
raises an alert, which is the intended behaviour for a given region, and the
transition reaching the coast is ordinary work rather than an obstacle: the
12 coastline edges on the rim are resampled at 30 m along the source
shoreline, and the coastline outside the hole is untouched.

**The site cannot move** — a fishery is given — so if the cost ever has to
come down, the levers are a coarser target (33 m would hold 4.5 s) or a
steeper gradation (shortening the 1,939 m transition, with the C4 cost
measured). Neither is needed to build this one.

## 7. What is still missing before this can be built

The review of revision 1 showed the generator is under-specified in ways that
would surface as failures during a run rather than as refusals before it.

1. **A feasibility and rejection contract.** Success must be conditional and
   the footprint explicit: a set of removable faces with immutable shared
   interface vertices, reported before filling. A failed QA run at the end is
   a diagnostic, not an algorithm; without a declared repair-or-reject policy
   there is no path from "it failed" to "a usable mesh".
2. **Interface segment survival.** `egfix` for every rim segment, verified
   after each cleanup stage — not only immediately after triangulation. No new
   vertex may split a frozen interface segment unless the retained element
   changes too, which enlarges the authorized footprint.
3. **Non-simple cuts.** Several components, pinch points, retained or physical
   islands, and coast-reaching cuts each need explicit handling or refusal.
4. **Seam quality as a whole-mesh property.** Valence and corner-angle budgets
   at seam vertices count retained neighbours; an OBC guard band is needed
   because an element incident to an OBC node can change orthogonality without
   any OBC coordinate moving.
5. **A named dt convention.** The recipe should state the metric and the
   Courant/safety convention. The present figure is an equilateral upper bound
   at Cr = 1 with no velocity allowance: a mesh diagnostic, not an operational
   time step.

## 8. Open questions

1. **Geometry from a file.** Needed before this is useful for real fishery
   boundaries. Shapefile or GeoJSON, with an optional attribute filter.
2. **Several regions at once.** `priority` is specified and validated but the
   generator's overlap handling is not designed. Overlapping cores with
   different targets need a single combined sizing field before cutting.
3. **Bathymetry inside the patch — and the base mesh does not have the
   depths this assumed.** Measured on `sample_repro_final.14`: min 2.00 m
   (the production floor is 3 m), max 735.3 m (the cap is 300 m), 1,160 nodes
   below 3 m, 2,258 edges above r = 0.2, worst 0.907. These are the SRTM15
   depths from notebook 325, not the M7001 production field. Inserting an
   M7001 patch here produces a **mixed-bathymetry** model. Either rebuild the
   baseline's depths once, globally, before freezing it, or declare the mixed
   operation deliberately. The default QA passing 21/21 does not contradict
   this: its depth floor is 2 m and it gates neither the cap nor the r-factor.

   The smoothing itself is a **constrained feasibility problem**, not a
   limiter run on a slice. With `rmax = 0.2`, `|h_i − h_j|/(h_i + h_j) ≤ 0.2`
   is exactly `2/3 ≤ h_i/h_j ≤ 3/2`, so a mutable node next to a frozen depth
   `H` must lie in `[H/1.5, 1.5H]`; a node adjacent to frozen 3 m and 12 m has
   an empty feasible set, and a k-edge path between frozen depths needs
   `H_max/H_min ≤ 1.5^k`. `rfactor_smooth()` moves both endpoints and has no
   fixed mask, so running it on the patch and restoring the rim afterwards
   re-introduces violations (3 and 12 become 6 and 9; restoring 3 gives
   r = 0.5). Five rim nodes in this selection are at 2 m, which conflicts with
   applying a 3 m floor to the rim at all. This needs a constrained solve with
   an explicit feasibility check, not a masked limiter.
4. **Repeated refinement.** Applying a second region to an already-patched
   mesh should work, but the frozen check then compares against the patched
   mesh, not the original. Whether to track a chain of base meshes is open.
5. **Does the fit pass re-run?** `coast_fit` moves boundary nodes; over a
   patch that touches the coast it should run again, but only on the patch.
