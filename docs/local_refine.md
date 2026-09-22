# Local refinement of an existing mesh — design

**Status: implemented, reviewed, and passing end to end.** The Futtsu recipe
builds a 30 m fishery into `sample_repro_final.14` with **QA 21/21**, 4,619
frozen nodes of which **none moved**, no retained face lost, no interface
segment split, and the open boundary untouched — and the generator now
*checks* all of that on the written file and refuses to report success
otherwise. This document is revised as the design moves, and records the
decisions and their reasons so that a later change is made knowingly rather
than by accident.

Revision history is the git history of this file.
`docs/local_refine_review.md` is an adversarial review of revision 1
(gpt-6-astra, 2026-09-22); every blocker it raised was reproduced
independently before revision 2 was written. Revision 4 replaced the design of
the generator with what was built. `docs/local_refine_implementation_review.md`
is an adversarial review of that implementation (gpt-6-astra, 2026-09-22):
sixteen findings and nine executable tests, **all nine of which reproduced**.
Revision 5 is what those fixed; §5.3 lists them.

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

`preserve` keeps the polyline geometrically identical — every original vertex
is kept and only interior points are added, so subdividing a segment cannot
move it — and is the strictest option. `resample` is the default
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
base polyline and verified against it — the Futtsu recipe raises it to 200 m,
because the source shoreline genuinely leaves the base polyline by up to
152 m over this stretch (280 m of a 6.3 km walk, median 7 m) where the base
mesh's 300–600 m spacing chords straight across an inlet. Refusing that at
100 m would only throw the detail away. **No mode touches the coastline outside
the hole**, and the two endpoints where the hole's coastline chain meets the
frozen coastline are fixed, so the polyline joins exactly.

Re-cutting the coastline from the original OSM data — as opposed to resampling
the shoreline polygon the generator was already given — is **not** part of
this operation.

## 5. The patch generator

`src/fvcom_mesh_tools/patch.py` (Apache-clean) and
`notebooks/420_local_refine.py` (the DistMesh call). Run it with

```bash
python notebooks/420_local_refine.py recipes/refine/futtsu_nori.yaml
python notebooks/421_local_refine_map.py outputs/refine_futtsu_nori
```

1. **Measure the ambient.** `ambient_size_field` is the base mesh's own edge
   length at each node, smoothed. Around Futtsu that is 445 m, not the 350 m
   the sizing recipe nominally asks for, and the transition has to reach what
   is actually there.
2. **Pre-flight.** `refine.preflight` on the core: depth, dry fraction, the
   dt alert. Refusal is reserved for the impossible.
3. **Select** (`select_patch`). Elements whose centroid is inside
   core ⊕ transition. The cut is **not** the footprint — whole triangles go,
   so it reaches 3,183 m against a requested 2,812 m, and the report says so.
   A vertex where the removed faces form two fans has four rim edges and no
   closed walk exists, so the cut grows until every rim node has exactly two.
   Refusals: reaching the open boundary or its guard band, emptying the mesh,
   severing the retained mesh.
4. **Rim** (`rim_constraints`). Walk each ring; keep every frozen node; replace
   each maximal run of free nodes with a re-cut coastline stretch. That every
   free rim node is on the physical boundary is structural, not lucky: an
   interface edge is shared with a retained face, so both its endpoints are
   frozen. Returns `pfix`, one `egfix` segment per rim edge, the base node id
   behind each fixed point, and the curve each new coastline point was cut
   from.
5. **Fill** (notebook). `om.generate_mesh(fd, fh, bbox, pfix=, egfix=,
   cleanup="none")` in **mesh CRS metres** — which is also what disables
   generate_mesh's internal tmerc sandwich, since it only engages for a bbox
   that looks like degrees. Then `collapse_thin_triangles` and
   `direct_smoother_lur`, both pfix-protected, by hand.
6. **Stitch** (`stitch_patch`). Retained nodes first and in base order, then
   the new ones; each fixed point matched to its patch vertex **exactly**, so
   a lost one is an error rather than a silent snap; depths from the base
   field (`depths_from_base`) over the **full** base triangulation.
7. **Repair the seam** (`improve_patch`). Flips restricted to patch faces,
   moves restricted to new interior nodes, slides restricted to the coastline
   curve the node was cut from. Every candidate scored and accepted only when
   it strictly improves the worst gate in its neighbourhood.
8. **Verify** (`verify_patch`) and the 21 QA gates.

### 5.1 Six things the implementation found

Each of these was a failure first and a fix second; each is measured.

| what looked right | what happened | the fix |
|---|---|---|
| use oceanmesh's default clean | `make_mesh_boundaries_traversable` takes neither `pfix` nor `egfix` and deleted **121 of 235** constrained rim points | `cleanup="none"`, then the pfix-protected stages by hand |
| also run `bound_connectivity` | its valence flips are blind to the sizing field and coarsened the 28.6 m core to **98.7 m** | not run |
| cut the coastline at `target_h_m` | a 30 m coastline 2.4 km out sat in a field asking for 400 m: min angle **0.8°**, quality **0.002** | cut it at the LOCAL size |
| size the patch as `min(target + g·d, ambient)` | the ramp reached 422 m at a rim whose base mesh is 811 m: C4 **0.649** against a 0.5 gate. Widening until the ramp catches up is not available — at 3.7 km the cut severs the mesh at the Futtsu spit | ramp from the target to the **local** base size over the declared width; the local slope then varies, and `effective_gradation` checks it against the 0.414 that C4 allows (Futtsu: 0.290) |
| the raw per-node ambient field | slope to **0.686**, above what C4 allows between neighbours; 26 % of edges above the recipe's own 0.165 | 20 Jacobi passes: p90 slope 0.243 → 0.094, median size 421 → 439 m |
| repair by Laplacian smoothing | the fill's own smoother has already put every interior node there — the offending node was **0.000 m** from its centroid — so the pass did nothing | search a ring of directions; 15 m off that centroid raised the neighbourhood health 0.664 → 0.705 |

Two more, both about what a check actually covers:

- **C4 is measured across an edge, and the far face need not touch the moved
  node.** Scoring only the moved node's own fan left two failures standing,
  each a patch element beside a retained element twice its area. The
  neighbourhood has to include everything the fan shares an edge with.
- **`preserve` was dropping vertices.** Re-walking a stretch by arc length
  replaced 18 base nodes with 13 and moved the coastline **74.5 m** — the
  opposite of the promise. Subdividing each original segment keeps every
  original vertex, and the departure went to 0.00 m.

### 5.2 What the contract is checked by

`verify_patch` reports, and gates on, seven separate things, because "the mesh
outside is unchanged" is seven claims wearing one coat:

| claim | Futtsu |
|---|---|
| frozen nodes did not move | 4,619 frozen, **0 moved**, max 0.0 m |
| frozen depths did not change | max 0.0 m |
| every retained face survives with the same vertices | 0 missing |
| nothing inverted | 0 |
| **interface segments are still shared, not split** | 0 split |
| the open boundary is the same list in the same order | unchanged |
| no duplicate or orphan nodes | 0 / 0 |

Interface survival is the one a plausible-looking mesh passes without it: a
fill may insert a vertex on a constrained edge, and on the coastline that is
harmless — the chord is where it was — but on an interface segment it leaves
a hanging node while every face and every coordinate still checks out.
`area_change_fraction` is reported rather than gated (+0.001 % here): it moves
legitimately when the coastline is resampled, and it is the one number that
notices a face quietly dropped. Every one of these is checked twice — once in
memory and once by re-reading the written fort.14, because the file is what
the model runs.

### 5.3 What the adversarial review changed

gpt-6-astra reviewed the implementation and supplied nine executable tests.
All nine reproduced before anything was changed; all nine are now regressions
in `tests/test_patch.py` (and one in `tests/test_fort14.py`). The findings
that changed behaviour:

| finding | what was wrong | what it is now |
|---|---|---|
| connectivity | `_require_connected` built a graph over NODES, so two triangles meeting at a single point counted as connected — the opposite of the edge connectivity the docstring promised | a graph over faces joined by shared edges. It immediately caught a real case: a cut at (400, 200) on the test grid leaves one triangle attached by its vertices alone, which the cut now absorbs |
| nested rings | `hole_polygon` unioned every face `polygonize` returned, filling the exclusions: four nested squares plus a disjoint one came back as area 101 where the domain is 57. An island inside the cut would have been handed to the filler as water | rings assembled by nesting parity, each shell carrying the holes whose immediate parent it is |
| a physical island | an island taken whole kept its base node ids, and stitching reads a non-negative base id as "frozen" and looks it up in a map of survivors, so an ordinary circle round a small island could not be stitched at all | those nodes are emitted as new points with their base coordinates |
| frozen means exact | a 5e-7 m change in a frozen coordinate or depth passed, because the check used a 1e-6 m tolerance | exact equality. The tolerance remains where it belongs, on matching a fill vertex to its fixed point |
| a face counted once | retained faces were compared as SETS, so a duplicated element — a second face laid on the first — passed every check | multisets, plus a gate on any edge carrying three faces |
| the written file | nothing re-read it, and `write_fort14` wrote depths at `.10e`, eleven significant figures: an interpolated 5.12345678912345 m came back 2.3e-11 m different | depths at `.17g`, and the frozen zone is verified again on the file that was written |
| flips | incidence and adjacency were computed once per sweep and only the two rewritten faces marked stale, so a later flip scored against a mesh that no longer existed; and valence was counted over the candidate's own fan, reading [5,5,7,7] where the mesh had [5,5,7,9] | both rebuilt after every accepted flip; valence tracked over the whole mesh. It is a cost rather than a veto, because forbidding every intermediate excess also blocks the sequences that end below the limit — and C5 is gated on the finished mesh |
| a zero-width transition | ambient already at the target gives width 0, and `d / 0` made the whole sizing field NaN — which DistMesh accepts | defined: target inside, base outside |
| a bbox | the driver guessed "circle" from the spread of vertex radii, and a square's corners are all equidistant from its centre, so a near-square bbox became a 257-point disc | the recipe's declared kind is carried on `RefineRegion`; polygon interiors are projected too |
| depths | `improve_patch` moved nodes after `stitch_patch` had evaluated the base field, and nobody re-evaluated it: a node moved from (0.2, 0.3) to (1, 1) on a 5 + x field kept 5.2 where 6.0 is right | `refresh_depths` after the repair |
| one shoreline for every stretch | the driver chose one ring by closest-pair distance to the whole free rim, contradicting the docstring that said the caller picks per stretch | every ring within 3 km is offered and the nearest to EACH stretch is used |
| QA | the driver wrote the mesh and stopped. That the worked case passed 21/21 was a fact about someone running QA afterwards, not about the generator | the 21 gates run on the written file and a failure ends the run |

And one the review noted in passing that was simply wrong in this document: a
30-30-120 cell's minimum altitude is **1/sqrt(3)** of an equilateral cell's
with the same side, not a half.

Two of the review's points stand unfixed and are listed in §7: there is no
caller-declared maximum for how far the cut may grow, and `effective_gradation`
measures the ramp term only, not the ambient term, so it is a diagnostic
rather than a C4 certificate.

### 5.4 The seed search, and why there is one

`improve_patch` is greedy: it accepts only moves that do not worsen the worst
gate in a neighbourhood, so it stops at a local optimum, and which one depends
on where the fill started. That is not a small effect. On this patch,
neighbouring configurations finished at 22.96°, 25.6°, 27.5° and 30.01°
against a 30° gate.

So the driver **searches**: for each seed it fills, repairs, stitches,
verifies, writes and runs the 21 gates, and it stops at the first mesh that
passes. Which seed was accepted is in the report. If none passes, the run
fails with the best attempt's failures named, and the right response is to
widen the transition or coarsen the target — not to lower the gate.

Two further things were needed before any seed passed:

- **The strict pass runs first, the soft pass only if it stalls.** The soft
  rule accepts a move that holds the worst margin and improves the saturated
  rest. It rescued a run stuck at 22.96°, and it cost a run that the strict
  rule took to 30.01°: wandering laterally changes which basin you end in.
- **A flip may never breach the valence gate, not even in passing.**
  Allowing an intermediate excess does buy reach, and it leaves nodes at 9
  that no later flip can bring down: with it, `preserve` failed C5 at every
  seed. A guarantee beats a heuristic. DistMesh can still hand over a
  valence-9 node of its own, and that is what the seed search is for.
- **The source shoreline is simplified to a quarter of the local element
  size before it is walked.** A coastline digitised at metres cannot be
  represented by 400 m elements; walking it at 400 m gives chords that turn
  sharply against each other, and every surviving QA failure was a coastal
  element 1.8–2.3 km out, where the transition is coarse and the coast is
  not. This is not detail thrown away — there is no room for it at that size
  — and what remains is still bounded by `coastline_tolerance_m`.

### Where the code goes

The package is Apache-2.0 and **must not import oceanmesh (GPL)**. The split:

- `src/fvcom_mesh_tools/refine.py` — the recipe, and the pre-flight refusals.
- `src/fvcom_mesh_tools/patch.py` — selection, rings, the coastline modes,
  sizing, stitching, seam repair, verification. License-clean, 47 tests.
- `notebooks/420_local_refine.py` — the DistMesh call, as notebook 325 does.
- `notebooks/421_local_refine_map.py` — the figures.

## 6. Worked case: Futtsu nori area — built

`recipes/refine/futtsu_nori.yaml`. The centre came from the data, not by eye:
the Futtsu tidal flat is the largest connected patch of water shallower than
1.5 m (T.P.) off Futtsu — 4.02 km², lon 139.768–139.821, lat 35.304–35.331 —
and the declared centre lies 654 m north of its edge.

| quantity | value |
|---|---|
| centre / radius | (139.7881, 35.3228) / 300 m |
| target | 30 m, **achieved 30.2 m** (median of 1,069 core edges; p90 31.2 m) |
| measured ambient around the site | 445 m |
| transition | 2,512 m |
| effective gradation | 0.290, against the 0.414 C4 allows |
| core depth (from the base mesh) | 2.5 – 6.16 m |
| **dt by minimum altitude** | preflight 3.34 s, **achieved 2.61 s** (base 11.86 s) |
| elements | 8,252 → 10,436 (+2,184); nodes 4,734 → 5,826 |
| elements removed / retained | 244 / 8,008 |
| rim: interface / coastline edges | 33 / 19; lengths 180 / 445 / 886 m |
| selection reach vs requested | 3,183 m vs 2,812 m |
| coastline nodes replaced / new | 18 / 12, on the source shoreline |
| coastline departure from that source | nodes 64 m, chords 118 m (tolerance 200 m) |
| **QA** | **21/21**, angles 30.01–119.26° — the base mesh's own range |
| frozen zone | 4,619 nodes, 0 moved, 0 depth change, 0 faces lost, 0 splits |
| water area | +0.001 % (the resampled coastline, not a lost or gained face) |
| DistMesh seed | **2**; seeds 0 and 1 failed 3 and 2 gates and were rejected. `preserve` passes at seed 0 |

`coastline: preserve` reaches 21/21 as well, with the coastline geometrically
identical. The two differ in what they buy: `resample` follows the source
shoreline up to 143 m closer than the base polyline, which is the point of
it, and pays for it with a harder mesh — it is the mode that needed the seed
search and the shoreline simplification.

**The cost is the time step, and it is real.** The base mesh allows 11.86 s by
minimum altitude; the patched mesh allows 2.61 s, and the binding element is
inside the core. The pre-flight alert predicted 3.34 s
from an equilateral cell over the deepest core sample; the achieved 2.61 s is
lower because real cells are not equilateral — which is exactly what the
alert says, though the factor it quotes for a 30-30-120 cell should be
1/sqrt(3), not a half. **The site cannot move** — a fishery is given — so
if the cost has to come down, the levers are a coarser target (40 m would hold
4.5 s) or a steeper gradation. Neither was needed to build this one.

## 7. What the review asked for, and where it stands

The review of revision 1 listed five things that would otherwise surface as
failures during a run rather than as refusals before it.

| asked for | state |
|---|---|
| **a feasibility and rejection contract** | done. `select_patch` refuses a cut that reaches the open boundary or its guard band, empties the mesh, severs it (by face adjacency), or pinches irreparably, and reports the removable face set and its real reach before anything is meshed. A fill that comes back poor is now rejected rather than written: the driver runs the 21 gates on the written file, tries the next seed, and fails the run when none passes |
| **interface segment survival** | done and gated. `verify_patch` checks every interface segment is an edge shared by two faces in the finished mesh, after every cleanup stage, not just after triangulation |
| **non-simple cuts** | done for the cases that arise. Several components, nesting and pinch points are handled: rings are walked from the rim graph, nesting parity decides both `ring_is_hole` and the polygon handed to the filler, a pinch grows the cut until the rim is a manifold, and a spike joined only at its vertices is absorbed. A physical island taken whole inside the cut now stitches, with a test |
| **seam quality as a whole-mesh property** | done. The repair scores C1, C2 and C4 over the fan **and everything it shares an edge with**, so a patch element is measured against its retained neighbour; valence counts incident faces exactly as the C5 gate does; the OBC has a guard band (500 m by default), not just non-intersection |
| **a named dt convention** | partly. `dt_s` is still an equilateral upper bound at Cr = 1 with no velocity allowance, but the alert now says so, the 30-30-120 factor is stated correctly as 1/sqrt(3), and the driver measures the **achieved** minimum-altitude step on the finished mesh and reports both: predicted 3.34 s, achieved 2.61 s |

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
   With `coastline: resample` the question largely dissolves — the new nodes
   are placed on the source shoreline to begin with — but the frozen
   coastline outside the hole keeps whatever offset the base fit left it.
6. **A global-rebuild option.** Asked for alongside the exterior-preserving
   patch and not yet written. It is the `sizing.py` route — add the region to
   the sizing recipe and rebuild — and its cost is precisely the thing this
   operation avoids: the ports and channels elsewhere come out subtly
   different, so a comparison against the base run cannot separate the
   fishery from everything else.
7. **The repair is greedy, and the driver searches around it.** See §5.4.
   Seeds 0 and 1 of the Futtsu recipe fail the gates and seed 2 passes; the
   run reports which it used. The search makes success checkable, not
   certain: when no seed passes, the honest response is to widen the
   transition or coarsen the target, never to lower the gate.
8. **The cut has no declared maximum.** `select_patch` reports how far it
   grew and how far it reached, and nothing compares that against a limit
   the caller stated, because there is no recipe key for one. The OBC guard
   is also node-to-node, so a long OBC segment can pass closer to the cut
   than the guard suggests (review finding 8).
9. **`effective_gradation` is a diagnostic, not a certificate.** The field is
   `target + (B(x) - target) * d(x) / W`, whose gradient has an ambient term
   that the reported number omits, and the `1 - 1/(1+g)^2` threshold assumes
   similar neighbouring triangles where C4 is an area ratio across an edge
   (review finding 13). C4 itself is gated on the finished mesh, which is
   where it belongs.
