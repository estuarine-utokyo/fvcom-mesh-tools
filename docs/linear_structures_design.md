# Breakwaters, piers and other linear structures: a design

Status: implemented and run in FVCOM on the Kimitsu port (§7-§17). §1-§5 are
the original design, §6 and §15.1 the owner's decisions; where they differ,
§15.1 (2026-09-24) supersedes §6 item 5.

## 1. What the owner saw, and why it happened

Inside the 30 m core of `recipes/refine/kimitsu_port_hires.yaml` there is a
small harbour -- an L-shaped breakwater, several piers, a basin -- that OSM
draws in detail and the delivered mesh ignores: the triangles run straight
across it (`outputs/refine_kimitsu_port_hires/structures_removed.png`).

The cause is `patch.filter_shoreline`, which is doing exactly what it was
written to do. It removes land narrower than `3 * h0` (90 m at h0 = 30 m) on the
grounds that an `h0` mesh cannot carry it **as an area**, and it closes water
narrower than the same width. Measured in the 900 m region (notebook 425, job
115443):

| filter | land removed | of which LINEAR (length > 5 x width, > h0) | water closed |
|---|---:|---:|---:|
| 1 x h0 (< 30 m) | 14,474 m2, 29 pieces | 5 | 801 m2 |
| 2 x h0 (< 60 m) | 25,150 m2, 26 pieces | 16 | 5,530 m2 |
| **3 x h0 (< 90 m), what ran** | **46,156 m2, 26 pieces** | **15** | **37,193 m2** |

The linear pieces are the structures: a pier 852 m long and 12.4 m wide, a
breakwater arm 359 m long and 60 m wide, and a dozen piers 40-80 m long and
2-9 m wide. The water it closed includes the harbour basin itself --
11,506 m2, 51.6 m wide -- and a 19,275 m2 basin 56.9 m wide.

The judgement was the wrong one for these features. "Too narrow to be an
area" was read as "too narrow to exist", and a breakwater is exactly the thing
that is too narrow to be an area and still has to be there. **A 12 m pier
cannot be meshed as land at 30 m, but it can be meshed as a line.**

## 2. What FVCOM needs for a wall -- an edge is not enough

The owner's proposal is to lay mesh edges along the structure. That is the
right geometry and not sufficient physics, and the reason is in the scheme:

* **Elevation and scalars live on nodes, in median-dual control volumes**
  (`tge.F`, the `XIJE/NIEC` construction: each face runs from an element
  centroid to an edge midpoint). The control volume of a node ON the line
  contains parts of the elements on BOTH sides of it. Water crossing the line
  at that node never crosses a control-volume face. An interior edge is not a
  face of anything a scalar sees.
* **Momentum lives on elements**, and the shared edge IS the face between the
  two elements on either side. Flux crosses it like any other.

So a line of edges does not block flow. What does is making the line part of
the **domain boundary**: the nodes along it are duplicated, the elements on one
side use one copy and those on the other side use the other, and the two sides
no longer share an edge. FVCOM then finds the wall itself, topologically --
`NBE = 0` on both sides, so `ISBCE = 1` -- exactly as it finds the coastline.

This is not a hypothesis. It is how the owner's own FVCOM tree already
describes `THIN_DAM` (`FVCOM/docs/sigma-z-dike-groin-thin-dam.md` §1): *"splits
the mesh (nodes duplicated, `N_DAM_MATCH`/`E_DAM_MATCH`) ... both sides own it
as a domain boundary (`NBE = 0` -> `ISBCE = 1`)"*, and `THIN_DAM`'s own
contribution is the re-stitching ABOVE a crest, which is the overtopping the
owner says a breakwater does not need.

**The consequence is the useful part:** a split mesh **without** `THIN_DAM` is
an impermeable, full-height wall, and it needs **no FVCOM change and no
rebuild**. `THIN_DAM` stays available for a structure that must overtop; it is
frozen in the owner's tree (same document) and is not part of this design.

## 3. Proposed design

### 3.1 Classify what the filter would have removed

Keep the width judgement, and split its output three ways instead of one:

| class | test | treatment |
|---|---|---|
| **area** | local width >= `k_area * h0` | land, as now -- the coastline goes round it |
| **wall** | local width < `k_area * h0` and length >= `L_min` | collapsed to its centreline, meshed as a split line |
| **drop** | length < `L_min` | removed, and listed in the report |

Width is measured **locally along the feature**, not per polygon, because real
structures are mixed: the L-shaped breakwater above has a 60 m body and a
10 m arm, and those are different classes.

The centreline comes from the medial axis of the thin part -- Voronoi vertices
of the densified outline that fall inside it, joined into a graph and pruned
of short spurs. That is numpy/scipy/shapely only (Apache-2.0). Where OSM maps a
structure as a LINE (`man_made=breakwater`, `man_made=pier` ways), the line is
the centreline and no medial axis is needed -- but those ways are **not in the
data on disk** (the Geofabrik free extract has no `man_made` layer); see §5.

### 3.2 Attach the wall to the coast

A pier's root touches land: its centreline is extended or trimmed to end
exactly on the coastline, and that point becomes a coastline node. The other
end is a **free tip**. A detached breakwater has two free tips.

### 3.3 Mesh the wall as a constraint

The centreline is resampled at the local size and handed to oceanmesh as
constrained edges -- `pfix` points and `egfix` segments in the INTERIOR of the
domain. oceanmesh already forces `egfix` through every retriangulation with
the CGAL constrained Delaunay binding (`mesh_generator.py:1077`); today the
driver uses it only for the rim. The mesher sees an ordinary interior line; it
knows nothing about walls.

### 3.4 Split the mesh along it

After meshing, a new function (Apache side, `patch.py` or its own module) cuts
the mesh along each wall:

* a node **inside** a wall: its fan of elements is divided by the two wall
  edges into two sectors, and one sector gets a new node at the same
  coordinates;
* a node where a wall **meets the coast**: already a boundary node; its fan
  is divided the same way, so the coastline continues on both sides of the
  pier;
* a **free tip**: its fan is one sector that wraps round the end of the wall,
  so it is **not** duplicated -- the boundary turns through 360 deg there;
* a node where walls **cross or branch**: one copy per sector.

The copies carry the same depth. They are recorded as pairs in the report
(and could be written in `THIN_DAM`'s `_dam_node.dat` format later, if a
structure ever needs to overtop).

### 3.5 What changes around it

* **QA**: the duplicate-node check (`duplicate_tol_m`) must accept the
  declared copies and still catch any other coincidence. The boundary checks
  must accept a zero-area island (a detached breakwater's ring).
* **The frozen contract**: unchanged. Walls exist only inside the hole, so
  every copy is a new node and every split element a patch element.
* **Bathymetry**: the ladder samples a copy once; both copies take the value.
  The r-factor limiter needs nothing -- there is no edge between copies.
* **Boundary lists** for fort.14: each wall contributes to the ring of the
  land it is attached to, or forms its own island ring when detached.

## 4. What must be tested before it is believed

Two things are not known and are the risk of the design:

1. **The free tip.** FVCOM's boundary treatment at a node whose boundary turns
   through 360 deg has not been exercised here. A cusp could give a degenerate
   boundary normal.
2. **Coincident nodes.** FVCOM builds its neighbour tables from connectivity,
   which the split handles, but anything in the tree that searches by
   coordinate would see two nodes at one point.

So the first deliverable is not the harbour. It is a **synthetic FVCOM test**:
a rectangular channel with a pier half-way across it and a detached
breakwater, tide at one end, dye released on one side. It passes if (a) the
run is stable, (b) no dye reaches the far side of the breakwater except round
its ends, and (c) the velocity normal to each wall is zero to round-off at
every output. Only after that, the Kimitsu harbour, finished and run against
the base as in §12 of `refine_coast_and_bathy_design.md`.

## 5. Decisions for the owner

1. **Representation.** Split mesh (recommended; exact geometry, no FVCOM
   change), or a thin land hole widened to a meshable width (no duplicate
   nodes, but the structure grows and its tip produces short edges), or
   `THIN_DAM` (needs a rebuild; frozen in the FVCOM tree; only worth it for
   overtopping).
2. **Thresholds.** `k_area` (proposed 2, i.e. an area feature needs two
   elements across; the filter used 3) and `L_min` (proposed `h0`: a
   structure shorter than one element is dropped and reported).
3. **Water.** The same filter closed a 51.6 m harbour basin at 3 x h0. With
   the structures kept as walls, should water be closed only below `2 * h0`
   (the basin would survive with about two elements across it), or not closed
   at all inside a harbour?
4. **Source.** The thin parts of the OSM land polygons (on disk now), and
   optionally OSM's `man_made=breakwater/pier` LINES, which would have to be
   extracted on a login node (ODbL, as the rest of OSM). Some breakwaters are
   mapped only as lines and are absent from the land polygons entirely.
5. **Volume.** A zero-width wall counts the structure's footprint as water on
   both sides: the 852 m x 12 m pier is about 10,000 m2. Acceptable, or should
   walls wider than some limit be kept as areas even below `k_area * h0`?

## 6. Decisions (owner, 2026-09-23)

| # | decision |
|---|---|
| 1 | **split mesh** |
| 2 | `k_area = 2`, `L_min = h0` |
| 3 | water is closed below **2 x h0** (it was 3) |
| 4 | OSM `man_made` lines: **not now**; later |
| 5 | every wall is **zero-width**, and the footprint it hands to the water is reported (superseded 2026-09-24, §15.1: land down to half an element stays land) |

On 5, what was being asked, said plainly: a real pier 852 m long and 12 m
wide covers about 10,000 m2 of sea. A zero-width wall lets the elements on
both sides reach its centreline, so that footprint becomes water in the model
-- at 5 m depth about 50,000 m3, a harbour's tidal prism enlarged by the
fraction of it the structures occupy. The alternative was to keep structures
above some width (say 20 m) as real land holes even below `2 * h0`, which is
geometrically exact and costs elements smaller than the target beside them,
and so the time step. The owner chose zero width throughout; the lost
footprint goes in the report.

Order of work, as §4 requires: the split function and its unit tests, then a
synthetic FVCOM case (a channel with a pier and a detached breakwater), then
the harbour.

## 7. The FVCOM test: a split wall stops water, an edge does not

Job 115449, notebook 427. The finished Kimitsu port case, which had already
integrated for 20 days, with ONE change: three paths of existing interior
edges split into walls -- a seal from coast to coast (24 edges), a pier
(8 edges, one free tip) and a detached breakwater (13 edges, two free tips).
45 nodes duplicated, 3 free tips, the domain cut in two. OBC ids, sponge,
tide, namelist and element order all kept; 64 ranks, 2 days. The control is
the same case unsplit, already run.

| | split | unsplit (same edges) |
|---|---:|---:|
| finished | exit 0, 97 records, TADA | exit 0, 97 records |
| sealed region (383 nodes): max \|zeta\| | **0.0 m** | 0.419 m |
| sealed region: M2 amplitude | **0.0 m** | 0.395 m |
| open water: M2 amplitude (median) | 0.3976 m | 0.3974 m |
| max speed at elements touching a free tip | 0.194 m/s | -- |
| 99th percentile speed, whole mesh | 0.245 m/s | 0.245 m/s |

* **A split wall is impermeable, to the last digit.** Behind the seal the
  level never left zero.
* **The same edges unsplit let the tide straight through**: 0.395 m in a
  region whose open-water neighbour has 0.397 m. §2's claim, that an edge is
  not a wall in FVCOM, is now measured in FVCOM rather than argued from
  `tge.F`.
* **Both risks of §4 are retired.** Coincident nodes did not trouble the
  build (METIS, neighbour tables, 64 ranks), and the free tips are not a hot
  spot -- the fastest element touching one runs below the mesh's own 99th
  percentile.
* The walls change the open-water tide by 0.19 mm.

Three things went wrong on the way, all in the test harness and none in the
wall: the synthetic channel of notebook 426 had a corner element with one
open and one solid side, which `tge.F` refuses; it then went unstable beside
its open boundary after 1.8 days with AND without walls, so it could not
separate a wall from a forcing problem and was set aside; and FVCOM reads
`INPUT_DIR` into 80 characters and truncates it silently.

## 8. The harbour, with walls (job 115511)

`recipes/refine/kimitsu_port_hires.yaml`, filtered at 2 x h0, walls extracted,
meshed as interior constraints and split. 24 wall pieces in the hole, 23 ends
rooted on the coast; 49 wall edges split, 47 nodes duplicated, 15 free tips,
one connected domain. Achieved 28.5 m median against 30 m, 99.9 % of the
water within 1.25x; frozen zone exact; coastline 0.3 m median from OSM.

QA 17/21 with 32 violations introduced, down from 78 at the first complete
run. The gates each earned:

| rule | why (measured) |
|---|---|
| walls kept 0.4 h clear of the coast except at roots | a wall hugging a quay made 3.1 / 4.8 deg elements and closed pockets off |
| a root is inserted INTO the coastline rim, judged at its foot | a root stopping short of the coast leaves a gap; judging at the wall end inserted a vertex twice |
| wall points within 0.25 h of another fixed point are that point | the mesher merges them, and the edge between becomes a self-edge |
| a detached single edge gets a midpoint | two free tips and nothing between: no second sector to split into |
| the shortest wall edge bordering water cut off from the OBC is withdrawn | withdrawing whole pieces removed an arm of the L-shaped breakwater |
| no wall meets another line at under 30 deg | elements between two lines cannot be wider than their angle |
| a free tip within 0.5 h of another line loses its last edge | dropping the acute crossing left a tip 4.4 m from the other wall |

**What remains is in the transition, not the harbour.** The worst element,
4.7 deg, is on the resolved COASTLINE 2.3 km east of the region, touching no
wall; the next wall-related ones (22-30 deg) are on the long pier where it runs
north into the transition. Both come from judging features at the declared
h0 across the whole hole (owner, 2026-09-23): features and walls that a 30 m
mesh can carry survive where the transition's elements are 150-400 m. The
earlier 3 x h0 port mesh had no introduced violations; the difference is
exactly the detail 2 x h0 lets through where elements are coarse.

## 9. Judged at the local size (owner, 2026-09-23; job 115531)

`patch.filter_shoreline_local`: the size field sampled over the patch
footprint, cut into octave bands, each band given the land filtered at its
lower bound, joined and cleaned once at h0. Walls shorter than the element at
their middle are dropped. Bands on this patch: 30, 60, 120 m (8.44 km2), 240 m
(34.71 km2), 480 m (44.03 km2), 960 m (0.13 km2).

| | one h0 everywhere (115511) | local size (115531) |
|---|---:|---:|
| minimum angle | 4.70 deg | **22.87 deg** |
| QA | 17/21 | **18/21** (C2 now passes) |
| violations introduced | 32 | **17** |
| walls split / free tips | 49 / 15 | 47 / 14 |
| achieved in the region | 28.5 m, 99.9 % | 28.5 m, 99.9 % |

The region is unchanged, as the rule promised; the transition lost the detail
its elements could not carry, and with it the 4.7 deg coastline element and
the walls running into coarse ground. What is left -- 14 elements between
22.9 and 29.0 deg, one of them the base's own element 2101 -- sits at wall
roots and tips in the harbour, where the repair may not move a wall node.

## 10. Wall nodes slide along their wall (job 115553)

Each wall piece's centreline joined the curves the repair may slide along. The
repair already confines a node to the span between two vertices and never
moves a vertex, so ends, corners and junctions stay put; a root has three
boundary neighbours and is not slid. The two sides of a wall are separate
boundaries, so a copy sliding on one side leaves the slit a slit -- which is
also why the resolution measure had to stop using matplotlib's trapezoid map:
two collinear boundary polylines whose vertices no longer match are refused
by it, as coincident nodes were before. It now asks an STRtree of the
triangles.

| | pinned (115531) | sliding (115553) |
|---|---:|---:|
| minimum angle | 22.87 deg | 22.85 deg |
| C1 violations | 14 | **8** (one is the base's element 2101) |
| C4 violations | 4 | **2** |
| violations introduced | 17 | **9** |
| achieved in the region | 28.5 m, 99.9 % | 28.6 m, 99.9 % |
| coastline from OSM | -- | median 0.2 m, max 37.4 m |

The seven patch C1 elements (22.9-28.8 deg) sit at wall roots and at two
spots on the resolved coastline; none is below 22.8 deg.

## 11. Roots that do not cut the coast short (job 115569)

Zooming on the remaining elements (notebook 430, the fill's saved constraints)
put most of them beside short rim edges: inserting a root split 30 m coastline
edges into 9.1-16.5 m pieces. A foot within half an element of either end of
its segment now takes that vertex (the root moves along the coast by at most
h/2), each wall is walked from its root rather than re-ended at it, and wall
edges under half an element are merged where the wall runs straight.

| | 115553 | 115569 |
|---|---:|---:|
| minimum angle | 22.85 deg | **24.84 deg** |
| C1 violations (incl. base element 2101) | 8 | **5** |
| C4 violations | 2 | **1** |
| violations introduced | 9 | **5** |
| water within 1.25x of target | 99.9 % | **100.0 %** |

What is left: two elements at wall roots in the harbour (24.8 and 26.8 deg),
where a short edge of the resolved coastline itself or a corner of the wall
sits next to the root, and two at a wall in the coarse transition south of
the region (25.3 and 28.8 deg).

## 12. Simplified centrelines and blunted roots (job 115573)

The zoom (notebook 430) on the four elements left by §11 found two causes.

- **A short rim edge at a root.** A wall continuing a pier that the width
  filter had cut short rooted at the cut end, which was 13.6 m wide against
  30 m elements (24.8 deg). A root now absorbs a neighbouring rim point
  closer than half an element: the two become one point at their midpoint.
  Only points the fill added may move, never a frozen one, and the hole is
  rebuilt afterwards.
- **Bends the size cannot resolve.** The medial axis keeps every corner, and
  `_subdivide` keeps every vertex it is given. So a 5.6 m bend left a 22 m
  edge beside 45 m ones at a root (26.8 deg), and a bend on the transition
  wall put a 123 deg element against a 160 m wall edge. Each piece is now
  simplified (Douglas-Peucker) at 0.2 of the local element before it is
  walked. Junctions between walls are kept whatever their bend.

| | 115569 | 115573 |
|---|---:|---:|
| minimum angle | 24.84 deg | **28.81 deg** |
| C1 violations (incl. base element 2101) | 5 | **2** |
| C4 violations | 1 | **0** |
| violations introduced | 5 | **1** |
| QA | 18/21 | **19/21** |
| water within 1.25x of target | 100.0 % | 100.0 % |

The five seeds gave 5, 3, 2, 1 and 2 introduced violations, and seed 3 was
accepted. The one left is at (393366, 3906973), 28.8 deg, where the
transition wall has 100-160 m edges. The other failing check, min_depth_clip,
is reported and not gated on this branch.

## 13. No node in a single element; the walled harbour in FVCOM (jobs 115612-115618)

The first 20-day M2 run on the walled mesh (115589) finished cleanly, but
three nodes kept an M2 amplitude of **exactly 0**, while every neighbour had
0.44 m. Each of them was in ONE element with two boundary sides: inside a
wall bend, or at a 44 deg corner of the resolved coastline. FVCOM never
updates such a node. Neither the base nor the unwalled port mesh had one,
and no QA gate looked for it.

- **QA:** a new FVCOM gate, `no_lone_corner_nodes`, fails any node in a
  single element that is not on the open boundary.
- **The geometric cause.** Once a wall is split, the water between two
  boundary lines that meet at under 60 deg holds one element if every angle
  is to stay at 30 or more. Bisecting that element leaves angles under 30:
  the patch went from 1 violation to 8-10. So:
  - a wall edge meeting another line at under 60 deg is dropped (was 30);
  - `patch.blunt_acute_corners` cuts off a resolved-coastline corner under
    60 deg on the water side. It walks one local element along each side
    and joins the two points reached with a chord; a first version cut at
    0.45 of the corner's own edges, left a 7.4 m chord and cut the external
    step from 2.46 to 1.37 s.
- **What is left:** a lone node the mesher still makes is opened by
  `walls.open_lone_corners`, before the repair so the new node is smoothed.
  `improve_patch` never flips an edge that would leave a node in one element.

| | 115573 | 115612 |
|---|---:|---:|
| violations introduced by the patch | 1 | **0** |
| QA | 19/21 | **20/22** (base element 2101 at 29.0 deg; min_depth reported only) |
| nodes in a single element | 3 | **0** |
| finished-mesh external step | 1.2 s | 1.0 s |

M2, 20 days, the walled mesh against the base (same forcing, DTE 1.0 s):

- **Both runs:** complete (TADA) and finite.
- **At the five gauges:** walled minus base is -0.2 to -0.4 mm in amplitude
  and +0.016 deg in phase.
- **In the harbour** (notebook 431, `outputs/m2wall_20260924_025411/`):
  - the amplitude is 0.438-0.443 m everywhere, with no frozen node;
  - the long breakwater separates a sheltered basin, whose amplitude is
    0.443 m against 0.439-0.441 m outside;
  - currents stop at the walls.

## 14. Drawn over the raw OSM: four defects (commit 289de6f)

Drawing the delivered coast and walls over the raw OSM land
(`notebooks/432_walls_vs_osm.py`) showed four things that the QA gates cannot
see, because each of them is a valid mesh:

| defect | cause | fix |
|---|---|---|
| a quay block standing in the harbour meshed as **water** | the rim re-draws the BASE coastline along the source; land the base never had -- here left detached when the filter made the pier joining it to the shore a wall -- has no stretch to be re-drawn from | `patch.island_rings`: every filtered land polygon inside the hole, with half an element to spare from the rim, becomes an island of the rim, resampled at the local size with its corners kept |
| a 15 m gap between the two arms of an L-shaped pier | a wall tip within half an element of another line lost its last edge | the tip **joins** the line: a vertex within half an element takes it, or the wall edge is split at the foot; the angles this makes go back through the 60 deg rule (§13) |
| a 43 m pier off the quay block dropped | the piece was judged on its length after the 0.4-element coast clearance (29 m < h0), although its end is rooted and regains the clearance | judged with the clearance added back for each end that will be rooted |
| a wall whose two sides no longer met | the repair slid the two copies of a wall node independently, up to 1.3 m apart across a bent wall | `walls.rejoin_copies` puts each group back at one position (the one with the best worst angle, never inverting) |

`notebooks/433_before_after.py` draws two refinements side by side over the
raw OSM, with the wall centrelines extracted from OSM and the pieces kept in
the hole, so a missing wall can be traced to the stage that lost it.

## 15. Narrow piers stay land; walls follow their pier (commit af263d2)

### 15.1 The owner's rules (2026-09-24)

The owner pointed out that the mesh size is a statement about the **water**:
a pier narrower than an element can still be meshed round as land, as long
as its length is resolvable, and making every pier narrower than two
elements a wall was not right. The owner also saw walls running at an angle
to piers that stand square to the quay.

| # | decision |
|---|---|
| 1 | land down to **0.5 x h** stays land (its outline is coastline); thinner land becomes a wall; shorter than h is dropped, as before |
| 2 | an end 0.5-0.75 x h wide closes onto **one point** |
| 3 | water narrower than **2 x h** is still closed |
| - | (after §15.4) the 0.5 x h land rule applies where the elements are up to 2 x h0; the coarse transition keeps 2 x h |

Why 0.5 x h: the only edge the pier's width sets is the one across its end,
and a triangle with a 0.5 h base and two sides of h has a 29 deg apex -- the
C1 limit. Below 0.75 x h the end therefore closes onto its midpoint.

### 15.2 What was built

- `filter_shoreline(..., land_width_factor=)` opens the land at its own
  threshold, separate from the water closing; `filter_shoreline_local`
  applies it in bands up to `land_width_max_band`.
- `patch._corner_walk` replaces the arc-length walk on the hires branch.
  The walk placed stations by arc length and cut every corner it passed,
  which on a pier under two elements wide cuts across the pier. Now:
  - every vertex of the simplified source is kept;
  - an edge under half an element on a straight run loses that vertex;
  - a STEP between two corners collapses to its midpoint (dropping one of
    the corners had cut a 14 m step in a quay into a 165 m chord across the
    water);
  - an end under 0.75 x h closes onto its midpoint **from half an element
    back**, so the pier's sides stay parallel (joining the root corners
    straight to the midpoint had made a 20 x 70 m pier a triangle).
- `walls.straighten_walls`: a wall whose middle (half an element or a
  quarter of its length off each end) lies within 0.1 h of a line becomes
  that line. The medial axis bends at a pier's root, where it branches into
  the corners of the junction with the quay.
- Rooting: a wall is rooted where it, **carried on straight, meets the
  coast**, not at the nearest point of the coast; and on a straight coast
  (bend < 15 deg) the coast vertex moves to the root rather than the root to
  the vertex. Taking the vertex had moved the wall sideways by up to half an
  element.
- A wall node with a copy no longer slides in the repair (§10): with straight
  walls the copies drifted up to 37 m apart before `rejoin_copies` pulled
  them together, bending the elements on one side. A free tip, which has no
  copy, still slides.

### 15.3 On the Kimitsu harbour

The structures the width filter took from the land, with each wall's angle to
its OSM footprint's long axis (`notebooks/433`):

| structure | length | width | before (§13) | after |
|---|---:|---:|---|---|
| south-west pier | 67 m | 30.4 m | wall, 2.3 deg | **land** |
| south pier | 72 m | 20.4 m | wall, 11.1 deg | **land**, end closed onto its midpoint |
| L-shaped pier | 245 m | 25.9 m | wall, 30.3 deg | **land** |
| quay extension | 113 m | 14.8 m | nothing (coast-hugging) | **land** |
| long breakwater | 687 m | 13.4 m | wall, 4.8 deg | land where 15 m or wider; its 2.2 m end a wall, 7.9 deg |
| thin piers | 34-35 m | 2.9-4.1 m | walls, 1.1-13 deg | walls, 0.1 and 2.1 deg |

In the run:

- **Straightening:** five walls were straightened, turning their
  end-to-end direction by 0.4-23.9 deg.
- **Rooting:** seven coast points moved to a root, and one tip joined a line.
- **Water closed:** the slips between the new land piers, all narrower than
  60 m, were closed, per rule 3.
- **Island:** the quay block of §14 is joined to the mainland now that the
  L-shaped pier is land, so no island was added.

### 15.4 Two things the new rules exposed

1. **The base's frozen coast.** The water a patch meshes is bounded by OSM
   land and by the frozen mesh's own coast, and the closing rule saw only the
   first. The base's land beside its frozen coast -- each retained coast
   edge buffered by two local elements, minus the base's water, 7.44 km2 here
   -- is now filtered together with the OSM land, so a gap between the two is
   closed by the same rule. It lies outside the hole.
2. **The coarse transition.** Applied in every band, the 0.5 x h rule kept a
   200 m spike of OSM land 3 km west of the port. The spike reached the
   frozen interface, whose edges are 440-600 m, and left elements of 17.8 and
   28.0 deg with two C4 jumps between them. The interface is open water on
   both sides, so no closing can remove that gap. The rule is therefore kept
   to bands 0-1 (elements up to 2 h0 = 60 m), where the structures the
   owner asked about are; the coarser bands keep 2 x h. The owner agreed
   (2026-09-24).
3. **A wall judged against the target.** The Futtsu coast recipe (§17) kept
   a 108 m stub wall among 150-250 m transition elements, which made an
   element of 27.3 deg. A wall piece clipped to the hole is now judged
   against the LOCAL element at its midpoint, like the width filter and the
   extraction, not against the target.

### 15.5 Result (job 115724) and FVCOM (jobs 115731-115735)

| | §13 (115612) | §15 (115724) |
|---|---:|---:|
| QA | 20/22 | 20/22 |
| violations introduced by the patch | 0 | 0 |
| nodes / elements | 6,058 / 10,958 | 6,018 / 10,862 |
| wall pieces / edges | 18 / 46 | 9 / 17 |
| water within 1.25x of target | 100.0 % | 100.0 % |
| finished-mesh external step | 1.0 s | 1.25 s |

The two gates that fail are the base's own element 2101 (29.0 deg) and
min_depth, which is reported and not gated on this branch.

M2, 20 days, against the base:

- **Both runs:** complete (TADA) and finite.
- **At the five gauges:** walled minus base is -0.2 to -0.4 mm in
  amplitude and +0.013 deg in phase.
- **In the harbour:** the amplitude is 0.437-0.442 m, with no frozen node.
  The long breakwater separates a sheltered basin at 0.442 m.

## 16. Figures: one colour for every solid boundary (commit 56a503c)

Mesh figures had drawn the coast blue and walls red in one notebook, black
and red in another, and not at all in a third. `plotting` now fixes the
colours for every mesh figure:

- `SOLID_BOUNDARY_COLOR` (black): every solid boundary edge -- coastline,
  quay, and a wall represented as a line alike, since a wall is boundary on
  both of its sides;
- `OPEN_BOUNDARY_COLOR` (red): the open boundary;
- `MESH_EDGE_COLOR`: thin interior edges.

`draw_mesh`, `boundary_segments` and `boundary_legend` apply them; notebooks
429, 431, 432 and 433 use them, and `notebooks/434_final_mesh.py` draws the
delivered mesh alone at `MESH_PNG_DPI` -- the whole patch, plus any close-ups
given in `FMESH_VIEWS` (`name:x0:x1:y0:y1`, plus-separated).

## 17. The other recipes, after §14-§16 (jobs 115717-115725)

Rebuilt to check that nothing else moved:

| recipe | branch | before §14 | first run (115717-115719) | after §15.4 item 3 (115723-115725) |
|---|---|---|---|---|
| `futtsu_coast_hires` | hires, coast in the region | 0 introduced | **3 introduced** (27.3 deg) | **0 introduced**, QA 20/22 |
| `futtsu_nori_hires` | hires, offshore region | -- | 0 introduced | 0 introduced, QA 20/22 |
| `futtsu_nori` | default (no hires) | 0 introduced, QA 20/21 | 0 introduced, QA 21/22 | (not affected) |
| `kimitsu_port_hires` | hires, port | -- | -- | 0 introduced, QA 20/22 |

The default branch shares only `improve_patch`'s new flip veto and the new
QA gate with this work, and passes both. On every hires recipe the two gates
that fail are the base's element 2101 and min_depth, which is reported only.

## 18. What is left

- **Depths.** The minimum depth and the r-factor smoothing are the next step
  (`refine_coast_and_bathy_design.md`). The FVCOM runs above were finished
  with `fmesh-refine-depths` at a 3 m floor and r <= 0.2 as a stand-in.
- **Structures hugging the coast.** A structure within 0.4 of an element of
  the coast is still not represented as a wall; since §15 one wider than
  half an element is land instead.
- **OSM `man_made` lines.** They are not used (§6 item 4).
