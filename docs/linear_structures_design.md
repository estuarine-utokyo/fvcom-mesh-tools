# Breakwaters, piers and other linear structures: a design

Status: design, for the owner's decision. Nothing here is implemented.

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
| 5 | every wall is **zero-width**, and the footprint it hands to the water is reported |

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
