# Refining the coastline and the bathymetry: design, revision 3, and what it delivered

A new **option** on the local refinement recipe. Today a refinement re-cuts the
mesh inside a declared region and inherits the base mesh's bathymetry node for
node (`depths_from_base`, owner 2026-09-22). The request is for the other half:

1. the declared region may **contain coastline**, re-cut at the region's size;
2. the depths inside come from an **M7001-based product**, not the base mesh;
3. the seabed is **smoothed last**, over the region *and its transition*;
4. everything outside region + transition keeps the base, unchanged;
5. every behaviour that exists today keeps working, as the default.

Revision 3 follows `refine_coast_and_bathy_design_review.md` (gpt-6-astra,
2026-09-23), which raised three P1 findings against revision 1, and the four
owner decisions taken after it. **Revision 1 was wrong in three places and one
of its measurements was wrong**; §0 says what changed and why, because a
design that quietly drops its errors teaches nobody. Revision 2 survived a day:
the owner's fourth decision -- use the depths as they are, and adjust the
minimum depth in a later step -- removed two of its stages.

## 0. What revision 1 got wrong

**0.1 The blend does not make the seam r-factor-safe. Claim withdrawn.**
Revision 1 argued that ramping the survey into the base field makes the seam
continuous "by construction" and so removes the feasibility problem
`depths_from_base` warns about. Continuity is not a bound on an edge
difference. The review's counterexample: base 3 m, survey 300 m, ramp over
330 m; a new node 30 m inside the interface takes 30 m, and its edge to the
frozen 3 m node has `r = 27/33 = 0.818`. The ramp is continuous everywhere.

The correct statement is quantitative, and it is useful because it can be
checked *before* meshing. A ramp spreading a source-base difference `D` over a
transition width `W` has depth slope `D/W`, so across an edge of length `h` at
local depth `H` the r-factor is about

```
    r  ~  (D / W) * h / (2H)          and     r <= rmax   iff   h <= 2 H W rmax / D
```

That is a per-region inequality in numbers the recipe already has, and §5 makes
it a pre-flight report.

**0.2 The feasibility condition was stated wrongly.** Revision 1 said two
frozen neighbours of a movable node must be within `R = (1+r)/(1-r)` of each
other. The correct condition is `B/A <= R²` -- the intervals `[A/R, AR]` and
`[B/R, BR]` intersect. The looser truth does not rescue the argument: new
connectivity can join retained nodes the base never joined, and the depth
bounds add obstructions the ratio test does not see (a frozen 1 m node beside a
movable node floored at 3 m is infeasible at `rmax = 0.2` whatever the blend
does). `limit_rfactor`'s own docstring carries the same error and is corrected
with this work.

**0.3 The blend weight cannot be read out of the size field.** Revision 1 set
`w = (h_ambient - h(x)) / (h_ambient - h_target)`. `patch.patch_sizing` returns
`min(base_size(x), region contributions) / distmesh_scale` with
`distmesh_scale = 1.2`, and it ramps to **the base mesh's own size at that
point**, which varies in space. With a nominal 400 m ambient and a 30 m target
that formula gives `w = 1.014` in the core and `w = 0.180` at the nominal outer
edge -- neither 1 nor 0. §4.1 replaces it with a geometric weight.

**0.4 One measurement was wrong: the Futtsu circle.** Revision 1 built the
300 m circle as `Point(lon, lat).buffer(300 / 111000 / cos(lat))`, which
buffers by a constant in degrees and so stretches latitude by `1/cos(35.32)`
= 1.226. The area came out 0.3466 km² against the true 0.2827 km², and the
sounding count with it. Re-measured in UTM 54N, as the driver works:

| region | area | soundings inside | equivalent spacing |
|---|---:|---:|---:|
| `futtsu_nori` circle r = 300 m | **0.2827 km²** | **8** | **188.0 m** |
| `futtsu_nori` (two-beds polygon) | 0.2880 km² | 9 | 178.9 m |
| `futtsu_nori_east` | 0.2640 km² | 16 | 128.5 m |
| `banzu_nori` | 1.1999 km² | 37 | 180.1 m |

The conclusion is unchanged; one of its numbers was not.

**0.5 A question I asked the owner was answerable here.** The JFA product's
CRS is recorded in `docs/DATA_INVENTORY.md` §1: EPSG:6677 (JGD2011 Japan Plane
Rectangular CS IX), columns `id, Y_easting_m, X_northing_m, elev_m, flag`, with
`.prj` files present. I should have read the inventory before asking.

## 1. The owner's decisions

| # | question | decision (2026-09-23) |
|---|---|---|
| 1 | the depth source | **M7001, T.P.-corrected, and accept the range it supports** |
| 2 | the coastline source | **OSM** |
| 3 | where M7001 has nothing | **use the grid product; where that is still not enough, extrapolate and report the area** |
| 4 | the depth floor | **none here. Use the depths as they are. The minimum-depth adjustment is the NEXT step** |
| 5 | where the r-factor smoothing lives | **the next step, with the minimum depth** |
| 6 | how the option is shaped | **one branch taken at the top. Its first fork is: keep the region's coastline, or resolve it. Resolving means following OSM faithfully -- nothing else** |

Decision 1 settles the review's P1-5 the way the review argued it: source
coarseness is a fact to report, not grounds to refuse the option the owner
asked for. There is **no source-spacing gate** and **no JFA source** in this
design.

Decision 2 is already the repository's resolved precedence.
`DATA_INVENTORY.md` records OSM land polygons (2026-03) as coastline
precedence #1, superseding MLIT C23, with the meshing shoreline being the
xcoast "true land" product. The option inherits a decision rather than making
one.

Decision 4 is the one that reshapes the pipeline, and it makes it **smaller**:
see §4. It also answers a question of mine that was badly posed -- I asked
which of three existing floors should be the *starting value*, when the
instruction is that this step applies no floor at all.

Decision 6 removes a question I asked three times and should not have asked
once. `coastline_tolerance_m` refuses a resample that departs from the base
polyline by more than the tolerance, and I kept asking how large to make it --
when the point of this branch is that the mesh's coastline SHOULD move onto
OSM. The base polyline runs 300-600 m between nodes and sits up to 68.5 m from
OSM at its segment midpoints; recovering that is the request, not a violation
to be bounded. So the tolerance is not applied on this branch, and the
departure is reported instead. It stays exactly as it is on every other
branch, which is what keeps requirement 5.

Decision 5 keeps the two together, and the reason is that a floor is what
makes the classic r-factor well defined on every edge. Smoothing here would
have needed two limiters -- `r <= rmax` on the wet edges and a `|dh| / L` slope
cap on the intertidal ones -- and an intertidal slope limit nobody has chosen.
Once the floor is declared, `limit_rfactor` applies unmodified to the whole
patch. So the original instruction "smooth last, including the transition"
is honoured by the next step, not abandoned: **this step's output is a mesh
carrying the source's own depths, and it is not the finished bathymetry.**

## 2. What the sources deliver here, reported and not gated

### 2.1 M7001's resolution, in region

Two T.P.-corrected M7001 products exist, and neither resolves a fishery.

| product | what it is | spacing |
|---|---|---|
| `M7001/TP/M7001_dem_tokyobay.nc` | the gridded DEM `interpolate_m7001_tp` reads | 181 x 222 m |
| `M7001/TP/M7001_TP.parquet` | the raw soundings, 3,950,314 rows | 128-188 m *in region* (§0.4) |

The raw soundings look dense until they are counted inside the region rather
than in a box around it: a nearest-neighbour statistic over a 3.6 x 4.4 km box
gives 38 m at Futtsu and 17.7 m at Banzu, but that is dominated by along-track
spacing and it is the wrong number to quote. `sqrt(area / N)` is an
**equivalent areal spacing** for an explicitly defined point population, not a
coverage certificate and not an information-resolution estimate; §5 records the
population (finite `z_tp`, `mark == 'L'` excluded -- those are chart-datum zero
lines with no depth), the CRS and the polygon identity alongside the number.

The values are also quantised: every sounding inside `futtsu_nori` is 3.08 m,
and the polygon's whole range is 3.08-4.08 m.

So the honest claim, and the one the reports will make: **a 30-40 m mesh
sampling M7001 does not acquire 30-40 m bathymetric information.** It can still
change the represented field -- the base's own discretisation, floor and
smoothing removed variation the source retains -- and that is what the option
delivers.

### 2.2 "Where M7001 has nothing" is not where these fisheries are

Measured over area samples inside each region, as the fraction each product
returns a finite value for:

| region | M7001 fine, 181 m | `tokyo_bay` 30 m grid | Kanto blend, 380 m |
|---|---|---|---|
| `futtsu_nori` circle | **100 %**, 2.47-4.15 m | 100 %, 2.65-4.50 m | 100 %, 2.70-4.33 m |
| `banzu_nori` | **99.9 %**, 0.50-7.35 m | 100 %, 0.07-7.08 m | 100 %, -0.40-7.02 m |
| `futtsu_nori` (two-beds) | 100 %, 2.52-4.06 m | 100 %, 2.18-4.61 m | 100 %, 2.32-4.41 m |
| `futtsu_nori_east` | 100 %, 2.64-4.51 m | 100 %, 2.62-5.05 m | 100 %, 2.62-5.05 m |

The M7001 grid is 39.9 % finite overall, but the missing 60 % is land and the
water outside the bay. **Inside these fisheries it is complete**, so the
fallback decision 3 asks for will not engage here -- and that is itself the
result, because it means the design must not *claim* M7001 coverage as
evidence of M7001 measurement.

**A finite grid value is not a sounding.** The 181 m product interpolates
across its own gaps, so "M7001 exists here" at grid level says nothing about
whether a survey point is nearby. §5 therefore reports, per node, the
**distance to the nearest actual M7001 sounding**. Where that exceeds the
product's own grid interval, the value is the product's interpolation and the
report says so. That is the measurable form of "where M7001 has nothing", and
it needs no threshold chosen in advance.

One further note carried from `DATA_INVENTORY.md` §2: M7001's marks start
>= 1 m below chart datum, so it has **no intertidal data** -- and these
fisheries are tidal flats. The gridded product fills that gap by interpolation
rather than leaving it empty, which is exactly why the sounding-distance
report exists.

## 3. Schema

The option is **one branch, taken at the top** (decision 6). A recipe either
declares `hires` or it does not; when it does not, nothing below is imported
and the recipe behaves exactly as it does today.

```yaml
# ------- NEW. Its presence IS the branch. Absent = everything as today -------
hires:
  coastline: resolve         # resolve | preserve   <- the first fork
  bathymetry: tokyo_bay      # the ladder of 3.1, or `base` to inherit as today
  scope: hole                # hole (the actual cut) | core
  blend: ramp                # ramp | none
```

`coastline: resolve` places the new boundary nodes on **OSM**, faithfully.
There is no departure veto on this branch: `coastline_tolerance_m` exists to
catch a coastline that moved when it should not have, and here it is meant to.
The departure is measured and reported (§5).

`coastline: preserve` keeps the base polyline and only subdivides it, which is
geometrically identical to the base -- the same guarantee the existing
`preserve` mode gives.

The top-level `coastline` and `coastline_tolerance_m` keys are **refused
alongside `hires`**. They steer the other branch, and a recipe that sets both
is a recipe asking for two different things about the same coastline.

`bathymetry: base` is offered because decision 6 makes the coastline fork the
FIRST one: a region whose coastline needs resolving does not necessarily need
its depths replaced, and the two are now independent.

Revision 1 proposed a `sources:` list with a generic `xyz_dir` loader. That is
withdrawn: the review is right that the measurements do not require it, the
owner scoped the work to M7001, and `CLAUDE.md` forbids abstractions built for
hypothetical features. `bathymetry` names **one policy**, not a user-assembled
stack, and the policy is the repository's own documented precedence.

`hmin_m`, `hmax_m` and `rfactor` are **gone** from revision 2's schema:
decisions 4 and 5 put them in the next step, and a key that this step does not
act on has no business being validated by it.

Validation, written as refusals and tested: unknown keys rejected; each value
in its enum; `scope: core` with `blend: ramp` refused as contradictory, since a
ramp not evaluated over the transition is not a ramp; `hires` together with
`coastline` or `coastline_tolerance_m` refused.

`coastline_h_m` (revision 1, per region) is **withdrawn**. `coastline_points`
requires a position-dependent size and its docstring records a 0.8 degree
minimum angle when 30 m coastal spacing met a transition asking for 400 m
triangles. The driver already cuts the coastline at the **local** size field,
which is the target inside the core and the ramp outside it. Nothing to add.

**One change `resolve` does require.** `_resample_on_source` simplifies the
source to `0.25 x the MEDIAN size over the whole substring` before walking it.
A stretch that is fine in the core and coarse in the transition therefore has
its core detail simplified away at the transition's scale (review P2-11). On
this branch "faithful to OSM" is the instruction, so the tolerance must be
pointwise -- `0.25 x the local size at each vertex` -- and the existing
behaviour must stay untouched on the other branch.

### 3.1 The ladder

Per node, in order, first finite value wins:

1. `M7001/TP/M7001_dem_tokyobay.nc` -- the survey authority;
2. `tokyo_bay/depth_0030-*.nc` -- 27.7 x 34.0 m, inner bay, 100 % finite over
   139.565-140.172 E / 35.101-35.856 N. This is decision 3's "grid data";
3. `tokyo_bay/kanto_M7001_srtm_15s.nc` -- the 380 m Kanto blend, for the bay
   mouth and the shelf, which is what `interpolate_m7001_tp` already falls back
   to internally;
4. **extrapolation** -- nearest finite value of the last product that had any,
   with the distance recorded. Decision 3 requires the extrapolated area to be
   reported, so it is reported as geometry, not as a count: the nodes, their
   convex extent and their distance to real data.

`interpolate_m7001_tp` currently **raises** when a point is covered by
neither of its two products. It is not modified -- existing callers depend on
that refusal -- and the ladder is a new composition in `bathy_patch.py` that
reuses its grid readers.

The raw sampler is that function's interpolation, **not `production_depths`**:
the latter floors and runs an unmasked Beckmann-Haidvogel smoother over the
whole connectivity, which would move frozen depths and pre-empt the next step.

Validation, written as refusals and tested: unknown keys rejected; `source`,
`scope`, `blend` in their enums; `scope: core` with `blend: ramp` refused as
contradictory, since a ramp not evaluated over the transition is not a ramp.

`coastline_h_m` (revision 1, per region) is **withdrawn**. `coastline_points`
requires a position-dependent size and its docstring records a 0.8 degree
minimum angle when 30 m coastal spacing met a transition asking for 400 m
triangles. The driver already cuts the coastline at the **local** size field,
which is the target inside the core and the ramp outside it. Nothing to add.

## 4. Pipeline

Only the depth stages change. Cut, rim, fill, repair and stitch are untouched.

```
  cut -> rim (coastline resampled on OSM at the LOCAL size) -> hole_polygon
      -> fill (DistMesh, notebook) -> clean -> stitch
      -> improve_patch                 <-- node coordinates become final HERE
      --------------------------------------------------------------------
      -> [1] sample the ladder at the FINAL coordinates, new nodes only
      -> [2] blend across the transition                          (4.1)
      -> [3] reports and gates, then write
      ====================================================================
         NEXT STEP: minimum depth, cap, and the r-factor smoothing that the
         original request asks to be done last, over the transition too
```

Stages 1-2 replace `refresh_depths` when `bathymetry` is declared and are not
imported when it is not. Stage 1 must follow `improve_patch`: the repair slides
boundary nodes along their coastline curve, and a depth sampled before the move
belongs to a coordinate the mesh no longer has (review finding 9, 2026-09-22).

**Decision 4 removes the clipping stage.** The written depths are the
source's, including the intertidal ones, and three things revision 2 said
about that were wrong (owner, 2026-09-23):

* **the case runs.** A node above T.P. zero is not a broken mesh, it is a
  tidal flat, and FVCOM integrates it with wetting and drying. The report
  states the requirement -- `WET_DRY_ON` and a `MIN_DEPTH` in the namelist --
  and the dry and intertidal area, rather than calling the case unrunnable;
* **the slope is computable.** What is undefined for `h_i + h_j <= 0` is the
  particular expression `|h_i - h_j| / (h_i + h_j)`, not the seabed gradient.
  The reports carry `|dh| / L` in m/m, which is defined on every edge, beside
  the classic r wherever the classic r means something;
* **the time step is computable.** `c = sqrt(g H)` is bounded by the DEEPEST
  water, not the shallowest, so an intertidal node cannot threaten the
  external step; a declared minimum depth covers the wave speed where one is
  needed. What must change is that `preflight` currently *refuses*
  non-positive depths -- that refusal is correct for a base mesh and wrong for
  this option, so it becomes a report of the dry fraction.

### 4.1 The blend weight, defined geometrically

Revision 1's size-field weight is withdrawn (§0.3). Let

* `d_core(x)` = distance to the declared region, 0 inside it (the union, when
  regions overlap, so overlap needs no separate rule);
* `d_int(x)` = distance to the **retained interface** -- the nodes the cut
  actually kept, from `select_patch`'s own rim, excluding the physical
  coastline.

```
    w(x) = d_int(x) / (d_int(x) + d_core(x))     w = 1 in the core, 0 at the interface
    depth(x) = w * source(x) + (1 - w) * base_interp(x)
```

This answers the review's P1-2 and P1-3 together. It is exactly 1 on the core
and exactly 0 on the retained interface **whatever shape the cut actually
took** -- and `select_patch` selects whole faces by centroid and then grows the
selection to repair pinches, so the cut is not the analytic buffer and a weight
keyed to the buffer would not vanish where it must. It needs no ambient size,
no `distmesh_scale`, and no new recipe parameter.

Excluding the coastline from `d_int` is not a detail: the free coastline inside
the cut is hole boundary, and forcing `w = 0` there would give the newly
refined shore its *base* depths, which is the opposite of the request.

Where `d_int + d_core = 0` -- a node on the core boundary that is also on the
interface, which happens when a transition is clipped by land -- `w = 1`.

`blend: none` writes the source throughout `scope`. It exists so §0.1 can be
measured rather than argued.

### 4.2 What the blend does to the next step's job

The blend no longer has an r-factor claim attached to it (§0.1, withdrawn), but
it still decides whether the next step can succeed, and that is forecastable
here. A ramp spreading a source-base difference `D` over a transition width `W`
has depth slope `D/W`, so across an edge of length `h` at local depth `H`,

```
    r ~ (D / W) * h / (2H)        and      r <= rmax   iff   h <= 2 H W rmax / D
```

Pre-flight reports that margin per region from numbers the recipe already has.
A region that fails it is telling the owner, before any meshing, that the next
step will have to move depths a long way -- or that its seam is infeasible,
which `limit_rfactor` can only discover afterwards.

The other correction from §0.2 belongs with it: the feasibility condition
between the two frozen neighbours of a movable node is `B/A <= R^2`, not `R`,
and `limit_rfactor`'s docstring carries the error and is corrected with this
work. A frozen-pair violation is **unaffected by any blend** -- both ends are
fixed -- so revision 1's test "`blend: none` produces a frozen-pair violation"
tested the wrong thing. The right expectation is *movable* violations, which is
what the reviewer measured on a two-triangle probe with fixed depths 1 and 4
(ratio 4 > R^2 = 2.25).

## 5. Reports and gates, labelled as what they are

The review is right that revision 1 called several things gates that had no
threshold.

**Gates** -- a failure refuses the candidate:

| # | gate |
|---|---|
| G1 | frozen depths equal the base **exactly**, compared to the base through the node map, after **every** stage including export -- not by re-reading a limiter statistic, which compares against its own input and would show zero after an earlier unmasked change |
| G2 | on `coastline: resolve`, every replaced stretch was matched to an OSM component. `_source_substring` returns `None` on failure and `_resample_on_source` then silently subdivides the base while `coastline_curve` returns the base as the reference -- so a fidelity check against that curve would certify zero departure without ever having used OSM (review P2-12). "Follow OSM faithfully" makes that fallback a refusal, not a fallback |
| G3 | the existing 21 QA checks, **with the minimum-depth check reported rather than gated**. `run_qa`'s `min_depth_m` defaults to 2.0 and decision 4 writes depths as they are -- 0.07 m at Banzu. That is a tidal flat under wetting and drying, not a defect, so the check reports it; it is re-gated in the next step, once a floor is declared |

**Reports** -- measured, printed, never a veto:

| what | why |
|---|---|
| per node: which rung of the ladder won, and the area fraction each won | decision 3, and `interpolate_m7001_tp` itself falls back internally, which a single source label would hide |
| the **extrapolated** nodes: geometry, extent and distance to real data | decision 3 asks for the area, so it is reported as area |
| per node: distance to the nearest real M7001 **sounding**, and the fraction of each region beyond the product's grid interval | §2.2. This is "where M7001 has nothing", measured rather than assumed |
| source spacing in each region, with population, CRS and polygon identity | §2.1; decision 1 says report, not gate |
| raw / blended depth statistics, area-weighted, per region | so "the survey was delivered" can be checked rather than assumed |
| max r on retained-to-new edges | the seam the next step inherits. Retained-to-retained edges are unchanged by construction and are not the interesting set |
| the feasibility margin `h <= 2 H W rmax / D` per region | §4.2 |
| coastline departure from the **matched OSM component**, over segment interiors | not only endpoint distance: a one-way nearest test accepts a shortcut that omits a narrow inlet |
| coastline departure from the **base** polyline | how far the shore moved. Reported on this branch, never a veto (decision 6) |
| depths at or below zero, and the area they cover | the wetting-and-drying requirement, and the next step's input |
| edge slope as `\|dh\| / L` on every edge, and the classic r where `h_i + h_j > 0` | the seam and the seabed gradient, in a measure that survives an intertidal node |

**The affected zone is the actual cut, not the analytic envelope.** `is_new`
means newly allocated, and `select_patch` grows its selection, so new nodes can
exist slightly outside `region.buffer(W)` and be blended while every retained
node passes G1. The reports name the actual cut and its enlargement over the
request; `hole_clearance` already measures the latter.

The driver writes **no sponge file** (it passes neither `sponge` nor
`write_empty_spg`), so there is nothing sponge-shaped for late smoothing to
invalidate; revision 1 said otherwise. `obc_depth_control=False` must stay on
the export call -- the helper's default is `True` and it rewrites OBC depths,
which are frozen.

## 6. Where the code goes

| file | change |
|---|---|
| `src/fvcom_mesh_tools/dem/tokyo_bay.py` (new) | the ladder of §3.1: `sample(lon, lat) -> (depth, rung, distance_to_data)`, plus `sounding_distance(lon, lat)` |
| `src/fvcom_mesh_tools/bathy_patch.py` (new) | `blend_weights(...)` (§4.1), `patch_depths(...)` = stages 1-2, and the reports of §5 |
| `src/fvcom_mesh_tools/patch.py` | a pointwise simplification tolerance for `resolve`, leaving the existing median rule on the other branch |
| `src/fvcom_mesh_tools/refine.py` | the `hires` schema and its refusals; `preflight` stops rejecting non-positive depths on this branch; `limit_rfactor` docstring corrected (§0.2) |
| `notebooks/420_local_refine.py` | **the top-level `if hires`** (decision 6), the coastline fork, the depth stage, QA minimum-depth reported not gated, and a report that says the case still needs the next step |
| `recipes/refine/futtsu_nori_bathy.yaml` (new) | the option on the existing Futtsu case |

Licence: numpy/scipy/shapely/netCDF4/pyproj only. DistMesh stays in the
notebook.

## 7. Backward compatibility, with a reference that exists

Revision 1 named "the committed output" as the oracle. `git ls-files outputs`
returns nothing, so that oracle does not exist. The reference is instead a
**controlled old-versus-new execution**: the recipe, the base revision, the
input file hashes, `LR_SEED` / `LR_MAX_ITER`, `FMESH_LAND` and the dependency
versions recorded, the run made before and after the change, and byte equality
required on `_grd.dat` and `_dep.dat`. One Futtsu recipe does not cover
`preserve` / `resample` / `spline`, multiple regions, or `rfactor_limit` off
and numeric; those keep their existing unit coverage.

## 8. Test plan

* `blend_weights`: 1 on the core, 0 on the retained interface for a cut that is
  **not** the analytic buffer, monotone between, 1 where `d_int + d_core = 0`,
  and not forced to 0 on the free coastline;
* the ladder: each rung wins where the one above it is NaN; extrapolation
  engages only when all three are NaN and records its distance; a point covered
  by none of them is reported, never silently zero;
* the sounding-distance report reproduces §2.1's counts for the four polygons,
  measured in UTM as the driver works -- not the degree-buffered values
  revision 1 quoted (§0.4);
* frozen depths bit-identical through stages 1-2 *and* through export,
  including when the source is negative at a node adjacent to the interface;
* a coastline stretch whose OSM match fails is refused, not silently
  subdivided (G2);
* `hires.coastline: preserve` leaves the base polyline geometrically identical,
  and `resolve` moves it onto OSM with no tolerance veto -- the two forks of
  decision 6, tested as two forks;
* the simplification tolerance is pointwise on `resolve`: a stretch fine in the
  core and coarse in the transition keeps its core detail, which the existing
  median-over-the-substring rule destroys;
* schema refusals: unknown key, `scope: core` with `blend: ramp`, `hires`
  beside `coastline` or `coastline_tolerance_m`, and any of revision 2's
  removed keys (`hmin_m`, `rfactor`) rejected rather than ignored;
* the R^2 feasibility condition, against the two-triangle probe in the review:
  fixed 1 and 4 must not converge, fixed 1 and 2 must (this guards the
  corrected docstring, and the next step);
* regression: the controlled run of §7;
* integration: `futtsu_nori_bathy.yaml`, reporting 188 m source spacing against
  a 30 m target, the ladder rungs used, and the depths below zero that the next
  step will have to deal with.

## 9. Open with the owner

Nothing. The six decisions of §1 settle the source, the fallback, the
shoreline, the floor, the smoothing and the shape of the option. The intertidal
question became a report rather than a decision once §2.2 made it measurable,
and the coastline tolerance became a non-question once decision 6 made the
branch explicit.


## 10. What it delivered

Implemented and run on OCTOPUS (job 115419, `recipes/refine/futtsu_coast_hires.yaml`,
a 700 m circle straddling the Futtsu shore at a 30 m target, five seeds).

| | delivered |
|---|---|
| mesh | NP 6,153 / NE 11,412, from a base of 3,210 / 5,645 |
| achieved resolution | **29.5 m** median cell against a 30 m target, **100 %** of the water within 1.25x |
| coastline | **0.00 m** from the curve each stretch was cut from, over 142 chords |
| bathymetry | M7001 2,500 / 30 m grid **589** / Kanto 0 / extrapolated **0** |
| depths | -2.44 to 28.58 m, 126 nodes at or above the datum, nothing floored or smoothed |
| frozen zone | coordinates, depths and retained faces identical; `frozen_exact: true` |
| size field | max slope **0.349** against a C4 reference of 0.414, **0.00 %** of samples above it |
| QA | **19/21, 0 introduced by the patch**; the two failures are the base's own |
| geometry | boundary **simple**, 0 crossing pairs, 0 duplicate nodes, 0 non-manifold edges |
| minimum angle | **28.992 deg**, at element 2101 -- the base's own, 19.7 km away |
| achieved dt | 2.70 s (predicted 3.46 s) |

The 589 nodes from the 30 m grid are decision 3 doing real work: M7001 has no
intertidal data and this region is a tidal flat. Nothing was extrapolated.

### 10.1 Four attempts at one defect

The first working mesh took four passes, and three of the four diagnoses were
wrong. They are recorded because the wrong ones cost the most time.

| # | diagnosis | outcome |
|---|---|---|
| 1 | shoreline detail near the core, at the 7.5 m simplification tolerance | **wrong** -- the offenders were 6 km away, in the coarse outer transition |
| 2 | follow the source only where the mesh is finer than the base polyline | **no effect** on the offenders (6.03 to 7.30 deg, same two elements) |
| 3 | the size field is a cliff outside the base mesh | **right** -- max slope 42 to 1.005 |
| 4 | the moved boundary crosses the frozen one | **right** -- and 0 introduced violations |

The measurement that found (3) was in every report from the first run and I
had not read it: every hires run recorded a size-field slope of 35-42 with 6-11
"partial support" samples, where the same patch on the default branch records
0.375 and none. `base_size_field` returns the field's MAXIMUM beyond the base
mesh, which is safe while the hole stays inside -- and `coastline: resolve`
takes it outside.

(4) was found by `notebooks/422_mesh_validity.py`, written for the purpose:
matplotlib's point locator refused a mesh `verify_patch` had passed, and the
difference was exactly one crossing pair, a retained boundary edge against a
resolved one. `verify_patch` does not test that the boundary is simple.

(2) is kept even though it did not fix this: it is what stops the other
failure, the one the offshore recipe showed -- resolving a coastline in a
400-1700 m field replaced 17 base nodes with 13, a coastline COARSER than the
base's.

### 10.2 What the offshore recipe is for

`recipes/refine/futtsu_nori_hires.yaml` asks for the same thing over 300 m of
open water 2 km from the shore. It is kept as the counter-example: the
coastline it touches lies at the outer edge of the transition, so `resolve`
has nothing to work with there and says so. **A coastline is cut at the local
size, so refining one needs the region to reach it.**

### 10.3 Still true, and still reported

* the depths are the source's. The case needs `WET_DRY_ON` and a `MIN_DEPTH`,
  and the minimum depth and the r-factor smoothing are the next step;
* 1,581 nodes are below run_qa's 2 m floor, which is reported and not gated;
* one stretch was kept on the base polyline to avoid a crossing, and its
  departure from OSM (up to 147 m) is in the report rather than hidden by the
  fact that it is 0 m from the curve it was cut from;
* M7001 still does not resolve 30 m. The mesh is refined and the seabed is
  resampled, and `sounding_distance` is what says which is which.
