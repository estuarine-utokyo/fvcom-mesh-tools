# Refining the coastline and the bathymetry: design, revision 2

A new **option** on the local refinement recipe. Today a refinement re-cuts the
mesh inside a declared region and inherits the base mesh's bathymetry node for
node (`depths_from_base`, owner 2026-09-22). The request is for the other half:

1. the declared region may **contain coastline**, re-cut at the region's size;
2. the depths inside come from an **M7001-based product**, not the base mesh;
3. the seabed is **smoothed last**, over the region *and its transition*;
4. everything outside region + transition keeps the base, unchanged;
5. every behaviour that exists today keeps working, as the default.

Revision 2 follows `refine_coast_and_bathy_design_review.md` (gpt-6-astra,
2026-09-23), which raised three P1 findings against revision 1, and three
owner decisions taken after it. **Revision 1 was wrong in three places and one
of its measurements was wrong**; §0 says what changed and why, because a
design that quietly drops its errors teaches nobody.

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

## 1. The owner's decisions, and what they settle

| # | question | decision (2026-09-23) |
|---|---|---|
| 1 | the 3 m depth floor | **not kept fixed; adjustable at the end** |
| 2 | the depth source | **M7001, T.P.-corrected, and accept the range it supports** |
| 3 | the coastline source | **OSM** |

Decision 2 settles the review's P1-5 the same way the review argued it: source
coarseness is a fact to report, not grounds to refuse the option the owner
asked for. There is **no source-spacing gate** and **no JFA source** in this
design. §2 reports what M7001 can and cannot deliver; it does not veto.

Decision 3 is already the repository's resolved precedence:
`DATA_INVENTORY.md` records OSM land polygons (2026-03) as coastline
precedence #1, superseding MLIT C23, with the meshing shoreline being the
xcoast "true land" product (land polygons minus inland water). So the option
inherits an existing decision rather than making a new one.

Decision 1 is the one with design consequences, and they are larger than a
constant: see §4.3.

## 2. What M7001 delivers here, reported and not gated

Two T.P.-corrected M7001 products exist, and neither resolves a fishery.

| product | what it is | spacing |
|---|---|---|
| `M7001/TP/M7001_dem_tokyobay.nc` | the gridded DEM `interpolate_m7001_tp` reads | 181 × 222 m |
| `M7001/TP/M7001_TP.parquet` | the raw soundings, 3,950,314 rows | 128-188 m *in region* (§0.4) |

The raw soundings look dense until they are counted inside the region rather
than in a box around it: a nearest-neighbour statistic over a 3.6 × 4.4 km box
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
smoothing removed variation that the source retains -- and that is what the
option delivers. It is not a new survey resolution and nothing in the output
will claim to be.

**One consequence of decision 2 needs the owner's eye** (§9.1): M7001 has
**no intertidal data** -- `DATA_INVENTORY.md` §2 records that its marks start
≥ 1 m below chart datum -- and these fisheries are tidal flats. 22 % of the
soundings inside `banzu_nori` are shallower than 3 m, minimum 1.14 m.

## 3. Schema

New optional block. **Absent, the recipe behaves as it does today.**

```yaml
# ------- existing keys, unchanged -------
coastline: resample            # this option needs it; `preserve` stays the default
coastline_tolerance_m: 120     # stays a per-recipe veto (review P2-12)
rfactor_limit: base

# ------- NEW: absent means "inherit the base depths", i.e. today -------
bathymetry:
  source: m7001tp              # the only value; dem.m7001.interpolate_m7001_tp
  scope: hole                  # hole (the actual cut) | core
  blend: ramp                  # ramp | none
  hmin_m: 3.0                  # applied in the final stage; see 4.3
  hmax_m: 300.0
  rfactor: 0.2                 # the final smoothing; must not be given with rfactor_limit
```

Revision 1 proposed a `sources:` list with a generic `xyz_dir` loader. That is
withdrawn: the review is right that §2 does not require it, the owner scoped
the work to M7001, and `CLAUDE.md` forbids abstractions built for hypothetical
features. One named source, one existing function.

The raw sampler is **`interpolate_m7001_tp`, not `production_depths`**: the
latter already floors and runs an unmasked Beckmann-Haidvogel smoother over the
whole connectivity, which would move frozen depths and pre-empt the final local
smoothing this option promises.

Validation, all of which must be written as refusals and tested: unknown keys
rejected; `source` in the enum; `scope`/`blend` in their enums; `hmin_m`,
`hmax_m` finite, positive, ordered; `rfactor` in `[0, 1)`; **`bathymetry.rfactor`
and `rfactor_limit` may not both be declared** -- the base's own worst r is a
statement about the base's bathymetry and means nothing once the depths come
from elsewhere. `scope: core` with `blend: ramp` is refused as contradictory:
a ramp that is not evaluated over the transition is not a ramp.

`coastline_h_m` (revision 1, per region) is **withdrawn**. `coastline_points`
requires a position-dependent size and its docstring records a 0.8° minimum
angle when 30 m coastal spacing met a transition asking for 400 m triangles.
The driver already cuts the coastline at the **local** size field, which is the
target inside the core and the ramp outside it. Nothing needs adding.

## 4. Pipeline

Only the depth stages change. Cut, rim, fill, repair and stitch are untouched.

```
  cut -> rim (coastline resampled on OSM at the LOCAL size) -> hole_polygon
      -> fill (DistMesh, notebook) -> clean -> stitch
      -> improve_patch                 <-- node coordinates become final HERE
      --------------------------------------------------------------------
      -> [1] sample M7001 at the FINAL coordinates, new nodes only
      -> [2] blend across the transition            (4.1)
      -> [3] floor at hmin_m, cap at hmax_m, MASKED to new nodes
      -> [4] r-factor smoothing over the cut, frozen depths held fixed
      -> [5] reports and gates, then write
```

Stages 1-4 replace `refresh_depths` + `limit_rfactor` when `bathymetry` is
declared and are not imported when it is not. Stage 1 must follow
`improve_patch`: the repair slides boundary nodes along their coastline curve,
and a depth sampled before the move belongs to a coordinate the mesh no longer
has (review finding 9, 2026-09-22).

### 4.1 The blend weight, defined geometrically

Revision 1's size-field weight is withdrawn (§0.3). Let

* `d_core(x)` = distance to the declared region, 0 inside it (the union, when
  regions overlap, so overlap needs no separate rule);
* `d_int(x)` = distance to the **retained interface** -- the nodes the cut
  actually kept, taken from `select_patch`'s own rim, excluding the physical
  coastline.

```
    w(x) = d_int(x) / (d_int(x) + d_core(x))          w = 1 in the core, 0 at the interface
    depth(x) = w * m7001(x) + (1 - w) * base_interp(x)
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

`blend: none` writes the source throughout `scope` and leaves stage 4 to cope.
It exists so §0.1 can be measured rather than argued.

### 4.2 The final smoothing, and what it can and cannot do

Stage 4 is `limit_rfactor(elements, depths, movable, rmax, depth_min, depth_max)`
-- the existing function. It moves only `movable` depths and it reports
`n_over_rmax_movable`, `n_over_rmax_frozen_pair` and `converged`.

Three corrections to revision 1 follow from §0.1-0.2:

* the driver must **reject a candidate that does not converge**, exactly as it
  already rejects one that breaks the frozen contract. Non-convergence is an
  infeasible seam, not a cosmetic residue;
* a frozen-pair violation is **unaffected by the blend** -- both ends are
  fixed -- so revision 1's test "`blend: none` produces a frozen-pair
  violation" was testing the wrong thing. The right expectation is *movable*
  violations and `converged = False`, which is what the reviewer measured on a
  two-triangle probe with fixed depths 1 and 4 (ratio 4 > R² = 2.25);
* the a-priori inequality of §0.1, `h <= 2 H W rmax / D`, is computed in
  pre-flight from the region's own numbers and reported, so an infeasible
  request is visible before the meshing rather than after it.

### 4.3 The floor is a knob at the end, so the depth stage must be separable

Decision 1 -- "not kept, but adjustable at the end" -- is not satisfied by
making `hmin_m` a recipe key. A refinement run costs a seed search over the
whole fill-repair-stitch loop; re-running it to try a different floor would be
absurd. So stages 1-2 are **cached** and stages 3-5 are **re-runnable**:

* the driver writes `<case>_dep_blended.npy` -- the blended field before any
  clipping or smoothing -- beside the case, with the node map;
* a new CLI, `fmesh-refine-depths <case> --hmin --hmax --rfactor`, reads that,
  applies stages 3-5 and rewrites `<case>_dep.dat`. The mesh is not touched, so
  the frozen contract is not at risk and no seed search happens.

This also settles what "floor then smooth" means: the floor is applied first
and the smoothing is bounded by it, as `production_depths` does, but both are
re-done together for each floor the owner tries. It is cheap -- a Gauss-Seidel
sweep over the edges of one patch.

Two consequences the review caught (P2-10) and that this must carry:

* `run_qa`'s `min_depth_m` defaults to **2.0**, and the project's own fixed
  setting in `DATA_INVENTORY.md` is a 2 m clip while the goto2023 production
  file used 3 m. A declared floor of 1 m would fail the existing gate and a
  declared 3 m is not certified by it. The QA call must be given the declared
  floor;
* pre-flight's `depth_of` must forecast **the clipped, blended field**, not raw
  survey depths: `preflight` requires strictly positive depths and raw
  positive-up land elevations would be rejected before the floor could act. The
  authoritative time step stays the one measured on the written depths.

## 5. Reports and gates, labelled as what they are

The review is right that revision 1 called several things gates that had no
threshold. Split:

**Gates** -- a failure refuses the candidate:

| # | gate |
|---|---|
| G1 | frozen depths equal the base **exactly**, compared to the base through the node map, after **every** stage including export -- not by re-reading `limit_rfactor`'s own statistic, which compares against its own input and would show zero after an earlier unmasked clip |
| G2 | `limit_rfactor` converged: no movable edge over `rmax` |
| G3 | no **new** frozen-pair edge over `rmax` that the base did not already have (the existing rule, kept) |
| G4 | every replaced coastline stretch was matched to an OSM component. `_source_substring` returns `None` on failure and `_resample_on_source` then silently subdivides the base while `coastline_curve` returns the base as the reference -- so a fidelity check against that curve would certify zero departure without ever using the requested shoreline (review P2-12). Under this option that fallback is a refusal |
| G5 | the existing 21 QA checks, with `min_depth_m` set to the declared floor |

**Reports** -- measured, printed, never a veto:

| what | why |
|---|---|
| source spacing in each region, with population, CRS and polygon identity | §2; decision 2 says report, not gate |
| which product won at each node, and the area fraction each won | `interpolate_m7001_tp` itself falls back from the fine grid to the Kanto/SRTM blend, and a single `m7001tp` label would hide it |
| fraction of the core clipped by `hmin_m` and by `hmax_m`, separately for stage 3 and for stage 4 | §2; and stage 4 can drive nodes onto a bound that stage 3 did not |
| raw / blended / clipped / final depth statistics, area-weighted | so "the survey was delivered" can be checked rather than assumed. A dense but datum-shifted survey floored to a flat plateau passes every numerical threshold above |
| max r on retained-to-new edges, before and after stage 4 | the seam. Retained-to-retained edges are unchanged by construction and are not the interesting set |
| coastline departure from the **matched OSM component**, over segment interiors | not only endpoint distance: a one-way nearest test accepts a shortcut that omits a narrow inlet |
| the a-priori feasibility margin `h <= 2 H W rmax / D` per region | §4.2 |
| achieved dt from the final written depths | a diagnostic at the gravity-wave convention, not a namelist |

**The affected zone is the actual cut, not the analytic envelope.** `is_new`
means newly allocated, and `select_patch` grows its selection, so new nodes can
exist slightly outside `region.buffer(W)` and be blended, clipped and smoothed
while every retained node passes G1. The reports name the actual cut and its
enlargement over the request; `hole_clearance` already measures the latter.

The driver writes **no sponge file** (it passes neither `sponge` nor
`write_empty_spg`), so there is nothing sponge-shaped for late smoothing to
invalidate; revision 1 said otherwise. `obc_depth_control=False` must stay on
the export call -- the helper's default is `True` and it rewrites OBC depths,
which are frozen.

## 6. Where the code goes

| file | change |
|---|---|
| `src/fvcom_mesh_tools/bathy_patch.py` (new) | `blend_weights(...)` (§4.1), `patch_depths(...)` = stages 1-3, `smooth_patch_depths(...)` = stage 4 |
| `src/fvcom_mesh_tools/cli/refine_depths.py` (new) | `fmesh-refine-depths`, §4.3 |
| `src/fvcom_mesh_tools/refine.py` | schema; `preflight` forecasts the clipped field; `limit_rfactor` docstring corrected (§0.2) |
| `notebooks/420_local_refine.py` | branch at the depth stage; reject non-convergence; QA floor |
| `recipes/refine/futtsu_nori_bathy.yaml` (new) | the option on the existing Futtsu case |

Licence: numpy/scipy/shapely/netCDF4/pyproj only. DistMesh stays in the
notebook.

## 7. Backward compatibility, with a reference that exists

Revision 1 named "the committed output" as the oracle. `git ls-files outputs`
returns nothing, so that oracle does not exist. The reference is instead a
**controlled old-versus-new execution**: the recipe, the base revision, the
input file hashes, `LR_SEED`/`LR_MAX_ITER`, `FMESH_LAND` and the dependency
versions recorded, the run made before and after the change, and byte equality
required on `_grd.dat` and `_dep.dat`. One Futtsu recipe does not cover
`preserve`/`resample`/`spline`, multiple regions, or `rfactor_limit` off and
numeric; those keep their existing unit coverage.

## 8. Test plan

* `blend_weights`: 1 on the core, 0 on the retained interface for a cut that is
  **not** the analytic buffer, monotone between, 1 where `d_int + d_core = 0`,
  and not forced to 0 on the free coastline;
* the R² feasibility condition, against the two-triangle probe in the review:
  fixed 1 and 4 must not converge, fixed 1 and 2 must;
* `blend: none` on a step field produces **movable** violations and
  `converged = False` (not frozen-pair violations, §4.2);
* frozen depths bit-identical through stages 1-4 *and* through export,
  including when the limiter hits its round cap and when `hmin_m` differs from
  the base minimum -- the unmasked-clip path of G1;
* a coastline stretch whose OSM match fails is refused, not silently
  subdivided (G4);
* schema refusals: unknown key, `rfactor` with `rfactor_limit`, `scope: core`
  with `blend: ramp`, unordered bounds;
* `fmesh-refine-depths` re-run with a different floor changes only the depth
  file and leaves the grid byte-identical;
* regression: the controlled run of §7;
* integration: `futtsu_nori_bathy.yaml`, reporting the 188 m source spacing
  against a 30 m target as the honest statement of what was delivered.

## 9. Open with the owner

1. **The intertidal gap.** Decision 2 scopes the source to M7001-T.P., and
   M7001 has no intertidal data (`DATA_INVENTORY.md`: marks start ≥ 1 m below
   chart datum) while these fisheries are tidal flats -- 22 % of the soundings
   inside `banzu_nori` are shallower than 3 m. The repository's own resolved
   precedence already covers this: "M7001 wins where both have data; the
   `tokyo_bay` 30 m grid wins in the intertidal/shallow gap and anywhere M7001
   has no soundings". That grid is `depth_0030-11+12+13+14+15.nc`, measured at
   **27.7 × 34.0 m**, covering 139.565-140.172 E / 35.101-35.856 N -- both
   fisheries, 100 % finite. It is a 2021 regrid from a prior study, not a
   survey authority. Is applying the documented precedence in scope, or is
   M7001 alone the instruction?
2. **The floor's range.** §4.3 makes it adjustable. The project's fixed
   setting is a 2 m clip, the goto2023 production file used 3 m, and `run_qa`
   defaults to 2 m. Which of those is the starting value, and is there a floor
   below which the run should be refused rather than reported?
3. **`coastline_tolerance_m`.** Kept as a per-recipe veto, per the review: a
   larger global default would weaken requirement 5 for every existing recipe.
   Each recipe using this option raises it deliberately. Confirm.
