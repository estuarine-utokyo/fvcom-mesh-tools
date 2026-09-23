# Review of the coastline and bathymetry refinement design

Reviewed against the current repository on 2026-09-23. This is a review of
`refine_coast_and_bathy_design.md`, not a claim that its proposed code exists.
The design needs revision before implementation, principally because its blend
does not establish feasibility, its spacing gate does not establish coverage,
and its proposed source policy goes beyond the requested M7001 option.

P1 blocks the design; P2 must change before code; P3 is worth considering.
Code references below describe existing behavior that the proposal must handle.

1. **P1 — Continuity does not remove the r-factor feasibility problem.**
   **Sections 4.1, 4.2 and 8.** Even an ideal continuous ramp with exactly zero
   weight on the interface does not bound differences across finite mesh
   edges. For example, let base depth be 3 m, survey depth 300 m, and weight
   increase linearly from zero at the interface to one 330 m inward. A new
   node 30 m inward has depth 30 m; its edge to a frozen 3 m node has
   `r = 27/33 = 0.81818`. The ramp is continuous everywhere. This is a
   mathematical counterexample, not a measurement on the production mesh.

   More fundamentally, a blend changes the initial guess, not the feasible
   set imposed by connectivity, fixed depths and depth bounds. With
   `R = (1+r)/(1-r)`, a movable node adjacent to fixed depths A and B needs
   the intersection of `[A/R, AR]` and `[B/R, BR]`. For B >= A the condition
   is **B/A <= R²**, not R as stated in the design and the limiter's current
   docstring. At r=0.2, a ratio of 2 is feasible; a ratio of 4 is not.
   More generally, fixed endpoints joined by a path of k edges need a ratio
   no greater than R^k. New connectivity can shorten such paths. Adjacent
   fixed rim nodes satisfying the base limit does not establish this property
   for all fixed neighbors or paths in the new mesh. Bounds add further
   obstructions: a fixed 1 m node cannot adjoin a movable node floored at
   3 m under r <= 0.2.

   The existing limiter was exercised on two triangles `[0,1,2]`, `[1,3,2]`,
   with nodes 0 and 3 fixed. Fixed depths 1 and 4 left **2 movable violations,
   0 frozen-pair violations**, and `converged=False`. Fixed depths 1 and 2
   converged. Fixed endpoints of 1 with a movable floor of 3 left **4 movable
   violations, 0 frozen-pair violations**. All fixed depths stayed exact.
   Thus the proposed `blend:none` test expecting a frozen-pair violation is
   the wrong test: blending new nodes cannot change a frozen-pair r-factor
   at all. Preserve the driver's rejection of both movable violations and
   newly introduced fixed-fixed violations. Describe blending as an initial
   field construction, measure its effect, and refuse infeasible/nonconverged
   candidates. Do not promise that it eliminates the problem.

2. **P1 — The proposed weight cannot be recovered from the actual size field.**
   **Sections 4.1 and 6.** `patch.patch_sizing` returns
   `min(B(x), contributions)/distmesh_scale`, with scale 1.2 by default.
   `B(x)` is spatially varying; the driver's per-region ambient median only
   sets transition width. There is no single ambient/target pair that inverts
   the combined field.

   Taking “the size field the fill was given” literally, with constant ambient
   400 m and target 30 m, yields w=**0.18018 at the nominal outer edge**, where
   the field is 400/1.2. In the core it yields w=375/370, greater than one.
   Removing the scaling alone does not fix spatial ambient variation or
   overlap: the minimum may come from another region with a different target
   and width. Clipping w into [0,1] conceals the ambiguity. If B equals the
   target, division is undefined; if B is already finer than the target,
   geometry still permits a bathymetry replacement but the size deficit no
   longer expresses that request.

   Define the bathymetry weight independently of this lossy scalar field,
   using explicit region membership and transition semantics, including
   overlaps and zero-width transitions. Reusing the existing geometric
   transition width need not introduce another user parameter. State how the
   core can retain w=1 when mesh size needs no refinement.

3. **P1 — The interface is not the transition's outer level set.**
   **Sections 4 and 4.1.** `select_patch` selects whole faces by centroid and
   grows the selection to repair pinches/spikes. Its actual interface can
   cross either side of the requested buffer boundary. Consequently, even a
   corrected geometric ramp can be nonzero at a retained interface node.
   Overriding that node with its base depth gives a mismatch between the
   prescribed interior field and its boundary trace. Evaluating weights only
   on new nodes does not make their values tend to zero at that interface.
   Linear finite-element interpolation on the stitched mesh is continuous
   because it shares nodes, irrespective of this blend; that fact says
   nothing about edge slopes.

   Also distinguish the **retained interface** from the **physical coast**:
   the free coastline inside the hole must normally take survey depths. It
   must not be forced to w=0 just because it belongs to the hole boundary.
   Frozen stretch endpoints need explicit handling, and survey/base blending
   on newly wetted land may require a base extension outside the old mesh.
   `depths_from_base` currently uses nearest-node fallback there; this is not
   a guaranteed continuous extension. Specify the weight against the actual
   interface, the coastal anchors, and the base support domain. Test final
   retained-to-new edges, not merely nominal buffer sample points.

4. **P2 — The density calculation is a useful warning, not a resolution gate;
   one supplied area is demonstrably wrong.** **Sections 2.3, 5/C23 and 8.**
   `sqrt(area/N)` has units of metres and is a valid *equivalent areal spacing*
   for an explicitly defined point population. It is not across-track
   spacing, a coverage certificate, or an information-resolution estimate.
   Arbitrarily many points on one track can make it pass while most of the
   polygon has no nearby observations. Conversely, a small polygon between
   well-supported grid samples can have N=0 yet be interpolable from nearby
   samples. Counting land area, duplicate tile points, contour vertices or
   invalid soundings changes the meaning again.

   The declared 300 m circle has area approximately **0.28274 km²**, not
   0.346 km². The driver's exact construction, `buffer(300, quad_segs=64)`,
   measured **282,714.9526 m²** in this review. With nine points that would
   imply **177.2365 m**, not 196 m. This does not validate the count of nine:
   the point selection must be repeated using the same projected geometry
   as the driver. The discrepancy is consistent with a projection/geometry
   error, but its cause was not reproduced. Do not make the four table values
   golden unit-test answers.

   Count unique, finite depth observations, explicitly excluding `mark == 'L'`
   with missing `z_tp`; do not interpret those chart-datum lines as T.P. zero
   soundings. Record CRS, polygon identity, wet-area mask, exclusions and
   effective product. Assess spatial support over area samples, including
   holes and poorly supported lobes, with local distances/interpolation
   support and explicit tolerances. For a regular product report its native
   grid intervals; count neither resampled output pixels nor mesh nodes as
   independent survey evidence. The quoted nearest-neighbor statistic can
   be dominated by along-track spacing, but asserting it *is* along-track
   spacing also needs directional evidence.

5. **P1 — Source coarseness does not authorize refusing the requested M7001
   option or substituting a different survey.** **Sections 2.3, 3 and 9.2.**
   A 40 m mesh sampling an approximately 180 m DEM does not acquire 40 m
   independent bathymetric detail. That warning is justified. The stronger
   claim that it contains “no information the base did not already have”
   does not follow: the base discretization, floor and smoothing can have
   removed variations still present in that DEM. Source resampling and
   geometric refinement can therefore change the represented field without
   claiming a new survey resolution.

   The owner's request is for depths from an **M7001-based product**. It does
   not require source spacing <= mesh target. A mandatory refusal at Banzu
   is a new scientific acceptance policy, and the proposed JFA-first Futtsu
   recipe is a different source choice. Neither follows from the density
   calculation. Separate honest source-resolution reporting from the owner's
   decision whether unsupported fine-scale information should block a run.
   Confirm authorization for JFA before making it the first implementation's
   purpose. The M7001-only option remains a coherent implementation of the
   stated request even if it cannot claim 40 m bathymetric information.

6. **P2 — The source stack can pass C23 while sampling the wrong or poorly
   supported product.** **Sections 3, 5 and 6.** “First source that covers a
   node” does not define coverage. A bounding box, convex hull, valid grid
   cell, and bounded-distance interpolation support are different. A convex
   hull over tiles can bridge missing tiles, land or channels. Nearest-neighbor
   extrapolation can let JFA claim the entire bay. Fine points concentrated
   in a small part of the core can dominate `spacing_in(geom)` while most
   nodes actually use M7001 fallback. The built-in sampler itself can fall
   back from M7001 to a Kanto/SRTM blend; a single `m7001tp` source label
   would hide that distinction. Hard first-wins boundaries inside the core
   introduce additional source seams that the outer ramp does not address.

   Define support, interpolation, extrapolation limits and missing-data
   behavior per supported product. Return actual product provenance and
   support diagnostics at sampled coordinates, and report area coverage by
   the source that actually wins. Test gaps, tile overlaps, coast/channel
   crossings and a fine-source footprint ending inside a core. Validate
   source coverage over the required transition as well as the core.

7. **P2 — The schema has unresolved semantics, and the generic loader is
   premature for the stated task.** **Sections 3 and 6.** A list itself is
   inexpensive, and if JFA plus M7001 is approved it has two concrete uses;
   it is not inherently a forbidden abstraction. But §2.3 does not require a
   generic `SourceStack` or `xyz_dir` reader. For the current request, a
   named M7001 product using `interpolate_m7001_tp` is enough. That existing
   interpolation function, **not `production_depths`**, is the reusable raw
   sampler: the latter already floors and smooths on connectivity with no
   fixed mask, contrary to the promised final local smoothing.

   Before code, specify required keys/defaults, unknown-key rejection,
   finite positive ordered bounds, r in [0,1), source identity, horizontal
   and vertical units/datum, path expansion, and the precedence/conflict
   rule between `bathymetry.rfactor` and `rfactor_limit: base|off|number`.
   Define all `scope`/`blend` combinations or remove unneeded combinations:
   `scope:core` cannot both leave the transition unsampled and perform the
   stated survey ramp throughout it. Explain whether the source is required
   where w=0. Section 8's `{sources: [m7001tp]}` is not the mapping schema
   in §3 and omits values called “declared, not inherited.” As written, an
   implementation cannot know whether to accept that test recipe.

8. **P2 — Clipping and the meaning of “frozen” leave depth-change paths that
   the limiter mask does not close.** **Sections 4 and 5/C22.** Stage 3 says
   “floor/cap” without a mask. Applying `np.clip` to the whole array changes
   retained values when the declared bounds differ from the base. The
   limiter's `max_frozen_depth_change_m` compares against its own input and
   would report zero after such an earlier change. C22 must compare to the
   original base through the node map after **all** stages, including export,
   not just reuse that limiter statistic. The existing `verify_patch` and
   native read-back check are the stronger controls; retain them.

   Separately, `is_new` means newly allocated, not inside the declared
   polygon-plus-buffer. Selection growth can remove old nodes outside that
   envelope, and new nodes there can be blended, clipped and smoothed while
   every retained node passes C22. This is already a distinction in the local
   refinement contract, but the new document conflates `hole` with
   `core + transition`. Explicitly adopt and report the actual cut as the
   affected zone, with a permitted enlargement, or enforce the owner's
   literal geographic freeze. Even w=0 outside the nominal buffer does not
   prevent the subsequent limiter from changing those new depths.

   I found no hidden mutation of retained depths in the current masked
   limiter. Export is another possible mutation path, but notebook 420
   already blocks it with `obc_depth_control=False` and checks the written
   native depths. Preserve that call explicitly; the export helper's default
   is True. A global `production_depths` call would reopen the path too.

9. **P2 — C22–C27 do not establish delivery of survey bathymetry.**
   **Section 5.** C24 has no acceptance threshold, C25 no numerical acceptance
   rule, C26 no tolerance or endpoint policy, and C27 no timestep criterion.
   They are not all gates as written. Reporting can be appropriate, especially
   for timestep, but call it reporting. Define the edge set in C25: the
   retained-retained interface segments are already unchanged; the useful
   test includes retained-to-new edges and all changed connectivity. “Same
   as the base” is undefined for edges the base did not have and is not
   implied by blending. The existing 21-check battery does not gate r-factor;
   retain the separate driver rejection discussed in finding 1.

   A dense but wrongly signed/datum-shifted survey can be floored to a smooth
   3 m plateau, keep frozen depths exact, pass the spacing and geometry
   checks, and yield a perfectly computable dt. C24 reports the clipping
   without refusing it. Similarly, final smoothing can erase most survey
   structure while satisfying every stated numerical threshold. Add product
   identity/datum verification and retain separate raw, blended, clipped and
   final statistics. Measure source-to-final changes and where they occur,
   preferably with area weighting rather than node weighting. Decide what
   changes constitute an unacceptable loss of the requested field, rather
   than treating r <= 0.2 as evidence of bathymetric fidelity. Stage-3 clipping
   fractions alone do not count additional bound hits caused by stage 4.

10. **P2 — “Only the depth stages change” conflicts with QA and preflight.**
    **Sections 4, 5 and 6.** Sampling after the final `improve_patch` pass and
    smoothing before acceptance is the right order. Current geometry repair
    does not depend on depths. Current `run_qa` calls, however, default to
    `min_depth_m=2.0`. A deliberately declared 1 m floor will still fail that
    gate if the rest of the driver is untouched. Conversely, a new 3 m
    requirement is not certified by the existing 2 m QA setting. Specify
    patch-specific depth acceptance while preserving inherited exceptions
    and the unchanged default path.

    Simply replacing preflight's `depth_of` with raw survey sampling is also
    wrong. It requires strictly positive depths, so valid positive-up JFA
    land elevations converted to negative depths in water selected by a
    different shoreline are rejected before the promised floor can operate.
    Even wet raw depths differ from the clipped/blended/final field used for
    dt. Use a forecast with the declared clipping and support policy, label
    smoothing uncertainty, and keep final dt authoritative. Preflight's
    degree-space geometry also differs from the driver's re-struck projected
    circles; use the same region definition for the new density/coverage
    gate, not another approximation.

    Existing notebook 420 already serializes after limiting, runs QA on that
    result, and computes achieved minimum-altitude dt from `written.depths`
    (lines 806–837). Preserve this order for each candidate and final output;
    compare the final native case's new depths too if it is the authoritative
    artifact. The OBC list/type is mapped from the base and needs no
    bathymetric rebuilding. This driver supplies neither `sponge` nor
    `write_empty_spg` to export, so it does **not** write a sponge file; there
    is no existing sponge calculation here for late smoothing to invalidate.
    Copying an external sponge file with old node IDs would need a separate
    mapping contract. A reported dt is still a diagnostic at the stated
    gravity-wave convention, not an updated FVCOM namelist or run-step choice.

11. **P2 — The new coastline size key can reintroduce an already measured
    geometry failure.** **Sections 1, 3 and 4.** `coastline_points` explicitly
    requires a *position-dependent* size, not a scalar region target. Its
    docstring records a minimum angle of 0.8 degrees when 30 m coastal
    spacing met a transition requesting 400 m triangles. The new
    `coastline_h_m` is not assigned a transition or overlap rule, and may
    conflict with the interior target even inside a core. Passing it as a
    constant over the whole free stretch could make every fill fail.

    Also, the current simplifier uses **0.25 times the median size over the
    source substring**, not a pointwise local tolerance. A mostly coarse
    transition can therefore erase fine core features before 30 m resampling
    begins. Sampling the simplified chord densely cannot recover them.
    Specify compatible boundary/interior sizing and preservation of core
    shoreline detail. C26 needs a scale-appropriate tolerance that tests the
    delivered segments against the original selected source, not merely a
    large tolerance against the base.

12. **P2 — Existing coastline resampling can silently preserve the base or
    omit requested source geometry.** **Sections 1, 2.1 and 5/C26.**
    `_source_substring` can return None; `_resample_on_source` then subdivides
    the base, and `coastline_curve` returns the base as the reference. A
    verification against that curve would certify zero departure without
    using the requested shoreline. `rim_constraints` also explicitly keeps
    an entirely free island's original ring instead of resampling it.
    Existing resampling replaces identified stretches; it is not an overlay
    that discovers new source islands/channels inside the hole.

    Make failure to match the required source a refusal for this option,
    or declare/report an explicit exception. State the supported shoreline
    topology and test an island, a missing source match and frozen anchors
    offset from the source. C26 must identify the matched source component
    and cover segment interiors; a one-way nearest-distance test alone can
    accept a shortcut omitting a narrow inlet. Fixed endpoints cannot be
    simultaneously exact base coordinates and on a displaced source, so
    specify their joining allowance separately. Raising the old tolerance
    does not settle any of these questions.

13. **P2 — The JFA schema's nodata convention is ambiguous, and the CRS
    question already has repository evidence.** **Sections 2.2, 3 and 9.5.**
    `docs/DATA_INVENTORY.md`, section 1, records EPSG:6677, T.P. signed
    elevations, the column order `(id, Y_easting, X_northing, elevation, flag)`
    and the existence of `.prj` files. I read `SHP/09MD5920.prj`: it declares
    JGD2011 Japan Zone 9. This is stronger evidence than “the footprint lands
    nearby”; horizontal agreement within 0.2 m does not establish a vertical
    datum or interchangeable CRSs.

    The first twelve records of `05MCSV_TP/09md5839_05g.txt` have plausible
    elevations from -1.10 to -1.20 m in column four and **-9999 in column
    five**. The schema's unqualified `nodata: -9999` and §9's “nodata in the
    fifth column” do not say whether those rows are valid. A generic filter
    that rejects any row containing the sentinel would discard these depth
    records and silently fall through to coarse M7001. Separate elevation
    nodata from flag semantics and confirm the product format. The repository
    documents the vertical convention, but external provenance may still be
    needed to confirm flag meaning and survey metadata. A rectangular mosaic
    extent and nominal lattice spacing alone do not prove valid coverage of
    the whole circle.

14. **P2 — The backward-compatibility oracle named in the design is absent.**
    **Sections 7 and 8.** `git ls-files outputs` returns no files in this
    checkout. Thus “against the committed output” does not identify a usable
    baseline. One Futtsu recipe also cannot certify preserve/resample/spline,
    polygon/multiple-region behavior and rfactor off/numeric defaults.
    Name a reproducible reference (baseline revision, external input hashes,
    seed sequence, relevant environment overrides and dependency versions),
    and compare old/new executions with the option absent. Keep byte equality
    for the claimed serialized files under that controlled setup, plus
    focused existing-mode coverage. These are proposed tests, not tests that
    can presently be claimed to pass for this unimplemented feature.

15. **P3 — Replace several owner questions with decisions the repository
    cannot answer; resolve the code questions locally.** **Section 9.**
    The five subjects are relevant, but they are not the five unresolved
    decisions in their current wording:

    - **Floor:** keep the question, including an explicit wet/dry and
      above-datum treatment if lowering it is contemplated. The 3 m baseline
      is already established in code and recipes. Even retaining 3 m while
      changing source changes the physics; a floor choice alone does not
      preserve a like-for-like seabed comparison.
    - **Banzu/source policy:** ask whether coarse M7001 resampling is wanted
      with disclosed support limits, or whether independent fine bathymetric
      information is a prerequisite. `bathymetry: {scope: core}` plus `warn`
      **does not inherit base depths**: it still samples the source in the
      core. Omitting the bathymetry block is the existing inherit option.
      “No suitable source found” is not evidence that none exists.
    - **Shoreline:** ask which physical shoreline convention should define
      the mesh. A zero T.P. elevation contour is an equipotential level
      intersection, not automatically the appropriate coastal boundary for
      this model or a vector shoreline's tide/date convention. A 3 m floor
      also destroys the claimed physical agreement at a zero-depth contour.
      Retaining the vector source with a clearly labeled contour diagnostic
      is reasonable, but “same survey” alone does not establish consistency.
    - **Extent and fidelity:** ask what actual cut enlargement and survey
      distortion from clipping/smoothing are acceptable. These affect the
      owner's unchanged-exterior and survey-delivery requirements directly.
      Keep a deliberate per-recipe coastal-departure veto as the current
      default; relaxing a global default would contradict requirement 5.
    - **Additional source authorization/provenance:** ask whether JFA is in
      scope and obtain only the missing format/provenance evidence. Do not
      ask the owner to rediscover the CRS recorded in the inventory and
      `.prj`. Remaining distribution provenance is separate from the package
      import policy; the proposed numpy/scipy/shapely approach introduces no
      GPL-import defect by itself.

## Verification and limits

Only this review file was added. No source, recipe, job script, output, or
design document was changed; the already staged design was left intact.
No dependencies were installed. No mesh generation, survey-wide scan, figure
batch, FVCOM run, or batch submission was performed. Therefore there are **no
new production node/element counts, QA results or achieved timesteps** to
report. The original sounding counts, grid-spacing measurements, full JFA
coverage and quantization claims were not independently reproduced; they
must not be presented as measurements from this review.

Read-only inspection used `rg`, `sed`, `cat`, and `head` on the requested
code/contracts, the data inventory, and the two small JFA metadata/data
samples named above. Exploratory searches for `patch_sizing.py`, `io/fvcom.py`
and `fvcom*` module paths failed because those paths do not exist; the relevant
definitions were then located in `patch.py` and `io/fvcom_native.py`.
`git ls-files outputs` produced no entries.

The following small in-memory probe ran successfully with
`PYTHONDONTWRITEBYTECODE=1` and
`/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python -`.
It completed in under one second and wrote no test file:

```python
import numpy as np
from shapely.geometry import Point
from fvcom_mesh_tools.refine import limit_rfactor

area = Point(0, 0).buffer(300, quad_segs=64).area
print(area, np.sqrt(area / 9))
tri = np.array([[0, 1, 2], [1, 3, 2]])
free = [False, True, True, False]
for depths, bounds in [
    ([1., 2., 2., 4.], {}),
    ([1., 2., 2., 2.], {}),
    ([1., 3., 3., 1.], {"depth_min": 3., "depth_max": 300.}),
]:
    result, info = limit_rfactor(tri, depths, free, .2, **bounds)
    print(result, {k: info[k] for k in (
        "converged", "n_over_rmax_movable",
        "n_over_rmax_frozen_pair", "max_frozen_depth_change_m")})
print((30. - 3.) / (30. + 3.))
print((400. - 400. / 1.2) / (400. - 30.))
```

These probes substantiate findings 1, 2 and the circle-area part of finding 4.
The other findings are traced design contradictions or explicitly identified
counterexamples/risks, not alleged failures of an implementation that does
not yet exist. A full density remeasurement would require a batch job under
the repository's execution rules; none was created because this task permits
only the review file to change.
