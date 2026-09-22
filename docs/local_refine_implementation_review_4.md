# Fourth review, completed at b28d0d5

Reviewed 2026-09-22. Read the three prior implementation reviews and the current
design, especially §§0, 5, 6 and 9. This report does not reopen their acknowledged
obligations. No production mesh or implementation was changed.

**The original 15 probes pass, but several fixes remain incomplete.** The most
important regression is rejection of legitimate inherited failures when QA
truncates its evidence. I found no new demonstrated blocker for the delivered
Futtsu mesh. Findings below concern the stated general-input/API contracts;
P1 denotes a blocker for the identified use, not evidence that Futtsu is wrong.

## 1. P1 for large inherited-failure sets — the new attribution rule rejects unchanged meshes

**CLAIM:** Conservative handling of missing evidence is appropriate, but the
driver itself discards evidence and then blames the patch for its absence.

**EVIDENCE:** `420_local_refine.py:651,696` caps offenders at 10,000;
`patch.py:1297-1305` turns excess counts into introduced violations. The probe
`test_driver_cap_rejects_large_unchanged_base` runs real default QA on a
73×73-node rectangular grid with vertical spacing 20 m and horizontal spacing
100 m. All **10,368 triangles** have inherited C1 failures; **10,000** are named.
With every triangle retained, attribution still returns **368 unattributed**
violations in one entry. A two-triangle example with a cap of one reproduces
the same mechanism. This is a false positive without any changed triangle.

**WHY IT MATTERS:** Every seed can fail on an exterior it is forbidden to change.
The policy does not make every base unusable: fully enumerated inherited C1,
C4 and C5 failures still attribute correctly. It does make success depend on
report truncation, contrary to §0.

**SUGGESTED FIX:** Separate complete gate attribution from display truncation.
Return all offender identities, or calculate introduced counts before truncating.
Keep missing-evidence rejection for genuinely unplaceable checks; do not restore
the empty-list acceptance hole. Increasing a fixed cap merely moves this failure.

## 2. P1 for the public read/control/export workflow — OBC type disappears on copying

**CLAIM:** Direct read→export preserves type, but a routine bathymetry-control
copy silently restores type 1.

**EVIDENCE:** `fvcom_native.py:575` attaches an undeclared `mesh.obc_type`;
`apply_obc_depth_control`, line 296, uses `dataclasses.replace`, which drops it.
`test_obc_type_survives_depth_control` sets type **2**, applies control, and
exports: the file contains **1**. The writer's internal control path is safe:
it captures the type at lines 327-328 before making its copy. The current 414
and 420 callers also explicitly pass the original type; they avoid this defect.

**WHY IT MATTERS:** A supported composition of public helpers changes boundary
semantics without a warning. This is incomplete propagation of the new metadata,
not a demonstrated regression in the present production driver.

**SUGGESTED FIX:** Declare the metadata on the mesh dataclass and preserve it in
all mesh transformations; test read→control→export as well as read→export.

## 3. P2 — conflict reporting still does not evaluate the sizing field

**CLAIM:** `_region_size` is not the field implemented by `patch_sizing`.

**EVIDENCE:** `patch.py:886-904` uses `target/(1-d/width)` without ambient size;
lines 873-881 implement `target + (base-target)*d/width`. In
`test_conflicts_matches_actual_field`, a 5 m box has a 1,000 m transition and a
12 m box is 480 m from it. On the actual 100 m grid the joint field at the
second box's centre remains **12 m**, exactly its isolated field. Nevertheless
`region_conflicts` says it gets **9.6153846 m** because of the first box, with
`by_core_overlap=False`. The first region's actual ramp is above 12 throughout
that box. The original transition-swallowing example still passes.

**WHY IT MATTERS:** The newly advertised field-based cost warning can invent
refinement which does not occur. This does not alter the meshing field itself.

**SUGGESTED FIX:** Pass the ambient/isolated/joint field into the diagnostic and
sample or bound their actual difference over each core. Label minimum sampled
size and affected area distinctly; geometric core overlap is not transition area.

## 4. P2 — gradient magnitude regresses narrow-footprint diagnostics

**CLAIM:** Missing derivative support is reported as zero slope.

**EVIDENCE:** `patch.py:948-965` masks outside samples before `np.gradient`, then
requires both components to be finite. For box `[0,100]×[0,20]`, spacing 10,
and `h=100+x`, **nine interior samples** produce **max_slope=0**, although the
true slope is **1**. The previous axis-difference implementation measured 1
on this example. See `test_thin_gradation_has_no_false_zero`.

**WHY IT MATTERS:** A narrow fishery can receive a reassuring zero and zero
fraction above the reference simply because no y derivative is available.
Final mesh C4 remains gated, so this is a diagnostic regression, not a C4 bypass.

**SUGGESTED FIX:** Sample a halo and evaluate derivatives before masking, or
adapt spacing/stencils to geometry. Report insufficient support as unknown,
with the number of valid gradient samples, never as zero.

## 5. P2 — malformed element attribution still crashes or falsely exonerates

**CLAIM:** The non-integer-ID fix covers the node lookup branch, not elements.

**EVIDENCE:** `patch.py:1270-1285`; parametrized
`test_malformed_element_is_not_inherited`: element ID **"0"** raises TypeError,
**0.5** returns no introduced violation, and **-1** also returns none.
`kind='edge'` likewise trusts its `elements` list without validating IDs.

**WHY IT MATTERS:** The claimed conservative API behavior is false. Current
`run_qa` emits valid integer element IDs, so these are malformed-input API
failures, not an established ordinary driver bypass.

**SUGGESTED FIX:** Require non-Boolean integer indices in range before claiming
inheritance; validate every incident-element ID. Unplaceable records remain
introduced/unknown, with a diagnostic instead of an exception.

## 6. P2 — finite-depth protection misses non-finite grid coordinates

**CLAIM:** A NaN in the grid still defeats the coordinate agreement check.

**EVIDENCE:** `fvcom_native.py:408-412,527-532`; `test_grid_nan_rejected` changes
one grid X to NaN and leaves the depth file finite and otherwise matching.
`read_fvcom_case` returns a mesh containing NaN instead of refusing it.

**WHY IT MATTERS:** The reader has accepted invalid geometry and failed to prove
the files correspond. Later triangulation/QA may reject it; I did not establish
successful end-to-end generation with NaN coordinates.

**SUGGESTED FIX:** Validate finiteness of both coordinate arrays before subtraction,
and require a finite nonnegative coordinate tolerance.

## 7. P2 API validation — the limiter can certify undefined r-factors

**CLAIM:** `converged=True` is sound for the tested finite positive domain, but
not for all accepted inputs.

**EVIDENCE:** `refine.py:592,609-620`; parametrized
`test_limiter_invalid_depth_refused` supplies one triangle with all depths equal
to 0, -1, NaN or infinity. None raises; each reports convergence. Zero/NaN/infinity
produce undefined ratios which compare false against the threshold.

**WHY IT MATTERS:** An invalid bathymetry array can receive a success certificate.
The native reader rejects non-finite depths, and production depths are positive;
this is not a demonstrated defect of the production-base run.

**SUGGESTED FIX:** Validate finite strictly positive depths, valid connectivity,
mask shape and consistent finite positive optional bounds. Explicitly reject
non-finite final ratios rather than treating them as within the limit.

## Checks that hold, and remaining scope

| Claim | What I checked / conclusion |
| --- | --- |
| Empty/truncated/unplaceable QA failures | Original empty-list and node-pair probes pass. Complete retained-element, seam-edge and C5 fan attribution controls pass. Missing evidence now blocks, correctly; finding 1 is the avoidable truncation regression. |
| Frozen-pair map direction | **Correct.** `node_map` maps base→output; lines 478-481 invert it and compare output-edge endpoints in base IDs. The new probe executes these actual driver statements with map `[2,-1,0,3,1]`: output edge `[2,0]` maps to existing base `[0,2]`, while `[0,1]` maps to absent `[2,4]`. Exactly one is new. The original actual-driver rejection probe also passes. |
| Limiter soundness on supported depths | Finite positive examples at rmax **0, 0.2, 0.999** converge, preserve frozen depths exactly and stay within bounds; independently measured edge ratios meet the limit. Original infeasible movable-edge and new frozen-diagonal examples are correctly nonconverged. Finite iteration exhaustion remains an honest refusal, not proof of infeasibility. No new bounds failure found. |
| Ceiling composition | Nested/swallowed cores, polygon holes, zero-width regions and reversed region order pass. Joint sizing equals the pointwise minimum of isolated fields and clips against ambient. Priority is explicitly irrelevant under the revised ceiling contract; this is not the old ignored-priority finding. The diagnostic fails separately in finding 3. |
| OBC types | Type **3** survives read→default export including the writer's default depth control; mixed types are refused. Original type-2 control passes. Copying metadata fails separately in finding 2. |
| Depth-file finiteness | New NaN, +infinity and -infinity depth probes are refused. Original dep-coordinate NaN probe passes. Grid-coordinate gap is finding 6. |
| File geometry | Projected EPSG:32654 GeoPackage with two features: attribute selection, reprojection and 20 m buffer pass; buffered 100 m square area is within 150 m² of `10000+8000+pi*400`. Ambiguous selection/out-of-range index refuse. A hole outside its shell refuses before buffering, independently of the original bow-tie probe. No local CRS/axis-order defect found. |
| Buffer limits | The implementation explicitly uses a local equirectangular approximation (111000 m/degree and centroid cosine), not an exact geodesic metre buffer. No global/polar/dateline guarantee was established. For the local Tokyo Bay scope the exercised approximation is reasonable; do not advertise arbitrary geographic extent as tested. |
| Gradation rotation | Four differently oriented affine fields have unit gradient magnitude as expected. Original diagonal ramp passes. Narrow support is finding 4. |
| Physical ramp | DTE **1.7, 2.3, 5.0 s** produce the declared 86400 s timescale within half one internal step, the inevitable integer-rounding error. Read-only FVCOM source confirms `tanh(IINT/IRAMP)` (`../FVCOM/src/fvcom.F:829-830`). Original DTE=1.5 probe passes. No new ramp defect found. |

**Native format assumptions:** This reader supports the project's ASCII,
row-ordered, metre-coordinate grid/depth files; it is not a general FVCOM format
or CRS autodetector. Ignoring node-row labels is consistent with the local
FVCOM reader, which reads `J,X1,Y1` but stores `XG2(I),YG2(I)`
(`../FVCOM/src/mod_input.F:4676-4692`). I initially tested duplicate labels as a
rejection requirement, then withdrew it after checking that source; the final
probe confirms row-order behavior. Connectivity is range checked; the dep column
wins over any grid depth column. Extra trailing dep/OBC rows are ignored by
slicing; strict row-count validation would make the accepted subset clearer.
The largest-loop/one-OBC assumptions remain the documented restricted topology,
not newly discovered support for multiple components. No new valid-format
counterexample to the project's ordinary input was established.

**Controlled M2 comparison:** For a verified base/refined pair, the current 414
preparer does impose the same DTE, ISPLIT, physical ramp, dates, physics, sigma
levels, harmonic anchors, along-arc forcing and sponge construction. Identical
ordered OBC coordinates are required. Its explicit type argument avoids finding
2. This is a controlled comparison of the complete refinement operation,
including declared interpolation/r-factor bathymetry changes; it is not a
pure-connectivity experiment or a convergence study. One shared step removes a
between-case time-step confound but cannot quantify temporal error in the tiny
reported M2 differences. A paired smaller-step run would address that.

The preparer alone does **not** certify an arbitrary supplied pair:
`414_refine_m2_prep.py:85-92` compares OBC coordinates, not OBC types, depths,
NEXT_OBC neighborhoods or exterior identity. For reuse, validate those from the
refinement map/report and compare post-control OBC depths and types before
staging. Different types are preserved, not normalized, so identical coordinates
alone are insufficient. I did not rerun FVCOM or independently validate the
historical M2 numbers. The design acknowledges that those integrations predate
the ramp correction; the size of its effect is not measured by this review.

## Single most likely silent wrong-mesh outcome at HEAD

**An unresolved requested core that nevertheless passes structural and quality
checks.** This is the already-recorded achieved-resolution acceptance gap, not
a new finding or a reopening of the old priority policy. HEAD now reports
per-region edge counts/median/p90/max, but only **after accepting a seed**
(`420_local_refine.py:720-750`); no rejection uses them. A thin/small core with
no edge midpoints inside yields `n_edges=0` and null size statistics without
changing success. All checks may say correct while the requested fishery is
not represented at its ceiling resolution. Add geometry-aware core coverage
and achieved-size acceptance before seed selection, with an explicit tolerance
policy. Do not treat the sampled sizing field itself as evidence of delivered
resolution. This is a risk ranking from the acceptance code, not a newly run
bad production mesh. I found no replacement ordinary-QA hole as direct as the
original `offenders=[]` defect.

## Verification and artifacts

The executable companion is
[`local_refine_implementation_review_4_probes.py`](local_refine_implementation_review_4_probes.py).
It contains the failing assertions as ordinary pytest code plus the passing
controls; it does not invoke meshing or FVCOM. From repository root:

```bash
MPLCONFIGDIR=/tmp/lr4-mpl \
/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python -m pytest -q \
  docs/local_refine_implementation_review_4_probes.py /tmp/test_lr4_review.py --tb=short
# 13 failed, 34 passed: 13 new failing assertions, 19 new passing controls,
# and all 15 original probes passed. Five expected invalid-arithmetic warnings.
```

The identical companion was run from `/tmp/test_lr4_followup.py` in **5.23 s**.
Initial probe harness mistakes in AST extraction were corrected before this
result; they are not findings. `/tmp/test_lr4_review.py` was left unchanged.
Full-suite result is recorded below. Only this report and its companion are
repository additions; no dependencies, commits, branches, output meshes or
sibling files were changed. No generation, finishing, figure batch or FVCOM
integration was run, and no batch submission was attempted. Consequently no new
production mesh QA/counts/dt claims are made here.

```bash
MPLCONFIGDIR=/tmp/lr4-mpl \
/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python -m pytest -q
# 638 passed, 5 existing rasterio PendingDeprecationWarnings in 68.60 s.

/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/ruff check \
  src/fvcom_mesh_tools/patch.py src/fvcom_mesh_tools/refine.py \
  src/fvcom_mesh_tools/sizing.py src/fvcom_mesh_tools/io/fvcom_native.py \
  notebooks/420_local_refine.py notebooks/383_m2_case_prep.py \
  notebooks/414_refine_m2_prep.py tests/test_patch.py tests/test_refine.py \
  tests/test_local_refine_driver.py
# All checks passed. The adversarial companion is not part of the normal suite.
```
