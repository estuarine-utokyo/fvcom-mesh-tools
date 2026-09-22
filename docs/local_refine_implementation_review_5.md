# Fifth review: verification of the fourth-review fixes

Reviewed 2026-09-22, against `95538f5`, starting with
`git diff b28d0d5..95538f5`. Scope is the seven fixes and the resolution gate,
including their immediate callers. Only this report and its
[pytest companion](local_refine_implementation_review_5_probes.py) were changed.

**The fixes are substantive, but the resolution gate still accepts spatially
unresolved requests.** This is reproduced using the actual saved Futtsu mesh,
without changing a mesh file. There are also remaining diagnostic and API
defects. I did **not** demonstrate corruption of any of the three delivered
Futtsu meshes under their original recipes.

## Disposition in the original numbering

| Fix | Verdict and reproduction |
| --- | --- |
| 1. Complete offender identities | The cap fix works in `run_qa`. New tests exercise caps **0, 1, 100**, element/node/edge offenders, C4 and optional implied-dt gating. Attribution is independent of the cap. The fallback has a separate API regression: it strips incident elements. See finding E. |
| 2. OBC type field | Type 3 survives depth control and native export/read. The exact `replace` defect is fixed. Explicit mesh reconstruction still loses it in `compact_nodes`; see B. |
| 3. Shared conflict expression | Independent linear-ramp arithmetic agrees with `patch_sizing`. The formula fix is real; sampled affected area can still be grossly wrong. See C. |
| 4. Masked gradation stencil | The outside-value control passes, so the production-edge slope spike is addressed. A transverse gradient in a narrow channel still reads as measured zero; see D. |
| 5. Malformed IDs | Boolean, string, fractional, negative, missing and out-of-range scalar element/incident IDs all block. An unhashable malformed ID crashes before validation; see G. No ordinary `run_qa` bypass found. |
| 6. Grid NaN/tolerance | NaN and both infinities in either grid coordinate are refused; NaN, infinity and negative tolerance are refused. No remaining defect found within this fix. See `fvcom_native.py:412-415,525-540` and the nine new parametrized reader controls. |
| 7. Limiter positive depths | The input check is real. Positive used nodes with an unused dry node still work and preserve the unused zero. Bounds can subsequently create zero depths and a false convergence certificate; see F. |
| New resolution gate | Empty and uniformly coarse requests reject; the actual 5% threshold works. Edge-count statistics do not prove spatial coverage. See A. |

File names below are relative to the repository. `patch.py`, `qa.py`,
`refine.py`, and `mesh_clean.py` mean files under `src/fvcom_mesh_tools/`;
`fort14.py` and `fvcom_native.py` are in its `io/` directory.

## A. P1 acceptance-contract gap: an unresolved lobe is hidden by a fine lobe

**Evidence:** `notebooks/420_local_refine.py:369-393` selects edge midpoints,
requires at least one inside the entire region, and gates only the median.
Lines 715-739 use that result in seed acceptance. Lines 791-795 repeat the
same test on the written mesh. Neither test checks coverage within the region.

**Reproduction on a real mesh:**
`test_defect_real_mesh_unsampled_lobe_is_not_hidden_by_fine_lobe` reads
`outputs/v4_futtsu_nori/TokyoBay_grd_futtsu_nori.14`. It declares a new region
consisting of the existing 300 m-radius core, a second 300 m-radius disc
5 km north, and a 10 m-wide connecting corridor. GEOS verifies that this is
one valid Polygon and that **0 m² lies outside the mesh**. No coordinates,
connectivity, depths, recipes, or output files are modified.

The polygon is **609,430.54 m²**. Its remote lobe is **282,714.95 m²**, over
46% of the region, and contains **zero edge midpoints**. That lobe correctly
fails when assessed alone. Joined to the fine lobe, it passes:

| Combined request | Measured value |
| --- | ---: |
| Target | 30 m |
| Counted edges | 1,171 |
| Median | 28.5277 m |
| p90 | 30.3081 m |
| Maximum | 219.3272 m |
| Missed regions | none |

An ordinary enlarged circle also passes:
`test_defect_real_mesh_expanded_request_is_not_delivered` buffers the original
core by 300 m. Its **2,313** counted edges have median **29.5365 m**, p90
**61.8272 m**, and maximum **145.7898 m**, against a 30 m request.

These execute the **actual driver function and actual tolerance assignment**
extracted with AST, not a reimplementation or a stub of the statistic. They
show a faulty acceptance decision on existing production geometry. They do
not show that rerunning DistMesh for these changed requests would generate
these candidates, nor that the original circle request was wrong.

The 1.05 factor (`420_local_refine.py:347-354`) is a modest empirical allowance,
not a coverage guarantee. The added threshold control accepts 31.49 m and
rejects 31.51 m for a 30 m request. Both bad examples above pass even at 1.00;
tightening the tolerance cannot solve this. Counting edges preferentially
weights finely meshed water. Even p90 misses the unsampled lobe.

**Repair direction:** retain a calibrated tolerance, but test spatial coverage
independently of edge density: subdivide/sample the requested water by area,
locate its containing cells, and enforce a documented local resolution and
coverage policy. A single global maximum would be a different, much stricter
contract; the delivered meshes already have isolated edges above target.

## B. P1 for the affected API composition: compaction still changes OBC type

**Evidence:** `fort14.py:55` now defaults the field to 1, correctly making
`dataclasses.replace` preserve it. But `mesh_clean.py:192-199` constructs a
new `Fort14Mesh` without the field. `fvcom_native.py:327-336` then exports
that new default. The new class documentation's assertion that every
transformation uses `replace` (`fort14.py:42-45`) is not true of this caller.

**Reproduction:** `test_defect_compaction_preserves_type` adds one unused node
to a type-3 mesh, compacts it, checks that the open-boundary node list is
preserved, and exports. The OBC file declares **1**, not **3**. This is a
legitimate cleanup operation; it does not require malformed connectivity or
an unsupported boundary type. The companion control verifies that direct
depth-control → export → read now preserves type 3.

This remains an **API metadata-propagation defect**, not a demonstrated
delivered-Futtsu corruption. The 420 native exporter explicitly uses the base
type (`420_local_refine.py:814-816`), so it avoids this composition. Carry the
field through explicit reconstructions as well as `replace`. The default 1
cannot distinguish an intentional type 1 from lost metadata.

## C. P2 diagnostic: a true 1% affected area becomes 0% or 100%

**Evidence:** `patch.py:923-935` sets sample spacing from area/400, lays out
an axis-aligned lattice, excludes boundary points, and adds one representative
point. `patch.py:1120-1138` converts the fraction of hit samples directly to
area, without an uncertainty qualifier. The driver presents this as percent
of area (`420_local_refine.py:238-242`).

**Reproduction:** both cases of `test_defect_thin_conflict_area` use a
10,000 × 10 m region at 30 m and a 100 × 10 m overlapping region at 5 m,
constant ambient 300 m, and zero transitions. True affected area is exactly
**1,000 m² / 100,000 m² = 1%**. Spacing is 15.81 m, so no lattice row lies
inside the thin region; only its representative point survives.

* Put the fine rectangle at x=0: `finer_than_declared` omits the coarse region.
* Put it at x=4,950: reported fraction is **100%**, area **100,000 m²**,
  `n_samples=1`.

The exact core-overlap branch, called without `base_size`, returns 1% in both
cases. Supplying the real field thus makes a previously exact core-area answer
worse, even though the pointwise formula is now correct. `overlapping_pairs`
still contains the correct 1,000 m², making the returned report contradictory.

This changes warnings, not the meshing field or QA acceptance. It is not
evidence of a wrong delivered Futtsu mesh. Use exact area where available;
otherwise use geometry-adapted, area-weighted sampling with an explicit
estimate/coverage qualification. A representative point is not an area sample
with weight equal to the entire polygon.

## D. P2 diagnostic: missing transverse support is still printed as flat

**Evidence:** `patch.py:1025-1044` accepts support on either axis, substitutes
zero for a missing component, and reports `measured=True`. The partial count
is returned, but `420_local_refine.py:305-312` only checks `measured` and prints
the numerical maximum/p99/fraction without qualifying partial support.

**Reproduction:** `test_defect_transverse_slope_not_reported_as_measured_flat`
uses a **1,000 × 40 m channel** and `h=100+y`. At the **default 25 m spacing**,
there are 39 interior samples, all partial. The result is `measured=True`,
**max=0, p99=0, fraction above reference=0**, while the gradient magnitude is
exactly **1**. At 5 m spacing, the same function measures 1.

`n_partial_samples` counts missing-axis support, **not one-sided stencils**:
the booleans returned at `patch.py:1021` only distinguish zero from nonzero
neighbour count. One-sided derivatives with both axes supported are exact
for affine fields; they need not be suspect merely because they are one-sided.
The companion outside-value test has zero partial samples, including its
one-sided boundary stencils, and measures 0.2 despite outside values of 10,000.

Thus the fix avoids the reported outside-mesh spike, but its lower bound can
still be arbitrarily far below the transverse truth. This is a realistic
narrow-channel geometry, not proof that the actual Futtsu hole has this defect.
Finished-mesh C4 remains independent. Distinguish partial bounds from measured
gradient magnitudes in both return data and driver text, or refine the lattice
until the missing direction is supported.

I also reran the current field diagnostic on **the actual Futtsu patch
footprint**, reconstructed as the union of the saved mesh's non-retained
triangles, using the actual native base and recorded region/transition width.
`test_saved_futtsu_field_on_actual_patch_footprint` measures:

| Spacing | Samples | Partial | Max slope | p99 |
| --- | ---: | ---: | ---: | ---: |
| 35.73596 m (default) | 27,035 | 7 | 0.381441 | 0.362128 |
| 17.86798 m | 108,084 | 5 | 0.384096 | 0.362446 |

Both have zero fraction above 0.414214. This supports the real-case fix for
the previously reported 24.3 spike; the narrow-channel counterexample does
not overturn that successful correction.

## E. P2 API regression: fallback identities discard incident elements

**Evidence:** `qa.py:161-165` copies only `kind` and `id` from a plain display
list. `patch.py:1397-1408` needs `elements` to place an edge offender.

**Reproduction:** `test_defect_fallback_keeps_incident_elements` constructs a
`QACheck` with one edge offender, ID `[0,2]`, incident elements `[0,1]`, and
both elements retained. Its derived identity loses `elements`; attribution
returns an introduced violation instead of `[]`. This is an API false positive
introduced by the fallback, even with **no truncation**.

I found **no `run_qa` gate that silently falls back to a truncated display
list and accepts the missing violations**. Helper-built lists carry complete
identities (`qa.py:441-480`); C4 explicitly retains incident elements
(`qa.py:865-884`). The fallback defect above concerns manually constructed
checks, not current C4 production attribution.

Two qualifications to “every offender's identity” remain:

* OBC duplicate occurrences are counted but never identified
  (`qa.py:763,771-784`). The new duplicate-walk control has one duplicate,
  no identities, and one **unattributed** violation. It fails closed; this
  is incomplete evidence, not a restored acceptance hole.
* `_offender_key` omits `segment` (`patch.py:1438-1439`). Two OBC segments
  containing the same invalid pair collide. The collision control executes
  `run_qa` on segments `[1,3]`, `[1,3]`: identities name segments **0,1**, but
  the decorated introduced records name **0,0**. Both violations still block,
  since attribution occurs before decoration (`patch.py:1392-1414`). This
  matters for report accuracy, not mesh acceptance in this reproduction.

Retain `elements` in fallback identities; include segment in decoration keys.
The complete-identity architecture itself is the right fix for the cap.

## F. P2 API validation: valid input depths can become certified zeros

**Evidence:** `refine.py:587-592` validates input depths only. Bounds at
lines 597-598 are unchecked; clipping at 615-620 may invalidate the depths.
Final comparison at 621-635 still treats NaN ratios as not over the limit.

**Reproduction:**
`test_defect_limiter_cannot_create_invalid_depth_and_certify_it` passes one
triangle, finite positive depths `[1,10,1]`, all nodes movable, rmax=0.2,
and `depth_max=0`. It returns **[0,0,0]**, **converged=True** after two rounds,
and **max_r_on_new_edges=NaN**. Reject the invalid bound or reject the final
invalid state; success is not valid. This is an invalid-parameter API defect,
not a failure for the production parameters.

The only non-test caller found by `rg -n 'limit_rfactor\(' --glob '*.py'` is
`420_local_refine.py:525-527`; it supplies the base's min/max depths. The saved
Futtsu artifacts and recipe use positive bathymetry. I found no production
caller passing a legitimate connected dry/mid-construction node and therefore
no demonstrated regression from the new positivity requirement. The new unused
dry-node test still converges and preserves the zero because validation is
restricted to connectivity-used nodes. Supporting connected dry bathymetry
would require an explicit wet/dry convention; this routine currently has none.

## G. P2 malformed-input API: hashing precedes ID validation

**Evidence:** `patch.py:1391` hashes display keys before `_element_id` at
1396; `_ident_key` at 1442-1448 returns a dict unchanged.

**Reproduction:** `test_defect_unhashable_malformed_id_is_conservative` uses
element ID `{"index": 0}`. It raises **TypeError: unhashable type: 'dict'**
instead of returning an unplaceable violation. This is strictly a malformed
API-input failure. `run_qa` does not emit this ID shape, and no bad mesh was
accepted. Scalar validation is fixed; either validate before keying or define
and enforce a narrower accepted record schema.

## Actual saved cases and limits of verification

`test_delivered_artifacts_measured_read_only` independently rereads the three
saved meshes, runs current QA with **max_offenders=0**, attributes against
each recorded retained-element count, executes the current resolution function,
and recomputes minimum-altitude implied dt. These are new measurements of
existing artifacts, **not new generation runs**:

| Saved output | Nodes | Elements | QA gates passed | Introduced | Implied dt (s) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `v4_futtsu_nori` | 4,413 | 8,039 | 20/21 | 0 | 2.4883003915 |
| `v4_futtsu_nori_polygon` | 5,991 | 10,755 | 21/21 | 0 | 2.2385176798 |
| `v4_futtsu_two_beds` | 6,179 | 11,131 | 21/21 | 0 | 2.9608912640 |

Original-request medians/p90/max (metres), all accepted:

* Circle: **28.527 / 30.189 / 47.388**, target 30, 1,166 edges.
* Polygon: **28.343 / 30.226 / 52.276**, target 30, 1,209 edges.
* Two-beds north: **29.451 / 30.578 / 37.870**, target 30, 1,155 edges.
* Two-beds east: **29.489 / 42.179 / 52.921**, target 45, 842 edges.
* Channel: **55.322 / 71.907 / 90.196**, target 60, 131 edges.

These numbers support a statistical resolution policy, not a literal maximum
edge-length ceiling. In particular, changing the gate to p90 ≤ 1.05×target
would reject the saved channel case. A revised policy needs measured spatial
coverage and a stated allowance for the normal meshing spread.

Beyond the requested fixes, the driver now stops if gradation is unmeasurable
(`420_local_refine.py:305-309`), and prioritizes number of missed regions before
number of introduced QA failures when selecting its best candidate
(`420_local_refine.py:735-747`). These are conservative behavior changes.
The gate operates on a serialized candidate (`420_local_refine.py:704-715`):
“fails instead of writing” should not be interpreted as “no provisional .14
file exists after failure.” I did not regenerate a failing case to test its
remaining filesystem artifacts, so this is a control-flow observation only.

No new mesh-generation or FVCOM jobs were run. Repository instructions require
batch execution for meshing, block `qsub`/`qstat`, and this task permits edits
only to the two review files. No job script was added, no dependency installed,
and no production output modified. The 32 old probes and 654 repository tests
were supplied as baseline information, not used as evidence that these fixes
are complete and not rerun here.

## Reproduction commands and results

```bash
source /octfs/work/G16445/v61021/miniforge3/etc/profile.d/conda.sh
conda activate oceanmesh-bench
MPLCONFIGDIR=/tmp/lr5-mpl XDG_CACHE_HOME=/tmp/lr5-cache \
  pytest -q -s docs/local_refine_implementation_review_5_probes.py --tb=short
ruff check docs/local_refine_implementation_review_5_probes.py
```

The companion deliberately uses ordinary failing assertions for the defects,
not xfails or assertions that bless the incorrect result. Artifact tests skip
explicitly if the local saved outputs are unavailable. The first run before
adding the lobe and tolerance controls was **8 failed, 29 passed in 36.22 s**;
it also emitted font-cache permission diagnostics. An initial Ruff import-order
failure was corrected in the allowed probe file.

The subsequent full run was **9 failed, 31 passed in 26.10 s**, with no skipped
artifact tests. The final added real-field control was then run separately:

```bash
MPLCONFIGDIR=/tmp/lr5-mpl XDG_CACHE_HOME=/tmp/lr5-cache \
  pytest -q -s docs/local_refine_implementation_review_5_probes.py \
  -k saved_futtsu_field --tb=short
# 1 passed, 40 deselected in 2.05 s
```

Together these exercise all **41** final probes: **9 defect failures and 32
passing controls/measurements**. The native base depth range in the last check
was **3–300 m**. Final Ruff check: **All checks passed!**

**Most likely current silent wrong-mesh outcome:** a requested region contains
a well-refined area and a substantial unresolved area, yet the fine area's
edge count dominates the median and the run reports that the resolution was
delivered. The current gate actually accepts that decision on the saved Futtsu
geometry above. I have not shown that the current mesher generates such a
candidate for these changed requests, or that the delivered original recipes
are wrong; the demonstrated defect is that the acceptance gate cannot detect it.
