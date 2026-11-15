# Review: new-region generalisation of fmesh-pipeline (2026-11-16)

Reviewer: gpt-6.1-sol (Codex), run by the coordinator. Target: commits `1988f22..14a7dd8`.
Sources: `scratch/reviews/newregion/r1_prompt.md`, `r1_final.md`.

## Round 1

### Prompt

````markdown
# Review request, round 1: new-region generalisation of fmesh-pipeline (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT modify files. You may run read-only
commands, python in memory, mocks and small synthetic meshes (keep each computation under a few minutes). Answer
in English as Markdown.

## Goal
Correctness and robustness of the target below under realistic use. Report every defect you can substantiate in
the target, of any severity.

## Target (concrete)
Commit range `1988f22..14a7dd8` (10 commits, e557279 ... 14a7dd8). Files:
`src/fvcom_mesh_tools/cli/pipeline.py`, `finishing.py`, `mesh_engine/polygon.py`, `obc_tools.py`,
`prep/shoreline.py`, `site_session.py`, `qa.py`, `autofinish/directives.py`, `autofinish/plan.py`,
`cli/meshplot.py`, `cli/prepshoreline.py`, `tests/test_new_region.py`, USER_GUIDE section 16, CHANGELOG.

The ONE question: can a user build a mesh for a NEW region (any UTM zone, a supplied domain polygon, a size
raster/regions, protected narrow channels, an open-boundary line) from one recipe with these changes, getting
either a correct mesh or a loud failure -- never a silently wrong mesh (wrong CRS, dropped protected water,
lost/altered open boundary, valence/lone-node passes that break conformity or the boundary) -- while the existing
Tokyo recipe (`tokyo_bay_v6`, no new keys) behaves exactly as before?

Specific points to check:
- CRS resolution order (`crs` > legacy `finish.utm_epsg` > UTM zone of prep.bbox/domain); no remaining silent
  EPSG:32654 default on any path the pipeline uses with a non-Tokyo recipe.
- Channel protection (`protect_polygons`, `on_dropped_protected`) really prevents finishing/siteops/polish from
  deleting elements inside protected water, and failure is loud.
- OBC carried through finish / obcfinal; `freeze_obc_band`, `obc.snap`, `finishing.auto.obc_lock: none` with
  `obc_perp_max_deg` guard; band restore checking C4.
- `max_valence` / `remove_lone_nodes` passes keep a valid conforming triangulation and do not move/remove OBC or
  protected nodes.
- `mesh_engine/polygon.py`: size field units (metres), `expected_drop` / `max_dropped_component_km2` logic,
  reproducibility with `seed`.
- `depths.product` = `external_dep_file` matching (`match_tol_m`) and `none`.
- Defaults: with none of the new keys, every code path is identical to 1988f22 (evidence of the author's md5
  regression: `scratch/regress_20261116/` -- base1b vs head*/new*; read, do not rerun the full Tokyo chain).

Real inputs: the VR-FVCOM recipe `$HOME/Github/VR-FVCOM/recipes/mesh/v3.yaml` (its `pipeline:` block, read with
`--section pipeline`), its outputs in `$HOME/Github/VR-FVCOM/scratch/mesh/v3/` (pipeline/ subdir, QA json), and
the Tokyo recipe `tokyo_bay_v6` in this repository.

## Out of scope -- do NOT report or investigate
- code the target only calls, and code that was not modified in the range;
- malformed or contrived input files, other platforms, injected faults, rare corner cases the real inputs do not
  produce (e.g. self-intersecting domain polygons, CRS strings that are not EPSG codes, empty recipes);
- performance, style and naming (unless it makes a realistic misuse likely);
- the Tokyo notebooks (`plotting.py` 32654 default kept deliberately for them; `fmesh-plot-mesh` now warns).

## What was done
Read `git log -p 1988f22..14a7dd8`. Intent: generalise the Tokyo-specific pipeline so a new region (first user:
Volta River mouth, Ghana, EPSG:32631, ~60-100 m channels in a 87 km river, lagoons, shelf to ~shelf edge) is
built fully automatically. VR v3 result (job 7591109): fmesh-mesh-qa 23/23 gates pass, NP 27,521 / NE 48,989,
1 component, 82 OBC nodes, min channel edge 59.9 m. Tokyo: md5-identical fort.14 at raw/finished/obc/siteops/
polished/final versus the old code plus the ledger fix (e557279), because the unpatched baseline crashes in
`finishing` on numpy int64.

## Please
1. A fresh audit of the Target above (not of code outside it).

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a realistic path
- minor: a clear robustness/clarity defect on a realistic path (an unusual or injected input does NOT qualify)
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or evidence, and a concrete fix. Then
`## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
````

### Answer (verbatim)

````markdown
Read-only review completed; no files changed. The stored Tokyo `base1b`, `head1`, `head2`, `new1`, and `new2` meshes have identical md5 values at every mesh stage. The stored Volta final mesh has NP 27,521 / NE 48,989 and 23/23 QA gates passing.

Verification used `git log -p`, source inspection, artifact comparisons, and Python `-B -c` in-memory reproductions. Four existing tests requiring no filesystem fixtures passed. Full generation was not rerun; importing OceanMesh encountered Matplotlib’s unwritable-cache error, so sizing probes used mocked generator imports.

1. **Major — Protection excludes nodes exactly on the protected shoreline.**  
   **Location:** `src/fvcom_mesh_tools/finishing.py:353`, `:381`.

   `contains_xy` excludes polygon boundaries. A narrow channel meshed with nodes on its two banks can therefore have no protected nodes, although every triangle lies inside protected water.

   **Reproduction:** A 60 m wide, 1,800 m long strip with 18 triangles, protection equal to its water polygon, and `on_dropped_protected="fail"` finished with **one triangle remaining**. Seventeen caps were deleted without an exception or component-drop record. Welding and optimization were mocked as no-ops.

   **Fix:** Classify protection using triangle overlap, including polygon boundaries and an appropriate geometric tolerance, rather than strict containment of vertices.

2. **Major — Lone-node removal ignores protection and can erase an inlet and its OBC.**  
   **Location:** `src/fvcom_mesh_tools/finishing.py:466`, `:478`; analogous unrestricted ear stripping in `mesh_engine/polygon.py:437`.

   `remove_lone_nodes` accepts `protect_polygons` but never uses it. Repeated ear deletion can consume an entire narrow arm; preserving connectivity at each deletion does not preserve protected water.

   **Reproduction:** On a fully protected 900 × 100 m strip, the default pass deleted **50,000 m²**, retreating its end from x=900 to x=400. A 300 × 100 m strip was deleted completely, including its open boundary, without raising.

   **Fix:** Protect affected elements and nodes before every deletion. Apply equivalent protection during polygon boundary repair, and fail when resolving a lone-node gate would require deleting protected water or OBC geometry.

3. **Major — R4 repair converts supplied open-boundary segments into land without requiring subsequent restoration.**  
   **Location:** `src/fvcom_mesh_tools/mesh_engine/polygon.py:473`, `:481`, `:509`; `cli/pipeline.py:159`.

   R4 repair advances the open-string endpoints and puts the excluded arc into the land string. Recipe validation permits `domain` → `build` → `qa` without an OBC restoration stage.

   **Evidence:** In the stored Volta raw mesh, **13 land-labelled edges totaling 8,614.8 m lie along the supplied OBC line**. Its open string covers arc positions 2,884.6–176,734.6 m of a 183,697.2 m line. The full Volta chain repairs this later; a permitted shorter recipe retains artificial walls across that open boundary.

   **Fix:** Preserve the supplied OBC extent, or fail when trimming occurs unless a required later restoration stage validates the complete line.

4. **Major — Freezing the OBC band also changes bathymetry at newly created nodes.**  
   **Location:** `src/fvcom_mesh_tools/finishing.py:263`, `:295`.

   Depths are replaced with random identity tags. For split nodes, restoration discards the optimizer’s interpolated depth and substitutes the nearest original node’s depth. This affects new nodes throughout the mesh, including outside the OBC band.

   **Reproduction:** Splitting an edge whose endpoint depths were 1 and 101 m should retain the interpolated midpoint depth of **51 m**. `_tag_obc_band` followed by `_restore_obc_band` restored it as **1 m**. `depths.product: none` leaves this alteration in place.

   **Fix:** Carry identity separately from bathymetry, or restore new-node depths using interpolation over the original triangulation.

5. **Major — `expected_drop` measures element count rather than water area.**  
   **Location:** `src/fvcom_mesh_tools/mesh_engine/polygon.py:344`.

   The unweighted centroid-membership mean lets finely meshed expected water authorize deletion of a much larger, coarsely meshed unexpected body.

   **Reproduction:** The actual component-classification routine classified a conforming **1.1 km²** component as `expected_drop=True` when the expected polygon covered only **0.1 km² (9.09%)**. That small region contained 40 of its 53 elements. As an unanchored component, it bypasses the 0.5 km² failure threshold.

   **Fix:** Calculate the fraction of component area covered by `expected_drop`, preferably using triangle–polygon intersection areas.

6. **Major — Valid raster coverage still produces a coarse sizing strip at its edges.**  
   **Location:** `src/fvcom_mesh_tools/mesh_engine/polygon.py:251`.

   Coverage validation accepts the raster’s full pixel footprint, but interpolation covers only pixel centres. The outer half-cell strip receives global `hmax`, overriding the supplied local size.

   **Reproduction:** A raster covering `[0,1000]²`, with centres at 50–950 m, western-column sizes of 100 m and other sizes of 400 m, passes coverage validation. The generator receives **100 m at x=50 but 400 m at x=0**, within the same western pixel footprint.

   **Fix:** Extend edge values through the pixel footprint, for example by clamping interpolation coordinates to the centre lattice.

7. **Major — Site fallback can accept a newly failing gate despite its acceptance contract.**  
   **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:888`, `:930`.

   A finite penalty of 1,000 does not enforce “no other gate starts failing.”

   **Reproduction:** `_gate_violations` scores 1,002 C1 offenders as **1,002**, while one newly failing non-tradeable gate scores **1,001**. Consequently, `s1 < s0` accepts such a result. On a sufficiently large mesh, reducing angle offenders can therefore pay for a new OBC-perpendicularity failure.

   **Fix:** Explicitly reject newly failing non-tradeable gates before comparing offender counts. Apply the configured OBC perpendicularity guard to fallback patches too.

8. **Major — The depth stage can discard the latest finishing stages.**  
   **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:1036`.

   Source selection ignores `polish_mesh`, `siteops_mesh`, and `finished_mesh`. A valid recipe containing `obc`, subsequent edits, and `depths`, but no `obcfinal` or `finishing`, selects the earlier OBC mesh.

   **Reproduction:** Using the stored Volta artifacts without the two optional finalization stages, a depth file matching the latest polished mesh was rejected: **27,521 nodes versus the selected OBC mesh’s 27,518**. `product: none` likewise reported retaining `vr_v3_obc.14`. When edits change connectivity without changing coordinates, external matching can succeed while silently discarding those edits.

   **Fix:** Select and propagate the latest produced mesh consistently, including for `product: none`.

9. **Major — Merging projected water files silently filters them with a geographic bbox.**  
   **Location:** `src/fvcom_mesh_tools/prep/shoreline.py:68`.

   `gpd.read_file(..., bbox=tuple(bbox))` applies longitude/latitude numbers in the source file’s CRS. Reprojection happens only after filtering.

   **Reproduction:** A valid in-memory EPSG:32631 water polygon at easting 250,000 m / northing 640,000 m reads normally without a filter. With the geographic bbox `(0,5,2,7)`, it returns zero features. Passing two such sources through `merge_water_files` produces an empty merged dataset.

   **Fix:** Transform the bbox into each source CRS before reading, or filter after reprojection to EPSG:4326.

10. **Minor — Prep does not consistently honor the resolved recipe CRS.**  
    **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:213`, `:262`.

    `_stage_prep` never forwards `--utm-epsg`; its perpendicular endpoint-extension call also omits the resolved CRS. Opening, skeletonization, and endpoint extension can therefore choose an automatic zone while later stages use the explicit recipe zone.

    **Reproduction:** With `crs: 32630` and a bbox spanning Ghana’s UTM-zone boundary, `_recipe_epsg` returns 32630, but the captured prep arguments contain no CRS override.

    **Fix:** Thread the resolved CRS into these prep operations while preserving the measured Tokyo behavior.

11. **Major — Persisted raw-mesh CRS survives a build-engine change and suppresses necessary projection.**  
    **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:575`, `:590`.

    Polygon builds persist `raw_mesh_crs`. Other build paths do not replace it, and provenance artifacts are merged on subsequent runs. Switching an existing output directory from `oceanmesh-polygon` to the geographic engine leaves the old projected-CRS marker attached to the new longitude/latitude mesh.

    **Reproduction:** After that sequence, `_stage_finish` passed coordinates such as **(0.5, 5.8)** unchanged to finishing while declaring **EPSG:32631**.

    **Fix:** Return an explicit source CRS from every build path and overwrite it on every build; transform according to that recorded source CRS.

## Verdict

VERDICT: FAIL (0 blocker, 10 major, 1 minor, 0 nit)
````

### Triage

Fixes: commit `32d211b`. Each finding has a regression test in `tests/test_new_region.py` (`test_r1_<N>_*`). All 12
tests fail on the reviewed code: they were run with `PYTHONPATH` set to a `git archive` of 14a7dd8. All 12 pass
on 32d211b. Behaviour of the reviewed code for #1 and #2:
`scratch/reviews/newregion/r1_evidence/old_behaviour.txt`.

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes. Reviewed code: a 30 m bank-to-bank strip whose protection is its own water polygon has 0 of 38 nodes protected, and cap deletion is allowed on 35 of its 36 elements. Test `test_r1_1` runs the full finish with the optimizer mocked. | yes | 32d211b: `protected_node_mask` (inside or within 1 m of the outline), used by finish, siteops and obc aftercare |
| 2 | major | yes. Reviewed code: on a protected 900 x 100 m strip, `remove_lone_nodes` deleted 10 elements (50,000 m2) and the end moved from x=900 to 400. On a 300 m strip with an open boundary it deleted all 6 elements and the open boundary. | yes | 32d211b: no cap deletion at all inside protected water; `remove_lone_nodes` never deletes an element with a protected or open-boundary node, and fails (policy `fail`) or reports the lone nodes left. Polygon engine: the boundary repair cannot leave lone nodes (FVCOM), so it reports the protected area it deletes and fails above `max_protected_ear_area_km2` (optional). VR v3 measured 2.37 km2. The region's coverage guard (all selected polygons >= 50 %) still applies. |
| 3 | major | yes. VR v3 raw mesh: the open string covers 2,885-176,735 m of the 183,697 m line, so 9,847 m are land. After the obc stage the final string covers 73-183,697 m (max offset 0.0 m). | yes | 32d211b: a polygon build with an open-boundary line requires obc/obcfinal (`_check_recipe`). The raw mesh records `raw_open_cover`. |
| 4 | major | yes. Reviewed code: an edge split between depths 1 and 101 m restored the new node as 1.0 m (`test_r1_4`). | yes | 32d211b: new-node depths are interpolated linearly on the original triangulation (51 m); the nearest node is used only outside it |
| 5 | major | yes. `test_r1_5`: a 2.25 km2 pond whose expected-drop polygon covers 4 % of its area (but most of its elements) was accepted by the reviewed code; it does not raise "dropped". | yes | 32d211b: area fraction from triangle-polygon intersection areas |
| 6 | major | yes. The reviewed code has no clamping. `test_r1_6`: on the reviewed code it fails on the missing helper; the behaviour was confirmed by the reviewer (100 m at x=50, hmax 400 m at x=0). | yes | 32d211b: `size_interpolator` clamps to the outermost cell centres |
| 7 | major | yes. `test_r1_7`, with mocked gate counts 1002 C1 -> 1000 C1 + 1 perpendicularity: the reviewed code kept the result. | yes | 32d211b: a newly failing or worse non-C1/C2/C4 gate refuses the fallback before totals are compared. Fallback patches go through `_execute_patches_guarded` (perpendicularity guard). |
| 8 | major | yes. `test_r1_8`: the reviewed code chose the obc mesh over the later polished mesh ("depths file has 9 nodes, the mesh 6"). | yes | 32d211b: `latest_mesh` artifact (the last mesh a stage wrote); the depth stage uses it and propagates it as `final_mesh`, including for `none` |
| 9 | major | yes. `test_r1_9`: an EPSG:32631 file filtered with lon/lat (0, 5, 2, 7) returned 0 features on the reviewed code. | yes | 32d211b: the bbox is transformed into each source CRS, then intersected in lon/lat |
| 10 | minor | yes. `test_r1_10`: no `--utm-epsg` in the prep argv with `crs: 32630`. | yes | 32d211b: an explicit `crs` is forwarded to prep and to the endpoint extension. Without `crs` the behaviour is unchanged (Tokyo). |
| 11 | major | yes. `test_r1_11`: the reviewed code's oceanmesh build returns no `raw_mesh_crs`, so a stale EPSG:32631 marker survives. | yes | 32d211b: every build records `raw_mesh_crs` (EPSG:4326 for the lon/lat engines); finish transforms from it |

No rebuttals.

Checks after the fixes:
- `pytest -q`: 1334 passed, 9 skipped. `ruff`: clean.
- Tokyo `tokyo_bay_v6` at 32d211b (job 7591382): raw, finished, obc, siteops, polished and final fort.14 are
  md5-identical to base1b. QA still has 5 gate failures (unchanged).
- VR v3 (VR-FVCOM job 7591383, `STAGES=domain:size:pipeline:qa:depth`): fmesh-mesh-qa 23/23, NP 27,521 /
  NE 48,989, one component, 82 open-boundary nodes, coverage guard passed, depth sha256 check passed.

## Round 2

### Prompt

````markdown
# Review request, round 2: new-region generalisation of fmesh-pipeline (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT modify files. You may run read-only
commands, python in memory, mocks and small synthetic meshes (keep each computation under a few minutes). Answer
in English as Markdown.

## Goal
Correctness and robustness of the target below under realistic use. Report every defect you can substantiate in
the target, of any severity.

## Target (concrete)
Commit range `1988f22..32d211b` (11 commits, e557279 ... 32d211b; 32d211b holds the round-1 fixes). Files:
`src/fvcom_mesh_tools/cli/pipeline.py`, `finishing.py`, `mesh_engine/polygon.py`, `obc_tools.py`,
`prep/shoreline.py`, `site_session.py`, `qa.py`, `autofinish/directives.py`, `autofinish/plan.py`,
`cli/meshplot.py`, `cli/prepshoreline.py`, `tests/test_new_region.py`, USER_GUIDE section 16, CHANGELOG.

The ONE question: can a user build a mesh for a NEW region (any UTM zone, a supplied domain polygon, a size
raster/regions, protected narrow channels, an open-boundary line) from one recipe with these changes, getting
either a correct mesh or a loud failure -- never a silently wrong mesh (wrong CRS, dropped protected water,
lost/altered open boundary, valence/lone-node passes that break conformity or the boundary) -- while the existing
Tokyo recipe (`tokyo_bay_v6`, no new keys) behaves exactly as before?

Specific points to check:
- CRS resolution order (`crs` > legacy `finish.utm_epsg` > UTM zone of prep.bbox/domain); no remaining silent
  EPSG:32654 default on any path the pipeline uses with a non-Tokyo recipe.
- Channel protection (`protect_polygons`, `on_dropped_protected`) really prevents finishing/siteops/polish from
  deleting elements inside protected water, and failure is loud.
- OBC carried through finish / obcfinal; `freeze_obc_band`, `obc.snap`, `finishing.auto.obc_lock: none` with
  `obc_perp_max_deg` guard; band restore checking C4.
- `max_valence` / `remove_lone_nodes` passes keep a valid conforming triangulation and do not move/remove OBC or
  protected nodes.
- `mesh_engine/polygon.py`: size field units (metres), `expected_drop` / `max_dropped_component_km2` logic,
  reproducibility with `seed`.
- `depths.product` = `external_dep_file` matching (`match_tol_m`) and `none`.
- Defaults: with none of the new keys, every code path is identical to 1988f22 (evidence of the author's md5
  regression: `scratch/regress_20261116/` -- base1b vs head*/new*; read, do not rerun the full Tokyo chain).

Real inputs: the VR-FVCOM recipe `$HOME/Github/VR-FVCOM/recipes/mesh/v3.yaml` (its `pipeline:` block, read with
`--section pipeline`), its outputs in `$HOME/Github/VR-FVCOM/scratch/mesh/v3/` (pipeline/ subdir, QA json), and
the Tokyo recipe `tokyo_bay_v6` in this repository.

## Out of scope -- do NOT report or investigate
- code the target only calls, and code that was not modified in the range;
- malformed or contrived input files, other platforms, injected faults, rare corner cases the real inputs do not
  produce (e.g. self-intersecting domain polygons, CRS strings that are not EPSG codes, empty recipes);
- performance, style and naming (unless it makes a realistic misuse likely);
- the Tokyo notebooks (`plotting.py` 32654 default kept deliberately for them; `fmesh-plot-mesh` now warns).

## What was done
Read `git log -p 1988f22..14a7dd8`. Intent: generalise the Tokyo-specific pipeline so a new region (first user:
Volta River mouth, Ghana, EPSG:32631, ~60-100 m channels in a 87 km river, lagoons, shelf to ~shelf edge) is
built fully automatically. VR v3 result (job 7591109): fmesh-mesh-qa 23/23 gates pass, NP 27,521 / NE 48,989,
1 component, 82 OBC nodes, min channel edge 59.9 m. Tokyo: md5-identical fort.14 at raw/finished/obc/siteops/
polished/final versus the old code plus the ledger fix (e557279), because the unpatched baseline crashes in
`finishing` on numpy int64.

## Previous rounds
Round 1 returned 10 major + 1 minor; all were verified as correct (each has a regression test `test_r1_<N>` in
`tests/test_new_region.py` that fails on 14a7dd8 and passes on 32d211b) and fixed in 32d211b. The triage table with
the fix for each is in `docs/newregion-review-20261116.md`, section "Round 1". In short:
1 `protected_node_mask` counts nodes inside or within 1 m of the protected outline; 2 no cap deletion inside
protected water, lone-node removal never deletes elements holding protected/OBC nodes, leftovers fail or are
reported; 3 a polygon build with an OBC line requires an `obc` or `obcfinal` stage, `raw_open_cover` recorded;
4 new-node depths interpolated on the original triangulation; 5 `expected_drop` uses intersection-area fraction;
6 `size_interpolator` clamps to the outer cell centres; 7 site fallback refused when any non-C1/C2/C4 gate newly
fails or worsens, patches pass the perpendicularity guard; 8 `latest_mesh` artifact used by depths (incl. `none`);
9 bbox transformed into each source CRS; 10 explicit `crs` forwarded to prep; 11 every build records
`raw_mesh_crs`. Known and deliberate: the polygon engine's boundary repair can delete protected water at
generation (single-element arms FVCOM cannot accept); the loss is logged and can be made fatal with
`build.max_protected_ear_area_km2` (VR v3 sets no threshold; VR's own coverage guard bounds it per water body).
Evidence: Tokyo job 7591382 md5-identical to base1b at all stages; VR v3 job 7591383 QA 23/23.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED / WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh audit of the Target above (not of code outside it).

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a realistic path
- minor: a clear robustness/clarity defect on a realistic path (an unusual or injected input does NOT qualify)
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or evidence, and a concrete fix. Then
`## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
````

### Answer (verbatim)

````markdown
Read-only review completed; no files changed.

The stored Tokyo `head1`, `head2`, `head3`, `new1`, and `new2` meshes match `base1b` at all six mesh stages. The stored Volta result identifies `32d211b`, has **27,521 nodes / 48,989 elements / 82 OBC nodes**, and passes **23/23 QA gates**. These results do not cover the defects below.

Round-1 status:

| Finding | Status | Reason |
|---|---|---|
| 1. Protection excludes shoreline nodes | PARTIAL | Outline nodes are now included, but vertex-only classification still misses protected water overlapping elements; finding 1 below. |
| 2. Lone-node removal erases protected water/OBC | PARTIAL | The finishing pass now preserves protected/OBC carriers. The mitigation for deliberately permitted generator repair remains incomplete; findings 4–5. |
| 3. Polygon build permits unrestored OBC trimming | RESOLVED | Validation now requires an `obc` or `obcfinal` stage, and raw coverage is recorded. |
| 4. Band restoration corrupts new-node depths | RESOLVED | New-node depths use interpolation over the original triangulation. |
| 5. Expected-drop classification counts elements | RESOLVED | Classification now uses triangle–polygon intersection area. |
| 6. Raster edge strip receives global maximum size | RESOLVED | Interpolation coordinates are clamped to cell centres. |
| 7. Fallback accepts newly failing gates | RESOLVED | Explicit rejection replaces the finite penalty; fallback patches receive the perpendicularity guard. |
| 8. Depths discard later mesh stages | RESOLVED | The runner records `latest_mesh`; both external depths and `none` use it. |
| 9. Projected water sources receive geographic filters | RESOLVED | Merge filtering transforms the bbox into each source CRS. |
| 10. Prep ignores resolved recipe CRS | PARTIAL | Explicit `crs` is forwarded; the supported legacy CRS override is still omitted; finding 8. |
| 11. Stale raw CRS survives engine changes | RESOLVED | Every build now records its source CRS. This introduces the figures regression in finding 7. |

Verification used source/diff inspection, artifact hashes, and in-memory synthetic meshes. The command `python -B -m pytest -q -s -p no:cacheprovider tests/test_new_region.py -k 'recipe_epsg or path_expands or cap_deletion_guard or r1_2 or r1_6 or site_fallback_not_needed or r1_7'` passed **8 tests**. Initial pytest capture and OceanMesh imports encountered unwritable temporary/cache directories. Polygon reproductions therefore used deterministic generator/clean mocks returning valid conforming meshes; full generation was not rerun.

1. **Major — Vertex-only protection still permits deletion of protected water.**

   **Location:** `src/fvcom_mesh_tools/finishing.py:367`, `:414`.

   Including outline nodes fixes bank-aligned protection, but an element overlapping a protected polygon remains unprotected when none of its vertices lies inside or within 1 m. This occurs with protected corridors narrower than the local mesh spacing.

   **Reproduction:** A conforming 4×4 node lattice with 300 m horizontal and 100 m vertical spacing, and a protected 900×60 m strip between y=20 and y=80, produces **zero protected nodes**. Calling `finish_constrained_mesh(..., on_dropped_protected="fail")`, with optimization mocked as a no-op, returns successfully after deleting two caps: **18 → 16 elements**, **9,000 m² of protected water lost**, and no component-drop record.

   **Fix:** Classify protected elements by positive triangle–polygon overlap. Carry that protection through deletion/component checks rather than relying exclusively on vertex membership.

2. **Major — Siteops protection detects disconnected components but permits direct erosion.**

   **Location:** `src/fvcom_mesh_tools/site_session.py:104`, `:225`; `src/fvcom_mesh_tools/cli/pipeline.py:709`.

   The new `protect_polygons` argument is used only after editing. The collapse operations receive OBC nodes as protected nodes, but receive no protection derived from the supplied polygons. A collapse can remove protected shoreline water while leaving one connected component, so the final component check succeeds.

   **Reproduction:** On a 4×4 lattice at 100 m spacing, move the lower-left corner from `(0,0)` to `(80,30)`, assign the top edge as OBC, and protect the entire original water polygon. One actual `apply_site_operators` pass performs `equalize` and `collapse`, returns successfully with **18 → 16 elements**, and loses **5,000 m²** of protected coverage.

   **Fix:** Apply protection before every destructive operator, including the short-edge sweep. Reject or revert edits that remove protected water; an after-edit component check is insufficient.

3. **Major — The polygon component filter does not enforce protected-water retention.**

   **Location:** `src/fvcom_mesh_tools/mesh_engine/polygon.py:359`, `:389`, `:408`, `:422`.

   `in_keep` is reporting metadata only. Both rejection lists ignore it, so an unexpected protected component below `max_dropped_component_km2` is dropped successfully. Neither `expected_drop` nor an explicit report policy is required.

   **Reproduction:** Supply a valid generated mesh containing a main component and a protected **100×100 m** component. With the default 0.5 km² threshold and no `expected_drop`, the function logs one dropped component holding protected water, removes it, and returns successfully. Its component record contains:

   ```text
   area_km2=0.01, in_keep=True, expected_drop=False, kept=False
   ```

   The protected-ear threshold can also be zero: that separate counter does not cover component deletion.

   **Fix:** Reject unexpected drops intersecting `keep_polygons`, independently of component size. Define explicitly how protection and `expected_drop` interact.

4. **Major — The new protected-repair threshold cannot be configured through the pipeline.**

   **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:447`; generator parameter at `src/fvcom_mesh_tools/mesh_engine/polygon.py:237`.

   Round 1 added `max_protected_ear_area_km2` to the generator, but `_stage_build_polygon` never forwards `build.max_protected_ear_area_km2`. A recipe setting it to zero silently receives the generator default, `None`.

   **Evidence:** Capturing the generator call from a recipe containing `build.max_protected_ear_area_km2: 0` confirms the keyword is absent. The stored Volta build reports **2,368,707.23 m²** of protected repair loss, demonstrating that this is a materially relevant guard.

   **Fix:** Forward the recipe setting to `generate_polygon_mesh` and add a pipeline-level regression that verifies excessive protected repair raises.

5. **Major — Protected-repair accounting can report zero despite substantial loss.**

   **Location:** `src/fvcom_mesh_tools/mesh_engine/polygon.py:514`.

   The new repair-loss counter counts entire deleted triangles only when their centroids lie inside protection. It misses overlapping triangles whose centroids fall outside, allowing the configured threshold to be bypassed even when calling the generator directly.

   **Reproduction:** On a conforming 5×5 node lattice at 300 m spacing, protect a **300×60 m** strip along the lower-left boundary. Boundary repair deletes its overlapping corner ear. With `max_protected_ear_area_km2=0`, generation succeeds and reports **0 m²** protected repair loss; triangle–polygon intersection measurement shows **16,200 m²** lost.

   **Fix:** Sum the intersection areas of deleted triangles with `keep_polygons`, using the same area-based principle applied by the round-1 `expected_drop` fix.

6. **Major — OBC-band restoration can introduce C4 failures while passing its guard.**

   **Location:** `src/fvcom_mesh_tools/finishing.py:332`, `:350`.

   The C4 guard compares only the maximum area-change value across affected pairs. An existing bad pair can therefore conceal newly failing pairs whose values remain below that maximum.

   **Reproduction:** Use a 4×4 lattice at 100 m spacing with the bottom edge as OBC. In the saved original mesh, offset nodes 5 and 6 by `(-40,-40)`; in the optimizer result, return node 5 to its nominal position. `_restore_obc_band` accepts restoring node 5:

   ```text
   restored=1, refused=0
   maximum C4: 0.8 → 0.8
   C4 violations: 5 → 6
   ```

   Four previously passing pairs become failing, with values approximately **0.8, 0.667, 0.667, and 0.571**.

   **Fix:** Compare C4 values pair by pair against `maximum(previous_value, threshold)`. Reject restoration when any previously passing pair starts failing or an existing violation worsens.

7. **Minor — Recording geographic raw CRS breaks raw-mesh figures, including Tokyo.**

   **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:1173`.

   The condition treats any truthy `raw_mesh_crs` as evidence of metric coordinates. Round 1 now records `"EPSG:4326"` for geographic build engines, so their raw meshes are incorrectly transformed from the recipe’s UTM CRS.

   **Reproduction:** Passing the actual Tokyo `head3` raw mesh and its recorded CRS through `_stage_figures`, with rendering mocked, transforms its first node from approximately:

   ```text
   (139.75119111, 35.45236670)
   ```

   to:

   ```text
   (136.512508, 0.000319759)
   ```

   The raw PNG is consequently generated for the wrong location. Mesh md5 equality does not catch this regression.

   **Fix:** Transform from the recorded source CRS, and skip reprojection when that CRS is already EPSG:4326.

8. **Minor — Prep still ignores a legacy CRS override.**

   **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:220`, `:275`.

   Both forwarding conditions check only the top-level `crs`. When the supported `finish.utm_epsg` alias resolves to a zone different from prep’s automatic choice, prep opening/skeletonization and endpoint extension use a different CRS from subsequent stages.

   **Reproduction:** With bbox `[-0.15, 5.0, 1.45, 6.5]` and `finish.utm_epsg: 32630`, `_recipe_epsg` returns **32630**, while the bbox’s automatic zone is **32631**. Captured prep arguments contain no `--utm-epsg`.

   **Fix:** Forward the resolved CRS for either explicit override source. Preserve the measured Tokyo behavior when its existing automatic and legacy choices agree.

## Verdict
VERDICT: FAIL (0 blocker, 6 major, 2 minor, 0 nit)
````

### Triage

Fixes:
- `670172c`: all 8 findings.
- `5c78274`: relative protected-drop tolerance. It was needed by VR v3; see #3.
- `448dbe5`: its test.

Each finding has a regression test `test_r2_<N>_*` in `tests/test_new_region.py`. All 8 fail on the reviewed code 32d211b (`PYTHONPATH` = `git archive 32d211b`) and pass on 448dbe5. Behaviour of the reviewed code for #1 and #2: `scratch/reviews/newregion/r1_evidence/r2_old_behaviour.txt`.

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes. Reviewed code, 300 x 100 m lattice, protected strip y 20-80 m: 0 protected nodes, finish 18 -> 17 elements, protected coverage 54,000 -> 45,000 m2, no exception. | yes | 670172c: `protected_element_mask` (protected vertex or positive overlap area; geometry tiled, 49k elements in 0.1 s); used for cap deletion, component drops (`keep_components_protected(protect_polygons=...)`) and lone-node removal. Fixed code: 18 -> 18 elements, 54,000 m2 kept. |
| 2 | major | yes. Reviewed code, the reviewer's lattice (corner moved to (80, 30)): `equalize` and `collapse`, 18 -> 16 elements, protected coverage 85,500 -> 85,000 m2. | yes | 670172c: before each collapse (short-edge sweep and site operators), a collapse with both ends on the mesh boundary next to protected water is refused (equalize keeps boundary apexes and split works on interior edges, so neither can change coverage). As a net, a protected-coverage loss over 1 m2 fails the stage. Fixed code: 85,500 m2 kept. |
| 3 | major | yes. `test_r2_3`: a 0.36 km2 protected pond was dropped silently by the reviewed code. | yes | 670172c + 5c78274. Rule: a dropped component of any size may not take protected water outside `expected_drop`. Allowed: max(`protected_drop_tol_m2` = 1 m2, `protected_drop_tol_frac` = 0.01 x its protected area). Else it fails (`protected_drop_policy: report` only logs). New `build.expected_drop_unresolved_m: w` adds the water the open-boundary component of the water opened at width w does not reach. VR v3 now sets 150 m (the creek resolve width): 16.1 km2. The relative 1 % came from VR job 7591568: a 4-element fragment in a neck had 20 m2 of its protected water just outside that geometry. |
| 4 | major | yes. `test_r2_4`: `build.max_protected_ear_area_km2: 0` did not reach the engine. | yes | 670172c: forwarded, together with `protected_drop_policy` / `_tol_m2` / `_tol_frac`; pipeline-level test raises. No threshold is set in VR v3 (owner to decide). |
| 5 | major | yes. `test_r2_5`: an ear overlapping a protected strip but with its centroid outside counted 0 m2 in the reviewed code. | yes | 670172c: intersection areas. Total loss and the part outside `expected_drop` are both reported; the threshold applies to the latter. VR v3: 2.21 km2 total, 0.63 km2 outside the expected drops. |
| 6 | major | yes. `test_r2_6`, reviewer's case: C4 violations 5 -> 6 after the restore in the reviewed code. | yes | 670172c: per-pair comparison against max(previous, 0.5) |
| 7 | minor | yes. `test_r2_7`: a lon/lat raw mesh was reprojected from UTM (139.75, 35.45) -> (136.51, 0.0003) in the reviewed code. | yes | 670172c: transform from the recorded CRS, none for EPSG:4326. The Tokyo raw PNG (job 7591693) is placed at 139.6-140.1 E, 35.0-35.7 N. |
| 8 | minor | yes. `test_r2_8`: no `--utm-epsg` with `finish.utm_epsg: 32630` over a zone-31 bbox. | yes | 670172c: `_prep_crs_override`: explicit `crs`, or the legacy key when it differs from prep's own zone. For Tokyo they agree (32654), so the argv is unchanged. |

No rebuttals.

Checks after the fixes (HEAD 448dbe5):
- `pytest -q`: 1342 passed, 9 skipped. `ruff`: clean.
- Tokyo `tokyo_bay_v6` (job 7591693, figures included): raw, finished, obc, siteops, polished and final fort.14 are md5-identical to base1b. QA still has 5 gate failures (unchanged). The raw figure is now correctly placed.
- VR v3 (VR-FVCOM 36ff15f, job 7591579):
  - fmesh-mesh-qa 23/23; NP 27,521 / NE 48,989; one component; 82 open-boundary nodes.
  - Coverage guard passed; depth sha256 check passed.
  - Min edge channel / lagoon 62.0 / 66.0 m; dt_E 9.72 s.
  - Siteops protected coverage unchanged (375,651,379 m2).
