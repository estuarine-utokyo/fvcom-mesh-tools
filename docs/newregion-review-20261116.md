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

## Round 3

### Prompt

````markdown
# Review request, round 3: new-region generalisation of fmesh-pipeline (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT modify files. You may run read-only
commands, python in memory, mocks and small synthetic meshes (keep each computation under a few minutes). Answer
in English as Markdown.

## Goal
Correctness and robustness of the target below under realistic use. Report every defect you can substantiate in
the target, of any severity.

## Target (concrete)
Commit range `1988f22..448dbe5` (32d211b holds the round-1 fixes; 670172c, 5c78274, 448dbe5 the round-2 fixes). Files:
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
Round 1 (10 major, 1 minor): fixed in 32d211b; round 2 judged 8 RESOLVED and 3 PARTIAL (#1, #2, #10), which are
covered by round-2 findings below.
Round 2 (6 major, 2 minor): all verified correct (tests `test_r2_1`..`test_r2_8` fail on 32d211b, pass on 448dbe5)
and fixed in 670172c (+ 5c78274, 448dbe5). Triage in `docs/newregion-review-20261116.md`, "Round 2":
1 `protected_element_mask`: protected if a protected vertex or positive overlap area (tiled geometry), used for cap
deletion, component drops and lone-node removal; 2 collapses with both ends on the boundary next to protected water
are refused before each collapse (short-edge sweep and site operators), and any protected-coverage loss > 1 m2 fails
the stage; 3 a dropped component of any size may not take protected water outside `expected_drop` beyond
max(1 m2, `protected_drop_tol_frac` (default 1 %) of its protected area) unless `protected_drop_policy: report`;
new `build.expected_drop_unresolved_m: w` adds water not reached by the open-boundary component of the water opened
at width w (VR v3 uses 150 m, VR-FVCOM commit 36ff15f); 4 `max_protected_ear_area_km2` and the drop policy and
tolerances are forwarded by the pipeline; 5 repair loss counted by intersection area (total and outside
expected_drop; the threshold applies to the latter); 6 C4 compared pair by pair in band restore; 7 raw mesh
transformed from its recorded CRS, none for EPSG:4326; 8 legacy `finish.utm_epsg` forwarded to prep when it differs
from prep's own zone.
Evidence: Tokyo job 7591693 (with figures) md5-identical to base1b at all stages, raw PNG now at 139.6-140.1 E;
VR v3 job 7591579 QA 23/23, coverage guard pass, generation repair loss 2.21 km2 (0.63 km2 outside expected drops;
the owner has not yet chosen a threshold, so none is set -- not a defect).

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
Read-only review completed; no files changed. The source and tests in the working tree match `448dbe5`.

The Tokyo `base1b`, `head1`–`head5`, `new1`, and `new2` artifacts are md5-identical at all six mesh stages. The saved Volta final mesh passes 23/23 gates, but the checks below demonstrate remaining silent failures.

The previous findings have these statuses:

| Finding | Status | Reason |
|---|---|---|
| R1 #1: Protection excludes shoreline nodes | RESOLVED | Boundary-inclusive node classification and element overlap now protect cap deletion. |
| R1 #2: Lone-node removal erases protected water/OBC | RESOLVED | Protected and OBC carriers are retained; unresolved protected lone nodes fail by default. Generator repair loss is measured and configurable. |
| R1 #3: Polygon build permits unrestored OBC trimming | RESOLVED | Recipe validation requires an OBC restoration stage and records raw coverage. Finding 4 below concerns subsequent fallback execution. |
| R1 #4: Band restoration corrupts new-node depths | RESOLVED | New-node depths use interpolation over the original triangulation. |
| R1 #5: Expected-drop classification counts elements | RESOLVED | The fraction uses intersection area. |
| R1 #6: Raster edges receive global maximum size | RESOLVED | Interpolation clamps coordinates to the cell-centre lattice. |
| R1 #7: Fallback accepts newly failing gates | PARTIAL | Explicit gate rejection replaces the finite penalty, but OBC preservation and the configured perpendicularity guard remain incomplete; findings 4–5. |
| R1 #8: Depths discard later mesh stages | RESOLVED | `latest_mesh` is maintained and used by both external depths and `none`. |
| R1 #9: Projected water sources receive geographic filters | RESOLVED | Filtering transforms the bbox into each source CRS. |
| R1 #10: Prep ignores resolved CRS | RESOLVED | Explicit and differing legacy overrides reach prep and endpoint extension. |
| R1 #11: Stale raw CRS survives engine changes | RESOLVED | Every build records its source CRS, which finish uses for transformation. |
| R2 #1: Vertex-only protection misses overlapping water | RESOLVED | Cap deletion, component checks and lone-node removal use element overlap. |
| R2 #2: Siteops permits direct protected-water erosion | RESOLVED | Protected boundary collapses are refused and coverage is checked. |
| R2 #3: Component filtering drops protected water | PARTIAL | The new measurement skips components containing at least half the elements; finding 1. |
| R2 #4: Repair threshold is not forwarded | RESOLVED | The pipeline forwards the threshold, policy and tolerances. |
| R2 #5: Repair accounting misses overlap | RESOLVED | Deleted-triangle intersection areas are counted, including the portion outside expected drops. |
| R2 #6: Band restoration hides new C4 failures | RESOLVED | C4 comparisons are now pair by pair. |
| R2 #7: Geographic raw figures are reprojected incorrectly | RESOLVED | Figures use the recorded source CRS and skip geographic reprojection. |
| R2 #8: Prep ignores legacy CRS override | RESOLVED | A legacy override differing from prep’s automatic zone is forwarded. |

1. **Major — The round-2 component fix skips a finely meshed protected component.**

   **Location:** [mesh_engine/polygon.py:377](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/polygon.py:377).

   Protection and expected-drop measurements run only when `sel.sum() < len(t) // 2`. With an explicit OBC, the component containing most elements need not be the component retained.

   **Reproduction:** An actual generator run used a 3,000 × 1,333.33 m sea rectangle and a separate protected 600 × 600 m pond. Sizes were approximately 333 m offshore and 60 m around the pond, with `seed=0`, `max_iter=30`, `ladder_k=1`, `clean_dj_cutoff=0.0005`, default protected-drop policy, and `max_protected_ear_area_km2=0`.

   The post-clean pond contained **127 of 177 elements** and **318,472.65 m²**. Generation dropped it and returned successfully. Its metadata incorrectly reported `in_keep=False`, `protected_m2=0`, and `protected_unexpected_m2=0`. A deterministic conforming-lattice reproduction likewise deleted the entire 0.36 km² pond.

   **Fix:** Determine the actual retained component IDs first. Measure every component that could be dropped; do not infer retention from element count. This defect was introduced by the round-2 protection fix.

2. **Major — Finishing protection does not cover welding, which can delete protected water and alter OBC nodes.**

   **Location:** [finishing.py:143](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/finishing.py:143), [finishing.py:186](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/finishing.py:186).

   The new guarded path protects cap deletion, but welding subsequently runs without that protection. Band tagging happens after welding, and there is no final protected-coverage or OBC-preservation check.

   **Reproduction:** Use boundary vertices `(0,0), (1,-1), (100,0), (100,100), (0,100)`, a centre node `(50,50)`, and five centre-fan triangles. Protect the complete polygon and set OBC nodes to `(1,-1), (100,0)`.

   With the actual welder and optimizer, default 2 m welding tolerance, `on_dropped_protected="fail"` and `freeze_obc_band=True`, finishing succeeds:

   ```text
   nodes/elements:       6/5 → 4/2
   protected coverage:  10,050 → 10,000 m²
   first OBC position:  (1,-1) → (0,0)
   caps removed:        0
   welded nodes:        1
   ```

   Only filesystem access was mocked. The saved Volta polish also reports one welded node and one dropped element.

   **Fix:** Extend the guarded path across welding and optimization. Refuse edits removing protected water or changing protected OBC geometry, and validate those invariants before returning.

3. **Major — Automated finishing bypasses protected-water coverage checks.**

   **Location:** [cli/pipeline.py:908](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:908), [cli/pipeline.py:949](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:949).

   `_execute_patches_guarded` checks perpendicularity only. It receives no protected geometry. The fallback’s later patch round has the same gap, despite its site-operator portion checking coverage.

   **Evidence:** I reran Volta’s `obcfinal` and automated finishing entirely in memory, mocking file writes. The resulting final coordinates and connectivity match the saved final mesh exactly.

   `obcfinal` preserves protected coverage. Subsequent automated finishing changes it from **375,653,269.04 to 375,649,415.40 m²**. Geometric differencing measures **4,827.59 m² lost** and **973.95 m² gained elsewhere**. The final mesh still passes all 23 gates.

   **Fix:** Pass protection through both patch-execution branches and fallback replanning. Reject trials that lose previously covered protected water, using geometric loss rather than allowing gains elsewhere to offset erosion.

4. **Major — Site fallback can accept complete loss of the OBC when protection is omitted.**

   **Location:** [cli/pipeline.py:1009](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:1009), [cli/pipeline.py:1024](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:1024).

   Without `protect_polygons`, fallback calls the site operators’ legacy component path, which clears boundary lists. The fallback neither restores nor validates the original OBC. OBC checks disappear from its failed-gate dictionary, allowing acceptance.

   **Reproduction:** Use a 5 × 5 lattice at 100 m spacing, node ID `5j+i`. Triangulate each cell with `(a,a+1,a+5)` and `(a+1,a+6,a+5)`, replacing the bottom-left pair with `(0,1,6),(0,6,5)` and the top-right pair with `(18,19,24),(18,24,23)`. Set OBC `[0,1,2,3,4]`, and move:

   ```text
   node 8 → (343.126008, 41.208301)
   node 7 → (214.748498, 22.858384)
   ```

   Calling the actual fallback with the rectangular shoreline and no protection returns:

   ```text
   before: OBC-perpendicularity=2, C1=7, C4=4
   after:  no failing gates
   outcome: kept
   open_boundaries: []
   ```

   **Fix:** Make fallback preserve boundaries independently of channel protection, and explicitly validate OBC existence, identity, order and boundary-edge adjacency before accepting.

5. **Major — The configured perpendicularity limit does not cover the complete fallback.**

   **Location:** [cli/pipeline.py:1009](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:1009), [cli/pipeline.py:1018](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:1018).

   Site operators run before `_execute_patches_guarded` establishes its reference deviations. The final acceptance check uses QA counts, whose perpendicularity threshold remains 20°, rather than the configured limit.

   **Reproduction:** Use the lattice from finding 4, protect the complete rectangle, set `obc_perp_max_deg=10`, and instead move:

   ```text
   node 17 → (232.394916693, 262.252995411)
   node 11 → (79.252582251, 161.262146908)
   ```

   Initially every OBC deviation is **0°**. The actual fallback accepts the result, removing all failing gates, while producing deviations:

   ```text
   [0°, 11.4891°, 13.8612°, 0°, 0°]
   ```

   These exceed the requested 10° guard while passing the QA threshold.

   **Fix:** Compare the complete fallback result against the original mesh’s per-node deviations and the configured limit. Reject or revert the whole trial when that comparison fails.

Verification included **17 passing filesystem-free tests**, actual small generator runs, in-memory finishing runs and artifact comparisons. Repeated polygon generation with the same seed produced identical points, cells and OBC strings. A protected fan test reduced maximum valence **10 → 5**, preserved coverage/OBC coordinates, and passed the conformity checks.

The full suite and Tokyo generation were not rerun. Initial OceanMesh import failed because Matplotlib required writable cache access; subsequent checks reused its existing cache with that access probe mocked. GDAL also emitted a PROJ database diagnostic; the reported runs completed.

## Verdict

VERDICT: FAIL (0 blocker, 5 major, 0 minor, 0 nit)
````

### Triage

Fixes, all marked `review newregion r3 #N`:
- `ef1e8b3`: the 5 findings.
- `60201fa`: follow-ups found while rebuilding VR under the new rules (#2).
- `9a3ec4d`: open-line coverage check in the obc stages (#2).

Regression tests `test_r3_1`..`test_r3_5` fail on the reviewed code 448dbe5 (`PYTHONPATH` = `git archive 448dbe5`) and pass on 9a3ec4d. The same holds for the follow-up tests: coincident OBC weld, lone node by flip and by split, open-line coverage.

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes. `test_r3_1`, the reviewer's case: a separate 0.36 km2 protected pond, meshed at 60 m so it holds most of the elements, was dropped and generation succeeded. | yes | ef1e8b3: the retained components are decided first (`anchor_components` inside `comp_rows`); every other component is measured |
| 2 | major | yes. `test_r3_2`, the reviewer's fan: the reviewed code welded the OBC node (1,-1) into (0,0) and lost 50 m2. | yes | ef1e8b3 + 60201fa + 9a3ec4d. Weld: never moves an OBC node except onto a coincident one (< 1 mm, a repair). A weld that loses protected water is refused. Optimizer: protected bank nodes go back where it lost water; a loss still above `protected_loss_tol_m2` (default 0), or a broken string, refuses the round. Finish validates the string and the protected loss before returning. Single-element nodes are fixed by an edge flip or split before any deletion. finish/polish keep the OBC string; its geometry is re-derived by the obc/obcfinal stages, which the recipe must have (r1 #3). Those stages now fail when the string no longer covers the supplied line as far as the raw mesh did, or a node is more than `max_move_m` off it. |
| 3 | major | yes. `test_r3_3`: on the reviewed code `_execute_patches_guarded` has no protection argument. The reviewer measured 4,827.6 m2 lost / 974 m2 gained on VR. | yes | ef1e8b3: `protected_loss_m2` gives the geometric loss over the changed triangles (no netting). Every patch (both branches) and the fallback patch round is reverted above `finishing.auto.protected_loss_tol_m2` (default 0). The site operators use the same geometric check. The stage logs water lost and gained. |
| 4 | major | yes. `test_r3_4`, the reviewer's lattice: the reviewed code returned `open_boundaries == []`. | yes | ef1e8b3: `apply_site_operators(keep_open_boundary=True)` in the fallback; the trial's open boundary is validated (existence, node count, positions on the original line in order, boundary edges) |
| 5 | major | yes. `test_r3_5`: the reviewed code accepted deviations 11.5 and 13.9 deg with `obc_perp_max_deg` 10. | yes | ef1e8b3: the whole trial is compared per node with max(original deviation, `obc_perp_max_deg`), with the reference taken before the site operators |

No rebuttals.

#### Trade-off of the strict no-loss rule on VR v3 (protected water 375 km2)

Each variant is a full VR chain on GENKAI (`scratch/mesh/r3_variants/`).

| setting (finish/polish tol, autofinish tol per patch) | QA gates | protected water lost | job |
|---|---|---|---|
| 0 / 0 m2 (default) | 20/23 (C1 18, C2 2, C4 7): the optimizer is refused in finish and polish (it would lose 107 / 145 m2 even after restoring 2,359 bank nodes) | 0 | 7592227 |
| 500 / 0 m2, bank nodes restored | 19/23: the restore undoes the optimizer's bank work | finish 107, polish 37 m2 | 7592101 |
| 400,000 / 0 m2 | 21/23 (C1 1 at 29.5 deg, C4 3 at 0.606) | finish 293,401, polish 3,397, autofinish 0 m2 | 7592277 |
| **400,000 / 1,500 m2 (VR v3)** | **23/23**; open line covered 73-183,697 of 183,697 m | finish 293,401, polish 3,397, autofinish 1,980 lost / 974 gained m2 (one patch reverted at 1,807 m2) | 7592343 |

The optimizer slides boundary nodes along the shoreline. Measured as geometric loss without netting, that moves 0.29 km2 (0.08 %) of the protected water out of the mesh in finish. The default stays 0 m2. VR v3 sets the tolerances explicitly in its recipe, with this table as the reason. The owner should review that choice.

A variant run while doing this also exposed a defect, now closed in 9a3ec4d. With polish run without `freeze_obc_band` (job 7592310), QA passed 23/23 while the eastern 26.7 km of the open line had become land. The obc stages now refuse that.

Checks after the fixes (HEAD 9a3ec4d):
- `pytest -q`: 1351 passed, 9 skipped. `ruff`: clean.
- Tokyo `tokyo_bay_v6` (job 7592397, figures included): md5-identical to base1b at all six stages; raw figure correctly placed.
- VR v3 (VR-FVCOM 63364eb, job 7592343):
  - QA 23/23; NP 27,521 / NE 48,989; one component; 82 open-boundary nodes covering the line.
  - Coverage guard passed; depth sha256 check passed.
  - Min edge channel / lagoon 62.0 / 66.0 m; dt_E 9.72 s.

## Round 4

### Prompt

````markdown
# Review request, round 4: new-region generalisation of fmesh-pipeline (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT modify files. You may run read-only
commands, python in memory, mocks and small synthetic meshes (keep each computation under a few minutes). Answer
in English as Markdown.

## Goal
Correctness and robustness of the target below under realistic use. Report every defect you can substantiate in
the target, of any severity.

## Target (concrete)
Commit range `1988f22..9a3ec4d` (32d211b round-1 fixes; 670172c, 5c78274, 448dbe5 round-2 fixes; ef1e8b3, 60201fa, 9a3ec4d round-3 fixes). Files:
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
Round 1 (10 major, 1 minor): fixed in 32d211b; round 2 judged 8 RESOLVED and 3 PARTIAL (#1, #2, #10), which are
covered by round-2 findings below.
Round 2 (6 major, 2 minor): all verified correct (tests `test_r2_1`..`test_r2_8` fail on 32d211b, pass on 448dbe5)
and fixed in 670172c (+ 5c78274, 448dbe5). Triage in `docs/newregion-review-20261116.md`, "Round 2":
1 `protected_element_mask`: protected if a protected vertex or positive overlap area (tiled geometry), used for cap
deletion, component drops and lone-node removal; 2 collapses with both ends on the boundary next to protected water
are refused before each collapse (short-edge sweep and site operators), and any protected-coverage loss > 1 m2 fails
the stage; 3 a dropped component of any size may not take protected water outside `expected_drop` beyond
max(1 m2, `protected_drop_tol_frac` (default 1 %) of its protected area) unless `protected_drop_policy: report`;
new `build.expected_drop_unresolved_m: w` adds water not reached by the open-boundary component of the water opened
at width w (VR v3 uses 150 m, VR-FVCOM commit 36ff15f); 4 `max_protected_ear_area_km2` and the drop policy and
tolerances are forwarded by the pipeline; 5 repair loss counted by intersection area (total and outside
expected_drop; the threshold applies to the latter); 6 C4 compared pair by pair in band restore; 7 raw mesh
transformed from its recorded CRS, none for EPSG:4326; 8 legacy `finish.utm_epsg` forwarded to prep when it differs
from prep's own zone.
Round 3 (5 major): all verified correct (tests `test_r3_1`..`test_r3_5` fail on 448dbe5, pass on 9a3ec4d), fixed in
ef1e8b3 (+ 60201fa, 9a3ec4d). Triage and a trade-off table in `docs/newregion-review-20261116.md`, "Round 3":
1 retained components decided first, every dropped component measured; 2 weld guard (OBC nodes only onto coincident
< 1 mm nodes, refused when protected water is lost), optimizer guard (protected bank nodes restored, round refused
above `protected_loss_tol_m2`), OBC string and protected loss validated before finish returns; lone nodes fixed by
flip/split before any deletion; obc/obcfinal check that the open string still covers the supplied line;
3 `protected_loss_m2` = geometric loss without netting; every autofinish patch (both branches and fallback round)
reverted above `finishing.auto.protected_loss_tol_m2`; 4 fallback keeps and validates the OBC; 5 whole fallback
trial compared per node against max(original deviation, `obc_perp_max_deg`).
Library defaults of all protected-loss tolerances are 0 m2. The VR v3 recipe (VR-FVCOM 7fb58cf) sets explicit
tolerances (finish/polish 400,000 m2 per step, autofinish 1,500 m2 per patch) because the optimizer slides bank nodes
along the shoreline; with 0 everywhere VR reaches 20/23. The tolerance VALUES are an owner decision pending, not a
defect; whether the code applies them correctly is in scope.
Evidence: Tokyo job 7592397 (with figures) md5-identical to base1b at all stages; VR v3 job 7592343 QA 23/23,
82 OBC nodes covering 73-183,697 m of the 183,697 m line, coverage guard pass.

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
Four major defects remain. No files were modified; the target source files match `9a3ec4d`.

Verification included 27 existing test functions, using in-memory I/O where necessary, plus focused synthetic-mesh probes. Initial cache/filesystem errors were resolved with narrower mocks; GDAL emitted a PROJ diagnostic during otherwise successful generator runs. Full regional generation was not rerun.

The stored Tokyo `head1`–`head8`, `new1`, and `new2` meshes are md5-identical to `base1b` at all six stages. The saved Volta result records NP **27,521**, NE **48,989**, **82 OBC nodes**, and **23/23** passing gates. Repeated small polygon generation with the same seed produced identical points, elements, and OBC strings.

Previous findings have these statuses:

| Finding | Status | Reason |
|---|---|---|
| R1 #1: Protection excludes shoreline nodes | RESOLVED | Boundary-inclusive node masks and triangle overlap protect cap deletion. |
| R1 #2: Lone-node removal erases protected water/OBC | RESOLVED | Flip/split repair preserves coverage; deletion guards protect water and OBC carriers. Polygon repair loss is measured and configurable. |
| R1 #3: Polygon build permits unrestored OBC trimming | PARTIAL | A restoration stage is required, but its coverage guard has findings 2–3 below. |
| R1 #4: Band restoration corrupts new-node depths | RESOLVED | New-node depths use interpolation over the original triangulation. |
| R1 #5: Expected-drop classification counts elements | RESOLVED | Classification uses intersection area. |
| R1 #6: Raster edges receive maximum size | RESOLVED | Interpolation clamps to the centre lattice. |
| R1 #7: Fallback accepts newly failing gates | RESOLVED | Explicit gate rejection and whole-trial perpendicularity checks replace the incomplete scoring contract. |
| R1 #8: Depths discard later stages | RESOLVED | Both depth products use the runner’s `latest_mesh`. |
| R1 #9: Geographic filters applied to projected water | RESOLVED | The bbox is transformed into the source CRS before filtering. |
| R1 #10: Prep ignores resolved CRS | RESOLVED | Explicit and differing legacy CRS overrides reach prep and endpoint extension. |
| R1 #11: Stale raw CRS survives engine changes | RESOLVED | Every build records its source CRS; finish transforms from it. |
| R2 #1: Vertex-only protection misses overlap | RESOLVED | Cap, component, and lone-node checks use element overlap. |
| R2 #2: Siteops permits direct protected erosion | RESOLVED | Boundary collapses are guarded and geometric loss is checked. |
| R2 #3: Component filtering drops protected water | RESOLVED | Actual retained components are determined before measuring every possible drop. |
| R2 #4: Repair threshold is not forwarded | RESOLVED | The pipeline forwards the threshold, policy, and tolerances. |
| R2 #5: Repair accounting misses overlap | RESOLVED | Intersection areas measure total and unexpected loss. |
| R2 #6: Band restoration hides C4 failures | RESOLVED | Restoration compares every affected element pair. |
| R2 #7: Geographic raw figures reprojected incorrectly | RESOLVED | Figures use the recorded source CRS. |
| R2 #8: Prep ignores legacy CRS override | RESOLVED | A differing legacy override is forwarded. |
| R3 #1: Fine protected component omitted from measurement | RESOLVED | Measurement follows the actual retained set. |
| R3 #2: Welding/optimization can lose protection or OBC | PARTIAL | Weld and optimizer guards work; deferred OBC geometry validation remains incomplete through findings 2–3. |
| R3 #3: Autofinish bypasses protected coverage checks | RESOLVED | Both execution branches and fallback patches enforce geometric loss tolerances. |
| R3 #4: Fallback loses OBC without protection | RESOLVED | Fallback preserves and validates the OBC independently of protection. The standalone stage still has finding 4. |
| R3 #5: Perpendicularity guard misses site operations | RESOLVED | The complete fallback trial is compared against the original deviations and configured limit. |

1. **Major — OBC aftercare deletes protected elements before checking protection.**

   **Location:** `src/fvcom_mesh_tools/obc_tools.py:460`, `:466`.

   The new protection branch checks components remaining after `_rm(mesh, ~bad)`. It never checks protected water inside the elements just deleted. Direct loss therefore escapes when no surviving component is dropped.

   **Reproduction:** Construct a conforming 10 × 3 lattice with nodes
   `(300000 + 100i, 640000 + 100j)`, triangulating each cell as
   `(a,a+1,a+10)` and `(a+1,a+11,a+10)`. Protect its complete 900 × 200 m water rectangle and supply its southern edge as the OBC line. Call `assign_west_south_obc` with EPSG:32631, `snap=False`, and `max_move_m=60`.

   The actual routine returns successfully:

   ```text
   NP / NE:                 30 / 36 → 27 / 32
   structural deletions:    4 elements
   protected water lost:    20,000 m²
   components_dropped:      []
   ```

   This is a gap in the newly added OBC protection path. The pipeline coverage guard can catch some accompanying OBC damage, but it does not measure protected loss and is bypassed on resumed runs as described below.

   **Fix:** Guard direct aftercare deletions using protected-element overlap, and validate geometric protected-water loss across the complete OBC operation. Refuse destructive repair or fail before returning the mesh.

2. **Major — Resumed runs silently disable the new OBC coverage guard.**

   **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:719`, `:1431`.

   `_check_open_cover` requires `artifacts["raw_open_cover"]`. This value is a dictionary, while provenance persists only string-valued artifacts. The next invocation loads only that filtered artifact mapping, so ordinary `--only obc` or `--only obcfinal` runs lack the reference and return from the guard without checking anything.

   **Evidence:** The real Volta provenance contains:

   ```text
   results.build.raw_open_cover.uncovered_m = 9847.255995723273
   artifacts.raw_open_cover                 = absent
   ```

   Calling the guard with the saved artifact mapping returns `None`; supplying the runtime coverage dictionary performs the check. Once another partial invocation rewrites `results`, even the previous build result is no longer retained there.

   This defect was introduced with the round-3 coverage guard.

   **Fix:** Persist and reload structured validation artifacts. For a polygon recipe with an OBC line, missing coverage metadata must trigger reconstruction or a loud failure, rather than the legacy bypass.

3. **Major — The OBC guard authorizes the raw mesh’s known artificial walls.**

   **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:735`.

   The allowed uncovered length is taken from the raw mesh *after* R4 trimming. That trimming is precisely the damage the required OBC stage is supposed to repair. Consequently, the guard permits both incomplete restoration and subsequent regression to that damaged extent.

   **Evidence:** For the saved Volta meshes:

   ```text
   supplied line:              183,697.17 m
   raw uncovered length:         9,847.26 m
   final uncovered length:          73.04 m
   ```

   Passing the unchanged saved raw mesh to the post-OBC guard succeeds, despite its known land-labelled sections along the supplied OBC line. A guard probe shortening the saved final open string by ten nodes also succeeds with **6,783.77 m uncovered**, because that remains below the raw allowance.

   This is a defect in the round-3 validation fix: “no worse than trimmed raw” does not establish preservation of the supplied open boundary.

   **Fix:** Validate against the intended wet extent of the supplied line, with an explicit endpoint tolerance. Persist successful restored coverage and prevent later stages from retreating from it.

4. **Major — Standalone siteops still clears the OBC when protection is absent.**

   **Location:** `src/fvcom_mesh_tools/cli/pipeline.py:802`; `src/fvcom_mesh_tools/site_session.py:297`.

   Round 3 added `keep_open_boundary=True` for the fallback, but `_stage_siteops` does not pass it. Without `protect_polygons`, standalone siteops still takes the legacy component path that clears boundary lists.

   **Reproduction:** Use the existing `_r3_lattice()` 5 × 5 fixture with OBC `[0,1,2,3,4]` and its rectangular shoreline. An actual one-pass `apply_site_operators` call with the standalone stage’s arguments produces:

   ```text
   elements:             32 → 32
   open strings:         one five-node string → []
   QA failing gates:     []
   OBC checks:           skipped
   ```

   `_check_recipe` accepts a new-region chain containing `obc`, `siteops`, `depths: {product: none}`, and `qa`, without `obcfinal`. The depth stage correctly propagates the latest mesh—now missing its open boundary—and QA accepts it.

   **Fix:** Enable OBC preservation for standalone new-region siteops independently of protection, and validate its result against the input OBC. Preserve the measured Tokyo behavior through an explicit legacy-path condition.

## Verdict

VERDICT: FAIL (0 blocker, 4 major, 0 minor, 0 nit)
````

### Triage

Fixes, marked `review newregion r4 #N`:
- `f076102`: the 4 findings.
- `c292eff`: follow-up to #4, found by the VR rebuild.

Regression tests `test_r4_1`..`test_r4_4` fail on the reviewed code 9a3ec4d (`PYTHONPATH` = `git archive 9a3ec4d`) and pass on c292eff.

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes. `test_r4_1`, the reviewer's 10 x 3 lattice: the reviewed code returned after deleting protected elements. | yes | f076102: an R4 / {1,2,1} / pinch element holding protected water is not deleted (the stage fails, naming it). The whole `assign_west_south_obc` may not lose protected water (geometric loss, `obc.protected_loss_tol_m2`, default 0). |
| 2 | major | yes. `test_r4_2`: the reviewed code returned `None` (no check) for a resumed run; `raw_open_cover` was not persisted. | yes | f076102: the guard needs only the domain water polygon (a path artifact) and fails loudly without it. `open_cover` / `raw_open_cover` dicts are persisted in the provenance and reloaded. A `--only obc` rerun checks again (test). |
| 3 | major | yes. `test_r4_3`: the reviewed code accepted a string 100 m short (allowance = the trimmed raw mesh's 9,847 m). | yes | f076102: the reference is the wet extent of the supplied line (line ∩ domain water outline). Each end must lie within `obc.end_tol_m` (default 250 m); every node within `max_move_m`. Each stage persists its coverage (`open_cover`), and a later stage may not retreat from it by more than 1 m. VR final: end gaps 73 / 0 m. |
| 4 | major | yes. `test_r4_4`: the reviewed `_stage_siteops` cleared the open boundary of a region recipe without protection. | yes | f076102: siteops of a region recipe (polygon build with a line, or protected water) keeps and validates the open boundary. Tokyo-style recipes take the legacy component step (Tokyo md5-identical). c292eff: in the boundary-edge check, coincident (< 1 mm) open nodes count as one. The OBC assignment leaves such a pair, which the polish weld merges later; the first VR rebuild (job 7592472) stopped on it. |

No rebuttals.

Checks after the fixes (HEAD c292eff):
- `pytest -q`: 1357 passed, 9 skipped. `ruff`: clean.
- Tokyo `tokyo_bay_v6` (job 7592486, figures included): md5-identical to base1b at all six stages.
- VR v3 (job 7592487, recipe tolerances unchanged):
  - QA 23/23; NP 27,521 / NE 48,989; 82 open-boundary nodes.
  - Open line covered 73-183,697 m of the wet 0-183,697 m at both obc and obcfinal; `open_cover` artifact persisted.
  - Coverage guard passed; depth check passed.
  - Min edge channel / lagoon 62.0 / 66.0 m; dt_E 9.72 s.

## Round 5

### Prompt

````markdown
# Review request, round 5: new-region generalisation of fmesh-pipeline (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT modify files. You may run read-only
commands, python in memory, mocks and small synthetic meshes (keep each computation under a few minutes). Answer
in English as Markdown.

## Goal
Correctness and robustness of the target below under realistic use. Report every defect you can substantiate in
the target, of any severity.

## Target (concrete)
Commit range `1988f22..c292eff` (32d211b round-1 fixes; 670172c, 5c78274, 448dbe5 round-2 fixes; ef1e8b3, 60201fa, 9a3ec4d round-3 fixes; f076102, c292eff round-4 fixes). Files:
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
Round 1 (10 major, 1 minor): fixed in 32d211b; round 2 judged 8 RESOLVED and 3 PARTIAL (#1, #2, #10), which are
covered by round-2 findings below.
Round 2 (6 major, 2 minor): all verified correct (tests `test_r2_1`..`test_r2_8` fail on 32d211b, pass on 448dbe5)
and fixed in 670172c (+ 5c78274, 448dbe5). Triage in `docs/newregion-review-20261116.md`, "Round 2":
1 `protected_element_mask`: protected if a protected vertex or positive overlap area (tiled geometry), used for cap
deletion, component drops and lone-node removal; 2 collapses with both ends on the boundary next to protected water
are refused before each collapse (short-edge sweep and site operators), and any protected-coverage loss > 1 m2 fails
the stage; 3 a dropped component of any size may not take protected water outside `expected_drop` beyond
max(1 m2, `protected_drop_tol_frac` (default 1 %) of its protected area) unless `protected_drop_policy: report`;
new `build.expected_drop_unresolved_m: w` adds water not reached by the open-boundary component of the water opened
at width w (VR v3 uses 150 m, VR-FVCOM commit 36ff15f); 4 `max_protected_ear_area_km2` and the drop policy and
tolerances are forwarded by the pipeline; 5 repair loss counted by intersection area (total and outside
expected_drop; the threshold applies to the latter); 6 C4 compared pair by pair in band restore; 7 raw mesh
transformed from its recorded CRS, none for EPSG:4326; 8 legacy `finish.utm_epsg` forwarded to prep when it differs
from prep's own zone.
Round 3 (5 major): all verified correct (tests `test_r3_1`..`test_r3_5` fail on 448dbe5, pass on 9a3ec4d), fixed in
ef1e8b3 (+ 60201fa, 9a3ec4d). Triage and a trade-off table in `docs/newregion-review-20261116.md`, "Round 3":
1 retained components decided first, every dropped component measured; 2 weld guard (OBC nodes only onto coincident
< 1 mm nodes, refused when protected water is lost), optimizer guard (protected bank nodes restored, round refused
above `protected_loss_tol_m2`), OBC string and protected loss validated before finish returns; lone nodes fixed by
flip/split before any deletion; obc/obcfinal check that the open string still covers the supplied line;
3 `protected_loss_m2` = geometric loss without netting; every autofinish patch (both branches and fallback round)
reverted above `finishing.auto.protected_loss_tol_m2`; 4 fallback keeps and validates the OBC; 5 whole fallback
trial compared per node against max(original deviation, `obc_perp_max_deg`).
Library defaults of all protected-loss tolerances are 0 m2. The VR v3 recipe (VR-FVCOM 7fb58cf) sets explicit
tolerances (finish/polish 400,000 m2 per step, autofinish 1,500 m2 per patch) because the optimizer slides bank nodes
along the shoreline; with 0 everywhere VR reaches 20/23. The tolerance VALUES are an owner decision pending, not a
defect; whether the code applies them correctly is in scope.
Round 4 (4 major): all verified correct (tests `test_r4_1`..`test_r4_4` fail on 9a3ec4d, pass on c292eff), fixed in
f076102 (+ c292eff). Triage in `docs/newregion-review-20261116.md`, "Round 4": 1 OBC aftercare refuses to delete an
R4 / {1,2,1} / pinch element holding protected water (stage fails naming the location), and the whole OBC assignment
may not lose protected water (`obc.protected_loss_tol_m2`, default 0); 2 the open-line guard needs only the domain
water polygon (persisted as a path; missing -> loud failure); `open_cover` / `raw_open_cover` persisted in provenance
and reloaded; 3 reference is the wet extent of the supplied line (its part on the domain water outline); each string
end within `obc.end_tol_m` (default 250 m) of it, every node within `max_move_m` of the line; each stage persists its
coverage and later stages may not retreat by more than 1 m; 4 for a region recipe (polygon build with a line, or
protected water) siteops keeps and validates the OBC; Tokyo-style recipes take the legacy path. Follow-up (c292eff):
coincident (< 1 mm) OBC nodes, which polish later welds, count as one node in the string validation.
Evidence: Tokyo job 7592486 (with figures) md5-identical to base1b at all stages; VR v3 job 7592487 QA 23/23, open
line covered 73-183,697 m of the wet 0-183,697 m at obc and obcfinal, coverage guard pass.

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
Three major findings remain. No files were modified.

Verification: 28 existing test functions passed using in-memory execution and mocks. The full filesystem-writing suite was not run. All six Tokyo meshes matched `base1b` across the 12 saved `head*`/`new*` runs. The saved c292eff VR result has NP 27,521 / NE 48,989, 82 OBC nodes, 23/23 passing gates, and persisted coverage of 73.04–183,697.17 m.

**Previous findings**

| Finding | Status | Reason |
|---|---|---|
| R1 #1: Shoreline nodes unprotected | RESOLVED | Boundary-inclusive node masks and triangle overlap cover them. |
| R1 #2: Lone-node removal loses protection/OBC | RESOLVED | Flip/split preserves coverage; deletion guards retain protected/OBC carriers. |
| R1 #3: Unrestored OBC trimming | PARTIAL | Restoration is required, but coverage can retreat locally; finding 2. |
| R1 #4: Band restoration corrupts depths | RESOLVED | New-node depths interpolate over the original triangulation. |
| R1 #5: Expected drops count elements | RESOLVED | Classification uses intersection area. |
| R1 #6: Coarse raster-edge strip | RESOLVED | Interpolation clamps to cell centres. |
| R1 #7: Fallback accepts new failed gates | RESOLVED | Explicit gate rejection and whole-trial perpendicularity checks apply. |
| R1 #8: Depths discard later stages | RESOLVED | Both depth products consume `latest_mesh`. |
| R1 #9: Projected water gets geographic filtering | RESOLVED | Filtering transforms the bbox into the source CRS. |
| R1 #10: Prep ignores resolved CRS | RESOLVED | Explicit and differing legacy overrides reach prep. |
| R1 #11: Stale raw CRS | RESOLVED | Every build records its source CRS; finish uses it. |
| R2 #1: Vertex-only protection misses overlap | RESOLVED | Cap, component and lone-node checks use element overlap. |
| R2 #2: Siteops erodes protected water | RESOLVED | Boundary collapses are guarded and geometric loss is checked. |
| R2 #3: Filtering drops protected components | RESOLVED | Every component outside the actual retained set is measured. Finding 1 concerns missing water within a retained component. |
| R2 #4: Repair threshold not forwarded | RESOLVED | Threshold, policy and tolerances reach the engine. |
| R2 #5: Repair accounting misses overlap | RESOLVED | The explicit repair loop measures intersection area and unexpected loss. |
| R2 #6: Band restore hides C4 failures | RESOLVED | Every affected adjacent-element pair is compared. |
| R2 #7: Geographic raw figures reprojected | RESOLVED | Figures honour the recorded source CRS. |
| R2 #8: Legacy CRS override omitted | RESOLVED | A differing legacy override reaches prep and endpoint extension. |
| R3 #1: Fine protected component skipped | RESOLVED | Measurement follows retained-component selection. |
| R3 #2: Welding/optimization loses protection or OBC | PARTIAL | Weld, optimizer and structural guards work; deferred coverage validation retains finding 2. |
| R3 #3: Autofinish bypasses protection | RESOLVED | Both patch branches and fallback enforce geometric loss tolerances. |
| R3 #4: Fallback loses unprotected OBC | RESOLVED | Fallback independently preserves and validates it. |
| R3 #5: Perpendicularity misses whole fallback | RESOLVED | The complete trial is compared per node against the original limit. |
| R4 #1: OBC aftercare deletes protected elements | RESOLVED | Direct deletions are checked and whole-assignment loss is bounded. |
| R4 #2: Resume disables coverage guard | RESOLVED | Structured coverage persists; missing domain metadata fails loudly. |
| R4 #3: Coverage authorizes artificial walls | PARTIAL | Wet-extent validation works, but gains at one end conceal retreat at the other; finding 2. |
| R4 #4: Standalone siteops clears OBC | PARTIAL | Fixed for polygon/protected recipes; non-Tokyo prep recipes remain affected; finding 3. |

**Findings**

1. **Major — Protected water can disappear inside a retained generator component without triggering protection.**

   **Location:** [mesh_engine/polygon.py:382](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/polygon.py:382), [polygon.py:447](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/polygon.py:447).

   Protection checks components selected for dropping and deletions made by the explicit boundary-repair loop. It never checks how much protected water the generator represented, or geometric loss caused by cleaning within a retained component.

   **Reproduction:** In EPSG:32631, use:

   ```python
   sea = box(300000, 640000, 301000, 640500)
   channel = box(300470, 640500, 300530, 641500)
   water = sea.union(channel)
   ```

   Build with background size 100 m, raster spacing 25 m, the sea’s southern edge as OBC, `keep_polygons=channel`, `max_iter=30`, `seed=0`, no `expected_drop`, and `max_protected_ear_area_km2=0`.

   Actual measurements:

   | Measurement | Result |
   |---|---:|
   | Supplied protected channel | 60,000 m² |
   | Represented immediately after generation | 2,979.25 m² |
   | Represented after OM2D cleaning | 2,553.29 m² |
   | Additional protected loss during cleaning | 425.95 m² |
   | Reported protected repair loss | **0 m²** |
   | Returned mesh | NP 47 / NE 71 |

   Both component records say `kept: true`; the build returns successfully with **95.74% of the protected channel absent**. The subsequent OBC assignment also returns, reporting zero additional protected loss. Separate QA reports three quality-gate failures, but does not identify the missing protected extent.

   **Fix:** Validate represented protected coverage against `water ∩ keep_polygons`, accounting for `expected_drop` and explicit coverage tolerances. Also measure geometric protected loss across OM2D cleaning and enforce the configured loss policy before returning. Component-drop checks alone cannot enforce this invariant.

2. **Major — The round-4 coverage guard nets endpoint retreat against extension at the opposite end.**

   **Location:** [cli/pipeline.py:772](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:772).

   The persisted record contains both endpoint gaps, but the retreat check compares only their sum. OBC reclassification can therefore turn previously open water into land while passing the advertised no-retreat guard.

   **Reproduction:** On a 400 m straight wet boundary with 100 m node spacing, call `_check_open_cover` for these consecutive boundary strings:

   | Stage | Covered interval | Endpoint gaps | Total uncovered |
   |---|---|---|---:|
   | `obc` | 100–400 m | 100 / 0 m | 100 m |
   | `obcfinal` | 0–300 m | 0 / 100 m | 100 m |

   Persist the first result as `open_cover` before checking the second. **Both calls succeed with default tolerances**, although the eastern 100 m of previously open boundary becomes land.

   This defect is introduced by the round-4 retreat check.

   **Fix:** Compare the two clipped coverage endpoints independently: the lower endpoint may not advance by more than 1 m, and the upper endpoint may not retreat by more than 1 m. Extension elsewhere must not cancel retreat.

3. **Major — New-region prep recipes still silently lose their OBC in standalone siteops.**

   **Location:** [cli/pipeline.py:850](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:850), [pipeline.py:714](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:714).

   The round-4 fix identifies a region only through the polygon engine or `protect_polygons`. The documented alternative using `prep`, regional water files and a non-Tokyo CRS takes the legacy branch when protection is absent. That branch clears boundary lists.

   **Reproduction:** Run `_stage_siteops` with:

   - `crs: 32631`;
   - a Ghana `prep.bbox` and `prep.obc_line`;
   - `build.engine: oceanmesh`;
   - `obc: {}` and `siteops: {passes: 1}`;
   - no `protect_polygons`.

   Using the existing 5 × 4 lattice fixture and five southern OBC nodes, with file I/O mocked and real site-operator execution:

   ```text
   _check_recipe: accepted
   _region_obc: False
   siteops: C1=0, C4=0; no edits
   OBC before: [[0, 1, 2, 3, 4]]
   OBC after: []
   stage: returned successfully
   ```

   With `depths.product: none`, this boundary-free siteops mesh becomes `final_mesh`.

   **Fix:** Preserve and validate the existing OBC for new-region prep recipes too. Restrict the compatibility branch to legacy Tokyo recipes, retaining the measured Tokyo behaviour without classifying every non-polygon recipe as legacy.

## Verdict

VERDICT: FAIL (0 blocker, 3 major, 0 minor, 0 nit)
````

### Triage

Fixes, marked `review newregion r5 #N`:
- `2ebabd7`: the 3 findings.
- docs commit: user guide and changelog.

Regression tests `test_r5_1`..`test_r5_3` fail on the reviewed code c292eff (`PYTHONPATH` = `git archive c292eff`) and pass on 2ebabd7.

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes. `test_r5_1`, the reviewer's 60 m channel at 100 m: the reviewed code returned. The fixed code measures the OM2D clean loss as 426 m2, matching the reviewer. | yes | 2ebabd7: the polygon build compares the protected water to carry (water ∩ `protect_polygons` − `expected_drop`) with what the final raw mesh represents (tiled triangle intersections), and measures the geometric protected loss of the OM2D clean. Strict defaults: `protected_coverage_tol_m2` / `_tol_frac` 0, `protected_clean_loss_tol_m2` 0, `protected_coverage_policy: fail`. Recorded as `raw_protected_coverage` (persisted). |
| 2 | major | yes. `test_r5_2`, the reviewer's 100-400 then 0-300 m case: the reviewed code passed both. | yes | 2ebabd7: each end is compared separately (lower may not advance, upper may not retreat, 1 m) |
| 3 | major | yes. `test_r5_3`: a prep-based Ghana recipe (crs 32631, no protection) lost its open boundary in siteops on the reviewed code. | yes | 2ebabd7: `_legacy_siteops`. Legacy means `siteops.legacy_component_path` if given; else no `crs`, EPSG 32654 (from `finish.utm_epsg` or the bbox), no polygon build, no protection. Every other recipe keeps and validates its open boundary. Tokyo v6 is legacy (test) and stays md5-identical. |

No rebuttals.

#### VR v3 numbers for finding 1 (the defaults are strict; VR sets explicit tolerances)

Measured in job 7592604 (`report` policy, build + obc only):
- protected water outside the expected drops: 377.11 km2;
- represented by the generated mesh: 373.64 km2, so 3.475 km2 (0.92 %) missing (bank chords, the boundary repair, narrow bits);
- removed by the OM2D clean: 0.544 km2.

The VR v3 recipe sets `protected_coverage_tol_frac: 0.01` and `protected_clean_loss_tol_m2: 600000` (owner to review). With the strict defaults the VR build would fail.

Checks after the fixes:
- `pytest -q`: 1360 passed, 9 skipped. `ruff`: clean.
- Tokyo (job 7592600, 2ebabd7, figures included): md5-identical to base1b at all six stages.
- VR v3 (job 7592660, 2ebabd7 + docs):
  - QA 23/23; NP 27,521 / NE 48,989; 82 open-boundary nodes.
  - Open line covered 73-183,697 m of the wet 0-183,697 m at obc and obcfinal.
  - Coverage guard passed; depth check passed.
  - Min edge channel / lagoon 62.0 / 66.0 m; dt_E 9.72 s.

## Round 6

### Prompt

```markdown
# Review request, round 6: new-region generalisation of fmesh-pipeline (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT modify files. You may run read-only
commands, python in memory, mocks and small synthetic meshes (keep each computation under a few minutes). Answer
in English as Markdown.

## Goal
Correctness and robustness of the target below under realistic use. Report every defect you can substantiate in
the target, of any severity.

## Target (concrete)
Commit range `1988f22..8f936bb` (32d211b round-1 fixes; 670172c, 5c78274, 448dbe5 round-2 fixes; ef1e8b3, 60201fa, 9a3ec4d round-3 fixes; f076102, c292eff round-4 fixes; 2ebabd7, 8f936bb round-5 fixes). Files:
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
Round 1 (10 major, 1 minor): fixed in 32d211b; round 2 judged 8 RESOLVED and 3 PARTIAL (#1, #2, #10), which are
covered by round-2 findings below.
Round 2 (6 major, 2 minor): all verified correct (tests `test_r2_1`..`test_r2_8` fail on 32d211b, pass on 448dbe5)
and fixed in 670172c (+ 5c78274, 448dbe5). Triage in `docs/newregion-review-20261116.md`, "Round 2":
1 `protected_element_mask`: protected if a protected vertex or positive overlap area (tiled geometry), used for cap
deletion, component drops and lone-node removal; 2 collapses with both ends on the boundary next to protected water
are refused before each collapse (short-edge sweep and site operators), and any protected-coverage loss > 1 m2 fails
the stage; 3 a dropped component of any size may not take protected water outside `expected_drop` beyond
max(1 m2, `protected_drop_tol_frac` (default 1 %) of its protected area) unless `protected_drop_policy: report`;
new `build.expected_drop_unresolved_m: w` adds water not reached by the open-boundary component of the water opened
at width w (VR v3 uses 150 m, VR-FVCOM commit 36ff15f); 4 `max_protected_ear_area_km2` and the drop policy and
tolerances are forwarded by the pipeline; 5 repair loss counted by intersection area (total and outside
expected_drop; the threshold applies to the latter); 6 C4 compared pair by pair in band restore; 7 raw mesh
transformed from its recorded CRS, none for EPSG:4326; 8 legacy `finish.utm_epsg` forwarded to prep when it differs
from prep's own zone.
Round 3 (5 major): all verified correct (tests `test_r3_1`..`test_r3_5` fail on 448dbe5, pass on 9a3ec4d), fixed in
ef1e8b3 (+ 60201fa, 9a3ec4d). Triage and a trade-off table in `docs/newregion-review-20261116.md`, "Round 3":
1 retained components decided first, every dropped component measured; 2 weld guard (OBC nodes only onto coincident
< 1 mm nodes, refused when protected water is lost), optimizer guard (protected bank nodes restored, round refused
above `protected_loss_tol_m2`), OBC string and protected loss validated before finish returns; lone nodes fixed by
flip/split before any deletion; obc/obcfinal check that the open string still covers the supplied line;
3 `protected_loss_m2` = geometric loss without netting; every autofinish patch (both branches and fallback round)
reverted above `finishing.auto.protected_loss_tol_m2`; 4 fallback keeps and validates the OBC; 5 whole fallback
trial compared per node against max(original deviation, `obc_perp_max_deg`).
Library defaults of all protected-loss tolerances are 0 m2. The VR v3 recipe (VR-FVCOM 7fb58cf) sets explicit
tolerances (finish/polish 400,000 m2 per step, autofinish 1,500 m2 per patch) because the optimizer slides bank nodes
along the shoreline; with 0 everywhere VR reaches 20/23. The tolerance VALUES are an owner decision pending, not a
defect; whether the code applies them correctly is in scope.
Round 4 (4 major): all verified correct (tests `test_r4_1`..`test_r4_4` fail on 9a3ec4d, pass on c292eff), fixed in
f076102 (+ c292eff). Triage in `docs/newregion-review-20261116.md`, "Round 4": 1 OBC aftercare refuses to delete an
R4 / {1,2,1} / pinch element holding protected water (stage fails naming the location), and the whole OBC assignment
may not lose protected water (`obc.protected_loss_tol_m2`, default 0); 2 the open-line guard needs only the domain
water polygon (persisted as a path; missing -> loud failure); `open_cover` / `raw_open_cover` persisted in provenance
and reloaded; 3 reference is the wet extent of the supplied line (its part on the domain water outline); each string
end within `obc.end_tol_m` (default 250 m) of it, every node within `max_move_m` of the line; each stage persists its
coverage and later stages may not retreat by more than 1 m; 4 for a region recipe (polygon build with a line, or
protected water) siteops keeps and validates the OBC; Tokyo-style recipes take the legacy path. Follow-up (c292eff):
coincident (< 1 mm) OBC nodes, which polish later welds, count as one node in the string validation.
Round 5 (3 major): all verified correct (tests `test_r5_1`..`test_r5_3` fail on c292eff, pass on 2ebabd7), fixed in
2ebabd7. Triage in `docs/newregion-review-20261116.md`, "Round 5": 1 the polygon build compares the protected water it
should carry (domain water ∩ protect_polygons minus expected_drop) with what the final raw mesh covers, and measures
the protected water removed by the OM2D clean; strict defaults (`protected_coverage_tol_m2`,
`protected_coverage_tol_frac`, `protected_clean_loss_tol_m2` = 0, `protected_coverage_policy: fail`); result persisted
as `raw_protected_coverage`; 2 each end of the open string compared separately against the earlier stage (1 m);
3 legacy siteops path only when `siteops.legacy_component_path` says so, or no `crs`, resolves to EPSG 32654, no
polygon build and no protection; every other recipe keeps and validates its OBC.
VR v3 (VR-FVCOM da4e566) sets `protected_coverage_tol_frac: 0.01` and `protected_clean_loss_tol_m2: 600000` after a
measurement (377.11 km2 protected, 373.64 km2 represented, 0.92 % missing, 0.544 km2 removed by clean). As before, the
tolerance VALUES are an owner decision pending, not a defect; whether the code applies them correctly is in scope.
Evidence: Tokyo job 7592600 (with figures) md5-identical to base1b at all stages; VR v3 job 7592660 QA 23/23, open
line covered 73-183,697 m of the wet 0-183,697 m at obc and obcfinal, coverage guard pass.

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
```

### Answer (verbatim)

One major finding remains. No files were modified by this review.

**Previous findings**

| Finding | Status | Reason |
|---|---|---|
| R1 #1: Shoreline nodes unprotected | RESOLVED | Boundary-inclusive masks and triangle overlap protect them. |
| R1 #2: Lone-node removal loses protection/OBC | RESOLVED | Flip/split repair preserves coverage; deletion guards retain protected/OBC carriers. |
| R1 #3: Unrestored polygon OBC trimming | RESOLVED | Restoration is required, and OBC stages validate the intended wet extent. |
| R1 #4: Band restoration corrupts depths | RESOLVED | New-node depths interpolate over the original triangulation. |
| R1 #5: Expected drops count elements | RESOLVED | Classification uses intersection area. |
| R1 #6: Coarse raster-edge strip | RESOLVED | Interpolation clamps to cell centres. |
| R1 #7: Fallback accepts newly failing gates | RESOLVED | Explicit gate rejection and whole-trial perpendicularity checks apply. |
| R1 #8: Depths discard later stages | RESOLVED | Both depth products consume `latest_mesh`. |
| R1 #9: Projected water receives geographic filtering | RESOLVED | Filtering transforms the bbox into the source CRS. |
| R1 #10: Prep ignores resolved CRS | RESOLVED | Explicit and differing legacy overrides reach prep. |
| R1 #11: Stale raw CRS | RESOLVED | Every build records its source CRS; finish uses it. |
| R2 #1: Vertex-only protection misses overlap | RESOLVED | Cap, component and lone-node checks include triangle overlap. |
| R2 #2: Siteops erodes protected water | RESOLVED | Boundary collapses are guarded and geometric loss is checked. |
| R2 #3: Filtering drops protected components | RESOLVED | Every component outside the actual retained set is measured. |
| R2 #4: Repair threshold not forwarded | RESOLVED | Threshold, policy and tolerances reach the engine. |
| R2 #5: Repair accounting misses overlap | RESOLVED | Intersection areas measure total and unexpected loss. |
| R2 #6: Band restore hides C4 failures | RESOLVED | Affected adjacent-element pairs are checked separately. |
| R2 #7: Geographic raw figures reprojected | RESOLVED | Figures honour the recorded source CRS. |
| R2 #8: Legacy CRS override omitted | RESOLVED | A differing override reaches prep and endpoint extension. |
| R3 #1: Fine protected component skipped | RESOLVED | Measurement follows retained-component selection. |
| R3 #2: Welding/optimization loses protection or OBC | RESOLVED | Weld/optimizer guards, structural validation and subsequent line checks address the reported cases. |
| R3 #3: Autofinish bypasses protection | RESOLVED | Both patch branches and fallback enforce geometric loss tolerances. |
| R3 #4: Fallback loses unprotected OBC | RESOLVED | Fallback independently preserves and validates it. |
| R3 #5: Perpendicularity misses whole fallback | RESOLVED | The complete trial is compared against the original per-node deviations. |
| R4 #1: OBC aftercare deletes protected elements | RESOLVED | Direct deletions are guarded; whole-assignment loss is bounded. |
| R4 #2: Resume disables coverage guard | RESOLVED | Structured coverage persists; missing domain metadata fails loudly. |
| R4 #3: Coverage authorizes artificial walls | RESOLVED | Validation uses the intended wet extent and independently checks endpoint retreat. |
| R4 #4: Standalone siteops clears OBC | PARTIAL | Preservation still excludes some new-region recipes; finding 1 below. |
| R5 #1: Protected water missing inside retained components | RESOLVED | Final represented coverage and cleaning loss are measured and enforced. |
| R5 #2: Endpoint retreat cancelled by opposite extension | RESOLVED | Each endpoint is compared separately. |
| R5 #3: New-region prep recipes lose OBC | PARTIAL | CRS inference into UTM zone 54 still selects the destructive compatibility branch. |

1. **Major — UTM zone 54 is still treated as proof that a recipe is Tokyo.**

   **Location:** [cli/pipeline.py:734](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:734), [pipeline.py:886](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:886).

   `_legacy_siteops()` returns true whenever a recipe has no explicit `crs`, resolves to EPSG:32654, uses a non-polygon engine and has no protection. That includes ordinary new-region prep recipes elsewhere in zone 54. Siteops then clears their boundary lists without validation.

   **Reproduction:** Use a prep/oceanmesh recipe for Sendai with `prep.bbox: [141.0, 38.2, 141.6, 38.8]`, a supplied OBC line, inferred CRS, no `protect_polygons`, and the stages `obc → siteops → depths: {product: none} → qa`, without `obcfinal`.

   I exercised `_stage_siteops()` using the existing `_r3_lattice()` fixture translated into EPSG:32654 near 141.2°E, 38.4°N, with southern OBC `[0,1,2,3,4]`. File I/O was mocked; site operators and QA were real:

   ```text
   _check_recipe:       accepted
   resolved CRS:        EPSG:32654
   _legacy_siteops:     True
   siteops edits:       none
   elements:            32 → 32
   open strings:        [[0,1,2,3,4]] → []
   QA failed gates:     []
   depths.product none: propagates the siteops mesh as final_mesh
   ```

   Thus a realistic new-region recipe silently produces an accepted mesh with its open boundary removed. This is an incomplete round-5 fix.

   **Fix:** Restrict automatic compatibility behaviour to identifiable shipped legacy Tokyo recipes, retaining explicit `siteops.legacy_component_path` as an override. Default other recipes to OBC preservation, including recipes that infer EPSG:32654. Add a same-zone, non-Tokyo regression test.

**Verification**

Executed **36 selected test functions from committed `8f936bb`** using Python `-B` and in-memory I/O where necessary; all passed. Repeated small polygon generation with the same seed produced identical points, elements and OBC strings. Initial Matplotlib cache imports failed under the read-only sandbox; an in-memory bootstrap allowed the checks to proceed. GDAL emitted a PROJ diagnostic during otherwise successful runs.

All six stored Tokyo stage meshes in **`head1`–`head11`, `new1`, and `new2`** are md5-identical to `base1b`. Saved VR artifacts confirm **NP 27,521 / NE 48,989 / 82 OBC nodes**, zero QA gate failures, and the reported coverage measurements. Full Tokyo and VR generation were not rerun. No additional regression introduced by the fixes was substantiated.

## Verdict

VERDICT: FAIL (0 blocker, 1 major, 0 minor, 0 nit)

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes: `_legacy_siteops` returned True for any recipe without `crs` resolving to EPSG:32654 with a non-polygon engine and no protection, including a Sendai bbox; test_r6_1 reproduces it (OBC cleared) | yes | 3a49832: only `siteops.legacy_component_path: true` selects the legacy path; the five shipped Tokyo recipes set it; test_r6_1. Tokyo v6 job 7592967 md5-identical to base1b at all six stages; ci_local job 7592971 PASS |

Note (not a review finding): GitHub CI run 43778047151 (pushed by the owner at 2ebabd7) failed 12 tests that need the oceanmesh fork; fixed in 3a49832 with `needs_oceanmesh` marks, and `tools/ci_local/run_ci_local.sh` now reproduces CI before every commit (CLAUDE.md).

## Round 7

### Prompt

```markdown
# Review request, round 7: new-region generalisation of fmesh-pipeline (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT modify files. You may run read-only
commands, python in memory, mocks and small synthetic meshes (keep each computation under a few minutes). Answer
in English as Markdown.

## Goal
Correctness and robustness of the target below under realistic use. Report every defect you can substantiate in
the target, of any severity.

## Target (concrete)
Commit range `1988f22..3a49832` (32d211b round-1 fixes; 670172c, 5c78274, 448dbe5 round-2 fixes; ef1e8b3, 60201fa, 9a3ec4d round-3 fixes; f076102, c292eff round-4 fixes; 2ebabd7, 8f936bb round-5 fixes; 3a49832 round-6 fix and CI test marks). Files:
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
Round 1 (10 major, 1 minor): fixed in 32d211b; round 2 judged 8 RESOLVED and 3 PARTIAL (#1, #2, #10), which are
covered by round-2 findings below.
Round 2 (6 major, 2 minor): all verified correct (tests `test_r2_1`..`test_r2_8` fail on 32d211b, pass on 448dbe5)
and fixed in 670172c (+ 5c78274, 448dbe5). Triage in `docs/newregion-review-20261116.md`, "Round 2":
1 `protected_element_mask`: protected if a protected vertex or positive overlap area (tiled geometry), used for cap
deletion, component drops and lone-node removal; 2 collapses with both ends on the boundary next to protected water
are refused before each collapse (short-edge sweep and site operators), and any protected-coverage loss > 1 m2 fails
the stage; 3 a dropped component of any size may not take protected water outside `expected_drop` beyond
max(1 m2, `protected_drop_tol_frac` (default 1 %) of its protected area) unless `protected_drop_policy: report`;
new `build.expected_drop_unresolved_m: w` adds water not reached by the open-boundary component of the water opened
at width w (VR v3 uses 150 m, VR-FVCOM commit 36ff15f); 4 `max_protected_ear_area_km2` and the drop policy and
tolerances are forwarded by the pipeline; 5 repair loss counted by intersection area (total and outside
expected_drop; the threshold applies to the latter); 6 C4 compared pair by pair in band restore; 7 raw mesh
transformed from its recorded CRS, none for EPSG:4326; 8 legacy `finish.utm_epsg` forwarded to prep when it differs
from prep's own zone.
Round 3 (5 major): all verified correct (tests `test_r3_1`..`test_r3_5` fail on 448dbe5, pass on 9a3ec4d), fixed in
ef1e8b3 (+ 60201fa, 9a3ec4d). Triage and a trade-off table in `docs/newregion-review-20261116.md`, "Round 3":
1 retained components decided first, every dropped component measured; 2 weld guard (OBC nodes only onto coincident
< 1 mm nodes, refused when protected water is lost), optimizer guard (protected bank nodes restored, round refused
above `protected_loss_tol_m2`), OBC string and protected loss validated before finish returns; lone nodes fixed by
flip/split before any deletion; obc/obcfinal check that the open string still covers the supplied line;
3 `protected_loss_m2` = geometric loss without netting; every autofinish patch (both branches and fallback round)
reverted above `finishing.auto.protected_loss_tol_m2`; 4 fallback keeps and validates the OBC; 5 whole fallback
trial compared per node against max(original deviation, `obc_perp_max_deg`).
Library defaults of all protected-loss tolerances are 0 m2. The VR v3 recipe (VR-FVCOM 7fb58cf) sets explicit
tolerances (finish/polish 400,000 m2 per step, autofinish 1,500 m2 per patch) because the optimizer slides bank nodes
along the shoreline; with 0 everywhere VR reaches 20/23. The tolerance VALUES are an owner decision pending, not a
defect; whether the code applies them correctly is in scope.
Round 4 (4 major): all verified correct (tests `test_r4_1`..`test_r4_4` fail on 9a3ec4d, pass on c292eff), fixed in
f076102 (+ c292eff). Triage in `docs/newregion-review-20261116.md`, "Round 4": 1 OBC aftercare refuses to delete an
R4 / {1,2,1} / pinch element holding protected water (stage fails naming the location), and the whole OBC assignment
may not lose protected water (`obc.protected_loss_tol_m2`, default 0); 2 the open-line guard needs only the domain
water polygon (persisted as a path; missing -> loud failure); `open_cover` / `raw_open_cover` persisted in provenance
and reloaded; 3 reference is the wet extent of the supplied line (its part on the domain water outline); each string
end within `obc.end_tol_m` (default 250 m) of it, every node within `max_move_m` of the line; each stage persists its
coverage and later stages may not retreat by more than 1 m; 4 for a region recipe (polygon build with a line, or
protected water) siteops keeps and validates the OBC; Tokyo-style recipes take the legacy path. Follow-up (c292eff):
coincident (< 1 mm) OBC nodes, which polish later welds, count as one node in the string validation.
Round 5 (3 major): all verified correct (tests `test_r5_1`..`test_r5_3` fail on c292eff, pass on 2ebabd7), fixed in
2ebabd7. Triage in `docs/newregion-review-20261116.md`, "Round 5": 1 the polygon build compares the protected water it
should carry (domain water ∩ protect_polygons minus expected_drop) with what the final raw mesh covers, and measures
the protected water removed by the OM2D clean; strict defaults (`protected_coverage_tol_m2`,
`protected_coverage_tol_frac`, `protected_clean_loss_tol_m2` = 0, `protected_coverage_policy: fail`); result persisted
as `raw_protected_coverage`; 2 each end of the open string compared separately against the earlier stage (1 m);
3 legacy siteops path only when `siteops.legacy_component_path` says so, or no `crs`, resolves to EPSG 32654, no
polygon build and no protection; every other recipe keeps and validates its OBC.
VR v3 (VR-FVCOM da4e566) sets `protected_coverage_tol_frac: 0.01` and `protected_clean_loss_tol_m2: 600000` after a
measurement (377.11 km2 protected, 373.64 km2 represented, 0.92 % missing, 0.544 km2 removed by clean). As before, the
tolerance VALUES are an owner decision pending, not a defect; whether the code applies them correctly is in scope.
Round 6 (1 major): verified correct and fixed in 3a49832: the legacy siteops path is selected only by the explicit
`siteops.legacy_component_path: true` (no CRS heuristic); the five shipped Tokyo recipes set it; test `test_r6_1`
(a zone-54 Sendai-like prep recipe keeps its OBC). The same commit adds `@pytest.mark.needs_oceanmesh` to 12 tests
that GitHub CI (no oceanmesh fork) failed, and `tools/ci_local/run_ci_local.sh` (CI reproduced locally; in scope
only as far as it reproduces `.github/workflows/ci.yml` faithfully).
Evidence: Tokyo job 7592967 md5-identical to base1b at all stages; ci_local job 7592971 PASS (CI mode 1329 passed,
36 skipped; lab mode 1361 passed); VR v3 job 7592660 QA 23/23 (unaffected: v3 has protection and a polygon build).

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
```

### Answer (verbatim)

Two findings remain: **one major and one minor**. No files were modified.

**Previous findings**

| Finding | Status | Reason |
|---|---|---|
| R1 #1: Shoreline nodes unprotected | RESOLVED | Protection includes outlines and triangle overlap. |
| R1 #2: Lone-node removal erases protection/OBC | RESOLVED | Flip/split repair preserves coverage; deletion guards protect carriers. |
| R1 #3: Polygon OBC trimming lacks restoration | PARTIAL | Restoration is required, but subsequent polish can undo it without validation; finding 1. |
| R1 #4: Band restoration corrupts depths | RESOLVED | New-node depths interpolate over the original triangulation. |
| R1 #5: Expected drops count elements | RESOLVED | Classification uses intersection area. |
| R1 #6: Raster-edge sizing becomes coarse | RESOLVED | Interpolation clamps to cell centres. |
| R1 #7: Fallback accepts newly failing gates | RESOLVED | Explicit gate rejection and perpendicularity guards apply. |
| R1 #8: Depths discard later stages | RESOLVED | Both products consume and propagate `latest_mesh`. |
| R1 #9: Projected water receives geographic filtering | RESOLVED | Filtering transforms the bbox into the source CRS. |
| R1 #10: Prep ignores resolved CRS | RESOLVED | Explicit and differing legacy overrides reach prep. |
| R1 #11: Stale raw CRS survives rebuild | RESOLVED | Every build records its source CRS. |
| R2 #1: Vertex-only protection misses overlap | RESOLVED | Cap, component and lone-node checks include overlap. |
| R2 #2: Siteops erodes protected water | RESOLVED | Boundary collapses are guarded; geometric loss is checked. |
| R2 #3: Filtering drops protected components | RESOLVED | Dropped protected area outside expected drops is bounded. |
| R2 #4: Repair controls not forwarded | RESOLVED | Thresholds, policy and tolerances reach the engine. |
| R2 #5: Repair accounting misses overlap | RESOLVED | Intersection areas measure total and unexpected loss. |
| R2 #6: Band restore hides C4 failures | RESOLVED | Affected adjacent-element pairs are checked individually. |
| R2 #7: Geographic raw figures reprojected | RESOLVED | Figures honour the recorded source CRS. |
| R2 #8: Legacy CRS override omitted | RESOLVED | A differing override reaches prep and endpoint extension. |
| R3 #1: Fine protected component skipped | RESOLVED | Measurement follows actual retained-component selection. |
| R3 #2: Welding/optimization loses protection/OBC | PARTIAL | The reported operator cases are guarded, but validation can precede a later boundary regression; finding 1. |
| R3 #3: Autofinish bypasses protection | RESOLVED | Both patch branches and fallback enforce geometric loss tolerances. |
| R3 #4: Fallback loses unprotected OBC | RESOLVED | Fallback preserves and validates it independently. |
| R3 #5: Perpendicularity misses whole fallback | RESOLVED | The complete trial receives the per-node comparison. |
| R4 #1: OBC aftercare deletes protected water | RESOLVED | Direct deletions are guarded; whole-assignment loss is bounded. |
| R4 #2: Resume disables coverage guard | RESOLVED | Coverage persists and missing domain metadata fails loudly. Persistence introduces finding 2. |
| R4 #3: Coverage authorizes artificial walls | RESOLVED | OBC-stage checks use intended wet extent and endpoint retreat. |
| R4 #4: Standalone siteops clears OBC | RESOLVED | Preservation is now the default. |
| R5 #1: Protected coverage missing within retained components | RESOLVED | Represented coverage and cleaning loss are measured and bounded. |
| R5 #2: Opposite extension cancels endpoint retreat | RESOLVED | Each endpoint is compared separately. |
| R5 #3: New-region prep recipes lose OBC | RESOLVED | Compatibility requires the explicit legacy flag. |
| R6 #1: Zone 54 mistaken for Tokyo | RESOLVED | CRS inference no longer selects legacy siteops; the Sendai regression passes. |

**Findings**

1. **Major — An accepted stage sequence can silently restore an earlier, shortened OBC.**

   **Locations:** [pipeline.py:164](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:164), [pipeline.py:949](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:949), [pipeline.py:974](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:974).

   `_check_recipe()` accepts `finish → obc → polish → depths → qa` without `siteops` or `obcfinal`. However, polish selects `finished_mesh` when `siteops_mesh` is absent, bypassing the validated `obc_mesh`. It then returns without checking coverage against the persisted OBC record.

   **Reproduction:** Using `_r3_lattice()`—25 nodes and 32 conforming triangles—supply these valid stage meshes:

   - Finished mesh: southern open string `[1,2,3]`, covering 100–300 m.
   - OBC mesh: restored string `[0,1,2,3,4]`, covering 0–400 m.
   - Persisted `open_cover`: 0–400 m.

   Run actual polish with `optimize_budget_s: 0`, `freeze_obc_band: true`, `max_valence: 8`, and `remove_lone_nodes: true`, followed by depths and pipeline QA. Only filesystem I/O was substituted:

   ```text
   recipe validation: accepted
   polish output OBC: 100–300 m
   depths.product none: promotes polished mesh to final_mesh
   pipeline QA: qa_failures = 0
   ```

   Calling the omitted coverage check rejects the result because both endpoints retreated by 100 m.

   The saved Volta artifacts demonstrate the realistic source-selection risk: `finished` covers approximately **2,885–176,735 m**, whereas `obc` restores **73–183,697 m**.

   **Fix:** Make polish consume the latest preceding mesh, including `obc_mesh`. Validate final OBC coverage after subsequent mesh edits, or require `obcfinal` whenever those edits cannot guarantee preservation. Preserve the explicit Tokyo compatibility path.

2. **Minor — Persisted coverage from an earlier build falsely rejects a legitimate rebuild.**

   **Locations:** [pipeline.py:785](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:785), [pipeline.py:794](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:794), [pipeline.py:1509](/home/pj24001722/ku40000343/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/pipeline.py:1509).

   Loading provenance retains `open_cover` even when domain/build are rerun. The first OBC stage of the new mesh is therefore compared against a different mesh generation. This is a regression introduced by coverage persistence.

   **Reproduction:** Run the existing `_domain()` geometry through actual `domain → build → obc`, with seed 0, 30 iterations and default endpoint tolerance:

   ```text
   background size 300 m: coverage 0–6000 m; PASS
   rebuild at 400 m:      coverage 0–5942.54 m; FAIL
   ```

   The failure reports an upper-end retreat from 6000 to 5943 m. Yet the **57.46 m** gap satisfies `end_tol_m: 250`; the identical 400 m build passes when only the previous generation’s coverage record is removed.

   **Fix:** Associate coverage records with a mesh-generation identity, or invalidate downstream coverage when rebuilding the raw mesh. Retain the persisted comparison when genuinely resuming that generation with `--only`.

**Verification**

Executed all **56 target test functions** using Python `-B`, in-memory vector/mesh I/O substitutions and a GDAL memory raster; all passed. Initial unadapted attempts encountered sandbox cache/temp-file restrictions. GDAL also emitted a PROJ diagnostic during successful checks.

All six stored Tokyo meshes in `head1`–`head12`, `new1`, and `new2` match `base1b` by md5. Saved Volta QA confirms **27,521 nodes, 48,989 elements, 82 OBC nodes and 23/23 gates passed**. Full Tokyo/Volta generation and full CI were not rerun. No additional defect in the round-6 compatibility fix or CI marks/script was substantiated.

## Verdict

VERDICT: FAIL (0 blocker, 1 major, 1 minor, 0 nit)

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes: `_check_recipe` accepted finish -> obc -> polish without obcfinal, and polish read `finished_mesh` when `siteops_mesh` was absent; test_r7_1 fails on a17ab77 | yes | 5f97c9d: obcfinal required after siteops/polish when the recipe has an open line; polish starts from `obc_mesh` without siteops; test_r7_1 |
| 2 | minor | yes: test_r7_2 (full rerun after a persisted wider coverage record) fails on a17ab77 (job 7593202) | yes | 5f97c9d: domain/build drop the persisted `open_cover`; `--only obc` still compares; test_r7_2 |

ci_local job 7593203 PASS before the commit.
