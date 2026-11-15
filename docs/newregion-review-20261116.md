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
