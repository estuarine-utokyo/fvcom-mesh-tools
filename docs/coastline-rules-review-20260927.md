# Review: coastline rules for local refinement (2026-09-27)

Reviewer: gpt-6-astra (Codex CLI), run by the codex-review loop. Scope:
`git diff 0c5d9b3..HEAD` at the start of round 1 (commit 59f3a6d).

## Round 1

### Answer (verbatim)

Read-only review completed; no files changed and no mesh generation or batch jobs ran. The full pytest invocation failed before collection because the sandbox has no writable temporary directory. A second attempt encountered Matplotlib’s writable-cache requirement. With `--capture=sys --noconftest -p no:cacheprovider`, **22 focused geometry tests passed**. The findings below also use direct in-memory reproductions and mocked driver blocks.

1. **Major — Default rim repair violates `hires.coastline: preserve`.**  
   **Location:** [420_local_refine.py:1215](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1215)

   The earlier blunting block explicitly discards changes under `preserve`, but this unconditional default repair applies them again. Free coastline points have `pfix_base == -1` even in preserve mode. Reproduction: the cyclic rim `[(0,0),(200,600),(0,900),(-200,600)]`, with all points free and `h=30`, loses its acute corner and **270 m² of water**. The hires branch also sets the departure tolerance to infinity, so that check cannot enforce preservation.

   **Fix:** Restrict geometry-changing rim repair and its retry to `resolve`, or implement a preservation mode that permits only changes proven to leave the boundary polyline identical.

2. **Major — The continuous-width rule can close a lagoon entrance containing resolvable water beyond it.**  
   **Location:** [patch.py:823](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:823)

   Touching one land component does not establish that a water feature is an unresolvable dead end. I reproduced this through `filter_shoreline_local`, using footprint `box(-2000,-2000,2000,2000)` and land equal to that footprint minus the union of:
   `box(-2100,-1500,-400,1500)`, `box(400,-1500,1500,1500)`, and the connecting channel `box(-500,-150,500,150)`.
   With constant `h=225`, `h0=30`, and spacing 30, disabling continuous width leaves **one water component**; enabling it creates **two**, closing the entrance while leaving the large lagoon wet.

   **Fix:** Check water connectivity and whether the closure isolates a basin that accommodates local elements. Preserve such entrances, consistent with the basin check already used by slit repair.

3. **Major — Refused wall pockets can permanently lose their walls.**  
   **Location:** [420_local_refine.py:689](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:689)

   Restoration tests only whether the pocket’s rounded representative point lies inside the remaining hole. A pocket crossing the rim can have substantial water inside while that point lies outside.

   Reproduction: land `box(-2000,-2000,2000,0)`, wall `[(0,0),(0,380),(400,380),(400,150),(40,150)]`, and `h=200`. Pocket closure reports its point at `(200,265)`. For hole `box(-100,100,100,500)`, `island_rings` refuses the pocket with **23,000 m² inside**, but `_lost` is empty. Wall length inside the hole falls from **440 m to 49.99 m**. Subsequent constraints contain neither the refused island nor the removed wall sections.

   **Fix:** Track pocket geometries and acceptance explicitly. Restore removed wall portions wherever the pocket remains water, using geometric intersection rather than representative-point containment.

4. **Minor — Short-edge removal can collapse a valid three-node ring into invalid constraints.**  
   **Location:** [patch.py:1571](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1571)

   The removal operation does not check whether the replacement edge already exists or whether the ring retains three distinct vertices. Calling `rim_repair` on `[(0,0),(10,0),(5,20)]`, cyclic edges, all-free points, and `h=30` raises `ValueError: A linearring requires at least 4 coordinates`. The same failure occurs with this triangle represented as an island inside a larger square.

   **Fix:** Reject vertex removal that creates a duplicate edge or undersized ring. Any deliberate whole-island removal needs a separate, explicit operation.

5. **Minor — A retreat-only or merge-only retry is neither searched nor restored.**  
   **Location:** [420_local_refine.py:2038](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:2038)

   `apply_rim_repair` has already mutated the rim, hole, and wall indices, but the retry condition considers only `n_points_removed` and `n_slits_closed`. It ignores `n_edges_merged` and `n_tips_stepped_back`.

   Executing the actual retry block with a mocked retreat-only repair produced **zero constraint-assembly calls, zero search calls, and unequal rim/constraint coordinates**. The state remains partially updated beside the original best candidate.

   **Fix:** Detect every geometry change, preferably by comparing the returned arrays. Search with the changed constraints or restore the snapshot on every unselected path, including exceptions.

6. **Minor — Empty land can turn open water into artificial land.**  
   **Location:** [patch.py:794](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:794)

   When the raster contains no land cells, the distance transform does not represent distance to actual coastline. Restricting disc centres to the footprint then permits false closures.

   Reproduction: `unresolvable_water(Polygon(), constant_2000, box(0,0,1000,1000), min_h=60, spacing=10)` returns approximately **999,954 m² of land** in an entirely land-free domain.

   **Fix:** Handle empty land explicitly and pad the distance-transform/centre domain sufficiently for discs outside the footprint. Clip the final result to the requested domain.

7. **Minor — Invalid continuous-width parameters can hang rather than fail.**  
   **Location:** [patch.py:801](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:801)

   `levels_ratio=1` leaves `lo` unchanged forever. A small in-memory invocation required interruption after 0.5 seconds. There is also no useful empty-footprint handling: `Polygon()` produces `ValueError: cannot convert float NaN to integer` during raster sizing.

   **Fix:** Validate finite positive spacing, radii, and sampled sizes; require finite `levels_ratio > 1`; and explicitly return an empty result or reject an empty footprint before allocating the raster.

8. **Minor — `seed_search` ignores its argument.**  
   **Location:** [420_local_refine.py:1960](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1960)

   The function iterates the global `seeds` instead of `seed_list`. Running the extracted function with global seeds `[0,1]` and calling `seed_search([17])` attempted **0 and 1**.

   **Fix:** Iterate `seed_list` and add a test where it differs from the global list.

9. **Minor — Provenance imports GPL code from inside the Apache package.**  
   **Location:** [provenance.py:68](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:68)

   `LIBRARIES` includes `oceanmesh`, and `_version` calls `importlib.import_module` for every entry. A mocked import recorder confirmed that ordinary `collect()` imports `oceanmesh`. This violates the explicit repository instruction against importing it from `fvcom_mesh_tools`.

   **Fix:** Obtain distribution versions through `importlib.metadata`, accept versions supplied by the notebook, or inspect the GPL dependency through a subprocess.

10. **Minor — Provenance omits inputs and settings that change the delivered mesh.**  
    **Location:** [420_local_refine.py:323](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:323)

    The record hashes the base depths but omits the three NetCDF products actually sampled for `bathymetry: tokyo_bay`. It also omits external region-geometry file hashes and effective `LR_ONLY_BELOW`/`LR_OBC_GUARD` settings. The former changes mesh repair; the latter changes selection. Consequently, different depth fields or mesh behavior can share the recorded provenance.

    **Fix:** Record resolved paths and hashes for consumed bathymetry and geometry inputs, plus all effective mesh-affecting environment settings.

11. **Minor — Untracked source files are reported as a clean checkout.**  
    **Location:** [provenance.py:59](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:59)

    `git status --untracked-files=no` excludes uncommitted new files. For example, an untracked experimental driver inside an otherwise clean repository can call `collect(code={"tool": __file__})` and receive the tracked commit with an empty `dirty` list, although that commit cannot reproduce the running program.

    **Fix:** Include relevant untracked source files and verify that each reported code entry is tracked at the recorded revision.

12. **Minor — The accepted-seed-only rebuild instructions are invalid for feedback retries.**  
    **Location:** [USER_GUIDE.md:305](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/docs/USER_GUIDE.md:305)

    The guide instructs users to rebuild with only the accepted seed. The retry geometry, however, depends on the best failing seed from the first search. A mock executing the actual retry block showed searches `[0,1]` and `[1]` both accepting seed 1 while using **different repair-focus coordinates**. Thus the accepted seed alone does not determine the repaired domain.

    **Fix:** Require replay of the original seed sequence and effective settings, or save and reload the selected repaired constraints. Record which search pass supplied the accepted candidate.

13. **Minor — Several tests pass without exercising their named operation or guard.**  
    **Locations:** [test_patch.py:2034](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/tests/test_patch.py:2034), [test_patch.py:2124](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/tests/test_patch.py:2124), [test_walls.py:407](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/tests/test_walls.py:407), [test_walls.py:445](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/tests/test_walls.py:445)

    The coast-movement refusal fixture has no edge below the **15 m** short-edge threshold; its purported short edge is approximately **40.11 m**. The cap-merge test actually reports **one point removal and zero merges**. The wall-width test uses a **40 m mouth with `h=30`**, so no closing chord is created. The fine-zone test’s pocket is only **70 m off the coast at `h=200`**, independently failing the 100 m clearance requirement.

    **Fix:** Construct fixtures that reach the intended branch, assert its operation counters, and verify that disabling the specific guard makes the corresponding test fail.

14. **Nit — The documented raster resolution is incorrectly fixed at 10 m.**  
    **Location:** [USER_GUIDE.md:255](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/docs/USER_GUIDE.md:255)

    The implementation passes `spacing=h0/3`, so 10 m applies only to a 30 m target.

    **Fix:** Document “one-third of the finest target size; 10 m for a 30 m target.”

## Verdict

VERDICT: FAIL (0 blocker, 3 major, 10 minor, 1 nit)


### Prompt

```markdown
# Review request, round 1: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (19 commits; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Please
1. A fresh, unrestricted audit of the scope above and everything it touches:
   correctness bugs, edge cases (empty geometry, MultiPolygon, degenerate
   rims, junctions, frozen points, index remapping), determinism and
   reproducibility, error handling, tests that do not test what they claim,
   and whether the docs and recipe comments match the code.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes: preserve path reaches apply_rim_repair (RIM_REPAIR was EXPERIMENTAL-only) | yes | RIM_REPAIR requires coastline resolve; test_rim_repair_is_off_where_the_coastline_is_preserved |
| 2 | major | yes: reviewer's lagoon, 1 -> 2 water components with continuous_width | yes | entrance check (_splits_water on the piece grown by 2 cells); test_continuous_width_leaves_a_lagoon_entrance_open |
| 3 | major | yes: pocket crossing hole, point outside, walls not restored | yes | restore judged on pocket polygon (wkt in report); test_a_refused_wall_pocket_gets_its_walls_back_even_off_its_point |
| 4 | minor | yes: triangle ring -> ValueError; cause was the angles operation (blunting) | yes | refuse removal creating doubled edge; blunt_acute_corners refuses a cut leaving < 3 points; test_rim_repair_keeps_a_triangle_ring_whole |
| 5 | minor | yes: counters ignored retreat; state left mutated | yes | change detected by array comparison; restore in finally; test_a_retry_that_moves_the_rim_and_does_no_better_is_put_back |
| 6 | minor | yes: empty land returned ~1 km2 | yes | empty land/footprint return empty; raster padded by the largest disc; test_unresolvable_water_needs_a_coast_and_sane_arguments |
| 7 | minor | yes: levels_ratio=1 loops | yes | argument validation; same test |
| 8 | minor | yes: seed_search iterated global seeds | yes | iterate seed_list; test_seed_search_tries_the_seeds_it_is_given |
| 9 | minor | yes: _version imported oceanmesh | yes | importlib.metadata; test_collect_imports_no_library |
| 10 | minor | yes: ladder NetCDFs, region files, LR_* not recorded | yes | recorded (files + environment) |
| 11 | minor | yes: --untracked-files=no | yes | untracked counted, path_tracked; test_git_state_counts_untracked_files... |
| 12 | minor | yes: retry geometry depends on first search | yes | search_pass recorded; USER_GUIDE §5 says replay provenance.seeds |
| 13 | minor | yes: fixtures missed their branch; merge never taken in 2,500 random caps | yes | fixtures rebuilt with branch asserts; merge_ends removed as dead code (proper fix) |
| 14 | nit | yes | yes | USER_GUIDE: one third of the finest target |


## Round 2

### Answer (verbatim)

Read-only review completed. No files changed, batch jobs submitted, or mesh generation performed.

The standard pytest command failed before collection because no writable temporary directory was available. A scoped run with `--capture=sys --noconftest -p no:cacheprovider` produced **203 passed, 34 failed, 50 errors**, with temporary-directory and Matplotlib-cache failures preventing suite verification. A focused selection passed **26 tests**; three driver regression tests also passed with writes mocked. The reproductions below ran in memory using the requested interpreter.

**Previous findings**

| Round 1 finding | Status | Reason |
|---|---|---|
| 1. Repair under preserve | RESOLVED | `RIM_REPAIR` now requires `hires.coastline == "resolve"`; regression test passes. |
| 2. Lagoon entrance closed | RESOLVED | The original lagoon fixture stays connected, and the entrance rejection counter is exercised. |
| 3. Refused pocket loses walls | RESOLVED | Restoration checks pocket geometry; the off-centre refusal regression passes. |
| 4. Triangle ring collapse | RESOLVED | Duplicate-edge removal and undersized blunting are refused; regression passes. |
| 5. Retreat-only retry | RESOLVED | Array comparison detects movement; `finally` restores unselected state. The mocked retry regression passes. |
| 6. Empty land becomes land | RESOLVED | Empty land returns an empty result. Finding 2 below is a separate regression for **all-land** input. |
| 7. Invalid continuous-width arguments | RESOLVED | The reported invalid spacing, radius and level-ratio cases are rejected; empty footprints are handled. |
| 8. Ignored seed argument | RESOLVED | `seed_search` iterates `seed_list`; the differing-argument regression passes. |
| 9. GPL import through provenance | RESOLVED | Versions use distribution metadata; the no-import test passes. |
| 10. Missing provenance inputs/settings | PARTIAL | Bathymetry and environment settings are recorded, but external geometry datasets remain incompletely hashed; see finding 6. |
| 11. Untracked source omitted | RESOLVED | Git status includes untracked files, and `path_tracked` is recorded. |
| 12. Accepted-seed-only replay | RESOLVED | The guide requires the original seed sequence for pass 2, and the selected pass is recorded. |
| 13. Ineffective fixtures | RESOLVED | Updated fixtures reach their guards. In-memory guard removal makes the coast-movement and pocket-width assertions fail. The merge operation was removed. |
| 14. Fixed 10 m raster documentation | RESOLVED | Documentation now states one-third of the finest target. |

**Findings**

1. **Major — Substantial resolved land can still be accepted as water.**  
   **Location:** [notebooks/420_local_refine.py:703](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:703)

   A polygon refused by `island_rings` for crossing the rim produces only a warning, regardless of how much resolvable land remains inside the meshed domain. Neither seed acceptance nor the final `ACCEPTED` condition gates this report.

   **Reproduction:** Use water `box(0,0,1000,1000)`, land `box(700,200,1200,800)`, and constant `h=30`. After `filter_shoreline_local` with continuous width enabled, `island_rings` returns **zero rings** and reports **180,000 m²** of land inside the water domain. That land survives filtering but is omitted from the constraints. Geometric mesh QA does not establish that it was represented correctly.

   **Fix:** Reconcile substantial crossing land with the rim, or fail before meshing when reconciliation would violate the frozen interface. Gate unresolved source-land overlap explicitly, distinguishing deliberately refused artificial wall pockets. This is residual behavior, not a regression introduced by round 1.

2. **Minor — The round 1 raster fix crashes when the padded raster contains only land.**  
   **Location:** [src/fvcom_mesh_tools/patch.py:823](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:823)

   The new unconditional reductions over `r[sea]` assume at least one water cell.

   **Reproduction:**
   ```python
   unresolvable_water(
       shapely.box(-1000, -1000, 1000, 1000),
       lambda q: np.full(len(q), 30.0),
       shapely.box(0, 0, 100, 100),
   )
   ```
   HEAD raises `ValueError: zero-size array to reduction operation minimum which has no identity`. Executing the pre-`7a6cffe` implementation with identical arguments returns an empty result.

   **Fix:** Return the empty report before these reductions when `water.any()` is false. Add an all-land regression alongside the empty-land test.

3. **Minor — Continuous width can invent an isolated island in open water.**  
   **Location:** [src/fvcom_mesh_tools/patch.py:860](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:860)

   The rule rejects pieces touching multiple land components but accepts pieces touching **zero** components. Such a piece is neither a coastal strip nor a dead end.

   **Reproduction:** Set footprint `box(-1000,-1000,1000,1000)` and land `box(-2000,-2000,2000,0)`. Use `h=1000` inside `abs(x)<200, 300<y<700`, and `h=30` elsewhere. With `radius_factor=.75`, `min_h=60`, and `spacing=10`, the function returns **129,583.45 m²** of artificial land containing `(0,500)`, separated from existing land by **320.01 m**.

   This uses an unusual discontinuous size field, but it is accepted by the API’s validation.

   **Fix:** Require exactly one adjacent land component; reject or report unattached pieces. Retain the connectivity check for attached pieces.

4. **Minor — Wall-pocket closure can erase an unsampled fine-resolution region.**  
   **Location:** [src/fvcom_mesh_tools/walls.py:685](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:685)

   “The finest element anywhere” is evaluated only at polygon vertices and one representative point. A finer region between those samples is missed.

   **Reproduction:** Use land `box(-2000,-2000,2000,0)` and wall:
   ```python
   [(0,0), (0,380), (400,380), (400,150), (40,150)]
   ```
   Set `h=30` within 40 m of `(100,265)` and `h=200` elsewhere. Calling `close_wall_pockets(..., min_h=60.000001)` closes **92,000 m²**, including the fine region, and removes the enclosing wall sections.

   **Fix:** Evaluate the size field over the pocket interior and edges at a controlled spacing, or use region geometry to establish a conservative minimum. Refuse closure where the fine-zone exclusion cannot be established.

5. **Minor — Rim repair leaves obsolete curves controlling subsequent boundary movement.**  
   **Location:** [notebooks/420_local_refine.py:1209](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1209), consumed at [notebooks/420_local_refine.py:1684](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1684)

   `apply_rim_repair` updates rim coordinates, edges and the hole, but leaves `rc["curves"]` unchanged. `improve_patch` consequently treats some newly created corners as sliding points on the old boundary.

   **Reproduction:** The existing pier-retreat fixture moves its left tip from `(190,297)` to `(190,270)`. For the valid local triangle fan:
   ```python
   nodes = [[190,270], [190,240], [210,270], [160,275]]
   elements = [[1,0,3], [0,2,3]]
   ```
   allowing only node 0 to slide makes one improvement round move that corner to **`(190,258)`** using the old curve. Using the repaired curve keeps it at **`(190,270)`**, because it is correctly recognized as a corner.

   **Fix:** Maintain separate source-reference and current sliding curves. Update sliding curves after geometry-changing repair, and include them in retry snapshot/restoration.

6. **Minor — Geometry provenance still omits consumed files.**  
   **Location:** [notebooks/420_local_refine.py:340](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:340)

   External region files are hashed as one pathname. For shapefiles, `_geometry_from_file` also consumes `.prj` and `.dbf`: changing the CRS changes placement, and changing an attribute can change which polygon `where:` selects. Neither change affects the recorded `.shp` hash.

   Conversely, `FMESH_LAND` is read through unrestricted `gpd.read_file`, but provenance always substitutes `.shp/.shx/.dbf/.prj` paths. A valid GeoJSON land input therefore hashes nonexistent substitute files instead of the consumed GeoJSON.

   **Fix:** Centralize dataset-file enumeration: hash the actual file for single-file formats and the consumed sidecars for shapefiles. Use it for both land and regions. This completes round 1 finding 10.

7. **Minor — Invalid falsy `hires.experimental` values bypass validation.**  
   **Location:** [src/fvcom_mesh_tools/refine.py:338](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/refine.py:338)

   `spec.get("experimental", []) or []` converts invalid values into a valid empty list before type checking.

   **Reproduction:** Each of `_hires({"experimental": value}, {})` for `False`, `0`, `""`, and `{}` succeeds and returns `experimental=[]`. This contradicts the subsequent list-only validation and silently accepts malformed recipes.

   **Fix:** Validate the supplied value before normalization. Handle `None` explicitly only if null is intentionally supported.

8. **Nit — The driver documents an ignored seed variable.**  
   **Location:** [notebooks/420_local_refine.py:17](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:17)

   The environment documentation advertises `LR_SEED`, default `0`. The search reads only `LR_SEEDS`, default `0,1,2,3,4`; setting the documented variable does not select the seed.

   **Fix:** Document `LR_SEEDS` and its comma-separated default. This is pre-existing.

## Verdict
VERDICT: FAIL (0 blocker, 1 major, 6 minor, 1 nit)


### Prompt

```markdown
# Review request, round 2: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 20 commits, the last, 7a6cffe, fixing round 1; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Round 1 (FAIL: 0 blocker, 3 major, 10 minor, 1 nit). All 14 were reproduced
and fixed in commit 7a6cffe (`git show 7a6cffe`); the triage table is in
`docs/coastline-rules-review-20260927.md`.
1. rim repair under preserve -> `RIM_REPAIR` requires coastline resolve.
2. lagoon entrance closed -> `unresolvable_water` leaves a piece whose
   closing splits the water (`_splits_water`, piece grown by 2 cells).
3. refused pocket lost walls -> restore judged on the pocket polygon (wkt).
4. triangle ring collapse -> removal refusing doubled edges; blunting
   refuses a cut leaving < 3 ring points.
5. retreat-only retry -> change by array comparison; restore in `finally`.
6./7. empty land/footprint, validation, raster padded by the largest disc.
8. `seed_search` uses `seed_list`.
9. provenance via `importlib.metadata` (no import).
10./11. provenance records ladder files, region files, LR_*/FMESH_*/DATA_DIR;
   untracked files, `path_tracked`.
12. `search_pass` recorded; USER_GUIDE §5 says replay `provenance.seeds`.
13. fixtures rebuilt; the merge-to-midpoint operation was removed as dead
   code (a random search of 2,500 caps never took it: whatever refused the
   removal refused the merge) rather than given a test that cannot reach it.
14. doc: raster of one third of the finest target.
After the fixes all six hires recipes were rebuilt on compute nodes and are
accepted (transition coasts move up to 122 m; Yokohama now through the
retry, search_pass 2). 885 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes: crossing land only warned; checked all six recipes -- their land in the hole is chord slivers, no piece holds a half-element disc | yes | run fails before meshing if land inside the hole holds a disc of 0.5 h (refused wall pockets excepted); land_meshed_as_water_m2 reported; test_land_wide_enough_for_an_element_may_not_be_meshed_as_water |
| 2 | minor | yes: all-land raster raised | yes | return empty when no water; test_unresolvable_water_on_all_land_and_on_detached_water |
| 3 | minor | yes: 129,583 m2 island from a jumping size field | yes | a piece touching no land is left open (n_detached_left_open); same test |
| 4 | minor | yes: fine region between vertices missed | yes | size sampled on an interior grid of a quarter element; test_close_wall_pockets_sees_a_fine_region_inside_the_pocket |
| 5 | minor | yes (reasoned from the slide logic and the reviewer's fan) | yes | points rim_repair makes or moves are pinned in the seam repair (RIM_PINNED, restored with the retry snapshot); test_rim_repair_pins_the_points_it_makes_or_moves |
| 6 | minor | yes: .shp only / substituted sidecars for GeoJSON | yes | provenance.dataset_files for land and region files; test_dataset_files_... |
| 7 | minor | yes: `or []` | yes | None -> [], anything else must be a list; test_hires_experimental_must_be_a_list |
| 8 | nit | yes | yes | driver header documents LR_SEEDS |
| (self) | minor | found while rebuilding the recipes: Yokohama logged 2,180 tip retreats, most by 0 m | yes | a tip steps back to just past the throat test (1.001 x gap_factor x h) and a step that gets nowhere is refused; test_rim_repair_steps_each_tip_back_once[1.0, 1.5] |

After the round-2 fixes all six hires recipes rebuild accepted: five
byte-identical to the round-1 rebuild, Yokohama (search_pass 2, 6 tip
retreats) moving by up to 22 m.


## Round 3

### Answer (verbatim)

Read-only review completed at `2d8d1d7`. No files changed, mesh generation performed, or batch jobs submitted.

Using the requested interpreter, the standard pytest invocation failed before collection because the sandbox has no writable temporary directory. With `--capture=sys --noconftest -p no:cacheprovider`, the expanded geometry selection passed **50 tests**. Three driver regressions also passed with filesystem writes mocked. The findings below use in-memory reproductions. `git diff --check` reported Markdown hard-break whitespace in the archived review; the working tree remains clean.

**Previous findings**

| Finding | Status | Reason |
|---|---|---|
| R1 #1: Repair under preserve | RESOLVED | `RIM_REPAIR` requires `resolve`; regression passes. |
| R1 #2: Lagoon entrance closed | RESOLVED | Original lagoon fixture remains connected; regression passes. |
| R1 #3: Refused pocket loses walls | RESOLVED | Restoration uses polygon overlap; off-centre regression passes. |
| R1 #4: Triangle ring collapse | RESOLVED | Undersized rings and duplicate-edge removals are refused. |
| R1 #5: Retreat-only retry | RESOLVED | Array comparison detects changes; mocked search/restoration regression passes. |
| R1 #6: Empty land becomes land | RESOLVED | Empty-land regression passes. |
| R1 #7: Invalid continuous-width arguments | RESOLVED | Reported invalid arguments are rejected; empty footprints are handled. |
| R1 #8: Ignored seed argument | RESOLVED | Mocked regression confirms the supplied seed list is used. |
| R1 #9: GPL import through provenance | RESOLVED | Distribution metadata replaces imports; regression passes. |
| R1 #10: Missing provenance inputs/settings | PARTIAL | Added inputs/settings are recorded, but dataset inventory still misses uppercase sidecars; finding 5. |
| R1 #11: Untracked source omitted | RESOLVED | Untracked files and `path_tracked` are recorded. |
| R1 #12: Accepted-seed-only replay | RESOLVED | Documentation requires the original sequence for pass 2; selected pass is recorded. |
| R1 #13: Ineffective fixtures | RESOLVED | Revised fixtures exercise their operations/guards; unreachable merge operation was removed. |
| R1 #14: Fixed 10 m documentation | RESOLVED | Documentation states one-third of the finest target. |
| R2 #1: Resolved land accepted as water | PARTIAL | Constant-size regression passes, but variable sizing and subsequent repairs bypass protection; findings 1–2. |
| R2 #2: All-land raster crashes | RESOLVED | All-land regression passes. |
| R2 #3: Detached artificial island | RESOLVED | Detached pieces are rejected; regression passes. |
| R2 #4: Unsampled fine region erased | PARTIAL | Original fixture passes, but the coarse sampling grid still misses smaller fine regions; finding 4. |
| R2 #5: Obsolete sliding curves | PARTIAL | Moved/created points are pinned, but unchanged points that become corners are missed; finding 3. |
| R2 #6: Geometry provenance omissions | PARTIAL | Lowercase sidecars and GeoJSON are handled; uppercase sidecars remain omitted; finding 5. |
| R2 #7: Falsy experimental values | RESOLVED | Direct calls reject `False`, `0`, `""`, and `{}`; `None` becomes `[]`. |
| R2 #8: Wrong seed-variable documentation | RESOLVED | Header documents `LR_SEEDS`. |
| Author’s zero-progress retreat finding | RESOLVED | Both threshold regressions pass; the longer retreat exposes the separate validation gap in finding 2. |

**Findings**

1. **Major — The land-overlap guard misses resolvable land when size varies across a polygon.**  
   Location: [420_local_refine.py:732](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:732).

   The guard samples `h_achieved` at one representative point and uses that radius everywhere in the polygon.

   **Reproduction:** Set the hole to `box(0,0,3000,1000)`, land to `box(-100,400,2500,500)`, and `h(x,y)=clip(30+0.2*x,30,400)`. The actual shoreline filter retains this land with `h0=30`, `land_width_factor=.5`, `land_width_max_band=1`, `keep_land=land`, and continuous width enabled. `island_rings` refuses it for crossing the rim.

   Executing the actual driver block then permits **250,000 m²** of land inside the hole: its representative point samples `h=280`, so erosion by 140 m is empty. Yet at `(100,450)`, `h=50`, and a **25 m radius disc fits entirely inside that land**. The field’s slope is only 0.2; this does not require a discontinuity.

   **Fix:** Test clearance against the size at potential disc centres throughout the overlap, using controlled/adaptive sampling or conservative bounds. A single representative-point size cannot implement the stated local-size criterion.

2. **Major — Later rim repairs invalidate the land check, and the report retains the pre-repair overlap.**  
   Locations: [420_local_refine.py:735](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:735), [420_local_refine.py:2105](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:2105).

   The overlap check runs before blunting, wall rooting, initial rim repair, and the feedback retry. None of those paths repeats it. The round-2 change making retreat respect `gap_factor=1.5` can now introduce land that fails the new guard’s own criterion.

   **Reproduction:** Use this cyclic, all-free water rim with constant `h=30`:

   ```python
   [(0,0), (190,0), (190,297), (310,297),
    (310,0), (600,0), (600,300), (0,300)]
   ```

   The pier land is `box(190,0,310,297)`. Execute the actual island/check block, default `apply_rim_repair()`, then the retry with `min_edge_factor=.75`, `gap_factor=1.5`, `rounds=1`, focused at `(190,269.97)` and `(310,269.97)`.

   The report still says **`land_meshed_as_water_m2=0`**, while the resulting hole overlaps **5,045.4 m²** of pier land. Eroding that overlap by `0.5*h=15` leaves **1,084.05 m²**, so it unequivocally fails the intended guard. Seed acceptance does not recheck it.

   **Fix:** Recompute and gate overlap after geometry changes, including retry candidates and the delivered boundary. Reject or restore repairs that violate the land contract, and report the selected geometry’s overlap.

3. **Minor — Pinning misses unchanged points that become corners when a neighbour is removed.**  
   Location: [420_local_refine.py:1242](/octfs/work/G16445/v610con-mesh-tools/notebooks/420_local_refine.py:1242).

   `RIM_PINNED` detects new coordinates, but a topology change can create a corner without moving its surviving vertex. That vertex remains eligible to slide along the obsolete source curve.

   **Reproduction:** With constant `h=30`, cyclic edges and all-free points, use:

   ```python
   [(0,0), (190,0), (190,277), (190,297), (202,297),
    (202,277), (202,0), (600,0), (600,500), (0,500)]
   ```

   Actual `apply_rim_repair()` removes `(190,297)` and returns **`RIM_PINNED=[]`**. Its surviving neighbour `(190,277)` is now a corner.

   For the valid fan with nodes `[(190,277),(190,250),(202,297),(160,280)]` and faces `[[1,0,3],[0,2,3]]`, allowing only node 0 to slide along the original pier curve moves it to **`(190,273.76)`** in one improvement round. Using the repaired curve keeps it at `(190,277)`.

   **Fix:** Maintain current sliding curves separately from source-reference curves, or also pin surviving vertices whose incident geometry changes. Include that state in retry restoration.

4. **Minor — Pocket sampling still misses fine regions smaller than the coarse sampling grid.**  
   Location: [walls.py:691](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:691).

   Grid spacing is chosen from the initial coarse samples. It therefore cannot establish that no finer region lies between those samples.

   **Reproduction:** Use land `box(-2000,-2000,2000,0)` and wall:

   ```python
   [(0,0), (0,380), (400,380), (400,150), (40,150)]
   ```

   Set `h=30` within 20 m of `(100,250)` and `h=200` elsewhere. With `min_h=60.000001`, `close_wall_pockets` closes **92,000 m²**, including the fine region. The 50 m sampling grid misses that region completely.

   **Fix:** Use declared fine-region geometry or conservative size bounds to establish the exclusion. For a sampled implementation, require a known minimum feature scale/variation bound and sample accordingly; quarter of the initially observed coarse size is insufficient.

5. **Minor — Uppercase shapefile sidecars remain absent from provenance.**  
   Location: [provenance.py:50](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:50).

   Extension recognition is case-insensitive, but sidecar lookup constructs only lowercase suffixes.

   **Reproduction:** For `LAND.SHP`, `LAND.SHX`, `LAND.DBF`, `LAND.PRJ`, and `LAND.CPG`, mocked filesystem inventory makes `dataset_files("…/LAND.SHP")` return **only `LAND.SHP`**. I also verified that GDAL reads an uppercase dataset through an in-memory ZIP, including its CRS and attributes. Changing the consumed uppercase `.PRJ` or `.DBF` therefore leaves the recorded dataset hash unchanged.

   **Fix:** Enumerate actual sibling filenames with case-insensitive sidecar matching, preserving their actual paths, and test uppercase/mixed-case datasets.

6. **Minor — Package provenance records the driver’s checkout rather than the imported package’s checkout.**  
   Location: [420_local_refine.py:336](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:336).

   The `"fvcom_mesh_tools"` entry uses the notebook’s `__file__`. Python can import the package from a different editable checkout or installed distribution; the driver neither records that origin nor checks that it matches.

   **Evidence:** Executing the extracted provenance assignment with a different driver-tree path while retaining the actual imported package attributes that other tree to `"fvcom_mesh_tools"`. The imported implementation remains in this repository. Running checkout A’s notebook with an environment installed from checkout B has this same mismatch.

   **Fix:** Record the driver and imported package separately, deriving the package entry from its actual module `__file__`. For an installation outside Git, record distribution identity and relevant source hashes.

7. **Minor — Geometry iterators are consumed twice, silently dropping land or walls.**  
   Locations: [patch.py:954](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:954), [walls.py:713](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:713).

   Both helpers accept iterable inputs but reuse them after consumption.

   **Reproductions:**

   - For land `[box(-500,-500,1500,500)]`, footprint `box(0,-1000,1000,1000)`, constant `h=120`, `h0=30`, and `keep_land` equal to that land polygon, `filter_shoreline_local` returns **2,000,000 m²** with a list but **1,000,000 m²** with `iter(list)`. The iterator case also incorrectly reports zero original land area.
   - For a single 500 m pier `[(0,0),(0,500)]`, coastal land below `y=0`, and `h=30`, `close_wall_pockets` reports zero closures in both cases. A list preserves **500 m** of wall; an iterator returns **zero walls**.

   **Fix:** Materialize iterable inputs once on entry and reuse the materialized geometry for every pass, report, and unchanged-result return.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 5 minor, 0 nit)


### Prompt

```markdown
# Review request, round 3: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 21 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Round 1 (FAIL 0/3/10/1): all fixed in 7a6cffe; your round-2 status table
had 13 RESOLVED and #10 PARTIAL.
Round 2 (FAIL 0/1/6/1): all fixed in 2d8d1d7 (`git show 2d8d1d7`; triage in
`docs/coastline-rules-review-20260927.md`).
1. crossing land only warned -> the driver fails before meshing if land
   inside the hole holds a disc of 0.5 h (refused wall pockets excepted);
   `land_meshed_as_water_m2` reported. (All six recipes' land in the hole
   is chord slivers, none holds such a disc; all still accepted.)
2. all-land raster -> empty result.
3. detached water -> left open (`n_detached_left_open`).
4. pocket size sampled on an interior grid of a quarter element.
5. points rim_repair makes or moves are pinned in the seam repair
   (`RIM_PINNED`, restored with the retry snapshot).
6. `provenance.dataset_files` (shapefile sidecars; GeoJSON itself) for land
   and region files -- completes round-1 #10.
7. `hires.experimental` must be a list (None -> []).
8. driver header documents `LR_SEEDS`.
Also fixed (found by us while rebuilding): `rim_repair`'s tip retreat put
the tip exactly at the throat distance, rounding left it inside, and it was
stepped back again by 0 m -- 2,180 times at Yokohama, and never out under
the retry's gap_factor 1.5. It now steps to 1.001 x gap_factor x h and a
no-progress step is refused.
All six hires recipes rebuild accepted on compute nodes. 896 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes: reviewer's sloping field; one point's h hid a 25 m disc | yes | patch.land_an_element_fits: disc centres searched inside the land (0.5 h_min in) on a grid of 0.25 h_min, each judged at its own size; test_land_an_element_fits_uses_the_size_where_the_disc_stands |
| 2 | major | yes: retreat at gap 1.5 handed 5,045 m2 of a 120 m pier to the water, report stale | yes | retreat refuses to hand over land a half-element disc fits in; the land check runs again after rim repair (fails the run) and on the retry (retry not taken); report updated; test_a_pier_tip_steps_back_only_over_land_no_element_fits_in, test_the_land_check_runs_again_after_a_rim_repair |
| 3 | minor | yes: survivor of a removal not pinned | yes | a point whose rim neighbours changed is pinned too; test_rim_repair_pins_a_survivor_that_becomes_a_corner |
| 4 | minor | yes: 20 m fine region missed at 50 m grid | yes | pocket grid spacing from min_h (the gate), not from the coarse size; test_close_wall_pockets_sees_a_fine_region_smaller_than_the_coarse_grid |
| 5 | minor | yes: uppercase sidecars | yes | case-insensitive sibling match; test_dataset_files_matches_sidecars_without_case |
| 6 | minor | yes: driver path recorded as the package | yes | package entry from sys.modules["fvcom_mesh_tools"].__file__, driver recorded separately |
| 7 | minor | yes: iterators spent | yes | materialised on entry; test_filter_shoreline_local_takes_an_iterator_of_land, test_close_wall_pockets_takes_an_iterator_of_walls |
| (self) | major | found while rebuilding: with the round-3 checks Yokohama failed -- its round-2 retry had stepped pier tips back over land an element fits in (375,761 3,925,679), which the new check refused | yes | when the retry fails the land check it is redone with `retreat_tips=False` (slits and short edges kept); Yokohama accepted that way; test_a_retry_that_wets_land_is_redone_without_stepping_tips_back |

After the round-3 fixes all six hires recipes rebuild accepted: five
byte-identical to the round-2 rebuild except Funabashi (5 m) and Futtsu
nori (22 m) from the wider pinning, Yokohama through the kept-tips retry.


## Round 4

### Answer (verbatim)

Read-only review completed at `fd2bbe1`; the working tree remains clean. No mesh generation or batch jobs ran.

The requested full pytest command failed before collection because no writable temporary directory exists. Focused runs using `--capture=sys --noconftest -p no:cacheprovider` reached **50 passed**, with one additional test blocked by Matplotlib’s writable-cache requirement. Three seed/retry regressions passed with filesystem writes mocked. Findings below use source inspection and small in-memory reproductions; production mesh QA was not rerun.

**Previous findings**

| Finding | Status | Reason |
|---|---|---|
| R1 #1: Repair under preserve | RESOLVED | `RIM_REPAIR` requires `resolve`. |
| R1 #2: Continuous-width lagoon closure | RESOLVED | The original entrance regression passes; connectivity is checked. |
| R1 #3: Refused pocket loses walls | RESOLVED | Restoration uses polygon overlap; regression passes. |
| R1 #4: Triangle ring collapse | RESOLVED | Duplicate edges and undersized rings are guarded. |
| R1 #5: Retreat-only retry | RESOLVED | Array comparison detects movement; mocked restoration regression passes. |
| R1 #6: Empty land becomes land | RESOLVED | Empty-land regression passes. |
| R1 #7: Invalid continuous-width arguments | RESOLVED | Reported invalid arguments are rejected; empty footprints are handled. |
| R1 #8: Ignored seed argument | RESOLVED | Mocked regression confirms the supplied seeds are used. |
| R1 #9: GPL import through provenance | RESOLVED | Versions use distribution metadata without importing libraries. |
| R1 #10: Missing provenance inputs/settings | RESOLVED | Bathymetry, geometry datasets and effective environment settings are recorded. |
| R1 #11: Untracked source omitted | RESOLVED | Untracked files and `path_tracked` are recorded. |
| R1 #12: Accepted-seed-only replay | RESOLVED | Documentation requires the original sequence for pass 2. |
| R1 #13: Ineffective fixtures | RESOLVED | The originally reported fixtures were corrected; finding 8 concerns new tests. |
| R1 #14: Fixed 10 m documentation | RESOLVED | Documentation specifies one-third of the finest target. |
| R2 #1: Resolvable land accepted as water | PARTIAL | The guard exists, but findings 1–3 leave acceptance gaps. |
| R2 #2: All-land raster crashes | RESOLVED | All-land regression passes. |
| R2 #3: Detached artificial island | RESOLVED | Detached pieces are rejected; regression passes. |
| R2 #4: Unsampled fine region erased | PARTIAL | Denser sampling fixes the original fixture, but finding 5 remains. |
| R2 #5: Obsolete sliding curves | RESOLVED | Repair-created, moved and adjacency-changed rim points are pinned. |
| R2 #6: Geometry provenance omissions | RESOLVED | Dataset enumeration handles single-file formats and uppercase sidecars. |
| R2 #7: Falsy experimental values | RESOLVED | Direct calls reject `False`, `0`, `""` and `{}`; `None` becomes `[]`. |
| R2 #8: Wrong seed-variable documentation | RESOLVED | Header documents `LR_SEEDS`. |
| R3 #1: Variable-size land guard | PARTIAL | Local sampling fixes the reported example but is not conservative; findings 1–2. |
| R3 #2: Repairs bypass land check | PARTIAL | Default repair/retry paths recheck, but disabling rim repair bypasses later validation; finding 3. |
| R3 #3: Unchanged new corners unpinned | RESOLVED | Incident-neighbour changes now cause pinning. |
| R3 #4: Pocket sampling grid | PARTIAL | The original fixture passes; a smooth-field counterexample remains, finding 5. |
| R3 #5: Uppercase sidecars | RESOLVED | Case-insensitive extension enumeration includes all five uppercase files. |
| R3 #6: Imported-package provenance | PARTIAL | Separate checkout origins are fixed; non-Git installations still lose identity, finding 6. |
| R3 #7: Exhausted iterators | RESOLVED | Both helpers materialize inputs; regressions pass. |
| Author: Zero-progress retreats | RESOLVED | Both threshold regressions pass. |
| Author: Retry without tip retreat | RESOLVED | Mocked regression confirms restoration and a second search with tips kept. |

**Findings**

1. **Major — The new land guard can miss a qualifying disc between samples.**

   Location: [patch.py:1797](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1797).

   An empty sample result is treated as proof that no qualifying centre exists. The representative-point fallback does not establish that either.

   **Reproduction:** With constant local size `60`, `h_min=30`, and this equilateral triangle:

   ```python
   r = 31.0
   land = Polygon([
       (0, 0), (2*np.sqrt(3)*r, 0), (np.sqrt(3)*r, 3*r)
   ])
   ```

   `land_an_element_fits` returns `[]`, although the in-centre has **31 m clearance**, exceeding the required **30 m**.

   I also extended this triangle below a hole’s boundary and ran the actual shoreline filter with continuous width enabled, followed by the actual driver island/check block. The filter retained it, `island_rings` refused it for crossing the rim, and the driver permitted **4,993.50 m²** of land inside the water domain.

   **Fix:** Use a conservative adaptive search with clearance/size bounds and an explicit tolerance. An unresolved cell must not certify absence. For constant sizes, polygon erosion provides a direct check.

2. **Major — The driver supplies a lower bound that the size field does not guarantee.**

   Locations: [420_local_refine.py:736](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:736), [420_local_refine.py:1299](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1299).

   `FINE_H / 2` comes from the smallest requested target. However, `patch_sizing` deliberately retains finer existing base resolution. Consequently, the guard can erase the entire candidate-centre domain before evaluating the actual field.

   **Reproduction:** For a 30 m requested target and a 20 m base field, the actual `patch_sizing` closure returns 20 m. With land `box(0,0,25,1000)`, the actual `wet_land_after_repair` returns:

   ```text
   ([], 25000.0)
   ```

   A 10 m radius disc fits comfortably. Calling the helper with the correct lower bound of 20 m detects it. The reproduction mocked only the base-field provider.

   **Fix:** Derive a guaranteed lower bound from both the target sizes and the base-field construction. Use that bound in every initial and retry check.

3. **Major — Switching off rim repair also switches off validation of earlier coastline changes.**

   Location: [420_local_refine.py:1302](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1302).

   The post-change land check is inside `if RIM_REPAIR`, but the earlier `blunt_acute_corners` call still changes resolved coastlines when `hires.rim_repair: false` or `LR_EXPERIMENTAL=none`. Its frozen-corner operation explicitly gives land to water.

   **Reproduction:** Scale the existing frozen-corner fixture by three:

   ```python
   p = 3*np.array([
       (0,0), (-31.5,-20.4), (-72.7,34.6), (-400,34.6),
       (-400,-400), (0,-400), (0,-64.9)
   ])
   base = [10, -1, -1, 11, 12, 13, 14]
   ```

   With cyclic edges and `h=30`, executing the actual blunting block adds **11,578.41 m²** of land to the water. The land helper finds qualifying centres, but with rim repair disabled no subsequent check runs; `land_meshed_as_water_m2` remains **0**.

   **Fix:** Run the land check unconditionally after coastline construction whenever `WET_LAND_SRC` exists. Rule switches should control repairs, not validation.

4. **Major — Slit repair can erase a basin that accommodates its local elements.**

   Location: [patch.py:1688](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1688).

   The basin test takes the minimum size only at vertices on its boundary path. A finer basin interior is therefore judged using coarser boundary sizes.

   **Reproduction:** Use cyclic, all-free points:

   ```python
   [(0,0), (500,0), (500,500), (260,500),
    (260,540), (290,540), (290,620), (210,620),
    (210,540), (240,540), (240,500), (0,500)]
   ```

   For `c=(250,580)`, define the continuous field:

   ```python
   h(q) = min(100, 30 + 0.4*max(norm(q-c)-5, 0))
   ```

   Default `rim_repair` closes one slit and removes **7,200 m²** of water, including `c`. Yet a **30 m radius disc** centred at `c` fits entirely in the original basin and `h(c)=30`. The same failure occurs with only the slit operation enabled.

   The land-overlap check cannot detect this: the error removes water rather than adding it.

   **Fix:** Search potential disc centres throughout the basin using their local sizes and conservative bounds. Preserve the basin whenever a qualifying disc exists or absence cannot be established.

5. **Minor — The finer pocket grid still does not enforce the fine-zone exclusion.**

   Location: [walls.py:693](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:693).

   Reducing the spacing fixes the previous example but leaves the same sampling assumption.

   **Reproduction:** Use coastal land `box(-2000,-2000,2000,0)` and:

   ```python
   wall = [(0,0), (0,250), (100,250), (100,150), (20,150)]
   h(q) = 59.9 + 0.1*norm(q-(30,180))
   min_h = 60.000001
   ```

   `close_wall_pockets` closes **10,000 m²**, including `(30,180)`, where `h=59.9` is below the exclusion threshold. This field is continuous with slope at most **0.1**; no discontinuous spike is needed.

   **Fix:** Use fine-zone geometry or conservative interpolation/variation bounds. Near-threshold samples require refinement or refusal, not unconditional closure.

6. **Minor — Non-Git package installations still have no recorded implementation identity.**

   Location: [provenance.py:125](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:125).

   The driver now passes the imported package’s correct path, but `collect` retains only `git_state(path)`. Outside Git that becomes `null`, discarding even the supplied path. `fvcom_mesh_tools` is also absent from `LIBRARIES`.

   **Evidence:** Simulating the documented non-Git return from `git_state` produces:

   ```json
   {"code": {"fvcom_mesh_tools": null}}
   ```

   Different installed implementations can therefore have identical package provenance.

   **Fix:** Always record the resolved module path and distribution identity; add source hashes when Git identity is unavailable.

7. **Minor — The new land guard treats invalid size values as successful absence checks.**

   Location: [patch.py:1803](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1803).

   **Reproduction:** For `box(0,0,100,100)` and `h_min=30`, size functions returning either `NaN` or positive infinity make `land_an_element_fits` return `[]`. The caller interprets this as permission to proceed.

   The new `factor` argument is also unvalidated: `factor=0` raises `ZeroDivisionError`, while `factor=-1` returns a qualifying point.

   **Fix:** Require finite positive `factor` and sampled sizes, validate the returned array shape, and raise a clear error instead of returning an empty result.

8. **Minor — New regression tests can pass without exercising their advertised protection.**

   Locations: [test_local_refine_driver.py:462](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/tests/test_local_refine_driver.py:462), [test_local_refine_driver.py:467](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/tests/test_local_refine_driver.py:467).

   **Evidence:** Replacing `rim_repair` with an identity stub makes `test_rim_repair_pins_a_survivor_that_becomes_a_corner` pass: its only assertion is conditional on removal having happened.

   Likewise, `test_the_land_check_runs_again_after_a_rim_repair` calls only the extracted checking function. It never executes the driver’s call site, so removing the post-repair invocation would not fail this test.

   **Fix:** Assert the prerequisite removal before checking pinning. Exercise the actual post-construction and retry blocks, including repair disabled and an introduced-land rejection.

9. **Nit — The guide describes a priority policy the implementation ignores.**

   Location: [USER_GUIDE.md:241](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/docs/USER_GUIDE.md:241).

   The table says `refine[].priority` determines which overlapping region wins. `patch_sizing` explicitly ignores priority and takes the finest requested contribution; its docstring and existing regression describe that policy.

   **Fix:** State that priority has no effect in local refinement and that the finest applicable size wins.

## Verdict

VERDICT: FAIL (0 blocker, 4 major, 4 minor, 1 nit)


### Prompt

```markdown
# Review request, round 4: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 22 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, the last commit round 3; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1 and 2: see your round-3 status table (all RESOLVED but R1#10,
R2#1, R2#4, R2#5, R2#6 PARTIAL).
Round 3 (FAIL 0/2/5/0): all fixed in fd2bbe1 (`git show fd2bbe1`; triage in
`docs/coastline-rules-review-20260927.md`).
1. `patch.land_an_element_fits`: disc centres searched inside the land
   (0.5 h_min in) on a grid of 0.25 h_min, each judged at its own size;
   used by the driver check.
2. retreat refuses to hand over land a half-element disc fits in; the land
   check runs again after rim repair (fails the run) and on the retry; a
   retry that fails it is redone with `retreat_tips=False` (we found that
   Yokohama's round-2 mesh had exactly this defect; it is now accepted
   through the kept-tips retry).
3. a rim point whose neighbours changed is pinned too.
4. pocket sampling grid at 0.25 * min_h.
5. shapefile sidecars matched case-insensitively.
6. provenance records the imported package (`sys.modules`) and the driver.
7. iterators materialised on entry.
All six hires recipes rebuild accepted. 905 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes: equilateral triangle, in-centre 31 m vs 30 m, returned [] | yes | absence is now proved: erosion by factor x a guaranteed lower bound of the size over the polygon (patch.size_lower_bound: grid readings minus slope x half-diagonal, floored); an unresolved polygon is reported; test_land_an_element_fits_proves_absence_and_does_not_sample_it |
| 2 | major | yes (reasoned): patch_sizing keeps a base finer than the target | yes | H_FLOOR = min(finest target, min base ambient size); used by every check; same test (20 m floor) |
| 3 | major | yes: the check sat inside `if RIM_REPAIR` | yes | check unconditional after the rim is final; test_the_driver_stops_when_the_final_rim_leaves_land_in_the_water[True/False] |
| 4 | major | yes: fine basin closed on vertex sizes | yes | slit closure judged on size_lower_bound over the dead end; test_rim_repair_keeps_a_fine_basin_behind_a_coarse_mouth |
| 5 | minor | yes: 59.9 m at the gate 60 m missed | yes | pocket gate and disc test on size_lower_bound; test_close_wall_pockets_refuses_near_the_fine_zone_threshold |
| 6 | minor | yes: git_state None discarded the path | yes | provenance.code_state: path, distribution version, git or a source hash; fvcom-mesh-tools in LIBRARIES; test_code_outside_git_is_identified_by_its_sources |
| 7 | minor | yes: NaN/inf sizes, factor 0/-1 | yes | validated (ValueError); test_land_an_element_fits_refuses_a_bad_size_or_factor |
| 8 | minor | yes: conditional assert; function-only test | yes | premise asserted; driver call site exercised with rim repair on and off |
| 9 | nit | yes | yes | USER_GUIDE: priority has no effect in local refinement |

The size bounds assume the sizing field changes by at most 1 m per metre
(`slope=1`); the fields here change by at most about 0.2 (smoothed base)
and the gradation (0.165) in the ramps. A field that jumps faster is
outside the contract of `size_lower_bound`, which says so.

After the round-4 fixes all six hires recipes rebuild accepted and byte-identical
to the round-3 rebuild.

## Round 5

### Answer (verbatim)

Read-only review completed at `f07c6d9`; the working tree remains clean. No mesh generation or batch jobs ran. Production mesh QA, counts and implied time steps were **not rerun**.

Verification used the requested interpreter:

- Full pytest invocation with `-p no:cacheprovider` failed before collection: no writable temporary directory.
- Focused pytest run with `--capture=sys --noconftest -p no:cacheprovider`: **43 passed**, three temporary-directory setup errors.
- Broader selection of tests requiring no filesystem writes: **191 passed**. Matplotlib’s cache-directory writability check was mocked; no cache was written.
- Five seed/retry and driver-land-check cases passed with filesystem writes mocked.
- Findings below were reproduced in memory.

**Previous findings**

| Finding | Status | Reason |
|---|---|---|
| R1 #1: Repair under preserve | RESOLVED | Repair and retry require `resolve`. |
| R1 #2: Lagoon entrance closure | RESOLVED | Connectivity protection remains; original regression passes. |
| R1 #3: Refused pocket loses walls | RESOLVED | Restoration uses pocket geometry rather than its representative point. |
| R1 #4: Triangle ring collapse | RESOLVED | Duplicate-edge and minimum-ring-size guards remain; regression passes. |
| R1 #5: Retreat-only retry | RESOLVED | Array comparison detects movement; mocked search/restoration regression passes. |
| R1 #6: Empty land becomes land | RESOLVED | Empty-land regression passes. |
| R1 #7: Invalid continuous-width arguments | RESOLVED | Reported parameter errors and empty footprints are handled. Finding 5 concerns the separate band filter. |
| R1 #8: Ignored seed argument | RESOLVED | Mocked regression confirms the supplied seeds are searched. |
| R1 #9: GPL import through provenance | RESOLVED | Distribution metadata replaces library imports. |
| R1 #10: Missing provenance inputs/settings | RESOLVED | Bathymetry, geometry datasets and effective environment settings are recorded. |
| R1 #11: Untracked source omitted | RESOLVED | Git status includes untracked files; tracking status is recorded. |
| R1 #12: Accepted-seed-only replay | RESOLVED | Pass-2 replay requires the recorded original seed sequence. |
| R1 #13: Ineffective fixtures | RESOLVED | Fixtures exercise their guards; the unreachable merge operation was removed. |
| R1 #14: Fixed 10 m documentation | RESOLVED | Guide specifies one-third of the finest target. |
| R2 #1: Resolvable land accepted as water | PARTIAL | Validation exists, but nested polygon collections bypass it; finding 4. |
| R2 #2: All-land raster crash | RESOLVED | Regression passes. |
| R2 #3: Detached artificial island | RESOLVED | Detached closure pieces are rejected; regression passes. |
| R2 #4: Fine region inside pocket erased | PARTIAL | Original fixture passes, but fine-zone exclusion still fails; finding 1. |
| R2 #5: Obsolete sliding curves | RESOLVED | Repair-created, moved and adjacency-changed rim points are pinned. |
| R2 #6: Geometry provenance omissions | RESOLVED | Dataset enumeration includes single-file formats and shapefile sidecars. |
| R2 #7: Falsy experimental values | RESOLVED | Schema explicitly rejects non-list values other than `None`. |
| R2 #8: Wrong seed-variable documentation | RESOLVED | Header documents `LR_SEEDS`. |
| R3 #1: Variable-size land guard | PARTIAL | Original example is fixed; the new bound’s guarantee remains incomplete, findings 1–2. |
| R3 #2: Repairs bypass land validation | RESOLVED | Final-rim validation runs with repair enabled or disabled; retry validation remains. |
| R3 #3: Unchanged new corners unpinned | RESOLVED | Changed adjacency causes pinning; strengthened regression passes. |
| R3 #4: Pocket sampling grid | PARTIAL | Denser sampling was replaced, but the replacement still misses a fine region; finding 1. |
| R3 #5: Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6: Imported-package provenance | RESOLVED | Actual imported paths and non-Git source identity are recorded. Finding 6 concerns failure handling introduced by that fallback. |
| R3 #7: Exhausted iterators | RESOLVED | Both helpers materialize inputs; regressions pass. |
| R4 #1: Land guard misses qualifying disc | PARTIAL | Triangle regression passes, but the claimed bound is not generally guaranteed; findings 1–2. |
| R4 #2: Incorrect size floor | RESOLVED | `H_FLOOR` includes the minimum base ambient size. |
| R4 #3: Repair switch disables validation | RESOLVED | Both driver call-site cases pass with writes mocked. |
| R4 #4: Fine basin behind coarse mouth erased | PARTIAL | Supplied basin regression passes; the driver does not satisfy the bound’s slope premise everywhere, finding 2. |
| R4 #5: Pocket fine-zone exclusion | PARTIAL | A continuous field within the stated slope contract still defeats it; finding 1. |
| R4 #6: Non-Git implementation identity | RESOLVED | Path, distribution and source hash are now retained. |
| R4 #7: Invalid sizes/factors accepted | RESOLVED | Reported land-helper cases now raise; regressions pass. |
| R4 #8: Tests omit their prerequisite/call site | RESOLVED | Removal is asserted and both driver branches are exercised. |
| R4 #9: Priority documentation | RESOLVED | Guide states that priority has no effect in local refinement. |
| Author: Zero-progress retreats | RESOLVED | Threshold regressions pass. |
| Author: Retry without tip retreat | RESOLVED | Mocked regression confirms restoration and a search with tips kept. |

**Findings**

1. **Minor — The “guaranteed” size bound still misses fine regions beside polygon edges.**

   Location: [patch.py:1809](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1809), used by [walls.py:692](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:692).

   The half-diagonal covering distance applies to the complete rectangular grid. Removing every grid point outside or on the polygon destroys that guarantee. Polygon vertices do not replace the missing samples along long edges.

   Reproduction:

   ```python
   wall = LineString([
       (0, 0), (0, 250), (100, 250), (100, 150), (20, 150)
   ])
   c = np.array([48.75, 150.1])
   size = lambda q: 59.9 + 0.8*np.linalg.norm(np.asarray(q)-c, axis=1)

   added, _, report = close_wall_pockets(
       [wall], box(-2000, -2000, 2000, 0), size,
       min_h=60.000001, size_floor=30.0,
   )
   ```

   This field has slope at most **0.8**, within the documented contract. Nevertheless:

   ```text
   reported lower bound: 61.2334453318
   actual size at c:     59.9
   pockets closed:       1
   area closed:          10,000 m²
   added.contains(c):    True
   ```

   A simpler direct counterexample is `box(0,0,100,1)`, `h_floor=20`, and `h(q)=30+||q-(50,0.5)||`: the returned bound is **76.466966**, exceeding the actual minimum **30**.

   **Fix:** Use samples with a proven covering distance, including boundary cells. For a globally defined Lipschitz field, retaining the complete bounding grid is one option. Otherwise establish coverage on intersecting cells or fall back to `h_floor`.

2. **Minor — The driver supplies a sizing field that violates the new bound’s slope contract.**

   Locations: [420_local_refine.py:629](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:629), [patch.py:2329](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2329).

   Outside the base mesh, `"nearest"` selects a nearest **node’s** ambient value. This is discontinuous across Voronoi boundaries whenever adjacent values differ. Twenty smoothing passes do not remove those discontinuities. Thus the driver cannot generally use `slope=1` to certify absence.

   Reproduction using the actual `base_size_field`, without mocking its arithmetic:

   ```python
   xy = np.array([
       (0,0), (100,0), (100,100), (0,100),
       (300,0), (300,100), (700,0), (700,100),
   ], dtype=float)
   tri = np.array([
       (0,1,3), (1,2,3), (1,4,2),
       (4,5,2), (4,6,5), (6,7,5),
   ])
   f = base_size_field(xy, tri, outside="nearest")
   f([[49.999999, -10], [50.000001, -10]])
   ```

   Result: **192.46580798** and **194.62454752** across **0.000002 m**, an apparent slope above **1,000,000**.

   A second reproduction using actual `patch_sizing` on a symmetric small grid returned a bound of **318.421719 m** where the field was **303.941149 m**. Even the complete bounding grid missed the smaller value, so fixing finding 1 alone does not fix this issue.

   **Fix:** Provide a proven Lipschitz extension, or compute bounds directly from the interpolation/Voronoi regions. Use the global floor where the slope premise cannot be established.

3. **Minor — The new lower bound weakens the wall-pocket minimum-area and coast-clearance guards.**

   Location: [walls.py:700](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:700).

   A lower bound is conservative for proving that no disc fits. It is not conservative for proving that an island is large enough or far enough from land.

   With constant `size=100`, `size_floor=30`, and `min_h=60.000001`, use:

   ```python
   wall = LineString([
       (0,0), (0,y+100), (width,y+100), (width,y), (20,y)
   ])
   land = box(-2000, -2000, 2000, 0)
   ```

   Actual results:

   | `width`, `y` | Accepted area | Coast clearance | Violation |
   |---|---:|---:|---|
   | `90`, `150` | 9,000 m² | 150 m | Below the 10,000 m² minimum |
   | `100`, `48` | 10,000 m² | 48 m | Below the 50 m clearance |
   | `90`, `48` | 9,000 m² | 48 m | Both |

   The 9,000 m² pocket also passes `island_rings`, producing one four-point island.

   **Fix:** Keep separate bounds for predicates with opposite directions. Use conservative upper estimates or direct checks for the minimum-area and clearance requirements; retain the lower bound for the fine-zone exclusion and disc-absence proof.

4. **Minor — Nested polygon collections bypass the land-in-water guard entirely.**

   Location: [patch.py:1835](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1835).

   The helper traverses only one level and skips a `MultiPolygon` inside a `GeometryCollection`. Such structures are ordinary outputs of `make_valid`.

   Reproduction:

   ```python
   raw = Polygon([
       (0,0), (100,0), (100,100), (0,100), (0,0),
       (-100,0), (-100,-100), (-200,-100),
       (-200,0), (-100,0), (0,0),
   ])
   land = shapely.make_valid(raw)
   land_an_element_fits(land, lambda q: np.full(len(q), 30.0), 30.0)
   ```

   `make_valid` returns a collection containing a `MultiPolygon` and a line. The two polygons total **20,000 m²**, each easily accommodating the required 15 m disc, but the guard returns **`[]`**.

   **Fix:** Recursively traverse polygonal members of collections. Ignore only non-area members, not nested polygon containers.

5. **Minor — The local shoreline filter silently accepts zero and negative sizing fields.**

   Location: [patch.py:951](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:951).

   `np.maximum(h, h0)` converts invalid non-positive sizes into valid finest-band sizes before validation. With continuous width disabled, nothing subsequently rejects them.

   Reproduction:

   ```python
   filter_shoreline_local(
       box(0,0,100,100),
       lambda q: np.full(len(q), -10.0),
       30.0,
       box(-50,-50,150,150),
   )
   ```

   Both `-10.0` and `0.0` return a successful **10,000 m²** result. NaN and infinity instead fail later with the misleading message `h0 must be finite and positive`, although `h0=30` is valid.

   **Fix:** Validate the sampled field’s shape, finiteness and positivity before assigning bands, consistently with `unresolvable_water`.

6. **Minor — Non-Git source hashing mishandles unreadable or changing files.**

   Location: [provenance.py:113](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:113).

   Each readable file is hashed twice. If the first read succeeds and the second fails, `.encode()` is called on `None`. If the first read fails, the file’s contents are silently omitted while a normal-looking aggregate hash is still reported.

   In-memory fault injection, with Git unavailable and one enumerated Python file:

   ```text
   file_sha256 results [valid_digest, None]
       → AttributeError: 'NoneType' object has no attribute 'encode'

   file_sha256 result None
       → source_sha256 contains an ordinary non-null digest
   ```

   The first case aborts provenance collection; the second presents an incomplete content identity without indicating incompleteness.

   **Fix:** Read each digest once. If any required source cannot be hashed, record a null/incomplete aggregate plus the failed paths, or raise a deliberate, descriptive error.

7. **Nit — The driver still documents the size-floor assumption that round 4 removed.**

   Location: [420_local_refine.py:741](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:741).

   The comment says `FINE_H / 2 is the smallest element anywhere`. A finer base disproves that statement—the reason `H_FLOOR` was introduced. The following call correctly uses `H_FLOOR`.

   **Fix:** Replace the stale explanation with the target/base minimum used by `H_FLOOR`.

## Verdict

VERDICT: FAIL (0 blocker, 0 major, 6 minor, 1 nit)


### Prompt

```markdown
# Review request, round 5: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 23 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, the last commit round 4; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-3: see your round-4 status table.
Round 4 (FAIL 0/4/4/1): all fixed in f07c6d9 (`git show f07c6d9`; triage in
`docs/coastline-rules-review-20260927.md`).
1. `patch.size_lower_bound` (grid readings minus slope x half-diagonal,
   floored); `land_an_element_fits` proves absence by erosion at that bound
   and reports any polygon it cannot clear. Assumption, now documented: the
   sizing field changes by at most `slope=1` m per metre (the real fields:
   about 0.2 for the smoothed base, 0.165 in the ramps).
2. `H_FLOOR` = min(finest target, finest base ambient size).
3. the land check runs on the final rim whatever the rules.
4. slit closure judged on the size bound over the dead end.
5. wall-pocket gate and disc test on the size bound.
6. `provenance.code_state`: path, distribution version, git or source hash.
7. NaN/inf/<=0 sizes and bad factors raise.
8. premise asserted; the driver call site is exercised with rim repair on
   and off.
9. USER_GUIDE: priority has no effect in local refinement.
All six hires recipes rebuild accepted, byte-identical to round 3.
914 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | minor | yes: reviewer's thin-box and wall-pocket reproductions (bound 76.47 > 30; pocket closed over a 59.9 m spot) | yes | `_covering_grid`: keep every grid point within half a diagonal of the geometry, not only those inside; tests `test_lower_bound_covers_the_edges_of_a_thin_box`, `test_pocket_over_a_fine_spot_beside_its_edge_stays_water` |
| 2 | minor | yes: reviewer's 8-node mesh, 192.47 vs 194.62 across 2e-6 m | yes | fixed differently from "make it Lipschitz": `base_size_field` and `patch_sizing` carry `size_bounds(geom, step)`, proved from element vertices (inside, the field is linear per element) and every node that can be nearest to a point of geom (outside); a region adds its target where its ramp reaches. No slope premise. `size_lower_bound` uses it when present. Tests over dense grids around the Voronoi jumps, with and without a region |
| 3 | minor | yes: widths/y (90,150), (100,48), (90,48) all closed | yes | `size_upper_bound`; the pocket's minimum-area and coast-clearance tests use it; parametrized test plus a pocket that still closes |
| 4 | minor | yes: make_valid collection, guard returned [] | yes | `_polygons` recurses collections; test |
| 5 | minor | yes: -10 and 0 returned 10,000 m2; NaN gave the h0 message | yes | sampled field validated before the h0 clip; h0/spacing validated; tests for -10, 0, NaN, inf and the normal path |
| 6 | minor | yes: code path read (double `file_sha256`, `None.encode`) | yes | each file read once; unreadable files listed in `source_unreadable` and `source_sha256` is None; rglob failure recorded; tests |
| 7 | nit | yes: comment read | yes | comment names H_FLOOR |

A first version of fix 2 lowered the bound to a region's target wherever its
ramp reached; all six hires recipes then failed the land check on land wider
than one target element. The bound now follows the ramp
(`target + (lo - target) * u(d_min)`).

After the round-5 fixes all six hires recipes rebuild accepted and
byte-identical to the round-4 rebuild. 931 tests pass.

## Round 6

### Answer (verbatim)

Reviewed `c88739e` read-only. No files changed; no mesh generation or batch jobs ran.

Using the specified Python, the full pytest command failed before collection because no writable temporary directory was available. An initial focused run also encountered cache/temp-file restrictions. After mocking Matplotlib cache setup and skipping tests requiring writes, **262 tests passed, 71 skipped**. Five additional driver regression cases passed with writes mocked. Both working-tree and staged diffs are empty. Mesh QA, counts and time steps were not independently rerun.

**Previous findings**

Statuses concern the original findings; newly discovered limitations are numbered below.

| Finding | Status | Reason |
|---|---|---|
| R1 #1 — Preserve mode | RESOLVED | Rim repair requires `resolve`. |
| R1 #2 — Lagoon entrance | RESOLVED | Connectivity guard and lagoon regression pass. |
| R1 #3 — Refused pocket walls | RESOLVED | Restoration uses pocket geometry; regression passes. |
| R1 #4 — Triangle collapse | RESOLVED | Minimum-ring and duplicate-edge guards remain; regression passes. |
| R1 #5 — Retreat-only retry | RESOLVED | Array comparison and snapshot restoration verified with mocks. |
| R1 #6 — Artificial land / raster padding | PARTIAL | Empty land is handled; variable-size padding remains insufficient—finding 2. |
| R1 #7 — Invalid continuous-width parameters | RESOLVED | Reported invalid parameters and empty footprints are handled in `unresolvable_water`. |
| R1 #8 — Ignored seed argument | RESOLVED | Mocked regression confirms supplied seeds are used. |
| R1 #9 — GPL provenance import | RESOLVED | Uses distribution metadata without importing libraries. |
| R1 #10 — Missing inputs/settings | RESOLVED | Bathymetry, region datasets and effective environment settings are recorded. |
| R1 #11 — Untracked source | RESOLVED | Git status includes untracked files; tracking status is recorded. |
| R1 #12 — Retry reproducibility | RESOLVED | Documentation requires the original seed sequence for pass 2. |
| R1 #13 — Ineffective fixtures | RESOLVED | Original fixtures were corrected; finding 5 concerns a new regression test. |
| R1 #14 — Raster documentation | RESOLVED | Documents one-third of the finest target. |
| R2 #1 — Resolvable land accepted as water | RESOLVED | Recursive land validation rejects it. |
| R2 #2 — All-land raster | RESOLVED | Regression passes. |
| R2 #3 — Detached artificial island | RESOLVED | Detached pieces are refused; regression passes. |
| R2 #4 — Fine region inside pocket | RESOLVED | Bound-based exclusion replaces sparse sampling. |
| R2 #5 — Obsolete sliding curves | RESOLVED | Repair-created and moved points are pinned. |
| R2 #6 — Geometry provenance | RESOLVED | Single-file datasets and shapefile sidecars are enumerated. |
| R2 #7 — Falsy experimental values | RESOLVED | Schema explicitly validates the list. |
| R2 #8 — Seed-variable documentation | RESOLVED | Header names `LR_SEEDS`. |
| R3 #1 — Variable-size land guard | RESOLVED | Structural field bounds replace the unsafe slope assumption. |
| R3 #2 — Repairs bypass validation | RESOLVED | Final-rim and retry checks remain; driver cases pass with mocked writes. |
| R3 #3 — Unchanged new corners | RESOLVED | Changed adjacency causes pinning; regression passes. |
| R3 #4 — Pocket sampling grid | RESOLVED | Replaced by conservative bounds. |
| R3 #5 — Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6 — Imported-package identity | RESOLVED | Actual imported package paths are recorded. |
| R3 #7 — Exhausted iterators | RESOLVED | Both helpers materialize their inputs; regressions pass. |
| R4 #1 — Missed qualifying disc | RESOLVED | Erosion uses the corrected lower bound. |
| R4 #2 — Incorrect size floor | RESOLVED | Floor includes the minimum base ambient size. |
| R4 #3 — Validation disabled with repair | RESOLVED | Both driver switch settings were verified with mocks. |
| R4 #4 — Fine basin behind coarse mouth | RESOLVED | Structural bounds and basin regression pass. |
| R4 #5 — Pocket fine-zone exclusion | RESOLVED | Corrected covering grid closes the reported gap. |
| R4 #6 — Non-Git identity | RESOLVED | Path, distribution and source hash are recorded; enumeration failures are covered by finding 1. |
| R4 #7 — Invalid sizes/factors | RESOLVED | Reported invalid cases raise; regressions pass. |
| R4 #8 — Missing test prerequisites | RESOLVED | Removal is asserted and driver call sites are exercised. |
| R4 #9 — Priority documentation | RESOLVED | Guide states that priority has no effect in local refinement. |
| R5 #1 — Grid coverage | RESOLVED | Boundary-adjacent grid samples restore the covering guarantee. |
| R5 #2 — Invalid slope premise | RESOLVED | Base and patch fields supply bounds from their structure. |
| R5 #3 — Pocket area/clearance | RESOLVED | Both predicates use an upper bound; regressions pass. |
| R5 #4 — Nested land collections | RESOLVED | Land validation recursively traverses polygons; regression passes. |
| R5 #5 — Invalid sampled sizes | RESOLVED | Shape, finiteness and positivity are checked before clipping. |
| R5 #6 — Unreadable source hashing | PARTIAL | Individual files are read once, but enumeration failures remain hidden—finding 1. |
| R5 #7 — Stale floor comment | RESOLVED | Comment now names `H_FLOOR`. |
| Author — Zero-progress retreats | RESOLVED | Threshold regressions pass. |
| Author — Retry without tip retreat | RESOLVED | Mocked regression verifies the second repair and search. |

**Findings**

1. **Minor — Directory-enumeration failures still produce a complete-looking source hash.**

   Location: [provenance.py:115](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:115).

   Python 3.12’s `Path.rglob()` suppresses directory-scanning errors internally. Consequently, the new surrounding `except OSError` does not detect unreadable directories.

   In-memory reproduction using the real source directory:

   ```python
   with patch("fvcom_mesh_tools.provenance.git_state", return_value=None), \
        patch("os.scandir", side_effect=PermissionError("injected")):
       result = code_state("no-such-dist", Path("src/fvcom_mesh_tools"))
   ```

   Result:

   ```text
   source_sha256 = e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
   source_unreadable = []
   ```

   This identifies an unreadable source tree as successfully hashed empty content. An unreadable subtree similarly disappears from the aggregate.

   **Fix:** Enumerate with explicit error reporting, such as `Path.walk(on_error=...)` or controlled `os.scandir()` recursion. Record failed directories and null the aggregate. Test enumeration failures, not only individual-file read failures.

2. **Minor — Continuous-width padding can exclude valid exterior disc centers and turn open water into land.**

   Location: [patch.py:803](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:803).

   Padding uses sizes sampled **inside** the footprint. Disc radii can increase outside it, so this does not establish that all contributing centers are included.

   Reproduction:

   ```python
   size = lambda q: np.minimum(
       200.0, 100.0 + 0.165*np.maximum(np.asarray(q)[:, 0], 0.0)
   )
   added, report = unresolvable_water(
       box(-1000, -4000, 0, 4000), size,
       box(0, -1000, 5, 1000),
       radius_factor=0.75, min_h=60.0, spacing=1.0,
   )
   ```

   Actual results:

   - **9,234.529 m²** becomes land; `(2, 0)` is filled.
   - No disc-center level is found.
   - A disc centered at `(86, 0)`, radius **85.6425 m**, has zero land overlap and contains `(2, 0)`.
   - Widening only the footprint from 5 m to 10 m makes the closure disappear.

   The field is bounded and has slope at most **0.165**. The missing coverage comes from padding, rather than a discontinuity or invalid field.

   **Fix:** Use a proven upper radius covering every potentially contributing exterior center. Where that cannot be established, conservatively leave affected boundary water open.

3. **Minor — `island_rings` crashes on nested polygon collections.**

   Location: [patch.py:1922](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1922).

   The recursive traversal added for land validation was not applied here. Passing the `make_valid()` geometry from `test_land_guard_reaches_polygons_nested_in_collections`, with water `box(-1000,-1000,1000,1000)` and constant size 30, raises:

   ```text
   AttributeError: 'MultiPolygon' object has no attribute 'exterior'
   ```

   A non-area member inside the water also reaches the same unchecked attribute access. Thus valid collection output from geometry repair cannot reach island construction.

   **Fix:** Iterate through `_polygons(land)` and ignore non-area members. Add a regression using the same nested collection as the land-validation test.

4. **Minor — An empty footprint crashes the local shoreline filter before its empty-domain behavior can apply.**

   Location: [patch.py:948](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:948).

   Reproduction:

   ```python
   filter_shoreline_local(
       box(0, 0, 100, 100),
       lambda q: np.full(len(q), 30.0),
       30.0, Polygon(),
   )
   ```

   This raises `ValueError: cannot convert float NaN to integer` when converting empty bounds into raster dimensions. The empty-footprint handling in `unresolvable_water` cannot help this caller.

   **Fix:** Handle emptiness before reading bounds. Consistent with the documented outside-footprint behavior, return the finest-band filtering with no local bands; alternatively reject empty footprints explicitly with a descriptive validation error.

5. **Minor — The new fine-spot pocket regression passes when its lower-bound protection is deliberately broken.**

   Location: [test_review_r5_bounds.py:38](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/tests/test_review_r5_bounds.py:38).

   This test still passes under:

   ```python
   with patch("fvcom_mesh_tools.patch.size_lower_bound", return_value=1000.0):
       test_pocket_over_a_fine_spot_beside_its_edge_stays_water()
   ```

   The fixture’s upper bound is **157.74284 m**, requiring **24,882.804 m²**, while its pocket has only **10,000 m²**. The newly added area guard therefore rejects it independently of fine-spot detection. This coverage defect was introduced by the round-5 fixes.

   **Fix:** Make the other predicates pass. Capping this fixture’s field at 90 m works: I verified that the real implementation leaves the pocket open, while the deliberately wrong lower bound closes it. Include that negative control.

6. **Nit — Driver comments still describe default rules as opt-in.**

   Locations: [420_local_refine.py:261](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:261), [420_local_refine.py:537](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:537), [420_local_refine.py:567](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:567).

   “None unless named” and the `keep_base_land`/`wall_pockets` “opt-in” comments contradict the defaults applied immediately below.

   **Fix:** Describe the four default-enabled rules and distinguish them from experimental overrides.

## Verdict

VERDICT: FAIL (0 blocker, 0 major, 5 minor, 1 nit)


### Prompt

```markdown
# Review request, round 6: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 24 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, f07c6d9 round 4, the last commit round 5; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-4: see your round-5 status table (all RESOLVED except the PARTIALs
that round-5 findings 1, 2 and 4 cover).
Round 5 (FAIL 0/0/6/1): all fixed in the last commit (`git show HEAD`;
triage in `docs/coastline-rules-review-20260927.md`).
1. `_covering_grid`: every grid point within half a diagonal of the
   geometry is read, not only those inside it.
2. Not "make the extension Lipschitz": `base_size_field` and
   `patch_sizing` fields carry `size_bounds(geom, step)`, proved without a
   slope premise -- inside the base mesh from the vertices of the elements
   the geometry meets (the field is linear per element), outside from every
   node that can be nearest to a point of the geometry; a region lowers the
   bound to `target + (lo - target) * u(d_min)` where the base bound is over
   its target. `size_lower_bound` uses a field's own `size_bounds` when it
   has one; the generic path keeps the documented slope premise.
   (A first version lowered the bound to the target wherever a ramp
   reached; all six recipes then failed the land check on land wider than
   one target element, which is why the ramp is followed.)
3. `size_upper_bound`; the wall-pocket area and coast-clearance tests use it.
4. `_polygons` recurses into collections.
5. `filter_shoreline_local` validates the sampled field (and h0, spacing).
6. `code_state` reads each file once; unreadable files go to
   `source_unreadable` and `source_sha256` is null.
7. stale comment fixed.
Regression tests: `tests/test_review_r5_bounds.py` (your reproductions).
All six hires recipes rebuild accepted, byte-identical to round 4.
931 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | minor | yes: `os.scandir` raising under `Path.rglob` hashed as empty (e3b0...) | yes | `Path.walk(on_error=...)`; failed directories go to `source_unreadable`, digest null; test injects a `PermissionError` from `os.scandir` |
| 2 | minor | yes: reviewer's field and footprint; 9,234 m2 closed | yes | the pad grows to a fixed point of `radius_factor * size_upper_bound(padded box) + 2 s` (overshoot 1.25); a pad that does not settle (growth not shrinking, 8 steps, or a raster over 5e7 cells) closes nothing and reports `unbounded_pad`. Tests: the reproduction, a dead end still closed, a diverging field. (A first cut without the divergence stop reached 378 GB in the test and was killed.) |
| 3 | minor | yes: `AttributeError: 'MultiPolygon' object has no attribute 'exterior'` | yes | `island_rings` iterates `_polygons(land)`; test |
| 4 | minor | yes: `cannot convert float NaN to integer` | yes | `filter_shoreline_local` refuses an empty footprint with a clear error; test |
| 5 | minor | yes: test passes with `size_lower_bound` patched to 1000 | yes | field capped at 90 m; negative control with the lower bound broken closes the pocket |
| 6 | nit | yes: comments read | yes | comments say the four rules are on by default |

The first rebuild with fix 2 moved every recipe's coast by 50-95 m and two
recipes failed QA: the raster origin moved with the pad and re-quantised the
land (lost cells went UP in some, which new centres cannot do). The origin
now stays where the footprint's own sizes put it and the pad grows by whole
cells, so a larger pad only adds centres (test
`test_a_larger_pad_only_adds_cells_around_the_same_raster`). After that all
six hires recipes rebuild accepted and byte-identical to round 5. 939 tests
pass.

## Round 7

### Answer (verbatim)

Reviewed `806d011` read-only. No files changed; no mesh generation or batch jobs ran. Git working-tree and staged diffs remain empty.

The full pytest invocation failed before collection because the sandbox has no writable temporary directory. An unadapted focused run produced 186 passes, 36 failures and 5 errors while encountering temporary-file/cache restrictions. With Matplotlib’s existing cache supplied in memory and write-dependent tests skipped, **269 tests passed, 72 skipped**. Another **9 driver cases passed with writes mocked**. Mesh QA, node/element counts and time steps were not independently rerun.

**Previous findings**

“PARTIAL” means the original reproduction is addressed but a substantiated variant remains.

| Finding | Status | Reason |
|---|---|---|
| R1 #1 — Preserve mode | RESOLVED | Geometry-changing rim repair requires `resolve`. |
| R1 #2 — Lagoon entrance | RESOLVED | Connectivity guard and lagoon regression pass. |
| R1 #3 — Refused pocket walls | RESOLVED | Restoration checks pocket geometry against remaining water. |
| R1 #4 — Ring collapse | PARTIAL | Three-node removal is blocked; four-node rings can still become invalid—finding 4. |
| R1 #5 — Retreat-only retry | RESOLVED | Array comparison detects changes; mocked search/restoration passes. |
| R1 #6 — Artificial land / padding | PARTIAL | Empty land and the smooth-field example are fixed; discontinuous fields still expose insufficient padding—finding 2. |
| R1 #7 — Invalid continuous-width parameters | RESOLVED | Reported invalid parameters and empty footprints are handled. |
| R1 #8 — Ignored seed argument | RESOLVED | Mocked regression confirms the supplied seed list is used. |
| R1 #9 — GPL provenance import | RESOLVED | Distribution metadata replaces imports. |
| R1 #10 — Missing provenance inputs/settings | RESOLVED | Bathymetry, geometry datasets and effective settings are recorded. |
| R1 #11 — Untracked source | RESOLVED | Git status includes untracked files; tracking status is recorded. |
| R1 #12 — Retry reproducibility | RESOLVED | Pass-2 replay requires the original seed sequence. |
| R1 #13 — Ineffective fixtures | RESOLVED | Original fixtures were corrected; unreachable merge operation removed. |
| R1 #14 — Raster documentation | RESOLVED | Guide specifies one-third of the finest target. |
| R2 #1 — Resolvable land accepted as water | RESOLVED | Recursive land validation rejects qualifying land in the hole. |
| R2 #2 — All-land raster | RESOLVED | Empty-water early return and regression pass. |
| R2 #3 — Detached artificial island | RESOLVED | Pieces adjacent to no land are refused. |
| R2 #4 — Fine region inside pocket | PARTIAL | Explicit certified floors work; the default floor still permits erasure—finding 3. |
| R2 #5 — Obsolete sliding curves | RESOLVED | Repair-created and moved points are pinned. |
| R2 #6 — Geometry provenance | RESOLVED | Single-file datasets and shapefile sidecars are enumerated. |
| R2 #7 — Falsy experimental values | RESOLVED | Schema explicitly validates the list. |
| R2 #8 — Seed-variable documentation | RESOLVED | Header names `LR_SEEDS`. |
| R3 #1 — Variable-size land guard | RESOLVED | Conservative bounds replace representative-point decisions. |
| R3 #2 — Repairs bypass validation | RESOLVED | Final-rim and retry checks remain; mocked driver cases pass. |
| R3 #3 — Unchanged new corners | RESOLVED | Changed adjacency causes pinning; regression passes. |
| R3 #4 — Pocket sampling grid | PARTIAL | Covering bounds are implemented, but an unjustified default floor overrides them—finding 3. |
| R3 #5 — Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6 — Imported-package identity | RESOLVED | Actual imported paths and non-Git source identity are recorded. |
| R3 #7 — Exhausted iterators | RESOLVED | Both helpers materialize their inputs; regressions pass. |
| R4 #1 — Missed qualifying disc | RESOLVED | Land validation uses conservative erosion with corrected bounds. |
| R4 #2 — Incorrect driver size floor | RESOLVED | `H_FLOOR` includes the minimum base ambient size. |
| R4 #3 — Validation disabled with repair | RESOLVED | Both driver switch settings pass mocked checks. |
| R4 #4 — Fine basin behind coarse mouth | PARTIAL | Explicit-floor regression passes; default `rim_repair` can still erase such a basin—finding 3. |
| R4 #5 — Pocket fine-zone exclusion | PARTIAL | Correct with a certified floor; default API path remains unsafe—finding 3. |
| R4 #6 — Non-Git identity | RESOLVED | Path, distribution and source hash are retained. |
| R4 #7 — Invalid sizes/factors | RESOLVED | Reported invalid cases raise; regressions pass. |
| R4 #8 — Missing test prerequisites | RESOLVED | Removal is asserted and driver call sites are exercised. |
| R4 #9 — Priority documentation | RESOLVED | Guide states that priority has no effect in local refinement. |
| R5 #1 — Grid coverage | RESOLVED | Boundary-adjacent samples provide the covering guarantee. |
| R5 #2 — Invalid slope premise for bounds | RESOLVED | Base and patch fields supply structural bounds. Finding 2 concerns the separate exterior-padding premise. |
| R5 #3 — Pocket area/clearance | RESOLVED | Both predicates use an upper bound; regressions pass. |
| R5 #4 — Nested land collections | RESOLVED | Land validation recursively traverses polygons. |
| R5 #5 — Invalid sampled sizes | RESOLVED | Shape, finiteness and positivity are checked before clipping. |
| R5 #6 — Unreadable source hashing | RESOLVED | Files are read once; scan failures now null the digest. |
| R5 #7 — Stale floor comment | RESOLVED | Comment names `H_FLOOR`. |
| R6 #1 — Hidden enumeration failures | RESOLVED | Injected `os.scandir` failure produces `source_sha256=None` and a recorded unreadable directory. |
| R6 #2 — Exterior disc centres | PARTIAL | Smooth-field reproduction passes, but the exterior slope premise is not established for discontinuous fields—finding 2. |
| R6 #3 — Nested island collections | RESOLVED | `island_rings` recursively visits polygons; regression passes. |
| R6 #4 — Empty filter footprint | RESOLVED | Explicit, descriptive rejection; regression passes. |
| R6 #5 — Ineffective fine-spot test | RESOLVED | Revised fixture and broken-bound negative control both pass. |
| R6 #6 — Opt-in comments | RESOLVED | Those comments now describe the defaults correctly. Finding 6 concerns a different stale statement. |

**Findings**

1. **Major — Fully free islands silently retain the base coastline under `resolve`.**

   Location: [patch.py:1206](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1206).

   The `is_free.all()` branch copies the base ring without consulting `coastline`, `shoreline` or local spacing. This pre-existing behavior remains reachable through the audited driver.

   **Reproduction:** Use a frozen outer square `box(-1000,-1000,1000,1000)`, an entirely free base-island ring around `box(-100,-100,100,100)`, source land `box(-50,-50,50,50)`, and constant `h=30`. Call `rim_constraints(..., coastline="resolve", shoreline=source.boundary)`.

   The returned island remains the **200 × 200 m base square**, omitting **30,000 m² of source water**. Subsequent actual calls to `island_rings` report no addition or refusal, and `land_an_element_fits(hole.intersection(source), ...)` returns `[]`. The validation checks land wrongly made water, so this opposite error escapes it. Geometry QA does not establish source-coastline fidelity.

   **Fix:** Resolve entirely free closed rings against the filtered source and resample them at local size. Preserve the base outline when the selected mode requires preservation. Add a regression covering a source island smaller than its base counterpart.

2. **Minor — The round-6 padding fix can still manufacture land with valid structural size bounds.**

   Location: [patch.py:805](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:805), [patch.py:825](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:825).

   Settling a bound over the current padded box does not establish the stated exterior slope premise. In particular, the driver’s `outside="nearest"` field can be discontinuous; supplying structural bounds does not make it Lipschitz.

   **Reproduction:** Set land to `box(-1000,-4000,0,4000)`, footprint to `box(0,-1000,5,1000)`, and:

   ```python
   def size(q):
       x = np.atleast_2d(q)[:, 0]
       return np.where(x < 70, 100., np.where(x < 120, 40., 200.))

   size.size_bounds = lambda g, step: (
       40., 100. if g.bounds[2] < 120 else 200.
   )
   ```

   These are valid conservative bounds. With `radius_factor=.75`, `min_h=60`, and `spacing=1`, the function adds **9,974.18 m² of land** and reports `unbounded_pad=False`.

   Yet the centre `(150,0)` has radius **150 m**, is **150 m** from land, and its disc covers `(2,0)`. Supplying the equally valid global upper bound of 200 enlarges the pad and produces **zero closure**.

   **Fix:** Use a certified global maximum for discontinuous fields—the base ambient maximum bounds the driver’s field—or establish an applicable exterior bound explicitly. Without that guarantee, refuse closure. Retain the aligned raster origin.

3. **Minor — Default size-floor estimates can override correct bounds and erase fine water.**

   Locations: [walls.py:690](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:690), [patch.py:1519](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1519).

   Both helpers default to half the smallest boundary-vertex size. That is not a guaranteed minimum over the interior. `size_lower_bound` then clamps its result upward to this estimate, defeating the conservative sampling.

   **Pocket reproduction:**

   ```python
   wall = LineString([(0,0), (0,500), (150,500), (150,200), (40,200)])
   land = box(-2000,-2000,2000,0)
   c = np.array([75., 350.])

   def size(q):
       return np.minimum(190., 30. + np.linalg.norm(np.atleast_2d(q)-c, axis=1))
   ```

   This field has slope at most 1. Calling `close_wall_pockets(..., min_h=60)` closes **45,000 m²**, including the point where `h=30`. Supplying `size_floor=30` closes nothing.

   I reproduced the same failure in slit repair with cyclic points:

   ```python
   [(0,0),(500,0),(500,500),(260,500),(260,540),(325,540),
    (325,840),(175,840),(175,540),(240,540),(240,500),(0,500)]
   ```

   Use the same field formula centred at `(250,690)`, all-free vertices, `operations=("slits",)`, and `retreat_tips=False`. The default closes one slit and removes the basin centre; `size_floor=30` preserves it.

   **Fix:** Require a certified floor for variable fields, or derive a conservative bound without clamping to a boundary estimate. Refuse closure when absence cannot be proved. The driver’s explicit `H_FLOOR` avoids this particular default-path defect.

4. **Minor — Short-edge repair can return a zero-area, invalid ring.**

   Location: [patch.py:1599](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1599).

   The round-1 guard preserves three vertices but does not ensure that they form a valid polygon. The intersection check also excludes adjacent edges that can overlap the replacement.

   **Reproduction:** Call default `rim_repair` with:

   ```python
   p = np.array([(0,0), (10,0), (20,0), (10,1)], float)
   e = np.array([(0,1), (1,2), (2,3), (3,0)])
   pfix_base = np.array([0, 1, 2, -1])
   # water = Polygon(p); constant h = 30
   ```

   The input polygon is valid with area **10 m²**. Repair removes `(10,1)` and returns three collinear points. `hole_polygon` then returns **area 0, `is_valid=False`**, without repair rejecting the operation.

   **Fix:** Validate each proposed affected ring before committing removal: positive area, simple boundary and no overlapping replacement edges. Reject a degenerate candidate while preserving the original constraints.

5. **Minor — Closed-ring arc selection can miss almost the entire stretch and depend on ring start position.**

   Locations: [patch.py:579](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:579), [patch.py:606](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:606).

   The 32 probes are uniformly spaced by vertex index. Uneven vertex density can hide the important portion of the stretch from the new arc-selection logic.

   **Reproduction:**

   ```python
   ring = LineString([(0,1000),(0,0),(1000,0),(1000,1000),(0,1000)])
   pts = np.vstack([
       np.column_stack([np.full(160, -20.), np.linspace(100, 101, 160)]),
       [(0,1000),(1000,1000),(1000,0),(-20,0)],
   ])
   ```

   `_source_substring(pts, ring)` selects the **100 m** arc for a **3,920.22 m** stretch; an input vertex is **1,345.36 m** from that arc. Rotating the same ring’s start coordinate to `(0,0)` selects the correct **3,900 m** arc, whose maximum input-vertex departure is **20 m**.

   **Fix:** Include every stretch vertex and samples spaced by arc length in the fit assessment. Make the selected geometric arc independent of the source ring’s arbitrary starting coordinate, with a regression for this uneven-density case.

6. **Nit — The driver still describes a shoreline signed-distance operation that it does not perform.**

   Location: [420_local_refine.py:475](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:475).

   The comment says OceanMesh’s processed shoreline signed-distance function becomes part of the fill domain. In fact, `_shl` is only inspected for logging; `shore` comes from `_filtered`, and `fd` at line 1352 uses the constructed hole. Its own docstring explains why the claimed intersection was abandoned.

   **Fix:** Describe the implemented filter → rim → hole path and remove the obsolete signed-distance claim.

## Verdict

VERDICT: FAIL (0 blocker, 1 major, 4 minor, 1 nit)


### Prompt

```markdown
# Review request, round 7: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 26 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, f07c6d9 round 4, c88739e round 5, the last commit round 6; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-5: see your round-6 status table (all RESOLVED except R1 #6 and
R5 #6, PARTIAL, which round-6 findings 2 and 1 cover).
Round 6 (FAIL 0/0/5/1): all fixed in the last commit (`git show HEAD`;
triage in `docs/coastline-rules-review-20260927.md`).
1. `code_state` walks with `Path.walk(on_error=...)`; unscannable
   directories go to `source_unreadable`, digest null.
2. `unresolvable_water`: the pad grows to a fixed point of
   `radius_factor * size_upper_bound(padded box) + 2 s` (overshoot 1.25);
   a pad that does not settle (growth not shrinking, 8 steps, raster over
   5e7 cells) closes nothing (`unbounded_pad`). Beyond the settled pad the
   premise is radius_factor x slope < 1, documented. The raster origin
   stays where the footprint's sizes put it and the pad grows by whole
   cells, so a larger pad only adds centres (a first version moved the
   origin, re-quantised the land and moved every recipe's coast 50-95 m).
3. `island_rings` iterates `_polygons`.
4. `filter_shoreline_local` refuses an empty footprint.
5. fine-spot pocket test capped at 90 m, negative control added.
6. driver comments.
Tests: `tests/test_review_r6.py`, `tests/test_review_r5_bounds.py`.
All six hires recipes rebuild accepted, byte-identical to round 5.
939 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | major | yes: code path read (`is_free.all()` copies the base ring whatever the mode); regression built on a 15 x 15 grid with a 200 m base island and a 100 m source island | yes (pre-existing; no hires recipe has an interior island) | under `resolve` with a shoreline, `rim_constraints` leaves a wholly free island out (`n_islands_left_to_source`) and the driver's `island_rings` adds the source's land inside the hole, which the land check then covers; `preserve`/`resample` still copy it. Tests for both |
| 2 | minor | yes: reviewer's stepped field with `size_bounds`, 9,974 m2 closed | yes | `patch_sizing` / `base_size_field` state a certified `size_min` / `size_max`; `unresolvable_water` checks annuli doubling out to `radius_factor * size_max` (a centre beyond ring R is R from the footprint, so the discs up to 2R must be under R); a field with `size_bounds` but no `size_max` closes nothing. Measured on the real 3,210-node base: 3.5-12 s per bound, settles after one doubling. Tests for both |
| 3 | minor | yes: pocket (45,000 m2 closed over h = 30) and slit reproductions | yes | `resolve_size_floor`: caller's floor, else the field's `size_min`, else an estimate that only sets the grid step; `size_lower_bound(certified=False)` and `land_an_element_fits(certified=False)` do not raise the bound to it (a bound at or below zero reads as "cannot prove absence"). Tests: both reproductions |
| 4 | minor | yes, but the degenerate ring came from `blunt_acute_corners` (the frozen-corner opening), not short-edge removal; verified by running each operation alone | yes | `_ring_is_polygon`: short-edge removal, frozen-corner opening and chord cuts each refuse a change that leaves their ring not a simple polygon of positive area. Tests: reproduction and the normal removal |
| 5 | minor | yes: 100 m arc chosen for a 3,920 m stretch; rotating the ring start changed the answer | yes | probes are every vertex plus 65 points evenly spaced by arc length; the arc choice is shorter vs longer (longer only when it fits twice as well), not direct vs wrap. Test on both ring starts; the existing Yokohama-arc test still passes |
| 6 | nit | yes: comment read | yes | comment describes filter -> rim -> hole; Shoreline for the log only |

All six hires recipes rebuild accepted and byte-identical to round 6. 948 tests pass.

## Round 8

### Answer (verbatim)

Reviewed `e88208f` read-only. No files changed, batch jobs submitted, or mesh generation run.

The full pytest command failed before collection because no temporary directory was writable. An adjusted run initially produced **204 passed, 37 failed, 15 errors** in filesystem/cache-dependent paths. With Matplotlib setup mocked in memory and temporary-file-dependent tests skipped, the scoped suite produced **231 passed, 34 skipped**, including all nine round-7 regressions. Another **26 in-memory hires-schema checks passed**. The author’s full-suite and compute-node results were not independently reproduced.

**Previous findings**

“RESOLVED” below refers to the reported defect; related new cases are identified separately.

| Finding | Status | Reason |
|---|---|---|
| R1 #1 — Preserve mode | RESOLVED | Rim repair requires `resolve`. |
| R1 #2 — Lagoon entrance | RESOLVED | Connectivity guard preserves the reproduced entrance. |
| R1 #3 — Refused pocket walls | RESOLVED | Restoration tests pocket geometry against remaining water. |
| R1 #4 — Ring collapse | RESOLVED | Original triangle and degenerate-ring cases are guarded. See finding 1 for a different operation. |
| R1 #5 — Retreat-only retry | RESOLVED | Array comparison detects changes; restoration uses `finally`. |
| R1 #6 — Artificial land/padding | RESOLVED | Empty-land handling and certified exterior bounds address the reported cases. |
| R1 #7 — Invalid continuous-width parameters | RESOLVED | Reported parameters and empty footprints are handled. |
| R1 #8 — Ignored seed argument | RESOLVED | Search iterates `seed_list`. |
| R1 #9 — GPL provenance import | RESOLVED | Versions use distribution metadata. |
| R1 #10 — Missing provenance inputs/settings | PARTIAL | Inputs/settings are recorded, but enumeration errors can omit sidecars—finding 5. |
| R1 #11 — Untracked source | RESOLVED | Untracked files and entry-point tracking status are recorded. |
| R1 #12 — Retry reproducibility | RESOLVED | Documentation requires the original seed sequence for pass 2. |
| R1 #13 — Ineffective fixtures | RESOLVED | Fixtures assert the relevant operations; unreachable merge code was removed. |
| R1 #14 — Raster documentation | RESOLVED | Resolution is documented relative to the finest target. |
| R2 #1 — Resolvable land accepted as water | RESOLVED | Recursive land validation checks qualifying land inside the hole. |
| R2 #2 — All-land raster | RESOLVED | Empty-water early return works. |
| R2 #3 — Detached artificial island | RESOLVED | Closure requires adjacent land. |
| R2 #4 — Fine region inside pocket | RESOLVED | Estimated floors no longer override conservative bounds. |
| R2 #5 — Obsolete sliding curves | RESOLVED | Repair-created and moved points are pinned. |
| R2 #6 — Geometry provenance | PARTIAL | Normal enumeration works; failed enumeration is silently incomplete—finding 5. |
| R2 #7 — Falsy experimental values | RESOLVED | Explicit list validation rejects them. |
| R2 #8 — Seed-variable documentation | RESOLVED | Header names `LR_SEEDS`. |
| R3 #1 — Variable-size land guard | RESOLVED | Decisions use conservative bounds. |
| R3 #2 — Repairs bypass validation | RESOLVED | Final-rim and retry land checks remain wired. |
| R3 #3 — Unchanged new corners | RESOLVED | Changed adjacency also causes pinning. |
| R3 #4 — Pocket sampling grid | RESOLVED | Covering samples and uncertified-floor handling address both gaps. |
| R3 #5 — Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6 — Imported-package identity | RESOLVED | Actual imported paths are recorded. |
| R3 #7 — Exhausted iterators | RESOLVED | Both helpers materialize their inputs. |
| R4 #1 — Missed qualifying disc | RESOLVED | Conservative erosion replaces sampled absence decisions. |
| R4 #2 — Incorrect driver size floor | RESOLVED | `H_FLOOR` includes the base ambient minimum. |
| R4 #3 — Validation disabled with repair | RESOLVED | Final land validation is outside the repair switch. |
| R4 #4 — Fine basin behind coarse mouth | RESOLVED | Default estimated floors no longer erase the reproduced basin. |
| R4 #5 — Pocket fine-zone exclusion | RESOLVED | Default estimated floors no longer override the fine-zone bound. |
| R4 #6 — Non-Git identity | RESOLVED | Path, distribution and source digest are retained. |
| R4 #7 — Invalid sizes/factors | RESOLVED | Reported invalid cases are rejected. |
| R4 #8 — Missing test prerequisites | RESOLVED | Tests assert prerequisites and exercise driver guards. |
| R4 #9 — Priority documentation | RESOLVED | Guide states that local refinement ignores priority. |
| R5 #1 — Grid coverage | RESOLVED | Boundary-adjacent samples provide coverage. |
| R5 #2 — Invalid slope premise | RESOLVED | Base and patch fields supply structural bounds. |
| R5 #3 — Pocket area/clearance | RESOLVED | Both predicates use upper bounds. |
| R5 #4 — Nested land collections | RESOLVED | Land validation recursively traverses polygons. |
| R5 #5 — Invalid sampled sizes | RESOLVED | Validation precedes clipping. |
| R5 #6 — Unreadable source hashing | RESOLVED | Failed source reads/scans invalidate the digest. |
| R5 #7 — Stale floor comment | RESOLVED | Comment names `H_FLOOR`. |
| R6 #1 — Hidden source-enumeration failures | RESOLVED | `Path.walk(on_error=...)` records failures. |
| R6 #2 — Exterior disc centres | RESOLVED | Certified global maxima bound the annulus search; unsupported bounded fields refuse closure. |
| R6 #3 — Nested island collections | RESOLVED | `island_rings` recursively visits polygons. |
| R6 #4 — Empty filter footprint | RESOLVED | Explicit descriptive rejection. |
| R6 #5 — Ineffective fine-spot test | RESOLVED | Revised fixture and negative control exercise the bound. |
| R6 #6 — Opt-in comments | RESOLVED | Comments describe the defaults correctly. |
| R7 #1 — Fully free islands | PARTIAL | Replacement works when source shoreline remains; filtered-away islands still abort upstream—finding 3. |
| R7 #2 — Discontinuous exterior sizes | RESOLVED | Both new padding regressions pass. |
| R7 #3 — Estimated default floors | RESOLVED | Pocket and slit reproductions now preserve fine water. |
| R7 #4 — Degenerate ring | RESOLVED | Original reproduction stays valid; ordinary removal still works. |
| R7 #5 — Closed-ring arc selection | PARTIAL | Uneven-density case is fixed; equal-length arcs remain start-dependent—finding 2. |
| R7 #6 — Obsolete driver comment | RESOLVED | Comment describes the implemented filter → rim → hole path. |

**Findings**

1. **Minor — Slit endpoint snapping can return a self-intersecting rim.**

   Location: [patch.py:1819](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1819).

   The slit path checks the original throat, then may replace it with a different chord by snapping to `other_end`. It commits that chord without the polygon-validity guard added to the other operations.

   Reproduced in memory with:

   ```python
   p = np.array([
       (0,0), (0,5), (20,5), (20,2), (10,2),
       (10,-2), (8,-1), (8,-20), (-100,-20), (-100,0)
   ], float)
   e = np.column_stack([np.arange(10), np.roll(np.arange(10), -1)])
   result = rim_repair(
       p, e, np.full(10, -1), Polygon(p),
       lambda q: np.full(len(q), 30.),
       protect=(5, 6), operations=("slits",), rounds=1,
       retreat_tips=False, focus=[[0, 0]], focus_factor=.001,
   )
   ```

   The input is valid, with area **2,243 m²**. Repair reports one closed slit, but `hole_polygon(result[0], result[1])` returns `is_valid=False`, with **`Self-intersection[8 -1.6]`**. The snapped chord crosses the retained edge beside the protected corner.

   **Fix:** Validate the proposed slit ring after endpoint snapping and before committing it. Try the unsnapped throat or refuse the operation when invalid. Also make `hole_polygon` explicitly reject invalid output. This is an existing slit-path gap left uncovered by round 7.

2. **Minor — Equal-length source arcs still depend on the ring’s starting vertex.**

   Location: [patch.py:612](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:612).

   When both arcs have equal length, `<=` designates the direct arc as “short.” The subsequent factor-of-two rule preserves that arbitrary choice unless the other arc fits dramatically better.

   Reproduction:

   ```python
   pts = np.array([(0,0), (60,40), (100,100)], float)
   a = LineString([(0,0), (100,0), (100,100), (0,100), (0,0)])
   b = LineString([(100,0), (100,100), (0,100), (0,0), (100,0)])
   ```

   `_source_substring(pts, a)` returns the **bottom/right** arc; the same call with `b` returns the **left/top** arc. Both are 200 m long; only the source ring’s starting coordinate changed.

   **Fix:** Handle equal or numerically indistinguishable lengths explicitly: choose by fit, then use a canonical geometric tie-break if necessary. Add rotated-start and reversed-orientation regressions. Round-7 finding 5 remains partial.

3. **Minor — The driver rejects a wholly free island that the filter correctly removes.**

   Location: [420_local_refine.py:616](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:616).

   The new island behavior permits an empty source to replace a free base island with water, but the driver exits before reaching it whenever filtering removes the last nearby shoreline.

   Using the round-7 island mesh, source land `box(695,695,705,705)`, constant `h=30`, and the driver’s fine-band land-width settings, the filter removes the **100 m²** source island. Executing the actual shoreline-update block then raises:

   > the h0 = 30 m filter left no shoreline near the free rim

   Calling the fixed `rim_constraints(..., coastline="resolve", shoreline=[])` directly succeeds: it reports one island left to the source and constructs the **640,000 m²** water hole.

   **Fix:** Permit empty source shoreline when the affected physical rings are wholly free interior islands. Pass an empty shoreline collection through to island replacement; retain appropriate refusal for anchored coastline stretches. Cover the earlier “no raw shoreline nearby” branch too.

4. **Minor — Empty land crashes `close_wall_pockets`.**

   Location: [walls.py:649](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:649), [walls.py:661](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:661).

   Reproduction:

   ```python
   wall = LineString([(0,0), (0,400), (110,400), (110,290), (20,290)])
   close_wall_pockets(
       [wall], Polygon(), lambda q: np.full(len(q), 100.),
       min_h=60.000001, size_floor=30.,
   )
   ```

   This raises **`AttributeError: 'NoneType' object has no attribute 'is_empty'`**. Unioning empty land yields an empty geometry collection whose boundary is `None`. With mainland `box(-2000,-2000,2000,0)`, the identical wall successfully closes **12,100 m²**.

   Empty filtered land is possible when all source structures become walls.

   **Fix:** Represent absent coastline with an empty line geometry and treat clearance from absent land as infinite. Merely filtering out `None` still leaves the later empty-land distance comparison returning NaN.

5. **Minor — Failed shapefile enumeration silently produces an incomplete provenance hash.**

   Location: [provenance.py:55](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:55).

   Injecting `PermissionError` into `Path.iterdir()` makes `dataset_files("/readable/land.shp")` return only the `.shp` path. `collect()` then produces a non-null aggregate hash, with no indication that `.dbf`, `.prj`, `.shx` or `.cpg` enumeration failed.

   Directory listing can fail while known files remain readable. Consequently, consumed attributes or CRS can change without changing the recorded fingerprint. This gap is in the sidecar-enumeration fix; the source-tree walker now handles the analogous failure correctly.

   **Fix:** Propagate enumeration failure or record an explicit incomplete-dataset status and null digest. Distinguish “sidecar absent” from “could not determine which sidecars exist.”

6. **Minor — Feedback retries overwrite first-pass rejected-mesh diagnostics.**

   Location: [420_local_refine.py:1984](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1984).

   Both searches reuse the same seeds and write `rejected_seed<seed>.npz`. Executing the actual rejection block twice with seed 3 and an in-memory save replacement produced:

   ```text
   /virtual/output/rejected_seed3.npz
   /virtual/output/rejected_seed3.npz
   ```

   Only the second attempt’s arrays remained. The attempts report can therefore retain two failures while the first failure’s diagnostic mesh has been lost.

   **Fix:** Include search pass and attempt identity in diagnostic filenames, and record the filename/pass alongside each attempt.

7. **Nit — The fill log prints stale constraint counts.**

   Location: [420_local_refine.py:1385](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1385).

   `rc["n_pfix"]` and `rc["n_egfix"]` describe the original `rim_constraints` result. Island insertion, wall construction and rim repair subsequently change the arrays without updating those counters.

   In the round-7 island fixture, the logged fixed-point count remains **32** after island insertion increases the rim to **44** points.

   **Fix:** Print `len(PFIX_ALL)` and `len(EGFIX_ALL)` for the constraints actually passed to the filler.

8. **Nit — A failed candidate is logged as “accepted” before the acceptance gate.**

   Location: [420_local_refine.py:2219](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:2219).

   When every seed introduces violations, the best failing candidate still reaches the unconditional `say(f"accepted seed {seed}")`. The driver later rejects it at line 2309. The `ACCEPTED` marker remains correctly gated, but the log gives a contradictory acceptance statement.

   **Fix:** Say “selected seed” at this point, or move “accepted seed” after all gates and successful marker creation.

## Verdict

VERDICT: FAIL (0 blocker, 0 major, 6 minor, 2 nit)


### Prompt

```markdown
# Review request, round 8: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 27 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, f07c6d9 round 4, c88739e round 5, 806d011 round 6, the last commit round 7; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-6: see your round-7 status table.
Round 7 (FAIL 0/1/4/1): all fixed in the last commit (`git show HEAD`;
triage in `docs/coastline-rules-review-20260927.md`).
1. `rim_constraints`: under `resolve` with a shoreline a wholly free island
   is left out (`n_islands_left_to_source`); the driver's `island_rings`
   adds the source's land inside the hole and the land check covers it.
   `preserve` / `resample` still copy it.
2. certified `size_min` / `size_max` on `patch_sizing` and
   `base_size_field`; `unresolvable_water` checks annuli doubling out to
   `radius_factor * size_max`; a field with `size_bounds` and no
   `size_max` closes nothing. A plain callable keeps the documented slope
   premise (radius_factor x slope < 1).
3. `resolve_size_floor`; `certified=False` bounds are never raised to an
   estimated floor.
4. `_ring_is_polygon` guards short-edge removal and both blunting paths
   (your reproduction's degenerate ring came from opening a frozen corner
   in `blunt_acute_corners`, run by the `angles` operation).
5. `_source_substring`: vertex + arc-length probes; shorter vs longer arc.
6. driver comment.
Tests: `tests/test_review_r7.py`. All six hires recipes rebuild accepted,
byte-identical to round 6. 948 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | minor | yes: reviewer's rim; hole self-intersects at (8, -1.6) | yes | the slit tries the snapped chord, then the unsnapped throat, and takes the first that leaves its ring a simple polygon and the whole rim a valid hole; else it is reported ("closing it would not leave a simple ring"). `hole_polygon` raises on an invalid result. Tests: reproduction, a bow tie |
| 2 | minor | yes: 200 m arcs, answer changed with the ring start | yes | equal-length arcs (within 1e-9 of the ring) are chosen by fit, then by the midpoint's coordinates; test with two starts and the reversed orientation |
| 3 | minor | yes: code path read (both SystemExit branches refuse an empty `shore`) | yes | the driver refuses an empty source only when some ring is ANCHORED (has frozen and free nodes); a rim of islands only passes an empty list, and `rim_constraints` leaves the islands to the source. Test on the round-7 island mesh |
| 4 | minor | yes: `AttributeError` on `Polygon()` | yes | empty land: no coast, no land faces, clearance infinite; test closes the pocket |
| 5 | minor | yes: `Path.iterdir` raising gave the .shp alone | yes | `dataset_files` ends the list with a `<stem>.<sidecars-unlisted>` marker that hashes to None, so the dataset digest is null and the record says why; tests for both paths |
| 6 | minor | yes: code read; both passes wrote `rejected_seed<k>.npz` | yes | `rejected_pass<P>_seed<k>.npz`, named in the attempt's `rejected_file`; `seed_search(..., search_pass=)` sets it. Test through the extracted `seed_search` |
| 7 | nit | yes | yes | the fill log prints `len(PFIX_ALL)`, `len(EGFIX_ALL)` |
| 8 | nit | yes | yes | "selected seed N (search pass P); the gates follow" |

The driver-test harness delimits the seed-search block by the log line, so it follows the new wording. All six hires recipes rebuild accepted and byte-identical to round 7. 956 tests pass.

## Round 9

### Answer (verbatim)

Reviewed `6fbe8e8` read-only. No files changed, jobs submitted, or mesh generation run. I found **three minor defects and two nits**. None is a demonstrated regression introduced by the last commit.

**Verification**

Using the requested Python interpreter:

- Full pytest invocation failed before collection because no temporary directory was writable.
- An adjusted scoped run produced **238 passed, 49 skipped, 11 failures**; all failures required temporary directories.
- With Matplotlib cache setup mocked and filesystem/mesher-dependent tests skipped: **285 passed, 98 skipped**.
- In-memory probes below reproduced all three functional findings.
- The author’s compute-node rebuilds and FVCOM results were not independently reproduced.

**Previous findings**

Statuses concern the reported behavior; related remaining gaps are identified explicitly.

| Finding | Status | Reason |
|---|---|---|
| R1 #1 — Preserve mode | RESOLVED | Geometry-changing rim repair requires `resolve`. |
| R1 #2 — Lagoon entrance | RESOLVED | Connectivity guard preserves the reproduced entrance. |
| R1 #3 — Refused pocket walls | RESOLVED | Restoration uses pocket geometry intersecting remaining water. |
| R1 #4 — Three-node ring collapse | RESOLVED | Removal guards retain a valid ring. |
| R1 #5 — Retreat-only retry | RESOLVED | Array comparison detects changes; `finally` restores unselected state. |
| R1 #6 — Empty land/padding | RESOLVED | Empty-land handling and exterior bounds cover the reported cases. |
| R1 #7 — Invalid continuous-width parameters | RESOLVED | Reported invalid parameters and empty footprints are handled. |
| R1 #8 — Ignored seed argument | RESOLVED | Search iterates `seed_list`. |
| R1 #9 — GPL provenance import | RESOLVED | Versions come from distribution metadata. |
| R1 #10 — Missing provenance inputs/settings | RESOLVED | Inputs/settings are recorded; failed sidecar enumeration now nulls the digest. |
| R1 #11 — Untracked source | RESOLVED | Untracked files and entry-point tracking status are recorded. |
| R1 #12 — Retry reproducibility | RESOLVED | Documentation requires the original seed sequence for pass 2. |
| R1 #13 — Ineffective fixtures | RESOLVED | Relevant prerequisites/counters are asserted; unreachable merge code was removed. |
| R1 #14 — Raster documentation | RESOLVED | Resolution is expressed relative to the finest target. |
| R2 #1 — Resolvable land accepted as water | RESOLVED | Recursive land validation checks qualifying land inside the hole. |
| R2 #2 — All-land raster | RESOLVED | Empty-water early return handles it. |
| R2 #3 — Detached artificial island | RESOLVED | Closure requires adjacent land. |
| R2 #4 — Fine region inside pocket | RESOLVED | Conservative bounds protect the reproduced fine region. |
| R2 #5 — Obsolete sliding curves | PARTIAL | `rim_repair` changes are pinned, but earlier blunting changes are not—finding 1. |
| R2 #6 — Geometry provenance | RESOLVED | Sidecars are enumerated; enumeration failure is explicit. |
| R2 #7 — Falsy experimental values | RESOLVED | Explicit list validation rejects them. |
| R2 #8 — Seed-variable documentation | RESOLVED | Header names `LR_SEEDS`. |
| R3 #1 — Variable-size land guard | RESOLVED | Decisions use conservative bounds. |
| R3 #2 — Repairs bypass validation | RESOLVED | Final-rim and retry land checks remain wired. |
| R3 #3 — Unchanged new corners | RESOLVED | Changed adjacency during `rim_repair` causes pinning. |
| R3 #4 — Pocket sampling grid | RESOLVED | Covering samples and uncertified-floor handling address the reported gaps. |
| R3 #5 — Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6 — Imported-package identity | RESOLVED | Actual imported paths are recorded. |
| R3 #7 — Exhausted iterators | RESOLVED | Both helpers materialize their inputs. |
| R4 #1 — Missed qualifying disc | RESOLVED | Conservative erosion replaces sampled absence decisions. |
| R4 #2 — Incorrect driver size floor | RESOLVED | `H_FLOOR` includes the base ambient minimum. |
| R4 #3 — Validation disabled with repair | RESOLVED | Final land validation is outside the repair switch. |
| R4 #4 — Fine basin behind coarse mouth | RESOLVED | Estimated floors no longer erase the reproduced basin. |
| R4 #5 — Pocket fine-zone exclusion | RESOLVED | Estimated floors no longer override the fine-zone bound. |
| R4 #6 — Non-Git identity | RESOLVED | Path, distribution and source digest are retained. |
| R4 #7 — Invalid sizes/factors | RESOLVED | Reported invalid cases are rejected. |
| R4 #8 — Missing test prerequisites | RESOLVED | Revised fixtures assert prerequisites and exercise guards. |
| R4 #9 — Priority documentation | RESOLVED | Guide states that local refinement ignores priority. |
| R5 #1 — Grid coverage | RESOLVED | Boundary-adjacent samples provide coverage. |
| R5 #2 — Invalid slope premise | RESOLVED | Base and patch fields supply structural bounds. |
| R5 #3 — Pocket area/clearance | RESOLVED | Both predicates use upper bounds. |
| R5 #4 — Nested land collections | RESOLVED | Land validation recursively traverses polygons. |
| R5 #5 — Invalid sampled sizes | RESOLVED | Validation precedes clipping. |
| R5 #6 — Unreadable source hashing | RESOLVED | Failed reads/scans invalidate the digest. |
| R5 #7 — Stale floor comment | RESOLVED | Comment names `H_FLOOR`. |
| R6 #1 — Hidden enumeration failures | RESOLVED | `Path.walk(on_error=...)` records failures. |
| R6 #2 — Exterior disc centres | RESOLVED | Certified maxima bound the search; unsupported bounded fields refuse closure. |
| R6 #3 — Nested island collections | RESOLVED | `island_rings` recursively visits polygons. |
| R6 #4 — Empty filter footprint | RESOLVED | Explicit descriptive rejection. |
| R6 #5 — Ineffective fine-spot test | RESOLVED | Revised fixture and negative control exercise the bound. |
| R6 #6 — Opt-in comments | RESOLVED | Reported comments now describe defaults correctly. |
| R7 #1 — Fully free islands | RESOLVED | Source replacement also permits islands removed entirely by filtering. |
| R7 #2 — Discontinuous exterior sizes | RESOLVED | Both padding regressions pass. |
| R7 #3 — Estimated default floors | RESOLVED | Pocket and slit regressions preserve fine water. |
| R7 #4 — Degenerate ring | RESOLVED | Original reproduction remains valid; ordinary removal still works. |
| R7 #5 — Closed-ring arc selection | RESOLVED | Uneven-density and equal-length cases are handled. |
| R7 #6 — Obsolete driver comment | RESOLVED | Reported comment describes filter → rim → hole. |
| R8 #1 — Snapped slit self-intersection | RESOLVED | Candidate validation and invalid-hole rejection pass the regression. |
| R8 #2 — Equal-length arcs | RESOLVED | Fit and geometric tie-breaking remove the reported start dependence. |
| R8 #3 — Filtered-away island | RESOLVED | Empty shoreline is permitted for islands-only rims. |
| R8 #4 — Empty pocket land | RESOLVED | Empty-land regression passes; mixed collections remain problematic—finding 3. |
| R8 #5 — Failed sidecar enumeration | RESOLVED | Injected `PermissionError` produces the marker and a null digest. |
| R8 #6 — Overwritten rejection diagnostics | RESOLVED | Filenames include search pass and are recorded as `rejected_file`. |
| R8 #7 — Stale fill counts | RESOLVED | Log counts the assembled constraint arrays. |
| R8 #8 — Premature “accepted” log | RESOLVED | Log says “selected seed”; acceptance remains gated. |

**Findings**

1. **Minor — Initial blunting changes can still slide along obsolete curves.**

   **Location:** [420_local_refine.py:781](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:781), [420_local_refine.py:1248](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1248).

   The first `blunt_acute_corners` call changes the rim without updating its sliding curves or recording protected corners. `RIM_PINNED` is subsequently initialized empty. `apply_rim_repair` compares against this already-blunted rim, so an unchanged repair pass cannot recover those earlier changes.

   **Reproduction:** With cyclic rim
   `[(0,0),(0,600),(-600,900),(-600,600)]`, base IDs `[-1,1,2,3]`, and constant `h=30`, initial blunting creates `(0,30)` and `(-21.213203,21.213203)`. Executing the extracted `apply_rim_repair` leaves **`RIM_PINNED=[]`**, with every change counter zero.

   For the adjacent fan:

   ```python
   nodes = [(0,30), (0,40), (-21.2132034356,21.2132034356), (-10,35)]
   elements = [[1,3,0], [0,3,2]]
   ```

   One `improve_patch` round, allowing only node 0 to slide, moves it to **`(0,26.485281)`** on the original curve. Supplying the blunted curve keeps it at **`(0,30)`**. Thus seam improvement can reshape the supposedly finished rim.

   **Fix:** Track changed corners from the first geometry edit onward, preserving their identity through wall construction and remapping; alternatively maintain current sliding curves separately from source-reference curves. Protection must also work with `hires.rim_repair: false`.

2. **Minor — Blunting can reverse an enclosed ring’s land/water role.**

   **Location:** [patch.py:2323](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2323).

   `_ring_is_polygon` checks only the edited ring. A blunting chord can leave another ring entirely outside its former enclosing shell. `hole_polygon` then interprets that former island as a water shell. The result is geometrically valid, so the new validity check does not detect the semantic inversion.

   **Reproduction:** Use these two cyclic rings:

   ```python
   outer = [(0,0), (200,600), (0,900), (-200,600)]
   island = [(-1,10), (1,10), (1,12), (-1,12)]
   ```

   Make the outer points free, the four island points frozen, and use constant `h=30`. Both default `rim_repair` and `operations=("angles",)` report one blunted corner.

   Before repair, the hole is a **179,996 m² Polygon**, excluding `(0,11)`. Afterwards it is a **valid 179,734 m² MultiPolygon containing `(0,11)`**: the island’s **4 m² of land became water**, despite all its frozen points surviving.

   This demonstrates a helper-level topology defect; acceptance by the complete driver was not demonstrated.

   **Fix:** Validate nesting and land/water roles of untouched rings before committing a repair. Refuse a cut that strands an enclosed ring. Polygon validity alone is insufficient.

3. **Minor — Nonempty mixed land collections still crash pocket closure.**

   **Location:** [walls.py:652](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/walls.py:652).

   A nonempty `GeometryCollection` containing polygons and lines also has `boundary=None`. The round-8 guard handles only empty geometry, so `coast=None` reaches the `g.is_empty` expression.

   **Reproduction:**

   ```python
   wall = LineString([(0,0), (0,400), (110,400), (110,290), (20,290)])
   land = GeometryCollection([
       box(-2000,-2000,2000,0),
       LineString([(1000,1000), (1100,1000)]),
   ])
   close_wall_pockets(
       [wall], land, lambda q: np.full(len(q), 100.),
       min_h=60., size_floor=30.,
   )
   ```

   This raises **`AttributeError: 'NoneType' object has no attribute 'is_empty'`**. I also reproduced it using polygon-and-line output directly from `shapely.make_valid` on a polygon with a retraced spike.

   **Fix:** Recursively extract and union polygonal land before computing its boundary and clearance, consistently with the other land helpers. Apply empty-land handling after extraction.

4. **Nit — Wall-component construction is dead code, and its docstring describes the abandoned behavior.**

   **Location:** [420_local_refine.py:1508](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1508).

   `walls_cut_off` constructs `wp` and unions wall components through `wfind`, then never uses the resulting components. Actual selection uses endpoints touching stranded water and, when coordinates are supplied, returns the shortest touching edge. The docstring still describes connected wall pieces.

   **Evidence:** No result from `wp` is consumed after `del wfind`; return selection uses only `we`, `stranded`, and `xy`.

   **Fix:** Remove the unused wall-component construction and document the current edge-selection behavior.

5. **Nit — Source-resolution comments and a helper retain the removed coarse-base policy.**

   **Location:** [patch.py:1303](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1303), [patch.py:1104](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1104).

   The comment says `resolve` follows the source in fine areas while retaining the base shape in coarse areas. `_resolve_stretch` explicitly implements source following over the whole stretch and explains why that earlier policy was removed. `_base_spacing`, associated with the former decision, has no callers in source, tests, or notebooks.

   **Fix:** Describe the delivered, simplified source polyline accurately and remove the unused helper.

## Verdict

VERDICT: FAIL (0 blocker, 0 major, 3 minor, 2 nit)


### Prompt

```markdown
# Review request, round 9: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 28 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, f07c6d9 round 4, c88739e round 5, 806d011 round 6, e88208f round 7, the last commit round 8; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-7: see your round-8 status table.
Round 8 (FAIL 0/0/6/2): all fixed in the last commit (`git show HEAD`;
triage in `docs/coastline-rules-review-20260927.md`).
1. slit closure: snapped chord, then the throat, whichever first leaves
   its ring a simple polygon and the whole rim a valid hole; else
   reported. `hole_polygon` raises on an invalid union.
2. equal-length arcs: by fit, then by midpoint coordinates.
3. driver: an empty source shoreline is refused only when some ring has
   an anchored stretch (`ANCHORED`); islands-only rims pass `[]`.
4. `close_wall_pockets` with empty land.
5. `dataset_files` appends a `<stem>.<sidecars-unlisted>` marker that
   hashes to None when the directory cannot be listed.
6. `rejected_pass<P>_seed<k>.npz`, recorded as `rejected_file`.
7. fill log counts; 8. "selected seed" log line.
Tests: `tests/test_review_r8.py`. All six hires recipes rebuild accepted,
byte-identical to round 7. 956 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | minor | yes: code read -- the first `blunt_acute_corners` changes the rim and `RIM_PINNED` starts empty; the reviewer's rim blunts to (0, 30) and (-21.2, 21.2) | yes | `_rim_changed` (made, moved, or neighbours changed) is shared by the first blunting and `apply_rim_repair`; the blunting's points are held by position (`PINNED_XY`) until the walls are rooted, then `RIM_PINNED` starts from them, whether or not rim repair runs. Test through the extracted driver functions |
| 2 | minor | yes: reviewer's outer + island; after blunting the hole contains (0, 11) | yes | `_rim_depths` (sorted nesting depth of every ring); `_ring_is_polygon(..., depths)` refuses an edit that changes it, for short edges, slits and both blunting paths. Tests: both call paths of the reproduction, and a plain blunt still cuts |
| 3 | minor | yes: reproduction raises `AttributeError` | yes | `close_wall_pockets` takes the polygons of the land (`_polygons`, recursive) before anything else; test against plain land |
| 4 | nit | yes: `wp`/`wfind` never used | yes | removed; docstring describes the touching / shortest edge selection |
| 5 | nit | yes: `_base_spacing` has no callers; the comment described the removed per-vertex policy | yes | removed; comment describes the delivered polyline |

All six hires recipes rebuild accepted with 0 introduced violations. Five are byte-identical to round 8; Funabashi moved by at most 8.79 m (Hausdorff), same seed 3 and dt 1.48 s -- the first blunting's point is now held in the seam repair. 960 tests pass.

## Round 10

### Answer (verbatim)

Found **three minor defects and one nit**. Two are gaps in the round-9 fixes; the wall-root collapse is pre-existing.

No files were modified. No mesh generation or batch jobs ran.

1. **Minor — Tip retreat still changes a ring’s land/water role.**

   **Location:** [patch.py:1731](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1731)

   Round 9 checks nesting for short-edge removal, slit closure and blunting, but `retreat()` commits its new coordinates without that check.

   **Reproduction:** Use these two cyclic rings:
   ```python
   outer = [(0,0), (190,0), (190,297), (200,297), (210,297),
            (210,0), (600,0), (600,300), (0,300)]
   lake = [(191,280), (193,280), (193,282), (191,282)]
   ```
   Concatenate their points and edges; freeze every point except outer point 2. Call `rim_repair` with constant `h=30`, `operations=("slits",)`, `rounds=1`, `focus=[[190,297]]`, and `focus_factor=.001`.

   The actual run reports one retreat. Nesting changes from **`(0,0)` to `(0,1)`**: `(192,281)` is water before repair and land afterwards. `hole_polygon` accepts both geometries as valid. The 4 m² water ring becomes an island.

   **Fix:** Before committing a retreat, validate the candidate ring and complete hole, and require unchanged nesting. Apply the guard consistently to every geometry-edit path.

2. **Minor — Wall-root folding loses the initial blunting’s pins.**

   **Locations:** [420_local_refine.py:943](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:943), [420_local_refine.py:1279](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1279)

   `PINNED_XY` stores coordinates before wall construction. `_blunt` subsequently merges and moves those points, but does not transfer their pin to the survivor. The coordinate lookup therefore misses the edited corner.

   **Reproduction:** Initial blunting of cyclic rim
   `[(0,0),(100,600),(0,900),(-100,600)]`, base IDs `[-1,1,2,3]`, and `h=30` creates pinned endpoints at approximately **`(±4.931970,29.591818)`**. Executing the driver’s extracted `_blunt` on either endpoint folds the cap into **`(0,29.591818)`**. That survivor is absent from `PINNED_XY`.

   If its wall is subsequently dropped and rim repair is disabled, it becomes eligible for sliding. On the explicit two-face fan:
   ```python
   nodes = [(100,600), (0,900), (-100,600), (0,29.59181771496431)]
   elements = [[0,1,3], [1,2,3]]
   ```
   one `improve_patch` round permitting only node 3 to slide moves it to **`(-42.879204,257.275223)`**. This reproduces with the original coastline stretch and a clipped wall curve `[(0,43),(0,200)]` supplied. Explicitly pinning the survivor produces zero moves.

   **Fix:** Transfer pin identity through root moves, folds and compaction, or use `_rim_changed` around wall construction and pin its surviving edits. Preserve those pins even when the associated wall is withdrawn.

3. **Minor — Wall-root folding can collapse a complete triangular ring.**

   **Location:** [420_local_refine.py:934](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:934)

   The driver’s `_blunt` has none of the minimum-ring-size or polygon-validity guards added to the package repair functions.

   **Reproduction:** Execute the actual AST-extracted `_blunt` with:
   ```python
   rim_xy = np.array([(0,0), (10,0), (0,10)], dtype=float)
   rim_eg = np.array([(0,1), (1,2), (2,0)])
   rim_base = np.full(3, -1)
   _segs, _folded = [], []
   # _blunt(0, 30.0)
   ```
   Starting from a valid **50 m²** polygon, it folds `(1,0)` and then `(2,0)`, leaving **no edges**. The subsequent `hole_polygon` raises `ValueError: the rim segments do not close a polygon`. The later `rim_repair` guards cannot rescue this earlier mutation.

   **Fix:** Validate each proposed fold before committing it. Reject folds that leave fewer than three distinct vertices, duplicate edges, invalid geometry or changed nesting.

4. **Nit — The wall-withdrawal call-site comment still describes removal of whole components.**

   **Location:** [420_local_refine.py:1674](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:1674)

   The comment says every touching connected wall piece is withdrawn. The call supplies `xy`, so `walls_cut_off` returns only the shortest touching edge. Round 9 corrected the helper’s docstring but left this contradictory explanation.

   **Fix:** Describe withdrawal of one shortest touching edge followed by another split/connectivity check.

**Status of every previous finding**

“RESOLVED” below refers to the reported defect and its reproduction; related uncovered paths are identified separately above.

| Finding | Status | Reason |
|---|---|---|
| R1 #1 — Preserve mode | RESOLVED | Rim repair requires `resolve`. |
| R1 #2 — Lagoon entrance | RESOLVED | Connectivity guard preserves the reported entrance. |
| R1 #3 — Refused pocket walls | RESOLVED | Restoration tests pocket geometry against remaining water. |
| R1 #4 — Three-node ring collapse | RESOLVED | Package repair guards the reported operation; driver folding has finding 3. |
| R1 #5 — Retreat-only retry | RESOLVED | Array comparison detects changes; `finally` restores unselected state. |
| R1 #6 — Empty land/padding | RESOLVED | Empty-land handling and exterior padding cover reported cases. |
| R1 #7 — Invalid continuous-width parameters | RESOLVED | Reported parameters and empty footprints are handled. |
| R1 #8 — Ignored seed argument | RESOLVED | Search iterates `seed_list`. |
| R1 #9 — GPL provenance import | RESOLVED | Versions use distribution metadata. |
| R1 #10 — Missing provenance inputs/settings | RESOLVED | Bathymetry, geometry inputs and relevant settings are recorded. |
| R1 #11 — Untracked source | RESOLVED | Untracked files and entry-point tracking status are recorded. |
| R1 #12 — Retry reproducibility | RESOLVED | Documentation requires the original seed sequence for pass 2. |
| R1 #13 — Ineffective fixtures | RESOLVED | Revised fixtures assert prerequisites; unreachable merge code was removed. |
| R1 #14 — Raster documentation | RESOLVED | Spacing is expressed relative to the finest target. |
| R2 #1 — Resolvable land accepted as water | RESOLVED | Land validation rejects qualifying land inside the hole. |
| R2 #2 — All-land raster | RESOLVED | Empty-water early return handles it. |
| R2 #3 — Detached artificial island | RESOLVED | Closure requires adjacent land. |
| R2 #4 — Fine region inside pocket | RESOLVED | Conservative bounds protect the reported fine region. |
| R2 #5 — Obsolete sliding curves | PARTIAL | Pinning works across rim repair, but wall folding loses pins: finding 2. |
| R2 #6 — Geometry provenance | RESOLVED | Sidecars are enumerated; enumeration failure invalidates the digest. |
| R2 #7 — Falsy experimental values | RESOLVED | Explicit list validation rejects them. |
| R2 #8 — Seed-variable documentation | RESOLVED | Header names `LR_SEEDS`. |
| R3 #1 — Variable-size land guard | RESOLVED | Decisions use conservative bounds. |
| R3 #2 — Repairs bypass validation | RESOLVED | Final-rim and retry land checks remain wired. |
| R3 #3 — Unchanged new corners | RESOLVED | Changed adjacency during rim repair causes pinning. |
| R3 #4 — Pocket sampling grid | RESOLVED | Covering samples and uncertified-floor handling address reported gaps. |
| R3 #5 — Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6 — Imported-package identity | RESOLVED | Actual imported paths are recorded. |
| R3 #7 — Exhausted iterators | RESOLVED | Both helpers materialize their inputs. |
| R4 #1 — Missed qualifying disc | RESOLVED | Conservative erosion replaces sampled absence decisions. |
| R4 #2 — Incorrect driver size floor | RESOLVED | `H_FLOOR` includes the base ambient minimum. |
| R4 #3 — Validation disabled with repair | RESOLVED | Final land validation is outside the repair switch. |
| R4 #4 — Fine basin behind coarse mouth | RESOLVED | Conservative bounds preserve the reproduced basin. |
| R4 #5 — Pocket fine-zone exclusion | RESOLVED | Conservative bounds protect the fine zone. |
| R4 #6 — Non-Git identity | RESOLVED | Path, distribution and source digest are recorded. |
| R4 #7 — Invalid sizes/factors | RESOLVED | Reported invalid inputs are rejected. |
| R4 #8 — Missing test prerequisites | RESOLVED | Revised fixtures exercise the intended guards. |
| R4 #9 — Priority documentation | RESOLVED | Guide states that local refinement ignores priority. |
| R5 #1 — Grid coverage | RESOLVED | Boundary-adjacent samples provide coverage. |
| R5 #2 — Invalid slope premise | RESOLVED | Base and patch fields supply structural bounds. |
| R5 #3 — Pocket area/clearance | RESOLVED | Both predicates use upper bounds. |
| R5 #4 — Nested land collections | RESOLVED | Land validation recursively visits polygons. |
| R5 #5 — Invalid sampled sizes | RESOLVED | Validation precedes clipping. |
| R5 #6 — Unreadable source hashing | RESOLVED | Read and scan failures invalidate the digest. |
| R5 #7 — Stale floor comment | RESOLVED | Comment names `H_FLOOR`. |
| R6 #1 — Hidden enumeration failures | RESOLVED | `Path.walk(on_error=...)` records failures. |
| R6 #2 — Exterior disc centres | RESOLVED | Certified maxima bound the search; unsupported fields refuse closure. |
| R6 #3 — Nested island collections | RESOLVED | `island_rings` recursively visits polygons. |
| R6 #4 — Empty filter footprint | RESOLVED | Explicit descriptive rejection. |
| R6 #5 — Ineffective fine-spot test | RESOLVED | Revised fixture and negative control exercise the bound. |
| R6 #6 — Opt-in comments | RESOLVED | Reported comments describe defaults correctly. |
| R7 #1 — Fully free islands | RESOLVED | Source replacement also permits filtered-away islands. |
| R7 #2 — Discontinuous exterior sizes | RESOLVED | Both exterior-padding regressions pass. |
| R7 #3 — Estimated default floors | RESOLVED | Estimated floors no longer override conservative lower bounds. |
| R7 #4 — Degenerate ring | RESOLVED | Reported removal/blunting paths retain valid polygons. |
| R7 #5 — Closed-ring arc selection | RESOLVED | Uneven-density and ring-start regressions pass. |
| R7 #6 — Obsolete driver comment | RESOLVED | Comment describes filter → rim → hole. |
| R8 #1 — Snapped slit self-intersection | RESOLVED | Candidate validation rejects the reported crossing. |
| R8 #2 — Equal-length arcs | RESOLVED | Fit and geometric tie-breaking remove start dependence. |
| R8 #3 — Filtered-away island | RESOLVED | Empty shoreline is permitted for islands-only rims. |
| R8 #4 — Empty pocket land | RESOLVED | Empty land and mixed collections are handled. |
| R8 #5 — Failed sidecar enumeration | RESOLVED | Failure marker forces a null digest. |
| R8 #6 — Overwritten rejection diagnostics | RESOLVED | Rejection filenames include search pass. |
| R8 #7 — Stale fill counts | RESOLVED | Log counts assembled constraints. |
| R8 #8 — Premature “accepted” log | RESOLVED | Log says “selected seed”; acceptance remains gated. |
| R9 #1 — Initial blunting pins | PARTIAL | Initial pins are created, but folding loses their identity: finding 2. |
| R9 #2 — Changed nesting | PARTIAL | Guards cover reported blunting paths but omit retreat: finding 1. |
| R9 #3 — Mixed pocket land | RESOLVED | Recursive polygon extraction handles the collection. |
| R9 #4 — Dead wall-component code/docstring | RESOLVED | Dead construction removed and docstring corrected; call-site nit remains. |
| R9 #5 — Obsolete resolution helper/comment | RESOLVED | Helper removed and comment corrected. |

**Verification**

Used the requested `oceanmesh-bench/bin/python`.

- Ordinary pytest failed before collection because temporary directories are unwritable.
- `--capture=sys -p no:cacheprovider` then encountered Matplotlib’s writable-cache requirement.
- An in-memory harness allowed reading existing Matplotlib caches and skipped filesystem-dependent and mesher tests: **756 passed, 200 skipped, 4 failed**. All four failures were OCSMesh imports requiring a writable temporary directory in `tests/test_mesh_clean.py`.
- Findings 1–3 were reproduced separately in memory, using package functions and AST-extracted driver functions.
- Final `git status --short` was empty.

No additional regression introduced specifically by the last commit was substantiated beyond its incomplete guards and pin propagation. Production mesh QA counts and implied time steps were not rerun.

## Verdict

VERDICT: FAIL (0 blocker, 0 major, 3 minor, 1 nit)


### Prompt

```markdown
# Review request, round 10: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 29 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, f07c6d9 round 4, c88739e round 5, 806d011 round 6, e88208f round 7, 6fbe8e8 round 8, the last commit round 9; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-8: see your round-9 status table.
Round 9 (FAIL 0/0/3/2): all fixed in the last commit (`git show HEAD`;
triage in `docs/coastline-rules-review-20260927.md`).
1. `_rim_changed` shared by the first `blunt_acute_corners` and
   `apply_rim_repair`; the blunting's points are held by position
   (`PINNED_XY`) until walls are rooted, then `RIM_PINNED` starts from
   them, with or without rim repair.
2. `_rim_depths`: short-edge, slit and both blunting edits are refused
   when the sorted nesting depths of the rim's rings change.
3. `close_wall_pockets` takes `_polygons` of the land first.
4. dead code in `walls_cut_off` removed, docstring updated.
5. `_base_spacing` removed; resolve comment corrected.
Tests: `tests/test_review_r9.py`. All six hires recipes rebuild accepted
(Funabashi moves 8.8 m from #1, same seed; the rest byte-identical).
960 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | minor | yes: code read -- `retreat` committed without the nesting check; reproduction turns a 4 m2 lake into land | yes | `retreat` builds the moved rim and refuses it unless its ring is a simple polygon and `_rim_depths` is unchanged. Test: the reproduction keeps the lake water and the depths |
| 2 | minor | yes: code read -- `_blunt` and the root move change coordinates that `PINNED_XY` holds by position | yes | `_pin_follows`: a pinned point that folds or moves passes its pin to where it lands. Test on the extracted `_blunt` |
| 3 | minor | yes: extracted `_blunt` folds a triangle to no edges | yes | each fold is built first and made only if the root's ring stays a simple polygon with unchanged nesting (`_ring_is_polygon`, `_rim_depths`, now imported by the driver). Tests: the triangle survives; a short edge on a square still folds |
| 4 | nit | yes | yes | call-site comment: the shortest touching edge, one at a time |

All six hires recipes rebuild accepted and byte-identical to round 9. 963 tests pass.

## Round 11

### Answer (verbatim)

Four minor defects remain. No regression unique to the last commit was substantiated; these are incomplete geometry and pinning safeguards.

1. **Minor — Wall-root movement still bypasses polygon and nesting validation.**

   **Location:** [420_local_refine.py:915](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:915)

   `_root` moves a vertex whenever its bend is below 15°. That permits a non-collinear move which changes the coastline. Unlike retreat and folding, this assignment has no candidate-ring guard.

   **Reproduction:** Execute the AST-extracted `_root` with cyclic rings:
   ```python
   outer = [(100,-5), (300,0), (300,300), (-100,300), (-100,0)]
   island = [(99,-4.95), (101,-4.95), (101,-4.9), (99,-4.9)]
   ```
   Set only outer vertex 0 free, use constant `h=30`, and call `_root([110,-4.75])`.

   Observed:
   - Vertex 0 moves to `(110,-4.75)`.
   - Nesting changes from `(0,1)` to `(0,0)`.
   - `hole_polygon` succeeds, but `(100,-4.925)` changes from land to water.

   **Fix:** Build the proposed coordinates first. Validate the edited polygon, inter-ring separation and each ring’s nesting before moving the vertex or transferring its pin.

2. **Minor — The new fold guard permits crossings between different rings.**

   **Locations:** [420_local_refine.py:953](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:953), [patch.py:2139](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2139)

   `_ring_is_polygon` checks the edited ring’s validity and the nesting-depth summary. It does not check whether that ring crosses another ring. Testing one vertex for containment can leave the summary unchanged despite a crossing.

   **Reproduction:** With all vertices free, execute extracted `_blunt(0,30)` on:
   ```python
   outer = [(0,0), (10,0), (300,0), (300,300), (0,300)]
   island = [(6,10), (6,20), (1,20), (1,10)]
   ```
   The fold `(1,0)` is committed, moving the root to `(5,0)`. Depths remain `(0,1)`, but the new coast crosses the island. Subsequent `hole_polygon` raises:
   ```text
   the assembled hole boundary is 1202.54 m against 1225.04 m
   of rim constraints; the rings touch or overlap
   ```

   **Fix:** Reject candidates whose ring boundaries intersect other rings. Validate the complete candidate rim before committing the fold, retaining the original geometry on refusal.

3. **Minor — Sorted nesting depths allow individual rings to exchange land/water roles.**

   **Location:** [patch.py:2108](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2108)

   Comparing sorted depths preserves only their counts. It does not establish the documented invariant that *every ring keeps its role*.

   **Reproduction:** With all vertices free, execute extracted `_blunt(0,30)` on:
   ```python
   outer = [(0,0), (10,10), (300,0), (300,300), (0,300)]
   island = [(1,20), (2,20), (2,21), (1,21)]
   lake = [(8,6), (9,6), (9,7), (8,7)]
   ```
   The fold is accepted. Before and after, `_rim_depths` returns `(0,0,1)`, and `hole_polygon` remains valid. Nevertheless:
   - `(1.5,20.5)` changes from land to water.
   - `(8.5,6.5)` changes from water to land.

   **Fix:** Match rings across an edit using surviving vertex identities and preserve their individual containment relationships. An unordered depth histogram is insufficient.

4. **Minor — Newly folded wall-root corners remain unpinned when neither predecessor was pinned.**

   **Location:** [420_local_refine.py:960](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/420_local_refine.py:960)

   `_pin_follows` successfully transfers existing pins, but creates none for previously unpinned geometry changed by wall construction. Subsequent `apply_rim_repair` compares against the already-folded rim, so an unchanged repair pass also creates no pin.

   **Reproduction:** Start with:
   ```python
   points = [(0,0), (10,0), (300,0), (300,300), (0,100)]
   base_ids = [-1, -1, 2, 3, 4]
   ```
   Use cyclic edges and constant `h=30`.

   Initial `blunt_acute_corners` reports zero edits. Extracted `_blunt(0,30)` folds `(1,0)` into `(5,0)`, leaving `PINNED_XY` empty. After compaction, ordinary `rim_repair` reports zero operations and `_rim_changed` returns no pins.

   On the compacted four-node mesh with faces `[[0,1,2],[0,2,3]]`, one `improve_patch` round allowing only that root to slide moves it from **`(5,0)` to `(10,0)`**. This uses the original coast curve and a clipped wall curve `[(5,13),(5,60)]`. Explicitly pinning the root produces zero moves.

   **Fix:** Record all surviving geometry and adjacency changes during wall construction, using `_rim_changed` around that stage, and include them in `RIM_PINNED`. Keep these pins when walls are subsequently withdrawn.

**Status of every previous finding**

“RESOLVED” refers to the original defect and reproduction. Broader remaining gaps are identified explicitly.

| Finding | Status | Reason |
|---|---|---|
| R1 #1 — Preserve mode | RESOLVED | Rim repair requires `resolve`. |
| R1 #2 — Lagoon entrance | RESOLVED | Connectivity guard protects the reported entrance. |
| R1 #3 — Refused-pocket walls | RESOLVED | Restoration checks pocket geometry against remaining water. |
| R1 #4 — Three-node collapse | RESOLVED | Reported removal and folding cases retain a polygon. |
| R1 #5 — Retreat-only retry | RESOLVED | Array comparison detects changes; unselected state is restored. |
| R1 #6 — Empty land/padding | RESOLVED | Empty-land handling and exterior padding cover reported cases. |
| R1 #7 — Invalid continuous-width parameters | RESOLVED | Reported invalid parameters and empty footprints are handled. |
| R1 #8 — Ignored seed argument | RESOLVED | Search iterates `seed_list`. |
| R1 #9 — GPL provenance import | RESOLVED | Versions come from distribution metadata. |
| R1 #10 — Missing provenance inputs/settings | RESOLVED | Bathymetry, region files and relevant settings are recorded. |
| R1 #11 — Untracked source | RESOLVED | Git status includes untracked files; entry-point tracking is recorded. |
| R1 #12 — Retry reproducibility | RESOLVED | Documentation requires the original seed sequence for pass 2. |
| R1 #13 — Ineffective fixtures | RESOLVED | Revised fixtures exercise prerequisites; unreachable merge code was removed. |
| R1 #14 — Raster documentation | RESOLVED | Resolution is expressed relative to the finest target. |
| R2 #1 — Resolvable land accepted as water | RESOLVED | Land validation rejects qualifying land inside the hole. |
| R2 #2 — All-land raster | RESOLVED | Empty-water early return handles it. |
| R2 #3 — Detached artificial island | RESOLVED | Closure requires adjacent land. |
| R2 #4 — Fine region inside pocket | RESOLVED | Conservative bounds protect the reported fine region. |
| R2 #5 — Obsolete sliding curves | PARTIAL | Existing pins survive, but newly folded roots can remain unpinned: finding 4. |
| R2 #6 — Geometry provenance | RESOLVED | Sidecars are enumerated; enumeration failure invalidates the digest. |
| R2 #7 — Falsy experimental values | RESOLVED | Explicit list validation rejects them. |
| R2 #8 — Seed-variable documentation | RESOLVED | Header names `LR_SEEDS`. |
| R3 #1 — Variable-size land guard | RESOLVED | Decisions use conservative bounds. |
| R3 #2 — Repairs bypass validation | RESOLVED | Final-rim and retry land checks remain wired. |
| R3 #3 — Unchanged new corners | RESOLVED | Changed adjacency during rim repair causes pinning. |
| R3 #4 — Pocket sampling grid | RESOLVED | Covering samples and uncertified-floor handling address reported gaps. |
| R3 #5 — Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6 — Imported-package identity | RESOLVED | Actual imported paths are recorded. |
| R3 #7 — Exhausted iterators | RESOLVED | Both helpers materialize their inputs. |
| R4 #1 — Missed qualifying disc | RESOLVED | Conservative erosion replaces sampled absence decisions. |
| R4 #2 — Incorrect driver size floor | RESOLVED | `H_FLOOR` includes the base ambient minimum. |
| R4 #3 — Validation disabled with repair | RESOLVED | Final land validation is outside the repair switch. |
| R4 #4 — Fine basin behind coarse mouth | RESOLVED | Conservative bounds protect the reproduced basin. |
| R4 #5 — Pocket fine-zone exclusion | RESOLVED | Conservative bounds protect the fine zone. |
| R4 #6 — Non-Git identity | RESOLVED | Path, distribution and source digest are recorded. |
| R4 #7 — Invalid sizes/factors | RESOLVED | Reported invalid inputs are rejected. |
| R4 #8 — Missing test prerequisites | RESOLVED | Revised fixtures exercise the intended guards. |
| R4 #9 — Priority documentation | RESOLVED | Guide states that local refinement ignores priority. |
| R5 #1 — Grid coverage | RESOLVED | Boundary-adjacent samples provide coverage. |
| R5 #2 — Invalid slope premise | RESOLVED | Base and patch fields supply structural bounds. |
| R5 #3 — Pocket area/clearance | RESOLVED | Both predicates use upper bounds. |
| R5 #4 — Nested land collections | RESOLVED | Land validation recursively visits polygons. |
| R5 #5 — Invalid sampled sizes | RESOLVED | Validation precedes clipping. |
| R5 #6 — Unreadable source hashing | RESOLVED | Read and scan failures invalidate the digest. |
| R5 #7 — Stale floor comment | RESOLVED | Comment names `H_FLOOR`. |
| R6 #1 — Hidden enumeration failures | RESOLVED | `Path.walk(on_error=...)` records failures. |
| R6 #2 — Exterior disc centres | RESOLVED | Padding accounts for exterior centres; unsupported bounds refuse closure. |
| R6 #3 — Nested island collections | RESOLVED | `island_rings` recursively visits polygons. |
| R6 #4 — Empty filter footprint | RESOLVED | Explicit descriptive rejection. |
| R6 #5 — Ineffective fine-spot test | RESOLVED | Revised fixture and negative control exercise the bound. |
| R6 #6 — Opt-in comments | RESOLVED | Reported comments describe defaults correctly. |
| R7 #1 — Fully free islands | RESOLVED | Source replacement permits filtered-away islands. |
| R7 #2 — Discontinuous exterior sizes | RESOLVED | Exterior-padding regressions pass. |
| R7 #3 — Estimated default floors | RESOLVED | Estimated floors no longer override conservative lower bounds. |
| R7 #4 — Degenerate ring | RESOLVED | Reported removal/blunting paths retain valid polygons. |
| R7 #5 — Closed-ring arc selection | RESOLVED | Uneven-density and ring-start regressions pass. |
| R7 #6 — Obsolete driver comment | RESOLVED | Comment describes filter → rim → hole. |
| R8 #1 — Snapped slit self-intersection | RESOLVED | Candidate validation rejects the reported crossing. |
| R8 #2 — Equal-length arcs | RESOLVED | Fit and geometric tie-breaking remove reported start dependence. |
| R8 #3 — Filtered-away island | RESOLVED | Empty shoreline is permitted for islands-only rims. |
| R8 #4 — Empty pocket land | RESOLVED | Empty land and mixed collections are handled. |
| R8 #5 — Failed sidecar enumeration | RESOLVED | Failure marker forces a null digest. |
| R8 #6 — Overwritten rejection diagnostics | RESOLVED | Rejection filenames include search pass. |
| R8 #7 — Stale fill counts | RESOLVED | Log counts assembled constraints. |
| R8 #8 — Premature “accepted” log | RESOLVED | Log says “selected seed”; acceptance remains gated. |
| R9 #1 — Initial blunting pins | RESOLVED | Initial pins are created and now survive root movement/folding. |
| R9 #2 — Changed nesting | PARTIAL | Reported cases are guarded, but root movement and depth-histogram comparison leave findings 1 and 3. |
| R9 #3 — Mixed pocket land | RESOLVED | Recursive polygon extraction handles the collection. |
| R9 #4 — Dead wall-component code/docstring | RESOLVED | Dead code removed; helper and call-site descriptions agree. |
| R9 #5 — Obsolete resolution helper/comment | RESOLVED | Helper removed and comment corrected. |
| R10 #1 — Retreat changes nesting | RESOLVED | The reported retreat is refused; its regression passes. |
| R10 #2 — Existing pins lost during rooting | RESOLVED | `_pin_follows` transfers existing pins; finding 4 concerns previously unpinned points. |
| R10 #3 — Fold collapses triangle | RESOLVED | Triangle-collapse regression passes; inter-ring validation remains incomplete as reported above. |
| R10 #4 — Withdrawal comment | RESOLVED | Comment describes withdrawing one shortest edge at a time. |

**Verification**

Used `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`.

- Ordinary `-m pytest -q tests/ -p no:cacheprovider` failed before collection because no writable temporary directory exists.
- A second full-suite attempt using system capture was interrupted without a test result.
- `tests/test_review_r10.py --capture=sys --noconftest -p no:cacheprovider`: **3 passed**.
- An in-memory pytest harness reading the existing Matplotlib cache ran the scoped package, driver and review tests: **290 passed, 75 skipped**. Filesystem-dependent tests were skipped.
- All four findings were reproduced in memory using package functions and AST-extracted driver functions.
- Final `git status --porcelain` was empty. No files changed, batch jobs submitted or mesh generation run. Production QA, node/element counts and implied time steps were not rerun.

## Verdict

VERDICT: FAIL (0 blocker, 0 major, 4 minor, 0 nit)


### Prompt

```markdown
# Review request, round 11: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 30 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, f07c6d9 round 4, c88739e round 5, 806d011 round 6, e88208f round 7, 6fbe8e8 round 8, 3903361 round 9, the last commit round 10; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-9: see your round-10 status table.
Round 10 (FAIL 0/0/3/1): all fixed in the last commit (`git show HEAD`;
triage in `docs/coastline-rules-review-20260927.md`).
1. `retreat` refuses a move whose ring is not a simple polygon or that
   changes `_rim_depths`.
2. driver `_pin_follows`: a pinned point that folds into a root or moves
   along the coast passes its pin to where it lands.
3. driver `_blunt` builds each fold first and makes it only if
   `_ring_is_polygon(..., _rim_depths(...))` holds.
4. withdrawal comment.
Every rim edit path is now: short-edge removal, slit closure, retreat,
both blunting paths, and the driver's root fold/move -- please check that
none is left without the ring/nesting guard.
Tests: `tests/test_review_r10.py`. All six hires recipes rebuild
accepted, byte-identical to round 9. 963 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | minor | yes: code read -- `_root` moved a vertex with no guard | yes | the move is built first and made only if `_rim_edit_ok`; test on the extracted `_root` |
| 2 | minor | yes: code read -- the guard checked one ring and a depth summary, not crossings | yes | `_rim_edit_ok(old, new)` replaces `_ring_is_polygon` / `_rim_depths` on every edit path (short edges, slits, retreat, both blunting paths, driver fold and root move): every ring a valid polygon of some area, no two rings intersecting |
| 3 | minor | yes: sorted depths are a histogram, and two rings can trade | yes | rings are matched to their old selves by shared vertices and each must sit inside the same (matched) rings as before. Tests: the island/lake swap, and a unit test of the guard |
| 4 | minor | yes: code read -- `_pin_follows` only moved existing pins | yes | `_pin`: a fold or a root move pins the point it makes and the rim points that now neighbour it (the corners wall construction made, as `_rim_changed` does for rim repair). Test on the extracted `_blunt` |

All six hires recipes rebuild accepted with 0 introduced violations; five byte-identical to round 10. Kimitsu: same boundary (Hausdorff 0.00), but with the wall-root corners pinned (#4) seed 0 introduced 3 violations (c1 at element 6353, two c4) and seed 1 is accepted, dt 1.84 s against 2.02 s. 968 tests pass.

## Round 12

### Answer (verbatim)

1. **Minor — Source islands lose their interior water rings.**

   **Location:** [patch.py:2078](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2078)

   `island_rings` resamples only `g.exterior`. Interior rings disappear when the driver adds the returned coastline to the hole. The resulting geometry remains valid, and the land-as-water check cannot detect water converted to land.

   **Reproduction:** Use constant `h=30`, water `box(0,0,1000,1000)`, and:
   ```python
   land = Polygon(
       box(100,100,900,900).exterior.coords,
       [box(300,300,700,700).exterior.coords],
   )
   ```
   Run `filter_shoreline_local` with `h0=30`, `land_width_factor=0.5`, `land_width_max_band=1`, and `continuous_width=True`; then run `island_rings` and assemble its returned rings as the driver does.

   Observed:
   - The filter preserves the interior ring, with zero land or water loss.
   - Island conversion returns one exterior ring and no refusal.
   - Expected water area: **520,000 m²**; assembled water area: **360,000 m²**.
   - `(500,500)` becomes land.
   - `land_an_element_fits(after.intersection(filtered), h, 30)` returns `[]`.

   **Fix:** Preserve interior rings and their nesting when constructing island constraints. If disconnected interior water is unsupported, reject or explicitly report that case instead of silently filling it.

2. **Minor — Refused edges remain permanently excluded after their surrounding geometry changes.**

   **Location:** [patch.py:1621](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1621)

   Each round initializes `done_refused` from all earlier refusals. An edge’s endpoints can remain unchanged while blunting changes its neighbors and makes removal permissible. That edge is never reconsidered, even with additional rounds.

   **Reproduction:** Use cyclic edges, all-free base IDs, `water=Polygon(p)`, constant `h=30`, and:
   ```python
   p = np.array([
       [36.1,23.5], [24.0,46.8], [5.1,76.2], [-7.4,19.5],
       [-9.2,23.5], [-74.3,-12.9], [4.9,-35.2], [58.2,-5.4],
   ])
   ```
   The input polygon is valid. Both `rounds=2` and `rounds=5` remove one point, blunt two corners, and leave the **4.38634 m** edge between `(-7.4,19.5)` and `(-9.2,23.5)`.

   An in-memory control changing only the initialization to `done_refused = set()` removes two points and leaves **zero edges below 15 m**, with either round count. No geometry guard was disabled.

   **Fix:** Invalidate refusal decisions when geometry changes—at least between rounds, preferably after edits affecting their prerequisites. Keep refusal history separately from the current eligibility cache.

3. **Minor — “Short edges left” reports edges that no longer exist.**

   **Location:** [patch.py:1894](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:1894)

   The final count and locations come from accumulated refusals, rather than the returned rim. Later blunting or slit operations can eliminate those edges.

   **Reproduction:**
   ```python
   p = np.array([(0.,0.), (10.,0.), (100.,20.)])
   e = np.array([(0,1), (1,2), (2,0)])
   result = rim_repair(
       p, e, np.full(3, -1), Polygon(p),
       lambda q: np.full(len(q), 30.),
   )
   ```
   The report says **three short edges remain**, including the original 10 m edge at `(5,0)`. Neither endpoint of that edge survives. The returned edge lengths are approximately **0.638, 42.195, 1.549, and 41.980 m**: only **two** are short.

   **Fix:** Compute remaining-edge diagnostics from the final coordinates, connectivity, and size field. Expose historical refusals under a separately named field if useful.

**Status of every previous finding**

“RESOLVED” below refers to the original defect and reproduction. The new findings above identify separate remaining gaps.

| Finding | Status | Reason |
|---|---|---|
| R1 #1 — Preserve mode | RESOLVED | Rim repair and retry require `resolve`. |
| R1 #2 — Lagoon entrance | RESOLVED | Connectivity guard protects the reported entrance. |
| R1 #3 — Refused-pocket walls | RESOLVED | Restoration checks pocket geometry against remaining water. |
| R1 #4 — Three-node collapse | RESOLVED | Removal and blunting reject degenerate candidate rings. |
| R1 #5 — Retreat-only retry | RESOLVED | Array comparison detects changes; unselected state is restored. |
| R1 #6 — Empty land/padding | RESOLVED | Empty-land handling and exterior padding cover the reported cases. |
| R1 #7 — Invalid continuous-width parameters | RESOLVED | Reported invalid parameters and empty footprints are handled. |
| R1 #8 — Ignored seed argument | RESOLVED | Search iterates `seed_list`. |
| R1 #9 — GPL provenance import | RESOLVED | Versions come from distribution metadata. |
| R1 #10 — Missing provenance inputs/settings | RESOLVED | Bathymetry, region files and relevant settings are recorded. |
| R1 #11 — Untracked source | RESOLVED | Git status includes untracked files; entry-point tracking is recorded. |
| R1 #12 — Retry reproducibility | RESOLVED | Documentation requires the original seed sequence for pass 2. |
| R1 #13 — Ineffective fixtures | RESOLVED | Revised fixtures exercise the guards; unreachable merge code was removed. |
| R1 #14 — Raster documentation | RESOLVED | Resolution is expressed relative to the finest target. |
| R2 #1 — Resolvable land accepted as water | RESOLVED | Land validation rejects qualifying land inside the hole. |
| R2 #2 — All-land raster | RESOLVED | Empty-water early return handles it. |
| R2 #3 — Detached artificial island | RESOLVED | Closure requires adjacent land. |
| R2 #4 — Fine region inside pocket | RESOLVED | Conservative bounds protect the reported fine region. |
| R2 #5 — Obsolete sliding curves | RESOLVED | Edited corners are pinned, including newly folded or moved roots. |
| R2 #6 — Geometry provenance | RESOLVED | Sidecars are enumerated; enumeration failure invalidates the digest. |
| R2 #7 — Falsy experimental values | RESOLVED | Explicit list validation rejects them. |
| R2 #8 — Seed-variable documentation | RESOLVED | Header names `LR_SEEDS`. |
| R3 #1 — Variable-size land guard | RESOLVED | Decisions use conservative bounds. |
| R3 #2 — Repairs bypass validation | RESOLVED | Final-rim and retry land checks remain wired. |
| R3 #3 — Unchanged new corners | RESOLVED | Changed adjacency causes pinning. |
| R3 #4 — Pocket sampling grid | RESOLVED | Covering samples and uncertified-floor handling address the reported gaps. |
| R3 #5 — Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6 — Imported-package identity | RESOLVED | Actual imported paths are recorded. |
| R3 #7 — Exhausted iterators | RESOLVED | Both helpers materialize their inputs. |
| R4 #1 — Missed qualifying disc | RESOLVED | Conservative erosion replaces sampled absence decisions. |
| R4 #2 — Incorrect driver size floor | RESOLVED | `H_FLOOR` includes the base ambient minimum. |
| R4 #3 — Validation disabled with repair | RESOLVED | Final land validation is outside the repair switch. |
| R4 #4 — Fine basin behind coarse mouth | RESOLVED | Conservative bounds protect the reproduced basin. |
| R4 #5 — Pocket fine-zone exclusion | RESOLVED | Conservative bounds protect the fine zone. |
| R4 #6 — Non-Git identity | RESOLVED | Path, distribution and source digest are recorded. |
| R4 #7 — Invalid sizes/factors | RESOLVED | Reported invalid inputs are rejected. |
| R4 #8 — Missing test prerequisites | RESOLVED | Revised fixtures exercise the intended guards. |
| R4 #9 — Priority documentation | RESOLVED | Guide states that local refinement ignores priority. |
| R5 #1 — Grid coverage | RESOLVED | Boundary-adjacent samples provide coverage. |
| R5 #2 — Invalid slope premise | RESOLVED | Base and patch fields supply structural bounds. |
| R5 #3 — Pocket area/clearance | RESOLVED | Both predicates use upper bounds. |
| R5 #4 — Nested land collections | RESOLVED | Land validation recursively visits polygons. |
| R5 #5 — Invalid sampled sizes | RESOLVED | Validation precedes clipping. |
| R5 #6 — Unreadable source hashing | RESOLVED | Read and scan failures invalidate the digest. |
| R5 #7 — Stale floor comment | RESOLVED | Comment names `H_FLOOR`. |
| R6 #1 — Hidden enumeration failures | RESOLVED | `Path.walk(on_error=...)` records failures. |
| R6 #2 — Exterior disc centres | RESOLVED | Padding accounts for exterior centres; unsupported bounds refuse closure. |
| R6 #3 — Nested island collections | RESOLVED | `island_rings` recursively visits polygons. Finding 1 concerns polygon interiors. |
| R6 #4 — Empty filter footprint | RESOLVED | Explicit descriptive rejection. |
| R6 #5 — Ineffective fine-spot test | RESOLVED | Revised fixture and negative control exercise the bound. |
| R6 #6 — Opt-in comments | RESOLVED | Reported comments describe defaults correctly. |
| R7 #1 — Fully free islands | RESOLVED | Source replacement handles the original smaller-island and absent-island cases. |
| R7 #2 — Discontinuous exterior sizes | RESOLVED | Exterior-padding regressions pass. |
| R7 #3 — Estimated default floors | RESOLVED | Estimated floors no longer override conservative lower bounds. |
| R7 #4 — Degenerate ring | RESOLVED | Reported removal/blunting paths retain valid polygons. |
| R7 #5 — Closed-ring arc selection | RESOLVED | Uneven-density and ring-start regressions pass. |
| R7 #6 — Obsolete driver comment | RESOLVED | Comment describes filter → rim → hole. |
| R8 #1 — Snapped slit self-intersection | RESOLVED | Candidate validation rejects the reported crossing. |
| R8 #2 — Equal-length arcs | RESOLVED | Fit and geometric tie-breaking remove reported start dependence. |
| R8 #3 — Filtered-away island | RESOLVED | Empty shoreline is permitted for islands-only rims. |
| R8 #4 — Empty pocket land | RESOLVED | Empty land and mixed collections are handled. |
| R8 #5 — Failed sidecar enumeration | RESOLVED | Failure marker forces a null digest. |
| R8 #6 — Overwritten rejection diagnostics | RESOLVED | Rejection filenames include search pass. |
| R8 #7 — Stale fill counts | RESOLVED | Log counts assembled constraints. |
| R8 #8 — Premature “accepted” log | RESOLVED | Log says “selected seed”; acceptance remains gated. |
| R9 #1 — Initial blunting pins | RESOLVED | Initial pins survive subsequent root movement and folding. |
| R9 #2 — Changed nesting | RESOLVED | Every listed rim-edit path checks matched containment relationships. |
| R9 #3 — Mixed pocket land | RESOLVED | Recursive polygon extraction handles the collection. |
| R9 #4 — Dead wall-component code/docstring | RESOLVED | Dead code removed; helper and call-site descriptions agree. |
| R9 #5 — Obsolete resolution helper/comment | RESOLVED | Helper removed and comment corrected. |
| R10 #1 — Retreat changes nesting | RESOLVED | Retreat uses `_rim_edit_ok`; regression passes. |
| R10 #2 — Existing pins lost during rooting | RESOLVED | `_pin` pins the destination after movement or folding. |
| R10 #3 — Fold collapses triangle | RESOLVED | Fold candidate validation rejects collapse; regression passes. |
| R10 #4 — Withdrawal comment | RESOLVED | Comment describes withdrawing one shortest edge at a time. |
| R11 #1 — Unguarded root movement | RESOLVED | Movement requires `_rim_edit_ok`; regression passes. |
| R11 #2 — Fold crosses another ring | RESOLVED | The guard checks every pair of candidate ring boundaries. |
| R11 #3 — Rings exchange roles | RESOLVED | Injective ring matching preserves individual containment relationships. |
| R11 #4 — Newly folded roots unpinned | RESOLVED | `_pin` pins the new root and its surviving neighbors; regression passes. |

All three new findings predate the latest commit. **No defect introduced specifically by the round-11 fixes was substantiated.**

**Verification**

Used `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`, with bytecode writing disabled.

- `-m pytest -q tests/ -p no:cacheprovider` failed before collection: no writable temporary directory.
- `-m pytest -q tests/test_review_r11.py --capture=sys --noconftest -p no:cacheprovider`: **5 passed**.
- An in-memory pytest harness using the existing Matplotlib cache ran the scoped package, driver and review tests: **297 passed, 98 skipped**. Filesystem-dependent and marked mesh-generation tests were skipped.
- All three findings were reproduced in memory. Finding 2 also had the single-change control described above.
- Final `git status --porcelain` was empty. No files changed, mesh generation ran, or batch jobs were submitted. Production QA, node/element counts and implied time steps were not rerun.

## Verdict
VERDICT: FAIL (0 blocker, 0 major, 3 minor, 0 nit)


### Prompt

```markdown
# Review request, round 12: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 31 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, f07c6d9 round 4, c88739e round 5, 806d011 round 6, e88208f round 7, 6fbe8e8 round 8, 3903361 round 9, 32cd93a round 10, the last commit round 11; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-10: see your round-11 status table.
Round 11 (FAIL 0/0/4/0): all fixed in the last commit (`git show HEAD`;
triage in `docs/coastline-rules-review-20260927.md`).
1-3. `patch._rim_edit_ok(old, new)` is now the one guard on every rim
   edit path (short edges, slits, retreat, both blunting paths, the
   driver's `_blunt` fold and `_root` move): every ring a valid polygon
   of positive area; no two ring boundaries intersecting; every ring
   matched to its old self by shared vertices (injective) and contained in
   the same matched rings as before. `_ring_is_polygon` and `_rim_depths`
   are gone.
4. driver `_pin`: a fold or a root move pins the point it makes and the
   rim points that now neighbour it.
Tests: `tests/test_review_r11.py`. All six hires recipes rebuild accepted
with 0 introduced; Kimitsu keeps its boundary but seed 0 now introduces 3
violations with the root corners pinned, so seed 1 is taken (dt 1.84 s,
was 2.02 s); the rest byte-identical. 968 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | minor | yes: code read -- `island_rings` resampled `g.exterior` only | yes | each lake of an added island becomes a ring too (depth two, water by the nesting parity `hole_polygon` uses), if its ring is valid, inside the island's and half an element clear of it; otherwise it is reported as a lake meshed as land. `n_lakes_added` in the report and the driver log. Tests: the reproduction (520,000 m2 of water, (500, 500) water), and a lake too close to its coast |
| 2 | minor | yes: code read -- `done_refused` was seeded from every earlier refusal | yes | refusals are remembered within one pass only; the history (`n_short_edges_refused`) is deduplicated. Test: the reproduction leaves no edge under half an element |
| 3 | minor | yes: code read -- the report counted refusals | yes | `short_edges_left` is read off the returned rim (edges under `min_edge_factor` of the local size, within the focus). Test: the count and lengths match the returned rim |

All six hires recipes rebuild accepted and byte-identical to round 11 (none has an island lake). 972 tests pass.

## Round 13

### Answer (verbatim)

1. **Minor — Independently resampled lakes can overlap, breaking hole assembly.**

   **Location:** [patch.py:2117](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2117)

   Each lake is checked against its island’s exterior, but not against the other resampled lakes. Valid, disjoint source holes can therefore produce incompatible constraints.

   **Reproduction:**
   ```python
   import numpy as np
   from shapely.geometry import Polygon, box
   from fvcom_mesh_tools.patch import island_rings
   from tests.test_review_r12 import _with_rings

   h = lambda q: np.full(len(q), 30.)
   water = box(0, 0, 1000, 1000)
   lake1 = Polygon([
       (200,200), (700,200), (700,440), (400,440),
       (400,450), (700,450), (700,700), (200,700),
   ])
   lake2 = box(580, 441, 620, 442)
   land = Polygon(
       box(100,100,900,900).exterior.coords,
       [lake1.exterior.coords, lake2.exterior.coords],
   )
   assert land.is_valid
   rings, report = island_rings(land, water, h)
   _with_rings(water, rings)
   ```

   The report claims **two lakes added, zero skipped**. Every returned ring is individually valid, but assembly raises:
   ```text
   ValueError: the assembled hole boundary is 9822.08 m against
   9871.1 m of rim constraints; the rings touch or overlap
   ```

   **Fix:** Validate relationships between all candidate rings after resampling, including intersections and containment. Refuse or adjust incompatible lake rings and report the decision before returning them.

   This defect is introduced by the round-12 lake addition.

2. **Minor — Refusing a lake reverses the role of islands inside it.**

   **Location:** [patch.py:2080](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2080), [patch.py:2124](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2124)

   After refusing a lake and declaring it land, `island_rings` still independently adds land polygons inside that lake. Without their enclosing water ring, these become **water** through nesting parity.

   **Reproduction:** With constant `h=30` and water `box(0,0,1000,1000)`, use:
   ```python
   outer = Polygon(
       box(100,100,900,900).exterior.coords,
       [box(105,300,700,700).exterior.coords],
   )
   land = MultiPolygon([outer, box(300,400,400,500)])
   ```

   This is valid source geometry. The function reports **two islands added, zero lakes added**, and says the lake is meshed as land. Assembling the returned rings instead makes the inner island’s **10,000 m² water**. The subsequent land check rejects this avoidable error.

   Using `box(300,400,310,410)` for the inner island produces **50 m² of land classified as water**, with `land_an_element_fits(...) == []`.

   **Fix:** Preserve a hierarchy of source rings and their land/water roles. When a lake is filled, suppress its now-redundant descendant island boundaries; do not append them as independent rings.

   The independent-polygon behavior predates round 12; the lake-refusal path remains incomplete.

3. **Minor — Lake clearance uses the smallest size anywhere on the ring instead of the local size.**

   **Location:** [patch.py:2112](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2112), [patch.py:2118](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/patch.py:2118)

   A fine-size corner can authorize insufficient clearance elsewhere on the lake. The existing island clearance check at lines 2094–2095 has the same problem.

   **Reproduction:**
   ```python
   h = lambda q: np.maximum(
       30., 30. + .1*((q[:,0]-150) + (q[:,1]-150))
   )
   land = Polygon(
       box(100,100,900,900).exterior.coords,
       [box(150,150,850,850).exterior.coords],
   )
   rings, report = island_rings(
       land, box(-1000,-1000,2000,2000), h
   )
   ```

   Observed: **one lake added, zero skipped**. The returned lake retains `(850,850)`, where `h=170 m`. Its distance from the island coast is **50 m**, below the required **85 m**. The test passes because another vertex has `h=30 m`, reducing the threshold everywhere to 15 m.

   **Fix:** Check clearance against the size along the relevant boundary, using conservative bounds between samples. A conservative whole-ring upper bound is also safe, though more restrictive.

   Round 12 introduces this defect for lakes by repeating the existing island test.

4. **Minor — Non-Git source hashes silently omit symlinked subpackages.**

   **Location:** [provenance.py:129](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:129)

   `Path.walk()` does not follow directory symlinks by default. It lists them as filenames, and the `.py` suffix filter discards them. Python can nevertheless import code through those directories.

   **Evidence:** An in-memory traversal mock represented:
   ```text
   pkg/
     __init__.py
     subpkg -> external Python package directory
   ```

   Changing the mocked contents of `subpkg/__init__.py` produced the **same `source_sha256`** both times. Only the top-level `__init__.py` was hashed, and `source_unreadable` remained `[]`.

   **Fix:** Either follow directory symlinks with cycle detection and hash their contents under logical package paths, or explicitly mark the source identity incomplete when encountering them.

   This predates the latest commit.

**Status of every previous finding**

Statuses distinguish the original reproductions from the additional cases above.

| Finding | Status | Reason |
|---|---|---|
| R1 #1 — Preserve mode | RESOLVED | Repair and retry require `resolve`. |
| R1 #2 — Lagoon entrance | RESOLVED | Connectivity guard protects the reported entrance. |
| R1 #3 — Refused-pocket walls | RESOLVED | Restoration checks the pocket polygon’s remaining water area. |
| R1 #4 — Three-node collapse | RESOLVED | Candidate guards reject degenerate rings. |
| R1 #5 — Retreat-only retry | RESOLVED | Array comparison detects movement; unselected state is restored. |
| R1 #6 — Empty land/padding | RESOLVED | Empty-land handling and exterior padding cover the reported cases. |
| R1 #7 — Invalid continuous-width parameters | RESOLVED | Reported invalid arguments and empty footprints are handled. |
| R1 #8 — Ignored seed argument | RESOLVED | Search iterates `seed_list`. |
| R1 #9 — GPL provenance import | RESOLVED | Versions use distribution metadata without importing GPL code. |
| R1 #10 — Missing provenance inputs/settings | RESOLVED | Bathymetry, region datasets and relevant settings are recorded. |
| R1 #11 — Untracked source | RESOLVED | Git status includes untracked files; entry-point tracking is recorded. |
| R1 #12 — Retry reproducibility | RESOLVED | Documentation requires the original seed sequence for pass 2. |
| R1 #13 — Ineffective fixtures | RESOLVED | Revised fixtures exercise prerequisites; unreachable merge code was removed. |
| R1 #14 — Raster documentation | RESOLVED | Resolution is stated relative to the finest target. |
| R2 #1 — Resolvable land accepted as water | RESOLVED | Land validation rejects qualifying land inside the hole. |
| R2 #2 — All-land raster | RESOLVED | Empty-water early return handles it. |
| R2 #3 — Detached artificial island | RESOLVED | Closure requires neighboring land. |
| R2 #4 — Fine region inside pocket | RESOLVED | Conservative bounds protect the reported fine region. |
| R2 #5 — Obsolete sliding curves | RESOLVED | Edited corners, folded roots and moved roots are pinned. |
| R2 #6 — Geometry provenance | RESOLVED | Sidecars are enumerated; enumeration failure invalidates the digest. |
| R2 #7 — Falsy experimental values | RESOLVED | Explicit list validation rejects them. |
| R2 #8 — Seed-variable documentation | RESOLVED | Header names `LR_SEEDS`. |
| R3 #1 — Variable-size land guard | RESOLVED | Decisions use conservative size bounds. |
| R3 #2 — Repairs bypass validation | RESOLVED | Final-rim and retry land checks remain wired. |
| R3 #3 — Unchanged new corners | RESOLVED | Changed adjacency causes pinning. |
| R3 #4 — Pocket sampling grid | RESOLVED | Covering samples and uncertified-floor handling address the reported gaps. |
| R3 #5 — Uppercase sidecars | RESOLVED | Extension matching is case-insensitive. |
| R3 #6 — Imported-package identity | RESOLVED | Actual imported package paths are recorded. |
| R3 #7 — Exhausted iterators | RESOLVED | Both helpers materialize their inputs. |
| R4 #1 — Missed qualifying disc | RESOLVED | Conservative erosion replaces sampled absence decisions. |
| R4 #2 — Incorrect driver size floor | RESOLVED | `H_FLOOR` includes the base ambient minimum. |
| R4 #3 — Validation disabled with repair | RESOLVED | Final land validation is outside the repair switch. |
| R4 #4 — Fine basin behind coarse mouth | RESOLVED | Conservative bounds protect the reported basin. |
| R4 #5 — Pocket fine-zone exclusion | RESOLVED | Conservative bounds protect the fine zone. |
| R4 #6 — Non-Git identity | PARTIAL | Ordinary source trees are identified; symlinked subpackages are omitted, finding 4. |
| R4 #7 — Invalid sizes/factors | RESOLVED | The reported invalid inputs are rejected. |
| R4 #8 — Missing test prerequisites | RESOLVED | Revised fixtures exercise the intended guards. |
| R4 #9 — Priority documentation | RESOLVED | Guide states that local refinement ignores priority. |
| R5 #1 — Grid coverage | RESOLVED | Boundary-adjacent samples provide coverage. |
| R5 #2 — Invalid slope premise | RESOLVED | Base and patch fields supply structural bounds. |
| R5 #3 — Pocket area/clearance | RESOLVED | Wall-pocket predicates use upper bounds; finding 3 concerns `island_rings`. |
| R5 #4 — Nested land collections | RESOLVED | Land validation recursively visits polygons. |
| R5 #5 — Invalid sampled sizes | RESOLVED | Validation precedes clipping. |
| R5 #6 — Unreadable source hashing | RESOLVED | File-read failures invalidate the digest. |
| R5 #7 — Stale floor comment | RESOLVED | Comment names `H_FLOOR`. |
| R6 #1 — Hidden enumeration failures | RESOLVED | `Path.walk(on_error=...)` records scan failures; finding 4 concerns skipped symlinks. |
| R6 #2 — Exterior disc centres | RESOLVED | Padding accounts for exterior centres; unsupported bounds refuse closure. |
| R6 #3 — Nested island collections | RESOLVED | Polygon traversal is recursive; nesting-policy gaps are reported separately. |
| R6 #4 — Empty filter footprint | RESOLVED | Explicit descriptive rejection. |
| R6 #5 — Ineffective fine-spot test | RESOLVED | Revised fixture and negative control exercise the bound. |
| R6 #6 — Opt-in comments | RESOLVED | Reported comments describe defaults correctly. |
| R7 #1 — Fully free islands | RESOLVED | Source replacement handles the original smaller-island and absent-island cases. |
| R7 #2 — Discontinuous exterior sizes | RESOLVED | Exterior-padding regressions pass. |
| R7 #3 — Estimated default floors | RESOLVED | Estimated floors no longer override conservative lower bounds. |
| R7 #4 — Degenerate ring | RESOLVED | Reported removal and blunting paths retain valid polygons. |
| R7 #5 — Closed-ring arc selection | RESOLVED | Uneven-density and ring-start regressions pass. |
| R7 #6 — Obsolete driver comment | RESOLVED | Comment describes filter → rim → hole. |
| R8 #1 — Snapped slit self-intersection | RESOLVED | Candidate validation rejects the reported crossing. |
| R8 #2 — Equal-length arcs | RESOLVED | Fit and geometric tie-breaking remove reported start dependence. |
| R8 #3 — Filtered-away island | RESOLVED | Empty shoreline is permitted for islands-only free rims. |
| R8 #4 — Empty pocket land | RESOLVED | Empty land and mixed collections are handled. |
| R8 #5 — Failed sidecar enumeration | RESOLVED | Failure marker forces a null digest. |
| R8 #6 — Overwritten rejection diagnostics | RESOLVED | Rejection filenames include search pass. |
| R8 #7 — Stale fill counts | RESOLVED | Log counts assembled constraints. |
| R8 #8 — Premature “accepted” log | RESOLVED | Log says “selected seed”; acceptance remains gated. |
| R9 #1 — Initial blunting pins | RESOLVED | Initial pins survive subsequent root movement and folding. |
| R9 #2 — Changed nesting | RESOLVED | Rim-edit paths check matched containment relationships. |
| R9 #3 — Mixed pocket land | RESOLVED | Recursive polygon extraction handles the collection. |
| R9 #4 — Dead wall-component code/docstring | RESOLVED | Dead code removed; helper and call-site descriptions agree. |
| R9 #5 — Obsolete resolution helper/comment | RESOLVED | Helper removed and comment corrected. |
| R10 #1 — Retreat changes nesting | RESOLVED | Retreat uses `_rim_edit_ok`; regression passes. |
| R10 #2 — Existing pins lost during rooting | RESOLVED | `_pin` pins the destination after movement or folding. |
| R10 #3 — Fold collapses triangle | RESOLVED | Fold candidate validation rejects collapse; regression passes. |
| R10 #4 — Withdrawal comment | RESOLVED | Comment describes withdrawing one shortest edge at a time. |
| R11 #1 — Unguarded root movement | RESOLVED | Movement requires `_rim_edit_ok`; regression passes. |
| R11 #2 — Fold crosses another ring | RESOLVED | Guard checks candidate ring boundaries against one another. |
| R11 #3 — Rings exchange roles | RESOLVED | Injective ring matching preserves individual containment relationships. |
| R11 #4 — Newly folded roots unpinned | RESOLVED | New root and surviving neighbors are pinned; regression passes. |
| R12 #1 — Source islands lose lakes | PARTIAL | Original single-lake reproduction passes; findings 1–3 expose remaining lake-handling gaps. |
| R12 #2 — Refused edges never reconsidered | RESOLVED | Refusal eligibility resets each pass; original regression passes. |
| R12 #3 — Stale remaining-edge diagnostics | RESOLVED | Diagnostics use returned edges within the requested focus; history is separate. |

**Verification**

Used `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python` with bytecode writing disabled.

- Full command, `-m pytest -q tests/ -p no:cacheprovider`, failed before collection because no writable temporary directory was available.
- Round-11/12 tests with `--capture=sys --noconftest -p no:cacheprovider`: **9 passed**.
- An in-memory pytest harness ran scoped package, driver and review tests using the existing Matplotlib cache: **299 passed, 75 skipped**. Filesystem-writing and mesh-generation tests were skipped.
- All four findings were reproduced in memory.
- `git diff --check` passed; final `git status --porcelain` was empty.

No files changed, batch jobs were submitted, or mesh generation ran. Production QA, node/element counts, implied time steps and the author’s six-recipe rebuild results were not rerun.

## Verdict
VERDICT: FAIL (0 blocker, 0 major, 4 minor, 0 nit)


### Prompt

```markdown
# Review request, round 13: coastline rules for local refinement (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (use `/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python`;
the tests run with `.../bin/python -m pytest -q tests/`). Do not submit
batch jobs and do not run mesh generation (`notebooks/420_local_refine.py`
needs the compute nodes). Answer in English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
`git diff 0c5d9b3..HEAD` (now 32 commits; 7a6cffe fixed round 1, 2d8d1d7 round 2, fd2bbe1 round 3, f07c6d9 round 4, c88739e round 5, 806d011 round 6, e88208f round 7, 6fbe8e8 round 8, 3903361 round 9, 32cd93a round 10, 18e7071 round 11, the last commit round 12; read `git log 0c5d9b3..HEAD` for the
reasons). The local-refinement driver `notebooks/420_local_refine.py`
refines a Tokyo Bay FVCOM mesh inside a region, re-cuts the coastline from
OSM (`hires.coastline: resolve`) and must pass QA with 0 violations
introduced. Changes in scope:

- `src/fvcom_mesh_tools/patch.py`
  - `rim_repair`: short edges (remove a free end, or merge a cap between
    two corners to its midpoint), slits (throat under one element; close a
    dead end no element fits in, otherwise step a pier tip back), angles
    (`blunt_acute_corners` again with wall roots protected), `focus=` for a
    QA-feedback retry; returns a remap for wall references.
  - `unresolvable_water` and `filter_shoreline_local(continuous_width=)`:
    coarse-zone water judged at the local size by a distance transform
    (radius 0.75 h, dead ends only, straits left open).
  - `filter_shoreline_local(keep_land=)`: in coarse bands a removed land
    piece is kept whole if most of it is base land.
  - `_source_substring`: the arc of a closed ring that fits the stretch.
  - `blunt_acute_corners(protect=)`; `island_rings` reports area inside.
  - Removal of `water_wedges`, `short_chords`, `seam_water`.
- `src/fvcom_mesh_tools/walls.py`: `close_wall_pockets`.
- `src/fvcom_mesh_tools/provenance.py` (new): commits, file hashes,
  library versions recorded in `report.json`.
- `src/fvcom_mesh_tools/refine.py`: `hires.rim_repair`,
  `continuous_width`, `keep_base_land`, `wall_pockets` (all default true),
  `hires.experimental` list.
- `notebooks/420_local_refine.py`: wiring; `LR_EXPERIMENTAL` override;
  `apply_rim_repair` / `assemble_constraints`; `seed_search` and the
  one-shot retry near the offenders with snapshot/restore; flat-face drop
  after the fill; wall-pocket restore when an island is refused;
  `rejected_seed<k>.npz`; provenance.
- `jobs/octopus/417_hires_refine.sh` (LR_EXPERIMENTAL documented).
- Tests: `tests/test_patch.py`, `test_walls.py`, `test_refine.py`,
  `test_provenance.py`, `test_local_refine_driver.py`,
  `test_review_hires_fixes.py`.
- Docs: `docs/USER_GUIDE.md` §4, §5 (reproducibility), §11; `CHANGELOG.md`;
  `recipes/refine/*.yaml` comments.

Evidence already gathered by the author (on compute nodes): all six hires
recipes accepted with the four rules on, byte-identical rebuilds from the
accepted seed, FVCOM smoke runs and 20-day M2 runs pass.

Out of scope: the rest of the package, unless these changes touch it.

## Previous rounds
Rounds 1-11: see your round-12 status table.
Round 12 (FAIL 0/0/3/0): all fixed in the last commit (`git show HEAD`;
triage in `docs/coastline-rules-review-20260927.md`).
1. `island_rings` adds each lake of an added island as a ring when valid,
   inside the island ring and half an element clear; otherwise reported
   ("a lake in an island ... meshed as land"). `n_lakes_added`.
2. `rim_repair` remembers refusals within one short-edge pass only;
   history deduplicated as `n_short_edges_refused`.
3. `short_edges_left` / `n_short_edges_left` are computed from the
   returned rim.
Tests: `tests/test_review_r12.py`. All six hires recipes rebuild accepted,
byte-identical to round 11. 972 tests pass.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the whole scope above and everything it
   touches.

## Severity
- blocker: produces wrong scientific results or loses data in normal use
- major: a failure or wrong result that can be accepted as success, in a
  realistic path
- minor: needs unusual input or an injected fault, or is a clear
  robustness/clarity defect
- nit: style, wording, dead code

## Required output
Numbered findings, each with severity, file:line, a reproduction or
evidence, and a concrete fix. Then `## Verdict` with exactly one line:
`VERDICT: PASS` (no finding of any severity) or
`VERDICT: FAIL (<n> blocker, <n> major, <n> minor, <n> nit)`.
```

### Triage

| id | severity | verified? (how) | correct? | action |
|---|---|---|---|---|
| 1 | minor | yes: code read -- lakes were checked against their island only | yes | every ring `island_rings` returns must not touch one already placed; a lake that would is filled and reported. Test: the two-lake reproduction assembles |
| 2 | minor | yes: code read -- polygons inside a refused lake were added on their own | yes | polygons are taken outside in (by area); one inside a lake filled as land is not added ("inside a lake meshed as land"). Tests: an island in a filled lake stays land; one in a kept lake is added |
| 3 | minor | yes: code read -- clearance used the finest size on the ring | yes, in substance; the proposed consequence (refuse) is not taken where it would do worse | `_clear_of` judges every edge at the size where it meets the gap (the foot of its shortest line to the other boundary) and returns the tightest edge. Three tiers: clear at the local size -> added; a gap under half the FINEST of those sizes carries no element -> left out / lake filled; in between, when an element fits on the piece, it is kept and named in `report["tight"]` (a driver WARNING), because leaving it out meshes land an element fits on as water (the driver then refuses the run) or fills a lake an element fits in as land. Found on the way: bounding the size along the edge (land, coarse nearest-node values), along the whole shortest line (ends on the coarse frozen interface) or at its middle (crosses the island from the far side) each refused a Funabashi island 159 m from the rim that has always meshed cleanly; at the foot the size is 367 m, so it is now kept and named. Tests: the reproduction's lake is named tight at > 150 m; a 5 m gap at 30 m is left out / filled |
| 4 | minor | yes: `Path.walk` does not descend into directory symlinks | yes | `code_state` walks with `follow_symlinks=True`, each real directory once. Test: a symlinked subpackage's change changes the digest; a cycle terminates |

All six hires recipes rebuild accepted and byte-identical to round 12; Funabashi's island is named in `tight`. 978 tests pass.
