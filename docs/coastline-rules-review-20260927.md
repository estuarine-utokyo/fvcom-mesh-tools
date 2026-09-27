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
