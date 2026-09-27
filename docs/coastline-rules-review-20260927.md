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
