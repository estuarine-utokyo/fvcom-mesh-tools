# Review: extending a base mesh outward (2026-10-01)

Reviewer gpt-6-astra (rounds 1-3) and gpt-6.1-sol (from round 4), Codex CLI, read-only; process in
`~/.claude/skills/codex-review`. Scope: the extension tool (notebooks
444-448, 453; `extend.py`, `extend_recipe.py`, `obc_design.py`,
`dem/sources.py`), its jobs, and the portability change (6d8b9a7, 6c068d2).

## Round 1

### Answer (verbatim)

Read-only audit; no files changed. The 16 existing pure geometry tests passed, scoped Ruff checks passed, and every OCTOPUS job script passed `bash -n`. Additional reproductions used small arrays, in-memory NetCDF and mocks. Production mesh generation, production QA and FVCOM integration were not run.

1. **Major — Masked NetCDF cells become valid bathymetry.** [dem/sources.py:89](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:89)  
   `np.asarray` discards the NetCDF mask while retaining its underlying fill values. An in-memory grid with a masked elevation of `-9999` returned depth `9999`, attributed to the first source; the valid fallback was never consulted. These values can survive clipping or contaminate smoothing. **Fix:** convert masked values to floating-point NaN before interpolation; test masked cells and fallback selection.

2. **Major — Failed final QA still produces a successful build exit.** [447_extend_merge.py:176](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:176)  
   A mocked execution with one failed gate printed `QA 21/22`, wrote the reports and returned normally. Notebook 445 consequently accepts the subprocess and prints `done`. This is a realistic path: the reviewed history documents builds with failed gates. **Fix:** preserve diagnostic reports, then exit nonzero whenever required gates fail. Test both successful and failed orchestration paths.

3. **Major — The final depths need not satisfy the reported r-factor.** [447_extend_merge.py:152](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:152)  
   The cap is reapplied after smoothing, but `r_after` is retained from before that cap. With fixed depth `100`, free depth `10`, cap `10` and `rmax=.2`, smoothing reports `.2`; the exported pair has r=`.8181818`. Separately, fixed depth `1` and free-node minimum `3` returns after 5,000 iterations with r=`.5`, which callers also accept. Notebook 453 repeats both problems. **Fix:** enforce bounds during smoothing, reject infeasible/nonconverged cases, and recompute r after every final transformation.

4. **Major — Constrained sizing bands can silently defeat the CFL floor.** [extend.py:88](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:88)  
   A synthetic lattice with floor `4000` and band target `1500` returned `1500` everywhere and merely reported `below_floor_fraction=1`. Notebook 446 accepts this report. Final QA does not receive `min_dt_s`, so it cannot enforce the recipe’s timestep promise after finishing or depth changes either. **Fix:** detect incompatible constraints explicitly and gate the final new elements against the requested timestep using final depths.

5. **Major — Boundary resampling does not guarantee the advertised spacing floor.** [obc_design.py:134](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:134)  
   `resample([[0,0],[8000,0]], 3000)` produces edges `[3000,3000,2000]`. Curved sections additionally shorten arc-length steps into chords, and variable spacing is checked only at the preceding node. Notebook 444 never checks the resulting edges against their depth-dependent floor. This contradicts USER_GUIDE §13. **Fix:** enforce the floor on actual output chords, account for depths along each edge, and redistribute the terminal remainder.

6. **Major — Missing sizing bathymetry is treated as zero depth.** [446_extend_generate.py:99](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:99)  
   Uncovered samples become zero through `nan_to_num`, removing the CFL constraint there. Notebook 444 does the same at line 100. Recipes explicitly allow different sizing and depth stacks, so later depth sampling can succeed while sizing was based on fictitious shallow water. **Fix:** reject uncovered wet sizing locations, or use an explicit conservative fallback and report its extent.

7. **Minor — The seaward-normal check can select an inland direction, and validation misses that crossing.** [obc_design.py:57](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:57), [444_design_obc.py:121](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:121)  
   For a 1-km-thick rectangular land strip, reversing polygon winding changed the returned normal from `180°` to `0°`: the latter crosses land but its distant test point has already emerged beyond it. Notebook 444 excludes the entire first and last edges from its land-intersection check; a synthetic 1-km crossing consequently measured zero. **Fix:** test the departure segment/local interior and validate the full boundary, exempting only endpoint contact.

8. **Major — A rejected boundary design overwrites the previous usable boundary.** [444_design_obc.py:149](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:149)  
   The CSV and report are written before `if bad: raise SystemExit`. A failed redesign therefore replaces an existing valid input with a rejected one. The extension recipe loader reads the CSV without inspecting its report’s `problems`. **Fix:** complete validation before publishing the CSV; preserve existing products on failure and publish successful products atomically.

9. **Minor — Export can violate the frozen-coordinate contract without detection.** [447_extend_merge.py:173](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:173)  
   Read-back checks only depths. The native writer uses eight decimal places: a base coordinate `0.123456789` passed the pre-export frozen check but read back as `0.12345679`. Higher-precision base depths also cannot round-trip through the six-decimal writer, although their mismatch is detected. **Fix:** use round-trip-safe serialization for frozen values and run the complete frozen-base verification on the exported case.

10. **Minor — Seam verification accepts overlapping elements on the same side.** [extend.py:158](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:158)  
    For the unit-square-style base used in the tests, adding a triangle sharing the east interface but with its third vertex inside the base passed `merge_outer` and `verify_frozen_base`. Counting one base and one outer owner does not establish a geometrically valid seam. **Fix:** require opposite-side incidence across every interface edge and reject overlap between the base and outer footprints.

11. **Minor — A one-edge land segment disappears.** [extend.py:196](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:196)  
    `land_segments([[0,1,2]], [[0,1,2]])` returns `[]`, although edge `2–0` is land. The accumulator holds starting vertices, so one valid land edge has `len(run)==1` and is discarded. **Fix:** emit any nonempty run with its terminating vertex; add one-edge and multiple-open-chain tests.

12. **Minor — Accepted short interfaces cannot reach generation.** [446_extend_generate.py:133](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:133)  
    Both ladders unconditionally use `skip_ends=2`. `build_obc_band` rejects chains with fewer than six nodes, while the recipe reader and merge API accept two-node boundaries/interfaces. **Fix:** adapt the ladder to chain length or reject unsupported lengths during preflight with a clear explanation.

13. **Minor — Cabinet Office cache entries leak between data roots.** [dem/sources.py:164](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:164)  
    The cache key contains only `(zone, area)`. A mocked first root containing depth `10`, followed by another root containing `20`, returned `10` for the second root without reading it. The registry keeps a process-wide `CaoNested` instance. **Fix:** include the resolved source identity and expected dimensions in the cache key.

14. **Minor — Invalid fine-grid samples erase valid coarse-grid coverage.** [dem/sources.py:198](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:198)  
    Coarse depth `10` followed by an overlapping fine grid containing NaN returned NaN. The finer interpolation overwrites `out` and updates `best` without checking validity. **Fix:** replace the coarser result only where the finer interpolation is finite; retain coarse coverage elsewhere.

15. **Minor — Degenerate M7001 windows abort the entire priority stack.** [dem/sources.py:119](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:119)  
    The minimum-point check occurs before duplicate removal. Three rows containing only two distinct locations reproduce `QhullError`; three or more collinear locations have the same problem. A following fallback source cannot run. **Fix:** validate distinct-point count and geometric rank, returning uncovered samples where interpolation is impossible.

16. **Major — Provenance omits the Cabinet Office depth data.** [dem/sources.py:132](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:132)  
    `CaoNested.files()` lists only calculation-area spreadsheets. Notebook 445 therefore hashes those tables but none of the `depth_*.dat` files actually used. Changing the default recipe’s principal bathymetry can leave its recorded source hashes unchanged. **Fix:** record and hash the exact depth files accessed by each generation/depth stage, alongside the tables.

17. **Minor — Dirty tracked code is incorrectly identified by its commit.** [provenance.py:126](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:126)  
    With `path_tracked=True` and `dirty=['extend.py']`, `code_state` returned `commit_identifies_code=True` and performed no source hashing. Different edits to that same file produce indistinguishable code identities. **Fix:** hash relevant working-tree sources whenever tracked or untracked changes can affect execution; do not mark the commit as sufficient.

18. **Minor — Recipe numeric validation accepts nonfinite and fractional control values.** [extend_recipe.py:73](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend_recipe.py:73)  
    Mocked recipes accepted `cfl_dt_s: .nan`, `max_iter: 0.5` and `gen_seed: 0.5`. The latter values are subsequently truncated with `int`, so effective execution differs from the recorded recipe. **Fix:** require finite real settings, integral iteration counts/seeds, valid seed ranges and finite ordered bounds.

19. **Minor — NaN spacing makes resampling loop indefinitely.** [obc_design.py:135](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:135)  
    For NaN, both the positivity check and termination comparison are false. A bounded callback reproduction observed repeated queries at `[nan,nan]`; execution otherwise keeps appending NaNs. **Fix:** reject nonfinite spacing and coordinates and verify that each iteration advances a finite distance.

20. **Major — Concurrent extension builds can write into the same output directory.** [445_extend_mesh.py:41](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:41)  
    Checking emptiness and then calling `mkdir(exist_ok=True)` does not reserve the directory. Two submissions using the same default recipe/output can both pass and overwrite each other’s generation files, meshes and reports. The refinement CLI already addresses this class of race with an exclusive reservation. **Fix:** reserve extension outputs atomically before launching either stage.

21. **Major — Re-depthing can overwrite its own source case or an existing experiment.** [453_redepth_extended.py:80](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:80)  
    There is no nonempty-output or source-alias check. Setting `OUTDIR=BUILT_DIR` and `--case-name` to the original case overwrites the input case. Reusing another output directory replaces its products and report without preserving the previous experiment. **Fix:** reject input/output aliases and reserve a fresh output directory before sampling or exporting.

22. **Major — Reused smoke directories can validate stale output.** [448_extend_smoke.py:74](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/448_extend_smoke.py:74)  
    Staging preserves existing `extended/output` files while replacing inputs. The checker does not bind history to those inputs or to the current invocation. In-memory old history with the requested dates and finite fields passed despite being unrelated to the staged mesh. **Fix:** require a fresh run directory and associate success with the staged-input identity; do not mix history from different attempts.

23. **Major — Smoke timestep selection precedes a depth change that can invalidate it.** [448_extend_smoke.py:62](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/448_extend_smoke.py:62)  
    `DTE` is selected before `apply_obc_depth_control`. On a nine-node synthetic mesh, the original allowance was `5.0482 s`, producing `DTE=5.0 s`; the adjusted mesh allowed only `3.0289 s`. The later manifest recomputes the allowance but does not reject the inconsistency. **Fix:** apply all depth transformations first, then select and verify the timestep.

24. **Minor — Relative smoke roots produce incorrect namelist paths.** [448_extend_smoke.py:72](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/448_extend_smoke.py:72)  
    `--root scratch/test` remains relative, so the namelist contains `scratch/test/extended/input/`. The job subsequently changes directory to `scratch/test/extended`; FVCOM resolves that path beneath the case directory again. **Fix:** resolve `--root` and `--case` immediately after parsing.

25. **Minor — The extension smoke job ignores the FVCOM executable override.** [448_extend_smoke.sh:25](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/448_extend_smoke.sh:25)  
    Unlike the other updated run jobs, this assignment ignores `FMESH_FVCOM` and always selects `$WORK_DIR/Github/FVCOM/src/fvcom`. An explicitly selected executable can therefore be silently replaced. **Fix:** use `${FMESH_FVCOM:-$WORK_DIR/Github/FVCOM/src/fvcom}` and record the chosen executable’s identity.

26. **Minor — The existing 383 job runs a case that preparation no longer creates.** [383_m2.sh:54](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/383_m2.sh:54)  
    Both its stale-output guard and run loop name `B_Adepth`; notebook 383 stages `B_m7001`. On a fresh run the third directory is missing. Moreover, the subshell is on the left of `||`, so `set -e` does not stop its commands after the failed `cd`. **Fix:** derive case names from the manifest and explicitly guard each directory change.

27. **Major — Existing package code still directly imports GPL oceanmesh.** [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332)  
    Direct imports also remain in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672` and `autofinish/directives.py:20`. These violate the expressly requested repository rule. `THIRD_PARTY_NOTICES.md:37` additionally endorses the conflicting arrangement. **Fix:** move those operations behind subprocess boundaries or into a separately licensed plugin, and reconcile the notices with the governing policy.

28. **Minor — Unvalidated case names can escape the protected output directory.** [extend_recipe.py:41](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend_recipe.py:41)  
    `case` is not restricted to a filename component. An absolute case prefix causes `outdir / f"{casename}_grd.dat"` to discard `outdir`; `../` can also escape it. Thus notebook 445’s empty-directory guard does not protect the actual destination. Notebook 453’s `--case-name` has the same issue. **Fix:** validate case names as nonempty basenames and assert that every resolved output stays inside its reserved directory.

29. **Minor — The native reader accepts scrambled open-boundary order.** [fvcom_native.py:566](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:566)  
    Validation compares sets rather than the ordered walk. On a nine-node square grid, the nonadjacent sequence `[0,2,1,5,8]` was accepted because it has the expected endpoints and node set. Extension generation uses this sequence directly for interface constraints and ladders. **Fix:** require a duplicate-free consecutive boundary walk, accepting either valid traversal direction.

## Verdict

VERDICT: FAIL (0 blocker, 13 major, 16 minor, 0 nit)

### Prompt

```markdown
# Review request, round 1: extending a base mesh outward (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (small synthetic inputs only; do not read the large data
under $DATA_DIR beyond listing it, and do not submit batch jobs). Answer in
English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
A tool that keeps a finished FVCOM base mesh exactly as it is and adds the
sea out to a new, designed open boundary (USER_GUIDE section 13). Read:

- `git show b0584f9 8e2739b 450ad44 d6d2a72 7044b6b 0f52d5b 69b50a4 d9e92fd b6d2ed8 765423c`
  (the extension tool and its documentation), and the current files:
  - `src/fvcom_mesh_tools/extend.py`, `extend_recipe.py`, `obc_design.py`,
    `dem/sources.py` (named bathymetry sources, priority stack, and the new
    `DATUM` registry / `non_tp_count` warning);
  - `notebooks/444_design_obc.py`, `445_extend_mesh.py`, `446_extend_generate.py`,
    `447_extend_merge.py`, `448_extend_smoke.py`, `453_redepth_extended.py`;
  - `recipes/extend/tokyo_bay_enshu.yaml`, `tokyo_bay_enshu_obc_design.yaml`;
  - `jobs/octopus/444_design_obc.sh`, `445_extend_mesh.sh`, `448_extend_smoke.sh`,
    `453_redepth_extended.sh`, `common.sh`;
  - tests: `tests/test_extend*.py`, `tests/test_obc_design*.py`,
    `tests/test_dem_sources.py` (whatever exists).
- Also in scope, just committed: the portability change --
  every job script and `common.sh` now take paths only from `$DATA_DIR` and
  `$WORK_DIR` (login profile), stop when they are unset, and derive the
  OCTOPUS FVCOM library directory as `FVCOM_LIBS` in `common.sh`; notebooks
  383/384/414 and `cli/refine_run.py` no longer fall back to `/octfs/...`.
  See commits 6d8b9a7 and 6c068d2 (`git log -5`).

Design intent:
- the base mesh's nodes, elements and depths are carried bit for bit
  (`verify_frozen_base`);
- the new part is generated with oceanmesh (subprocess; GPL code must never
  be imported by `fvcom_mesh_tools`), with fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Please
1. A fresh, unrestricted audit of the scope above and everything it touches:
   correctness (geometry, indexing, orientation, frozen-base contract,
   depth handling, units/datums), failure modes that could be accepted as
   success, reproducibility (seeds, ordering, MPI-independence where
   relevant), error handling, the job scripts under `set -euo pipefail`,
   the licence rules (GPL not imported; Cabinet Office data never written
   out as-is), tests (missing happy-path and error-case coverage), and
   whether the documentation (USER_GUIDE section 13, CHANGELOG) matches
   the code.

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

| id | sev | verified? (how) | correct? | action |
|---|---|---|---|---|
| F1 | major | code read; test with a masked grid | yes | 662b142 |
| F2 | major | code read (447 returned after QA failures) | yes | 3e61359 |
| F3 | major | code read; tests (cap after smoothing, infeasible bounds) | yes | 960b232 |
| F4 | major | code read; **the current Enshu mesh fails the new gate** (new elements 9.7 s vs base 14.4 s) | yes | report 960b232; 446 gate 65119e6; 447 gate 3e61359 |
| F5 | major | code read; test (8 km at 3 km gave 3000/3000/2000) | yes | 65119e6 |
| F6 | major | code read (nan_to_num -> 0) | yes | 65119e6 |
| F7 | minor | test with a 1 km strip | yes | 65119e6 |
| F8 | major | code read (CSV written before the check) | yes | 65119e6 |
| F9 | minor | code read (read-back checked depths only) | yes | 3e61359 |
| F10 | minor | test (third vertex inside the base) | yes | 960b232 |
| F11 | minor | test (one-edge run dropped) | yes | 960b232 |
| F12 | minor | code read (skip_ends=2 vs 6-node minimum) | yes | 65119e6 |
| F13 | minor | test (two data roots) | yes | 662b142 |
| F14 | minor | test (NaN fine grid over a coarse value) | yes | 662b142 |
| F15 | minor | test (two distinct / collinear points) | yes | 662b142 |
| F16 | major | code read (files() listed tables only) | yes | 662b142 |
| F17 | minor | test with a mocked dirty tracked tree | yes | 6ca8c51 |
| F18 | minor | tests (NaN, 0.5) | yes | 3e61359 |
| F19 | minor | test (NaN spacing) | yes | 65119e6 |
| F20 | major | code read (check-then-mkdir race) | yes | 3e61359 |
| F21 | major | code read (no alias / non-empty check) | yes | 3e61359 |
| F22 | major | code read (reused root kept old output) | yes | 3e61359 |
| F23 | major | code read (DTE before depth control) | yes | 3e61359 |
| F24 | minor | code read (relative root) | yes | 3e61359 |
| F25 | minor | code read (FMESH_FVCOM ignored) | yes | 3e61359 |
| F26 | minor | code read (B_Adepth vs B_m7001; cd in a subshell left of ||) | yes | 6ca8c51 |
| F27 | major | grep: four modules import oceanmesh; THIRD_PARTY_NOTICES endorses it, against CLAUDE.md | yes, pre-existing | **owner decision needed** (subprocess vs separate GPL plugin) |
| F28 | minor | tests (../, absolute, a/b) | yes | 3e61359 |
| F29 | minor | test (3x3 grid, [0, 2, 1]) | yes | 6ca8c51 |

Tests after the fixes: 1096 passed.


## Round 2

### Answer (verbatim)

Several fixes remain partial, and the fixes introduce additional defects. No files were changed; the pre-existing untracked PNG remains untouched. No production data or batch jobs were used.

Verification: 22 pure geometry tests passed with `--noconftest`, bytecode/cache writes disabled, and Numba JIT disabled. Scoped Ruff checks and all OCTOPUS shell syntax checks passed. Ordinary pytest collection failed because its dependency-import hook required a writable Matplotlib cache. Additional checks used small arrays and in-memory mocks. Production generation, production QA and FVCOM integration were not run.

**Status of previous findings**

References below point to the numbered findings in this review.

| Previous | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked grid values become NaN; an in-memory priority-stack test selected the fallback correctly. |
| F2 | RESOLVED | Failed final gates now cause a nonzero exit. The resulting provenance regression is finding 18. |
| F3 | RESOLVED | Bounds participate in smoothing; nonconvergence raises; callers recompute the final r-factor. Additional limiter defects are findings 14–15. |
| F4 | RESOLVED | Generation checks the OBC band, and 447 gates final new-element allowance against the base. The reported Enshu failure is not counted as a defect. |
| F5 | PARTIAL | Constant-spacing terminal remainder is fixed; variable spacing and chord floors remain incorrect—finding 4. |
| F6 | PARTIAL | Missing wet lattice and design-point samples are rejected; generated-node samples still receive fabricated depths—finding 19. |
| F7 | PARTIAL | Full-line validation catches crossings, but normal selection still skips thin land—finding 5. The new validation also has a numerical regression—finding 6. |
| F8 | PARTIAL | Rejected designs preserve the CSV, but concurrent publication can report success with another writer’s CSV—finding 2. |
| F9 | PARTIAL | 447 detects frozen-coordinate changes after export; serialization remains lossy, and 453 does not check coordinates—finding 8. |
| F10 | PARTIAL | Same-side seam incidence is rejected; overlap away from the interface remains unchecked—finding 9. |
| F11 | RESOLVED | One-edge land runs are retained; regression test passed. |
| F12 | RESOLVED | Unsupported chains now receive an explicit six-node requirement before mesh generation. |
| F13 | RESOLVED | Cache keys include resolved roots and dimensions; mocked distinct roots returned distinct depths. |
| F14 | RESOLVED | Only finite fine-grid interpolations replace coarse values. |
| F15 | RESOLVED | Distinct-point count and geometric rank are checked before triangulation. |
| F16 | RESOLVED | Provenance inventory now includes CAO depth files, as well as tables. |
| F17 | PARTIAL | Ordinary dirty paths trigger hashing; incoming renames can still be misidentified as committed code—finding 10. |
| F18 | PARTIAL | Settings and depth controls are validated; geographic bounds remain unchecked—finding 11. |
| F19 | PARTIAL | NaN spacing is rejected, but finite spacing can still cause nontermination—finding 12. |
| F20 | RESOLVED | 445 uses an exclusive `.reserved` creation. |
| F21 | PARTIAL | Existing output and source aliases are rejected, but the output is not reserved atomically—finding 1. |
| F22 | PARTIAL | Sequential reuse is rejected; concurrent staging remains possible—finding 1. Refusing reuse also destroys the old success marker—finding 21. |
| F23 | RESOLVED | Smoke depth control now precedes timestep selection. |
| F24 | RESOLVED | Smoke case and root paths are resolved before staging. |
| F25 | RESOLVED | The smoke job honors `FMESH_FVCOM` and logs its hash. |
| F26 | RESOLVED | Case names match preparation; directory changes are explicitly guarded. |
| F27 | NOT RESOLVED | The acknowledged GPL-import policy conflict remains. Finding 22 records it as previously known, not new. |
| F28 | RESOLVED | Recipe and re-depth case names are restricted to filename components. |
| F29 | PARTIAL | Scrambled walks are rejected, but a valid reversed all-open walk is now rejected—finding 13. |

**Findings**

1. **Major — Re-depth and smoke output directories still permit concurrent writers.**  
   [453_redepth_extended.py:54](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:54), [448_extend_smoke.py:59](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/448_extend_smoke.py:59). Both check absence/emptiness, perform intervening work, then create directories with `exist_ok=True`. Executing each actual guard twice against a mocked absent destination allowed both callers through. Concurrent re-depth jobs can overwrite cases/reports; concurrent smoke jobs can mix staged inputs and history. **Fix:** acquire an exclusive reservation before processing inputs, following 445’s pattern, and retain ownership through the run. Residual F21/F22.

2. **Major — The new atomic boundary publication uses shared temporary filenames.**  
   [444_design_obc.py:170](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:170). A deterministic two-thread reproduction of the actual publication loop produced: writer A succeeded, writer B raised `FileNotFoundError`, the published CSV contained B’s coordinates, and the JSON described A. Both writers use `<destination>.tmp`, so A can rename B’s payload. **Fix:** use unique temporary files and exclusive publication ownership; bind the report to the CSV with a content hash. Introduced by the F8 fix.

3. **Major — Re-depthing bypasses the extension’s final acceptance checks.**  
   [453_redepth_extended.py:95](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:95). The script neither runs final QA nor checks the new-element timestep allowance. In a mocked six-node depth/export/report execution, changing the two new depths from 10 to 15 m satisfied `r=0.2` and wrote `redepth.json`, while allowance fell from **100.964 s** in the base to **82.437 s** in the extension. Its input check also omits base connectivity. **Fix:** verify the complete frozen-base contract and repeat the applicable final acceptance checks after re-depthing. Deliberately rejected sensitivity cases should be explicitly identified as such.

4. **Minor — Resampling still violates its spacing-floor contract.**  
   [obc_design.py:167](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:167), [444_design_obc.py:152](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:152). On a 10,900 m straight line, spacing `1000` before x=5200 and `2500` thereafter produced an edge from x=4360 to x=5450: **1090 m against a 2500 m floor**. Stretching positions changes where the spacing function is evaluated. Separately, a quarter-circle reproduction returned chord/floor **0.998972**, which 444’s `0.995` threshold accepts despite the documented “never below” rule. **Fix:** re-evaluate spacing after redistribution and solve against actual chord lengths, lengthening/coarsening until the floor holds. Residual F5; redistribution is a fix regression.

5. **Minor — Coast-normal selection still jumps over thin land.**  
   [obc_design.py:59](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:59). The first probe is 60 m away at the default chord. For `box(-50000, 0, 50000, 50)` and a query south of it, reversing polygon winding changed the returned normal from **180° to 0°**, crossing the entire strip. Full-line validation subsequently rejects the design instead of selecting the valid direction. **Fix:** inspect the departure segment/local polygon interior continuously, rather than 25 separated points. Residual F7.

6. **Minor — The new land-crossing gate rejects harmless endpoint roundoff.**  
   [444_design_obc.py:150](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:150). A rectangular coast rotated 13° and translated to metric coordinates `(400000, 3900000)` produced a valid outward normal with an intersection length of **5.82×10⁻¹¹ m** at its snapped endpoint. Testing `crosses_land_m > 0` rejects it. **Fix:** distinguish endpoint contact/numerical residue from an actual crossing using a documented geometric tolerance. Introduced by the F7 fix.

7. **Minor — The published boundary is not the geometry that was validated.**  
   [444_design_obc.py:170](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:170). Validation uses full-precision metric nodes, but CSV coordinates are rounded to six decimal places. A synthetic coast endpoint at lon/lat `(139.00000049, 35.00000049)` passed before serialization; its serialized/reprojected departure crossed **0.05345 m** of land. **Fix:** serialize with adequate precision and validate the reloaded coordinates before publication.

8. **Minor — Native export still cannot preserve arbitrary valid base coordinates and depths.**  
   [fvcom_native.py:45](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:45), [453_redepth_extended.py:97](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:97). An intercepted native write changed `1000.123456789` to `1000.12345679`. The new 447 check detects this but consequently rejects otherwise valid high-precision bases; 453 checks only depths and can silently change coordinates. Depths likewise retain six-decimal serialization. **Fix:** use round-trip-safe float serialization and verify coordinates, connectivity and depths after export in both paths. Residual F9.

9. **Minor — Opposite-side incidence does not establish a nonoverlapping extension.**  
   [extend.py:165](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:165). A connected synthetic outer mesh shared the interface correctly but wrapped around an endpoint and overlapped **100,000 m²** of the base elsewhere. `verify_frozen_base` returned success. This reproduction establishes the verifier’s gap, not a production QA pass. **Fix:** additionally check outer-element intersection with the base footprint, allowing only the prescribed shared boundary. Residual F10.

10. **Minor — Dirty-code detection misses renames into the inspected source tree.**  
    [provenance.py:132](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:132). With porcelain entry `other/new.py -> src/fvcom_mesh_tools/new.py`, the function treats the entire rename expression as one path. A mocked Git state returned `commit_identifies_code=True` with no source hash. **Fix:** parse NUL-delimited porcelain records correctly, including both rename paths, before deciding whether the commit identifies the sources. Residual F17.

11. **Minor — Recipe geographic bounds remain unvalidated.**  
    [extend_recipe.py:75](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend_recipe.py:75). Mocked file preconditions allowed both `[NaN, 33, 141, 36]` and reversed bounds `[141, 36, 137, 33]` through `load_extend_recipe`. Only the number of entries is checked. **Fix:** require finite numeric coordinates, valid geographic ranges, and strictly ordered minima/maxima before reserving outputs or reading land. Residual F18.

12. **Minor — Finite spacing can still make resampling loop indefinitely.**  
    [obc_design.py:162](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:162). On a 100,000,000 m line, a callback returning `99,999,999` initially and `1e-12` subsequently repeatedly queried x=`99,999,999`: floating-point addition no longer advanced. A bounded callback stopped the reproduction after 16 calls. **Fix:** require each proposed position to be finite and strictly greater than the preceding position. Residual F19.

13. **Minor — The stricter native reader rejects a valid reverse traversal.**  
    [fvcom_native.py:564](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:564). With all four square boundary nodes open, mocked reads accepted `[0,1,2,3]` but rejected `[0,3,2,1]`. Both neighbours of the first node belong to the open-node set, so the membership-based direction heuristic selects the wrong traversal. **Fix:** choose direction from the supplied second node, then validate consecutive adjacency and uniqueness. Introduced by the F29 fix.

14. **Minor — The limiter rejects convergence on its last allowed iteration.**  
    [extend.py:267](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:267). With depths `[10,100]`, both nodes free, `rmax=.2`, `hmin=3`, and `max_iter=1`, the correction reaches `[44,66]`, exactly `r=.2`; the function nevertheless raises “limit … not reached”. **Fix:** check the recomputed final r-factor before raising on exhaustion. Introduced by the F3 fix.

15. **Minor — Nonfinite depths can be returned as successfully limited.**  
    [extend.py:249](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:249). Calling the public limiter with depths `[10, NaN]`, one fixed and one free node, returned `(array([10, NaN]), 0, NaN)` without error because NaN fails the `r > limit` comparison. Normal pipeline sampling guards this case, but the limiter itself does not. **Fix:** validate finite positive depths and valid bounds/controls, and explicitly reject nonfinite computed ratios.

16. **Minor — Ladder direction depends on an undocumented input orientation.**  
    [446_extend_generate.py:164](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:164). The interface is always reversed and the new boundary never reversed before constructing left-side ladders. Reversing a valid six-node boundary at latitude 35 changed the guide latitude from **35.011261** to **34.988739**. The loader accepts either order, while the ladder filter checks land rather than membership in the extension domain; guides can therefore be outside the domain or all discarded. **Fix:** determine the extension-facing side geometrically and orient each ladder accordingly.

17. **Minor — Later sizing bands silently override earlier constrained-band targets.**  
    [extend.py:93](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:93). On a 5×5 lattice, a 1000 m band at y=0 followed by a 5000 m band at y=1000, with gradation `.2`, returned **4800 m throughout the first band**. The report contains no conflict indication. These constraints are incompatible with the promised exact band sizes. **Fix:** check joint feasibility of all band constraints and reject conflicts, or explicitly document/report a supported priority policy.

18. **Minor — Failed builds and re-depth variants lack sufficient provenance.**  
    [445_extend_mesh.py:68](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:68), [453_redepth_extended.py:102](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:102). Now that 447 exits nonzero on rejection, `check=True` prevents 445 from reaching its provenance collection and `report.json`. Re-depth reports separately record source names but omit input/source/code hashes and the effective depth controls. **Fix:** capture input identity and effective settings before processing, and write a status-bearing report on both success and failure; apply the same provenance scheme to re-depth variants. The failed-build portion is introduced by the F2 fix.

19. **Minor — Generated nodes without sizing bathymetry still receive an invented 2 m depth.**  
    [446_extend_generate.py:284](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:284). The new coverage check concerns lattice points only. The later sample at actual generated nodes still runs through `nan_to_num(dn, nan=2.0)` and feeds those artificial depths into finishing, whose operations use depth-dependent timestep constraints. Coverage at lattice points does not prove coverage at every generated node. **Fix:** reject uncovered generated wet nodes, or use an explicit conservative fallback with reported coverage. Residual F6.

20. **Minor — Fillet radii are not validated.**  
    [obc_design.py:101](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:101). A radius of `-1000` for `[[0,0],[10000,0],[10000,10000]]` returned a path extending to x=12000 and y=−2000 instead of refusing an invalid radius. Radius zero silently retains a sharp corner, contrary to the rounded-corner contract. **Fix:** require finite positive radii and valid, nonzero adjacent segments; validate arc sampling parameters as well.

21. **Minor — Refusing smoke-directory reuse destroys the previous success marker.**  
    [448_extend_smoke.sh:20](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/448_extend_smoke.sh:20). The shell removes `SMOKE_OK` before Python checks whether the root is already populated. Resubmitting a completed run now deletes its success marker and then refuses to stage anything. **Fix:** reserve/validate the destination before touching markers, and leave an existing completed run unchanged when rejecting reuse. This interaction was introduced by the F22 guard.

22. **Major — Previously known F27 remains open: package code directly imports oceanmesh.**  
    [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332). Direct imports also remain in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672`, and `autofinish/directives.py:20`; `THIRD_PARTY_NOTICES.md:37` still endorses that arrangement despite the repository’s explicit policy. **Fix:** implement the owner-selected subprocess or separately licensed plugin architecture and reconcile the notices. **Previously reported, not a new finding; included once in the open-defect verdict.**

## Verdict

VERDICT: FAIL (0 blocker, 4 major, 18 minor, 0 nit)

### Prompt

```markdown
# Review request, round 2: extending a base mesh outward (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (small synthetic inputs only; do not read the large data
under $DATA_DIR beyond listing it, and do not submit batch jobs). Answer in
English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
A tool that keeps a finished FVCOM base mesh exactly as it is and adds the
sea out to a new, designed open boundary (USER_GUIDE section 13). Read:

- `git show b0584f9 8e2739b 450ad44 d6d2a72 7044b6b 0f52d5b 69b50a4 d9e92fd b6d2ed8 765423c`
  (the extension tool and its documentation), and the current files:
  - `src/fvcom_mesh_tools/extend.py`, `extend_recipe.py`, `obc_design.py`,
    `dem/sources.py` (named bathymetry sources, priority stack, and the new
    `DATUM` registry / `non_tp_count` warning);
  - `notebooks/444_design_obc.py`, `445_extend_mesh.py`, `446_extend_generate.py`,
    `447_extend_merge.py`, `448_extend_smoke.py`, `453_redepth_extended.py`;
  - `recipes/extend/tokyo_bay_enshu.yaml`, `tokyo_bay_enshu_obc_design.yaml`;
  - `jobs/octopus/444_design_obc.sh`, `445_extend_mesh.sh`, `448_extend_smoke.sh`,
    `453_redepth_extended.sh`, `common.sh`;
  - tests: `tests/test_extend*.py`, `tests/test_obc_design*.py`,
    `tests/test_dem_sources.py` (whatever exists).
- Also in scope, just committed: the portability change --
  every job script and `common.sh` now take paths only from `$DATA_DIR` and
  `$WORK_DIR` (login profile), stop when they are unset, and derive the
  OCTOPUS FVCOM library directory as `FVCOM_LIBS` in `common.sh`; notebooks
  383/384/414 and `cli/refine_run.py` no longer fall back to `/octfs/...`.
  See commits 6d8b9a7 and 6c068d2 (`git log -5`).

Design intent:
- the base mesh's nodes, elements and depths are carried bit for bit
  (`verify_frozen_base`);
- the new part is generated with oceanmesh (subprocess; GPL code must never
  be imported by `fvcom_mesh_tools`), with fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Round 1 found 29 (13 major, 16 minor); all verified correct. Fixes:
662b142 (F1, F13-F16), 960b232 (F3, F4 report, F10, F11), 65119e6 (F4 446
gate, F5-F8, F12, F19), 3e61359 (F2, F4 447 gate, F9, F18, F20-F25, F28),
6ca8c51 (F17, F26, F29). Read `git log ddeb8bf -8` and the diffs. Notes:
- F4: the time-step gate in 447 compares the new elements' allowance with
  the base's; the CURRENT Enshu mesh fails it (9.7 s vs 14.4 s). That is a
  mesh/recipe matter to be decided by the owner, not a code defect left open.
- F27 (oceanmesh imported by mesh_engine/oceanmesh.py, multiscale.py,
  mesh_clean.py, autofinish/directives.py; THIRD_PARTY_NOTICES endorses it):
  confirmed; the remedy (subprocess boundary vs a separate GPL plugin
  package) is an architectural/licensing decision put to the owner; not yet
  changed. Report it as still open; do not count it as a new finding.
- The round-1 record is docs/extend-tools-review-20261001.md.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the scope and everything it touches.

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

| id | sev | verified? (how) | correct? | action |
|---|---|---|---|---|
| R2-1 | major | reproduced the check-then-create gap; outdir.reserve | yes | b5d0994 |
| R2-2 | major | code read (shared .tmp) | yes | b5d0994 |
| R2-3 | major | code read (no QA/dt gate in 453) | yes | b5d0994 |
| R2-4 | minor | test (10.9 km, 1000/2500 m; quarter circle) | yes | b5d0994 |
| R2-5 | minor | test (50 m strip) | yes | b5d0994 |
| R2-6 | minor | code read (> 0 test) | yes | b5d0994 |
| R2-7 | minor | code read (6 decimals published) | yes | b5d0994 |
| R2-8 | minor | test (1000.123456789) | yes | b5d0994 |
| R2-9 | minor | test (overlap away from the seam) | yes | b5d0994 |
| R2-10 | minor | test (porcelain rename) | yes | b5d0994 |
| R2-11 | minor | tests (NaN, reversed, lat -95) | yes | b5d0994 |
| R2-12 | minor | code read (non-advancing step) | yes | b5d0994 |
| R2-13 | minor | test (four open nodes, reversed) | yes | b5d0994 |
| R2-14 | minor | test ([10, 100], max_iter 1) | yes | b5d0994 |
| R2-15 | minor | tests (NaN, rmax 1.5) | yes | b5d0994 |
| R2-16 | minor | code read (fixed flip flags) | yes | b5d0994 |
| R2-17 | minor | test (two bands 1 km apart) | yes | b5d0994 |
| R2-18 | minor | code read (check=True; no 453 provenance) | yes | b5d0994 |
| R2-19 | minor | code read (nan -> 2 m) | yes | b5d0994 |
| R2-20 | minor | tests (-1000, 0, NaN) | yes | b5d0994 |
| R2-21 | minor | code read (rm before the check) | yes | b5d0994 |
| R2-22 (F27) | major | as round 1 | yes, pre-existing | owner decision pending |

Tests after the fixes: 1110 passed. An integration run (444 design, 445
build, 453 re-depth into scratch) checks the notebook changes on real data.


## Round 3

### Answer (verbatim)

Read-only review completed. No files changed; the existing untracked PNG was left untouched. No large datasets, production meshes or batch jobs were used.

Verification:

- Pure geometry tests: **28 passed**, using `pytest --noconftest -s -q -p no:cacheprovider tests/test_extend.py tests/test_obc_design.py` with bytecode writing disabled.
- Scoped `ruff check --no-cache`: passed.
- `bash -n` on every `jobs/octopus/*.sh`: passed.
- Additional reproductions used small arrays, in-memory streams and mocked notebook execution.
- Initial pytest attempts failed before running tests: output capture required a writable temporary directory; subsequently, the optional-engine hook imported Matplotlib, which required a writable cache. Disabling capture and that hook allowed the pure tests to run.

The documented real-data bathymetry/floor conflict is **not** counted as a defect below.

### Previous findings: round 1

“RESOLVED” refers to the reported defect; separate regressions are identified in the findings.

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked elevations become NaN rather than valid fill-value depths. A separate interpolation regression remains, finding 7. |
| F2 | RESOLVED | Failed required QA gates produce a nonzero build exit. |
| F3 | RESOLVED | Bounds are applied during limiting; nonconvergence fails; final exported depths are checked again. |
| F4 | PARTIAL | Raw band conflicts and final new-element timestep allowance are checked, but tolerated smoothing can bypass the band-floor gate: finding 1. |
| F5 | PARTIAL | Resampling checks actual chords and both endpoints; publication still accepts below-floor edges: finding 2. |
| F6 | RESOLVED | Missing sizing coverage is rejected at boundary samples, wet lattice points and generated nodes. |
| F7 | RESOLVED | Continuous departure segments and the complete boundary are checked against land. |
| F8 | RESOLVED | Rejected geometric designs leave the existing CSV intact. Publication has a separate failure mode, finding 5. |
| F9 | RESOLVED | Native serialization round-trips doubles; both pipelines verify the exported frozen base. The supplementary fort.14 writer has a separate defect, finding 8. |
| F10 | RESOLVED | Opposite-side seam incidence and outer-versus-base overlap are checked. |
| F11 | RESOLVED | Single-edge land runs are retained. |
| F12 | RESOLVED | Generation explicitly rejects constrained lines with fewer than six nodes. |
| F13 | RESOLVED | CAO cache keys include the resolved data root and dimensions. |
| F14 | RESOLVED | Invalid fine-grid interpolation no longer overwrites valid coarse coverage. |
| F15 | RESOLVED | Insufficient distinct points and collinear windows return uncovered samples. |
| F16 | RESOLVED | CAO provenance now includes depth files as well as area tables. |
| F17 | RESOLVED | Dirty relevant sources are hashed; rename paths are parsed correctly. |
| F18 | RESOLVED | Extension settings, seeds, depth controls and geographic bounds receive the reported validation. |
| F19 | RESOLVED | Nonfinite spacing and nonadvancing steps are rejected. |
| F20 | RESOLVED | Extension output reservation uses exclusive marker creation. |
| F21 | RESOLVED | Re-depth refuses source/output overlap and reserves its destination atomically. |
| F22 | RESOLVED | Smoke staging atomically reserves a fresh root, excluding previous history. |
| F23 | RESOLVED | Notebook 448 applies depth control before selecting its timestep. Notebook 414 retains the analogous problem, finding 10. |
| F24 | RESOLVED | Smoke case and root paths are resolved before use. |
| F25 | RESOLVED | The smoke job respects and logs `FMESH_FVCOM`. |
| F26 | RESOLVED | The 383 job uses `B_m7001` and explicitly guards directory changes. |
| F27 | NOT RESOLVED | Direct package imports remain; owner decision pending. Finding 14, counted once. |
| F28 | RESOLVED | Extension and re-depth case names are restricted to filename components. |
| F29 | RESOLVED | Native reading requires a unique consecutive boundary walk. |

### Previous findings: round 2

| ID | Status | Reason |
|---|---|---|
| R2-1 | RESOLVED | Both re-depth and smoke use exclusive directory reservation. |
| R2-2 | RESOLVED | Unique temporaries and a publication lock remove the reported shared-temporary collision; the report includes the CSV hash. Separate publication fault: finding 5. |
| R2-3 | RESOLVED | Re-depth verifies connectivity and the exported frozen base, checks overlap, runs QA and compares timestep allowances. Overrides are explicitly labelled. |
| R2-4 | PARTIAL | The resampler’s chord/endpoint defects are fixed; the permissive publication gate remains: finding 2. |
| R2-5 | RESOLVED | Continuous segment testing catches the 50 m land strip. |
| R2-6 | RESOLVED | The land-crossing check allows 1 mm of numerical residue. |
| R2-7 | RESOLVED | Validation uses coordinates reconstructed from the published nine-decimal values. |
| R2-8 | RESOLVED | Native coordinates and depths use round-trip-safe serialization. |
| R2-9 | RESOLVED | Positive-area overlap with the base is checked beyond the interface. |
| R2-10 | RESOLVED | NUL-separated porcelain parsing includes both rename paths. |
| R2-11 | RESOLVED | Geographic bounds must be finite, ordered and within the accepted ranges. |
| R2-12 | RESOLVED | Every resampling step must advance finitely. |
| R2-13 | RESOLVED | Boundary traversal direction follows the supplied second node. |
| R2-14 | RESOLVED | Convergence on the final permitted iteration is accepted. |
| R2-15 | RESOLVED | Nonfinite depths on limited edges and invalid controls are rejected. |
| R2-16 | PARTIAL | Orientation is inferred geometrically, but an unsuccessful left probe is treated as evidence for the right: finding 4. |
| R2-17 | PARTIAL | Lists receive the deviation check, with the explicitly accepted 5% policy; iterators bypass it: finding 3. |
| R2-18 | PARTIAL | Ordinary stage failures and completed variants receive provenance, but several early failures still do not: finding 6. |
| R2-19 | RESOLVED | Uncovered generated nodes are rejected. |
| R2-20 | RESOLVED | Nonpositive/nonfinite radii are rejected. A distinct reversal defect remains, finding 13. |
| R2-21 | RESOLVED | Rejecting smoke-root reuse no longer removes `SMOKE_OK`. |
| R2-22 | NOT RESOLVED | Same outstanding F27; finding 14. |

### Findings

1. **Major — Tolerated band smoothing can bypass the timestep-floor gate.**  
   [extend.py:94](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:94), [446_extend_generate.py:131](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:131).

   **Evidence:** On a 5×5 lattice spaced 100 m apart, with gradation `0.2`, an interface band of 1,000 m and an open-boundary band `[1000,1040,1040,1040,1040]`, setting the floor equal to the latter targets produced `[1000,1020,1040,1040,1040]`. The deviation was accepted at **1.923%**, and `band_1_below_floor_cells` remained **0**. Generation therefore accepts a field below its floor. The counter examines the original targets, before smoothing. The later base-relative timestep check does not enforce the recipe’s specified floor.

   **Fix:** Recompute band-floor violations against the final composed field, after all bands. Apply the intended interface exception explicitly. This regression is exposed by `6360e8d`; it is independent of the documented deep-water design conflict.

2. **Minor — Boundary publication still accepts spacing below the floor.**  
   [444_design_obc.py:161](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:161).

   **Evidence:** The acceptance threshold remains `0.995`, although resampling now handles chord shortening. A UTM54 segment from `(400000,3900000)` to `(410000,3900000)`, with floor 3,001 m for `402000 < x < 402999.99998` and 3,000 m elsewhere, resampled successfully. Nine-decimal geographic publication moved an endpoint across that floor transition. Rechecking the published geometry gave edge/floor **0.99966675**, which 444 accepts.

   **Fix:** Enforce the floor after serialization, allowing only a documented numerical distance tolerance. Recompute/coarsen the boundary when publication changes its required spacing. Residual F5/R2-4.

3. **Minor — Iterator-valued bands bypass the new conflict check.**  
   [extend.py:106](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:106).

   **Evidence:** `bands` is traversed twice. Passing `iter([b0,b1])` exhausts it during composition, so validation never executes. The round-2 incompatible-band example—1,000 m and 5,000 m bands separated by 1,000 m—raises for a list but returns **4,800 m throughout the first band** for an iterator, without deviation diagnostics.

   **Fix:** Materialize and validate the bands once before either traversal. Introduced by the R2-17 fix.

4. **Minor — Ladder orientation assumes the right side is valid without checking it.**  
   [446_extend_generate.py:180](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:180).

   **Evidence:** Executing the actual selector with a six-node line, 1,000 m spacing and a 100 m-wide extension immediately on its left returned `False`: every 250 m left probe overshot the domain. The right side was never tested. The subsequent ladder filter checks land, not domain membership, so wrong-side guides can survive.

   **Fix:** Test both sides at adaptive local distances, reject ambiguous/unresolvable geometry, and check actual guide points and edges against the extension domain. Residual R2-16.

5. **Minor — Boundary publication can leave mismatched products after failure.**  
   [444_design_obc.py:190](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:190).

   **Evidence:** An in-memory execution of the actual publication block, failing the second `os.replace`, left **new CSV + old JSON**, then released the lock. The loader does not verify the report’s CSV hash. The PNG is also published outside the lock, allowing concurrent designs to leave a figure from another publication.

   **Fix:** Publish a complete versioned artifact set through one atomic pointer/directory switch, and have consumers verify its identity. Include the figure in the publication protocol and clean abandoned temporaries. Additional failure mode in the R2-2 implementation.

6. **Minor — Failed-run provenance remains incomplete and can itself abort reporting.**  
   [445_extend_mesh.py:90](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:90), [453_redepth_extended.py:81](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:81).

   **Evidence:** Mocking a failed 446 subprocess and unavailable oceanmesh caused 445 to pass `None` into `code_state`, raising `TypeError` with **zero report writes**. Mocked execution of 453 with two uncovered new nodes likewise exited with **zero report writes**. Missing input data, limiter failures and malformed stage JSON also precede report publication. Input/code identity is still collected after processing rather than captured before execution.

   **Fix:** Write an initial status/provenance record before processing; represent unavailable dependencies explicitly; finalize failure information through guarded exception handling. Snapshot or detect changes to execution inputs. Residual R2-18.

7. **Minor — The masked-grid fix discards valid samples with zero-weight masked neighbours.**  
   [dem/sources.py:90](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:90).

   **Evidence:** An in-memory 2×2 elevation grid `[[−10,−20],[−30,masked]]`, queried exactly at the valid `−10` vertex, returned **NaN**, not depth **10**. The interpolator propagates the neighbouring NaN even though its interpolation weight is zero. The priority stack consequently substitutes a fallback or reports missing coverage.

   **Fix:** Ignore zero-weight contributors while requiring every positive-weight contributor to be valid; preserve exact valid vertex/edge samples. Regression from the F1 masking fix.

8. **Minor — The supplementary fort.14 export still violates its round-trip contract.**  
   [fort14.py:182](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:182).

   **Evidence:** Actual writer/reader execution through in-memory streams changed coordinate `0.12345678912345678` to `0.123456789123457`. Coordinates still use `.15f`. Notebook 447 exports this supplementary product but performs its frozen-base read-back only on the native case.

   **Fix:** Use round-trip-safe coordinate formatting in fort.14 too, and verify the advertised frozen contract for that product. Pre-existing; native serialization fixes do not cover it.

9. **Minor — QA can approve a self-overlapping triangulation.**  
   [qa.py:1110](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:1110).

   **Evidence:** Construct eight nodes: a centre and seven radius-1,000 m vertices ordered at angles `4πk/7`; connect the centre to consecutive vertices. With positive depths, metric coordinates and no open boundary, this connected mesh passed **16/16 applicable gates**. Its summed element area was **3,412,247.69 m²**, versus union area **2,101,798.05 m²**. The new extension overlap check only compares outer elements with the base.

   **Fix:** Add geometric self-intersection/positive-area element-overlap checks, permitting legitimate shared edges and vertices. This reproduction establishes a QA gap, not a production extension failure.

10. **Minor — Notebook 414 still selects its timestep before depth control.**  
    [414_refine_m2_prep.py:125](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/414_refine_m2_prep.py:125).

    **Evidence:** A nine-node synthetic mesh gave an initial allowance of **15.14456 s**, selecting **15 s**. Applying the actual OBC depth-control function reduced the allowance to **4.78913 s**. Unlike corrected notebook 448, 414 applies that transformation later, at line 196. Its subsequent `mesh_metrics` check rejects the automatically selected timestep after case files have already been written.

    **Fix:** Transform both meshes before computing their common allowance, selecting the timestep or writing products. Pre-existing in the explicitly scoped portability notebook.

11. **Minor — Design recipes accept negative CFL controls and silently remove their constraint.**  
    [444_design_obc.py:105](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:105).

    **Evidence:** The actual `spacing()` function, with depth 4,000 m, `cfl_dt_s=-18`, `cfl_cr=0.9` and `min_m=3000`, returned **3,000 m** without rejection. A positive 18 s control requires **3,961.82 m**. The publication check reuses the same invalid controls, so it cannot detect this error.

    **Fix:** Validate design-recipe numeric controls as finite and positive before geometry or data processing, matching extension-recipe validation.

12. **Minor — Short boundary designs crash while constructing their report.**  
    [444_design_obc.py:146](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:146).

    **Evidence:** Resampling an 8,000 m line at 3,000 m produces three valid nodes. For that result, `nodes[2:-2]` is empty and the report raises `ValueError: min() iterable argument is empty`. Two- and four-node results have the same issue.

    **Fix:** Reject unsupported boundary lengths explicitly before processing, or represent the unavailable interior-distance statistic as null. Do not let optional reporting determine geometric validity.

13. **Minor — Filleting silently preserves a sharp reversal.**  
    [obc_design.py:105](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:105).

    **Evidence:** `fillet([[0,0],[10000,0],[0,0]], [1000])` returns the original backtracking polyline unchanged. The same branch handles straight continuation and a 180° reversal, despite their different geometry and the rounded-corner contract.

    **Fix:** Treat straight continuation separately; reject reversals and overlapping/backtracking paths that cannot support the requested tangent fillet. Pre-existing, separate from the corrected radius validation.

14. **Major — Previously reported F27 remains open: package code directly imports oceanmesh.**  
    [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

    **Evidence:** Direct imports also remain in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672` and `autofinish/directives.py:20`. `THIRD_PARTY_NOTICES.md:37` still endorses the arrangement contrary to the repository’s explicit package-import policy.

    **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary and reconcile the notices. **Previously reported, not new; counted once.**

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 12 minor, 0 nit)

### Prompt

```markdown
# Review request, round 3: extending a base mesh outward (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (small synthetic inputs only; do not read the large data
under $DATA_DIR beyond listing it, and do not submit batch jobs). Answer in
English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
A tool that keeps a finished FVCOM base mesh exactly as it is and adds the
sea out to a new, designed open boundary (USER_GUIDE section 13). Read:

- `git show b0584f9 8e2739b 450ad44 d6d2a72 7044b6b 0f52d5b 69b50a4 d9e92fd b6d2ed8 765423c`
  (the extension tool and its documentation), and the current files:
  - `src/fvcom_mesh_tools/extend.py`, `extend_recipe.py`, `obc_design.py`,
    `dem/sources.py` (named bathymetry sources, priority stack, and the new
    `DATUM` registry / `non_tp_count` warning);
  - `notebooks/444_design_obc.py`, `445_extend_mesh.py`, `446_extend_generate.py`,
    `447_extend_merge.py`, `448_extend_smoke.py`, `453_redepth_extended.py`;
  - `recipes/extend/tokyo_bay_enshu.yaml`, `tokyo_bay_enshu_obc_design.yaml`;
  - `jobs/octopus/444_design_obc.sh`, `445_extend_mesh.sh`, `448_extend_smoke.sh`,
    `453_redepth_extended.sh`, `common.sh`;
  - tests: `tests/test_extend*.py`, `tests/test_obc_design*.py`,
    `tests/test_dem_sources.py` (whatever exists).
- Also in scope, just committed: the portability change --
  every job script and `common.sh` now take paths only from `$DATA_DIR` and
  `$WORK_DIR` (login profile), stop when they are unset, and derive the
  OCTOPUS FVCOM library directory as `FVCOM_LIBS` in `common.sh`; notebooks
  383/384/414 and `cli/refine_run.py` no longer fall back to `/octfs/...`.
  See commits 6d8b9a7 and 6c068d2 (`git log -5`).

Design intent:
- the base mesh's nodes, elements and depths are carried bit for bit
  (`verify_frozen_base`);
- the new part is generated with oceanmesh (subprocess; GPL code must never
  be imported by `fvcom_mesh_tools`), with fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Round 1: 29 findings, all fixed except F27. Round 2: 22 findings (21 new +
F27), all verified; fixed in b5d0994 and 6360e8d (BAND_TOLERANCE 5 %: the
real interface band departs 0.9 % from its sizes because the base spacing
varies faster than the gradation allows -- smoothing, not a band conflict).
Read `git log a2afa1c..HEAD` and `git show b5d0994 6360e8d`, and the record
docs/extend-tools-review-20261001.md.
Notes:
- Integration run on real data (2026-10-01): 444 redesign passed (191 nodes,
  every edge >= 3000 m, ends 90.0 deg); 453 re-depth passed with
  --allow-failing-gates; 445 now stops in 446 because 11,388 open-boundary
  band cells are below the time-step floor computed from the raw depths.
  That is the same design issue as round-1 F4 (deep water near the
  boundary). The owner will cap the maximum depth and smooth the slopes in
  the next stage; the sizing floor will then use the capped depth. Not a
  code defect to report.
- F27 (oceanmesh imported in package code): owner decision pending; report
  as still open, not as new.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the scope and everything it touches.

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

| id | sev | verified? (how) | correct? | action |
|---|---|---|---|---|
| R3-1 | major | test (5x5 lattice, smoothed band) | yes | fixed; floor findings now warnings (owner), 7007e7e |
| R3-2 | minor | code read (0.995) | yes | fixed, 7007e7e |
| R3-3 | minor | test (iterator) | yes | fixed, 7007e7e |
| R3-4 | minor | code read (right side untested) | yes | fixed, 7007e7e |
| R3-5 | minor | code read (order, PNG outside lock) | yes | fixed, 7007e7e |
| R3-6 | minor | code read (None into code_state; early exits) | yes | fixed, 7007e7e |
| R3-7 | minor | test (2x2 grid) | yes | fixed, 7007e7e |
| R3-8 | minor | code read (.15f) | yes | fixed, 7007e7e |
| R3-9 | minor | test (fan winding twice) | yes | fixed (new QA gate), 7007e7e |
| R3-10 | minor | code read (414 order) | yes | fixed, 7007e7e |
| R3-11 | minor | code read | yes | fixed, 7007e7e |
| R3-12 | minor | code read (nodes[2:-2] empty) | yes | fixed, 7007e7e |
| R3-13 | minor | test (reversal) | yes | fixed, 7007e7e |
| R3-14 (F27) | major | as round 1 | yes, pre-existing | owner decision pending |

Owner decision (2026-10-01): meshes are made from the real depths; the
time step is settled in the depth stage (maximum depth, smoothing). The
band-floor check (446) and the new-element time-step gate (447, 453)
therefore report warnings and the values, and do not fail the build.
Tests after the fixes: 1115 passed.

## Round 4

The first attempt (gpt-6-astra) stopped at the Codex usage limit before its
answer. On the owner's instruction (2026-10-01) the round was run again with
gpt-6.1-sol, with the same prompt.

### Answer (verbatim)

Read-only review completed. No files changed; the pre-existing untracked PNG remains untouched. No shared datasets or batch jobs were used.

Verification:

- `pytest --noconftest -s -q -p no:cacheprovider tests/test_extend.py tests/test_obc_design.py`: **30 passed**.
- Two additional datum/interpolation tests: **2 passed**.
- Scoped `ruff check --no-cache`: passed.
- `bash -n` on every OCTOPUS shell script: passed.
- Additional checks used small arrays, mocked notebook execution and in-memory file operations. Production generation, production QA and FVCOM integration were **not run**.

The owner-approved timestep warnings are not counted as defects.

## Previous findings

“RESOLVED” refers to the reported defect. Separate regressions are listed below.

### Round 1

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked NetCDF elevations become NaN and permit fallback. |
| F2 | RESOLVED | Failed required QA gates cause a nonzero build exit. |
| F3 | RESOLVED | Limiting enforces bounds, rejects nonconvergence and checks final rounded depths. |
| F4 | RESOLVED | Final band-floor counts and timestep allowances are reported under the owner’s warning policy. |
| F5 | RESOLVED | Resampling checks chords and both endpoints; publication permits only numerical tolerance. |
| F6 | RESOLVED | Missing sizing coverage is rejected at boundary, lattice and generated-node samples. |
| F7 | RESOLVED | Continuous departure segments and the full boundary receive land-crossing checks. |
| F8 | RESOLVED | Rejected designs preserve the existing CSV. Publication failure remains separately problematic: finding 6. |
| F9 | RESOLVED | Native values round-trip, and exported frozen-base checks run. |
| F10 | RESOLVED | Opposite-side seam incidence and outer-versus-base overlap are checked. |
| F11 | RESOLVED | Single-edge land segments are retained. |
| F12 | RESOLVED | Generation explicitly rejects constrained lines shorter than six nodes. |
| F13 | RESOLVED | CAO cache keys include resolved data root and dimensions. |
| F14 | RESOLVED | Invalid fine-grid results retain valid coarse coverage. |
| F15 | RESOLVED | Insufficient distinct points and collinear windows return uncovered samples. |
| F16 | RESOLVED | CAO source inventories include depth files. |
| F17 | RESOLVED | Dirty relevant sources are hashed, including renamed paths. |
| F18 | RESOLVED | Extension numeric controls, seeds, bounds and geographic limits are validated. |
| F19 | RESOLVED | Nonfinite spacing and nonadvancing resampling steps are rejected. |
| F20 | RESOLVED | Extension output directories receive exclusive reservations. |
| F21 | RESOLVED | Re-depth rejects source/output overlap and reserves its destination. |
| F22 | RESOLVED | Smoke staging reserves a fresh root, excluding stale history. |
| F23 | RESOLVED | Smoke depth control precedes timestep selection. |
| F24 | RESOLVED | Smoke case and root paths are resolved. |
| F25 | RESOLVED | The smoke job respects `FMESH_FVCOM`; relative overrides have a separate defect: finding 12. |
| F26 | RESOLVED | The 383 case names match preparation, and directory changes are guarded. |
| F27 | NOT RESOLVED | Direct package imports remain; owner decision pending. Finding 1, counted once. |
| F28 | RESOLVED | Extension and re-depth case names are restricted to filename components. |
| F29 | RESOLVED | Native reading validates unique, consecutive boundary traversal. |

### Round 2

| ID | Status | Reason |
|---|---|---|
| R2-1 | RESOLVED | Re-depth and smoke use exclusive directory reservations. |
| R2-2 | RESOLVED | Unique temporaries and a publication lock prevent competing-writer payload swaps. |
| R2-3 | RESOLVED | Re-depth repeats frozen-base, overlap and QA checks; sensitivity overrides are labelled. |
| R2-4 | RESOLVED | Actual chord lengths and endpoint floors are checked, including publication precision. |
| R2-5 | RESOLVED | Continuous departure testing catches thin land strips. |
| R2-6 | RESOLVED | Land-crossing validation tolerates numerical endpoint residue. |
| R2-7 | RESOLVED | Validation uses reconstructed published coordinates. |
| R2-8 | RESOLVED | Native coordinates and depths round-trip. |
| R2-9 | RESOLVED | Positive-area overlap with the base is checked beyond the seam. |
| R2-10 | RESOLVED | Porcelain parsing includes both rename paths. |
| R2-11 | RESOLVED | Geographic bounds are finite, ordered and range-checked. |
| R2-12 | RESOLVED | Resampling explicitly requires finite forward progress. |
| R2-13 | RESOLVED | Boundary direction follows the supplied second node. |
| R2-14 | RESOLVED | Convergence on the final permitted iteration is accepted. |
| R2-15 | RESOLVED | Nonfinite limited depths and invalid limiter controls are rejected. |
| R2-16 | PARTIAL | Both sides are now probed, but actual guides can still leave the domain: finding 5. |
| R2-17 | RESOLVED | Bands are materialized and checked against the accepted 5% deviation policy. |
| R2-18 | PARTIAL | Stage failures generally receive provenance; early failures and changing inputs remain: findings 7–8. |
| R2-19 | RESOLVED | Uncovered generated nodes are rejected. |
| R2-20 | RESOLVED | Invalid radii and reversals are rejected. |
| R2-21 | RESOLVED | Rejecting smoke-root reuse preserves the previous success marker. |
| R2-22 | NOT RESOLVED | Same outstanding F27; finding 1. |

### Round 3

| ID | Status | Reason |
|---|---|---|
| R3-1 | RESOLVED | Band-floor violations are counted on the final composed field; warnings follow owner policy. |
| R3-2 | RESOLVED | Publication uses a numerical tolerance instead of the previous 0.5% allowance. |
| R3-3 | RESOLVED | Materializing bands prevents iterator exhaustion. |
| R3-4 | PARTIAL | Both-side probing fixes the original selector error, but guide containment remains unchecked: finding 5. |
| R3-5 | PARTIAL | Hash verification prevents accepting mismatched products, but publication still destroys the previous coherent set: finding 6. |
| R3-6 | PARTIAL | Missing oceanmesh no longer passes `None` to `code_state`; early failures and input identity remain incomplete: findings 7–8. |
| R3-7 | RESOLVED | Valid NetCDF vertex/edge samples survive zero-weight masked neighbours. Separate shape regression: finding 3. |
| R3-8 | RESOLVED | fort.14 coordinates round-trip; 447 verifies the supplementary export. |
| R3-9 | PARTIAL | The original fan is detected, but the global tolerance can hide substantial local overlap: finding 4. |
| R3-10 | RESOLVED | Notebook 414 selects its timestep after depth control. |
| R3-11 | RESOLVED | Boundary spacing controls must be finite and positive. |
| R3-12 | RESOLVED | Short designs receive an explicit rejection before interior statistics. |
| R3-13 | RESOLVED | Filleting rejects reversals separately from straight continuation. |
| R3-14 | NOT RESOLVED | Same outstanding F27; finding 1. |

## Findings

1. **Major — Previously reported F27 remains open: package code directly imports oceanmesh.**  
   [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

   **Evidence:** Direct imports also remain in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672` and `autofinish/directives.py:20`. `THIRD_PARTY_NOTICES.md:37` endorses the arrangement contrary to the explicit repository policy.

   **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary, and reconcile the notices. **Previously reported; owner decision pending; counted once.**

2. **Major — M2 convergence can be falsely accepted because half-window amplitudes are interpolated incorrectly.**  
   [384_m2_analysis.py:264](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/384_m2_analysis.py:264).

   **Reproduction:** Fit unit-amplitude M2 signals at three nodes, with first-half phases `[-10,0,10]°` and second-half phases `[-20,0,20]°`. At equal station weights, the current calculation reports amplitude change approximately **0 m**, satisfying the **2 mm** convergence criterion. Interpolating the fitted coefficients gives station amplitudes **0.98987184 → 0.95979508 m**, a **30.08 mm** change. Arithmetic phase interpolation also turns `[359,1]° → [1,3]°` into a reported **−178°** change instead of **2°**.

   **Fix:** Retain each half’s fitted coefficients and interpolate them before deriving station amplitude and phase, as the full-window calculation already does. Pre-existing.

3. **Minor — The new bilinear interpolator loses the query shape.**  
   [dem/sources.py:64](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:64).

   **Reproduction:** Direct `Grid.depth()` sampling of a constant 2×2 grid with `(2,2)` longitude/latitude arrays returns shape **`(4,)`**. Previously it returned `(2,2)`. `shape_out` is captured after flattening. The priority-stack wrapper masks this regression by supplying one-dimensional queries.

   **Fix:** Capture and validate the original query shape before flattening, then restore it. Introduced by the R3-7 fix.

4. **Minor — The new overlap gate dilutes local defects against the entire mesh area.**  
   [qa.py:670](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:670).

   **Reproduction:** Combine the seven-element, twice-winding fan at radius 1 m with a triangle of area **5×10⁹ m²**. The actual `no_element_overlap` gate reports **PASS**, despite **1.31045 m²** of overlap in the fan—approximately **38%** of its summed area. This establishes a gate defect; the synthetic mesh does not pass every other gate.

   **Fix:** Check candidate element pairs with a tolerance tied to local element geometry and numerical precision. Keep the global area comparison as an additional diagnostic. Residual R3-9.

5. **Minor — Correctly oriented ladders can still place fixed constraints outside the extension domain.**  
   [446_extend_generate.py:153](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:153).

   **Reproduction:** Execute the actual selector and ladder on a six-node line with 1,000 m spacing and a 100 m-wide domain immediately to its left. The selector returns **True**. The ladder places both inner points **1,250 m** left of the line, where the synthetic signed distance is **+1,150 m**, and retains both because the land geometry is empty. Short side probes therefore do not validate the actual constraints.

   **Fix:** Check every guide point and guide segment against the extension domain; adapt the offset or reject unsupported geometry explicitly. Residual R2-16/R3-4.

6. **Minor — Failed publication still replaces the previous report and figure, leaving the old boundary unusable.**  
   [444_design_obc.py:232](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:232).

   **Reproduction:** In-memory execution of the actual publication block, failing the CSV `os.replace`, leaves **old CSV + new JSON + new PNG**, then releases the lock. The loader correctly rejects the hash mismatch, but the previous coherent report and figure have already been overwritten. A failed redesign consequently disables the previously usable boundary.

   **Fix:** Publish a complete versioned artifact set through one atomic pointer switch, retaining the previous version on failure. Residual R3-5.

7. **Minor — Provenance failures occur before the failure-report handler is registered.**  
   [445_extend_mesh.py:71](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:71), [453_redepth_extended.py:72](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:72).

   **Reproduction:** Injecting an exception into `collect()` during the actual 445 initialization produces **zero report writes and zero registered handlers**, after output reservation. Notebook 453 has the same ordering. In 445, unset `DATA_DIR` also exits after reservation but before handler registration.

   **Fix:** Initialize and persist the run record immediately after reservation, then guard preflight and provenance collection with structured failure reporting. Represent incomplete provenance explicitly. Residual R3-6.

8. **Minor — Recorded provenance is not bound to the inputs actually consumed by the stages.**  
   [445_extend_mesh.py:93](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:93).

   **Evidence:** Both stages reload the live recipe; the parent records its earlier settings and hashes without checking for subsequent changes. A mocked execution with stage settings changed from `fin_seed: 42` to `99` still produces **`status: "ok"`**, parent settings **42**, stage settings **99**, and the original recipe hash. Base files are omitted from initial provenance and hashed only after processing at line 109. Re-depth likewise omits the base files from its provenance collection.

   **Fix:** Run from immutable recipe/input snapshots where practical; otherwise verify identities before consumption and at completion. Capture base identity initially and reject or clearly label changed inputs. Residual R3-6’s input-identity requirement.

9. **Minor — CAO interpolation still discards valid grid-centre samples.**  
   [dem/sources.py:229](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:229).

   **Reproduction:** An in-memory CAO grid `[[10,20],[30,NaN]]`, queried exactly at its valid upper-left centre, returns **NaN** because the zero-weight NaN term propagates. Separately, an entirely valid constant grid queried at its easternmost centre returns **NaN** because `fi < nx-1` excludes that centre. These cases unnecessarily select a fallback or lose coverage.

   **Fix:** Apply the same zero-weight-aware interpolation used for NetCDF grids, with inclusive centre bounds and clipped cell indices. Pre-existing; the R3-7 remedy covers only `Grid`.

10. **Minor — Boundary design can publish a self-crossing boundary as valid.**  
    [444_design_obc.py:163](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:163).

    **Reproduction:** With a straight northern coast, use metric vertices  
    `[(0,0),(0,-10000),(100000,-40000),(0,-40000),(100000,-10000),(100000,0)]`.  
    The actual acceptance block reports **90° at both ends**, **zero land crossing**, sufficient spacing and **no problems**, although `LineString.is_simple` is false. Filleting with 5,000 m radii and resampling at 3,000 m also succeeds: **73 nodes**, minimum edge approximately **3,000 m**, still self-crossing. Generation later rejects the domain.

    **Fix:** Validate simplicity, repeated vertices and the intended closed domain before publishing any products. Pre-existing.

11. **Minor — Notebook 414 accepts a negative manual timestep.**  
    [414_refine_m2_prep.py:131](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/414_refine_m2_prep.py:131).

    **Reproduction:** Execute the actual override guards with `dte=-1`, allowance `100`, interval `1800` and `ISPLIT=10`. They accept and print **`EXTSTEP_SECONDS = -1.0`**, because the negative internal step still divides the interval and does not exceed the allowance. The later mesh-margin comparison also does not reject it.

    **Fix:** Require a supplied timestep to be finite and strictly positive before divisibility or allowance checks. Pre-existing.

12. **Minor — Relative FVCOM executable overrides change meaning after the job enters the run directory.**  
    [448_extend_smoke.sh:45](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/448_extend_smoke.sh:45).

    **Evidence:** `FMESH_FVCOM=../FVCOM/src/fvcom` identifies the sibling executable from the submission repository and is hashed there at line 26. The MPI invocation follows `cd "$CASE_DIR"`, so the same string now identifies a different path beneath the run root. Other portable FVCOM run jobs use the same pattern.

    **Fix:** Resolve and validate the chosen executable to an absolute path before changing directories. This affects the executable overrides added by the portability work.

13. **Minor — Smoke staging accepts paths that FVCOM truncates.**  
    [448_extend_smoke.py:98](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/448_extend_smoke.py:98), [383_m2_case_prep.py:155](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:155).

    **Evidence:** The inspected FVCOM source declares `INPUT_DIR` and `OUTPUT_DIR` as `CHARACTER(LEN=80)` in `mod_main.F:183`. From this repository, `--root scratch/smoke_round4` produces an input-directory string **86 characters** long. Staging writes it and emits `STAGED`; FVCOM subsequently reads a truncated path. The refinement workflow already recognizes this limit, but extension smoke staging does not.

    **Fix:** Validate the final namelist directory strings against the solver’s capacity before reservation or export, and require a shorter root. Pre-existing, exposed by portable/custom paths.

14. **Minor — The smoke checker accepts history fields with no spatial data.**  
    [cli/check_run.py:180](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:180).

    **Reproduction:** In-memory history with two valid half-hour timestamps and `zeta`, `ua`, `va` each shaped **`(2,0)`**, together with matching dates and a `TADA` log, returns **`ok: True`** and no reasons. Empty arrays satisfy the finiteness check. Notebook 448’s job can therefore accept structurally invalid history.

    **Fix:** Validate variable dimensions and nonzero spatial sizes, and compare node/element counts with the staged mesh. Pre-existing robustness defect under injected corrupt output.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 12 minor, 0 nit)

### Prompt

```markdown
# Review request, round 4: extending a base mesh outward (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (small synthetic inputs only; do not read the large data
under $DATA_DIR beyond listing it, and do not submit batch jobs). Answer in
English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
A tool that keeps a finished FVCOM base mesh exactly as it is and adds the
sea out to a new, designed open boundary (USER_GUIDE section 13). Read:

- `git show b0584f9 8e2739b 450ad44 d6d2a72 7044b6b 0f52d5b 69b50a4 d9e92fd b6d2ed8 765423c`
  (the extension tool and its documentation), and the current files:
  - `src/fvcom_mesh_tools/extend.py`, `extend_recipe.py`, `obc_design.py`,
    `dem/sources.py` (named bathymetry sources, priority stack, and the new
    `DATUM` registry / `non_tp_count` warning);
  - `notebooks/444_design_obc.py`, `445_extend_mesh.py`, `446_extend_generate.py`,
    `447_extend_merge.py`, `448_extend_smoke.py`, `453_redepth_extended.py`;
  - `recipes/extend/tokyo_bay_enshu.yaml`, `tokyo_bay_enshu_obc_design.yaml`;
  - `jobs/octopus/444_design_obc.sh`, `445_extend_mesh.sh`, `448_extend_smoke.sh`,
    `453_redepth_extended.sh`, `common.sh`;
  - tests: `tests/test_extend*.py`, `tests/test_obc_design*.py`,
    `tests/test_dem_sources.py` (whatever exists).
- Also in scope, just committed: the portability change --
  every job script and `common.sh` now take paths only from `$DATA_DIR` and
  `$WORK_DIR` (login profile), stop when they are unset, and derive the
  OCTOPUS FVCOM library directory as `FVCOM_LIBS` in `common.sh`; notebooks
  383/384/414 and `cli/refine_run.py` no longer fall back to `/octfs/...`.
  See commits 6d8b9a7 and 6c068d2 (`git log -5`).

Design intent:
- the base mesh's nodes, elements and depths are carried bit for bit
  (`verify_frozen_base`);
- the new part is generated with oceanmesh (subprocess; GPL code must never
  be imported by `fvcom_mesh_tools`), with fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-3 and their triage are in docs/extend-tools-review-20261001.md.
Round 3 (13 findings + F27) was fixed in 7007e7e; read `git show 7007e7e`.
Owner decision (2026-10-01): meshes are made from the real depths and the
time step is settled in the depth stage, so the band-floor check (446) and
the new-element time-step gate (447, 453) REPORT warnings and values and do
not fail the build. Please do not report that as a defect.
Integration run on real data after 7007e7e: the full 445 build completed,
QA 23/23 (with the new no_element_overlap gate), NP 14,740 / NE 27,135 (the
same counts as the existing mesh), warnings as designed.
F27 (oceanmesh imported in package code): owner decision pending; report it
as still open, not as new.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the scope and everything it touches.

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
| F1 (F27) | major | as round 1 | yes, pre-existing | owner decision pending |
| F2 | major | numeric check of the reviewer's case (0.9899 -> 0.9598 m; 2 deg not -178) | yes | fixed, b3457b4 |
| F3 | minor | code read (shape taken after ravel) | yes | fixed, b3457b4; test |
| F4 | minor | test (1 m fan beside a 5e9 m2 element) | yes | fixed (pairwise gate), b3457b4; test |
| F5 | minor | code read (guides tested against land only) | yes | fixed, b3457b4; real build |
| F6 | minor | fault injection (CSV replace fails) | yes | fixed (stage, keep .prev, roll back), b3457b4 |
| F7 | minor | run without DATA_DIR (no report before the fix) | yes | fixed, b3457b4 |
| F8 | minor | code read (live recipe re-read by stages) | yes | fixed (hash again at the end), b3457b4; test |
| F9 | minor | test (valid centre beside NaN; last column) | yes | fixed, b3457b4; test. Also found: projection round-off (3e-11 cells) gave a NaN corner weight 1e-21; indices within 1e-9 of a grid line are snapped |
| F10 | minor | shapely is_simple on the reviewer's vertices | yes | fixed, b3457b4 |
| F11 | minor | call with -1, NaN, 0 | yes | fixed, b3457b4 |
| F12 | minor | code read (cd before mpiexec) | yes | fixed (fmesh_fvcom in common.sh, all 9 jobs), b3457b4 |
| F13 | minor | path length of the reviewer's root (86) | yes | fixed (383/414/448, before writing), b3457b4 |
| F14 | minor | code read (empty arrays are finite) | yes | fixed, b3457b4; tests |

## Round 5

### Answer (verbatim)

Read-only review; no files changed. The pre-existing untracked PNG remains untouched. No production mesh, bathymetry data or batch jobs were used.

Verification: **30 pure extension/OBC tests passed**, four in-memory DEM regression checks passed, mocked recipe checks passed, scoped Ruff passed, and every OCTOPUS script passed `bash -n`. Initial pytest attempts failed because capture and the repository’s conftest required writable temporary/cache directories; the successful run used `--noconftest -s -p no:cacheprovider` with bytecode and plugin autoload disabled.

**Previous findings — round 1**

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked grid cells become NaN, permitting fallback. |
| F2 | RESOLVED | Required QA failures produce a nonzero build exit. |
| F3 | RESOLVED | Bounds are enforced during limiting; convergence and final rounded r-factor are checked. |
| F4 | RESOLVED | Floor violations and final timestep allowances are reported under the owner’s warning policy. |
| F5 | RESOLVED | Resampling checks actual chords and both endpoint floors; publication checks serialized geometry. |
| F6 | RESOLVED | Missing sizing coverage is rejected at boundary, lattice and generated-node samples. |
| F7 | RESOLVED | Departure segments and the full boundary receive land-crossing checks. |
| F8 | RESOLVED | Rejected designs preserve the published CSV. Publication recovery has a separate residual below. |
| F9 | RESOLVED | Native serialization round-trips and exported frozen-base verification runs. |
| F10 | RESOLVED | Opposite-side seam incidence and outer-versus-base overlap are checked. |
| F11 | RESOLVED | Single-edge land runs are retained. |
| F12 | RESOLVED | Generation explicitly rejects constrained lines shorter than six nodes. |
| F13 | RESOLVED | CAO cache keys include resolved root and dimensions. |
| F14 | RESOLVED | Invalid fine-grid samples retain coarse coverage. |
| F15 | RESOLVED | Insufficient distinct points and collinear M7001 windows return uncovered samples. |
| F16 | RESOLVED | CAO inventories include depth files. |
| F17 | RESOLVED | Relevant dirty sources are hashed, including rename paths. |
| F18 | RESOLVED | Recipe numbers, integer controls, seeds and bounds are validated. |
| F19 | RESOLVED | Resampling rejects nonfinite spacing and lack of forward progress. |
| F20 | RESOLVED | Extension outputs receive exclusive reservations. |
| F21 | RESOLVED | Re-depth rejects overlapping destinations and reserves its output. |
| F22 | RESOLVED | Smoke staging reserves a fresh root. |
| F23 | RESOLVED | Depth control precedes timestep selection. |
| F24 | RESOLVED | Smoke paths are resolved before use. |
| F25 | RESOLVED | The executable override is respected and resolved before changing directories. |
| F26 | RESOLVED | Case names match preparation; directory changes are guarded. |
| F27 | NOT RESOLVED | Direct GPL imports remain; owner decision pending. Counted once below. |
| F28 | RESOLVED | Case prefixes are restricted to filename components. |
| F29 | RESOLVED | Native reading validates unique, consecutive boundary traversal. |

**Previous findings — round 2**

| ID | Status | Reason |
|---|---|---|
| R2-1 | RESOLVED | Re-depth and smoke use exclusive output reservations. |
| R2-2 | RESOLVED | Unique temporaries and a publication lock prevent competing-writer swaps. |
| R2-3 | RESOLVED | Re-depth repeats frozen-base, overlap and QA checks; sensitivity overrides are labelled. |
| R2-4 | RESOLVED | Chord lengths, endpoint floors and publication precision are checked. |
| R2-5 | RESOLVED | Continuous departure testing catches thin land strips. |
| R2-6 | RESOLVED | Land-crossing checks tolerate only numerical endpoint residue. |
| R2-7 | RESOLVED | Validation uses reconstructed published coordinates. |
| R2-8 | RESOLVED | Native coordinates and depths round-trip. |
| R2-9 | RESOLVED | Overlap with the base is checked beyond the seam. |
| R2-10 | RESOLVED | Porcelain parsing includes both rename paths. |
| R2-11 | RESOLVED | Geographic bounds receive finite/range/order validation. |
| R2-12 | RESOLVED | Resampling requires finite forward progress. |
| R2-13 | RESOLVED | Traversal direction follows the supplied second node. |
| R2-14 | RESOLVED | Convergence on the last allowed iteration is accepted. |
| R2-15 | RESOLVED | Limited-edge depths and limiter controls are validated. |
| R2-16 | PARTIAL | Orientation and guide-point containment are checked; guide segments can still cross land between samples: finding 4. |
| R2-17 | RESOLVED | Bands are materialized and checked against the accepted deviation policy. |
| R2-18 | PARTIAL | Reporting and input hashing improved; findings 3 and 8 remain. |
| R2-19 | RESOLVED | Uncovered generated nodes are rejected. |
| R2-20 | RESOLVED | Invalid radii, zero sides and reversals are rejected. |
| R2-21 | RESOLVED | Rejecting root reuse preserves the previous success marker. |
| R2-22 | NOT RESOLVED | Same outstanding F27. |

**Previous findings — round 3**

| ID | Status | Reason |
|---|---|---|
| R3-1 | RESOLVED | Violations are counted on the composed field and reported under owner policy. |
| R3-2 | RESOLVED | Publication allows only numerical spacing tolerance. |
| R3-3 | RESOLVED | Materializing bands prevents iterator exhaustion. |
| R3-4 | PARTIAL | Both sides are probed; segment containment remains incomplete: finding 4. |
| R3-5 | PARTIAL | Ordinary replacement failures roll back; rollback failures destroy unrecovered backups: finding 5. |
| R3-6 | PARTIAL | Missing dependencies and early provenance failures are handled; findings 3 and 8 remain. |
| R3-7 | RESOLVED | Zero-weight masked neighbours no longer invalidate Grid samples. |
| R3-8 | RESOLVED | fort.14 values round-trip and the supplementary export is verified. |
| R3-9 | RESOLVED | Pairwise overlap checking removes global-area dilution. |
| R3-10 | RESOLVED | Notebook 414 selects its timestep after depth control. |
| R3-11 | RESOLVED | Boundary spacing controls must be finite and positive. |
| R3-12 | RESOLVED | Short designs receive explicit rejection. |
| R3-13 | RESOLVED | Filleting rejects reversals. |
| R3-14 | NOT RESOLVED | Same outstanding F27. |

**Previous findings — round 4**

| ID | Status | Reason |
|---|---|---|
| R4-F1 | NOT RESOLVED | Same outstanding F27. |
| R4-F2 | RESOLVED | Half-window amplitudes/phases now derive from interpolated harmonic coefficients. |
| R4-F3 | RESOLVED | `_bilinear` preserves query shape. The separate CAO implementation still flattens queries: finding 9. |
| R4-F4 | RESOLVED | Overlap tolerance is tied to the smaller meeting element. |
| R4-F5 | PARTIAL | Guide points are checked, but nine segment samples do not establish containment: finding 4. |
| R4-F6 | PARTIAL | A single replacement failure restores the old set; a rollback failure remains destructive: finding 5. |
| R4-F7 | PARTIAL | Handler registration precedes provenance, but follows a fallible print: finding 8. |
| R4-F8 | PARTIAL | Persistent changes after initial hashing are detected; changes between recipe loading and initial hashing escape: finding 3. |
| R4-F9 | RESOLVED | CAO preserves valid centres beside NaNs and includes the last row/column. |
| R4-F10 | RESOLVED | Boundary simplicity and repeated adjacent nodes are checked. |
| R4-F11 | RESOLVED | Manual timestep overrides must be finite and positive. |
| R4-F12 | RESOLVED | All nine run jobs resolve/check the executable before changing directories. |
| R4-F13 | PARTIAL | ASCII path lengths are checked; multibyte paths bypass the solver’s byte capacity: finding 7. |
| R4-F14 | PARTIAL | Empty/inconsistent spatial arrays are rejected; an unreadable named grid silently disables mesh-size checking: finding 6. |

**Findings**

1. **Major — Previously reported F27 remains open: package code directly imports oceanmesh.**  
   [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

   **Evidence:** Imports remain there and in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672` and `autofinish/directives.py:20`. These conflict with the explicit repository policy; `THIRD_PARTY_NOTICES.md` still endorses the conflicting arrangement.

   **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary and reconcile the notices. **Previously reported, owner decision pending; counted once.**

2. **Major — A finishing flip can corrupt a valid triangulation while reporting successful repair.**  
   [algorithms/obc_finish.py:265](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/algorithms/obc_finish.py:265).

   **Reproduction:** Call `flip_c4_edges` with default controls on:

   ```python
   nodes = np.array([[3951, 982], [2396, 1743],
                     [926, 262], [2471, 635]], float)
   elements = np.array([[3, 2, 0], [1, 3, 0]])
   ```

   The original triangles do not overlap. The function reports `fixed: [(0,3)]` and produces `[[1,2,3],[1,2,0]]`, introducing **869,917.5 m² of overlap** and the same footprint change. Independently reversing each proposed triangle to make its area positive hides an illegal flip across a concave quadrilateral.

   **Fix:** Require a strictly convex quadrilateral and a valid crossing of its diagonals before flipping. Preserve the original patch footprint and refresh affected edge ownership after accepted flips. Pre-existing; reached by 447’s finishing chain. Final overlap QA protects the build from accepting this particular corruption.

3. **Minor — Recipe loading precedes its initial hash, leaving a provenance race.**  
   [445_extend_mesh.py:45](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:45), [453_redepth_extended.py:92](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:92).

   **Reproduction:** Execute the actual 445 driver with virtual file I/O. Load `fin_seed: 42`, then change the live recipe to `99` during source preflight, before `collect`. Both stages consume `99`; both hash comparisons see `99`. The driver exits successfully with **`status: "ok"`**, parent settings **42**, stage settings **99**, and **`inputs_changed_during_build: []`**. Re-depth likewise loads its controls before hashing the recipe.

   **Fix:** Parse and hash the same recipe bytes, retaining the recipe’s original location for relative-path resolution. Compare that digest with the live file before launching work and at completion. Residual R4-F8.

4. **Minor — Ladder segments can cross land between the nine containment samples.**  
   [446_extend_generate.py:165](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:165).

   **Reproduction:** Execute the actual `ladder` function on a six-node straight boundary with 1,000 m spacing. Place a 10 m-wide island at **35%** along its inner guide segment. Both guide-point clearance checks pass; samples at 10%, 20%, …, 90% miss the island. The result is **`keep=[True,True]`, `join=[True]`**, although the joined `LineString` intersects land.

   **Fix:** Check the whole segment against land geometry. Establish containment in the generation domain geometrically, or use an adaptive check with a defensible distance bound. Residual R4-F5; the new segment check is insufficient.

5. **Minor — A rollback failure deletes backups that have not been restored.**  
   [444_design_obc.py:254](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:254).

   **Reproduction:** Execute the actual publication block with virtual files. Fail CSV replacement, then fail restoration of the previous JSON. The old CSV remains, but the new JSON/PNG remain published; cleanup **deletes the old PNG backup** and releases the lock. The loader rejects the mixed CSV/report set. A single CSV-replacement failure correctly restores all three files.

   **Fix:** Attempt every restoration independently, remove a backup from recovery tracking only after restoration succeeds, and retain all unrecovered backups plus a recovery marker when rollback fails. Introduced by the R4-F6 rollback implementation.

6. **Minor — An unreadable named staged grid silently disables history mesh-size validation.**  
   [cli/check_run.py:169](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:169).

   **Reproduction:** An in-memory history stack with valid timestamps, finite fields and a `TADA` log returns **`ok: True`** when the namelist names `INPUT_DIR` and `GRID_FILE`, but opening that grid raises `FileNotFoundError`. `_grid_counts` returns `None`, so the new comparison is skipped.

   **Fix:** When a grid is named, inability to read and validate its counts must add a failure reason. Preserve optional checking only for namelists that genuinely omit the grid specification. Gap in the R4-F14 fix.

7. **Minor — The FVCOM path guard counts characters instead of encoded bytes.**  
   [383_m2_case_prep.py:157](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:157).

   **Reproduction:** The actual guard accepts:

   ```python
   Path("/tmp") / ("漢" * 26) / "extended/input"
   ```

   Including the trailing slash, the namelist directory contains **47 Unicode characters but 99 UTF-8 bytes**. FVCOM declares `INPUT_DIR` and `OUTPUT_DIR` as default `CHARACTER(LEN=80)` in `mod_main.F:183–184`, so this path exceeds their capacity.

   **Fix:** Check the byte length using the encoding actually used to write the namelist, or explicitly restrict these directory strings to ASCII. Applies to 383, 414 and 448. Introduced by the R4-F13 guard.

8. **Minor — The first post-reservation print still precedes failure-handler registration.**  
   [445_extend_mesh.py:51](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:51).

   **Reproduction:** Execute the actual initialization block and inject `BrokenPipeError` into the first `say` call. The output has been reserved, but **zero exit handlers are registered** and `STATE` has not been initialized. Consequently no failure report is written.

   **Fix:** Initialize state and register the handler immediately after `reserve`, before logging or any other fallible operation. Residual R4-F7.

9. **Minor — Direct CAO sampling still loses multidimensional query shape.**  
   [dem/sources.py:244](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:244).

   **Reproduction:** Direct `CaoNested.depth` sampling with `(2,2)` longitude/latitude arrays returns shape **`(4,)`**. An empty mocked area inventory reproduces this without accessing data; populated inventories follow the same flattened return path. Grid and M7001 sampling preserve the query shape.

   **Fix:** Capture and validate query shapes before flattening, then reshape the result. Pre-existing and separate from the resolved `_bilinear` regression.

10. **Minor — Mesh QA passes NaN and infinite depths.**  
    [qa.py:1000](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:1000).

    **Reproduction:** On the pristine 16-node/18-element test grid, replace one interior depth with NaN or positive infinity. `run_qa(..., channel_check=False)` reports **zero failed gates** in both cases. With NaN, the depth gate even reports **`min = nan m, violations = 0`**. Only comparing depths with the minimum does not establish finiteness.

    **Fix:** Add a required finite-depth check before depth/CFL calculations. Keep nonfinite depths explicitly failing even when informational timestep checks are disabled. Pre-existing; native extension input reading guards this case, but general mesh QA does not.

**Numbering discrepancy**

A repository-local **amplifier** is confirmed, though I did not establish the source of the initial 446 numbering variation. [mesh_clean_phase_h.py:3870](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_clean_phase_h.py:3870) processes failing elements in ID order and distributes seeded random trials in that order. On an identical small mesh with only element rows reversed, seed 42 and a one-pass, 20-trial budget produced **200.927 m maximum geometric difference**, with **3 versus 5** reported fixes. The subsequent vertex-removal chain also uses numbering-dependent candidate ordering.

This explains how numbering changes can affect finishing, but does **not** prove which operation removed the two production nodes. To isolate that, canonicalize outer node/element ordering before finishing and record mesh fingerprints/counts after each finishing stage; retain the frozen base’s ordering during merge.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 8 minor, 0 nit)

### Prompt

```markdown
# Review request, round 5: extending a base mesh outward (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (small synthetic inputs only; do not read the large data
under $DATA_DIR beyond listing it, and do not submit batch jobs). Answer in
English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
A tool that keeps a finished FVCOM base mesh exactly as it is and adds the
sea out to a new, designed open boundary (USER_GUIDE section 13). Read:

- `git show b0584f9 8e2739b 450ad44 d6d2a72 7044b6b 0f52d5b 69b50a4 d9e92fd b6d2ed8 765423c`
  (the extension tool and its documentation), and the current files:
  - `src/fvcom_mesh_tools/extend.py`, `extend_recipe.py`, `obc_design.py`,
    `dem/sources.py` (named bathymetry sources, priority stack, and the new
    `DATUM` registry / `non_tp_count` warning);
  - `notebooks/444_design_obc.py`, `445_extend_mesh.py`, `446_extend_generate.py`,
    `447_extend_merge.py`, `448_extend_smoke.py`, `453_redepth_extended.py`;
  - `recipes/extend/tokyo_bay_enshu.yaml`, `tokyo_bay_enshu_obc_design.yaml`;
  - `jobs/octopus/444_design_obc.sh`, `445_extend_mesh.sh`, `448_extend_smoke.sh`,
    `453_redepth_extended.sh`, `common.sh`;
  - tests: `tests/test_extend*.py`, `tests/test_obc_design*.py`,
    `tests/test_dem_sources.py` (whatever exists).
- Also in scope, just committed: the portability change --
  every job script and `common.sh` now take paths only from `$DATA_DIR` and
  `$WORK_DIR` (login profile), stop when they are unset, and derive the
  OCTOPUS FVCOM library directory as `FVCOM_LIBS` in `common.sh`; notebooks
  383/384/414 and `cli/refine_run.py` no longer fall back to `/octfs/...`.
  See commits 6d8b9a7 and 6c068d2 (`git log -5`).

Design intent:
- the base mesh's nodes, elements and depths are carried bit for bit
  (`verify_frozen_base`);
- the new part is generated with oceanmesh (subprocess; GPL code must never
  be imported by `fvcom_mesh_tools`), with fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-4 and their triage are in docs/extend-tools-review-20261001.md.
Round 4 (your previous answer; 14 findings) was fixed in b3457b4; read
`git show b3457b4`. Summary of what was done per finding:
- F1 (= F27, oceanmesh imported in package code): owner decision pending;
  report it as still open, not as new.
- F2 384: half-window amplitude/phase now derived from interpolated
  harmonic coefficients, like the full window.
- F3 `_bilinear` keeps the query shape (new `_bilinear_index` core).
- F4 `no_element_overlap`: pairwise check of meeting elements (STRtree,
  `touches` excluded) against the smaller element's area; the global
  summed-minus-union area is kept as a diagnostic.
- F5 446 `ladder`: guide points kept only with sdf < -h/4 (inside the
  domain) and off land; guide segments joined only when 9 interior samples
  are inside the domain.
- F6 444: all products staged to temporaries first; published files kept
  as `.prev`, replaced, and rolled back if a replacement fails. A process
  death half-way is left to the hash check of the loader (documented).
- F7 445/453: exit handler registered right after `reserve`; a failure
  before provenance writes `"provenance": null, "provenance_complete": false`.
- F8 445/453: recipe, open boundary (+ its report if present), base and
  (453) the built case are hashed at the start and again at the end
  (`provenance.changed_files`); a change fails the build, and in 453 no
  `--allow-failing-gates` accepts it. We did not copy inputs to snapshots
  because the recipe's relative paths and the loader's sidecar check bind
  to their locations; detection at the end is the chosen contract.
- F9 CAO uses the same zero-weight-aware interpolation with inclusive
  bounds. We also found that the projection round-off (3e-11 cells) put a
  1e-21 weight on a NaN corner, so indices within 1e-9 of a grid line are
  snapped (`_snap`), in both the NetCDF and CAO paths.
- F10 444: `LineString.is_simple` and repeated nodes are checked.
- F11 414: `--dte` must be finite and > 0 (0 no longer means "auto").
- F12 `common.sh` `fmesh_fvcom` resolves and checks the executable; all
  nine FVCOM job scripts use it.
- F13 `check_fvcom_dirs` (383) refuses run directories over 80 characters;
  called by `namelist()`, and before any write in 383, 414 and 448 (448:
  before `reserve`).
- F14 `check_run`: (time, space) with space > 0, the same in every stack,
  and equal to the staged grid's counts when the namelist names
  INPUT_DIR/GRID_FILE.
Verification after b3457b4: full test suite 1121 passed (batch
job 122652); on real data 444 published the Enshu boundary (191 nodes, no
problems), the full 445 build passed QA 23/23 (pairwise overlap gate
included) with ladders fully kept and joined (interface 9/9, 8 segments;
open boundary 189/189, 188 segments), and 453 re-depth completed.
Observation we are checking ourselves (job 122678): the 445 build gave
NP 14,738 / NE 27,133, the previous build (same sizing field, bit-equal to
1e-12) 14,740 / 27,135. The 446 DistMesh output had the same counts and the
same node set to 0.2 mm but a different node numbering, and 447's finishing
then diverged by two nodes. So generation looks run-to-run nondeterministic
at round-off level; if you can locate the cause in this repository's code
(not the oceanmesh fork), report it.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths and the time step is settled in the depth stage, so the band-floor
check (446) and the new-element time-step gate (447, 453) REPORT warnings
and do not fail the build. Please do not report that as a defect.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the scope and everything it touches.

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
| F1 (F27) | major | as round 1 | yes, pre-existing | owner decision pending |
| F2 | major | reviewer's 4-node case run: flip reported fixed, overlap | yes | fixed (`_flip_pair`, all three flip sites), e2918dc; tests |
| F3 | minor | code read (load before hash) | yes | fixed (digest of parsed bytes, passed to stages), e2918dc; test |
| F4 | minor | code read (samples only) | yes | fixed (whole-segment test), e2918dc |
| F5 | minor | fault injection (CSV replace, then JSON restore fail) | yes | fixed (.RECOVER marker), e2918dc; test |
| F6 | minor | code read (None skips the check) | yes | fixed, e2918dc; test |
| F7 | minor | guard called on a 26-kanji path | yes | fixed (ASCII and bytes), e2918dc |
| F8 | minor | code read (say before register) | yes | fixed, e2918dc |
| F9 | minor | code read (ravel, no reshape) | yes | fixed, e2918dc; test |
| F10 | minor | code read (NaN < min is False) | yes | fixed, e2918dc; tests |

Numbering discrepancy (raised by us in the prompt): resolved. Two 446 runs of
the same code in one job (122678) gave bit-identical meshes, identical also
to the 445 build of job 122652, so generation is deterministic. The change
from the round-3 build came from the round-4 interpolation rewrite (F9):
the sizing depths moved by up to 1.8e-12 m, enough to change DistMesh's
node numbering, which the order-dependent finishing then turned into two
nodes. Same inputs and code give the same mesh; a code change may move a
mesh at round-off level.

## Round 6

### Answer (verbatim)

Read-only review completed. **11 findings remain: 6 major and 5 minor.** Finding 1 is the previously reported, pending F27. Finding 7 is a defect in the new recovery cleanup.

No files were modified, shared datasets were read, or batch jobs were submitted. The existing untracked `recipes/extend/tokyo_bay_enshu_obc.png` was left untouched.

Verification in the `fvcom-mesh-tools` conda environment:

- `pytest -q -s --noconftest -p no:cacheprovider tests/test_extend.py tests/test_obc_design.py tests/test_obc_finish.py`, with bytecode and plugin autoload disabled: **32 passed**.
- Five selected in-memory DEM tests covering datum registration, masked neighbours, interpolation shape and Cao interpolation: **5 passed**.
- Scoped `ruff check … --no-cache`: **passed**.
- `bash -n` invoked independently for every OCTOPUS job script: **56 scripts, no failures**.
- Additional reproductions used synthetic geometry, diskless NetCDF, virtual file operations and shell mocks.
- Normal pytest startup failed because the sandbox has no writable temporary directory; the capture-disabled retry encountered Matplotlib’s writable-cache requirement through `conftest.py`. The full suite, production meshes and FVCOM integration were **not run**.

**Previous findings — round 1**

Statuses refer to the originally reported defect; separately discovered gaps are identified below.

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked grid cells become NaN and permit fallback. |
| F2 | RESOLVED | Required QA failures produce a nonzero build exit. |
| F3 | RESOLVED | Limiting enforces bounds, convergence and the final rounded r-factor. |
| F4 | RESOLVED | Floor violations and timestep allowances follow the owner’s explicit warning policy. |
| F5 | RESOLVED | Actual chords, endpoint floors and serialized geometry are checked. |
| F6 | RESOLVED | Missing boundary, wet-lattice and generated-node sizing coverage is rejected. |
| F7 | RESOLVED | Continuous departures and the entire boundary are checked against land. |
| F8 | RESOLVED | Validation precedes publication; rejected designs preserve the CSV. Recovery remains partial separately. |
| F9 | RESOLVED | Serialization round-trips and exported frozen-base verification runs. |
| F10 | RESOLVED | Opposite-side seam incidence and overlap beyond the seam are checked. |
| F11 | RESOLVED | Single-edge land runs are retained. |
| F12 | RESOLVED | Unsupported constrained lines shorter than six nodes are explicitly rejected. |
| F13 | RESOLVED | Cao cache keys include the resolved root and dimensions. |
| F14 | RESOLVED | Invalid fine samples retain valid coarse coverage. |
| F15 | RESOLVED | Insufficient distinct points and collinear M7001 windows return uncovered samples. |
| F16 | RESOLVED | Cao inventories include the depth files. |
| F17 | RESOLVED | Relevant dirty sources and rename paths are hashed. |
| F18 | RESOLVED | Numeric settings, integer controls, seeds and bounds are validated. |
| F19 | RESOLVED | Resampling rejects nonfinite spacing and failure to progress. |
| F20 | RESOLVED | Extension outputs receive exclusive reservations. |
| F21 | RESOLVED | Re-depth rejects overlapping destinations and reserves its output. |
| F22 | RESOLVED | Smoke staging reserves a fresh root. |
| F23 | RESOLVED | Depth control precedes timestep selection. |
| F24 | RESOLVED | Smoke paths are resolved before use. |
| F25 | RESOLVED | Executable overrides are respected and resolved before directory changes. |
| F26 | RESOLVED | Case names match preparation and directory changes are guarded. |
| F27 | NOT RESOLVED | Package imports still conflict with repository policy; owner decision pending. Finding 1. |
| F28 | RESOLVED | Case prefixes are restricted to filename components. |
| F29 | RESOLVED | Native reading validates unique, consecutive boundary traversal. |

**Previous findings — round 2**

| ID | Status | Reason |
|---|---|---|
| R2-1 | RESOLVED | Re-depth and smoke reserve their outputs exclusively. |
| R2-2 | RESOLVED | Unique temporaries and a publication lock prevent competing-writer swaps. |
| R2-3 | RESOLVED | Re-depth repeats frozen-base, overlap and QA checks; overrides are labelled. |
| R2-4 | RESOLVED | Chords, endpoint floors and publication precision are checked. |
| R2-5 | RESOLVED | Continuous departure checks catch thin land strips. |
| R2-6 | RESOLVED | Land-crossing tolerance is limited to numerical endpoint residue. |
| R2-7 | RESOLVED | Checks use reconstructed published coordinates. |
| R2-8 | RESOLVED | Native coordinates and depths round-trip. |
| R2-9 | RESOLVED | Outer-versus-base overlap is checked beyond the seam. |
| R2-10 | RESOLVED | Porcelain parsing includes both rename paths. |
| R2-11 | RESOLVED | Geographic bounds receive finite, range and order validation. |
| R2-12 | RESOLVED | Resampling requires finite forward progress. |
| R2-13 | RESOLVED | Traversal direction follows the supplied second node. |
| R2-14 | RESOLVED | Convergence on the final permitted iteration is accepted. |
| R2-15 | RESOLVED | Limited-edge depths and limiter controls are validated. |
| R2-16 | RESOLVED | Both-side selection, guide points and complete guide segments are checked. |
| R2-17 | RESOLVED | Bands are materialized and checked against the accepted deviation policy. |
| R2-18 | RESOLVED | Failure reporting starts before fallible logging; parsed recipe identities bind the stages. |
| R2-19 | RESOLVED | Uncovered generated nodes are rejected. |
| R2-20 | RESOLVED | Invalid radii, zero sides and reversals are rejected. |
| R2-21 | RESOLVED | Rejecting smoke-root reuse preserves previous success markers. |
| R2-22 | NOT RESOLVED | Same pending F27; counted once as finding 1. |

**Previous findings — round 3**

| ID | Status | Reason |
|---|---|---|
| R3-1 | RESOLVED | Violations are counted on the composed field and reported under owner policy. |
| R3-2 | RESOLVED | Publication permits numerical spacing tolerance only. |
| R3-3 | RESOLVED | Materializing bands prevents iterator exhaustion. |
| R3-4 | RESOLVED | Both-side probing and complete guide-segment checks address the reported failures. |
| R3-5 | PARTIAL | Ordinary rollback and recovery work; cleanup can still destroy an unrestored backup. Finding 7. |
| R3-6 | RESOLVED | Missing dependencies, early failures and the reported recipe-identity race are handled. |
| R3-7 | RESOLVED | Zero-weight masked neighbours preserve valid grid samples. |
| R3-8 | RESOLVED | fort.14 values round-trip and supplementary export is verified. |
| R3-9 | RESOLVED | Pairwise overlap checks remove global-area dilution. |
| R3-10 | RESOLVED | Notebook 414 selects its timestep after depth control. |
| R3-11 | RESOLVED | Boundary spacing controls must be finite and positive. |
| R3-12 | RESOLVED | Short designs receive explicit rejection. |
| R3-13 | RESOLVED | Filleting rejects reversals. |
| R3-14 | NOT RESOLVED | Same pending F27. |

**Previous findings — round 4**

| ID | Status | Reason |
|---|---|---|
| R4-F1 | NOT RESOLVED | Same pending F27. |
| R4-F2 | RESOLVED | Half-window amplitudes/phases derive from interpolated harmonic coefficients. |
| R4-F3 | RESOLVED | `_bilinear` preserves query shape; Cao now does too. |
| R4-F4 | RESOLVED | Overlap tolerance is tied to the smaller meeting element. |
| R4-F5 | RESOLVED | Complete guide segments are tested against land and constrained lines. |
| R4-F6 | PARTIAL | Recovery survives the previously injected faults, but another cleanup failure loses a backup. Finding 7. |
| R4-F7 | RESOLVED | State initialization and handler registration precede the first log call. |
| R4-F8 | RESOLVED | Parsed recipe/OBC bytes are hashed and checked by the stages. Scientific-source coverage has a separate gap: finding 8. |
| R4-F9 | RESOLVED | Cao retains valid centres beside missing cells and includes final rows/columns. |
| R4-F10 | RESOLVED | Boundary simplicity and repeated adjacent nodes are checked. |
| R4-F11 | RESOLVED | Manual timestep overrides must be finite and positive. |
| R4-F12 | RESOLVED | Run jobs resolve/check the executable before directory changes. |
| R4-F13 | RESOLVED | The guard requires ASCII and checks the 80-byte limit. Relocated cases have a separate gap: finding 11. |
| R4-F14 | RESOLVED | Spatial dimensions are checked, and unreadable named grids fail validation. |

**Previous findings — round 5**

| ID | Status | Reason |
|---|---|---|
| R5-F1 | NOT RESOLVED | Same pending F27. |
| R5-F2 | RESOLVED | All three flip sites use the convexity guard; rewritten ownership rows are skipped. |
| R5-F3 | RESOLVED | Hashes identify parsed bytes; drivers compare them and stages enforce expected identities. |
| R5-F4 | RESOLVED | Complete ladder segments receive geometric land/constrained-line checks. |
| R5-F5 | PARTIAL | Restoration attempts are independent, but cleanup can still delete an unrecovered backup. Finding 7. |
| R5-F6 | RESOLVED | An injected unreadable named grid returned `ok=False`. |
| R5-F7 | RESOLVED | Non-ASCII paths and paths over 80 bytes are refused. |
| R5-F8 | RESOLVED | Handler registration precedes the first post-reservation log call. |
| R5-F9 | RESOLVED | Cao restores query shape; its regression test passed. |
| R5-F10 | RESOLVED | Synthetic NaN and infinite depths each fail `min_depth_clip`. |

The generation-nondeterminism hypothesis is **WITHDRAWN** based on the supplied repeated-run evidence. Element-order dependence alone does not establish nondeterminism. The owner-authorized band-floor and timestep warnings are not findings.

**Findings**

1. **Major — Previously reported F27 remains open: package code imports oceanmesh.**  
   [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

   **Evidence:** Direct imports remain here, at line 165, in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672` and `autofinish/directives.py:20`. They conflict with the explicit repository rule prohibiting GPL imports within `fvcom_mesh_tools`. `THIRD_PARTY_NOTICES.md` still endorses the conflicting arrangement.

   **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary and reconcile the notices. **Previously reported; owner decision pending; counted once.**

2. **Major — Clipping land creates artificial coastlines that pass boundary-design validation.**  
   [444_design_obc.py:55](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:55).

   **Reproduction:** Use true land `box(-50000, 0, 50000, 40000)`, clipped to `box(-30000, -10000, 30000, 10000)`. Call `coast_normal` near `(0,11000)` and `(20000,11000)`, then design the rounded boundary through `(0,40000)` and `(20000,40000)`, with 5-km corner radii and 3-km spacing.

   The actual design functions and 444’s acceptance block produce **26 nodes, length 75.464 km, both endpoint angles 90°, zero reported land crossing and `problems=[]`**. Against the unclipped land, **the entire 75.464 km intersects land**. The clip edge is treated as a coastline, and checks use the same incomplete geometry. Neither design nor generation requires the footprint to remain inside the land window.

   **Fix:** Preserve genuine coastline segments, reject endpoint/chord neighbourhoods affected by clipping, and require complete land coverage of the design/generation footprint. Validate against coastline geometry without artificial clip edges.

3. **Major — Run validation ignores the namelist’s `OUTPUT_DIR` and can accept stale histories.**  
   [cli/check_run.py:185](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:185).

   **Reproduction:** Give the namelist `OUTPUT_DIR='current_output/'`, leave that directory empty, and provide a valid older history in `run/output/`. With a matching tiny grid, finite two-record NetCDF and a `TADA` log, the actual checker returned **`ok=True`, no reasons**, and searched only **`/virtual/run/output`**.

   Thus it can approve a run whose configured output contains no history. Conversely, valid runs using another output directory are rejected.

   **Fix:** Resolve `OUTPUT_DIR` relative to the run directory or as an absolute path, inspect that directory, and use a documented default only when the setting is absent.

4. **Major — The renumber benchmark accepts a zero-exit early STOP as a successful timing.**  
   [416_renumber_benchmark.sh:99](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/416_renumber_benchmark.sh:99).

   **Reproduction:** Execute the actual benchmark loop with mocked `mpiexec` printing `STOP: unable to initialize` and returning zero. It creates no histories and never prints `TADA`. The script’s error-pattern search misses that message; the loop emits successful **`BENCH … seconds=0.01`** records and exits **0**.

   The job checks neither simulation completion nor output health, so a failed initialization can appear exceptionally fast.

   **Fix:** Validate each invocation’s specific log and histories with `check_run` before accepting its timing. Require completed dates, finite fields and successful completion; propagate validation failures.

5. **Major — Speed verification exits successfully after an identity mismatch or failed QA.**  
   [407_speed_verify.sh:44](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/407_speed_verify.sh:44), [line 46](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/407_speed_verify.sh:46).

   **Reproduction:** Execute the actual identity-check tail with `cmp` returning 1 and the diagnostic Python command reporting a 10-m displacement successfully. The script prints **`DIFFERS from the reference`**, reaches its final timestamp and exits **0**. The preceding QA pipeline explicitly suppresses its failure with `|| true`.

   This contradicts the job’s stated requirement to produce an identical mesh.

   **Fix:** Preserve diagnostics but exit nonzero after a mismatch. Require QA success as well before reporting successful verification.

6. **Major — A relative speed-verification reference can become the generated file itself.**  
   [407_speed_verify.sh:23](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/407_speed_verify.sh:23).

   **Reproduction:** Set:

   ```bash
   FMESH_REF=outputs/sample_repro/sample_repro_final.14
   ```

   `REF` retains this relative string. After the job changes into its scratch directory, both `cmp` operands refer to the newly generated file there. Executing the actual tail with those operands prints **`IDENTICAL to the reference`** and exits **0**, without inspecting the repository reference.

   **Fix:** Resolve and verify `FMESH_REF` as an absolute existing file before changing directories or creating the work area.

7. **Minor — The new rollback cleanup can again delete an unrestored backup.**  
   [444_design_obc.py:270](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:270).

   **Reproduction:** Execute the actual publication block using virtual files, injecting three failures:

   - CSV publication fails.
   - Restoration of the old JSON fails.
   - Removal of the unused CSV backup fails once.

   Line 271 raises before `kept.clear()` and before writing `.RECOVER`. The `finally` block then deletes the unrestored JSON backup. The resulting set contains **old CSV, new JSON, old PNG, no old JSON backup and no recovery marker**; the publication lock is released.

   **Fix:** Immediately remove unrestored backups from ordinary cleanup tracking. Write recovery information before unrelated cleanup, and make cleanup operations independent so one failure cannot destroy recovery state. **Defect in the e2918dc recovery fix.**

8. **Minor — Completion-time provenance checks omit scientific source files.**  
   [445_extend_mesh.py:126](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:126), [453_redepth_extended.py:167](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:167).

   **Reproduction:** Execute 445’s actual orchestration/report block with virtual I/O. Its initial provenance hashes a bathymetry file containing `depth A`. Change that file persistently to `depth B` during the mocked generation stage. Both stages succeed.

   The driver reports **`status="ok"` and `inputs_changed_during_build=[]`**. Rechecking all captured provenance files correctly identifies **`bathymetry_srtm15plus`** as changed. The final check receives only `INPUTS`, excluding bathymetry and OSM data; 453 likewise excludes bathymetry.

   **Fix:** Recheck every captured scientific input, including dataset sidecars and inventories, or consume immutable snapshots. Reject changed scientific sources consistently with the existing moving-input policy.

9. **Minor — Native export accepts invalid connectivity and nonfinite mesh values.**  
   [io/fvcom_native.py:66](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:66).

   **Reproduction:** Export a three-node triangle using mocked file writers:

   ```python
   nodes = [[0, 0], [1000, 0], [0, 1000]]
   elements = [[-1, 0, 1]]
   depths = [5, 5, 5]
   ```

   NumPy interprets `-1` as the final node; orientation and unique-count checks pass. `export_fvcom_case` returns successfully after writing four products, including invalid grid connectivity **`1 0 1 2`**.

   Separate exports with a NaN coordinate or NaN depth also succeeded and wrote literal **`nan`** values. Extension read-back/QA protects its completed build, but the public exporter still accepts these invalid inputs.

   **Fix:** Before any writes, validate array shapes, integral connectivity within `[0, NP)`, and finite coordinates/depths.

10. **Minor — ASCII-valid directory names can corrupt the generated Fortran namelist.**  
    [383_m2_case_prep.py:178](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:178), [line 204](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:204).

    **Reproduction:** The actual path guard accepts `/tmp/run'case/input`. The namelist function emits:

    ```fortran
    INPUT_DIR= '/tmp/run'case/input/',
    ```

    The embedded apostrophe breaks the string. With a literal backslash followed by `t` in the directory name, `re.subn` interprets the replacement escape and writes a **tab instead of those pathname characters**.

    **Fix:** Serialize Fortran strings correctly, doubling embedded apostrophes and rejecting unsupported control characters. Use a callable regex replacement so pathname backslashes remain literal. This function is also used by 414 and 448.

11. **Minor — Relocated smoke and benchmark cases bypass the 80-byte directory guard.**  
    [416_renumber_benchmark.sh:57](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/416_renumber_benchmark.sh:57), [423_m2_smoke.sh:37](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/423_m2_smoke.sh:37).

    **Reproduction:** Use a run root consisting of `/tmp/` followed by 55 ASCII characters. The original refined output directory, including its trailing slash, is **76 bytes** and passes the actual guard.

    These jobs relocate it beneath `bench/` or `smoke/` and replace the namelist path directly. The new directory is **82 bytes**. The actual guard rejects that relocated path, but neither job invokes it before writing or launching FVCOM.

    **Fix:** Validate every relocated input/output directory before writing the copied namelist and launching the solver. Reuse the common namelist path validation.

## Verdict

VERDICT: FAIL (0 blocker, 6 major, 5 minor, 0 nit)

### Prompt

```markdown
# Review request, round 6: extending a base mesh outward (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (small synthetic inputs only; do not read the large data
under $DATA_DIR beyond listing it, and do not submit batch jobs). Answer in
English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
A tool that keeps a finished FVCOM base mesh exactly as it is and adds the
sea out to a new, designed open boundary (USER_GUIDE section 13). Read:

- `git show b0584f9 8e2739b 450ad44 d6d2a72 7044b6b 0f52d5b 69b50a4 d9e92fd b6d2ed8 765423c`
  (the extension tool and its documentation), and the current files:
  - `src/fvcom_mesh_tools/extend.py`, `extend_recipe.py`, `obc_design.py`,
    `dem/sources.py` (named bathymetry sources, priority stack, and the new
    `DATUM` registry / `non_tp_count` warning);
  - `notebooks/444_design_obc.py`, `445_extend_mesh.py`, `446_extend_generate.py`,
    `447_extend_merge.py`, `448_extend_smoke.py`, `453_redepth_extended.py`;
  - `recipes/extend/tokyo_bay_enshu.yaml`, `tokyo_bay_enshu_obc_design.yaml`;
  - `jobs/octopus/444_design_obc.sh`, `445_extend_mesh.sh`, `448_extend_smoke.sh`,
    `453_redepth_extended.sh`, `common.sh`;
  - tests: `tests/test_extend*.py`, `tests/test_obc_design*.py`,
    `tests/test_dem_sources.py` (whatever exists).
- Also in scope, just committed: the portability change --
  every job script and `common.sh` now take paths only from `$DATA_DIR` and
  `$WORK_DIR` (login profile), stop when they are unset, and derive the
  OCTOPUS FVCOM library directory as `FVCOM_LIBS` in `common.sh`; notebooks
  383/384/414 and `cli/refine_run.py` no longer fall back to `/octfs/...`.
  See commits 6d8b9a7 and 6c068d2 (`git log -5`).

Design intent:
- the base mesh's nodes, elements and depths are carried bit for bit
  (`verify_frozen_base`);
- the new part is generated with oceanmesh (subprocess; GPL code must never
  be imported by `fvcom_mesh_tools`), with fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-5 and their triage are in docs/extend-tools-review-20261001.md.
Round 5 (your previous answer; 10 findings) was fixed in e2918dc; read
`git show e2918dc`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 `algorithms/obc_finish.py`: new `_flip_pair` requires a strictly
  convex quadrilateral (both diagonals cross inside each other); used by
  all three flip sites (`flip_for_obc_perp`, `fix_r4`, `flip_c4_edges`).
  `flip_c4_edges` skips an edge whose rows an earlier flip rewrote.
- F3 `load_extend_recipe` hashes the bytes it parses (`recipe_sha256`,
  `open_boundary_sha256`); 445/453 compare them with the provenance hash;
  445 passes both to its stages via env, and 446/447 call `check_expected`.
- F4 446 `ladder`: each joined segment is tested whole against land
  (prepared `land_all`) and both constrained lines, plus the sdf samples.
- F5 444: restorations attempted one by one; unrestored backups kept with
  `<csv>.RECOVER`; 444 and the recipe loader refuse while it exists.
- F6 `check_run`: a named but unreadable grid is a failure reason.
- F7 `check_fvcom_dirs`: ASCII only, and <= 80 bytes.
- F8 445: STATE initialised and handler registered before the first log line
  (only a `def` lies between `reserve` and `atexit.register`).
- F9 `CaoNested.depth` keeps the query shape.
- F10 `min_depth_clip` fails non-finite depths (kept as one gate so the
  documented gate count is unchanged).
Numbering discrepancy: resolved by us. Two 446 runs of the same code gave
bit-identical meshes (job 122678), so generation is deterministic; the
change came from the round-4 interpolation rewrite moving sizing depths by
<= 1.8e-12 m. Your note that finishing is element-order dependent is
correct but, with deterministic generation, does not break reproducibility;
we do not treat it as a defect unless you show nondeterminism.
Verification after e2918dc: full test suite 1129 passed (batch
job 122691); on real data 444 published the Enshu boundary with no
problems; the full 445 build passed QA 23/23 with ladders fully kept and
joined, NP 14,738 / NE 27,133, and its grd/dep files are bit-identical to
the round-4 build (so no production flip was across a concave
quadrilateral); 453 re-depth completed.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step gate (447,
453) REPORT warnings and do not fail the build. Not a defect.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the scope and everything it touches.

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
| F1 (F27) | major | as round 1 | yes, pre-existing | owner decision pending |
| F2 | major | code read (land clipped to bbox, checks use the same land) | yes | fixed (window checks in 444 and 446), 3455b8e |
| F3 | major | code read (glob of run/output only) | yes | fixed, 3455b8e; test |
| F4 | major | code read (no completion check) | yes | fixed (check_run per run), 3455b8e |
| F5 | major | mock of the job tail (differs / QA fail -> exit 1) | yes | fixed, 3455b8e |
| F6 | major | code read (REF relative after cd) | yes | fixed, 3455b8e |
| F7 | minor | fault injection (3 failures) | yes | fixed, 3455b8e |
| F8 | minor | code read (INPUTS only) | yes | fixed (all provenance files), 3455b8e |
| F9 | minor | code read (no index/finite check) | yes | fixed, 3455b8e; tests |
| F10 | minor | namelist() with an apostrophe and a backslash | yes | fixed (io/fvcom_namelist), 3455b8e; tests |
| F11 | minor | code read (relocated path unchecked) | yes | fixed (relocate_case), 3455b8e; test |

The reviewer withdrew the generation-nondeterminism question on the
repeated-run evidence.

## Round 7

### Answer (verbatim)

Read-only review completed: **12 findings remain—2 major and 10 minor.** Finding 1 is the previously reported F27, still awaiting the owner’s decision.

No files were modified, shared datasets were read, or batch jobs were submitted. The existing untracked `recipes/extend/tokyo_bay_enshu_obc.png` was left untouched.

Verification:

- Capture-disabled pytest with `--noconftest -p no:cacheprovider`: **30 geometry tests passed**, plus **10 selected flip, DEM, recipe and namelist tests passed**.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` for every OCTOPUS shell script: **56 scripts, no failures**.
- Additional reproductions used small synthetic inputs, diskless NetCDF, mocked filesystem operations and shell commands.
- Ordinary pytest startup failed because capture requires writable temporary storage; its capture-disabled retry encountered Matplotlib’s writable-cache requirement through `conftest.py`. Bypassing that hook enabled the tests above.
- The full suite, production mesh generation and FVCOM integration were **not run**. No production QA counts or timestep measurements were independently obtained.

## Previous findings

Statuses assess the original defect and any remaining gap in its stated guarantee. References below point to this review’s numbered findings.

### Round 1

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked grid cells become NaN and permit fallback. |
| F2 | RESOLVED | Required QA failures produce a nonzero build exit. |
| F3 | RESOLVED | Limiting checks bounds, convergence and the final rounded r-factor. |
| F4 | RESOLVED | Floor and timestep warnings follow the explicit owner policy. |
| F5 | RESOLVED | Actual chords, endpoint floors and serialized geometry are checked. |
| F6 | RESOLVED | Missing boundary, lattice and generated-node sizing coverage is rejected. |
| F7 | RESOLVED | Continuous departures and the complete boundary are checked against land. |
| F8 | RESOLVED | Rejected designs preserve the existing CSV. Rollback gaps are tracked separately. |
| F9 | RESOLVED | Values round-trip and exported frozen-base verification runs. |
| F10 | RESOLVED | Opposite-side seam incidence and overlap beyond the seam are checked. |
| F11 | RESOLVED | Single-edge land runs are retained. |
| F12 | RESOLVED | Unsupported constrained lines shorter than six nodes are rejected explicitly. |
| F13 | RESOLVED | Cao cache keys include the resolved data root and dimensions. |
| F14 | RESOLVED | Invalid fine samples retain valid coarse coverage. |
| F15 | RESOLVED | Degenerate M7001 point sets return uncovered samples. Finding 2 is separate. |
| F16 | RESOLVED | Cao provenance includes existing depth files. Inventory changes remain separate. |
| F17 | RESOLVED | Relevant dirty sources and both rename paths are recorded. Finding 10 concerns later changes. |
| F18 | RESOLVED | Numeric settings, integer controls, seeds and bounds are validated. |
| F19 | RESOLVED | Resampling rejects nonfinite spacing and failure to progress. |
| F20 | RESOLVED | Extension outputs receive exclusive reservations. |
| F21 | RESOLVED | Re-depth rejects overlapping destinations and reserves its output. |
| F22 | RESOLVED | Smoke staging reserves a fresh root. |
| F23 | RESOLVED | Depth control precedes timestep selection. |
| F24 | RESOLVED | Smoke paths are resolved before use. |
| F25 | RESOLVED | Executable overrides are respected and resolved before directory changes. |
| F26 | RESOLVED | Case names match preparation and directory changes are guarded. |
| F27 | NOT RESOLVED | GPL imports still conflict with repository policy; owner decision pending. Finding 1. |
| F28 | RESOLVED | Extension case prefixes are restricted to filename components. |
| F29 | RESOLVED | Native reading validates unique, consecutive boundary traversal. |

### Round 2

| ID | Status | Reason |
|---|---|---|
| R2-1 | RESOLVED | Re-depth and smoke reserve outputs exclusively. |
| R2-2 | RESOLVED | Unique temporaries and a publication lock prevent competing-writer swaps. |
| R2-3 | RESOLVED | Re-depth repeats frozen-base, overlap and QA checks; overrides are labelled. |
| R2-4 | RESOLVED | Chords, endpoint floors and publication precision are checked. |
| R2-5 | RESOLVED | Continuous departure checks catch thin land strips. |
| R2-6 | RESOLVED | Land-crossing tolerance is limited to numerical endpoint residue. |
| R2-7 | RESOLVED | Checks use reconstructed published coordinates. |
| R2-8 | RESOLVED | Native coordinates and depths round-trip. |
| R2-9 | RESOLVED | Outer-versus-base overlap is checked beyond the seam. |
| R2-10 | RESOLVED | Porcelain parsing includes both rename paths. |
| R2-11 | RESOLVED | Geographic bounds receive finite, range and order validation. |
| R2-12 | RESOLVED | Resampling requires finite forward progress. |
| R2-13 | RESOLVED | Traversal direction follows the supplied second node. |
| R2-14 | RESOLVED | Convergence on the final permitted iteration is accepted. |
| R2-15 | RESOLVED | Limited-edge depths and limiter controls are validated. |
| R2-16 | RESOLVED | Both sides, guide points and complete guide segments are checked. |
| R2-17 | RESOLVED | Bands are materialized and follow the accepted deviation policy. |
| R2-18 | PARTIAL | Early failures retain provenance; consumption and late-reporting gaps remain. Findings 8–10, 12. |
| R2-19 | RESOLVED | Uncovered generated nodes are rejected. |
| R2-20 | RESOLVED | Invalid radii, zero sides and reversals are rejected. |
| R2-21 | RESOLVED | Refusing smoke-root reuse preserves previous success markers. |
| R2-22 | NOT RESOLVED | Same pending F27; counted once. |

### Round 3

| ID | Status | Reason |
|---|---|---|
| R3-1 | RESOLVED | Violations are counted on the composed field under the owner’s warning policy. |
| R3-2 | RESOLVED | Publication permits numerical spacing tolerance only. |
| R3-3 | RESOLVED | Materializing bands prevents iterator exhaustion. |
| R3-4 | RESOLVED | Both-side probing and complete guide-segment checks address the reported failures. |
| R3-5 | PARTIAL | Rollback preserves backups, but failed recovery-marker publication permits their later destruction. Finding 11. |
| R3-6 | PARTIAL | Early failure handling works; consumed-input identity and late failure status remain incomplete. Findings 8–10, 12. |
| R3-7 | RESOLVED | Zero-weight masked neighbours preserve valid samples. |
| R3-8 | RESOLVED | fort.14 values round-trip and supplementary export is verified. |
| R3-9 | RESOLVED | Pairwise overlap checks remove global-area dilution. |
| R3-10 | RESOLVED | Notebook 414 selects its timestep after depth control. |
| R3-11 | RESOLVED | Boundary spacing controls must be finite and positive. |
| R3-12 | RESOLVED | Short designs receive explicit rejection. |
| R3-13 | RESOLVED | Filleting rejects reversals. |
| R3-14 | NOT RESOLVED | Same pending F27. |

### Round 4

| ID | Status | Reason |
|---|---|---|
| R4-F1 | NOT RESOLVED | Same pending F27. |
| R4-F2 | RESOLVED | Half-window amplitudes/phases derive from interpolated harmonic coefficients. |
| R4-F3 | RESOLVED | Bilinear interpolation preserves query shape. |
| R4-F4 | RESOLVED | Overlap tolerance uses the smaller meeting element. |
| R4-F5 | RESOLVED | Complete guide segments are tested against land and constrained lines. |
| R4-F6 | PARTIAL | Ordinary restoration and the previous cleanup fault are handled; finding 11 remains. |
| R4-F7 | RESOLVED | State initialization and handler registration precede fallible logging. |
| R4-F8 | PARTIAL | Recipe identities are enforced, but actual OBC consumption remains unbound. Finding 8. |
| R4-F9 | RESOLVED | Cao retains valid centres beside missing cells and includes final rows/columns. |
| R4-F10 | RESOLVED | Boundary simplicity and repeated adjacent nodes are checked. |
| R4-F11 | RESOLVED | Manual timestep overrides must be finite and positive. |
| R4-F12 | RESOLVED | Run jobs resolve/check the executable before directory changes. |
| R4-F13 | RESOLVED | The guard requires printable ASCII and checks the 80-byte limit, including relocated cases. |
| R4-F14 | RESOLVED | Spatial dimensions are checked and unreadable named grids fail validation. |

### Round 5

| ID | Status | Reason |
|---|---|---|
| R5-F1 | NOT RESOLVED | Same pending F27. |
| R5-F2 | RESOLVED | All flip sites use the convexity guard and skip stale ownership rows. |
| R5-F3 | PARTIAL | YAML hashes identify parsed bytes; the boundary is still read separately after verification. Finding 8. |
| R5-F4 | RESOLVED | Complete ladder segments receive land/constrained-line checks. |
| R5-F5 | PARTIAL | Restoration and cleanup are independent, but marker-write failure remains unsafe. Finding 11. |
| R5-F6 | RESOLVED | An unreadable named grid fails run validation. |
| R5-F7 | RESOLVED | Non-ASCII paths and paths over 80 bytes are refused. |
| R5-F8 | RESOLVED | Handler registration precedes the first post-reservation log call. |
| R5-F9 | RESOLVED | Cao preserves query shape; its selected regression test passed. |
| R5-F10 | RESOLVED | Nonfinite depths fail the minimum-depth QA gate. |

### Round 6

| ID | Status | Reason |
|---|---|---|
| R6-F1 | NOT RESOLVED | Same pending F27. |
| R6-F2 | RESOLVED | 444 checks the buffered design inside the window; 446 checks sea inside its inset. |
| R6-F3 | PARTIAL | Ordinary relative/absolute `OUTPUT_DIR` works; valid escaped quotes select the wrong directory. Finding 3. |
| R6-F4 | RESOLVED | Each benchmark run is validated before timing acceptance, using saved conda Python without `LD_LIBRARY_PATH`. |
| R6-F5 | RESOLVED | QA failure and identity mismatch now produce failure exits. |
| R6-F6 | PARTIAL | Relative references resolve correctly, but self-reference rejection follows destructive staging. Finding 7. |
| R6-F7 | PARTIAL | The previous cleanup fault preserves backups; recovery-marker failure leaves another destructive path. Finding 11. |
| R6-F8 | PARTIAL | Existing provenance files are rehashed, but dynamically discovered inventories are not recollected. Finding 9. |
| R6-F9 | PARTIAL | Wrapper connectivity/finiteness checks work; boundary dtype and public-writer validation remain incomplete. Finding 4. |
| R6-F10 | RESOLVED | The writer doubles apostrophes and preserves literal backslashes. Reader incompatibility is finding 3. |
| R6-F11 | RESOLVED | Relocated directory lengths are checked. The new helper’s source deletion is finding 6. |

The generation-nondeterminism hypothesis remains **WITHDRAWN**, based on the supplied repeated-run evidence. Owner-authorized band-floor and new-element timestep warnings are **not defects**.

## Findings

1. **Major — Previously reported F27 remains open: package code imports oceanmesh.**  
   [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

   **Evidence:** Direct imports remain here and at line 165, in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672`, and `autofinish/directives.py:20`. They conflict with the explicit `CLAUDE.md` policy prohibiting GPL imports inside `fvcom_mesh_tools`. `THIRD_PARTY_NOTICES.md` still describes the conflicting arrangement.

   **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary and reconcile the notices. **Previously reported; owner decision pending; counted once.**

2. **Major — M7001 depth and source attribution depend on the other query points.**  
   [dem/sources.py:162](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:162).

   **Reproduction:** Mock the parquet reader with four `N` soundings at `(139.7,34.7)`, `(140.3,34.7)`, `(139.7,35.3)`, `(140.3,35.3)`, all `z_tp=-10`, and add a constant 100-m fallback source. Actual `sample()` results:

   ```text
   query [(140,35)]:
       depth=[100], which=[fallback]

   query [(140,35), (139.7,34.7), (140.3,35.3)]:
       depth=[10,10,10], which=[m7001,m7001,m7001]
   ```

   The unchanged centre lies inside the soundings’ hull. The query-derived 0.25° window excludes all soundings for the first call, while the second call includes them. The priority stack’s uncovered subset also changes this window, so earlier-source coverage can alter M7001 results at other points. Scalar design queries and batch sizing/depth queries can therefore disagree.

   **Fix:** Define interpolation and coverage from a stable dataset triangulation, or a deterministic per-location neighbourhood independent of query batching. Add batch-invariance tests. Pre-existing, newly identified.

3. **Minor — The new `OUTPUT_DIR` handling misreads valid Fortran escaped quotes.**  
   [cli/check_run.py:75](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:75), [cli/check_run.py:168](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:168).

   **Reproduction:** The new writer correctly emits:

   ```fortran
   OUTPUT_DIR = '/virtual/current''case/',
   ```

   `_nml_value()` returns `/virtual/current`. With no history in the actual apostrophe-containing directory, a valid two-record diskless NetCDF in the truncated directory, matching grid counts and a `TADA` log, the checker returned **`ok=True`, `output_dir='/virtual/current'`, `reasons=[]`**. Quoted `INPUT_DIR` values are likewise truncated.

   **Fix:** Use a namelist reader that handles doubled quote delimiters and unescapes them. Test writer/reader round trips for both quote styles. This is an interaction between the R6-F3 and R6-F10 fixes.

4. **Minor — Export validation still silently changes fractional boundary IDs, and standalone writers bypass it.**  
   [io/fvcom_native.py:323](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:323), [io/fvcom_native.py:87](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:87).

   **Reproduction:** For a valid three-node CCW mesh, an open boundary `[0.9,1.9]` passes `_check_exportable()`. Actual `export_fvcom_case(..., obc_depth_control=False)` succeeds and writes node IDs **1 and 2**, silently truncating the supplied indices.

   Separately, actual public `write_grd()` accepts connectivity `[[-1,0,1]]` on the same mesh: NumPy treats `-1` as the final node during validation, and the writer emits **`1 0 1 2`**, including invalid FVCOM node zero. The new wrapper guard is never called by this public writer.

   **Fix:** Centralize validation across public writers. Require boundary arrays to be one-dimensional with integral, finite, in-range indices before any cast; apply connectivity and finiteness checks to standalone exports too. R6-F9 is incomplete.

5. **Minor — Optional export controls accept nonfinite values and are validated after other files are overwritten.**  
   [io/fvcom_native.py:158](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:158), [io/fvcom_native.py:183](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:183), [io/fvcom_native.py:369](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:369).

   **Reproduction:** Actual export with `cor=[nan,35,35]` and `sponge=[(0.9,nan,inf)]` succeeds, emitting a `nan` Coriolis value and sponge row **`1 nan inf`**.

   A separate call with wrong-length `cor=[35]` raises only after `_grd.dat`, `_dep.dat` and `_obc.dat` have been written. Against an existing destination, invalid optional input therefore leaves a partially replaced case.

   **Fix:** Prevalidate all optional controls before opening outputs: shape, finiteness, integral node IDs and valid radius/damping ranges. Stage the complete output set before replacing an existing case. Pre-existing, newly identified.

6. **Minor — The new `relocate_case()` helper can delete its source.**  
   [io/fvcom_namelist.py:62](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_namelist.py:62).

   **Reproduction:** With filesystem operations mocked, calling `relocate_case('/tmp/staged', '/tmp/staged')` records **`shutil.rmtree('/tmp/staged')`**, then `copytree()` fails because the source was deleted. A destination symlink resolving to the source has the same result. A destination that is an ancestor of the source also deletes it.

   **Fix:** Reject equal or overlapping resolved source/destination paths before mutation. Render and validate the rewritten namelist before replacing the destination, and publish a staged copy safely. **Introduced by the new R6-F11 helper.**

7. **Minor — Speed verification rejects a self-reference after deleting it.**  
   [407_speed_verify.sh:34](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/407_speed_verify.sh:34), [407_speed_verify.sh:52](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/407_speed_verify.sh:52).

   **Reproduction:** Run the actual staging and guard statements with shell commands mocked and `FMESH_REF` equal to the existing `$WORK/outputs/sample_repro/sample_repro_final.14`. They first execute **`rm -f` on the reference**, and only later reject it with exit 2. Preceding `rsync` can also replace its contents.

   **Fix:** Calculate the produced path and compare resolved references before any staging, copying, removal or generation. This is an incomplete R6-F6 fix.

8. **Minor — The boundary bytes actually used by generation are not bound to the verified digest.**  
   [extend_recipe.py:80](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend_recipe.py:80), [446_extend_generate.py:71](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:71).

   **Reproduction:** Execute the actual recipe-loading, expected-identity check and OBC assignment with mocked reads. The loader validates CSV A and hashes A; `check_expected()` accepts A. The later `read_open_boundary()` returns CSV B, with latitudes changed from **34.0 to 34.1**. Generation’s OBC then contains B while `open_boundary_sha256` still identifies A. Restoring A before completion also defeats the final file check.

   The loader itself validates and hashes the CSV in separate reads, unlike its YAML handling.

   **Fix:** Read the boundary once, hash those exact bytes, parse them and return the verified coordinates for stage consumption. Alternatively use immutable input snapshots. R4-F8/R5-F3 remain partial.

9. **Minor — Completion checks miss newly discovered scientific input files.**  
   [provenance.py:263](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/provenance.py:263), [dem/sources.py:193](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:193).

   **Reproduction:** Collect Cao provenance when the mocked inventory contains an unchanged area table. Add a depth file for an area already described by that table before sampling. `CaoNested.files()` now includes the new depth file, but actual `changed_files(PROV, list(PROV["files"]))` returns **`[]`**: it rehashes only the originally captured paths.

   Sampling discovers depth files dynamically, so a stage can consume the new file without recording it or failing the completion check. The same inventory problem applies to newly appearing dataset sidecars.

   **Fix:** Recollect scientific inventories at completion and compare path membership as well as content. Prefer recording exact consumed files or using snapshots. The expanded R6-F8 check remains incomplete.

10. **Minor — A build can accept changed execution code while reporting the initial code identity.**  
    [445_extend_mesh.py:120](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:120), [445_extend_mesh.py:128](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:128).

    **Evidence/reproduction:** The driver captures code identity before launching subprocesses, but both stages execute scripts and import packages from the live repositories. Completion checks only `PROV["files"]`; no code identity is rechecked.

    Mocked execution of the actual driver tail with successful stages, unchanged data and simulated live code changing from A to B produced **`status="ok"`, `inputs_changed=[]`, provenance code=A**. A real edit between stages can make the second stage execute B. Re-depth has the same omission.

    **Fix:** Execute an immutable code snapshot, or verify relevant package, driver, stage and oceanmesh identities at stage boundaries and completion. This concerns changes after capture, independently of the resolved initial dirty-code recording. Pre-existing, newly identified.

11. **Minor — Failure to write the recovery marker releases the lock and lets the next publication destroy the retained backup.**  
    [444_design_obc.py:280](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:280), [444_design_obc.py:297](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:297).

    **Reproduction:** Execute the actual publication/rollback block against a virtual filesystem. Inject CSV replacement failure, JSON restoration failure, then `.RECOVER` write failure. The revised cleanup preserves the old JSON backup, but **both the lock and recovery marker are absent** afterward. The published JSON remains new.

    Execute the next publication’s backup loop: `shutil.copy2()` overwrites that retained `.json.prev` with the current new JSON, destroying the old report. Observed backup contents changed from **old JSON to new JSON**.

    **Fix:** Retain the exclusive lock when recovery cannot be recorded, or persist recovery intent before publication begins. Refuse to overwrite unresolved backups. The previously reported cleanup fault is fixed, but the recovery guarantee remains partial.

12. **Minor — Re-depth leaves an overall success report when figure publication fails.**  
    [453_redepth_extended.py:193](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:193), [453_redepth_extended.py:229](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:229).

    **Reproduction:** Execute the actual report-write, completion-state assignment and figure-save statements with `fig.savefig()` raising `OSError`. The command fails without its PNG, but **`redepth.json` retains `status="ok"` and `STATE["done"]` is true**, disabling the failure handler. Matplotlib import or plotting errors take the same path.

    **Fix:** Finalize overall success only after required products are published. If plotting is optional, record its failure separately and define the command’s exit/status consistently. Pre-existing, newly identified.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 10 minor, 0 nit)

### Prompt

```markdown
# Review request, round 7: extending a base mesh outward (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (small synthetic inputs only; do not read the large data
under $DATA_DIR beyond listing it, and do not submit batch jobs). Answer in
English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
A tool that keeps a finished FVCOM base mesh exactly as it is and adds the
sea out to a new, designed open boundary (USER_GUIDE section 13). Read:

- `git show b0584f9 8e2739b 450ad44 d6d2a72 7044b6b 0f52d5b 69b50a4 d9e92fd b6d2ed8 765423c`
  (the extension tool and its documentation), and the current files:
  - `src/fvcom_mesh_tools/extend.py`, `extend_recipe.py`, `obc_design.py`,
    `dem/sources.py` (named bathymetry sources, priority stack, and the new
    `DATUM` registry / `non_tp_count` warning);
  - `notebooks/444_design_obc.py`, `445_extend_mesh.py`, `446_extend_generate.py`,
    `447_extend_merge.py`, `448_extend_smoke.py`, `453_redepth_extended.py`;
  - `recipes/extend/tokyo_bay_enshu.yaml`, `tokyo_bay_enshu_obc_design.yaml`;
  - `jobs/octopus/444_design_obc.sh`, `445_extend_mesh.sh`, `448_extend_smoke.sh`,
    `453_redepth_extended.sh`, `common.sh`;
  - tests: `tests/test_extend*.py`, `tests/test_obc_design*.py`,
    `tests/test_dem_sources.py` (whatever exists).
- Also in scope, just committed: the portability change --
  every job script and `common.sh` now take paths only from `$DATA_DIR` and
  `$WORK_DIR` (login profile), stop when they are unset, and derive the
  OCTOPUS FVCOM library directory as `FVCOM_LIBS` in `common.sh`; notebooks
  383/384/414 and `cli/refine_run.py` no longer fall back to `/octfs/...`.
  See commits 6d8b9a7 and 6c068d2 (`git log -5`).

Design intent:
- the base mesh's nodes, elements and depths are carried bit for bit
  (`verify_frozen_base`);
- the new part is generated with oceanmesh (subprocess; GPL code must never
  be imported by `fvcom_mesh_tools`), with fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-6 and their triage are in docs/extend-tools-review-20261001.md.
Round 6 (your previous answer; 11 findings) was fixed in 3455b8e; read
`git show 3455b8e`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 444: the boundary buffered by (chord_m + 5 km) must lie within the
  land window (metric transform of the segmentized bbox); 446: the sea to
  mesh (domain polygon minus land) must stay 0.02 deg inside land.bbox.
- F3 check_run looks in the namelist's OUTPUT_DIR (relative to the run or
  absolute), else <run>/output; `info["output_dir"]` reports it.
- F4 416: each run is validated with `fvcom_mesh_tools.cli.check_run`
  (conda python kept as PYBIN, LD_LIBRARY_PATH unset) before its time.
- F5/F6 407: REF resolved with `realpath -e` before cd and must not equal
  the produced file; QA failure or identity mismatch exits 1.
- F7 444: marker written and unrestored backups popped from cleanup before
  any other fallible step; each cleanup unlink independent.
- F8 445/453: completion check over every provenance file.
- F9 `export_fvcom_case` validates indices, dtype and finiteness first.
- F10/F11 new `io/fvcom_namelist.py` (fortran_string, set_value with a
  callable, check_fvcom_dirs, relocate_case), used by 383/414/448 and the
  416/423 jobs.
Verification after 3455b8e: full test suite 1138 passed (batch
job 122760); on real data 444 published the Enshu boundary with no
problems (window check passed), the full 445 build passed QA 23/23 (window
check passed), its grd file bit-identical to rounds 4-5, and 453 completed.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step gate (447,
453) REPORT warnings and do not fail the build. Not a defect.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the scope and everything it touches.

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
| F1 (F27) | major | as round 1 | yes, pre-existing | owner decision pending |
| F2 | major | test: centre of four soundings alone vs in a batch | yes | fixed (whole-dataset triangulation), 227536f; test |
| F3 | minor | round-trip test of fortran_string | yes | fixed, 227536f; test |
| F4 | minor | code read (cast before check; write_grd bypass) | yes | fixed, 227536f; test |
| F5 | minor | code read (cor checked after writes) | yes | fixed (prevalidate, staged export), 227536f; test |
| F6 | minor | code read (rmtree before copy) | yes | fixed, 227536f; test |
| F7 | minor | code read (rm -f before guard) | yes | fixed; mock of the job head, 227536f |
| F8 | minor | code read (two reads of the CSV) | yes | fixed (open_boundary_lonlat), 227536f |
| F9 | minor | code read (paths fixed at collect) | yes | fixed (changed_inventory), 227536f; test |
| F10 | minor | code read (no code recheck) | yes | fixed, 227536f |
| F11 | minor | fault injection (marker write fails, then republish) | yes | fixed, 227536f |
| F12 | minor | code read (report before figure) | yes | fixed, 227536f |

## Round 8

### Answer (verbatim)

The review found **one previously reported major issue and 12 additional minor defects**. No files were modified.

Statuses below assess the original reproductions. New defects in related code are listed separately.

**Previous findings — round 1**

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked grid cells become NaN, allowing fallback. |
| F2 | RESOLVED | Required QA failures produce failure exits. |
| F3 | RESOLVED | Limiting checks bounds, convergence and rounded r-factor. |
| F4 | RESOLVED | Warning behavior follows the explicit owner policy. |
| F5 | RESOLVED | Actual chords, endpoint floors and published coordinates are checked. |
| F6 | RESOLVED | Missing boundary, lattice and generated-node coverage is rejected. |
| F7 | RESOLVED | Continuous departures and complete boundaries receive land checks. |
| F8 | RESOLVED | Rejected designs preserve existing CSVs. |
| F9 | RESOLVED | Native values round-trip; exported base verification runs. Finding 12 concerns the verifier’s bit comparison. |
| F10 | RESOLVED | Seam incidence and overlap beyond the seam are checked. |
| F11 | RESOLVED | Single-edge land runs are retained. |
| F12 | RESOLVED | Unsupported constrained lines shorter than six nodes are explicitly rejected. |
| F13 | RESOLVED | Cao cache keys distinguish data roots and dimensions. Finding 9 concerns subsequent file changes. |
| F14 | RESOLVED | Invalid fine samples retain valid coarse coverage. |
| F15 | RESOLVED | Degenerate M7001 point sets return uncovered samples. |
| F16 | RESOLVED | Cao provenance includes existing depth files; inventories are recollected. |
| F17 | RESOLVED | Relevant dirty sources and both rename paths are recorded. |
| F18 | RESOLVED | Numeric controls, integer settings, seeds and bounds are validated. |
| F19 | RESOLVED | Resampling rejects nonfinite spacing and failure to progress. |
| F20 | RESOLVED | Extension outputs receive exclusive reservations. |
| F21 | RESOLVED | Re-depth rejects overlapping destinations and reserves outputs. |
| F22 | RESOLVED | Smoke staging reserves a fresh root. |
| F23 | RESOLVED | Depth control precedes timestep selection. |
| F24 | RESOLVED | Smoke paths are resolved before use. |
| F25 | RESOLVED | Executable overrides are respected and resolved before directory changes. |
| F26 | RESOLVED | Job 383 uses the prepared case names and guards directory changes. Finding 13 concerns job 413. |
| F27 | NOT RESOLVED | GPL imports remain; owner decision pending. Counted once as finding 1. |
| F28 | RESOLVED | Case prefixes are restricted to filename components. |
| F29 | RESOLVED | Native reading validates unique, consecutive boundary traversal. |

**Previous findings — round 2**

| ID | Status | Reason |
|---|---|---|
| R2-1 | RESOLVED | Re-depth and smoke reserve outputs exclusively. |
| R2-2 | RESOLVED | Unique temporaries and publication locking prevent competing swaps. |
| R2-3 | RESOLVED | Re-depth repeats frozen-base, overlap and QA checks; overrides are labelled. |
| R2-4 | RESOLVED | Chords, endpoint floors and publication precision are checked. |
| R2-5 | RESOLVED | Continuous departure checks catch thin land strips. |
| R2-6 | RESOLVED | Crossing tolerance is limited to numerical endpoint residue. |
| R2-7 | RESOLVED | Checks use reconstructed published coordinates. |
| R2-8 | RESOLVED | Native coordinates and depths round-trip. |
| R2-9 | RESOLVED | Outer/base overlap is checked beyond the seam. |
| R2-10 | RESOLVED | Porcelain parsing includes both rename paths. |
| R2-11 | RESOLVED | Geographic bounds receive finite, range and order checks. |
| R2-12 | RESOLVED | Resampling requires finite forward progress. |
| R2-13 | RESOLVED | Traversal direction follows the supplied second node. |
| R2-14 | RESOLVED | Final-iteration convergence is accepted. |
| R2-15 | RESOLVED | Limited-edge depths and limiter controls are validated. |
| R2-16 | RESOLVED | Both sides and complete guide segments are checked. |
| R2-17 | RESOLVED | Bands are materialized and follow the owner-approved deviation policy. |
| R2-18 | RESOLVED | Early provenance, consumed-input identity, completion inventories/code and late reporting are now handled. |
| R2-19 | RESOLVED | Uncovered generated nodes are rejected. |
| R2-20 | RESOLVED | Invalid radii, zero sides and reversals are rejected. |
| R2-21 | RESOLVED | Refusing smoke-root reuse preserves previous success markers. |
| R2-22 | NOT RESOLVED | Same pending F27; counted once. |

**Previous findings — round 3**

| ID | Status | Reason |
|---|---|---|
| R3-1 | RESOLVED | Composed-field violations are reported under the owner’s warning policy. |
| R3-2 | RESOLVED | Publication permits numerical spacing tolerance only. |
| R3-3 | RESOLVED | Materializing bands prevents iterator exhaustion. |
| R3-4 | RESOLVED | Both-side probing and complete guide checks cover the reproductions. |
| R3-5 | RESOLVED | Failed marker publication retains the lock; existing backups block publication. |
| R3-6 | RESOLVED | Early handling, consumed-input identity and late failure reporting are fixed. |
| R3-7 | RESOLVED | Zero-weight masked neighbours preserve valid samples. |
| R3-8 | RESOLVED | fort.14 values round-trip; supplementary export is verified. |
| R3-9 | RESOLVED | Pairwise overlap checks avoid global-area dilution. |
| R3-10 | RESOLVED | Notebook 414 selects timestep after depth control. |
| R3-11 | RESOLVED | Spacing controls must be finite and positive. |
| R3-12 | RESOLVED | Short designs receive explicit rejection. |
| R3-13 | RESOLVED | Filleting rejects reversals. |
| R3-14 | NOT RESOLVED | Same pending F27. |

**Previous findings — round 4**

| ID | Status | Reason |
|---|---|---|
| R4-F1 | NOT RESOLVED | Same pending F27. |
| R4-F2 | RESOLVED | Half-window amplitudes/phases derive from interpolated harmonic coefficients. |
| R4-F3 | RESOLVED | Bilinear interpolation preserves query shape. |
| R4-F4 | RESOLVED | Overlap tolerance uses the smaller meeting element. |
| R4-F5 | RESOLVED | Complete guide segments receive land/constrained-line checks. |
| R4-F6 | RESOLVED | Restoration, recovery locking and backup preservation address the reported faults. |
| R4-F7 | RESOLVED | State and failure handlers precede fallible logging. |
| R4-F8 | RESOLVED | CSV coordinates and their identity now come from the same bytes. |
| R4-F9 | RESOLVED | Cao retains valid centres beside missing cells and includes final rows/columns. |
| R4-F10 | RESOLVED | Boundary simplicity and repeated adjacent nodes are checked. |
| R4-F11 | RESOLVED | Timestep overrides must be finite and positive. |
| R4-F12 | RESOLVED | Run jobs resolve/check executables before directory changes. |
| R4-F13 | RESOLVED | Printable ASCII and the 80-byte directory limit are enforced. |
| R4-F14 | RESOLVED | History dimensions are checked; unreadable named grids fail validation. |

**Previous findings — round 5**

| ID | Status | Reason |
|---|---|---|
| R5-F1 | NOT RESOLVED | Same pending F27. |
| R5-F2 | RESOLVED | Flip sites apply convexity guards and skip stale ownership rows. |
| R5-F3 | RESOLVED | YAML and CSV identities identify the bytes actually parsed. |
| R5-F4 | RESOLVED | Complete ladder segments receive geometry checks. |
| R5-F5 | RESOLVED | Marker failure retains recovery protection; backups cannot be overwritten. |
| R5-F6 | RESOLVED | An unreadable named grid fails run validation. |
| R5-F7 | RESOLVED | Non-ASCII and overlength paths are refused. |
| R5-F8 | RESOLVED | Handler registration precedes initial post-reservation logging. |
| R5-F9 | RESOLVED | Cao preserves query shape. |
| R5-F10 | RESOLVED | Nonfinite depths fail minimum-depth QA. |

**Previous findings — round 6**

| ID | Status | Reason |
|---|---|---|
| R6-F1 | NOT RESOLVED | Same pending F27. |
| R6-F2 | RESOLVED | Buffered design and sea geometry are checked against their windows. |
| R6-F3 | RESOLVED | Ordinary paths and doubled quoted delimiters are parsed correctly. Finding 8 is a different parsing defect. |
| R6-F4 | RESOLVED | Every benchmark run is validated before timing acceptance. |
| R6-F5 | RESOLVED | QA failures and identity mismatches produce failure exits. |
| R6-F6 | RESOLVED | Reference/work-area checks now precede staging. |
| R6-F7 | RESOLVED | Recovery-marker failure retains locking and protects backups. |
| R6-F8 | RESOLVED | Completion checks rehash files and recollect inventories. |
| R6-F9 | RESOLVED | Original fractional-index and writer-bypass reproductions are rejected. Finding 4 concerns floating connectivity. |
| R6-F10 | RESOLVED | Writer/reader agreement covers doubled delimiters and literal backslashes. |
| R6-F11 | RESOLVED | Relocation checks directory lengths and source/destination overlap. |

**Previous findings — round 7**

| ID | Status | Reason |
|---|---|---|
| R7-F1 | NOT RESOLVED | Same pending F27. |
| R7-F2 | RESOLVED | M7001 uses a cached whole-dataset triangulation; query composition no longer changes results. |
| R7-F3 | RESOLVED | `_nml_value` decodes doubled delimiters. |
| R7-F4 | RESOLVED | Original fractional, nonfinite and out-of-range index cases are rejected across writers. |
| R7-F5 | PARTIAL | Invalid optional inputs are rejected before publication; publication failures can still leave a mixed case. Finding 2. |
| R7-F6 | RESOLVED | Resolved source overlaps are refused; rendering and copying precede replacement. Finding 3 concerns destination recovery. |
| R7-F7 | RESOLVED | Job 407 checks its work area before staging. |
| R7-F8 | RESOLVED | CSV coordinates/hash share one read; 446 consumes those coordinates. |
| R7-F9 | RESOLVED | Completion inventories are recollected. |
| R7-F10 | RESOLVED | Completion code collection detects changes. Finding 7 concerns excessive comparison scope. |
| R7-F11 | RESOLVED | Lock survives marker failure; existing `.prev` blocks publication. |
| R7-F12 | RESOLVED | Figure generation precedes the success report. |

The generation-nondeterminism hypothesis remains **WITHDRAWN**, based on the supplied repeated-run evidence. Owner-authorized band-floor and new-element timestep warnings are not findings.

**Findings**

1. **Major — Previously reported F27 remains open: package code imports oceanmesh.**  
   [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

   **Evidence:** Direct imports remain here and at line 165, in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672`, and `autofinish/directives.py:20`. This conflicts with the explicit repository policy prohibiting GPL imports inside `fvcom_mesh_tools`. The extension subprocess does not eliminate those other imports.

   **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary and reconcile `THIRD_PARTY_NOTICES.md`. **Existing finding; owner decision pending; counted once.**

2. **Minor — Export staging does not protect an existing case during publication failure.**  
   [io/fvcom_native.py:435](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:435).

   **Reproduction:** With an existing case, mock the second `Path.replace` to raise `OSError`. The in-memory filesystem then contains the **new grid, old depths and old OBC**. The original grid has been overwritten, and `finally` deletes the remaining staged files.

   **Fix:** Preserve originals and restore them on publication failure, retaining recovery artifacts if restoration fails; alternatively publish a complete versioned directory through one atomic switch. The new staging fix resolves optional-input failures but leaves this publication gap.

3. **Minor — Relocation destroys the previous destination if final publication fails.**  
   [io/fvcom_namelist.py:86](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_namelist.py:86).

   **Reproduction:** Mock `work.rename(dst)` to raise after a successful staged copy. `shutil.rmtree(dst)` has already deleted the existing destination; the `finally` block also removes the staged replacement. The in-memory reproduction confirmed that the old destination no longer exists.

   **Fix:** Reject existing destinations, or move the old destination into a recoverable backup and restore it if publication fails. Source-overlap protection does not protect existing destination data.

4. **Minor — The new validator accepts floating connectivity that export subsequently cannot index.**  
   [io/fvcom_native.py:335](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:335).

   **Reproduction:** A valid CCW triangle with `elements=np.array([[0., 1., 2.]])` passes `_check_exportable`. `_validate_for_export` then raises `IndexError: arrays used as indices must be of integer (or boolean) type`, because `_signed_areas` indexes with the original floating array.

   **Fix:** Either explicitly require integer connectivity or propagate the normalized integer connectivity returned by `_indices` into every subsequent operation. This is a regression from the former explicit dtype rejection.

5. **Minor — Sponge prevalidation consumes an iterator, then export successfully writes an empty sponge.**  
   [io/fvcom_native.py:411](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:411).

   **Reproduction:** Pass `sponge=iter([(0, 500., .001)])`. `_check_sponge` consumes it, but export discards the normalized rows and passes the exhausted iterator to `write_spg` at line 429. The mocked export completed with zero sponge rows.

   **Fix:** Normalize once and pass the returned rows onward. Although the annotation specifies `Sequence`, an unsupported iterator must be rejected explicitly if it is not supported; accepting it and silently dropping its contents is unsafe. This double consumption was introduced by the fix.

6. **Minor — OBC types silently truncate fractional values.**  
   [io/fvcom_native.py:133](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:133).

   **Reproduction:** `write_obc(..., obc_type=[1.9])` succeeds and writes type `1` for both nodes in a two-node segment. It silently changes the requested boundary condition. Scalar `True` also follows the integer branch and is formatted as `True`.

   **Fix:** Validate finite, whole, non-boolean types before conversion, then enforce the FVCOM range. Pre-existing.

7. **Minor — Completion code checks reject unchanged code when unrelated repository files change.**  
   [445_extend_mesh.py:135](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:135); also `453_redepth_extended.py:177`.

   **Reproduction:** Keep the commit, package, driver, version and source bytes unchanged; add `custom_result/t_grd.dat` to the repository-wide dirty inventory. Both snapshots report `commit_identifies_code=True`, but their code dictionaries compare unequal. A permitted output directory outside ignored `outputs/` can therefore make the build reject its own products as a code change. An unrelated documentation edit has the same effect.

   **Fix:** Compare relevant source identities, rather than the complete repository dirty inventory retained by `provenance.code_state`. Keep unrelated dirty paths as diagnostics. The completion comparison introduced this failure path.

8. **Minor — Namelist parsing can select assignments embedded inside another quoted value and accept stale history.**  
   [cli/check_run.py:79](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:79).

   **Reproduction:** Use this valid namelist fragment:

   ```fortran
   &NML_CASE
     CASE_TITLE = "OUTPUT_DIR = '/virtual/stale/'",
   /
   &NML_IO
     OUTPUT_DIR = '/virtual/current/',
   /
   ```

   `_nml_value(..., "OUTPUT_DIR")` returns `/virtual/stale/`. In the complete mocked checker reproduction, the current directory had no history, while matching synthetic stale history produced `ok=True` with no reasons.

   **Fix:** Tokenize namelist syntax while respecting quoted values and group boundaries. Search actual assignments in the appropriate group; reject ambiguous duplicates. Pre-existing parser limitation survives the delimiter fix.

9. **Minor — Cao caches stale depths indefinitely after a file changes.**  
   [dem/sources.py:244](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:244).

   **Reproduction:** Mock one depth file first containing four 10-m values, then four 20-m values. Two calls using the same Cao instance/root/area/shape both return 10 m; `read_bytes` is called only once. Module-level source instances persist across API calls.

   **Fix:** Include the resolved file identity and change metadata in the cache key, as M7001 does, or give the cache an explicit run lifetime. Otherwise a later run can use old samples while recording the updated file’s provenance. Pre-existing; root isolation is already fixed.

10. **Minor — Grid sampling assumes latitude/longitude dimension order and can silently transpose depths.**  
    [dem/sources.py:143](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:143).

    **Reproduction:** A synthetic NetCDF variable with dimensions `("lon", "lat")`, two coordinates on each axis, and elevations `[[-10,-20],[-30,-40]]` returns depth **30 m** at the first longitude/second latitude; the correct depth is **20 m**. Equal axis lengths conceal the shape error.

    **Fix:** Inspect variable dimension names, validate coordinate dimensions and shape, and transpose supported layouts or reject unsupported ones clearly. Pre-existing; no claim that the current production datasets use this layout.

11. **Minor — Sponge serialization can turn validated positive radii into zero.**  
    [io/fvcom_native.py:187](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:187).

    **Reproduction:** `write_spg(..., [(0, 1e-5, .001)])` passes validation but writes:

    ```text
    Sponge Node Number = 1
    1 0.0000 0.001000
    ```

    FVCOM’s `mod_setup.F:551` divides distance by `R_SPG`. Fixed decimal formatting also rounds sufficiently small positive damping values to zero.

    **Fix:** Use round-trip-safe numeric formatting for radius and damping, or validate their serialized values before publication. Pre-existing; the new positivity check does not ensure positive exported radii.

12. **Minor — Frozen-base verification checks numerical equality rather than bit equality.**  
    [extend.py:174](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:174).

    **Reproduction:** In a valid synthetic merged mesh, change a base coordinate from `-0.0` to `+0.0`. `verify_frozen_base` passes despite different floating-point bit patterns. `np.array_equal` also does not enforce matching dtypes.

    **Fix:** Check dtype and compare the preserved arrays’ bit representations, including depths at line 178. This substantiates a verifier gap against the explicit “bit for bit” contract; it does not demonstrate a production mesh alteration.

13. **Minor — The M2 dependency chain’s analysis job checks nonexistent case names.**  
    [jobs/octopus/413_m2_analysis.sh:21](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/413_m2_analysis.sh:21).

    **Evidence:** Jobs 410/411 prepare and run `A`, `B_own`, and `B_m7001`. Job 413 instead requires `base/RUN_OK` and `refined/RUN_OK`. A shell fault injection supplying the manifest and all three successful chain markers exits **2**, reporting the missing `base/RUN_OK`.

    **Fix:** Read required case names from the manifest, or use the chain’s actual three names. Pre-existing and separate from the repaired job-383 finding.

**Verification and limits**

- `pytest --noconftest -s -q -p no:cacheprovider tests/test_extend.py tests/test_obc_design.py`: **30 passed**.
- Direct diskless execution of zero-argument tests across the reviewed modules: **45 passed**, including those geometry tests.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` on all **56** OCTOPUS scripts: **passed**.
- Small synthetic and mocked reproductions established the findings above.
- Ordinary pytest startup failed because writable temporary storage was unavailable. Disabling capture then encountered Matplotlib’s cache requirement through `conftest.py`; bypassing that hook enabled geometry tests.
- The full suite, production mesh generation and FVCOM integration were **not run**. No production node/element counts, QA gates or timestep measurements were independently obtained. The supplied 1144-test and real-data results remain prior verification evidence.

## Verdict

VERDICT: FAIL (0 blocker, 1 major, 12 minor, 0 nit)

### Prompt

```markdown
# Review request, round 8: extending a base mesh outward (fvcom-mesh-tools)

Read-only review of the git repository at the current directory. Do NOT
modify files. You may run read-only commands, python in memory, mocks and
fault injections (small synthetic inputs only; do not read the large data
under $DATA_DIR beyond listing it, and do not submit batch jobs). Answer in
English as Markdown.

## Goal
World-class correctness and robustness. Report every defect you can
substantiate, of any severity, in or outside the change, including
pre-existing ones.

## What was done
A tool that keeps a finished FVCOM base mesh exactly as it is and adds the
sea out to a new, designed open boundary (USER_GUIDE section 13). Read:

- `git show b0584f9 8e2739b 450ad44 d6d2a72 7044b6b 0f52d5b 69b50a4 d9e92fd b6d2ed8 765423c`
  (the extension tool and its documentation), and the current files:
  - `src/fvcom_mesh_tools/extend.py`, `extend_recipe.py`, `obc_design.py`,
    `dem/sources.py` (named bathymetry sources, priority stack, and the new
    `DATUM` registry / `non_tp_count` warning);
  - `notebooks/444_design_obc.py`, `445_extend_mesh.py`, `446_extend_generate.py`,
    `447_extend_merge.py`, `448_extend_smoke.py`, `453_redepth_extended.py`;
  - `recipes/extend/tokyo_bay_enshu.yaml`, `tokyo_bay_enshu_obc_design.yaml`;
  - `jobs/octopus/444_design_obc.sh`, `445_extend_mesh.sh`, `448_extend_smoke.sh`,
    `453_redepth_extended.sh`, `common.sh`;
  - tests: `tests/test_extend*.py`, `tests/test_obc_design*.py`,
    `tests/test_dem_sources.py` (whatever exists).
- Also in scope, just committed: the portability change --
  every job script and `common.sh` now take paths only from `$DATA_DIR` and
  `$WORK_DIR` (login profile), stop when they are unset, and derive the
  OCTOPUS FVCOM library directory as `FVCOM_LIBS` in `common.sh`; notebooks
  383/384/414 and `cli/refine_run.py` no longer fall back to `/octfs/...`.
  See commits 6d8b9a7 and 6c068d2 (`git log -5`).

Design intent:
- the base mesh's nodes, elements and depths are carried bit for bit
  (`verify_frozen_base`);
- the new part is generated with oceanmesh (subprocess; GPL code must never
  be imported by `fvcom_mesh_tools`), with fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-7 and their triage are in docs/extend-tools-review-20261001.md.
Round 7 (your previous answer; 12 findings) was fixed in 227536f; read
`git show 227536f`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 `M7001Points`: one `LinearNDInterpolator` over the whole N/M dataset,
  cached per (path, size, mtime); no query window.
- F3 `_nml_value` reads doubled delimiters in quoted values.
- F4 `_indices` (whole, finite, in range, checked before the cast) for
  elements, boundaries and sponge nodes; every public writer checks.
- F5 `export_fvcom_case`: cor/sponge validated first; files written to a
  temporary directory in OUTDIR and moved into place.
- F6 `relocate_case`: overlap refused (resolved paths), namelist rendered
  first, copy built in a sibling temporary directory.
- F7 407: work-area check before any staging.
- F8 loader reads the CSV once; `open_boundary_lonlat` and
  `open_boundary_sha256` come from the same bytes; 446 uses the former.
- F9 `provenance.changed_inventory`: datasets listed again at completion.
- F10 445/453: `collect(code=...)` again at completion; must match.
- F11 444: lock held until the marker is written; an existing `.prev`
  blocks publication.
- F12 453: success report written after the figure.
Verification after 227536f: full test suite 1144 passed (batch
job 122846); on real data 444 published the Enshu boundary, the full 445
build passed QA 23/23 with status ok and no inputs or code changed, its grd
file bit-identical to rounds 4-6, and 453 completed with status ok. The
whole-dataset M7001 triangulation (3.79 M points) builds in 55 s with 3.5
GB peak memory, once per process (job 122941); later queries are instant.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step gate (447,
453) REPORT warnings and do not fail the build. Not a defect.

## Please
1. Status of every previous finding: RESOLVED / PARTIAL / NOT RESOLVED /
   WITHDRAWN, with reasons.
2. Defects introduced by the fixes.
3. A fresh, unrestricted audit of the scope and everything it touches.

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
| F1 (F27) | major | as round 1 | yes, pre-existing | owner decision pending |
| F2 | minor | test (second move fails) | yes | fixed, 80a803c; test |
| F3 | minor | test (rename into place fails) | yes | fixed, 80a803c; test |
| F4 | minor | code read (float elements index later) | yes | fixed, 80a803c; test |
| F5 | minor | code read (iterator spent) | yes | fixed, 80a803c; test |
| F6 | minor | code read (int(1.9)) | yes | fixed, 80a803c; test |
| F7 | minor | code read (repo-wide dirty list compared) | yes | fixed (code_identity), 80a803c; test |
| F8 | minor | test (KEY = inside a quoted value) | yes | fixed, 80a803c; test |
| F9 | minor | test (file changed, cache served old) | yes | fixed, 80a803c; test |
| F10 | minor | test ((lon, lat) variable) | yes | fixed, 80a803c; test. Real grids are (lat, lon) |
| F11 | minor | code read (.4f) | yes | fixed (round-trip format), 80a803c; test |
| F12 | minor | test (-0.0 vs +0.0) | yes | fixed, 80a803c; test |
| F13 | minor | mock of the job's check | yes | fixed (cases from the manifest), 80a803c |
