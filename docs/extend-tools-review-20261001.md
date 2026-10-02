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

## Round 9

### Answer (verbatim)

Found **one previously reported major issue, one additional major issue, and 12 minor defects**. No files were modified.

Statuses assess the reported defects. Related defects found in this audit are listed separately; incomplete recovery guarantees are marked PARTIAL.

**Previous findings — round 1**

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked cells become NaN, enabling source fallback. |
| F2 | RESOLVED | Required QA failures produce nonzero exits. |
| F3 | RESOLVED | Limiting checks bounds, convergence and final rounded r-factor. |
| F4 | RESOLVED | Warning behavior implements the explicit owner policy. |
| F5 | RESOLVED | Chords, endpoint floors and published coordinates are checked. |
| F6 | RESOLVED | Missing boundary, lattice and generated-node coverage is rejected. |
| F7 | RESOLVED | Continuous departures and complete boundaries receive land checks. |
| F8 | RESOLVED | Rejected designs preserve existing CSVs. |
| F9 | RESOLVED | Native values round-trip; exported-base verification runs. |
| F10 | RESOLVED | Opposite seam incidence and global base overlap are checked. |
| F11 | RESOLVED | Single-edge land runs are retained. |
| F12 | RESOLVED | Unsupported constrained lines shorter than six nodes are rejected. |
| F13 | RESOLVED | CAO cache distinguishes roots and dimensions; file changes now invalidate it. |
| F14 | RESOLVED | Invalid fine samples retain valid coarse coverage. |
| F15 | RESOLVED | Degenerate M7001 point sets return uncovered samples. |
| F16 | RESOLVED | Depth files enter provenance; inventories are recollected. |
| F17 | RESOLVED | Relevant dirty sources and both rename paths are recorded. |
| F18 | RESOLVED | Recipe numeric controls, counts, seeds and bounds are validated. |
| F19 | RESOLVED | Resampling rejects nonfinite spacing and failure to progress. |
| F20 | RESOLVED | Extension outputs receive exclusive reservations. Finding 2 concerns older M2 workflows. |
| F21 | RESOLVED | Re-depth rejects overlapping destinations and reserves outputs. |
| F22 | RESOLVED | Smoke staging reserves a fresh root. |
| F23 | RESOLVED | Depth control precedes timestep selection. |
| F24 | RESOLVED | Smoke paths are resolved before use. |
| F25 | RESOLVED | Executable overrides are respected and resolved before directory changes. |
| F26 | RESOLVED | Job 383 uses the prepared names and checks directory changes. |
| F27 | NOT RESOLVED | GPL imports remain; owner decision pending. Finding 1, counted once. |
| F28 | RESOLVED | Recipe and re-depth case prefixes are restricted to filename components. |
| F29 | RESOLVED | Native reading validates unique, consecutive boundary traversal. |

**Previous findings — round 2**

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Re-depth and smoke reserve outputs exclusively. |
| F2 | RESOLVED | Unique temporaries and publication locking prevent competing swaps. |
| F3 | RESOLVED | Re-depth repeats frozen-base, overlap and QA checks; overrides are labelled. |
| F4 | RESOLVED | Chords, endpoint floors and publication precision are checked. |
| F5 | RESOLVED | Continuous departure checks catch thin land strips. |
| F6 | RESOLVED | Crossing tolerance covers numerical endpoint residue only. |
| F7 | RESOLVED | Checks use reconstructed published coordinates. |
| F8 | RESOLVED | Native coordinates and depths round-trip. |
| F9 | RESOLVED | Outer/base overlap is checked beyond the seam. |
| F10 | RESOLVED | Porcelain parsing includes both rename paths. |
| F11 | RESOLVED | Geographic bounds receive finite, range and ordering checks. |
| F12 | RESOLVED | Resampling requires finite forward progress. |
| F13 | RESOLVED | Traversal direction follows the supplied second node. |
| F14 | RESOLVED | Final-iteration convergence is accepted. |
| F15 | RESOLVED | Limited-edge depths and limiter controls are validated. |
| F16 | RESOLVED | Both sides and complete guide segments are checked. |
| F17 | RESOLVED | Bands are materialized and follow the approved deviation policy. |
| F18 | RESOLVED | Early provenance, consumed-input identity, completion checks and late reporting are handled. |
| F19 | RESOLVED | Uncovered generated nodes are rejected. |
| F20 | RESOLVED | Invalid radii, zero sides and reversals are rejected. |
| F21 | RESOLVED | Refusing smoke-root reuse preserves previous success markers. |
| F22 | NOT RESOLVED | Same pending F27. |

**Previous findings — round 3**

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Composed-field violations follow the owner’s warning policy. |
| F2 | RESOLVED | Publication permits numerical spacing tolerance only. |
| F3 | RESOLVED | Materializing bands prevents iterator exhaustion. |
| F4 | RESOLVED | Both-side probing and complete guide checks cover the reproduction. |
| F5 | RESOLVED | Marker failure retains locking; existing backups block publication. |
| F6 | RESOLVED | Early handling, consumed-input identity and late failure reporting are fixed. |
| F7 | RESOLVED | Zero-weight masked neighbours preserve valid samples. |
| F8 | RESOLVED | fort.14 numeric values round-trip; supplementary export is verified. |
| F9 | RESOLVED | Pairwise overlap checks avoid global-area dilution. |
| F10 | RESOLVED | Notebook 414 selects timestep after depth control. |
| F11 | RESOLVED | Design spacing controls must be finite and positive. |
| F12 | RESOLVED | Short designs receive explicit rejection. |
| F13 | RESOLVED | Filleting rejects reversals. |
| F14 | NOT RESOLVED | Same pending F27. |

**Previous findings — round 4**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Same pending F27. |
| F2 | RESOLVED | Half-window fits interpolate harmonic coefficients before deriving amplitude/phase. |
| F3 | RESOLVED | Bilinear interpolation preserves query shape. |
| F4 | RESOLVED | Overlap tolerance uses the smaller meeting element. |
| F5 | RESOLVED | Complete guide segments receive geometry checks. |
| F6 | RESOLVED | Boundary restoration, recovery locking and backup preservation cover the faults. |
| F7 | RESOLVED | State and failure handlers precede fallible logging. |
| F8 | RESOLVED | CSV coordinates and identity come from the same bytes. |
| F9 | RESOLVED | CAO retains valid centres beside missing cells, including final rows/columns. |
| F10 | RESOLVED | Boundary simplicity and repeated adjacent nodes are checked. |
| F11 | RESOLVED | Timestep overrides must be finite and positive. |
| F12 | RESOLVED | Executables are resolved and checked before directory changes. |
| F13 | RESOLVED | Printable ASCII and the 80-byte directory limit are enforced. |
| F14 | RESOLVED | History dimensions are checked; unreadable named grids fail validation. |

**Previous findings — round 5**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Same pending F27. |
| F2 | RESOLVED | Flip sites apply convexity guards and skip stale ownership rows. |
| F3 | RESOLVED | Recipe identity hashes the bytes actually parsed. |
| F4 | RESOLVED | Complete ladder segments receive geometry checks. |
| F5 | RESOLVED | Marker failure retains recovery protection; backups cannot be overwritten. |
| F6 | RESOLVED | An unreadable named grid fails run validation. |
| F7 | RESOLVED | Non-ASCII and overlength paths are refused. |
| F8 | RESOLVED | Handler registration precedes initial post-reservation logging. |
| F9 | RESOLVED | CAO preserves query shape. |
| F10 | RESOLVED | Nonfinite depths fail minimum-depth QA. |

**Previous findings — round 6**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Same pending F27. |
| F2 | RESOLVED | Buffered design and sea geometry are checked against their windows. |
| F3 | RESOLVED | Ordinary paths and doubled quoted delimiters are parsed correctly. |
| F4 | RESOLVED | Benchmark runs are validated before timing acceptance. |
| F5 | RESOLVED | QA failures and identity mismatches produce failure exits. |
| F6 | RESOLVED | Reference/work-area checks precede staging. |
| F7 | RESOLVED | Recovery-marker failure retains locking and protects backups. |
| F8 | RESOLVED | Completion checks rehash files and recollect inventories. |
| F9 | RESOLVED | Original fractional-index and writer-bypass cases are rejected. |
| F10 | RESOLVED | Writer/reader agreement covers doubled delimiters and literal backslashes. |
| F11 | RESOLVED | Relocation checks directory lengths and source/destination overlap. |

**Previous findings — round 7**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Same pending F27. |
| F2 | RESOLVED | M7001 uses whole-dataset triangulation; query composition no longer changes results. |
| F3 | RESOLVED | `_nml_value` decodes doubled delimiters. |
| F4 | RESOLVED | Writers reject the original fractional, nonfinite and out-of-range indices. |
| F5 | PARTIAL | Prevalidation and ordinary publication rollback work; interrupted recovery remains destructive. Finding 3. |
| F6 | RESOLVED | Source overlap is refused; rendering/copying precede destination replacement. |
| F7 | RESOLVED | Job 407 checks its work area before staging. |
| F8 | RESOLVED | CSV coordinates/hash share one read; generation consumes those coordinates. |
| F9 | RESOLVED | Completion inventories are recollected. |
| F10 | RESOLVED | Completion detects code changes while ignoring unrelated repository dirt. |
| F11 | RESOLVED | Lock survives marker failure; existing `.prev` blocks publication. |
| F12 | RESOLVED | Figure generation precedes success reporting. |

**Previous findings — round 8**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Same pending F27. |
| F2 | PARTIAL | Ordinary move failures restore previous files; interruption during restoration deletes recovery copies. Finding 3. |
| F3 | PARTIAL | Ordinary destination-move failures restore the old case; interrupted restoration can erase it. Finding 3. |
| F4 | RESOLVED | Export requires integer element dtype. |
| F5 | RESOLVED | Sponge iterators are normalized once and the resulting rows are written. |
| F6 | RESOLVED | OBC types must be whole, non-boolean numbers. |
| F7 | RESOLVED | Comparison excludes repository-wide dirt but retains relevant code identity. |
| F8 | RESOLVED | Quoted assignments are ignored; conflicting duplicate assignments fail validation. |
| F9 | RESOLVED | CAO cache keys include resolved path, size, mtime_ns and shape. |
| F10 | RESOLVED | Grid dimension order is checked; `(lon, lat)` data are transposed. |
| F11 | RESOLVED | Sponge radius and damping use round-trip formatting. |
| F12 | RESOLVED | Coordinates/depths are compared by dtype and bytes. Finding 6 concerns merge compatibility. |
| F13 | RESOLVED | Job 413 obtains case names from manifest `runs`. |

The generation-nondeterminism hypothesis remains **WITHDRAWN** based on the supplied repeated-run evidence. Owner-authorized band-floor and new-element timestep warnings are not defects.

**Findings**

1. **Major — Known F27 remains open: package code imports oceanmesh.**  
   [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

   **Evidence:** Direct imports remain here and at line 165, in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672`, and `autofinish/directives.py:20`. They conflict with the explicit prohibition in `CLAUDE.md:32`. The extension subprocess does not remove these other imports.

   **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary and reconcile `THIRD_PARTY_NOTICES.md`. **Previously reported; owner decision pending; counted once.**

2. **Major — Older M2 workflows allow concurrent jobs to overwrite the same run.**  
   [383_m2.sh:19](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/383_m2.sh:19), [410_m2_chain.sh:17](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/410_m2_chain.sh:17).

   **Evidence:** Job 383 always uses `scratch/m2_383`. Two submissions can both pass the history-file checks before either produces history, then prepare and integrate into the same directories. Their logs are also opened with truncation. Chains 410 and 415 use one-second timestamps and `mkdir -p`, allowing identical roots. Job 412 has no per-case execution lock and removes existing history at line 48.

   This permits one job to overwrite another’s inputs/output or validate output produced by the other job. The extension output reservation does not protect these workflows.

   **Fix:** Reserve run roots atomically before staging/submission, use unique roots, and acquire a per-case execution lock before invalidating markers or deleting history.

3. **Minor — Interrupted rollback deletes the recovery copies. Introduced by `80a803c`.**  
   [fvcom_native.py:476](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:476), [fvcom_namelist.py:97](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_namelist.py:97).

   **Reproduction:** In-memory filesystem mocks injected an `OSError` during publication, followed by `KeyboardInterrupt` during restoration. Export left the new grid alongside old depths/OBC and deleted `.previous`, losing the original grid. Relocation deleted its temporary directory containing the entire previous destination.

   Both restoration handlers catch only `OSError`; their outer `finally` blocks still clean up after another exception. Ordinary `OSError` rollback restored the previous contents successfully.

   **Fix:** Preserve recovery storage by default once replacement begins. Enable cleanup only after publication or restoration is confirmed complete; retain and identify backups on every exceptional restoration exit, including `BaseException`.

4. **Minor — Merge silently truncates malformed connectivity and boundary indices.**  
   [extend.py:140](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:140).

   **Reproduction:** Merging adjacent squares with outer connectivity `[[0., 1.9, 2.], [0., 2., 3.]]` succeeds. The resulting outer elements are `[[1,4,5],[1,5,2]]`, and `verify_frozen_base` passes. The invalid `1.9` has become node 1 before subsequent validation can detect it.

   The interface arrays and `outer_open` are similarly cast without validation.

   **Fix:** Validate connectivity dtype, shape and range before conversion; validate boundary/interface indices as finite whole numbers before casting. Reuse the checked index normalization used by native writers.

5. **Minor — QA reports success for fractional boundary IDs.**  
   [qa.py:536](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:536).

   **Reproduction:** On the pristine fixture from `tests/test_qa.py`, set `open_boundaries=[np.array([4.9,8.9])]`. `run_qa(..., channel_check=False)` reports **23/23 gates passed**, with 16 nodes, 18 elements and implied dt **142.784312 s**. The input still contains fractional IDs; QA evaluated their truncated substitutes.

   The integrity gate at lines 571–574 also casts before checking.

   **Fix:** Validate index shapes, integrality, finiteness and range before any cast. Return a failed integrity gate for malformed connectivity or boundary indices.

6. **Minor — Merge violates the frozen-base contract for float32 bases. Exposed by the stronger verifier.**  
   [extend.py:155](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:155).

   **Reproduction:** A four-node base with float32 coordinates and depths is accepted by `merge_outer`. Appending the outer mesh changes both arrays to float64; `verify_frozen_base` then raises `base node coordinates changed`.

   The API promises unchanged base nodes/depths and declares no float64-only precondition.

   **Fix:** Preserve the base coordinate/depth dtype when allocating appended arrays, with an explicit policy for representing new values. Alternatively, document and enforce a float64 precondition at merge entry.

7. **Minor — Namelist replacement silently deletes neighbouring assignments.**  
   [fvcom_namelist.py:50](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_namelist.py:50).

   **Reproduction:** Given the valid namelist line:

   ```fortran
   INPUT_DIR='old/', NC_OUT_INTERVAL='seconds = 1800.',
   ```

   `set_value(text, "INPUT_DIR", "'/new/input/'")` produces only:

   ```fortran
   INPUT_DIR= '/new/input/',
   ```

   `relocate_case` uses this helper and can therefore return success after removing other model settings.

   **Fix:** Replace only the selected assignment’s value using quote/comment-aware parsing. Until supported, reject lines containing neighbouring assignments before writing the relocated case.

8. **Minor — Sponge validation rejects a numeric NumPy row array.**  
   [fvcom_native.py:376](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:376).

   **Reproduction:** `_check_sponge(mesh, np.array([[1,1000,.001]]))` raises “The truth value of an array with more than one element is ambiguous.” The rows themselves satisfy the subsequent numeric validation.

   **Fix:** Use an explicit `sponge is None` check, then materialize the iterable with `list(sponge)`.

9. **Minor — Corner trimming can remove mutually supporting triangles and then crash.**  
   [extend.py:374](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:374).

   **Reproduction:** `trim_lone_corners([[0,1,2],[0,2,3]], [True,True])` marks both triangles for removal because each sees the shared opposite edge. The next iteration calls `t.max()` on an empty array and raises a zero-size reduction error.

   The simultaneous removals also invalidate the comment’s assertion that removal leaves no new lone node.

   **Fix:** Assess removals against the surviving topology, preventing mutually dependent deletions. Handle empty input/result explicitly before reductions.

10. **Minor — The limiter skips requested free-node bounds when no edge touches a free node.**  
    [extend.py:305](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:305).

    **Reproduction:**

    ```python
    rfactor_smooth_free(
        [5,5,1], np.array([0]), np.array([1]),
        [False,False,True], rmax=.2, hmin=3, hmax=4,
    )
    ```

    returns depths `[5,5,1]`, zero iterations and r=0. The free node remains below `hmin`, contrary to the documented bounds guarantee.

    **Fix:** Apply and validate free-node depth bounds before the no-live-edge return. Only the edge iteration should be skipped.

11. **Minor — The fort.14 reader rejects a valid single-triangle mesh.**  
    [fort14.py:120](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:120).

    **Reproduction:** A synthetic standard file with `NE=1`, `NP=3`, three node rows, one triangle and complete boundary records raises:

    ```text
    element block shape (5,) does not match expected (1, 5)
    ```

    `np.loadtxt` squeezes the single element row.

    **Fix:** Read tabular blocks with `ndmin=2`; add a single-triangle round-trip case.

12. **Minor — The fort.14 reader ignores record IDs and the element type field.**  
    [fort14.py:117](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:117), [fort14.py:125](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:125).

    **Reproduction:** An in-memory two-element/four-node file containing duplicate node ID 1, duplicate element ID 1, and an element row `1 4 1 2 3` is accepted as an ordinary triangular mesh. Those metadata columns are discarded.

    Thus malformed identifiers can silently change the interpretation of connectivity/boundary references, and a non-triangle declaration is treated as a triangle.

    **Fix:** Validate unique node/element IDs and triangle type 3. Either enforce canonical ordering or map explicit node IDs to array positions before translating references.

13. **Minor — M2 chain wrappers do not protect comma-containing configured paths.**  
    [410_m2_chain.sh:25](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/410_m2_chain.sh:25), [415_refine_m2_chain.sh:39](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/415_refine_m2_chain.sh:39).

    **Evidence:** With `WORK_DIR=/scratch/user,project`, the constructed argument starts:

    ```text
    FMESH_RUN_ROOT=/scratch/user,project/scratch/m2_...
    ```

    `qsub -v` uses commas as variable separators, so this does not transmit the intended root. Shell quoting does not escape that delimiter. Candidate/base/refined paths have the same problem. No jobs were submitted.

    **Fix:** Validate every `-v` value before the first submission, using the comma guard already present in `refine_workflow.sh:51`.

14. **Minor — Relative `FVCOM_LIBS` overrides change meaning after entering the case directory. Introduced by the portability change.**  
    [common.sh:38](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/common.sh:38).

    **Reproduction:** Executing the assignment in memory with `FVCOM_LIBS=vendor` leaves `vendor` unchanged. Job 448 builds `LD_LIBRARY_PATH=vendor/lib:vendor/lib64:...`, then changes into the case directory before launching MPI. Those library paths now refer to directories under the case, rather than under the submission directory.

    Executable overrides already receive the necessary absolute-path treatment; the new library override does not.

    **Fix:** Resolve and validate `FVCOM_LIBS` to an absolute directory in `common.sh`, before any caller changes directory.

**Verification**

- In the `fvcom-mesh-tools` environment, `pytest --noconftest -s -q -p no:cacheprovider tests/test_extend.py tests/test_obc_design.py`, with bytecode disabled: **31 passed**.
- Direct invocation of 49 tests requiring no fixtures, including those geometry tests: **49 passed**.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` on all 56 OCTOPUS shell scripts: **passed**.
- In-memory fault injections and synthetic reproductions produced the failures above. Additional mocked checks confirmed dimension-order handling, CAO file-change invalidation and quoted/duplicate namelist parsing.
- The full suite and production mesh/FVCOM runs were not repeated. No shared data were read and no batch jobs were submitted. Git status remained unchanged.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 12 minor, 0 nit)

### Prompt

```markdown
# Review request, round 9: extending a base mesh outward (fvcom-mesh-tools)

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
Rounds 1-8 and their triage are in docs/extend-tools-review-20261001.md.
Round 8 (your previous answer; 13 findings) was fixed in 80a803c; read
`git show 80a803c`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 `export_fvcom_case`: replaced files are copied to `<stage>/.previous`
  first and moved back if a later move fails; if that fails, the stage is
  kept and named in the error.
- F3 `relocate_case`: the old destination is renamed aside and renamed back
  on failure; kept and named if that fails.
- F4 elements require an integer dtype.
- F5 sponge rows normalised once and passed to `write_spg`.
- F6 `write_obc`: types must be whole, non-bool numbers.
- F7 `provenance.code_identity` drops the repository-wide `git.dirty` list
  from the comparison (dirt under the code's path still changes
  `commit_identifies_code`/`dirty_under_path`/`source_sha256`).
- F8 `_nml_value` searches a copy with quoted text blanked, and refuses a
  key set twice with different values (check_run reports it as a reason).
- F9 CAO cache key: resolved file path, size, mtime_ns, shape.
- F10 `Grid.depth`: (lat, lon) as is, (lon, lat) transposed, else refused.
- F11 `write_spg` uses the round-trip `_num` format.
- F12 `verify_frozen_base`: dtype + bytes for coordinates and depths.
- F13 413 reads the case list from `manifest.json` `runs`.
Verification after 80a803c: full test suite 1152 passed (batch
job 122957); on real data 444 published the Enshu boundary, the full 445
build passed QA 23/23 with status ok and nothing changed during the build
(the review record was edited meanwhile, as F7 intends), grd and dep files
bit-identical to rounds 4-7, and 453 completed with status ok.
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
| F2 | major | code read (fixed root, 1-s stamps, no lock) | yes | fixed (unique roots, per-case lock), f05220a; lock mocked |
| F3 | minor | test (KeyboardInterrupt during restore) | yes | fixed, f05220a; test |
| F4 | minor | test (1.9 in outer elements) | yes | fixed, f05220a; test |
| F5 | minor | test (fractional OBC ids) | yes | fixed, f05220a; test |
| F6 | minor | test (float32 base) | yes | fixed (float64 required), f05220a; test |
| F7 | minor | test (two assignments on a line) | yes | fixed, f05220a; test; 383 namelist output checked (whitespace only) |
| F8 | minor | code read (array truthiness) | yes | fixed, f05220a |
| F9 | minor | test (two triangles of a square) | yes | fixed, f05220a; test |
| F10 | minor | test (no live edge) | yes | fixed, f05220a; test |
| F11 | minor | test (NE=1) | yes | fixed, f05220a; test |
| F12 | minor | test (duplicate ids, type 4) | yes | fixed, f05220a; tests |
| F13 | minor | code read | yes | fixed, f05220a |
| F14 | minor | evaluation of the assignment | yes | fixed, f05220a |

## Round 10

### Answer (verbatim)

Found **4 major, 8 minor and 2 nit findings**, including the previously reported GPL issue. No files were modified.

Statuses below distinguish fixes to the original reproductions from remaining gaps. Repeated reports of the GPL issue count as one finding.

**Previous findings — round 1**

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked bathymetry becomes NaN and permits fallback. |
| F2 | RESOLVED | Required QA failures cause nonzero exits. |
| F3 | RESOLVED | The original cap/r-factor conflict and nonconvergence are rejected; final r-factor is recomputed. Finding 12 is a separate bounds issue. |
| F4 | RESOLVED | Warning behavior implements the explicit owner decision. |
| F5 | RESOLVED | Actual chords, endpoint floors and published coordinates are checked. |
| F6 | RESOLVED | Missing sizing coverage at boundaries, lattice samples and generated nodes is rejected. |
| F7 | RESOLVED | Continuous departure and full-boundary land checks cover the reproduction. |
| F8 | RESOLVED | Validation precedes publication; rejected designs preserve existing CSVs. |
| F9 | RESOLVED | Native values round-trip, and exported-base verification runs. |
| F10 | RESOLVED | Seam-side incidence and outer/base overlap are checked. |
| F11 | RESOLVED | Single-edge land runs are retained. |
| F12 | RESOLVED | Unsupported constrained lines shorter than six nodes are rejected explicitly. |
| F13 | RESOLVED | CAO cache keys distinguish roots, dimensions and file changes. |
| F14 | RESOLVED | Invalid fine samples retain valid coarse coverage. |
| F15 | RESOLVED | Degenerate M7001 point sets return uncovered samples. |
| F16 | RESOLVED | Depth files enter provenance; completion inventories are recollected. |
| F17 | RESOLVED | Relevant dirty sources and both rename paths are recorded. |
| F18 | RESOLVED | Numeric controls, counts, seeds and bounds are validated. |
| F19 | RESOLVED | Resampling rejects nonfinite spacing and failure to progress. |
| F20 | RESOLVED | Extension output directories are reserved exclusively. |
| F21 | RESOLVED | Re-depth rejects overlapping destinations and reserves outputs. |
| F22 | RESOLVED | Extension smoke staging reserves a fresh root. |
| F23 | RESOLVED | Depth control precedes timestep selection. |
| F24 | RESOLVED | Smoke paths are resolved before use. |
| F25 | RESOLVED | Executable overrides are respected and resolved before directory changes. |
| F26 | RESOLVED | Job 383 uses prepared case names and guards directory changes. |
| F27 | NOT RESOLVED | Direct GPL imports remain; owner decision pending. Finding 1. |
| F28 | RESOLVED | Recipe and re-depth case names are restricted to filename components. |
| F29 | RESOLVED | Native reading validates unique, consecutive boundary traversal. |

**Previous findings — round 2**

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Re-depth and extension smoke reserve outputs exclusively. |
| F2 | RESOLVED | Unique temporaries and publication locking prevent competing boundary swaps. |
| F3 | RESOLVED | Re-depth repeats frozen-base, overlap and QA checks; overrides are identified. |
| F4 | RESOLVED | Chords, endpoint floors and publication precision are checked. |
| F5 | RESOLVED | Continuous departure checks catch thin land strips. |
| F6 | RESOLVED | Crossing tolerance accommodates numerical endpoint residue. |
| F7 | RESOLVED | Checks use reconstructed published coordinates. |
| F8 | RESOLVED | Native coordinates and depths round-trip. |
| F9 | RESOLVED | Outer/base overlap is checked beyond the seam. |
| F10 | RESOLVED | Porcelain parsing includes both rename paths. |
| F11 | RESOLVED | Geographic bounds receive finite, range and ordering checks. |
| F12 | RESOLVED | Resampling requires finite forward progress. |
| F13 | RESOLVED | Traversal follows the supplied second node. |
| F14 | RESOLVED | Final-iteration convergence is accepted. |
| F15 | RESOLVED | Limited-edge depths and limiter controls are validated. |
| F16 | RESOLVED | Both sides and complete guide segments are checked. |
| F17 | RESOLVED | Bands are materialized and follow the approved warning policy. |
| F18 | RESOLVED | Early provenance, consumed-input identity, completion checks and late reporting are handled. |
| F19 | RESOLVED | Uncovered generated nodes are rejected. |
| F20 | RESOLVED | Invalid radii, zero sides and reversals are rejected. |
| F21 | RESOLVED | Refusing extension smoke-root reuse preserves previous markers. |
| F22 | NOT RESOLVED | Alias of pending F27; finding 1. |

**Previous findings — round 3**

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Composed-field violations follow the explicit owner warning policy. |
| F2 | RESOLVED | Publication allows only numerical spacing tolerance. |
| F3 | RESOLVED | Materializing bands prevents iterator exhaustion. |
| F4 | RESOLVED | Both-side probing and complete guide checks cover the reproduction. |
| F5 | PARTIAL | Marker-failure protection works, but interrupted recovery can still delete backups; finding 5. |
| F6 | RESOLVED | Early handling, consumed-input identity and late failure reporting are fixed. |
| F7 | RESOLVED | Zero-weight masked neighbours preserve valid samples. |
| F8 | RESOLVED | fort.14 values round-trip; supplementary export is verified. |
| F9 | RESOLVED | Pairwise overlap checks avoid global-area dilution. |
| F10 | RESOLVED | Notebook 414 selects its timestep after depth control. |
| F11 | RESOLVED | Design spacing controls must be finite and positive. |
| F12 | RESOLVED | Unsupported short designs receive explicit rejection. |
| F13 | RESOLVED | Filleting rejects reversals. |
| F14 | NOT RESOLVED | Alias of pending F27; finding 1. |

**Previous findings — round 4**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Alias of pending F27; finding 1. |
| F2 | RESOLVED | Half-window fitting interpolates harmonic coefficients before deriving amplitude/phase. |
| F3 | RESOLVED | Bilinear interpolation preserves query shape. |
| F4 | RESOLVED | Overlap tolerance uses the smaller meeting element. |
| F5 | RESOLVED | Complete guide segments receive geometry checks. |
| F6 | PARTIAL | Ordinary restoration works; interrupted restoration remains destructive in 444. Finding 5. |
| F7 | RESOLVED | State and failure handlers precede fallible logging. |
| F8 | RESOLVED | CSV coordinates and identity come from the same bytes. |
| F9 | RESOLVED | CAO retains valid centres beside missing cells, including final rows/columns. |
| F10 | RESOLVED | Boundary simplicity and repeated adjacent nodes are checked. |
| F11 | RESOLVED | Timestep overrides must be finite and positive. |
| F12 | RESOLVED | Executables are resolved and checked before directory changes. |
| F13 | RESOLVED | Printable ASCII and the 80-byte directory limit are enforced. |
| F14 | RESOLVED | Required history dimensions are checked; unreadable named grids fail validation. Finding 4 concerns unchecked fields. |

**Previous findings — round 5**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Alias of pending F27; finding 1. |
| F2 | RESOLVED | Flip sites apply convexity guards and skip stale ownership rows. |
| F3 | RESOLVED | Recipe identity hashes the bytes actually parsed. |
| F4 | RESOLVED | Complete ladder segments receive geometry checks. |
| F5 | PARTIAL | Marker failures preserve recovery protection; other exceptional recovery exits still delete backups. Finding 5. |
| F6 | RESOLVED | An unreadable named grid fails run validation. |
| F7 | RESOLVED | Non-ASCII and overlength run paths are refused. |
| F8 | RESOLVED | Handlers precede initial post-reservation logging. |
| F9 | RESOLVED | CAO preserves query shape. |
| F10 | RESOLVED | Nonfinite depths fail minimum-depth QA. |

**Previous findings — round 6**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Alias of pending F27; finding 1. |
| F2 | RESOLVED | Buffered design and sea geometry are checked against their windows. |
| F3 | RESOLVED | Ordinary paths and doubled quoted delimiters parse correctly. |
| F4 | RESOLVED | Benchmark runs undergo validation before timing acceptance. Finding 4 concerns incomplete field validation. |
| F5 | RESOLVED | QA failures and identity mismatches cause failure exits. |
| F6 | RESOLVED | Reference/work-area checks precede staging. |
| F7 | PARTIAL | Recovery-marker failure retains locking, but interrupted restoration bypasses that protection. Finding 5. |
| F8 | RESOLVED | Completion checks rehash files and recollect inventories. |
| F9 | RESOLVED | The original native fractional-index and writer-bypass cases are rejected. |
| F10 | RESOLVED | Writer/reader agreement covers doubled delimiters and literal backslashes. |
| F11 | RESOLVED | Relocation checks directory lengths and source/destination overlap. |

**Previous findings — round 7**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Alias of pending F27; finding 1. |
| F2 | RESOLVED | Whole-dataset M7001 triangulation removes query-composition dependence. |
| F3 | RESOLVED | `_nml_value` decodes doubled delimiters. |
| F4 | RESOLVED | Native writers reject the original fractional, nonfinite and out-of-range indices. |
| F5 | PARTIAL | Prevalidation and ordinary rollback work; interruption after replacement can escape tracking. Finding 6. |
| F6 | RESOLVED | Source overlap is refused; rendering and copying precede replacement. Remaining recovery defects are covered by finding 6. |
| F7 | RESOLVED | Job 407 checks its work area before staging. |
| F8 | RESOLVED | CSV coordinates/hash share one read; generation consumes those coordinates. |
| F9 | RESOLVED | Completion inventories are recollected. |
| F10 | RESOLVED | Completion detects relevant code changes while ignoring unrelated dirt. |
| F11 | PARTIAL | Marker-failure locking and existing-backup refusal work; interrupted restoration still destroys recovery copies. Finding 5. |
| F12 | RESOLVED | Figure generation precedes success reporting. |

**Previous findings — round 8**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Alias of pending F27; finding 1. |
| F2 | PARTIAL | Ordinary rollback and the original interrupted-restoration case are fixed; post-replacement interruption still loses backups. Finding 6. |
| F3 | PARTIAL | Original restoration interruption is fixed; interruption after moving the old case aside remains destructive. Finding 6. |
| F4 | RESOLVED | Native export requires integer element dtype. |
| F5 | RESOLVED | Sponge iterators are normalized once and those rows are written. |
| F6 | RESOLVED | OBC types must be whole, non-boolean numbers. |
| F7 | RESOLVED | Comparison excludes repository-wide dirt while retaining relevant code identity. |
| F8 | RESOLVED | Quoted assignments are ignored; conflicting duplicates fail validation. |
| F9 | RESOLVED | CAO cache keys include resolved path, size, mtime and shape. |
| F10 | RESOLVED | Grid dimension order is checked; `(lon, lat)` data are transposed. |
| F11 | RESOLVED | Sponge values use round-trip formatting. |
| F12 | RESOLVED | Coordinates/depths are compared by dtype and bytes. Finding 10 concerns connectivity. |
| F13 | RESOLVED | Job 413 reads case names from the manifest. |

**Previous findings — round 9**

| ID | Status | Reason |
|---|---|---|
| F1 | NOT RESOLVED | Alias of pending F27; finding 1. |
| F2 | PARTIAL | 383/410/415 scratch-root reservation and the 412 execution lock are fixed. Analysis destinations still collide; finding 2. Other workflow roots remain unsafe; finding 3. |
| F3 | PARTIAL | Interruption during restoration now preserves backups. State/tracking updates still occur after destructive renames; finding 6. |
| F4 | RESOLVED | Merge validates outer/interface/boundary indices before conversion. |
| F5 | PARTIAL | Fractional IDs now fail the integrity gate, but OBC conversion still precedes validation and shapes are unchecked; finding 8. |
| F6 | RESOLVED | Merge documents and enforces float64 base coordinates/depths. |
| F7 | RESOLVED | Only the selected value is replaced; neighbouring assignments survive. The replacement introduced the comment regression in finding 7. |
| F8 | RESOLVED | Explicit `is None` permits NumPy sponge rows. |
| F9 | RESOLVED | Sequential edge counts prevent mutually dependent deletions; empty meshes are guarded. Finding 11 concerns reporting at the iteration limit. |
| F10 | RESOLVED | Free-node bounds apply before returning when there are no live edges. |
| F11 | RESOLVED | `ndmin=2` supports a single-triangle file. |
| F12 | RESOLVED | Canonical node/element IDs, triangle type and element references are checked. Boundary references remain unchecked; finding 9. |
| F13 | RESOLVED | Both chain wrappers reject commas in transmitted values. |
| F14 | RESOLVED | `FVCOM_LIBS` becomes absolute before directory changes; existence need not be required for unrelated jobs. |

The generation-nondeterminism hypothesis remains **WITHDRAWN**, consistent with the supplied repeated-run evidence. Owner-authorized band-floor and new-element timestep warnings are not defects.

**Findings**

1. **Major — Known F27 remains open: package code directly imports oceanmesh.**  
   [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

   **Evidence:** Direct imports also remain at line 165, in `mesh_engine/multiscale.py:52–53,139,145`, `mesh_clean.py:1672`, and `autofinish/directives.py:20`. They violate the explicit repository policy in `CLAUDE.md:32–34`. The extension subprocess does not remove these imports.

   **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary and reconcile the third-party notices. **Previously reported; owner decision pending; counted once.**

2. **Major — Unique M2 scratch roots still publish into shared analysis destinations.**  
   [383_m2.sh:77](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/383_m2.sh:77), [384_m2_analysis.py:464](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/384_m2_analysis.py:464), [415_refine_m2_chain.sh:62](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/415_refine_m2_chain.sh:62).

   **Reproduction:** Executing notebook 384’s argument handling with analysis mocked shows that roots `/scratch/m2_383.123` and `/scratch/m2_383.124` both select `/repo/outputs/m2_383` and its `m2_383.png`. Job 383 supplies no `--output`. JSON, CSV and figure publication can overwrite or mix two experiments.

   Chain 415 now reserves unique scratch roots, but its analysis directory remains `outputs/m2r_$STAMP`, with one-second timestamp resolution.

   **Fix:** Derive and exclusively reserve each analysis destination from the unique run identifier; pass it explicitly from the submitting job.

3. **Major — Refinement workflows can still assign different experiments the same scratch root.**  
   [refine_workflow.sh:44](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/refine_workflow.sh:44), [421_finish_and_run.sh:21](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/421_finish_and_run.sh:21), [423_m2_smoke.sh:30](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/423_m2_smoke.sh:30).

   **Reproduction:** Evaluating the root-selection block with the supplied work-directory length, timestamp `20261002_120000`, and names `kimitsu_port_hires` and `kisarazu_hires` produces the same root for both:
   ```
   /octfs/work/G16445/v61021/scratch/m2_20261002_120000
   ```
   The length fallback discards the recipe name, and the root is not reserved.

   Jobs 421 and 423 invalidate shared markers and stage/replace cases without a lock. The new 412 lock therefore does not protect preparation or smoke execution. Concurrent workflows can overwrite inputs/history and accept another experiment’s output.

   **Fix:** Reserve a unique short root atomically. Protect every stage that mutates an existing run with locks acquired before marker invalidation, staging or history removal; coordinate those locks with 412.

4. **Major — Smoke validation accepts nonfinite three-dimensional model fields.**  
   [check_run.py:261](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:261), [448_extend_smoke.sh:56](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/448_extend_smoke.sh:56).

   **Reproduction:** An in-memory history with four correctly spaced records, matching grid dimensions, finite `zeta/ua/va`, all-NaN `u`, and all-infinite `w` returned:
   ```python
   {"n_records": 4, "ok": True, "reasons": []}
   ```
   With a successful exit and clean log, job 448 consequently writes `SMOKE_OK`. A failed 3D solution can coexist with apparently healthy surface/barotropic output.

   **Fix:** Check finiteness and record dimensions of emitted numeric model fields, including 3D velocities and enabled scalars/turbulence fields. Process records incrementally to bound memory, and test corrupted 3D output with healthy barotropic fields.

5. **Minor — Interrupted boundary rollback deletes the original backups and releases the lock.**  
   [444_design_obc.py:277](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:277).

   **Reproduction:** Executing the actual publication block against an in-memory filesystem: publish the new JSON, fail PNG publication with `OSError`, then interrupt JSON restoration with `KeyboardInterrupt`. The final state contains old CSV, new JSON and old PNG; all `.prev` backups and the lock are deleted.

   Restoration catches only `OSError`, while `finally` unconditionally removes remaining backups unless the marker-specific path changed the state.

   **Fix:** Preserve backups and locking by default once publication starts. Allow cleanup only after publication or restoration is confirmed complete, and retain recovery information on every exceptional restoration exit.

6. **Minor — Export and relocation record destructive renames too late.**  
   [fvcom_native.py:470](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:470), [fvcom_namelist.py:122](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_namelist.py:122).

   **Reproduction:** Filesystem mocks perform a rename and then raise `KeyboardInterrupt`, modelling interruption after the filesystem operation but before Python bookkeeping:

   - Export replaces the grid before adding it to `moved`. Rollback sees no moved grid, declares restoration complete, deletes its backup, and leaves new grid with old depth/OBC files.
   - Relocation moves the old destination into `tmp/previous` before changing `state` from `"staging"`. Cleanup deletes the temporary directory and the entire old destination.

   The round-9 fix protects interruptions during restoration but leaves these publication windows. The relocation state assignment is newly introduced by `f05220a`.

   **Fix:** Record recovery intent before destructive operations; reconcile filesystem state during recovery. Preserve backups whenever the outcome is uncertain. Add post-effect interruption tests alongside ordinary rename-failure tests.

7. **Minor — The new namelist replacement parser treats comments as assignments. Introduced by `f05220a`.**  
   [fvcom_namelist.py:47](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_namelist.py:47).

   **Reproduction:** Calling `set_value(..., "INPUT_DIR", "'/new/'")` gives:

   | Input | Result |
   |---|---|
   | Commented example followed by active `INPUT_DIR` | Rejects valid text: “expected once, found 2” |
   | Comment containing an unmatched apostrophe, followed by active `INPUT_DIR` | “expected once, found 0” |
   | Only `! INPUT_DIR='old/',` | Returns success after modifying only the comment |

   `_mask_quoted` neither masks Fortran `!` comments nor prevents quotes inside comments from affecting subsequent lines.

   **Fix:** Scan strings and comments together, masking comments without interpreting their quotes. Match only active assignments and preserve surrounding text.

8. **Minor — QA converts boundary IDs before validating them and does not validate index-array shapes.**  
   [qa.py:536](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:536), [qa.py:574](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:574).

   **Reproduction:** A three-node/one-element mesh with OBC `["bad", "1"]` raises `ValueError` during the initial int64 conversion, before producing the failed integrity report. On the pristine 16-node/18-element fixture, OBC `[[4,8]]` passes value checks and later raises `IndexError`; four-column elements likewise reach geometry with an invalid shape.

   **Fix:** Validate connectivity as `(NE,3)` and boundary chains as one-dimensional arrays, including dtype, finiteness, integrality and range, before conversion or geometry. Return failed integrity checks for malformed inputs.

9. **Minor — fort.14 boundary indices remain unchecked and can be silently truncated.**  
   [fort14.py:93](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:93), [fort14.py:210](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:210).

   **Reproduction:** A synthetic single-triangle file with external boundary IDs `0` and `4` for `NP=3` is accepted, returning `[-1,3]`. Writing OBC indices `[1.9,2.9]` succeeds and emits external IDs `2,3`, silently changing the boundary.

   **Fix:** Validate writer boundaries as finite, whole, in-range indices before opening the destination. Validate reader boundary IDs against `1..NP`; also reject invalid counts and inconsistent totals.

10. **Minor — Frozen-base verification permits connectivity dtype and byte changes.**  
    [extend.py:176](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:176), [extend.py:200](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:200).

    **Reproduction:** Merge a four-node/two-element base with int32 connectivity into an adjacent square. The result has six nodes/four elements, and base connectivity changes from int32/24 bytes to int64/48 bytes. `verify_frozen_base` succeeds and reports one valid interface edge.

    Connectivity uses numerical equality, whereas the documented contract requires elements bit for bit.

    **Fix:** Compare connectivity with the existing dtype-and-bytes helper. Preserve the base integer dtype with overflow checks, or document and enforce an int64 precondition at merge entry.

11. **Minor — Corner trimming reports stale surviving lone nodes at the iteration limit.**  
    [extend.py:409](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:409).

    **Reproduction:**
    ```python
    trim_lone_corners([[0,1,2], [0,2,3]], [True,True], max_rounds=1)
    ```
    returns `[[0,2,3]]`, one dropped element and `lone_nodes_left=[3]`. Actual surviving lone nodes are `[0,2,3]`. The report was calculated before the final removal.

    **Fix:** Recompute surviving lone nodes from the returned mesh and report whether the iteration limit prevented further trimming.

12. **Minor — Rounding after depth limiting can violate recipe depth bounds.**  
    [447_extend_merge.py:157](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:157), [453_redepth_extended.py:128](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:128).

    **Reproduction:** With one fixed and one free node, `rmax=.2`:

    - `hmin=3.0000004` yields free depth `3.0000004`; rounding changes it to `3.0`, below the floor.
    - `hmax=3.0000006` yields free depth `3.0000006`; rounding changes it to `3.000001`, above the cap.

    These finite recipe bounds are accepted. Both results pass the subsequent r-factor check; final QA does not enforce these recipe-specific bounds.

    **Fix:** Remove unnecessary six-decimal rounding now that serialization round-trips, or quantize inward and revalidate both bounds and r-factor on final depths.

13. **Nit — M2 instructions describe obsolete roots, case names and figure paths.**  
    [m2_383.md:21](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/docs/m2_383.md:21).

    **Evidence:** Documentation still describes fixed `scratch/m2_383`, `B_Adepth`, mandatory archival before resubmission, and `outputs/figures/383_m2_comparison.png`. Current job 383 uses a reserved per-job root, stages `B_m7001`, and notebook 384 defaults its figure to `outputs/m2_383/m2_383.png`.

    **Fix:** Update the instructions to the current root reservation, case names and analysis destinations, incorporating finding 2’s fix.

14. **Nit — fort.14 writer documentation describes formatting it no longer uses.**  
    [fort14.py:171](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:171).

    **Evidence:** The docstring promises coordinates written with 15 decimal digits. Line 199 actually uses `float(...).__repr__`, producing shortest round-trip representations.

    **Fix:** Describe the current round-trip formatting rather than a fixed decimal count.

**Verification**

- In the project conda environment, with bytecode and pytest cache disabled:
  ```
  python -B -m pytest --noconftest -s -q -p no:cacheprovider \
      tests/test_extend.py tests/test_obc_design.py
  ```
  **34 passed.**
- Direct invocation of fixture-free tests across extension, design, DEM sources, recipes and namelist handling: **43 passed**, overlapping the pytest selection.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` on all **56** OCTOPUS scripts: **passed**.
- Synthetic arrays and in-memory filesystem/NetCDF mocks reproduced the findings above.
- A broader direct-test attempt importing M2 notebook tests stopped because Matplotlib required a writable configuration/cache location.
- Full-suite and production mesh/FVCOM runs were not repeated. No shared-data files were read and no batch jobs were submitted. Git status remained unchanged, including the pre-existing untracked boundary PNG.

## Verdict

VERDICT: FAIL (0 blocker, 4 major, 8 minor, 2 nit)

### Prompt

```markdown
# Review request, round 10: extending a base mesh outward (fvcom-mesh-tools)

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
Rounds 1-9 and their triage are in docs/extend-tools-review-20261001.md.
Round 9 (your previous answer; 14 findings) was fixed in f05220a; read
`git show f05220a`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 410/415: `mktemp -d` run roots; 383: FMESH_RUN_ROOT or
  `m2_383.<job>`, created with `mkdir` (refused if present); 412: per-case
  lock `mkdir <case>/.running` before RUN_OK is invalidated, removed on exit.
- F3 export/relocate: cleanup only in states staging/done/restored.
- F4 `merge_outer` validates indices with `_indices` before use.
- F5 `node_index_valid` counts bad ids before any cast; elements must be
  integer dtype.
- F6 `merge_outer` requires float64 base nodes and depths.
- F7 `set_value` replaces only the value span (quoted text masked).
- F8 `_check_sponge`: `is None`.
- F9 `trim_lone_corners`: edge counts decremented as removals are chosen;
  empty mesh guarded.
- F10 `rfactor_smooth_free`: bounds applied before the no-live-edge return.
- F11/F12 `read_fort14`: `ndmin=2`, ids 1..N in order, type 3, references
  in range.
- F13 comma guard in 410/415.
- F14 `FVCOM_LIBS` made absolute with `realpath -m` (existence is not
  required there: every job sources common.sh, most never run FVCOM).
Verification after f05220a: full test suite 1161 passed and 1
failed (batch job 123040): `test_r1_failed_prerequisite_invalidates_old_marker`
cut 412 at the old `rm -f "${FMESH_RUN_ROOT` text, which the lock replaced;
the test now starts the 412 block at `CASE_DIR=` and passes (commit after
f05220a). On real data 444 published the Enshu boundary, the full 445 build
passed QA 23/23 with status ok and nothing changed during the build, grd
and dep files bit-identical to rounds 4-8 (trim_lone_corners dropped the
same 2 elements), and 453 completed with status ok.
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
| F2 | major | code read (default outputs/m2_383) | yes | fixed (reserved outputs/<root name>), f236d34 |
| F3 | major | code read (length fallback, no locks in 421/423) | yes | fixed (mktemp root, .staging lock), f236d34; locks mocked |
| F4 | major | test (NaN u, inf w) | yes | fixed, f236d34; test; run on two real histories |
| F5 | minor | fault injection (interrupt during restore) | yes | fixed (state-based cleanup), f236d34 |
| F6 | minor | tests (interrupt after the move) | yes | fixed (record before move), f236d34; tests |
| F7 | minor | tests (comments) | yes | fixed, f236d34; test |
| F8 | minor | tests (string ids, 2-D chain) | yes | fixed, f236d34; tests |
| F9 | minor | tests (ids 0/4, NETA mismatch, 1.9) | yes | fixed, f236d34; test |
| F10 | minor | test (int32 elements) | yes | fixed (int64 required, bytes compared), f236d34; test |
| F11 | minor | test (max_rounds=1) | yes | fixed, f236d34; test |
| F12 | minor | test (3.0000004 floor) | yes | fixed (round_depths_inside), f236d34; test |
| F13 | nit | doc read | yes | fixed, f236d34 |
| F14 | nit | doc read | yes | fixed, f236d34 |

## Round 11

### Answer (verbatim)

Read-only review completed: **2 major and 6 minor findings**, including the previously known GPL policy issue. No files were modified, shared datasets were read, or batch jobs were submitted. The existing untracked PNG remains untouched.

Verification:

- **50 existing tests passed**, executed directly with `python -B` in memory, selecting test functions requiring no filesystem-writing fixtures.
- Additional probes used small synthetic meshes, diskless NetCDF, mocked filesystem operations and shell variable parsing.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n`: **56 OCTOPUS scripts passed**.
- A suspected fractional fort.14 parsing defect was **not reproduced**: the installed NumPy rejected the malformed integers.
- The full suite, production generation and FVCOM integration were **not run**. Production measurements supplied in the request were not independently repeated.

## Previous findings

IDs below retain their numbering within each previous round. Grouped IDs share the stated status and reasons.

The GPL finding is **NOT RESOLVED** in every round: **R1-F27, R2-F22, R3-F14, and R4–R10-F1**. Direct imports remain; owner decision pending. It is counted once as finding 1 below.

### Round 1

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Masked elevations become NaN, allowing fallback. |
| F2 | RESOLVED | Failed required QA gates produce a failure exit. |
| F3 | RESOLVED | Bounds participate in limiting; convergence and final rounded r-factor are checked. |
| F4 | RESOLVED | Floor violations and final timestep allowances follow the explicit owner warning policy. |
| F5 | RESOLVED | Resampling and publication check actual chord lengths against both endpoint floors. |
| F6 | RESOLVED | Missing sizing coverage is rejected at design, wet-lattice and generated-node samples. |
| F7 | RESOLVED | Normal departure segments and the complete boundary receive land-crossing checks. |
| F8 | RESOLVED | Rejected designs preserve published products; exceptional recovery retains protection. |
| F9 | RESOLVED | Native serialization round-trips doubles; exported frozen-base verification runs. |
| F10 | RESOLVED | Opposite-side seam incidence and outer/base overlap are checked. |
| F11 | RESOLVED | Single-edge land runs are retained. |
| F12 | RESOLVED | Unsupported constrained chains receive an explicit six-node requirement. |
| F13 | RESOLVED | CAO cache identity includes resolved source path and dimensions. |
| F14 | RESOLVED | Invalid fine-grid interpolation preserves valid coarse coverage. |
| F15 | RESOLVED | Insufficient distinct points and collinear soundings return uncovered samples. |
| F16 | RESOLVED | CAO provenance inventories include depth files and tables. |
| F17 | RESOLVED | Relevant dirty sources are hashed, including incoming renames. |
| F18 | RESOLVED | Settings, seeds, depth controls and geographic bounds receive numeric validation. |
| F19 | RESOLVED | Nonfinite spacing and nonadvancing resampling steps are rejected. |
| F20 | RESOLVED | Extension output reservation is exclusive. |
| F21 | RESOLVED | Re-depth rejects source/output overlap and reserves its destination. |
| F22 | RESOLVED | Extension smoke staging reserves a fresh root and excludes previous history. |
| F23 | RESOLVED | Depth control precedes smoke timestep selection. |
| F24 | RESOLVED | Smoke paths are resolved before staging. |
| F25 | RESOLVED | The smoke job honors and identifies `FMESH_FVCOM`. |
| F26 | RESOLVED | Job 383 uses the prepared case names and guards directory changes. |
| F28 | RESOLVED | Extension and re-depth case names are restricted to filename components. |
| F29 | RESOLVED | Native reading requires a unique consecutive boundary walk. |

### Round 2

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Re-depth and extension smoke use exclusive reservations. |
| F2 | RESOLVED | Unique temporaries and publication locking prevent competing payload swaps. |
| F3 | RESOLVED | Re-depth repeats frozen-base, overlap and QA checks; sensitivity overrides are labelled. |
| F4 | RESOLVED | Resampling checks chords and both endpoint floors. |
| F5 | RESOLVED | Continuous departure testing detects thin land strips. |
| F6 | RESOLVED | Land-crossing validation permits only numerical endpoint residue. |
| F7 | RESOLVED | Validation uses coordinates reconstructed from the published CSV precision. |
| F8 | RESOLVED | Native coordinates and depths use round-trip-safe serialization. |
| F9 | RESOLVED | Outer/base overlap is checked beyond the seam. |
| F10 | RESOLVED | Porcelain parsing includes both rename paths. |
| F11 | RESOLVED | Geographic bounds are finite, ordered and range-checked. |
| F12 | RESOLVED | Resampling requires finite forward progress. |
| F13 | RESOLVED | Boundary traversal direction follows the supplied second node. |
| F14 | RESOLVED | Final-iteration convergence is accepted. |
| F15 | RESOLVED | Limited depths and limiter controls are validated. |
| F16 | RESOLVED | Both sides are probed and complete guide segments receive geometry checks. |
| F17 | RESOLVED | Bands are materialized and checked under the accepted deviation policy. |
| F18 | RESOLVED | Failure reporting, consumed-input identity and completion inventories/code are handled. |
| F19 | RESOLVED | Uncovered generated nodes are rejected. |
| F20 | RESOLVED | Invalid radii, zero sides and reversals are rejected. |
| F21 | RESOLVED | Refusing smoke-root reuse preserves previous success markers. |

### Round 3

| ID | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Final composed-field violations are reported under owner policy. |
| F2 | RESOLVED | Publication permits numerical spacing tolerance only. |
| F3 | RESOLVED | Materializing bands prevents iterator exhaustion. |
| F4 | RESOLVED | Both-side probing and complete guide checks cover the reported defect. |
| F5 | RESOLVED | Recovery failures retain backups and locking; existing backups block publication. |
| F6 | RESOLVED | Early handlers, input identity and late failure reporting are present. |
| F7 | RESOLVED | Zero-weight masked neighbours preserve valid samples. |
| F8 | RESOLVED | fort.14 values round-trip; supplementary export is checked. |
| F9 | RESOLVED | Pairwise overlap checks avoid global-area dilution. |
| F10 | RESOLVED | Notebook 414 selects its timestep after depth control. |
| F11 | RESOLVED | Design spacing controls must be finite and positive. |
| F12 | RESOLVED | Short designs receive explicit rejection before interior statistics. |
| F13 | RESOLVED | Filleting rejects reversals. |

### Round 4

| ID | Status | Reason |
|---|---|---|
| F2 | RESOLVED | Half-window amplitude and phase derive from interpolated harmonic coefficients. |
| F3 | RESOLVED | Bilinear interpolation preserves query shape. |
| F4 | RESOLVED | Overlap tolerance uses the smaller meeting element. |
| F5 | RESOLVED | Complete guide segments receive land and constrained-line checks. |
| F6 | RESOLVED | Publication rollback preserves the previous coherent set and unresolved recovery storage. |
| F7 | RESOLVED | State initialization and handler registration precede fallible logging. |
| F8 | RESOLVED | Parsed YAML/CSV identities are enforced; completion checks cover inputs and inventories. |
| F9 | RESOLVED | CAO retains valid centres beside missing cells, including final rows/columns. |
| F10 | RESOLVED | Boundary simplicity and repeated adjacent nodes are checked. |
| F11 | RESOLVED | Manual timesteps must be finite and positive. |
| F12 | RESOLVED | Executables are resolved and checked before directory changes. |
| F13 | RESOLVED | Printable ASCII and the 80-byte directory limit are enforced. |
| F14 | RESOLVED | Required barotropic fields have nonempty shape checks; unreadable named grids fail. Finding 4 concerns additional fields. |

### Round 5

| ID | Status | Reason |
|---|---|---|
| F2 | RESOLVED | Flip sites use convexity guards and skip stale ownership rows. |
| F3 | RESOLVED | Recipe and boundary hashes identify the bytes parsed. |
| F4 | RESOLVED | Complete ladder segments receive intersection checks. |
| F5 | RESOLVED | Exceptional recovery retains backups and publication protection. |
| F6 | RESOLVED | Unreadable named grids fail run validation. |
| F7 | RESOLVED | Non-ASCII and overlength directories are refused. |
| F8 | RESOLVED | Failure handlers precede initial post-reservation logging. |
| F9 | RESOLVED | CAO preserves query shape. |
| F10 | RESOLVED | Nonfinite depths fail minimum-depth QA. |

### Round 6

| ID | Status | Reason |
|---|---|---|
| F2 | RESOLVED | Buffered design and generation sea are checked against their land windows. |
| F3 | RESOLVED | The checker honors ordinary and quoted `OUTPUT_DIR` values. |
| F4 | RESOLVED | Benchmark runs undergo completion validation before timing acceptance. |
| F5 | RESOLVED | QA failures and identity mismatches cause failure exits. |
| F6 | RESOLVED | Reference/work-area checks precede destructive staging. |
| F7 | RESOLVED | Unknown recovery outcomes retain backup storage and locking. |
| F8 | RESOLVED | Completion checks cover all recorded scientific inputs and recollect inventories. |
| F9 | RESOLVED | Native writers reject the reported invalid connectivity and nonfinite values. |
| F10 | RESOLVED | Namelist writing preserves doubled delimiters and literal backslashes. |
| F11 | RESOLVED | Relocation checks directory lengths and source/destination overlap. |

### Round 7

| ID | Status | Reason |
|---|---|---|
| F2 | RESOLVED | Whole-dataset M7001 triangulation removes query-composition dependence. |
| F3 | RESOLVED | Namelist reading decodes doubled quoted delimiters. |
| F4 | RESOLVED | Native writers reject fractional, nonfinite and out-of-range indices. |
| F5 | RESOLVED | Prevalidation, staged publication and tracked rollback protect existing case files. |
| F6 | RESOLVED | Relocation refuses source overlap and prepares the copy before replacement. |
| F7 | RESOLVED | Speed verification checks its work area before staging. |
| F8 | RESOLVED | Generation consumes coordinates parsed from the verified CSV bytes. |
| F9 | RESOLVED | Completion inventories are recollected. |
| F10 | RESOLVED | Relevant execution-code identity is compared at completion. |
| F11 | RESOLVED | Marker failure retains locking; existing backups cannot be overwritten. |
| F12 | RESOLVED | Figure publication precedes the re-depth success report. |

### Round 8

| ID | Status | Reason |
|---|---|---|
| F2 | RESOLVED | Existing exports are restored on publication failure; interrupted restoration preserves recovery copies. |
| F3 | RESOLVED | The previous relocation destination is preserved/restored on the reported fault paths. Finding 6 covers an initially absent destination. |
| F4 | RESOLVED | Native export requires integer connectivity. |
| F5 | RESOLVED | Sponge iterators are normalized once and those rows are written. |
| F6 | RESOLVED | OBC types must be whole, non-boolean numbers. |
| F7 | RESOLVED | Code comparison excludes repository-wide unrelated dirt. |
| F8 | RESOLVED | Assignments inside quoted values are ignored; conflicting duplicates fail. |
| F9 | RESOLVED | CAO cache keys include resolved path, size, modification time and shape. |
| F10 | RESOLVED | Grid dimension order is checked and transposed when necessary. |
| F11 | RESOLVED | Sponge values use round-trip-safe formatting. |
| F12 | RESOLVED | Frozen coordinates, depths and connectivity are checked by dtype and bytes. |
| F13 | RESOLVED | Job 413 obtains case names from the manifest. |

### Round 9

| ID | Status | Reason |
|---|---|---|
| F2 | PARTIAL | Unique roots and execution locks fix the original collisions; shared finishing products remain exposed—finding 2. |
| F3 | RESOLVED | Unknown rollback outcomes preserve backups; destructive moves are tracked before execution. |
| F4 | RESOLVED | Merge validates indices before conversion. |
| F5 | PARTIAL | Fractional IDs fail, but accepted float/unsigned IDs can bypass OBC QA—finding 3. |
| F6 | RESOLVED | Merge requires float64 base coordinates and depths. |
| F7 | RESOLVED | Replacement preserves neighbouring assignments and ignores comments. |
| F8 | RESOLVED | Explicit `is None` permits NumPy sponge rows. |
| F9 | RESOLVED | Sequential edge counts prevent mutually supporting deletions; empty meshes are guarded. |
| F10 | RESOLVED | Free-node bounds apply before no-live-edge returns. |
| F11 | RESOLVED | `ndmin=2` supports single-triangle fort.14 files. |
| F12 | RESOLVED | Canonical record IDs, triangle type and node references are checked. |
| F13 | RESOLVED | Chain wrappers reject commas in transmitted values. |
| F14 | RESOLVED | `FVCOM_LIBS` becomes absolute before case-directory changes. |

### Round 10

| ID | Status | Reason |
|---|---|---|
| F2 | RESOLVED | Analysis destinations are reserved and derive from unique run names in the affected paths. |
| F3 | PARTIAL | Scratch roots are unique and root/case locks coordinate 421/423/412; the shared finishing destination remains—finding 2. |
| F4 | PARTIAL | Written nonfinite 3-D fields are rejected; empty additional fields still pass—finding 4. |
| F5 | RESOLVED | Boundary cleanup occurs only after a known outcome; interrupted restoration retains protection. |
| F6 | RESOLVED | Moves are tracked before execution and previous destinations survive the reported interruptions. Finding 6 is the fresh-destination residual. |
| F7 | RESOLVED | `set_value` masks comments, including quotes inside comments. |
| F8 | PARTIAL | Malformed shapes/dtypes fail, but safe-cast rejection silently removes otherwise accepted OBC arrays—finding 3. |
| F9 | RESOLVED | Boundary IDs and totals are checked; invalid indices fail before opening the writer. Finding 5 concerns omitted shape validation. |
| F10 | RESOLVED | Merge requires int64 base connectivity; verification compares dtype and bytes. |
| F11 | RESOLVED | Trimming reports the returned mesh and iteration-limit status. |
| F12 | RESOLVED | Rounded depths remain inside inward-rounded recipe bounds. |
| F13 | RESOLVED | M2 documentation describes the current roots, cases and output locations. |
| F14 | RESOLVED | fort.14 documentation describes its current numeric serialization. |

The generation-nondeterminism hypothesis remains **WITHDRAWN**, consistent with the supplied repeated-run evidence. Owner-authorized timestep warnings are not defects.

## Findings

1. **Major — Known F27 remains open: package code directly imports oceanmesh.**  
   [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332).

   **Evidence:** Direct imports also remain at line 165, in `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672` and `autofinish/directives.py:20`. They conflict with the explicit repository policy. The extension’s subprocess boundary does not address these other imports.

   **Fix:** Implement the owner-selected subprocess or separately licensed plugin boundary and reconcile `THIRD_PARTY_NOTICES.md`. **Previously reported; owner decision pending; counted once.**

2. **Major — Different run-root locks do not protect shared finishing products.**  
   [421_finish_and_run.sh:42](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/421_finish_and_run.sh:42), [finish_depths.py:217](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/finish_depths.py:217).

   **Evidence:** Two `FMESH_SKIP_REFINE=1` workflows can reuse the same accepted refinement with different depth controls. They now receive distinct run roots and therefore distinct `.staging` locks, but both finish into `<refinement>/fvcom_finished` and subsequently stage that shared case.

   An in-memory interleaving of the actual exporter let both exports return successfully. Caller B requested **20 m** depths but subsequently read caller A’s **5 m** depths from its returned paths. The overwritten case passed QA with **zero failed gates**. Staging does not verify that this case belongs to its finishing invocation.

   **Fix:** Finish into a directory owned by the unique run root using `--outdir`, then stage from that directory. Alternatively, lock the shared finishing destination through both publication and staging.

3. **Minor — QA silently treats accepted float and uint64 OBC IDs as no open boundary. Introduced by `f236d34`.**  
   [qa.py:540](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:540).

   **Reproduction:** On the existing pristine fixture, supply the duplicate chain `[4, 8, 4]`:

   | dtype | Reported OBC nodes | Failed gates |
   |---|---:|---:|
   | int64 | 2 | 1 |
   | float64 | 0 | 0 |
   | uint64 | 0 | 0 |

   `casting="safe"` rejects float64 and uint64 conversion regardless of their actual values. The exception handler substitutes an empty OBC set, while `_bad_ids` subsequently accepts these integral, in-range arrays. Ordering, perpendicularity, boundary placement and reachability checks are skipped.

   **Fix:** Validate first, then normalize accepted arrays to int64. A normalization failure must fail integrity validation rather than mean “no open boundary.”

4. **Minor — Empty additional history fields pass the strengthened smoke check.**  
   [check_run.py:292](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:292).

   **Reproduction:** A diskless history with valid times and finite `zeta`, `ua`, `va`, but `u(time=2, siglay=0, nele=1)`, returned **`ok=True`, `reasons=[]`**. The equivalent nonempty NaN `u` correctly failed.

   Every empty record satisfies `np.isfinite(...).all()` vacuously. Thus the new checks can approve a 3-D output field containing no spatial solution.

   **Fix:** Reject zero-length non-time dimensions in written record fields before finiteness checks. Validate known spatial dimensions against the staged mesh where applicable.

5. **Minor — fort.14 shape errors still truncate an existing destination before failing.**  
   [fort14.py:214](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:214).

   **Reproduction:** With valid three-node coordinates/depths and integer connectivity `[[0, 1, 2, 0]]`, `_indices(..., ndim=2)` passes. The writer opens the destination with `"w"`, writes the header and nodes, then raises **`ValueError: too many values to unpack (expected 3)`**.

   A mocked destination initially containing `OLD CASE` was replaced by the partial file. Coordinate/depth shape errors have the same late-failure risk.

   **Fix:** Validate connectivity shape exactly `(NE, 3)`, coordinates `(NP, 2)` and depths `(NP,)` before opening the destination. Stage the complete file before replacing an existing file.

6. **Minor — Interrupted relocation falsely declares restoration when the destination was initially absent.**  
   [fvcom_namelist.py:136](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_namelist.py:136).

   **Reproduction:** Start without `dst`; let `work.rename(dst)` complete, then inject `KeyboardInterrupt`. No `previous` directory exists, so rollback does nothing, sets `state="restored"` and deletes staging storage. The function raises while the newly published destination remains.

   The initially absent destination has not been restored to its original state, despite the explicit state declaration.

   **Fix:** Track publication into a previously absent destination and remove that newly installed copy during rollback. Declare restoration complete only after restoring either the previous copy or previous absence.

7. **Minor — Extension acceptance does not detect missing water or unsupported artificial islands.**  
   [447_extend_merge.py:195](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:195), also the corresponding acceptance path in notebook 453.

   **Reproduction:** From a small regular ocean mesh, remove the two triangles covering an interior 1 km square, entirely outside the frozen base. Derive land boundaries as the pipeline does.

   The resulting mesh had **49 nodes / 70 elements**, passed **23/23 QA gates**, passed the frozen contract for **14 base nodes / 12 base elements / 6 interface edges**, and passed `check_no_overlap`. Its footprint nevertheless omitted **1,000,000 m²** of intended water, represented as an artificial island.

   The acceptance checks validate surviving geometry and base preservation, without comparing new coverage or closed boundary loops against the intended sea/land geometry. This reproduction exercises the acceptance checks; the full generation/finishing pipeline was not run.

   **Fix:** Carry the intended domain into final acceptance. Check missing/excess coverage and require new island boundaries to correspond to supplied land, with explicit tolerances and reporting for permitted coastline/resolution changes.

8. **Minor — Base paths containing spaces break finishing/staging portability.**  
   [421_finish_and_run.sh:56](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/421_finish_and_run.sh:56).

   **Reproduction:** The Python helper prints three paths separated by spaces, and `read` splits them on whitespace. For paths beneath `/work/a b/`, the assignments become:

   ```text
   B_GRD=/work/a
   B_DEP=b/base_grd.dat
   B_OBC=/work/a b/base_dep.dat /work/a b/base_obc.dat
   ```

   A valid accepted refinement therefore cannot be staged when its recorded base paths contain spaces. Quoting later shell expansions cannot repair the earlier split.

   **Fix:** Transfer paths as NUL-delimited entries into an array, checking the producer’s exit status, or perform validation and copying directly in Python.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 6 minor, 0 nit)

### Prompt

```markdown
# Review request, round 11: extending a base mesh outward (fvcom-mesh-tools)

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
Rounds 1-10 and their triage are in docs/extend-tools-review-20261001.md.
Round 10 (your previous answer; 14 findings) was fixed in f236d34; read
`git show f236d34`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 384 `--output` defaults to outputs/<run root name> and is taken with
  `outdir.reserve`; 415 passes that name; 410's message matches 413.
- F3 refine_workflow: `mktemp -d` root (template shortened if the 80-byte
  rule needs it); 421/423: root lock `<root>/.staging` taken with mkdir
  before invalidation, refused while any `<root>/*/.running` exists; 412
  refuses while `.staging` exists (each checks after taking its own lock,
  so a race makes both refuse rather than both proceed).
- F4 check_run: every float variable with a leading time dimension is
  checked record by record (besides the zeta/ua/va shape checks).
- F5/F6 444, export_fvcom_case, relocate_case: explicit states; cleanup only
  in staging/done/restored; each move recorded before it is made.
- F7 `set_value` blanks `!` comments.
- F8 QA `_bad_ids`: shape (1-D chains, (NE, 3) elements) and dtype first.
- F9 `read_fort14`: NOPE/NETA/NBOU/NVEL non-negative and consistent; ids in
  1..NP. `write_fort14` validates indices first. (Weir boundaries, which
  carry node pairs, were never read correctly; they now fail the NVEL check
  rather than being misread.)
- F10 `merge_outer` requires int64 elements; `verify_frozen_base` compares
  elements by dtype and bytes.
- F11 trim report from the returned mesh, plus `round_limit_reached`.
- F12 `round_depths_inside` (bounds moved inward to the 1e-6 grid).
- F13/F14 docs.
Verification after f236d34: full test suite 1172 passed (batch
job 123146); on real data 444 published the Enshu boundary, the full 445
build passed QA 23/23 with status ok and nothing changed during the build,
grd and dep files bit-identical to rounds 4-9, 453 completed with status
ok, and the stricter check_run accepted two existing real FVCOM smoke
histories (97 records each, 3-D fields included).
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
| F2 | major | code read (shared fvcom_finished/) | yes | fixed (finish into the run root), 88bd449 |
| F3 | minor | test (float/uint64 ids) | yes | fixed, 88bd449; test |
| F4 | minor | test (empty siglay) | yes | fixed, 88bd449; test |
| F5 | minor | test (4-column elements over an existing file) | yes | fixed (shape checks, staged write), 88bd449; test. Also found: mkstemp files were published 0600 (here and in 444); now umask mode |
| F6 | minor | test (interrupt, no previous copy) | yes | fixed, 88bd449; test |
| F7 | minor | test (1 km hole) | yes | fixed for 447 (check_island_holes); 453 keeps the 447-accepted geometry, so the check belongs to 447 |
| F8 | minor | mock with a spaced path | yes | fixed, 88bd449 |

## Round 12

### Answer (verbatim)

The review found **3 major, 3 minor and 1 nit**. No files were modified. The pending GPL finding is counted once as an existing finding.

## Previous findings

IDs below use `R<round>-F<finding>`. Grouped IDs share the stated status and resolution reason.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4-F1 through R11-F1 | **NOT RESOLVED** | Package code still imports oceanmesh directly. Owner decision remains pending; finding 1 below. |
| Generation-nondeterminism hypothesis | **WITHDRAWN** | The supplied repeated production runs contradict it; this review found no new supporting evidence. |
| R1-F1, F13–F16 | **RESOLVED** | Masked samples fall through correctly; cache identity, degenerate soundings, coarse-grid preservation and source inventories address the reported cases. |
| R1-F2–F4 | **RESOLVED** | Required QA failures fail the build; depth bounds and final r-factor are checked. Timestep warnings follow the explicit owner policy. |
| R1-F5–F8, F12, F19 | **RESOLVED** | Actual chords, coverage, continuous land intersections, protected publication, minimum chain length and resampling progress are checked. |
| R1-F9–F11, F29 | **RESOLVED** | Native serialization preserves doubles; frozen-base, seam, overlap and boundary-walk checks address the reported defects. Single-edge land runs are retained. |
| R1-F17–F18, F28 | **RESOLVED** | Relevant source identity, recipe controls and case names are validated. |
| R1-F20–F26 | **RESOLVED** | Output reservations, source/output separation, smoke staging, depth-before-timestep ordering, absolute paths, executable selection and case-directory handling are present. |
| R2-F1–F2, F21 | **RESOLVED** | Exclusive reservations and publication locking protect the reported competing-write and reused-root paths. |
| R2-F3 | **PARTIAL** | Re-depth repeats frozen-base, overlap and QA checks, but can bypass the subsequently added island acceptance gate; finding 2. |
| R2-F4–F7, F13, F16–F17, F19–F20 | **RESOLVED** | Chords, departures, published coordinates, traversal direction, guides, bands, coverage and fillet controls address the reported geometry cases. |
| R2-F8–F9 | **RESOLVED** | Native doubles round-trip and overlap is checked beyond the interface. |
| R2-F10–F12, F14–F15, F18 | **RESOLVED** | Rename identity, numeric bounds, resampling progress, final-iteration convergence, limiter inputs and provenance/failure handling are corrected. |
| R3-F1–F4, F11–F13 | **RESOLVED** | Composed-field warnings follow owner policy; spacing, materialized bands, guides, short chains and reversals are checked. |
| R3-F5–F6 | **RESOLVED** | Exceptional recovery retains protection; early handlers and completion identity checks are present. |
| R3-F7–F10 | **RESOLVED** | Zero-weight missing corners preserve samples; numeric fort.14 round-trips, local overlap checks and depth-before-timestep ordering address the original reproductions. |
| R4-F2–F3, F9 | **RESOLVED** | Harmonic-coefficient interpolation and shape-preserving interpolation address the reported cases, including valid edge centres. |
| R4-F4–F5, F10 | **RESOLVED** | Local overlap tolerance, complete-guide checks and boundary simplicity checks are present. |
| R4-F6–F8 | **RESOLVED** | Coherent rollback, early handlers and parsed-input/completion identity checks address the reported failures. |
| R4-F11–F14 | **RESOLVED** | Manual timestep, executable, FVCOM directory and required barotropic-field validation address the reported cases. |
| R5-F2, F4 | **RESOLVED** | Convexity/stale-ownership guards and complete ladder intersection checks are present. |
| R5-F3, F5, F8 | **RESOLVED** | Parsed-byte identity and protected failure handling address the reported cases. |
| R5-F6–F7, F9–F10 | **RESOLVED** | Unreadable grids, invalid directories, lost query shape and nonfinite QA depths are handled. |
| R6-F2, F5, F8 | **RESOLVED** | Land-window checks, failure exits and completion inventories address the reported paths. |
| R6-F3–F4 | **RESOLVED** | The checker honors `OUTPUT_DIR`; benchmark timings now require completion validation. Finding 3 concerns a separate concurrency defect. |
| R6-F6–F7, F11 | **RESOLVED** | Destructive staging prechecks, protected uncertain recovery and relocation path checks are present. |
| R6-F9–F10 | **RESOLVED** | Native writers reject the reported invalid values; namelist writing preserves quoting and literal backslashes. |
| R7-F2 | **RESOLVED** | Whole-dataset M7001 triangulation removes query-composition dependence. |
| R7-F3–F4 | **RESOLVED** | Quoted delimiters are decoded and native indices are validated before conversion. |
| R7-F5–F7, F11 | **RESOLVED** | Export/relocation staging, overlap refusal, work-area prechecks and protected recovery address the reported failures. |
| R7-F8–F10, F12 | **RESOLVED** | Verified CSV coordinates are consumed; inventories/code are rechecked; figure publication precedes re-depth success reporting. |
| R8-F2–F6, F11–F12 | **RESOLVED** | Export and relocation rollback, connectivity/type validation, sponge materialization, numeric serialization and frozen-byte checks address the reported cases. |
| R8-F7–F10, F13 | **RESOLVED** | Relevant code comparison, quoted-assignment handling, cache identity, grid dimension order and manifest case selection are corrected. |
| R9-F2–F3, F13–F14 | **RESOLVED** | Unique workflow roots, locks, isolated finishing, pretracked recovery, comma rejection and absolute library paths address the reported paths. |
| R9-F4–F6 | **RESOLVED** | Merge indices are checked before conversion; accepted OBC IDs are normalized after validation; base coordinates/depths must be float64. |
| R9-F7–F8 | **RESOLVED** | Namelist replacement preserves neighbours and comments; NumPy sponge arrays are accepted explicitly. |
| R9-F9–F12 | **RESOLVED** | Sequential trimming, free-depth bounds, single-triangle reading and canonical fort.14 references address the reported cases. |
| R10-F2–F3 | **RESOLVED** | Analysis destinations are reserved; affected workflow roots/locks and 421’s finishing destination are isolated. |
| R10-F4 | **PARTIAL** | Normal time-leading 3D fields are checked; nonleading-time fields still bypass validation, including NaN fields; finding 4. |
| R10-F5–F7 | **RESOLVED** | Cleanup follows known outcomes, destructive moves are pretracked and namelist comments are masked. |
| R10-F8 | **RESOLVED** | QA validates shapes/dtypes and normalizes accepted boundary arrays after validation. |
| R10-F9–F12 | **RESOLVED** | Reported boundary-index, frozen-connectivity, trimming-report and depth-rounding cases are corrected. |
| R10-F13–F14 | **RESOLVED** | The reported M2 workflow and fort.14 serialization documentation was updated. |

Round 11 specifically:

| Finding | Status | Reason |
|---|---|---|
| F1 | **NOT RESOLVED** | Pending F27; finding 1 below. |
| F2 | **RESOLVED** | 421 finishes under its run root and requires exactly one resulting case. |
| F3 | **RESOLVED** | Valid float/unsigned OBC IDs are validated and normalized; invalid IDs fail integrity validation. |
| F4 | **PARTIAL** | Time-leading empty fields now fail, but fields with a nonleading time dimension still escape the check; finding 4. |
| F5 | **PARTIAL** | The reported nonempty shape-error truncation is fixed by prevalidation and atomic publication. Empty malformed connectivity still bypasses shape validation and publishes an unreadable file; finding 6. |
| F6 | **RESOLVED** | Relocation rollback restores an initially absent destination to absence. |
| F7 | **PARTIAL** | 447 rejects unsupported islands, but 453 does not establish that its source passed that gate; finding 2. |
| F8 | **RESOLVED** | Python checks/copies the base files, preserving spaced paths. |

No additional regression was substantiated as newly introduced by `88bd449`, the portability change or the DATUM change. Findings 2, 4 and 6 expose remaining gaps in the fixes. The authorized band-floor and timestep warning policy is not reported as a defect.

## Findings

1. **Major — Direct GPL imports remain inside the package. Existing F27, still open.**

   **Locations:** [mesh_engine/oceanmesh.py:165](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:165), [mesh_engine/multiscale.py:52](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/multiscale.py:52), [mesh_clean.py:1672](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_clean.py:1672), [autofinish/directives.py:20](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/autofinish/directives.py:20).

   **Evidence:** Searching package source finds eight direct oceanmesh import statements, including generation, multiscale generation, smoothing and remeshing paths. These remain contrary to the repository’s explicit subprocess-only policy.

   **Fix:** Complete the owner’s pending decision and move these execution paths behind external subprocess boundaries. Reconcile the third-party notices with the resulting arrangement.

2. **Major — Re-depth can accept geometry rejected by the island gate.**

   **Locations:** [447_extend_merge.py:182](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:182), [453_redepth_extended.py:108](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:108), [453_redepth_extended.py:150](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:150).

   **Evidence:** 447 exports the native case before calling `check_island_holes` at line 198. A rejected build therefore leaves readable case files. 453 neither requires a successful source-build report nor repeats the island check.

   The synthetic missing-water mesh reproduced the previous result:

   - **49 nodes / 70 elements**, with **1,000,000 m²** missing water.
   - Frozen contract passed: **14 base nodes / 12 base elements / 6 interface edges**.
   - Overlap check passed and QA passed **23/23**.
   - `check_island_holes` rejected the unsupported hole.
   - Executing 453’s actual acceptance code with mocked file I/O reached **`status: "ok"`**, with no problems and no sensitivity override.

   Thus “453 only reads geometry accepted by 447” is an unenforced assumption.

   **Fix:** Require successful build acceptance tied to the exact input case files, or run the island/land gate when accepting an unverified source. A depth-only operation must not promote geometry rejected by the build.

3. **Major — Concurrent benchmarks share mutable cases and can accept the wrong experiment.**

   **Locations:** [416_renumber_benchmark.sh:51](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/416_renumber_benchmark.sh:51), [416_renumber_benchmark.sh:94](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/416_renumber_benchmark.sh:94), [416_renumber_benchmark.sh:110](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/416_renumber_benchmark.sh:110).

   **Evidence:** Every invocation uses `$FMESH_RUN_ROOT/bench`, without a reservation or execution lock. It replaces the same cases, deletes the same histories and writes the same per-repeat logs.

   Executing the actual staging block in memory demonstrated this interleaving:

   - A requests two days and stages `END_DATE=2021-01-03`.
   - B requests one day and replaces those same cases with `END_DATE=2021-01-02`.
   - A subsequently runs/checks B’s namelist.
   - The actual checker accepts B’s one-day history: **`ok=True`, `reasons=[]`**.

   A’s benchmark can therefore report a two-day experiment using a one-day run. Overlapping execution can also delete another invocation’s output.

   **Fix:** Give every benchmark invocation a unique working directory, or hold an exclusive lock throughout staging and all repeats. Coordinate source copying with the staging lock and bind validation to the requested experiment configuration.

4. **Minor — Scientific fields are ignored when their time dimension is not first.**

   **Location:** [check_run.py:285](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:285).

   **Evidence:** The additional-field check skips every variable whose first dimension is not `"time"`.

   Diskless NetCDF probes with valid times and finite `zeta/ua/va` both returned **`ok=True`, `n_records=2`, `reasons=[]`**:

   - `u(siglay=2, time=2, nele=1)` containing only NaNs.
   - `u(siglay=0, time=2, nele=1)` containing no spatial solution.

   This requires a nonstandard history layout, so it is a robustness defect rather than a demonstrated normal FVCOM-output failure.

   **Fix:** Detect the time dimension wherever it occurs and slice records along that axis. Alternatively, explicitly reject unsupported layouts for scientific fields instead of silently skipping them.

5. **Minor — A zero-duration smoke test is staged successfully and can pass completion validation.**

   **Locations:** [448_extend_smoke.py:55](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/448_extend_smoke.py:55), [448_extend_smoke.py:77](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/448_extend_smoke.py:77), [check_run.py:200](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:200).

   **Evidence:** `--days` is not required to be positive. Executing 448’s staging code in memory with `--days 0` produced `START_DATE == END_DATE` and wrote **`STAGED`**.

   Separately, the actual checker accepted that interval with a finite initial record and a clean completion log: **`ok=True`, `n_records=1`, `reasons=[]`**. Subsecond positive durations can also round to the same formatted timestamp. Jobs 416 and 423 have the same unchecked duration calculation.

   **Fix:** Validate finite positive durations before reserving/staging. Require the formatted end time to exceed the start and allow at least one integration step. Reject nonpositive integration intervals in smoke completion validation.

6. **Minor — fort.14 validation still permits successful publication of unreadable files.**

   **Locations:** [fort14.py:216](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:216), [fort14.py:224](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:224), [fort14.py:277](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:277).

   **Evidence:** Actual writer/reader calls with mocked filesystem operations showed:

   - A valid triangle with land boundary type `0.5` publishes normally; reading fails with `invalid literal for int(): '0.5'`.
   - Connectivity shaped `(0, 4)` skips validation because `n_elements == 0`, publishes normally, and then fails reading with an element-block shape mismatch.

   Both calls reached `os.replace`. Atomic writing therefore still replaces an existing destination with these invalid products.

   **Fix:** Validate connectivity shape and dtype independently of its row count. Reject empty meshes if unsupported by the reader. Validate land boundary types as whole, non-boolean integers and serialize their normalized integer values before publication.

7. **Nit — Timestep documentation and option help contradict the warning policy.**

   **Locations:** [USER_GUIDE.md:681](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/docs/USER_GUIDE.md:681), [453_redepth_extended.py:57](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:57), [extend.py:14](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:14).

   **Evidence:** The guide promises that the extension does not limit the timestep; module text repeats that requirement. The re-depth override help describes a failing “time-step gate,” although timestep differences only populate warnings and never require that override.

   **Fix:** Describe sizing as a target and the final timestep comparison as advisory. Restrict override help to the QA failures it actually permits. This finding concerns wording, not the authorized warning policy.

## Verification

- Disk-free pytest runs: **79 passed**—38 extension/OBC tests, 40 selected QA/DEM/checker cases, and one recipe test.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n`: **56 job scripts passed**.
- In-memory probes exercised the findings above and verified source priority, uncovered samples, DATUM counts and non-T.P. warnings.
- Initial pytest attempts failed because temporary-file/capture or matplotlib-cache creation was blocked. A broader selection produced **40 passed / 4 fixture errors** for the same filesystem restriction; the disk-free selection then passed.
- No production mesh generation, full suite, FVCOM integration or batch submission was performed. The owner’s production measurements were not independently repeated.
- Tracked files remained unchanged; the pre-existing untracked boundary PNG remained untouched.

## Verdict

VERDICT: FAIL (0 blocker, 3 major, 3 minor, 1 nit)

### Prompt

```markdown
# Review request, round 12: extending a base mesh outward (fvcom-mesh-tools)

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
Rounds 1-11 and their triage are in docs/extend-tools-review-20261001.md.
Round 11 (your previous answer; 8 findings) was fixed in 88bd449; read
`git show 88bd449`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 421: `finish_depths --outdir <root>/finished`; exactly one case there.
- F3 QA `_as_ids`: validate, then normalise; invalid ones fail the gate.
- F4 check_run: empty dimension in a record field fails; node/nele match.
- F5 `write_fort14`: shapes checked first; written via a temporary file and
  `os.replace`; mkstemp files (also in 444) given the umask mode.
- F6 `relocate_case`: rollback restores "no destination" too.
- F7 new `extend.check_island_holes` in 447: a closed boundary loop touching
  new nodes, without open-boundary nodes, must intersect the supplied land
  (land_with_base) over a positive area. Not in 453: it re-depths the mesh
  447 accepted and does not change its nodes or elements (it verifies the
  frozen base and reads the built case as is). A coverage-area comparison
  against the intended sea was not added: the coast is resolved at the
  mesh size, so area differences there are expected; the hole check is the
  gate.
- F8 421: base files checked and copied in Python.
Verification after 88bd449: full test suite 1179 passed (batch
job 123197); on real data 444 published the Enshu boundary (files now
0644), the full 445 build passed QA 23/23 with status ok, the island check
found 9 islands in the new part, each on supplied land, grd and dep files
bit-identical to rounds 4-10, 453 completed with status ok, and check_run
still accepts the two real FVCOM smoke histories.
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
| F2 | major | code read (export before the check; 453 assumed acceptance) | yes; our round-11 rebuttal withdrawn | fixed, 08e8fd9 |
| F3 | major | code read (shared bench/) | yes | fixed (mktemp -d), 08e8fd9 |
| F4 | minor | test (time second) | yes | fixed (refused), 08e8fd9; test |
| F5 | minor | test (0 days); first attempt compared with NC_FIRST_OUT and broke an earlier test, corrected to START_DATE | yes | fixed (end_after), 08e8fd9; tests |
| F6 | minor | tests (type 0.5, (0, 4)) | yes | fixed, 08e8fd9; test |
| F7 | nit | doc read | yes | fixed, 08e8fd9 |

## Round 13

### Answer (verbatim)

Read-only review completed. No files changed; the pre-existing untracked boundary PNG remains untouched. Production generation, FVCOM integration and batch jobs were not run.

## Previous findings

Grouped ranges below include every numbered finding from rounds 1–12. Statuses concern the reported defects; distinct remaining gaps are numbered separately below.

| Previous finding(s) | Status | Reason |
|---|---|---|
| R1-F27, R2-F22, R3-F14, R4–R12-F1 | **NOT RESOLVED** | Direct oceanmesh imports remain. Owner decision pending; counted once below. |
| Generation-nondeterminism hypothesis | **WITHDRAWN** | Supplied repeated-run evidence contradicts it; no new supporting evidence found. |
| R1-F1, F13–F16 | **RESOLVED** | Masks, source-root/cache identity, degenerate soundings, fallback preservation and depth-file inventories are handled. |
| R1-F2–F4 | **RESOLVED** | QA failures fail execution; depth bounds and final r-factor are checked. Timestep warnings follow the owner’s explicit policy. |
| R1-F5–F8, F12, F19 | **RESOLVED** | Chord spacing, coverage, continuous departure/crossing checks, protected publication, chain length and resampling progress address the reported cases. |
| R1-F9–F11, F29 | **RESOLVED** | Serialization, frozen-base checks, seam incidence, single-edge land runs and ordered boundary validation are corrected. |
| R1-F17–F18, F28 | **RESOLVED** | Relevant code identity, recipe numeric controls and case names are validated. |
| R1-F20–F26 | **RESOLVED** | Reservations, source/output separation, fresh smoke roots, depth-before-timestep ordering, absolute paths, executable overrides and case names are corrected. |
| R2-F1–F3, F18, F21 | **RESOLVED** | Reservation/publication protection, re-depth QA, failure reports and smoke-root protection address the original reproductions. |
| R2-F4–F9, F13, F16–F17, F19–F20 | **RESOLVED** | Reported chord, departure, precision, overlap, traversal, guide, band, coverage and fillet cases are handled. |
| R2-F10–F12, F14–F15 | **RESOLVED** | Rename identity, numeric validation, resampling progress and limiter validation/convergence are corrected. |
| R3-F1–F4, F11–F13 | **RESOLVED** | Composed-field diagnostics, spacing controls, reusable bands, side selection, short chains and reversals are handled. |
| R3-F5–F10 | **RESOLVED** | Publication recovery, early failure handling, interpolation, fort.14 precision, overlap and timestep ordering address the reported cases. |
| R4-F2–F5, F9–F10 | **RESOLVED** | Reported interpolation, local overlap, guide-placement and boundary-simplicity cases are corrected. |
| R4-F6–F8, F11–F14 | **RESOLVED** | Rollback, failure reporting, input rechecks, timestep/executable/path validation and nonempty barotropic arrays address the original cases. |
| R5-F2–F10 | **RESOLVED** | Flip guards, parsed-byte identity, complete guide intersections, recovery protection, grid/path/query validation, early handlers and finite depths are corrected. |
| R6-F2–F6, F8, F11 | **RESOLVED** | Land windows, output-directory handling, completion checks, failure exits, destructive-staging prechecks, input inventories and relocation validation are present. |
| R6-F7, F9–F10 | **RESOLVED** | Reported recovery, native-value and Fortran-string cases are corrected. |
| R7-F2–F4, F8–F10, F12 | **RESOLVED** | Whole-source triangulation, quote/index validation, verified CSV consumption, inventories/code rechecks and figure-before-success ordering are present. |
| R7-F5–F7, F11 | **RESOLVED** | Staged export/relocation, overlap refusal, work-area prechecks and protected recovery address the reported failures. |
| R8-F2–F13 | **RESOLVED** | Reported rollback, connectivity/type, sponge, identity, namelist, cache, dimension-order, serialization and manifest-selection cases are corrected. |
| R9-F2–F14 | **RESOLVED** | Original root collisions, execution/recovery races, index/dtype problems, namelist/sponge handling, trimming, free-depth bounds, reading and library-path cases are corrected. |
| R10-F2–F14 | **RESOLVED** | Original destination/finishing collisions, additional-field validation, recovery tracking, comment parsing, QA normalization, serialization and documentation cases are corrected. |
| R11-F2–F8 | **RESOLVED** | Finishing is isolated; IDs are normalized; empty record dimensions and malformed connectivity fail; relocation restores absence; island checks run in both acceptance paths; spaced base paths are handled. |
| R12-F2 | **PARTIAL** | Static status/hash rejection and the island gate are enforced, but acceptance is not bound to the subsequently consumed source snapshot; finding 6. |
| R12-F3 | **RESOLVED** | Separate benchmarks no longer share their destinations. Findings 2–3 concern remaining source-copy coordination and a new path-length regression. |
| R12-F4 | **RESOLVED** | The reported additional floating fields with time second now fail. Required barotropic fields have a separate bypass; finding 4. |
| R12-F5 | **PARTIAL** | Zero/subsecond formatted intervals are rejected, but initial-only history can still certify a short positive run; finding 5. |
| R12-F6 | **RESOLVED** | Empty/malformed connectivity and fractional land-boundary types are rejected before publication. Finding 8 concerns another unchecked input. |
| R12-F7 | **RESOLVED** | Guide, module text and option help now describe timestep differences as advisory. |

## Findings

1. **Major — Direct GPL imports remain inside the package. Existing F27, still open.**

   **Locations:** [mesh_engine/oceanmesh.py:165](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:165), `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672`, `autofinish/directives.py:20`.

   **Evidence:** Package source still contains direct oceanmesh imports in generation, multiscale, smoothing and remeshing paths, contrary to the repository’s explicit subprocess-only rule.

   **Fix:** Complete the owner’s pending decision, move these operations behind subprocess boundaries or into a separately licensed plugin, and reconcile the notices. This is not a new finding.

2. **Major — Benchmark source copying can mix two experiments.**

   **Locations:** [416_renumber_benchmark.sh:61](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/416_renumber_benchmark.sh:61), [414_refine_m2_prep.sh:23](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/414_refine_m2_prep.sh:23).

   **Evidence:** 416 copies `base` and `refined` separately without coordinating with the source root’s staging lock. 414 likewise does not participate in that lock.

   Executing 416’s actual staging block with an in-memory relocation mock reproduced:

   - Copy experiment A’s base at **5 m**.
   - Restage the source pair as experiment B at **20 m**.
   - Copy B’s refined case.
   - Staging finishes normally with **A/base + B/refined**.

   Both small synthetic cases independently passed **23/23 QA gates**, each with **16 nodes / 18 elements**; implied dt was **142.7843 s** and **71.3922 s** respectively. Completion checks do not establish that the pair belongs to one experiment. FVCOM execution was not performed.

   **Fix:** Hold the root staging lock across both source copies, and make every writer—including 414—honor that protocol and active-run locks. Record the copied pair’s input identity.

3. **Minor — Unique benchmark directories break the default chain’s FVCOM paths. Introduced by `08e8fd9`.**

   **Location:** [416_renumber_benchmark.sh:54](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/416_renumber_benchmark.sh:54).

   **Reproduction:** With the root produced by 415:

   ```text
   /octfs/work/G16445/v61021/scratch/m2r_20261002_120000_abcd
   ```

   The old `bench/base/output/` path is **77 bytes**. The new `bench_abcd/base/output/` path is **82 bytes**. Actual `check_fvcom_dirs` accepts the former and rejects the latter. The original staged `refined/output/` is only **74 bytes**, so the source is valid.

   **Fix:** Reserve a shorter unique benchmark root directly beneath the scratch directory, and validate its longest final FVCOM path before staging.

4. **Minor — Required barotropic fields bypass dimension-order validation.**

   **Location:** [check_run.py:266](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:266).

   **Reproduction:** Diskless history with three valid timestamps, finite `zeta(node,time)`, `ua(nele,time)` and `va(nele,time)`, each shaped `(3,3)`, returned:

   ```text
   ok=True, n_records=3, reasons=[]
   ```

   Shape checks cannot distinguish these square transposed arrays. The strengthened additional-field loop explicitly excludes these three variables.

   **Fix:** Require their dimensions to be exactly `("time","node")` and `("time","nele")`, or explicitly normalize supported layouts before validation.

5. **Minor — A short positive smoke interval can pass with no integration output.**

   **Locations:** [check_run.py:233](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:233), [check_run.py:333](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:333).

   **Reproduction:** With start `00:00`, end `00:10`, a 1,800-second output interval, a clean `TADA` log and only the finite **initial** record at `00:00`, the actual checker returned **`ok=True`, `n_records=1`, `reasons=[]`**. The end exceeds the start, but the whole requested run fits inside the completion tolerance.

   A second probe omitted `START_DATE`, set `NC_FIRST_OUT == END_DATE`, and likewise passed one record.

   **Fix:** Require a valid start and positive interval, and require output evidence after the start. Preserve legitimate final-only output where `NC_FIRST_OUT == END_DATE > START_DATE`. Staging should also require at least one integration step.

6. **Minor — Re-depth checks one source version and records/consumes another. Introduced by `08e8fd9`.**

   **Locations:** [453_redepth_extended.py:104](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:104), [453_redepth_extended.py:120](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:120).

   **Reproduction:** Executing the actual acceptance/provenance block in memory, replace the source report and three case files once—after their acceptance checks, before `collect`.

   The block accepted report A with `status="ok"`, then captured report B with `status="failed"` and B’s case files. The final actual `changed_files` check returned **`[]`**, because it compared B against B. The recipe digest guard does not bind the accepted report or case hashes to that snapshot.

   **Fix:** Parse/hash the report from one byte read and verify its digest against provenance; compare the provenance case hashes with the accepted product hashes before reading the mesh. Prefer a locked or immutable source snapshot.

7. **Minor — Newly added land sidecars escape re-depth’s input-change check. Introduced by `08e8fd9`.**

   **Location:** [453_redepth_extended.py:204](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:204).

   **Evidence:** Land’s shapefile inventory is captured initially, but the final `changed_inventory` call relists only bathymetry.

   With the actual provenance helpers, adding `land_with_base.prj` after capture produced **`changed=[]`** under 453’s check. Including a fresh land inventory correctly returned **`["land"]`**. A sidecar appearing before `gpd.read_file` can affect the consumed dataset without appearing in its provenance.

   **Fix:** Relist `dataset_files(LAND)` in the final inventory comparison, alongside bathymetry.

8. **Minor — A multiline title still publishes an unreadable fort.14.**

   **Location:** [fort14.py:259](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fort14.py:259).

   **Reproduction:** A valid three-node triangle with `title="title\nextra line"` passed `write_fort14` validation and reached the mocked `os.replace`. Actual read-back of the emitted text failed:

   ```text
   ValueError: invalid literal for int() with base 10: 'extra'
   ```

   Atomic publication can therefore replace an existing usable destination with an unreadable product.

   **Fix:** Validate the title as a single text line before staging; reject embedded line separators.

9. **Minor — Public sizing helpers accept invalid controls that yield negative sizes or fail to terminate.**

   **Location:** [extend.py:82](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:82).

   **Reproduction:** On a **2×2** lattice, positive ambient sizes and a band filled with `-10` returned **`h=-10` everywhere**, reporting deviation `-0.0`.

   Separately, `grade=-0.2` repeatedly lowers values around negative-cost cycles in `_limit`. A bounded fault injection stopped after 100 heap pops; popped values fell from **1000** to **−1800.1429**, with processing continuing.

   Recipe validation protects the normal driver, but these exported helpers lack equivalent checks.

   **Fix:** Validate finite nonnegative gradation, matching finite lattice arrays, positive ambient/band sizes and valid floor values before entering `_limit`. Preserve NaN solely as the intentional off-band sentinel.

## Verification

- Disk-free pytest selection: **50 passed in 1.81 s**, covering extension/OBC geometry, selected DEM interpolation/datum cases, recipe completeness and namelist controls.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n`: **56 scripts passed**.
- Findings reproduced using small arrays, diskless NetCDF, in-memory snapshots and bounded fault injections.
- Initial direct pytest attempts ran no tests because matplotlib required a writable cache. The successful selection disabled the optional oceanmesh import.
- Full-suite and production measurements supplied by the owner were not independently repeated. No large shared datasets were read.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 7 minor, 0 nit)

### Prompt

```markdown
# Review request, round 13: extending a base mesh outward (fvcom-mesh-tools)

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
Rounds 1-12 and their triage are in docs/extend-tools-review-20261001.md.
Round 12 (your previous answer; 7 findings) was fixed in 08e8fd9; read
`git show 08e8fd9`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 447: island check before any file is written. 453: requires the
  source build's report.json status "ok" and the case files' sha256 equal
  to its products_sha256; runs check_island_holes on
  <build>/generate/land_with_base.shp (which is now in its provenance).
  Our round-11 claim that 453 only reads accepted geometry was right to be
  rejected; it is now enforced.
- F3 416: `BENCH=$(mktemp -d "$RUN_ROOT/bench_XXXX")`.
- F4 check_run: a float field with time not first is a failure reason.
- F5 `io.fvcom_namelist.end_after(start, days)` used by 416/423/448 (448
  before `reserve`); check_run: END_DATE must be after START_DATE (not
  NC_FIRST_OUT, which may equal END_DATE legitimately).
- F6 `write_fort14`: (NE, 3) integer elements for any NE, empty mesh
  refused, land boundary types whole numbers.
- F7 wording in USER_GUIDE, extend.py, 447, 453.
Verification after 08e8fd9: full test suite 1186 passed (batch
job 123303); on real data 444 published the Enshu boundary, the full 445
build passed QA 23/23 with status ok and 9 new islands each on land, grd
and dep bit-identical to rounds 4-11; 453 accepted that build (report ok,
hashes matched), found the same 9 islands, and completed with status ok;
check_run still accepts the two real FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F2 | major | code read (copies outside the lock; 414 outside the protocol) | yes | fixed (fmesh_stage_lock), 92ea195; lock mocked |
| F3 | minor | path length arithmetic | yes | fixed (b416_XXXX under scratch), 92ea195 |
| F4 | minor | test (zeta(node, time)) | yes | fixed, 92ea195; test; fixtures moved to FVCOM's layout |
| F5 | minor | test (10-minute run, initial record only) | yes | fixed, 92ea195; test |
| F6 | minor | code read (report read twice) | yes | fixed, 92ea195 |
| F7 | minor | code read (land not relisted) | yes | fixed, 92ea195 |
| F8 | minor | test (two-line title) | yes | fixed, 92ea195; test |
| F9 | minor | tests (negative grade, negative band, NaN ambient) | yes | fixed, 92ea195; tests |

## Round 14

### Answer (verbatim)

Read-only review completed. **3 major, 5 minor and 1 nit findings remain.** The pending GPL finding is counted once as previously known. No files were modified.

## Previous findings

Grouped ranges below cover every numbered finding from rounds 1–13. “RESOLVED” refers to the reported defect; separate remaining gaps appear below.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R13-F1 | **NOT RESOLVED** | Direct oceanmesh imports remain. Owner decision pending; finding 1 below. |
| Generation-nondeterminism hypothesis | **WITHDRAWN** | Supplied repeated production runs contradict it; no new supporting evidence found. |
| R1-F1, F13–F16 | **RESOLVED** | Mask handling, source/cache identity, fallback preservation, degenerate soundings and source inventories address the reported cases. |
| R1-F2–F4 | **RESOLVED** | Required QA failures fail execution; depth bounds and final r-factor are checked. Timestep warnings follow the owner’s explicit policy. |
| R1-F5–F8, F12, F19 | **RESOLVED** | Actual chords, coverage, continuous land checks, protected publication, chain length and resampling progress are checked. |
| R1-F9–F11, F29 | **RESOLVED** | Serialization, frozen-base verification, seam/overlap checks, single-edge land runs and boundary traversal are corrected. |
| R1-F17–F18, F28 | **RESOLVED** | Relevant code identity, recipe controls and case names are validated. |
| R1-F20–F26 | **RESOLVED** | Reservations, source/output separation, fresh smoke roots, depth-before-timestep ordering, absolute paths, executable overrides and case names address the reported defects. |
| R2-F1–F3, F18, F21 | **RESOLVED** | Reservation/publication protection, re-depth acceptance, failure reporting and reuse refusal address the original reproductions. |
| R2-F4–F9, F13, F16–F17, F19–F20 | **RESOLVED** | Reported spacing, departure, precision, overlap, traversal, guide, band, coverage and fillet cases are handled. |
| R2-F10–F12, F14–F15 | **RESOLVED** | Rename identity, geographic controls, resampling progress and limiter validation/convergence are corrected. |
| R3-F1–F4, F11–F13 | **RESOLVED** | Composed-field diagnostics, spacing controls, reusable bands, side selection, short chains and reversals are handled. |
| R3-F5–F10 | **RESOLVED** | Publication recovery, early failure handling, interpolation, fort.14 precision, overlap and timestep ordering address the reported cases. |
| R4-F2–F5, F9–F10 | **RESOLVED** | Reported interpolation, local overlap, guide-placement and boundary-simplicity cases are corrected. |
| R4-F6–F8, F11–F14 | **RESOLVED** | Rollback, failure reporting, input rechecks, timestep/executable/path validation and required-field checks address the original cases. |
| R5-F2–F10 | **RESOLVED** | Flip guards, parsed-byte identity, complete guide intersections, protected recovery, grid/path/query validation, early handlers and finite depths are corrected. |
| R6-F2–F11 | **RESOLVED** | Land windows, output-directory handling, completion checks, failure exits, staging prechecks, recovery, inventories, native values, Fortran strings and relocation validation address the reported cases. |
| R7-F2–F12 | **RESOLVED** | Whole-source triangulation, quoting/index checks, staged export/relocation, work-area checks, verified CSV consumption, identity rechecks and publication ordering are corrected. |
| R8-F2–F13 | **RESOLVED** | Reported rollback, connectivity/type, sponge, identity, namelist, cache, dimension-order, serialization and manifest-selection defects are corrected. |
| R9-F2–F14 | **RESOLVED** | Original root collisions, recovery races, index/dtype defects, namelist/sponge handling, trimming, depth bounds, reading and library paths are corrected. Finding 2 concerns another staging entry point. |
| R10-F2–F14 | **RESOLVED** | Original destination/finishing collisions, field validation, recovery tracking, comment parsing, QA normalization, serialization and documentation cases are corrected. |
| R11-F2–F8 | **RESOLVED** | Finishing is isolated; IDs are normalized; malformed/empty fields and connectivity fail; relocation restores absence; re-depth acceptance includes the island gate; spaced paths work. |
| R12-F2 | **RESOLVED** | 453 reads acceptance once and binds the consumed report and case hashes to accepted provenance. |
| R12-F3 | **RESOLVED** | Benchmark destinations are unique. Finding 6 concerns the new destination’s parent directory. |
| R12-F4 | **RESOLVED** | Additional floating fields with a nonleading time dimension are rejected. |
| R12-F5 | **RESOLVED** | Invalid formatted durations fail; initial-only output cannot certify completion. |
| R12-F6–F7 | **RESOLVED** | Empty/malformed fort.14 connectivity and fractional boundary types fail; timestep documentation follows owner policy. |
| R13-F2 | **RESOLVED** | 414 participates in the staging protocol; 416 holds the lock across both copies. |
| R13-F3 | **RESOLVED** | Benchmark copies use a shorter path; relocation checks FVCOM’s directory limit. |
| R13-F4–F5 | **RESOLVED** | Required barotropic dimensions are exact; START_DATE is required and output must advance beyond it. |
| R13-F6–F7 | **RESOLVED** | Acceptance is read once and hash-bound; land datasets are relisted at completion. |
| R13-F8 | **RESOLVED** | fort.14 titles must occupy one line. |
| R13-F9 | **PARTIAL** | Negative gradation, invalid ambient sizes and invalid bands are rejected. `graded_up` still accepts nonfinite lattice coordinates; finding 5. The fix also introduced finding 4. |

## Findings

1. **Major — Known F27: package code still imports oceanmesh directly.**  
   **Location:** `src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332`; also `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672`, `autofinish/directives.py:20`; `THIRD_PARTY_NOTICES.md:37`.

   **Evidence:** These execution paths contain direct oceanmesh imports. The third-party notices explicitly endorse this arrangement, contrary to the repository’s stated subprocess separation policy. The extension subprocess path does not remove the other imports.

   **Fix:** Complete the pending owner decision, move these operations behind subprocess boundaries or into a separately licensed plugin, and reconcile the notices. **Previously reported; not new.**

2. **Major — Job 411 can overwrite inputs while job 412 is running.**  
   **Location:** `jobs/octopus/411_m2_prep.sh:25`, `:33`; `notebooks/383_m2_case_prep.py:321`.

   **Evidence:** 411 checks only for existing `m2_*.nc` histories and then prepares cases without acquiring `.staging` or checking `.running`. In a realistic interleaving, 412 acquires `<case>/.running` and starts initialization before its first history exists; 411’s guard passes and 383 overwrites grid, depth, forcing and namelist files. Two simultaneous 411 preparations are likewise uncoordinated.

   An in-memory shell harness executing the actual 411 body with filesystem operations mocked reached preparation successfully with **zero staging-lock calls**. Job 412’s protection at lines 21–25 depends on writers participating in the protocol.

   **Fix:** Acquire `fmesh_stage_lock "$RUN_ROOT"` before inspecting or writing the cases, and hold it through preparation.

3. **Major — Job 414 preserves obsolete staging and smoke-success markers when replacing cases.**  
   **Location:** `jobs/octopus/414_refine_m2_prep.sh:23–33`; consumer: `jobs/octopus/412_m2_run.sh:37`.

   **Reproduction:** Stage and smoke-test a root through 421/423, then rerun 414 with different base/refined inputs before the long integrations start. The production case output directories can still be empty, so 414 permits replacement. It leaves `STAGED` and `SMOKE_OK` intact. Consequently, 412 with `FMESH_REQUIRE_SMOKE=1` accepts the replacement cases using the previous cases’ smoke result.

   The mocked actual 414 body completed preparation with one lock call and the pre-existing smoke marker unchanged. Preparation also writes the two cases sequentially, so a mid-preparation failure can retain old acceptance markers over a mixed state.

   **Fix:** After acquiring the lock and passing reuse guards, invalidate staging/smoke/run acceptance markers before the first write. Publish staging acceptance only after complete preparation; bind smoke acceptance to the staged inputs.

4. **Minor — The round 13 gradation checks reject valid NumPy scalar values.**  
   **Location:** `src/fvcom_mesh_tools/extend.py:55`, `:89`.

   **Reproduction:** Both calls below fail with “the gradation must be finite and non-negative,” despite a valid finite positive value:

   ```python
   graded_up(values, x, y, np.float32(0.2))
   compose_sizing(values, x, y, grade=np.float32(0.2))
   ```

   `np.int64(1)` has the same type-check problem. These values supported the arithmetic before the new `isinstance(grade, (int, float))` checks.

   **Fix:** Accept real numeric scalars, including NumPy scalars, normalize them to `float`, and then validate finiteness and sign. Explicitly reject booleans if intended.

5. **Minor — Sizing input validation remains incomplete and inconsistently normalizes inputs.**  
   **Location:** `src/fvcom_mesh_tools/extend.py:53–57`, `:91–94`, `:110–112`.

   **Reproduction:** With:

   ```python
   values = np.array([[1., 100.], [1., 100.]])
   x = np.array([[0., np.nan], [0., np.nan]])
   y = np.array([[0., 0.], [1., 1.]])
   ```

   `graded_up(values, x, y, .2)` returns the unchanged field successfully. A NaN size also survives in its returned field.

   Separately, `compose_sizing` validates converted coordinate arrays but passes the original objects to `_limit`: valid nested-list coordinates raise `TypeError`. Matching one-dimensional inputs pass validation and then raise an internal unpacking error.

   **Fix:** Normalize arrays once and share explicit matching, finite, two-dimensional lattice validation between both helpers. Reject NaN/+infinity size bounds in `graded_up`, while preserving the intentional negative-infinity sentinel used for off-band lower bounds.

6. **Minor — The new benchmark scratch location assumes its parent already exists.**  
   **Location:** `jobs/octopus/416_renumber_benchmark.sh:54`.

   **Evidence:** `mktemp -d "$WORK_DIR/scratch/b416_XXXX"` creates only its final directory. Neither 416 nor `common.sh` creates `$WORK_DIR/scratch`. A valid staged root elsewhere, with a configured WORK_DIR containing the conda installation but no scratch directory, therefore fails with `ENOENT`. The previous location beneath the existing run root did not have this prerequisite.

   **Fix:** Create and validate `$WORK_DIR/scratch` before calling `mktemp`.

7. **Minor — Run validation truncates fractional timestamps before checking cadence.**  
   **Location:** `src/fvcom_mesh_tools/cli/check_run.py:332`.

   **Reproduction:** An in-memory history with correctly shaped finite fields and timestamps at **0, 0.25, 0.5, 0.75 and 1 second**, matching `NC_OUT_INTERVAL='seconds = 0.25'`, is rejected solely because “the history times do not increase record by record.”

   Conversion to `datetime64[s]` turns several distinct records into identical timestamps. The parser accepts fractional timestamps and intervals, so this rejection contradicts accepted inputs.

   **Fix:** Calculate differences without discarding fractional seconds, using datetime subtraction or microsecond-resolution NumPy timestamps converted to seconds.

8. **Minor — Pure test selections unnecessarily initialize the optional oceanmesh backend.**  
   **Location:** `tests/conftest.py:29–30`.

   **Evidence:** Collection always calls `_importable("oceanmesh")`, even when no selected test has `needs_oceanmesh`. In this read-only environment, that import initializes matplotlib and aborts otherwise independent test collection because its cache cannot be created. An in-memory fault injection likewise showed an optional-backend initialization `OSError` escaping collection containing only unmarked tests.

   **Fix:** First check whether any selected test requires oceanmesh. Probe/import the backend only when such tests exist.

9. **Nit — Job 401’s usage and experiment description refer to different jobs/configurations.**  
   **Location:** `jobs/octopus/401_verify_best.sh:12–20`.

   **Evidence:** Its usage instructs submitting `395_target_mode.sh`. Its header describes `field`, `achieved` and `ach1w`, while the actual invocations at lines 71–72 run `best` and `legacy`.

   **Fix:** Update the submission command and header to describe job 401’s actual two variants.

## Verification

- **47 selected tests passed**, covering `test_extend.py`, `test_obc_design.py`, `test_dem_sources.py` and `test_extend_recipe.py`, using disk-free tests with:
  `pytest.main(["-q", "-s", "--noconftest", "-p", "no:cacheprovider", *selected_nodes])`.
- **All 56 OCTOPUS shell scripts passed `bash -n`.**
- **Focused Ruff checks passed** with `ruff check --no-cache` on the reviewed extension modules, notebooks and affected CLI modules.
- In-memory sizing, history-validation, optional-import fault injection and shell-body mocks reproduced the findings described above.
- Initial pytest attempts failed because ordinary capture required writable temporary storage, and then because unconditional optional-backend initialization required a writable matplotlib cache. The successful selection avoided both.
- No production mesh generation, real-data depth processing or FVCOM integrations were run. The supplied production measurements were not independently repeated.
- Git status remained unchanged, including the pre-existing untracked `recipes/extend/tokyo_bay_enshu_obc.png`. No batch jobs were submitted.

## Verdict

VERDICT: FAIL (0 blocker, 3 major, 5 minor, 1 nit)

### Prompt

```markdown
# Review request, round 14: extending a base mesh outward (fvcom-mesh-tools)

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
Rounds 1-13 and their triage are in docs/extend-tools-review-20261001.md.
Round 13 (your previous answer; 9 findings) was fixed in 92ea195; read
`git show 92ea195`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 `common.sh` `fmesh_stage_lock`/`fmesh_stage_unlock` (same `.staging`
  protocol as 421/423): 416 holds it across both copies; 414 takes it.
- F3 416: `mktemp -d "$WORK_DIR/scratch/b416_XXXX"`; relocate_case checks
  the 80-byte limit.
- F4 check_run: zeta (time, node), ua/va (time, nele) exactly.
- F5 check_run: START_DATE required; last record must be after it.
- F6 453: one read of report.json; provenance hashes of the report and of
  the case (combined as `collect` combines them) must equal the accepted
  ones.
- F7 453: `dataset_files(LAND)` relisted at completion.
- F8 `write_fort14`: one-line title.
- F9 `compose_sizing`/`graded_up` input checks.
Verification after 92ea195: full test suite 1192 passed (batch
job 123355); on real data 444, the full 445 build (QA 23/23, status ok, 9
islands on land, grd and dep bit-identical to rounds 4-12) and 453 (status
ok) all passed, and the stricter check_run (exact zeta/ua/va dimensions,
START_DATE, output after it) still accepts the two real FVCOM smoke
histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F2 | major | code read (411 without the lock) | yes | fixed, ea04210 |
| F3 | major | code read (old markers kept) | yes | fixed, ea04210 |
| F4 | minor | test (np.float32 grade) | yes, introduced in round 13 | fixed (_grade), ea04210; test |
| F5 | minor | tests (NaN x, lists, 1-D) | yes | fixed (_lattice), ea04210; test |
| F6 | minor | code read | yes | fixed, ea04210 |
| F7 | minor | test (0.25 s records) | yes | fixed, ea04210; test |
| F8 | minor | code read | yes | fixed, ea04210 |
| F9 | nit | doc read | yes | fixed, ea04210 |

## Round 15

### Answer (verbatim)

Review baseline: `554241b`, including `ea04210`. Concurrent uncommitted relicensing edits appeared during the review; F27 remains open under your explicit instructions. I modified no files.

The statuses below distinguish resolution of the original reproduction from separate defects found this round.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R14-F1 | NOT RESOLVED | Direct oceanmesh imports remain; owner decision pending. Finding 1 below counts this once. |
| Generation-nondeterminism hypothesis | WITHDRAWN | Supplied repeated production runs contradict it; no new supporting evidence found. |
| R1-F1, F13–F16 | RESOLVED | Mask handling, source/cache identity, fallback preservation, degenerate soundings and source inventories address the reported cases. |
| R1-F2–F4 | RESOLVED | Required QA failures fail execution; depth bounds and final r-factor are checked. Timestep warnings follow the owner’s policy. |
| R1-F5–F8, F12, F19 | RESOLVED | Chords, coverage, continuous land checks, protected publication, chain length and resampling progress address the original cases. |
| R1-F9–F11, F29 | RESOLVED | Serialization, frozen-base verification, seam/overlap checks, single-edge land runs and boundary traversal are corrected. |
| R1-F17–F18, F28 | RESOLVED | Relevant code identity, recipe controls and case names are validated. Finding 7 concerns the separate recipe `name`. |
| R1-F20–F26 | RESOLVED | Reservations, source/output separation, fresh smoke roots, depth-before-timestep ordering, absolute paths, executable overrides and case names address the reported defects. |
| R2-F1–F21 | RESOLVED | Publication/reservation protection, re-depth acceptance, geometry and spacing checks, rename identity, numeric validation, limiter checks, failure reporting and source coverage address the original reproductions. |
| R3-F1–F13 | RESOLVED | Composed-field diagnostics, spacing controls, reusable bands, ladder-side selection, recovery, interpolation, serialization, overlap and timestep ordering are corrected. |
| R4-F2–F14 | RESOLVED | Interpolation, overlap tolerances, guide placement, recovery, input rechecks, boundary simplicity and runtime/path/array validation address the reported cases. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identity, continuous guide intersections, protected recovery, grid/path/query validation, early handlers and finite depths are corrected. |
| R6-F2–F11 | RESOLVED | Land windows, output handling, completion checks, failure exits, staging prechecks, recovery, inventories, native values, Fortran strings and relocation checks address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source triangulation, quoting/index checks, staged publication, work-area checks, verified CSV consumption, identity rechecks and publication ordering are corrected. |
| R8-F2–F13 | RESOLVED | Reported rollback, connectivity/type, sponge, identity, namelist, cache, dimension-order, serialization and manifest-selection defects are corrected. |
| R9-F2–F14 | RESOLVED | Original root collisions, recovery races, index/dtype defects, namelist/sponge handling, trimming, depth bounds, reading and library paths are corrected. |
| R10-F2–F14 | RESOLVED | Destination/finishing collisions, field validation, recovery tracking, comment parsing, QA normalization, serialization and documentation address the original cases. |
| R11-F2–F8 | RESOLVED | Finishing is isolated; malformed IDs/fields/connectivity fail; relocation restores absence; both acceptance paths check water holes; spaced paths work. Finding 2 concerns covered land, a different gap. |
| R12-F2–F7 | RESOLVED | Re-depth acceptance is hash-bound; benchmark destinations are unique; record dimensions and duration/advancement are checked; malformed fort.14 inputs fail; documentation follows owner policy. |
| R13-F2–F8 | RESOLVED | Staging locks, shorter benchmark paths, required dimensions/START_DATE, output advancement, acceptance hashes, land inventories and single-line titles address the reported cases. |
| R13-F9 | RESOLVED | The original invalid sizing controls and lattice cases are now rejected. Finding 4 concerns a newly accepted scalar type. |
| R14-F2 | RESOLVED | 411 takes the staging lock before inspection or writing. Finding 3 concerns invalidation before a refusal. |
| R14-F3 | RESOLVED | 414 invalidates old acceptance markers under the lock before replacement. Finding 3 concerns the no-replacement path. |
| R14-F4 | RESOLVED | Valid NumPy real scalars are accepted. The expanded type check introduces finding 4. |
| R14-F5 | RESOLVED | Arrays are normalized to matching 2-D shapes; nonfinite lattice coordinates and invalid size values are rejected. |
| R14-F6 | RESOLVED | 416 creates `$WORK_DIR/scratch` before `mktemp`. |
| R14-F7 | RESOLVED | Cadence uses microseconds; the synthetic 0.25-second history passes. |
| R14-F8 | PARTIAL | Explicit node selection avoids the backend probe, but `-k`/`-m` deselection happens too late; finding 5. |
| R14-F9 | RESOLVED | 401’s header describes its actual variants and submission command. |

1. **Major — existing F27 remains open: package code imports oceanmesh directly.**

   **Location:** [mesh_engine/oceanmesh.py:332](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_engine/oceanmesh.py:332); also `mesh_engine/multiscale.py:52`, `mesh_clean.py:1672`, and `autofinish/directives.py:20`.

   **Evidence:** `rg -n 'import oceanmesh|from oceanmesh' src/fvcom_mesh_tools` still identifies direct imports inside the package. This remains contrary to the subprocess-only contract specified in this review. It is the previously reported finding, not a new licensing conclusion.

   **Fix:** Complete the owner-selected resolution. Under the stated contract, move these calls behind subprocess interfaces and reconcile the package documentation accordingly. Concurrent relicensing edits are not treated here as an approved change to that contract.

2. **Major — the extension can cover real land while passing every acceptance gate.**

   **Location:** [446_extend_generate.py:303](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:303), [447_extend_merge.py:183](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:183), [extend.py:314](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:314); repeated in `453_redepth_extended.py:148`. The optional QA land check is at `qa.py:1133`.

   **Evidence:** CDT repair constrains the interface, open boundary and ladders, then keeps triangles according to their **centroids**. `check_island_holes` checks holes that exist; it cannot detect an island whose hole has disappeared. Neither extension acceptance path supplies `land_solid_shp` to QA, and that optional gate also checks centroids.

   A synthetic reproduction using the existing `_pristine()` QA fixture, with `SPACING=5000`, and land polygon `box(8400, 5100, 9900, 6100)` produced:

   ```text
   mesh:                         16 nodes, 18 elements
   mesh intersection with land:  1,500,000 m²
   nodes inside land:            0
   centroids inside land:        0
   check_island_holes:           {'n_new_islands': 0}
   ordinary QA:                  23/23 passed
   QA with optional land input:  24/24 passed
   ```

   Thus the actual acceptance checks certify a mesh covering the entire island. This demonstrates the acceptance gap; it does not establish that the supplied production mesh contains it.

   **Fix:** Before publication, check intersections between new element polygons and true land polygons, with an explicit approved coastline tolerance. Exclude the preserved base footprint and intentionally removed features explicitly. Preserve required shoreline topology during repair; adding the existing centroid gate alone is insufficient.

3. **Minor — rejected preparation deletes valid acceptance markers. Introduced by `ea04210`.**

   **Location:** [411_m2_prep.sh:28](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/411_m2_prep.sh:28), [414_refine_m2_prep.sh:26](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/414_refine_m2_prep.sh:26).

   **Reproduction:** Start with an accepted run root containing histories and `STAGED`, `SMOKE_OK` and case `RUN_OK` markers. Both scripts acquire the lock, remove the markers, then detect the existing histories and exit 2 without replacing the cases.

   I exercised the extracted shell bodies with mocked locking, `rm`, and `compgen`: both removed the markers before the reuse refusal. Consequently, `413_m2_analysis.sh:27` rejects the unchanged accepted run because `RUN_OK` has disappeared.

   **Fix:** Under the lock, perform the history/reuse and other no-write preflight checks first. Invalidate acceptance immediately before the first case replacement.

4. **Minor — complex NumPy gradation silently becomes a real gradation. Introduced by `ea04210`.**

   **Location:** [extend.py:56](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:56).

   **Reproduction:** `_grade` accepts every `np.number`, including `np.complex128`. With a 100-m lattice and values `[[1,100],[1,100]]`:

   ```python
   grade = np.complex128(0.2 + 7j)
   graded_up(values, x, y, grade)
   # [[80, 100], [80, 100]]

   compose_sizing(values, x, y, grade=grade)
   # [[1, 21], [1, 21]]
   ```

   Both succeed after `float(grade)` discards the imaginary component, emitting a `ComplexWarning`. The documented real-number validation is bypassed.

   **Fix:** Accept real scalar types explicitly, such as `numbers.Real` or NumPy integer/floating scalars, and reject complex scalars before conversion. Cover both public sizing functions.

5. **Minor — deselected backend tests still trigger the oceanmesh import. R14-F8 remains partial.**

   **Location:** [tests/conftest.py:29](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/tests/conftest.py:29).

   **Reproduction:** Run:

   ```text
   pytest -q -s -p no:cacheprovider \
     tests/test_mesh_engine_oceanmesh_safety.py \
     -k test_repair_flipped_elements_alias_is_public
   ```

   The repository collection hook examines marked items before pytest’s selection hooks remove them. A pure selected test therefore still causes `_importable("oceanmesh")`.

   In an actual in-memory pytest run, a `tryfirst` test plugin replaced `_importable` with a function raising `RuntimeError("UNSELECTED OCEANMESH PROBE")`. The invocation above reached that function and exited 3 with a collection `INTERNALERROR`, although the backend tests were excluded.

   **Fix:** Probe after selection hooks have completed, preferably using a post-yield collection-hook wrapper, and inspect only surviving items. Verify both `-k` selection and `-m 'not needs_oceanmesh'`.

6. **Minor — QA reports a single isolated element as having no isolated elements. Pre-existing.**

   **Location:** [qa.py:648](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:648).

   **Reproduction:** Run metric QA on one triangle with nodes `[[1000,1000],[2000,1000],[2000,2000]]`, connectivity `[[0,1,2]]`, and depths of 10 m. The `no_isolated_elements` check passes with `isolated = 0`, although the element has no neighbours.

   The `ne > 1` condition suppresses the offending element. Other gates reject this mesh, so this is an incorrect gate result and diagnostic rather than an overall acceptance bypass.

   **Fix:** Compute `np.where(deg == 0)[0]` for every nonempty mesh, including a single element.

7. **Minor — an unchecked recipe name can move the default output outside the repository. Pre-existing.**

   **Location:** [extend_recipe.py:66](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend_recipe.py:66), [445_extend_mesh.py:48](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/445_extend_mesh.py:48).

   **Reproduction:** In a mocked, otherwise valid recipe, set:

   ```yaml
   name: x/../../../escaped-review-output
   ```

   `load_extend_recipe` accepts it. With no explicit output argument, the driver’s default resolves to:

   ```text
   /octfs/work/G16445/v61021/Github/escaped-review-output
   ```

   This is outside the repository and its documented `outputs/` directory. The reservation protects existing nonempty directories, but does not prevent this unintended destination.

   **Fix:** Validate `name` as a single safe path component, or derive a validated slug separately. Assert that the automatically constructed destination remains beneath `REPO/outputs` before reservation.

Verification used the prescribed conda environment:

- Read-only pytest selection across `test_extend.py`, `test_extend_recipe.py`, `test_obc_design.py`, `test_dem_sources.py` and `test_check_run.py`: **52 passed**, using `pytest.main(["-q", "-s", "-p", "no:cacheprovider", *nodes])` with bytecode and plugin autoload disabled.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` on all **56** OCTOPUS shell scripts: **passed**.
- Diskless synthetic 0.25-second FVCOM history: `check_run` accepted all five records with no rejection reasons.
- Intentional faults produced the refusal and collection failures described above.
- No production mesh generation, FVCOM execution, shared-data reads, batch submission or full-suite rerun was performed. The production verification in your prompt remains supplied evidence.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 5 minor, 0 nit)

### Prompt

```markdown
# Review request, round 15: extending a base mesh outward (fvcom-mesh-tools)

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
Rounds 1-14 and their triage are in docs/extend-tools-review-20261001.md.
Round 14 (your previous answer; 9 findings) was fixed in ea04210; read
`git show ea04210`. Per finding:
- F1 (= F27): owner decision pending; report it as still open, not as new.
- F2 411: `fmesh_stage_lock` before inspecting/writing; old markers removed.
- F3 414: STAGED, SMOKE_OK and */RUN_OK removed under the lock before the
  cases are replaced (415's chain does not use STAGED; 421 publishes it).
- F4/F5 `extend._grade` and `extend._lattice`, used by `graded_up` and
  `compose_sizing`.
- F6 416 `mkdir -p $WORK_DIR/scratch`.
- F7 check_run cadence in microseconds.
- F8 conftest: oceanmesh probed only if a selected test is marked.
- F9 401 header.
Verification after ea04210: full test suite 1194 passed (batch
job 123456); on real data 444, the full 445 build (QA 23/23, status ok, 9
islands on land, grd and dep bit-identical to rounds 4-13) and 453 (status
ok) passed, and check_run accepts the two real FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 (F27) | major | as round 1 | yes | resolved by the owner's decision (2026-10-02): the package is relicensed to GPL-3.0-or-later, e37a433 |
| F2 | major | measured on the real build: 2,445 covered land pieces, all rocks and islets, at most 1.37 times the local element area | partly | fixed for land the mesh could resolve (check_land_cover, ratio 2), 5c394d8; rebutted for land smaller than the elements, which the owner's resolution principle drops (the reviewer's 1.5 km2 island under 5 km elements) |
| F3 | minor | code read (markers dropped before the reuse refusal) | yes, introduced in round 14 | fixed, 5c394d8 |
| F4 | minor | code read (np.complex128 is np.number) | yes, introduced in round 14 | fixed, 5c394d8; test |
| F5 | minor | code read (hook ran before deselection) | yes | fixed (hookwrapper), 5c394d8; -k and -m runs checked |
| F6 | minor | test (one triangle) | yes | fixed, 5c394d8; test |
| F7 | minor | test (name with ..) | yes | fixed, 5c394d8; test |

## Round 16

### Answer (verbatim)

Six findings remain: **2 major, 3 minor, 1 nit**. The repository was not modified.

Status refers to each original finding and reproduction. Ranges below include every finding from rounds 1–15.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R15-F1 | RESOLVED | `LICENSE`, `NOTICE`, project metadata and current policy permit oceanmesh imports under GPL-3.0-or-later. The licensing findings below concern different dependencies. |
| Generation-nondeterminism hypothesis | WITHDRAWN | Supplied repeated production measurements contradict it; no new supporting evidence found. |
| R1-F1–F26, F28–F29 | RESOLVED | Source handling, numeric validation, serialization, seam checks, reservations, provenance and smoke preparation address the original cases. F4 is resolved under the owner’s warning policy; mandatory timestep rejection is no longer required. |
| R2-F1–F21 | RESOLVED | Atomic reservation/publication, re-depth checks, interpolation, geometry checks and failure reporting address the original reproductions. |
| R3-F1–F13 | RESOLVED | Final-field diagnostics, resampling, reusable bands, ladder selection, publication recovery and serialization are corrected. F1 follows the approved warning policy. |
| R4-F2–F14 | RESOLVED | M2 interpolation, local overlap tolerance, guide placement, recovery, input identity, boundary simplicity and runtime validation are corrected. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identity, continuous intersections, recovery protection, early handlers and finite-depth checks address the reported cases. |
| R6-F2–F11 | RESOLVED | Land-window checks, output-directory handling, completion checks, failure exits, provenance inventories and path guards address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source triangulation, index/control validation, protected relocation, verified CSV consumption, identity rechecks and publication ordering are corrected. |
| R8-F2–F13 | RESOLVED | Rollback, connectivity/types, sponge handling, code identity, namelist parsing, caches, dimension order and manifest selection are corrected. |
| R9-F2–F14 | RESOLVED | Original concurrency, recovery, index/dtype, parsing, trimming, depth-bound, reader and library-path defects are corrected. |
| R10-F2–F14 | RESOLVED | Unique destinations, field checks, recovery tracking, parsing, QA normalization, serialization and documentation address the original cases. |
| R11-F2–F8 | RESOLVED | Finishing products are isolated; malformed arrays fail; relocation restores absence; water-hole checks and spaced paths are handled. |
| R12-F2–F7 | RESOLVED | Re-depth requires hash-bound acceptance; benchmarks are isolated; record dimensions, duration and advancement are checked; serialization and documentation are corrected. |
| R13-F2–F9 | RESOLVED | Locks, shorter paths, required dimensions/start time, output advancement, acceptance identity, inventories, titles and sizing validation address the original cases. |
| R14-F2–F9 | RESOLVED | Staging locks, marker invalidation, real NumPy scalars, lattice normalization, scratch creation, fractional cadence, post-selection backend probing and job documentation are corrected. |
| R15-F2 | PARTIAL | The small-island reproduction is withdrawn under the owner’s resolution principle. The new gate rejects the ordinary resolvable-land example, but its median calculation has finding 3 below. |
| R15-F3 | RESOLVED | Both preparation scripts check for existing histories before removing acceptance markers. |
| R15-F4 | RESOLVED | Complex NumPy gradation is rejected before conversion. |
| R15-F5 | RESOLVED | Both `-k` and `-m` fault probes pass without probing the deselected backend. |
| R15-F6 | RESOLVED | QA counts a single element as isolated. |
| R15-F7 | RESOLVED | Recipe `name` is validated; the default destination is checked beneath `outputs/`. |

1. **Major — permitted OCSMesh imports load GPL-incompatible Triangle.**

   **Location:** [mesh_compose/convert.py:28](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/mesh_compose/convert.py:28), [environment.yml:29](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/environment.yml:29), [THIRD_PARTY_NOTICES.md:50](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/THIRD_PARTY_NOTICES.md:50). The generation adapter also imports OCSMesh at `mesh_engine/ocsmesh.py:86`.

   **Evidence:** In the prescribed environment, converting a three-node mesh reaches this import chain:

   ```text
   fort14_to_meshdata
     → ocsmesh.internal
     → ocsmesh package initialization
     → ocsmesh.engines.factory
     → ocsmesh.engines.triangle
     → import triangle
   ```

   An in-memory import sentinel raised `FORBIDDEN_TRIANGLE_IMPORT` at that final import. No Triangle engine selection was necessary.

   The [Python wrapper](https://raw.githubusercontent.com/drufat/triangle/master/README.rst) wraps Shewchuk’s native Triangle library. That library restricts sale and inclusion in commercial products, contrary to the repository’s new GPL-compatible-import policy. Its wrapper’s LGPL metadata does not remove the native library’s restrictions. [Triangle author’s licensing statement](https://www.cs.cmu.edu/~quake/triangle.html).

   **Fix:** Remove Triangle from the default environment and require an OCSMesh installation whose permitted paths do not import it. Verify the transitive imports and correct the unconditional redistribution statement. Any retained Triangle functionality needs a separately handled external interface.

2. **Major — the environment-build job imports the explicitly excluded JIGSAW core.**

   **Location:** [build_env.sh:51](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/build_env.sh:51), [environment.yml:34](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/environment.yml:34), [THIRD_PARTY_NOTICES.md:31](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/THIRD_PARTY_NOTICES.md:31).

   **Evidence:** The job explicitly imports `jigsawpy`. Importing the installed module reaches:

   ```text
   jigsawpy/__init__.py:58 → import libsaw
   jigsawpy/libsaw.py:76   → ctypes.cdll.LoadLibrary("libjigsaw.so")
   ```

   A mocked native-library loader reproduced that call without loading the core. Thus this smoke import contradicts the stated “never imported” policy.

   The notice also incorrectly identifies the wrappers as LGPL-3.0. Both the installed Python header and the [upstream jigsaw-python license](https://raw.githubusercontent.com/dengwirda/jigsaw-python/master/LICENSE.md) carry restrictive JIGSAW terms.

   **Fix:** Remove `jigsawpy` from this import smoke test and the default environment, retaining it only as a separately handled optional installation if needed. Correct the wrapper-license description and document its actual handling.

3. **Minor — boundary-only contacts can make the land-cover gate accept resolvable covered land. Introduced by `5c394d8`.**

   **Location:** [extend.py:383](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:383), particularly line 389.

   **Evidence:** `predicate="intersects"` includes elements that touch a land piece only along an edge or at a vertex. Their full areas participate in the median despite covering none of the land.

   A conforming, nonoverlapping synthetic mesh with **10 nodes and 12 elements** produced:

   | Measurement | Result |
   |---|---:|
   | Triangular land area | 1,732,050.81 m² |
   | Land covered | 100% |
   | Elements covering positive land area | 3 |
   | Area of each covering element | 577,350.27 m² |
   | Land/local-element ratio | **3.00** |
   | Additional boundary-only contacts | 9 |
   | Reported ratio | **1.64948 — accepted** |

   The nine surrounding elements have areas of 952,627.94 or 1,147,483.66 m² and inflate the median. This violates the configured ratio of 2 even under the approved resolution principle.

   This is a gate-level reproduction: ordinary metric QA passed **16/17** gates; its remaining failure was numerical roundoff at an exact 30° angle. It does not establish a production-build acceptance bypass.

   **Fix:** Compute the median from elements with positive-area land overlap, excluding boundary-only contacts with an explicit numerical tolerance.

4. **Minor — invalid acceptance-check controls silently disable checks.**

   **Location:** [extend.py:375](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:375), line 391; also [extend.py:310](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:310).

   **Reproduction:** On the existing 6×6-node synthetic grid, a 9,000,000 m² land square normally fails with ratio **18.0**. However:

   ```python
   check_land_cover(mesh, land, 0, ratio=np.nan)
   # accepted; max_covered_area_ratio=18.0, ratio_limit=nan

   check_land_cover(mesh, land, mesh.n_elements + 1)
   # accepted; n_land_pieces_covered=0

   check_land_cover(mesh, land, -1)
   # accepted; n_land_pieces_covered=0
   ```

   Separately, an entirely duplicated outer triangle normally fails `check_no_overlap`, but `rel_tol=np.nan` accepts it.

   **Fix:** Before slicing or comparing, validate the base-element count as an integer in `[0, NE]`, the land ratio as finite and positive, and the overlap tolerance as finite and nonnegative. Invalid controls should raise. The land-check cases are introduced by `5c394d8`; the overlap case is pre-existing.

5. **Minor — frozen-base verification can certify an unmerged base by accepting an incomplete interface. Pre-existing.**

   **Location:** [extend.py:264](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:264).

   **Reproduction:** For the two-triangle square base with open boundary `[1, 2]`:

   ```python
   verify_frozen_base(base, base, [])
   # succeeds; n_interface_edges=-1

   verify_frozen_base(base, base, [1])
   # succeeds; n_interface_edges=0
   ```

   Its open-boundary edge remains exterior, contrary to the function’s stated contract. On a correctly merged example, passing `[1.9, 2.9]` also succeeds because the indices are silently truncated.

   **Fix:** Validate indices before casting and require the supplied interface to describe the base’s actual open-boundary edges, or derive those edges directly from `base.open_boundaries`. Preserve the legitimate empty-interface case only for bases without an open boundary.

6. **Nit — a current comment still states the superseded oceanmesh import prohibition.**

   **Location:** [446_extend_generate.py:12](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:12).

   **Evidence:** The comment says oceanmesh is imported in a subprocess and “never from the fvcom_mesh_tools package.” Current package code imports it, and the approved relicensing permits that.

   **Fix:** Describe this script’s retained subprocess organization and, if useful, its historical origin, following the updated comment in notebook 420.

Verification used the prescribed conda environment:

- Read-only pytest selection: **82 passed**, with bytecode, plugin autoload and cache writing disabled.
- Backend-selection fault probes: **1 passed/2 deselected** with `-k`; **2 passed/1 deselected** with `-m`.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n`: **56 scripts passed**.
- Synthetic geometry, invalid-control and import/native-load probes reproduced the findings above. Initial pytest capture failed because writable temporary files were unavailable; rerunning with `-s` passed. Dependency probes required cache/JIT mocks before reaching the intended sentinels.
- No production mesh generation, FVCOM execution, large shared-data reads or batch submissions occurred. Tracked files remained unchanged; the pre-existing untracked boundary PNG remained untouched.

## Verdict

VERDICT: FAIL (0 blocker, 2 major, 3 minor, 1 nit)

### Prompt

```markdown
# Review request, round 16: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-15 and their triage are in docs/extend-tools-review-20261001.md.

F27 (GPL imports) is RESOLVED by the owner's decision of 2026-10-02: the
package is relicensed from Apache-2.0 to GPL-3.0-or-later by its sole
copyright holder (commit e37a433: LICENSE is the GNU GPL v3 text; NOTICE,
pyproject `license`, THIRD_PARTY_NOTICES.md, CLAUDE.md's licence policy and
the docs updated). Importing oceanmesh (GPL-3.0-or-later) is therefore
allowed; the review contract's "subprocess-only" rule no longer applies.
Please check the relicensing for consistency (any place still claiming
Apache-2.0 or forbidding the import, any GPL-incompatible import).

Round 15 (your previous answer; 7 findings) was fixed in 5c394d8; read
`git show 5c394d8`. Per finding:
- F1: resolved as above.
- F2 new `extend.check_land_cover`, run in 447 before writing and in 453:
  a land piece (true land, base footprint removed) more than half under the
  new elements and larger than LAND_COVER_RATIO = 2 times the median area
  of those elements fails. REBUTTAL for the rest: land smaller than the
  elements over it is dropped on purpose (the owner's resolution principle:
  what the element size cannot carry goes); your 1.5 km2 island under 5 km
  elements is such a case. On the real build 2,445 covered pieces, all rocks
  and islets, have a ratio of at most 1.37. Coastal elements reaching onto
  land through the coast's approximation are not pieces covered by half.
- F3 411/414 remove the markers after the reuse checks, before writing.
- F4 `_grade`: int, float, np.integer, np.floating only.
- F5 conftest: hookwrapper(trylast); probes after deselection.
- F6 `no_isolated_elements` for any NE.
- F7 recipe `name` checked as one component; 445 default output must be
  under outputs/.
Verification after 5c394d8: full test suite 1198 passed (batch
job 123753); on real data 444, the full 445 build (QA 23/23, status ok, 9
islands on land, land cover: 2,445 covered pieces, maximum ratio 1.37,
grd bit-identical to rounds 4-14) and 453 (status ok, same land-cover
result) passed; check_run accepts the two real FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | major | code read (ocsmesh imports triangle at load) | yes | owner decision (2026-10-02): optional, private use; out of the default environment, documented, c76c0c6 |
| F2 | major | code read (build_env imports jigsawpy; wrapper licence misstated) | yes | same decision, c76c0c6 |
| F3 | minor | code read (edge contacts in the median) | yes, introduced in round 15 | fixed, 19f94e1; test |
| F4 | minor | tests (NaN ratio, out-of-range counts, NaN tolerance) | yes | fixed, 19f94e1; test |
| F5 | minor | test (empty and partial interfaces) | yes | fixed, 19f94e1; test |
| F6 | nit | doc read | yes | fixed, 19f94e1 |

We also checked a suspicion of our own: the autofinish directives call
`oceanmesh.remesh_patch(..., target_h=...)`; the function oceanmesh exports
is the DistMesh one in `mesh_merge.py`, which takes `target_h` and uses
neither JIGSAW nor Triangle, so the call is right.

## Round 17

### Answer (verbatim)

Found **8 minor defects and 1 nit**. No files were modified. No production mesh generation, FVCOM execution, shared-data reads or batch submissions occurred.

Status below refers to the original findings and reproductions; grouped ranges include every finding.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R15-F1 | RESOLVED | The approved GPL-3.0-or-later relicensing permits oceanmesh imports. |
| Generation-nondeterminism hypothesis | WITHDRAWN | The supplied repeated production measurements contradict it. |
| R1-F1–F26, F28–F29 | RESOLVED | Current source handling, validation, serialization, seam checks, reservations and smoke preparation address the original cases. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Reservation, publication, acceptance, interpolation, geometry and provenance fixes address the original cases. |
| R3-F1–F13 | RESOLVED | Final-field diagnostics, resampling, bands, ladders, recovery and serialization are corrected. F1 follows the approved warning policy. |
| R4-F2–F14 | RESOLVED | Interpolation, local overlap checks, guide placement, recovery, identity, boundary and runtime checks address the original cases. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte hashes, continuous intersections, recovery protection, early handlers and finite-depth checks are present. |
| R6-F2–F11 | RESOLVED | Window, output, completion, failure-exit, provenance, native-value and path checks address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source triangulation, validation, protected relocation, verified CSV consumption and identity checks are corrected. |
| R8-F2–F13 | RESOLVED | Rollback, types, sponge handling, identity, parsing, caches, axis order and manifest selection are corrected. |
| R9-F2–F14 | RESOLVED | Original concurrency, recovery, index, dtype, parsing, trimming, depth-bound and library-path defects are corrected. |
| R10-F2–F14 | RESOLVED | Destination isolation, field validation, recovery tracking, parsing, serialization and the originally reported documentation defects are corrected. |
| R11-F2–F8 | RESOLVED | Finishing products are isolated; malformed arrays fail; relocation, water-hole and spaced-path handling are corrected. |
| R12-F2–F7 | RESOLVED | Accepted-build hashes, isolated benchmarks, record dimensions, duration, advancement and serialization address the original cases. |
| R13-F2–F9 | RESOLVED | Locks, paths, history checks, acceptance identity, inventories, titles and sizing validation are corrected. |
| R14-F2–F9 | RESOLVED | Locks, marker ordering, scalar/lattice validation, scratch creation, fractional cadence, backend selection and the reported job documentation are corrected. |
| R15-F2 | RESOLVED | The gate now excludes boundary-only contacts from its median. The small-island objection remains withdrawn under the owner’s resolution principle. |
| R15-F3–F7 | RESOLVED | Marker preservation, complex-scalar rejection, post-selection probing, isolated-element detection and recipe-name validation are corrected. |
| R16-F1–F2 | PARTIAL | Default dependencies and smoke imports are corrected; optional private-use handling is documented. Contradictory statements remain in the notices—finding 1. |
| R16-F3 | RESOLVED | The median uses elements with positive land overlap, subject to numerical tolerance. |
| R16-F4 | RESOLVED | The originally reported NaN controls and out-of-range element counts are rejected. |
| R16-F5 | RESOLVED | Whole indices and equality with the base’s open-boundary edge set prevent the original incomplete-interface bypass. |
| R16-F6 | RESOLVED | Notebook 446’s comment reflects the approved relicensing. |

No new runtime defect was substantiated in the `19f94e1` changes. The remaining notices inconsistency makes `c76c0c6` incomplete.

1. **Minor — Third-party notices contradict the approved private-backend policy.**

   **Location:** [THIRD_PARTY_NOTICES.md:62](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/THIRD_PARTY_NOTICES.md:62), also line 14.

   **Evidence:** The redistribution summary still includes OCSMesh among components that may ship together with this package, requiring only GPL compliance and attribution. Lines 29–31 expressly forbid that arrangement. The introduction also says the JIGSAW core is “never imported,” despite the newly permitted private-use paths.

   **Fix:** Remove OCSMesh from the permitted-shipping row, explicitly exclude the optional private backends from redistribution, and qualify the import prohibition with the approved private-use exception. This finding concerns documentary consistency, not the existence of those backends.

2. **Minor — One job still takes its scientific input path from `$HOME`.**

   **Location:** [422_base_rfactor.sh:15](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/422_base_rfactor.sh:15).

   **Evidence:** After sourcing `common.sh`, it sets:

   ```bash
   G=$HOME/Github/TB-FVCOM/input/goto2023/grid
   ```

   Changing `$WORK_DIR` does not change this selection. With distinct home and work directories, the job either fails to find the intended repository or reads a different home-directory checkout. This contradicts the portability requirement covering every job script.

   **Fix:** Use `G="$WORK_DIR/Github/TB-FVCOM/input/goto2023/grid"`.

3. **Minor — Boundary publication can overwrite its own design recipe.**

   **Location:** [444_design_obc.py:38](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:38), publication at line 281.

   **Reproduction/evidence:** Passing the same path for the design and output produces `design_path == out_csv`; no guard rejects it. An in-memory execution of the actual publication loop replaced the original YAML contents with the accepted CSV. Successful cleanup then removes the backup. A sidecar destination can likewise alias the recipe—for example, a JSON-formatted design at `boundary.json` with output `boundary.csv`.

   **Fix:** Before processing, reject any resolved CSV, report or figure destination that aliases the design input, including existing hard-link aliases. Require distinct product destinations.

4. **Minor — Nonfinite merge tolerance silently bypasses interface matching.**

   **Location:** [extend.py:223](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:223).

   **Reproduction:** Using the two-square synthetic base/outer example, displacing both outer interface nodes by `(5000, 5000)` gives a **7,071.068 m** mismatch. `merge_outer(..., tol_m=np.nan)` accepts it, substitutes the base coordinates, and the resulting mesh passes `verify_frozen_base`. Infinite tolerance also accepts it.

   **Fix:** Validate `tol_m` as a finite, nonnegative real before comparing or merging. Add NaN, infinity and negative-tolerance regressions.

5. **Minor — Invalid QA thresholds can produce a passing report and CLI exit zero.**

   **Location:** [qa.py:1039](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:1039), with the same problem in other threshold comparisons.

   **Reproduction:** A synthetic **7-node, 6-element** hexagonal fan with depths of 10 m fails `min_depth_m=11`, but passes `min_depth_m=np.nan`. With mesh reads and report writes mocked in memory:

   ```text
   fmesh-mesh-qa … --min-depth 11  → exit 1
   fmesh-mesh-qa … --min-depth nan → exit 0
   ```

   Similarly, `max_valence=5` fails while `max_valence=np.nan` passes.

   **Fix:** Validate all QA thresholds centrally for finite values, appropriate ranges and integral counts before evaluating the mesh. This concerns invalid explicitly supplied controls, not the approved extension timestep-warning policy.

6. **Minor — A missing explicitly requested land dataset silently removes the land gate.**

   **Location:** [qa.py:1134](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:1134).

   **Reproduction:** On the same synthetic fan, a mocked solid-land polygon covering the mesh produces **6 land violations** and one failed gate. Making only the requested shapefile’s existence check return false produces **zero failed gates**, with `land_overlap` absent from the report.

   **Fix:** When `land_solid_shp` is supplied, require a readable dataset and fail explicitly if it is missing. Reserve omission of the gate for `None`. This is an API reproduction; it does not demonstrate a default 445 acceptance bypass.

7. **Minor — QA crashes on nonfinite coordinates instead of reporting invalid geometry.**

   **Location:** [qa.py:625](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:625).

   **Reproduction:** The valid synthetic fan passes QA. Replacing its centre coordinates with NaN causes:

   ```text
   GEOSException:
   IllegalArgumentException: Points of LinearRing do not form a closed linestring
   ```

   Coordinates reach geometry calculations without the validation that indices and depths receive.

   **Fix:** Add an early coordinate shape/finiteness gate and return a failed QA report before projection or geometry construction. Native readers already reject this input; direct `run_qa` callers remain affected.

8. **Minor — Frozen-base verification reports a negative interface-edge count.**

   **Location:** [extend.py:302](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:302).

   **Reproduction:** For a valid closed two-triangle base with no open boundary:

   ```python
   verify_frozen_base(base, base, [])
   # n_interface_edges: -1
   ```

   The legitimate empty-interface case succeeds, but its diagnostic count is wrong.

   **Fix:** Report the number of verified interface edges, such as `len(want)`, rather than `len(ib) - 1`.

9. **Nit — The user guide still instructs users to edit machine-specific job paths.**

   **Location:** [USER_GUIDE.md:125](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/docs/USER_GUIDE.md:125).

   **Evidence:** The OCTOPUS section says job paths belong to this account, lists hard-coded FVCOM/scratch locations, and tells another user to edit them. Current scripts derive these paths from the environment.

   **Fix:** Document `$DATA_DIR`, `$WORK_DIR`, `FMESH_FVCOM` and `FVCOM_LIBS`; retain the accounting-group customization separately.

Verification:

- Read-only pytest selection across extension, OBC, bathymetry, native I/O, QA and provenance tests: **81 passed**.
- Initial scoped pytest run: **51 passed, 16 setup errors**, all caused by unavailable writable temporary directories.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n`: **56 scripts passed**.
- Synthetic probes and in-memory fault injections reproduced findings 3–8.
- `git diff --exit-code`: clean; the pre-existing untracked boundary PNG remained untouched.

## Verdict

VERDICT: FAIL (0 blocker, 0 major, 8 minor, 1 nit)

### Prompt

```markdown
# Review request, round 17: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-16 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later since e37a433 (owner's decision); F27 is
resolved (you confirmed this in round 16).

Round 16 (your previous answer; 6 findings) was fixed in 19f94e1 and
c76c0c6; read both. Per finding:
- F1/F2 (Triangle via OCSMesh, JIGSAW via jigsawpy): owner's decision
  (2026-10-02): these stay OPTIONAL, PRIVATE-USE backends. environment.yml
  no longer installs triangle, ocsmesh or jigsawpy (it says how to install
  them privately); build_env.sh no longer imports them; the pyproject
  extras say private use only and "all" holds only GPL-compatible extras;
  THIRD_PARTY_NOTICES.md has a section for them; CLAUDE.md states the rule
  (no new code may depend on them). The OCSMesh paths (--engine ocsmesh,
  fmesh-mesh-combine overlap/neighbor, --repair-skewed-elements) import
  OCSMesh lazily and remain. Please do not re-report their existence; do
  report any place where the default paths still reach Triangle or JIGSAW,
  or where the documents say otherwise.
- F3 `check_land_cover`: median over elements holding positive land area.
- F4 control validation in `check_land_cover` and `check_no_overlap`.
- F5 `verify_frozen_base`: interface edges must equal the base's open
  boundary edges; indices checked whole.
- F6 446 comment.
Our own check: `oceanmesh.remesh_patch` (exported) is the DistMesh one in
`mesh_merge.py` (takes target_h; no JIGSAW/Triangle), so the autofinish
directives call is correct.
Verification after 19f94e1: full test suite 1200 passed (batch job
123886); real-data 444/445/453 passed as before (QA 23/23, land cover max
ratio 1.37, grd bit-identical). After c76c0c6 the needs_ocsmesh tests run
when OCSMesh is installed and are skipped otherwise.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | minor | doc read | yes | fixed, 87f3a06 |
| F2 | minor | code read ($HOME in 422) | yes | fixed, 87f3a06 |
| F3 | minor | run 444 with the design as output (refused after the fix) | yes | fixed, 87f3a06 |
| F4 | minor | test (NaN, inf, negative tol_m) | yes | fixed, 87f3a06; test |
| F5 | minor | test (NaN thresholds) | yes | fixed, 87f3a06; test |
| F6 | minor | test (missing land file) | yes | fixed, 87f3a06; test |
| F7 | minor | test (NaN node) | yes | fixed, 87f3a06; test |
| F8 | minor | test (closed base) | yes | fixed, 87f3a06; test |
| F9 | nit | doc read | yes | fixed, 87f3a06 |

## Round 18

### Answer (verbatim)

Found **1 major and 5 minor defects**. No files were modified. Production mesh generation and FVCOM execution were not run; no shared-data reads or batch submissions occurred.

Previous-finding statuses below cover every numbered finding; ranges refer to the original reproductions.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R15-F1 | RESOLVED | The approved GPL-3.0-or-later relicensing permits oceanmesh imports. |
| Generation-nondeterminism hypothesis | WITHDRAWN | The supplied repeated production measurements contradict it. |
| R1-F1–F26, F28–F29 | RESOLVED | Source handling, validation, serialization, seam checks, reservations and smoke preparation address the original cases. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Exclusive reservations, publication recovery, acceptance checks, interpolation, geometry and provenance fixes address the original cases. |
| R3-F1–F13 | RESOLVED | Final-field diagnostics, resampling, bands, ladders, recovery and serialization are corrected. F1 follows the approved warning policy. |
| R4-F2–F14 | RESOLVED | Harmonic coefficients are interpolated correctly; local overlap, guide placement, publication, identity, boundary and runtime checks address the original cases. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte hashes, continuous intersections, recovery protection, early handlers and finite-depth checks are present. |
| R6-F2–F11 | RESOLVED | Window, output, completion, failure-exit, provenance, native-value and path checks address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source triangulation, validation, protected relocation, verified CSV consumption and identity checks are corrected. |
| R8-F2–F13 | RESOLVED | Rollback, types, sponge handling, identity, parsing, caches, axis order and manifest selection are corrected. |
| R9-F2–F14 | RESOLVED | Original concurrency, recovery, index, dtype, parsing, trimming, depth-bound and library-path defects are corrected. |
| R10-F2–F14 | RESOLVED | Destination isolation, field validation, recovery tracking, parsing, serialization and documentation fixes address the original cases. |
| R11-F2–F8 | RESOLVED | Finishing products are isolated; malformed indices fail; relocation, unsupported water holes and spaced paths are handled. |
| R12-F2–F7 | RESOLVED | Accepted-build hashes, isolated benchmarks, record dimensions, duration, advancement and serialization address the original cases. |
| R13-F2–F9 | RESOLVED | Locks, paths, history checks, acceptance identity, inventories, titles and sizing validation are corrected. |
| R14-F2–F9 | RESOLVED | Locks, marker ordering, scalar/lattice validation, scratch creation, fractional cadence, backend selection and documentation are corrected. |
| R15-F2 | RESOLVED | Resolvable covered land is checked. The small-island objection remains withdrawn under the owner’s resolution principle. |
| R15-F3–F7 | RESOLVED | Marker preservation, complex-scalar rejection, backend probing, isolated-element detection and recipe-name validation are corrected. |
| R16-F1–F2 | RESOLVED | The approved private-use policy is reflected in the default environment, imports and corrected notices. |
| R16-F3–F6 | RESOLVED | Positive-area contacts, control validation, complete-interface verification and the relicensing comment address the original cases. |
| R17-F1 | RESOLVED | Notices now consistently exclude the private backends from redistribution. |
| R17-F2 | RESOLVED | Job 422 obtains the grid from `$WORK_DIR`. |
| R17-F3 | PARTIAL | Success products cannot overwrite the design, but rejection reports remain unprotected and product destinations need not be distinct—findings 2 and 3. |
| R17-F4 | RESOLVED | Nonfinite and negative merge tolerances are rejected. |
| R17-F5 | RESOLVED | The original invalid QA thresholds are rejected centrally. |
| R17-F6 | RESOLVED | A named missing land dataset raises explicitly. |
| R17-F7 | RESOLVED | The original NaN-coordinate reproduction returns failed QA before geometry. Finding 4 is a separate regression. |
| R17-F8 | RESOLVED | The diagnostic counts verified interface edges, including zero for a closed base. |
| R17-F9 | RESOLVED | USER_GUIDE §2.4 documents environment-derived paths. |

1. **Major — The M2 smoke test can retain a boundary condition that ignores its forcing.**

   **Location:** [448_extend_smoke.py:93](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/448_extend_smoke.py:93).

   **Evidence:** The export inherits `mesh.obc_type`. An in-memory execution of the staging code using a **16-node, 18-element** synthetic case with type 3 wrote:

   ```text
   OBC Node Number = 2
   1 5 3
   2 9 3
   ```

   It also wrote M2 forcing and a manifest claiming “uniform M2 on the open boundary.” Type 3 is the zero-elevation clamp documented in `docs/fvcom_source_constraints.md`; it does not apply that tide. The completion checker requires advancing, finite output, so a quiescent, unforced integration can receive `SMOKE_OK`. This is a realistic path for a finished case using a supported non-tidal OBC type.

   **Fix:** Explicitly export the smoke boundary with the intended tidal type, normally `obc_type=1`, and record and verify that type in the manifest and staged file.

2. **Minor — A rejected design can still overwrite its own recipe. Residual R17-F3.**

   **Location:** [444_design_obc.py:204](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:204).

   **Reproduction:** Use a JSON-formatted design at `boundary.rejected.json` and request `boundary.csv`. JSON is valid YAML. The new guard accepts these paths because it checks only the CSV, `.json` and `.png` success products. Executing the actual rejection branch against an in-memory filesystem replaced the design with:

   ```json
   {"problems": ["injected rejection"]}
   ```

   An existing rejection-report hard link or symlink to a normal YAML design has the same risk; `write_text` follows it.

   **Fix:** Include the rejection destination in the design-alias checks, including existing file identities, before processing or writing any report.

3. **Minor — Colliding success-product paths produce a successful but unusable publication. Residual R17-F3.**

   **Location:** [444_design_obc.py:41](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:41), publication at [line 286](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:286).

   **Reproduction:** Request output `boundary.json`. Its CSV destination and report destination are then identical. Executing the actual publication loop with mocked files first published the JSON report, then replaced it with CSV, without an exception. Loading this boundary subsequently treats that same file as its JSON sidecar and raises `JSONDecodeError`. An output ending in `.png` similarly collides with the figure.

   **Fix:** Require distinct resolved product destinations before staging. Requiring a `.csv` output suffix would also prevent these direct collisions.

4. **Minor — The new coordinate check can prevent invalid-index QA from returning a report. Introduced by `87f3a06`.**

   **Location:** [qa.py:652](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:652).

   **Reproduction:** In the synthetic pristine mesh, set one connectivity reference to `999`, make coordinates an object array, and set one coordinate to `"bad"`.

   ```text
   Before 87f3a06: returned QA report with one failed gate
   Current: ValueError: could not convert string to float: 'bad'
   ```

   The unchecked conversion now occurs before the existing invalid-index return. A malformed coordinate array therefore suppresses the diagnostic report.

   **Fix:** Catch conversion failures and report invalid coordinates through the early integrity gate. Validate coordinate shape and real numeric dtype before conversion.

5. **Minor — QA accepts malformed depth arrays and emits invalid timestep diagnostics. Pre-existing.**

   **Location:** [qa.py:1084](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:1084), timestep calculation at [line 1163](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:1163).

   **Reproduction:** On the **16-node, 18-element** pristine mesh, `run_qa(..., coords="metric", channel_check=False)`:

   - accepts 17 finite depths and reports `passed=True`;
   - accepts depths shaped `(16, 1)` and reports `passed=True`;
   - crashes with `IndexError` for 15 depths.

   With depths `np.arange(5., 21.)`, the vector reports worst timestep element **16**. Reshaping the same values to `(16, 1)` broadcasts the timestep calculation into an element-by-element matrix and reports worst element **288**, although there are only 18 elements.

   **Fix:** Require depths shaped exactly `(NP,)` before geometry or diagnostics, and return failed integrity QA for a mismatch.

6. **Minor — Direct use of the merge stage can overwrite the adopted base case. Pre-existing.**

   **Location:** [447_extend_merge.py:59](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:59), export at [line 193](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:193).

   **Evidence:** The script advertises a standalone invocation but does not protect its output from input aliases or require the driver’s reservation. A valid recipe may use `case == base_case`; supplying the base directory as `OUTDIR` then replaces the source files.

   Executing the actual export statement and native exporter against an in-memory filesystem changed the base grid from **4 nodes / 2 elements** to **6 nodes / 4 elements**, removed the original backup copies, and still passed `verify_frozen_base` against the previously loaded base object. Later QA failure cannot restore those inputs. The normal 445 driver prevents this through its reservation; direct 447 invocation remains exposed.

   **Fix:** Reject output products that alias any input. Require authenticated ownership of the driver’s reserved output, or reserve a fresh destination for standalone execution.

Finding 4 is the substantiated regression introduced by the latest fixes. Findings 2 and 3 are incomplete aspects of the design-publication fix; findings 1, 5 and 6 are pre-existing.

Verification:

- Read-only `pytest.main` selection across extension, recipe, OBC design, bathymetry, QA and provenance tests, with `--capture=sys -p no:cacheprovider`: **89 passed, 57 deselected**. Deselected tests required writable temporary directories.
- The initial pytest attempt failed before collection because default file-descriptor capture required a writable temporary directory.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n`: **56 scripts passed**.
- Small synthetic arrays and in-memory publication/export probes reproduced the findings above.
- `git diff --exit-code`: **clean**. The pre-existing untracked boundary PNG remained untouched.

## Verdict

VERDICT: FAIL (0 blocker, 1 major, 5 minor, 0 nit)

### Prompt

```markdown
# Review request, round 18: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-17 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 17 (your previous answer; 9 findings, no major) was fixed in
87f3a06; read it. Per finding:
- F1 THIRD_PARTY_NOTICES introduction and redistribution summary.
- F2 422 uses $WORK_DIR.
- F3 444 refuses products that are the design (path or hard link).
- F4 `merge_outer` validates tol_m.
- F5 `run_qa` validates thresholds (`_check_qa_controls`).
- F6 `run_qa` raises on a named but missing land_solid_shp.
- F7 non-finite coordinates fail `node_index_valid` (count unchanged) and
  end QA before geometry.
- F8 `n_interface_edges` = edges verified.
- F9 USER_GUIDE §2.4.
Verification after 87f3a06: full test suite 1203 passed (batch
job 124167); on real data 444, 445 (QA 23/23, status ok, land cover max
ratio 1.37, grd bit-identical to rounds 4-16) and 453 (status ok) passed;
check_run accepts the two real FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | major | code read (export inherits the case's OBC type) | yes | fixed (type 1, read back, manifest), 371de18 |
| F2 | minor | run 444 with the design as the rejection report (refused after the fix) | yes | fixed, 371de18 |
| F3 | minor | run 444 with a .json output (refused after the fix) | yes | fixed, 371de18 |
| F4 | minor | test (object coordinates) | yes, introduced in round 17 | fixed, 371de18; test |
| F5 | minor | tests (17, 15 and (16, 1) depths) | yes | fixed, 371de18; tests |
| F6 | minor | code read | yes | fixed, 371de18 |

## Round 19

### Answer (verbatim)

Found **3 major and 4 minor defects**. No files were modified. Production mesh generation and FVCOM execution were not run; no shared-data reads or batch submissions occurred.

Previous-finding statuses cover every numbered finding. Ranges refer to the original reproductions.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R15-F1 | RESOLVED | The approved GPL relicensing permits oceanmesh imports. |
| Generation-nondeterminism hypothesis | WITHDRAWN | The supplied repeated production measurements contradict it. |
| R1-F1–F26, F28–F29 | RESOLVED | Original source, geometry, serialization, validation, reservation and smoke reproductions are addressed. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Reservations, publication recovery, interpolation, seam checks and provenance address the original cases. |
| R3-F1–F13 | RESOLVED | Final-field reporting, chord spacing, band iteration, ladders, serialization and timestep ordering are corrected. F1 follows the warning policy. |
| R4-F2–F14 | RESOLVED | Harmonic coefficients, interpolation shape, local overlap, publication recovery, input identity and runtime checks are corrected. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte hashes, continuous intersections, recovery protection and finite-depth checks are present. |
| R6-F2–F11 | RESOLVED | Window, output, completion, failure-exit, provenance, native-value and path checks address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source triangulation, export validation, protected relocation, verified CSV consumption and completion identity checks are corrected. |
| R8-F2–F13 | RESOLVED | Rollback, type validation, sponge handling, parsing, caches, axis order and manifest selection are corrected. |
| R9-F2–F14 | RESOLVED | Original concurrency, recovery, index, dtype, trimming, depth-bound and library-path defects are corrected. |
| R10-F2–F14 | RESOLVED | Destination isolation, field validation, recovery tracking, parsing, serialization and documentation address the original cases. |
| R11-F2–F8 | RESOLVED | Finishing products are isolated; malformed indices fail; relocation, unsupported water holes and spaced paths are handled. |
| R12-F2–F7 | RESOLVED | Accepted-build hashes, isolated benchmarks, record dimensions, duration, advancement and serialization address the original cases. |
| R13-F2–F9 | RESOLVED | Locks, paths, history checks, acceptance identity, inventories, titles and sizing validation are corrected. |
| R14-F2–F9 | RESOLVED | Locks, markers, lattice validation, scratch creation, fractional cadence, backend selection and documentation are corrected. |
| R15-F2 | PARTIAL | Resolvable covered islands are checked; resolvable mainland portions can still escape—finding 1. |
| R15 small-island objection | WITHDRAWN | The owner’s resolution principle permits dropping unresolvable islets. Finding 1 concerns a resolvable peninsula. |
| R15-F3–F7 | RESOLVED | Marker preservation, scalar validation, backend probing, isolated-element detection and recipe-name validation are corrected. |
| R16-F1–F2 | RESOLVED | The approved private-use policy is reflected in the environment and notices. |
| R16-F3–F6 | RESOLVED | Positive-area contacts, control validation, complete-interface verification and the relicensing comment address the original cases. |
| R17-F1–F9 | RESOLVED | Notices, environment paths, design protection, QA controls, missing-land handling, coordinate checks and interface counts are corrected. |
| R18-F1 | RESOLVED | 448 explicitly exports type 1, reads it back and records both boundary types. |
| R18-F2 | RESOLVED | The rejection destination is checked against the design, including file identity. |
| R18-F3 | RESOLVED | Requiring `.csv` prevents the original success-product collisions. |
| R18-F4 | RESOLVED | The object-coordinate reproduction returns failed integrity QA without conversion failure. |
| R18-F5 | RESOLVED | Depths must have numeric dtype and exactly `(NP,)` shape. |
| R18-F6 | PARTIAL | Base-directory aliases are refused, but an existing reservation bypasses output ownership—finding 2. |

No additional behavioral regression unique to `371de18` was substantiated. Its reservation protection remains incomplete; the other findings are pre-existing.

1. **Major — Land-cover acceptance ignores resolvable portions of mainland.**

   **Location:** [extend.py:427](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:427).

   **Reproduction:** On the pristine **16-node, 18-element** grid, each element has area **0.5 km²**. Supply a peninsula head `box(1000, 750, 2750, 2500)`, area **3.0625 km²**. Supplied alone, the gate rejects it as **6.1 times** the element area.

   Connect that same head through a narrow corridor to a large mainland polygon outside the mesh. The mesh now covers **3.1625 km²** of land, but the gate returns:

   ```text
   n_land_pieces_covered = 0
   max_covered_area_ratio = 0.0
   QA passed = True
   ```

   The `covered <= 0.5 * piece.area` condition measures the entire connected mainland polygon. Adding land outside the mesh therefore disables detection of the unchanged, resolvable covered head. Both 447 and 453 use this gate.

   **Fix:** Add a local mainland-overlap check based on covered patches, local element size and permitted coastline approximation. Do not use the entire connected mainland’s area to exempt an interior covered peninsula. Test that attaching exterior mainland does not change acceptance.

2. **Major — Standalone extension stages do not enforce output ownership.**

   **Location:** [447_extend_merge.py:69](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:69); related path at [446_extend_generate.py:53](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:53).

   **Evidence:** `.reserved` remains after successful builds. Executing the actual new guard against an earlier output containing that marker bypassed `reserve()` entirely. It checks neither ownership nor completion.

   An in-memory native-export reproduction then replaced the existing case’s depths from **5 m to 6 m**, retaining the old reservation marker. A standalone rerun—or another merge while 445 owns the directory—can replace case files and reports. The normal 445 reservation does not protect against this bypass. Standalone 446 likewise accepts a populated destination and overwrites generation products.

   **Fix:** Pass and validate a reservation token between 445 and its stages. Standalone invocations must reserve fresh destinations. Refuse completed outputs and acquire an exclusive merge-stage lock to prevent simultaneous merges.

3. **Major — The refinement M2 experiment can still ignore its tidal forcing.**

   **Location:** [414_refine_m2_prep.py:210](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/414_refine_m2_prep.py:210).

   **Evidence:** Unlike corrected 448, 414 exports each input’s inherited boundary type. A supported type-3 synthetic case produces:

   ```text
   OBC Node Number = 2
   1 5 3
   2 9 3
   ```

   The subsequent code writes M2 elevation forcing. Type 3 clamps elevation to zero, so that forcing is ignored. Complete, finite output can still be analyzed; an unforced zero response also satisfies the half-window convergence criterion. Different inherited types additionally undermine the claim that the two experiments differ only in the refinement.

   **Fix:** Stage both M2 cases with the same explicitly selected tidal boundary type, normally 1. Read back and record the staged types, as 448 now does.

4. **Minor — QA uses extra coordinate columns to approve invalid planar geometry.**

   **Location:** [qa.py:264](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:264); acceptance at [qa.py:656](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:656).

   **Reproduction:**

   ```python
   mesh = _pristine()
   mesh.nodes[:, 0] *= 0.1
   mesh.nodes = np.c_[mesh.nodes, 10 * mesh.nodes[:, 0]]
   run_qa(mesh, coords="metric", channel_check=False)
   ```

   Before adding the third column, all **18 elements** fail the minimum-angle gate at **5.71°**, and the reported timestep is **14.2784 s**. Afterward, QA reports **passed=True**, with timestep **142.7843 s**.

   The integrity check permits `(NP, >=2)`, but `_metric_nodes` retains every column. Edge lengths and angles use three dimensions, while signed areas and overlap use x/y.

   **Fix:** Normalize geometry to `nodes[:, :2]` throughout QA, or require exactly two coordinate columns. Test that an additional column cannot change planar QA or timestep results.

5. **Minor — A rejection report can corrupt the previous boundary through an alias.**

   **Location:** [444_design_obc.py:212](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:212).

   **Reproduction:** Make `boundary.rejected.json` a hard link or symlink to the existing `boundary.csv`, with a separate design recipe. The product guards accept this because they compare products only against the design.

   Executing the actual rejection branch with an in-memory hard-link model replaced the previous CSV with:

   ```json
   {"problems": ["injected rejection"]}
   ```

   It then reported that `boundary.csv` was “left as it was.”

   **Fix:** Publish rejection reports using a unique temporary file and `os.replace`, and reject rejection destinations that alias success products.

6. **Minor — The native grid reader accepts invalid record identifiers.**

   **Location:** [fvcom_native.py:553](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:553).

   **Reproduction:** Replace every node-record identifier in the pristine native grid with `999`, leaving coordinates, connectivity and depth rows unchanged. `read_fvcom_case` accepts it as the same **16-node, 18-element** mesh because `read_grd` discards the identifiers.

   This is not a valid FVCOM grid: the production reader uses the first node identifier `1` to detect the transition from connectivity to node records. The fort.14 reader already rejects comparable identifier corruption.

   **Fix:** Validate cell identifiers as `1..NE` and node identifiers as `1..NP` in order before returning the mesh. Preserve support for the documented trailing columns.

7. **Minor — Smoke completion accepts history from a different mesh with matching counts.**

   **Location:** [check_run.py:226](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:226).

   **Reproduction:** An in-memory NetCDF history had valid times, `TADA`, finite `zeta/ua/va`, and the staged **3-node, 1-element** counts. The staged coordinates were near `(0, 0)` with **5 m** depths; history coordinates were translated near `(100000, 90000)` with **100 m** depths.

   ```text
   ok = True
   reasons = []
   ```

   `_grid_counts` binds history only to header counts. Fresh smoke directories prevent the original ordinary stale-output case, but injected or misplaced foreign histories remain indistinguishable.

   **Fix:** Compare history connectivity, coordinates and bathymetry against the staged inputs, allowing FVCOM’s orientation reversal and documented numeric tolerances. Require these identity checks before issuing success.

Verification:

- Read-only `pytest.main` selection across extension, recipes, OBC design, DEM sources, QA and provenance: **93 passed, 57 deselected**. Deselected tests required writable fixtures.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n`: **all 56 OCTOPUS scripts passed**.
- Small synthetic arrays, in-memory file mocks and an in-memory NetCDF reproduced the findings. Initial history probes omitted grid-read mocks and failed closed; completing the mocks reproduced finding 7.
- `git diff --exit-code`: **clean**. The pre-existing untracked boundary PNG remained untouched.

## Verdict

VERDICT: FAIL (0 blocker, 3 major, 4 minor, 0 nit)

### Prompt

```markdown
# Review request, round 19: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-18 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 18 (your previous answer; 6 findings) was fixed in 371de18; read it.
Per finding:
- F1 448: `export_fvcom_case(..., obc_type=1)`, written _obc.dat read back
  (types must all be 1), manifest records obc_type and case_obc_type.
- F2/F3 444: output must end in .csv; CSV, .json, .png and .rejected.json
  checked against the design by path, resolved path and samefile.
- F4/F5 `run_qa`: coordinates (NP, >=2) numeric checked without
  conversion; depths exactly (NP,) numeric; failures go to node_index_valid
  and end QA.
- F6 447: OUTDIR may not be or contain the base directory; without the
  reservation marker it reserves OUTDIR itself.
Verification after 371de18: full test suite 1207 passed (batch
job 124205); on real data 444, 445 (QA 23/23, status ok, grd bit-identical
to rounds 4-17) and 453 (status ok) passed; check_run accepts the two real
FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | major | test (head joined to a mainland) | yes | fixed (eroded covered patches; measured on the real build first), 94364a4; test |
| F2 | major | code read (marker accepted without ownership) | yes | fixed (reservation token, outdir.claim), 94364a4; test |
| F3 | major | code read (414 inherits the OBC type) | yes | fixed, 94364a4 |
| F4 | minor | test (extra column) | yes | fixed, 94364a4; test |
| F5 | minor | code read (write_text through a link) | yes | fixed, 94364a4 |
| F6 | minor | test (node id 999) | yes | fixed, 94364a4; test |
| F7 | minor | test (translated history) | yes | fixed, 94364a4; test; a real history passes |

## Round 20

### Answer (verbatim)

Found **3 major defects, 5 minor defects and 1 nit**. No files were modified.

Previous-finding statuses refer to the original reproductions. The inclusive ranges below cover every numbered finding from rounds 1–19.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R15-F1 | RESOLVED | The approved GPL relicensing permits oceanmesh imports. |
| R1-F1–F26, F28–F29 | RESOLVED | Source handling, validation, serialization, seam checks, provenance, reservations and smoke preparation address the original cases. F4 follows the approved timestep-warning policy. |
| R2-F1–F21 | RESOLVED | Reservation, publication recovery, re-depth acceptance, interpolation, geometry and provenance fixes address the reported cases. |
| R3-F1–F13 | RESOLVED | Final-field diagnostics, spacing, band iteration, ladders, recovery, serialization and timestep ordering are corrected. F1 follows the warning policy. |
| R4-F2–F14 | RESOLVED | Harmonic interpolation, query shape, local overlap, guide placement, recovery, input identity and runtime checks are corrected. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identity, continuous intersections, protected recovery, early handlers and finite-depth checks address the original cases. |
| R6-F2–F11 | RESOLVED | Window, output, completion, failure-exit, inventory, native-value, quoting and path checks address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source triangulation, validation, protected relocation, verified CSV consumption and identity checks are corrected. |
| R8-F2–F13 | RESOLVED | Rollback, type normalization, sponge precision, parsing, caches, axis order, serialization and manifest selection are corrected. |
| R9-F2–F14 | RESOLVED | Original concurrency, recovery, index, dtype, parsing, trimming, depth-bound and library-path defects are corrected. |
| R10-F2–F14 | RESOLVED | Destination isolation, field validation, recovery tracking, parsing, QA normalization, serialization and the reported documentation defects are corrected. |
| R11-F2–F8 | RESOLVED | Finishing products are isolated; malformed arrays fail; relocation, unsupported water holes and spaced paths are handled. |
| R12-F2–F7 | RESOLVED | Accepted-build hashes, isolated benchmarks, record dimensions, positive duration, advancement and serialization address the original cases. |
| R13-F2–F9 | RESOLVED | Locks, paths, history validation, acceptance identity, land inventories, titles and sizing controls are corrected. |
| R14-F2–F9 | RESOLVED | Locks, marker ordering, scalar/lattice validation, scratch creation, fractional cadence, backend selection and job 401 documentation are corrected. Finding 9 below concerns other scripts. |
| R15-F2 | RESOLVED | The original covered-land gap is addressed under the approved resolution principle. Finding 6 concerns the new erosion implementation. |
| R15-F3–F7 | RESOLVED | Marker preservation, complex-scalar rejection, backend probing, isolated-element reporting and recipe-name validation are corrected. |
| R16-F1–F2 | RESOLVED | Environment and notices follow the approved optional private-backend policy. |
| R16-F3–F6 | RESOLVED | Positive-area ownership, control validation, complete-interface verification and the licensing comment address the original cases. |
| R17-F1–F9 | RESOLVED | Notices, environment paths, design protection, QA controls, missing-land handling, coordinate validation and interface counts are corrected. Later fixes close the publication aliases. |
| R18-F1–F6 | RESOLVED | Explicit boundary type, protected rejection/product paths, QA array validation and token-based output ownership address the original cases. |
| R19-F1 | PARTIAL | Connected-patch erosion fixes mainland dilution, but its distance uses an area-derived proxy rather than the specified median edge length—finding 6. |
| R19-F2 | RESOLVED | Driver tokens bind stages to their reservation; standalone stages reserve fresh outputs; 447 refuses finished builds. |
| R19-F3 | RESOLVED | Both 414 cases explicitly use boundary type 1, with read-back verification and recording. |
| R19-F4 | RESOLVED | QA geometry uses only the first two coordinate columns. |
| R19-F5 | RESOLVED | Unique temporary rejection reports and replacement protect aliased existing products. |
| R19-F6 | RESOLVED | Native grid record identifiers must be consecutive and ordered. Finding 7 concerns separate format checks. |
| R19-F7 | PARTIAL | Coordinates and connectivity are checked for the first history stack; later stacks escape that check—finding 2. |
| Generation-nondeterminism hypothesis | WITHDRAWN | Supplied repeated production measurements contradict it. |
| R15 small-island objection | WITHDRAWN | The owner permits dropping land below mesh resolution. No finding below challenges that policy. |

1. **Major — Standalone merge can accept an extension for a different designed boundary.**

   **Location:** [447_extend_merge.py:77](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:77), especially the boundary selection at line 99; [446_extend_generate.py:368](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:368).

   **Evidence:** 447 loads the current recipe but neither checks the generated mesh’s new boundary against its prescribed boundary nor consumes an identity-bound generation manifest. `check_expected` checks the current recipe against driver expectations; it does not establish which recipe produced `GEN`.

   I executed 447 in memory with mocked I/O, finishing and depth sampling, retaining the real merge, frozen-base checks, acceptance checks and QA. A generated boundary with **2 nodes** was accepted against a recipe prescribing **6 different boundary points**. The stage exited successfully with **16 nodes, 18 elements and QA 23/23**. No production mesh generation occurred.

   A realistic trigger is running 446, editing the boundary CSV or selecting another recipe, then running standalone 447 against the old generation directory.

   **Fix:** Record recipe, boundary, base and generated-artifact identities in `generate.json`; validate them in 447 before finishing. Also compare the generated new boundary’s ordered coordinates with the recipe’s boundary.

2. **Major — Later history stacks can belong to another mesh and still pass smoke completion. Residual R19-F7.**

   **Location:** [check_run.py:299](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:299).

   **Reproduction:** Two diskless NetCDF stacks had matching counts, finite fields and contiguous six-hour timestamps covering the requested day. Stack 0001 matched the staged grid. Stack 0002 had its x coordinates shifted by **100,000 m**. The complete `check_run` returned **`ok=True`, 5 records, no reasons**.

   The condition `not stamps[:-len(times)]` restricts `_grid_identity` to the first stack. Mixed or stale later histories therefore remain acceptable.

   **Fix:** Validate coordinates and connectivity in every history stack. Parse the staged grid once and reuse it.

3. **Major — Smoke completion accepts history from a different bathymetry case. Pre-existing.**

   **Location:** [check_run.py:221](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:221), [check_run.py:135](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:135).

   **Reproduction:** A virtual namelist named both the staged grid and depth file. The staged depths were **5 m**; a diskless history carried identical coordinates/connectivity but **`h=500 m`** at every node. All timestamps and required fields were valid. `check_run` returned **`ok=True`, no reasons**, and never read the staged depth file.

   This matters particularly for 453 experiments: cases intentionally share geometry while differing scientifically in depths. The existing geometry check cannot distinguish their histories.

   **Fix:** Parse `DEPTH_FILE`, validate its coordinates against the staged grid, and compare each history’s static bathymetry with the staged depths, accounting for declared output precision and explicitly applied model depth controls.

4. **Minor — The new history identity check discards NetCDF masks. Introduced by `94364a4`.**

   **Location:** [check_run.py:145](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:145), also line 152.

   **Reproduction:** In the complete diskless-history fixture, x had `_FillValue=1000` and a masked first node whose underlying value matched the staged x coordinate. `check_run` returned **`ok=True`, no reasons**.

   `np.asarray` removes the mask before comparison. The same conversion is used for `nv`, so missing connectivity can likewise be treated as valid underlying values.

   **Fix:** Reject masked coordinates/connectivity before conversion; validate dimensions, finite coordinates and integral connectivity explicitly.

5. **Minor — The 1 mm identity tolerance rejects valid float32 FVCOM histories. Introduced by `94364a4`.**

   **Location:** [check_run.py:147](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:147).

   **Reproduction:** Correct UTM coordinates near `(380000.123456789, 3900000.123456789)`, stored as NetCDF float32, differed from their staged doubles by up to **0.123456789 m**. `_grid_identity` rejected them as another mesh.

   Float32 output is supported by FVCOM: [mod_nctools.F:3964](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_nctools.F:3964) converts double variables to `NF90_FLOAT` for `DOUBLE_PRECISION && SINGLE_OUTPUT`. Thus a supported executable configuration can fail the smoke check solely because of serialization precision.

   **Fix:** Compare against staged values rounded to the history variable’s storage precision, or use a justified dtype-aware rounding tolerance. Retain mask, finite-value and connectivity checks.

6. **Minor — Land erosion uses an equilateral-area proxy instead of the promised median edge. Introduced by `94364a4`.**

   **Location:** [extend.py:422](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:422).

   **Reproduction:** A **36-node, 50-element** grid of 1,000 m squares split into right triangles covered `box(0, 2000, 5000, 3050)`, a **1,050 m wide** land strip.

   The actual median triangle-edge length is **1,000 m**. Erosion by the specified half-edge distance leaves **200,000 m²**, requiring rejection. The implementation instead computes an equivalent equilateral edge of **1,074.57 m**; its larger erosion removes the strip and returns `max_area_left_after_erosion_m2=0`.

   This is a gate-level synthetic reproduction, not a production build result.

   **Fix:** Compute actual edge lengths from the positive-area owning elements and use their median. Keep connected-patch erosion and the approved resolution principle.

7. **Minor — Native grid parsing still accepts malformed or surplus records. Pre-existing.**

   **Location:** [fvcom_native.py:548](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:548), also lines 560–565.

   **Reproduction:** With `Path.open` mocked to `StringIO`, actual `read_grd`:

   - Accepted a three-node grid with an additional fourth coordinate record, silently discarding it.
   - Accepted uniformly truncated node records containing only `NODE# X`, returning nodes of shape **`(3, 1)`** despite its `(N,2)` contract.

   The surplus-record case is incompatible with FVCOM’s coordinate reader, which reads through EOF and rejects more nodes than declared.

   **Fix:** Require the declared record count and minimum fields per record; validate final array shapes before returning. Preserve the documented allowance for trailing columns within records.

8. **Minor — Notebook 383 inherits experiment dates from the external template. Pre-existing portability defect.**

   **Location:** [383_m2_case_prep.py:165](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:165).

   **Reproduction:** Its declared experiment starts on **2021-01-01** and lasts **20 days**, but `namelist` only replaces `END_DATE`. Using the existing small template with `START_DATE` and `NC_FIRST_OUT` changed to 2012 produced:

   ```text
   START_DATE   = 2012-01-01 00:00:00
   NC_FIRST_OUT = 2012-01-01 00:00:00
   END_DATE     = 2021-01-21 00:00:00
   ```

   The current local template happens to match the constants. A different FVCOM checkout selected through `WORK_DIR` need not, causing the staged experiment to disagree with its declared timing.

   **Fix:** Explicitly set `START_DATE=START` and `NC_FIRST_OUT=START`, together with the required date format and timezone, and verify the emitted interval.

9. **Nit — Three sweep-script headers describe and submit the wrong experiment. Pre-existing.**

   **Location:** [397_finish_sweep.sh:12](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/397_finish_sweep.sh:12), [398_obc_sweep.sh:12](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/398_obc_sweep.sh:12), [399_obc_sweep2.sh:12](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/399_obc_sweep2.sh:12).

   **Evidence:** All three describe `field`, `achieved` and `ach1w` variants and instruct submission of `395_target_mode.sh`. Their actual variants respectively adjust finishing/generation controls, OBC controls, and OBC skip/curvature controls.

   **Fix:** Update each experiment description, variant list and `qsub` example to match that script.

Verification: targeted `pytest -q -s -p no:cacheprovider` selections covering extension, boundary design, recipe and DEM tests passed **53 tests**. Scoped `ruff check --no-cache` passed. `bash -n` passed for **all 56 job scripts**. An earlier broader pytest selection produced **48 passes, 16 deselections and 4 setup errors** because fixtures required a writable temporary directory unavailable in this read-only sandbox. The reproductions above used small synthetic inputs and mocks.

Production generation and FVCOM execution were not run. No large shared-data files were read and no batch jobs were submitted. The supplied real-run results were not independently rerun.

## Verdict
VERDICT: FAIL (0 blocker, 3 major, 5 minor, 1 nit)

### Prompt

```markdown
# Review request, round 20: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-19 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 19 (your previous answer; 7 findings) was fixed in 94364a4; read it.
Per finding:
- F1 `check_land_cover` rewritten: each connected patch of covered land is
  eroded by LAND_COVER_ERODE (0.5) times the median edge of the elements
  holding positive area of it; any area left fails. Measured on the real
  build before choosing: 3,681 patches, none survives 0.5. The rule keeps
  the owner's resolution principle (land narrower than an element is
  dropped on purpose).
- F2 `outdir.reserve(path, token)`, `outdir.claim(path, owner)`: 445 passes
  FMESH_OUT_TOKEN; 446 (owner = its parent) and 447 accept only that
  reservation; alone they reserve fresh; 447 refuses a finished build.
- F3 414: obc_type=1 for both cases, read back, recorded.
- F4 QA geometry on nodes[:, :2].
- F5 444 rejection report via mkstemp + os.replace; products may not be
  one file.
- F6 `read_grd` checks ids.
- F7 `check_run` compares history x, y and nv with the staged grid.
Verification after 94364a4: full test suite 1211 passed (batch
job 124336); on real data 444, 445 (stages accepted 445's reservation
token; QA 23/23, status ok, land cover 3,681 patches with nothing left
after erosion, grd bit-identical to rounds 4-18) and 453 (status ok)
passed; check_run, now tied to the staged grid, accepts the two real FVCOM
smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | major | code read (no identity link between 446 and 447) | yes | fixed, e661cf9 |
| F2 | major | test (second stack on another mesh) | yes, introduced in round 19 | fixed, e661cf9; test |
| F3 | major | test (h = 500 m); a real history's h equals its depth file exactly | yes | fixed, e661cf9; test |
| F4 | minor | test (masked x) | yes, introduced in round 19 | fixed, e661cf9; test |
| F5 | minor | test (float32 UTM) | yes, introduced in round 19 | fixed, e661cf9; test |
| F6 | minor | test (1,050 m strip) | yes, introduced in round 19 | fixed, e661cf9; test |
| F7 | minor | tests (surplus, short rows); every real grid in the tree still reads | yes | fixed, e661cf9; tests |
| F8 | minor | namelist printed | yes | fixed, e661cf9 |
| F9 | nit | job read | yes | fixed, e661cf9 |

## Round 21

### Answer (verbatim)

Found **1 major, 9 minor and 1 nit**. No files were modified.

Previous-finding statuses refer to their original reproductions. These inclusive ranges cover every numbered finding from rounds 1–20.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R15-F1 | RESOLVED | The approved GPL relicensing addresses the import restriction. |
| R1-F1–F26, F28–F29 | RESOLVED | Source handling, validation, serialization, seams, provenance, reservations and smoke preparation address the original cases. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Reservation, recovery, re-depth acceptance, interpolation, geometry and provenance checks address the reported cases. |
| R3-F1–F13 | RESOLVED | Final-field diagnostics, spacing, bands, ladders, recovery, serialization and timestep ordering are corrected. F1 follows the warning policy. |
| R4-F2–F14 | RESOLVED | Interpolation, query shape, overlap, guide placement, recovery, identity and runtime checks address the original cases. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identity, continuous intersections, protected recovery, handlers and finite-depth checks are corrected. |
| R6-F2–F11 | RESOLVED | Window, output, completion, failure-exit, inventory, native-value, quoting and path checks address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source interpolation, validation, protected relocation, CSV consumption and identity checks are corrected. |
| R8-F2–F13 | RESOLVED | Rollback, types, sponge precision, parsing, caches, axis order, serialization and manifest selection address the reported cases. |
| R9-F2–F14 | RESOLVED | Original concurrency, recovery, index, dtype, parsing, trimming, depth-bound and library-path defects are corrected. |
| R10-F2–F14 | RESOLVED | Isolation, field validation, recovery tracking, parsing, QA normalization, serialization and documentation are corrected. |
| R11-F2–F8 | RESOLVED | Finishing is isolated; malformed arrays fail; relocation, water-hole acceptance and spaced paths are handled. |
| R12-F2–F7 | RESOLVED | Accepted-build hashes, benchmark isolation, dimensions, duration, advancement and serialization address the original cases. |
| R13-F2–F9 | RESOLVED | Locks, paths, history validation, acceptance identity, land inventories, titles and sizing controls are corrected. |
| R14-F2–F9 | RESOLVED | Locks, markers, scalar/lattice validation, scratch creation, fractional cadence, backend selection and documentation are corrected. |
| R15-F2–F7 | RESOLVED | Covered-land checking now uses the approved resolution principle and actual edges; marker, scalar, backend, isolation and recipe-name fixes remain present. |
| R16-F1–F2 | RESOLVED | Environment and notices follow the approved optional-backend policy. |
| R16-F3–F6 | RESOLVED | Ownership, controls, complete-interface verification and the licensing comment address the original cases. |
| R17-F1–F9 | RESOLVED | Notices, environment paths, design protection, QA controls, land handling, coordinates and interface counts are corrected; later fixes close publication aliases. |
| R18-F1–F6 | RESOLVED | Explicit boundary types, protected products, QA arrays and token-based ownership address the original cases. |
| R19-F1–F7 | RESOLVED | Actual-edge erosion, reservation tokens, boundary types, coordinate slicing, atomic rejection, ordered IDs and every-stack identity checks address all seven cases. |
| R20-F1 | RESOLVED | Generation identities and ordered boundary comparison reject the original mismatched-design case. Findings 2–3 identify additional handoff gaps. |
| R20-F2–F5 | RESOLVED | Every stack checks coordinates, depths and integral connectivity, rejects masks, and accommodates storage precision. |
| R20-F6–F7 | RESOLVED | Erosion uses actual median edge lengths; GRD parsing enforces exact record counts and minimum fields. |
| R20-F8 | PARTIAL | Experiment dates are now explicit, but the required date format and timezone remain inherited; finding 6. |
| R20-F9 | RESOLVED | The three misleading job headers were corrected. |
| Generation-nondeterminism hypothesis; R15 small-island objection | WITHDRAWN | Supplied production measurements contradict nondeterminism; the owner explicitly permits dropping land below mesh resolution. |

1. **Major — Smoke completion accepts arbitrarily overlong, potentially stale history.**

   **Location:** [check_run.py:426](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:426).

   **Reproduction:** A diskless synthetic run declared START_DATE January 1, END_DATE January 2 and daily output. Its history instead contained January 1–5. With matching staged coordinates, depths and connectivity, finite fields and a successful log, `check_run` returned `ok=True`, five records and no reasons. The correctly bounded January 1–2 baseline also passed.

   The end check rejects only premature termination. Reusing a directory after shortening a run can therefore accept later stacks left by the previous run, even with the new mesh-identity checks.

   **Fix:** Reject history beyond END_DATE plus the permitted timing tolerance. Bind output acceptance to the current invocation, or require fresh history files when staging a reused directory.

2. **Minor — The new generation manifest can hash a replacement base rather than the base consumed. Introduced by e661cf9.**

   **Location:** [446_extend_generate.py:373](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:373); corresponding check at [447_extend_merge.py:83](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:83).

   **Evidence:** 446 reads the base at line 65, but hashes its files after meshing. Executing the actual manifest assignment against virtual files, with the GRD contents replaced between consumption and hashing, recorded the replacement’s digest rather than the consumed file’s digest. Standalone 447 compares against that same replacement.

   Driver 445’s initial provenance check catches a persistent change; this finding concerns standalone stages and the read/hash race.

   **Fix:** Capture input identities before consumption, bind parsing to those bytes, and recheck before publishing. Apply the same binding to 447’s base reads and identity checks.

3. **Minor — Generated land geometry is omitted from the generation identity.**

   **Location:** [447_extend_merge.py:91](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:91); manifest at [446_extend_generate.py:371](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/446_extend_generate.py:371).

   **Reproduction:** I executed 447 with synthetic geometry and mocked I/O, finishing and sampling, retaining the real merge, frozen-base checks, land/overlap acceptance and QA. The case had **30 nodes, 40 elements and constant 10 m depth**.

   With the original generated land, a **3,000,000 m²** covered-land rectangle produced a resolvable-patch rejection. Replacing only the supplied `land_with_base.shp` geometry with the base footprint left every checked manifest identity unchanged and produced successful acceptance with **QA 23/23** and zero reported covered-land patches.

   The land artifact controls both finishing and the acceptance reference. Its substitution consequently escapes the new identity check. Driver provenance identifies the original land source, not this intermediate dataset.

   **Fix:** Include the generated shapefile and all sidecars in the generation manifest’s inventory and digests; verify them before consumption and publication.

4. **Minor — Staged-depth validation ignores its header and accepts files FVCOM rejects. Introduced by e661cf9.**

   **Location:** [check_run.py:146](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:146).

   **Reproduction:** In an otherwise valid four-node synthetic run, replacing `Node Number = 4` with either `Node Number = 400` or `garbage` still produced `ok=True`. The parser discards the first line without validating it.

   FVCOM’s `READ_DEPTH` requires the header and checks its declared count against the grid.

   **Fix:** Use a strict shared DEP reader that validates the header, declared count, exact rows, coordinates and finite depths.

5. **Minor — Native readers silently discard surplus records and invent missing OBC types.**

   **Location:** [fvcom_native.py:592](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:592), [fvcom_native.py:616](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:616), [fvcom_native.py:622](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:622).

   **Reproduction:** Mocked file reads showed that a DEP declaring four nodes but containing five depth records returns four depths, silently ignoring the fifth. An OBC declaring two records containing only `1 2` and `2 3` returns node IDs `[1, 2]` and synthesizes types `[1, 1]`. `read_fvcom_case` accepted the combined malformed inputs.

   Both readers slice to the declared count, hiding trailing records. FVCOM requires explicit OBC types and rejects excess records.

   **Fix:** Enforce exact record counts, required fields and explicit valid OBC types. Reject missing types rather than assigning scientific boundary conditions implicitly.

6. **Minor — Experiment date format and timezone still depend on the external template. Residual R20-F8.**

   **Location:** [383_m2_case_prep.py:165](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:165).

   **Reproduction:** Supplying the small external template with `DATE_FORMAT='DMY'` and `TIMEZONE='EST'` produced:

   ```text
   DATE_FORMAT  = DMY
   TIMEZONE     = EST
   START_DATE   = 2021-01-01 00:00:00
   END_DATE     = 2021-01-21 00:00:00
   NC_FIRST_OUT = 2021-01-01 00:00:00
   ```

   The dates are emitted in YMD format while FVCOM is instructed to parse DMY. Independently, retaining EST shifts the experiment’s UTC interpretation. The inherited staging paths used by 414/448 share this dependency.

   **Fix:** Explicitly set `DATE_FORMAT='YMD'` and `TIMEZONE='UTC'`; verify staging against a template containing different valid settings.

7. **Minor — Duplicate YAML keys silently override scientific controls.**

   **Location:** [extend_recipe.py:56](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend_recipe.py:56); [444_design_obc.py:63](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/444_design_obc.py:63).

   **Reproduction:** Loading the actual extension recipe with duplicate entries accepted `rfactor: 0.2` followed by `rfactor: 0.8` as **0.8**, and `cfl_dt_s: 18` followed by `cfl_dt_s: 1` as **1**. Both loaders use ordinary `yaml.safe_load`, which silently keeps the last value.

   **Fix:** Use a shared safe loader that rejects duplicate keys in every mapping and reports their locations. The repository already has a duplicate-rejecting loader pattern in sizing code.

8. **Minor — Bathymetry sampling discards coordinate masks.**

   **Location:** [sources.py:131](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:131); [sources.py:353](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:353).

   **Reproduction:** A diskless 2×2 NetCDF grid had longitude `[139, 139.1]`, with the first coordinate masked through `_FillValue`. Despite the unknown coordinate, interpolation at `(139.05, 35.05)` returned **15 m**. Separately, passing a masked query longitude through `sample` to a synthetic source returned depth **10 m**, source index **0**, rather than rejecting or marking the query uncovered.

   `np.asarray` exposes masked arrays’ underlying coordinate values.

   **Fix:** Reject masked/nonfinite grid axes before conversion. Reject masked query coordinates or propagate their masks into uncovered output.

9. **Minor — Invalid timing controls escape the structured failure verdict.**

   **Location:** [check_run.py:269](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:269); [check_run.py:279](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:279).

   **Reproduction:** Mocked namelists with `END_DATE='garbage'` or an invalid calendar date raised uncaught `ValueError`. `_interval('seconds=1e999', '')` raised uncaught `OverflowError`; a zero interval returned a zero timedelta.

   Date parsing occurs outside a failure handler, and interval handling catches only `ValueError`.

   **Fix:** Convert date and interval parsing errors into failure reasons. Require finite, positive intervals and valid positive cycle controls before constructing timedeltas.

10. **Minor — Submission wrappers schedule chains without requiring DATA_DIR. Portability gap.**

    **Location:** [410_m2_chain.sh:20](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/410_m2_chain.sh:20); equivalent submission paths in `415_refine_m2_chain.sh` and `refine_workflow.sh`.

    **Reproduction:** Running 410 with DATA_DIR unset and WORK_DIR set, while mocking `qsub`, directory creation and temporary-directory allocation, exited successfully and attempted **five submissions**. No real jobs or directories were created.

    WORK_DIR is checked before submission; DATA_DIR is deferred to child scripts. Missing required configuration therefore schedules a chain that fails later instead of stopping at its entry point.

    **Fix:** Validate both variables at the beginning of each submission wrapper, before allocating directories or calling `qsub`.

11. **Nit — Unused job-ID parser remains in the submission script.**

    **Location:** [410_m2_chain.sh:31](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/410_m2_chain.sh:31).

    **Evidence:** `id()` is defined but never called; each submission instead extracts its ID using separate `grep`/`head` commands.

    **Fix:** Remove the dead function or use one shared parser consistently.

Verification: the read-only-compatible pytest selection passed **58 tests** using `pytest.main(['-q', '-s', '-p', 'no:cacheprovider', ...])`. An initial broader selection passed 58 tests but had **four setup errors** because indirect fixtures required writable temporary directories. Scoped `ruff check --no-cache` passed, and `bash -n` passed for all **56** OCTOPUS shell scripts. The reproductions above used small synthetic inputs and mocks. Production mesh generation and FVCOM execution were not rerun; no batch jobs were submitted.

## Verdict

VERDICT: FAIL (0 blocker, 1 major, 9 minor, 1 nit)

### Prompt

```markdown
# Review request, round 21: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-20 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 20 (your previous answer; 9 findings) was fixed in e661cf9; read it.
Per finding:
- F1 446 writes generate.json "inputs" (recipe, open boundary, base files,
  outer_utm.14 digests); 447 compares them with its own before finishing,
  and compares the new open boundary's nodes with the recipe's (1 mm,
  either direction).
- F2-F5 check_run: `_staged_mesh` (grid + DEPTH_FILE, checked against each
  other) and `_grid_identity` on every stack: x, y, h compared at the
  storage precision (float32: 2 ULP, at least 1 mm), masked values refused,
  nv integral (3, nele) and equal up to node order. A real FVCOM history
  has h equal to its depth file exactly.
- F6 erosion by the median of the real edge lengths.
- F7 `read_grd`: exact record count; cell rows >= 4 fields, node rows >= 3.
- F8 383 sets START_DATE and NC_FIRST_OUT.
- F9 job headers.
Verification after e661cf9: full test suite 1216 passed (batch
job 124412); on real data 444, 445 (447 accepted the generation identity
and the boundary match; QA 23/23, status ok, land cover nothing left after
erosion, grd bit-identical to rounds 4-19) and 453 (status ok) passed;
check_run, now checking x, y, h and nv on every stack, accepts the two real
FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | major | test (Jan 1-5 for a 1-day run) | yes | fixed, da3cb29; test |
| F2 | minor | code read (hashed after use) | yes, introduced in round 20 | fixed, da3cb29 |
| F3 | minor | code read (land not in the manifest) | yes | fixed, da3cb29 |
| F4 | minor | code read (header discarded) | yes, introduced in round 20 | fixed, da3cb29 |
| F5 | minor | code read; every real OBC and DEP file checked first | yes | fixed, da3cb29 |
| F6 | minor | namelist printed | yes | fixed, da3cb29 |
| F7 | minor | test (repeated key) | yes | fixed (yaml_strict), da3cb29; test |
| F8 | minor | tests (masked axis, masked query) | yes | fixed, da3cb29; tests |
| F9 | minor | tests (garbage date, 0 and 1e999 intervals) | yes | fixed, da3cb29; tests |
| F10 | minor | run with DATA_DIR unset | yes | fixed, da3cb29 |
| F11 | nit | code read | yes | fixed, da3cb29 |

## Round 22

### Answer (verbatim)

Found **1 major and 8 minor defects**. No files were modified.

Statuses below cover every numbered finding from rounds 1–21. Earlier findings are assessed against their original reproductions; additional gaps are reported separately.

**Previous-finding status**

| Findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; F1 in each round R4–R15 | RESOLVED | Approved GPL relicensing permits the imports. |
| R1-F1–F26, F28–F29 | RESOLVED | Source handling, depth limiting, spacing, seam checks, serialization, provenance and reservations address the original cases. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Reservations, publication recovery, re-depth acceptance, interpolation and geometry checks address the original cases. |
| R3-F1–F13 | RESOLVED | Final-field reporting, spacing, ladders, serialization and timestep ordering are corrected. F1 follows the warning policy. |
| R4-F2–F14 | RESOLVED | Interpolation, overlap, guide placement, recovery, code identity and runtime checks address the reported cases. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identity, continuous intersections, recovery protection and finite-depth checks remain present. |
| R6-F2–F11 | RESOLVED | Window, completion, failure-exit, inventory, native-value, quoting and path checks remain present. |
| R7-F2–F12 | RESOLVED | Whole-source interpolation, validation, protected relocation, verified CSV consumption and identity checks remain present. |
| R8-F2–F13 | RESOLVED | Rollback, type normalization, sponge precision, parsing, caches, axis order and manifest selection address the original cases. |
| R9-F2–F14 | RESOLVED | Concurrency, recovery, index validation, dtype enforcement, trimming, parsing and absolute library paths address the original cases. |
| R10-F2–F14 | RESOLVED | Isolation, field validation, recovery tracking, parsing, QA normalization and serialization remain corrected. |
| R11-F2–F8 | RESOLVED | Finishing isolation, malformed-array handling, relocation and rejection of holes containing no land address the original cases. Finding 6 concerns holes containing insufficient land. |
| R12-F2–F7 | RESOLVED | Accepted-build product hashes, benchmark isolation, dimensions, duration, advancement and serialization address the original cases. |
| R13-F2–F9 | RESOLVED | Locks, paths, history validation, acceptance hashes, inventories, titles and sizing controls remain corrected. |
| R14-F2–F9 | RESOLVED | Locks, marker ordering, scalar/lattice validation, scratch creation, fractional cadence and backend selection remain corrected. |
| R15-F2–F7 | RESOLVED | Land-cover checking follows the approved resolution principle and uses actual edges; the remaining original fixes are present. |
| R16-F1–F6 | RESOLVED | Notices and environment follow the approved backend policy; ownership, controls and complete-interface verification remain corrected. |
| R17-F1–F9 | RESOLVED | Environment paths, product protection, QA controls, land handling, coordinates and interface counts remain corrected. |
| R18-F1–F6 | RESOLVED | Explicit boundary types, protected products, QA array checks and reservation ownership address the original cases. |
| R19-F1–F7 | RESOLVED | Actual-edge erosion, reservation tokens, boundary types, coordinate slicing, atomic rejection, ordered IDs and every-stack identity checks address the original cases. |
| R20-F1–F7, F9 | RESOLVED | Generation identities, ordered boundary comparison, history identities, strict GRD records, actual-edge erosion and corrected job headers address the original cases. |
| R20-F8 | RESOLVED | 383 now explicitly sets dates, `DATE_FORMAT='YMD'` and `TIMEZONE='UTC'`. |
| R21-F1 | RESOLVED | History beyond END_DATE plus tolerance is rejected. |
| R21-F2 | RESOLVED | Both stages hash the base before reading and recheck before publishing their acceptance reports. |
| R21-F3 | PARTIAL | Generated land and sidecars are checked initially, but are not rechecked before merge acceptance; finding 3. |
| R21-F4–F5 | RESOLVED | Strict DEP reading, exact record counts and explicit OBC types address the original malformed files. Finding 8 concerns another required field. |
| R21-F6 | RESOLVED | Date format and timezone are explicit. |
| R21-F7 | RESOLVED | Duplicate mapping keys are rejected. Finding 5 is a compatibility regression introduced by that fix. |
| R21-F8 | RESOLVED | Masked grid axes are rejected; masked query coordinates become NaN. |
| R21-F9 | PARTIAL | Original END_DATE and ordinary interval errors become reasons, but START_DATE can escape validation and cycle controls remain incompletely checked; findings 2 and 4. |
| R21-F10–F11 | RESOLVED | Submission wrappers check both environment variables first; 410’s unused function is removed. |
| Generation-nondeterminism hypothesis; R15 small-island objection | WITHDRAWN | Production measurements contradict nondeterminism; the owner permits dropping land below mesh resolution. |

**Findings**

1. **Major — Re-depth accepts a build for a different prescribed boundary. Pre-existing.**

   **Location:** [453_redepth_extended.py:112](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:112), [453_redepth_extended.py:145](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:145).

   **Reproduction:** Executed the complete script with small synthetic geometry and mocked filesystem, sampling and plotting, retaining the real frozen-base, depth-limiter, geometry and QA checks. The accepted build came from recipe A; supplied recipe B prescribed its boundary **1,000 m farther east**, with the same base and case name.

   The script wrote `status="ok"` for **30 nodes, 40 elements, QA 23/23**, retaining A’s boundary. Base and new timestep allowances were both **100.96 s**.

   Product hashes establish which build files were read, but no check establishes that their geometry matches the supplied recipe.

   **Fix:** Verify the recipe’s base identity and ordered prescribed boundary against the accepted build and loaded mesh. Make any supported recipe overrides explicit.

2. **Minor — A valid NC_FIRST_OUT hides an invalid START_DATE. Introduced by da3cb29.**

   **Location:** [check_run.py:288](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:288), [check_run.py:328](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:328).

   **Reproduction:** A diskless two-record history spanning January 1–2, with finite fields and a successful log, used:

   ```fortran
   START_DATE='garbage',
   NC_FIRST_OUT='2020-01-01 00:00:00',
   END_DATE='2020-01-02 00:00:00',
   NC_OUT_INTERVAL='days=1',
   ```

   `check_run` returned **`ok=True`, `reasons=[]`**.

   `_date("START_DATE")` is skipped when NC_FIRST_OUT succeeds. The later parsing exception is suppressed under the incorrect assumption that it was already reported.

   **Fix:** Parse and validate START_DATE independently once; reuse that result for output-start fallback and advancement checks.

3. **Minor — Merge accepts generated land changed after its initial identity check. Residual R21-F3.**

   **Location:** [447_extend_merge.py:86](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:86), [447_extend_merge.py:95](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:95), [447_extend_merge.py:272](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:272).

   **Reproduction:** Executed the complete merge stage with synthetic geometry, retaining real repair, depth limiting, frozen-base verification, land checks, overlap checks and QA. After the manifest comparison, the mocked land read persistently replaced the original land with the base footprint.

   The original land contained a **5,400,000 m²** resolvable covered patch, which the real land-cover check rejected. The replacement produced successful merge acceptance with **30 nodes, 40 elements, QA 23/23**, and zero covered patches. Its current digests differed from the manifest.

   Only the base is rehashed before acceptance. Driver 445 checks raw inputs, rather than these generated artifacts.

   **Fix:** Recheck generated-land inventory/digests and `outer_utm.14` before acceptance. Consume immutable snapshots or bind parsing to the verified contents.

4. **Minor — Cycle controls can crash validation or produce invalid accepted intervals. Residual R21-F9.**

   **Location:** [check_run.py:224](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:224), [check_run.py:232](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:232).

   **Reproduction:** Direct in-memory calls produced:

   | Input | Result |
   |---|---|
   | `cycles=1`, `EXTSTEP_SECONDS=1`, `ISPLIT=1e999` | Uncaught `OverflowError` |
   | `cycles=1`, `EXTSTEP_SECONDS=-1`, `ISPLIT=-10` | Accepted 10-second interval |
   | `cycles=1`, `EXTSTEP_SECONDS=1`, `ISPLIT=1.9` | Silently truncated to a 1-second interval |
   | `seconds=1e-9` | Returned zero timedelta despite the positive-interval check |

   **Fix:** Validate each cycle control before conversion: finite positive step and positive integral ISPLIT. Convert parsing/conversion failures to failure reasons, and require a nonzero resulting timedelta.

5. **Minor — The duplicate-key fix rejects previously supported YAML merge keys. Introduced by da3cb29.**

   **Location:** [yaml_strict.py:19](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/yaml_strict.py:19).

   **Reproduction:**

   ```yaml
   depths:
     <<: {min_m: 3, max_m: null, rfactor: 0.2}
   ```

   `yaml.safe_load` returns the expected depth mapping. `load_unique` raises:

   ```text
   ConstructorError: could not determine a constructor for
   tag:yaml.org,2002:merge
   ```

   This mapping contains no duplicate or unknown depth keys. Replacing a recipe’s ordinary depths mapping with it therefore breaks a previously valid recipe.

   **Fix:** Preserve SafeLoader’s merge handling while detecting duplicate explicit keys. Test merged mappings, aliases and duplicate-key rejection.

6. **Minor — A tiny land patch legitimizes a much larger missing-water hole. Pre-existing.**

   **Location:** [extend.py:343](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:343).

   **Reproduction:** On a synthetic 1 km grid, removed eight triangles to create a **2 km × 2 km** hole and placed a **100 m × 100 m** islet inside it. After compacting the orphan node:

   - Mesh: **80 nodes, 120 elements**.
   - `check_island_holes`: accepted one island.
   - Land-cover check: zero covered patches.
   - QA: **23/23**.
   - Missing sea: **3,990,000 m²**.

   The hole check requires only positive land intersection, with no limit on excluded water.

   **Fix:** Check holes against the reference wet domain or apply a resolution-aware bound to water excluded around their land. This should still permit dropping sub-resolution islets.

7. **Minor — QA accepts masked mesh coordinates, depths and indices as known values. Pre-existing.**

   **Location:** [qa.py:633](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:633), [qa.py:658](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:658).

   **Reproduction:** Starting from `tests/test_qa.py::_pristine`, separately masked one coordinate, depth, connectivity entry or open-boundary ID while preserving its underlying value. Every variant reported **QA 23/23** on the **16-node, 18-element** fixture.

   `np.asarray` discards masks before integrity validation. Native export validation also accepted the masked-coordinate, masked-depth and masked-boundary variants.

   **Fix:** Reject masked required mesh values before conversion in QA and writer validation. Share the policy across coordinates, depths, connectivity and boundaries.

8. **Minor — The strict OBC reader never validates its first integer field. Pre-existing.**

   **Location:** [fvcom_native.py:627](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:627).

   **Reproduction:** Mocked reads of:

   ```text
   OBC Node Number = 2
   garbage 5 1
   2 9 1
   ```

   returned IDs `[4, 8]` and types `[1, 1]`. Replacing `garbage` with `1.5` also passed.

   FVCOM’s `READ_COLDSTART_OBC_GRID` reads all three fields as integers. Its initial scan skips the malformed first row, leaving fewer data rows than declared.

   **Fix:** Parse and validate the counter field as an integer. Its value need not be consecutive because FVCOM ignores the counter.

9. **Minor — Re-depth does not bind its land reference to the accepted build. Pre-existing omission exposed by the new manifest.**

   **Location:** [453_redepth_extended.py:112](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:112), [453_redepth_extended.py:118](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:118).

   **Evidence:** Accepted-build validation checks only GRD, DEP and OBC product hashes. It never compares the current generated land with `_built["generate"]["inputs"]["land_sha256"]`.

   Instead, it captures the current land as fresh provenance and later checks against that snapshot. A land file or sidecar replaced **before invocation** therefore escapes the accepted-build identity check. In the complete mocked 453 execution, differing accepted/current land digests did not prevent `status="ok"`.

   **Fix:** Compare the complete current land inventory and digests with the accepted build’s recorded generation identity before reading it; retain the existing checks for changes during the run.

**Verification**

- Read-only pytest selection: **63 passed in 0.97 s**, using `pytest.main(['-q', '-s', '-p', 'no:cacheprovider', …])` with bytecode disabled and filesystem fixtures excluded.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n`: **all 56 OCTOPUS scripts passed**.
- Reproductions used small synthetic arrays, diskless NetCDF and mocks.
- Full-suite tests, production meshing and FVCOM integration were **not rerun**. No shared data were read and no batch jobs were submitted.
- Git status remained unchanged.

## Verdict
VERDICT: FAIL (0 blocker, 1 major, 8 minor, 0 nit)

### Prompt

```markdown
# Review request, round 22: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-21 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 21 (your previous answer; 11 findings) was fixed in da3cb29; read it.
Per finding:
- F1 check_run: last record after END_DATE + tolerance fails.
- F2 446/447 hash the base before reading it and recheck before
  publishing.
- F3 generate.json "inputs" has land_sha256 over every shapefile sidecar;
  447 checks it.
- F4 check_run uses `read_dep` (strict).
- F5 `read_dep`/`read_obc`: exact counts; OBC types explicit 1-10.
- F6 383 sets DATE_FORMAT and TIMEZONE.
- F7 `yaml_strict.load_unique` in extend_recipe, base_recipe and 444.
- F8 `Grid.depth` refuses masked axes; `sample` fills masked queries NaN.
- F9 check_run date/interval errors become reasons.
- F10 410/415/refine_workflow check DATA_DIR and WORK_DIR first.
- F11 dead id() removed.
Verification after da3cb29: full test suite 1220 passed (batch
job 124529); on real data 444, 445 (QA 23/23, status ok, grd bit-identical
to rounds 4-20) and 453 (status ok) passed; check_run accepts the two real
FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | major | code read (no boundary/base link in 453) | yes | fixed, b27aa4a |
| F2 | minor | test (garbage START_DATE with NC_FIRST_OUT) | yes, introduced in round 21 | fixed, b27aa4a; test |
| F3 | minor | code read | yes | fixed, b27aa4a |
| F4 | minor | tests (ISPLIT 1e999, -10, 1.9; 1e-9 s) | yes | fixed, b27aa4a; test |
| F5 | minor | test (merge key) | yes, introduced in round 21 | fixed, b27aa4a; test |
| F6 | minor | measured on the real build's nine island holes; test (2 km hole, 100 m islet) | yes | fixed (water farther than half an edge from land < one element), b27aa4a; tests |
| F7 | minor | tests (masked nodes, depths, elements, OBC) | yes | fixed, b27aa4a; tests |
| F8 | minor | test ("garbage", "1.5" counters) | yes | fixed, b27aa4a; test |
| F9 | minor | code read | yes | fixed with F1, b27aa4a |

## Round 23

### Answer (verbatim)

Found **7 minor defects**. No files were modified.

**Previous-finding status**

The grouped ranges below cover every numbered finding from rounds 1–21. Statuses assess the reported defects; related new gaps are identified separately.

| Findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; F1 in each round R4–R15 | RESOLVED | Approved GPL relicensing resolves the import-policy conflict. |
| R1-F1–F26, F28–F29 | RESOLVED | Source handling, depth limiting, spacing, geometry, serialization, provenance and reservations address the original cases. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Reservations, recovery, re-depth validation, interpolation and geometry checks address the original cases. |
| R3-F1–F13 | RESOLVED | Final-field diagnostics, spacing, ladders, serialization and timestep ordering are corrected. F1 follows the approved warning policy. |
| R4-F2–F14 | RESOLVED | Interpolation, overlap, guide containment, recovery, identity and runtime checks remain corrected. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identity, continuous intersections, recovery and finite-depth checks remain present. |
| R6-F2–F11 | RESOLVED | Window, completion, failure-exit, inventory, native-value, quoting and path checks address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source interpolation, normalization, protected relocation, verified CSV consumption and provenance checks remain present. |
| R8-F2–F13 | RESOLVED | Rollback, types, sponge precision, parsing, caches, axis order and manifest selection remain corrected. |
| R9-F2–F14 | RESOLVED | Concurrency, recovery, indices, dtypes, trimming, parsing and library paths address the original cases. |
| R10-F2–F14 | RESOLVED | Isolation, field validation, recovery tracking, parsing, QA normalization and serialization remain corrected. |
| R11-F2–F6, F8 | RESOLVED | Finishing isolation, malformed-array handling, relocation and path checks remain corrected. |
| R11-F7 | PARTIAL | Large unsupported holes are rejected, but the new gate again accepts holes containing no land: finding 2. |
| R12-F2–F7 | RESOLVED | Accepted-product hashes, benchmark isolation, dimensions, duration, advancement and serialization address the original cases. |
| R13-F2–F9 | RESOLVED | Locks, paths, history validation, accepted-product binding, inventories, titles and sizing controls remain corrected. |
| R14-F2–F9 | RESOLVED | Locks, marker ordering, scalar/lattice validation, scratch creation, fractional cadence and backend selection remain corrected. |
| R15-F2–F7 | RESOLVED | Land coverage follows the accepted resolution principle and uses actual edges; the remaining fixes remain present. |
| R16-F1–F6 | RESOLVED | Notices/environment follow the accepted backend policy; ownership, controls and complete-interface checks remain corrected. |
| R17-F1–F9 | RESOLVED | Environment paths, product protection, QA controls, land handling, coordinates and interface counts remain corrected. |
| R18-F1–F6 | RESOLVED | Boundary types, protected products, QA array checks and reservation ownership address the original cases. |
| R19-F1–F7 | RESOLVED | Actual-edge erosion, ownership tokens, boundary types, coordinate slicing, atomic rejection, ordered IDs and every-stack identity checks remain present. |
| R20-F1–F9 | RESOLVED | Generation/history identities, strict GRD records, erosion, explicit dates/format/timezone and job headers remain corrected. |
| R21-F1–F6, F8–F11 | RESOLVED | Original history, base/land identity, native-record, date, DEM-mask, timing-control, wrapper and dead-code cases are addressed. Finding 5 is a separate date-relationship gap. |
| R21-F7 | PARTIAL | Ordinary duplicate keys are rejected; duplicates inside inline merged mappings bypass validation: finding 1. |
| Generation-nondeterminism hypothesis; R15 small-island objection | WITHDRAWN | Supplied production measurements contradict nondeterminism; the owner permits dropping land below mesh resolution. |

Round 22 individually:

| Finding | Status | Reason |
|---|---|---|
| F1 | RESOLVED | Generation boundary/base identities and ordered boundary coordinates are checked. |
| F2 | RESOLVED | `START_DATE` is parsed independently. |
| F3 | RESOLVED | 447 rechecks land inventory and `outer_utm.14` before accepting the merge. |
| F4 | RESOLVED | Cycle controls and intervals that round to zero are rejected. |
| F5 | RESOLVED | Legal YAML merges work again. The fix introduces finding 1. |
| F6 | RESOLVED | The original large-hole/tiny-islet case is rejected. The replacement gate introduces finding 2. |
| F7 | RESOLVED | Masked required mesh arrays are rejected by QA and native writers. Finding 6 concerns optional physics arrays. |
| F8 | PARTIAL | Counters undergo integer conversion, but Python-only integer spellings still pass: finding 7. |
| F9 | PARTIAL | Land replacement before invocation is rejected; replacement between acceptance and provenance capture escapes: finding 3. |

**Findings**

1. **Minor — Inline YAML merges bypass duplicate-key rejection. Introduced by b27aa4a.**  
   [yaml_strict.py:33](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/yaml_strict.py:33)

   `flatten_mapping()` recursively processes inline merge mappings without invoking `_mapping()` on those mappings. Their duplicate keys therefore escape the new pre-flattening check.

   Reproduction:

   ```python
   load_unique(
       "depths:\n"
       "  <<: {min_m: 3, max_m: null, rfactor: 0.2, rfactor: 0.8}\n"
   )
   ```

   This returns `rfactor: 0.8` without an error. Nested inline merges also bypass the check. A conflicting scientific control is silently overwritten despite `load_unique()` promising duplicate rejection in every mapping.

   **Fix:** Validate the original mapping-node graph recursively before flattening, tracking visited aliases. Preserve legitimate explicit overrides of merged values.

2. **Minor — The replacement island gate accepts an entire missing water element as an island. Introduced by b27aa4a.**  
   [extend.py:384](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:384)

   The old positive-land-intersection requirement disappeared. A hole now passes whenever its unsupported water area is no greater than the local median element area—even when it contains no land. Equality also contradicts the specified “below” threshold.

   In-memory reproduction: triangulate a regular 9×9 node grid at 1 km spacing, remove one interior triangle, assign depth 10 m, derive land boundaries, and supply empty land geometry. With valid outer-boundary topology and `n_base_nodes=27`, the **81-node, 127-element** mesh returns:

   ```text
   n_new_islands = 1
   max_water_in_hole_per_element = 1.0
   QA = 23/23
   ```

   The supposed island is **500,000 m² of missing open water**. Smaller-than-median missing triangles also pass, so changing `>` alone is insufficient.

   **Fix:** Retain a positive-area land-intersection requirement alongside the resolution-aware water limit. Enforce the specified strict limit with an explicit numerical tolerance.

3. **Minor — 453 does not bind accepted land identity to the provenance snapshot it consumes. Introduced by the new identity check.**  
   [453_redepth_extended.py:130](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:130), [453_redepth_extended.py:141](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:141)

   Land is compared with the accepted generation digest, then inventoried and hashed again by `collect()`. Unlike the recipe, report and built case, the collected land digest is never compared back to the accepted digest.

   Fault injection replaced land A with land B immediately after the initial acceptance hash. Provenance recorded B, consumption used B, and final change checks compared B with B. On a **54-node, 80-element** synthetic extension, the real geometry, frozen-base, r-factor and QA checks yielded:

   ```text
   status = "ok"
   QA = 23/23
   inputs_changed_during_run = []
   ```

   B added a harmless patch outside the mesh, isolating the identity defect from geometry failures. Filesystem/export operations and plotting were mocked; no products were written.

   **Fix:** Compare the collected land inventory and individual digests with `_gin["land_sha256"]` before consuming land. Bind collected base digests similarly, or consume verified immutable snapshots.

4. **Minor — Masked history timestamps become valid timestamps in `check_run`. Pre-existing.**  
   [check_run.py:371](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:371)

   `netCDF4.chartostring()` discards the mask on the character array while retaining its underlying characters. The checker consequently uses timestamps containing missing characters to establish cadence and coverage.

   A diskless NetCDF reproduction used two otherwise valid `Times` strings for January 1–2, 2020, with `_FillValue=b"2"`. Five timestamp characters were masked. With finite `zeta`, `ua`, `va`, matching dates and a completion log, the actual checker returned:

   ```text
   ok = True
   reasons = []
   ```

   **Fix:** Check the raw `Times` array for masked characters before calling `chartostring()` and reject any masked required timestamp character.

5. **Minor — An impossible `NC_FIRST_OUT` can pass run validation. Pre-existing.**  
   [check_run.py:302](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:302), [check_run.py:455](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:455)

   The checker validates `END_DATE > START_DATE`, but never requires `NC_FIRST_OUT` to fall inside that interval. Its one-output-interval tolerance can hide an invalid first-output setting.

   Diskless histories with finite fields and a completion log passed both cases:

   | Integration window | Output interval | `NC_FIRST_OUT` | History timestamps |
   |---|---|---|---|
   | Jan 1–2, 2020 | 1 day | Dec 31, 2019 | Jan 1 and Jan 2 |
   | Jan 1–2, 2020 | 1 day | Jan 3, 2020 | Jan 2 only |

   Both returned `ok=True, reasons=[]`. FVCOM’s cold-start initialization explicitly rejects first output before the start or after the end, so these configurations cannot establish a valid smoke run.

   **Fix:** Validate `START_DATE <= NC_FIRST_OUT <= END_DATE` independently of history tolerance. Continue allowing equality with `END_DATE`.

6. **Minor — Native Coriolis and sponge writers silently discard masks. Pre-existing.**  
   [fvcom_native.py:374](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:374), [fvcom_native.py:388](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:388)

   The round-22 mask protection covers mesh arrays, but `_check_cor()` and `_check_sponge()` still convert optional masked inputs with `np.asarray()`, exposing their underlying data.

   With a valid three-node mesh and writer output captured in memory:

   ```python
   cor = np.ma.array([35., 36., 37.], mask=[True, False, False])
   sponge = np.ma.array([[0., 100., .001]], mask=[[False, True, False]])
   ```

   `write_cor()` writes the missing value as `35.000000`; `write_spg()` writes the missing radius as `100.0`. Masked sponge node IDs and damping values also pass.

   **Fix:** Reject masks before converting Coriolis arrays or sponge rows. Preserve support for iterators by consuming them once and checking each resulting row.

7. **Minor — OBC integer validation accepts spellings FVCOM cannot read. Residual round-22 F8 defect.**  
   [fvcom_native.py:639](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:639)

   Python’s `int()` accepts underscores and Unicode digits, which are outside FVCOM’s list-directed integer input syntax.

   Actual reader reproduction using an in-memory file:

   ```text
   OBC Node Number = 2
   1_0 5 1
   2 9 1
   ```

   `read_obc(..., with_types=True)` accepts this and returns node indices `[4, 8]` with types `[1, 1]`. Replacing the first counter with full-width `１` also succeeds. FVCOM reads these fields as Fortran integers; the malformed row cannot be consumed as the declared OBC record.

   **Fix:** Require an ASCII integer token such as `[+-]?[0-9]+` before conversion. Apply the same lexical validation to node IDs and types.

**Verification**

- Read-only pytest selection: **106 passed, 114 deselected, 7.19 s**. Used `pytest.main()` with `-q -s -p no:cacheprovider`, disabled plugin autoload/bytecode, and an in-memory collection plugin excluding filesystem-writing fixtures. Covered extension, recipe, OBC design, DEM sources, run checking, QA and native I/O tests.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` on all **56 OCTOPUS shell scripts**: **passed**.
- Executed the small in-memory reproductions above. The initial 453 plotting setup failed because Matplotlib required a writable cache; the successful acceptance reproduction stubbed plotting.
- No full-suite, real-data meshing, production re-depth or FVCOM run was performed. No batch jobs were submitted or shared data read.
- Git status remained unchanged; the pre-existing untracked `recipes/extend/tokyo_bay_enshu_obc.png` was untouched.

## Verdict
VERDICT: FAIL (0 blocker, 0 major, 7 minor, 0 nit)

### Prompt

```markdown
# Review request, round 23: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-22 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 22 (your previous answer; 9 findings) was fixed in b27aa4a; read it.
Per finding:
- F1/F9 453: the accepted build's report["generate"]["inputs"] must match
  the recipe's open-boundary digest, the current base digests and the
  current land inventory digests; the build's open boundary must equal the
  recipe's coordinates (1 mm, either direction).
- F2 START_DATE parsed independently.
- F3 447 rechecks land and outer_utm.14 before merge.json.
- F4 cycle controls validated; zero-rounding intervals refused.
- F5 yaml_strict: duplicates among written keys; merge keys flattened.
- F6 `check_island_holes`: water inside a hole farther than half the local
  median edge from land must be below the local median element area. Real
  build: nine holes, worst 43,000 m2 vs ~170,000 m2 elements.
- F7 masked values refused in QA (node_index_valid) and in the writers.
- F8 OBC counter must be an integer.
Verification after b27aa4a: full test suite 1228 passed (batch
job 124612); on real data 444, 445 (QA 23/23, status ok, nine islands,
worst water-in-hole 0.27 of a local element, grd bit-identical to rounds
4-21) and 453 (accepted the build under the new identity checks; status ok)
passed; check_run accepts the two real FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | minor | test (inline and aliased merges) | yes, introduced in round 22 | fixed, 9e3cd06; test |
| F2 | minor | test (one lost triangle) | yes, introduced in round 22 | fixed, 9e3cd06; test |
| F3 | minor | code read | yes | fixed, 9e3cd06 |
| F4 | minor | test (masked Times) | yes | fixed, 9e3cd06; test |
| F5 | minor | tests (NC_FIRST_OUT before and after the run) | yes | fixed, 9e3cd06; test |
| F6 | minor | tests (masked cor, sponge) | yes | fixed, 9e3cd06; test |
| F7 | minor | tests ("1_0", full-width digit) | yes | fixed, 9e3cd06; test |

## Round 24

### Answer (verbatim)

Six substantiated **minor** findings remain. One is a regression introduced by `9e3cd06`. No files were modified.

## Previous-finding status

Ranges are inclusive; the status applies to every finding in each range. **RESOLVED** refers to the reported case, with newly discovered gaps listed separately below. These assessments use current code, recorded triage and the read-only checks described below; production verification was not repeated.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; F1 in R4–R15 | RESOLVED | Approved GPL relicensing resolves the original import-policy conflict. |
| R1-F1–F26, F28–F29 | RESOLVED | Source handling, limiting, spacing, geometry, serialization, provenance and reservations address the original cases. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Reservations, publication recovery, re-depth validation, interpolation and geometry checks address the reported cases. |
| R3-F1–F13 | RESOLVED | Final-field reporting, spacing, ladders, serialization and timestep ordering are corrected. F1 follows the approved warning policy. |
| R4-F2–F14 | RESOLVED | Interpolation, overlap, guide containment, recovery, identity and runtime fixes remain present. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identity, continuous intersections, recovery and finite-depth checks remain present. |
| R6-F2–F11 | RESOLVED | Window, completion, failure-exit, inventory, native-value, quoting and path checks address the original cases. |
| R7-F2–F12 | RESOLVED | Whole-source interpolation, normalization, protected relocation, CSV consumption and provenance checks remain present. |
| R8-F2–F13 | RESOLVED | Rollback, types, sponge precision, parsing, caches, axis order and manifest selection remain corrected. |
| R9-F2–F14 | RESOLVED | Concurrency, recovery, indices, dtypes, trimming, parsing and absolute library paths address the original cases. |
| R10-F2–F14 | RESOLVED | Isolation, field validation, recovery tracking, parsing, QA normalization and serialization remain corrected. |
| R11-F2–F6, F8 | RESOLVED | Finishing isolation, malformed-array handling, relocation and path checks remain corrected. |
| R11-F7 | PARTIAL | Ordinary unsupported holes are rejected, but holes bounded entirely by base nodes escape checking: finding 3. |
| R12-F2–F7 | RESOLVED | Accepted-product hashes, benchmark isolation, dimensions, duration, advancement and serialization address the original cases. |
| R13-F2–F9 | RESOLVED | Locks, paths, history validation, accepted-product binding, inventories, titles and sizing controls remain corrected. |
| R14-F2–F9 | RESOLVED | Locks, marker ordering, scalar/lattice validation, scratch creation, fractional cadence and backend selection remain corrected. |
| R15-F2–F7 | RESOLVED | Land coverage follows the accepted resolution principle and uses actual edges; the remaining fixes remain present. |
| R16-F1–F6 | RESOLVED | Notices/environment follow the approved optional-backend policy; ownership, controls and complete-interface checks remain corrected. |
| R17-F1–F9 | RESOLVED | Environment paths, product protection, QA controls, land handling, coordinates and interface counts remain corrected. |
| R18-F1–F6 | RESOLVED | Boundary types, protected products, QA array checks and reservation ownership address the original cases. |
| R19-F1–F7 | RESOLVED | Actual-edge erosion, reservation tokens, boundary types, coordinate slicing, atomic rejection, ordered IDs and every-stack identity checks remain present. |
| R20-F1–F9 | RESOLVED | Generation/history identities, GRD records, erosion, explicit dates/format/timezone and job headers address the original cases. |
| R21-F1–F6, F8–F11 | RESOLVED | Original history, base/land identity, native-record, date, DEM-mask, timing-control, wrapper and dead-code cases are addressed. |
| R21-F7 | PARTIAL | Ordinary and merged-map duplicates are rejected, but repeated written merge keys silently override: finding 2. |
| R22-F1 | RESOLVED | Generation boundary/base identities and ordered boundary coordinates are checked. |
| R22-F2 | RESOLVED | `START_DATE` is parsed independently. |
| R22-F3 | RESOLVED | 447 rechecks land inventory and generated mesh identity before acceptance. |
| R22-F4 | RESOLVED | Invalid cycle controls and intervals rounding to zero are rejected. |
| R22-F5 | PARTIAL | Simple legal merges work; reusing an overriding merged alias is now rejected: finding 1. |
| R22-F6 | RESOLVED | The reported large-hole/tiny-islet case is rejected. |
| R22-F7 | RESOLVED | QA and native writers reject masked required mesh arrays. Finding 4 concerns masks discarded before those checks. |
| R22-F8 | RESOLVED | All three OBC row fields receive ASCII integer validation. Finding 5 concerns other numeric fields. |
| R22-F9 | RESOLVED | Current and collected land identities are bound to the accepted build. |
| R23-F1 | RESOLVED | Original inline, aliased and sequence-merge duplicate cases are rejected. The fix introduces finding 1. |
| R23-F2 | RESOLVED | The reported missing-element hole now fails the positive-land and strict water-area requirements. Finding 3 bypasses those requirements entirely. |
| R23-F3 | RESOLVED | 453 compares collected land and base digests with accepted generation identities before consumption; confirmed by code inspection. |
| R23-F4 | RESOLVED | Masked `Times` characters are rejected; reproduced in memory. |
| R23-F5 | RESOLVED | First output outside `[START_DATE, END_DATE]` is rejected; reproduced in memory. |
| R23-F6 | RESOLVED | Masked Coriolis and sponge fields are rejected; reproduced in memory. |
| R23-F7 | RESOLVED | Underscore and Unicode spellings in OBC row fields are rejected; reproduced in memory. |
| Generation-nondeterminism hypothesis; R15 small-island objection | WITHDRAWN | Supplied production measurements contradict the hypothesis; the owner permits dropping land below mesh resolution. |

## Findings

1. **Minor — Reusing a legal overriding YAML merge alias raises a false duplicate error. Introduced by `9e3cd06`.**

   **Location:** [yaml_strict.py:50](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/yaml_strict.py:50).

   **Reproduction:**

   ```python
   load_unique(
       "x: &x {a: 1}\n"
       "y: &y {<<: *x, a: 2}\n"
       "z: {<<: *y}\n"
   )
   ```

   `yaml.safe_load()` correctly gives `a: 2` in both `y` and `z`. The current loader instead raises:

   ```text
   ValueError: duplicate key 'a' at line 2, column 16
   ```

   `flatten_mapping()` mutates the aliased mapping to contain inherited and explicit entries. A later `_check_keys()` call uses a fresh visited set and checks that already-flattened mapping, mistaking a legitimate override for a written duplicate.

   **Fix:** Validate the original composed mapping graph before any flattening, using loader-wide visitation or cached original-key validation. Preserve explicit overrides of merged values.

2. **Minor — Repeated written YAML merge keys bypass duplicate rejection. Pre-existing residual gap.**

   **Location:** [yaml_strict.py:31](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/yaml_strict.py:31).

   **Reproduction:**

   ```python
   load_unique(
       "depths:\n"
       "  <<: {rfactor: 0.2}\n"
       "  <<: {rfactor: 0.8}\n"
   )
   ```

   This returns `{"depths": {"rfactor": 0.8}}` without an error. Merge keys take the branch that bypasses `seen`, allowing conflicting scientific settings to be silently overwritten.

   **Fix:** Reject a second written `<<` key in the same mapping. Continue allowing one merge key containing a sequence of mappings and legitimate explicit overrides.

3. **Minor — New water holes bounded entirely by base nodes are exempt from island validation. Pre-existing.**

   **Location:** [extend.py:373](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:373).

   The condition `max(loop) < n_base_nodes` assumes that any such loop existed in the base. A new hole between nonconsecutive nodes of a concave interface can also satisfy it.

   **In-memory reproduction:** Used a 9×9 node grid at 1 km spacing, with standard two-triangle cells and the two lone outer corners retriangulated. The base contains cells with `i >= 6`, plus:

   ```text
   (3,2), (4,2), (5,2), (3,3), (5,3), (5,5)
   ```

   After placing base nodes first, removed the outer triangle originally indexed `[40, 41, 50]`, now `[19, 20, 24]`. All three nodes are below the base-node count of 38. Its removal preserves every constrained interface edge. The final OBC contains only new nodes.

   Actual results:

   ```text
   base:   38 nodes, 44 elements
   merged: 81 nodes, 127 elements
   frozen-base verification: passed, 18 interface edges
   overlap check: passed
   island check: n_new_islands=0, max_water_in_hole_per_element=0
   land-cover check: passed
   QA: 26/26
   ```

   With empty land geometry, the omitted triangle is **500,000 m² of missing open water**, yet the island check never examines it.

   **Fix:** Identify preserved base holes from the base’s actual boundary topology or geometry. Exempt only those unchanged holes; validate every newly created hole regardless of its node IDs.

4. **Minor — Merge and frozen-base verification discard masked depths, promoting unknown values to accepted depths. Pre-existing.**

   **Location:** [extend.py:239](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:239), [extend.py:259](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:259), [extend.py:526](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:526).

   **Reproduction:** Merged two adjacent square meshes, each containing four nodes and two triangles. The base depths were:

   ```python
   np.ma.array([5., 6., 7., 8.], mask=[True, False, False, False])
   ```

   After filling new depths, the merged depths were an ordinary array:

   ```text
   [5., 6., 7., 8., 8., 8.]
   ```

   `verify_frozen_base()` passed, and native export validation accepted the result. The unknown first depth had become a known 5 m value. The verifier also discards the original mask through `np.ascontiguousarray()`.

   Separately, `rfactor_smooth_free()` accepted a masked fixed depth of 10 m and used it to smooth its free neighbour. Thus downstream mask guards cannot protect values whose masks have already disappeared.

   **Fix:** Reject masked required inputs before conversion in merge, frozen-base verification and limiting. Apply the policy consistently to required coordinates, depths, connectivity and limiter inputs.

5. **Minor — Native headers, GRD fields and DEP fields still accept Python-only numeric spellings. Pre-existing.**

   **Location:** [fvcom_native.py:543](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:543), [fvcom_native.py:572](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:572), [fvcom_native.py:607](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:607).

   **Actual reader reproductions using mocked in-memory files:**

   | Input alteration | Accepted result |
   |---|---|
   | OBC header `OBC Node Number = ２` | Two OBC records accepted |
   | GRD cell counter `１` | Cell accepted |
   | GRD connectivity token `0_1` | Interpreted as node 1 |
   | DEP depth token `1_0` | Interpreted as 10 m |

   Python `int()` and `float()` accept these spellings. FVCOM consumes the corresponding fields through list-directed Fortran reads, as shown in [mod_input.F:4470](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_input.F:4470), [mod_input.F:4606](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_input.F:4606) and [mod_input.F:4913](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_input.F:4913). These spellings cannot establish valid native records.

   **Fix:** Share ASCII integer and supported Fortran-real lexical validation across all consumed numeric fields and headers before conversion. The OBC row-field fix alone is insufficient.

6. **Minor — `check_run` accepts timestamps on an unrelated dimension. Pre-existing.**

   **Location:** [check_run.py:376](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/cli/check_run.py:376).

   **Reproduction:** Created a diskless NetCDF history with:

   ```text
   time=2, stamp=2, DateStrLen=26
   Times(stamp, DateStrLen)
   zeta(time, node)
   ua(time, nele)
   va(time, nele)
   ```

   `Times` held valid January 1–2, 2020 strings. The fields were finite, the mocked namelist specified matching dates and cadence, and the log contained `TADA!`.

   The actual checker returned:

   ```text
   ok=True
   reasons=[]
   n_records=2
   ```

   Equal lengths conceal the fact that timestamps are not indexed by the fields’ record dimension. The checker validates field dimension names but never validates `Times` dimensions.

   **Fix:** Require the supported two-dimensional character layout with `time` first, and verify its record count against the `time` dimension before decoding.

## Verification

- Read-only pytest selection: **108 passed, 116 deselected, 6.73 s**. Used `pytest.main()` with `-q -s -p no:cacheprovider`, disabled bytecode/plugin autoload, and an in-memory collection plugin excluding filesystem-writing fixtures and optional-backend tests.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` on all **56 OCTOPUS scripts**: **passed**.
- Executed the synthetic reproductions above. The initial hole fixture failed two unrelated QA checks; correcting its outer-corner triangulation and OBC endpoints produced the reported **26/26** result.
- Portability changes and datum reporting were inspected; no additional defect was substantiated there.
- Full-suite, production mesh, production re-depth and FVCOM integration runs were **not performed**. No shared-data files were read and no batch jobs were submitted.
- Git status remained unchanged; the pre-existing untracked OBC PNG was untouched.

## Verdict
VERDICT: FAIL (0 blocker, 0 major, 6 minor, 0 nit)

### Prompt

```markdown
# Review request, round 24: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-23 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 23 (your previous answer; 7 minor findings, no major) was fixed in
9e3cd06; read it. Per finding:
- F1 `yaml_strict._check_keys` recurses into merged mappings (inline,
  aliased, sequences) before flattening; visited nodes tracked.
- F2 `check_island_holes` requires positive land in the hole AND water
  strictly below one local element.
- F3 453 compares the provenance's combined land and base digests with the
  accepted generation record.
- F4 masked Times refused.
- F5 NC_FIRST_OUT in [START_DATE, END_DATE].
- F6 masked cor and sponge refused.
- F7 OBC fields must match [+-]?[0-9]+.
Verification after 9e3cd06: full test suite 1232 passed (batch
job 124625); on real data 444, 445 (QA 23/23, status ok, nine islands all
holding land, grd bit-identical to rounds 4-22) and 453 (status ok, land
and base bound to the accepted build) passed; check_run accepts the two
real FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | minor | test (overriding alias reused) | yes, introduced in round 23 | fixed (graph check before construction), 0e6658a |
| F2 | minor | test (two << keys) | yes | fixed, 0e6658a |
| F3 | minor | test (hole among base-id nodes) | yes | fixed (base holes from the base's own loops), 0e6658a; test |
| F4 | minor | test (masked depths) | yes | fixed, 0e6658a; test |
| F5 | minor | tests (full-width digits, 1_0); every grid in use still reads | yes | fixed, 0e6658a; tests |
| F6 | minor | test (Times on another dimension) | yes | fixed, 0e6658a; test |

## Round 25

### Answer (verbatim)

Seven substantiated findings remain: **1 major and 6 minor**. The negative-count bypass in finding 3 was introduced by `0e6658a`; the other defects predate that commit. No files were modified.

## Previous-finding status

Ranges are inclusive: each status applies to every finding named. Assessments use current code, recorded triage and the synthetic checks below. Production verification was not repeated.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; F1 in R4–R15 | RESOLVED | Approved GPL relicensing resolves the original import-policy conflict. |
| R1-F1–F26, F28–F29 | RESOLVED | Source handling, limiting, spacing, geometry, serialization, provenance and reservations address the reported cases. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Reservation, recovery, re-depth, interpolation and geometry fixes remain present. |
| R3-F1–F13 | RESOLVED | Final-field reporting, spacing, ladders, serialization and timestep ordering are corrected. F1 follows the approved warning policy. |
| R4-F2–F14 | RESOLVED | Coefficient interpolation, overlap, containment, recovery, identity and runtime checks address the original cases. Finding 1 below concerns a separate convergence criterion. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identity, intersections, recovery and finite-depth checks remain corrected. |
| R6-F2–F11 | RESOLVED | Land-window, completion, failure-exit, inventory, native-value, quoting and path checks remain corrected. |
| R7-F2–F12 | RESOLVED | Whole-source interpolation, validation, protected relocation, verified CSV consumption and provenance checks remain present. |
| R8-F2–F13 | RESOLVED | Rollback, normalization, sponge precision, parsing, caches, axis order and manifest selection address the reported cases. |
| R9-F2–F14 | RESOLVED | Concurrency, recovery, merge indices, dtypes, trimming, parsing and library paths remain corrected. Finding 5 concerns the limiter’s separate index inputs. |
| R10-F2–F14 | RESOLVED | Isolation, field validation, recovery tracking, parsing, QA normalization and serialization remain corrected. |
| R11-F2–F6, F8 | RESOLVED | Finishing isolation, malformed-array, relocation and path fixes remain present. |
| R11-F7 | PARTIAL | Notebook paths reject the reported unsupported holes; the exported island checker retains a bypass when the element count is omitted. See finding 3. |
| R12-F2–F7 | RESOLVED | Accepted-product hashes, benchmark isolation, dimensions, duration, advancement and serialization remain corrected. |
| R13-F2–F9 | RESOLVED | Locks, paths, history validation, acceptance hashes, inventories, titles and sizing controls remain corrected. |
| R14-F2–F9 | RESOLVED | Locks, marker ordering, scalar/lattice validation, scratch creation, fractional cadence and backend selection remain corrected. |
| R15-F2–F7 | RESOLVED | Land coverage follows the approved resolution principle; the remaining reported cases are addressed. |
| R16-F1–F6 | RESOLVED | Notices and environment follow the approved backend policy; ownership, controls and complete-interface checks remain corrected. |
| R17-F1–F9 | RESOLVED | Environment paths, product protection, QA controls, land handling, coordinates and interface counts remain corrected. |
| R18-F1–F6 | RESOLVED | Boundary types, protected products, QA array validation and reservation ownership remain corrected. |
| R19-F1–F7 | RESOLVED | Actual-edge erosion, ownership tokens, boundary types, coordinate slicing, atomic rejection, ordered IDs and history identities remain corrected. |
| R20-F1–F9 | RESOLVED | Generation/history identities, strict GRD records, erosion, explicit dates/format/timezone and job headers address the reported cases. |
| R21-F1–F11 | RESOLVED | History limits, consumed-input identities, native records, dates, duplicate keys, DEM masks, timing controls, wrapper environments and dead code are addressed. |
| R22-F1–F9 | RESOLVED | Boundary/base and land identities, independent date parsing, timing validation, legal merges, land-hole checks, QA masks and OBC lexical checks address the original cases. |
| R23-F1–F7 | RESOLVED | Merged-map duplicates, positive-land requirements, accepted provenance, masked Times, output-date bounds, optional physics masks and ASCII integers are corrected. |
| R24-F1 | RESOLVED | Graph validation before construction permits reuse of legal overriding aliases; reproduced in memory. |
| R24-F2 | RESOLVED | A second written `<<` is rejected; reproduced in memory. |
| R24-F3 | PARTIAL | 447/453 pass the correct element count and the original explicit-count case fails correctly. The optional fallback and unvalidated count still bypass checking; finding 3. |
| R24-F4 | PARTIAL | Masked coordinates, depths, connectivity and limiter inputs are rejected. Masked interface and boundary indices still lose their masks; finding 2. |
| R24-F5 | RESOLVED | Native ASCII integers and Fortran reals are enforced. D exponents were accepted; underscore and Unicode numeric spellings were rejected in synthetic probes. |
| R24-F6 | RESOLVED | Times must use `(time, DateStrLen)` character dimensions; confirmed with diskless netCDF fixtures. |
| Generation-nondeterminism hypothesis; R15 small-island objection | WITHDRAWN | Supplied production measurements contradict nondeterminism; the owner permits omission of land below mesh resolution. |

## Findings

1. **Major — The spinup gate accepts substantial phase drift as convergence.**

   **Location:** [384_m2_analysis.py:306](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/384_m2_analysis.py:306).

   The gate checks only the difference between half-window amplitudes. It computes phase changes at line 293, but excludes them from `spinup.converged`; the report then tells readers to trust constants where that flag is true.

   **Reproduction:** Execute the actual `analyze()` function with mocked reads/writes and a three-node mesh carrying this unit-amplitude signal:
   ```python
   P = 44714.
   t = np.arange(960) * P / 48
   z = np.cos(2*np.pi*t/P - np.deg2rad(60 + 60*t/(20*P)))
   ```
   Analysis completed successfully and reported:

   - Half-window amplitude change: **0.000535 m**, below the 0.002 m tolerance.
   - Half-window phase change: **29.762°**.
   - `spinup.converged: true`.
   - Full-window amplitude: **0.958903 m**, although the signal’s amplitude is 1 m.

   Phase evolution during a transient can therefore produce accepted, biased harmonic constants.

   **Fix:** Gate agreement of the interpolated complex harmonic coefficients, for example `abs(h2 - h1) < TOL_M`, and report that measure. This checks amplitude and phase together without unstable phase thresholds at negligible amplitude.

2. **Minor — Merge and frozen-base verification still discard masked interface and boundary indices.**

   **Locations:** [extend.py:231](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:231), [extend.py:289](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:289), [fvcom_native.py:324](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/io/fvcom_native.py:324).

   The new `_no_masks` guards omit `interface_outer`, `interface_base`, `outer_open` and the original base’s open-boundary chains. `_indices()` and `_edges()` convert these with `np.asarray`, exposing masked backing values as ordinary indices.

   **Reproduction:** Merge two adjacent 1 km squares. Mask the first entry of each interface/open-boundary array separately, using valid backing indices. Every case is accepted. The result has ordinary open-boundary indices `[4, 5]`, `verify_frozen_base()` reports one valid interface edge, and `_check_exportable()` accepts the result after filling new depths. A masked original base OBC chain also passes verification.

   **Fix:** Reject active masks inside `_indices()` before conversion, and validate the original base boundary chains before `_edges()` strips their representation.

3. **Minor — The island checker retains its original default bypass and adds an unchecked-count bypass.**

   **Locations:** [extend.py:388](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:388), [extend.py:396](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:396).

   Omitting `n_base_elements` still exempts every loop whose node IDs are below `n_base_nodes`, including newly created holes. When supplied, the count is used directly as a Python slice: a negative count can classify new holes as base holes. The latter bypass is introduced by `0e6658a`.

   **Reproduction:** Use the new `test_a_new_hole_among_base_nodes_is_checked` fixture: a 6×6 lattice at 1 km spacing, with triangle 14 removed and empty land.

   ```python
   check_island_holes(mesh, Polygon(), 18)      # succeeds: n_new_islands == 0
   check_island_holes(mesh, Polygon(), 18, -1)  # succeeds: n_new_islands == 0
   check_island_holes(mesh, Polygon(), 18, 10)  # rejects missing sea
   ```

   The last call correctly identifies **500,000 m² of missing open water**. The current test explicitly preserves acceptance by the unsafe default. The 447/453 calls supply correct counts and avoid these bypasses.

   **Fix:** Require trustworthy base topology—an element count or explicit base-hole loops—and validate counts as integers within mesh bounds before slicing. Remove the node-ID-only exemption.

4. **Minor — Sizing and boundary-design helpers silently use masked backing data.**

   **Locations:** [extend.py:71](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:71), [extend.py:134](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:134), [obc_design.py:143](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:143).

   These conversions discard masks before validating required values.

   **Reproduction:** On a 2×2, 1 km lattice, supply ambient sizes of 4,000 m and a fully masked band whose backing data are 1,500 m. `compose_sizing(..., grade=.1)` returns an ordinary 1,500 m field and certifies `band_0_cells == 4` with zero deviation. All four constraints were masked.

   Similarly, `resample()` accepts a line with its endpoint’s x coordinate masked and returns ordinary coordinates extending to the hidden 10,000 m endpoint.

   **Fix:** Check masks before converting sizing values, lattice coordinates, bands and design coordinates. Required coordinates should reject masks; optional missing constraints need an explicit documented representation.

5. **Minor — The r-factor limiter accepts negative node indices and smooths the wrong node.**

   **Location:** [extend.py:556](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:556).

   Edge endpoints are converted but never checked for valid mesh indices. NumPy interprets negative indices from the end.

   **Reproduction:**
   ```python
   rfactor_smooth_free(
       np.array([100., 10.]),
       np.array([-1]), np.array([0]),
       np.array([False, True]),
       rmax=.2, hmin=1,
   )
   ```
   This succeeds, returning depths approximately `[100., 66.666667]`, one iteration and maximum r approximately 0.2. An invalid endpoint silently becomes node 1.

   **Fix:** Validate matching one-dimensional endpoint arrays, integral indices in `[0, len(h0))`, and matching depth/free-mask shapes before indexing. Reuse the checked index helper after fixing its mask handling.

6. **Minor — Strict YAML rejects the legal unquoted scalar key `=`.**

   **Location:** [yaml_strict.py:53](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/yaml_strict.py:53).

   The graph checker constructs keys before SafeLoader performs its special handling of `tag:yaml.org,2002:value`.

   **Reproduction:**
   ```python
   yaml.safe_load("=: 1\n")  # {'=': 1}
   load_unique("=: 1\n")     # ConstructorError
   ```
   The same failure occurs with the parent version of `yaml_strict.py`; it is pre-existing, not introduced by round 24.

   **Fix:** Treat the value-tagged scalar key as `"="` during duplicate checking, consistently with SafeLoader. Also detect duplicates between quoted and unquoted spellings.

7. **Minor — An accepted recipe depth floor can conflict with final QA’s hardcoded 2 m floor.**

   **Locations:** [extend_recipe.py:155](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend_recipe.py:155), [447_extend_merge.py:256](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:256), [453_redepth_extended.py:225](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:225).

   The recipe permits any positive `min_m`. Both final stages omit `min_depth_m` when calling QA, whose default is 2 m.

   **Reproduction:** The mocked recipe loader accepts `depths.min_m = 1`. On the pristine 4×4 QA fixture with uniform 1.5 m depths:

   - Default QA fails only `min_depth_clip`, with **16 violations**.
   - `run_qa(mesh, min_depth_m=1)` passes **23/23 gates**.

   Thus supported recipe settings can produce depths that satisfy the requested bounds and r-factor but fail a different, implicit final threshold.

   **Fix:** Make the accepted recipe and QA policies consistent. Apply configured bounds explicitly to new nodes and document the frozen-base depth policy. If 2 m is mandatory universally, reject incompatible recipes before generation.

## Verification

- In-memory pytest collection across extension, recipe, OBC design, DEM sources, native I/O, run checking and QA tests: **110 passed, 121 deselected**. Tests requiring filesystem fixtures or optional backends were excluded; cache writing and bytecode generation were disabled.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` on all **56** OCTOPUS shell scripts: **passed**.
- Findings above were exercised with small synthetic arrays, mocked I/O or diskless netCDF. Two initial probe setup errors—an incorrect import and `stringtochar` incompatibility—were corrected.
- Git status remained unchanged, including the pre-existing untracked OBC PNG. No shared-data reads or batch submissions occurred. Production mesh generation, implied timesteps and FVCOM smoke runs were **not rerun**.

## Verdict
VERDICT: FAIL (0 blocker, 1 major, 6 minor, 0 nit)

### Prompt

```markdown
# Review request, round 25: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-24 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 24 (your previous answer; 6 minor findings) was fixed in 0e6658a;
read it. Per finding:
- F1/F2 `yaml_strict`: `_UniqueLoader.construct_document` walks the
  composed node graph once (visited by id) before construction; in each
  mapping, duplicates among written keys and a second '<<' are errors.
  SafeLoader's own merge handling then constructs the document.
- F3 `check_island_holes(..., n_base_elements)`: base holes are the closed
  boundary loops of the base elements; 447 and 453 pass it.
- F4 `_no_masks` in merge_outer, verify_frozen_base, rfactor_smooth_free.
- F5 `_fint` / `_freal` (ASCII integer; Fortran real with D exponents)
  for native headers, GRD and DEP fields; OBC uses the same integer pattern.
- F6 Times must be (time, DateStrLen) characters.
Verification after 0e6658a: full test suite 1239 passed (batch
job 124631); on real data 444, 445 (QA 23/23, status ok, nine new islands
checked, the base's own holes left alone, grd bit-identical to rounds
4-23) and 453 (status ok) passed; check_run accepts the two real FVCOM
smoke histories; every native grid in use reads under the lexical checks.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | major | arithmetic (30 deg drift: |h2-h1| = 0.52 m vs 2 mm) | yes | fixed (complex-coefficient gate), 02b4ca2 |
| F2 | minor | test (masked interface index) | yes | fixed, 02b4ca2; test |
| F3 | minor | tests (default and negative count) | yes | fixed (count required and checked), 02b4ca2; tests |
| F4 | minor | tests (masked band, masked line) | yes | fixed, 02b4ca2; tests |
| F5 | minor | test (edge end -1) | yes | fixed, 02b4ca2; test |
| F6 | minor | test ('=' key) | yes | fixed, 02b4ca2; test |
| F7 | minor | code read (QA default 2 m) | yes | fixed (recipe floor, or the base's own minimum if lower), 02b4ca2 |

## Round 26

### Answer (verbatim)

Seven findings: **0 blocker, 0 major, 7 minor, 0 nit**. Finding 1 is introduced by `02b4ca2`; finding 2 is an incomplete round 25 fix. The others predate that commit. No files were modified.

Previous-finding statuses below use current code, recorded triage and small synthetic checks. Ranges are inclusive; the status applies to every named finding. Production verification was not repeated.

| Previous findings | Status | Reason |
|---|---|---|
| R1-F27; R2-F22; R3-F14; R4–R15-F1 | RESOLVED | Approved GPL relicensing resolves the original import-policy conflict. |
| R1-F1–F26, F28–F29 | RESOLVED | Interpolation, geometry, serialization, provenance, reservations and failure handling address the reported cases. F4 follows the approved warning policy. |
| R2-F1–F21 | RESOLVED | Atomic reservations, publication recovery, re-depth acceptance and geometry corrections remain present. |
| R3-F1–F13 | RESOLVED | Final-field reporting, chord spacing, ladders, serialization and timestep ordering remain corrected. F1 follows the approved warning policy. |
| R4-F2–F14 | RESOLVED | Coefficient interpolation, overlap, containment, recovery, identities and runtime checks remain corrected. |
| R5-F2–F10 | RESOLVED | Flip guards, parsed-byte identities, intersection checks, recovery and finite-depth checks remain present. |
| R6-F2–F11 | RESOLVED | Land-window, completion, failure-exit, inventory, native-value, quoting and path corrections remain present. |
| R7-F2–F12 | RESOLVED | Whole-source interpolation, validation, relocation protection and input/code identity checks remain present. |
| R8-F2–F13 | RESOLVED | Rollback, normalization, sponge precision, parsing, cache identity, axis order and manifest selection remain corrected. |
| R9-F2–F14 | RESOLVED | Concurrency, recovery, merge indices, dtypes, trimming, parsing and library paths address the original cases. |
| R10-F2–F14 | RESOLVED | Isolation, field validation, recovery tracking, parsing, QA normalization and serialization remain corrected. |
| R11-F2–F8 | RESOLVED | Finishing isolation, malformed-array and relocation fixes remain present. F7’s omitted-count island-check bypass is now removed. |
| R12-F2–F7 | RESOLVED | Accepted-product identities, benchmark isolation, dimensions, duration, advancement and serialization remain corrected. |
| R13-F2–F9 | RESOLVED | Locks, paths, history validation, acceptance hashes, inventories, titles and sizing controls remain corrected. |
| R14-F2–F9 | RESOLVED | Locks, marker ordering, lattice validation, scratch creation, fractional cadence and backend selection remain corrected. |
| R15-F2–F7 | RESOLVED | Land checks follow the approved resolution policy; remaining reported cases remain corrected. |
| R16-F1–F6 | RESOLVED | Optional-backend handling follows the owner’s decision; ownership, controls and complete-interface checks remain corrected. |
| R17-F1–F9 | RESOLVED | Environment paths, product protection, QA controls, land handling, coordinates and interface counts remain corrected. |
| R18-F1–F6 | RESOLVED | Boundary types, protected products, QA array validation and reservation ownership remain corrected. |
| R19-F1–F7 | RESOLVED | Actual-edge erosion, ownership tokens, boundary types, coordinate slicing, atomic rejection, ordered IDs and history identities remain corrected. |
| R20-F1–F9 | RESOLVED | Generation/history identities, strict GRD records, erosion, explicit dates/format/timezone and job headers remain corrected. |
| R21-F1–F11 | RESOLVED | History limits, consumed-input identities, native records, dates, duplicate keys, DEM masks, timing controls, wrapper environments and dead code remain addressed. |
| R22-F1–F9 | RESOLVED | Boundary/base and land identities, date parsing, timing validation, YAML merges, land-hole checks, QA masks and OBC lexical checks remain corrected. |
| R23-F1–F7 | RESOLVED | Merged-map duplicates, positive land in holes, accepted provenance, masked Times, output-date bounds, optional physics masks and ASCII integers remain corrected. |
| R24-F1–F6 | RESOLVED | Graph-first YAML checking, merge-key rejection, required/validated base-element counts, mesh/index mask checks, native lexical checks and Times dimensions address the original cases. |
| R25-F1 | RESOLVED | Actual analysis code reports coefficient change **0.509616 m** and `converged: false` for the drifting-phase reproduction; the steady signal passes. |
| R25-F2 | RESOLVED | `_indices` rejects active masks; original base OBC chains also receive mask checks. Finding 3 concerns separate fractional indices. |
| R25-F3 | RESOLVED | The element count is required and validated; the node-ID fallback is gone. |
| R25-F4 | PARTIAL | Masked sizing values, bands, floor, line coordinates and targets are rejected, but `band_field` still discards lattice masks: finding 2. |
| R25-F5 | RESOLVED | Depth/free-mask shapes, boolean flags, endpoint shapes and endpoint indices are checked. |
| R25-F6 | RESOLVED | Plain `=` loads as a string; duplicate quoted/unquoted spellings are rejected. |
| R25-F7 | RESOLVED | The original implicit 2 m floor conflict is corrected. The replacement introduces finding 1. |
| Generation-nondeterminism hypothesis; R15 objection to omitted land below mesh resolution | WITHDRAWN | Supplied repeated-run measurements and the owner’s resolution policy remain controlling. |

1. **Minor — The new QA floor can silently accept nonpositive frozen-base depths.**

   **Locations:** [447_extend_merge.py:259](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/447_extend_merge.py:259), [453_redepth_extended.py:228](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/453_redepth_extended.py:228).

   `min(recipe_floor, base_minimum)` automatically lowers the depth gate to zero or a negative value when the base contains such a depth. There is no independent check of the adopted base’s positive-depth sea policy. The limiter checks only free nodes and endpoints of edges touching them, so a bad depth inside the base escapes that check.

   **Reproduction:** Split the pristine 4×4 QA fixture into a frozen eastern base and western extension, and set one base node away from the interface to zero depth. Actual checks returned:

   - Base: **12 nodes, 12 elements**; merged mesh: **16 nodes, 18 elements**.
   - Limiter: **0 iterations**, maximum free-edge r-factor **0**.
   - Frozen-base verification: passed, including **3 interface edges**.
   - Island and overlap checks: passed.
   - QA with the new floor calculation: **23/23**.

   A separate pristine-fixture probe also passed **23/23** with a base depth of **−1 m**. These are injected malformed-base cases, hence minor.

   **Fix:** Validate the base against the supported depth policy before using its minimum. For these positive-depth sea tools, require finite, strictly positive base depths. If nonpositive depths are deliberately supported, require an explicit policy and corresponding QA instead of silently lowering the gate.

2. **Minor — `band_field` still converts masked lattice coordinates into certified constraints.**

   **Location:** [extend.py:100](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:100).

   The added mask guard covers `line_xy` and `targets`, but omits `x` and `y`. Shapely’s conversion strips their masks. The returned ordinary array subsequently passes `compose_sizing`’s new guards.

   **Reproduction:** On a 2×2 lattice at 1 km spacing, fully mask `x` or `y`, then call:

   ```python
   b = band_field(masked_x, y, [[0, 0], [1000, 0]], [1500, 1500], 2000)
   h, report = compose_sizing(
       np.full((2, 2), 4000.), x, y, grade=.1, bands=[b]
   )
   ```

   The result is an ordinary **1500 m** field, with `band_0_cells: 4` and `band_0_max_rel_deviation: 0`. Unknown coordinates have become verified constraints.

   **Fix:** Reject masks on both lattice coordinate arrays before constructing Shapely points; also validate their matching lattice shape and finiteness there.

3. **Minor — Frozen-base verification truncates fractional original OBC indices.**

   **Location:** [extend.py:295](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:295).

   `interface_base` goes through `_indices`, but the base’s original chains go directly through `_edges`, which casts to `int64`. The new mask check does not validate their values.

   **Reproduction:** Merge two adjacent squares with interface `[1, 2]`, then replace the original base OBC chain with `[1.9, 2.9]`. `verify_frozen_base(merged, base, [1, 2])` succeeds and reports one verified interface edge.

   Thus it certifies the interface after interpreting invalid original indices as different, valid ones.

   **Fix:** Pass every original base OBC chain through `_indices(c, nb, ...)` before computing its edge set.

4. **Minor — Masked mutability can authorize deletion of an unknown/frozen element.**

   **Location:** [extend.py:640](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:640).

   `trim_lone_corners` converts `mutable` with `np.asarray(..., bool)`, exposing masked backing flags as deletion permissions.

   **Reproduction:** Use the existing six-triangle fan plus cape triangle `[1, 7, 2]`. Mark the six fan triangles immutable and mask the cape’s mutability flag, whose backing value is `True`. The function drops the cape and returns:

   ```text
   n_elements_dropped: 1
   lone_nodes_left: []
   ```

   The unknown permission was treated as affirmative permission.

   **Fix:** Reject active masks on connectivity, mutability and protected-node inputs before conversion. Validate a boolean mutability vector matching the element count.

5. **Minor — Depth rounding accepts unknown depths and nonfinite bounds.**

   **Location:** [extend.py:538](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/extend.py:538).

   `round_depths_inside` has no mask or finite-bound validation. Its ordered-bound check does not catch NaN.

   **Reproduction:**

   ```python
   round_depths_inside(np.ma.array([999.], mask=True), 3., None)
   # array([999.])

   round_depths_inside([5.], np.nan, None)
   # array([nan])

   round_depths_inside([5.], 1., np.nan)
   # array([nan])
   ```

   The exported helper therefore returns either hidden depths or a result that cannot satisfy its stated interval. Recipe validation prevents the NaN-bound cases in current notebook calls.

   **Fix:** Reject active depth masks and nonfinite depth values; validate finite ordered bounds before computing the rounding interval.

6. **Minor — Ray intersection accepts masked coordinates and returns NaNs as a successful intersection.**

   **Location:** [obc_design.py:70](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/obc_design.py:70).

   Conversion discards masks, while NaN parameters bypass the `s <= 0 or t <= 0` rejection.

   **Reproduction:**

   ```python
   ray_intersection(
       np.ma.array([0., 0.], mask=[True, False]),
       [1, 0], [1, 1], [0, -1]
   )
   # array([1., 0.])

   ray_intersection([np.nan, 0], [1, 0], [1, 1], [0, -1])
   # array([nan, nan])
   ```

   Both return normally despite lacking a known valid ray intersection. Notebook 444’s later fillet validation catches nonfinite vertices, but callers of this exported helper receive a false success.

   **Fix:** Reject masks, require finite `(2,)` vectors and nonzero directions, and require finite, strictly positive solved parameters.

7. **Minor — Datum accounting silently uses masked or invalid source indices.**

   **Location:** [dem/sources.py:386](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/dem/sources.py:386).

   `non_tp_count` strips masks and does not validate attribution indices.

   **Reproduction:**

   ```python
   non_tp_count(
       ["m7001", "srtm15plus"],
       np.ma.array([1, 1], mask=True)
   )
   # (2, ["srtm15plus"])
   ```

   Two unknown attributions are reported as definite MSL-source depths. Out-of-range indices are instead silently omitted from the count. Current `sample` returns valid unmasked indices, so this concerns unusual inputs to the exported accounting helper.

   **Fix:** Reject active masks and validate integral indices in `{-1, 0, …, len(names)-1}` before counting.

Verification:

- Filesystem-free pytest selections, using `--capture=sys -p no:cacheprovider` and excluding filesystem fixtures/optional backends: **116 passed, 135 deselected** across two runs.
- Scoped `ruff check --no-cache`: **passed**.
- `bash -n` on all **56** OCTOPUS shell scripts: **passed**.
- Findings reproduced with small arrays and in-memory mocks. The spinup reproduction executed the actual analysis function bodies with mocked I/O and plotting.
- Initial pytest capture and direct Matplotlib import failed because the sandbox has no writable temporary/cache directory; corrected approaches required no writes. A fractional fort.14 hypothesis was rejected by the reader and is not reported.
- Git status remained unchanged, including the pre-existing untracked OBC PNG. No shared-data contents were read and no jobs were submitted. Production mesh counts, implied timesteps and FVCOM smoke results were **not rerun**.

## Verdict
VERDICT: FAIL (0 blocker, 0 major, 7 minor, 0 nit)

### Prompt

```markdown
# Review request, round 26: extending a base mesh outward (fvcom-mesh-tools)

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
- the new part is generated with oceanmesh (run by 445 as a subprocess
  stage; the package may import oceanmesh since the relicensing), with
  fixed points/edges and ladders on
  both constrained lines, `cleanup="none"`, a constrained-Delaunay repair,
  flat-element removal; then finishing, coast fit, merge, a repair limited
  to the new part and kept off the open boundary, depths from the recipe's
  source stack, an r-factor limit with base depths held, export and QA;
- the open boundary is designed orthogonal to the coast at both ends, with
  straight legs and filleted corners, spacing never below the CFL floor.

Out of scope: the oceanmesh fork itself; the tide tools (notebooks 449-454,
`tide_models.py`), reviewed separately.

## Previous rounds
Rounds 1-25 and their triage are in docs/extend-tools-review-20261001.md.
The package is GPL-3.0-or-later (e37a433); OCSMesh/Triangle/JIGSAW are
optional private-use backends outside the default environment (c76c0c6;
owner's decision) -- do not re-report their existence, only inconsistencies.

Round 25 (your previous answer; 7 findings) was fixed in 02b4ca2; read it.
Per finding:
- F1 384: convergence on max |h2 - h1| of the interpolated complex
  coefficients (reported as half_window_coefficient_change_m).
- F2 `_indices` refuses masks; base OBC chains checked in
  verify_frozen_base.
- F3 `check_island_holes(mesh, land, n_base_elements)`: required, an
  integer in [0, NE]; base holes only from the base elements' loops.
- F4 masks refused in `_lattice`, compose_sizing bands/floor, band_field,
  fillet, resample.
- F5 limiter: shapes and edge indices checked.
- F6 '=' key handled as '='.
- F7 QA in 447/453 with min_depth_m = min(recipe min_m, base minimum).
Verification after 02b4ca2: full test suite 1242 passed (batch
job 124632); on real data 444, 445 (QA 23/23 with the recipe's depth
floor, status ok, grd bit-identical to rounds 4-24) and 453 (status ok)
passed; check_run accepts the two real FVCOM smoke histories.
Owner decision (2026-10-01), unchanged: meshes are made from the real
depths; the band-floor check (446) and the new-element time-step comparison
(447, 453) REPORT warnings and do not fail the build. Not a defect.

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
| F1 | minor | code read (min() follows a bad base depth) | yes, introduced in round 25 | fixed, 3df171d |
| F2 | minor | test (masked x) | yes | fixed, 3df171d; test |
| F3 | minor | test (base chain 1.9, 2.9) | yes | fixed, 3df171d; test |
| F4 | minor | tests (masked, non-boolean mutable) | yes | fixed, 3df171d; tests |
| F5 | minor | tests (masked depth, NaN bounds) | yes | fixed, 3df171d; tests |
| F6 | minor | tests (masked, NaN rays) | yes | fixed, 3df171d; tests |
| F7 | minor | tests (masked, out-of-range indices) | yes | fixed, 3df171d; tests |
