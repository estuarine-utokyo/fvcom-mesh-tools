# Review: extending a base mesh outward (2026-10-01)

Reviewer gpt-6-astra (Codex CLI), read-only; process in
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
