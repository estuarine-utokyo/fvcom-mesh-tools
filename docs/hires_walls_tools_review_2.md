**Second adversarial review — `77234ad` and `18d7c85`, 2026-09-24.**
The work is substantially improved, but **not yet fit for the documented
unattended acceptance chain**: failed attempts can leave or even create success
markers, and the run checker accepts missing required fields and seriously
incomplete output. Geometry fixes passed the additional small checks below;
the saved products do not establish these failure paths are safe. CI now
credibly represents the package tests and selected lint scope: I reproduced
811 local passes and 787 passes/19 skips without the laboratory imports/data
environment. Its setup is consistent with the action documentation, but a
fresh GitHub environment solve was not run here. Green CI still does not
establish workflow acceptance or full generator integration.

**Confirmed findings, ordered by severity**

- **R1 — Major: stale success markers survive failed attempts (F5 residual).**
  `jobs/octopus/421_finish_and_run.sh:23`,
  `jobs/octopus/423_m2_smoke.sh:28`,
  `jobs/octopus/412_m2_run.sh:22`,
  `src/fvcom_mesh_tools/cli/check_run.py:129`.
  Prerequisite checks precede removal of the stage's old marker. The three
  `test_r1_failed_prerequisite_*` cases execute those actual shell statements:
  each exits **2**, retaining respectively `STAGED`, `SMOKE_OK`, or `RUN_OK`.
  The checker probe first writes a valid marker, replaces the log with an
  early STOP, then gets exit **1** with the old marker still present.
  A concrete reused-directory chain is: 421 removes `STAGED` and fails
  finishing; 423 refuses the missing `STAGED` but retains old `SMOKE_OK`;
  412's prerequisite therefore passes. Likewise 413 only tests existence of
  both `RUN_OK` files (`413_m2_analysis.sh:22`). Fresh workflow directories
  reduce exposure, but the individually documented jobs can reuse directories.
  **Fix:** invalidate the attempt's own and downstream markers before any
  failure path; bind markers and prerequisite checks to an attempt ID and
  input identity. Publish markers atomically. Moving `rm` alone does not
  cover an environment-setup failure before it, or concurrent attempts.

- **R2 — Major: 412 creates `RUN_OK` despite a failed MPI exit (F5 regression).**
  `jobs/octopus/412_m2_run.sh:46`, `:57`.
  Solver failure sets `status=1`, but the checker still receives `--marker`
  unconditionally. `test_r2_*` executes the real script tail with mocked
  environment commands and an MPI stub returning **7**, a TADA log, and a
  tiny complete NetCDF. The real checker writes `RUN_OK`; the job then exits
  **1**. This also models old complete output surviving a failed rerun:
  412 does not clean its output directory. NQSV releases 413, whose marker
  check accepts that failed job. No real solver or scheduler was invoked.
  **Fix:** run diagnostics regardless, but publish `RUN_OK` only after both
  the process status and checker pass; isolate each attempt's output.
  423 correctly includes the solver exit in its final `fail` decision.

- **R3 — Major: required output fields are optional in the checker (F6 residual).**
  `src/fvcom_mesh_tools/cli/check_run.py:90`.
  `test_r3_*` supplies valid Times through END_DATE and a TADA log, with
  respectively no fields, only `zeta`, or `zeta` plus `ua`. All three return
  **ok=True**. `if var in ds.variables` skips the missing-field failure.
  The existing successful fixture in `tests/test_check_run.py:22` itself
  creates only `zeta`, so it endorses the weaker behavior. This contradicts
  USER_GUIDE §9's requirement for finite `zeta`, `ua`, and `va` in every record.
  **Fix:** identify the required history-output stream and require all three
  fields with nonempty, compatible time dimensions, then check every record.
  Keep ancillary NetCDF files separate from that required stream.

- **R4 — Major: a missing-output gap enlarges the completion tolerance (F6 residual).**
  `src/fvcom_mesh_tools/cli/check_run.py:87`, `:101`.
  `test_r4_*` specifies hourly output and END_DATE **2020-01-03**, but supplies
  only **Jan 1 and Jan 2** records, with all three finite fields and TADA.
  The checker returns **ok=True**, accepting an output end 24 requested
  intervals short. Its tolerance comes from the last observed gap, not the
  configured cadence; a damaged/incomplete series can thus relax its own gate.
  Across files, the final file with two records also determines that tolerance,
  independently of which file supplied the latest timestamp.
  **Fix:** derive a bounded tolerance from the namelist's output schedule,
  verify time ordering/coverage for the required history stream, and reject
  gaps rather than using them to justify missing completion.

- **R5 — Minor: batch latitude still changes extrapolated depths and distances (F9 residual).**
  `src/fvcom_mesh_tools/dem/tokyo_bay.py:192`, `:226`.
  Search-box completeness is fixed, but the metric still uses the batch's
  mean latitude. `test_r5_depth_*` uses two finite cells near `(139,35)`:
  querying that location alone returns **20 m**; adding a query at latitude
  **35.6** changes the same location to **10 m**. Both are Tokyo Bay-scale
  latitudes; no global/polar counterexample is needed. The sounding probe
  likewise changes **909.258769 m** to **905.912725 m** for the same point and
  same sounding. These are synthetic counterexamples, not demonstrated
  corruption of the saved meshes.
  **Fix:** use a fixed product projection/reference latitude, or a distance
  metric defined per query independently of the other queries. Cover both
  public `sample` and `sounding_distance`, not just `_nearest_finite` with
  a fixed `lat0`.

- **R6 — Minor: the guide overstates CRS detection (F13 residual).**
  `docs/USER_GUIDE.md:133`, `notebooks/420_local_refine.py:125`.
  The guide says the two limits are checked and refused. The CRS test can
  only establish numerical plausibility after assuming EPSG:32654, not the
  source CRS. Executing the actual preflight with a tiny triangle around
  `(135,34.5)` projected into **EPSG:32653** accepts it and interprets its
  centre as **141.0005 degrees east**. Lon/lat input and zero/two open arcs
  were correctly rejected in the same check. This does not request lifting
  either owner-imposed restriction.
  **Fix:** call this a plausibility check in the guide and make the user's
  responsibility to establish the base CRS explicit; if automatic refusal
  of other CRSs is required, require trustworthy CRS metadata/declaration.
  Naked coordinate arrays cannot identify their UTM zone. Reproduced with
  an inline AST extraction and pyproj; no persistent failing probe asserts
  an impossible inference from the coordinates alone.

**Disposition of every original finding**

| Finding | Disposition | Evidence / boundary of the fix |
|---|---|---|
| F1 | Fixed | `patch.py:1223` rebuilds adjacency after each cut and remaps only live points. Existing two-corner test passes. Additional concave quadrilateral `[(0,0),(100,50),(0,100),(30,50)]` blunts three consecutive acute corners, retains a valid rim and valid indices; repeating with the fourth point frozen retains its exact coordinate and ID. |
| F2 | Fixed | `patch.py:1055` uses exact closing-point equality. The translated 100 x 30 m ring regression passes. |
| F3 | Fixed | `dem/tokyo_bay.py:91` fills the masked array before ndarray conversion. Real finite-fill NetCDF regression passes. |
| F4 | Fixed | `420_local_refine.py:600` records `acute_corners_kept` and suppresses assignment of the proposed cuts in preserve mode. The actual-block geometric-identity test passes. |
| F5 | Partially fixed | Fresh successful stages now publish markers and downstream stages check them. R1/R2 defeat their meaning on failure/reuse. The structural marker tests do not exercise these paths. |
| F6 | Partially fixed | Missing TADA, fatal strings, unreadable NetCDF, nonfinite existing fields and substantially short regular output are rejected by the passing checker tests. R3/R4 remain. |
| F7 | Fixed | `cli/finish_depths.py:190` rejects nonconvergence before writing; diagnostic override still returns 3. `cli/refine_depths.py:132` evaluates the final capped field. Existing infeasible frozen-depth probe passes. Job 421 invokes the finisher as an ordinary command under `set -e`, so exit 3 stops staging. Stale downstream markers are R1, not an exit-code failure in the finisher. |
| F8 | Fixed | `cli/finish_depths.py:117` requires finite `0 <= r < 1`, and positive rounds when enabled; public equal smoother validates inputs. All four invalid-r CLI regressions pass. |
| F9 | Partially fixed | Widening proves completeness for the chosen distance metric in both searches. R5 demonstrates that the public calls still change that metric with batching. |
| F10 | Fixed for the reported overwrite case | `cli/refine_run.py:73`, job 417 and workflow reject nonempty directories, including pre-report failures. Regression passes. Atomic reservation of an initially empty directory remains unimplemented; no concurrent overwrite was reproduced here. |
| F11 | Fixed | `420_local_refine.py:380`, `:398`, `:562`, `:623` decouple filtering/islands/walls from `shore`. Full filter-block checks described below pass; existing island insertion regression passes. |
| F12 | Fixed for the delivered mesh | `walls.py:623`, `:644`, `cli/meshqa.py:110` persist, validate and load explicit pairs. Existing CLI duplicate-exemption/opt-out/hash tests pass; copied/renumbered checks below pass. |
| F13 | Partially fixed | Restrictions are documented and the single-arc guard runs before meshing. CRS plausibility rejects obvious lon/lat inputs, but the stronger promise to reject other CRSs is R6. |
| F14 | Fixed for the reported path issue | `refine_workflow.sh:45` rejects commas in recipe/output/run-root/view values before submission; the comma-recipe regression passes. It does not validate literally every transported value: numeric settings and `LR_SEEDS` still rely on their documented numeric/colon syntax. |

For F11 I executed the complete filter block, including initialization of
`_shl` and `_land_filtered`, with default preserve, default resample, and
hires preserve. All short-circuit without requiring `_keep`, `shore`, or
land. Offshore resolve with no free coastline, a 1 km square water footprint
and a 200 x 200 m island sets `_land_filtered=True`, keeps **40,000 m²** of
land, and leaves `shore=None`. The real filtering/wall helpers ran; only the
shapefile writer and oceanmesh Shoreline constructor were stubbed. This is
stronger than the existing guard-name test, but is not a full generator run.

For F12 a byte-identical renamed copy validates with `--wall-pairs` pointing
to its original sidecar. A genuinely renumbered six-node mesh rejects the
old sidecar with CLI exit **2**. `--no-wall-pairs` bypasses it and reports
ordinary QA failures, including duplicates; mapping the pair through the
inverse permutation and writing a fresh sidecar validates again. A copy
without its sidecar loses exemptions; a renamed copy needs a renamed sidecar
or `--wall-pairs`. The native-case renumber CLI does not propagate this
fort.14 metadata automatically. Hash refusal is the correct behavior for an
old node map, not a regression that should be fixed by relaxing the hash.

**CI and documentation assessment**

The `environment-file`, matching environment name, extra `pre-commit` spec
and `bash -el {0}` shell agree with the
[setup-micromamba v2 instructions](https://github.com/mamba-org/setup-micromamba/blob/v2/README.md).
The YAML is a normal environment file, not a lockfile, so additional package
specs are supported. It provides the scientific stack and setuptools/wheel
needed by the no-deps/no-build-isolation editable install. I found no concrete
CI setup defect. Fresh conda-forge resolution, remote hook installation and
the actual GitHub runner remain unverified; the unpinned environment also
means a future solve need not match this machine.

Installed Ruff is **0.16.7**; the cached pre-commit executable is **0.6.9**.
Both pass `ruff check --no-cache .`. The newer Ruff also passes the explicit
`notebooks/42*.py notebooks/43*.py` invocation. The pinned
[hook definition](https://github.com/astral-sh/ruff-pre-commit/blob/v0.6.9/.pre-commit-hooks.yaml)
uses `--force-exclude`, so pre-commit does not accidentally reintroduce all
historical notebooks through explicit filename arguments. The dedicated CI
lint invocation does include the maintained scripts. Version drift is a
maintenance risk, not a demonstrated current incompatibility. I did not run
mutating pre-commit hooks because this review must leave existing files intact.

Excluding historical research scripts and intentionally failing review
probes is reasonable, as is limiting default pytest discovery to `tests`.
It should not be read as full workflow coverage: the maintained generator
is linted but not executed end-to-end by GitHub; jobs and the active
`414_refine_m2_prep.py`/`384_m2_analysis.py` path are outside that notebook
lint selection. The oceanmesh-dependent behavior is explicitly skipped.
Current tests check marker strings, and isolate the F11 island block; they
do not establish failure/retry semantics or complete driver integration.
CI lists **13 of 19** console scripts despite its "every" comment. I checked
all **19** parser entry points with oceanmesh/xcoast imports blocked and
DATA_DIR unset: all passed. The six omissions are a coverage gap, not a
currently broken CLI.

USER_GUIDE now correctly documents preserve behavior, nonconvergence and
the diagnostic override, the required base CRS/one arc, nonempty-output
refusal and persisted wall pairs. Remaining acceptance claims in §3/§6/§8
and §9 are too strong because of R1-R4; §2.5 needs R6's qualification.
Also clarify that §8's r-factor gate covers edges with a movable endpoint:
`cli/refine_depths.py:142` excludes frozen/frozen violations
from convergence. Document the allowed output-time tolerance and the
copy/renumber procedure for sidecars. Rechecking raw mesh QA reproduces its
wall exemption, not necessarily exit 0: inherited angle and raw-depth
failures are still reported. The production-depth bit-for-bit reproduction
claim in §8 was not rerun in this review.

**Verification and remaining limits**

Commands ran in `oceanmesh-bench`; bytecode and pytest/Ruff caches were
disabled to avoid changes to existing repository files:

```bash
. /octfs/work/G16445/v61021/miniforge3/etc/profile.d/conda.sh
conda activate oceanmesh-bench
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider
ruff check --no-cache .
ruff check --no-cache notebooks/42*.py notebooks/43*.py
/octfs/home/v61021/.cache/pre-commit/repoixdw2wm9/py_env-python3.12/bin/ruff check --no-cache .
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider docs/hires_walls_tools_review_2_probes.py --tb=short
ruff check --no-cache docs/hires_walls_tools_review_2_probes.py
for s in jobs/octopus/{412_m2_run,413_m2_analysis,417_hires_refine,421_finish_and_run,423_m2_smoke,refine_workflow}.sh; do
    bash -n "$s" || exit
done
```

Results: **811 passed**, seven rasterio deprecation warnings, **49.52 s**;
all listed lint/syntax checks passed. New probes: **11 failed**, all at their
intended behavioral assertions (R1-R5), approximately **4 s**. These are
regression expectations, not xfails. The initial ten-probe run also failed
as intended; the output-cadence probe brought the total to eleven.

CI emulation used `pytest.main(['-q', '-p', 'no:cacheprovider'])` in a fresh
Python process after removing DATA_DIR and installing a MetaPathFinder that
raises ModuleNotFoundError for oceanmesh/xcoast and their submodules:
**787 passed, 19 skipped**, seven warnings, **38.48 s**. This emulates missing
imports in that process, not a fresh conda solve or every child process's
environment. Separate parser checks used the same import block and invoked
each `project.scripts` target with `sys.argv=[name, '--help']`.

Read-only checks of the current saved `.14`, report, marker and sidecar:

| `outputs/refine_...` | Nodes / elements | Saved gates passed | Introduced | Wall pairs | ACCEPTED |
|---|---:|---:|---:|---:|---|
| kimitsu_port_hires | 6,018 / 10,862 | 20/22 | 0 | 17 | yes |
| futtsu_coast_hires | 6,156 / 11,410 | 20/22 | 0 | 1 | yes |
| futtsu_nori_hires | 4,458 / 8,120 | 20/22 | 0 | 2 | yes |
| futtsu_nori | 4,413 / 8,039 | 21/22 | 0 | 0 | yes |

Counts come from reading the meshes; QA counts come from saved reports.
Every available sidecar's SHA-256/node-count check passed. No mesh, depth
finish, implied-dt calculation, FVCOM integration, figure batch, qsub or
qstat was run. No new dependency was installed. The only repository files
created are this report and `docs/hires_walls_tools_review_2_probes.py`.
No confirmed finding needs a batch job to reproduce.

Unconfirmed risks remain separate from the findings: concurrent output
reservation (F10), newly created frozen/frozen edges inheriting the
finisher's exemption, and full real-OSM/oceanmesh behavior for offshore
islands or land filtered entirely away. No new failing mesh established any
of these. The owner's successful production runs are valuable evidence for
the ordinary path, but do not cover the reproduced acceptance failures.
