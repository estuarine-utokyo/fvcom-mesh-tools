# Third adversarial review

**Verdict: not yet fit for unattended acceptance as documented.** Reviewed
`6f9e574` (including `2acacdb` and `f36e9fc`), 2026-09-25. The earlier
failure/retry bugs are repaired, and the ordinary suite passes. However,
`fmesh-check-run` still certifies incomplete history, and concurrent
refinements can bypass the promised overwrite protection. The successful
production chain reported by the owner does not exercise these cases.
Supervised use with isolated directories and independent history checks is
reasonable; these findings do not establish damage to the saved meshes.

Only this report and `docs/hires_walls_tools_review_3_probes.py` were created.
All eight final probes fail at their intended behavioral assertions. They
use tiny synthetic inputs; none runs a mesher, FVCOM, or the scheduler.

**Confirmed findings, ordered by severity**

1. **T1 — Major: restart output can substitute for history.**
   `src/fvcom_mesh_tools/cli/check_run.py:98`, `:109`, `:121`.
   `test_t1_restart_is_not_history` supplies only
   `m2_restart_0001.nc`, one finite record at END_DATE: **ok=True,
   n_records=1**. `test_t1_restart_must_not_fill_history_gap` removes an
   entire middle history stack but supplies those times in restart files:
   **ok=True, n_records=25**. Presence of `zeta` does not identify the stream.
   This is realistic field content: the local FVCOM
   `src/mod_ncdio.F:1589` onward constructs restart output from time, zeta,
   velocity and averaged-velocity objects. Its default restart filename is
   set in `src/mod_input.F:1471`. These sibling sources were read only.
   **Fix:** identify the configured history stream explicitly, validate its
   required variables and coverage independently, and exclude restart and
   average streams from history timestamps. An explicit history pattern or
   case identity is preferable to treating every file with zeta as history.

2. **T2 — Major: missing initial stacks and reversed records pass.**
   `src/fvcom_mesh_tools/cli/check_run.py:109`, `:121`, `:131`.
   Both `test_t2_*` cases declare START_DATE/NC_FIRST_OUT at hour 0,
   END_DATE at hour 24, hourly output, and `NC_OUTPUT_STACK=12`.
   Removing the first 12-record stack returns **ok=True, n_records=13**.
   Reversing the first stack's records also returns **ok=True,
   n_records=25**. Checking only the latest time misses initial loss;
   `sorted(set(stamps))` erases evidence of order and overlap.
   **Fix:** derive expected first output from the configured output schedule
   and startup mode; check leading and trailing coverage and record order
   within each stack. Reconcile legitimate boundary duplicates explicitly
   rather than silently sorting/deduplicating arbitrary output. Ordinary
   three-stack output passed the positive control, so splitting itself is
   not broken.

3. **T3 — Major: valid interval syntax silently disables the gap gate.**
   `src/fvcom_mesh_tools/cli/check_run.py:75`, `:87`, `:125`.
   `test_t3_*` supplies only hours 0 and 24 for hourly output. Both
   `NC_OUT_INTERVAL = 'cycles = 360'` with a 10-second internal step and
   `NC_OUT_INTERVAL = "seconds = 3600."` return **ok=True,
   n_records=2**, without `output_interval_s`. The regex understands
   neither cycles nor double-quoted strings; it treats an unsupported
   declaration as an absent one and skips gap checking. FVCOM's
   `src/mod_set_time.F:1466` and `:1786` explicitly support cycles and
   multiply them by the internal timestep.
   **Fix:** parse the supported namelist syntax and convert cycles using the
   actual integration step. Fail closed on a present but unrecognized or
   invalid interval. Whitespace around `=` already works; positive controls
   also passed for seconds, fractional days, and the checker's hours/minutes
   forms. Those latter controls test the checker, not FVCOM unit support.

4. **T4 — Major: concurrent output reservation is now reproduced.**
   `src/fvcom_mesh_tools/cli/refine_run.py:73`, `:90`;
   `notebooks/420_local_refine.py:90`;
   `jobs/octopus/417_hires_refine.sh:38`.
   `test_t4_*` runs two real wrapper invocations against the same initially
   absent output. A barrier holds the stub generators before their first
   write. **Both return 0**; the second stub overwrites the first stub's
   `report.json` (`writes=[False, True]`). Only the expensive generator is
   replaced: admission through the wrapper is real. The actual generator's
   `mkdir(exist_ok=True)` provides no later exclusion. The shell guards have
   the same check-before-use structure.
   **Fix:** acquire an exclusive reservation before starting any generator,
   enforce it at the common entry used by CLI and jobs, and bind it to the
   attempt. Handle failures and stale reservations explicitly. This is a
   reproduced admission/overwrite race, not a claim that two real mesh jobs
   were run or that a particular corrupted mesh was produced.

5. **T5 — Minor: a new frozen/frozen edge is still called inherited.**
   `src/fvcom_mesh_tools/cli/refine_depths.py:135`, `:142`;
   `src/fvcom_mesh_tools/cli/finish_depths.py:184`;
   `docs/USER_GUIDE.md` §8.
   `test_t5_*` supplies a four-node, two-triangle square with depths
   `[3, 4.5, 6.75, 4.5]`. With original diagonal 1–3, every edge has
   **r <= 0.2**. Changing only connectivity to diagonal 0–2 gives
   **r=0.3846153846**, one over-limit frozen pair, but the shared finisher
   still returns **converged_at_write_precision=True**. The nodes and depths
   remain unchanged. This confirms the numerical exemption; it does not
   demonstrate that DistMesh generated this topology in production.
   The guide accurately describes the movable-end gate, but its next claim
   that such an edge “is the base's own” is false without connectivity
   provenance. The non-ladder generator already performs that provenance
   check (`420_local_refine.py:1574`); the later depth finisher does not.
   **Fix:** distinguish inherited edges from new connections through
   `node_map` and base connectivity, reject new unfixable violations, and
   retain the intended exemption for inherited ones. At minimum, document
   that convergence is not a guarantee against newly introduced fixed-pair
   violations. The probe expresses that stronger acceptance requirement.

**Disposition of the requested earlier findings**

| Finding | Disposition | Evidence and limits |
|---|---|---|
| R1 | Fixed for the reproduced sequential failure/retry cases | Existing prerequisite/recheck tests pass. Invalidation now precedes common.sh. Additional actual-invalidation-block tests abort safely for unset and empty FMESH_RUN_ROOT and remove all relevant old markers before a simulated SIGKILL after partial output. No scheduler was invoked. Concurrent shared-directory attempts remain outside this guarantee. |
| R2 | Fixed | Existing nonzero-MPI regression passes. Additional execution of the actual 412 tail with an initially empty output directory, stub MPI producing valid output, and the real checker returns 0 and writes RUN_OK. Unmatched `rm -f .../*.nc` is safe. Marker publication is conditional on solver status. |
| R3 | Partially fixed | Missing ua/va/Times and nonfinite fields are rejected for recognized history, but stream recognition is wrong: T1. |
| R4 | Partially fixed | The observed gap no longer enlarges tolerance; the old regression passes. T2/T3 still defeat schedule coverage. |
| R5 | Fixed within the documented Tokyo Bay scope | Existing batch-invariance regressions pass. 100 individual cropped searches matched an exhaustive spherical tree; see units/completeness below. |
| R6 | Fixed for its reported CRS claim | §2.5 now explicitly calls the check plausibility-only, names the neighboring-zone limitation and assigns CRS verification to the user. Single-arc restriction remains explicit. |
| F5 | Partially fixed overall | Sequential marker handoff is repaired, including interruption after partial output, but checker false positives can still create success markers for incomplete products. Shared attempts also lack isolation. |
| F6 | Partially fixed | Ordinary failure checks pass; T1–T3 remain. |
| F9 | Fixed within scope | Per-query sphere metric removes batch-latitude dependence; cropped-search controls pass. No antimeridian claim is made. |
| F13 | Fixed | The guide now accurately states both base restrictions and the limit of CRS detection. |

The interruption checks ran the real initial `rm` blocks of 421, 423 and
412, wrote a partial file, then killed that child shell with SIGKILL.
Each returned **-9**, retained the partial file, and left its own/downstream
markers absent. This models the relevant interruption point, not NQSV's
entire termination protocol. The six unset/empty-root cases returned
nonzero before a stub `rm` could be called: `${FMESH_RUN_ROOT:?...}` guards
the whole command expansion. Marker writes are still non-atomic and not
bound to attempt identity; do not interpret this as a concurrency guarantee.

**Geometry, remaining risks, and the plotting change**

- **DEM completeness and units:** `_xyz` returns Earth-centered coordinates
  in metres using the radius implied by `_M_PER_DEG`; cKDTree returns chord
  distance in metres, not an angular distance. A 0.01-degree eastward step
  at 35 N measured **909.258768 m**. Chord and surface distance differ, but
  no materially wrong Tokyo Bay result was reproduced. With a seeded sparse
  51×51 grid over 139.4–140.4 E / 34.8–35.8 N, all **100** separate cropped
  queries matched exhaustive nearest elevations and distances within
  **1e-6 m**. `_pad_needed` uses the most extreme box latitude and converts
  metres back to degree margins; its conservatism here covers the small
  chord/arc difference. `_nearest_finite` eventually uses the whole grid.
  The soundings search uses the same margin and has a large-distance cutoff;
  it is not a global completeness proof. No new in-scope F9 failure found.
- **Frozen/frozen risk:** now a confirmed local numerical/documentation
  limitation (T5), with a synthetic connectivity witness. Whether the full
  hires generator produces an offending new retained-node edge remains
  unconfirmed. No new mesh, QA count or implied dt was measured.
- **Offshore islands with real OSM:** leave documented as an integration
  risk, not a new finding. The corrected filter/island guards remain
  independent of a free base coastline, and the existing helper/block tests
  pass. Neither those tests nor the owner's port run establishes a complete
  offshore real-OSM/oceanmesh run. I did not run that integration or invent
  a failed geometry. No confirmed finding here needs a batch job. A future
  offshore-island acceptance run needs an island-containing recipe; rerunning
  the existing offshore fishery alone is not evidence that an island was
  included.
- **`6f9e574`: no confirmed regression.** `mesh_alpha` is keyword-only and
  defaults to 1, preserving existing callers. The new reader returns
  zero-based OBC indices as `draw_mesh` expects. A tiny artist-level check
  with alpha 0.25 found **two solid edges and one open edge**, with opaque
  boundary collections. Existing plotting tests pass. No figure batch or
  visual assessment of dense M2 maps was run. Selecting the first OBC file
  assumes the staged input contains the intended single case; this was not
  reproduced as a failure.

**What users need to know beyond the guide**

Use a unique output/run directory per attempt and never start simultaneous
writers to it. `RUN_OK` currently does not prove complete, ordered history:
check the configured first output, all history stacks and their cadence
independently, excluding restart/average files. Use the checker's supported
single-quoted cadence form until T3 is fixed. A rerun through 412 deletes
**all `output/*.nc`**, so preserve needed history/restart products elsewhere
first. Frozen-pair violations are exempt even when the patch created their
connectivity; inspect that report before treating a converged product as
globally r-limited. Real offshore-island integration remains unverified.

**Verification**

Commands run from the repository root in `oceanmesh-bench`:

```bash
. /octfs/work/G16445/v61021/miniforge3/etc/profile.d/conda.sh
conda activate oceanmesh-bench
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider
ruff check --no-cache .
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider \
  docs/hires_walls_tools_review_3_probes.py --tb=short
ruff check --no-cache --no-force-exclude docs/hires_walls_tools_review_3_probes.py
```

Suite: **825 passed, seven rasterio deprecation warnings, 58.86 s**.
Review probes: **eight intentional failures**, covering T1–T5; no xfails.
Repository and explicit probe lint passed. Additional checks described above
ran through an inline `PYTHONDONTWRITEBYTECODE=1 python -` script with
temporary files. Matplotlib emitted cache/fontconfig permission warnings,
then its artist assertions passed; no image was saved. An initial candidate
using a D exponent was discarded after source inspection showed FVCOM's
interval tokenizer does not accept it; the final syntax probe uses ordinary
double quotes instead. An exploratory search for `mod_value.F` failed
because the parser actually lives in `mod_utils.F`; that source was then read.

No existing file was modified, dependency added, batch job submitted, mesh
generated, full depth product finished, or FVCOM integration run. No qsub
or qstat was attempted. The owner-reported GitHub CI and jobs 115839–115845
were not rerun here. The outstanding work is T1–T5 and full real-OSM offshore
integration; successful routine runs do not close those acceptance gaps.
