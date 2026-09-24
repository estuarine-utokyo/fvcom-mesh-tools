# Fourth adversarial review

**Verdict: the original T1–T5 counterexamples are repaired, but fix U1 before
handing users the documented general-purpose run-acceptance gate.** One major
false-positive remains: commented namelist settings can override the actual
output schedule. Three minor findings concern a short-run regression and
unusual provenance/reservation states. The findings are diminishing in scope,
but are not all nits or documentation. After these bounded fixes and their
regressions, I would stop broad review iteration and move to supervised user
validation; this review does not justify another mesh-algorithm redesign.

Reviewed `6ca5e80`, 2026-09-25. Only this report and
`docs/hires_walls_tools_review_4_probes.py` were created. Existing files were
not edited. All four new probes fail at their intended behavioral checks.

**Confirmed findings, prioritised**

1. **U1 — Major: comments can make incomplete history pass.**
   `src/fvcom_mesh_tools/cli/check_run.py:52`, `:118`, `:174`.
   `test_u1_commented_interval_cannot_hide_missing_records` puts
   `! NC_OUT_INTERVAL = 'days = 1',` before an active hourly declaration.
   History contains only hours 0 and 24. The checker returns **ok=True,
   n_records=2, output_interval_s=86400**, although the active interval is
   3,600 seconds. Without the comment the same gap is rejected. `_nml_value`
   searches raw text and takes its first match, including Fortran comments;
   the same parser handles the date and integration settings. This is an
   acceptance defect, not merely failure to understand malformed input.
   **Suggested fix:** parse the relevant namelist groups and active values,
   at least handling comments outside quoted strings before matching. Do not
   strip `!` inside strings. Test commented dates/intervals and quoted text.
   The older regex also searched comments: this is a newly confirmed gap,
   not a claim that comment handling first broke in this commit.

2. **U2 — Minor: a valid one-record history crashes the checker.**
   `src/fvcom_mesh_tools/cli/check_run.py:171`, `:174`.
   `test_u2_valid_single_record_history_returns_a_verdict` sets first output
   equal to END_DATE, 24 hours after START_DATE, with an hourly interval and
   one finite record at that time. It raises **ValueError: zero-size array
   to reduction operation maximum which has no identity**, rather than
   returning success. `np.diff` is empty and `steps.max()` is unguarded.
   This is a regression in the rewritten history loop. It fails closed,
   so does not create a false success marker.
   **Suggested fix:** check cadence only when at least two timestamps exist;
   still apply first/last coverage checks to one timestamp. Cover both a
   complete single-record schedule and an incomplete one.

3. **U3 — Minor: the base pathname is treated as mesh identity.**
   `src/fvcom_mesh_tools/cli/refine_depths.py:58`, `:65`, `:71`, `:175`;
   `notebooks/420_local_refine.py:283`.
   `test_u3_replaced_base_cannot_excuse_a_new_frozen_edge` uses four nodes,
   two triangles, depths `[3, 4.5, 6.75, 4.5]` and the T5 diagonal flip.
   With the actual original base, the new frozen diagonal has
   **r=0.3846153846**, **n_over_frozen_pair_new=1**, and convergence is false.
   Replacing only the external `base.14` at the recorded path with the flipped
   topology changes the verdict to **converged=True**, **new=0, inherited=1**.
   The refinement and node map did not change. Thus the T5 fix works only
   while that external file retains its identity. No real saved mesh was
   shown to suffer this replacement.
   **Suggested fix:** persist the base connectivity/mapped inherited edges,
   bound to the refinement, or record and verify a base fingerprint before
   granting exemptions. Treat legacy reports without identity evidence
   conservatively. At minimum require immutable base files explicitly in
   the guide. This is a reproducible provenance limitation with a specific
   external-mutation prerequisite, hence minor rather than a routine-path
   major defect.

4. **U4 — Minor: a stale delegated token can fall back to a new claim.**
   `notebooks/420_local_refine.py:95`, `:98`, `:100`.
   `test_u4_stale_delegated_token_must_not_reclaim_nonempty_output` executes
   the actual reservation block with `LR_RESERVATION` set, `.reserved`
   absent, and an existing `report.json`. The block **admits the caller and
   creates a new reservation**, instead of rejecting the missing delegated
   claim. It does not check that the directory is empty. This establishes
   admission to an old output, not an executed meshing overwrite. The
   wrapper and job 417 separately reject this nonempty directory, so this
   finding is limited to direct driver invocation or an altered handoff.
   **Suggested fix:** when a token is supplied, require an existing exact
   match and fail otherwise; only token-free standalone entry may reserve
   a fresh directory, and it should enforce nonempty-output refusal too.

**Disposition of T1–T5 and requested edge cases**

| Finding | Disposition | Evidence |
|---|---|---|
| T1 | Fixed | Both restart-substitution regressions pass. Case `m2.a+[b]` works with inference from its namelist name and with explicit `casename`; the regex escapes the name. |
| T2 | Fixed for the reported loss/order cases; short-stream regression U2 remains | Missing initial stacks and reversed records are rejected by the existing regressions. Delayed first output and hot-start controls below pass. |
| T3 | Partially fixed overall | Both quoted cadence forms and cycles regressions pass; an unreadable quoted interval fails closed. U1 still reads the wrong effective schedule from a valid commented namelist. |
| T4 | Fixed for competing fresh wrapper/job claims; partially hardened direct entry | Existing concurrent-wrapper probe passes. Matching token is admitted; stale token with an existing different reservation is refused. Missing reservation takes the U4 fallback. |
| T5 | Fixed with the original readable base; partial provenance guarantee | New frozen diagonal fails. Both native grd and fort.14 produce the expected five mapped edges. A moved-away base returns None and grants no exemptions. A replaced readable base produces U3. |

The full staged namelist read was
`/octfs/work/G16445/v61021/scratch/m2_20260924_181844/base/m2_run.nml`.
Its actual START_DATE and NC_FIRST_OUT are both 2021-01-01; it does **not**
already have delayed first output. In a temporary copy, changing NC_FIRST_OUT
to January 20 while retaining START_DATE January 1 and END_DATE January 21
passes with **49 synthetic half-hour records**, January 20–21. A second copy
with `STARTUP_TYPE='hotstart'`, START_DATE January 20 and a startup filename
also passes with those 49 records. These are checker controls using a real
namelist structure, not FVCOM integrations or validation of a restart file.
Read-only FVCOM source inspection (`FVCOM/src/mod_set_time.F:86`, `:224`)
confirms ordinary hotstart uses START_DATE and the declared NC_FIRST_OUT;
there is no basis here to demand pre-restart history. Crashrestart/forecast
modes are not established by this control.

A crashed reservation is recoverable as documented: §6, lines 308–312,
says nonempty outputs must be moved aside and `.reserved` persists. A
temporary wrapper control returned **2** for a directory containing only
the abandoned reservation; after moving the **whole directory**, its fresh
dry run returned **0**. There is no automatic expiry, and none is promised.
Confirm the old writer is stopped before recovery; preserve its directory
for diagnosis. Removing just the claim file is not the documented recovery.

**Final handoff assessment**

- **(a) Defects to fix:** U1 is the remaining blocker to trusting the guide's
  run-acceptance claim for user-edited namelists. U2 is a small normal-input
  regression worth fixing in the same change. U3/U4 are bounded hardening
  fixes; immutable bases and the documented wrapper/job entry points avoid
  their demonstrated triggers, but unrestricted provenance/direct-entry
  claims require the fixes.
- **(b) Limits to document:** until U3 is fixed, archive the exact base and
  preserve its contents and numbering; a pathname alone is not provenance.
  A moved base loses inherited-edge exemptions (already described in §8).
  Reservation recovery is already adequately described; an explicit
  “move the whole stopped attempt, not just `.reserved`” sentence would
  improve it. Keep the documented EPSG:32654, single-OBC-arc and Tokyo Bay
  staging limits. Do not advertise the hot-start checker control as a full
  hot-start integration test.
- **(c) Risks to accept:** full real-OSM offshore-island integration remains
  unverified here, as in review 3. Shared downstream run directories still
  require one writer per attempt; refinement reservation is not a lock for
  the entire FVCOM workflow. The owner-reported successful jobs and saved
  meshes are useful routine-path evidence, not proof of all geometry or
  interruption cases. No additional reproduced meshing/DEM/wall regression
  was found; the earlier regression files pass in the full suite.

**Verification**

Commands from the repository root, in `oceanmesh-bench` (caches/bytecode
disabled to respect the read-only review):

```bash
. /octfs/work/G16445/v61021/miniforge3/etc/profile.d/conda.sh
conda activate oceanmesh-bench
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider
ruff check --no-cache .
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider \
  docs/hires_walls_tools_review_4_probes.py --tb=short
ruff check --no-cache --no-force-exclude docs/hires_walls_tools_review_4_probes.py
```

Results: **834 passed**, seven rasterio deprecation warnings, **103.60 s**;
repository Ruff passed; explicit probe Ruff passed (docs is normally
excluded). New probes: **four intentional failures in 0.97 s**, no xfails.
The additional controls above ran with inline Python, temporary fixtures,
and `runpy.run_path` to reuse the new probe helpers. Existing tracked-file
diff remained empty.

A read-only `check_run` of the existing
`/octfs/work/G16445/v61021/scratch/m2_20260924_181844/smoke/base` also passed:
**97 records**, 1,800-second interval, January 1–3, 2021, reaching END_DATE.
No output or marker there was rewritten.

No dependencies added, mesh generated, production depth product finished,
FVCOM solver run, figure batch run, or scheduler invoked. Mesh QA gates,
production node/element counts and implied dt were **not remeasured**.
All new findings reproduce without a batch job, so no qsub command is needed
to close them. This review does not rerun the owner's CI or job 417.
