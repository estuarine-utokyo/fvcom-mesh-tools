# Hires, walls, depth finishing and tooling review

Reviewed `8b07e32^..86ac83d` on 2026-09-24. **Not fit for users as documented
yet.** The four saved meshes preserve the frozen coordinates/depths and have
zero lone-corner nodes, but small counterexamples expose coastline corruption,
ignored land, and workflow stages that continue after failed acceptance.
No source, tests, notebooks, jobs or existing documents were changed. The only
new files are this review and `hires_walls_tools_review_probes.py`.

**Confirmed major findings**

- **F1 — Adjacent blunted corners leave deleted node references.**
  `src/fvcom_mesh_tools/patch.py:1263`, `:1271`, `:1281`.
  Each cut uses the original adjacency. A later cut can remove an endpoint
  used by an earlier cut's newly queued edges; the conflict test checks only
  the removed sets. Probe `test_f1_*` blunts the two 45-degree corners of
  `[(0,0),(100,0),(50,50)]` at h=10 and returns edges `[-1,1]` and `[-1,3]`.
  NumPy then interprets these as the last point; `hole_polygon` constructs a
  self-touching rim instead of the intended two cuts. This is before meshing,
  not a DistMesh seed failure.
  **Fix:** apply cuts against updated topology, or reserve all cut endpoints
  and reject conflicting edits; assert valid indices and a simple rim before
  passing constraints to the mesher.

- **F2 — UTM coordinates cause the island walk to erase a real corner.**
  `src/fvcom_mesh_tools/patch.py:1052`.
  Default `np.allclose` uses relative tolerance: near northing 3,909,000,
  points 30 m apart can count as the repeated closing point. `_ring_at_size`
  already removes the actual closing point before calling this routine.
  Probe `test_f2_*` translates the same 100 x 30 m rectangle from the origin
  to `(392000,3909000)` at h=10. Its returned area changes from **3,000 to
  1,500 m²**, despite no width rule justifying that loss.
  **Fix:** use exact equality for a repeated source vertex, or an explicitly
  small absolute tolerance with `rtol=0`. Test translated and rotated rings.

- **F3 — NetCDF masks are discarded before the depth ladder sees them.**
  `src/fvcom_mesh_tools/dem/tokyo_bay.py:88`.
  `np.asarray(masked_array)` removes the mask; the subsequent `np.ma.filled`
  cannot recover it. Probe `test_f3_*` reads a real tiny NetCDF with finite
  `_FillValue=-9999`: the returned elevation remains -9999, so sampling would
  treat missing survey data as a 9,999 m depth and suppress fallback.
  The installed grid30 product declares finite `_FillValue=9.96921e36`;
  M7001 and Kanto declare NaN. The three grid30 rows inspected had no masked
  cells, so this does **not** establish corruption in the saved meshes.
  **Fix:** fill the masked array before converting it to an ndarray; cover
  finite fill values, `missing_value`, interpolation and fallback in tests.

- **F4 — `hires.coastline: preserve` still changes the coastline.**
  `notebooks/420_local_refine.py:556`.
  Blunting is guarded by `HIRES is not None`, not by `coastline == "resolve"`.
  Probe `test_f4_*` executes this actual driver block with `preserve` and an
  acute free coast: it removes **14.3674 m²** of water. This contradicts the
  independent preserve/resolve decision and the geometric identity promised
  in `docs/refine_coast_and_bathy_design.md` §3.
  **Fix:** restrict geometric blunting to resolve mode. Preserve mode must
  repair topology without moving its coastline, or reject an unsatisfiable
  quality target explicitly.

- **F5 — The workflow runs downstream stages on rejected products.**
  `jobs/octopus/refine_workflow.sh:48`, `:61`;
  `jobs/octopus/421_finish_and_run.sh:21`;
  `notebooks/420_local_refine.py:1887`, `:1903`.
  Precise path: a candidate can preserve the frozen contract and meet
  resolution while failing QA; the generator writes `node_map.npy`, the
  FVCOM case and report, then exits nonzero for introduced violations.
  NQSV's after-termination dependency still releases job 421, which checks
  neither the report nor a success marker. The finisher accepts the available
  files and staging proceeds. Similarly, optional 20-day runs follow a failed
  smoke job: their namelists were staged before smoke ran, so mere existence
  cannot certify smoke success. No scheduler submission was needed to
  establish these paths.
  **Fix:** write an atomic, run-specific success marker only after each
  stage's acceptance checks; require it in every dependent stage. Keep failed
  diagnostic artefacts distinct from accepted products.

- **F6 — Smoke can pass without completion or readable, finite output.**
  `jobs/octopus/423_m2_smoke.sh:60–77` versus `docs/USER_GUIDE.md` §9.
  Probe `test_f6_*` executes the unmodified validation loop with a stub solver
  that exits zero after `STOP: integration ended early`, plus a 26-byte
  placeholder `.nc`. Both cases pass and the script exits **0**. It checks
  selected error strings and file existence/size, but never checks `TADA`,
  reads NetCDF variables, counts records or checks the final simulation time.
  **Fix:** require completion, the requested end time/record coverage and
  finite required fields. Make this validation determine the success marker
  from F5. The probe invokes no real MPI, FVCOM or scheduler.

- **F7 — An unconverged depth product is exported as successful.**
  `src/fvcom_mesh_tools/cli/finish_depths.py:170–191`.
  Probe `test_f7_*` covers two cases: zero sweeps leave r=**0.941748**; default
  2,000 sweeps cannot reconcile a frozen 1 m node with new nodes floored at
  3 m (r=**0.5**). Both print NOT CONVERGED, write an `rfac0p2` product and
  return **0**. Job 421 therefore stages it as run-ready. The diagnostic is
  honest, but its process status and downstream behavior defeat it.
  **Fix:** return nonzero and withhold the accepted case on nonconvergence;
  permit diagnostic output only through an explicit override. Check the
  final capped field: patch-only capping can create new violations against
  uncapped frozen nodes even after the smoother itself converged.

- **F11 — Resolve ignores new islands when the hole has no base coastline.**
  `notebooks/420_local_refine.py:335`, `:355`, `:532`.
  An offshore hole with no free base coastline leaves `shore=None` and
  `_shl=None`. Consequently neither filtering/wall extraction nor island
  insertion runs, even when OSM contains resolvable land inside that hole.
  Probe `test_f11_*` shows that the helper accepts a 200 x 200 m island inside
  a 1 km water square at h=30, but the actual driver's island block skips it
  entirely. There is no skipped-feature report. This contradicts the §14
  fix for land the base never represented.
  **Fix:** discover and filter land from the hole footprint independently of
  whether a base coastline stretch exists; then add its islands and walls.

**Confirmed minor findings**

- **F8 — Invalid r-factors silently disable or trivialize smoothing.**
  `src/fvcom_mesh_tools/cli/finish_depths.py:142`;
  `src/fvcom_mesh_tools/refine.py:665`.
  All four `test_f8_*` cases return 0: `nan` and `-0.2` become “no smoothing”,
  while `1.2` and `inf` claim convergence without changing a steep field.
  **Fix:** require finite `0 <= r < 1`, with exactly zero meaning disabled;
  validate rounds and validate the public equal smoother consistently with
  `limit_rfactor`.

- **F9 — “Nearest finite cell” depends on the other query points.**
  `src/fvcom_mesh_tools/dem/tokyo_bay.py:135–139`.
  The first search box containing *any* finite cell terminates the search,
  without proving that a nearer cell lies nowhere outside it. In
  `test_f9_*`, querying `(0,0)` alone chooses elevation -10 at **2,982.58 m**;
  adding a second query makes the first choose its actual nearest elevation
  -20 at **2,331 m**. The same location thus receives different extrapolated
  depths depending on batching.
  **Fix:** expand until each nearest-distance circle is contained in the
  search region, or query a cached complete spatial index. Apply the same
  completeness reasoning to `sounding_distance`'s cropped search.

- **F10 — Refusal to overwrite only protects directories with a report.**
  `src/fvcom_mesh_tools/cli/refine_run.py:68`;
  `jobs/octopus/417_hires_refine.sh:43`;
  `jobs/octopus/refine_workflow.sh:32`.
  `test_f10_*` supplies a nonempty output directory with constraints but no
  report; the wrapper invokes the generator and returns 0. Pre-seed failures
  can leave exactly this state, since constraints/shapefiles precede the
  first report. The guide's final “refuses to overwrite” assurance is broader.
  **Fix:** refuse existing nonempty outputs, and reserve the output directory
  atomically to avoid concurrent-run races.

- **F12 — Standalone QA cannot reproduce the delivered wall exemption.**
  `notebooks/420_local_refine.py:1806`, `:1816`;
  `src/fvcom_mesh_tools/qa.py:665`;
  `src/fvcom_mesh_tools/cli/meshqa.py:99`.
  The saved Kimitsu QA excuses its wall copies, but running ordinary
  `run_qa` on the delivered `.14` reports **17 duplicate pairs**. `copy_of`
  is stripped from the persisted report, and the documented QA CLI has no
  mechanism to load the declared pairs. Rechecking the delivered mesh thus
  changes its verdict. **Fix:** persist explicit, mesh-bound wall metadata
  and let the CLI validate/load it; never infer permission from coincident
  coordinates alone.

- **F13 — The new guide overstates supported base meshes.**
  `docs/USER_GUIDE.md` §2.1/§4 versus
  `notebooks/420_local_refine.py:79`, `:1672`.
  “Any FVCOM grd/dep/obc works” omits two code requirements: coordinates are
  assumed to be **EPSG:32654**, and serialization requires **exactly one**
  open-boundary arc. A different CRS is interpreted incorrectly; zero or
  multiple arcs are rejected only after meshing. These implementation limits
  predate parts of this range; the new general-purpose claim exposes them.
  **Fix:** document and preflight both restrictions, or add explicit CRS and
  multiple-arc support before making the general claim.

- **F14 — Commas in recipe/output paths are not safe in `qsub -v`.**
  `jobs/octopus/refine_workflow.sh:44–49`.
  For a valid filename `a,b.yaml`, expansion yields
  `FMESH_RECIPE=recipes/refine/a,b.yaml,LR_SEEDS=0:1:2:3:4`. NQSV parses the
  comma as another variable, despite shell quoting; derived output paths
  have the same problem. Seeds already avoid this through colon separation.
  **Fix:** reject comma-containing transported values before any submission,
  or transport them through a job configuration file. This is a precise
  argument-construction finding; qsub was not invoked.

**Verification and bounds of the review**

Ran in `oceanmesh-bench`, with bytecode/cache writes disabled where applicable:

```bash
. /octfs/work/G16445/v61021/miniforge3/etc/profile.d/conda.sh
conda activate oceanmesh-bench
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider
ruff check --no-cache notebooks/420_local_refine.py \
  src/fvcom_mesh_tools/{patch,walls,bathy_patch,refine,qa,plotting}.py \
  src/fvcom_mesh_tools/dem/tokyo_bay.py \
  src/fvcom_mesh_tools/cli/{finish_depths,refine_depths,refine_run,plot_views}.py
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider \
  docs/hires_walls_tools_review_probes.py --tb=short
ruff check --no-cache docs/hires_walls_tools_review_probes.py
for script in jobs/octopus/{417_hires_refine,421_finish_and_run,423_m2_smoke,refine_workflow}.sh; do
  bash -n "$script" || exit
done
```

Existing suite: **786 passed, 7 rasterio deprecation warnings, 103.21 s**.
Scoped lint and shell syntax: **passed**. Review probes: **14 failed in
1.89 s**, all at the intended behavioral assertions (10 findings, including
parameterized cases). Probe lint initially found import ordering; corrected
in the new probe file and rerun successfully. These tests intentionally fail;
they are regression expectations, not xfails or tests of the buggy result.

Read each saved report, QA JSON, `.14`, constraints archive, and available
wall-stage archive. Independently compared mapped retained coordinates,
depths and the ordered OBC against the production base files: exact in all
four. Direct incidence counts found zero lone nodes in all four. Counts and
raw-depth implied dt below are from the saved artefacts, **not new runs**:

| Output | Nodes / elements | Saved QA | Introduced | Raw implied dt (s) |
|---|---:|---:|---:|---:|
| kimitsu_port_hires | 6,018 / 10,862 | 20/22 | 0 | 2.0563 |
| futtsu_coast_hires | 6,156 / 11,410 | 20/22 | 0 | 2.5909 |
| futtsu_nori_hires | 4,458 / 8,120 | 20/22 | 0 | 2.9668 |
| futtsu_nori | 4,413 / 8,039 | 21/22 | 0 | 2.4883 |

The saved verification reports also record zero missing retained faces,
split interface edges and nonmanifold edges. Kimitsu records 17 wall edges,
17 copies, seven free tips and one component. These raw dt values must not
be substituted for a finished-depth FVCOM step.

No additional frozen-node movement or copy-map corruption was established
in the split/open-lone/rejoin path. Retained faces precede patch faces;
opening requires both faces mutable; added nodes receive new identity
entries in `copy_of`; copied wall nodes are pinned during repair. Folded
rim mappings are propagated in reverse fold order. These protections and
the saved results are positive evidence, not proof for all geometries.
The lone-node gate uses actual incidence and its complete offender identities
survive display truncation. No demonstrated new misattribution defect was
found there. No direct GPL import was found in the scoped package modules;
the new wrapper uses a subprocess as required.

Compared `smooth_rfactor_equal` with the actual sibling TB-FVCOM algorithm.
On a five-node, three-triangle synthetic field `[3,12,120,400,650]`, the
all-movable outputs were bit-identical after **32 sweeps**. The equal method
floors, smooths uncapped depths and then caps, matching that algorithm's
order; no volume-conservation claim remains. **The full production-depth
reproduction claim was not rerun**: whole-mesh finishing belongs in a batch
job. Patch-only post-cap convergence and frozen-pair attribution deserve
separate acceptance tests; fixed-pair violations are currently reported but
excluded from the convergence Boolean without checking base connectivity.
No new mesh exhibiting that latter case was generated, so it remains a risk,
not an additional confirmed mesh finding.

The tests that passed do not cover adjacent simultaneous blunts, translated
closed rings, finite NetCDF fill values, offshore land discovery, or failed
workflow handoffs. Existing island/wall helper tests cannot establish driver
integration. No refinement, figure batch, full depth finish, FVCOM run,
qsub or qstat was executed. No dependencies were added.

After fixes, regeneration still needs the delegating session. An isolated
Kimitsu run, preserving existing outputs, can be submitted with:

```bash
qsub -v "FMESH_RECIPE=recipes/refine/kimitsu_port_hires.yaml,LR_OUT=outputs/review_hires_kimitsu,LR_SEEDS=0:1:2:3:4" jobs/octopus/417_hires_refine.sh
```

Use a fresh output name if that directory already exists; repeat for the
other three recipes before acceptance. The workflow's current 80-character
guard is conservative for its hardcoded root and input/output suffixes; no
off-by-one failure was established. A changed account root needs validation
of the actual generated namelist paths, not reliance on this constant.
