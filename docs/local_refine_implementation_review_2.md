# Second adversarial review of local refinement

Reviewed 2026-09-22 at **8caca27**, including bb3dd35. Read the first
`docs/local_refine_implementation_review.md` before reviewing the changes.
This is a review only; no implementation changes.

**The general-input contract is still not enforced.** The most consequential
remaining gaps are coverage of the intended water domain and complete coastline
provenance. The first review's exact-depth, face-multiplicity, adjacency, corner
pinning and explicit-geometry fixes are substantial; their original reproducers
are not presented here as new findings.

The new executable appendix produces **9 failures and 8 passing controls**.
Two failures concern the same artifact-publication finding. Existing repository
tests independently pass: **593 passed, 5 warnings in 43.50 s**. No DistMesh
production run or production QA rerun was performed. Severity below is about
impact on the affected input, not a claim that the measured Futtsu result is bad.
P1 means a contract blocker; P2 means a restricted-input failure or diagnostic /
operational defect. File references are relative to this repository.

## 1. P1 — An unchanged whole island is compared with another coastline

**CLAIM:** Giving all-free island points ID -1 fixes stitching, but the driver
cannot reliably process a cut containing both a whole island and an anchored
mainland stretch. The island is missing from the curves used for both sliding
and departure verification.

**EVIDENCE:** `patch.py:626-637` emits the whole ring without adding a curve or
`curve_of_pfix` entries; only the anchored branch adds a curve at line 653.
The driver makes every new boundary node slidable and passes the collected
curves (`420_local_refine.py:318-339`). It then measures **all** new boundary
edges against their union (`371-395`). `test_island_coastline_reference_covers_unchanged_island`
uses a 15x15 regular grid, removes the cell [600,700] x [600,700] as an island,
and cuts the box [-1,850] x [400,900]. Selection removes 83 faces, retains 307,
and adds no repair faces. The only curve is mainland x=0, y=400..900.
The unchanged island vertices measure **600, 700, 700, 600 m** from it. Thus the
driver's departure metric is at least **700 m**, despite zero island movement;
a 200 m tolerance rejects this identity geometry. Separate identity stitching
and frozen-map verification on this same input pass.

**WHY IT MATTERS:** The island fix is incomplete precisely on a common compound
input. Every seed can be rejected for a false geometric departure. With no
anchored stretch, `slide_on` is empty and the metric instead returns zero; adding
an unrelated mainland stretch changes the result. Nearest-curve binding can
also assign island nodes to the mainland; I did not demonstrate an accepted
island move onto that curve and do not claim one here.

**SUGGESTED FIX:** Register closed island curves, including their closing segment,
and carry explicit node/edge-to-curve provenance through stitching. In preserve
mode pin original island vertices. Measure each boundary component against its
own reference, including islands even when there are no anchored stretches.
Do not infer ownership from a nearest curve.

## 2. P1 — The verifier still accepts a missing piece of the patch

**CLAIM:** The new multiset and edge-count gates certify retained-face integrity,
not coverage of the intended final hole. Removing an interior patch face passes.
This is the unaddressed coverage part of first-review finding 2, now reproduced
against the new gates with a different mutation.

**EVIDENCE:** `verify_patch`, `patch.py:1031-1105`, compares retained multiplicities
and requires interface edges to have two faces, but does not constrain newly
created boundary edges inside the hole. `test_missing_patch_face_rejected`
removes face 101 from an identity stitch with a 250 m cut. All its vertices are
non-frozen and still used by other faces. Result: **ok=True**, 72 frozen nodes
exact, zero missing retained faces, zero interface splits, zero extra faces,
zero non-manifold edges, zero inverted elements, zero duplicates, zero orphans;
**area_change_fraction=-0.0078125**. One 5,000 m2 water triangle disappeared.

A separate small QA probe rebuilt boundary rings before and after that deletion:
both meshes passed **15/15 applicable gates**, with no OBC (the OBC checks are
skipped). This is not a claim of a 21/21 driver run: that driver requires an OBC.
It shows why ordinary mesh QA cannot distinguish an accidental new land hole
from a legitimate one without the intended domain.

**WHY IT MATTERS:** Centroid clipping and cleanup can leave gaps, and the final
boundary reconstruction can label the gap as another island. Neither retained
identity nor quality gates prove that the requested water remains water.
Reporting ungated area change is not sufficient.

**SUGGESTED FIX:** Compare the union of final patch faces against the **final rim
domain**, allowing only explicitly declared coastline changes. Check missing
and excess area separately, unexpected boundary components, crossings and vertex
fans. Preserve the intended coastline-component identities during serialization.

## 3. P1 — Per-stretch source selection still chooses a wrong bank on closest contact

**CLAIM:** Moving source selection inside `_source_substring` fixes the global
one-ring mistake, but `min(lines, key=here.distance)` still chooses by a single
closest pair, not correspondence of the stretch. This is a remaining portion
of first-review finding 11, with a new per-stretch counterexample.

**EVIDENCE:** `patch.py:539-553`. In
`test_source_assignment_considers_whole_stretch`, the base is [(0,0),(5,0),(10,0)].
One candidate follows it at y=0.1; the other departs diagonally from (0,0) to
(10,-8). The diagonal wins because it touches the first endpoint. The chosen
substring is **[(0,0),(6.09756098,-4.87804878)]**, although the other candidate
stays 0.1 from the entire stretch. These two candidate lines do not intersect
over their stated extents. The length fallback does not reject the chosen piece.

**WHY IT MATTERS:** Offering more rings within 3 km can introduce a wrong winner.
A sufficiently permissive coastline tolerance accepts a different shore instead
of identifying the mismatch. The per-stretch change alone is not a bank-selection
algorithm. The closed-ring seam fallback also remains at lines 550-553.

**SUGGESTED FIX:** Require compatible projections at both anchors and a bounded
fit along the whole base stretch; preserve source component identity when known.
Compare both cyclic routes for closed rings. Record fallback and ambiguity
explicitly rather than silently substituting preservation.

## 4. P2 — The best failed candidate can disagree with the mesh on disk

**CLAIM:** When no candidate passes QA, the retained best candidate, node map,
QA and achieved statistics can describe a different mesh from the `.14` file.
There is also a latent path-argument error in the new serialization helper.

**EVIDENCE:** `420_local_refine.py:510-527`: publication tests
`seed != seeds[0]`, rather than whether the last file written belongs to `best`.
`test_failed_best_seed_matches_file` executes the actual seed-search source with
lightweight deterministic mesher/QA substitutes: seed 0 has one failed gate,
seed 1 has two. `reports['seed']` and `node_map.npy` are **0**, but the file
contains **1**. The cached `written` and `qa` remain from seed 0. The later report
and dt calculations use those cached objects. The run eventually exits nonzero;
**this is not a demonstrated false successful exit**. If a candidate passes QA,
the loop breaks immediately, and I found no such mismatch in that path.

Separately, `serialise(candidate, out, path)` writes `path` at line 461 but reads
**global `out14`** at line 470. `test_serialise_reads_its_path_argument` executes
that function with candidate.14 and records a read of stale.14. Current callers
pass `out14`, so this second bug is latent until using per-attempt paths.

**WHY IT MATTERS:** A failed-run artifact can be investigated or resubmitted with
the wrong node map and QA evidence. An existing older `report.json` / node map
can also outlive an early failure, because reports are only saved at the end.

**SUGGESTED FIX:** Write each attempt to a distinct temporary artifact, verify
that path, then publish the chosen mesh, its map and its report together.
Always serialize the chosen candidate, or compare against the last written
candidate identity. Use `read_fort14(path)`. Clearly identify failed artifacts.

## 5. P2 — Candidate-specific exceptions bypass the seed search and its report

**CLAIM:** Several ordinary fill failures abort the whole search instead of
recording a failed attempt and trying the next seed.

**EVIDENCE:** `420_local_refine.py:271-307,488-497` has no exception boundary around
`attempt(seed)`. `stitch_patch` raises for a missing or merged constrained point
(`patch.py:918-925`), and the backend can raise its own fixed-point failure.
`test_seed_local_failure_does_not_abort_search` runs the actual search with seed
0 raising a lost-constraint `ValueError`; seed 1 is never called. With every
attempt returning None, lines 515-517 also exit saying “see attempts in the
report”, although `report.json` is not written until line 558.

**WHY IT MATTERS:** A bad first seed defeats the purpose of searching. Batch time
is lost, and failure diagnostics may refer to a nonexistent or previous report.

**SUGGESTED FIX:** Catch narrowly defined candidate-generation / candidate-contract
failures, append their diagnostics, and continue. Keep configuration and coding
errors fatal. Persist the attempts report on all exits, including total failure;
validate unsupported OBC counts before expensive meshing.

## 6. P2 — Parity is correct for strict nesting, but touching rings lose constraints

**CLAIM:** Correct nesting parity still does not validate the geometric rim.
Two rings touching along an edge with distinct IDs are silently unioned into a
different boundary.

**EVIDENCE:** `hole_polygon`, `patch.py:702-725`, calls `union_all` without checking
ring intersections or equality of constraints and result boundary.
`test_geometrically_touching_rings_rejected` supplies squares [0,2] x [0,2] and
[2,4] x [0,2] with distinct node IDs. It returns a valid rectangle; **2 units of
constraint geometry are absent from its boundary**. The test expecting rejection
fails. A passing control checks **five nesting levels**, shuffled ring order,
a disjoint shell, and a hole belonging to that smaller shell: correct area **72**
and no lost boundary constraints.

**WHY IT MATTERS:** A valid returned Polygon does not certify the original rim.
Such a touching-ID input is not a conforming base triangulation; this is a
validation gap, not a claim that `select_patch` creates it from a valid base.
Re-cut curves are another possible source of geometric contact. A point contact
also needs explicit treatment because degree checks operate on IDs, not geometry.

**SUGGESTED FIX:** Reject geometrically overlapping/touching rings unless a
supported contact topology is explicitly constructed. Check each ring's validity
and prove the returned domain boundary represents the constraints.

## 7. P2 — A valid sole retained face is always absorbed as a “spike”

**CLAIM:** Zero retained edge-neighbours does not imply a defective remnant when
the retained domain consists of exactly one triangle.

**EVIDENCE:** `_isolated_faces`, `patch.py:95-113`, marks that triangle isolated;
`select_patch:176-198` absorbs it and then refuses the empty mesh.
`test_valid_one_face_remainder_is_not_a_spike` takes one of the two triangles
of a 100 m square. The other is a valid edge-connected retained mesh, the rim
already has degree two, and the interface is one edge. Selection nevertheless
raises **“the footprint removes the entire mesh”**.

**WHY IT MATTERS:** This is a false rejection on a small but valid input, and a
semantic change from merely repairing vertex-attached spikes. It is not evidence
of production-scale runaway growth.

**SUGGESTED FIX:** Distinguish a sole connected retained face from isolated faces
coexisting with another retained component. If a minimum retained-domain size is
required, state and validate it separately instead of calling it non-manifold.

## 8. P2 — Clipped interpolation still does not bound the new-edge r-factor

**CLAIM:** The depth fix refreshes positions correctly, but the categorical
“refinement ... never worsen[s]” r-factor statement remains false for changed
connectivity. This is the remaining mathematical part of first-review finding 9.

**EVIDENCE:** `refine.py:342-363` proves a same-base-triangle bound and then applies
it to refinement generally. `test_depth_interpolation_does_not_guarantee_new_edge_rfactor`
uses a 3x2 grid with x-column depths **2, 3, 4.5**. Every old edge has r <= **0.2**.
New points (0.1,0.5) and (1.9,0.5) interpolate to **2.1, 4.35**, both inside the
base range, both inside the base domain, with outside count **0**. An edge between
them has r = **0.3488372093**. The segment is entirely in the water rectangle and
can be an edge of a new triangulation; no extrapolation or clipping is involved.
This is an interpolation/edge counterexample, not a measured DistMesh output.

**WHY IT MATTERS:** The global clip does not establish seabed-slope safety. The
current `run_qa` does not contain an r-factor gate, so its pass count does not
repair this claim. A new edge can connect values from different old triangles.

**SUGGESTED FIX:** Keep the same-triangle qualification, remove the global
promise, and report/check actual final-edge r-factor if that is a required
simulation constraint. Retain the interpolation source and outside counts.

## Remaining first-review findings, ranked by practical cost

These are residual obligations, not newly discovered copies of the old report.
Each row states the still-open claim, current evidence, impact and remedy.

| Priority / first finding | CLAIM | EVIDENCE | WHY IT MATTERS | SUGGESTED FIX |
| --- | --- | --- | --- | --- |
| P1 / 8 | No caller-enforceable maximum cut or correct segment-distance OBC guard | `patch.py:124-128,178-186,217-226,234-252`: no maximum envelope; all incident faces added; KD-tree uses OBC nodes | Can redefine what is frozen or approach the middle of a long OBC segment; cost is remeshing outside authorization | Bound the union of selected faces against an explicit allowed envelope; measure cut-to-OBC-polyline distance |
| P1 / 10 (tolerance portion) | Preserving corners does not fix resample/spline tolerance to the base | `patch.py:405-407` still tests only interior output vertices; driver `371-395` uses nine samples/edge against source curves, not the base | Wrong coastline can pass even though corner-sliding regression is fixed | Store base and final reference polylines per stretch; check full geometry against the declared metric |
| P1 / 5 (metadata portion) | Boundary type/grouping outside the patch can still change | Driver `431-460` refuses multiple OBCs but rebuilds every land ring as 20/21; `477-479` compares only OBC output vs constructed list and number of land lists | Simulation boundary semantics can change outside the cut | Map untouched base metadata verbatim; splice affected stretches; compare complete output lists/types against expected metadata |
| P2 / 12 (priority portion) | Declared overlapping-region priority is ignored | `refine.py:135-139`; driver `213-215,227-229`; `patch.py:816-826` takes minima without priority | Over-refinement adds mesh and time-step cost that the recipe did not request | Implement precedence or refuse unsupported overlaps |
| P2 / 13 | Effective-gradation / hard-C4-impossibility claim remains unsound | `patch.py:830-860`: reports `(ambient-target)/width`, omits the ambient spatial derivative, labels it a hard limit | Misleads transition choices and feasibility assessment; final QA limits downstream harm | Measure the actual field gradient; label as diagnostic, gate actual mesh C4 |
| P2 / 16 | Achieved resolution still absent from the driver's report | Driver `550-555` reports global dt and counts but no core edge statistics; hmin remains target/1.2 at `247` | A passing coarse patch can miss the requested refinement; calibration is empirical | Measure core/transition edge distributions and report target error on each accepted result |
| Low / 3 (literal bits / coordinate generality) | `.17g` depth output does not make all values bit-identical | `patch.py:1026-1027` uses numerical equality (signed zeros compare equal); `io/fort14.py:174` still uses `.15f` for coordinates | Mostly a literal-contract / local-near-zero-coordinate concern, not demonstrated UTM corruption here | Define numerical versus bit identity; use round-trip coordinate formatting if supporting arbitrary coordinates |

Disposition of **all 16** first-review findings: **1** strict nesting fixed,
invalid contacts still open (6 above); **2** multiplicity fixed, coverage open
(2); **3** depth serialization and numerical exact gate fixed, narrow residual
above; **4** stale adjacency / subset valence fixed in examined code; **5** zero /
multiple OBC explicitly refused, metadata still open; **6** island stitch fixed,
provenance incomplete (1); **7** edge connectivity fixed, singleton absorption
regression (7); **8** maximum extent / guard still open, last-round recheck fixed;
**9** depth refresh fixed, r-factor statement open (8); **10** correct-curve
preserve pinning fixed, tolerance/provenance still open; **11** global source
selection fixed, correspondence/cyclic seam still open (3); **12** zero width
fixed, priorities still open; **13** diagnostic claim still open; **14** explicit
shape projection fixed; **15** final QA is now called and failed QA causes a
nonzero exit, but failure-path artifact handling is broken (4–5); **16** global
achieved dt added and altitude comparison corrected (`refine.py:293`), achieved
resolution remains unmeasured by this driver.

## Claims I could not break, and limits of those checks

| Requested claim | What was actually checked / conclusion |
| --- | --- |
| 1 — face connectivity, spike loop | Face graph really joins shared edges. Growth is monotone (`removed |= add`), so no oscillation. For a finite mesh each changing iteration consumes another face, additionally bounded by `max_repair_rounds`. Isolated-face removal alone cannot recursively isolate edge-neighbours: an isolated face has none. Pinch repair can still grow the cut without a spatial bound. A one-round-limit control matched the default result and had no residual spike; post-loop recheck tests both defects correctly. Singleton exception: finding 7. |
| 2 — nesting parity | Could not break strict valid nesting: five levels, shuffled rings, disjoint components, nonlargest immediate parent all passed. Contact validation fails separately (6). |
| 3 — island IDs / stitching | Compound island/mainland identity stitch passed, correct frozen set, free base map entries -1. Duplicate constraint rows are **rejected**, not silently duplicated: `stitch_patch` requires unique nearest patch indices. Two rings sharing a node ID also violate `boundary_rings` degree two. No valid conforming shared-node-island counterexample was found. Downstream reference handling fails (1). |
| 4 — exact frozen / multiset / edge gates | One-ULP frozen-depth mutation rejected in new control; existing regression suite covers coordinate mutation and duplicate retained face. Multiset and >2 edge gates themselves were not defeated. Their conjunction still admits missing patch coverage (2). OBC coordinate equality does not validate metadata. |
| 5 — `.17g`, depth clip | Values immediately below/above 2, 1e-100, 1e100 and a long mantissa round-tripped exactly. Affine interpolation and two outside fallback points gave the expected values and outside count 2. No real extrapolation bug hidden by clipping was demonstrated. The clip is global and silent, so it should report any substantial overshoot instead of assuming every overshoot is roundoff. |
| 5b — curve-vertex pin and span clamp | Could not defeat preserve mode when supplied the correct full curve: original corner and rectangle-slide tests pass. `_vertex_span` pins curve vertices; straight-segment candidates remain within their original vertex bracket. Positive incident-face area prevents passing an adjacent collinear boundary node on a valid fan. This does not fix absent/wrong provenance (1) or resampled chords not originally on the unsimplified source. |
| 6 — repaired adjacency, global valence, monotonicity | A new instrumented 5x5 perturbed-grid test checked global worst margin and global valence for every candidate accepted by `_better`, including soft comparisons, during one round. No decrease beyond 1e-12 and no valence >8. By inspection, moves score all incident triangles and every edge-adjacent area ratio; flips score every face incident to their four vertices and rebuild adjacency. Thus all changing angle/area terms are represented: local hard monotonicity implies global hard monotonicity for these terms, not just the neighbourhood. No counterexample found on valid topology. Soft and cleanup comparisons permit a 1e-12 numerical decrease; literal exact monotonicity is not claimed. Valence is a separate veto; an already excessive immutable valence may remain. Cleanup cannot cycle indefinitely: it is round-bounded and each accepted cleanup flip strictly reduces total positive excess valence, without creating new excess at receiving vertices. |
| 7 — width zero | No NaN at interior, boundary, exterior or out-of-base samples. Interior target is capped by the base size (no coarsening), consistent with current sizing policy. |
| 8 — declared geometry | Existing bbox test and new projected polygon-with-hole test pass. Parser accepts Polygon, not MultiPolygon, and explicitly rejects the latter; that is not a projection regression. |
| 9 — driver integration / state | Refresh is after repair and before verification. Unsupported OBC count really is refused, but late. QA uses reread coordinates/depths. Each successful passing seed ends the loop. New per-stretch source selection is insufficient (3); retry/publication has failures (4–5). No mutation leak of `sel`, `hole`, prepared geometry, `rc`, or old depths was established. `improve_patch` copies arrays; `refresh_depths` returns a copy. Read-only inspection of sibling `oceanmesh/mesh_generator.py:1864-1872` shows pfix copied by `_unpack_pfix`; `mesh_improve.py:57-64` copies collapse inputs and uses pfix for lookup. No GPL library was imported for review probes and no real multiseed meshing run was done. |

The score proof assumes a valid planar triangulation and truthful ownership
masks. It is not a proof of geometric coverage, vertex manifoldness, or a global
soft-score ascent: saturated neighbourhood sums can count shared area terms
differently. Neither of those is the claimed hard-margin guarantee.

Fourth-column compatibility audit: the central reader uses whitespace numeric
`np.loadtxt` (`io/fort14.py:103-110`), not fixed-width or mandatory-exponent parsing.
Searches of package I/O and notebook depth consumers found no demonstrated
`.17g` incompatibility. FVCOM native export still formats depth with `.6f`
(`io/fvcom_native.py:46,103`), so native files do not promise exact fort.14 depth
identity; this is existing downstream quantization, not a regression introduced
by `.17g`. Full unit tests cover the installed package's readers/writers; no
external solver execution was performed.

## Missing tests that would have caught these defects

Add the appendix failures at the contract boundary they exercise. In particular,
use an island **and** an anchored coast, not only an interior island cut; test
removed patch coverage, not only retained faces; and execute the driver's seed
state machine with deterministic substitutes instead of depending on lucky
production seeds. For serialization, vary the destination path. Test no passing
seed, best-first failed seed, exception on first seed, no contract-valid seed,
and pre-existing output artifacts. Require non-vacuous accepted-operation
instrumentation for repair. Compare every constrained rim segment with the
assembled domain boundary under strict nesting and geometric contacts.

Further valuable cases not fully executed here: actual two-curve/island repair
with provenance carried through the backend; closed-source seam rotation;
full-polyline tolerance under simplification; outside-depth counts around a
concave coast; complete OBC/land metadata round trips; maximum allowed cut and
long-segment OBC guards. These should precede additional successful-site seeds.

## Verification and reproducibility

Only the new review document is a repository change. No output mesh was modified,
no dependencies installed, no sibling repository edited, and no commit/branch
created. Small synthetic/unit probes ran on the login node. No batch job was
submitted or bypass attempted. No production node/element counts, achieved dt,
or QA metrics were recomputed; the numbers in this report are the synthetic
measurements identified above.

Commands, from repository root:

```bash
MPLCONFIGDIR=/tmp/lr2-mpl \
/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python \
-m pytest -q /tmp/test_local_refine_review2.py --tb=short
# 9 failed, 8 passed in 3.11 s (the intentional regression assertions below)

MPLCONFIGDIR=/tmp/lr2-mpl \
/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python -m pytest -q
# 593 passed, 5 rasterio PendingDeprecationWarnings in 43.50 s
```

Lint verification:

```bash
/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/ruff check .
# Failed: 1164 existing errors across the repository, principally unrelated notebooks.

/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/ruff check \
  src/fvcom_mesh_tools/patch.py src/fvcom_mesh_tools/refine.py \
  src/fvcom_mesh_tools/io/fort14.py tests/test_patch.py \
  notebooks/420_local_refine.py --output-format concise
# All checks passed. No lint fixes applied.
```

An earlier 14-test version produced 8 failed / 6 passed; the final version adds
two passing controls and the serialization-path failure. The final rerun also moved the wrong-bank
fixture below the base line so the two offered source lines are disjoint. Initial exploratory
circular island/mainland cuts were refused for splitting the retained domain;
the reported box fixture was then chosen and passes selection. One exploratory
Shapely MultiLineString construction needed `list(array)` rather than an ndarray;
that harness error is not counted as a product finding.

## Executable appendix

Save this block as `/tmp/test_local_refine_review2.py` and run the first command
above. Tests are ordinary assertions, not xfails. Driver tests extract actual
source and replace expensive meshing / I/O observation where stated; they do
not import the batch driver or GPL backend. Nine failures are expected at 8caca27.

```python
from pathlib import Path
from types import SimpleNamespace
import ast
import runpy
import numpy as np
import pytest
import shapely
from pyproj import Transformer
import fvcom_mesh_tools.patch as patch
from fvcom_mesh_tools.refine import depths_from_base, RefineRegion
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14, write_fort14

ROOT = Path('/octfs/work/G16445/v61021/Github/fvcom-mesh-tools')
DRIVER = (ROOT / 'notebooks/420_local_refine.py').read_text()
helpers = runpy.run_path(str(ROOT / 'tests/test_patch.py'))
grid = helpers['grid_mesh']


def extract_function(name, env):
    fn = next(n for n in ast.parse(DRIVER).body
              if isinstance(n, ast.FunctionDef) and n.name == name)
    exec(compile(ast.Module(body=[fn], type_ignores=[]), '<driver>', 'exec'), env)
    return env[name]


def identity(radius=250):
    xy, tri = grid()
    dep = np.ones(len(xy)) * 8
    sel, rc, pn, pt = helpers['replay_patch'](
        xy, tri, shapely.Point(400, 400).buffer(radius))
    out, faces, depth, nm, _ = patch.stitch_patch(
        xy, tri, dep, sel, pn, pt, rc['pfix'], rc['pfix_base'])
    return xy, tri, dep, sel, out, faces, depth, nm


def test_missing_patch_face_rejected():
    xy, tri, dep, sel, out, faces, depth, nm = identity()
    k = next(i for i in range(len(sel.retained), len(faces))
             if not np.isin(faces[i], nm[sel.frozen_nodes]).any())
    ver = patch.verify_patch(xy, dep, tri, sel, out,
                            np.delete(faces, k, axis=0), depth, nm)
    assert not ver['ok'], ver


def island_selection():
    xy, tri = grid(15, 15)
    c = xy[tri].mean(1)
    tri = tri[~((c[:, 0] > 600) & (c[:, 0] < 700)
                & (c[:, 1] > 600) & (c[:, 1] < 700))]
    sel = patch.select_patch(xy, tri, shapely.box(-1, 400, 850, 900))
    rc = patch.rim_constraints(xy, sel, size=1e9, coastline='preserve')
    return xy, tri, sel, rc


def test_island_coastline_reference_covers_unchanged_island():
    xy, _, sel, rc = island_selection()
    island = next(r for r, hole in zip(sel.rings, sel.ring_is_hole) if hole)
    # The driver's curve union and metric, evaluated even at unchanged nodes.
    curves = shapely.MultiLineString(rc['curves'])
    departure = shapely.distance(shapely.points(xy[island]), curves).max()
    assert departure == 0, departure


def ring_input(boxes):
    xy = np.vstack([np.array([[a,b],[c,b],[c,d],[a,d]], float)
                    for a,b,c,d in boxes])
    edges = np.array([[i, i//4*4+(i+1)%4] for i in range(len(xy))])
    return xy, edges


def test_geometrically_touching_rings_rejected():
    xy, edges = ring_input([(0,0,2,2), (2,0,4,2)])
    with pytest.raises(ValueError):
        patch.hole_polygon(xy, edges)


def seed_env(tmp_path, *, exceptional=False):
    def attempt(seed):
        if exceptional and seed == 0:
            raise ValueError('one constrained point absent from filled patch')
        return (np.array([seed]), None, None, np.array([seed])), {'seed': seed}
    def serialise(candidate, out, path):
        seed = int(candidate[0][0])
        path.write_text(str(seed))
        return SimpleNamespace(seed=seed), None
    def qa(mesh, **kwargs):
        return SimpleNamespace(n_gate_total=21, n_gate_failed=mesh.seed + 1,
                               checks=[])
    return dict(Path=Path, np=np, OUT=tmp_path,
                cfg={'base_mesh': 'base.14'}, recipe=Path('recipe.yaml'),
                os=SimpleNamespace(environ={'LR_SEEDS': '0,1'}), reports={},
                say=lambda *a: None, attempt=attempt, serialise=serialise,
                run_qa=qa,
                read_fort14=lambda p: SimpleNamespace(seed=int(p.read_text())))


def run_seed_selection(env):
    start = DRIVER.index('out14 = OUT /')
    end = DRIVER.index('say(f"accepted seed {seed}")', start)
    exec(compile(DRIVER[start:end], '<actual seed search>', 'exec'), env)


def test_failed_best_seed_matches_file(tmp_path):
    env = seed_env(tmp_path)
    run_seed_selection(env)
    assert int(env['out14'].read_text()) == env['reports']['seed']


def test_seed_local_failure_does_not_abort_search(tmp_path):
    env = seed_env(tmp_path, exceptional=True)
    run_seed_selection(env)
    assert env['reports']['seed'] == 1


def test_valid_one_face_remainder_is_not_a_spike():
    xy, tri = grid(2, 2)
    sel = patch.select_patch(xy, tri,
        shapely.Point(xy[tri[0]].mean(0)).buffer(1))
    assert len(sel.retained) == 1


def test_source_assignment_considers_whole_stretch():
    pts = np.array([[0,0],[5,0],[10,0]], float)
    correct = shapely.LineString([(0,.1),(10,.1)])
    wrong = shapely.LineString([(0,0),(10,-8)])
    chosen = patch._source_substring(pts, [correct, wrong])
    assert np.allclose(chosen[:, 1], .1), chosen


def test_depth_interpolation_does_not_guarantee_new_edge_rfactor():
    xy, tri = grid(3, 2, 1)
    dep = 2 * 1.5 ** xy[:, 0]
    edges, _ = patch._edge_table(tri)
    before = (np.abs(dep[edges[:,0]] - dep[edges[:,1]])
              / dep[edges].sum(1)).max()
    q = np.array([[.1,.5],[1.9,.5]])
    new, outside = depths_from_base(xy, tri, dep, q)
    assert outside == 0
    after = abs(new[0] - new[1]) / new.sum()
    assert after <= before, (before, after, new)


def test_five_levels_disjoint_and_nonlargest_parent():
    xy, edges = ring_input([(2,2,8,8), (20,0,24,4), (0,0,10,10),
                           (3,3,7,7), (21,1,23,3), (1,1,9,9), (4,4,6,6)])
    out = patch.hole_polygon(xy, edges)
    assert out.is_valid
    assert out.area == pytest.approx(100-64+36-16+4+16-4)
    assert shapely.MultiLineString(list(xy[edges])).difference(out.boundary).is_empty


def test_last_repair_round_rechecked():
    xy, tri = grid()
    foot = shapely.Point(400,200).buffer(150)
    a = patch.select_patch(xy, tri, foot)
    b = patch.select_patch(xy, tri, foot, max_repair_rounds=1)
    assert np.array_equal(a.removed, b.removed)
    assert not patch._isolated_faces(tri, b.removed).any()


def test_island_identity_and_free_map():
    xy, tri, sel, rc = island_selection()
    taken = np.unique(tri[sel.removed]); local = np.full(len(xy), -1)
    local[taken] = np.arange(len(taken))
    out, faces, dep, nm, _ = patch.stitch_patch(
        xy, tri, np.ones(len(xy)), sel, xy[taken], local[tri[sel.removed]],
        rc['pfix'], rc['pfix_base'])
    assert patch.verify_patch(xy, np.ones(len(xy)), tri, sel,
                             out, faces, dep, nm)['ok']
    assert (nm[sel.free_nodes] == -1).all()


def test_global_margin_and_valence_on_accepted_operations(monkeypatch):
    xy, tri = grid(5, 5)
    xy[12] += [30, -20]
    original_score, original_better = patch._scores, patch._better
    records = {}; accepted = []
    def score(x, t, faces, adj, *gates):
        result = original_score(x, t, faces, adj, *gates)
        global_score = original_score(x, t, np.arange(len(t)), adj, *gates)[0]
        records[id(result)] = (result, global_score, np.bincount(t.ravel()).max())
        return result
    def better(after, before, soft):
        ok = original_better(after, before, soft)
        if ok:
            accepted.append(1)
            assert records[id(after)][1] >= records[id(before)][1] - 1e-12
            assert records[id(after)][2] <= 8
        return ok
    monkeypatch.setattr(patch, '_scores', score)
    monkeypatch.setattr(patch, '_better', better)
    patch.improve_patch(xy, tri, np.arange(len(xy)) == 12,
                        np.ones(len(tri), bool), rounds=1, soft=True)
    assert accepted


def test_width_zero_depth_roundtrip_and_exact_gate(tmp_path):
    xy, tri, dep, sel, out, faces, depth, nm = identity()
    field = patch.patch_sizing(xy, tri,
        [(shapely.box(300,300,500,500), 30, 0)], distmesh_scale=1)
    assert field([[400,400]])[0] == 30
    assert np.isfinite(field([[200,200], [300,400], [900,900]])).all()
    values = np.resize([np.nextafter(2.,0.), np.nextafter(2.,3.),
                        1e-100, 1e100, 5.12345678912345], len(xy))
    mesh = Fort14Mesh(title='precision', nodes=xy, depths=values, elements=tri,
                      open_boundaries=[], land_boundaries=[])
    path = tmp_path / 'depth.14'; write_fort14(mesh, path)
    assert np.array_equal(read_fort14(path).depths, values)
    depth[nm[sel.frozen_nodes[0]]] = np.nextafter(depth[0], np.inf)
    assert not patch.verify_patch(xy, dep, tri, sel, out, faces, depth, nm)['ok']


def test_projection_keeps_polygon_holes():
    env = dict(np=np, shapely=shapely,
               to_m=Transformer.from_crs(4326,32654,always_xy=True))
    fn = extract_function('region_in_metres', env)
    region = SimpleNamespace(kind='polygon', geometry=shapely.Polygon(
        [(139.7,35.3),(139.8,35.3),(139.8,35.4),(139.7,35.4)],
        [[(139.72,35.32),(139.74,35.32),(139.74,35.34),(139.72,35.34)]]))
    projected = fn(region)
    assert projected.is_valid and len(projected.interiors) == 1
    assert len(projected.exterior.coords) == 5


def test_duplicate_constraint_rows_refused_not_duplicated():
    xy, tri, sel, rc = island_selection()
    taken = np.unique(tri[sel.removed]); local = np.full(len(xy), -1)
    local[taken] = np.arange(len(taken))
    row = np.flatnonzero(rc['pfix_base'] == -1)[-1]
    with pytest.raises(ValueError, match='two constrained points'):
        patch.stitch_patch(xy, tri, np.ones(len(xy)), sel,
            xy[taken], local[tri[sel.removed]],
            np.vstack([rc['pfix'], rc['pfix'][row]]),
            np.append(rc['pfix_base'], -1))


def test_depth_outside_fallback_is_counted():
    xy, tri = grid(3, 3, 1)
    dep = 5 + xy[:, 0] + 2 * xy[:, 1]
    q = np.array([[.2,.3], [1.7,1.8], [-.1,-.1], [3,3]])
    result, n = depths_from_base(xy, tri, dep, q)
    assert n == 2
    assert np.allclose(result, [5.8,10.3,5,11])


def test_serialise_reads_its_path_argument(tmp_path):
    xy = np.array([[0,0],[100,0],[100,100],[0,100]], float)
    tri = np.array([[0,1,2],[0,2,3]])
    base = Fort14Mesh(title='base', nodes=xy, depths=np.ones(4)*8,
                     elements=tri, open_boundaries=[np.array([0,1])],
                     land_boundaries=[(20,np.array([1,2,3,0]))])
    seen = []
    def reader(path):
        seen.append(path)
        return read_fort14(tmp_path / 'candidate.14')
    env = dict(np=np, shapely=shapely, base=base, sel=None,
               recipe=Path('recipe.yaml'), Fort14Mesh=Fort14Mesh,
               om=SimpleNamespace(boundary_loops=lambda _: [np.arange(4)]),
               say=lambda *a: None, write_fort14=write_fort14,
               read_fort14=reader, verify_patch=lambda *a, **kw: {'ok': True},
               out14=tmp_path / 'stale.14')
    serialise = extract_function('serialise', env)
    candidate = (xy, tri, base.depths, np.arange(4))
    serialise(candidate, {}, tmp_path / 'candidate.14')
    assert seen == [tmp_path / 'candidate.14'], seen
```
