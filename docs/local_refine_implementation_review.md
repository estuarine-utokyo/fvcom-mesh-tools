# Adversarial review of local refinement, revision 4

Reviewed 2026-09-22. Scope: all of `patch.py`, `refine.py`, the 420 driver,
`test_patch.py`, the design and Futtsu recipe. File references below are relative
to this repository. This is a review, not an implementation change.

**Recommendation: do not treat this generator as enforcing the stated contract
on general inputs.** The successful Futtsu build does not exercise the failures
below. Nine executable regression tests at the end all fail on this checkout.
P1 denotes a release blocker for the affected supported input or contract;
P2 denotes a misleading diagnostic or an additional robustness limitation.

## 1. P1 — Polygon construction fills every nested exclusion

**CLAIM:** Correct nesting flags are discarded before meshing; `egfix` is not
necessarily the boundary of the domain given to DistMesh.

**EVIDENCE:** `patch.py:628-632` polygonizes all edges and unions **all** resulting
faces. For outer square area 100 and inner square area 16, `boundary_rings`
returns `[False, True]`, but `hole_polygon` returns area 100, not 84. Four nested
squares plus a disjoint square give correct flags `[False, True, False, True,
False]`, but area **101 instead of 57**. See `test_nested_domain_preserves_parity`.
The inner constraints become interior lines, absent from `hole.boundary` used
by `420_local_refine.py:205-213` and its centroid clipping at lines 251-254.

**WHY IT MATTERS:** Land islands can be filled with water; a retained island,
if allowed through selection, can be overlaid with new faces. A valid Shapely
polygon is not evidence that it represents the intended water domain.

**SUGGESTED FIX:** Build shells with their immediate odd-depth holes, and union
only even-depth water components. Validate that the resulting boundary equals
the emitted constraints geometrically. Reject crossings/touches that cannot be
represented as the intended manifold domain; do not silently repair them by
noding and unioning every face.

## 2. P1 — Verification accepts extra frozen faces and non-manifold edges

**CLAIM:** Retained-face *presence* is checked, but retained-face multiplicity,
exterior occupancy and general mesh conformity are not.

**EVIDENCE:** `patch.py:901-903` reduces both face lists to sets of sorted node
IDs. Lines 929-935 check edge count only on designated interface segments.
Appending a duplicate of the first retained face to the identity-fill fixture
returns **`ok=True`**, no missing faces, no inverted elements, no duplicate nodes,
no split interfaces, and `area_change_fraction=0.0078125`. See the first appendix
test. Its positive-area duplicate creates edges with three incident faces away
from the interface. Area is not gated (`patch.py:960-962`).

**WHY IT MATTERS:** The frozen exterior can acquire overlapping elements while
all advertised checks pass. Interface edge multiplicity two also does not prove
opposite-side incidence or exclude a separate hanging vertex on that segment.
The verifier does not compare against the actual hole domain or check extra
boundary loops caused by missing patch faces.

**SUGGESTED FIX:** Compare retained oriented-face multisets exactly under a valid,
injective map; reject extra faces in the frozen exterior. Check all edge counts,
vertex fans, opposing incidence, geometric crossings/overlap, and coverage of the
intended final hole. Gate area/coverage against that hole, which already accounts
for legitimate coastline changes, instead of accepting arbitrary area change.

## 3. P1 — “Bit-for-bit” is neither checked nor preserved on disk

**CLAIM:** Frozen floating-point values can change without rejection, including
changes introduced after verification by the writer.

**EVIDENCE:** `patch.py:898-899,950,960` accepts coordinate and depth differences
up to `1e-6`. Editing each by `5e-7` passes. The driver verifies at
`420_local_refine.py:347`, then writes at line 384 without re-reading.
`io/fort14.py:168` formats depths as `.10e`: depth `5.12345678912345` round-trips
with a change of **-2.3450574815342407e-11 m**. Both cases have failing tests below.
This is a general precision failure; it does not assert that the specific Futtsu
input contains more precision than that formatter retains.

**WHY IT MATTERS:** This violates the explicit identity contract even where the
physical magnitude is negligible. A successful in-memory check cannot certify
the delivered fort.14.

**SUGGESTED FIX:** Separate the numerical tolerance for matching fill vertices
from exact comparison of copied frozen data. Use round-trip float formatting for
all values, and re-read and verify the written mesh, boundaries included. If
literal float-bit identity includes signed zero, compare bit representations.

## 4. P1 — Accepted flips can violate max_valence; scoring uses stale topology

**CLAIM:** The repair's valence guard is incorrect after an earlier flip in the
same sweep, and its monotonicity argument uses incorrect adjacency.

**EVIDENCE:** `patch.py:1136-1140` computes incidence and adjacency once per flip
sweep. Lines 1166-1167 mark only the two rewritten faces as touched. Later flips
with different face IDs can share their vertices or neighbours. Thus their
`faces` union at lines 1152-1153 omits newly incident faces, and `before` at
1154 uses stale adjacency. `_valence_ok` at lines 1292-1298 counts only that subset.

On the existing deterministic `perturbed_patch()` fixture, one flip-only sweep
accepts a flip whose subset counts are **[5,5,7,7]**, while true counts are
**[5,5,7,9]**, with `max_valence=8`. The appendix checks every accepted flip and
fails with `9 <= 8`. This sweep happens to end at maximum valence 8; the bug is
an accepted intermediate violation, not a claim that this fixture's final mesh
has valence 9.

Instrumenting `_health` on the same sweep found 79 calls with stale adjacency
on the scored faces, eight with a different score; one was **0.1986335967 instead
of 0.2806914797**. The post-flip score uses fresh adjacency, so it is compared
against a different, incorrect pre-flip objective. I did not establish a decrease
in the *global* minimum health on this fixture.

**WHY IT MATTERS:** The advertised acceptance predicates are false. Later sweeps
need not undo a forbidden valence increase. “Strictly improving” cannot be
inferred from these comparisons.

**SUGGESTED FIX:** Update/recompute incidence and adjacency after accepted flips,
or invalidate every affected vertex neighbourhood before the next candidate.
Count valence globally or maintain exact counts with flip deltas. Compare before
and after on the same complete affected neighbourhood. Test every accepted
operation, not just final minimum angle.

## 5. P1 — Open/land boundary identity is not verified, and the writer assumes one OBC

**CLAIM:** `open_boundary_unchanged` checks old node coordinates, not the output
boundary lists; output boundary metadata is rebuilt after verification.

**EVIDENCE:** `verify_patch` has no argument for new open/land boundary lists.
`patch.py:910-914` merely evaluates the old list through the map. In the driver,
`420_local_refine.py:368-381` indexes `obc[0]`, tests a **set**, and writes one
new segment `ring[:stop+1]`, not the list `obc` it just mapped. Zero OBCs cause
an IndexError at line 368. Multiple disjoint OBCs generally fail the set test;
adjacent segments can be merged into one. All mainland/island land-boundary
types and grouping are reconstructed as 20/21 at lines 374-375, regardless of
base metadata outside the patch.

**WHY IT MATTERS:** Boundary segmentation/order/type is part of the simulation
input, not merely geometry. The current verifier cannot notice its corruption.
Even if the one-OBC Futtsu case preserves its order, that does not implement the
list-of-segments contract exposed by fort.14.

**SUGGESTED FIX:** Carry every OBC segment through the map verbatim. Preserve
unmodified land-boundary metadata and splice only the affected runs. Validate
actual output lists, ordering, types and boundary-edge incidence after boundary
construction and after serialization. Explicitly support or preflight-refuse
zero/multiple OBCs before meshing.

## 6. P1 — A physical island inside the cut cannot be stitched

**CLAIM:** The all-free ring branch gives free nodes IDs which stitching assumes
must already have survived as frozen nodes.

**EVIDENCE:** `patch.py:564-567` keeps the base IDs of an entirely free physical
island ring. `stitch_patch` maps only frozen nodes (`799-801`), interprets every
nonnegative `pfix_base` as frozen (`818-820`), then raises at 821. The appendix
constructs a regular grid with one central cell removed as land and selects the
surrounding water. The retained exterior remains connected. Four island rim
nodes are free; even an exact identity fill fails with “selection and the rim
disagree”.

**WHY IT MATTERS:** A circle enclosing a small physical island is an ordinary
valid refinement input. This fails independently of the polygon-union bug.
The comment's “every maximal run ... bounded by two frozen nodes” is false for
an entire physical boundary component, even though the narrower free-rim-node
claim is true.

**SUGGESTED FIX:** Distinguish frozen references from preserved-but-free base
nodes. Either carry the latter into the node map explicitly or emit them as new
points with appropriate provenance. Implement refinement of closed unanchored
coastline rings instead of preserving arbitrarily coarse island edges.

## 7. P1 — Retained connectivity is tested on nodes, not face adjacency

**CLAIM:** The code accepts a cut that leaves retained pieces joined at only one
vertex, despite explicitly promising edge connectivity.

**EVIDENCE:** `_require_connected`, `patch.py:237-242`, creates a graph of **node
IDs** linked by triangle edges. Two triangles sharing one node are connected in
that graph. The appendix's six-node, four-triangle semicircular fan removes its
two middle triangles. The remaining two share only vertex 0. Selection succeeds,
but the test expects the documented refusal. The removed rim has degree two,
so rim repair does not catch this case.

**WHY IT MATTERS:** The retained domain has been severed in the sense the caller
explicitly prohibits. A manifold removed rim alone does not prove that the
retained domain is a manifold with boundary.

**SUGGESTED FIX:** Use a graph whose vertices are retained faces and whose links
are shared edges. Validate retained boundary vertex fans separately. Compare
component structure with the original if initially disconnected domains are to
be supported.

## 8. P1 — Cut growth has no caller-enforceable spatial maximum

**CLAIM:** Selection can redefine the frozen zone beyond the requested envelope;
its report is not a bound, and its OBC guard measures the wrong geometry.

**EVIDENCE:** `patch.py:136-159` selects by centroid, then takes all faces at bad
vertices. There is no allowed-envelope/max-overreach argument or recipe key.
The only reach reported (`202-219`) is distance from the footprint centroid,
which is not departure from an arbitrary/multi-component footprint. The driver
prints this and immediately proceeds (`420_local_refine.py:163-190`). Its guard
at `patch.py:189` is **node-to-node** distance, not distance from the cut to OBC
segments. A long OBC segment can pass close to a removed node while both endpoint
distances exceed the guard.

**WHY IT MATTERS:** A distant whole triangle or repair fan can be remeshed beyond
what the caller permits, yet verification passes because it trusts the enlarged
selection's definition of frozen. Reporting growth cannot enforce the design's
“any enlargement must stay inside a declared maximum”. The OBC guard can also
underestimate proximity along long segments.

**SUGGESTED FIX:** Require an explicit maximum allowed cut/enlargement; compare
the union of removed triangles with it after repair. Report actual excess area
and maximum departure from the allowed geometry. Use distance between cut
geometry and complete OBC polylines, including the planned re-cut shoreline.

Termination itself is bounded: each successful repair iteration adds faces and
there is a 50-round cap. I found no infinite loop for a valid finite mesh. A
repair resolving the last pinch on iteration 50 still raises because the loop
never performs the next success check; check the final rim before declaring
failure at the cap.

## 9. P1 — Repair moves nodes but leaves their depths at the old positions

**CLAIM:** The delivered new-node bathymetry is not necessarily the base field
evaluated at the delivered coordinates.

**EVIDENCE:** Depth interpolation occurs at `patch.py:838-843`. The driver then
replaces `nodes,elements` at `420_local_refine.py:301-303`, but never recomputes
`depths`. `_before_repair` at line 293 is unused. Verification only compares
frozen depths. A small fan with corners `(0,0),(2,0),(2,2),(0,2)` and interior
point `(0.2,0.3)` moves that point to `(1,1)` in one repair round. For affine base
depth `5+x`, the attached depth remains **5.2 instead of 6.0**.

**WHY IT MATTERS:** The bathymetry contract and any achieved dt estimate based on
those depths are wrong, even though frozen-zone verification passes.
Furthermore, `refine.py:331-339` proves an r-factor bound only for points in the
same old triangle; remeshing can create edges across several old triangles, so
its global “never worsen” claim does not follow.

**SUGGESTED FIX:** Re-evaluate base depths for every moved non-frozen node after
repair; report extrapolation counts again. Measure achieved r-factor and dt on
the finished output. Do not claim same-triangle interpolation proves an r-factor
bound for arbitrary new connectivity.

## 10. P1 — `preserve` can lose its corners during repair; tolerance changes reference

**CLAIM:** The driver does not preserve the strict polyline invariant after rim
construction, and the final tolerance is not the promised tolerance to the base.

**EVIDENCE:** `_subdivide` keeps original vertices, but `rim_constraints` marks
all intermediate results, including those original free corners, as new with
base ID -1 (`patch.py:586-588`). The driver marks all new boundary nodes slidable
(`289-291`). `_candidates` projects to the **whole** curve (`1241-1246`), with no
constraint to retain a corner or stay between neighbouring arc-length stations.
The code itself acknowledges corner cutting at driver lines 335-338 but only
refuses it above the ordinary tolerance, which is 200 m in this recipe.

Initially `coastline_points` checks only output **vertices** against the base
polyline (`patch.py:365-371`). A measured example uses the base V-shaped chain
`[(0,0),(1,1),(2,0)]`, straight source `[(0,0),(2,0)]`, `size=1e9`, and
`tolerance_m=0.1`. `coastline_points(..., mode='resample')` accepts the two
endpoints, whose connecting segment has **Hausdorff departure 1 m** from the
base, ten times the tolerance. There are no intermediate vertices to check.
After repair the driver checks
nine samples per edge against the **source curve**, not the base (`315-333`).

**WHY IT MATTERS:** In preserve mode a surviving corner may slide away and be
chorded off while every new node is on the original polyline. In other modes,
staying near the source does not establish the specified maximum departure from
the base. Coarse sampling can miss a peak between samples.

**SUGGESTED FIX:** Pin original corners in preserve mode and restrict inserted
points to their original segments with ordered stations. Retain curve provenance
through stitching instead of recovering it by nearest curve. Check final full
polylines against the declared reference with a specified one-/two-sided metric
and a conservative tolerance algorithm.

## 11. P1 — One nearest shoreline is used for every coastline stretch

**CLAIM:** The driver contradicts its own per-stretch source-ring requirement.

**EVIDENCE:** `420_local_refine.py:173-179` makes one MultiPoint from **all** free
rim nodes and chooses one exterior ring by minimum distance. That is the closest
pair distance, not a fit to each stretch. It passes that one line to all runs in
`rim_constraints`. Land-polygon interior rings are excluded. `_source_substring`
(`patch.py:495-506`) does not wrap a closed ring across its coordinate seam; it
falls back to the base when the non-wrapped route exceeds three times base length.

**WHY IT MATTERS:** Two patches on different islands/banks can project to the
same source ring. A tolerance check can accept a wrong bank at a narrow strait;
fallback can silently turn requested resampling into preservation. Translating
the starting vertex of a closed source ring can change the chosen route.

**SUGGESTED FIX:** Associate each individual free run with its actual source
component, including interior rings, before resampling. Choose the appropriate
cyclic substring using the existing stretch, and report/refuse ambiguity or
fallback rather than silently changing modes.

## 12. P1 — Accepted zero-width transitions produce NaN; priorities are ignored

**CLAIM:** The field is undefined for a permitted no-coarsening input, and recipe
priority has no effect on overlapping regions.

**EVIDENCE:** `refine.py:102-103` derives width 0 when ambient <= target.
`patch.py:723` computes `d/width`, including `0/0` throughout the core; the
appendix obtains `[nan]` on a uniform grid with target 200 m. `effective_gradation`
can nevertheless report within-limit because `d < 0` selects no samples.

`RefineRegion` accepts priority at `refine.py:122-126`, but the driver passes only
`(geometry,target,width)` to sizing (`187-203`), and `patch.py:720-725` takes the
minimum across regions. A high-priority 60 m core overlapped by a low-priority
30 m core receives 30 m. The design acknowledges overlap handling as unfinished,
but the parser/driver accept and execute it instead of refusing it.

**WHY IT MATTERS:** Valid recipes can poison DistMesh sizing or silently pay for
resolution the declared precedence was meant to avoid. The driver's `hmin =
min(target)/1.2` at line 221 is also not a true field floor when local base sizes
are finer than every target.

**SUGGESTED FIX:** Define width-zero behavior explicitly (normally retain the
ambient where no refinement is needed), validate finite positive field results,
and derive a real conservative seeding minimum. Implement the priority rule
across cores/transitions or reject unsupported overlaps during preflight.

## 13. P2 — “Effective gradation” omits the ambient gradient and is not a C4 certificate

**CLAIM:** The reported slope is not the gradient of the implemented field;
`within_c4` and the claimed hard impossibility threshold are misleading.

**EVIDENCE:** In a transition, `patch.py:724` implements
`h = target + (B(x)-target)*d(x)/W`. Its gradient is
`(B-target)*grad(d)/W + (d/W)*grad(B)`. `effective_gradation` at lines 747-755
measures only `(B-target)/W` at base nodes. Twenty Jacobi passes (`668-673`) do
not impose a spatial Lipschitz bound, particularly across short edges or at
extrapolation outside the base triangulation. On a falling non-monotone base
field, the second term can reverse the ramp; on a rising one it can steepen it.

The formula `1-1/(1+g)^2` at lines 733-738 assumes similar neighbouring triangles;
C4 is actually an area/altitude ratio across a shared edge. This contradicts the
qualification already made in design section 3.1. Even a reported OVER only
prints a message in the driver (`158-160`), it does not refuse.

**WHY IT MATTERS:** A report can understate the true gradient, while another
report can falsely declare a feasible mesh impossible. Large cores are not
inherently a division problem: distance is zero inside; width zero and spatially
varying ambient are the actual field-definition issues.

**SUGGESTED FIX:** Measure/bound the actual interpolated field gradient, including
its ambient term and exterior policy. Treat it as a meshing diagnostic, not a
necessary/sufficient C4 test. Gate C4 on actual final neighbouring triangle areas.

## 14. P1 — Geometry projection turns boxes into circles and discards polygon holes

**CLAIM:** The meshed core can differ substantially from the declared geometry.

**EVIDENCE:** `420_local_refine.py:94-99` guesses “circle” from the coefficient
of variation of **vertex radii**, then reconstructs a circle. A nearly square
bbox `[139.78,35.32,139.79,35.32816]` becomes a **257-coordinate circle instead of
a five-coordinate projected box** (appendix test). A square's vertices are all
approximately equidistant from its centroid, so this is not a rare numerical
accident. Both branches ignore `g.interiors`, although the shared geometry parser
supports GeoJSON holes. The preflight uses the original geometry while selection
and sizing use the altered one.

**WHY IT MATTERS:** This changes which water is refined, where the frozen zone
starts, and the meaning of the preflight measurements.

**SUGGESTED FIX:** Preserve geometry type/source parameters in `RefineRegion`;
project all rings for polygons and boxes, and reconstruct a circle only when the
recipe explicitly declared a circle, using its declared metric radius.

## 15. P1 — The executable never runs the claimed final quality gates

**CLAIM:** The driver writes a mesh after frozen verification even if seam repair
has stopped short of the quality requirements.

**EVIDENCE:** After `improve_patch`, `420_local_refine.py:333-384` only gates
coastline departure, `verify_patch['ok']`, and its boundary-loop assumptions.
There is no invocation of the 21-gate QA suite, no angle/area/valence rejection,
and no achieved minimum-altitude dt calculation anywhere in this driver.
`improve_patch` deliberately returns at a local optimum or round limit.

**WHY IT MATTERS:** Independent QA of the worked case does not make QA part of
the generator. Another fill can produce a successful exit and written artifact
with failed C1/C2/C4/C5 or other solver gates. The wrong valence guard makes this
particularly material.

**SUGGESTED FIX:** Run the established whole-mesh QA on the finished, serialized
mesh and record/reject failures before reporting successful completion. Include
actual achieved resolution and minimum-altitude dt in the report. Preserve a
failed candidate as an explicitly failed diagnostic artifact if useful.

## 16. P2 — The scale and dt estimates are empirical, not bounds on every generated mesh

**CLAIM:** The driver has no achieved-resolution check to support transferring
the 1.2 calibration to new patch geometries, and core-only preflight cannot bound
the global time step.

**EVIDENCE:** The driver divides sizing and hmin by 1.2 (`51,203,221`). In the
local OceanMesh checkout, `../oceanmesh/oceanmesh/mesh_generator.py:1629` computes
`L0 = hbars * L0mult * median(L)/median(hbars)`: uniform scaling of `fh` cancels
from this force expression. Scaling still changes seeding/hmin, so it can affect
achieved resolution, but it is not an algebraic edge-length guarantee. The
existing test at `tests/test_patch.py:341` checks only division by 1.2.

`refine.py:264-279` samples core depths and assumes equilateral target-sized
cells. It ignores possibly deeper transition cells and retained-domain limits.
The statement that a 30-30-120 cell “halves” the equilateral altitude depends on
which edge is held fixed: with two sides h, altitude is h/2 versus sqrt(3)h/2,
a factor 1/sqrt(3), not 1/2. Also, sampling can miss depth extrema. This is a
conditional estimate, not a bound on the global output step.

**WHY IT MATTERS:** These diagnostics can understate cost or misstate achievable
resolution for inputs unlike Futtsu, even after correcting the mesh defects.

**SUGGESTED FIX:** Keep calibration explicit and empirical; measure core edge
statistics on every result. Report achieved whole-mesh and per-zone altitude dt
with the controlling cell and depth convention. Fix the altitude comparison and
label the sampled preflight estimate accordingly.

## Specific attacks that did not establish additional defects

- For a valid manifold base triangulation, the narrow claim **free rim node =>
  physical boundary** is sound. A nonphysical removed-rim edge has a retained
  incident face, making both endpoints frozen. There is no valid counterexample
  to that implication in this review. The erroneous extrapolation is that every
  free run has frozen anchors; finding 6 supplies the counterexample.
- `boundary_rings` immediate-parent parity works for strict 3+ level nesting and
  disjoint simple rings (four levels were checked). Shared-node touching rings
  have degree four and are rejected/repaired. Geometric touching/crossing with
  distinct node IDs or a self-intersecting ring is not rejected by its degree
  test or checked for polygon validity. Such inputs are not valid conforming
  base meshes; they need explicit rejection rather than a nesting interpretation.
  Re-cut curves can independently introduce geometrical crossings (finding 1).
- `_convex_quad` tests one diagonal's side signs. For two correctly oriented,
  nonoverlapping incident triangles in a valid planar mesh, the other diagonal's
  side test follows from the input, so this is sufficient in that context. It
  does not validate already-overlapping input, nor check an existing replacement
  edge. Do not report that as a demonstrated fresh inversion on valid input.
- `_health` rejects nonpositive area in moved incident faces. With correct masks,
  no frozen node is moved, and flips cannot rewrite retained face IDs. I found
  no direct off-by-one in the zero-based stitch map. Nearest pfix matching is
  tolerance-based, not literally exact, but frozen coordinates are copied from
  the base; the all-free-ring interpretation is the concrete stitching failure.
- The 1 mm duplicate threshold rejects distinct points closer than 1 mm; it
  is not a general hanging-node/intersection test. The recipe has no minimum
  target bound compatible with that threshold. Document/configure this limit;
  do not confuse it with the exact frozen-identity check.

## Missing tests that explain these escapes

`tests/test_patch.py:151` exercises only two-level parity, not polygon construction
with nested rings. The area test at 283 uses one simple hole. There is no complete
physical-island identity fill, face-multiset mutation, written-file identity,
multiple/zero OBC, or recipe-to-driver geometry test. The new-depth test at 480
ends by asserting only **frozen** depths (line 494); it would not catch new-depth
corruption or stale depths after moves.

Most critically, `test_sliding_stays_on_the_boundary_chain` (548-561) builds a
`ring` but **never passes `slide_on=ring`**. `improve_patch` disables sliding when
no lines are supplied, so the `if moved.any()` assertion is vacuous. It covers
neither actual sliding nor preservation of original corners. The repair tests
check final minimum angle but not each operation's health and full valence.

## Verification and limits

Only small synthetic/unit probes ran on the login node; no mesh generation,
batch submission, production QA rerun, or modification of existing outputs.
No dependencies were installed. No commits or branches were created.

The appendix was saved temporarily as `/tmp/test_local_refine_review.py` and run:

```bash
MPLCONFIGDIR=/tmp/lr-review-mpl \
/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python \
-m pytest -q /tmp/test_local_refine_review.py --tb=no
```

Result: **9 failed, 1 warning in 1.95 s**, as expected for regression tests asserting
the missing behavior. The warning is `invalid value encountered in divide` at
`patch.py:723`. All nine failures are attributable to the defects above. Earlier
round-trip probe setup attempts failed first for omitted required boundary
arguments, then because the existing reader cannot handle a one-element block;
the reported precision measurement uses the normal multi-element grid and is
independent of those setup failures. The repository's full test suite and ruff
were not rerun: this change only adds this report.

## Runnable failing tests

Save the following block as `/tmp/test_local_refine_review.py` and run the command
above from this repository. The absolute ROOT is intentional for this workspace.
These are unmarked failures, not xfails; no production code is patched. The AST
extraction tests the driver's actual projection function without importing the
batch meshing driver or GPL dependencies.

```python
from pathlib import Path
import ast
import runpy
import numpy as np
import pytest
import shapely
from pyproj import Transformer
import fvcom_mesh_tools.patch as patch
from fvcom_mesh_tools.refine import RefineRegion
from fvcom_mesh_tools.io.fort14 import Fort14Mesh, read_fort14, write_fort14

ROOT = Path('/octfs/work/G16445/v61021/Github/fvcom-mesh-tools')
helpers = runpy.run_path(str(ROOT / 'tests/test_patch.py'))
grid = helpers['grid_mesh']


def identity():
    xy, tri = grid()
    depth = np.full(len(xy), 8.0)
    sel, rc, pn, pt = helpers['replay_patch'](
        xy, tri, shapely.Point(400, 400).buffer(150))
    out, faces, dep, mapping, _ = patch.stitch_patch(
        xy, tri, depth, sel, pn, pt, rc['pfix'], rc['pfix_base'])
    return xy, depth, tri, sel, out, faces, dep, mapping


def test_duplicate_retained_face_is_rejected():
    xy, depth, tri, sel, out, faces, dep, mapping = identity()
    faces = np.vstack([faces, faces[0]])
    assert not patch.verify_patch(
        xy, depth, tri, sel, out, faces, dep, mapping)['ok']


def test_frozen_means_exact():
    xy, depth, tri, sel, out, faces, dep, mapping = identity()
    v = mapping[sel.frozen_nodes[0]]
    out[v, 0] += 5e-7
    dep[v] += 5e-7
    assert not patch.verify_patch(
        xy, depth, tri, sel, out, faces, dep, mapping)['ok']


def test_nested_domain_preserves_parity():
    xy = np.vstack([np.array([[a,a],[b,a],[b,b],[a,b]], float)
                    for a,b in [(0,10),(1,9),(2,8),(3,7),(20,21)]])
    edges = np.array([[i, i//4*4+(i+1)%4] for i in range(len(xy))])
    assert patch.boundary_rings(xy, edges)[1] == [False, True, False, True, False]
    assert patch.hole_polygon(xy, edges).area == pytest.approx(57)


def test_physical_island_identity_fill():
    xy, tri = grid()
    c = xy[tri].mean(axis=1)
    tri = tri[~((c[:,0]>300)&(c[:,0]<400)&(c[:,1]>300)&(c[:,1]<400))]
    sel = patch.select_patch(xy, tri, shapely.Point(350,350).buffer(220))
    rc = patch.rim_constraints(xy, sel, size=1e9, coastline='preserve')
    taken = np.unique(tri[sel.removed])
    local = np.full(len(xy), -1)
    local[taken] = np.arange(len(taken))
    patch.stitch_patch(xy, tri, np.ones(len(xy)), sel, xy[taken],
                       local[tri[sel.removed]], rc['pfix'], rc['pfix_base'])


def test_zero_width_is_finite():
    xy, tri = grid()
    field = patch.patch_sizing(xy, tri,
        [(shapely.Point(400,400).buffer(100), 200, 0)], distmesh_scale=1)
    assert np.isfinite(field([[400,400]])).all()


def test_every_accepted_flip_respects_global_valence(monkeypatch):
    xy, tri, _ = helpers['perturbed_patch']()
    original = patch._valence_ok
    accepted = []
    def checked(tri, faces, verts, limit):
        result = original(tri, faces, verts, limit)
        if result:
            accepted.append(max(int((tri == v).sum()) for v in verts))
        return result
    monkeypatch.setattr(patch, '_valence_ok', checked)
    patch.improve_patch(xy, tri, np.zeros(len(xy), bool),
                        np.ones(len(tri), bool), rounds=1, max_valence=8)
    assert max(accepted) <= 8


def test_frozen_depths_round_trip(tmp_path):
    xy, tri = grid()
    depth = np.full(len(xy), 5.12345678912345)
    mesh = Fort14Mesh(title='test', nodes=xy, elements=tri, depths=depth,
                     open_boundaries=[], land_boundaries=[])
    path = tmp_path / 'mesh.14'
    write_fort14(mesh, path)
    assert np.array_equal(read_fort14(path).depths, depth)


def test_bbox_does_not_become_circle():
    tree = ast.parse((ROOT/'notebooks/420_local_refine.py').read_text())
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
              and n.name == 'region_in_metres')
    env = dict(np=np, shapely=shapely,
               to_m=Transformer.from_crs(4326, 32654, always_xy=True))
    exec(compile(ast.Module(body=[fn], type_ignores=[]), '<driver>', 'exec'), env)
    region = RefineRegion(dict(name='square', target_h_m=30,
                              geometry={'bbox':[139.78,35.32,139.79,35.32816]}))
    projected = env['region_in_metres'](region)
    assert len(projected.exterior.coords) == 5


def test_retained_faces_must_be_edge_connected():
    xy = np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0]], float)
    tri = np.array([[0,1,2],[0,2,3],[0,3,4],[0,4,5]])
    with pytest.raises(ValueError, match='splits the retained mesh'):
        patch.select_patch(xy, tri, shapely.Point(0, 2/3).buffer(.4))
```
