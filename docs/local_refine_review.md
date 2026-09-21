# Design review: local refinement

Reviewed 2026-09-22. Scope: the proposed operation in `local_refine.md`, especially
section 5, rather than an audit of a completed generator. Only this report was
written. Findings below are ordered by the cost of ignoring them. Measurements
use the actual file named by `recipes/refine/futtsu_nori.yaml`; analytical
estimates are explicitly distinguished from measurements.

**Recommendation:** do not implement section 5 literally. The first recipe is
ineligible under its own coastline rule, its advertised timestep uses the wrong
geometric measure, and its base file does not carry the assumed production
bathymetry. A local operation remains defensible when preserving the exterior
is a hard requirement, but it needs an explicit constrained interface, a
topological footprint, and a feasibility/rejection contract. For this first
fishery, a global sizing-region rebuild is the lower-risk route if exterior
identity can be relaxed. That rebuild still needs correct timestep and
bathymetry treatment.

## 1. Blocker: the worked example already reaches the coast

The document and recipe establish only that the **core** is clear of land.
Section 5 requires the **hole** to be clear. These are substantially different
tests:

```
W = (350 - 30) / 0.165 = 1939.394 m
outer radius = 300 + W = 2239.394 m
```

Read-only selection on the supplied base, in EPSG:32654, gives:

| Measured quantity | Result |
|---|---:|
| Base nodes / elements | 4,734 / 8,252 |
| Centre to closest mesh boundary segment | 1,014.459 m |
| Elements selected by the proposed centroid rule | 173 |
| Area of those elements | 14.998625 km² |
| Edges on the selected submesh boundary | 47 |
| Interface edges shared with retained elements | 35 |
| Selected physical boundary edges | 12 |
| Selected land-boundary / OBC nodes | 13 / 0 |

Thus `touch_coast: false` must reject this recipe. This is established from the
mesh's own physical boundary, without relying on a possibly different land
dataset. The 654 m distance from the tidal-flat edge is not distance from the
coastline and does not establish eligibility.

Even the optimistic smooth-disc clearance permits only
`W < 1014.459 - 300 = 714.459 m`. Reaching 350 m in that space requires
`g > 320/714.459 = 0.44789`, before any element-selection buffer. Merely changing
`g` to 0.33 is insufficient. Raising the target at the original gradation
would require roughly `350 - 0.165*714.459 = 232.11 m`, defeating the objective.

`preflight()` samples only core points for land, accepts `land=None` despite its
docstring requiring land when coast touching is forbidden, and has no mesh or
OBC argument. Its success therefore cannot certify the documented selection
constraints. A 64-by-64 lattice also cannot prove that a thin land sliver does
not intersect a core. Use geometric intersection for that condition, then
inspect actual selected faces and protected node stars.

**Required decision:** move the example, explicitly support a coastline patch,
or use a global rebuild. Do not silently turn on coastline modification after
selection discovers the conflict.

## 2. Blocker: 30 m does not deliver the stated timestep

The certified 11.9 s figure uses **minimum triangle altitude / wave speed**;
`outputs/verify_409.115302/summary.txt` states this explicitly.
`coast_fit._implied_dt()` guards both shortest edge and minimum altitude.
`refine.preflight()` uses the requested edge length, and `qa.py`'s timestep
check uses shortest edge. These are different quantities.

At the document's maximum core depth, 4.15 m:

```
c = sqrt(9.81 * 4.15)
30/c                         = 4.701784 s  (document/preflight)
(sqrt(3)/2 * 30)/c            = 4.071865 s  (equilateral altitude)
4.5*c / (sqrt(3)/2)           = 33.154344 m (required equilateral edge)
```

The ideal equilateral 30 m patch already misses 4.5 s. At the recipe comment's
4.11 m, the corresponding numbers are approximately 4.72 s and 4.09 s; that
discrepancy between recipe and document does not rescue it. At 3 m the ideal
altitude estimate is about 4.79 s, so the deeper part matters. A nominal 30 m
statistical target might coexist with larger triangles over the deeper water;
the proposed constant core field does not promise that adaptation.

Passing C1 does not fix this: a 30°–75°–75° triangle with shortest edge 30 m
has minimum altitude `30*sin(75°)`, and a 30°–30°–120° triangle with two
30 m sides has minimum altitude 15 m. Both meet C1/C2; the latter allows
only 2.35 s at 4.15 m. A quality gate is not a timestep guarantee.

Check actual final altitudes and maximum vertex depths throughout the patch,
seam, and unchanged base. Sampling only the core misses a deeper transition
cell and depth changes introduced by smoothing. Give the recipe a named
timestep metric and Courant/safety convention; the current estimate implicitly
uses Cr=1, without velocity or operational safety allowance. This is a mesh
diagnostic, not a proof of an FVCOM operational timestep.

The cost statement also overpromises: `11.9/4.071865 = 2.92` times as many
external steps in the ideal example, before the increased element count.
Runtime is not determined entirely by core depth. The size distribution,
element shape, added work, transition depths, and solver configuration matter.
The prose's “largest target” should be “smallest permissible target.”

## 3. Blocker: the base file contradicts the bathymetry premise

Actual depths in the recipe's fort.14 are:

| Measured quantity | Result |
|---|---:|
| Minimum / maximum depth | 2.0 / 735.315384 m |
| Nodes below 3 m | 1,160 |
| Maximum edge r-factor | 0.906991 |
| Edges with r > 0.2 + 1e-9 | 2,258 |
| Depths of existing nodes inside the 300 m core | 3.649794, 5.980503 m |
| Selected-element vertex depth range | 2.0–32.000004 m |
| Interface node depth range | 2.0–23.094765 m |
| Interface nodes below 3 m | 5 |

These measurements do not identify the depths' source, but they disprove that
this file already implements the claimed 3 m floor, 300 m cap, and r<=0.2
recipe. A regional M7001 insertion here creates a mixed-bathymetry model.
Holding the exterior fixed cannot simultaneously turn this into a globally
production-compliant depth field.

The current default QA run **does pass 21/21** on this base. That is not a
contradiction: its default minimum depth is 2 m, and it has no r-factor or
300 m maximum-depth gate. “21 gates” is not the complete acceptance policy.

Choose a baseline whose depths already have the intended provenance, or make
an explicit one-time global depth conversion before establishing the frozen
baseline. Alternatively declare the mixed-source depth operation intentionally
and assess it as such. Include the depth artifact/checksum in provenance;
changing it implicitly would compromise the claimed controlled experiment.

## 4. Blocker: fixed vertices are not a conforming interface

Section 5 treats a list of `pfix` points as a boundary-edge constraint. The
local oceanmesh implementation distinguishes them:

* `mesh_generator.py:1077` constructs segment constraints only for `egfix`.
* At retriangulation, it uses CDT when those constraints exist and ordinary
  Delaunay otherwise, then removes triangles by the domain predicate.
* `pfix` coordinates are restored and forces on them are suppressed. This
  does not require the prescribed edges between them to appear.

On a concave or stair-stepped rim, an unconstrained triangulation can connect
across an indentation. Centroid-based domain rejection does not prove that
surviving triangles lie entirely in the hole, nor that all required boundary
segments survive. Deduplicating coordinates cannot repair a missing interface
edge, an overlap, a T-junction, or a crack. These failures can remain even
when the global mesh graph is connected.

Use the union of removed triangles as the meshing domain, extract its exact
boundary segments, and constrain every interface segment. The fork already
has `pfix+egfix`, used by notebook 325 for the OBC. Verify segment survival
after every cleanup stage, not just immediately after triangulation. No new
vertex may split a frozen interface segment unless the retained element is
also changed, which expands the authorized footprint.

Several hundred fixed points are not inherently fatal. The problem is the
shape and spacing constraints they impose and the absence of exterior forces
and topology in a patch-only solve. In this particular selection there are
only 47 selected-boundary edges, not hundreds. Expect density adjustment,
retriangulation, and stalled quality near awkward corners; do not promise a
quality result from convergence of the force iteration.

The statement “two rim nodes closer than the local target size collapse” is
not a mathematical spacing threshold. In the current fork, initial points
prepend all pfix, deduplication uses a machine-scale tolerance
`1024*eps*extent`, and the named error detects a non-injective nearest-vertex
mapping after pruning. Lost points or topology are also causes. Closeness to
the target alone does not prove that distinct fixed vertices must merge.
Also require a small correspondence distance: nearest-neighbour uniqueness
alone does not prove that a supposedly surviving vertex is the right one.

The original DistMesh formulation describes a relative sizing function and
fixed node positions, not an exact edge-length or segment-preservation
guarantee. See [Persson's thesis](https://persson.berkeley.edu/thesis/persson-thesis.pdf).
The fork-specific conclusions above come from its local source, not an
assumption that it behaves exactly like the original MATLAB implementation.

## 5. High cost: the rim does not match a 350 m ambient field

The measured interface edge lengths are **286.70 / 485.74 / 819.25 m**
(minimum / median / maximum). “Already equals the local base-mesh edge length”
is tautological; it does not establish compatibility with the proposed
clipped sizing field. The recipe supplies no ambient field, and the preflight
accepts only a scalar ambient size.

With a field clipped to 350 m, an immutable 819 m rim edge is already 2.34
times the nominal local target. For a patch-side triangle on that edge to
have all angles >=30°, its opposite vertex must be at least
`819.25*tan(30°)/2 = 236.50 m` from the edge; at least one other side must be
`819.25217/sqrt(3) = 472.995 m` or longer. It cannot be an approximately
equilateral 350 m rim cell. A field target is soft enough to permit this,
but the claimed matching argument is false.

The QA constraints depend on the frozen exterior:

* C4 is `(A_large-A_small)/A_large <= 0.5`, hence an area ratio <=2.
  Across a shared rim edge this is an **altitude ratio <=2**, not a universal
  edge-length ratio <=sqrt(2). The latter assumes geometrically similar
  triangles. A gradient bound alone cannot certify C4.
* A seam vertex's final valence counts retained and replacement neighbours.
  Patch-only valence <=8 is insufficient. With six retained neighbours and
  both rim neighbours already counted, only two new neighbours are available.
* Patch corner angle sectors are fixed by the cut. If their total angle is
  theta, the number of incident new triangles k must satisfy
  `30*k <= theta <= 130*k`, as well as the full-mesh valence budget.
* Merely avoiding intersection with an OBC arc can still modify an element
  incident to an OBC node. That can change orthogonality while every OBC
  coordinate remains fixed. Protect the OBC node stars/needed guard band.

Define ambient reconstruction and interface compatibility from actual edge
lengths, retained triangle areas, angles, and node stars. Reserve a collar in
which the old triangulation can be retained or used as a starting condition.
Allow a deterministic, bounded footprint expansion when the interface is
incompatible, or reject. A final failed QA run is a diagnostic; without a
repair/reselection policy it is not an algorithm for obtaining a usable mesh.

## 6. High cost: selection and the frozen contract disagree

Centroid selection removes whole triangles that straddle the analytic
footprint. In the measured selection, **22 selected vertices lie outside the
2,239.394 m disc**, and cut-boundary vertices extend to **2,614.767 m** from
the centre. Conversely, centroid-outside triangles can intersect the core
and leave coarse intrusions for general small or narrow regions.

A centroid cut can have multiple components, pinched vertices, holes, or
retained triangle islands. Its boundary need not be one closed chain:

* Separate face components require separate fills or an explicitly merged
  footprint. Components touching only at a vertex do not make a manifold
  patch domain.
* A retained mesh island inside the selected area is a hole whose nodes and
  faces must survive. A physical island is also a hole, with a land-boundary
  identity. Both require all boundary loops and containment/orientation data.
* A coast-reaching selection has an interface chain joined to a physical
  boundary chain. It is not a water-only closed interface. Coast replacement
  needs defined endpoints, protected endpoint stars, and boundary-list edits.
* Buffering a polygon with holes can shrink or erase those holes. Decide
  whether a hole means “not core” or “protected from all modification”; these
  are different contracts.

The exterior-identity requirement is achievable for an accepted patch, but
**not as simultaneous guarantees of arbitrary target feasibility, an exact
analytic affected region, and immutable surroundings**. Make success
conditional. Define the footprint as an explicit set of removable faces,
with shared interface vertices immutable. Report its actual geometry and
departure from the requested envelope before filling. Any enlargement must
respect a declared maximum change region; never silently relax the promise.

Specify frozen coordinates, depths, retained face connectivity/orientation,
boundary membership/order/type, and node identity under an explicit ID map.
Incident replacement faces at interface nodes may change; retained faces may
not. Renumbering is compatible with topological identity under a map, not
literal equality of the serialized connectivity table. For downstream
node-indexed forcing/restart data, a map or stable IDs is essential. Exterior
mesh identity also does not imply unchanged simulated flow there.

`frozen_changes()` is insufficient even as the eventual verifier:

* It requires equal-shaped arrays in row correspondence, whereas this
  operation inserts/deletes nodes and renumbers them.
* It does not inspect connectivity, depths, boundaries, node loss, or extra
  exterior faces. Its caller can mark away a violation with the affected mask.
* It ignores columns after xy and has no CRS check. Its tolerance is called
  metres regardless of input units, and 1e-6 tolerance is not “identical.”
* A frozen coordinate replaced with NaN currently returns `ok=True`, because
  `NaN > tol` is false. This was reproduced with a one-node input.

Verification should use a base hash plus immutable face/node ID sets computed
before meshing, exact copies of frozen payloads, and explicit old-to-new maps.
Compare retained oriented triangles after mapping, boundary sequences, depths,
and coordinates. Check finite values, injective maps, every interface edge's
two-face incidence, all other edge incidences, no unreferenced or duplicate
vertices, and geometric coverage of exactly the removed-face union (no gaps
or overlaps). Re-read the output and repeat the invariants. Run full-mesh QA,
the connectivity/width checks, explicit r-factor, and altitude timestep checks
after all finishing passes. This is substantially stronger than a count of
moved coordinates.

## 7. High cost: patch bathymetry is a constrained solve, not a sliced limiter

The proposed idea can be sound **if every seam edge is included, the rim is
fixed in depth, and the constraints are feasible**. Define mutable depth
nodes separately from geometry nodes. All nodes incident to retained faces,
including the rim, should retain their depths if the exterior model is frozen.
If the rim is instead mutable, include every rim-to-exterior edge in the
limiter and admit that retained elements' depth fields changed.

For positive depths and rmax=0.2:

```
|hi-hj|/(hi+hj) <= 0.2  <=>  2/3 <= hi/hj <= 3/2
```

Thus a mutable node next to a frozen depth H must lie in
`[H/1.5, 1.5*H]`, intersected with [3,300]. Constraints from multiple frozen
neighbours must have a common solution. Example: a new node adjacent to
fixed depths 3 and 12 m needs both `[2,4.5]` and `[8,18]`, an empty
intersection. More generally a k-edge path between fixed depths requires
`Hmax/Hmin <= 1.5**k`. Remeshing can create short graph paths that invalidate
otherwise compatible fixed data. With log depths, these are difference
constraints `|log hi-log hj| <= log(1.5)`; feasibility can be checked rather
than discovered through 2,000 ineffective iterations.

The existing `rfactor_smooth()` moves both endpoints of a bad edge, averages
incident corrections, and has no fixed mask. Running it on the patch and
restoring rim depths afterward reintroduces violations. For example, depths
3 and 12 are moved to 6 and 9 on an isolated edge (r=0.2); restoring the fixed
3 produces r=0.5. Simply dropping edges with fixed endpoints misses the seam.
Masking half the correction is also a different algorithm requiring explicit
convergence and residual checks. The current function returns its residual
at the iteration limit; callers must reject an excessive residual.

Use a constrained objective that stays close to the sampled M7001 depths,
fixed rim values, explicit bounds, and all relevant graph-edge inequalities.
Report maximum r separately in the patch, across the interface, and globally;
report depth adjustments and integrated volume change. Recheck timestep after
depth changes. A uniform floor/cap does not itself worsen a valid positive
depth ratio, but applying different depth policies on opposite sides does not
give that guarantee. In this actual base, five rim nodes are 2 m: preserving
them conflicts immediately with applying a 3 m minimum to the entire rim.

Frozen-boundary smoothing is not generally identical to a global production
smoothing run and should not be described as such. The production survey grid
is approximately 180 m: 30 m triangles give about six edges per source-grid
spacing, not new 30 m bathymetric information. Refinement can still help flow
or aquaculture representations, but depth fidelity alone is not a justification.

## 8. Medium/high cost: target and recipe semantics need decisions now

Keep an achieved-resolution quantity in the user interface, but replace
“achieved edge length approximately target” with a measurable acceptance rule:
which edges count (e.g. edges wholly inside the wet core), a statistic and
tolerance, a maximum unresolved fraction, and separately a hard timestep
test. Median alone can hide a coarse strip through the fishery. Publish the
distribution and spatial coverage. Do not treat target_h_m as a minimum
edge/altitude or as an exact length for every edge.

The factor 1.2 is calibration, not a backend law. This fork sets `L0mult=1.1`
when pfix is present (`mesh_generator.py:981`), and computes
`L0 = hbars*L0mult*median(L)/median(hbars)` (`:1629`). Patch extent, constrained
edge population, initialization, and density control influence the result.
Notebook 325 scales distance-derived size and slope, not every field term.
A patch must calibrate and report its achieved statistics; blindly feeding 30
to the backend, or unconditionally dividing all terms by 1.2, does not implement
the promised contract.

Decide the following schema/manifest semantics before adding generators:

| Missing or ambiguous contract | Concrete decision needed |
|---|---|
| Version and baseline identity | Schema version; base mesh and depth hashes; backend/version/seed and resolved settings in result manifest. |
| CRS and distances | Read/declare mesh CRS; transform geometry and buffer in metres. EPSG:32654 is specific to this case, not every permitted lon/lat polygon. |
| Ambient size | Saved field with provenance, or documented reconstruction from the base; a coastal scalar does not describe this rim. |
| Shoreline source | Reference the processed generator coastline, including channel edits, rather than infer it from raw OSM or a filename elsewhere. |
| Coast policy | Separate “core may overlap land,” “transition may reach coast,” and “coast may be replaced.” A boolean currently conflates these. |
| Footprint limit | Authorize/reject cut expansion, exclusions, protected channels/islands, and OBC star protection. |
| Resolution/CFL | Target statistics/tolerances; achieved versus field gradation; altitude metric and Cr convention. |
| Depth policy | Source/datum, floor/cap/rmax, immutable rim, residual and feasibility behavior. |
| Transition | Define extra width behavior and whether zero is allowed for a no-op. The loader currently rejects zero even when required width is zero. |
| Unsupported input | Reject multiple regions/coast touching until implemented; priority validation alone is not support. |

Some of these can be fixed documented defaults or generated provenance rather
than user knobs. Also decide whether target>=ambient is a no-op, coarsening, or
an error: the current width helper returns zero for both equality and coarsening.
Do not inherit sizing.py's coarsening semantics accidentally. Validate direct
API arguments as well as YAML: finite positive gradation/depth arrays and
valid shapes are part of meaningful preflight results.

## 9. Medium cost: the transition and count arithmetic is inconsistent

For the stated **350 m** ambient and g=0.165, the correct width table is:

| Target (m) | W (m) | W / 300 | Diameter / target |
|---:|---:|---:|---:|
| 15 | 2,030.303 | 6.768 | 40 |
| 30 | 1,939.394 | 6.465 | 20 |
| 50 | 1,818.182 | 6.061 | 12 |
| 100 | 1,515.152 | 5.051 | 6 |

The document's 15, 50, and 100 m widths instead approximately use **330 m**
ambient; its 30 m row uses 350 m. “Elements across” is undefined. If it means
core diameter divided by target, use the last column. If it means radial
transition size intervals, the continuous estimate is
`integral(dr/h(r)) = log(350/target)/0.165`, giving 19.090, 14.889, 11.793,
and 7.593 respectively. Neither interpretation reproduces the printed column;
triangle row counts additionally depend on orientation/altitude.

`preflight()` buffers lon/lat by `W/111000` before metric conversion. At
35.3228° this supplies only approximately `W*cos(latitude) = 1582 m` of
east-west transition, versus 1939 m north-south. Project first, then buffer.
The reported affected area is 13.250335 km²; a true circular footprint has
`pi*(300+1939.394)^2 = approximately 15.755 km²`.

The midpoint-size count is not an integration of inverse squared size.
For a circular wet domain, equilateral cells, a=√3/4, core radius R, target h,
ambient H, and W=(H-h)/g, the actual continuous estimate is:

```
Ntransition = (2*pi/a) * integral[0,W] (R+s)/(h+g*s)^2 ds
            = (2*pi/a) * [log(H/h)/g^2
                         + (R-h/g)/g * (1/h-1/H)]
```

With R=300, h=30, H=350, g=0.165:

* Core: 725.520 elements.
* Transition: **1,626.135**, versus the reported 829.569.
* Uniform ambient elements replaced in the full disc: 297.012.
* Net added: **2,054.643 = 24.90% of 8,252**, versus 1,305 =15.82%.

These are corrected **idealized estimates**, not predictions for the
coast-clipped, variable-density real patch. Actual centroid selection removes
173 elements, not 250 or 297. Integrate over the actual wet footprint and
spatial field when estimating a real candidate. Enlarging the core also
increases transition circumference and thus transition cost even if W stays
constant; “comparatively cheap” needs quantified marginal costs.

## Answers to all five section 7 questions

1. **Geometry files:** useful, but not the first blocker. Support CRS, layer,
   deterministic feature/attribute selection, multiple selected features,
   Polygon/MultiPolygon and holes, and a content hash. Resolve paths against
   the recipe. A shapefile needs its sidecars. Define whether holes are hard
   exclusions. File loading does not solve topology at the cut.
2. **Several regions:** combine transition footprints before selecting faces
   and solve connected patches jointly. Define priority over cores versus
   transitions; a high-priority coarse transition must not silently erase a
   fine core unless that is explicitly intended. Gradation can propagate fine
   sizes beyond either nominal footprint. A minimum-size envelope is a
   sensible refinement policy, but differs from unconditional “higher wins.”
   Reject unsupported multi-region inputs initially. Check permutation
   invariance and disjoint, nested, overlapping-transition, and conflicting
   core cases when implemented.
3. **Bathymetry:** conditionally sound as the constrained solve in finding 7;
   unsound as a call to the existing unrestricted limiter followed by restoring
   frozen depths. Fix the baseline mismatch first.
4. **Repeated refinement:** compare against the immediate input with a manifest
   linking its hash to ancestors and storing ID maps/changed-face sets. To
   claim equality to the original outside all patches, compose those maps and
   compare against the original there. Reject/coarsen/no-op behavior for a
   previously finer region must be explicit. Sequential patches can be order
   dependent; do not promise equivalence to a joint run.
5. **Coast fit:** unnecessary for a truly interior patch. For a coast patch,
   operate on the stitched mesh with only authorized coast nodes movable,
   fixed interface/endpoints and protected exterior stars. A cropped patch
   alone loses exterior C4 and timestep context. Existing `fixed` handling
   expands to neighbours by default, so it can freeze the endpoint collar;
   account for that in the generator. Fitting cannot guarantee exact coastline
   placement because guards may refuse moves. Define whether depths are
   resampled after coordinates move; then recheck r-factor, timestep, topology,
   and all final invariants. Never run an unrestricted global fit afterward.

## A minimum defensible implementation and stop/go test

Before optimizing a fill backend, implement candidate extraction and an
interface feasibility report. Restrict the first supported case to one
interior manifold patch; reject coasts, OBC stars, nonmanifold cuts, and
unsupported holes explicitly. Save an immutable baseline snapshot and
topological footprint. Feed exact interface segments to a constrained fill,
preserve an exterior-aware collar, stitch by IDs, solve constrained depths,
and verify the full mesh. Fail without publishing a replacement on any
invariant failure. Bound retries and footprint expansion.

The acceptance suite should include a concave rim where point constraints
alone lose an edge, alternating selected triangles at a pinched vertex, a
retained interior island, an OBC-adjacent face whose centroid misses the arc,
a full-valence rim node, a C4 mismatch to an exterior triangle, a nonfeasible
fixed-depth example, and insertion/renumbering in the frozen verifier. For the
first production candidate, report actual rim distributions, core resolution
coverage, selected area/count, seam quality, depths/r-factor, both timestep
measures, and preserved exterior identity. Existing specification unit tests
do not demonstrate any of these generator properties.

If that engineering effort is justified only by “the patch ought to be
cheaper,” choose the existing sizing-region rebuild. It already supports the
field edit and avoids an artificial immutable interface, though it still
requires full QA and calibration. Local refinement is justified by a concrete
need for exterior mesh identity, not by the size of the 300 m core: the
required influence radius is kilometres. Neither route restores 4.5 s for
uniform equilateral 30 m cells at 4.15 m depth.

## Verification record and limitations

Read the requested design, recipe, specification tests, sizing documentation
and implementation, notebook 325, coast-fit documentation and implementation,
M7001 module, and QA implementation. Also read local oceanmesh generator
source as a read-only dependency, the fort.14 reader, and the supplied mesh
summary. No oceanmesh import or mesh generation was performed.

All Python commands used the existing environment:

```bash
. /octfs/work/G16445/v61021/miniforge3/etc/profile.d/conda.sh
conda activate oceanmesh-bench
PYTHONDONTWRITEBYTECODE=1 pytest -q -p no:cacheprovider tests/test_refine.py
```

Result: **24 passed in 0.87 s**. No new dependencies or test files.

Read-only inline Python diagnostics used `read_fort14`, a pyproj
4326→32654 transform of the declared centre, Euclidean centroid distances,
`coast_fit.boundary_edges`, `dem.m7001.node_edges`, and
`coast_fit._implied_dt`. The essential selection calculation was:

```python
m = read_fort14('outputs/verify_409.115302/fit/sample_repro_final.14')
p, t = m.nodes, m.elements
c = np.array(Transformer.from_crs(4326, 32654, always_xy=True)
             .transform(139.7881, 35.3228))
selected = np.linalg.norm(p[t].mean(axis=1) - c, axis=1) < 300 + 320/.165
cut = boundary_edges(t[selected])
physical_edges = set(map(tuple, np.sort(boundary_edges(t), axis=1)))
seam = np.array([e for e in cut if tuple(sorted(e)) not in physical_edges])
ei, ej = node_edges(t)
r = abs(m.depths[ei]-m.depths[ej]) / (m.depths[ei]+m.depths[ej])
```

`run_qa(m, coords='metric')` returned **passed=True, 0 failed, 21 gates**.
Measured base timestep minima were **16.597904 s by edge** and
**11.855275 s by altitude**, reinforcing the metric mismatch. A separate
`preflight()` call with constant 4.15 m depth and ambient 350 m reproduced
the advertised 4.701784 s, 725.447 core, 829.569 transition, 249.799 replaced,
and 1,305.217 added estimates; this reproduces the estimator, not M7001 data.

Two diagnostic attempts failed and were corrected: NumPy rejected `np.cross`
on 2-D vectors, replaced with the explicit 2-D determinant for selected area;
QA rejected `coords='projected'`, rerun using the documented `'metric'` value.
These were review-command errors, not findings against the proposed generator.

No patch was generated, no figure batch or smoothing run was launched, and no
NQSV command was submitted. The reported production M7001 core depth range
and tidal-flat provenance were not independently regenerated. The measured
coast conflict is with the existing mesh shoreline; a processed-land check is
still required for any revised candidate. There are no measured final-patch
QA gates, node/element counts, or timestep because the generator does not yet
exist. The final-patch behavior discussed above is a set of design obligations
and analytical counterexamples, not a claimed failed meshing experiment.
