# Third review: is local refinement finished?

Reviewed at **30e0e70**, 2026-09-22. Read both implementation reviews and
`docs/local_refine.md` revision 5. This report reviews readiness, not another
round of the already-open defects. No generator or mesh files were changed.

**The Futtsu demonstration is finished; a generally usable generator is not.**
The current driver enforces substantially more than a plausible-looking mesh:
it checks the intended boundary, frozen geometry/depths, retained faces and
interface, then checks the serialized candidate and runs QA. I found no new
structural run-stopper in the delivered mesh. That does not qualify its known
SRTM15 bathymetry for the intended production model.

The next work should make acceptance structural and establish the intended
production baseline, then support actual user geometries. Another seed sweep
or another repair heuristic is not the highest-value next step.

## 1. Contract: correct current call sequence, incomplete public guarantee

### P1 for a general-purpose API — success can mean incomplete verification

**CLAIM:** `verify_patch['ok']` does not imply that the patch contract was
checked. The present driver supplies the missing evidence correctly; a new
caller or refactor can omit it and silently obtain success. This is the API
readiness consequence of the coverage fix, not a claim that that fix fails
when used as intended.

**EVIDENCE:** `src/fvcom_mesh_tools/patch.py:1069` defaults
`expected_boundary=None`; lines 1162–1165 condition the boundary comparison on
its presence; lines 1196–1199 do not require `boundary_checked` for `ok`.
The appendix's intact identity-fill test returns **ok=True,
boundary_checked=False**. Its control supplies the boundary, accepts the
intact fill and rejects deletion of an interior face with **three unexpected
boundary edges**.

`notebooks/420_local_refine.py:310` constructs the boundary before repair;
lines 415–418 and 479–485 pass it into both verifications. Those calls hold.
Moving the construction after sliding ordinarily raises on a missing coordinate
(`patch.py:1040–1042`), rather than silently disabling coverage. The silent
failure is omitting the evidence, not the timing error itself.

**WHY IT MATTERS:** An optional check plus a Boolean success result is too easy
to mistake for a certificate. `boundary_checked=False` is useful diagnostics,
but is not a safe default for a method advertised as verifying the contract.
The mapping's lifetime also depends on knowing that repair preserves node IDs.

**SUGGESTED FIX:** Make the invariant structural:

> Every candidate eligible for acceptance carries the intended boundary and
> its persistent vertex identities, derived from the declared cut/rim, and
> final success requires checking that evidence on the delivered candidate.

Have stitching return a candidate containing `pfix_to_output`, expected boundary
edges and curve ownership alongside arrays and the base node map. Stitching
already has `idx` and `patch_map` (`patch.py:941–960`); their composition gives
that map without another coordinate lookup. Repair should preserve those IDs
or explicitly remap the evidence. Require this contract object in verification;
missing evidence must raise or make overall acceptance false. A separate
`verify_frozen_zone` can remain a partial diagnostic. Keep intended evidence
independent of the candidate's observed boundary: deriving both from the output
would make the coverage check tautological. Carry the same evidence through
serialization, including the boundary semantics whose preservation is already
an open item.

Boundary-edge equality is valuable under valid planar-triangulation assumptions;
it is not, alone, a general geometric proof against crossings, overlapping
components or incorrect coastline location. Do not rename it an unconditional
domain certificate. The delivered mesh's independent overlap check below passed.

## 2. What fails first away from Futtsu?

This is a ranked forecast from the implementation and existing synthetic tests,
not seven new production runs. All rows assume the driver's current metric
frame, **EPSG:32654** (`420_local_refine.py:53,93–94`). A base in another CRS first
needs an explicit CRS contract; fort.14 coordinates alone do not supply one.

| Rank | Change | First practical limitation / readiness |
| --- | --- | --- |
| 1 | Several open boundaries | Certain refusal in `serialise` at lines 439–445, after expensive filling/repair. No silent merge in the current driver. Move this refusal to preflight immediately; implement multiple arcs only when a needed base requires them. |
| 2 | Narrow channel | The buffered cut may span both banks and disconnect the retained domain, which selection correctly refuses. If it survives, bank correspondence, coastline tolerance relative to channel width and too few cells across the channel become the risks. A 200 m tolerance is not a channel-preservation specification. Global connectivity and angle gates do not certify the intended hydraulic cross-section. Channel w/h is advisory by default (`qa.py:895–917`). |
| 3 | Region touching/enclosing an island | Basic nesting, stitching and island-curve registration now have passing regressions. But a wholly free island retains its original vertices without local-size subdivision (`patch.py:640–658`): a coarse island rim can defeat fine surrounding resolution/QA at every seed. A partially touched island also exercises the already-open cyclic-source/provenance and tolerance issues. I found no reason to reopen the fixed island-ID or missing-curve bugs. |
| 4 | Two regions | Disjoint regions already feed one union cut and one minimum sizing field; two regions are not categorically unsupported. Joint growth may sever retained water. Overlap with conflicting priorities silently follows finest-wins, the recorded open defect. Preflight-refuse unsupported precedence now; do not spend seed-search time on the wrong field. |
| 5 | Much larger region | More fine cells and a larger cut cost memory and repair time, and eventually hit OBC/connectivity refusals. Core cell demand scales approximately as area/target². Widely separated regions also enlarge the common fill bounding box at the finest hmin. No reason to expect seed retries to solve those resource or feasibility limits. No scaling benchmark was run. |
| 6 | Polygon fishery | Inline GeoJSON Polygon, including holes, already parses and projects; file loading is the missing convenience. An ordinary compact wet polygon is the lowest-risk extension. Thin, highly concave or fragmented features need explicit sample/coverage and resolution checks: preflight uses a fixed bounding-box lattice (`refine.py:252–259`), not a geometry-adaptive sample. Polygon holes exclude the core target, not necessarily all transition refinement. MultiPolygon is not supported by the parser. |
| 7 | Much deeper site | Depth alone does not break the geometric construction. Its first consequence is global time-step cost, approximately proportional to sqrt(depth) at fixed cell geometry, and the suitability of the inherited depth field. Transition depths can control the step. The design intentionally makes dt advisory; preserve that decision and report the achieved cost. |

Ranks describe likely failure of reuse, not universal severity. A long thin
polygon is effectively the channel case; a distant two-region request is also
the large-bounding-box case. No real polygon/channel/island multiseed build was
performed in this review.

## 3. Seed search: reasonable engineering with a limited acceptance predicate

**First passing seed is defensible.** This is randomized numerical construction
with deterministic acceptance checks, not repeated statistical testing where
trying more seeds invalidates a significance threshold. The present state
machine retries candidate `ValueError`s, records attempts, checks disk QA and
rewrites the selected candidate. Its regression tests pass. No additional
seed-state defect was established here.

A passing mesh can still be wrong when the predicates do not express the
request: incorrect bank/coastline geometry, unresolved core size, ignored region
priority, unsuitable depth field, or changed boundary semantics. These are the
already-recorded open requirements, not stochastic failures. Validating a rim
against itself can also consistently certify the wrong intended domain. More
seeds cannot repair a deterministic specification error. A seed can pass close
to a numerical gate; the Futtsu global minimum angle 30.008° is not evidence of
robust margin, and much of that limiting exterior is immutable.

Do **not** replace this with an unspecified “best of N.” First define hard
acceptance, including per-region achieved-size statistics and the intended
physical/bathymetric checks. Keep dt advisory as requested by the owner. Use a
fixed, recorded seed budget and stop at the first acceptable candidate. If a
user later wants optimization, compare *passing* candidates using a declared
objective, for example element count at required resolution, achieved dt and
mutable-zone quality margin. Optimizing dt alone rewards coarsening the
fishery; optimizing global minimum angle can be meaningless when the base binds
it. Ranking failed candidates by failed-gate count is diagnostic only, not a
notion of scientific closeness to acceptance.

Report exhaustion as “no acceptable candidate found within this budget,” not
proof that the target is impossible. Record backend/version and settings for
reproducibility as well as seed. Full best-of-N optimization is not worth doing
before the missing acceptance requirements.

## 4. The delivered artefact

Read `outputs/refine_futtsu_nori/sample_repro_final_futtsu_nori.14`, its map and
report, and the report's named base. Independent measurements:

| Check | Result |
| --- | --- |
| Counts | 4,734 → 5,826 nodes; 8,252 → 10,436 elements |
| Frozen copied coordinates/depths | All 4,619 mapped survivors exactly equal |
| Actual open-boundary list | One segment; exactly the mapped base list, including order |
| Boundary listing | All 1,230 topological boundary edges listed once; zero unlisted, non-boundary or repeated edges (closing island loops implicitly) |
| Triangle overlap | Zero positive-area intersecting pairs above 1e-6 m²; maximum intersection area 0 |
| Core edge lengths, midpoint within 300 m | 1,069 edges; min 21.91, median 30.226, p90 31.178, max 40.275 m |
| Re-evaluated new-node depths | All 1,207 match exactly; **six final nodes outside the base domain** use nearest-node fallback |
| Minimum element area | 232.45 m²; no tiny-element collapse |
| Total water area | 1,354,417,069.13 → 1,354,436,414.81 m²; +19,345.68 m², consistent with reported +0.001428% |
| Minimum-altitude dt | 11.8553 → 2.60802 s, reproducing the reported cost |
| Rerun default QA | **21/21 passed** |

The six outside-domain depths are a real qualification to “same interpolated
seabed,” already allowed by `depths_from_base`'s documented fallback. They are
not stale depths or a new demonstrated blocker. A final outside count is more
useful than adding the stitch and refresh counts, which need not be disjoint.

Two additional advisories in the rerun: **8 elements with estimated w/h < 1**
and **787 with w/h < 3**; **82 non-Delaunay internal edges (0.55%)**. Neither is
an established run-stopper. The channel metric depends on boundary polylines,
which are rebuilt here and known to be defective in the base, so raw before/
after counts would not isolate physical channel change. Inspect any relevant
channel against actual bank geometry before elevating this heuristic to a veto.
Default QA's informational dt is **3.01 s**, because it uses shortest edge
(`qa.py:952–960`); use the driver's **2.608 s altitude convention**, not the
larger advisory number.

**I found nothing else measured here that justifies stopping this artefact on
mesh structure alone.** Your local/global r-factor and production depth-floor/
cap measurements already prevent calling it the specified production-depth
mesh. Those measurements are accepted as supplied, not independently repeated
here. A controlled experiment on the same SRTM15 baseline is a different,
legitimate purpose if explicitly chosen.

The step ratio is 4.55 and the element-count ratio 1.265; their product is about
**5.75 times element-external-step work**. This is a workload estimate, not a
measured FVCOM runtime multiplier. It is a more useful budget warning than the
step ratio alone.

Limits: no FVCOM startup/integration, forcing/restart transfer, independent
survey validation, or full shoreline-overlap audit was run. The driver's QA
call does not enable the optional solid-land check (`qa.py:459–480,971`). Thus
21/21 is mesh QA, not validation of the entire model configuration or hydraulic
fidelity. Node/element renumbering also requires dependent model inputs to be
prepared for this mesh; the node map is not a remapped restart.

## 5. What to do next, in order

These are dispositions of known obligations, not additional findings.

1. **Decide the baseline and purpose.** For production, establish the intended
   M7001 depth field, floor/cap and slope treatment globally before freezing;
   then rerun the patch and measure it. For a refinement-only experiment,
   explicitly retain SRTM15 in both runs. Do not quietly mix the two. This costs
   users more scientifically than any meshing convenience on the list.
2. **Close acceptance structurally.** Carry contract evidence with the candidate
   as above; add per-region achieved-size reporting/acceptance with an explicit
   tolerance policy. Put unsupported OBC/priority refusals before meshing. These
   are small changes relative to another scientific run on the wrong mesh.
3. **Add the caller-declared maximum cut and segment-distance OBC guard.** This
   protects the comparison's frozen domain and limits unexpected work. It is a
   higher priority than improving a gradation diagnostic.
4. **Preserve/splice land-boundary semantics for general inputs.** Boundary
   types/grouping are model inputs. Do not blindly splice the known malformed
   base lists: validate or repair the baseline once under an explicit policy,
   then preserve valid untouched metadata. For the present ordinary 20/21
   Futtsu boundary, this is not an additional demonstrated blocker.
5. **Deliver geometry-from-file with CRS/filter validation and one actual
   polygon case.** High practical value: manually transcribing real fishery
   polygons is costly and error-prone. Before advertising channel/island
   support, add corresponding end-to-end cases with physical geometry checks;
   the old tolerance/provenance items matter more there than more seeds.
6. **Resolve priority semantics only to the extent needed.** Reject conflicting
   overlapping priorities now, or explicitly expose finest-wins as the policy.
   Implement full precedence when an actual user needs it. Two disjoint regions
   do not justify designing a complicated overlap hierarchy.
7. **Correct `effective_gradation`'s label cheaply; defer a certificate.** Remove
   the `within_c4` implication in output/docstrings. Actual C4 is already gated.
   A rigorous ambient-gradient bound is low value compared with delivered
   resolution, domain and depth checks; it would not certify C4 anyway.
8. **Repeated refinement: a sequential integration test and provenance first.**
   Check two disjoint patches and then overlapping ones against each immediate
   base, retaining parent hashes/maps. No need for a general ancestry database
   until original-baseline identity across many edits is a real requirement.
9. **Defer a separate global-rebuild product.** The sizing route already exists.
   Use it when the base itself needs rebuilding or patch feasibility fails and
   the owner accepts global changes. Building another interface now does not
   finish the exterior-preserving contract and complicates controlled comparisons.

A sensible stopping criterion is therefore not “works for every geometry.” It
is a declared supported envelope, mandatory acceptance evidence, early refusal
outside that envelope, one real noncircular use case, and an artefact validated
against the intended bathymetry/model inputs. Do that before extending the
repair search.

## Verification and runnable evidence

Repository suite:

```bash
/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python -m pytest -q
# 604 passed, 5 rasterio PendingDeprecationWarnings in 76.19 s
```

The focused contract test below was run as `/tmp/test_lr3_contract.py`:

```bash
MPLCONFIGDIR=/tmp/lr3-mpl \
/octfs/work/G16445/v61021/miniforge3/envs/oceanmesh-bench/bin/python \
-m pytest -q /tmp/test_lr3_contract.py --tb=short
# 1 intentionally failing API-contract assertion, 1 passing control in 2.25 s
```

No DistMesh generation, repair sweep, batch submission, installation, commit,
branch or sibling-repository edit. Only short read-only artefact probes and
unit tests ran on the login node. Lint was not rerun for this Markdown-only
change. One measurement invocation emitted a harmless Matplotlib cache-path
warning; subsequent probes used writable `MPLCONFIGDIR`.

### Contract test

```python
from pathlib import Path
import runpy
import numpy as np
import shapely
import fvcom_mesh_tools.patch as p
H = runpy.run_path(str(Path.cwd() / 'tests/test_patch.py'))

def test_unchecked_boundary_cannot_certify_contract():
    xy, tri = H['grid_mesh']()
    dep = np.full(len(xy), 8.0)
    sel, rc, pn, pt = H['replay_patch'](xy, tri, shapely.Point(400,400).buffer(250))
    nodes, faces, depths, nm, _ = p.stitch_patch(xy, tri, dep, sel, pn, pt, rc['pfix'], rc['pfix_base'])
    # Even an intact candidate must not be certified without the domain evidence.
    report = p.verify_patch(xy, dep, tri, sel, nodes, faces, depths, nm)
    assert report['boundary_checked'] is False
    assert not report['ok'], report

def test_complete_evidence_rejects_missing_face():
    xy, tri = H['grid_mesh']()
    dep = np.full(len(xy), 8.0)
    sel, rc, pn, pt = H['replay_patch'](xy, tri, shapely.Point(400,400).buffer(250))
    nodes, faces, depths, nm, _ = p.stitch_patch(xy, tri, dep, sel, pn, pt, rc['pfix'], rc['pfix_base'])
    boundary = p.boundary_after_patch(tri, sel, rc, nodes, nm)
    assert p.verify_patch(xy,dep,tri,sel,nodes,faces,depths,nm,expected_boundary=boundary)['ok']
    k = next(i for i in range(len(sel.retained),len(faces)) if not np.isin(faces[i], nm[sel.frozen_nodes]).any())
    result=p.verify_patch(xy,dep,tri,sel,nodes,np.delete(faces,k,axis=0),depths,nm,expected_boundary=boundary)
    assert not result['ok']
    assert result['n_unexpected_boundary_edges']==3
```

### Artefact measurement script

Run from the repository root with the environment's Python and
`MPLCONFIGDIR=/tmp/lr3-mpl OPENBLAS_NUM_THREADS=1`. This is a read-only probe,
not a production generation job.

```python
import json
from pathlib import Path
import numpy as np
import shapely
from pyproj import Transformer
from fvcom_mesh_tools.io.fort14 import read_fort14
from fvcom_mesh_tools.refine import depths_from_base
root = Path.cwd()
r = json.loads((root/'outputs/refine_futtsu_nori/report.json').read_text())
b = read_fort14(r['base_mesh'])
m = read_fort14(root/'outputs/refine_futtsu_nori/sample_repro_final_futtsu_nori.14')
nm = np.load(root/'outputs/refine_futtsu_nori/node_map.npy')
keep = nm >= 0
print('frozen', keep.sum(), np.array_equal(b.nodes[keep],m.nodes[nm[keep]]), np.array_equal(b.depths[keep],m.depths[nm[keep]]))
print('obc_exact',len(b.open_boundaries),len(m.open_boundaries),all(np.array_equal(nm[x],y) for x,y in zip(b.open_boundaries,m.open_boundaries)))
polys = shapely.polygons(m.nodes[m.elements])
tree = shapely.STRtree(polys)
i,j=tree.query(polys,predicate='intersects'); sel=i<j; i,j=i[sel],j[sel]
a=shapely.area(shapely.intersection(polys[i],polys[j]))
print('positive_overlap_pairs_gt_1e-6_m2',int((a>1e-6).sum()),'max_overlap_m2',a.max())
e,c=np.unique(np.sort(np.vstack([m.elements[:,[0,1]],m.elements[:,[1,2]],m.elements[:,[2,0]]]),axis=1),axis=0,return_counts=True)
bound={tuple(x) for x in e[c==1]}
listed=[]
for x in m.open_boundaries:
    listed.extend(tuple(sorted(z)) for z in zip(x[:-1],x[1:]))
for typ,x in m.land_boundaries:
    listed.extend(tuple(sorted(z)) for z in zip(x[:-1],x[1:]))
    if typ==21 and x[-1]!=x[0]: listed.append(tuple(sorted((x[-1],x[0]))))
print('boundary_edges',len(bound),'listed_unique',len(set(listed)),'unlisted',len(bound-set(listed)),'nonboundary_listed',len(set(listed)-bound),'repeated',len(listed)-len(set(listed)))
cx,cy=Transformer.from_crs(4326,32654,always_xy=True).transform(139.7881,35.3228)
mid=m.nodes[e].mean(axis=1); core=np.linalg.norm(mid-[cx,cy],axis=1)<300
length=np.linalg.norm(m.nodes[e[:,0]]-m.nodes[e[:,1]],axis=1)
print('core_midpoint_edges',core.sum(),'length_min_median_p90_max',np.quantile(length[core],[0,.5,.9,1]).tolist())
new=np.setdiff1d(np.arange(m.n_nodes),nm[keep]); dep,out=depths_from_base(b.nodes,b.elements,b.depths,m.nodes[new])
print('new_depth_final_outside',out,'max_abs_depth_error',np.max(abs(dep-m.depths[new])))
for label,mesh in [('base',b),('patch',m)]:
    pts=mesh.nodes[mesh.elements]; area=shapely.area(shapely.polygons(pts)); side=np.stack([np.linalg.norm(pts[:,(k+1)%3]-pts[:,k],axis=1) for k in range(3)],axis=1); dt=2*area/side.max(axis=1)/np.sqrt(9.81*mesh.depths[mesh.elements].max(axis=1))
    print(label,'nodes',mesh.n_nodes,'faces',mesh.n_elements,'dt',dt.min(),'min_area',area.min(),'area',area.sum(),'depth_range',[mesh.depths.min(),mesh.depths.max()])

from fvcom_mesh_tools.qa import run_qa
qa = run_qa(m)
print('QA', qa.n_gate_total, qa.n_gate_failed)
for check in qa.checks:
    if not check.gate:
        print(check.check_id, check.observed)
```
