# Local refinement of an existing mesh — design

**Status: partially implemented.** The specification and its pre-flight checks
are in `src/fvcom_mesh_tools/refine.py` with tests; the patch generator is
designed here and not yet written. This document is revised as the design
moves, and records the decisions and their reasons so that a later change is
made knowingly rather than by accident.

Revision history is the git history of this file.

## 1. What problem this solves

A user wants one part of the domain at a much finer resolution than the rest —
a fishery, a port basin, a discharge site — and wants **the rest of the mesh
left alone**.

The existing `regions:` block in a sizing recipe (`src/fvcom_mesh_tools/
sizing.py`, `docs/sizing_regions.md`) does not do this. It edits the **sizing
field**, after which the whole mesh is regenerated. DistMesh is a global
relaxation: perturbing the field anywhere moves nodes everywhere. The result
is a different mesh, not the old mesh with a refined patch.

So there are two operations, sharing one vocabulary for geometry and targets:

| | outside the region | cost | status |
|---|---|---|---|
| **A. sizing region** | changes | one full rebuild | implemented (`sizing.py`) |
| **B. local refinement** | **unchanged, and checked** | one patch | this document |

B is what the rest of this document specifies.

## 2. The three zones and the contract

| zone | definition | guarantee |
|---|---|---|
| **core** | inside the declared geometry | achieved edge length ≈ `target_h_m` |
| **transition** | the annulus the gradation needs to reach the ambient size | remeshed; changes |
| **frozen** | everything else | node coordinates and connectivity **identical** to the base mesh |

The frozen guarantee is the point of the exercise, so it is **checked, not
asserted**: `refine.frozen_changes()` counts the nodes outside core+transition
whose coordinates moved, and the run fails if that count is not zero.

## 3. Two facts the caller cannot see, and must not have to guess

### 3.1 The transition width is not a free choice

With a linear size gradation `g`, the edge length grows as
`h(r) = h_target + g · (r − r_core)`, so reaching the ambient size takes

```
W = (h_ambient − h_target) / g
```

For the certified Tokyo Bay chain (`g = 0.165`, ambient ≈ 350 m):

| target | elements across | transition W | W / core radius (300 m) |
|---|---|---|---|
| 15 m | 20.2 | 1,909 m | 6.4 |
| **30 m** | **15.7** | **1,939 m** | **6.5** |
| 50 m | 12.4 | 1,697 m | 5.7 |
| 100 m | 7.8 | 1,394 m | 4.6 |

**The transition is several times the size of the region itself.** The user
declares the target; the width follows. An explicit `transition_m` is checked
against `W` and rejected when the gradation cannot cover it, with the three
ways out named in the message (widen it, coarsen the target, steepen the
gradation).

Steepening the gradation shortens the transition but spends the C4 margin
(adjacent element area change ≤ 0.5, i.e. a linear ratio ≤ √2 ≈ 1.41 per
edge). `g = 0.33` halves the width; beyond that, expect QA failures.

### 3.2 The time step is the cost that bites

FVCOM integrates with **one global external step**, so the smallest element
anywhere sets it for the whole run:

```
dt = h / sqrt(g · H)
```

| depth | h = 30 m | h = 50 m | h = 100 m | h = 350 m |
|---|---|---|---|---|
| 2 m | 6.77 s | 11.29 s | 22.58 s | 79.0 s |
| 3 m | 5.53 s | 9.22 s | 18.43 s | 64.5 s |
| 5 m | 4.28 s | 7.14 s | 14.28 s | 50.0 s |
| 10 m | 3.03 s | 5.05 s | 10.10 s | 35.3 s |

The certified mesh allows 11.9 s. A 30 m region therefore costs a factor of
2–4 in run time, **depending entirely on the depth of that region**: shallow
water is cheap to refine, deep water is not.

`dt_floor_s` is **mandatory** in the recipe, and `refine.preflight()` refuses a
target the water cannot carry, naming the largest target that water does
permit. A failed pre-flight costs a second; a failure discovered after meshing
costs the run.

### 3.3 Node count, for completeness

Core elements ≈ `area / (√3/4 · h²)`. A 300 m disc at 30 m is ~725 elements;
its transition adds ~830. Against a base mesh of 8,252 elements that is
**+16 %** for one small fishery, and the transition is the larger half.
Enlarging the core is comparatively cheap — the transition width does not grow
with it.

## 4. Specification

A refinement recipe is YAML, version-controlled, and lives in
`recipes/refine/`. Environment variables (`SR_*`) remain for sweeps; a durable
specification belongs in a file.

```yaml
base_mesh: ../../outputs/verify_409.115302/fit/sample_repro_final.14
dt_floor_s: 4.5
gradation: 0.165

refine:
  - name: futtsu_nori
    geometry:
      circle: {center: [139.7881, 35.3228], radius_m: 300}
    target_h_m: 30        # achieved edge length, as elsewhere in the chain
    touch_coast: false
    priority: 0
```

### Keys

| key | meaning |
|---|---|
| `base_mesh` | the mesh to patch; relative paths resolve against the recipe |
| `dt_floor_s` | mandatory; the external time step the result must still allow |
| `gradation` | size growth per metre through the transition |
| `refine[].name` | unique, non-empty |
| `refine[].geometry` | `circle`, `bbox` or GeoJSON `Polygon`, all lon/lat |
| `refine[].target_h_m` | **achieved** edge length, matching `SR_H_TARGET=achieved` |
| `refine[].transition_m` | optional; derived when absent, checked when present |
| `refine[].touch_coast` | may the core overlap land? default `false` |
| `refine[].priority` | higher wins an overlap; ties take the smaller target |

Unknown keys, non-finite values and invalid geometry are errors.

### Geometry forms

`circle: {center: [lon, lat], radius_m: N}` is the most direct way to say
"this fishery, roughly here" and is the recommended form for exploration.
`bbox` and GeoJSON `Polygon` cover the rest. For production work a real
boundary — a fishery-right polygon read from GeoJSON or a shapefile — is more
reproducible than typed coordinates; **reading geometry from a file is not yet
implemented** and is the first extension to make.

### `touch_coast`

When false, a core overlapping land is refused. When true, the coastline
inside the patch is re-cut at the target size. That improves shoreline
fidelity where it matters most, but it changes the land-boundary node list, so
the frozen guarantee then covers only the water boundary away from the patch.
**Default false**, because it keeps the first runs simple to verify.

## 5. The patch generator (designed, not yet written)

1. **Select.** Build `core ⊕ transition` in the mesh CRS (EPSG:32654). Mark
   every element whose centroid falls inside it.
2. **Cut.** Delete those elements. What remains is the base mesh with a hole;
   the hole's rim is a closed chain of **existing** nodes.
3. **Fill.** Generate a triangulation inside the hole with
   `h(r) = target + g · (r − r_core)`, clipped above by the ambient field, and
   with the rim nodes as fixed points (`pfix`). The rim spacing already equals
   the local base-mesh edge length, so the ambient end of the gradation
   matches by construction.
4. **Stitch.** Concatenate, deduplicate the rim nodes, renumber.
5. **Verify.** `fmesh-mesh-qa` (21 gates), `frozen_changes()` = 0, the 342
   connectivity comparator, the 364 width gate, and the implied dt against
   `dt_floor_s`.

### Constraints the generator must respect

- The hole must not reach the **open boundary**: OBC nodes are an input.
  Refuse if `core ⊕ transition` intersects the OBC arc.
- With `touch_coast: false` the hole must not reach the coastline either, so
  the land boundary is untouched.
- Every rim node must survive as a distinct vertex. Two rim nodes closer than
  the local target size collapse onto one and DistMesh fails; the oceanmesh
  fork now reports which points collapsed (`generate_mesh`, commit `529a462`),
  so this surfaces as a named error rather than a shape mismatch.

### Where the code goes

The package is Apache-2.0 and **must not import oceanmesh (GPL)**. The split:

- `src/fvcom_mesh_tools/refine.py` — specification, pre-flight, selection,
  stitching, verification. License-clean.
- `notebooks/4xx_local_refine.py` — the DistMesh call, as notebook 325 already
  does. Notebooks are not part of the distributed package.

## 6. Worked case: Futtsu nori area

`recipes/refine/futtsu_nori.yaml`. The centre came from the data, not by eye:
the Futtsu tidal flat is the largest connected patch of water shallower than
1.5 m (T.P.) off Futtsu — 4.02 km², lon 139.768–139.821, lat 35.304–35.331 —
and the declared centre lies 654 m north of its edge.

| quantity | value |
|---|---|
| centre / radius | (139.7881, 35.3228) / 300 m |
| target | 30 m |
| core depth | 3.00 – 4.15 m |
| **dt** | **4.70 s** (floor 4.5 s) |
| transition | 1,939 m (required: 1,939 m) |
| elements core / transition / replaced | 725 / 830 / 250 |
| **elements added** | **+1,305** (base 8,252) |

Of 1,200 candidate cells 650–750 m north of the flat with the core clear of
land, 112 also met the 4.5 s floor. The chosen point is the one nearest the
middle of that band in longitude; the shallowest candidates reach 4.88 s but
sit on the flat itself.

## 7. Open questions

1. **Geometry from a file.** Needed before this is useful for real fishery
   boundaries. Shapefile or GeoJSON, with an optional attribute filter.
2. **Several regions at once.** `priority` is specified and validated but the
   generator's overlap handling is not designed. Overlapping cores with
   different targets need a single combined sizing field before cutting.
3. **Bathymetry inside the patch.** The base mesh's depths come from the
   production M7001 recipe (`dem/m7001.py`). New nodes must be given depths
   the same way, then the r-factor smoothing re-run — but re-running it
   globally would move frozen nodes' depths. Proposal: interpolate M7001 at
   the new nodes, apply the 3 m floor and 300 m cap, and run the r-factor
   limiter **only over the patch and its rim**, holding frozen depths fixed.
   Not yet decided.
4. **Repeated refinement.** Applying a second region to an already-patched
   mesh should work, but the frozen check then compares against the patched
   mesh, not the original. Whether to track a chain of base meshes is open.
5. **Does the fit pass re-run?** `coast_fit` moves boundary nodes; over a
   patch that touches the coast it should run again, but only on the patch.
