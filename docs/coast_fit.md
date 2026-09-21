# Fitting the mesh shoreline to the coastline

## The problem

The mesh boundary produced by DistMesh is not on the coastline, and the gap is
large enough to see by eye.

DistMesh moves points along bar forces and then projects a point back to the
domain boundary **only when that point has stepped outside it**.  A point that
stops just inside is left where it is.  Which triangles survive is decided
afterwards, by deleting the ones whose centroid falls on land.  Nothing in that
sequence asks a boundary node to sit *on* the coast, so the mesh shoreline
wanders across the real one in both directions.

Measured on the certified Tokyo Bay mesh (4,734 nodes, 1,219 land boundary
nodes), distance from each boundary node to the nearest point of the coastline:

| reference | median | p90 | p99 | max | nodes on the land side |
|---|---|---|---|---|---|
| the polygon handed to oceanmesh | 10.1 m | 84.0 m | 198.4 m | 565.6 m | 388 (32 %) |
| raw OSM land | 26.4 m | 150.6 m | 275.1 m | 565.6 m | 568 (47 %) |

117 of 1,219 nodes sit further than a quarter of their local edge length from
the coastline they were meant to follow (313 against raw OSM).

Two things it is *not*:

* **Not shoreline smoothing.**  `Shoreline()` densifies at `0.5 * h0` and then
  applies a 5-point boxcar.  Measured separately, that filter displaces the
  shoreline by a median of 1.8 m and a p90 of 14.7 m -- it shows in the p99
  tail, not in the bulk.
* **Not our preprocessing.**  The gap against the polygon oceanmesh was
  actually given is already 10 / 84 / 566 m.  The extra gap against raw OSM is
  the recipe edits, the normalize pass and the waterway carving, which are
  deliberate.

## The pass

`fvcom_mesh_tools.coast_fit.fit_boundary_to_coast` relaxes the land boundary
nodes onto the reference polygon:

1. Each sweep moves every movable boundary node `relax` (0.8) of the way to
   its nearest point on the coastline, so the mesh eases towards the coast
   instead of snapping to it.
2. A node never ends further than `max_move_frac` (0.5) local edge lengths
   from where it started.
3. A move is applied only if it survives every guard; otherwise that one node
   is put back and the sweep continues.

The guards, in the order they are tested:

| guard | rejects a move that would |
|---|---|
| orientation | flip an incident element or collapse it to zero area |
| C1 / C2 | take an incident element below 30° or above 130° |
| C4 | push the area change across a shared edge past 0.5 |
| time step | take an incident element below `dt_floor_s` |
| wet centroids | put an incident element's centroid on land |

Each bound is waived for an element that already violated it, provided the
move improves that element -- so the pass can help a bad element but never
creates a new violation.

`dt_floor_s` defaults to the mesh's own current minimum, i.e. **the fit is not
allowed to cost any time step at all**.  This matters: without the guard the
same fit costs 5 % of dt (16.34 → 15.44 s) because pulling a boundary node onto
the coast can flatten one triangle, and dt is set by the single worst one.
Pass `dt_floor_s=0` to switch the guard off.

No node is inserted or deleted and no edge is flipped.  Connectivity, the
one-element-wide ledger and the land-breach census are therefore unchanged by
construction, and only the coordinates differ.

### Which coastline to fit to

The reference must be **the polygon the mesh generator was given**
(`land_channel_adj.shp`), not the raw source data.  The difference between the
two is preprocessing we chose on purpose -- channel widening, closed dead ends,
filled data cracks -- and fitting to raw OSM would quietly undo it.

## Measured result

On the passing `achieved` mesh, 6.5 s of wall time:

| | before | after |
|---|---|---|
| \|offset\| median | 10.1 m | **0.0 m** |
| \|offset\| p90 | 84.0 m | **11.8 m** |
| boundary nodes moved | -- | 1,085 / 1,219 |
| QA gates | 21/21 PASS | 21/21 PASS |
| min / max interior angle | 30.01° / 116.00° | 30.01° / 119.26° |
| implied dt | 16.34 s | 16.38 s |
| nodes / elements | 4,734 / 8,252 | unchanged |

## Limits

The worst node stays 566 m off.  Its offset is larger than the move budget, so
no amount of relaxation reaches it: closing that kind of gap needs a node
inserted on the boundary edge, which changes the node count and is a separate
decision.

Fidelity finer than the local edge length is out of reach by construction.  A
cape narrower than `h` cannot be represented by a boundary that has no node
inside it; that is a sizing question, not a fitting one.

## Usage

```python
from fvcom_mesh_tools.coast_fit import fit_boundary_to_coast

res = fit_boundary_to_coast(mesh.nodes, mesh.elements, land_utm,
                            fixed=obc_node_ids, depths=mesh.depths)
mesh.nodes[:, :2] = res.nodes[:, :2]
print(res.summary())
```

In the sample chain the pass runs at the end of `notebooks/331_finish2.py` and
is on by default; `SR_COAST_FIT=off` skips it.  `notebooks/404_coast_fit.py`
measures and maps the remaining gap, and `notebooks/406_coast_fit_compare.py`
draws the worst sites before and after.
