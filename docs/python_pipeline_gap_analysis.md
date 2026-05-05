# Python pipeline gap analysis (Phase 4)

This note records what is missing from the Python mesh-generation pipeline
relative to the legacy OceanMesh2D (MATLAB) workflow used to produce
``data/mesh/reference/tokyo_bay/tb_futtsu20220311.14``. Numbers come from
the PoC #5 / PoC #6 runs; see ``outputs/05_*`` and ``outputs/06_*``.

The reference workflow is OceanMesh2D + manual post-edits in MATLAB.
The Python pipeline currently exercised is:

```
DEM (NetCDF) -> ocsmesh.Geom(zmax=0)
              -> ocsmesh.Hfun(hmin, hmax)
              -> ocsmesh.MeshDriver(engine="gmsh").run()
              -> ocsmesh.EuclideanMesh2D.write(format="grd")
```

## 1. Headline parity numbers

Reference is the legacy fort.14; "OCSMesh minimal" is PoC #5 with
``hmin=200 m``, ``hmax=5000 m``, ``zmax=0``, gmsh engine.

| Metric                                | Reference | OCSMesh minimal |
| ------------------------------------- | --------- | --------------- |
| NP                                    | 95,551    | 15,757 (×0.165) |
| NE                                    | 182,603   | 19,711 (×0.108) |
| Mesh-generation wall time             | (legacy)  | 26.75 s         |
| Edge length, p50 (m)                  | 277       | 111             |
| Edge length, p95 (m)                  | 23,169    | 562             |
| Edge length, max (m)                  | 107,889   | 8,717           |
| Triangle alpha-quality, mean          | 0.979     | 0.725           |
| Fraction alpha < 0.3                  | 0.00 %    | 5.87 %          |
| Triangle min interior angle, p50 (deg)| 54.3      | 33.9            |
| Fraction min angle < 20 deg           | 0.00 %    | 19.59 %         |
| Open-boundary segments / nodes        | 1 / 193   | 0 / 0           |
| Land-boundary segments / nodes        | 54 / 8,464| 0 / 0           |
| Land-boundary ibtypes                 | {20, 21}  | (none)          |

Domain extents differ on purpose: the reference covers the Pacific halo
out to roughly (117 E, 18 N) — (165 E, 62 N), while the minimal pipeline
meshes only the DEM bounding box (Tokyo Bay proper). Edge-length stats
are therefore *not* directly comparable in scale, but the spread (p50 vs
p95 vs max) is — and the reference uses graded sizing across four orders
of magnitude while the minimal Hfun is uniform.

## 2. Capability gaps

Each gap is rated S / M / L for implementation effort within Python.

### 2.1 Boundary classification (S–M)

Status: **missing entirely.** OCSMesh's ``meshdata_to_grd`` writes a
fort.14 with ``0 ! total number of open boundaries`` and ``0 ! total
number of land boundaries`` whenever ``mesh.boundaries`` is empty, which
is the default after ``MeshDriver.run()``.

What we actually need:

- An *open* boundary at the seaward arc.
- *Land* boundaries for every coast segment (ibtype 0 by ADCIRC
  convention; the reference uses 20 / 21 — non-zero values are an FVCOM
  / project-specific convention worth documenting separately).
- A way to tell the two apart automatically (e.g. seaward = farthest
  from coastline shapefile, or seaward-side-of-DEM mask).

OCSMesh ships ``Mesh.boundaries`` and a ``Boundaries`` helper, but they
require either a manual partition or a polygonal seam. The boundary
walker we already have in ``oceanmesh-tools`` covers most of this and
can be ported.

### 2.2 Depth interpolation (S)

Status: **silently zero.** With our minimal driver, the GRD writer falls
through to ``np.zeros(len(coords))`` because ``mesh.values`` is None.

Fix: call ``mesh.interpolate(raster_collection, method=...)`` after
``driver.run()`` and before ``mesh.write(...)``. Single line; only listed
here so it is not forgotten when the pipeline matures.

### 2.3 Mesh quality (M, partially closed)

The OceanMesh2D reference has 0.00 % triangles below the 20-degree
min-angle threshold and a mean alpha of 0.979. The raw gmsh output has
20.75 % below 20 degrees and mean alpha 0.718.

PoCs #8 / #9 / #10 closed part of this gap; the rest is sizing-driven.

| Pass on PoC #7 mesh             | mean alpha | frac<20° |
| ------------------------------- | ---------- | -------- |
| raw gmsh + buildmesh            | 0.7180     | 20.75 %  |
| Laplacian only (PoC #8)         | 0.7169     | 21.44 %  |
| Edge-swap only (PoC #9)         | 0.7367     | 17.74 %  |
| Swap + smooth combo (PoC #10)   | **0.7498** | **17.12 %** |
| Reference                       | 0.979      | 0.00 %   |

Findings:

- Pure Laplacian smoothing barely moves the metric on this mesh -
  slivers are *topologically* trapped; node moves alone cannot fix
  three-near-collinear vertices.
- Edge swap (Lawson / min-angle flip) is monotonically helpful and
  cheap (0.45 s on 19 k triangles).
- The swap+smooth combination plateaus at ~17 % bad triangles after
  3-4 rounds. The plateau is set by the *initial size function*: with
  a coarse uniform Hfun(hmin=200, hmax=5000) the bay edge nodes are
  forced into thin triangles that no in-place rearrangement can fix.

Closing the residual gap requires either:

1. Adaptive sizing during generation (coastline-distance Hfun, see 2.4
   below).
2. Local refinement of bad-quality regions (split + retriangulate any
   triangle below threshold). Not implemented.
3. Constrained Delaunay with explicit feature edges. Not implemented.

The combo loop is exposed as ``fmesh-buildmesh --quality-pass N``;
default 0 (off) since the user usually wants the option to inspect
the raw output first.

### 2.4 Coastline-aware sizing (DONE, PoC #12)

The reference uses graded sizing driven by both DEM gradient and
distance-to-shoreline (typical OceanMesh2D recipe). The minimal
``Hfun(raster, hmin=200, hmax=5000)`` was the dominant residual error
source: even after a 6-round swap+smooth quality pass, ~17 % of
triangles stayed below the 20-degree min-angle threshold.

PoC #12 plumbs ``Hfun.add_feature`` with the MLIT C23 Tokyo Bay
coastline shapefile (1,325 LineStrings clipped to the DEM bbox -> 571
features). Results on the same DEM / hmin / hmax:

| Config            | NP     | NE     | alpha mean | frac<20° | wall |
| ----------------- | ------ | ------ | ---------- | -------- | ---- |
| baseline (no coast) | 15,757 | 19,711 | 0.718      | 20.75 %  | 40 s |
| coast only          | 18,953 | 26,103 | 0.825      | **3.03 %** | 32 s |
| coast + 6-round qp  | 18,953 | 26,103 | **0.848**  | **2.70 %** | 33 s |
| reference           | 95,551 | 182,603| 0.979      | 0.00 %   | -    |

So coastline-aware sizing alone cuts ``frac<20deg`` by **86 %**
(20.75 -> 3.03 %); the combined pipeline gets to **2.70 %**. Mesh size
grows by ~20 % (15.7k -> 18.9k nodes) - a cheap cost for the gain.

CLI surface:
``fmesh-buildmesh DEM out.14 --coastline coast.shp --coast-target-size METRES --coast-expansion-rate RATE``.

Multiple ``--coastline`` flags can be combined; inputs are reprojected
to EPSG:4326 and clipped to the DEM bbox before being fed to OCSMesh.

### 2.5 Open-boundary perpendicularity (DONE for fixes, partial for generation)

PoC #4 / ``fmesh-perpfix`` already corrects edge perpendicularity at the
open boundary in a finished mesh. Generating a mesh that is born
perpendicular requires geometry-side support: either inserting a
"perpendicular spine" at the open arc during meshing, or running
``fmesh-perpfix`` as a post-step. Recommended near-term: do the latter.
That avoids modifying the mesher.

### 2.6 Channel widening / river inflow nodes (DONE, PoC #16)

Reference includes hand-edited ibtype=21 segments around river mouths
(Edo, Ara, Sumida, Tama, Tsurumi). ``fmesh-buildmesh`` now reads a
CSV/GeoJSON/shapefile of river-mouth points via
``--river-inflow-points`` and snaps each point to the nearest
land-boundary node (skipping nodes already inside any ibtype=21
segment), then splits the parent land segment into prefix / river /
suffix. Tunables: ``--river-segment-nodes`` (nodes per inflow
segment, default 5), ``--river-ibtype`` (default 21),
``--river-snap-tol-m`` (default unbounded). PoC #16 on Tokyo Bay
snaps 5 rivers at 310-702 m and produces 5 ibtype=21 segments
totalling 25 nodes, with no quality regression
(``frac<20deg`` 1.13 %, ``alpha`` mean 0.847).

### 2.7 Mesh-generation determinism / reproducibility (S)

gmsh produces a different mesh every time unless seeded. OCSMesh does
not currently expose the seed through ``MeshDriver``. Worth filing a
follow-up to set ``Mesh.Algorithm`` / ``Mesh.RandomSeed`` explicitly
once we move beyond PoCs.

### 2.8 Cross-basin portability (DONE, PoC #17)

Validated on Osaka Bay using the SRTM15+ global DEM (subset to bbox)
and the GSHHS-f L1 global coastline. The same flag set as PoC #16
(no Tokyo-Bay-specific tuning) produced a clean fort.14:

* NP=3,412, NE=5,451, ``flipped=0``,
* alpha mean 0.898, ``frac<20deg`` 0.18 %,
* 3 independent open arcs detected (Akashi Strait + Tomogashima
  + a small Awaji-south passage) - the bbox classifier handles
  multi-strait geometries without merging unrelated arcs,
* 4 ibtype=21 river segments (Yodo, Yamato, Mukou, Kanzaki).

Two operational notes from PoC #17:

* SRTM15+ ships without an embedded CRS, so ocsmesh.Raster cannot
  open it directly. The PoC subsets the bbox and writes a CF-tagged
  EPSG:4326 GeoTIFF before passing it to ``fmesh-buildmesh``. A
  reusable ``fmesh-subset-dem`` helper would be a natural follow-up.
* GSHHS-f L1 has far fewer line strings than MLIT C23 (8 vs. 571 in
  the bbox), so ``Hfun.add_feature`` runs in ~1 s and the resulting
  mesh is much sparser. Quality is *better* (fewer enforced kinks
  to honour); fidelity is lower (no fine-scale Tokyo-Bay-style
  estuary detail).

## 3. License & dependency notes

- ``ocsmesh`` itself is CC0-1.0; importing it from Apache-2.0 code is
  fine.
- ``ocsmesh`` pulls ``gmsh`` (GPL-2.0+) at runtime when the gmsh engine
  is selected. We use it as an external tool via ``ocsmesh``, not by
  linking to ``libgmsh`` directly, but downstream redistributors should
  be aware of the GPL footprint.
- ``triangle`` (the alternative engine) is released under a custom
  research-only license and is **not** OSI-approved; avoid making it a
  default.

## 4. Recommended next iteration

In order of return-on-effort:

1. **Boundary classification + depth interpolation** (DONE, PoCs #5-#7).
   ``fmesh-buildmesh`` produces a parseable, FVCOM-shaped fort.14 in
   one shot; ``fmesh-perpfix`` runs automatically as a post-step.
2. **Quality post-processing loop** (PARTIAL, PoCs #8-#10). Edge-swap
   plus damped Laplacian smoothing is wired in as
   ``fmesh-buildmesh --quality-pass N``. On the baseline mesh it cuts
   ``frac<20deg`` from 20.75 % to 17.12 % (plateau set by sizing).
3. **Coastline-aware sizing** (DONE, PoC #12). Driving Hfun with the
   MLIT C23 coastline shapefile slashes ``frac<20deg`` to 3.03 %
   alone, or 2.70 % with the quality pass on top. Available as
   ``fmesh-buildmesh --coastline path.shp ...``.
4. **Local refinement of bad-quality regions** (DONE, PoC #15).
   Longest-edge bisection (Rivara-style) for triangles below the
   user-supplied min-angle threshold. Available as
   ``fmesh-buildmesh --refine-min-angle DEG``. Cuts ``frac<20deg``
   from 2.84 % to 1.82 % at +6.7 % mesh-size cost, no quality
   regression. Centroid insertion was tried first but empirically
   *increased* sliver count (it splits a sliver into one wedge that
   keeps the long base). Iteration includes regression-rollback so
   the worst-case is "no change", not "worse than before".
5. **Open-segment merging** (DONE, PoC #14).
   ``--open-merge-coast-gap NODES`` bridges short coast intrusions
   between two open arcs. With ``gap=50`` the Tokyo Bay output
   collapses from 3 open segments to 1 contiguous arc of 778 nodes.
6. **Island merging / sliver-island filtering** (DONE, PoCs #13-#14).
   ``--min-polygon-area-m2`` drops detached single-pixel water
   bodies; ``--min-island-area-m2`` fills holes (islands) below a
   metric area threshold. With ``min_island=1e5`` the Tokyo Bay
   output goes from 165 to 23 land segments with no quality
   regression and a slight wall-time reduction (gmsh has fewer
   features to honour).
7. **River-inflow segments** (DONE, PoC #16).
   ``fmesh-buildmesh --river-inflow-points`` reads a points file
   (CSV/GeoJSON/shapefile) and produces ibtype=21 land segments at
   each river mouth. Nodes already inside an existing ibtype=21
   segment are skipped, so re-running the pipeline is idempotent.
   On Tokyo Bay this matches the reference structure (5 rivers, one
   open arc, 26 ibtype=20 + 5 ibtype=21 land segments).
8. **Cross-basin validation** (DONE, PoC #17). Same flag set as
   PoC #16 ports cleanly to Osaka Bay using SRTM15+ + GSHHS-f L1
   inputs. Three open arcs (Akashi + Tomogashima + Awaji-south)
   are detected without merging; quality is good (alpha 0.90,
   ``frac<20deg`` 0.18 %).
9. Item 2.7 (reproducibility) and a reusable DEM subset helper
   (currently inlined in the PoC) are the remaining minor levers.
