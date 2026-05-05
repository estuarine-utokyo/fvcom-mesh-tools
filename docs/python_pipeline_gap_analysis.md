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

### 2.3 Mesh quality (M)

The OceanMesh2D reference has 0.00 % triangles below the 20-degree
min-angle threshold and a mean alpha of 0.979. The raw gmsh output has
19.59 % below 20 degrees and mean alpha 0.725. OceanMesh2D applies
several iterations of:

- bound-preserving Laplacian smoothing,
- edge swapping ("flip"), and
- removing slivers / collapsed triangles.

OCSMesh exposes some quality post-processing under ``ocsmesh.utils`` and
``MeshDriver`` has ``cleanup_isolates`` / ``cleanup_duplicates`` but no
full equivalent. Closing this gap means either (a) calling
``meshkernel`` orthogonalize + smoothing (PoC #3 showed this *worsens*
the open-boundary perpendicularity but improves global quality), (b)
porting the OceanMesh2D quality loop, or (c) integrating ``triangle``
or another smoother. **Largest single block of work** on the roadmap.

### 2.4 Coastline-aware sizing (M)

The reference uses graded sizing driven by both DEM gradient and
distance-to-shoreline (typical OceanMesh2D recipe). Our minimal Hfun is
just ``Hfun(raster, hmin=200, hmax=5000)`` — no shoreline shapefile
input, no ``add_feature_size`` calls.

OCSMesh has ``Hfun.add_feature_size`` and ``Hfun.add_courant_size``;
plumbing the existing
``data/coastline/tokyo_bay/coastline.shp`` (or whatever resolved by
``oceanmesh-tools scan``) into Hfun is straightforward but needs a small
glue layer because OCSMesh's API is geometry-first while OceanMesh2D's
is raster-first.

### 2.5 Open-boundary perpendicularity (DONE for fixes, partial for generation)

PoC #4 / ``fmesh-perpfix`` already corrects edge perpendicularity at the
open boundary in a finished mesh. Generating a mesh that is born
perpendicular requires geometry-side support: either inserting a
"perpendicular spine" at the open arc during meshing, or running
``fmesh-perpfix`` as a post-step. Recommended near-term: do the latter.
That avoids modifying the mesher.

### 2.6 Channel widening / river inflow nodes (M–L)

Reference includes hand-edited refinement around river mouths (Sumida,
Tama, Edo). The Python pipeline has nothing equivalent. Likely
implementation: shapefile-driven local sizing + ``omesh14-edit-bdy``-
style post-edits to mark ibtype=21 inflow nodes.

### 2.7 Mesh-generation determinism / reproducibility (S)

gmsh produces a different mesh every time unless seeded. OCSMesh does
not currently expose the seed through ``MeshDriver``. Worth filing a
follow-up to set ``Mesh.Algorithm`` / ``Mesh.RandomSeed`` explicitly
once we move beyond PoCs.

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

1. **Boundary classification + depth interpolation** (small wins,
   prerequisite for everything else). Output a fort.14 the FVCOM
   harness can actually run.
2. **Quality post-processing loop** (port OceanMesh2D's
   smooth + swap + sliver-removal, or chain MeshKernel). Without this,
   the Python output will not pass FVCOM stability checks.
3. **Coastline-aware sizing** (close the resolution-distribution gap).
4. **Run ``fmesh-perpfix`` automatically** as a post-step in any wrapper
   we publish.
5. Items 2.6 (river channels) and 2.7 (reproducibility) deferred until
   1–4 are stable.
