# Regional resolution inputs

Set `SR_SIZING=recipes/sizing/tokyo_bay.yaml` for notebook 325. YAML supports
comments (including the disabled 30 m example) and uses the existing pipeline
PyYAML dependency. Geometry uses the existing Shapely vector extra. No new
dependencies are added. Recipe values override SR_H0, SR_MAXEL, SR_GRADE, SR_DT
and SR_CRMIN. Other legacy controls remain active, including SR_CFL_DT and
SR_CFL_CRMAX for the graded depth-dilated guard. In the present certified chain
DT/CRMIN are legacy knobs; the recipe CFL pair diagnoses regions, and does not
replace that separately calibrated guard.

Required keys: `coastal_target_m`, `max_edge_length_m`, positive `gradation`,
`cfl: {dt_s: positive, cr: positive}`, and `regions` (possibly empty).
Optional `hmin_m` is an absolute scalar lower bound, independent of the coastal
target: leaving it unset permits a 30 m region inside a 290 m coastal field.
Each region requires a unique nonempty `name`, `target_h_m`, and `geometry`.
Geometry is `{bbox: [west, south, east, north]}` or a GeoJSON Polygon with closed
rings and optional holes, in lon/lat degrees. Optional `transition_m` is a
nonnegative linear blend distance outside the polygon; `priority` is numeric.
Unknown keys, invalid polygons, nonfinite values and negative sizes fail.

Higher priority wins overlap; equal priority chooses the smaller requested
target (including transition footprints), independent of list order. Each
blend uses the original ambient field. Targets replace sizes, so they may
coarsen as well as refine. Gradation can subsequently lower any target.
Regions run after legacy guards and corridors; boundary constraints themselves
are unchanged, so refinement near pinned boundaries needs a mesh QA run.
With no SR_SIZING the new code never reads or transforms the sizing field.
The example's empty region list also skips the new limiter.

`apply_sizing_regions` accepts finite positive 2-D metre values and matching
2-D longitude/latitude arrays. It copies inputs. The caller converts any
OceanMesh degree-sized values. Local distances use 111000 m/degree and cosine
of mean latitude; use local domains away from poles/dateline. The limiter is
an eight-neighbour multi-source Dijkstra lower envelope, enforcing the slope
on graph edges, not a continuous Euclidean gradient bound. This avoids a GPL
OceanMesh import in the Apache package. Runtime is O(N log N).

Reports include polygon area (including portions outside the lattice), sample
count, target, achieved minimum/median, hmin yield fraction, CFL floor range,
target-below-floor fraction, and minimum permitted external dt at configured
Cr and at Cr=1, regionally and globally. Statistics use lattice sample centres,
not fractional cell coverage. Unsampled and dry-only dt statistics are null.
Depth must be finite, nonnegative and positive-down. CFL never clips targets.
A 30 m target at 10 m depth permits about 3.03 s at Cr=1, or 1.36 s at Cr=.45.
FVCOM uses one global external step, so a small region can limit the entire run.
Regions smaller than lattice spacing may have no samples: inspect sample_count
and refine the input lattice if necessary.

Run `python notebooks/391_sizing_report.py mesh.14 recipes/sizing/tokyo_bay.yaml`
for per-element CSV and a global dt on stderr. Mesh coordinates must be lon/lat,
depths positive-down, with land boundary segments. Measured dt uses minimum
triangle altitude and maximum vertex depth. Constraint attribution reconstructs
coast + grade*distance, maxel, CFL, regions, and an element-adjacency gradation
limit; it is explicitly an estimate, not proof of historical binding. A mesh
and recipe cannot recover the builder's shoreline, feature field, mouth zones,
dilated-depth guard or boundary-corridor history. Exact historical attribution
requires those intermediate fields. Use NQSV for large-mesh diagnostics; this
script does not submit jobs. No mesh build is needed for the unit tests.

The supplied OCTOPUS wrapper runs the default certified mesh diagnostic:
`qsub jobs/octopus/391_sizing_report.sh`. It writes CSV and diagnostics to
job-specific files in `logs/`. Override `SIZING_MESH` / `SR_SIZING` through the
batch environment for another case.
