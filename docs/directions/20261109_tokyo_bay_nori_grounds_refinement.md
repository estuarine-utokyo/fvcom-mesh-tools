# Direction: refine the Tokyo Bay production mesh over the nori fishing grounds

- Date: 2026-11-09
- Requested by: the owner (TB-FVCOM dye session, 2026-11-09)
- Owning repository: `fvcom-mesh-tools` (the mesh and its QA). The model inputs on the new mesh
  are TB-FVCOM work (hydro first, then dye); see the companion direction
  `TB-FVCOM/dye/docs/directions/20261109_nori_fertilization_dye.md`.
- Related: `docs/USER_GUIDE.md` (local refinement), `docs/local_refine.md`,
  `docs/refine_coast_and_bathy_design.md`, and the demonstration recipes
  `recipes/refine/futtsu_nori*.yaml`, `banzu_nori.yaml`, `futtsu_two_beds.yaml`.

## Goal

A locally refined version of the TB-FVCOM production mesh in which the nori fishing grounds
(区画漁業権 areas used for nori culture) are resolved finely enough to:

- inject a fertiliser-like tracer inside each ground (dye area sources), and
- resolve the concentration field over the grounds.

Everything outside the regions and their transitions must stay bit-for-bit the production mesh.

## Base

The production case of the `current` hydro baseline, exactly as `banzu_nori.yaml` declares it:
`$WORK_DIR/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat`, depths
`TokyoBay_dep_m7001tp_rfac0p2_cap300.dat`, OBC `TokyoBay_obc.dat`. Write paths from
`$WORK_DIR` (the existing recipes still say `~/Github/...`; do not copy that pattern into new
recipes).

## Inputs needed from the owner (stop and ask before meshing)

1. **The fishery-right polygons** (区画漁業権, the nori grounds): file, source, licence and
   coordinate system. The existing geometries are DEMONSTRATION rectangles, not
   fishery-right boundaries. Put the real file under `$DATA_DIR` (shared input data), not in
   this repository, unless its licence allows it here.
2. **Which grounds** to refine (for example Futtsu, Banzu/Kisarazu, Funabashi/Sanbanze, others
   in Chiba), or all of them.
3. **Target resolution** inside the grounds (the demonstrations use 30 m), and the acceptable
   growth of the mesh and of the FVCOM time step (the demonstrations show a strip like Banzu
   roughly doubling the node count).

## What to do

1. Write a recipe `recipes/refine/tokyo_bay_nori_grounds.yaml` with one region per ground,
   read from the polygon file with `geometry: {file: ..., where: {...}, buffer_m: ...}`.
   Keep `coastline: preserve` and `rfactor_limit: base`, as in the demonstrations, unless the
   owner decides otherwise.
2. Run the workflow on OCTOPUS (`jobs/octopus/refine_workflow.sh`, inside batch jobs, never
   on a login node; `octopus` skill and watcher after every qsub). That covers QA,
   `finish_depths`, and the 2-day FVCOM smoke of base and refined mesh.
3. Deliver with the mesh:
   - `TokyoBay_grd.dat`, the depth file and `TokyoBay_obc.dat` of the refined mesh;
   - the frozen-zone bit-identity report;
   - the CFL/time-step estimate;
   - **the old-to-new node map** (`node_map`, base node -> refined node, -1 where a base node
     was removed). Every node-indexed TB-FVCOM input must be transferred with it: river and
     sewer nodes in the river maps, MPOS and analysis stations (for example the dye bar nodes
     542, 1556, 1573, 2107, 2316, 2663), OBC nodes, mud IC, and the groundwater channel lists.
     Report every such node that falls inside a refined region, since those must be re-chosen,
     not mapped.
4. Report the result to the owner: node and element counts, time step, smoke-run outcome, and
   any ground the tool could not mesh as asked. The adoption of the mesh into TB-FVCOM is a
   separate, owner-approved step.

## Constraints

- Commit on the default branch; never push. Stage files explicitly.
- Generated meshes and run output stay in the repository's ignored output directories.
- Paths from `$DATA_DIR` / `$WORK_DIR`.

## Done when

- The refined mesh and its node map are delivered with the reports above, the smoke run passes,
  and the owner has the numbers to decide on adoption.
