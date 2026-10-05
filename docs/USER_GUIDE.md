# Local refinement: user guide

This guide is for someone who wants a finer mesh in one part of an existing
FVCOM model -- a port, a fishery, a stretch of coast -- without rebuilding the
rest. It says what the tool does, how to run it, how to check what it made,
how to finish the depths, how to test the result in FVCOM, and how to improve
the tool when a new mesh shows it something it cannot yet do.

The design records behind it are `docs/local_refine.md` (the refinement),
`docs/refine_coast_and_bathy_design.md` (coastline and depths) and
`docs/linear_structures_design.md` (breakwaters and piers). Read them when a
rule below needs its reason.

---

## 1. What the tool does

You declare one or more **regions** and a target element size for each. The
tool cuts a hole in the base mesh around them and fills it again:

| zone | what happens |
|---|---|
| **region (core)** | meshed at the target size |
| **transition** | grades from the target up to the base's own size; its width is worked out for you |
| **frozen** | everything else -- **bit-for-bit unchanged**: nodes, depths, elements, open boundary. The tool checks this on the written file and refuses to report success otherwise |

With the `hires` option, the region and its transition also take their
**coastline from OSM** and their **depths from the M7001 survey** (then the
30 m grid, the Kanto blend, then an extrapolation whose area is reported).

What the mesh size can carry is decided by the size itself:

| feature | rule | becomes |
|---|---|---|
| land at least half an element wide (0.5 h), longer than h | resolvable | **land** -- its outline is coastline |
| land thinner than 0.5 h, longer than h (a thin pier, a breakwater) | too thin to mesh round | a **wall**: a line of mesh edges that water cannot cross |
| anything shorter than one element | not resolvable | dropped, and reported |
| water narrower than two elements (2 h) | cannot hold an element with neighbours | closed (made land) |
| a pier end 0.5-0.75 h wide | its end edge would make angles under 30 deg | the end closes onto one point |

`h` is the **local** element size: the target inside the region, the growing
transition size outside it. The half-element land rule is applied where the
elements are up to twice the target; further out the coarse transition keeps
the two-element rule.

A **wall** is a zero-width structure. In FVCOM an ordinary interior edge does
not stop water, so the tool splits the mesh along the wall: the nodes on it
are duplicated so that each side is boundary, exactly like a coastline
(`docs/linear_structures_design.md` §2, tested in FVCOM in §7).

Every result passes through the project's 22-check FVCOM acceptance gate
(`fmesh-mesh-qa`). The number that matters is **violations introduced by the
patch**: the base mesh may carry its own (the Tokyo Bay base has one element
at 29.0 deg), and the refinement must add none.

---

## 2. Requirements

### 2.1 Repositories

| what | where | access | needed for |
|---|---|---|---|
| **fvcom-mesh-tools** (this repository, GPL-3.0-or-later; Apache-2.0 before 2026-10-02) | https://github.com/estuarine-utokyo/fvcom-mesh-tools | public | everything |
| **oceanmesh, the laboratory's fork** (GPL-3.0) | https://github.com/estuarine-utokyo/oceanmesh, branch `main` (tested at `76903e3`) | public | the fill in `fmesh-refine` |
| **xcoast** | https://github.com/estuarine-utokyo/xcoast | public | making the OSM land polygons for a new area (§2.3); land in some figures |
| **a base FVCOM model** | e.g. `TB-FVCOM` (`input/goto2023/grid/`), the laboratory's model repository | laboratory only | the recipe's `base_mesh`, `base_depth`, `base_obc` -- an FVCOM grd/dep/obc **in UTM zone 54N metres (EPSG:32654) with exactly one open boundary arc** (§2.5) |
| **FVCOM** | the laboratory's `FVCOM` repository, branch `uk-fabm/v5.1.0-dev`, built | laboratory only | the FVCOM tests (§9) only |

**The oceanmesh fork is required.** The oceanmesh on PyPI or conda-forge, and
the upstream CHLNDDEV repository, will not work. The refinement hands the
mesher fixed points and fixed edges (`pfix` / `egfix`: the coastline rim and
the walls), and only the fork forces those through a constrained Delaunay
triangulation (its README §6.3). The generator also looks for the fork at
**`~/Github/oceanmesh`** first (`notebooks/420_local_refine.py`), so clone it
there.

### 2.2 Environment

Python 3.12 in a conda environment built from this repository's
`environment.yml`, conda-forge only, and named after the repository
(`fvcom-mesh-tools`). The three repositories above are
installed editable, without letting pip fetch anything:

```bash
# on a login node (compute nodes have no network)
mamba env create -n fvcom-mesh-tools -f environment.yml
mamba activate fvcom-mesh-tools
git clone https://github.com/estuarine-utokyo/oceanmesh.git ~/Github/oceanmesh
git clone https://github.com/estuarine-utokyo/xcoast.git ~/Github/xcoast
(cd ~/Github/oceanmesh && pip install -e . --no-deps --no-build-isolation)  # compiles C++
(cd ~/Github/xcoast && pip install -e . --no-deps --no-build-isolation)
pip install -e . --no-deps --no-build-isolation                          # this repository
```

- **Compiling oceanmesh:** it compiles C++ (CGAL). On OCTOPUS,
  `qsub jobs/octopus/build_env.sh` does the compile step as a batch job; on
  GENKAI, `pjsub jobs/genkai/build_env.sh`.
- **Other packages:** add them with `mamba install -c conda-forge ...`,
  never with pip. The `pip install -e` lines above are the only use of pip.
- **Editable install:** `fmesh-refine` needs this repository installed
  editable, because it runs `notebooks/420_local_refine.py` from the
  checkout.

### 2.3 Data

Everything below is read from `DATA_DIR` (on OCTOPUS
`/octfs/work/G16445/share/Data`, the laboratory's shared data area):

| data | path under `DATA_DIR` | used for | without it |
|---|---|---|---|
| OSM land polygons | `geodata/OSM/coastmask_cache/<area>/land.shp`, made by xcoast (OSM, ODbL) | the coastline, the width filter, walls | `fmesh-refine` stops; build the cache for a new area with xcoast, or pass `--land` |
| M7001 (JHA) on T.P. | `geodata/bathymetry/M7001/TP/M7001_dem_tokyobay.nc`, `M7001_TP.parquet` | `hires.bathymetry: tokyo_bay` | use `bathymetry: base` (inherit the base's depths) |
| Tokyo Bay 30 m grid, Kanto blend | `geodata/bathymetry/tokyo_bay/...` | the same depth ladder | as above |

M7001 is licensed survey data and is not public. The depth ladder
(`fvcom_mesh_tools.dem.tokyo_bay`) covers Tokyo Bay only. Another area needs
`bathymetry: base`, or a ladder of its own.

### 2.4 OCTOPUS

- **Login nodes** have the network and no heavy computing: install, clone
  and submit there.
- **Compute nodes** run everything else, through `qsub`. The job scripts are
  in `jobs/octopus/`. GENKAI (`pjsub`) has its own site layer in
  `jobs/genkai/` (`build_env.sh`; `extend_check.sh`, the end-to-end check of
  the extension tools). Both source `jobs/common_core.sh`, the
  machine-independent set-up; only job ids, conda location and library
  defaults differ.
- **Paths come from the environment, set once in the login profile:**
  `$DATA_DIR` (input data) and `$WORK_DIR` (repositories, conda, scratch,
  runs); a job stops if either is unset. The FVCOM executable is
  `$WORK_DIR/Github/FVCOM/src/fvcom` unless `FMESH_FVCOM` names another,
  and its libraries are `FVCOM_LIBS` (default in `jobs/octopus/common.sh`).
  The only account-specific setting left in the scripts is the accounting
  group, `#PBS --group=G16445`, which another group changes.
- **Monitoring:** start a monitor after every `qsub`, and read the log when
  the job ends.

### 2.5 What the base mesh must be

Two limits of the generator:

- **Coordinates in UTM zone 54N metres (EPSG:32654).** The generator
  projects OSM, the regions and the depth products into that frame. A base
  in another CRS must be reprojected first. **Making sure of this is the
  user's responsibility.** The files carry bare numbers, so the check before
  meshing can only refuse coordinates that are implausible read as
  EPSG:32654 (degrees, or a centre far from 136-146 E). A base in a
  neighbouring UTM zone looks plausible and is not caught.
- **Exactly one open-boundary arc.** The boundary lists of the refined mesh
  are rebuilt around a single arc. This one is checked before meshing, and
  a base with zero or several arcs is refused.

Lifting either is future work, not a setting.

## 3. Quick start (OCTOPUS)

1. **Write a recipe.** Copy one and change the region (§4):

   ```bash
   cp recipes/refine/kimitsu_port_hires.yaml recipes/refine/my_port.yaml
   ```

2. **Submit the whole workflow** from the repository root on a login node:

   ```bash
   FMESH_VIEWS="port:391900:394100:3908350:3910600" \
     bash jobs/octopus/refine_workflow.sh recipes/refine/my_port.yaml
   ```

   It submits a chain of jobs and prints their ids and log files:

   | stage | job | what it makes |
   |---|---|---|
   | refine | `417_hires_refine.sh` | the mesh, QA and report in `outputs/refine_my_port/` |
   | figures | `418` + `notebooks/434_final_mesh.py` | `final_mesh_*.png` |
   | depths | `421_finish_and_run.sh` | the finished depths (`fmesh-finish-depths`) and an M2 test pair against the base |
   | smoke | `423_m2_smoke.sh` | a 2-day FVCOM run of the base and the refined mesh |
   | M2 (optional) | `412` x 2 + `413` | 20-day M2 runs and the tide-gauge comparison; set `FMESH_M2=1` |

   The depth product is set by `FMESH_HMIN` (3), `FMESH_HMAX` (300) and
   `FMESH_RFACTOR` (0.2).

3. **Monitor it** with the command the script prints. When it ends, read the
   logs in order (§6) and look at the figures (§7).

   NQSV starts a job when the one before it **ends, whatever the outcome**.
   So each stage writes a marker only when it passed, and the next stage
   refuses to start without it:

   | marker | written by | only when |
   |---|---|---|
   | `ACCEPTED` | the refinement | every gate passed (see below) |
   | `STAGED` | the depths stage | the depth product converged |
   | `SMOKE_OK` | the smoke test | both runs pass `fmesh-check-run` |
   | `RUN_OK` | each M2 run | the run passes `fmesh-check-run` |

   `ACCEPTED` needs the frozen contract, the resolution, the written case and
   0 QA violations introduced. A stage that stops with "no ... marker" means
   the stage before it failed; read that stage's log.

To run one step by hand instead, the same commands exist on their own (§10).
`fmesh-refine` runs the generator; run it inside a batch job, not on a login
node.

---

## 4. Writing a recipe

A minimal `hires` recipe:

```yaml
base_mesh:  ~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat
base_depth: ~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_dep_m7001tp_rfac0p2_cap300.dat
base_obc:   ~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_obc.dat

dt_expected_s: 4.5      # the external step the model runs at (advisory)
gradation: 0.165        # how fast element size may grow away from the region

hires:
  coastline: resolve    # resolve (OSM) | preserve (keep the base coastline)
  bathymetry: tokyo_bay # tokyo_bay (M7001 ladder) | base (inherit)
  scope: hole           # hole (region + transition) | core
  blend: ramp           # ramp | none -- how new depths meet the frozen ones
  # the coastline rules below are all on; any one can be switched off:
  # rim_repair: false | continuous_width: false | keep_base_land: false |
  # wall_pockets: false

refine:
  - name: my_port
    geometry:
      circle: {center: [139.8228, 35.3230], radius_m: 900}
    target_h_m: 30
    priority: 0
```

| key | required | meaning |
|---|---|---|
| `base_mesh` | yes | the base FVCOM `_grd.dat` (or a fort.14) |
| `base_depth` | with a `.dat` base | which depth file the model runs with -- the grd's own column may be a different product |
| `base_obc` | no | the open-boundary node list |
| `dt_expected_s` | yes | the model's external step; a region too fine for it raises an **alert**, not an error |
| `gradation` | yes | the size growth rate outside the region; it sets the transition width |
| `hires` | no | present = the coastline/bathymetry branch above; absent = the default branch, which keeps the base coastline and depths (`coastline: preserve/resample/spline`, `rfactor_limit`, `coastline_tolerance_m` apply there). `hires.coastline: resolve` follows OSM (sharp corners cut, islands and walls added); `preserve` keeps the base coastline exactly |
| `hires.rim_repair`, `hires.continuous_width`, `hires.keep_base_land`, `hires.wall_pockets` | no | the coastline rules (below); each `true` by default, `false` switches it off |
| `refine` | yes | one or more regions |
| `refine[].geometry` | yes | `circle: {center: [lon, lat], radius_m}`, a `bbox`, a GeoJSON `Polygon`, or `{file: area.geojson, where: {...}, buffer_m: 25}` |
| `refine[].target_h_m` | yes | the target element size, m |
| `refine[].priority` | no | no effect in local refinement: where regions overlap the finest size wins (`patch_sizing`); kept for sizing recipes |

Unknown keys are errors, and relative paths resolve against the recipe's own
directory. Before meshing, the run checks each region (depth, dryness, time
step) and prints what it found.

**Coastline rules (all on by default).** Four rules turn the OSM coastline
into one a mesh of the local size can carry. Each was written for a port
where the earlier rules failed, each was tried on every recipe before it
was made the default, and each can be switched off in a recipe
(`hires.<rule>: false`), which restores the behaviour before it exactly.

| rule | what it does | since |
|---|---|---|
| `continuous_width` | in the coarse zone (elements over twice the target), water is judged at the local element size itself, not at the lower bound of an octave band: water narrower than 1.5 elements that ends at one body of land is closed (a strait between two bodies is left). Distance transform on a raster of one third of the finest target (10 m for a 30 m target), smoothed; it only adds land to what the band filter keeps. The fix at the source for the band seams (§11) | 2026-09-27 |
| `keep_base_land` | in the coarse zone the base decides what the width filter would change: a piece of land the filter would remove stays land -- the whole piece -- if the base mesh has most of it as land, and a piece of water it would fill stays water if the base has most of it as water (the transition serves efficiency, not detail; owner). Water was added when a 240 m band closed an Odaiba channel that the project's own base draws, and the seam with the frozen channel beyond failed QA | 2026-09-27, water 2026-09-28 |
| `wall_pockets` | water shut in by walls alone, off the coast, in the coarse zone, narrower than two elements and at least one element in area, becomes an island; its walls come back if the rim refuses the island | 2026-09-27 |
| `rim_repair` | before the fill, the finished coastline is checked against the local size and repaired, in order, up to two rounds: (1) a point beside an edge under half an element is removed, if that crosses nothing, moves little and sharpens no water angle below 60 deg -- otherwise the edge is reported; (2) a point within one element of a coast it is not next to is a throat -- the dead end beyond becomes land if no element fits in it, otherwise (a pier tip nearly touching the quay across) the tip steps back to one element; (3) water corners under 60 deg are cut again. Frozen points and wall roots never move. If every seed fails, the rim is repaired once more near the best seed's offenders and the seeds are tried again | 2026-09-26 |

With all four on, all seven recipes are accepted and none names a rule.
Three more rules (`water_wedges`, `short_chords`, `seam_water`) worked round
the band seams before `continuous_width` fixed them, and were removed on
2026-09-27 (commit 4a8c503 has them). `LR_EXPERIMENTAL` (job 417) overrides
the set for a trial run.

**Choosing a target.** The external time step scales with the smallest element
over the square root of the depth. The run tells you the step the target
allows. For example, 30 m elements over 8.65 m of water allow 2.8 s against a
model running at 4.5 s, so the refined model costs 1.6x the steps.

---

## 5. What a run does

In order:

1. select the hole and read the base;
2. filter the OSM shoreline at the local size, and extract walls from what
   the filter removes;
3. build the coastline rim: add islands, cut sharp corners, root walls on the
   coast;
4. fill the hole with DistMesh (oceanmesh), for **several random seeds**;
5. for each seed:
   - stitch the fill to the frozen mesh and split it along the walls;
   - repair it;
   - take its depths;
   - verify the frozen contract;
   - write and run QA;
6. keep the seed with the fewest violations introduced.

**Islands and lakes (step 3).** Whether a feature is resolved at all is
decided once, by the shoreline filter at the local element size (step 2):
what the size cannot carry is not resolved, and in the transition zone
efficiency comes before detail. Step 3 then delivers what the filter kept.
Each island's outline, and each lake's, is either resampled at the local
size or kept as the source draws it, and one outline's choice can crowd a
neighbour. `island_rings` searches these choices -- every combination for
up to eight lakes per island, a bounded search beyond, and a repeat of the
pass preferring source outlines next to anything refused -- and keeps the
result that loses least area. It is a best effort, not a proof of the
optimum: whatever it still cannot deliver is named in `report.json`
(`islands_added.skipped`, and `tight` for a feature kept closer than half an
element, which the QA gate then judges), and land left as water stops the run.

A port-sized region takes 5-10 minutes on one core.

**Reproducibility.** The seeds are fixed (`LR_SEEDS`, default 0-4), and the
oceanmesh fork seeds its random generator with them, so the same inputs give
the same mesh: rebuilding the five port recipes on 2026-09-25 gave
byte-identical fort.14 files, and a run of seed 3 alone gave the same mesh as
a search over 0-4 that chose 3. "The same inputs" means more than the seed.
`report.json` records it all under `provenance`:
- the commit of this repository and of the oceanmesh fork, with any
  uncommitted or untracked files (the run also warns about them);
- SHA-256 of the recipe, base mesh, depth and OBC files, the OSM
  shapefile, the depth products the ladder read, and any region file;
- the Python and library versions (from package metadata; nothing is
  imported for it), the opt-in rules used, and every `LR_*` / `FMESH_*`
  setting and `DATA_DIR`.

To remake a mesh, check out both commits, use the same inputs and settings,
and set `LR_SEEDS` to the list the report records under
`provenance.seeds`. The accepted `seed` alone is enough only when
`search_pass` is 1: a mesh from the retry near the offenders (`search_pass`
2) depends on which seed failed best in the first search, so the whole list
is needed. Identical results across machines or library versions are not
guaranteed: floating-point differences can move a node.

---

## 6. Outputs and how to read them

`outputs/refine_<name>/`:

| file | what it is |
|---|---|
| `<base>_<name>.14` | the refined mesh (fort.14) |
| `<base>_<name>_qa.json` | the acceptance gate, check by check, with offender coordinates |
| `report.json` | everything the run decided and measured |
| `fvcom/` | the FVCOM case (grd, dep, obc, cor), depths as the source gives them |
| `fvcom_finished/` | after `fmesh-finish-depths`: the case with finished depths, plus `<case>_dep_<variant>.dat` |
| `node_map.npy` | base node -> refined node, the frozen contract's evidence |
| `<base>_<name>_walls.json` | the node pairs a wall duplicated on purpose, bound to the `.14` by its SHA-256; `fmesh-mesh-qa` reads it, so the delivered mesh gets the same verdict when checked again. Copy or rename it together with the mesh, or point to it with `--wall-pairs`. A renumbered mesh needs the pairs mapped through the new numbering and written again (`walls.write_wall_pairs`); the old file is refused |
| `ACCEPTED` | written last, only when every gate passed (§3) |
| `shoreline_filtered.shp` | the OSM land after the width filter |
| `fill_constraints.npz`, `walls_stages.npz` | the rim and walls the mesher was given, for inspection |
| `final_mesh_*.png` | the figures (§7) |

The run log ends with lines like these:

```
[lr]     seed 0: QA 20/22, 0 introduced by the patch
[lr] selected seed 0 (search pass 1); the gates follow
[lr] achieved in my_port: median cell 28.5 m against a 30 m target, 100.0 % of the water within 1.25x
```

- **`0 introduced`** is the goal.
- **`selected seed`** is the best candidate, not yet the result: the gates
  after it (resolution, depths, the written file's QA) decide, and only a
  run that passes them all writes `ACCEPTED`.
- **The two failed checks.** On the Tokyo Bay base they are the base's own
  element 2101 and `min_depth_clip`. On the hires branch the depths are the
  source's until they are finished, so `min_depth_clip` is reported and not
  gated.
- **`no seed produced a mesh that keeps the frozen-zone contract`** means no
  seed passed. `report.json` holds each attempt and why it failed.
- **A run that exits with violations introduced** still leaves its files for
  diagnosis, but no `ACCEPTED`, so the chain goes no further.
- **An output directory that is not empty** is refused (`fmesh-refine`, job
  417, `refine_workflow.sh`). Move it aside first.
- **One run per output directory.** A run claims its directory by creating
  `.reserved` atomically, so of two runs started together, one is refused.
  The file stays after the run. If a run crashed, recover like this: make
  sure it has stopped, then move the WHOLE directory aside (keep it for
  diagnosis) and start again. Do not delete only `.reserved`.

---

## 7. Looking at the mesh

Every mesh figure uses **fixed colours**:

- **black**: every solid boundary -- coastline, quay, and a wall
  represented as a line;
- **red**: the open boundary;
- **thin grey**: interior edges.

| figure | how | use it to |
|---|---|---|
| final mesh | `fmesh-plot-views <out dir> --view name:x0:x1:y0:y1` (or `final_mesh_*.png`) | see the delivered mesh, whole and in close-ups |
| mesh over the raw OSM | `notebooks/433_before_after.py <before dir> <after dir>` (job 418 with `FMESH_SCRIPT=433_before_after.py`, `FMESH_ARG2`) | check that every pier, quay and breakwater is where OSM has it; compare two runs; it also prints each wall's angle to its pier |
| QA map | `notebooks/429_wall_qa_map.py <out dir>` | find where the QA offenders are |
| zoom on constraints | `jobs/octopus/427_constraint_zoom.sh` with `FMESH_POINTS=x:y:half+...` | see the fixed points, constrained edges and angles around one offender |

A valid mesh can still be wrong about the port. The QA gate cannot see a
quay block meshed as water, or a breakwater a few degrees off. So always look
at the mesh over the raw OSM.

---

## 8. Finishing the depths

The run leaves the depths as the source gives them -- a tidal flat can be
-4.7 m. `fmesh-finish-depths` makes a run-ready depth file the way
`TB-FVCOM/input/goto2023/grid` makes its variants:

```bash
# a refinement: only the patch's depths move; the frozen base depths are kept
fmesh-finish-depths outputs/refine_my_port --hmin 3 --hmax 300 --rfactor 0.2

# any FVCOM case, every node, from a chosen depth file
fmesh-finish-depths ~/Github/TB-FVCOM/input/goto2023/grid/TokyoBay_grd.dat \
    --dep TokyoBay_dep_m7001tp_raw.dat --hmin 3 --hmax 300 --rfactor 0.2 --outdir out/
```

| option | meaning | in the file name |
|---|---|---|
| `--hmin` | minimum depth, m | `min3m` |
| `--hmax` | maximum depth, m (default 300) | `cap300` |
| `--rfactor` | limit on \|h_i - h_j\| / (h_i + h_j) over every edge (default 0.2; 0 = no smoothing) -- the sigma-coordinate stability condition | `rfac0p2` |
| `--whole-mesh` | on a refinement, move the frozen depths too | |
| `--method` | `equal` (default): TB-FVCOM's smoother -- floor, smooth the uncapped field, then cap; `limit`: the refinement's limiter | |
| `--allow-unconverged` | write the product even if the limit was not reached (diagnosis only) | |

`--rfactor` must be 0 or a number in (0, 1). If the limit is not reached --
too few `--rounds`, or a frozen depth the floor cannot meet -- the command
writes nothing and exits 3. So a chain never stages an unfinished product.

The limit is judged on every edge with at least one end the command may
move, and on every edge between two frozen nodes that the patch CREATED
(`n_over_frozen_pair_new`). The patch may connect two retained nodes that
were never neighbours. Only an edge the base mesh already had (read from
the base named in `report.json`) is excused as the base's own
(`n_over_frozen_pair_inherited`). The base is identified by the SHA-256
recorded in `report.json` when the refinement ran. When that base cannot be
read, has changed since, or the report is older and has no hash, no frozen
pair is excused. So keep the base files as they were: do not overwrite a
base that a refinement used.

The output is `<case>_dep_min3m_rfac0p2_cap300.dat`, with a JSON report of
what each step moved. From `TokyoBay_dep_m7001tp_raw.dat` it reproduces
`TokyoBay_dep_m7001tp_rfac0p2_cap300.dat` at every node.

---

## 9. Testing in FVCOM

`jobs/octopus/421_finish_and_run.sh` stages an M2 test pair: the base and the
refined mesh, with the same forcing, sponge, sigma levels and external step
(the finer mesh's). The staging is written for the Tokyo Bay goto2023 base.
Then:

- **`423_m2_smoke.sh`** runs both for 2 days. Each run is judged by
  `fmesh-check-run`, which requires all of:
  - `TADA` in the log, and no fatal message;
  - output that reaches the namelist's `END_DATE`;
  - finite `zeta`, `ua` and `va` in every record.

  A zero exit code alone is not success: some FVCOM STOP paths return 0.

  In detail, the history output is `<casename>_0001.nc`, `_0002.nc`, ...
  numbered without a gap; restart and other files are not history. It must:
  - carry `zeta`, `ua`, `va` and `Times`;
  - have times that increase record by record, with no gap over 1.5 times
    the namelist's `NC_OUT_INTERVAL` (seconds, minutes, hours, days, or
    `cycles` of the internal step). A declaration that cannot be read fails;
  - start at `NC_FIRST_OUT` (else `START_DATE`) and reach `END_DATE`, each
    within one interval (one hour when none is declared).

  Namelist comments (`!` outside quotes) are ignored. A hot start is judged
  from its own START_DATE and NC_FIRST_OUT. That is a check of the checker,
  not a tested hot-start workflow.

  `412_m2_run.sh` deletes the case's `output/*.nc` before it runs, so an
  earlier attempt's output is never judged as this one's. Copy anything you
  need from there first.

  A stage removes its own marker and every later one before it starts, and
  writes its own only on success. `RUN_OK` needs the solver's exit status
  and the check to pass. A marker left over from an earlier attempt in the
  same directory therefore never lets the next stage start.
- **`412_m2_run.sh` x 2 + `413_m2_analysis.sh`** run 20 days of M2 and
  compare the amplitude and phase at the tide gauges.

A refinement should change the far field by a fraction of a millimetre (the
Kimitsu port: -0.2 to -0.4 mm, +0.013 deg). Inside the region, look at the
M2 amplitude map (`notebooks/431_harbour_currents.py`): a node stuck at
exactly zero amplitude is a defect.

---

## 10. Command reference

| command | what it does |
|---|---|
| `fmesh-refine RECIPE [--out] [--seeds] [--land]` | the refinement (runs `notebooks/420_local_refine.py`) |
| `fmesh-finish-depths SOURCE --hmin --hmax --rfactor` | the depth product (§8) |
| `fmesh-plot-views MESH [--view ...]` | mesh figures in the fixed colours |
| `fmesh-mesh-qa MESH` | the 22-check FVCOM acceptance gate, on any fort.14 (reads `<stem>_walls.json` when present) |
| `fmesh-check-run RUN_DIR` | did an FVCOM run finish with usable output (§9) |
| `fmesh-refine-depths` | the earlier depth finisher (floor, cap, then the refinement's limiter); kept for old scripts |
| `jobs/octopus/refine_workflow.sh RECIPE` | the whole chain as batch jobs (§3) |

The generator lives in `notebooks/`, not in the package. This split dates
from when the package was Apache-2.0 and could not import the GPL-3.0
oceanmesh; since 2026-10-02 the package is GPL-3.0-or-later and may, but
the split is kept. Its parts -- selection, rim, walls, stitching,
verification, QA -- are in the package and have unit tests.

---

## 11. Known limits

- Tested on one base mesh (goto2023): the Kimitsu port, the Futtsu coast, an
  offshore fishery, the default branch, the Tokyo port at Odaiba
  (`recipes/refine/tokyo_odaiba_hires.yaml`, with the Daiba islands) and
  Funabashi port (`recipes/refine/funabashi_port_hires.yaml`) and Yokohama
  inner harbour (`recipes/refine/yokohama_port_hires.yaml`, a 2.3 km
  circle). Odaiba
  needed two coastline rules Kimitsu had not shown -- a new place will find
  new cases (§13).
- A structure hugging the coast within 0.4 element, thinner than half an
  element, is not represented.
- **Band seams.** The coastline filter judges width in octave bands of the
  element size, each at its LOWER bound, so water two lower bounds wide --
  only one element where the elements are twice the bound -- survives. Where the next band closes a
  channel that this band keeps, the join is a straight cut across the
  channel, and its corners can be acute (Funabashi, 1.4 km west of the
  region). The fix at the source -- judging each channel at its own local
  size -- is `continuous_width`, on by default since 2026-09-27 (§4).
  Judged at two whole elements it closed 1.26 km2 of the Odaiba port; at
  1.5 elements (the bands' average), leaving straits open, every recipe
  is accepted with it. Quarter-octave bands were tried first and moved the
  Kimitsu transition coast by 211 m.
- OSM `man_made` lines (breakwaters mapped only as lines) are not used.
- The result depends on the DistMesh seed; the best of several is kept.
- The M2 staging in job 421 is specific to the Tokyo Bay goto2023 base.
- The base must be in EPSG:32654 with one open-boundary arc (§2.5).
- With `hires.coastline: preserve` the coastline is kept exactly, including
  corners sharper than 60 deg. These are reported (`acute_corners_kept` in
  `report.json`), not cut. A node left in one element at such a corner is
  opened by the repair, and the QA gate says if that failed.

---

## 12. Building the base mesh

The local refinement needs a whole-bay base. The Tokyo Bay base this project
built -- `TokyoBayTool`, which follows the goto2023 mesh (SMS and hand work)
closely -- is made from a recipe and the raw data alone:

```bash
# on a login node, from the repository root
qsub -v FMESH_RECIPE=recipes/base/tokyo_bay_tool.yaml jobs/octopus/440_base_mesh.sh
# -> outputs/base_tokyo_bay_tool/TokyoBayTool{.14,_grd,_dep,_obc,_cor}.dat
```

It takes about six minutes on 16 cores. The steps
(`notebooks/440_base_mesh.py`):

1. **Land.** OSM land polygons minus the inland water connected to the sea,
   for the recipe's window, from `DATA_DIR` (`geodata/OSM/...`). A missing
   source stops the run; it is never replaced by a download.
2. **Generation.** `notebooks/325_sample_repro.py`: oceanmesh DistMesh with
   the open boundary constrained, a two-zone Courant sizing on SRTM15, the
   hand-drawn geometry corrections of `recipes/edits/sample_repro/`, and the
   OSM waterways.
3. **Finishing.** `notebooks/331_finish2.py`: narrow-channel policy,
   open-boundary finishing, the coastline fit.
4. **Depths and the FVCOM case.** `notebooks/422_tool_base.py`: M7001 on
   T.P., 3 m floor, r-factor 0.2, 300 m cap -- the production recipe of the
   hydro baseline -- then QA.

**The recipe** (`recipes/base/tokyo_bay_tool.yaml`) names everything:

| key | what |
|---|---|
| `open_boundary` | the open boundary, an input: `lon,lat` of its nodes in order (`tokyo_bay_obc.csv`, goto2023's 13 nodes). Mesh nodes are constrained onto it |
| `domain` | the closure of the meshing domain around the open boundary, its bbox, and where the coast size at the southern closure is read (`tokyo_bay_domain.json`) |
| `land` | the OSM window and the smallest inland water kept |
| `edits` | the directory of hand-drawn corrections, applied in file-name order |
| `settings` | every generation and finishing setting, written out -- a default changed in the code must not change the mesh -- including the seeds (`SR_GEN_SEED`, `SR_FIN_SEED`) |
| `depths` | the depth product (`m7001_production`) |
| `reference` | the mesh the recipe reproduces: its fort.14 hash, node and element counts, open-boundary node count, wet area, and the relative tolerances |

**What "reproduce" means.** The environment is not pinned: it follows the
latest conda-forge releases, and the oceanmesh fork evolves
upward-compatibly -- a new feature must leave the existing examples
substantively reproducible. A mesh may therefore change slightly (and may
get better). The build **reproduces** the reference when every QA gate
passes, the open boundary keeps its node count, and the node count, element
count and wet area stay within the recipe's relative tolerances (5 %, 5 %,
1 %). Otherwise `440_base_mesh.py` exits with status 3. Byte identity is
reported separately.

The seeds are fixed and the thread count is fixed at 16. `report.json`
records:
- the commit of this repository and of the oceanmesh fork, with any changes
  not committed, and the version of every library;
- the SHA-256 of the recipe, its files and every raw input;
- the effective settings and seeds;
- the hashes of the products and the comparison with the reference
  (`reproduction`).

The build of 2026-09-28 (`aa6e1c6`, environment `oceanmesh-bench`)
reproduced all six files of the 2026-09-22 base byte for byte.

**What another user needs.**
- The same `DATA_DIR` files. `report.json` lists their hashes, and a
  different OSM extract gives a different coastline.
- The oceanmesh fork (estuarine-utokyo/oceanmesh), current `main`.
- The `fvcom-mesh-tools` environment, built from `environment.yml` (§1).

M7001 is licensed (§2.3). Without it the mesh can be rebuilt, but not its
depths.

**Another bay** needs its own recipe: an open-boundary file, a domain file,
a land window, and edits of its own. The two-zone sizing seam in notebook
325 (`SR_ZW_LAT`, `SR_ZE_LAT`, and the seam longitudes in the code) is Tokyo
Bay's.

## 13. Extending a base mesh outward

A wider mesh keeps a finished base mesh **as it is** and adds the sea out to
a new open boundary (owner, 2026-09-29). The base's open boundary becomes an
interior line; the base's nodes, elements and depths come through bit for
bit, which the build checks.

```bash
# 1. the open boundary (an input): design it, check it, write its nodes
qsub -v FMESH_DESIGN=recipes/extend/tokyo_bay_enshu_obc_design.yaml,\
FMESH_OBC_CSV=recipes/extend/tokyo_bay_enshu_obc.csv jobs/octopus/444_design_obc.sh
# 2. the mesh (about 25 min)
qsub -v FMESH_RECIPE=recipes/extend/tokyo_bay_enshu.yaml jobs/octopus/445_extend_mesh.sh
# 3. does it run? (2 days, uniform M2 on the open boundary; stability only)
qsub -v FMESH_CASE=outputs/extend_tokyo_bay_enshu/TokyoBayEnshu,\
FMESH_RUN_ROOT=$WORK_DIR/scratch/smoke_<stamp>,WORK_DIR=$WORK_DIR jobs/octopus/448_extend_smoke.sh
```

**The open boundary** (`notebooks/444_design_obc.py`, `obc_design.py`):
orthogonal to the coast at both ends (the coast direction is the chord
3 km either side of the end point; ends sit on long straight coasts),
straight sides, corners rounded with circular arcs, and a node spacing
never below the time-step floor `dt*sqrt(g*H)/Cr`. The design YAML lists the
sides as legs (`bearing`, and `until_lat`, `until_lon` or `until:
end_normal`). The report beside the CSV gives the angles at the ends, the
spacing, the depths and whether the line crosses land.

**The recipe** (`recipes/extend/*.yaml`, `extend_recipe.py`) names the base
case, the open boundary, the land window, the bathymetry sources, every
sizing setting and the depth rules:

| key | what |
|---|---|
| `base`, `base_case` | the finished FVCOM case adopted unchanged |
| `open_boundary` | the new boundary's nodes (from 444) |
| `bathymetry.sizing`, `bathymetry.depths` | source lists in priority order (below) |
| `settings` | `coast_h_m`, `max_edge_m`, `gradation`, `cfl_dt_s`, `cfl_cr`, the band half-widths on the two constrained lines, `lattice_m`, `dm_scale`, the seeds |
| `depths` | `min_m`, `max_m` (null = no cap), `rfactor` for the new nodes |

**Bathymetry sources** (`dem/sources.py`), selectable by name; each point
takes the first source that covers it:

| name | datum | what |
|---|---|---|
| `cao_shutochokka_2025` | T.P. | Cabinet Office nested grids (10-2430 m), finest first |
| `m7001` | T.P. | M7001 soundings and low-tide line, linear between points |
| `m7001_tokyobay` | T.P. | M7001 gridded at ~180 m, Tokyo Bay only |
| `srtm15_kanto`, `srtm15plus`, `gebco_2024` | mean sea level | global grids, to fill what the survey products do not cover |

`dem/sources.py` keeps each source's datum (`DATUM`); `sample` warns when
points take their depth from a source not on T.P., and 447/453 report the
count (`nodes_not_on_tp`). The Cabinet Office grids may not be redistributed
as they are; they are read in place (see their README in `$DATA_DIR`). For dredged pits M7001 is the
authority (the Cabinet Office grids miss or misplace them).

**What a build does** (`notebooks/445` runs `446` then `447`):

1. **Generation** (`446`): the base is land for this stage. The sizing is
   the coast-distance field limited to the gradation, raised to the
   time-step floor (a graded dilation that aims to keep the extension from
   limiting dt; bands along the constrained lines and the final depths can
   still undercut it, which is reported as a warning, not a failure),
   and set to each constrained line's own spacing on a band along it. The
   base's open boundary and the new one are fixed points and edges, each
   with a *ladder* (a second fixed line one local size inside, as the base's
   own open boundary has). oceanmesh's default clean is off -- it deletes
   fixed nodes -- and instead free nodes that DistMesh projected onto a fixed
   line are dropped and the rest re-triangulated with the same constrained
   Delaunay; flat elements are removed.
2. **Finishing and merge** (`447`): the open-boundary finishing chain and
   the coastline fit on the new part; the merge onto the base through the 13
   shared nodes, with the frozen-base check; a repair restricted to the new
   part (cape tips one element wide are dropped, `improve_patch` moves only
   new interior nodes and flips only new elements, never on the open
   boundary; the perpendicularity pass); depths of the new nodes from
   `bathymetry.depths`, floored, optionally capped, and r-factor limited
   with the base depths held fixed; the FVCOM case and QA.

The first Tokyo Bay extension (`tokyo_bay_enshu`): on OCTOPUS (2026-09-29)
14,740 nodes, 27,135 elements, QA 22/22, all new depths from the Cabinet Office
grids, and a two-day FVCOM smoke run finite (max |zeta| 0.47 m, max speed
0.56 m/s at the Uraga strait). The current build (23 QA gates; GENKAI,
2026-10-05, reproduced bit for bit by three jobs) has 14,673 nodes, 27,011
elements and QA 23/23, and the two-day smoke run on it, with the FVCOM rebuilt
on 2026-10-05 (120 ranks, external step 5.625 s, 26 s wall), is finite again
(max |zeta| 0.469 m, max depth-mean speed 0.563 m/s). Its real forcing will
come from JCOPE-T DA re-extracted over the wider domain.

## 14. Astronomical tide on the open boundary

A tide-only run checks a mesh against tide gauges before any other forcing
exists: no wind, rivers or density, only the astronomical tide on the open
boundary (and, optionally, the tidal potential inside the domain).

```bash
# a 200-day run of one case, NAO.99Jb on its open boundary (about 40 min)
qsub -v WORK_DIR=$WORK_DIR,DATA_DIR=$DATA_DIR,\
FMESH_CASE=outputs/extend_tokyo_bay_enshu/TokyoBayEnshu,FMESH_RUN_ROOT=$WORK_DIR/scratch/tide_<stamp>/enshu,\
FMESH_EQUI=1,FMESH_FVCOM=$WORK_DIR/local/fvcom-equi/bin/fvcom jobs/octopus/449_tide_nao_run.sh
# harmonic constants against the gauges (several runs side by side)
qsub -v DATA_DIR=$DATA_DIR,FMESH_RUNS=a=<run>+b=<run>,FMESH_OUT=<dir> jobs/octopus/450_tide_gauge_compare.sh
# maps against the tide model that forces the run; time series at gauges
qsub -v DATA_DIR=$DATA_DIR,FMESH_RUN=<run>,FMESH_OUT=<dir> jobs/octopus/451_tide_model_map.sh
qsub -v DATA_DIR=$DATA_DIR,FMESH_RUNS=a=<run>+b=<run>,FMESH_OUT=<png> jobs/octopus/452_tide_timeseries.sh
```

`qsub -v` cannot pass commas (and passes spaces unreliably): lists in these
jobs are joined by `+`.

**The forcing** (`tide_models.py`, `notebooks/449`). The harmonic constants
of a tide model (NAO.99Jb for now, `$DATA_DIR/tides/models`) are sampled at
every open-boundary node and converted to FVCOM's spectral form
`A cos(2 pi t/T - phi)`, with `t` counted from the file's Time Origin:

- `phi = G - V0 - u` and `A = f * amplitude`, where `G` is the Greenwich
  phase lag, `V0` the astronomical argument at the Time Origin, and `f`, `u`
  the nodal factor and angle, frozen at the middle of the run (FVCOM cannot
  vary them);
- the astronomy is utide's, so a harmonic analysis with utide (450, 451)
  uses the same convention. A test builds an FVCOM series from converted
  constants and analyses it back to the constants it came from.

The usual errors in this step are each much larger than anything else here:
leaving out `V0` (a different phase error for every constituent: 306 deg for
M2 and 295 deg for O1 at 2021-01-01 00 UT), mixing JST constants with a UTC
run (nine hours, about 260 deg of M2), and leaving out `f` and `u` (K1 about
10 %, O1 about 18 % in amplitude in some years).

**The tidal potential** (`--equilibrium`, `FMESH_EQUI=1`). FVCOM's
`-DEQUI_TIDE` adds the equilibrium tide to the pressure gradient. Its own
astronomy for real dates is a monthly approximation; the owner's FVCOM now
also reads `f` and `V0 + u` from each component line of the tide file, so
that the potential uses exactly the open boundary's astronomy
(`tide_models.EQUILIBRIUM`: Cartwright-Tayler amplitudes, elasticity factor
`1 + k - h` after Wahr). That FVCOM is built out of tree by FVCOM's
`octopus/build_fvcom_equi.sh` into `$WORK_DIR/local/fvcom-equi/bin/fvcom`.
On the Enshu mesh (about 300 km across, down to 5.7 km deep) it raises M2 by
0.6-1.0 cm (about 2 %) at the Tokyo Bay gauges and S2 by 0.2-0.4 cm, both
towards the observations; diurnal constituents change by less than 0.1 cm.
On the base mesh every change is under 0.4 cm. It is not essential at these
sizes, but it is right and costs nothing: keep it on for wide meshes.
Without f and V0+u in the file FVCOM uses its own monthly astronomy; that
path used to crash (SIGSEGV, fixed in FVCOM 30cc7c96), and in the channel of
notebook 454 it now agrees with the utide path to 2 mm for a 1.4 cm
potential response. The utide path stays the default here: it is the same
astronomy as the open boundary.

**First results** (2021-01-16 to 07-20, 185 days, 7 gauges, NAO.99Jb,
8 constituents): rms vector difference M2 1.7 cm on the Enshu mesh with
the potential (base 2.3 cm), S2 0.6 cm (1.0 cm). The diurnal constituents
come out 4-9 % too large on the Enshu mesh (K1 +2.2 cm); NAO.99Jb itself is
3-5 % above the 2021 gauges for K1.

**Which tide model, and where the excess comes from** (2026-09-30). Along
the base mesh's open boundary, inside the Enshu run, K1 was 1.06 times
NAO.99Jb. It did not move with the tidal potential (1.060), the bottom
friction (Cd minimum 0.0015 to 0.006: 1.062 to 1.058) or the bathymetry
(the new nodes from M7001 instead of the Cabinet Office grids, 1.4 % mean
change: 1.061). The three tide models then gave (rms vector difference at
the 7 gauges, cm; `--tide-model`, all with the tidal potential):

| forcing | mesh | M2 | S2 | K1 | O1 | Q1 |
|---|---|---|---|---|---|---|
| NAO.99Jb | Enshu | **1.7** | **0.6** | 2.3 | 1.4 | 0.6 |
| TPXO10-atlas-v2 | Enshu | 2.1 | 1.2 | 1.5 | 0.8 | 0.3 |
| FES2022b | Enshu | 2.2 | 1.2 | 1.5 | **0.7** | **0.2** |
| NAO.99Jb | base | 2.3 | 1.0 | 1.0 | 0.8 | 0.3 |
| TPXO10-atlas-v2 | base | 2.4 | 1.6 | 1.3 | 0.7 | 0.2 |
| FES2022b | base | **1.7** | 1.0 | **0.7** | **0.6** | **0.2** |

- **FES2022b is the best forcing on the base mesh**, whose open boundary is
  at the bay mouth: its constants there are close to the gauges
  (mean amplitude error M2 +0.1 cm, K1 +0.4 cm).
- **On the Enshu mesh every forcing comes out 3-5 % high at the bay
  mouth.** Forced by FES2022b or TPXO10 (which agree on the outer
  boundary), the run carries M2 0.362 m and K1 0.241 m along the base
  mesh's open boundary, where FES2022b has 0.351 m and 0.231 m (ratios
  1.03 and 1.045); the gauges then see M2 +1.2 cm and K1 +1.5 cm. NAO.99Jb
  looks better for M2 only because its outer-boundary M2 is lower. So the
  wide mesh itself amplifies the tide between the outer boundary and the
  bay mouth by a few per cent, for all constituents; not the bathymetry
  source, not the friction range tried.
- **Neither the open-boundary depth control nor self-attraction and loading
  explains it.** Keeping the outer boundary depths as built
  (`--no-obc-depth-control`; FVCOM would change them by up to 830 m) gives
  the same constants to the last digit: the elevation is prescribed there.
  Self-attraction and loading in the scalar approximation (`--sal-beta 0.1`,
  the owner's FVCOM reads a `SAL Beta` line) slows the tide and makes the
  excess larger (M2 at the gauges +2.4 cm instead of +1.2 cm on the Enshu
  mesh, +1.0 instead of +0.1 cm on the base mesh): leave it off.
- **Nor does a topographic (internal-tide) drag.** The owner's FVCOM takes
  `BOTTOM_ITD_COEFFICIENT` / `BOTTOM_ITD_MIN_DEPTH` (449
  `--itd-coefficient`): a linear bottom stress `C H |grad H|^2 u` below
  200 m. With C = 5e-4 and 2e-3 1/s the ratio to FES2022b at the bay mouth
  goes from 1.031 to 1.029 (M2) and 1.045 to 1.041 (K1): nearly nothing.
- Where it arises, by gauge (model / observed amplitude, FES2022b forcing):
  K1 is already 1.06 at Mera, on the open coast, and 1.05-1.08 in the bay;
  M2 is 0.98 at Mera but 1.02-1.07 in the bay. So the diurnal excess builds
  up outside, before the bay mouth, and the semidiurnal one inside the
  Uraga channel region, where the wide run's M2 field differs in shape
  from the base run's (the base run, forced at the bay mouth by FES2022b,
  is within 1 % in the bay).
- **A Flather boundary** (owner's FVCOM; 449 `--flather` with
  `--tide-model=tpxo10`). FVCOM's radiation boundaries (GWI, BKI, ORE)
  take no tide, so a Flather boundary was added: the tide file's optional
  `UnAmplitude` / `UnPhase` sections give the outward normal depth-mean
  velocity (TPXO transport over the model depth, on the outward normal
  FVCOM uses), and FVCOM solves `eta = eta_T + sqrt(D/g)(u_n - u_T)`
  implicitly with the boundary node's continuity, `u_n` being the flux that
  leaves the node's control volume. Two earlier versions failed and are
  recorded in FVCOM's history: an explicit one (u_n from the surrounding
  elements) and an implicit one that read the continuity value from ELF,
  which FVCOM leaves unchanged at open-boundary nodes (the interior flux is
  in `XFLUX_OBCN`). Found with **notebook 454**, an idealised channel
  (200 km x 20 km x 100 m, closed at one end, M2) against linear theory:
  clamped 0.35 cm rms, Flather with the theoretical velocity 0.15 cm rms,
  and Flather with `u_T = 0` stable.
- **Flather on the Tokyo Bay meshes** (TPXO10 elevation and transport,
  tidal potential, 200 days): stable on both. On the Enshu mesh the bay-mouth
  ratio to FES2022b moves from 1.03 (M2) and 1.045 (K1) to 0.984 and 1.027;
  at the gauges M2 improves (rms 2.1 to 1.7 cm, amplitude +1.2 to -0.7 cm)
  and K1 (1.5 to 1.2 cm), while S2 (1.2 to 1.5 cm) and O1 (0.8 to 1.3 cm)
  get worse. On the base mesh it changes little. So the boundary condition
  accounts for most of the wide mesh's few-per-cent excess; what remains
  is of the size of the differences between the tide models themselves.
  The best fit overall is still the base mesh forced by FES2022b (clamped).
  FES2022 offers no transports here, so FES + Flather needs another
  velocity source.

### 14.1 Which constituents to use

The open-boundary file may carry any number of constituents, and the tide
models offer many (NAO.99Jb 16, FES2022 34, TPXO10-atlas about 15). The
Japan Meteorological Agency predicts with 40 at a gauge. For the open
boundary of a regional model, **use the diurnal and semidiurnal
constituents; leave out the following**:

| constituents | why they are left out |
|---|---|
| **Sa, Ssa** (annual, semiannual) | Mostly not astronomical: seasonal heating, steric height, winds and currents. The ocean model or the reanalysis that gives the non-tidal sea level (JCOPE-T DA) carries them; adding them as tide counts them twice. |
| **Mm, Mf, Msf, Mtm, Msqm** (long period) | Small (a few cm at most) and close to equilibrium; a low-pass filtered reanalysis keeps them in its non-tidal part, so they would be counted twice. Their periods (9-32 days) also need long records to be told apart in a validation. |
| **S1** | Radiational: driven by the daily cycle of air pressure and wind, not by gravity. It belongs to the meteorological forcing. (S2 also has a radiational part, but a small one, and the tide models include it; keep S2.) |
| **M3, M4, M6, M8, MN4, MS4, MKS2, N4, S4** (shallow water, overtides) | Generated inside the bay by the model's own nonlinearity (advection, friction, the finite depth). At an open boundary on the shelf they are millimetres; a gauge in the inner bay needs them, which is why the JMA list is long, but the model makes them. Forcing them at the boundary adds a small, model-inconsistent signal. |

**Measured (2026-10-01, FES2022b, clamped, tidal potential, 185 days):**
the 9 minor diurnal/semidiurnal constituents (2N2 EPS2 J1 L2 LAMBDA2 MU2
NU2 R2 T2) on top of the major 8 lower the tide-series error at the gauges
(450 `tide_rms`: model minus the observation's own tide, diurnal and
faster) from 3.00 to 2.77 cm on the base mesh and from 3.33 to 3.09 cm on
the Enshu mesh; S2 improves most (rms 1.0 to 0.6 cm), since T2 and R2 no
longer fold into it. Use them.

**Worth adding** when a tide model has them: 2N2, MU2, NU2, L2, T2,
LAMBDA2, EPS2, R2 (semidiurnal) and J1, OO1 (diurnal). Each is a few
millimetres to 1 cm, but together they can improve a time series by 1-2 cm.
Mind two points: M1 is defined differently in different tables, so check
its convention before using it; and the tidal potential
(`tide_models.EQUILIBRIUM`) has coefficients only for the major eight, so a
run with `--equilibrium` needs the table extended (or the minor ones given a
zero potential) first.

## 15. Improving the tool with AI

Every new mesh is a new coastline, and it will find a case the tool has not
met: a pier at an unusual angle, a basin narrower than expected, a feature
that touches the frozen boundary. This tool was built, and is meant to be
kept up, with an AI coding assistant (Claude Code; a second model for review).
What worked:

1. **Show, do not describe.**
   - Give the assistant the recipe, the output directory and the figure you
     are looking at, and say what is wrong in the figure's own terms: "the
     pier at (392584, 3908790) is drawn at 11 deg to the quay".
   - Ask it to draw the mesh over the raw OSM (§7). Most defects are
     visible there and invisible to QA.
2. **Ask for the cause before the fix.** Ask the assistant to measure first:
   zoom on the offender, read `report.json`, compare runs. It should explain
   the mechanism in the code. A fix aimed at a symptom tends to move the
   problem somewhere else.
3. **Fix with a rule, not a case.** The fix should be a statement about the
   geometry ("a sector under 60 deg holds one element"), never a threshold
   tuned to one port. Thresholds and policies -- what is land, what is a wall,
   what is closed -- are the owner's decisions. The assistant should propose
   options with their consequences and let you choose.
4. **Test, then run everything.**
   - Each fix comes with a unit test that reproduces the case.
   - Then rebuild **every** recipe in `recipes/refine/` and compare "violations
     introduced" with the previous results. A fix that helps one mesh has
     broken another more than once.
   - Compare the coastlines too, not only QA. A rule added for one port can
     move another port's coastline and still pass every gate: the Odaiba
     step rule moved the Futtsu coast up to 47 m off OSM with 0 violations.
     Keep the previous output (move it aside, e.g. `refine_x.before`) and run

     ```bash
     python notebooks/435_compare_boundaries.py outputs/refine_x.before outputs/refine_x
     ```

     It prints the largest movement of the solid boundary. For each place
     that changed, it also prints the distance to the filtered OSM coast
     before and after. "identical" is the expected answer for recipes the
     change was not aimed at. A movement elsewhere is either intended, and
     then stated in the commit message, or a regression to fix.
   - For a change that matters, run the FVCOM smoke test too.
5. **Get a second opinion** on larger changes: an adversarial review by
   another model, given the design and the diff (for example
   `codex exec -m gpt-6-astra -s workspace-write ...`, with the repository's
   `AGENTS.md`). Reproduce each finding before acting on it.
6. **Write it down.** Record each decision and its measured effect in the
   design log (`docs/*_design.md`), and commit each verified change on its
   own. The next person -- or the next assistant session -- starts from that
   record.

The tool's own guards make this safe to iterate:

- the frozen contract is verified on the written file;
- the QA gate counts only what the patch introduced;
- a run refuses to overwrite an earlier result.

Keep the human decisions with the human: accepting a mesh, choosing
thresholds, and pushing to the shared repository.
