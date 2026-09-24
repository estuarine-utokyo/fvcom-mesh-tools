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
| **fvcom-mesh-tools** (this repository, Apache-2.0) | https://github.com/estuarine-utokyo/fvcom-mesh-tools | public | everything |
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
`environment.yml`, conda-forge only. The three repositories above are
installed editable, without letting pip fetch anything:

```bash
# on a login node (compute nodes have no network)
mamba env create -n oceanmesh-bench -f environment.yml
mamba activate oceanmesh-bench
git clone https://github.com/estuarine-utokyo/oceanmesh.git ~/Github/oceanmesh
git clone https://github.com/estuarine-utokyo/xcoast.git ~/Github/xcoast
(cd ~/Github/oceanmesh && pip install -e . --no-deps --no-build-isolation)  # compiles C++
(cd ~/Github/xcoast && pip install -e . --no-deps --no-build-isolation)
pip install -e . --no-deps --no-build-isolation                          # this repository
```

- **Compiling oceanmesh:** it compiles C++ (CGAL). On OCTOPUS,
  `qsub jobs/octopus/build_env.sh` does the compile step as a batch job.
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
  in `jobs/octopus/`.
- **Paths in the job scripts are this account's.** The accounting group is
  `G16445`; the FVCOM binary is
  `/octfs/work/G16445/v61021/Github/FVCOM/src/fvcom`; the scratch directory is
  under `/octfs/work/G16445/v61021/scratch`. Another user edits them.
- **Monitoring:** start a monitor after every `qsub`, and read the log when
  the job ends.

### 2.5 What the base mesh must be

Two limits of the generator, checked before any meshing and refused with a
reason:

- **Coordinates in UTM zone 54N metres (EPSG:32654).** The generator
  projects OSM, the regions and the depth products into that frame. A base
  in another CRS must be reprojected first.
- **Exactly one open-boundary arc.** The boundary lists of the refined mesh
  are rebuilt around a single arc.

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
| `refine` | yes | one or more regions |
| `refine[].geometry` | yes | `circle: {center: [lon, lat], radius_m}`, a `bbox`, a GeoJSON `Polygon`, or `{file: area.geojson, where: {...}, buffer_m: 25}` |
| `refine[].target_h_m` | yes | the target element size, m |
| `refine[].priority` | no | which region wins where two overlap |

Unknown keys are errors, and relative paths resolve against the recipe's own
directory. Before meshing, the run checks each region (depth, dryness, time
step) and prints what it found.

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

A port-sized region takes 5-10 minutes on one core.

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
| `<base>_<name>_walls.json` | the node pairs a wall duplicated on purpose, bound to the `.14` by its SHA-256; `fmesh-mesh-qa` reads it, so the delivered mesh gets the same verdict when checked again |
| `ACCEPTED` | written last, only when every gate passed (§3) |
| `shoreline_filtered.shp` | the OSM land after the width filter |
| `fill_constraints.npz`, `walls_stages.npz` | the rim and walls the mesher was given, for inspection |
| `final_mesh_*.png` | the figures (§7) |

The run log ends with lines like these:

```
[lr]     seed 0: QA 20/22, 0 introduced by the patch
[lr] accepted seed 0
[lr] achieved in my_port: median cell 28.5 m against a 30 m target, 100.0 % of the water within 1.25x
```

- **`0 introduced`** is the goal.
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

The generator lives in `notebooks/`, not in the package, because it imports
the GPL-3.0 oceanmesh; the Apache-2.0 package may not (`CLAUDE.md`). Its
parts -- selection, rim, walls, stitching, verification, QA -- are in the
package and have unit tests.

---

## 11. Known limits

- Tested on the Futtsu and Kimitsu coast of Tokyo Bay, one base mesh: a port,
  a coastal region, an offshore fishery, and the default branch. A new place
  will find new cases (§12).
- A structure hugging the coast within 0.4 element, thinner than half an
  element, is not represented.
- OSM `man_made` lines (breakwaters mapped only as lines) are not used.
- The result depends on the DistMesh seed; the best of several is kept.
- The M2 staging in job 421 is specific to the Tokyo Bay goto2023 base.
- The base must be in EPSG:32654 with one open-boundary arc (§2.5).
- With `hires.coastline: preserve` the coastline is kept exactly, including
  corners sharper than 60 deg. These are reported (`acute_corners_kept` in
  `report.json`), not cut. A node left in one element at such a corner is
  opened by the repair, and the QA gate says if that failed.

---

## 12. Improving the tool with AI

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
