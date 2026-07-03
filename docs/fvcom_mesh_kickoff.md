# fvcom-mesh-tools — AI-Assisted Mesh Generation Kickoff

> Purpose: give the AI agent the **full picture up front** (scope, data, fixed
> settings, reconnaissance, constraints, workflow, done-criteria). Register the
> concise completion condition separately with `/goal` (§11), and use `/loop` with
> the mesh-QA validator (§9) to iterate until acceptance criteria pass.
> Keep this file in the repo (or as a skill under `.claude/skills/`) so it loads as
> durable context rather than being re-typed each session.
>
> Language policy: chat and progress reports are in Japanese; **all code, tracked
> documents, and figure labels are in English.**

---

## 1. Goal

Build an **AI-assisted toolchain that produces FVCOM unstructured triangular meshes
from prose intent**, and that can **revise an existing FVCOM mesh** from prose.

Two intertwined objectives:

- **A general-purpose, config-driven meshing toolkit** (`oceanmesh` +
  `oceanmesh-tools` + `xcoast` + `fvcom-mesh-tools`).
- **An interactive construction workflow**: the user directs meshing in prose while
  looking at figures, and the agent edits the mesh recipe and regenerates.

Working philosophy is **phased**:

- **Phase 0 — reconnaissance.** Survey the existing toolkit and reference scripts;
  plan reuse/port before writing new code (§6).
- **Phase 1 — get a mesh that FVCOM runs on.** Produce a valid prototype mesh for
  Tokyo Bay. Correctness and runnability first (§7).
- **Phase 2 — interactive refinement.** Improve the mesh through chat, and **promote
  each ad-hoc fix into a reusable, backward-compatible tool feature**, generalizing
  the toolkit through the refinement process.

## 2. Scope, tools, and modification policy

Primary tools:

- `oceanmesh` (Python) — sizing functions and DistMesh-based generation; the active
  generator for this Python pre/post workflow.
- `oceanmesh-tools`, `xcoast` (land/coast rendering + OSM tools), `fvcom-mesh-tools`
  (FVCOM I/O + QA).
- `oceanmesh2d` (MATLAB) — higher-feature reference implementation.

Modification policy:

- `oceanmesh`, `oceanmesh-tools`, `xcoast`, `fvcom-mesh-tools` may be improved freely
  **provided backward compatibility is preserved** (public signatures, recipe schema,
  output formats). Document every change.
- Features present in `oceanmesh2d` but missing in `oceanmesh` may be ported in.
- Prefer additive, config-driven changes over silent behavioral changes. When a
  Phase-2 interactive fix proves useful, generalize it into the toolkit.

## 3. Operating principles

1. **Prose → recipe → mesh → figures → QA → iterate.** The agent emits an editable,
   version-controlled **mesh recipe** (`recipe.yaml`) capturing every decision (CRS,
   datum, source precedence, region polygons, sizing functions, boundary treatment,
   quality targets). The recipe — not ad-hoc reasoning — drives generation, so runs
   are reproducible and auditable.
2. **Deterministic generation.** Same recipe + inputs → same mesh (fix seeds; record
   tool versions).
3. **Close the loop.** After generation, render figures (§10) and run the QA gate (§9),
   report a pass/fail table in Japanese, and iterate until all criteria pass or an
   unavoidable violation is surfaced for a human decision.
4. **Provenance.** With each mesh, write a record: datasets + versions, CRS, vertical
   datum, precedence applied, tool versions, and the exact `recipe.yaml` used.
5. **Approval gates.** Ask for explicit approval before destructive/irreversible
   actions (overwriting a production mesh, `git push`, `gh pr create/merge`).
6. **State assumptions inline.** On ambiguity, choose a sensible default from
   ocean-modeling knowledge, record it in the recipe, and flag it — do not stall.

## 4. Input data

Deliverable: **audit all inputs first** into `docs/DATA_INVENTORY.md`. Per dataset
record CRS/EPSG, vertical datum, resolution, coverage, source/authority, license,
known issues, and precedence.

Bathymetry — `$DATA_DIR/bathymetry/`
- `M7001`, `tokyo_bay` — primary inner-bay sources (document how they overlap and which
  wins where).
- `mesh500` — baseline for the Japanese coast **outside** Tokyo Bay.
- `Futtsu_JFA_2023` — high-detail tidal-flat / shallows bathymetry **north of Cape
  Futtsu**; highest precedence within its footprint.
- `GEBCO` — low accuracy near Japan; **fallback only**.

Coastline / rivers
- `$DATA_DIR/OSM/` — key source for coastline and river-bank geometry.
- `$DATA_DIR/coastline/` — use `tokyo_bay/` as the base (from Japanese government data).
  `GSHHS` is low accuracy near Japan; **fallback only**.

Precedence (agent applies, and records the resolved precedence):
- Bathymetry (high → low, within footprint): `Futtsu_JFA_2023` → `M7001` / `tokyo_bay`
  → `mesh500` → `GEBCO`.
- Coastline (high → low): `coastline/tokyo_bay` → `OSM` → `GSHHS`.

## 5. Fixed project settings

- **Horizontal CRS: UTM Zone 54N.** Reproject all inputs to it. Spherical is out of
  scope for now.
- **Vertical datum: T.P. (Tokyo Peil).** If a source's datum is unclear, treat mean sea
  level as 0 and record the assumption.
- **Minimum depth: clip to 2 m.** Depths shallower than 2 m are set to 2 m.
- **No wetting/drying this time.** The domain is bounded by the coastline and river
  boundaries; land inundation (e.g. tsunami runup) is not modeled.

## 6. Phase 0 — toolkit & reference reconnaissance (do this first)

Before writing new code, survey the current state and prefer reuse/port over reinvention.
Deliverable: `docs/TOOLCHAIN_SURVEY.md`.

1. **Inventory the four toolkits** (`oceanmesh`, `oceanmesh-tools`, `xcoast`,
   `fvcom-mesh-tools`): existing modules, public APIs, recipe/config schema, FVCOM I/O,
   QA, and figure utilities. Record what already exists vs. what must be added.
2. **Study the reference MATLAB scripts** in `Github/OceanMesh2D/Tokyo_Bay/` (local to
   the working machine). Extract the intended Tokyo Bay workflow: domain/bbox, shoreline
   source, edge-length (sizing) functions, boundary handling, and bathymetry steps. Map
   each step to its Python `oceanmesh` equivalent for reuse.
3. **Anchor on the OceanMesh2D pipeline** the reference scripts follow:
   `geodata` (geospatial input) → `edgefx` (mesh size functions) → `meshgen` (generate)
   → `msh` (store / read / write / inspect / visualize + model inputs). The recipe (§3)
   mirrors this decomposition; the Python `oceanmesh` API mirrors it too.
4. **Reuse proven example patterns** rather than ad-hoc logic: the OceanMesh2D *JBAY*
   example (high-resolution coastal mesh with CFL-limiting) for shallow/flat resolution;
   the *GBAY* example's polyline/thalweg channel size function for narrow rivers/canals
   (§7.4).
5. **Log gaps** between the MATLAB reference and the Python toolkit as Phase-2
   generalization candidates (§1).
6. **Derive the authoritative mesh constraint set from the local FVCOM source**, not
   from third-party summaries. Inspect how boundary nodes/cells are classified and how
   the OBC is applied — e.g. the node/cell boundary markers `ISONB` / `ISBCE` and the
   element-neighbour array `NBE`, set in `TRIANGLE_GRID_EDGE` (tge.F) / `SETUP_DOMAIN`,
   plus the open-boundary setup in `mod_obcs.F` (names may vary by FVCOM version).
   Confirm in particular the §7.2 rule that **no element carries both an open-boundary
   and a solid-boundary edge**, then encode every verified rule into the QA checker (§9)
   so acceptance is checked against the source, not against assumptions.

## 7. Phase-1 constraints (FVCOM-runnable mesh)

### 7.1 Mesh geometry quality (acceptance thresholds)
The mesh must be a valid (non-overlapping, no degenerate/sliver) Delaunay
triangulation satisfying:

| Metric                                   | Requirement |
|------------------------------------------|-------------|
| Minimum interior angle                   | ≥ 30.0°     |
| Maximum interior angle                   | ≤ 130.0°    |
| Element area change (adjacent gradation) | ≤ 0.5       |
| Connecting elements per node             | ≤ 8         |

Elements should be near-equilateral in general (open-boundary elements excepted).

### 7.2 Open boundary (verify against FVCOM source — see §6, item 6)
- Open-boundary elements should be **as close to right-angled as possible, with one
  interior edge normal to the open boundary** (FVCOM manual §6, Fig. 6.2; SMS guide).
  FVCOM runs without this, but the radiation/clamped OBC otherwise produces increased
  numerical noise from wave reflection.
- **No single triangle may straddle both boundary types**: an element must not have an
  open-boundary edge and a solid (land) boundary edge at the same time. At an
  open/land junction the corner node is shared, but the open edge and the land edge
  belong to different elements. (Confirm the exact rule in the source; it also avoids
  the corner element-overlap issue noted in the FVCOM grid docs.)
- FVCOM does not distinguish open vs land boundary in the grid file itself; the
  `_obc.dat` node list designates which boundary nodes are open (absent → treated as
  land). Emit the ordered open-boundary node list accordingly.
- Place the open boundary in relatively deep, smooth water; straight or a smooth arc;
  roughly uniform along-boundary resolution.

### 7.3 Boundary conforming
- Coastline and river-bank polygons are **forced boundaries**.
- **Resample/simplify** raw coastline (e.g. OSM) to the local target element size before
  use, to avoid over-dense boundaries and micro-elements.

### 7.4 Narrow channels (rivers / canals)
- Configurable **minimum element size for narrow channels** and a configurable
  **inclusion threshold**: important channels are widened to the specified size and
  meshed; unimportant small channels are excluded (state the decision criterion).
- Ensure flow can pass: at least **1 element across (preferably 2–3)** along the channel;
  run a **wet-domain connectivity check** to detect and fix pinch-offs / blocked channels.
- Reuse the thalweg/polyline channel size function pattern (§6.4).

### 7.5 Bathymetry (depth field is independent of mesh geometry)
- Interpolate node depths from the merged, datum-harmonized (T.P.) sources.
- Clip to minimum depth 2 m (§5).
- **rx0 / max-slope (≤ 0.1) smoothing is a SEPARATE depth-field post-process, decoupled
  from mesh generation.** It is applied to the depth field, not by moving nodes, and is
  **optional** — often unnecessary given the sigma-z + Delft internal-PGF correction.
  Do not apply it in Phase 1 unless requested; keep it available for Phase 2.
- Any depth-dependent (CFL / wavelength) sizing in §8 reads the **real (unsmoothed)**
  bathymetry.

## 8. Sizing / resolution model (prose → sizing functions)

Translate prose ("raise the Banzu tidal-flat resolution to ~50 m") into region
polygons + sizing-function overrides, combining:

- Distance-to-shore sizing (finer near coast).
- Feature size (medial-axis) sizing to resolve narrow features.
- Wavelength / CFL sizing from **real** local depth (report the implied external-mode
  Δt so cost is visible).
- Bathymetric-gradient sizing (optional).
- User region overrides from prose (via §10 reference grid or named regions), blended
  into the background size field under the gradation limit (§7.1 area-change ≤ 0.5).

Ambiguous spatial terms are resolved with ocean-modeling judgment and recorded in the
recipe.

## 9. QA gate / validator (use with `/loop`)

A single command (e.g. `fvcom-mesh-tools qa <mesh>`) that:

1. Checks every §7.1 threshold; reports a Japanese pass/fail table with worst offenders.
2. Verifies Delaunay validity, no overlaps/slivers, boundary conformity.
3. Verifies open-boundary perpendicularity, node ordering, and that **no element has
   both an open-boundary and a solid-boundary edge** (§7.2).
4. Runs wet-domain connectivity / channel-flow checks (§7.4).
5. Reports min-depth compliance and the implied stable Δt.
6. Reports node/element counts.

Wire this as the `/loop` validator so the agent iterates until it passes.

## 10. Visualization & reference grid (instruction-addressing)

Every mesh iteration produces:

- **Whole-domain figure**: the full mesh with **land rendered by `xcoast`**, plus a
  **rectangular reference grid overlay** defined in UTM 54N (true squares; default
  **5 km × 5 km**, configurable). Columns labeled A, B, C… (west→east), rows 1, 2, 3…
  (north→south); a cell is e.g. "C4". Named-region aliases (e.g. `Banzu`, `Futtsu`,
  `Arakawa-mouth`) are also drawn and are usable interchangeably with grid cells.
- **Regional zoom panels** per grid cell / named region.

The agent parses region references from chat ("C4", "C4–D5", "Banzu") back into UTM
polygons for the recipe, so instructions like *"raise resolution in C4 to 50 m"* are
unambiguous. The reference-grid definition (cell size, labels, aliases) lives in the
recipe.

## 11. `/goal` registration (concise completion condition)

> `/goal` Produce/revise the FVCOM mesh for `<case>` per fvcom_mesh_kickoff.md so that
> all §9 QA checks pass (all §7.1 thresholds met, open boundary perpendicular and
> ordered, no blocked channels, depths interpolated to T.P. and clipped to 2 m),
> emitting the §12 deliverables plus recipe + provenance and the §10 figures. Not done
> until the QA report is all-pass or remaining violations are surfaced for my decision.

## 12. Outputs / deliverables

- FVCOM grid inputs: `<case>_grd.dat`, `<case>_dep.dat`, `<case>_obc.dat`
  (+ hooks for `<case>_cor.dat`, `<case>_spg.dat`; `<case>_sigma.dat` is out of scope
  for the mesh tool unless requested).
- Interoperability export (e.g. SMS `.2dm`).
- Figures from §10 (whole-domain + reference grid, and regional zooms).
- `recipe.yaml`, the QA report, the provenance record, `docs/DATA_INVENTORY.md`, and
  `docs/TOOLCHAIN_SURVEY.md`.

## 13. Definition of Done (Phase 1)

- [ ] `TOOLCHAIN_SURVEY.md` and `DATA_INVENTORY.md` written; CRS = UTM 54N, datum = T.P.
- [ ] Mesh constraint set verified against the local FVCOM source (§6, item 6) and encoded in the QA checker.
- [ ] Mesh passes all §7.1 thresholds (QA table all-pass).
- [ ] Open boundary: near-right-angled elements normal to the boundary, in smooth/deep water, ordered node list emitted; **no element carries both an open-boundary and a solid-boundary edge**.
- [ ] Coastline/river boundaries conformed; raw coastline resampled to target size.
- [ ] Narrow channels: important ones meshed at min size and flow-connected; unimportant ones excluded per stated criterion.
- [ ] Depths interpolated to T.P. and clipped to 2 m; no rx0 smoothing applied in Phase 1.
- [ ] §12 deliverables + `recipe.yaml` + provenance + §10 figures produced; run reproducible.

## Appendix — Example prose instructions

Phase-2 interactive intents the toolkit must eventually handle:

- Extend the Arakawa and Edogawa rivers ~2 km further upstream in the existing mesh.
- Raise resolution around the Banzu tidal flats (Kisarazu) to ~50 m — or "raise C4 to 50 m".
- Reposition the open boundary along the ~50 m isobath near the bay mouth and enforce
  perpendicular boundary elements.
- Align mesh edges to the rim of a dredged pit so a sigma-z staircase conforms.
- Apply rx0 ≤ 0.1 depth smoothing over a specified region (leaving mesh geometry unchanged).
