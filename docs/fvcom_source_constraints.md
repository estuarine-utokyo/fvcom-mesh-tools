# FVCOM mesh constraint set — derived from the local Fortran source

Verified 2026-07-03 against the production build. This document is the canonical
reference the QA checker (`fvcom_mesh_kickoff.md` §9) must encode; acceptance is
checked against the source, not against manual summaries.

## Authoritative source and version

- Production binary: TB-FVCOM `ersem/bin/fvcom_online`, built from
  `~/Github/FVCOM/src` = **FVCOM 5.1** (per `TB-FVCOM/CLAUDE.md:203-204`,
  `FVCOM/src/mod_utils.F`, `fvcom.F`).
- Production compile flags (`TB-FVCOM/ersem/make.inc`): **`-DGCN`** (so `*_gcn.F`
  files are active, not `*_gcy.F`), **CARTESIAN** (coords in metres),
  `-DMULTIPROCESSOR` + `-DMETIS_5`, `-DPROJ`, `-DTVD`, `-DDOUBLE_PRECISION`,
  `-DLIMITED_NO`; **no `-DPLBC`** (this keeps the critical stop in rule R4 active),
  no `-DWET_DRY`, no `-DTHIN_DAM`.
- `uk-fvcom` (FVCOM 5.0 + FABM) has byte-identical `tge.F` and equivalent readers
  (`mod_input.F` offsets ≈ −35 lines); every constraint below holds in both trees.
  Line numbers cited are FVCOM 5.1.

## Boundary classification (tge.F, TRIANGLE_GRID_EDGE)

Value meanings (`tge.F:73-86`):

| Marker | 0 | 1 | 2 | 3 |
|---|---|---|---|---|
| `ISONB(node)` | interior | solid boundary | open boundary | — |
| `ISBCE(elem)` | interior | solid boundary | open boundary | two solid boundary edges |
| `ISBC(edge)` | interior | boundary | — | — |

Assignment order matters:

1. `NBE(I,j)` = element sharing edge j (0 if none) — `tge.F:141-172`.
   `NBE(I,1)` is the neighbour opposite node 1 (across edge N2–N3).
2. Solid marking is purely topological: any zero-neighbour edge ⇒ `ISBCE(I)=1`
   and both edge nodes `ISONB=1` — `tge.F:204-217`.
3. OBC list then **overwrites**: `ISONB(I_OBC_N(:)) = 2` — `tge.F:540-542`.
   An open/land junction node (in the OBC list) therefore ends up `ISONB=2`.
4. `SUM(ISONB(NV(I,1:3))) == 4` ⇒ `ISBCE(I)=2` — `tge.F:558-581` (see R4).
5. Element with ≥2 zero-neighbour edges, not already open ⇒ `ISBCE(I)=3` —
   `tge.F:583-588`.
6. Open elements get porosity `EPOR=0` (`tge.F:857-860`) — their momentum is
   replaced by ghost/OBC treatment.

## Open-boundary mechanics (mod_obcs.F, mod_input.F)

- `_obc.dat` is read without any adjacency or order validation
  (`mod_input.F:4113-4234`); adjacency is reconstructed from mesh topology in
  `SETUP_OBC` (`mod_obcs.F:892-1197`), **not** from file order. Consecutive file
  lines need NOT be spatially adjacent.
- **Every OBC node must be mesh-adjacent to ≥1 other OBC node**, else
  `'NO ADJACENT NODE FOUND FOR BOUNDARY NODE'` → `PSTOP`
  (`mod_obcs.F:937-941`). OBC nodes must form connected chains.
- Junction handling: a segment endpoint is the OBC node with exactly one OBC
  neighbour (`NADJN_OBC==1`); its flux goes wholly to its single adjacent boundary
  cell (`mod_obcs.F:1033-1036`). Junctions are safe **provided the junction node is
  in the OBC list**.
- `NEXT_OBC` (first interior neighbour by max inward-normal dot product,
  `mod_obcs.F:1068-1095`) assumes each OBC node has a non-OBC neighbour; a
  one-node-wide neck leaves it 0 **silently** — QA should flag that geometry.
- OBC types 1–10 (`mod_obcs.F:30-39`): odd = elevation only, even = + nonlinear
  flux; {1,2} tidal (ASL), {3,4} clamped-zero, {5,6} gravity-wave radiation,
  {7,8} Blumberg–Khanta, {9,10} Orlanski. Reader enforces 1 ≤ type ≤ 10.

## R4 — the mixed-boundary-element rule (CONFIRMED, hard stop)

`tge.F:558-581`: for each element, let S = `SUM(ISONB(NV(I,1:3)))`.

- **S > 4** (node flags {2,2,1}=5 or {2,2,2}=6) → `PSTOP` with message
  "IT HAS EITHER TWO SIDES OF OPEN BOUNDARY OR ONE OPEN BOUNDARY AND ONE SOLID
  BOUNDARY". Active in production (`#if !defined(PLBC)` and PLBC is unset).
- **This is STRONGER than the prose rule** in the kickoff §7.2: an element whose
  open edge (two ISONB=2 endpoints) coexists with *any* third boundary node —
  solid (1) **or** open (2) — is fatal, even if that third node does not form a
  solid edge of that element.
- **S == 4** → `ISBCE=2`. Caveat: {2,2,0} (canonical open element) and {2,1,1}
  (one OBC node + two solid nodes, no open edge) **both** classify as open and get
  `EPOR=0`. FVCOM does not catch the {2,1,1} mis-tag; the QA checker must
  independently require every `ISBCE=2` element to contain an actual external edge
  whose two endpoints are consecutive OBC nodes.
- Why mixed cells cannot work downstream (GCN path): the solid-edge tangent
  `ALPHA(I)` is computed from "the one" `NBE==0` edge (`shape_coef_gcn.F:221-278`)
  — ill-defined with two zero edges; ghost-cell velocity branches solid-reflect vs
  open-copy on `ISBCE` (`ghostuv.F:52,82-84`); `bcond_gcn.F:381-453` zeroes
  `ISBCE==3` velocities and rotates `ISBCE==1`. (GCY variant differs at corner
  cells — irrelevant unless the build flag changes.)

## Hard startup checks (grid geometry error stops)

| # | Check | Location | Action |
|---|---|---|---|
| 1 | Element with no neighbours (isolated element) | `tge.F:191-198` | PSTOP |
| 2 | Mixed/double boundary element, S > 4 (R4) | `tge.F:558-581` | PSTOP |
| 3 | Non-manifold boundary at a node (ISONB ordering fails) | `tge.F:352-355` | FATAL_ERROR |
| 4 | Interior node element-fan does not close | `tge.F:317-327` | PSTOP |
| 5 | `Node/Cell Number =` headers missing in `_grd.dat` | `mod_input.F:3989-4007` | FATAL_ERROR |
| 6 | Connectivity row count ≠ NGL; max node index ≠ MGL | `mod_input.F:4081-4100` | FATAL_ERROR |
| 7 | CCW winding — tested on **element #1 only** (`IS_TRI_CW`, incl. collinearity tolerance) | `mod_input.F:4397-4407, 5249-5295` | FATAL_ERROR |
| 8 | Coordinate/depth/coriolis row counts ≠ MGL | `mod_input.F:4379-4501, 4630-4665` | FATAL_ERROR |
| 9 | OBC node id ∉ [1, MGL] or type ∉ [1, 10]; row count mismatch | `mod_input.F:4197-4224` | FATAL_ERROR |
| 10 | Sponge node id/range/count | `mod_input.F:4694-4789` | FATAL_ERROR |
| 11 | OBC node with no adjacent OBC node | `mod_obcs.F:937-941` | PSTOP |
| 12 | Degenerate solid-boundary cell (shape determinant < 1e-6) | `shape_coef_gcn.F:314-321` | PSTOP |

Warnings only (run continues): triangle area < 1e-6 m² and nodal control-volume
area < 1e-6 m² (`cell_area.F:144, 327`); depth/coriolis coordinates mismatching the
grid file (`mod_input.F:4528-4531, 4658-4661`).

## What FVCOM does NOT check (QA must add)

- Duplicate/coincident nodes — pass silently.
- Unreferenced (isolated) nodes — only the max-index test exists (#6); a mid-range
  orphan node is undetected.
- CW / zero / negative-area elements beyond element #1 — winding is checked once;
  `cell_area` uses `ABS`, so flipped elements are invisible.
- Node valence: >50 incident elements overruns a scratch array
  (`CELLS(MT,50)`, `tge.F:142`; the `'bad'` print at `tge.F:152` does not stop) —
  treat 50 as an undocumented hard cap; the project's quality cap is 8.
- All angle / skewness / aspect-ratio / adjacent-area-ratio metrics — the manual
  C1–C5 criteria are advisory to the executable; the mesh runs (badly) without them.

## File formats as parsed (mod_input.F, free-format, 1-indexed)

- `casename_grd.dat`: header contains `Node Number = <int>` and
  `Cell Number = <int>`; connectivity rows `CELL# N1 N2 N3`; then coordinate rows
  `NODE# X Y` (row order = node id). **On read, columns are swapped:
  `NVG = (N1, N3, N2)`** — the file must be CCW (SMS convention); FVCOM converts to
  its internal CW convention itself. Under CARTESIAN, X/Y are metres.
- `casename_dep.dat`: header `Node Number = <int>` (= MGL); rows `X Y H` in node
  order; H positive-down.
- `casename_obc.dat`: header `OBC Node Number = <int>` (0 allowed = no OBC); rows
  `OBCNODE# GLOBALNODE# TYPE` (sequence counter ignored; type 1–10).
- `casename_cor.dat`: header `Node Number = <int>`; rows `X Y COR`.
- `casename_spg.dat`: header `Sponge Node Number = <int>` (0 allowed); rows
  `GLOBALNODE# RADIUS DAMPING`.

## Minimum constraint set for the QA gate

FVCOM-enforced (mesh rejected at startup otherwise):

1. CCW node ordering in `_grd.dat` (check ALL elements, not just #1).
2. Connectivity row count = NGL; node indices dense with max = MGL.
3. No isolated elements (every element shares ≥1 edge).
4. **R4**: no element with `SUM(ISONB) > 4` — no open edge sharing a triangle with
   any other boundary node.
5. Every OBC node mesh-adjacent to ≥1 other OBC node (connected chains).
6. OBC ids ∈ [1, MGL], types ∈ [1, 10]; dep/cor/spg node counts = MGL.
7. Manifold boundary topology at every node (no bowtie/pinch nodes).

QA-added (model runs but should not):

8. No duplicate/coincident nodes; no unreferenced nodes.
9. No zero/negative-area elements anywhere; area ≥ 1e-6 m².
10. Valence ≤ 8 (project) — never approach the code's implicit 50 cap.
11. `ISBCE=2` elements must have a real open external edge (catch the {2,1,1}
    mis-classification).
12. No one-node-wide necks at the OBC (`NEXT_OBC=0` hazard).
13. Manual quality criteria C1–C5 (min angle ≥ 30°, max angle ≤ 130°,
    area ratio ≤ 0.5, valence ≤ 8) — advisory to FVCOM, mandatory for this project.
