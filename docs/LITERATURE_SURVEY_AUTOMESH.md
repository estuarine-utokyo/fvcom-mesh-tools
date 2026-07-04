# Automated Coastal Ocean Mesh Generation: Literature Survey (2020-2026)

Compiled 2026-07-05 (web research agent; load-bearing citations
Crossref-verified). Context: fully automated FVCOM mesh pipeline for
Tokyo Bay (~300 m), reference = hand-finished OceanMesh2D+SMS mesh
(goto2023).

## 1. SCHISM ecosystem — the "AI-assisted" papers

Headline: the recent SCHISM mesh-automation papers are NOT machine
learning; "data-driven" means DEM-driven deterministic automation.

### 1.1 RiverMapper / RiverMeshTools (the key paper)

Ye, F., Cui, L., Zhang, Y.J., Wang, Z., Moghimi, S., Myers, E.,
Seroka, G., Zundel, A., Mani, S., Kelley, J.G.W. (2023). "A parallel
Python-based tool for meshing watershed rivers at continental scale."
*Environmental Modelling & Software* 166, 105731.
DOI: 10.1016/j.envsoft.2023.105731.
Code: github.com/schism-dev/RiverMeshTools — **Apache-2.0**, active.

- pyDEM: depression filling, D8 flow accumulation -> thalweg network.
- RiverMapper: bank detection by perpendicular elevation scans,
  longitudinal river arcs (bank + interior), width-scaled arc count
  and spacing, confluence snapping (~25% of local width). Output: SMS
  .map or OCSMesh shapefiles. ~10 min for ~30,000 rivers.
- **Pseudo-channels** (`i_pseudo_channel`): channels below resolvable
  width become FIXED-WIDTH templates with a FIXED number of
  cross-channel element rows ("artificial dredging in the mesh") —
  the published, operational answer to sub-resolution channels:
  widen to a prescribed template, let depth carry conveyance.
- Caveat for FVCOM: SCHISM tolerates elongated/skewed elements;
  FVCOM's 30-degree gate requires capping aspect ratio (2-3 rows at
  ~300/150 m).

### 1.2 OCSMesh (NOAA)

Mani, S., Calzada, J.R., Moghimi, S., Zhang, Y.J., Myers, E.,
Pe'eri, S. (2021). NOAA Tech. Memo. NOS CS 47.
DOI: 10.25923/csba-m072. Code: github.com/noaa-ocs-modeling/OCSMesh —
**CC0-1.0 (public domain)**, very active (v2.x, gmsh default engine).

- Channels: `utils.get_polygon_channels` — erode domain polygon by
  -width/2, re-buffer, difference; residuals refined. Also ingests
  RiverMapper arcs (`add_feature`).
- Open boundaries: `Mesh.boundaries.auto_generate(threshold=0)` —
  classification only; line placement is manual/domain-given.
- Cleanup: deterministic targeted passes (isolates, pinched nodes,
  folded boundary elements, slivers); smoothing delegated to engine.
- CC0 = liftable verbatim into our Apache-2.0 package.

### 1.3 End-to-end automation (2025-2026)

- Cassalho, F., Mani, S., Moghimi, S., Ye, F., Zhang, Y.J. (2026).
  "OCSMesh and an automated creek-to-ocean mesh generation workflow."
  *Ocean Modelling* 203, 102774 — pyDEM -> RiverMapper -> OCSMesh v2,
  no SMS; reproduced the hand-built STOFS-3D-Atlantic mesh in <6 h
  (NOAA webinar 2025-08-05); drives automated STOFS-3D-Pacific.
- Calzada, J.R., Zhang, Y.J., Ye, F., Cui, L. (2025). "Geomesh."
  *EMS* 192, 106587. JIGSAW-based, DEM-only input; no public code.
- Operational status: STOFS-3D-Atlantic v3.1 still "minor manual
  edits" in SMS; the fully automatic path is being phased in.

## 2. Other generators (summary)

| Tool | Sub-res channels | OBC | Cleanup | License |
|---|---|---|---|---|
| OceanMesh2D v6 | thalweg size fn; v6 high-fidelity shoreline constraints (pfix/egfix, MATLAB CDT) | auto-classify | msh.clean modes, in-loop | GPL-3.0 |
| oceanmesh (Py) | none | none | traversable + delete_boundary_faces + Laplacian | GPL-3.0 |
| JIGSAW-GEO | none | none | hill-climb; RAD2=1.05 -> min angle ~28.4 deg | non-OSI |
| MeshKernelPy | depth-based refine | none | richest ops: fractional-area sliver deletion, valency flips, orthogonalization | LGPL core/MIT wrapper |
| ADMESH+ 2024 | **automatic medial-axis width classification; 1D/2D split; constrained centerlines in DistMesh** | n/r | DistMesh q~0.97 | GPL-3.0 (MATLAB) |

Key existence proof: Pringle et al. 2021 (*GMD* 14:1125) built global
ADCIRC meshes fully scripted with OceanMesh2D v3 — "No
post-processing hand edits of any mesh were necessary for numerical
stability" (ADCIRC tolerance, not FVCOM's gates).

ADMESH+ (Kang & Kubatko 2024, *GMD* 17:1603-1625,
DOI: 10.5194/gmd-17-1603-2024): water-mask medial axis with width =
sum of two bank distances; narrow regions -> 1D centerlines FIXED as
constraints in 2D DistMesh — independently validates our
water-skeleton seeding (PoC #91); reimplement (GPL, MATLAB), do not
port.

Folklore rule in print (Remacle & Lambrechts, *CAD* 2018): no island
and no channel smaller than the local mesh size may remain in the
domain.

## 3. ML / RL / LLM meshing (2023-2026)

Nothing production-usable: RL triangle improvers (Thacher/Persson
2025, arXiv:2504.03610; GNNRL-Smoothing 2025) have no code and no
hard bounds; GMSNet (2024) is speed-only; UM2N (NeurIPS 2024) is
r-adaptivity, not generation; LLM-CFD agents orchestrate classical
meshers without quality gating. Survey arXiv:2512.23719 confirms: no
coastal applications, no hard-angle-bound repair. Deterministic
gated passes remain the right architecture; hard-gated auto-finishing
is open, publishable territory.

## 4. FVCOM specifically

- **No FVCOM-targeted automated mesh generator has ever been
  published** (incl. Chen lab). Manual v3.1.6: "no automatic mesh
  generator is supplied"; SMS walkthrough with hand node-dragging.
- The C1-C5 gates ARE the FVCOM manual's SMS dialog settings (Ch. 20
  pp. 374-376): min angle 30.0, max angle 130.0, area change 0.5,
  max 8 elements/node; PLUS the OBC-triangle interior edge normal to
  the open boundary (documented anti-noise requirement — validates
  our perpendicularity fixer).
- NOAA operational FVCOM meshes (LMHOFS, NGOFS2): hand-built in SMS.
- No published "widen the channel, conserve the cross-section"
  operator for FVCOM; nearest print: subgrid-channel literature
  (Neal et al. 2012 *WRR* 48:W11506) and RiverMapper's pseudo-
  channels. Publication opportunity for our pipeline.

## 5. Recommendations

Sub-resolution channels/harbors:
1. Pseudo-channel templates (RiverMapper, Apache-2.0) — fixed-width,
   fixed-row channel meshes; cap rows for the 30-degree gate.
2. Medial-axis width classification (ADMESH+ algorithm,
   reimplemented) — decides WHICH features get templates.
3. OCSMesh channel detection + feature arcs (CC0, lift verbatim).

Auto-finishing:
1. OceanMesh2D bounded cleanup recipe incl. in-generation cleanup
   cadence and mean-quality-minus-3-sigma termination (reimplement).
2. MeshKernelPy local operators (fractional-area sliver deletion,
   valency flips) as an optional backend.
3. No ML adoption now; our gated deterministic operator suite is the
   contribution.

Open-boundary note: no surveyed tool automates OBC line placement
(classification only) — our geometric OBC construction has no
published competitor.

Strategic: SCHISM solved channel automation by RELAXING quality;
FVCOM cannot. Combining SCHISM-ecosystem channel automation
(Apache/CC0) + OM2D-style bounded cleanup + MeshKernel operators
under FVCOM's strict gates is unpublished territory.
