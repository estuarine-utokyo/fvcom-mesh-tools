# Channel resolution policy: design proposal (not implemented)

Status: PROPOSAL, 2026-09-20. Drafted by gpt-6-astra (Codex CLI) at the
owner's request, reviewed by Claude; nothing in it is implemented.

Question: given a target grid size near the coast, which rivers, port
canals and navigation channels should be resolved (refined or widened
beyond reality) and which should be filled as land? The aim is a small
set of options, with numeric arguments, that makes the decision nearly
automatic and reproducible.

Read together with `DESIGN_HISTORY.md` (why earlier approaches failed),
`STAGE2_DESIGN.md` (detector/planner/executor structure) and
`LITERATURE_SURVEY_AUTOMESH.md` (prior art). Formulas, defaults and
tolerances below are engineering proposals that still need site and
model verification.

---

**Recommendation.** Add a versioned `waterways` policy block to the existing recipe. Decide **which water bodies and connections must survive**, then select `keep`, `refine`, `widen`, or `fill` subject to numerical and physical acceptance. Preserve `review` as an explicit outcome when the available data cannot justify automation. A narrow channel’s width alone should not decide its fate.

**Repository findings and one measurement discrepancy**

- `src/fvcom_mesh_tools/waterways.py` already provides detection, keep/close classification, arc widening, bridge handling, achieved-width checks, and unresolved-water normalization. Extend these mechanisms.
- Current classification uses resolution-relative proxies: `min_basin_cells=6`, `min_resolve_width_frac=0.2`, canal extent/width tests, and through/anchor connectivity. These are useful geometric heuristics, but they do not measure hydraulic importance.
- `channel_arcs.py` provides constrained, axis-based carving and barrier protection. Retain it; evaluate interacting corridors together so two individually acceptable widenings cannot consume their shared wall.
- `channel_policy.py` and `notebooks/331_finish2.py` can delete small basins during finishing; 331 uses `min_basin_elements=25`, strict one-wide detection, and at most eight pruning rounds. A new protected-feature decision must survive this stage.
- `notebooks/325_sample_repro.py` uses conservative settings distinct from library defaults: `SR_WIDEN_FACTOR=0.875`, `SR_ATTAIN_BAR=1.5`, normalization on, forced rows/refinement off. Policy scale is `1.2*H0`, hence 348 m at default `H0=290`.
- `recipes/*.yaml` and `src/fvcom_mesh_tools/cli/pipeline.py` supply the extension point, but the current recipe CLI is not a configuration-complete replacement for the certified notebook chain.
- The four JSON edits distinguish data correction (`edit_001`), domain-edge correction (`edit_003`), purpose-driven harbour removal (`edit_004`), and channel representation (`edit_005`). Preserve that distinction.
- **The supplied “every gate passes” summary conflicts with the stored [job 115210 log](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/logs/385_no_edits.115210.log:624).** It records 21/21 mesh-QA gates passing and zero critical connectivity failures, but two confirmed one-wide cells at OW05 and notebook 364 `GATE FAIL`: 21 stray elements and 71 wall crossings. Its job script tolerates checker failures with `|| true`.
- Recorded no-edits results are 3,406 nodes, 5,863 elements, implied dt 16.35 s; the comparison’s certified mesh has 3,393 nodes, 5,849 elements, implied dt 16.32 s. These are existing-log measurements, not new runs. Any later adjudication of those flags remains unverified.

**Decision criteria and available evidence**

| Category | Criterion and decision use | Available now; missing evidence |
|---|---|---|
| Numerical | Natural width profile \(W(s)\), target spacing \(h(s)\), achieved cross-channel rows. Keep naturally resolved water; otherwise compare refinement and widening. | Land polygons plus centerlines provide geometric estimates; OSM centerlines alone do not establish width or navigability. |
| Numerical | External-mode timestep, advection, cost, and grading halo. Reject refinement that violates the declared timestep budget. | DEM plus candidate mesh support the repository’s implied-dt calculation; velocity, water-level envelope, solver safety factor, and actual runtime require model configuration/results. |
| Numerical | No bank-to-bank one-cell choke on retained routes; acceptable angle, area change, valence, boundary topology. | Actual generated mesh, QA, and notebook 346; requested width is insufficient evidence. |
| Numerical | Remaining wall thickness between **all** neighboring carved corridors; no accidental basin merger. | Land polygons and corridor geometry; narrow wall accuracy, culverts, gates, and bridge status may require survey/as-built information. |
| Physical | Basin area \(A_b\), volume \(V_b\), tidal prism \(P\), and exchange. A tiny entrance can serve a large, important basin. | Area is geometric; volume needs adequate submerged DEM coverage. Tidal range/phase and actual exchange need forcing, observations, or simulation. |
| Physical | River freshwater and constituent loads; preserve the route or explicitly relocate the source. | River list identifies known sources; discharge \(Q(t)\), concentrations, temperature, and salinity are external unless actually supplied in that list. |
| Physical | Through-path connectivity, alternate-route length, and flow partition. Preserve meaningful shortcuts even if another graph path exists. | Geometry supplies connectivity/detour proxies; hydraulic significance needs depths, forcing, resistance, and model results. |
| Physical | Stratification, salt intrusion, two-layer exchange, sill control. Prefer true geometry and adequate vertical resolution. | DEM may identify sills; density profiles, seasonal discharge, turbulence settings, and observations are required. |
| Purpose | Hydrodynamics only: omit a peripheral basin only when its influence on the modeled receiving water is acceptable. | User supplies the receiving-water objectives and error tolerance. “No port water-quality objective” does not mean “unimportant.” |
| Purpose | Port water quality: retain the basin, relevant entrances, storage, sources, and exchange paths. | User identifies receptors/processes; loading and reaction data are additional inputs. |
| Purpose | Dye/tracer transport: retain source-to-receptor paths, travel-time controls, storage, and flow splits. | Source/receptor locations and release scenario are required; a connected graph alone does not ensure correct transport. |

Useful screening quantities are \(P\approx A_b\Delta\eta_b\) and \(\tau\approx V_b/(Q_r+\epsilon P/T)\), where \(T\) is the tidal period and \(\epsilon\) represents effective replacement after return flow. Neither basin tidal range nor \(\epsilon\) is known from OSM. These are ranking estimates, not flushing predictions; avoid double-counting river flow in a measured exchange estimate.

**What widening changes, and compensation**

At unchanged depth and roughness, widening generally increases cross-sectional area, storage, and conveyance. It can increase basin tidal amplitude/prism, alter flow partition and frictional losses, and change residence time. The residence-time response is not necessarily monotonic: exchange, storage, recirculation, and return flow all change.

Use surveyed or credible DEM-derived sections at stations along the original axis. Never interpret elevations on newly carved land as channel bathymetry. Let \(A_0,R_0,n_0\) be original section area, hydraulic radius, and Manning roughness at a specified reference level; primed quantities describe the achieved mesh geometry.

- **Preserve area:** choose the engineered section so \(A'=A_0\). For rectangular sections, \(D'=D_0W_0/W'\). This preserves section volume per unit length and mean velocity at a given discharge.
- **Preserve conveyance:** under the steady uniform-flow approximation, \(K=AR^{2/3}/n\), \(Q=K\sqrt{|S|}\). Set \(n'=n_0(A'/A_0)(R'/R_0)^{2/3}\).
- **Preserve both:** first match area through depth, then match \(K\) through roughness. For a wide rectangle widened by factor \(r\), \(D'=D_0/r\) and \(n'\approx n_0r^{-2/3}\).
- **Depth held fixed:** preserving \(K\) instead requires \(n'\approx r n_0\) for a wide rectangle, while artificial storage remains. This is a different approximation and must be declared.
- For varying sections, check integrated resistance: \(\Delta H\approx Q|Q|\int ds/K(s)^2\). Abrupt transitions, bends, bridges, and entrance losses need separate treatment.
- Apply compensation using **achieved**, not requested, sections; repeat after final bathymetry interpolation, clipping, and finishing. Map roughness into the actual FVCOM friction law; Manning values cannot simply be written as drag coefficients.

These formulas are proposed engineering approximations, not implemented or validated repository behavior. Effective roughness must remain within site-supported bounds; otherwise reject the surrogate.

**What compensation cannot recover:** matching \(A\) and \(K\) at one level does not preserve tidal storage over all levels, wetting/drying, wave speed, lateral jets, bend circulation, dispersion, bed stress, stratification, or habitat geometry. In a rectangular channel, \(dA/d\eta=W\); widening changes this even when reference-level area matches. Depth reduction can destroy a sill or salt wedge. No scalar roughness correction fixes that.

For port water quality or stratified exchange, prefer local refinement. Widening is acceptable only as a tested surrogate for the stated outputs. If neither works within budget, report that the chosen model resolution is inadequate; do not silently fill the basin.

**Proposed options**

The following are **new-schema defaults**, not changes to the certified defaults. Absent `waterways`, preserve existing behavior. Version 1 initially selects `legacy_sr`; `purpose_v1` is an explicit opt-in.

| Key under `waterways` | Units / default | Meaning |
|---|---|---|
| `version`, `profile` | `1`, `legacy_sr` | Freeze interpretation; alternative `purpose_v1`. |
| `purpose` | enum, `hydrodynamics` | Also `port_water_quality`, `tracer`; scoped overrides identify relevant basins/routes. |
| `action` | enum, `auto` | `auto`, `keep`, `refine`, `widen`, `fill`, `review`; explicit actions still face acceptance. |
| `representation` | enum, `refine_then_widen` | Also `refine_only`, `widen_only`; used for unresolved required features. |
| `target_h_m` | m, inherit local build sizing | Desired spacing; distinguish nominal sizing from measured edge lengths. |
| `min_rows` | count, `2` | Minimum realized transverse resolution; a topology floor, not a transport-convergence guarantee. |
| `min_implied_dt_s` | s, `null` | Required for automatic refinement/widening; null means review when timestep feasibility matters. |
| `max_width_factor` | ratio, `2.0` | Stationwise \(W'/W_0\) ceiling for automatic widening; provisional conservative bound. |
| `min_land_gap_m` | m, `150` | Retain existing guard initially; check final shared walls. This is a compatibility default, not a universal physical rule. |
| `compensation` | enum, `area_conveyance` | Also `area`, `conveyance`, `none`; non-default approximations require explicit evidence. |
| `max_effect_fraction` | fraction, `0.10` | Provisional maximum change in relevant storage, prism, exchange, travel/residence time, or receptor response against a declared reference. |
| `overrides` | list, `[]` | Stable feature selectors, local values, rationale, and evidence references. |

Existing model settings remain authoritative for minimum wet depth, friction law/bounds, vertical grid, grading, and runtime timestep. Do not duplicate them in this block. Missing physical data are `unknown`, never zero. `max_effect_fraction` is an engineering starting point requiring owner agreement, not a literature-derived universal tolerance.

**Evaluation order**

1. Validate units, coordinate systems, source coverage, feature identities, and conflicting overrides. Apply explicit data/domain corrections before measuring features.
2. Build a basin–channel graph from corrected **physical** geometry, retaining original widths and provenance. Use existing detectors as candidate generators, including water without OSM centerlines.
3. Mark required basins/routes from purpose regions, river sources, and through-path obligations. River and through-path preservation are standing invariants; exceptions need explicit source relocation or connectivity evidence.
4. Keep naturally resolved features. For unresolved, unprotected appendages, permit `fill` only with bounded physical effect; otherwise return `review`. Fill the designated basin/entrance complex, not just its throat.
5. For required unresolved water, try the declared representation order: refine at true width; then, if allowed, widen about the physical axis with compensation. Check timestep, width-factor ceiling, depth, and neighboring walls before meshing.
6. Normalize approved excluded water to land before final corridor realization, following 325’s fixed two-pass architecture. Preserve the physical-axis reference so normalization cannot redefine “natural width.”
7. Generate and finish once with bounded existing operators. Pass protected basin/route IDs into finishing; a failed protected channel becomes `review`, not an automatic deletion.
8. Evaluate final geometry and physical acceptance. Persist the decision and reasons; no unbounded regeneration-until-pass loop.

Overrides follow **global → basin/region → feature** precedence. Conflicting equal-priority overrides fail validation. Select by OSM IDs plus a geometry fingerprint, or a versioned local basin/arc ID; never by transient element number.

Illustrative fragments to merge into an existing build recipe; IDs and evidence paths are placeholders, not existing artifacts:
```yaml
# Port A: internal water quality is an objective.
waterways:
  version: 1
  profile: purpose_v1
  purpose: hydrodynamics
  target_h_m: 300
  min_implied_dt_s: 15
  overrides:
    - feature: basin:port_a
      purpose: port_water_quality
      representation: refine_only
      target_h_m: 50
      min_rows: 4
      evidence: studies/port_a_requirements.yaml
      reason: Preserve basin storage and stratified entrance exchange.
```
For an illustrative 240 m entrance, 50 m sizing requests roughly 4.8 spacings across; actual rows and transport convergence still need checking. This purpose override protects the basin and its necessary entrance network.

```yaml
# Port B: no internal receptor; omission has supporting evidence.
waterways:
  version: 1
  profile: purpose_v1
  purpose: hydrodynamics
  target_h_m: 300
  min_implied_dt_s: 15
  overrides:
    - feature: basin:port_b
      action: fill
      evidence: studies/port_b_omission_comparison.yaml
      reason: No river/source or through route; receiving-bay effect below 10%.
```
If Port B’s evidence is absent, filling is proposed for review, not automatically accepted. Conversely, an essential 400 m canal at \(h=300\) m could request 600 m width: factor 1.5. At original depth 6 m, area compensation gives 4 m; roughness and all-level exchange checks remain necessary.

**Relationship to the edit ledger**

Keep `recipes/edits/sample_repro/*.json` as reproducible geometry instructions. Add decision provenance linking policy overrides to generated or retained edits. Data corrections remain explicit evidence-backed inputs; purpose decisions can eventually generate deterministic `land_patch` or arc edits using existing operations.

Preserve filename order during migration, especially edit_004 before edit_005. Later replace implicit order with explicit dependencies. Never retire an edit merely because mesh QA passes without it. Record intended component connections explicitly; retain 342’s `allow_connect` safeguard, preferably scoped to named component pairs.

**Human judgment and printed evidence**

Human judgment remains in selecting objectives/receptors, accepting data corrections, identifying operational structures, bounding hydraulic uncertainty, and approving surrogate-channel assumptions. Automation should repeat those choices consistently.

Print a compact summary followed by ranked site records:

1. **Blocking:** lost source/receptor routes, protected basin deletion, unauthorized component mergers, wall failures, timestep/depth conflicts.
2. **High physical impact:** large widening factors, uncertain depth/roughness, stratified entrances, substantial lost storage, and dominant artificial exchange routes.
3. **Lower priority:** small peripheral closures with good evidence; ambiguous OSM artifacts.

Each record should show stable ID, atlas reference/coordinates, purpose, proposed action/reason, alternatives rejected, physical/requested/achieved width profiles, row count, wall minimum, depth/roughness changes, basin area/volume, source loads, timestep/cost estimate, uncertainties, and exact override syntax.

Provide a paired physical-versus-engineered map, section plots, and graph changes. Print missing evidence explicitly and conclude with separate statuses: `geometry_pass`, `physics_screen_pass`, `model_validation_pass`, `review_required`.

**Acceptance and detecting a bad choice**

- Retain all applicable `fmesh-mesh-qa` structural/OBC gates: positive elements, manifold boundaries, C1 ≥30°, C2 ≤130°, C4 normalized area difference ≤0.5, and valence ≤8.
- Gate the declared implied timestep. Current `qa.py` uses \(\min_e[L_{\min,e}/\sqrt{gH_{\max,e}}]\), without a safety factor, and gates it only when `min_dt_s` is supplied. Do not confuse this with an approved FVCOM runtime timestep.
- Use 342 for retained-route continuity and unauthorized connections, 346 for confirmed one-wide chokes, and 364 for stray resolution/wall integrity. Generalize their sample-specific references to the declared policy domain while retaining legacy checks.
- Audit **inside** intended widening tubes as well as outside them: edit-aware exclusions must not exempt an accidental connection from physical review.
- Require achieved rows, shared-wall clearance, compensated section area/conveyance, and acceptable storage curves over the operating water-level range. Unknown quantities cannot earn physical acceptance.
- For important widened/closed features, compare against a locally refined physical-geometry reference or observations using identical forcing: tidal amplitude/phase, prism, entrance discharge, freshwater balance, flow partition, and dye washout/travel time.
- Test relevant low/high discharge and tidal conditions; add salinity profiles and two-layer exchange for stratified sites. Water-quality objectives also need concentration/process sensitivity.
- Diagnose artificial dominance with \(f_j=P_j/\sum_kP_k\), where \(P_j\) is consistently measured incoming volume through entrance \(j\). Flag a surrogate supplying over half the basin’s inflow; dominance alone is not failure if the real entrance also dominates.
- Reject automatic acceptance when dominance is newly created, reference response changes exceed tolerance, or plausible width/depth/roughness uncertainty changes the selected action. Report residence-time distributions or tracer mass decay where a single flushing time hides recirculation.

**Migration preserving Tokyo Bay reproducibility**

1. Freeze the certified chain’s code/backend versions, source hashes, seeds, all resolved `SR_*` values, edit revisions/order, and intermediate geometry—not just final meshes. Reconcile the 115210 checker discrepancy separately.
2. Introduce a typed configuration adapter with no algorithm changes. Map `SR_H0`, `SR_GRADE`, `SR_DT`, `SR_CRMIN`, and CFL/OBC switches into build configuration; preserve their distinct meanings.
3. Map `SR_WATERWAYS`, `SR_NORMALIZE`, `SR_WIDEN_FACTOR`, `SR_ATTAIN_BAR`, `SR_FORCE2ROWS`, `SR_CH_REFINE`, `SR_CH_POLICY`, and `SR_EDITS_EXCLUDE` exactly into `legacy_sr`. Freeze scattered detector/finishing thresholds in that profile.
4. Keep environment aliases temporarily, but print the effective configuration and reject conflicting YAML/environment values. Establish the recipe wrapper by invoking the same certified stage functions/order.
5. Add feature IDs, purpose metadata, evidence reports, and physical measurements in shadow mode. Confirm unchanged geometry/connectivity and mesh output; document any platform-dependent numerical differences.
6. Enable `purpose_v1` feature by feature in separate output directories. Carry protection into finishing and compare decisions against the four human edits.
7. Add depth/friction compensation only with explicit model/export support; depth-field design is outside the existing `STAGE2_DESIGN.md` scope and needs a separate implementation stage.
8. Preserve `legacy_sr` permanently as a regression profile. Deprecate environment switches only after equivalent recipe runs pass the complete measured acceptance suite.

This follows `docs/DESIGN_HISTORY.md`: fix geometry upstream, constrain generation, avoid broad morphology that erases important water, and avoid repair cascades. It also follows the detector/planner/executor/reporting structure in `docs/STAGE2_DESIGN.md`; `docs/phase_h_user_guide.md` and `docs/detector_repair_matrix.md` describe repair capabilities, not physical-policy authority.

**External prior art: verified versus modeling judgment**

| System | Relevant practice and proposed lesson |
|---|---|
| SCHISM/RiverMapper | Verified documentation exposes fixed-width pseudo-channels and cross-channel row counts. Adopt explicit channel templates and network objects; this does **not** establish that widening automatically preserves hydraulics. [RiverMapper documentation](https://schism-dev.github.io/schism/v5.11/mesh-generation/meshing-for-compound-floods/generate-river-map.html) |
| ADMESH+ | The published method classifies narrow regions and generates coupled 1D–2D representations. Adopt width classification; do not interpret its 1D treatment as permission for one-cell-wide FVCOM channels. [Kang and Kubatko, 2024](https://gmd.copernicus.org/articles/17/1603/2024/) |
| OceanMesh2D | Automated spatial sizing provides the resolution machinery; application-specific keep/fill decisions remain a separate concern. Repository history documents additional constraints and cleanup lessons specific to this port. [OceanMesh2D paper](https://gmd.copernicus.org/articles/12/1847/2019/gmd-12-1847-2019.pdf) |
| Delft3D/D-Flow FM | D-Flow FM supports 1D and 2D domain parts: an alternative when explicit narrow-channel geometry is unaffordable. That capability is not established for this FVCOM pipeline. [Technical reference](https://content.oss.deltares.nl/dhydro/D-Flow_FM_Technical_Reference_Manual.pdf) |
| TELEMAC | Its older official manual recommends resolving the low-water bed with several points, supporting explicit transverse-resolution checks. Those recommendations are not interchangeable with FVCOM triangle gates. [TELEMAC-2D manual](https://www.opentelemac.org/downloads/Archive/v6p0/telemac2d_user_manual_v6p0.pdf) |
| SMS | Repository `STAGE2_DESIGN.md` records manual redraw, refinement, smoothing, and topology edits. From general modeling knowledge, the important transferable practice is preserving expert intent in feature-level instructions; current SMS-specific automation capabilities were not verified here. |

The compensation formulas, proposed tolerances, and purpose-dependent validation strategy are engineering synthesis requiring site/model verification. `docs/LITERATURE_SURVEY_AUTOMESH.md` is useful context, but its broader operational and novelty claims were not independently revalidated.

**Verification scope:** read-only inspection using `cat`, `sed`, and `rg`, plus primary-source web checks. No files were created or modified; no tests, mesh generation, or hydrodynamic runs were performed. Outstanding work is schema implementation, physical-data acquisition, FVCOM friction/export verification, model sensitivity testing, and reconciliation of the stored acceptance results.