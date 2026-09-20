# Fixed-h channel representation: one-wide vs widened (study)

Status: STUDY, 2026-09-20. Drafted by gpt-6-astra (Codex CLI) at the
owner's request, reviewed by Claude. Nothing here is implemented and no
simulation was run for it.

The question, from the owner: the coastal target size h is FIXED (cost).
Under that constraint, is it better to widen a sub-h channel until it
has two element rows -- which makes the modelled channel several times
wider than reality -- or to keep the true width with a one-element-wide
channel, which this project currently forbids? Filling it as land,
anisotropic elements, partial retention and a river boundary condition
are the other candidates. "Use a finer mesh" is not an answer here.

This study supersedes the corresponding parts of
`channel_resolution_policy_design.md` (whose `min_rows: 2` /
`max_width_factor: 2.0` defaults it contradicts).

---

**Scope and evidence.** This is a read-only design study; no files changed and no simulations ran. “Source finding” below means this checkout’s implementation; formulas, screening thresholds and proposed tests are engineering analysis requiring validation. The coastal target \(h\) stays fixed in every production candidate. Smaller transverse edges inherent in true-width strips are accounted for explicitly in timestep and cost; they are not a change to the coastal sizing field.

**1. What one-element-wide water actually means in this FVCOM build**

The checked configuration selects **GCN, without ghost cells**, has `DYE_RELEASE`, and disables `WET_DRY` and GOTM: [make.inc:350](/octfs/work/G16445/v61021/Github/FVCOM/src/make.inc:350), [make.inc:240](/octfs/work/G16445/v61021/Github/FVCOM/src/make.inc:240), [make.inc:362](/octfs/work/G16445/v61021/Github/FVCOM/src/make.inc:362), [make.inc:438](/octfs/work/G16445/v61021/Github/FVCOM/src/make.inc:438). Conclusions about its bank treatment should not be transferred automatically to ghost-cell builds.

| Source finding | Consequence for a strip between two banks |
|---|---|
| Momentum uses triangle-centred velocities and fluxes across triangle edges: [advave_edge_gcn.F:149](/octfs/work/G16445/v61021/Github/FVCOM/src/advave_edge_gcn.F:149). Nodal control-volume faces connect centroids to edge midpoints: [tge.F:614](/octfs/work/G16445/v61021/Github/FVCOM/src/tge.F:614). Elevation advances by `ELF=ELRK-DTK*XFLUX/ART1`: [extel_edge.F:997](/octfs/work/G16445/v61021/Github/FVCOM/src/extel_edge.F:997). | Bank nodes still have water-filled tracer control elements (TCEs). Having no interior nodes does **not** eliminate water storage or longitudinal transport. This staggering also agrees with the [FVCOM user manual, p. 46](https://etchellsfleet27.com/wp-content/uploads/2020/06/FVCOM_User_Manual_v3.1.6.pdf). |
| For one solid edge, the external velocity is projected onto the bank tangent; the internal-layer velocities receive the same projection: [bcond_gcn.F:393–434](/octfs/work/G16445/v61021/Github/FVCOM/src/bcond_gcn.F:393), [bcond_gcn.F:490–535](/octfs/work/G16445/v61021/Github/FVCOM/src/bcond_gcn.F:490). | This is an element-velocity constraint, not nodal no-slip. A straight strip can carry nonzero along-channel flow. Every bank-adjacent triangle loses its normal component; a strip composed entirely of such triangles cannot represent an independent interior lateral circulation. |
| Two solid edges produce `ISBCE=3`: [tge.F:584](/octfs/work/G16445/v61021/Github/FVCOM/src/tge.F:584). Both horizontal velocity components are zeroed in external and internal modes: [bcond_gcn.F:381](/octfs/work/G16445/v61021/Github/FVCOM/src/bcond_gcn.F:381), [bcond_gcn.F:485](/octfs/work/G16445/v61021/Github/FVCOM/src/bcond_gcn.F:485). | A corner/cap triangle can be stagnant by construction. That is different from an alternating strip whose triangles each have one bank edge. |
| **GCN sets all momentum reconstruction coefficients `A1U/A2U` to zero in solid-boundary triangles**: [shape_coef_gcn.F:213–221](/octfs/work/G16445/v61021/Github/FVCOM/src/shape_coef_gcn.F:213), [shape_coef_gcn.F:312–331](/octfs/work/G16445/v61021/Github/FVCOM/src/shape_coef_gcn.F:312). This routine is selected at [grid_metrics.F:71](/octfs/work/G16445/v61021/Github/FVCOM/src/grid_metrics.F:71). | In the canonical strip, reconstruction is piecewise constant throughout. This is a concrete loss of spatial accuracy, including along-channel accuracy, rather than simply too few samples across the width. Near-collinear boundary-neighbour centroids also trigger the explicit determinant failure at lines 314–320. |
| Momentum advection is upwind, using reconstructed velocities; viscous stress uses those reconstruction gradients: [adv_uv_edge_gcn.F:338](/octfs/work/G16445/v61021/Github/FVCOM/src/adv_uv_edge_gcn.F:338), [387–401](/octfs/work/G16445/v61021/Github/FVCOM/src/adv_uv_edge_gcn.F:387), [508–550](/octfs/work/G16445/v61021/Github/FVCOM/src/adv_uv_edge_gcn.F:508). External-mode counterparts are [advave_edge_gcn.F:340](/octfs/work/G16445/v61021/Github/FVCOM/src/advave_edge_gcn.F:340) and [447–475](/octfs/work/G16445/v61021/Github/FVCOM/src/advave_edge_gcn.F:447). | Where both triangles are solid-boundary cells, these viscous stresses vanish because both gradients vanish, even with constant mixing. Momentum advection reduces to first-order reconstruction. Merely increasing the viscosity coefficient cannot recover the missing stress there. |
| Edge boundary flags suppress external momentum flux; internal advection is suppressed but its viscous boundary factor differs: the flux lines just cited. | Do not describe this implementation as a generic, resolved lateral no-slip boundary layer. The whole strip may lack the interior stencil needed to develop a credible cross-channel shear profile. |
| Bottom drag survives: logarithmic `CBC`, its minimum, and bottom-layer quadratic stress are computed at [brough.F:88–127](/octfs/work/G16445/v61021/Github/FVCOM/src/brough.F:88), [187–214](/octfs/work/G16445/v61021/Github/FVCOM/src/brough.F:187). | One row does **not** remove bottom friction or vertical shear. A straight, well-mixed, bottom-friction-dominated connector can still have useful bulk hydraulics. Missing bank drag, lateral shear and bend losses remain physical/model errors. |
| The internal step runs `ISPLIT` external steps, then layer momentum; transport adjustment reconciles the modes: [internal_step.F:933](/octfs/work/G16445/v61021/Github/FVCOM/src/internal_step.F:933), [1212](/octfs/work/G16445/v61021/Github/FVCOM/src/internal_step.F:1212), [adjust2d3d.F:141–176](/octfs/work/G16445/v61021/Github/FVCOM/src/adjust2d3d.F:141). Vertical continuity uses layer fluxes and nodal volumes: [vertvl_edge.F:126](/octfs/work/G16445/v61021/Github/FVCOM/src/vertvl_edge.F:126), [598–627](/octfs/work/G16445/v61021/Github/FVCOM/src/vertvl_edge.F:598). | One row is not intrinsically a mode-split failure. Short transverse dimensions still constrain external waves, internal advection and transport consistency. Stability and conservation must be measured separately. |
| Dye uses nodal reconstruction, diffusion and antisymmetric TCE flux deposits: [mod_dye.F:1213](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_dye.F:1213), [1369–1452](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_dye.F:1369). Nodal mixing has its own boundary-aware stencil: [viscofh.F:58–119](/octfs/work/G16445/v61021/Github/FVCOM/src/viscofh.F:58). | Scalars are **not** automatically frozen or subject to the same zero-gradient rule as momentum. Two bank-node chains can transport dye, but unresolved lateral mixing, first-order momentum and return-flow errors can badly bias dispersion and residence time. |

**Interpretation:** a straight one-wide strip can approximate a hydraulic connector; it cannot be presumed to reproduce a channel velocity field. Parallel-bank projection permits an approximately plug-like lateral flow, with vertical structure retained. At bends, alternating bank tangents constrain successive velocities differently, potentially introducing artificial losses and suppressing secondary circulation. Their magnitude requires testing.

The project’s strict detector flags three land-boundary vertices with **at most one** boundary edge: [channel_policy.py:199–203](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/channel_policy.py:199). That is not FVCOM’s two-solid-edge zero-velocity condition. The existing blanket prohibition therefore combines genuine stencil concerns with a broader geometric rule.

Without wetting–drying, require \(D(s)+\eta_{\min}(s)>\text{MIN_DEPTH}\) everywhere, with uncertainty margin. Compensated shallow channels cannot be allowed to dry and then trusted to close correctly. Wet/dry-specific masks in [mod_dye.F:1357](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_dye.F:1357) are inactive here.

**2. Quantifying widening and compensation**

The following are derived screening equations, not FVCOM validation results. Assume a rectangular channel of fixed length \(L\), initially \(A_0=W_0D_0\), widened by \(r=n_{\rm rows}h/W_0\). Set \(d=D'/D_0\), \(c=C_d'/C_{d0}\), and \(a=A'/A_0=rd\).

For well-mixed, bottom-drag-dominated flow, \(Q=K_C\sqrt{|S|}\), with \(K_C=A\sqrt{gD/C_d}\). Thus \(k=K_C'/K_{C0}=rd^{3/2}/\sqrt c\). This treats effective bottom velocity as proportional to section-mean velocity; stratification or changed vertical profiles invalidate a simple coefficient mapping.

| Representation | \(d\) | \(c\) | Section area / mean channel storage \(a\) | Conveyance \(k\) | Channel travel time \(V_c/Q\), same head |
|---|---:|---:|---:|---:|---:|
| No compensation | 1 | 1 | \(r\) | \(r\) | 1 |
| Area only, fixed effective drag | \(1/r\) | 1 | 1 | \(r^{-1/2}\) | \(r^{1/2}\) |
| Conveyance only, fixed depth | 1 | \(r^2\) | \(r\) | 1 | \(r\) |
| Area and conveyance | \(1/r\) | \(1/r\) | 1 | 1 | 1 |

“Same roughness length” is not “same \(C_d\)”: this build recalculates drag from bottom-layer height, uses a 3 m depth floor in the original law, and applies `CBCMIN` ([brough.F:115–127](/octfs/work/G16445/v61021/Github/FVCOM/src/brough.F:115)). It also has a user-defined drag path ([brough.F:154–157](/octfs/work/G16445/v61021/Github/FVCOM/src/brough.F:154)). Compensation must be checked against the **realised** drag and bottom velocity.

For Manning conveyance, use \(K_M=AR^{2/3}/n\), \(R=WD/(W+2D)\); exactly, \(n'/n_0=a(R'/R_0)^{2/3}\) preserves \(K_M\). In the wide-channel limit, area-only gives \(k=r^{-2/3}\), and preserving both requires \(n'/n_0=r^{-2/3}\). Manning \(n\), roughness length and FVCOM drag are not interchangeable.

Storage and prism require separate accounting:

- At mean level, \(V_c'/V_c=a\). For channel-plus-basin volume, \(V'/V=1+f_V(a-1)\), where \(f_V=V_c/(V_b+V_c)\).
- At any common elevation \(\eta\), \(A'(\eta)-A_0(\eta)=W_0[(rd-1)D_0+(r-1)\eta]\). Even area compensation leaves relative section error \((r-1)\eta/(D_0+\eta)\).
- **Every widening has \(dV_c'/d\eta=rWL\)**. Channel prism at the same tidal range increases by \(r\), regardless of depth or roughness compensation.
- If basin and channel share the same tidal range, total prism ratio is \(1+f_S(r-1)\), \(f_S=WL/(B+WL)\), where \(B\) is basin surface area. Actual basin amplitude changes, so this is a storage screen, not a prediction.
- For \(\tau\approx V/(\epsilon Q_{\rm exchange})\), \(\tau'/\tau=[1+f_V(a-1)](\epsilon/\epsilon')/(Q_{\rm exchange}'/Q_{\rm exchange})\). Under equal head and replacement efficiency, substitute \(k\) for the discharge ratio. For basin-only residence, the volume factor is 1.

Useful **10% error crossings** under those assumptions:

| Quantity | First failure as \(r\) increases |
|---|---|
| Uncompensated section, channel storage, conveyance; any widening’s channel storage slope | \(r>1.10\) |
| Uncompensated basin residence, equal head | \(r>1/0.9=1.111\) |
| Area-only channel travel time, fixed \(C_d\) | \(r>1.1^2=1.21\) |
| Area-only conveyance, fixed \(C_d\) | \(r>1/0.9^2=1.235\) |
| Area-only Manning conveyance | \(r>0.9^{-3/2}=1.171\) |
| Total storage-slope/prism screen | \(r>1+0.10/f_S\) |
| Area-compensated section at low tide \(-a_\eta\) | \((r-1)a_\eta/(D_0-a_\eta)>0.10\) |
| Permanently wet area-compensated channel | Failure when \(D_0/r-a_\eta\le D_{\min}\) |

Example: \(h=300\) m, \(W_0=150\) m, \(D_0=6\) m, two-row width 600 m gives \(r=4\). Uncompensated area and nominal conveyance are **300% too large**. Area-only yields \(D'=1.5\) m and half the original conveyance at fixed drag. Matching both needs \(C_d'=C_{d0}/4\); if \(C_{d0}=0.003\), this means 0.00075. The working M2 preparation specifies minimum drag 0.003 ([383_m2_case_prep.py:152–153](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:152)); that configuration cannot realise this correction through the original roughness law.

At a 0.5 m low tide, this compensated example has section area \(600\) versus \(825\ {\rm m^2}\): **27% too small**, despite exact mean-level area.

Matching area and conveyance preserves reference-level mean velocity at given \(Q\), nominal hydraulic resistance, channel volume and lumped inertia. It does **not** preserve storage slope, all-level resistance, local wave speed \(\sqrt{gD}\), bottom stress, vertical mixing, lateral jets or tracer replacement efficiency. Conveyance-only also changes momentum flux \(\rho Q^2/A\) by \(1/r\) at fixed \(Q\).

**3. Computable conditions that change the preferred representation**

Geometry alone cannot uniquely determine flow error. The same channel can be negligible under one forcing and essential under another. Available data support the following reproducible screening; absent depth, discharge or roughness evidence remains unknown.

| Compute from available inputs | What it establishes; proposed flip criterion |
|---|---|
| Width profile from polygon intersections normal to OSM axes; \(W_0/h\), \(r=2h/W_0\), uncertainty intervals | If uncompensated section/conveyance must be within 10%, require \(r\le1.1\). Thus genuine sub-\(h\) channels, where \(r>2\), cannot justify uncompensated two-row widening through geometry alone. |
| DEM depth, \(D_0/h\), \(D_0/W_0=(D_0/h)/(W_0/h)\), sills | Depth uncertainty directly affects area and resistance. Large \(D/W\) weakens bottom-only drag assumptions because bank effects matter. Area compensation fails if the low-water/depth bound above fails. |
| Channel length, \(L/h\), bend radius \(R_b\), bank-direction change over one element | A long, straight, slowly varying strip is a candidate bulk connector. As a provisional screen require \(h/R_b\le0.1\) and fractional width/depth change per \(h\le0.1\); larger changes make the GCN projection and first-order reconstruction more consequential. |
| Basin surface area \(B\), DEM volume \(V_b\), \(f_S,f_V\) | These determine artificial storage bounds. For widening, test \((r-1)f_S\le0.1\). For filling, rank omitted storage/prism relative to the receiving bay; a small entrance does not imply a small basin effect. |
| Polygon connectivity graph, alternate route, river-list intersections | A through path affects flow partition even without a large terminal basin. A dead end can still dominate tidal storage or contain a tracer receptor. Neither may be filled solely because \(W/h\) is small. |
| \(\Lambda=C_dL/D\), with an explicit uncertain drag bracket | Measures integrated frictional head loss, \(\Delta H_f\approx\Lambda U^2/g\). Large \(\Lambda\) favours assessing bulk resistance; small \(\Lambda\) makes entrance inertia, jets and bends relatively more important. |
| Actual candidate triangle altitudes, aspect ratios, determinant conditioning and cost | True-width one-wide and anisotropic candidates can reduce timestep despite fixed nominal \(h\). Screen using transverse altitude and wave speed, then measured stability; edge length alone is insufficient for very skew triangles. |

For quantitative tidal screening, solve a cheap **lumped channel–basin model**, explicitly an approximation:
\[
I\dot Q+JQ|Q|=\eta_m-\eta_b,\qquad B\dot\eta_b=Q,\qquad
I=L/(gA),\quad J=C_dL/(gDA^2).
\]
With equivalent linear resistance \(R=(8/3\pi)J\widehat Q\),
\[
\widehat\eta_b/\widehat\eta_m=[1-\omega^2IB+i\omega RB]^{-1}.
\]
Use the working M2 frequency and declared amplitude scenarios; iterate \(\widehat Q\). Add distributed channel storage when \(WL/B\) is appreciable. This uses geometry plus declared forcing/drag assumptions, not unavailable observations.

The response flips around \(\omega RB\sim1\), or near \(\omega^2IB\sim1\): resistance/inertia changes then strongly affect amplitude and phase. Well below both, basin tide nearly follows the bay; widening cannot multiply its prism by \(r\), although channel storage still increases. Well above the friction threshold, exchange is restricted and altered conveyance becomes decisive.

For one-wide screening, estimate longitudinal donor-cell diffusion as \(K_{\rm num}\sim U\ell/2\), with actual along-flow cell scale \(\ell\). This is an order-of-magnitude modified-equation estimate, **not an FVCOM error bound**. Test whether \(K_{\rm num}/(UL_g)=\ell/(2L_g)\) exceeds 0.1 for the acceleration or tracer-gradient length \(L_g\). Short jets/fronts can fail while section-integrated volume remains useful.

Objectives change the decision: volume needs correct storage and \(Q\); momentum additionally needs \(Q^2/A\), direction and mouth geometry; tracers need storage, source routes, return flow and mixing. Passing a prism test does not establish a credible residence time.

A river BC is a strong candidate when the omitted reach supplies known freshwater/load and its tidal storage and delay are negligible for the objective. FVCOM injects nodal/edge volume ([extel_edge.F:523–535](/octfs/work/G16445/v61021/Github/FVCOM/src/extel_edge.F:523)) and momentum using discharge, source area, direction and vertical distribution ([adv_uv_edge_gcn.F:1492–1523](/octfs/work/G16445/v61021/Github/FVCOM/src/adv_uv_edge_gcn.F:1492)). A river list alone does not provide those time series.

Do **not** substitute a river node for a reversibly exchanging back basin: river velocity handling clips the incoming-direction component ([bcond_gcn.F:403–407](/octfs/work/G16445/v61021/Github/FVCOM/src/bcond_gcn.F:403)). Elevation OBCs diagnose boundary flux from continuity ([mod_obcs.F:1238–1255](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_obcs.F:1238)); they do not automatically supply the missing basin impedance. A bidirectional storage/flux coupling would need separate implementation and verification.

Partial retention is preferable when the seaward mouth/jet matters but the omitted upstream reach contributes less than the accepted storage, resistance, delay and load errors. Preserve or relocate river inputs explicitly. Closing just the neck of an important tidal basin is not a harmless truncation.

**4. Experiment that would decide this**

**Geometry and forcing, proposed:** rectangular main bay \(6\times3\) km, depth 10 m; channel \(L=1.8\) km, \(W_0=150\) m, \(D_0=6\) m; square back basin of area 1.44 km², depth 5 m. Use identical smooth mouth transitions defined in physical coordinates. Production \(h=300\) m; channel longitudinal spacing approximately \(h\). Force the offshore edge with M2 amplitude 0.5 m, uniform phase, initially no wind, river or density gradient, zero Coriolis, ten sigma layers and native closure.

Use a fixed physical sponge at the offshore boundary, remote from the channel, identical across variants. Reuse spectral-file construction ([383_m2_case_prep.py:111–121](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:111)), sponge-writing machinery ([283–297](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:283)), and M2/M4/M6 least-squares analysis ([384_m2_analysis.py:35–54](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/384_m2_analysis.py:35)).

| Run per geometry | Purpose |
|---|---|
| F8: true width, eight transverse intervals; reference bay spacing 75 m | Fully resolved physical-geometry reference, outside the production cost constraint. |
| F12: twelve intervals; reference bay spacing 50 m | Check that reference prism/phase are within 2%, dye time within 5%; otherwise reference remains unqualified. |
| O1: true width, one transverse interval, alternating diagonals | Test the forbidden bank-to-bank strip directly. |
| W0: 600 m width, two transverse intervals, original depth/drag | Measure uncompensated widening bias. |
| WA: same width, depth 1.5 m, original effective drag | Isolate area compensation. |
| WK: same width/depth as W0, drag multiplied by 16 | Isolate conveyance compensation. |
| WAK: depth 1.5 m, effective drag divided by four | Test simultaneous reference-level area/resistance preservation. |
| X: fill channel and back basin | Measure omission’s effect on the main bay; back-basin exchange is lost, not a successful zero-error result. |
| A2: true width, two transverse intervals of 75 m, longitudinal spacing 300 m | Test anisotropy with unchanged coastal \(h\), including its extra timestep cost. |

Use controlled user-defined drag for the analytical comparison, initially \(C_{d0}=0.003\); separately verify whether the production roughness law can realise the successful surrogate. Do not silently clip WAK’s requested drag.

Run these nine variants for back-basin areas **0.36, 1.44 and 5.76 km²**: **27 simulations** spanning storage/loading. Add F8/O1/A2 for a smooth 90° bend: **three more**. Add three timestep-halving and three mixing-sensitivity reruns selected from the central case: **36 simulations**, provisionally one job each, plus one preparation and one analysis job: **38 jobs**.

Ramp over two M2 periods, spin up for ten periods, and require successive-cycle discharge/amplitude changes below 1%. Then initialise \(C=1\) in the back basin and \(C=0\) elsewhere through a staged restart; switch off releases and decay. Dye startup choices and decay controls exist ([mod_dye.F:384](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_dye.F:384), [5017–5027](/octfs/work/G16445/v61021/Github/FVCOM/src/mod_dye.F:5017)); the spatial restart preparation must be added and checked. Run at least 40 further periods, extending until an e-folding is observed or reporting a censored lower bound.

Measure:

- Incoming and outgoing section-integrated volumes per cycle, net discharge and \(V_b^{\max}-V_b^{\min}\); distinguish exchange prism from net transport.
- Complex M2 discharge and basin elevation, phase lag, M4 distortion and channel head loss.
- Main-bay complex M2 elevation/current errors on common physical sampling locations, including mouth jets.
- Basin dye inventory, concentration and first e-folding time; also export/return fraction and the full decay curve.
- Water and dye budget residuals, positivity, uniform-tracer preservation in a diagnostic rerun, bottom stress, cross-channel velocity structure, and 2-D versus depth-integrated 3-D transport.
- Triangle quality, node/element counts, shortest altitude, implied and actual stable timestep, and core-hours.

Proposed acceptance: exchange/prism and dye e-folding within **10%**, basin phase within **5°**, main-bay complex M2 elevation error below **2% of forcing amplitude**, and normalized water/dye budget residuals below **0.1%**. Near-zero quantities need absolute tolerances. A numerically stable run alone passes none of the physical criteria.

Use a shared conservative timestep initially, then benchmark each candidate’s largest verified timestep. A2 is not free: its transverse spacing is half O1’s. If “fixed cost” also forbids the resulting timestep reduction, it is infeasible regardless of its accuracy.

Runtime anchor: the existing 11-day, \(DTE=5\) s M2 cases took **138–148 s on eight ranks**, from [383_m2.115182.log:98–193](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/logs/383_m2.115182.log:98), with duration/timestep at [383_m2_case_prep.py:40–45](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/notebooks/383_m2_case_prep.py:40). Budget roughly **10–60 minutes per coarse run and 1–4 hours per reference**, 15–50 aggregate eight-rank job-hours before long tracer extensions; these are planning estimates to replace after a two-period pilot.

Follow the existing OCTOPUS runner’s environment and binary hashing ([383_m2.sh:31–58](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/jobs/octopus/383_m2.sh:31)). A future driver could be submitted with `qsub jobs/octopus/channel_fixed_h_study.sh`; **that proposed script does not exist and was not created or submitted in this study**.

Evidence favouring O1 would be accurate exchange/phase and dye despite its stencil limitations, at lower physical bias/cost than widening. WAK wins only if it meets those metrics across basin sizes and does not require unacceptable drag/depth changes. A2 wins if it improves momentum/tracer behaviour enough to justify measured cost. X wins only for a main-bay-only objective when omission passes its receiving-water criteria.

**5. Policy after that analysis**

**Replace the universal two-row prohibition with a representation-specific acceptance rule.** At fixed \(h\), true-width one-wide water is often the better *candidate* for a straight hydraulic connector than widening it several-fold. It is not yet a validated default for channel velocity fields or tracer residence. Large uncompensated widening is the least defensible automatic choice; filling and river substitution can be better when their omitted processes are demonstrably irrelevant to the objective.

The proposal currently specifies `refine_then_widen`, universal `min_rows: 2`, and `max_width_factor: 2.0` ([channel_resolution_policy_design.md:75–87](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/docs/channel_resolution_policy_design.md:75)). Those defaults do not answer this fixed-\(h\) problem.

| Proposed option mapping—not implemented | Numeric argument and decision |
|---|---|
| Keep `target_h_m`, `purpose`, `action`, `overrides`, `max_effect_fraction` | Lock the supplied \(h\); use `max_effect_fraction: 0.10` for exchange/storage/residence, with separate phase/main-bay tolerances above. |
| Add `representation: fixed_h_compare` | Compare true-width strip, anisotropic, widening and omission/BC candidates. Remove refinement from this profile’s action order. |
| `action: keep`, new `representation: one_row`, `min_rows: 1` | Permit only a protected, tested connector; screen bank-direction/section change at 0.1 per element and require the timestep/cost budget. Retain explicit zero-velocity-corner and stencil-conditioning checks. |
| New `representation: anisotropic`, `min_rows: 2`, initial `max_aspect_ratio: 4` | Trial bound matching A2; validate conditioning, timestep and transport. The current 30° minimum-angle gate requires an explicit scoped exception, not a claim of passing unchanged QA ([qa.py:464–465](/octfs/work/G16445/v61021/Github/fvcom-mesh-tools/src/fvcom_mesh_tools/qa.py:464)). |
| `action: widen`, `compensation: none` | With a 10% hydraulic tolerance, set `max_width_factor: 1.10`, not 2.0. Most genuinely unresolved channels will fail this screen. |
| `action: widen`, `compensation: area_conveyance` | No universal acceptable \(r\). Enforce wet-depth, realised-drag, all-level section and \((r-1)f_S\le0.10\) screens, followed by model validation; otherwise `review`. |
| `action: fill` | Require receiving-water effect below tolerance, no required through route/receptor, and explicit treatment of every river/load. |
| New `representation: partial` or `river_bc` | Record retained length/source location; require omitted storage, delay and resistance effects below 10% for the objective. Preserve discharge and constituent loads; specify momentum area/direction if relevant. |
| Keep `min_implied_dt_s` and `review` | Add a measured core-hour budget. If no candidate passes, return `review` with the quantified failure; do not silently widen or delete. |

Finishing must preserve approved one-row features instead of rediscovering and removing them through the existing strict flag. Two nominal rows must also be checked for actual interior momentum stencils; row count alone does not establish them.

**Verification and open work:** inspected guidance, policy, solver kernels, build flags, M2 preparation/analysis and stored logs using `rg`, `nl`, `sed` and `cat`; cross-checked staggering against the FVCOM manual. Some exploratory filename lookups failed and were resolved to the actual source files. No tests, mesh generation or FVCOM runs were performed; no new counts, QA gates or timestep measurements are claimed. The empirical ranking, spatial dye staging, production drag feasibility and fixed-cost acceptance remain to be established by the proposed experiment.