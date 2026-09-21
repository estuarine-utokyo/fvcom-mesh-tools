# One-element-wide channels

`SR_ONE_WIDE=forbid|allow` overrides `one_wide` in the recipe selected by
`SR_SIZING`. Omission means `forbid`, preserving the certified chain and its
existing connectivity exceptions. No new option or environment variable is required.

`allow` permits a channel to remain one element wide only when its policy
`kind` is `port` or `dead-end`. Records classified as `canal` or `through`
can carry connections and retain the existing two-row widening, attainability
checks, branch width treatment and minimum refinement rows. The classification
comes from `detect_waterways`; it is consumed directly by
`apply_waterway_policy`, never re-derived there. The rule applies to each
classified network, including its branches.

Kept ports and dead ends may use one row. Their width-only closure and
attainability rejection are disabled.
The existing keep/close/ignore decisions, basin/extent selection and barrier
guards still apply; permission for one row does not turn a close into a keep.
Global feature sizing uses one row; optional channel refinement uses the
per-record minimum (one for ports/dead ends, two otherwise). The 1.2 mesh-edge
calibration, size/CFL floors, manual edits and quality repair remain.

`SR_WIDEN_FACTOR` and `SR_ATTAIN_BAR` supply the existing canal/through-treatment
baseline; ports and dead ends override these with one-row settings. `SR_FS`
cannot override allow's global feature sizing. `SR_CH_REFINE` retains its
on/off meaning. Allow adds no fixed points and disables forced ladders even
if `SR_FORCE2ROWS` is set. Forbid retains all legacy tuning behavior.
The failed kept-bank constraints, thinning and continuity gates have been removed.

In allow mode width-only normalization and finishing throat pruning and
choke splitting/widening remain disabled. Canal/through protection is established
by the geometry policy before meshing. DistMesh does not guarantee realized
row counts; QA, notebook 342 connectivity, notebook 346 one-wide flags and
implied dt must still be checked before certification. No QA gates change.

Generation writes `outputs/sample_repro/channel_policy.json`. Finishing reuses
that setting and rejects explicit conflicts; legacy meshes without metadata
use the environment/recipe/default. Regenerate when changing policy.

Run `qsub jobs/octopus/402_one_wide_continuity.sh` from the repository root to
compare allow, achieved forbid and field forbid. Results go to a new
`outputs/one_wide_402.<jobid>/` directory. The historical job filename is retained.
