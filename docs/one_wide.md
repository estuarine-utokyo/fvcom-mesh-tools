# One-element-wide channels

`SR_ONE_WIDE=forbid|allow` overrides `one_wide` in the recipe selected by
`SR_SIZING`. Omission means `forbid`, reproducing the certified chain including
its existing connectivity exceptions (the name is a policy, not a guarantee
that the ledger is empty). `recipes/sizing/tokyo_bay.yaml` declares that default.

`allow` permits one row without asserting that it is harmless. There is no
`allow_if_harmless`: severance (342), topology, QA angles/areas, and implied dt
are not combined into a new automatic acceptance test. Existing quality and
CFL processing and all reports remain active. Inspect those reports before
certifying a new mesh.

The switch controls feature sizing (3 -> 1), channel refinement minimum rows
(2 -> 1), geometry detection width floors, branch width floors, widening
(2 -> 1), attainability rejection, forced ladders, thin-stub and residual-stub
closure, duplicate width checks, width-based normalization, finishing throat
pruning, and finishing choke splitting/widening. In allow mode normalization
returns no fills: width alone is insufficient evidence to declare water land.
The 1.2 mesh-edge calibration remains; it is not a two-row multiplier.
Basin/extent selection, barrier guards, manual geometry/mesh edits, minimum
size and CFL floors, and quality repair remain. These can change or exclude
water, and explicit manual edits can still widen a channel. DistMesh does not
guarantee a particular realized row count.

`SR_FS`, `SR_WIDEN_FACTOR`, `SR_ATTAIN_BAR`, and `SR_FORCE2ROWS` retain their
legacy behavior under forbid and cannot override allow. `SR_CH_REFINE` retains
its on/off meaning and uses the selected minimum row count. Regional sizing
can independently request finer cells. QA and notebook 346 keep reporting
one-wide cells in both modes.

Generation writes `outputs/sample_repro/channel_policy.json` alongside the
UTM mesh. Finishing reuses that setting and rejects explicit conflicts; legacy
meshes without metadata use the environment/recipe/default. Regenerate to
change policy consistently across stages. No change to the QA gates is made.
