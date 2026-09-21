# One-element-wide channels

`SR_ONE_WIDE=forbid|natural|allow` overrides `one_wide` in the recipe selected
by `SR_SIZING`. Omission means `forbid`, preserving the certified chain and its
existing connectivity exceptions.

## Two independent decisions

A channel policy answers two questions, and it is worth keeping them apart:

1. **Which channels are kept?** `detect_waterways` decides, from the shoreline
   geometry alone, whether a sub-resolution waterway is a resolve target at
   all. Bars: `min_canal_width_frac`, `min_resolve_width_frac`, and the
   normalize pass that fills unresolvable water.
2. **How wide is a kept channel carved?** `apply_waterway_policy` either pushes
   the banks into land until two standard rows fit, or cuts the channel at its
   natural width.

| setting | which channels are kept | how wide they are carved |
|---|---|---|
| `forbid` | strict bars, normalize runs | widened until two rows fit |
| `natural` | **same as forbid** | kept `port` / `dead-end` at natural width, one row |
| `allow` | bars lowered, normalize skipped | same as `natural` |

`relaxes_selection()` and `permits_one_row()` in `fvcom_mesh_tools.one_wide`
name the two axes; every caller states which one it cares about.

## Why `natural` exists

`allow` moves both decisions at once, and the measured cost is the selection
change, not the single rows. On Tokyo Bay (job 115255):

| | `forbid` | `allow` |
|---|---|---|
| kept channels | 29 | 51 |
| confirmed one-wide cells | 14 | 92 |
| real chokes | 6 | 16 |
| QA | 21/21 PASS | 4 failures |
| implied dt | 16.34 s | 15.40 s |

Of the ten elements behind those four failures, only **two** lie inside a
confirmed one-wide site: the failures come from meshing far harder geometry,
not from the rows. So "permit one row" was never available on its own.

The reason to want it on its own is shoreline fidelity. Widening pushes
**1,702 ha** of Tokyo Bay coastline into the water, which makes a channel look
considerably larger than it is. That shows up directly in the coastline
measurements (`docs/coast_fit.md`): the fitted mesh sits a median of 10 m from
the polygon it was given but 26 m from raw OSM, and the difference is our own
preprocessing, of which widening is the largest part.

## What each setting does in detail

`natural` and `allow` permit one row only when the policy `kind` is `port` or
`dead-end`. Records classified as `canal` or `through` can carry connections
and retain two-row widening, attainability checks, branch width treatment and
minimum refinement rows. The classification comes from `detect_waterways` and
is consumed directly by `apply_waterway_policy`, never re-derived there. The
rule applies to each classified network, including its branches.

For a kept port or dead end, width-only closure and attainability rejection are
disabled. Permission for one row never turns a `close` into a `keep`.

Global feature sizing differs: `allow` drops the size field to one row, while
`natural` leaves it at `SR_FS` exactly as `forbid` does — `natural` changes the
carve, not the field. Optional channel refinement uses the per-record minimum
(one for ports and dead ends, two otherwise) under both.

`SR_WIDEN_FACTOR` and `SR_ATTAIN_BAR` supply the canal/through baseline; ports
and dead ends override these with one-row settings. `SR_FS` cannot override
`allow`'s global feature sizing. `SR_CH_REFINE` retains its on/off meaning.
`allow` adds no fixed points and disables forced ladders even if
`SR_FORCE2ROWS` is set. The 1.2 mesh-edge calibration, size/CFL floors, manual
edits and quality repair remain under every setting.

Width-only normalization, finishing throat pruning and choke splitting are
disabled whenever one row is permitted. Canal/through protection is established
by the geometry policy before meshing. DistMesh does not guarantee realized row
counts; QA, notebook 342 connectivity, notebook 346 one-wide flags and implied
dt must still be checked before certification. **No QA gate changes.**

## Metadata

Generation writes `outputs/sample_repro/channel_policy.json`. Finishing reuses
that setting and rejects explicit conflicts; legacy meshes without metadata use
the environment/recipe/default. Regenerate when changing policy.

## Jobs

- `qsub jobs/octopus/402_one_wide_continuity.sh` — `allow` against achieved and
  field `forbid`. Results in `outputs/one_wide_402.<jobid>/`.
- `qsub jobs/octopus/408_natural.sh` — all three settings under one coastline
  fit. Results in `outputs/natural_408.<jobid>/`.
