"""Single policy switch for width-based channel enforcement.

Three settings, which differ along two independent axes -- WHICH channels are
kept, and HOW WIDE they are carved:

``forbid``   the certified chain: only channels that can be made two standard
             rows wide are kept, and their banks are pushed into land until
             they are.  Connectivity exceptions preserved.
``natural``  the same selection as ``forbid``, but a kept ``port`` or
             ``dead-end`` is carved at its NATURAL width instead of being
             widened.  Chosen when shoreline fidelity matters more than a
             uniform two rows: widening moved 1,702 ha of Tokyo Bay coastline
             into the water, which is most of the mesh-vs-OSM gap that the
             coastline fit cannot close (see docs/coast_fit.md).
``allow``    one row permitted AND the keep bar lowered, so many more narrow
             channels are kept.  Measured on Tokyo Bay: 51 kept channels
             instead of 29, one-wide cells 14 -> 92, and four QA failures --
             the extra geometry, not the single rows, is what fails.

Environment overrides the recipe.  Canal and through records retain the legacy
row baseline under every setting.
"""
from __future__ import annotations

import os

#: Every accepted value, in order of how much they depart from the certified chain.
MODES = ("forbid", "natural", "allow")

#: Record kinds a single row is the honest representation of.  A through route
#: or a canal is a corridor the flow has to negotiate; a port basin or a dead
#: end simply is one cell wide.
ONE_ROW_KINDS = ("port", "dead-end")


def parse_one_wide(value="forbid"):
    if not isinstance(value, str) or value not in MODES:
        raise ValueError(f"one_wide must be one of {MODES}")
    return value


def configured_one_wide(recipe=None, *, environ=None):
    env = os.environ if environ is None else environ
    return parse_one_wide(env.get("SR_ONE_WIDE", (recipe or {}).get("one_wide", "forbid")))


def relaxes_selection(one_wide) -> bool:
    """Does this setting change WHICH channels are kept?

    Only ``allow`` does.  ``natural`` keeps the ``forbid`` selection, so the
    normalize pass and the detection width bars stay exactly as certified.
    """
    return parse_one_wide(one_wide) == "allow"


def permits_one_row(one_wide, kind=None) -> bool:
    """Does this setting carve a kept record at one row / natural width?

    ``kind`` is a policy record kind; omit it to ask about the setting alone.
    """
    mode = parse_one_wide(one_wide)
    if mode == "forbid":
        return False
    return kind is None or kind in ONE_ROW_KINDS


def generation_options(one_wide, *, environ=None):
    """Resolve global sizing and the baseline for canal/through treatment."""
    env = os.environ if environ is None else environ
    mode = parse_one_wide(one_wide)
    # `natural` leaves the SIZING FIELD alone: it changes how wide the channel
    # is cut, not how many rows the field asks for.  Only `allow` drops the
    # field to a single row.
    loose = mode == "allow"
    return dict(
        feature_rows=1.0 if loose else float(env.get("SR_FS", 3.0)),
        min_rows=1 if loose else 2,
        widen_factor=float(env.get("SR_WIDEN_FACTOR", "0.875")),
        attain_bar_h=float(env.get("SR_ATTAIN_BAR", "1.5")),
        force_two_rows=False if loose else env.get("SR_FORCE2ROWS", "off") == "on",
    )


def finishing_one_wide(metadata_path, recipe=None, *, environ=None):
    """Reuse generation's policy; refuse an explicitly conflicting finish setting."""
    import json
    from pathlib import Path

    env = os.environ if environ is None else environ
    path = Path(metadata_path).resolve()
    requested = configured_one_wide(recipe, environ=env)
    if not path.exists():
        return requested  # legacy meshes have no sidecar
    recorded = parse_one_wide(json.loads(path.read_text())["one_wide"])
    explicit = "SR_ONE_WIDE" in env or "one_wide" in (recipe or {})
    if explicit and requested != recorded:
        raise ValueError("one_wide differs from generation metadata; regenerate with this policy")
    return recorded
