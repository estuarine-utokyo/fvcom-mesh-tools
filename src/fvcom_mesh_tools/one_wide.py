"""Single policy switch for width-based channel enforcement.

``forbid`` preserves the certified chain (including its connectivity exceptions).
``allow`` permits one row only for port and dead-end policy records. QA still applies.
Environment overrides the recipe. Canal and through records retain the legacy row baseline.
"""
from __future__ import annotations

import os


def parse_one_wide(value="forbid"):
    if not isinstance(value, str) or value not in ("forbid", "allow"):
        raise ValueError("one_wide must be 'forbid' or 'allow'")
    return value


def configured_one_wide(recipe=None, *, environ=None):
    env = os.environ if environ is None else environ
    return parse_one_wide(env.get("SR_ONE_WIDE", (recipe or {}).get("one_wide", "forbid")))


def generation_options(one_wide, *, environ=None):
    """Resolve global sizing and the baseline for canal/through treatment."""
    env = os.environ if environ is None else environ
    allow = parse_one_wide(one_wide) == "allow"
    return dict(
        feature_rows=1.0 if allow else float(env.get("SR_FS", 3.0)),
        min_rows=1 if allow else 2,
        widen_factor=float(env.get("SR_WIDEN_FACTOR", "0.875")),
        attain_bar_h=float(env.get("SR_ATTAIN_BAR", "1.5")),
        force_two_rows=False if allow else env.get("SR_FORCE2ROWS", "off") == "on",
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
