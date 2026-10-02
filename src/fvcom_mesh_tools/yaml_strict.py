"""YAML for recipes: a key given twice is an error, not a silent override.

``yaml.safe_load`` keeps the last of two equal keys, so ``rfactor: 0.2``
followed by ``rfactor: 0.8`` became 0.8 with nothing said (review of the
extend tools, round 21 F7).
"""

from __future__ import annotations

import yaml

__all__ = ["load_unique"]


class _UniqueLoader(yaml.SafeLoader):
    pass


def _mapping(loader, node):
    out = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=True)
        if key in out:
            mark = key_node.start_mark
            raise ValueError(f"duplicate key {key!r} at line {mark.line + 1}, "
                             f"column {mark.column + 1}")
        out[key] = loader.construct_object(value_node, deep=True)
    return out


_UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def load_unique(text: str):
    """``yaml.safe_load`` that refuses a key repeated in any mapping."""
    return yaml.load(text, Loader=_UniqueLoader)  # noqa: S506 -- a SafeLoader
