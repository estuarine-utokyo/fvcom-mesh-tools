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
    # Duplicates are looked for among the keys written in this mapping; merge
    # keys (<<) are then expanded as SafeLoader does, an explicit key
    # overriding a merged one (review round 22 F5).
    seen = set()
    for key_node, _value in node.value:
        if key_node.tag == "tag:yaml.org,2002:merge":
            continue
        key = loader.construct_object(key_node, deep=True)
        if key in seen:
            mark = key_node.start_mark
            raise ValueError(f"duplicate key {key!r} at line {mark.line + 1}, "
                             f"column {mark.column + 1}")
        seen.add(key)
    loader.flatten_mapping(node)
    return loader.construct_mapping(node, deep=True)


_UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def load_unique(text: str):
    """``yaml.safe_load`` that refuses a key repeated in any mapping."""
    return yaml.load(text, Loader=_UniqueLoader)  # noqa: S506 -- a SafeLoader
