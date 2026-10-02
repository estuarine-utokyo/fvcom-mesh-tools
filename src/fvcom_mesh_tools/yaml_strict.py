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


_MERGE = "tag:yaml.org,2002:merge"


def _check_keys(loader, node, visited) -> None:
    """Refuse a key written twice in ``node``, and in every mapping merged
    into it (inline or by alias), which flattening never hands to the
    mapping constructor (review rounds 22 F5, 23 F1)."""
    if id(node) in visited:
        return
    visited.add(id(node))
    seen = set()
    for key_node, value_node in node.value:
        if key_node.tag == _MERGE:
            parts = (value_node.value if isinstance(value_node, yaml.SequenceNode)
                     else [value_node])
            for part in parts:
                if isinstance(part, yaml.MappingNode):
                    _check_keys(loader, part, visited)
            continue
        key = loader.construct_object(key_node, deep=True)
        if key in seen:
            mark = key_node.start_mark
            raise ValueError(f"duplicate key {key!r} at line {mark.line + 1}, "
                             f"column {mark.column + 1}")
        seen.add(key)


def _mapping(loader, node):
    # Duplicates are looked for among the keys written in this mapping and
    # in what it merges; merge keys (<<) are then expanded as SafeLoader
    # does, an explicit key overriding a merged one.
    _check_keys(loader, node, set())
    loader.flatten_mapping(node)
    return loader.construct_mapping(node, deep=True)


_UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _mapping)


def load_unique(text: str):
    """``yaml.safe_load`` that refuses a key repeated in any mapping."""
    return yaml.load(text, Loader=_UniqueLoader)  # noqa: S506 -- a SafeLoader
