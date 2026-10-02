"""YAML for recipes: a key given twice is an error, not a silent override.

``yaml.safe_load`` keeps the last of two equal keys, so ``rfactor: 0.2``
followed by ``rfactor: 0.8`` became 0.8 with nothing said (review of the
extend tools, round 21 F7).
"""

from __future__ import annotations

import yaml

__all__ = ["load_unique"]


_MERGE = "tag:yaml.org,2002:merge"


class _UniqueLoader(yaml.SafeLoader):
    """SafeLoader that checks the whole composed document first.

    Every mapping in the node graph is checked once, before any merge key is
    expanded (expansion rewrites aliased mappings in place): a key written
    twice, or two merge keys in one mapping, is an error; an explicit key
    overriding a merged one is not (review rounds 21 F7 - 24 F2).
    """

    def construct_document(self, node):
        _check_graph(self, node)
        return super().construct_document(node)


def _check_graph(loader, root) -> None:
    visited, stack = set(), [root]
    while stack:
        node = stack.pop()
        if id(node) in visited:
            continue
        visited.add(id(node))
        if isinstance(node, yaml.SequenceNode):
            stack.extend(node.value)
        elif isinstance(node, yaml.MappingNode):
            seen, merges = set(), 0
            for key_node, value_node in node.value:
                stack.append(value_node)
                mark = key_node.start_mark
                where = f"line {mark.line + 1}, column {mark.column + 1}"
                if key_node.tag == _MERGE:
                    merges += 1
                    if merges > 1:
                        raise ValueError(f"a second merge key '<<' at {where}")
                    continue
                stack.append(key_node)
                # the plain key '=' carries the value tag, which SafeLoader
                # reads as the string '=' (round 25 F6)
                if key_node.tag == "tag:yaml.org,2002:value":
                    key = "="
                else:
                    key = loader.construct_object(key_node, deep=True)
                if key in seen:
                    raise ValueError(f"duplicate key {key!r} at {where}")
                seen.add(key)


def load_unique(text: str):
    """``yaml.safe_load`` that refuses a key repeated in any mapping."""
    return yaml.load(text, Loader=_UniqueLoader)  # noqa: S506 -- a SafeLoader
