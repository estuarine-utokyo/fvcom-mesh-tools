"""Shadow package for tools/ci_local/run_ci_local.sh.

It makes ``import oceanmesh`` fail as it does on GitHub's CI.
"""

raise ModuleNotFoundError(
    "No module named 'oceanmesh' (hidden by tools/ci_local to reproduce GitHub CI)",
    name="oceanmesh")
