"""Tests for ``fvcom_mesh_tools.cli.buildmesh`` argparse-level behaviour.

These tests intentionally do not invoke a mesher — they only exercise
the CLI front matter (argparse, deprecation warnings, validation
errors) up to the point a real DEM would be needed. Running an
end-to-end build is the job of the PoC notebooks.
"""
from __future__ import annotations

import warnings

import pytest

from fvcom_mesh_tools.cli import buildmesh


def test_engine_ocsmesh_emits_deprecation_warning(capsys) -> None:
    """``--engine ocsmesh`` must fire a DeprecationWarning AND a stderr
    notice. The DEM does not exist, so the call returns 2 (validation
    failure); the warning fires before the validation, which is what
    we want — even a misuse with a missing DEM still surfaces the
    deprecation to the user."""
    argv = ["nonexistent.nc", "out.14", "--engine", "ocsmesh"]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        rc = buildmesh.main(argv)

    msgs = [str(w.message) for w in caught
            if issubclass(w.category, DeprecationWarning)]
    assert any("--engine ocsmesh" in m and "DEPRECATED" in m for m in msgs), (
        f"expected DeprecationWarning mentioning '--engine ocsmesh' and "
        f"'DEPRECATED', got: {msgs!r}"
    )

    err = capsys.readouterr().err
    assert "DEPRECATED" in err
    assert "engine_complementarity" in err  # cross-link is present

    # DEM does not exist -> validation error 2; we only care that the
    # deprecation fired before that.
    assert rc == 2


def test_engine_oceanmesh_does_not_emit_deprecation(capsys) -> None:
    """``--engine oceanmesh`` (the default) must not emit a
    DeprecationWarning."""
    argv = ["nonexistent.nc", "out.14", "--engine", "oceanmesh"]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        rc = buildmesh.main(argv)

    deprecations = [w for w in caught
                    if issubclass(w.category, DeprecationWarning)]
    assert not deprecations, (
        f"oceanmesh path should not deprecate; got {deprecations!r}"
    )

    err = capsys.readouterr().err
    assert "DEPRECATED" not in err
    assert rc == 2  # missing DEM


def test_engine_default_is_oceanmesh(capsys) -> None:
    """Passing no --engine should behave exactly as ``--engine oceanmesh``."""
    argv = ["nonexistent.nc", "out.14"]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", DeprecationWarning)
        rc = buildmesh.main(argv)

    deprecations = [w for w in caught
                    if issubclass(w.category, DeprecationWarning)]
    assert not deprecations
    assert capsys.readouterr().err  # missing-DEM message; nonzero
    assert rc == 2


def test_engine_invalid_value_rejected() -> None:
    """argparse choices guard the ``--engine`` value space."""
    argv = ["nonexistent.nc", "out.14", "--engine", "jigsaw"]
    with pytest.raises(SystemExit):
        buildmesh.main(argv)
