# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working in this repository.

## Project overview

`fvcom-mesh-tools` is a Python toolkit for FVCOM unstructured mesh generation,
repair, and visualization. It wraps multiple mesh backends (OCSMesh,
MeshKernelPy, stompy, gmsh, JIGSAW) behind a common interface rather than
reimplementing mesh algorithms from scratch. The repo is sibling to
`oceanmesh-tools`, which is scoped to OceanMesh2D (MATLAB) input scanning and
`fort.14` post-processing.

## Key commands

```bash
# Editable install with dev/test deps
pip install -e ".[dev]"

# Tests
pytest -q

# Lint
ruff check .
```

A `Makefile` exposes `install`, `test`, `lint`, `format`, `clean` targets.

## Before every commit: reproduce GitHub CI locally (owner rule, 2026-11-16)

A plain `pytest -q` in the laboratory environment is NOT enough. GitHub's CI has no oceanmesh fork
and no `WORK_DIR`/`DATA_DIR`; on 2026-11-16 twelve new tests that generate meshes with oceanmesh
passed locally (1360 passed) and failed on GitHub (CI run 43778047151) because they lacked
`@pytest.mark.needs_oceanmesh`. So, before committing any change under `src/`, `tests/`,
`notebooks/`, `tools/`, `pyproject.toml` or `environment.yml`:

1. Run `tools/ci_local/run_ci_local.sh` -- on GENKAI as a batch job, `pjsub jobs/genkai/ci_local.sh`
   from the repository root (~2 min; log in `logs/ci_local.<jobid>.log`). It runs every CI step:
   ruff (repository and notebooks 42*/43*), pytest in **CI mode** (oceanmesh hidden by
   `tools/ci_local/hide_oceanmesh/`, `WORK_DIR`/`DATA_DIR` unset), pytest in lab mode (full
   environment), `pre-commit run --all-files`, and the CLI `--help` smoke checks of `ci.yml`.
2. Commit only when the log ends with `=== ci_local: PASS`, and quote the job id in the commit
   message or the report. A failing step is fixed first, never committed "to fix later".
3. A new test that imports or runs oceanmesh (directly or through `mesh_engine.polygon`, the
   pipeline build stage, ...) carries `@pytest.mark.needs_oceanmesh`; CI mode is what proves it.
4. When `.github/workflows/ci.yml` gains a step, add it to `run_ci_local.sh` in the same commit.
5. After a push (owner only), check the run with `gh run list --limit 3`; a red run is reported.

## License policy (do not violate)

- The project is **GPL-3.0-or-later** (relicensed from Apache-2.0 on
  2026-10-02 by its sole copyright holder, because it imports `oceanmesh`).
- Import only packages whose licenses are compatible with GPL-3.0
  (permissive ones such as MIT/BSD/CC0, LGPL, GPL-3.0, GPL-2.0-or-later).
  `oceanmesh` (GPL-3.0-or-later) may be imported. Do not import or bundle
  code under GPL-2.0-only or with non-free/extra restrictions; the one
  exception is the existing optional OCSMesh paths below, and no new code
  may depend on them.
- Triangle (loaded by OCSMesh at import) and JIGSAW (`jigsawpy`) restrict
  commercial distribution and are not GPL-compatible. They are not in the
  default environment; the OCSMesh paths that need them stay optional, are
  imported lazily, and must never be bundled or redistributed with the
  package (owner, 2026-10-02; THIRD_PARTY_NOTICES.md).
- New source files need no license header; `LICENSE`, `NOTICE` and the
  `license` field in `pyproject.toml` state the license.
- See `THIRD_PARTY_NOTICES.md` for the canonical handling rules.

## Code organization

- `src/fvcom_mesh_tools/` — package source (src layout)
- `src/fvcom_mesh_tools/backend.py` — `Backend` protocol that backend wrappers
  implement
- `tests/` — pytest tests
- `tests/fixtures/` — small fort.14 / NetCDF fixtures (must be small; large
  data lives outside the repo)

## Conventions

- All file paths in code: absolute (`Path.resolve()`).
- fort.14 node IDs are 1-indexed.
- Documentation, code comments, and commit messages: English (per the user's
  global rule).
- Add tests for both happy paths and error cases.
- Avoid pre-designing abstractions for hypothetical features. Add Protocol
  methods only when a concrete implementation needs them.

## Environment

- Python 3.12+.
- Install scientific deps via `mamba install -c conda-forge ...`, not
  `pip install`, for compiled libraries (numpy, netCDF4, hdf5, geopandas,
  rasterio, etc.).
- ruff: `line-length = 100`, `target-version = py312`.

## Status

Pre-alpha. Public API and CLI are unstable. Do not tag releases until the
backend abstraction settles.
