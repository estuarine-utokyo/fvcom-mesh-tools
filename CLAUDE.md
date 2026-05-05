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

## License policy (do not violate)

- The project is **Apache-2.0**.
- **Do not `import` GPL-licensed packages** (`oceanmesh`, `gmsh`) from
  `fvcom_mesh_tools`. They must be invoked as a subprocess, or kept in a
  separate plugin package with its own GPL license.
- JIGSAW (`jigsawpy`) carries a non-OSI license that restricts commercial
  distribution. Keep it as an optional extra; do not bundle it.
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
