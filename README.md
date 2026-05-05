# fvcom-mesh-tools

> ⚠️ **Pre-alpha** — under active development. APIs and CLIs are unstable.

Python toolkit for FVCOM unstructured mesh generation, repair, and visualization.

`fvcom-mesh-tools` provides a unified Python interface for building
high-quality FVCOM-ready unstructured meshes (`fort.14`), with a focus on:

- **Open-boundary edge orthogonalization** — enforce edges perpendicular to open boundaries
- **River-channel connectivity repair** — fix narrow channels where flow does not connect
- **Mesh quality inspection and visualization** — element quality, boundary classification, fort.14 plots

The package wraps several mature mesh tools behind a common backend interface
rather than reimplementing meshing algorithms from scratch.

## Backend strategy

| Role | Backend | License | How used |
|------|---------|---------|----------|
| Mesh generation (preferred) | [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) | CC0-1.0 | imported |
| Orthogonalization / smoothing | [MeshKernelPy](https://github.com/Deltares/MeshKernelPy) | MIT | imported |
| Grid utilities | [stompy](https://github.com/rustychris/stompy) | MIT | imported |
| Geometry-aware sizing (optional) | [JIGSAW](https://github.com/dengwirda/jigsaw-python) | custom (non-OSI) | imported, optional extra |
| External mesher | [gmsh](https://gmsh.info/) | GPL-2.0+ | invoked as subprocess |

GPL-licensed backends are deliberately invoked as subprocesses to keep this
project's distribution under Apache-2.0. See `THIRD_PARTY_NOTICES.md`.

## Installation

Python ≥3.12. conda-forge is the recommended channel for the scientific stack:

```bash
mamba env create -f environment.yml
mamba activate fvcom-mesh
pip install -e .
```

Optional backends are pulled in via extras:

```bash
pip install -e ".[ocsmesh,meshkernel,test]"
```

## Development

```bash
make install   # editable install with dev deps
make test      # pytest -q
make lint      # ruff check
```

## Status and scope

This project covers the FVCOM-specific mesh tooling that was originally
considered for [`oceanmesh-tools`](https://github.com/estuarine-utokyo/oceanmesh-tools).
`oceanmesh-tools` continues to focus on OceanMesh2D (MATLAB) input scanning
and `fort.14` post-processing.

## License

Apache License 2.0. See `LICENSE` and `NOTICE`. Third-party backend licenses
are documented in `THIRD_PARTY_NOTICES.md`.
