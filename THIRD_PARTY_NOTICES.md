# Third-Party Software Notices

`fvcom-mesh-tools` is licensed under the Apache License 2.0 (see `LICENSE`).
This document records the licensing terms of third-party software that
`fvcom-mesh-tools` may invoke or import as a backend.

The project intentionally avoids importing GPL-licensed code into the
`fvcom-mesh-tools` distribution. GPL-licensed backends, when used, must be
called as a subprocess or maintained as a separate optional plugin package
under its own GPL license.

## Permissive backends (safe to `import`)

| Backend | License | Notes |
|---------|---------|-------|
| [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) | CC0-1.0 (public domain dedication) | NOAA Coastal Survey; ADCIRC/SCHISM-oriented Python mesh generator; preferred backend |
| [MeshKernelPy](https://github.com/Deltares/MeshKernelPy) | MIT | Deltares; unstructured-grid orthogonalization and smoothing |
| [stompy](https://github.com/rustychris/stompy) | MIT | UnstructuredGrid utilities. Not on PyPI; install from git |
| [PyFVCOM](https://github.com/pwcazenave/PyFVCOM) | MIT | FVCOM postprocessing helpers |

## Copyleft backends (do NOT `import` from `fvcom_mesh_tools`)

| Backend | License | Handling |
|---------|---------|----------|
| [gmsh](https://gmsh.info/) | GPL-2.0-or-later (no Python API exception in `LICENSE.txt`) | Invoke as a subprocess; never `import gmsh` from this package |
| [oceanmesh](https://github.com/CHLNDDEV/oceanmesh) (Python) | GPL-3.0-or-later | Invoke as a subprocess, or use as a separately licensed plugin |

## Restricted-license backends

| Backend | License | Handling |
|---------|---------|----------|
| [JIGSAW / jigsawpy](https://github.com/dengwirda/jigsaw-python) | Custom (non-OSI; commercial distribution by arrangement only) | Optional extra; do not bundle. Document in user-facing install instructions |

## Installation hints

PyPI extras (for permissive and copyleft-when-used-by-subprocess backends):

```bash
pip install "fvcom-mesh-tools[ocsmesh,meshkernel]"
```

conda-forge is the preferred channel for the scientific stack:

```bash
mamba install -c conda-forge ocsmesh meshkernel
```

stompy and jigsawpy:

```bash
pip install git+https://github.com/rustychris/stompy.git
pip install jigsawpy   # NOTE: review the JIGSAW license terms first
```

If you redistribute `fvcom-mesh-tools` together with any of these backends,
you must comply with each backend's license terms in addition to Apache-2.0.
