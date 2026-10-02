# Third-Party Software Notices

`fvcom-mesh-tools` is licensed under the GNU General Public License,
version 3 or (at your option) any later version (GPL-3.0-or-later; see
`LICENSE` and `NOTICE`). Until 2026-10-02 it was Apache-2.0; its sole
copyright holder relicensed it because the package imports the
GPL-3.0-or-later `oceanmesh`. This document records the licensing terms of
third-party software that `fvcom-mesh-tools` may invoke or import as a
backend, and how each combines with the GPL.

Any backend whose license is compatible with GPL-3.0 may be imported. The
compatible licenses in use are permissive ones (CC0, MIT), GPL-3.0 itself
and GPL-2.0-or-later. A component whose terms add restrictions the GPL does
not allow (JIGSAW's core) is never imported or bundled.

## Backends imported as Python modules

| Backend | License | Notes |
|---------|---------|-------|
| [oceanmesh](https://github.com/CHLNDDEV/oceanmesh) (the laboratory's fork: https://github.com/estuarine-utokyo/oceanmesh) | GPL-3.0-or-later | The DistMesh generator, smoothing and patch remeshing. Default `fmesh-buildmesh --engine`. Same license as this package. |
| [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) | CC0-1.0 (public domain dedication) | NOAA Coastal Survey; used by `--engine ocsmesh` (gmsh-driven generation) and by `fmesh-mesh-combine` (`ops.combine_mesh`). |
| [MeshKernelPy](https://github.com/Deltares/MeshKernelPy) | MIT | Deltares; orthogonalization and smoothing. Optional. |
| [stompy](https://github.com/rustychris/stompy) | MIT | UnstructuredGrid utilities. Not on PyPI; install from git. |
| [PyFVCOM](https://github.com/pwcazenave/PyFVCOM) | MIT | FVCOM postprocessing helpers. |

## Backends invoked as external tools, or not used directly

| Backend | License | Handling |
|---------|---------|----------|
| [gmsh](https://gmsh.info/) | GPL-2.0-or-later | OCSMesh's `MeshDriver(engine="gmsh")` runs gmsh; compatible with this package's GPL-3.0-or-later in any case. |
| [JIGSAW / jigsawpy](https://github.com/dengwirda/jigsaw-python) | LGPL-3.0 (`jigsawpy` wrappers); the JIGSAW core C++ has its own license that restricts commercial distribution | Pulled in transitively as an OCSMesh dependency on conda-forge. Never imported by this package and never bundled with it: its core's terms are not GPL-compatible. |

## Installation hints

Most dependencies live on conda-forge; `oceanmesh` (the laboratory's fork)
is the exception. Use the project's `environment.yml` and
`docs/USER_GUIDE.md` §2:

```bash
mamba env create -f environment.yml
mamba activate fvcom-mesh-tools
pip install --no-deps -e <path to the oceanmesh fork>
pip install --no-deps -e .             # this package, editable
```

## Redistribution summary

| You ship | Required to comply with |
|----------|------------------------|
| `fvcom-mesh-tools`, alone or with `oceanmesh`, OCSMesh, gmsh, MeshKernelPy, stompy, PyFVCOM | GPL-3.0-or-later (with the permissive components' attribution notices) |
| anything that also bundles the JIGSAW core | not permitted under the GPL; ship JIGSAW separately, if at all |
