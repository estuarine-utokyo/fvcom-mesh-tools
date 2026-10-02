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

## Backends imported as Python modules (GPL-compatible)

| Backend | License | Notes |
|---------|---------|-------|
| [oceanmesh](https://github.com/CHLNDDEV/oceanmesh) (the laboratory's fork: https://github.com/estuarine-utokyo/oceanmesh) | GPL-3.0-or-later | The DistMesh generator, smoothing and patch remeshing (`oceanmesh.remesh_patch`, the DistMesh one in `mesh_merge`). Default `fmesh-buildmesh --engine`. Same license as this package. |
| [MeshKernelPy](https://github.com/Deltares/MeshKernelPy) | MIT | Deltares; orthogonalization and smoothing. Optional. |
| [stompy](https://github.com/rustychris/stompy) | MIT | UnstructuredGrid utilities. Not on PyPI; install from git. |
| [PyFVCOM](https://github.com/pwcazenave/PyFVCOM) | MIT | FVCOM postprocessing helpers. |

## Optional backends that are not GPL-compatible (private use only)

These are **not** in the default environment (owner's decision,
2026-10-02). They serve only optional paths, which import them lazily. They
may be installed for private use (`mamba install -c conda-forge ocsmesh
jigsawpy triangle`), but must **never be redistributed together with this
package**.

| Backend | License | Used by | Why private only |
|---------|---------|---------|------------------|
| [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) | CC0-1.0 itself | `fmesh-buildmesh --engine ocsmesh` (deprecated), `fmesh-mesh-combine`, `fmesh-meshclean --repair-skewed-elements` | Importing it loads [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) at package initialisation (`ocsmesh.engines.triangle`). |
| [Triangle](https://www.cs.cmu.edu/~quake/triangle.html) (via the `triangle` Python wrapper) | Shewchuk's terms: no sale or inclusion in commercial products without permission (the wrapper's own LGPL does not lift them) | reached through OCSMesh; the fork's separate `oceanmesh/remesh_patch.py` jigsaw path, which this package does not call | Restricts commercial distribution; not GPL-compatible. |
| [JIGSAW / jigsawpy](https://github.com/dengwirda/jigsaw-python) | JIGSAW's own license (wrappers and core) restricting commercial use | OCSMesh's optional engine; the fork's separate `remesh_patch.py` jigsaw path | Restricts commercial distribution; not GPL-compatible. |

## Backends invoked as external tools

| Backend | License | Handling |
|---------|---------|----------|
| [gmsh](https://gmsh.info/) | GPL-2.0-or-later | OCSMesh's `MeshDriver(engine="gmsh")` runs gmsh; compatible with this package's GPL-3.0-or-later in any case. |

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
