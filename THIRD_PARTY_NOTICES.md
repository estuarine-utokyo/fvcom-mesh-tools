# Third-Party Software Notices

`fvcom-mesh-tools` itself is licensed under the Apache License 2.0
(see `LICENSE`). This document records the licensing terms of
third-party software that `fvcom-mesh-tools` may invoke or import as
a backend, and the obligations that arise when redistributing
combined works.

The toolkit's design splits backends into three classes:

* **Permissive backends** can be `import`ed freely; they impose no
  obligations on downstream redistribution beyond their own
  attribution requirements.
* **Copyleft backends imported as Python modules** make any
  redistributed combined work subject to the backend's copyleft
  license. We currently use one such backend (``oceanmesh``,
  GPL-3.0-or-later) by deliberate choice, traded against its much
  better mesh quality.
* **Copyleft backends invoked as external tools** (subprocess /
  shell-out, no symbol linkage) do **not** propagate copyleft into
  the combined work.

## Permissive backends (safe to `import`)

| Backend | License | Notes |
|---------|---------|-------|
| [OCSMesh](https://github.com/noaa-ocs-modeling/OCSMesh) | CC0-1.0 (public domain dedication) | NOAA Coastal Survey; ADCIRC/SCHISM-oriented Python mesh generator. Used both by `--engine ocsmesh` (gmsh-driven generation) and by `fmesh-mesh-combine` (overlap / neighbor strategies via `ops.combine_mesh`). |
| [MeshKernelPy](https://github.com/Deltares/MeshKernelPy) | MIT | Deltares; unstructured-grid orthogonalization and smoothing. Optional; not required by the default pipeline. |
| [stompy](https://github.com/rustychris/stompy) | MIT | UnstructuredGrid utilities. Not on PyPI; install from git. |
| [PyFVCOM](https://github.com/pwcazenave/PyFVCOM) | MIT | FVCOM postprocessing helpers. |
| [JIGSAW / jigsawpy](https://github.com/dengwirda/jigsaw-python) | LGPL-3.0 (`jigsawpy` Python wrappers); JIGSAW core C++ has its own custom license — review for redistribution | Pulled in transitively as an OCSMesh dependency on conda-forge. We do not import it directly. |

## Copyleft backends imported as Python modules

| Backend | License | Handling |
|---------|---------|----------|
| [oceanmesh](https://github.com/CHLNDDEV/oceanmesh) (Python; the OceanMesh2D port) | GPL-3.0-or-later | Imported by `src/fvcom_mesh_tools/mesh_engine/oceanmesh.py` and used as the default `fmesh-buildmesh --engine`. **Distributing this toolkit together with `oceanmesh` makes the combined work GPL-3.0-or-later.** Source distributions of `fvcom-mesh-tools` itself remain Apache-2.0. The `--engine ocsmesh` path is fully GPL-free at link time. |

## Copyleft backends invoked as external tools

| Backend | License | Handling |
|---------|---------|----------|
| [gmsh](https://gmsh.info/) | GPL-2.0-or-later | OCSMesh's `MeshDriver(engine="gmsh")` shells out to gmsh; we do not link `libgmsh` directly. The combined work is therefore not propagating gmsh's GPL into our Python code. Downstream users who *do* link gmsh (e.g. by using `MeshDriver(engine="jigsaw")`-vs.-gmsh in their own code) take on the GPL obligations. |

## Installation hints

Most dependencies live on conda-forge; `oceanmesh` is the exception.
Use the project's `environment.yml`:

```bash
mamba env create -f environment.yml
mamba activate fvcom-mesh
pip install --no-deps oceanmesh        # GPL-3.0-or-later (PyPI)
pip install --no-deps -e .             # this package, editable
```

If you must avoid the GPL footprint entirely, omit the
`pip install oceanmesh` step and run with `fmesh-buildmesh
--engine ocsmesh`. The default `--engine oceanmesh` will fail with an
`ImportError` in that case.

## Combined-work obligations summary

| You ship | Required to comply with |
|----------|------------------------|
| `fvcom-mesh-tools` source only | Apache-2.0 |
| `fvcom-mesh-tools` + OCSMesh + gmsh | Apache-2.0 + CC0 + GPL-2.0+ runtime obligations |
| `fvcom-mesh-tools` + `oceanmesh` | GPL-3.0-or-later on the combined work |

If in doubt, prefer the OCSMesh path, or ship `oceanmesh` as a
separately-installed tool that the user wires up themselves.
