"""Backend protocol for mesh-generation/operation engines (OCSMesh, MeshKernel, gmsh, ...)."""

from typing import Protocol, runtime_checkable


@runtime_checkable
class Backend(Protocol):
    name: str

    def is_available(self) -> bool:
        ...
