"""fvcom-mesh-tools: utilities for FVCOM unstructured mesh generation, repair, and visualization."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("fvcom-mesh-tools")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
