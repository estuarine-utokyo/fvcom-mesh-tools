"""Stage-2 automated finishing (docs/STAGE2_DESIGN.md)."""

from .detect import detect_violations
from .directives import apply_directives
from .execute import execute_patches
from .plan import plan_patches, write_ledger

__all__ = [
    "apply_directives",
    "detect_violations",
    "execute_patches",
    "plan_patches",
    "write_ledger",
]
