"""Stage-2 automated finishing (docs/STAGE2_DESIGN.md)."""

from .detect import detect_violations
from .execute import execute_patches
from .plan import plan_patches, write_ledger

__all__ = ["detect_violations", "execute_patches", "plan_patches", "write_ledger"]
