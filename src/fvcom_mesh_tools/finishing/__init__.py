"""Stage-2 automated finishing (docs/STAGE2_DESIGN.md)."""

from .detect import detect_violations
from .plan import plan_patches, write_ledger

__all__ = ["detect_violations", "plan_patches", "write_ledger"]
