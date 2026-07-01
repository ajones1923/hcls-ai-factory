"""
System / cost-governance endpoint (PRD §5.1.3 NFR-COST). Exposes the per-tier token +
USD cost ledger the model router maintains, so a run's spend is observable and the cap
is auditable. Offline runs show counted tokens at zero cost (never billed).
"""
from __future__ import annotations

from fastapi import APIRouter

from src.utils.cost import get_ledger

router = APIRouter(prefix="/system", tags=["system"])


@router.get("/usage")
def usage():
    return get_ledger().totals()
