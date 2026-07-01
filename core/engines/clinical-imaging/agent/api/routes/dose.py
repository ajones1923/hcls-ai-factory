"""Dose intelligence endpoints.

Provides radiation dose tracking, cumulative dose monitoring,
DRL comparison, and population-level dose analytics.
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Any, Dict, List, Optional

from src.dose_intelligence import DoseIntelligenceEngine, DoseRecord

router = APIRouter(prefix="/dose", tags=["Dose Intelligence"])

# Singleton engine instance
_engine = DoseIntelligenceEngine()


# ── Request / Response Models ──────────────────────────────────────────


class DoseRecordRequest(BaseModel):
    """Request body for recording a dose."""

    patient_id: str
    study_date: str
    modality: str
    protocol: str
    body_region: str
    effective_dose_msv: float
    dlp_mgy_cm: Optional[float] = None
    ctdi_vol_mgy: Optional[float] = None
    dap_mgy_cm2: Optional[float] = None
    num_exposures: int = 1
    scanner_model: Optional[str] = None
    pediatric: bool = False


class CumulativeDoseResponse(BaseModel):
    """Cumulative dose summary response."""

    patient_id: str
    total_effective_dose_msv: float
    study_count: int
    date_range: Dict[str, str] = {}
    by_modality: Dict[str, float] = {}
    by_body_region: Dict[str, float] = {}
    alert_level: str = "normal"
    alert_message: Optional[str] = None


class DRLComparisonRequest(BaseModel):
    """Request body for DRL comparison."""

    patient_id: str = "anonymous"
    study_date: str = "2026-01-01"
    modality: str = "CT"
    protocol: str
    body_region: str = "chest"
    effective_dose_msv: float
    pediatric: bool = False


class DRLComparisonResponse(BaseModel):
    """DRL comparison result response."""

    protocol: str
    patient_dose_msv: float
    drl_msv: float
    achievable_dose_msv: float
    ratio: float
    status: str
    optimization_suggestions: List[str] = []


# ── Endpoints ──────────────────────────────────────────────────────────


@router.post("/record")
def record_dose(req: DoseRecordRequest):
    """Record a radiation dose for a patient study.

    Stores the dose record in the in-memory registry for cumulative
    tracking and population analytics.
    """
    record = DoseRecord(**req.model_dump())
    _engine.record_dose(record)
    return {
        "status": "recorded",
        "patient_id": req.patient_id,
        "protocol": req.protocol,
        "effective_dose_msv": req.effective_dose_msv,
        "registry_size": _engine.registry_size,
    }


@router.get("/cumulative/{patient_id}", response_model=CumulativeDoseResponse)
def get_cumulative_dose(patient_id: str, period_days: int = 365):
    """Get cumulative radiation dose for a patient.

    Parameters
    ----------
    patient_id : str
        Patient identifier.
    period_days : int
        Lookback period in days (default 365).

    Returns
    -------
    CumulativeDoseResponse
        Cumulative dose with alert level and breakdown.
    """
    result = _engine.get_cumulative_dose(patient_id, period_days=period_days)
    return CumulativeDoseResponse(**result.model_dump())


@router.post("/compare-drl", response_model=DRLComparisonResponse)
def compare_to_drl(req: DRLComparisonRequest):
    """Compare a study dose to national Diagnostic Reference Levels.

    Returns the comparison status (below_achievable, below_drl,
    above_drl, significantly_above) with optimization suggestions.
    """
    record = DoseRecord(
        patient_id=req.patient_id,
        study_date=req.study_date,
        modality=req.modality,
        protocol=req.protocol,
        body_region=req.body_region,
        effective_dose_msv=req.effective_dose_msv,
        pediatric=req.pediatric,
    )
    result = _engine.compare_to_drl(record)
    return DRLComparisonResponse(**result.model_dump())


@router.get("/population")
def get_population_summary():
    """Get population-level dose statistics across all recorded studies.

    Returns aggregate statistics including total records, unique patients,
    per-protocol and per-modality breakdowns, and alert distribution.
    If the registry is empty, demo data is generated automatically.
    """
    if _engine.registry_size == 0:
        _engine.generate_demo_data(n_records=200)
    return _engine.get_population_dose_summary()
