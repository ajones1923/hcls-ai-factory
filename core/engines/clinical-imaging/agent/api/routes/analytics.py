"""Analytics endpoints for population-level imaging intelligence.

GPU-accelerated via RAPIDS cuDF with pandas fallback.

Author: Adam Jones
Date: April 2026
"""

from typing import Any, Dict

from fastapi import APIRouter, HTTPException, Query
from loguru import logger

from src.analytics_engine import (
    CohortCriteria,
    CohortResult,
    ImagingAnalyticsEngine,
    PopulationSummary,
    TrendResult,
)

router = APIRouter(prefix="/api/analytics", tags=["Analytics"])

# ── Module-level engine (initialized on first import) ──
_engine: ImagingAnalyticsEngine | None = None


def _get_engine() -> ImagingAnalyticsEngine:
    """Lazy-initialize the analytics engine singleton."""
    global _engine
    if _engine is None:
        _engine = ImagingAnalyticsEngine()
    return _engine


# ═══════════════════════════════════════════════════════════════════════
# ENDPOINTS
# ═══════════════════════════════════════════════════════════════════════


@router.get("/population", response_model=PopulationSummary)
async def population_summary():
    """Get population-level imaging statistics.

    Returns aggregate counts, distributions, finding prevalence,
    and critical finding rates across all registered studies.
    """
    engine = _get_engine()
    try:
        return engine.population_summary()
    except Exception as e:
        logger.error(f"Population summary failed: {e}")
        raise HTTPException(status_code=500, detail="Analytics computation failed")


@router.post("/cohort", response_model=CohortResult)
async def cohort_query(criteria: CohortCriteria):
    """Find imaging studies matching criteria.

    Supports filtering by modality, body region, severity,
    finding type, date range, patient age, and finding size.
    """
    engine = _get_engine()
    try:
        return engine.cohort_query(criteria)
    except Exception as e:
        logger.error(f"Cohort query failed: {e}")
        raise HTTPException(status_code=500, detail="Cohort query failed")


@router.get("/trends/{metric}", response_model=TrendResult)
async def temporal_trends(
    metric: str,
    granularity: str = Query("monthly", regex="^(weekly|monthly|quarterly)$"),
):
    """Get temporal trends for a metric.

    Supported metrics: critical_findings, ct_volume, total_volume.
    Granularity: weekly, monthly, quarterly.
    """
    engine = _get_engine()
    try:
        return engine.temporal_trends(metric, granularity)
    except Exception as e:
        logger.error(f"Trend computation failed: {e}")
        raise HTTPException(status_code=500, detail="Trend computation failed")


@router.get("/correlation")
async def finding_correlation(
    finding_a: str = Query(..., min_length=1, description="First finding type"),
    finding_b: str = Query(..., min_length=1, description="Second finding type"),
) -> Dict[str, Any]:
    """Get co-occurrence rate between two finding types.

    Computes patient-level co-occurrence: percentage of patients
    with finding_a who also have finding_b.
    """
    engine = _get_engine()
    try:
        return engine.finding_correlation(finding_a, finding_b)
    except Exception as e:
        logger.error(f"Correlation computation failed: {e}")
        raise HTTPException(status_code=500, detail="Correlation computation failed")


@router.get("/severity-by-modality")
async def severity_by_modality() -> Dict[str, Dict[str, int]]:
    """Get severity cross-tabulated by modality.

    Returns a nested dict: {modality: {severity: count, ...}, ...}
    """
    engine = _get_engine()
    try:
        return engine.severity_by_modality()
    except Exception as e:
        logger.error(f"Severity-by-modality failed: {e}")
        raise HTTPException(status_code=500, detail="Cross-tabulation failed")


@router.get("/concordance")
async def ai_concordance() -> Dict[str, Any]:
    """Compare AI workflow results vs radiologist reports.

    Returns agreement rate and list of disagreements for studies
    with both AI and radiologist severity annotations.
    """
    engine = _get_engine()
    try:
        return engine.ai_concordance()
    except Exception as e:
        logger.error(f"Concordance computation failed: {e}")
        raise HTTPException(status_code=500, detail="Concordance computation failed")


@router.post("/generate-demo-data")
async def generate_demo(
    n_studies: int = Query(500, ge=1, le=50000, description="Number of synthetic studies"),
) -> Dict[str, Any]:
    """Generate synthetic demo data for analytics.

    Creates n_studies with realistic distributions of modality,
    body region, severity, and finding types.
    """
    engine = _get_engine()
    try:
        count = engine.generate_demo_data(n_studies=n_studies)
        return {
            "generated": count,
            "total_studies": engine.study_count,
            "backend": engine.backend_name,
        }
    except Exception as e:
        logger.error(f"Demo data generation failed: {e}")
        raise HTTPException(status_code=500, detail="Demo data generation failed")


@router.get("/status")
async def analytics_status() -> Dict[str, Any]:
    """Return analytics engine status."""
    engine = _get_engine()
    return {
        "backend": engine.backend_name,
        "rapids_available": engine._rapids_available,
        "registered_studies": engine.study_count,
    }
