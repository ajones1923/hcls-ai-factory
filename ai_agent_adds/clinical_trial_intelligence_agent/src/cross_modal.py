"""Cross-agent integration for the Clinical Trial Intelligence Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into a unified trial assessment.

Supported cross-agent queries:
  - query_oncology_agent()   -- molecular matches for precision oncology
  - query_pgx_agent()        -- pharmacogenomic metabolism screening
  - query_cardiology_agent() -- cardiac safety assessment
  - query_biomarker_agent()  -- biomarker enrichment strategies
  - integrate_cross_agent_results() -- unified assessment

All functions degrade gracefully: if an agent is unavailable, a warning
is logged and a default response is returned.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

from config.settings import settings

logger = logging.getLogger(__name__)


# ===================================================================
# CROSS-AGENT QUERY FUNCTIONS
# ===================================================================


def query_oncology_agent(
    patient_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Oncology Intelligence Agent for molecular trial matches.

    Sends patient biomarker and genomic data to the oncology agent to
    identify precision medicine trials targeting the patient's molecular
    profile.

    Args:
        patient_profile: Patient data including diagnosis, biomarkers,
            and genomic variants.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``molecular_matches``, and ``recommendations``.
        Returns a default response if the agent is unavailable.
    """
    try:
        import requests

        biomarkers = patient_profile.get("biomarkers", [])
        genomic_variants = patient_profile.get("genomic_variants", [])
        diagnosis = patient_profile.get("diagnosis", "")

        response = requests.post(
            f"{settings.ONCOLOGY_AGENT_URL}/api/query",
            json={
                "question": f"Find molecular-targeted trials for {diagnosis}",
                "patient_context": {
                    "biomarkers": biomarkers,
                    "genomic_variants": genomic_variants,
                    "diagnosis": diagnosis,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "oncology",
            "molecular_matches": data.get("matches", []),
            "recommendations": data.get("recommendations", []),
            "confidence": data.get("confidence", 0.0),
        }

    except ImportError:
        logger.warning("requests library not available for oncology agent query")
        return _unavailable_response("oncology")
    except Exception as exc:
        logger.warning("Oncology agent query failed: %s", exc)
        return _unavailable_response("oncology")


def query_pgx_agent(
    patient_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Pharmacogenomics Intelligence Agent for metabolism screening.

    Checks the patient's pharmacogenomic profile to identify potential
    drug-gene interactions that could affect trial drug metabolism.

    Args:
        patient_profile: Patient data including genomic variants and
            current medications.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``pgx_results``, and ``warnings``.
    """
    try:
        import requests

        genomic_variants = patient_profile.get("genomic_variants", [])
        medications = patient_profile.get("medications", [])

        response = requests.post(
            f"{settings.PGX_AGENT_URL}/api/query",
            json={
                "question": "Screen for pharmacogenomic interactions",
                "patient_context": {
                    "genomic_variants": genomic_variants,
                    "medications": medications,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "pharmacogenomics",
            "pgx_results": data.get("pgx_results", []),
            "warnings": data.get("warnings", []),
            "metabolizer_status": data.get("metabolizer_status", {}),
        }

    except ImportError:
        logger.warning("requests library not available for PGx agent query")
        return _unavailable_response("pharmacogenomics")
    except Exception as exc:
        logger.warning("PGx agent query failed: %s", exc)
        return _unavailable_response("pharmacogenomics")


def query_cardiology_agent(
    drug: str,
    patient_profile: Optional[Dict[str, Any]] = None,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Cardiology Intelligence Agent for cardiac safety assessment.

    Checks whether a trial drug has known cardiovascular safety concerns
    (QT prolongation, cardiomyopathy risk, etc.).

    Args:
        drug: Drug name to assess.
        patient_profile: Optional patient cardiac history.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cardiac_safety``, and ``risk_flags``.
    """
    try:
        import requests

        query_data: Dict[str, Any] = {
            "question": f"Assess cardiac safety profile of {drug}",
        }
        if patient_profile:
            query_data["patient_context"] = patient_profile

        response = requests.post(
            f"{settings.CARDIOLOGY_AGENT_URL}/api/query",
            json=query_data,
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "cardiology",
            "cardiac_safety": data.get("safety_assessment", {}),
            "risk_flags": data.get("risk_flags", []),
            "qt_risk": data.get("qt_prolongation_risk", "unknown"),
        }

    except ImportError:
        logger.warning("requests library not available for cardiology agent query")
        return _unavailable_response("cardiology")
    except Exception as exc:
        logger.warning("Cardiology agent query failed: %s", exc)
        return _unavailable_response("cardiology")


def query_biomarker_agent(
    biomarkers: List[str],
    therapeutic_area: str = "",
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Biomarker Intelligence Agent for enrichment strategies.

    Identifies biomarker-driven enrichment strategies for clinical trial
    design and patient stratification.

    Args:
        biomarkers: List of biomarker names/results.
        therapeutic_area: Target therapeutic area.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``enrichment_strategies``, and
        ``stratification_recommendations``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.BIOMARKER_AGENT_URL}/api/query",
            json={
                "question": f"Biomarker enrichment strategies for {therapeutic_area}",
                "biomarkers": biomarkers,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "biomarker",
            "enrichment_strategies": data.get("strategies", []),
            "stratification_recommendations": data.get("stratification", []),
        }

    except ImportError:
        logger.warning("requests library not available for biomarker agent query")
        return _unavailable_response("biomarker")
    except Exception as exc:
        logger.warning("Biomarker agent query failed: %s", exc)
        return _unavailable_response("biomarker")


# ===================================================================
# INTEGRATION FUNCTION
# ===================================================================


def integrate_cross_agent_results(
    results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Integrate results from multiple cross-agent queries into a unified assessment.

    Combines molecular matches, PGx warnings, cardiac safety flags, and
    biomarker strategies into a single assessment suitable for trial
    eligibility and safety determination.

    Args:
        results: List of cross-agent result dicts (from the query_* functions).

    Returns:
        Unified assessment dict with:
          - ``agents_consulted``: List of agent names queried.
          - ``agents_available``: List of agents that responded successfully.
          - ``combined_warnings``: Aggregated warnings from all agents.
          - ``combined_recommendations``: Aggregated recommendations.
          - ``safety_flags``: Combined safety concerns.
          - ``overall_assessment``: Summary assessment text.
    """
    agents_consulted: List[str] = []
    agents_available: List[str] = []
    combined_warnings: List[str] = []
    combined_recommendations: List[str] = []
    safety_flags: List[str] = []

    for result in results:
        agent = result.get("agent", "unknown")
        agents_consulted.append(agent)

        if result.get("status") == "success":
            agents_available.append(agent)

            # Collect warnings
            warnings = result.get("warnings", [])
            combined_warnings.extend(
                f"[{agent}] {w}" for w in warnings
            )

            # Collect recommendations
            recs = result.get("recommendations", [])
            combined_recommendations.extend(
                f"[{agent}] {r}" for r in recs
            )

            # Collect safety flags
            risk_flags = result.get("risk_flags", [])
            safety_flags.extend(
                f"[{agent}] {f}" for f in risk_flags
            )

            # PGx-specific warnings
            pgx_results = result.get("pgx_results", [])
            for pgx in pgx_results:
                if isinstance(pgx, dict) and pgx.get("impact", "").lower() in ("high", "critical"):
                    combined_warnings.append(
                        f"[pharmacogenomics] {pgx.get('gene', '')}: {pgx.get('recommendation', '')}"
                    )

    # Generate overall assessment
    if not agents_available:
        overall = "No cross-agent data available. Proceeding with trial agent data only."
    elif safety_flags:
        overall = (
            f"Cross-agent analysis identified {len(safety_flags)} safety concern(s). "
            f"Review recommended before proceeding."
        )
    elif combined_warnings:
        overall = (
            f"Cross-agent analysis completed with {len(combined_warnings)} warning(s). "
            f"All flagged items should be reviewed."
        )
    else:
        overall = (
            f"Cross-agent analysis completed successfully. "
            f"{len(agents_available)} agent(s) consulted with no safety concerns."
        )

    return {
        "agents_consulted": agents_consulted,
        "agents_available": agents_available,
        "combined_warnings": combined_warnings,
        "combined_recommendations": combined_recommendations,
        "safety_flags": safety_flags,
        "overall_assessment": overall,
    }


# ===================================================================
# HELPERS
# ===================================================================


def _unavailable_response(agent_name: str) -> Dict[str, Any]:
    """Return a standard unavailable response for a cross-agent query.

    Args:
        agent_name: Name of the unavailable agent.

    Returns:
        Dict with ``status`` set to ``"unavailable"``.
    """
    return {
        "status": "unavailable",
        "agent": agent_name,
        "message": f"{agent_name} agent is not currently available",
        "recommendations": [],
        "warnings": [],
    }
