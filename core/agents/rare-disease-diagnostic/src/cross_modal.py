"""Cross-agent integration for the Rare Disease Diagnostic Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into a unified rare disease assessment.

Supported cross-agent queries:
  - query_cardiology_agent()   -- cardiac genetics (channelopathies, cardiomyopathies)
  - query_biomarker_agent()    -- metabolic profiles for IEMs
  - query_pgx_agent()          -- post-diagnosis pharmacogenomic dosing
  - query_imaging_agent()      -- phenotype correlation with imaging data

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


def query_cardiology_agent(
    gene: str,
    phenotypes: List[str] = None,
    patient_profile: Optional[Dict[str, Any]] = None,
    timeout: float = None,
) -> Dict[str, Any]:
    """Query the Cardiology Intelligence Agent for cardiac genetics assessment.

    Checks for inherited cardiac conditions (channelopathies, cardiomyopathies)
    that overlap with rare disease differential diagnoses.

    Args:
        gene: Gene symbol to assess for cardiac involvement.
        phenotypes: List of HPO phenotype terms.
        patient_profile: Optional patient cardiac history.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cardiac_findings``, and ``recommendations``.
    """
    timeout = timeout or settings.CROSS_AGENT_TIMEOUT
    phenotypes = phenotypes or []

    try:
        import requests

        query_data: Dict[str, Any] = {
            "question": f"Assess cardiac genetics for gene {gene} with phenotypes: {', '.join(phenotypes[:5])}",
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
            "cardiac_findings": data.get("findings", []),
            "recommendations": data.get("recommendations", []),
            "cardiac_genes": data.get("cardiac_genes", []),
        }

    except ImportError:
        logger.warning("requests library not available for cardiology agent query")
        return _unavailable_response("cardiology")
    except Exception as exc:
        logger.warning("Cardiology agent query failed: %s", exc)
        return _unavailable_response("cardiology")


def query_biomarker_agent(
    metabolites: List[str] = None,
    disease_category: str = "",
    patient_profile: Optional[Dict[str, Any]] = None,
    timeout: float = None,
) -> Dict[str, Any]:
    """Query the Biomarker Intelligence Agent for metabolic profiles.

    Identifies metabolic biomarker patterns characteristic of inborn
    errors of metabolism (IEMs) and lysosomal storage disorders.

    Args:
        metabolites: List of metabolite names or results.
        disease_category: Disease category (e.g., "lysosomal", "metabolic").
        patient_profile: Optional patient metabolic data.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``metabolic_profile``, and ``biomarker_patterns``.
    """
    timeout = timeout or settings.CROSS_AGENT_TIMEOUT
    metabolites = metabolites or []

    try:
        import requests

        response = requests.post(
            f"{settings.BIOMARKER_AGENT_URL}/api/query",
            json={
                "question": f"Metabolic biomarker profile for {disease_category}",
                "biomarkers": metabolites,
                "disease_category": disease_category,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "biomarker",
            "metabolic_profile": data.get("profile", {}),
            "biomarker_patterns": data.get("patterns", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for biomarker agent query")
        return _unavailable_response("biomarker")
    except Exception as exc:
        logger.warning("Biomarker agent query failed: %s", exc)
        return _unavailable_response("biomarker")


def query_pgx_agent(
    gene: str,
    medications: List[str] = None,
    patient_profile: Optional[Dict[str, Any]] = None,
    timeout: float = None,
) -> Dict[str, Any]:
    """Query the Pharmacogenomics Agent for post-diagnosis dosing guidance.

    Checks pharmacogenomic interactions relevant to rare disease treatment
    (enzyme replacement therapy dosing, substrate reduction therapy, etc.).

    Args:
        gene: Gene symbol related to the diagnosis.
        medications: Current or proposed medications.
        patient_profile: Optional patient genomic data.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``pgx_results``, and ``dosing_recommendations``.
    """
    timeout = timeout or settings.CROSS_AGENT_TIMEOUT
    medications = medications or []

    try:
        import requests

        response = requests.post(
            f"{settings.PGX_AGENT_URL}/api/query",
            json={
                "question": f"Pharmacogenomic dosing for {gene}-related therapy",
                "patient_context": {
                    "gene": gene,
                    "medications": medications,
                    **(patient_profile or {}),
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
            "dosing_recommendations": data.get("dosing", []),
            "warnings": data.get("warnings", []),
        }

    except ImportError:
        logger.warning("requests library not available for PGx agent query")
        return _unavailable_response("pharmacogenomics")
    except Exception as exc:
        logger.warning("PGx agent query failed: %s", exc)
        return _unavailable_response("pharmacogenomics")


def query_imaging_agent(
    phenotypes: List[str] = None,
    imaging_findings: List[str] = None,
    timeout: float = None,
) -> Dict[str, Any]:
    """Query the Imaging Agent for phenotype correlation with imaging data.

    Correlates clinical phenotypes with radiological and imaging findings
    to support rare disease differential diagnosis.

    Args:
        phenotypes: List of HPO phenotype terms.
        imaging_findings: List of imaging findings.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``correlations``, and ``imaging_recommendations``.
    """
    timeout = timeout or settings.CROSS_AGENT_TIMEOUT
    phenotypes = phenotypes or []
    imaging_findings = imaging_findings or []

    try:
        import requests

        response = requests.post(
            f"{settings.GENOMICS_AGENT_URL}/api/query",
            json={
                "question": "Correlate imaging findings with phenotype",
                "phenotypes": phenotypes,
                "imaging_findings": imaging_findings,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "imaging",
            "correlations": data.get("correlations", []),
            "imaging_recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for imaging agent query")
        return _unavailable_response("imaging")
    except Exception as exc:
        logger.warning("Imaging agent query failed: %s", exc)
        return _unavailable_response("imaging")


# ===================================================================
# INTEGRATION FUNCTION
# ===================================================================


def integrate_cross_agent_results(
    results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Integrate results from multiple cross-agent queries into a unified assessment.

    Args:
        results: List of cross-agent result dicts (from the query_* functions).

    Returns:
        Unified assessment dict with agents_consulted, combined_findings,
        combined_warnings, and overall_assessment.
    """
    agents_consulted: List[str] = []
    agents_available: List[str] = []
    combined_findings: List[str] = []
    combined_warnings: List[str] = []
    combined_recommendations: List[str] = []

    for result in results:
        agent = result.get("agent", "unknown")
        agents_consulted.append(agent)

        if result.get("status") == "success":
            agents_available.append(agent)

            # Collect findings
            for key in ("cardiac_findings", "correlations", "biomarker_patterns"):
                items = result.get(key, [])
                combined_findings.extend(f"[{agent}] {item}" for item in items)

            # Collect warnings
            warnings = result.get("warnings", [])
            combined_warnings.extend(f"[{agent}] {w}" for w in warnings)

            # Collect recommendations
            recs = result.get("recommendations", result.get("dosing_recommendations", []))
            combined_recommendations.extend(f"[{agent}] {r}" for r in recs)

    # Generate overall assessment
    if not agents_available:
        overall = "No cross-agent data available. Proceeding with rare disease agent data only."
    elif combined_warnings:
        overall = (
            f"Cross-agent analysis completed with {len(combined_warnings)} warning(s). "
            f"All flagged items should be reviewed by clinical geneticist."
        )
    else:
        overall = (
            f"Cross-agent analysis completed successfully. "
            f"{len(agents_available)} agent(s) consulted with no concerns."
        )

    return {
        "agents_consulted": agents_consulted,
        "agents_available": agents_available,
        "combined_findings": combined_findings,
        "combined_warnings": combined_warnings,
        "combined_recommendations": combined_recommendations,
        "overall_assessment": overall,
    }


# ===================================================================
# HELPERS
# ===================================================================


def _unavailable_response(agent_name: str) -> Dict[str, Any]:
    """Return a standard unavailable response for a cross-agent query."""
    return {
        "status": "unavailable",
        "agent": agent_name,
        "message": f"{agent_name} agent is not currently available",
        "recommendations": [],
        "warnings": [],
    }
