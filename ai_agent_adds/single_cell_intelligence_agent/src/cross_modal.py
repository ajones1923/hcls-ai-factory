"""Cross-agent integration for the Single-Cell Intelligence Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into unified single-cell analyses.

Supported cross-agent queries:
  - query_oncology_agent()       -- tumor molecular profiling correlation
  - query_cart_agent()           -- CAR-T target validation enrichment
  - query_biomarker_agent()      -- biomarker panel enrichment
  - query_drug_discovery_agent() -- drug candidate integration
  - query_imaging_agent()        -- spatial-imaging correlation
  - integrate_cross_agent_results() -- unified assessment

All functions degrade gracefully: if an agent is unavailable, a warning
is logged and a default response is returned.

Author: Adam Jones
Date: March 2026
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

from config.settings import settings

logger = logging.getLogger(__name__)


# ===================================================================
# CROSS-AGENT QUERY FUNCTIONS
# ===================================================================


def query_oncology_agent(
    tumor_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Oncology Intelligence Agent for tumor profiling correlation.

    Cross-references single-cell TME profiles with bulk tumor molecular
    data, clinical trial eligibility, and targeted therapy recommendations.

    Args:
        tumor_data: Tumor characteristics including cancer type, molecular
            markers, and TME classification from single-cell analysis.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``molecular_matches``, and ``recommendations``.
    """
    try:
        import requests

        cancer_type = tumor_data.get("cancer_type", "")
        tme_class = tumor_data.get("tme_class", "")

        response = requests.post(
            f"{settings.ONCOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Correlate single-cell TME profile with molecular "
                    f"targets for {cancer_type} ({tme_class})"
                ),
                "patient_context": tumor_data,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "oncology",
            "molecular_matches": data.get("matches", []),
            "trial_recommendations": data.get("recommendations", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for oncology agent query")
        return _unavailable_response("oncology")
    except Exception as exc:
        logger.warning("Oncology agent query failed: %s", exc)
        return _unavailable_response("oncology")


def query_cart_agent(
    target_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the CAR-T Intelligence Agent for target validation enrichment.

    Cross-references single-cell expression profiles with CAR-T target
    databases, safety profiles, and clinical trial data.

    Args:
        target_data: Target antigen data including expression levels,
            normal tissue expression, and heterogeneity scores.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``target_validation``, and ``safety_profile``.
    """
    try:
        import requests

        target_gene = target_data.get("target_gene", "")

        response = requests.post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Validate CAR-T target {target_gene} using single-cell "
                    f"expression data and assess on-target off-tumor risk"
                ),
                "patient_context": target_data,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "cart",
            "target_validation": data.get("validation", {}),
            "safety_profile": data.get("safety", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for CAR-T agent query")
        return _unavailable_response("cart")
    except Exception as exc:
        logger.warning("CAR-T agent query failed: %s", exc)
        return _unavailable_response("cart")


def query_biomarker_agent(
    biomarker_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Biomarker Intelligence Agent for panel enrichment.

    Cross-references single-cell-derived biomarker candidates with
    established biomarker databases and clinical assay availability.

    Args:
        biomarker_data: Biomarker candidates including gene symbols,
            cell type specificity, and expression metrics.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``biomarker_validation``, and ``panel_recommendations``.
    """
    try:
        import requests

        genes = biomarker_data.get("genes", [])

        response = requests.post(
            f"{settings.BIOMARKER_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Validate single-cell biomarker candidates: {', '.join(genes[:10])}"
                ),
                "biomarkers": genes,
                "patient_context": biomarker_data,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "biomarker",
            "biomarker_validation": data.get("validation", {}),
            "panel_recommendations": data.get("panel_recommendations", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for biomarker agent query")
        return _unavailable_response("biomarker")
    except Exception as exc:
        logger.warning("Biomarker agent query failed: %s", exc)
        return _unavailable_response("biomarker")


def query_drug_discovery_agent(
    drug_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Drug Discovery pipeline for compound integration.

    Cross-references single-cell drug response predictions with
    the BioNeMo/MolMIM/DiffDock drug discovery pipeline for structural
    insights and docking validation.

    Args:
        drug_data: Drug response data including predicted sensitivities,
            target genes, and resistance mechanisms.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``docking_results``, and ``compound_suggestions``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.GENOMICS_AGENT_URL}/api/query",
            json={
                "question": "Validate drug candidates from single-cell predictions",
                "patient_context": drug_data,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "drug_discovery",
            "docking_results": data.get("docking", {}),
            "compound_suggestions": data.get("compounds", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for drug discovery query")
        return _unavailable_response("drug_discovery")
    except Exception as exc:
        logger.warning("Drug discovery agent query failed: %s", exc)
        return _unavailable_response("drug_discovery")


def query_imaging_agent(
    spatial_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Imaging Intelligence Agent for spatial-imaging correlation.

    Cross-references spatial transcriptomics niches with histopathology
    imaging features for multi-modal tissue characterization.

    Args:
        spatial_data: Spatial niche data including cell type composition,
            niche labels, and platform metadata.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``imaging_correlation``, and ``recommendations``.
    """
    try:
        import requests

        spatial_data.get("niches", [])
        platform = spatial_data.get("platform", "unknown")

        response = requests.post(
            f"{settings.GENOMICS_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Correlate spatial transcriptomics niches ({platform}) "
                    f"with histopathology imaging features"
                ),
                "patient_context": spatial_data,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "imaging",
            "imaging_correlation": data.get("correlation", {}),
            "recommendations": data.get("recommendations", []),
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

    Combines oncology correlations, CAR-T validation, biomarker enrichment,
    drug discovery insights, and imaging correlations into a single
    multi-modal assessment.

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

    # Generate overall assessment
    if not agents_available:
        overall = "No cross-agent data available. Proceeding with single-cell agent data only."
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
    """Return a standard unavailable response for a cross-agent query."""
    return {
        "status": "unavailable",
        "agent": agent_name,
        "message": f"{agent_name} agent is not currently available",
        "recommendations": [],
        "warnings": [],
    }
