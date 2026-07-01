"""Cross-agent integration for the CAR-T Intelligence Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into unified CAR-T therapy assessments.

Supported cross-agent queries:
  - query_biomarker_agent()   -- target antigen expression data
  - query_oncology_agent()    -- tumor profile and disease context
  - query_single_cell_agent() -- TME profiling for target validation
  - query_cardiology_agent()  -- baseline cardiac assessment pre-lymphodepletion
  - query_trial_agent()       -- CAR-T clinical trial matching
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


def query_biomarker_agent(
    target_antigens: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Precision Biomarker Agent for target antigen expression data.

    Retrieves patient-specific target expression levels (e.g. CD19, CD22,
    BCMA) to inform CAR-T construct selection and predict antigen escape
    risk. Expression heterogeneity data guides dual-targeting decisions.

    Args:
        target_antigens: Dict containing target gene symbols, patient ID,
            and sample type (e.g. bone marrow, peripheral blood).
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``expression_data``, and ``recommendations``.
    """
    try:
        import requests

        antigens = target_antigens.get("antigens", [])

        response = requests.post(
            f"{settings.BIOMARKER_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Provide expression levels and heterogeneity data for "
                    f"CAR-T target antigens: {', '.join(antigens[:10])}"
                ),
                "patient_context": target_antigens,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "biomarker",
            "expression_data": data.get("expression", {}),
            "heterogeneity": data.get("heterogeneity", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for biomarker agent query")
        return _unavailable_response("biomarker")
    except Exception as exc:
        logger.warning("Biomarker agent query failed: %s", exc)
        return _unavailable_response("biomarker")


def query_oncology_agent(
    patient_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Oncology Intelligence Agent for tumor profile and disease context.

    Retrieves tumor molecular profile, prior treatment lines, disease burden,
    and staging data to contextualize CAR-T eligibility and predict response.
    Critical for bridging therapy selection and lymphodepletion planning.

    Args:
        patient_profile: Dict containing cancer type, stage, prior therapies,
            molecular markers, and performance status.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``tumor_profile``, ``prior_therapies``, and
        ``recommendations``.
    """
    try:
        import requests

        cancer_type = patient_profile.get("cancer_type", "")
        stage = patient_profile.get("stage", "")

        response = requests.post(
            f"{settings.ONCOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Provide tumor molecular profile and treatment history "
                    f"for CAR-T candidacy assessment: {cancer_type} stage {stage}"
                ),
                "patient_context": patient_profile,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "oncology",
            "tumor_profile": data.get("tumor_profile", {}),
            "prior_therapies": data.get("prior_therapies", []),
            "disease_burden": data.get("disease_burden", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for oncology agent query")
        return _unavailable_response("oncology")
    except Exception as exc:
        logger.warning("Oncology agent query failed: %s", exc)
        return _unavailable_response("oncology")


def query_single_cell_agent(
    tumor_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Single-Cell Intelligence Agent for TME profiling.

    Retrieves tumor microenvironment composition at single-cell resolution
    to validate CAR-T target expression uniformity, identify immunosuppressive
    niches, and predict T-cell exhaustion risk in the TME.

    Args:
        tumor_data: Dict containing tumor sample ID, cancer type, biopsy
            site, and target antigens for TME-level validation.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``tme_profile``, ``target_validation``, and
        ``recommendations``.
    """
    try:
        import requests

        cancer_type = tumor_data.get("cancer_type", "")
        targets = tumor_data.get("target_antigens", [])

        response = requests.post(
            f"{settings.SINGLE_CELL_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Profile TME composition and validate CAR-T target "
                    f"expression ({', '.join(targets[:5])}) for {cancer_type}"
                ),
                "patient_context": tumor_data,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "single_cell",
            "tme_profile": data.get("tme_profile", {}),
            "target_validation": data.get("target_validation", {}),
            "exhaustion_risk": data.get("exhaustion_risk", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for single-cell agent query")
        return _unavailable_response("single_cell")
    except Exception as exc:
        logger.warning("Single-cell agent query failed: %s", exc)
        return _unavailable_response("single_cell")


def query_cardiology_agent(
    patient_id: str,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Cardiology Intelligence Agent for baseline cardiac assessment.

    Retrieves cardiac function data (LVEF, troponin, BNP, ECG findings)
    required before lymphodepletion chemotherapy. Cytokine release syndrome
    (CRS) carries significant cardiovascular risk; baseline assessment
    guides CRS management and tocilizumab timing.

    Args:
        patient_id: Patient identifier for cardiac history lookup.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cardiac_assessment``, ``risk_flags``, and
        ``recommendations``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.CARDIOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    "Provide baseline cardiac assessment for CAR-T "
                    "lymphodepletion candidacy including LVEF, troponin, "
                    "BNP, and arrhythmia history"
                ),
                "patient_context": {"patient_id": patient_id},
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "cardiology",
            "cardiac_assessment": data.get("assessment", {}),
            "risk_flags": data.get("risk_flags", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for cardiology agent query")
        return _unavailable_response("cardiology")
    except Exception as exc:
        logger.warning("Cardiology agent query failed: %s", exc)
        return _unavailable_response("cardiology")


def query_trial_agent(
    cart_product: Dict[str, Any],
    patient_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Clinical Trial Intelligence Agent for CAR-T trial matching.

    Matches the patient to active CAR-T clinical trials based on product
    type (e.g. axi-cel, liso-cel, cilta-cel), disease indication, prior
    lines of therapy, and molecular eligibility criteria.

    Args:
        cart_product: Dict containing CAR-T product name, target antigen,
            construct generation, and costimulatory domain.
        patient_profile: Dict containing diagnosis, stage, prior therapies,
            and performance status.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``matched_trials``, and ``recommendations``.
    """
    try:
        import requests

        product_name = cart_product.get("product_name", "")
        cancer_type = patient_profile.get("cancer_type", "")

        response = requests.post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Match patient to CAR-T clinical trials for "
                    f"{product_name} in {cancer_type}"
                ),
                "patient_context": {
                    **patient_profile,
                    "cart_product": cart_product,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "trial",
            "matched_trials": data.get("trials", []),
            "eligibility_summary": data.get("eligibility_summary", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for trial agent query")
        return _unavailable_response("trial")
    except Exception as exc:
        logger.warning("Trial agent query failed: %s", exc)
        return _unavailable_response("trial")


# ===================================================================
# INTEGRATION FUNCTION
# ===================================================================


def integrate_cross_agent_results(
    results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Integrate results from multiple cross-agent queries into a unified assessment.

    Combines biomarker expression, oncology context, TME profiling, cardiac
    risk, and trial matching into a single CAR-T therapy assessment.

    Args:
        results: List of cross-agent result dicts (from the query_* functions).

    Returns:
        Unified assessment dict with:
          - ``agents_consulted``: List of agent names queried.
          - ``agents_available``: List of agents that responded successfully.
          - ``combined_warnings``: Aggregated warnings from all agents.
          - ``combined_recommendations``: Aggregated recommendations.
          - ``safety_flags``: Combined safety concerns (cardiac, CRS, neurotox).
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
        overall = (
            "No cross-agent data available. Proceeding with "
            "CAR-T agent data only."
        )
    elif safety_flags:
        overall = (
            f"Cross-agent analysis identified {len(safety_flags)} safety "
            f"concern(s) relevant to CAR-T therapy. Cardiac and CRS risk "
            f"factors must be reviewed before lymphodepletion."
        )
    elif combined_warnings:
        overall = (
            f"Cross-agent analysis completed with {len(combined_warnings)} "
            f"warning(s). All flagged items should be reviewed before "
            f"CAR-T infusion."
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
