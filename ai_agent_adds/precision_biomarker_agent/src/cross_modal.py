"""Cross-agent integration for the Precision Biomarker Intelligence Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into unified biomarker assessments.

Supported cross-agent queries:
  - query_oncology_agent()  -- correlate biomarkers with tumor molecular profile
  - query_cart_agent()      -- validate CAR-T target suitability from biomarker data
  - query_pgx_agent()       -- check pharmacogenomic implications of biomarker findings
  - query_trial_agent()     -- match biomarker-driven clinical trials
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
    biomarker_panel: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Oncology Intelligence Agent to correlate biomarkers with tumor profile.

    Cross-references biomarker panel results (e.g. HER2, PD-L1, MSI status,
    TMB) with tumor molecular profiling data to identify actionable targets,
    validate companion diagnostic findings, and refine treatment selection.

    Args:
        biomarker_panel: Dict containing biomarker names, values, reference
            ranges, cancer type, and sample metadata.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``molecular_correlations``, and ``recommendations``.
    """
    try:
        import requests

        markers = biomarker_panel.get("markers", [])
        cancer_type = biomarker_panel.get("cancer_type", "")

        response = requests.post(
            f"{settings.ONCOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Correlate biomarker panel findings with tumor molecular "
                    f"profile for {cancer_type}: "
                    f"{', '.join(m.get('name', '') for m in markers[:10])}"
                ),
                "patient_context": biomarker_panel,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "oncology",
            "molecular_correlations": data.get("correlations", {}),
            "actionable_targets": data.get("actionable_targets", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for oncology agent query")
        return _unavailable_response("oncology")
    except Exception as exc:
        logger.warning("Oncology agent query failed: %s", exc)
        return _unavailable_response("oncology")


def query_cart_agent(
    target_antigens: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the CAR-T Intelligence Agent to validate target suitability.

    Cross-references biomarker-identified surface antigen expression data
    with CAR-T target databases to assess therapeutic candidacy. Evaluates
    antigen density thresholds, heterogeneity patterns, and on-target
    off-tumor risk from normal tissue expression profiles.

    Args:
        target_antigens: Dict containing antigen names, expression levels,
            tumor vs. normal differential, and patient diagnosis.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``target_validation``, ``safety_profile``,
        and ``recommendations``.
    """
    try:
        import requests

        antigens = target_antigens.get("antigens", [])

        response = requests.post(
            f"{settings.CART_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Validate CAR-T target suitability from biomarker "
                    f"expression data: {', '.join(antigens[:10])}"
                ),
                "patient_context": target_antigens,
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


def query_pgx_agent(
    drug_list: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Pharmacogenomics Intelligence Agent for drug-gene interactions.

    Cross-references biomarker-guided therapy recommendations with PGx
    data to identify drug-gene interactions, dosing adjustments, and
    HLA-mediated hypersensitivity risks before treatment initiation.

    Args:
        drug_list: Dict containing drug names, gene panel results,
            patient metabolizer status, and indication context.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``pgx_interactions``, ``dosing_adjustments``,
        and ``recommendations``.
    """
    try:
        import requests

        drugs = drug_list.get("drugs", [])

        response = requests.post(
            f"{settings.PGX_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Check pharmacogenomic implications for biomarker-guided "
                    f"therapy: {', '.join(drugs[:10])}"
                ),
                "patient_context": drug_list,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "pgx",
            "pgx_interactions": data.get("interactions", []),
            "dosing_adjustments": data.get("dosing_adjustments", []),
            "hla_risks": data.get("hla_risks", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for PGx agent query")
        return _unavailable_response("pgx")
    except Exception as exc:
        logger.warning("PGx agent query failed: %s", exc)
        return _unavailable_response("pgx")


def query_trial_agent(
    biomarker_panel: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Clinical Trial Intelligence Agent for biomarker-driven trials.

    Matches the patient to clinical trials that use specific biomarkers
    as inclusion/exclusion criteria (e.g. HER2+, MSI-H, TMB-high, PD-L1
    CPS >= 10) for biomarker-stratified enrollment.

    Args:
        biomarker_panel: Dict containing biomarker names, values, cancer
            type, and patient eligibility parameters.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``matched_trials``, and ``recommendations``.
    """
    try:
        import requests

        markers = biomarker_panel.get("markers", [])
        cancer_type = biomarker_panel.get("cancer_type", "")

        response = requests.post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Match patient to biomarker-driven clinical trials for "
                    f"{cancer_type} with markers: "
                    f"{', '.join(m.get('name', '') for m in markers[:10])}"
                ),
                "patient_context": biomarker_panel,
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

    Combines oncology correlations, CAR-T validation, PGx interactions,
    and trial matching into a single biomarker-driven assessment.

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
        overall = (
            "No cross-agent data available. Proceeding with "
            "biomarker agent data only."
        )
    elif safety_flags:
        overall = (
            f"Cross-agent analysis identified {len(safety_flags)} safety "
            f"concern(s). PGx interactions and drug contraindications must "
            f"be reviewed before treatment initiation."
        )
    elif combined_warnings:
        overall = (
            f"Cross-agent analysis completed with {len(combined_warnings)} "
            f"warning(s). All flagged items should be reviewed."
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
