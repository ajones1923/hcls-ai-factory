"""Cross-agent integration for the Precision Autoimmune Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into a unified autoimmune assessment.

Supported cross-agent queries:
  - query_oncology_agent()    -- immune-related adverse events from immunotherapy
  - query_cardiology_agent()  -- myocarditis risk from checkpoint inhibitors
  - query_neurology_agent()   -- autoimmune encephalitis risk from immunotherapy
  - query_pgx_agent()         -- drug-gene interactions for biologics
  - query_biomarker_agent()   -- autoantibody-disease activity correlation
  - query_trial_agent()       -- autoimmune/immunotherapy trial matching
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
    patient_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Oncology Intelligence Agent for immune-related adverse events.

    Assesses risk of immune-related adverse events (irAEs) from cancer
    immunotherapy in patients with pre-existing autoimmune conditions.
    Evaluates checkpoint inhibitor safety, autoimmune flare risk during
    immunotherapy, and recommends monitoring protocols.

    Args:
        patient_profile: Patient data including autoimmune diagnoses,
            current immunosuppressive therapy, and cancer type/stage.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``irae_risk_assessment``, ``monitoring_plan``,
        and ``recommendations``.
    """
    try:
        import requests

        diagnoses = patient_profile.get("autoimmune_diagnoses", [])
        cancer_type = patient_profile.get("cancer_type", "")

        response = requests.post(
            f"{settings.ONCOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Assess immune-related adverse event risk for patient with "
                    f"autoimmune conditions {diagnoses} receiving immunotherapy "
                    f"for {cancer_type}"
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
            "irae_risk_assessment": data.get("assessment", {}),
            "monitoring_plan": data.get("monitoring_plan", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for oncology agent query")
        return _unavailable_response("oncology")
    except Exception as exc:
        logger.warning("Oncology agent query failed: %s", exc)
        return _unavailable_response("oncology")


def query_cardiology_agent(
    patient_id: str,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Cardiology Intelligence Agent for myocarditis risk assessment.

    Evaluates myocarditis and pericarditis risk in autoimmune patients,
    particularly those receiving checkpoint inhibitors or biologic therapies.
    Assesses cardiac involvement in systemic autoimmune diseases such as
    lupus myocarditis, rheumatoid pericarditis, and scleroderma cardiomyopathy.

    Args:
        patient_id: Patient identifier for cardiac risk lookup.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cardiac_risk``, ``myocarditis_risk_score``,
        and ``recommendations``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.CARDIOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    "Assess myocarditis and cardiac risk from checkpoint inhibitors "
                    "and autoimmune cardiac involvement"
                ),
                "patient_context": {
                    "patient_id": patient_id,
                    "assessment_type": "autoimmune_cardiac_risk",
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "cardiology",
            "cardiac_risk": data.get("assessment", {}),
            "myocarditis_risk_score": data.get("risk_score", None),
            "risk_flags": data.get("risk_flags", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for cardiology agent query")
        return _unavailable_response("cardiology")
    except Exception as exc:
        logger.warning("Cardiology agent query failed: %s", exc)
        return _unavailable_response("cardiology")


def query_neurology_agent(
    patient_id: str,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Neurology Intelligence Agent for autoimmune encephalitis risk.

    Evaluates risk of autoimmune encephalitis and neurological immune-related
    adverse events in patients receiving immunotherapy. Assesses CNS
    demyelination risk, neuromyelitis optica spectrum disorder overlap,
    and neuropsychiatric lupus manifestations.

    Args:
        patient_id: Patient identifier for neurological risk assessment.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``encephalitis_risk``, ``neuro_assessment``,
        and ``recommendations``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.NEUROLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    "Assess autoimmune encephalitis risk and neurological "
                    "immune-related adverse events from immunotherapy"
                ),
                "patient_context": {
                    "patient_id": patient_id,
                    "assessment_type": "autoimmune_encephalitis_risk",
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "neurology",
            "encephalitis_risk": data.get("risk_assessment", {}),
            "neuro_assessment": data.get("assessment", {}),
            "risk_flags": data.get("risk_flags", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for neurology agent query")
        return _unavailable_response("neurology")
    except Exception as exc:
        logger.warning("Neurology agent query failed: %s", exc)
        return _unavailable_response("neurology")


def query_pgx_agent(
    drug_list: List[str],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the PGx Intelligence Agent for drug-gene interactions with biologics.

    Evaluates pharmacogenomic interactions for biologic therapies including
    TNF inhibitors, IL-6 blockers, JAK inhibitors, and anti-CD20 agents.
    Identifies patients at risk for adverse reactions or reduced efficacy
    based on genetic polymorphisms (e.g., TPMT, NUDT15 for thiopurines;
    CYP metabolism for JAK inhibitors).

    Args:
        drug_list: List of drug names (biologics, DMARDs, immunosuppressants)
            to evaluate for pharmacogenomic interactions.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``pgx_interactions``, ``dose_adjustments``,
        and ``recommendations``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.PGX_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Evaluate drug-gene interactions for autoimmune biologics: "
                    f"{', '.join(drug_list)}"
                ),
                "drugs": drug_list,
                "context": "autoimmune_biologics",
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "pgx",
            "pgx_interactions": data.get("interactions", []),
            "dose_adjustments": data.get("dose_adjustments", []),
            "risk_flags": data.get("risk_flags", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for PGx agent query")
        return _unavailable_response("pgx")
    except Exception as exc:
        logger.warning("PGx agent query failed: %s", exc)
        return _unavailable_response("pgx")


def query_biomarker_agent(
    autoantibody_panel: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Biomarker Intelligence Agent to correlate autoantibodies with disease activity.

    Maps autoantibody profiles to disease activity indices, flare prediction,
    and treatment response monitoring. Evaluates ANA patterns, anti-dsDNA
    titers, ANCA specificity, anti-CCP levels, and organ-specific
    autoantibodies (anti-TPO, anti-GAD, anti-PLA2R) in clinical context.

    Args:
        autoantibody_panel: Dict of autoantibody results including antibody
            names, titers, and reference ranges.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``antibody_correlation``, ``disease_activity``,
        and ``recommendations``.
    """
    try:
        import requests

        antibodies = autoantibody_panel.get("antibodies", {})

        response = requests.post(
            f"{settings.BIOMARKER_AGENT_URL}/api/query",
            json={
                "question": (
                    "Correlate autoantibody panel with disease activity "
                    "and treatment response"
                ),
                "biomarkers": list(antibodies.keys()) if isinstance(antibodies, dict) else [],
                "patient_context": autoantibody_panel,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "biomarker",
            "antibody_correlation": data.get("analysis", {}),
            "disease_activity": data.get("disease_activity", {}),
            "panel_recommendations": data.get("panel_recommendations", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for biomarker agent query")
        return _unavailable_response("biomarker")
    except Exception as exc:
        logger.warning("Biomarker agent query failed: %s", exc)
        return _unavailable_response("biomarker")


def query_trial_agent(
    autoimmune_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Clinical Trial Intelligence Agent for autoimmune/immunotherapy trials.

    Matches patients to active clinical trials for autoimmune diseases
    and immunotherapy. Considers disease type, severity, prior treatment
    failures, biomarker eligibility, and geographic location for trial
    matching across conditions including RA, SLE, vasculitis, IBD,
    and emerging CAR-T therapy for refractory autoimmune disease.

    Args:
        autoimmune_profile: Patient autoimmune profile including diagnoses,
            disease activity scores, prior treatments, and biomarkers.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``matching_trials``, ``eligibility_assessment``,
        and ``recommendations``.
    """
    try:
        import requests

        diagnoses = autoimmune_profile.get("diagnoses", [])
        prior_therapies = autoimmune_profile.get("prior_therapies", [])

        response = requests.post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Match autoimmune patient to clinical trials for "
                    f"{', '.join(diagnoses) if diagnoses else 'autoimmune disease'} "
                    f"with prior treatment failures: "
                    f"{', '.join(prior_therapies) if prior_therapies else 'none'}"
                ),
                "patient_context": autoimmune_profile,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "trial",
            "matching_trials": data.get("matches", []),
            "eligibility_assessment": data.get("eligibility", {}),
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

    Combines oncology irAE risk, cardiac safety, neurological risk,
    pharmacogenomic interactions, biomarker correlations, and trial
    matches into a single autoimmune assessment.

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
            "No cross-agent data available. Proceeding with autoimmune "
            "agent data only."
        )
    elif safety_flags:
        overall = (
            f"Cross-agent analysis identified {len(safety_flags)} safety concern(s). "
            f"Review recommended before proceeding with treatment changes."
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
