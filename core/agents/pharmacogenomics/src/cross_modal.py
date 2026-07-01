"""Cross-agent integration for the Pharmacogenomics Intelligence Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into unified pharmacogenomic assessments.

Supported cross-agent queries:
  - query_oncology_agent()   -- get planned therapy for drug-gene interaction context
  - query_cardiology_agent() -- assess cardiac drug interactions
  - query_neurology_agent()  -- assess neurotoxic drug interactions
  - query_trial_agent()      -- match PGx-guided clinical trials
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
    """Query the Oncology Intelligence Agent for planned therapy context.

    Retrieves the patient's planned or active chemotherapy regimen to
    contextualize drug-gene interaction analysis. Critical for identifying
    high-risk PGx pairs (e.g. DPYD/fluoropyrimidines, UGT1A1/irinotecan,
    TPMT/thiopurines) before treatment initiation.

    Args:
        patient_profile: Dict containing cancer type, planned regimen,
            prior therapies, and current medications.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``planned_therapies``, ``drug_gene_context``,
        and ``recommendations``.
    """
    try:
        import requests

        cancer_type = patient_profile.get("cancer_type", "")
        regimen = patient_profile.get("planned_regimen", "")

        response = requests.post(
            f"{settings.ONCOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Provide planned therapy details for PGx drug-gene "
                    f"interaction assessment: {cancer_type}, regimen {regimen}"
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
            "planned_therapies": data.get("planned_therapies", []),
            "drug_gene_context": data.get("drug_gene_context", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for oncology agent query")
        return _unavailable_response("oncology")
    except Exception as exc:
        logger.warning("Oncology agent query failed: %s", exc)
        return _unavailable_response("oncology")


def query_cardiology_agent(
    drug_list: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Cardiology Intelligence Agent for cardiac drug interaction assessment.

    Evaluates cardiotoxicity risk for prescribed drugs in the context of
    pharmacogenomic findings. Identifies QT-prolongation risk (e.g.
    CYP2D6 poor metabolizers on tamoxifen/ondansetron), cardiomyopathy
    risk with anthracyclines, and arrhythmia risk with tyrosine kinase
    inhibitors.

    Args:
        drug_list: Dict containing drug names, PGx metabolizer phenotypes,
            cardiac history, and concurrent medications.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cardiac_risk_assessment``, ``risk_flags``,
        and ``recommendations``.
    """
    try:
        import requests

        drugs = drug_list.get("drugs", [])

        response = requests.post(
            f"{settings.CARDIOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Assess cardiac drug interaction risk for PGx-flagged "
                    f"medications: {', '.join(drugs[:10])}"
                ),
                "patient_context": drug_list,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "cardiology",
            "cardiac_risk_assessment": data.get("assessment", {}),
            "qt_prolongation_risk": data.get("qt_risk", {}),
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
    drug_list: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Neurology Intelligence Agent for neurotoxic drug interaction assessment.

    Evaluates neurotoxicity risk for prescribed drugs in the context of
    pharmacogenomic findings. Identifies peripheral neuropathy risk (e.g.
    CYP2C19 poor metabolizers on voriconazole), CNS toxicity with
    methotrexate, and seizure threshold alterations with enzyme-inhibiting
    drug combinations.

    Args:
        drug_list: Dict containing drug names, PGx metabolizer phenotypes,
            neurological history, and concurrent medications.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``neurotoxicity_assessment``, ``risk_flags``,
        and ``recommendations``.
    """
    try:
        import requests

        drugs = drug_list.get("drugs", [])

        response = requests.post(
            f"{settings.NEUROLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Assess neurotoxicity risk for PGx-flagged "
                    f"medications: {', '.join(drugs[:10])}"
                ),
                "patient_context": drug_list,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "neurology",
            "neurotoxicity_assessment": data.get("assessment", {}),
            "neuropathy_risk": data.get("neuropathy_risk", {}),
            "risk_flags": data.get("risk_flags", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for neurology agent query")
        return _unavailable_response("neurology")
    except Exception as exc:
        logger.warning("Neurology agent query failed: %s", exc)
        return _unavailable_response("neurology")


def query_trial_agent(
    pgx_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Clinical Trial Intelligence Agent for PGx-guided trials.

    Matches the patient to clinical trials that incorporate pharmacogenomic
    stratification, genotype-guided dosing arms, or PGx biomarker-based
    inclusion criteria (e.g. CYP2D6 ultrarapid metabolizer exclusion,
    DPYD-deficiency adapted dosing protocols).

    Args:
        pgx_profile: Dict containing gene-drug pairs, metabolizer
            phenotypes, diplotypes, and patient eligibility parameters.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``matched_trials``, and ``recommendations``.
    """
    try:
        import requests

        gene_drug_pairs = pgx_profile.get("gene_drug_pairs", [])

        response = requests.post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Match patient to PGx-guided clinical trials based on "
                    f"pharmacogenomic profile: "
                    f"{', '.join(str(p) for p in gene_drug_pairs[:10])}"
                ),
                "patient_context": pgx_profile,
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

    Combines oncology therapy context, cardiac risk, neurotoxicity risk,
    and trial matching into a single PGx-guided assessment.

    Args:
        results: List of cross-agent result dicts (from the query_* functions).

    Returns:
        Unified assessment dict with:
          - ``agents_consulted``: List of agent names queried.
          - ``agents_available``: List of agents that responded successfully.
          - ``combined_warnings``: Aggregated warnings from all agents.
          - ``combined_recommendations``: Aggregated recommendations.
          - ``safety_flags``: Combined safety concerns (cardiac, neuro, drug-gene).
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
            "PGx agent data only."
        )
    elif safety_flags:
        overall = (
            f"Cross-agent analysis identified {len(safety_flags)} safety "
            f"concern(s). Cardiac and neurotoxicity risk factors must be "
            f"reviewed before prescribing PGx-flagged medications."
        )
    elif combined_warnings:
        overall = (
            f"Cross-agent analysis completed with {len(combined_warnings)} "
            f"warning(s). All flagged drug-gene interactions should be "
            f"reviewed by clinical pharmacology."
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
