"""Cross-agent integration for the Neurology Intelligence Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into a unified neurological assessment.

Supported cross-agent queries:
  - query_imaging_agent()      -- neuroimaging correlation
  - query_cardiology_agent()   -- cardiac/stroke overlap assessment
  - query_biomarker_agent()    -- fluid biomarker enrichment
  - query_oncology_agent()     -- neuro-oncology cross-referencing
  - query_rare_disease_agent() -- neurogenetic rare disease correlation
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


def query_imaging_agent(
    imaging_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Imaging Intelligence Agent for neuroimaging correlation.

    Sends neuroimaging findings to the imaging agent for advanced
    analysis, pattern recognition, and multi-modal correlation.

    Args:
        imaging_data: Imaging findings including modality, sequences, and
            preliminary interpretations.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``imaging_analysis``, and ``recommendations``.
    """
    try:
        import requests

        modality = imaging_data.get("modality", "")
        findings = imaging_data.get("findings", "")

        response = requests.post(
            f"{settings.IMAGING_AGENT_URL}/api/query",
            json={
                "question": f"Analyze neuroimaging findings: {findings}",
                "patient_context": {
                    "modality": modality,
                    "findings": findings,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "imaging",
            "imaging_analysis": data.get("analysis", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for imaging agent query")
        return _unavailable_response("imaging")
    except Exception as exc:
        logger.warning("Imaging agent query failed: %s", exc)
        return _unavailable_response("imaging")


def query_cardiology_agent(
    patient_profile: Dict[str, Any],
    clinical_context: str = "",
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Cardiology Intelligence Agent for cardiac/stroke overlap.

    Assesses cardioembolic risk factors, atrial fibrillation screening,
    and cardiac evaluation needs for stroke and cerebrovascular patients.

    Args:
        patient_profile: Patient data including cardiac history, ECG findings.
        clinical_context: Additional clinical context (e.g., stroke subtype).
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cardiac_assessment``, and ``risk_flags``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.CARDIOLOGY_AGENT_URL}/api/query",
            json={
                "question": f"Assess cardioembolic risk and cardiac evaluation needs: {clinical_context}",
                "patient_context": patient_profile,
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
            "af_screening": data.get("af_screening", {}),
        }

    except ImportError:
        logger.warning("requests library not available for cardiology agent query")
        return _unavailable_response("cardiology")
    except Exception as exc:
        logger.warning("Cardiology agent query failed: %s", exc)
        return _unavailable_response("cardiology")


def query_biomarker_agent(
    biomarkers: List[str],
    clinical_context: str = "",
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Biomarker Intelligence Agent for fluid biomarker enrichment.

    Identifies relevant CSF, blood, and genetic biomarkers for
    neurological conditions including Alzheimer's ATN panel, NfL,
    and autoimmune antibody panels.

    Args:
        biomarkers: List of biomarker names/results.
        clinical_context: Clinical context for biomarker interpretation.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``biomarker_analysis``, and ``panel_recommendations``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.BIOMARKER_AGENT_URL}/api/query",
            json={
                "question": f"Analyze neurological biomarkers: {clinical_context}",
                "biomarkers": biomarkers,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "biomarker",
            "biomarker_analysis": data.get("analysis", {}),
            "panel_recommendations": data.get("panel_recommendations", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for biomarker agent query")
        return _unavailable_response("biomarker")
    except Exception as exc:
        logger.warning("Biomarker agent query failed: %s", exc)
        return _unavailable_response("biomarker")


def query_oncology_agent(
    tumor_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Oncology Intelligence Agent for neuro-oncology assessment.

    Cross-references brain tumor molecular profiles with targeted therapy
    trials and precision oncology knowledge.

    Args:
        tumor_data: Tumor characteristics including WHO grade, molecular
            markers (IDH, MGMT, 1p19q), and histology.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``molecular_matches``, and ``trial_recommendations``.
    """
    try:
        import requests

        response = requests.post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            json={
                "question": "Find targeted therapy trials for brain tumor",
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


def query_rare_disease_agent(
    patient_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Rare Disease Diagnostic Agent for neurogenetic correlation.

    Cross-references neurological presentations with rare disease
    databases (OMIM, Orphanet) for genetic epilepsy, leukodystrophies,
    neurometabolic disorders, and neuromuscular genetic conditions.

    Args:
        patient_profile: Patient data including phenotype, genetic variants,
            and family history.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``genetic_matches``, and ``diagnostic_recommendations``.
    """
    try:
        import requests

        phenotypes = patient_profile.get("phenotypes", [])
        variants = patient_profile.get("genomic_variants", [])

        response = requests.post(
            f"{settings.RARE_DISEASE_AGENT_URL}/api/query",
            json={
                "question": "Evaluate neurogenetic differential diagnosis",
                "patient_context": {
                    "phenotypes": phenotypes,
                    "genomic_variants": variants,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "rare_disease",
            "genetic_matches": data.get("matches", []),
            "diagnostic_recommendations": data.get("recommendations", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for rare disease agent query")
        return _unavailable_response("rare_disease")
    except Exception as exc:
        logger.warning("Rare disease agent query failed: %s", exc)
        return _unavailable_response("rare_disease")


# ===================================================================
# INTEGRATION FUNCTION
# ===================================================================


def integrate_cross_agent_results(
    results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Integrate results from multiple cross-agent queries into a unified assessment.

    Combines imaging analysis, cardiac assessment, biomarker data,
    oncology matches, and rare disease correlations into a single
    neurological assessment.

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
        overall = "No cross-agent data available. Proceeding with neurology agent data only."
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
