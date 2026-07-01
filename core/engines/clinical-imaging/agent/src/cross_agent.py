"""Cross-agent integration for the Imaging Intelligence Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into imaging workflow decisions, with
particular focus on the pediatric oncology pathway where coordinated
imaging across specialties (oncology staging, cardiac monitoring,
CNS surveillance) is critical.

Supported cross-agent queries:
  - query_oncology_agent()    -- staging requirements and response criteria
  - query_trial_agent()       -- trial-specific imaging requirements
  - query_cardiology_agent()  -- cardiac imaging for cardiotoxicity monitoring
  - query_neurology_agent()   -- brain MRI for CNS involvement
  - integrate_cross_agent_results() -- unified assessment

All functions degrade gracefully: if an agent is unavailable, a warning
is logged and a default response is returned so that the imaging agent
can continue with locally available data.

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
    cancer_type: str,
    patient_context: Optional[Dict[str, Any]] = None,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Oncology Intelligence Agent for staging and response criteria.

    Retrieves imaging-relevant staging requirements for a given cancer
    type, including:
      - Initial staging scan protocols (CT, PET/CT, MRI)
      - RECIST 1.1 / iRECIST / Lugano / RANO response criteria
      - Recommended imaging intervals during treatment
      - Modality-specific measurement requirements

    For pediatric oncology, this includes COG-specific staging protocols
    and age-appropriate imaging considerations (radiation dose reduction,
    sedation needs).

    Args:
        cancer_type: Cancer diagnosis (e.g., "ALL", "neuroblastoma",
            "Wilms tumor", "osteosarcoma", "medulloblastoma").
        patient_context: Optional patient data for protocol customization
            (age, body habitus, contrast allergy history).
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``staging_protocol``, ``response_criteria``,
        ``imaging_intervals``, and ``recommendations``.
    """
    try:
        import requests

        context = patient_context or {}

        response = requests.post(
            f"{settings.ONCOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Provide imaging staging requirements and response criteria "
                    f"for {cancer_type}. Include required modalities, measurement "
                    f"criteria, and recommended imaging intervals."
                ),
                "patient_context": {
                    "cancer_type": cancer_type,
                    **context,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "oncology",
            "staging_protocol": data.get("staging_protocol", {}),
            "response_criteria": data.get("response_criteria", {}),
            "imaging_intervals": data.get("imaging_intervals", []),
            "modality_requirements": data.get("modality_requirements", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for oncology agent query")
        return _unavailable_response("oncology")
    except Exception as exc:
        logger.warning("Oncology agent query failed: %s", exc)
        return _unavailable_response("oncology")


def query_trial_agent(
    patient_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Clinical Trial Intelligence Agent for imaging requirements.

    Retrieves trial-protocol-mandated imaging requirements including:
      - Required imaging modalities and sequences
      - Scan timing relative to treatment cycles
      - Measurement methodology (e.g., central review vs. local)
      - Quality assurance criteria (slice thickness, contrast timing)
      - Biomarker imaging endpoints (e.g., SUV thresholds, ADC values)

    Args:
        patient_profile: Patient data including enrolled trial ID,
            treatment arm, and imaging history.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``imaging_requirements``, ``scan_schedule``,
        ``quality_criteria``, and ``recommendations``.
    """
    try:
        import requests

        trial_id = patient_profile.get("trial_id", "")

        response = requests.post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Retrieve imaging requirements for clinical trial {trial_id}. "
                    f"Include required modalities, scan schedule, measurement "
                    f"methodology, and quality assurance criteria."
                ),
                "patient_context": patient_profile,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "clinical_trial",
            "imaging_requirements": data.get("imaging_requirements", {}),
            "scan_schedule": data.get("scan_schedule", []),
            "quality_criteria": data.get("quality_criteria", {}),
            "measurement_methodology": data.get("measurement_methodology", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for trial agent query")
        return _unavailable_response("clinical_trial")
    except Exception as exc:
        logger.warning("Trial agent query failed: %s", exc)
        return _unavailable_response("clinical_trial")


def query_cardiology_agent(
    patient_id: str,
    imaging_context: Optional[Dict[str, Any]] = None,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Cardiology Intelligence Agent for cardiac imaging coordination.

    Coordinates cardiac imaging for cardiotoxicity monitoring in patients
    receiving anthracycline-based chemotherapy.  Retrieves:
      - Current cardiac function status (LVEF, GLS baseline)
      - Recommended cardiac imaging protocol and timing
      - Cardiotoxicity risk level for imaging prioritization
      - Integration with serial biomarker data (troponin, BNP trends)

    This is critical for the pediatric oncology pathway where baseline
    and serial echocardiograms must be coordinated with staging scans
    to minimize patient visits and sedation events.

    Args:
        patient_id: Patient identifier for cross-referencing cardiac data.
        imaging_context: Optional context including current imaging schedule,
            available modalities, and scheduling constraints.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cardiac_status``, ``imaging_protocol``,
        ``risk_level``, and ``coordination_recommendations``.
    """
    try:
        import requests

        context = imaging_context or {}

        response = requests.post(
            f"{settings.CARDIOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Provide cardiac imaging requirements for patient {patient_id} "
                    f"receiving cardiotoxic chemotherapy. Include current cardiac "
                    f"function, recommended protocol, and monitoring schedule."
                ),
                "patient_context": {
                    "patient_id": patient_id,
                    "request_type": "cardiotoxicity_monitoring",
                    **context,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "cardiology",
            "cardiac_status": data.get("cardiac_status", {}),
            "imaging_protocol": data.get("imaging_protocol", {}),
            "risk_level": data.get("risk_level", "unknown"),
            "coordination_recommendations": data.get("coordination_recommendations", []),
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
    neurological_context: Optional[Dict[str, Any]] = None,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Neurology Intelligence Agent for CNS imaging coordination.

    Coordinates brain and spine MRI for patients with potential CNS
    involvement, including:
      - CNS staging for leukemia/lymphoma (diagnostic LP correlation)
      - Brain metastasis surveillance protocols
      - Posterior fossa tumor monitoring (medulloblastoma, ependymoma)
      - Neurotoxicity assessment (methotrexate leukoencephalopathy,
        radiation necrosis vs. tumor recurrence)

    For pediatric oncology, coordinates with neurology on age-appropriate
    MRI protocols, sedation planning, and developmental monitoring.

    Args:
        patient_id: Patient identifier for cross-referencing neurological data.
        neurological_context: Optional context including diagnosis, prior
            CNS imaging, neurological symptoms, and intrathecal therapy history.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cns_assessment``, ``mri_protocol``,
        ``surveillance_schedule``, and ``recommendations``.
    """
    try:
        import requests

        context = neurological_context or {}

        response = requests.post(
            f"{settings.NEUROLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Provide CNS imaging requirements for patient {patient_id}. "
                    f"Include brain/spine MRI protocol, surveillance schedule, "
                    f"and neurotoxicity assessment needs."
                ),
                "patient_context": {
                    "patient_id": patient_id,
                    "request_type": "cns_imaging_coordination",
                    **context,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "neurology",
            "cns_assessment": data.get("cns_assessment", {}),
            "mri_protocol": data.get("mri_protocol", {}),
            "surveillance_schedule": data.get("surveillance_schedule", []),
            "neurotoxicity_risk": data.get("neurotoxicity_risk", {}),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for neurology agent query")
        return _unavailable_response("neurology")
    except Exception as exc:
        logger.warning("Neurology agent query failed: %s", exc)
        return _unavailable_response("neurology")


# ===================================================================
# INTEGRATION FUNCTION
# ===================================================================


def integrate_cross_agent_results(
    results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Integrate results from multiple cross-agent queries into a unified assessment.

    Combines oncology staging requirements, trial imaging mandates,
    cardiac monitoring needs, and CNS surveillance into a single
    coordinated imaging plan for the pediatric oncology pathway.

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

            # Check cardiac risk level
            risk_level = result.get("risk_level", "")
            if risk_level in ("high", "critical"):
                safety_flags.append(
                    f"[{agent}] Cardiac risk level: {risk_level}"
                )

    # Generate overall assessment
    if not agents_available:
        overall = (
            "No cross-agent data available. Proceeding with "
            "imaging agent data only."
        )
    elif safety_flags:
        overall = (
            f"Cross-agent analysis identified {len(safety_flags)} safety concern(s). "
            f"Imaging protocol may require modification. Review recommended."
        )
    elif combined_warnings:
        overall = (
            f"Cross-agent analysis completed with {len(combined_warnings)} warning(s). "
            f"All flagged items should be reviewed by the imaging team."
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
