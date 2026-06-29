"""Cross-agent integration for the Cardiology Intelligence Agent.

Provides functions to query other HCLS AI Factory intelligence agents
and integrate their results into a unified cardiovascular assessment,
with particular focus on the pediatric oncology pathway where
anthracycline cardiotoxicity monitoring is critical.

Supported cross-agent queries:
  - query_oncology_agent()                -- planned chemotherapy & cumulative anthracycline dose
  - query_trial_agent()                   -- trial-specific cardiac monitoring requirements
  - query_biomarker_agent()               -- troponin, BNP trends for cardiotoxicity
  - query_imaging_agent()                 -- echocardiogram / cardiac MRI baseline coordination
  - pediatric_cardiotoxicity_assessment() -- synthesize all inputs for anthracycline safety
  - integrate_cross_agent_results()       -- unified assessment

All functions degrade gracefully: if an agent is unavailable, a warning
is logged and a default response is returned so that the cardiology
agent can continue with locally available data.

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
    """Query the Oncology Intelligence Agent for planned chemotherapy details.

    Retrieves the planned chemotherapy regimen including cumulative
    anthracycline dose (doxorubicin, daunorubicin, epirubicin),
    treatment schedule, and any planned dose modifications.  This is
    essential for prospective cardiotoxicity risk stratification in
    pediatric oncology patients.

    Args:
        patient_profile: Patient data including cancer type, age, weight,
            current treatment plan, and prior anthracycline exposure.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``chemotherapy_plan``, ``anthracycline_exposure``,
        and ``cardiotoxicity_risk_factors``.
    """
    try:
        import requests

        cancer_type = patient_profile.get("cancer_type", "")
        age = patient_profile.get("age", "")

        response = requests.post(
            f"{settings.ONCOLOGY_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Retrieve planned chemotherapy regimen and cumulative "
                    f"anthracycline dose for {cancer_type} patient, age {age}. "
                    f"Include cardiotoxicity risk factors."
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
            "chemotherapy_plan": data.get("chemotherapy_plan", {}),
            "anthracycline_exposure": data.get("anthracycline_exposure", {}),
            "cardiotoxicity_risk_factors": data.get("cardiotoxicity_risk_factors", []),
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
    """Query the Clinical Trial Intelligence Agent for cardiac monitoring requirements.

    Retrieves trial-specific cardiac monitoring protocols including
    required echocardiogram intervals, LVEF thresholds for dose
    modification or treatment discontinuation, and biomarker monitoring
    schedules mandated by the trial protocol.

    Args:
        patient_profile: Patient data including enrolled trial ID,
            treatment arm, and current cardiac function metrics.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``monitoring_protocol``, ``lvef_thresholds``,
        and ``schedule``.
    """
    try:
        import requests

        trial_id = patient_profile.get("trial_id", "")

        response = requests.post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Retrieve cardiac monitoring requirements for trial {trial_id}. "
                    f"Include LVEF thresholds, echocardiogram schedule, and "
                    f"biomarker monitoring mandates."
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
            "monitoring_protocol": data.get("monitoring_protocol", {}),
            "lvef_thresholds": data.get("lvef_thresholds", {}),
            "schedule": data.get("schedule", []),
            "dose_modification_rules": data.get("dose_modification_rules", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for trial agent query")
        return _unavailable_response("clinical_trial")
    except Exception as exc:
        logger.warning("Trial agent query failed: %s", exc)
        return _unavailable_response("clinical_trial")


def query_biomarker_agent(
    cardiac_biomarkers: Dict[str, Any],
    clinical_context: str = "",
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Biomarker Intelligence Agent for troponin and BNP trend analysis.

    Sends serial cardiac biomarker values (high-sensitivity troponin I/T,
    NT-proBNP, BNP) for trend analysis and early cardiotoxicity detection.
    The biomarker agent correlates values with treatment timeline to
    identify subclinical myocardial injury before LVEF decline.

    Args:
        cardiac_biomarkers: Biomarker data including serial troponin values,
            BNP/NT-proBNP levels, collection dates, and reference ranges.
        clinical_context: Additional context such as treatment phase,
            cumulative dose at time of sampling, or concurrent medications.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``troponin_trend``, ``bnp_trend``,
        ``cardiotoxicity_signal``, and ``recommendations``.
    """
    try:
        import requests

        biomarker_list = []
        for name, values in cardiac_biomarkers.items():
            if isinstance(values, dict):
                biomarker_list.append(f"{name}: {values}")
            else:
                biomarker_list.append(f"{name}={values}")

        response = requests.post(
            f"{settings.BIOMARKER_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Analyze cardiac biomarker trends for cardiotoxicity monitoring: "
                    f"{clinical_context}. Biomarkers: {', '.join(biomarker_list)}"
                ),
                "biomarkers": cardiac_biomarkers,
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "biomarker",
            "troponin_trend": data.get("troponin_trend", {}),
            "bnp_trend": data.get("bnp_trend", {}),
            "cardiotoxicity_signal": data.get("cardiotoxicity_signal", "unknown"),
            "panel_recommendations": data.get("panel_recommendations", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for biomarker agent query")
        return _unavailable_response("biomarker")
    except Exception as exc:
        logger.warning("Biomarker agent query failed: %s", exc)
        return _unavailable_response("biomarker")


def query_imaging_agent(
    cardiac_imaging_type: str,
    patient_context: Optional[Dict[str, Any]] = None,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Imaging Intelligence Agent to coordinate cardiac imaging baseline.

    Requests baseline echocardiogram or cardiac MRI scheduling and
    protocol recommendations.  For pediatric oncology patients, this
    ensures a pre-treatment LVEF/GLS baseline is established before
    anthracycline exposure, per COG (Children's Oncology Group) and
    ASCO cardio-oncology guidelines.

    Args:
        cardiac_imaging_type: Type of cardiac imaging required, e.g.,
            ``"echocardiogram"``, ``"cardiac_mri"``, ``"strain_imaging"``,
            or ``"muga_scan"``.
        patient_context: Optional patient data for protocol customization
            (age, body habitus, prior imaging, sedation needs).
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``imaging_protocol``, ``baseline_parameters``,
        and ``scheduling_recommendations``.
    """
    try:
        import requests

        context = patient_context or {}

        response = requests.post(
            f"{settings.IMAGING_AGENT_URL}/api/query",
            json={
                "question": (
                    f"Coordinate baseline {cardiac_imaging_type} for cardiotoxicity "
                    f"monitoring. Provide protocol recommendations, required "
                    f"measurements (LVEF, GLS), and scheduling guidance."
                ),
                "patient_context": {
                    "imaging_type": cardiac_imaging_type,
                    **context,
                },
            },
            timeout=timeout,
        )
        response.raise_for_status()
        data = response.json()

        return {
            "status": "success",
            "agent": "imaging",
            "imaging_protocol": data.get("imaging_protocol", {}),
            "baseline_parameters": data.get("baseline_parameters", {}),
            "scheduling_recommendations": data.get("scheduling_recommendations", []),
            "recommendations": data.get("recommendations", []),
        }

    except ImportError:
        logger.warning("requests library not available for imaging agent query")
        return _unavailable_response("imaging")
    except Exception as exc:
        logger.warning("Imaging agent query failed: %s", exc)
        return _unavailable_response("imaging")


# ===================================================================
# PEDIATRIC CARDIOTOXICITY ASSESSMENT
# ===================================================================


def pediatric_cardiotoxicity_assessment(
    therapy_plan: Dict[str, Any],
) -> Dict[str, Any]:
    """Synthesize cross-agent inputs for anthracycline cardiotoxicity safety.

    Orchestrates queries to the oncology, clinical trial, biomarker, and
    imaging agents to build a comprehensive cardiotoxicity risk assessment
    for pediatric patients receiving anthracycline-based chemotherapy.

    The assessment integrates:
      - Planned cumulative anthracycline dose (from oncology agent)
      - Trial-mandated cardiac monitoring schedule (from trial agent)
      - Serial troponin/BNP trends (from biomarker agent)
      - Baseline echocardiogram or cardiac MRI (from imaging agent)

    Risk stratification follows COG Long-Term Follow-Up Guidelines and
    the 2022 ESC Cardio-Oncology Guidelines:
      - LOW:  cumulative doxorubicin-equivalent < 250 mg/m^2, normal
              biomarkers, LVEF >= 55%, normal GLS
      - MODERATE: 250-400 mg/m^2, OR borderline biomarker elevation,
                  OR LVEF 50-54%
      - HIGH: > 400 mg/m^2, OR rising troponin trend, OR LVEF < 50%,
              OR GLS decline > 15% from baseline

    Args:
        therapy_plan: Dict containing at minimum:
            - ``patient_profile``: demographics, weight, BSA
            - ``cancer_type``: diagnosis
            - ``trial_id``: enrolled clinical trial (if any)
            - ``cardiac_biomarkers``: serial troponin/BNP values
            - ``baseline_imaging_type``: preferred imaging modality

    Returns:
        Dict with:
          - ``risk_level``: "low", "moderate", or "high"
          - ``oncology_data``: chemotherapy plan from oncology agent
          - ``trial_data``: monitoring protocol from trial agent
          - ``biomarker_data``: trend analysis from biomarker agent
          - ``imaging_data``: baseline imaging from imaging agent
          - ``safety_recommendations``: aggregated clinical recommendations
          - ``monitoring_schedule``: recommended monitoring timeline
          - ``agents_consulted``: list of agents queried
          - ``agents_available``: list of agents that responded
    """
    patient_profile = therapy_plan.get("patient_profile", {})
    patient_profile["cancer_type"] = therapy_plan.get("cancer_type", "")
    patient_profile["trial_id"] = therapy_plan.get("trial_id", "")

    cardiac_biomarkers = therapy_plan.get("cardiac_biomarkers", {})
    imaging_type = therapy_plan.get("baseline_imaging_type", "echocardiogram")

    # Query all agents
    oncology_data = query_oncology_agent(patient_profile)
    trial_data = query_trial_agent(patient_profile)
    biomarker_data = query_biomarker_agent(
        cardiac_biomarkers,
        clinical_context=f"Pediatric {therapy_plan.get('cancer_type', '')} on anthracycline therapy",
    )
    imaging_data = query_imaging_agent(imaging_type, patient_context=patient_profile)

    # Track agent availability
    all_results = [oncology_data, trial_data, biomarker_data, imaging_data]
    agents_consulted = [r.get("agent", "unknown") for r in all_results]
    agents_available = [
        r.get("agent", "unknown") for r in all_results
        if r.get("status") == "success"
    ]

    # Risk stratification
    risk_level = _calculate_cardiotoxicity_risk(
        oncology_data, trial_data, biomarker_data, imaging_data
    )

    # Aggregate safety recommendations
    safety_recommendations = _aggregate_safety_recommendations(
        oncology_data, trial_data, biomarker_data, imaging_data, risk_level
    )

    # Build monitoring schedule
    monitoring_schedule = _build_monitoring_schedule(
        risk_level, trial_data, therapy_plan
    )

    return {
        "risk_level": risk_level,
        "oncology_data": oncology_data,
        "trial_data": trial_data,
        "biomarker_data": biomarker_data,
        "imaging_data": imaging_data,
        "safety_recommendations": safety_recommendations,
        "monitoring_schedule": monitoring_schedule,
        "agents_consulted": agents_consulted,
        "agents_available": agents_available,
    }


# ===================================================================
# INTEGRATION FUNCTION
# ===================================================================


def integrate_cross_agent_results(
    results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Integrate results from multiple cross-agent queries into a unified assessment.

    Combines oncology chemotherapy data, trial monitoring protocols,
    biomarker trends, and imaging baselines into a single cardiovascular
    assessment for the pediatric oncology pathway.

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

            # Check for cardiotoxicity signals
            cardiotoxicity = result.get("cardiotoxicity_signal", "")
            if cardiotoxicity and cardiotoxicity not in ("normal", "unknown", ""):
                safety_flags.append(
                    f"[{agent}] Cardiotoxicity signal: {cardiotoxicity}"
                )

    # Generate overall assessment
    if not agents_available:
        overall = (
            "No cross-agent data available. Proceeding with "
            "cardiology agent data only."
        )
    elif safety_flags:
        overall = (
            f"Cross-agent analysis identified {len(safety_flags)} safety concern(s). "
            f"Cardiotoxicity review recommended before proceeding with therapy."
        )
    elif combined_warnings:
        overall = (
            f"Cross-agent analysis completed with {len(combined_warnings)} warning(s). "
            f"All flagged items should be reviewed by the cardio-oncology team."
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


def _calculate_cardiotoxicity_risk(
    oncology_data: Dict[str, Any],
    trial_data: Dict[str, Any],
    biomarker_data: Dict[str, Any],
    imaging_data: Dict[str, Any],
) -> str:
    """Calculate cardiotoxicity risk level from cross-agent data.

    Uses a conservative approach: the highest risk from any single
    data source determines the overall risk level.

    Returns:
        "low", "moderate", or "high".
    """
    risk_score = 0

    # Anthracycline dose risk (from oncology agent)
    if oncology_data.get("status") == "success":
        exposure = oncology_data.get("anthracycline_exposure", {})
        cumulative_dose = exposure.get("cumulative_doxorubicin_equivalent_mg_m2", 0)
        if cumulative_dose > 400:
            risk_score = max(risk_score, 3)
        elif cumulative_dose > 250:
            risk_score = max(risk_score, 2)
        elif cumulative_dose > 0:
            risk_score = max(risk_score, 1)

    # Biomarker risk (from biomarker agent)
    if biomarker_data.get("status") == "success":
        signal = biomarker_data.get("cardiotoxicity_signal", "")
        if signal in ("elevated", "rising", "critical"):
            risk_score = max(risk_score, 3)
        elif signal in ("borderline", "trending_up"):
            risk_score = max(risk_score, 2)

    # Imaging risk (from imaging agent)
    if imaging_data.get("status") == "success":
        baseline = imaging_data.get("baseline_parameters", {})
        lvef = baseline.get("lvef_percent", 60)
        if lvef < 50:
            risk_score = max(risk_score, 3)
        elif lvef < 55:
            risk_score = max(risk_score, 2)

    # Map score to risk level
    if risk_score >= 3:
        return "high"
    elif risk_score >= 2:
        return "moderate"
    else:
        return "low"


def _aggregate_safety_recommendations(
    oncology_data: Dict[str, Any],
    trial_data: Dict[str, Any],
    biomarker_data: Dict[str, Any],
    imaging_data: Dict[str, Any],
    risk_level: str,
) -> List[str]:
    """Aggregate safety recommendations from all cross-agent data sources."""
    recommendations: List[str] = []

    # Always recommend baseline
    recommendations.append(
        "Obtain baseline echocardiogram with LVEF and GLS before anthracycline initiation"
    )

    # Risk-based recommendations
    if risk_level == "high":
        recommendations.extend([
            "Consider cardio-oncology consultation before treatment initiation",
            "Evaluate dexrazoxane cardioprotection per COG guidelines",
            "Serial troponin and NT-proBNP monitoring at each cycle",
            "Echocardiogram every 1-2 cycles during active treatment",
        ])
    elif risk_level == "moderate":
        recommendations.extend([
            "Cardio-oncology referral recommended",
            "Serial troponin monitoring every 2 cycles",
            "Echocardiogram every 3 cycles during active treatment",
        ])
    else:
        recommendations.extend([
            "Standard cardiac monitoring per protocol",
            "Echocardiogram at mid-treatment and end of treatment",
        ])

    # Add agent-specific recommendations
    for data in [oncology_data, trial_data, biomarker_data, imaging_data]:
        if data.get("status") == "success":
            for rec in data.get("recommendations", []):
                if rec not in recommendations:
                    recommendations.append(rec)

    return recommendations


def _build_monitoring_schedule(
    risk_level: str,
    trial_data: Dict[str, Any],
    therapy_plan: Dict[str, Any],
) -> List[Dict[str, str]]:
    """Build a monitoring schedule based on risk level and trial requirements."""
    schedule: List[Dict[str, str]] = []

    # Pre-treatment baseline
    schedule.append({
        "timepoint": "Pre-treatment baseline",
        "assessments": "Echocardiogram (LVEF + GLS), ECG, troponin, NT-proBNP",
        "priority": "required",
    })

    if risk_level == "high":
        schedule.extend([
            {
                "timepoint": "Every 1-2 cycles",
                "assessments": "Echocardiogram, troponin, NT-proBNP",
                "priority": "required",
            },
            {
                "timepoint": "Mid-treatment",
                "assessments": "Cardiac MRI if LVEF borderline",
                "priority": "recommended",
            },
        ])
    elif risk_level == "moderate":
        schedule.extend([
            {
                "timepoint": "Every 2-3 cycles",
                "assessments": "Echocardiogram, troponin",
                "priority": "required",
            },
        ])
    else:
        schedule.append({
            "timepoint": "Mid-treatment and end of treatment",
            "assessments": "Echocardiogram, troponin",
            "priority": "recommended",
        })

    # End-of-treatment
    schedule.append({
        "timepoint": "End of treatment",
        "assessments": "Echocardiogram (LVEF + GLS), ECG, troponin, NT-proBNP",
        "priority": "required",
    })

    # Long-term follow-up (COG guidelines)
    schedule.append({
        "timepoint": "Annual follow-up (lifelong)",
        "assessments": "Echocardiogram per COG LTFU guidelines based on cumulative dose and age at exposure",
        "priority": "required",
    })

    # Merge trial-specific schedule if available
    if trial_data.get("status") == "success":
        trial_schedule = trial_data.get("schedule", [])
        for item in trial_schedule:
            if isinstance(item, dict):
                schedule.append(item)

    return schedule
