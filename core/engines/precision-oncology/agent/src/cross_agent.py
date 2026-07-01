"""Cross-agent integration for the Precision Oncology Intelligence Agent.

Provides HTTP-based query functions to consult other HCLS AI Factory
intelligence agents and integrate their results into a unified pediatric
oncology assessment.  The module is designed around the clinical workflow
for pediatric precision oncology where multiple disciplines (immunotherapy,
pharmacogenomics, cardiology, neurology, imaging, single-cell analysis)
must converge to produce a safe, individualised treatment plan.

Supported cross-agent queries:
  - query_cart_agent()         -- CAR-T eligibility for hematologic malignancies
  - query_biomarker_agent()    -- biomarker panel enrichment / risk stratification
  - query_trial_agent()        -- precision medicine trial matching
  - query_cardiology_agent()   -- anthracycline cardiotoxicity risk
  - query_neurology_agent()    -- methotrexate / vincristine neurotoxicity risk
  - query_pgx_agent()          -- pharmacogenomic interaction screening
  - query_imaging_agent()      -- staging imaging protocol recommendation
  - query_single_cell_agent()  -- TME profiling for immunotherapy assessment
  - integrate_cross_agent_results() -- unified multi-agent assessment

All functions degrade gracefully: if an agent is unavailable, a warning
is logged and a default response with status="unavailable" is returned.

Author: Adam Jones
Date: March 2026
License: Apache 2.0
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List

import httpx

from config.settings import settings

logger = logging.getLogger(__name__)


# ===================================================================
# HELPERS
# ===================================================================


def _unavailable_response(agent_name: str) -> Dict[str, Any]:
    """Return a standard unavailable response for a cross-agent query.

    Args:
        agent_name: Name of the unavailable agent.

    Returns:
        Dict with ``status`` set to ``"unavailable"`` and empty result lists.
    """
    return {
        "status": "unavailable",
        "agent": agent_name,
        "message": f"{agent_name} agent is not currently available",
        "recommendations": [],
        "warnings": [],
    }


def _post(
    url: str,
    payload: Dict[str, Any],
    timeout: float,
) -> Dict[str, Any]:
    """Issue a POST request via httpx and return the JSON body.

    Args:
        url: Fully-qualified endpoint URL.
        payload: JSON-serialisable request body.
        timeout: Read/connect timeout in seconds.

    Returns:
        Parsed JSON response dict.

    Raises:
        httpx.HTTPStatusError: If the response status code >= 400.
        httpx.RequestError: On connection / timeout failures.
    """
    with httpx.Client(timeout=timeout) as client:
        resp = client.post(url, json=payload)
        resp.raise_for_status()
        return resp.json()


# ===================================================================
# CROSS-AGENT QUERY FUNCTIONS
# ===================================================================


def query_cart_agent(
    target_antigens: List[str],
    patient_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the CAR-T Intelligence Agent for immunotherapy eligibility.

    Clinical rationale: Pediatric patients with relapsed/refractory
    B-cell acute lymphoblastic leukemia (B-ALL) may be eligible for
    CAR-T therapy targeting CD19 or CD22.  Early cross-referencing with
    the CAR-T agent enables timely referral before the patient
    progresses beyond eligibility criteria.  The agent evaluates antigen
    expression levels, prior therapy lines, and disease status to
    determine candidacy for approved products (e.g. tisagenlecleucel)
    or investigational constructs.

    Args:
        target_antigens: Antigen targets to evaluate (e.g. ["CD19", "CD22"]).
        patient_data: Patient demographics, diagnosis, prior therapies,
            disease status, and antigen expression data.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``eligibility``, ``recommended_constructs``,
        ``warnings``, and ``recommendations``.
    """
    try:
        diagnosis = patient_data.get("diagnosis", "")
        data = _post(
            f"{settings.CART_AGENT_URL}/api/query",
            payload={
                "question": (
                    f"Evaluate CAR-T eligibility for {diagnosis} "
                    f"targeting {', '.join(target_antigens)}"
                ),
                "patient_context": {
                    "target_antigens": target_antigens,
                    "diagnosis": diagnosis,
                    "prior_therapies": patient_data.get("prior_therapies", []),
                    "disease_status": patient_data.get("disease_status", ""),
                    "age": patient_data.get("age"),
                    "antigen_expression": patient_data.get("antigen_expression", {}),
                },
            },
            timeout=timeout,
        )
        return {
            "status": "success",
            "agent": "cart",
            "eligibility": data.get("eligibility", {}),
            "recommended_constructs": data.get("recommended_constructs", []),
            "warnings": data.get("warnings", []),
            "recommendations": data.get("recommendations", []),
            "confidence": data.get("confidence", 0.0),
        }
    except Exception as exc:
        logger.warning("CAR-T agent query failed: %s", exc)
        return _unavailable_response("cart")


def query_biomarker_agent(
    tumor_profile: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Biomarker Intelligence Agent for risk-stratification enrichment.

    Clinical rationale: Pediatric oncology risk stratification relies on
    an expanding panel of molecular biomarkers (e.g. MYCN amplification
    in neuroblastoma, MRD in ALL, 1p/19q co-deletion in CNS tumours).
    The biomarker agent enriches the tumour profile with validated
    prognostic and predictive markers, enabling refined risk-group
    assignment (low / intermediate / high) that directly impacts
    therapy intensity decisions.

    Args:
        tumor_profile: Tumour genomic and phenotypic data including
            cancer type, known biomarkers, histology, and grade.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``enriched_biomarkers``,
        ``risk_stratification``, ``prognostic_markers``, and
        ``recommendations``.
    """
    try:
        cancer_type = tumor_profile.get("cancer_type", "")
        data = _post(
            f"{settings.BIOMARKER_AGENT_URL}/api/query",
            payload={
                "question": (
                    f"Enrich biomarker panel for risk stratification "
                    f"in pediatric {cancer_type}"
                ),
                "tumor_context": {
                    "cancer_type": cancer_type,
                    "known_biomarkers": tumor_profile.get("biomarkers", []),
                    "histology": tumor_profile.get("histology", ""),
                    "grade": tumor_profile.get("grade", ""),
                    "genomic_variants": tumor_profile.get("genomic_variants", []),
                },
            },
            timeout=timeout,
        )
        return {
            "status": "success",
            "agent": "biomarker",
            "enriched_biomarkers": data.get("enriched_biomarkers", []),
            "risk_stratification": data.get("risk_stratification", {}),
            "prognostic_markers": data.get("prognostic_markers", []),
            "recommendations": data.get("recommendations", []),
            "warnings": data.get("warnings", []),
        }
    except Exception as exc:
        logger.warning("Biomarker agent query failed: %s", exc)
        return _unavailable_response("biomarker")


def query_trial_agent(
    patient_profile: Dict[str, Any],
    mutations: List[Dict[str, Any]],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Clinical Trial Intelligence Agent for precision-medicine trial matches.

    Clinical rationale: Pediatric oncology patients -- especially those
    with relapsed/refractory disease -- benefit disproportionately from
    access to molecularly-guided clinical trials.  This function sends
    the patient's mutation profile and clinical context to the trial
    agent, which matches against its curated ClinicalTrials.gov index
    to surface open, enrolment-active studies with relevant molecular
    eligibility criteria.

    Args:
        patient_profile: Patient demographics, diagnosis, disease status,
            and prior therapies.
        mutations: List of mutation dicts with gene, variant, consequence,
            and actionability fields.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``matched_trials``, ``recommendations``,
        and ``warnings``.
    """
    try:
        diagnosis = patient_profile.get("diagnosis", "")
        gene_list = [m.get("gene", "") for m in mutations if m.get("gene")]
        data = _post(
            f"{settings.TRIAL_AGENT_URL}/api/query",
            payload={
                "question": (
                    f"Find precision medicine trials for pediatric {diagnosis} "
                    f"with mutations in {', '.join(gene_list)}"
                ),
                "patient_context": {
                    "diagnosis": diagnosis,
                    "age": patient_profile.get("age"),
                    "disease_status": patient_profile.get("disease_status", ""),
                    "prior_therapies": patient_profile.get("prior_therapies", []),
                    "mutations": mutations,
                },
            },
            timeout=timeout,
        )
        return {
            "status": "success",
            "agent": "clinical_trial",
            "matched_trials": data.get("matched_trials", []),
            "recommendations": data.get("recommendations", []),
            "warnings": data.get("warnings", []),
            "match_count": data.get("match_count", 0),
        }
    except Exception as exc:
        logger.warning("Clinical Trial agent query failed: %s", exc)
        return _unavailable_response("clinical_trial")


def query_cardiology_agent(
    therapy_plan: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Cardiology Intelligence Agent for anthracycline cardiotoxicity risk.

    Clinical rationale: Anthracyclines (doxorubicin, daunorubicin) are
    backbone agents in many pediatric oncology protocols (e.g. ALL,
    AML, Ewing sarcoma, osteosarcoma) but carry significant cumulative
    dose-dependent cardiotoxicity risk.  Children are especially
    vulnerable because they have decades of potential cardiac exposure
    ahead.  This function sends the planned regimen to the cardiology
    agent to obtain a risk assessment, recommended cumulative dose
    limits, and echocardiographic monitoring schedules.

    Args:
        therapy_plan: Planned therapy details including drug names,
            cumulative doses, schedule, and patient cardiac history.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``cardiac_risk``, ``risk_flags``,
        ``monitoring_schedule``, and ``recommendations``.
    """
    try:
        drugs = therapy_plan.get("drugs", [])
        drug_names = ", ".join(
            d.get("name", d) if isinstance(d, dict) else str(d)
            for d in drugs
        )
        data = _post(
            f"{settings.CARDIOLOGY_AGENT_URL}/api/query",
            payload={
                "question": (
                    f"Assess cardiotoxicity risk for pediatric regimen "
                    f"containing {drug_names}"
                ),
                "therapy_context": {
                    "drugs": drugs,
                    "cumulative_doses": therapy_plan.get("cumulative_doses", {}),
                    "schedule": therapy_plan.get("schedule", ""),
                    "cardiac_history": therapy_plan.get("cardiac_history", {}),
                    "patient_age": therapy_plan.get("patient_age"),
                },
            },
            timeout=timeout,
        )
        return {
            "status": "success",
            "agent": "cardiology",
            "cardiac_risk": data.get("cardiac_risk", {}),
            "risk_flags": data.get("risk_flags", []),
            "monitoring_schedule": data.get("monitoring_schedule", {}),
            "recommendations": data.get("recommendations", []),
            "warnings": data.get("warnings", []),
        }
    except Exception as exc:
        logger.warning("Cardiology agent query failed: %s", exc)
        return _unavailable_response("cardiology")


def query_neurology_agent(
    therapy_plan: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Neurology Intelligence Agent for CNS toxicity risk assessment.

    Clinical rationale: Methotrexate (intrathecal and high-dose IV) and
    vincristine are cornerstones of pediatric ALL and CNS-tumour
    protocols.  Methotrexate can cause acute and delayed
    leukoencephalopathy; vincristine is a dose-limiting peripheral
    neurotoxin.  In young children, neurotoxicity has lifelong
    developmental consequences.  This function requests a
    neurotoxicity risk profile so the oncology team can adjust dosing,
    select neuroprotective strategies, or plan neurocognitive
    monitoring.

    Args:
        therapy_plan: Planned therapy details including drug names,
            doses, routes of administration, and patient neurological
            history.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``neurotoxicity_risk``, ``risk_flags``,
        ``monitoring_recommendations``, and ``recommendations``.
    """
    try:
        drugs = therapy_plan.get("drugs", [])
        drug_names = ", ".join(
            d.get("name", d) if isinstance(d, dict) else str(d)
            for d in drugs
        )
        data = _post(
            f"{settings.NEUROLOGY_AGENT_URL}/api/query",
            payload={
                "question": (
                    f"Assess neurotoxicity risk for pediatric regimen "
                    f"containing {drug_names}"
                ),
                "therapy_context": {
                    "drugs": drugs,
                    "doses": therapy_plan.get("doses", {}),
                    "routes": therapy_plan.get("routes", {}),
                    "neurological_history": therapy_plan.get("neurological_history", {}),
                    "patient_age": therapy_plan.get("patient_age"),
                    "cns_involvement": therapy_plan.get("cns_involvement", False),
                },
            },
            timeout=timeout,
        )
        return {
            "status": "success",
            "agent": "neurology",
            "neurotoxicity_risk": data.get("neurotoxicity_risk", {}),
            "risk_flags": data.get("risk_flags", []),
            "monitoring_recommendations": data.get("monitoring_recommendations", []),
            "recommendations": data.get("recommendations", []),
            "warnings": data.get("warnings", []),
        }
    except Exception as exc:
        logger.warning("Neurology agent query failed: %s", exc)
        return _unavailable_response("neurology")


def query_pgx_agent(
    drug_list: List[str],
    patient_id: str,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Pharmacogenomics Intelligence Agent for drug-gene interactions.

    Clinical rationale: Pharmacogenomic variation significantly affects
    drug metabolism in pediatric patients.  For example, TPMT/NUDT15
    polymorphisms dictate thiopurine (6-MP) dosing in ALL; UGT1A1
    variants affect irinotecan toxicity; DPYD deficiency is
    life-threatening with fluoropyrimidines.  Screening the planned
    drug list against the patient's PGx profile prevents serious
    adverse drug reactions and enables dose individualisation before
    therapy begins.

    Args:
        drug_list: List of drug names planned for the patient.
        patient_id: Patient identifier for PGx profile lookup.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``pgx_results``, ``metabolizer_status``,
        ``dose_adjustments``, and ``warnings``.
    """
    try:
        data = _post(
            f"{settings.PGX_AGENT_URL}/api/query",
            payload={
                "question": (
                    f"Screen pharmacogenomic interactions for "
                    f"{', '.join(drug_list)}"
                ),
                "patient_context": {
                    "patient_id": patient_id,
                    "drugs": drug_list,
                },
            },
            timeout=timeout,
        )
        return {
            "status": "success",
            "agent": "pharmacogenomics",
            "pgx_results": data.get("pgx_results", []),
            "metabolizer_status": data.get("metabolizer_status", {}),
            "dose_adjustments": data.get("dose_adjustments", []),
            "warnings": data.get("warnings", []),
            "recommendations": data.get("recommendations", []),
        }
    except Exception as exc:
        logger.warning("PGx agent query failed: %s", exc)
        return _unavailable_response("pharmacogenomics")


def query_imaging_agent(
    cancer_type: str,
    stage: str,
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Imaging Intelligence Agent for staging protocol recommendations.

    Clinical rationale: Accurate staging drives risk-group assignment and
    therapy selection in pediatric oncology.  Different tumour types
    require distinct imaging protocols (e.g. MIBG + CT for
    neuroblastoma, MRI brain/spine for medulloblastoma, PET-CT for
    Hodgkin lymphoma).  The imaging agent recommends an evidence-based
    staging protocol based on tumour type and current stage, including
    modality, timing, and response-assessment criteria (e.g. RECIST,
    RANO, Lugano).

    Args:
        cancer_type: Cancer type or histological diagnosis.
        stage: Current or suspected disease stage.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``imaging_protocol``,
        ``recommended_modalities``, ``response_criteria``, and
        ``recommendations``.
    """
    try:
        data = _post(
            f"{settings.IMAGING_AGENT_URL}/api/query",
            payload={
                "question": (
                    f"Recommend staging imaging protocol for "
                    f"pediatric {cancer_type} stage {stage}"
                ),
                "clinical_context": {
                    "cancer_type": cancer_type,
                    "stage": stage,
                    "population": "pediatric",
                },
            },
            timeout=timeout,
        )
        return {
            "status": "success",
            "agent": "imaging",
            "imaging_protocol": data.get("imaging_protocol", {}),
            "recommended_modalities": data.get("recommended_modalities", []),
            "response_criteria": data.get("response_criteria", {}),
            "recommendations": data.get("recommendations", []),
            "warnings": data.get("warnings", []),
        }
    except Exception as exc:
        logger.warning("Imaging agent query failed: %s", exc)
        return _unavailable_response("imaging")


def query_single_cell_agent(
    tumor_data: Dict[str, Any],
    timeout: float = settings.CROSS_AGENT_TIMEOUT,
) -> Dict[str, Any]:
    """Query the Single-Cell Analysis Agent for tumour microenvironment profiling.

    Clinical rationale: The tumour microenvironment (TME) composition --
    particularly the abundance and functional state of tumour-infiltrating
    lymphocytes (TILs), myeloid-derived suppressor cells, and
    cancer-associated fibroblasts -- is a key determinant of
    immunotherapy response.  In pediatric solid tumours, which are
    often immunologically "cold", single-cell RNA-seq and spatial
    transcriptomics data from the single-cell agent can reveal whether
    checkpoint inhibitors, bispecific antibodies, or cellular therapies
    are likely to benefit the patient.

    Args:
        tumor_data: Tumour data including cancer type, sample ID,
            available -omics data types, and known immune markers.
        timeout: Request timeout in seconds.

    Returns:
        Dict with ``status``, ``tme_profile``, ``immune_infiltration``,
        ``immunotherapy_likelihood``, and ``recommendations``.
    """
    try:
        cancer_type = tumor_data.get("cancer_type", "")
        data = _post(
            f"{settings.SINGLE_CELL_AGENT_URL}/api/query",
            payload={
                "question": (
                    f"Profile tumour microenvironment for immunotherapy "
                    f"assessment in pediatric {cancer_type}"
                ),
                "tumor_context": {
                    "cancer_type": cancer_type,
                    "sample_id": tumor_data.get("sample_id", ""),
                    "available_omics": tumor_data.get("available_omics", []),
                    "immune_markers": tumor_data.get("immune_markers", {}),
                },
            },
            timeout=timeout,
        )
        return {
            "status": "success",
            "agent": "single_cell",
            "tme_profile": data.get("tme_profile", {}),
            "immune_infiltration": data.get("immune_infiltration", {}),
            "immunotherapy_likelihood": data.get("immunotherapy_likelihood", {}),
            "recommendations": data.get("recommendations", []),
            "warnings": data.get("warnings", []),
        }
    except Exception as exc:
        logger.warning("Single-cell agent query failed: %s", exc)
        return _unavailable_response("single_cell")


# ===================================================================
# INTEGRATION FUNCTION
# ===================================================================


def integrate_cross_agent_results(
    results: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Integrate results from multiple cross-agent queries into a unified assessment.

    Combines eligibility, risk assessments, trial matches, PGx warnings,
    cardiac/neurotoxicity flags, imaging protocols, and TME profiles into
    a single assessment suitable for the pediatric oncology tumour board.

    Args:
        results: List of cross-agent result dicts (from the query_* functions).

    Returns:
        Unified assessment dict with:
          - ``agents_consulted``: List of agent names queried.
          - ``agents_available``: List of agents that responded successfully.
          - ``combined_warnings``: Aggregated warnings from all agents.
          - ``combined_recommendations``: Aggregated recommendations.
          - ``safety_flags``: Combined safety concerns.
          - ``trial_matches``: Aggregated trial matches.
          - ``overall_assessment``: Summary assessment text.
    """
    agents_consulted: List[str] = []
    agents_available: List[str] = []
    combined_warnings: List[str] = []
    combined_recommendations: List[str] = []
    safety_flags: List[str] = []
    trial_matches: List[Dict[str, Any]] = []

    for result in results:
        agent = result.get("agent", "unknown")
        agents_consulted.append(agent)

        if result.get("status") != "success":
            continue

        agents_available.append(agent)

        # -- Warnings --
        for w in result.get("warnings", []):
            combined_warnings.append(f"[{agent}] {w}")

        # -- Recommendations --
        for r in result.get("recommendations", []):
            combined_recommendations.append(f"[{agent}] {r}")

        # -- Risk / safety flags --
        for f in result.get("risk_flags", []):
            safety_flags.append(f"[{agent}] {f}")

        # -- PGx high-impact items --
        for pgx in result.get("pgx_results", []):
            if isinstance(pgx, dict) and pgx.get("impact", "").lower() in (
                "high",
                "critical",
            ):
                combined_warnings.append(
                    f"[pharmacogenomics] {pgx.get('gene', '')}: "
                    f"{pgx.get('recommendation', '')}"
                )

        # -- Dose adjustments as warnings --
        for adj in result.get("dose_adjustments", []):
            if isinstance(adj, dict):
                combined_warnings.append(
                    f"[pharmacogenomics] Dose adjustment for "
                    f"{adj.get('drug', '')}: {adj.get('recommendation', '')}"
                )

        # -- Trial matches --
        trial_matches.extend(result.get("matched_trials", []))

    # -- Overall assessment --
    if not agents_available:
        overall = (
            "No cross-agent data available. Proceeding with oncology "
            "agent data only."
        )
    elif safety_flags:
        overall = (
            f"Cross-agent analysis identified {len(safety_flags)} safety "
            f"concern(s) from {len(agents_available)} agent(s). "
            f"Review by the tumour board is recommended before finalising "
            f"the treatment plan."
        )
    elif combined_warnings:
        overall = (
            f"Cross-agent analysis completed with {len(combined_warnings)} "
            f"warning(s) from {len(agents_available)} agent(s). All flagged "
            f"items should be reviewed."
        )
    else:
        overall = (
            f"Cross-agent analysis completed successfully. "
            f"{len(agents_available)} agent(s) consulted with no safety "
            f"concerns identified."
        )

    return {
        "agents_consulted": agents_consulted,
        "agents_available": agents_available,
        "combined_warnings": combined_warnings,
        "combined_recommendations": combined_recommendations,
        "safety_flags": safety_flags,
        "trial_matches": trial_matches,
        "overall_assessment": overall,
    }
