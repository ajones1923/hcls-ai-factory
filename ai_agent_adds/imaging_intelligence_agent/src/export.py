"""Export agent responses to Markdown, JSON, PDF, and FHIR R4 formats."""

import json
import uuid
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from loguru import logger

from src.models import AgentResponse, FindingSeverity, WorkflowResult


# ═══════════════════════════════════════════════════════════════════════
# FHIR R4 CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

FHIR_LOINC_SYSTEM = "http://loinc.org"
FHIR_SNOMED_SYSTEM = "http://snomed.info/sct"
FHIR_DICOM_SYSTEM = "http://dicom.nema.org/resources/ontology/DCM"
FHIR_INTERP_SYSTEM = (
    "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation"
)

# SNOMED CT codes for common imaging findings
SNOMED_FINDING_CODES: Dict[str, str] = {
    "hemorrhage": "50960005",
    "nodule": "416940007",
    "consolidation": "95436008",
    "pneumothorax": "36118008",
    "effusion": "60046008",
    "fracture": "125605004",
    "cardiomegaly": "8186001",
    "mass": "4147007",
    "edema": "267038008",
    "normal": "17621005",
}

# Severity -> FHIR Observation Interpretation code
SEVERITY_INTERPRETATION: Dict[str, str] = {
    "critical": "HH",
    "urgent": "H",
    "significant": "A",
    "routine": "N",
    "normal": "N",
}

# Modality keyword -> DICOM modality code
MODALITY_DICOM_CODES: Dict[str, str] = {
    "ct": "CT",
    "mri": "MR",
    "xray": "DX",
    "cxr": "CR",
    "ultrasound": "US",
    "pet": "PT",
    "pet_ct": "PT",
    "mammography": "MG",
    "fluoroscopy": "RF",
}


def export_markdown(response: AgentResponse) -> str:
    """Export response as Markdown string."""
    # (Similar to agent.generate_report but more detailed)
    md = [
        "# Imaging Intelligence Report\n",
        f"**Query:** {response.question}\n",
        f"**Timestamp:** {response.timestamp}\n",
        f"\n## Analysis\n\n{response.answer}\n",
        f"\n## Evidence ({response.evidence.hit_count} items)\n",
    ]

    for collection, hits in response.evidence.hits_by_collection().items():
        md.append(f"\n### {collection} ({len(hits)} results)\n")
        for hit in hits[:5]:
            md.append(f"- [{hit.id}] (score: {hit.score:.3f}) {hit.text[:200]}...\n")

    if response.workflow_results:
        md.append("\n## Workflow Results\n")
        for wr in response.workflow_results:
            md.append(f"\n### {wr.workflow_name}\n")
            md.append(f"- **Status:** {wr.status.value if hasattr(wr.status, 'value') else wr.status}\n")
            md.append(f"- **Severity:** {wr.severity.value if hasattr(wr.severity, 'value') else wr.severity}\n")
            if wr.classification:
                md.append(f"- **Classification:** {wr.classification}\n")
            if wr.findings:
                md.append("- **Findings:**\n")
                for finding in wr.findings:
                    desc = finding.get("description", str(finding))
                    md.append(f"  - {desc}\n")
            if wr.measurements:
                md.append("- **Measurements:**\n")
                for key, value in wr.measurements.items():
                    md.append(f"  - {key}: {value}\n")

    if response.nim_services_used:
        md.append("\n## NVIDIA NIM Services Used\n")
        for nim in response.nim_services_used:
            md.append(f"- {nim}\n")

    md.append("\n---\n*Research use only.*\n")
    return "".join(md)


def export_json(response: AgentResponse) -> str:
    """Export response as JSON string."""
    return response.model_dump_json(indent=2)


def export_pdf(response: AgentResponse, output_path: str) -> str:
    """Export response as PDF file."""
    try:
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer

        doc = SimpleDocTemplate(output_path, pagesize=letter)
        styles = getSampleStyleSheet()
        story = []

        story.append(Paragraph("Imaging Intelligence Report", styles["Title"]))
        story.append(Spacer(1, 12))
        story.append(Paragraph(f"<b>Query:</b> {response.question}", styles["Normal"]))
        story.append(Paragraph(f"<b>Timestamp:</b> {response.timestamp}", styles["Normal"]))
        story.append(Spacer(1, 12))
        story.append(Paragraph("Analysis", styles["Heading2"]))

        # Split answer into paragraphs
        for para in response.answer.split("\n\n"):
            if para.strip():
                story.append(Paragraph(para.strip(), styles["Normal"]))
                story.append(Spacer(1, 6))

        story.append(Spacer(1, 12))
        story.append(Paragraph(f"Evidence: {response.evidence.hit_count} items from {response.evidence.total_collections_searched} collections", styles["Normal"]))
        story.append(Spacer(1, 12))
        story.append(Paragraph("<i>Research use only. All findings require clinician review.</i>", styles["Normal"]))

        doc.build(story)
        logger.info(f"PDF exported to {output_path}")
        return output_path
    except ImportError:
        logger.error("reportlab not installed - cannot export PDF")
        return ""


# ═══════════════════════════════════════════════════════════════════════
# FHIR R4 EXPORT
# ═══════════════════════════════════════════════════════════════════════


def _make_fullurl() -> str:
    """Generate a FHIR-compliant urn:uuid for bundle entry fullUrl."""
    return f"urn:uuid:{uuid.uuid4()}"


def _snomed_coding(category: str) -> Dict[str, Any]:
    """Return a SNOMED CT coding dict for a finding category.

    Falls back to the generic 'Clinical finding' code 404684003
    if the category is not in the mapping.
    """
    code = SNOMED_FINDING_CODES.get(category.lower(), "404684003")
    display = category.replace("_", " ").title()
    return {
        "system": FHIR_SNOMED_SYSTEM,
        "code": code,
        "display": display,
    }


def _severity_to_interpretation(severity: str) -> Dict[str, Any]:
    """Map FindingSeverity string to FHIR v3-ObservationInterpretation coding."""
    code = SEVERITY_INTERPRETATION.get(severity.lower(), "N")
    display_map = {
        "HH": "Critical high",
        "H": "High",
        "A": "Abnormal",
        "N": "Normal",
    }
    return {
        "coding": [
            {
                "system": FHIR_INTERP_SYSTEM,
                "code": code,
                "display": display_map.get(code, "Normal"),
            }
        ]
    }


def _detect_modality(response: AgentResponse) -> Optional[str]:
    """Try to detect imaging modality from workflow names or question text."""
    combined = response.question.lower()
    for wr in response.workflow_results:
        combined += " " + wr.workflow_name.lower()
    for keyword, dicom_code in MODALITY_DICOM_CODES.items():
        if keyword in combined:
            return dicom_code
    return None


def _build_observation(
    finding: Dict[str, Any],
    measurements: Dict[str, float],
    report_ref: str,
    severity_str: str,
) -> Dict[str, Any]:
    """Build a FHIR Observation resource from a workflow finding."""
    obs_id = _make_fullurl()
    category_val = finding.get("category", "unknown")
    description = finding.get("description", "")
    severity = finding.get("severity", severity_str)

    # Build components for measurements
    components: List[Dict[str, Any]] = []
    for mkey, mval in measurements.items():
        # Determine unit from key suffix
        unit = ""
        unit_code = ""
        if mkey.endswith("_ml"):
            unit = "mL"
            unit_code = "mL"
        elif mkey.endswith("_mm"):
            unit = "mm"
            unit_code = "mm"
        elif mkey.endswith("_cm"):
            unit = "cm"
            unit_code = "cm"
        elif mkey.endswith("_hu"):
            unit = "HU"
            unit_code = "[hnsf'U]"

        component: Dict[str, Any] = {
            "code": {
                "text": mkey.replace("_", " "),
            },
        }
        if unit:
            component["valueQuantity"] = {
                "value": mval,
                "unit": unit,
                "system": "http://unitsofmeasure.org",
                "code": unit_code,
            }
        else:
            component["valueQuantity"] = {
                "value": mval,
            }
        components.append(component)

    observation: Dict[str, Any] = {
        "fullUrl": obs_id,
        "resource": {
            "resourceType": "Observation",
            "status": "final",
            "code": {
                "coding": [_snomed_coding(category_val)],
                "text": category_val,
            },
            "valueString": description,
            "interpretation": [_severity_to_interpretation(severity)],
        },
    }

    if components:
        observation["resource"]["component"] = components

    return observation


def export_fhir(
    response: AgentResponse,
    patient_id: str = "anonymous",
    practitioner_id: str = "AI-system",
) -> str:
    """Export response as a FHIR R4 DiagnosticReport Bundle (JSON string).

    Converts an AgentResponse (with optional WorkflowResults) into a
    valid FHIR R4 Bundle of type "collection" containing:
      - DiagnosticReport with LOINC category and SNOMED conclusionCodes
      - Observation resources for each workflow finding
      - ImagingStudy stub with DICOM modality
      - Patient stub

    Args:
        response: The AgentResponse to export.
        patient_id: Patient identifier (default "anonymous").
        practitioner_id: Performer identifier (default "AI-system").

    Returns:
        A JSON string containing the FHIR R4 Bundle.
    """
    now_iso = datetime.now(timezone.utc).isoformat()
    entries: List[Dict[str, Any]] = []

    # --- Patient resource stub ---
    patient_fullurl = _make_fullurl()
    entries.append({
        "fullUrl": patient_fullurl,
        "resource": {
            "resourceType": "Patient",
            "id": patient_id,
            "identifier": [
                {
                    "system": "urn:oid:imaging-intelligence-agent",
                    "value": patient_id,
                }
            ],
        },
    })

    # --- ImagingStudy stub ---
    imaging_study_fullurl = _make_fullurl()
    modality_code = _detect_modality(response)
    modality_list = []
    if modality_code:
        modality_list.append({
            "system": FHIR_DICOM_SYSTEM,
            "code": modality_code,
        })
    entries.append({
        "fullUrl": imaging_study_fullurl,
        "resource": {
            "resourceType": "ImagingStudy",
            "id": str(uuid.uuid4()),
            "status": "available",
            "subject": {"reference": patient_fullurl},
            "modality": modality_list,
            "description": response.question,
        },
    })

    # --- Observation resources for workflow findings ---
    observation_refs: List[str] = []
    all_finding_categories: List[str] = []

    for wr in response.workflow_results:
        severity_str = (
            wr.severity.value if hasattr(wr.severity, "value") else str(wr.severity)
        )
        if wr.findings:
            for finding in wr.findings:
                obs = _build_observation(
                    finding, wr.measurements, "", severity_str
                )
                entries.append(obs)
                observation_refs.append(obs["fullUrl"])
                cat = finding.get("category", "")
                if cat:
                    all_finding_categories.append(cat)
        else:
            # No findings -- create a single observation for the workflow itself
            obs = _build_observation(
                {"category": "normal", "description": f"{wr.workflow_name}: No findings"},
                wr.measurements,
                "",
                severity_str,
            )
            entries.append(obs)
            observation_refs.append(obs["fullUrl"])

    # --- DiagnosticReport ---
    report_fullurl = _make_fullurl()

    # conclusionCode: SNOMED codes for all unique finding categories
    conclusion_codes = []
    seen_codes = set()
    for cat in all_finding_categories:
        coding = _snomed_coding(cat)
        if coding["code"] not in seen_codes:
            seen_codes.add(coding["code"])
            conclusion_codes.append({"coding": [coding]})

    # Build extensions for cross-modal results
    extensions: List[Dict[str, Any]] = []
    if response.cross_modal:
        extensions.append({
            "url": "urn:imaging-intelligence:cross-modal-result",
            "valueString": response.cross_modal.enrichment_summary
            or response.cross_modal.trigger_reason,
        })

    report_resource: Dict[str, Any] = {
        "resourceType": "DiagnosticReport",
        "id": str(uuid.uuid4()),
        "status": "final",
        "category": [
            {
                "coding": [
                    {
                        "system": FHIR_LOINC_SYSTEM,
                        "code": "LP29684-5",
                        "display": "Radiology",
                    }
                ]
            }
        ],
        "code": {
            "coding": [
                {
                    "system": FHIR_LOINC_SYSTEM,
                    "code": "18748-4",
                    "display": "Diagnostic imaging study",
                }
            ],
            "text": "Imaging Intelligence Agent Report",
        },
        "subject": {"reference": patient_fullurl},
        "effectiveDateTime": response.timestamp,
        "issued": now_iso,
        "performer": [
            {
                "reference": f"Practitioner/{practitioner_id}",
                "display": practitioner_id,
            }
        ],
        "imagingStudy": [{"reference": imaging_study_fullurl}],
        "result": [{"reference": ref} for ref in observation_refs],
        "conclusion": response.answer,
    }

    if conclusion_codes:
        report_resource["conclusionCode"] = conclusion_codes

    if extensions:
        report_resource["extension"] = extensions

    report_entry = {
        "fullUrl": report_fullurl,
        "resource": report_resource,
    }
    entries.append(report_entry)

    # --- Assemble Bundle ---
    bundle: Dict[str, Any] = {
        "resourceType": "Bundle",
        "id": str(uuid.uuid4()),
        "type": "collection",
        "timestamp": now_iso,
        "entry": entries,
    }

    return json.dumps(bundle, indent=2)
