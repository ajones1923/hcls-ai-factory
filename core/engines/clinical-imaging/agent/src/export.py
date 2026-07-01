"""Export agent responses to Markdown, JSON, PDF, FHIR R4, and DICOM SR formats."""

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

# SNOMED CT codes for common imaging findings (100+ codes)
SNOMED_FINDING_CODES: Dict[str, str] = {
    # --- General / Legacy ---
    "normal": "17621005",

    # --- Chest / Lung findings ---
    "pulmonary nodule": "427359005",
    "nodule": "427359005",
    "lung nodule": "427359005",
    "lung mass": "309067003",
    "ground-glass opacity": "713171002",
    "ground glass opacity": "713171002",
    "ggo": "713171002",
    "consolidation": "45199005",
    "pleural effusion": "60046008",
    "effusion": "60046008",
    "pneumothorax": "36118008",
    "pulmonary embolism": "59282003",
    "pe": "59282003",
    "atelectasis": "46621007",
    "emphysema": "87433001",
    "bronchiectasis": "12295008",
    "pulmonary fibrosis": "51615001",
    "lung fibrosis": "51615001",
    "mediastinal lymphadenopathy": "274744001",
    "lymphadenopathy": "274744001",
    "cardiomegaly": "8186001",
    "enlarged heart": "8186001",
    "pericardial effusion": "373945007",

    # --- Head / Brain findings ---
    "intracranial hemorrhage": "1386000",
    "hemorrhage": "1386000",
    "ich": "1386000",
    "subdural hematoma": "35486000",
    "sdh": "35486000",
    "epidural hematoma": "75063001",
    "edh": "75063001",
    "subarachnoid hemorrhage": "21454007",
    "sah": "21454007",
    "cerebral infarction": "432504007",
    "cerebral infarct": "432504007",
    "stroke": "432504007",
    "midline shift": "89694002",
    "hydrocephalus": "230745008",
    "brain tumor": "126952004",
    "brain mass": "126952004",
    "intracranial mass": "126952004",
    "white matter lesion": "38341003",
    "white matter disease": "38341003",
    "leukoaraiosis": "38341003",

    # --- Abdominal findings ---
    "hepatic lesion": "300332007",
    "liver lesion": "300332007",
    "hepatocellular carcinoma": "25370001",
    "hcc": "25370001",
    "pancreatic mass": "126859007",
    "pancreatic lesion": "126859007",
    "renal mass": "126877002",
    "renal lesion": "126877002",
    "kidney mass": "126877002",
    "splenomegaly": "16294009",
    "enlarged spleen": "16294009",
    "ascites": "389026000",
    "bowel obstruction": "81060008",
    "small bowel obstruction": "81060008",
    "sbo": "81060008",
    "appendicitis": "74400008",
    "cholecystitis": "76581006",
    "abdominal aortic aneurysm": "233985008",
    "aaa": "233985008",

    # --- Breast findings ---
    "breast mass": "290077003",
    "breast lesion": "290077003",
    "breast calcification": "129748006",
    "breast calcifications": "129748006",
    "architectural distortion": "129752003",
    "breast density": "129747001",
    "dense breast": "129747001",

    # --- Thyroid findings ---
    "thyroid nodule": "237495002",
    "thyroid lesion": "237495002",
    "thyroid calcification": "396353007",

    # --- Musculoskeletal ---
    "fracture": "125605004",
    "bone fracture": "125605004",
    "degenerative joint disease": "396275006",
    "osteoarthritis": "396275006",
    "djd": "396275006",
    "disc herniation": "76107001",
    "herniated disc": "76107001",
    "disk herniation": "76107001",
    "compression fracture": "207957008",
    "vertebral compression fracture": "207957008",

    # --- Vascular ---
    "stenosis": "415582006",
    "vascular stenosis": "415582006",
    "aneurysm": "432119003",
    "dissection": "308546005",
    "vascular dissection": "308546005",
    "aortic dissection": "308546005",
    "thrombosis": "64572001",
    "thrombus": "64572001",
    "dvt": "64572001",

    # --- General findings ---
    "mass": "4147007",
    "lesion": "4147007",
    "mass nos": "4147007",
    "calcification": "50960005",
    "calcification nos": "50960005",
    "edema": "423666004",
    "inflammation": "257552002",
    "inflammatory": "257552002",
    "necrosis": "6574001",
    "necrotic": "6574001",
    "fibrosis": "263756000",
}

# Severity -> FHIR Observation Interpretation code
SEVERITY_INTERPRETATION: Dict[str, str] = {
    "critical": "HH",
    "urgent": "H",
    "significant": "A",
    "routine": "N",
    "normal": "N",
}

# ═══════════════════════════════════════════════════════════════════════
# NVIDIA BRANDING CONSTANTS
# ═══════════════════════════════════════════════════════════════════════

NVIDIA_GREEN = (118, 185, 0)  # RGB
NVIDIA_DARK = (30, 30, 30)
HEADER_HEIGHT = 50
PAGE_MARGIN = 40

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
        f"# Imaging Intelligence Report\n",
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
        md.append(f"\n## Workflow Results\n")
        for wr in response.workflow_results:
            md.append(f"\n### {wr.workflow_name}\n")
            md.append(f"- **Status:** {wr.status.value if hasattr(wr.status, 'value') else wr.status}\n")
            md.append(f"- **Severity:** {wr.severity.value if hasattr(wr.severity, 'value') else wr.severity}\n")
            if wr.classification:
                md.append(f"- **Classification:** {wr.classification}\n")
            if wr.findings:
                md.append(f"- **Findings:**\n")
                for finding in wr.findings:
                    desc = finding.get("description", str(finding))
                    md.append(f"  - {desc}\n")
            if wr.measurements:
                md.append(f"- **Measurements:**\n")
                for key, value in wr.measurements.items():
                    md.append(f"  - {key}: {value}\n")

    if response.nim_services_used:
        md.append(f"\n## NVIDIA NIM Services Used\n")
        for nim in response.nim_services_used:
            md.append(f"- {nim}\n")

    md.append(f"\n---\n*Research use only.*\n")
    return "".join(md)


def export_json(response: AgentResponse) -> str:
    """Export response as JSON string."""
    return response.model_dump_json(indent=2)


def _clean_markdown(text: str) -> str:
    """Strip common markdown formatting for PDF rendering."""
    import re
    # Bold **text** or __text__ -> <b>text</b>
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'__(.+?)__', r'<b>\1</b>', text)
    # Italic *text* or _text_ -> <i>text</i>
    text = re.sub(r'\*(.+?)\*', r'<i>\1</i>', text)
    text = re.sub(r'(?<!\w)_(.+?)_(?!\w)', r'<i>\1</i>', text)
    # Inline code `text` -> text
    text = re.sub(r'`(.+?)`', r'\1', text)
    # Headings ### text -> text
    text = re.sub(r'^#{1,6}\s+', '', text, flags=re.MULTILINE)
    # Bullet points - item -> item
    text = re.sub(r'^\s*[-*]\s+', '  \u2022 ', text, flags=re.MULTILINE)
    return text


def _severity_color(severity_str: str):
    """Return a reportlab color for a finding severity level."""
    from reportlab.lib import colors
    severity_map = {
        "critical": colors.Color(0.8, 0.0, 0.0),    # red
        "urgent": colors.Color(0.9, 0.4, 0.0),      # orange
        "significant": colors.Color(0.8, 0.6, 0.0),  # amber
        "routine": colors.Color(0.2, 0.6, 0.2),      # green
        "normal": colors.Color(0.3, 0.3, 0.3),       # gray
    }
    return severity_map.get(severity_str.lower(), colors.Color(0.3, 0.3, 0.3))


def export_pdf(response: AgentResponse, output_path: str) -> str:
    """Export response as NVIDIA-branded PDF file.

    Generates a professional PDF report with NVIDIA green header bar,
    structured evidence tables, workflow results, and NIM service listing.
    Style matches the Precision Oncology Agent PDF output.

    Args:
        response: The AgentResponse to render.
        output_path: File path for the generated PDF.

    Returns:
        The output file path on success, or empty string on failure.
    """
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        )
        from reportlab.platypus.flowables import HRFlowable
    except ImportError:
        logger.error("reportlab not installed - cannot export PDF. "
                      "Install with: pip install reportlab")
        return ""

    # --- Brand colors ---
    nvidia_green = colors.Color(
        NVIDIA_GREEN[0] / 255, NVIDIA_GREEN[1] / 255, NVIDIA_GREEN[2] / 255
    )
    nvidia_dark = colors.Color(
        NVIDIA_DARK[0] / 255, NVIDIA_DARK[1] / 255, NVIDIA_DARK[2] / 255
    )

    # --- Document setup ---
    doc = SimpleDocTemplate(
        output_path,
        pagesize=letter,
        leftMargin=PAGE_MARGIN,
        rightMargin=PAGE_MARGIN,
        topMargin=PAGE_MARGIN,
        bottomMargin=PAGE_MARGIN,
    )

    styles = getSampleStyleSheet()

    # --- Custom paragraph styles ---
    styles.add(ParagraphStyle(
        "NVTitle",
        parent=styles["Title"],
        textColor=colors.white,
        fontSize=20,
        leading=24,
        spaceAfter=6,
    ))
    styles.add(ParagraphStyle(
        "NVHeading",
        parent=styles["Heading2"],
        textColor=nvidia_dark,
        fontSize=14,
        leading=18,
        spaceBefore=16,
        spaceAfter=8,
    ))
    styles.add(ParagraphStyle(
        "NVBody",
        parent=styles["BodyText"],
        fontSize=10,
        leading=13,
    ))
    styles.add(ParagraphStyle(
        "NVCode",
        parent=styles["BodyText"],
        fontName="Courier",
        fontSize=8,
        leading=10,
        backColor=colors.Color(0.95, 0.95, 0.95),
    ))
    styles.add(ParagraphStyle(
        "NVDisclaimer",
        parent=styles["BodyText"],
        fontSize=7,
        leading=9,
        textColor=colors.gray,
    ))

    elements = []

    # ── 1. Green Header Bar ──────────────────────────────────────────
    header_data = [[Paragraph("Imaging Intelligence Report", styles["NVTitle"])]]
    header_table = Table(
        header_data,
        colWidths=[7.3 * inch],
        rowHeights=[HEADER_HEIGHT],
    )
    header_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), nvidia_green),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 12),
    ]))
    elements.append(header_table)
    elements.append(Spacer(1, 12))

    # ── 2. Report Metadata ───────────────────────────────────────────
    meta_lines = [
        f"<b>Query:</b> {response.question}",
        f"<b>Timestamp:</b> {response.timestamp}",
        f"<b>Collections Searched:</b> {response.evidence.total_collections_searched}",
        f"<b>Total Evidence Items:</b> {response.evidence.hit_count}",
    ]
    for ml in meta_lines:
        elements.append(Paragraph(ml, styles["NVBody"]))
    elements.append(Spacer(1, 8))

    # ── 3. Analysis Section ──────────────────────────────────────────
    elements.append(Paragraph("Analysis", styles["NVHeading"]))
    for para in response.answer.split("\n\n"):
        cleaned = para.strip()
        if cleaned:
            cleaned = _clean_markdown(cleaned)
            elements.append(Paragraph(cleaned, styles["NVBody"]))
            elements.append(Spacer(1, 6))
    elements.append(Spacer(1, 8))

    # ── 4. Evidence Table ────────────────────────────────────────────
    hits_by_coll = response.evidence.hits_by_collection()
    if response.evidence.hit_count > 0 and hits_by_coll:
        elements.append(Paragraph("Evidence", styles["NVHeading"]))
        table_data = [["Collection", "ID", "Score", "Summary"]]
        for collection, hits in hits_by_coll.items():
            for hit in hits[:10]:
                summary_text = hit.text[:150] + "..." if len(hit.text) > 150 else hit.text
                summary_text = _clean_markdown(summary_text)
                table_data.append([
                    str(collection),
                    str(hit.id),
                    f"{hit.score:.3f}",
                    Paragraph(summary_text, styles["NVBody"]),
                ])

        col_widths = [1.2 * inch, 1.0 * inch, 0.7 * inch, 4.4 * inch]
        t = Table(table_data, colWidths=col_widths, repeatRows=1)
        t.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), nvidia_green),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
            ("FONTSIZE", (0, 0), (-1, 0), 9),
            ("FONTSIZE", (0, 1), (-1, -1), 8),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.whitesmoke, colors.white]),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("LEFTPADDING", (0, 0), (-1, -1), 4),
            ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ]))
        elements.append(t)
        elements.append(Spacer(1, 10))

    # ── 5. Workflow Results ──────────────────────────────────────────
    if response.workflow_results:
        elements.append(Paragraph("Workflow Results", styles["NVHeading"]))
        for wr in response.workflow_results:
            severity_str = (
                wr.severity.value if hasattr(wr.severity, "value") else str(wr.severity)
            )
            status_str = (
                wr.status.value if hasattr(wr.status, "value") else str(wr.status)
            )
            elements.append(Paragraph(
                f"<b>{wr.workflow_name}</b> &mdash; Status: {status_str} | "
                f"Severity: {severity_str}",
                styles["NVBody"],
            ))
            if wr.classification:
                elements.append(Paragraph(
                    f"<b>Classification:</b> {wr.classification}",
                    styles["NVBody"],
                ))

            # Findings table with severity color coding
            if wr.findings:
                findings_data = [["Finding", "Category", "Severity"]]
                for finding in wr.findings:
                    desc = finding.get("description", str(finding))
                    cat = finding.get("category", "")
                    sev = finding.get("severity", severity_str)
                    findings_data.append([
                        Paragraph(_clean_markdown(desc), styles["NVBody"]),
                        str(cat),
                        str(sev),
                    ])

                ft = Table(findings_data, colWidths=[4.0 * inch, 1.5 * inch, 1.8 * inch], repeatRows=1)
                sev_colors = []
                for row_idx in range(1, len(findings_data)):
                    sev_val = str(findings_data[row_idx][2]).lower()
                    sev_colors.append(
                        ("TEXTCOLOR", (2, row_idx), (2, row_idx), _severity_color(sev_val))
                    )
                ft.setStyle(TableStyle([
                    ("BACKGROUND", (0, 0), (-1, 0), nvidia_green),
                    ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                    ("FONTSIZE", (0, 0), (-1, 0), 9),
                    ("FONTSIZE", (0, 1), (-1, -1), 8),
                    ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
                    ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.whitesmoke, colors.white]),
                    ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                    ("LEFTPADDING", (0, 0), (-1, -1), 4),
                    ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ] + sev_colors))
                elements.append(ft)
                elements.append(Spacer(1, 6))

            # Measurements as key-value pairs
            if wr.measurements:
                elements.append(Paragraph("<b>Measurements:</b>", styles["NVBody"]))
                for key, value in wr.measurements.items():
                    elements.append(Paragraph(
                        f"&nbsp;&nbsp;&bull; {key}: {value}",
                        styles["NVBody"],
                    ))
            elements.append(Spacer(1, 8))

    # ── 6. NIM Services Section ──────────────────────────────────────
    if response.nim_services_used:
        elements.append(Paragraph("NVIDIA NIM Services Used", styles["NVHeading"]))
        for nim in response.nim_services_used:
            elements.append(Paragraph(f"&bull; {nim}", styles["NVBody"]))
        elements.append(Spacer(1, 10))

    # ── 7. Footer ────────────────────────────────────────────────────
    elements.append(Spacer(1, 20))
    elements.append(HRFlowable(width="100%", color=colors.grey, thickness=0.5))
    elements.append(Spacer(1, 4))
    elements.append(Paragraph(
        "Research use only. Not for clinical decision-making. "
        "All findings require radiologist review.",
        styles["NVDisclaimer"],
    ))

    doc.build(elements)
    logger.info(f"PDF exported to {output_path}")
    return output_path


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
                    "system": "urn:oid:imaging-intelligence-engine",
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


# ═══════════════════════════════════════════════════════════════════════
# DICOM STRUCTURED REPORT EXPORT
# ═══════════════════════════════════════════════════════════════════════


def export_dicom_sr(
    workflow_result: WorkflowResult,
    output_path: str = "output.dcm",
    **kwargs,
) -> Dict[str, Any]:
    """Export workflow result as DICOM Structured Report.

    Delegates to DICOMSRExporter to generate a TID 1500 Measurement
    Report. The resulting .dcm file is viewable in OHIF and retrievable
    via standard DICOM C-FIND queries.

    Args:
        workflow_result: The WorkflowResult to export.
        output_path: File path for the generated .dcm SR.
        **kwargs: Additional arguments passed to
            DICOMSRExporter.export_workflow_result() (e.g., patient_id,
            patient_name, source_dicom_path, study_date,
            accession_number).

    Returns:
        Dict with keys: success, output_path, sop_instance_uid, error,
        content_items_count.
    """
    from src.export_dicom_sr import DICOMSRExporter, DICOMSRResult

    exporter = DICOMSRExporter()
    result: DICOMSRResult = exporter.export_workflow_result(
        workflow_result=workflow_result,
        output_path=output_path,
        **kwargs,
    )
    logger.info(
        f"DICOM SR export {'succeeded' if result.success else 'failed'}: "
        f"{result.output_path or result.error}"
    )
    return result.model_dump()
