"""
Export module for Oncology Intelligence Agent.

Generates Markdown, JSON, PDF, and FHIR R4 exports for oncology reports
and Molecular Tumor Board (MTB) packets.

Author: Adam Jones
Date: February 2026
"""

import json
import logging
import uuid
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Default brand colors — can be overridden via OncoSettings
# (PDF_BRAND_COLOR_R/G/B environment variables with ONCO_ prefix)
NVIDIA_GREEN = (118, 185, 0)  # RGB for NVIDIA brand green
NVIDIA_DARK = (30, 30, 30)
HEADER_HEIGHT = 50
PAGE_MARGIN = 40


def _get_brand_color():
    """Get brand color from settings if available, else use default."""
    try:
        from config.settings import settings as _s
        return (_s.PDF_BRAND_COLOR_R, _s.PDF_BRAND_COLOR_G, _s.PDF_BRAND_COLOR_B)
    except Exception:
        return NVIDIA_GREEN

EVIDENCE_LEVEL_LABELS = {
    "level_1": "Level 1 — FDA-approved / Standard of Care",
    "level_2": "Level 2 — Clinical Evidence / Consensus",
    "level_3": "Level 3 — Case Reports / Early Trials",
    "level_4": "Level 4 — Preclinical / Biological Rationale",
    "level_R": "Level R — Resistance Evidence",
}

FHIR_LOINC_CODES = {
    "genomic_report": "81247-9",       # Master HL7 genomic report
    "gene_studied": "48018-6",         # Gene studied
    "variant": "69548-6",              # Genetic variant assessment
    "therapeutic_implication": "51969-4",  # Genetic analysis summary
    "tumor_mutation_burden": "94076-7",   # TMB
    "microsatellite_instability": "81695-9",  # MSI
}

FHIR_SNOMED_CANCER_CODES = {
    "nsclc": ("254637007", "Non-small cell lung cancer"),
    "breast": ("254837009", "Malignant neoplasm of breast"),
    "colorectal": ("363406005", "Malignant tumor of colon"),
    "melanoma": ("372244006", "Malignant melanoma"),
    "pancreatic": ("363418001", "Malignant tumor of pancreas"),
    "ovarian": ("363443007", "Malignant tumor of ovary"),
    "prostate": ("399068003", "Malignant tumor of prostate"),
    "glioblastoma": ("393563007", "Glioblastoma multiforme"),
    "aml": ("91861009", "Acute myeloid leukemia"),
    "cml": ("92818009", "Chronic myeloid leukemia"),
    "bladder": ("399326009", "Malignant neoplasm of urinary bladder"),
    "renal": ("363518003", "Malignant tumor of kidney"),
    "gastric": ("363349007", "Malignant tumor of stomach"),
    "hepatocellular": ("25370001", "Hepatocellular carcinoma"),
    "esophageal": ("363402007", "Malignant tumor of esophagus"),
    "thyroid": ("363478007", "Malignant tumor of thyroid gland"),
    "head_and_neck": ("255055008", "Malignant neoplasm of head and neck"),
    "sclc": ("254632001", "Small cell carcinoma of lung"),
    "cholangiocarcinoma": ("312104005", "Cholangiocarcinoma"),
    "endometrial": ("254878006", "Malignant neoplasm of endometrium"),
    "cervical": ("363354003", "Malignant tumor of cervix"),
    "sarcoma": ("424413001", "Sarcoma"),
    "mesothelioma": ("62061000", "Mesothelioma"),
}


# ---------------------------------------------------------------------------
# Helper: normalise input to a dict
# ---------------------------------------------------------------------------

def _normalise_input(mtb_packet_or_response: Any) -> Dict:
    """Accept an MTBPacket model, a dict, or a string and return a dict."""
    if isinstance(mtb_packet_or_response, dict):
        return mtb_packet_or_response
    if isinstance(mtb_packet_or_response, str):
        try:
            return json.loads(mtb_packet_or_response)
        except json.JSONDecodeError:
            return {"raw_text": mtb_packet_or_response}
    # Pydantic model or dataclass — attempt .dict() / .model_dump()
    for method in ("model_dump", "dict", "__dict__"):
        fn = getattr(mtb_packet_or_response, method, None)
        if callable(fn):
            return fn()
        if isinstance(fn, dict):
            return fn
    return {"raw_text": str(mtb_packet_or_response)}


def _safe_get(data: Dict, *keys, default=""):
    """Nested dict safe-get."""
    current = data
    for k in keys:
        if isinstance(current, dict):
            current = current.get(k, default)
        else:
            return default
    return current


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ===================================================================
# 1. Markdown Export
# ===================================================================

def export_markdown(
    mtb_packet_or_response: Any,
    title: Optional[str] = None,
) -> str:
    """
    Generate a formatted Markdown report suitable for Molecular Tumor Board
    review or clinical documentation.

    Sections:
        - Header / patient summary
        - Variant table
        - Evidence summary
        - Therapy ranking
        - Clinical trial matches
        - Open questions / follow-up

    Parameters
    ----------
    mtb_packet_or_response : dict | MTBPacket | str
        The MTB packet or agent response to render.
    title : str, optional
        Override report title.

    Returns
    -------
    str
        Markdown-formatted report.
    """
    data = _normalise_input(mtb_packet_or_response)

    # If this is a plain-text response, wrap it minimally
    if "raw_text" in data and len(data) == 1:
        return f"# Oncology Intelligence Report\n\n{data['raw_text']}\n"

    report_title = title or data.get("title", "Oncology Intelligence — Molecular Tumor Board Report")
    lines: List[str] = []

    # --- Header ---
    lines.append(f"# {report_title}")
    lines.append("")
    lines.append(f"**Generated:** {_timestamp()}")
    lines.append("**Pipeline:** HCLS AI Factory — Oncology Intelligence Agent")
    if data.get("patient_id"):
        lines.append(f"**Patient ID:** {data['patient_id']}")
    if data.get("cancer_type"):
        lines.append(f"**Cancer Type:** {data['cancer_type']}")
    if data.get("sample_id"):
        lines.append(f"**Sample ID:** {data['sample_id']}")
    lines.append("")

    # --- Summary ---
    summary = data.get("summary") or data.get("clinical_summary", "")
    if summary:
        lines.append("## Clinical Summary")
        lines.append("")
        lines.append(str(summary))
        lines.append("")

    # --- Variant Table ---
    variants = data.get("variants") or data.get("somatic_variants") or []
    if variants:
        lines.append("## Somatic Variant Profile")
        lines.append("")
        lines.append("| Gene | Variant | Type | VAF | Consequence | Tier |")
        lines.append("|------|---------|------|-----|-------------|------|")
        for v in variants:
            gene = v.get("gene", v.get("gene_symbol", ""))
            variant_name = v.get("variant_name", v.get("hgvs", v.get("variant", "")))
            vtype = v.get("variant_type", v.get("type", ""))
            vaf = v.get("vaf", v.get("allele_frequency", ""))
            if isinstance(vaf, float):
                vaf = f"{vaf:.2%}"
            consequence = v.get("consequence", v.get("effect", ""))
            tier = v.get("tier", v.get("evidence_level", ""))
            lines.append(f"| {gene} | {variant_name} | {vtype} | {vaf} | {consequence} | {tier} |")
        lines.append("")

    # --- Biomarkers ---
    biomarkers = data.get("biomarkers", {})
    if biomarkers:
        lines.append("## Biomarker Summary")
        lines.append("")
        tmb = biomarkers.get("tmb") or biomarkers.get("tumor_mutation_burden")
        msi = biomarkers.get("msi") or biomarkers.get("microsatellite_instability")
        pdl1 = biomarkers.get("pdl1") or biomarkers.get("pd_l1_expression")
        if tmb is not None:
            lines.append(f"- **Tumor Mutation Burden (TMB):** {tmb} mut/Mb")
        if msi is not None:
            lines.append(f"- **Microsatellite Instability (MSI):** {msi}")
        if pdl1 is not None:
            lines.append(f"- **PD-L1 Expression:** {pdl1}")
        for k, v in biomarkers.items():
            if k not in ("tmb", "tumor_mutation_burden", "msi",
                         "microsatellite_instability", "pdl1", "pd_l1_expression"):
                lines.append(f"- **{k.replace('_', ' ').title()}:** {v}")
        lines.append("")

    # --- Evidence ---
    evidence_items = data.get("evidence") or data.get("evidence_items") or []
    if evidence_items:
        lines.append("## Evidence Summary")
        lines.append("")
        for idx, ev in enumerate(evidence_items, 1):
            ev_gene = ev.get("gene", "")
            ev_level = ev.get("evidence_level", ev.get("level", ""))
            ev_label = EVIDENCE_LEVEL_LABELS.get(ev_level, ev_level)
            ev_source = ev.get("source", ev.get("reference", ""))
            ev_text = ev.get("summary", ev.get("text", ev.get("description", "")))
            lines.append(f"### {idx}. {ev_gene} — {ev_label}")
            lines.append("")
            if ev_text:
                lines.append(ev_text)
            if ev_source:
                lines.append(f"\n*Source:* {ev_source}")
            lines.append("")

    # --- Therapy Ranking ---
    therapies = data.get("therapies") or data.get("therapy_ranking") or data.get("recommendations") or []
    if therapies:
        lines.append("## Therapy Ranking")
        lines.append("")
        lines.append("| Rank | Therapy | Target(s) | Evidence | Line | Notes |")
        lines.append("|------|---------|-----------|----------|------|-------|")
        for idx, tx in enumerate(therapies, 1):
            name = tx.get("name", tx.get("drug", tx.get("therapy", "")))
            targets = tx.get("targets", tx.get("target", ""))
            if isinstance(targets, list):
                targets = ", ".join(targets)
            ev_level = tx.get("evidence_level", tx.get("level", ""))
            line = tx.get("line_of_therapy", tx.get("line", ""))
            notes = tx.get("notes", tx.get("rationale", ""))
            lines.append(f"| {idx} | {name} | {targets} | {ev_level} | {line} | {notes} |")
        lines.append("")

    # --- Clinical Trial Matches ---
    trials = data.get("clinical_trials") or data.get("trial_matches") or []
    if trials:
        lines.append("## Clinical Trial Matches")
        lines.append("")
        for trial in trials:
            nct = trial.get("nct_id", trial.get("id", ""))
            trial_title = trial.get("title", "")
            phase = trial.get("phase", "")
            status = trial.get("status", "")
            match_rationale = trial.get("match_rationale", trial.get("rationale", ""))
            lines.append(f"- **{nct}** — {trial_title}")
            if phase or status:
                lines.append(f"  - Phase: {phase} | Status: {status}")
            if match_rationale:
                lines.append(f"  - *Match rationale:* {match_rationale}")
        lines.append("")

    # --- Pathway context ---
    pathways = data.get("pathways") or data.get("pathway_context") or []
    if pathways:
        lines.append("## Pathway Context")
        lines.append("")
        for pw in pathways:
            pw_name = pw.get("name", pw.get("pathway", ""))
            pw_desc = pw.get("description", pw.get("summary", ""))
            lines.append(f"- **{pw_name}:** {pw_desc}")
        lines.append("")

    # --- Resistance Mechanisms ---
    resistance = data.get("resistance_mechanisms") or data.get("resistance") or []
    if resistance:
        lines.append("## Known Resistance Mechanisms")
        lines.append("")
        for rm in resistance:
            rm_name = rm.get("mechanism", rm.get("name", ""))
            rm_drug = rm.get("drug", "")
            rm_desc = rm.get("description", rm.get("summary", ""))
            lines.append(f"- **{rm_name}** (affects {rm_drug}): {rm_desc}")
        lines.append("")

    # --- Open Questions ---
    questions = data.get("open_questions") or data.get("follow_up") or []
    if questions:
        lines.append("## Open Questions / Follow-Up")
        lines.append("")
        for q in questions:
            if isinstance(q, str):
                lines.append(f"- {q}")
            elif isinstance(q, dict):
                lines.append(f"- {q.get('question', q.get('text', str(q)))}")
        lines.append("")

    # --- Disclaimer ---
    lines.append("---")
    lines.append("")
    lines.append("*This report is generated by the HCLS AI Factory Oncology Intelligence Agent "
                 "and is intended for research use only. Clinical decisions should be made "
                 "in consultation with qualified healthcare professionals.*")
    lines.append("")

    return "\n".join(lines)


# ===================================================================
# 2. JSON Export
# ===================================================================

def export_json(mtb_packet_or_response: Any) -> dict:
    """
    Export MTB packet or agent response as structured JSON.

    Parameters
    ----------
    mtb_packet_or_response : dict | MTBPacket | str
        The MTB packet or agent response to export.

    Returns
    -------
    dict
        Structured JSON-serialisable dictionary with standard fields.
    """
    data = _normalise_input(mtb_packet_or_response)

    export = {
        "meta": {
            "format": "hcls-ai-factory-oncology-report",
            "version": "1.0.0",
            "generated_at": _timestamp(),
            "pipeline": "Oncology Intelligence Agent",
            "author": "HCLS AI Factory",
        },
        "patient_id": data.get("patient_id"),
        "cancer_type": data.get("cancer_type"),
        "sample_id": data.get("sample_id"),
        "clinical_summary": data.get("summary") or data.get("clinical_summary"),
        "variants": data.get("variants") or data.get("somatic_variants") or [],
        "biomarkers": data.get("biomarkers", {}),
        "evidence": data.get("evidence") or data.get("evidence_items") or [],
        "therapy_ranking": data.get("therapies") or data.get("therapy_ranking") or data.get("recommendations") or [],
        "clinical_trials": data.get("clinical_trials") or data.get("trial_matches") or [],
        "pathways": data.get("pathways") or data.get("pathway_context") or [],
        "resistance_mechanisms": data.get("resistance_mechanisms") or data.get("resistance") or [],
        "open_questions": data.get("open_questions") or data.get("follow_up") or [],
    }

    # Remove None values at top level
    export = {k: v for k, v in export.items() if v is not None}

    logger.info("Exported JSON report with %d variants, %d therapies",
                len(export.get("variants", [])),
                len(export.get("therapy_ranking", [])))
    return export


# ===================================================================
# 3. PDF Export
# ===================================================================

def export_pdf(
    mtb_packet_or_response: Any,
    output_path: str,
) -> str:
    """
    Generate an NVIDIA-themed PDF report via ReportLab.

    Features:
        - NVIDIA green header bar with white title text
        - Structured tables for variants, therapies, trials
        - Evidence summary sections
        - Footer with timestamp and disclaimer

    Parameters
    ----------
    mtb_packet_or_response : dict | MTBPacket | str
        The MTB packet or agent response to render.
    output_path : str
        File path for the generated PDF.

    Returns
    -------
    str
        The output file path on success.

    Raises
    ------
    ImportError
        If ReportLab is not installed.
    """
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.platypus import (
            SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
            HRFlowable, PageBreak,
        )
    except ImportError:
        raise ImportError(
            "ReportLab is required for PDF export. "
            "Install with: pip install reportlab"
        )

    data = _normalise_input(mtb_packet_or_response)

    brand_rgb = _get_brand_color()
    nvidia_green = colors.Color(
        brand_rgb[0] / 255, brand_rgb[1] / 255, brand_rgb[2] / 255
    )
    nvidia_dark = colors.Color(
        NVIDIA_DARK[0] / 255, NVIDIA_DARK[1] / 255, NVIDIA_DARK[2] / 255
    )

    doc = SimpleDocTemplate(
        output_path,
        pagesize=letter,
        leftMargin=PAGE_MARGIN,
        rightMargin=PAGE_MARGIN,
        topMargin=PAGE_MARGIN,
        bottomMargin=PAGE_MARGIN,
    )

    styles = getSampleStyleSheet()

    # Custom styles
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
        "NVDisclaimer",
        parent=styles["BodyText"],
        fontSize=7,
        leading=9,
        textColor=colors.gray,
    ))

    elements: List = []

    # --- Green Header Bar ---
    report_title = data.get("title", "Oncology Intelligence — Molecular Tumor Board Report")
    header_data = [[Paragraph(report_title, styles["NVTitle"])]]
    header_table = Table(header_data, colWidths=[7.3 * inch], rowHeights=[HEADER_HEIGHT])
    header_table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), nvidia_green),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 12),
    ]))
    elements.append(header_table)
    elements.append(Spacer(1, 12))

    # --- Patient / Meta Info ---
    meta_lines = [f"<b>Generated:</b> {_timestamp()}"]
    if data.get("patient_id"):
        meta_lines.append(f"<b>Patient ID:</b> {data['patient_id']}")
    if data.get("cancer_type"):
        meta_lines.append(f"<b>Cancer Type:</b> {data['cancer_type']}")
    if data.get("sample_id"):
        meta_lines.append(f"<b>Sample ID:</b> {data['sample_id']}")
    for ml in meta_lines:
        elements.append(Paragraph(ml, styles["NVBody"]))
    elements.append(Spacer(1, 8))

    # --- Summary ---
    summary = data.get("summary") or data.get("clinical_summary", "")
    if summary:
        elements.append(Paragraph("Clinical Summary", styles["NVHeading"]))
        elements.append(Paragraph(str(summary), styles["NVBody"]))
        elements.append(Spacer(1, 8))

    # --- Variant Table ---
    variants = data.get("variants") or data.get("somatic_variants") or []
    if variants:
        elements.append(Paragraph("Somatic Variant Profile", styles["NVHeading"]))
        table_data = [["Gene", "Variant", "Type", "VAF", "Consequence", "Tier"]]
        for v in variants:
            gene = str(v.get("gene", v.get("gene_symbol", "")))
            variant_name = str(v.get("variant_name", v.get("hgvs", v.get("variant", ""))))
            vtype = str(v.get("variant_type", v.get("type", "")))
            vaf = v.get("vaf", v.get("allele_frequency", ""))
            if isinstance(vaf, float):
                vaf = f"{vaf:.2%}"
            consequence = str(v.get("consequence", v.get("effect", "")))
            tier = str(v.get("tier", v.get("evidence_level", "")))
            table_data.append([gene, variant_name, vtype, str(vaf), consequence, tier])

        t = Table(table_data, repeatRows=1)
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

    # --- Therapy Ranking Table ---
    therapies = (data.get("therapies") or data.get("therapy_ranking")
                 or data.get("recommendations") or [])
    if therapies:
        elements.append(Paragraph("Therapy Ranking", styles["NVHeading"]))
        table_data = [["Rank", "Therapy", "Target(s)", "Evidence", "Line", "Notes"]]
        for idx, tx in enumerate(therapies, 1):
            name = str(tx.get("name", tx.get("drug", tx.get("therapy", ""))))
            targets = tx.get("targets", tx.get("target", ""))
            if isinstance(targets, list):
                targets = ", ".join(targets)
            ev_level = str(tx.get("evidence_level", tx.get("level", "")))
            line = str(tx.get("line_of_therapy", tx.get("line", "")))
            notes = str(tx.get("notes", tx.get("rationale", "")))
            table_data.append([str(idx), name, str(targets), ev_level, line, notes])

        t = Table(table_data, repeatRows=1)
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

    # --- Clinical Trial Matches ---
    trials = data.get("clinical_trials") or data.get("trial_matches") or []
    if trials:
        elements.append(Paragraph("Clinical Trial Matches", styles["NVHeading"]))
        for trial in trials:
            nct = trial.get("nct_id", trial.get("id", ""))
            trial_title = trial.get("title", "")
            phase = trial.get("phase", "")
            status = trial.get("status", "")
            rationale = trial.get("match_rationale", trial.get("rationale", ""))
            trial_text = f"<b>{nct}</b> — {trial_title}"
            if phase or status:
                trial_text += f"<br/>Phase: {phase} | Status: {status}"
            if rationale:
                trial_text += f"<br/><i>Match: {rationale}</i>"
            elements.append(Paragraph(trial_text, styles["NVBody"]))
            elements.append(Spacer(1, 4))

    # --- Disclaimer Footer ---
    elements.append(Spacer(1, 20))
    elements.append(HRFlowable(width="100%", color=colors.grey, thickness=0.5))
    elements.append(Spacer(1, 4))
    elements.append(Paragraph(
        "This report is generated by the HCLS AI Factory Oncology Intelligence Agent "
        "and is intended for research use only. Clinical decisions should be made "
        "in consultation with qualified healthcare professionals.",
        styles["NVDisclaimer"],
    ))

    doc.build(elements)
    logger.info("PDF report written to %s", output_path)
    return output_path


# ===================================================================
# 4. FHIR R4 Export
# ===================================================================

def export_fhir_r4(
    mtb_packet: Any,
    patient_id: str,
) -> dict:
    """
    Export an MTB packet as a FHIR R4 Bundle.

    Resources created:
        - Patient
        - DiagnosticReport (master genomic report)
        - Observation (one per variant, plus TMB / MSI if present)

    Coding systems used:
        - SNOMED CT for cancer type
        - LOINC for genomic observations
        - HGVS for variant descriptions

    Parameters
    ----------
    mtb_packet : dict | MTBPacket | str
        The MTB packet to convert.
    patient_id : str
        Patient identifier for the FHIR Patient resource.

    Returns
    -------
    dict
        FHIR R4 Bundle (type=collection).
    """
    data = _normalise_input(mtb_packet)
    bundle_id = str(uuid.uuid4())
    timestamp = datetime.now(timezone.utc).isoformat()

    entries: List[Dict] = []

    # --- Patient Resource ---
    patient_resource_id = str(uuid.uuid4())
    patient_resource = {
        "resourceType": "Patient",
        "id": patient_resource_id,
        "identifier": [{
            "system": "urn:hcls-ai-factory:patient",
            "value": patient_id,
        }],
        "active": True,
    }
    entries.append({
        "fullUrl": f"urn:uuid:{patient_resource_id}",
        "resource": patient_resource,
    })

    # --- Observation Resources (one per variant) ---
    observation_ids: List[str] = []
    variants = data.get("variants") or data.get("somatic_variants") or []

    for variant in variants:
        obs_id = str(uuid.uuid4())
        observation_ids.append(obs_id)

        gene = variant.get("gene", variant.get("gene_symbol", "Unknown"))
        hgvs = variant.get("variant_name", variant.get("hgvs", variant.get("variant", "")))
        vaf = variant.get("vaf", variant.get("allele_frequency"))
        consequence = variant.get("consequence", variant.get("effect", ""))

        observation = {
            "resourceType": "Observation",
            "id": obs_id,
            "status": "final",
            "category": [{
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                    "code": "laboratory",
                    "display": "Laboratory",
                }],
            }],
            "code": {
                "coding": [{
                    "system": "http://loinc.org",
                    "code": FHIR_LOINC_CODES["variant"],
                    "display": "Genetic variant assessment",
                }],
            },
            "subject": {
                "reference": f"urn:uuid:{patient_resource_id}",
            },
            "effectiveDateTime": timestamp,
            "valueCodeableConcept": {
                "coding": [{
                    "system": "http://varnomen.hgvs.org",
                    "code": hgvs,
                    "display": f"{gene} {hgvs}",
                }],
                "text": f"{gene} {hgvs}",
            },
            "component": [],
        }

        # Gene studied component
        observation["component"].append({
            "code": {
                "coding": [{
                    "system": "http://loinc.org",
                    "code": FHIR_LOINC_CODES["gene_studied"],
                    "display": "Gene studied [ID]",
                }],
            },
            "valueCodeableConcept": {
                "coding": [{
                    "system": "http://www.genenames.org/geneId",
                    "display": gene,
                }],
                "text": gene,
            },
        })

        # VAF component
        if vaf is not None:
            vaf_value = vaf if isinstance(vaf, (int, float)) else None
            if vaf_value is not None:
                observation["component"].append({
                    "code": {
                        "coding": [{
                            "system": "http://loinc.org",
                            "code": "81258-6",
                            "display": "Allelic frequency [NFr]",
                        }],
                    },
                    "valueQuantity": {
                        "value": vaf_value,
                        "unit": "relative frequency",
                        "system": "http://unitsofmeasure.org",
                        "code": "1",
                    },
                })

        # Consequence component
        if consequence:
            observation["component"].append({
                "code": {
                    "coding": [{
                        "system": "http://loinc.org",
                        "code": "48006-1",
                        "display": "Molecular consequence type",
                    }],
                },
                "valueCodeableConcept": {
                    "text": consequence,
                },
            })

        entries.append({
            "fullUrl": f"urn:uuid:{obs_id}",
            "resource": observation,
        })

    # --- TMB Observation ---
    biomarkers = data.get("biomarkers", {})
    tmb = biomarkers.get("tmb") or biomarkers.get("tumor_mutation_burden")
    if tmb is not None:
        tmb_id = str(uuid.uuid4())
        observation_ids.append(tmb_id)
        tmb_obs = {
            "resourceType": "Observation",
            "id": tmb_id,
            "status": "final",
            "category": [{
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                    "code": "laboratory",
                    "display": "Laboratory",
                }],
            }],
            "code": {
                "coding": [{
                    "system": "http://loinc.org",
                    "code": FHIR_LOINC_CODES["tumor_mutation_burden"],
                    "display": "Tumor mutation burden",
                }],
            },
            "subject": {"reference": f"urn:uuid:{patient_resource_id}"},
            "effectiveDateTime": timestamp,
            "valueQuantity": {
                "value": float(tmb) if not isinstance(tmb, str) else 0,
                "unit": "mutations/megabase",
                "system": "http://unitsofmeasure.org",
                "code": "{mutations}/Mb",
            },
        }
        entries.append({
            "fullUrl": f"urn:uuid:{tmb_id}",
            "resource": tmb_obs,
        })

    # --- MSI Observation ---
    msi = biomarkers.get("msi") or biomarkers.get("microsatellite_instability")
    if msi is not None:
        msi_id = str(uuid.uuid4())
        observation_ids.append(msi_id)
        msi_obs = {
            "resourceType": "Observation",
            "id": msi_id,
            "status": "final",
            "category": [{
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/observation-category",
                    "code": "laboratory",
                    "display": "Laboratory",
                }],
            }],
            "code": {
                "coding": [{
                    "system": "http://loinc.org",
                    "code": FHIR_LOINC_CODES["microsatellite_instability"],
                    "display": "Microsatellite instability [Interpretation]",
                }],
            },
            "subject": {"reference": f"urn:uuid:{patient_resource_id}"},
            "effectiveDateTime": timestamp,
            "valueCodeableConcept": {
                "text": str(msi),
            },
        }
        entries.append({
            "fullUrl": f"urn:uuid:{msi_id}",
            "resource": msi_obs,
        })

    # --- Specimen Resource ---
    specimen_id = str(uuid.uuid4())
    specimen_resource = {
        "resourceType": "Specimen",
        "id": specimen_id,
        "subject": {"reference": f"urn:uuid:{patient_resource_id}"},
        "type": {
            "coding": [{
                "system": "http://snomed.info/sct",
                "code": "119376003",
                "display": "Tissue specimen",
            }],
            "text": data.get("sample_type", "Tumor tissue"),
        },
        "collection": {
            "collectedDateTime": timestamp,
        },
    }
    if data.get("sample_id"):
        specimen_resource["identifier"] = [{
            "system": "urn:hcls-ai-factory:specimen",
            "value": data["sample_id"],
        }]
    entries.append({
        "fullUrl": f"urn:uuid:{specimen_id}",
        "resource": specimen_resource,
    })

    # --- Condition Resource (cancer diagnosis) ---
    cancer_type_raw = data.get("cancer_type", "").lower().strip()
    cancer_coding_cond = FHIR_SNOMED_CANCER_CODES.get(cancer_type_raw)
    if cancer_coding_cond:
        condition_id = str(uuid.uuid4())
        condition_resource = {
            "resourceType": "Condition",
            "id": condition_id,
            "clinicalStatus": {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/condition-clinical",
                    "code": "active",
                    "display": "Active",
                }],
            },
            "verificationStatus": {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/condition-ver-status",
                    "code": "confirmed",
                    "display": "Confirmed",
                }],
            },
            "code": {
                "coding": [{
                    "system": "http://snomed.info/sct",
                    "code": cancer_coding_cond[0],
                    "display": cancer_coding_cond[1],
                }],
                "text": data.get("cancer_type", ""),
            },
            "subject": {"reference": f"urn:uuid:{patient_resource_id}"},
            "recordedDate": timestamp,
        }
        stage = data.get("stage")
        if stage:
            condition_resource["stage"] = [{
                "summary": {"text": f"Stage {stage}"},
            }]
        entries.append({
            "fullUrl": f"urn:uuid:{condition_id}",
            "resource": condition_resource,
        })

    # --- MedicationRequest Resources (therapy recommendations) ---
    therapies_data = (data.get("therapies") or data.get("therapy_ranking")
                      or data.get("recommendations") or [])
    for tx in therapies_data[:10]:  # cap to top 10 recommendations
        med_id = str(uuid.uuid4())
        drug_name = tx.get("name", tx.get("drug_name", tx.get("drug", "")))
        if not drug_name:
            continue
        med_request = {
            "resourceType": "MedicationRequest",
            "id": med_id,
            "status": "draft",
            "intent": "proposal",
            "medicationCodeableConcept": {
                "text": drug_name,
            },
            "subject": {"reference": f"urn:uuid:{patient_resource_id}"},
            "authoredOn": timestamp,
            "note": [],
        }
        ev_level = tx.get("evidence_level", tx.get("level", ""))
        if ev_level:
            med_request["note"].append({"text": f"Evidence level: {ev_level}"})
        guideline = tx.get("guideline_recommendation", tx.get("rationale", ""))
        if guideline:
            med_request["note"].append({"text": guideline})
        entries.append({
            "fullUrl": f"urn:uuid:{med_id}",
            "resource": med_request,
        })

    # --- DiagnosticReport Resource ---
    report_id = str(uuid.uuid4())
    cancer_type = data.get("cancer_type", "").lower().strip()
    cancer_coding = FHIR_SNOMED_CANCER_CODES.get(cancer_type)

    diagnostic_report = {
        "resourceType": "DiagnosticReport",
        "id": report_id,
        "status": "final",
        "category": [{
            "coding": [{
                "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                "code": "GE",
                "display": "Genetics",
            }],
        }],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": FHIR_LOINC_CODES["genomic_report"],
                "display": "Master HL7 genetic variant reporting panel",
            }],
        },
        "subject": {"reference": f"urn:uuid:{patient_resource_id}"},
        "specimen": [{"reference": f"urn:uuid:{specimen_id}"}],
        "effectiveDateTime": timestamp,
        "issued": timestamp,
        "result": [
            {"reference": f"urn:uuid:{obs_id}"} for obs_id in observation_ids
        ],
        "conclusion": data.get("summary") or data.get("clinical_summary", ""),
    }

    # Add cancer type coding if we have a SNOMED match
    if cancer_coding:
        diagnostic_report["conclusionCode"] = [{
            "coding": [{
                "system": "http://snomed.info/sct",
                "code": cancer_coding[0],
                "display": cancer_coding[1],
            }],
        }]

    entries.append({
        "fullUrl": f"urn:uuid:{report_id}",
        "resource": diagnostic_report,
    })

    # --- Build Bundle ---
    bundle = {
        "resourceType": "Bundle",
        "id": bundle_id,
        "type": "collection",
        "timestamp": timestamp,
        "meta": {
            "profile": [
                "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report",
            ],
        },
        "entry": entries,
    }

    logger.info(
        "Exported FHIR R4 Bundle with %d entries (patient=%s)",
        len(entries), patient_id,
    )
    return bundle


# ═══════════════════════════════════════════════════════════════════════════
#  Convenience wrappers for reports router
# ═══════════════════════════════════════════════════════════════════════════


def case_to_markdown(case_data: Any) -> str:
    """Convert a case dict/snapshot to Markdown report."""
    return export_markdown(case_data)


def markdown_to_pdf(markdown_text: str) -> bytes:
    """Convert a Markdown string to PDF bytes via ReportLab."""
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.units import inch
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
    from reportlab.lib.styles import getSampleStyleSheet
    import io as _io

    buf = _io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=letter,
                            leftMargin=0.75 * inch, rightMargin=0.75 * inch,
                            topMargin=0.75 * inch, bottomMargin=0.75 * inch)
    styles = getSampleStyleSheet()
    story = []
    for line in markdown_text.split("\n"):
        stripped = line.strip()
        if stripped.startswith("# "):
            story.append(Paragraph(stripped[2:], styles["Heading1"]))
        elif stripped.startswith("## "):
            story.append(Paragraph(stripped[3:], styles["Heading2"]))
        elif stripped.startswith("### "):
            story.append(Paragraph(stripped[4:], styles["Heading3"]))
        elif stripped.startswith("**") and stripped.endswith("**"):
            story.append(Paragraph(f"<b>{stripped[2:-2]}</b>", styles["Normal"]))
        elif stripped == "---":
            story.append(Spacer(1, 12))
        elif stripped:
            story.append(Paragraph(stripped, styles["Normal"]))
        else:
            story.append(Spacer(1, 6))
    doc.build(story)
    return buf.getvalue()


def case_to_fhir_bundle(case_data: Any) -> dict:
    """Convert a case dict/snapshot to a FHIR R4 Bundle."""
    data = _normalise_input(case_data)
    patient_id = data.get("patient_id", "unknown")
    return export_fhir_r4(case_data, patient_id=patient_id)
