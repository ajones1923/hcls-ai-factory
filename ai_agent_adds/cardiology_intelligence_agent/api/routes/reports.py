"""Cardiology report generation and export routes.

Provides endpoints for generating structured cardiology reports in
multiple formats: Markdown, JSON, PDF, and FHIR R4 DiagnosticReport.
Supports clinical summaries, risk assessment reports, GDMT optimization
reports, and workflow-specific outputs.

Author: Adam Jones
Date: March 2026
"""

import base64
import io
import json
import uuid
from datetime import datetime, timezone
from typing import Optional

from fastapi import APIRouter, HTTPException, Request
from loguru import logger
from pydantic import BaseModel, Field

try:
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.colors import HexColor
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        HRFlowable,
    )
    HAS_REPORTLAB = True
except ImportError:
    HAS_REPORTLAB = False

router = APIRouter(prefix="/v1/reports", tags=["reports"])


# =====================================================================
# Schemas
# =====================================================================

class ReportRequest(BaseModel):
    """Request to generate a cardiology report."""
    report_type: str = Field(
        ...,
        description="Type: clinical_summary | risk_assessment | gdmt_report | workflow_report | imaging_report",
    )
    format: str = Field("markdown", pattern="^(markdown|json|pdf|fhir)$")
    patient_id: Optional[str] = None
    title: Optional[str] = None
    data: dict = Field(default={}, description="Report payload (risk results, workflow output, etc.)")
    include_evidence: bool = True
    include_recommendations: bool = True


class ReportResponse(BaseModel):
    report_id: str
    report_type: str
    format: str
    generated_at: str
    title: str
    content: str  # Markdown/JSON string or base64 for PDF
    metadata: dict = {}


# =====================================================================
# Report Templates
# =====================================================================

def _generate_markdown_header(title: str, patient_id: Optional[str] = None) -> str:
    """Standard markdown report header."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        f"# {title}",
        "",
        f"**Generated:** {now}",
        "**Agent:** Cardiology Intelligence Agent v1.0.0",
    ]
    if patient_id:
        lines.append(f"**Patient ID:** {patient_id}")
    lines.extend(["", "---", ""])
    return "\n".join(lines)


def _risk_assessment_markdown(data: dict) -> str:
    """Format risk calculator results as markdown."""
    sections = []
    results = data.get("results", [])
    if not results and "calculator" in data:
        results = [data]

    for r in results:
        calc = r.get("calculator", "Unknown")
        score = r.get("score", "N/A")
        category = r.get("risk_category", "N/A")
        interp = r.get("interpretation", "")
        recs = r.get("recommendations", [])

        section = [
            f"## {calc}",
            "",
            "| Metric | Value |",
            "|--------|-------|",
            f"| Score | **{score}** |",
            f"| Risk Category | **{category}** |",
            f"| Interpretation | {interp} |",
            "",
        ]
        if recs:
            section.append("### Recommendations")
            for rec in recs:
                section.append(f"- {rec}")
            section.append("")

        details = r.get("details", {})
        if details:
            section.append("### Details")
            for k, v in details.items():
                section.append(f"- **{k.replace('_', ' ').title()}:** {v}")
            section.append("")

        sections.append("\n".join(section))

    return "\n".join(sections)


def _gdmt_report_markdown(data: dict) -> str:
    """Format GDMT optimization output as markdown."""
    lines = [
        f"## Heart Failure Phenotype: {data.get('hf_phenotype', 'Unknown')}",
        f"- **LVEF:** {data.get('lvef', 'N/A')}%",
        f"- **NYHA Class:** {data.get('nyha_class', 'N/A')}",
        "",
        "## Four Pillars Status",
        "",
        "| Pillar | Status |",
        "|--------|--------|",
    ]
    for pillar, status in data.get("four_pillars_status", {}).items():
        lines.append(f"| {pillar} | **{status}** |")

    lines.extend(["", "## Optimization Recommendations", ""])
    for rec in data.get("recommendations", []):
        if isinstance(rec, dict):
            lines.append(f"### {rec.get('medication_class', 'Unknown')}")
            lines.append(f"- **Current:** {rec.get('current_status', 'N/A')}")
            lines.append(f"- **Action:** {rec.get('recommendation', 'N/A')}")
            if rec.get("target_dose"):
                lines.append(f"- **Target:** {rec['target_dose']}")
            if rec.get("evidence_level"):
                lines.append(f"- **Evidence:** {rec['evidence_level']}")
            if rec.get("guideline_ref"):
                lines.append(f"- **Reference:** {rec['guideline_ref']}")
            if rec.get("caution"):
                lines.append(f"- **Caution:** {rec['caution']}")
            lines.append("")

    if data.get("summary"):
        lines.extend(["## Summary", "", data["summary"], ""])

    return "\n".join(lines)


def _generate_fhir_diagnostic_report(data: dict, title: str, patient_id: Optional[str]) -> dict:
    """Generate a FHIR R4 DiagnosticReport resource."""
    now = datetime.now(timezone.utc).isoformat()
    report = {
        "resourceType": "DiagnosticReport",
        "id": str(uuid.uuid4()),
        "status": "final",
        "category": [{
            "coding": [{
                "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                "code": "CUS",
                "display": "Cardiac Ultrasound",
            }]
        }],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "18842-5",
                "display": "Cardiology Assessment",
            }],
            "text": title,
        },
        "effectiveDateTime": now,
        "issued": now,
        "conclusion": data.get("summary", data.get("interpretation", "")),
        "presentedForm": [{
            "contentType": "text/markdown",
            "data": title,
        }],
    }
    if patient_id:
        report["subject"] = {"reference": f"Patient/{patient_id}"}

    # Add contained observations for risk scores
    observations = []
    results = data.get("results", [])
    if "calculator" in data:
        results = [data]
    for r in results:
        obs = {
            "resourceType": "Observation",
            "id": str(uuid.uuid4()),
            "status": "final",
            "code": {"text": r.get("calculator", "Risk Score")},
            "valueQuantity": {
                "value": r.get("score", 0),
                "unit": "score",
            },
            "interpretation": [{
                "text": r.get("risk_category", ""),
            }],
        }
        observations.append(obs)

    if observations:
        report["contained"] = observations
        report["result"] = [{"reference": f"#{obs['id']}"} for obs in observations]

    return report


# =====================================================================
# PDF Generation (ReportLab)
# =====================================================================

# NVIDIA dark-theme palette
_NAVY = HexColor("#1B2333") if HAS_REPORTLAB else None
_TEAL = HexColor("#1AAFCC") if HAS_REPORTLAB else None
_GREEN = HexColor("#76B900") if HAS_REPORTLAB else None
_WHITE = HexColor("#FFFFFF") if HAS_REPORTLAB else None
_LIGHT_GRAY = HexColor("#D0D0D0") if HAS_REPORTLAB else None
_DARK_ROW = HexColor("#232F3E") if HAS_REPORTLAB else None


def _generate_pdf(
    title: str,
    report_type: str,
    patient_id: Optional[str],
    data: dict,
) -> str:
    """Build a styled PDF using ReportLab and return its base64-encoded bytes."""
    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf,
        pagesize=letter,
        topMargin=0.6 * inch,
        bottomMargin=0.5 * inch,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
    )

    # -- Custom styles --
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        "CardioTitle",
        parent=styles["Title"],
        fontSize=20,
        textColor=_TEAL,
        spaceAfter=6,
    ))
    styles.add(ParagraphStyle(
        "CardioHeading",
        parent=styles["Heading2"],
        fontSize=14,
        textColor=_GREEN,
        spaceBefore=12,
        spaceAfter=4,
    ))
    styles.add(ParagraphStyle(
        "CardioBody",
        parent=styles["Normal"],
        fontSize=10,
        textColor=_LIGHT_GRAY,
        leading=14,
    ))
    styles.add(ParagraphStyle(
        "CardioMeta",
        parent=styles["Normal"],
        fontSize=8,
        textColor=_LIGHT_GRAY,
        spaceAfter=2,
    ))

    story: list = []

    # -- Title & metadata --
    story.append(Paragraph(title, styles["CardioTitle"]))
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    story.append(Paragraph(f"Generated: {now} | Cardiology Intelligence Agent v1.0.0", styles["CardioMeta"]))
    if patient_id:
        story.append(Paragraph(f"Patient ID: {patient_id}", styles["CardioMeta"]))
    story.append(HRFlowable(width="100%", thickness=1, color=_TEAL, spaceAfter=10))

    # -- Report body by type --
    if report_type == "risk_assessment":
        results = data.get("results", [])
        if not results and "calculator" in data:
            results = [data]
        for r in results:
            story.append(Paragraph(r.get("calculator", "Risk Score"), styles["CardioHeading"]))
            table_data = [
                ["Metric", "Value"],
                ["Score", str(r.get("score", "N/A"))],
                ["Risk Category", str(r.get("risk_category", "N/A"))],
                ["Interpretation", str(r.get("interpretation", ""))],
            ]
            t = Table(table_data, colWidths=[2.5 * inch, 4 * inch])
            t.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), _TEAL),
                ("TEXTCOLOR", (0, 0), (-1, 0), _WHITE),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [_NAVY, _DARK_ROW]),
                ("TEXTCOLOR", (0, 1), (-1, -1), _LIGHT_GRAY),
                ("GRID", (0, 0), (-1, -1), 0.5, _TEAL),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]))
            story.append(t)
            story.append(Spacer(1, 8))
            for rec in r.get("recommendations", []):
                story.append(Paragraph(f"&bull; {rec}", styles["CardioBody"]))

    elif report_type == "gdmt_report":
        story.append(Paragraph(
            f"HF Phenotype: {data.get('hf_phenotype', 'Unknown')} | "
            f"LVEF: {data.get('lvef', 'N/A')}% | "
            f"NYHA: {data.get('nyha_class', 'N/A')}",
            styles["CardioBody"],
        ))
        story.append(Spacer(1, 6))
        pillars = data.get("four_pillars_status", {})
        if pillars:
            story.append(Paragraph("Four Pillars Status", styles["CardioHeading"]))
            table_data = [["Pillar", "Status"]]
            for pillar, status in pillars.items():
                table_data.append([pillar, status])
            t = Table(table_data, colWidths=[3.25 * inch, 3.25 * inch])
            t.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), _TEAL),
                ("TEXTCOLOR", (0, 0), (-1, 0), _WHITE),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [_NAVY, _DARK_ROW]),
                ("TEXTCOLOR", (0, 1), (-1, -1), _LIGHT_GRAY),
                ("GRID", (0, 0), (-1, -1), 0.5, _TEAL),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 4),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
            ]))
            story.append(t)
        for rec in data.get("recommendations", []):
            if isinstance(rec, dict):
                story.append(Paragraph(rec.get("medication_class", ""), styles["CardioHeading"]))
                story.append(Paragraph(
                    f"Current: {rec.get('current_status', 'N/A')} | "
                    f"Action: {rec.get('recommendation', 'N/A')}",
                    styles["CardioBody"],
                ))

    else:
        # Generic data rendering
        for key, value in data.items():
            story.append(Paragraph(key.replace("_", " ").title(), styles["CardioHeading"]))
            if isinstance(value, list):
                for item in value:
                    story.append(Paragraph(f"&bull; {item}", styles["CardioBody"]))
            elif isinstance(value, dict):
                for k, v in value.items():
                    story.append(Paragraph(f"<b>{k}:</b> {v}", styles["CardioBody"]))
            else:
                story.append(Paragraph(str(value), styles["CardioBody"]))

    # -- Footer --
    story.append(Spacer(1, 20))
    story.append(HRFlowable(width="100%", thickness=0.5, color=_GREEN, spaceAfter=4))
    story.append(Paragraph(
        "Cardiology Intelligence Agent | HCLS AI Factory | NVIDIA DGX Spark",
        styles["CardioMeta"],
    ))

    def _on_page(canvas, doc_obj):
        """Draw NVIDIA dark background on every page."""
        canvas.saveState()
        canvas.setFillColor(_NAVY)
        canvas.rect(0, 0, letter[0], letter[1], fill=1, stroke=0)
        canvas.restoreState()

    doc.build(story, onFirstPage=_on_page, onLaterPages=_on_page)
    pdf_bytes = buf.getvalue()
    buf.close()
    return base64.b64encode(pdf_bytes).decode("ascii")


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/generate", response_model=ReportResponse)
async def generate_report(request: ReportRequest, req: Request):
    """Generate a formatted cardiology report."""
    report_id = str(uuid.uuid4())[:12]
    now = datetime.now(timezone.utc).isoformat()
    title = request.title or f"Cardiology {request.report_type.replace('_', ' ').title()}"

    try:
        if request.format == "fhir":
            fhir_resource = _generate_fhir_diagnostic_report(
                request.data, title, request.patient_id,
            )
            content = json.dumps(fhir_resource, indent=2)

        elif request.format == "json":
            content = json.dumps({
                "report_id": report_id,
                "title": title,
                "type": request.report_type,
                "generated": now,
                "patient_id": request.patient_id,
                "data": request.data,
            }, indent=2)

        elif request.format == "pdf":
            if not HAS_REPORTLAB:
                raise HTTPException(
                    status_code=501,
                    detail="PDF generation requires reportlab. Install with: pip install reportlab",
                )
            content = _generate_pdf(title, request.report_type, request.patient_id, request.data)

        else:  # markdown
            header = _generate_markdown_header(title, request.patient_id)
            if request.report_type == "risk_assessment":
                body = _risk_assessment_markdown(request.data)
            elif request.report_type == "gdmt_report":
                body = _gdmt_report_markdown(request.data)
            else:
                # Generic markdown body
                body_lines = []
                for key, value in request.data.items():
                    body_lines.append(f"## {key.replace('_', ' ').title()}")
                    if isinstance(value, list):
                        for item in value:
                            body_lines.append(f"- {item}")
                    elif isinstance(value, dict):
                        for k, v in value.items():
                            body_lines.append(f"- **{k}:** {v}")
                    else:
                        body_lines.append(str(value))
                    body_lines.append("")
                body = "\n".join(body_lines)
            content = header + body

        with req.app.state.metrics_lock:
            req.app.state.metrics["report_requests_total"] += 1

        return ReportResponse(
            report_id=report_id,
            report_type=request.report_type,
            format=request.format,
            generated_at=now,
            title=title,
            content=content,
            metadata={
                "agent": "cardiology-intelligence-agent",
                "version": "1.0.0",
                "data_keys": list(request.data.keys()),
            },
        )

    except Exception as exc:
        logger.error(f"Report generation failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@router.get("/formats")
async def list_formats():
    """List supported report export formats."""
    return {
        "formats": [
            {"id": "markdown", "name": "Markdown", "extension": ".md", "mime": "text/markdown", "description": "Human-readable clinical report"},
            {"id": "json", "name": "JSON", "extension": ".json", "mime": "application/json", "description": "Structured data export"},
            {"id": "pdf", "name": "PDF", "extension": ".pdf", "mime": "application/pdf", "description": "Printable clinical report"},
            {"id": "fhir", "name": "FHIR R4", "extension": ".json", "mime": "application/fhir+json", "description": "HL7 FHIR R4 DiagnosticReport resource"},
        ],
        "report_types": [
            "clinical_summary",
            "risk_assessment",
            "gdmt_report",
            "workflow_report",
            "imaging_report",
        ],
    }
