"""Clinical trial report generation and export routes.

Provides endpoints for generating structured clinical trial reports in
multiple formats: Markdown, JSON, PDF, and FHIR R4 DiagnosticReport.
Supports protocol summaries, matching reports, safety reports, regulatory
summaries, and competitive landscape reports.

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
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import (
        Paragraph,
        SimpleDocTemplate,
        Spacer,
        Table,
        TableStyle,
    )
    HAS_REPORTLAB = True
except ImportError:
    HAS_REPORTLAB = False

# -- NVIDIA Dark Theme colours --
_NAVY = colors.HexColor("#1B2333") if HAS_REPORTLAB else None
_TEAL = colors.HexColor("#1AAFCC") if HAS_REPORTLAB else None
_GREEN = colors.HexColor("#76B900") if HAS_REPORTLAB else None

router = APIRouter(prefix="/v1/reports", tags=["reports"])


# =====================================================================
# Schemas
# =====================================================================

class ReportRequest(BaseModel):
    """Request to generate a clinical trial report."""
    report_type: str = Field(
        ...,
        description=(
            "Type: protocol_summary | matching_report | safety_report | "
            "regulatory_summary | competitive_landscape | eligibility_analysis | "
            "diversity_report | dct_assessment"
        ),
    )
    format: str = Field("markdown", pattern="^(markdown|json|pdf|fhir)$")
    patient_id: Optional[str] = None
    trial_id: Optional[str] = None
    title: Optional[str] = None
    data: dict = Field(default={}, description="Report payload (matching results, safety data, etc.)")
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

def _generate_markdown_header(title: str, trial_id: Optional[str] = None, patient_id: Optional[str] = None) -> str:
    """Standard markdown report header."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        f"# {title}",
        "",
        f"**Generated:** {now}",
        "**Agent:** Clinical Trial Intelligence Agent v1.0.0",
    ]
    if trial_id:
        lines.append(f"**Trial ID:** {trial_id}")
    if patient_id:
        lines.append(f"**Patient ID:** {patient_id}")
    lines.extend(["", "---", ""])
    return "\n".join(lines)


def _protocol_summary_markdown(data: dict) -> str:
    """Format protocol optimization results as markdown."""
    lines = []
    complexity = data.get("complexity_score", "N/A")
    percentile = data.get("percentile_rank", "N/A")

    lines.extend([
        "## Protocol Complexity Assessment",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Complexity Score | **{complexity}** |",
        f"| Percentile Rank | **{percentile}** |",
        "",
    ])

    recs = data.get("optimization_recommendations", [])
    if recs:
        lines.append("## Optimization Recommendations")
        lines.append("")
        for rec in recs:
            lines.append(f"- {rec}")
        lines.append("")

    risks = data.get("risk_factors", [])
    if risks:
        lines.append("## Risk Factors")
        lines.append("")
        for risk in risks:
            if isinstance(risk, dict):
                lines.append(f"- **{risk.get('factor', 'Unknown')}:** {risk.get('description', '')}")
            else:
                lines.append(f"- {risk}")
        lines.append("")

    return "\n".join(lines)


def _matching_report_markdown(data: dict) -> str:
    """Format patient-trial matching results as markdown."""
    lines = []
    matches = data.get("matches", [])

    lines.extend([
        "## Patient-Trial Matching Results",
        "",
        f"**Total Screened:** {data.get('total_screened', 'N/A')}",
        f"**Matches Found:** {len(matches)}",
        "",
    ])

    if matches:
        lines.extend([
            "| Trial ID | Title | Phase | Score | Confidence |",
            "|----------|-------|-------|-------|------------|",
        ])
        for m in matches:
            if isinstance(m, dict):
                lines.append(
                    f"| {m.get('trial_id', 'N/A')} | "
                    f"{m.get('trial_title', 'N/A')[:50]} | "
                    f"{m.get('phase', 'N/A')} | "
                    f"**{m.get('overall_score', 'N/A')}** | "
                    f"{m.get('confidence', 'N/A')} |"
                )
        lines.append("")

    return "\n".join(lines)


def _safety_report_markdown(data: dict) -> str:
    """Format safety signal detection results as markdown."""
    lines = [
        "## Safety Signal Analysis",
        "",
    ]

    severity = data.get("severity_distribution", {})
    if severity:
        lines.extend([
            "### Severity Distribution",
            "",
            "| Severity | Count |",
            "|----------|-------|",
        ])
        for sev, count in severity.items():
            lines.append(f"| {sev.title()} | {count} |")
        lines.append("")

    signals = data.get("signals_detected", [])
    if signals:
        lines.append("### Detected Signals")
        lines.append("")
        for sig in signals:
            if isinstance(sig, dict):
                lines.append(f"- **{sig.get('event_type', 'Unknown')}** (PRR: {sig.get('prr', 'N/A')})")
            else:
                lines.append(f"- {sig}")
        lines.append("")

    recs = data.get("recommendations", [])
    if recs:
        lines.append("### Recommendations")
        lines.append("")
        for rec in recs:
            lines.append(f"- {rec}")
        lines.append("")

    return "\n".join(lines)


def _generate_fhir_research_study(data: dict, title: str, trial_id: Optional[str]) -> dict:
    """Generate a FHIR R4 ResearchStudy-like resource."""
    now = datetime.now(timezone.utc).isoformat()
    resource = {
        "resourceType": "ResearchStudy",
        "id": str(uuid.uuid4()),
        "status": "active",
        "title": title,
        "description": data.get("summary", ""),
        "period": {"start": now},
        "meta": {
            "lastUpdated": now,
            "source": "clinical-trial-intelligence-agent",
        },
    }
    if trial_id:
        resource["identifier"] = [{"system": "https://clinicaltrials.gov", "value": trial_id}]
    return resource


# =====================================================================
# PDF Generation (ReportLab)
# =====================================================================

def _generate_pdf(
    title: str,
    report_type: str,
    data: dict,
    trial_id: Optional[str] = None,
    patient_id: Optional[str] = None,
) -> str:
    """Generate a branded PDF report and return it as a base64-encoded string."""
    buf = io.BytesIO()
    doc = SimpleDocTemplate(
        buf,
        pagesize=letter,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
    )

    styles = getSampleStyleSheet()

    # Custom styles using NVIDIA dark theme colours
    title_style = ParagraphStyle(
        "NVTitle",
        parent=styles["Title"],
        textColor=_NAVY,
        fontSize=20,
        spaceAfter=12,
    )
    heading_style = ParagraphStyle(
        "NVHeading",
        parent=styles["Heading2"],
        textColor=_TEAL,
        fontSize=14,
        spaceAfter=8,
    )
    body_style = ParagraphStyle(
        "NVBody",
        parent=styles["BodyText"],
        fontSize=10,
        leading=14,
        spaceAfter=6,
    )

    elements: list = []

    # Title
    elements.append(Paragraph(title, title_style))
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    meta_parts = [f"Generated: {now}", "Agent: Clinical Trial Intelligence Agent v1.0.0"]
    if trial_id:
        meta_parts.append(f"Trial ID: {trial_id}")
    if patient_id:
        meta_parts.append(f"Patient ID: {patient_id}")
    elements.append(Paragraph(" | ".join(meta_parts), body_style))
    elements.append(Spacer(1, 0.25 * inch))

    # Report-type specific sections
    if report_type == "protocol_summary":
        _pdf_protocol_summary(elements, data, heading_style, body_style)
    elif report_type == "matching_report":
        _pdf_matching_report(elements, data, heading_style, body_style)
    elif report_type == "safety_report":
        _pdf_safety_report(elements, data, heading_style, body_style)
    else:
        _pdf_generic(elements, data, heading_style, body_style)

    doc.build(elements)
    pdf_bytes = buf.getvalue()
    buf.close()
    return base64.b64encode(pdf_bytes).decode("utf-8")


def _nv_table_style() -> TableStyle:
    """Return a TableStyle using the NVIDIA dark theme."""
    return TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), _NAVY),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, 0), 10),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
        ("TOPPADDING", (0, 0), (-1, 0), 8),
        ("BACKGROUND", (0, 1), (-1, -1), colors.HexColor("#F4F6F8")),
        ("GRID", (0, 0), (-1, -1), 0.5, _TEAL),
        ("FONTSIZE", (0, 1), (-1, -1), 9),
        ("TOPPADDING", (0, 1), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 1), (-1, -1), 5),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.HexColor("#F4F6F8"), colors.white]),
    ])


def _pdf_protocol_summary(elements, data, heading_style, body_style):
    complexity = data.get("complexity_score", "N/A")
    percentile = data.get("percentile_rank", "N/A")

    elements.append(Paragraph("Protocol Complexity Assessment", heading_style))
    table_data = [
        ["Metric", "Value"],
        ["Complexity Score", str(complexity)],
        ["Percentile Rank", str(percentile)],
    ]
    t = Table(table_data, colWidths=[3 * inch, 3 * inch])
    t.setStyle(_nv_table_style())
    elements.append(t)
    elements.append(Spacer(1, 0.2 * inch))

    recs = data.get("optimization_recommendations", [])
    if recs:
        elements.append(Paragraph("Optimization Recommendations", heading_style))
        for rec in recs:
            elements.append(Paragraph(f"\u2022 {rec}", body_style))
        elements.append(Spacer(1, 0.15 * inch))

    risks = data.get("risk_factors", [])
    if risks:
        elements.append(Paragraph("Risk Factors", heading_style))
        for risk in risks:
            if isinstance(risk, dict):
                elements.append(Paragraph(
                    f"\u2022 <b>{risk.get('factor', 'Unknown')}:</b> {risk.get('description', '')}",
                    body_style,
                ))
            else:
                elements.append(Paragraph(f"\u2022 {risk}", body_style))


def _pdf_matching_report(elements, data, heading_style, body_style):
    matches = data.get("matches", [])
    elements.append(Paragraph("Patient-Trial Matching Results", heading_style))
    elements.append(Paragraph(
        f"Total Screened: {data.get('total_screened', 'N/A')} | "
        f"Matches Found: {len(matches)}",
        body_style,
    ))
    elements.append(Spacer(1, 0.1 * inch))

    if matches:
        table_data = [["Trial ID", "Title", "Phase", "Score", "Confidence"]]
        for m in matches:
            if isinstance(m, dict):
                table_data.append([
                    m.get("trial_id", "N/A"),
                    m.get("trial_title", "N/A")[:45],
                    m.get("phase", "N/A"),
                    str(m.get("overall_score", "N/A")),
                    str(m.get("confidence", "N/A")),
                ])
        if len(table_data) > 1:
            t = Table(table_data, colWidths=[1.1 * inch, 2.4 * inch, 0.9 * inch, 0.8 * inch, 0.9 * inch])
            t.setStyle(_nv_table_style())
            elements.append(t)


def _pdf_safety_report(elements, data, heading_style, body_style):
    elements.append(Paragraph("Safety Signal Analysis", heading_style))

    severity = data.get("severity_distribution", {})
    if severity:
        elements.append(Paragraph("Severity Distribution", heading_style))
        table_data = [["Severity", "Count"]]
        for sev, count in severity.items():
            table_data.append([sev.title(), str(count)])
        t = Table(table_data, colWidths=[3 * inch, 3 * inch])
        t.setStyle(_nv_table_style())
        elements.append(t)
        elements.append(Spacer(1, 0.15 * inch))

    signals = data.get("signals_detected", [])
    if signals:
        elements.append(Paragraph("Detected Signals", heading_style))
        for sig in signals:
            if isinstance(sig, dict):
                elements.append(Paragraph(
                    f"\u2022 <b>{sig.get('event_type', 'Unknown')}</b> (PRR: {sig.get('prr', 'N/A')})",
                    body_style,
                ))
            else:
                elements.append(Paragraph(f"\u2022 {sig}", body_style))

    recs = data.get("recommendations", [])
    if recs:
        elements.append(Spacer(1, 0.1 * inch))
        elements.append(Paragraph("Recommendations", heading_style))
        for rec in recs:
            elements.append(Paragraph(f"\u2022 {rec}", body_style))


def _pdf_generic(elements, data, heading_style, body_style):
    for key, value in data.items():
        elements.append(Paragraph(key.replace("_", " ").title(), heading_style))
        if isinstance(value, list):
            for item in value:
                elements.append(Paragraph(f"\u2022 {item}", body_style))
        elif isinstance(value, dict):
            for k, v in value.items():
                elements.append(Paragraph(f"<b>{k}:</b> {v}", body_style))
        else:
            elements.append(Paragraph(str(value), body_style))
        elements.append(Spacer(1, 0.1 * inch))


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/generate", response_model=ReportResponse)
async def generate_report(request: ReportRequest, req: Request):
    """Generate a formatted clinical trial report."""
    report_id = str(uuid.uuid4())[:12]
    now = datetime.now(timezone.utc).isoformat()
    title = request.title or f"Clinical Trial {request.report_type.replace('_', ' ').title()}"

    try:
        if request.format == "fhir":
            fhir_resource = _generate_fhir_research_study(
                request.data, title, request.trial_id,
            )
            content = json.dumps(fhir_resource, indent=2)

        elif request.format == "json":
            content = json.dumps({
                "report_id": report_id,
                "title": title,
                "type": request.report_type,
                "generated": now,
                "trial_id": request.trial_id,
                "patient_id": request.patient_id,
                "data": request.data,
            }, indent=2)

        elif request.format == "pdf":
            if not HAS_REPORTLAB:
                raise HTTPException(
                    status_code=501,
                    detail="PDF generation requires the reportlab package. Install with: pip install reportlab",
                )
            content = _generate_pdf(title, request.report_type, request.data,
                                    request.trial_id, request.patient_id)

        else:  # markdown
            header = _generate_markdown_header(title, request.trial_id, request.patient_id)
            if request.report_type == "protocol_summary":
                body = _protocol_summary_markdown(request.data)
            elif request.report_type == "matching_report":
                body = _matching_report_markdown(request.data)
            elif request.report_type == "safety_report":
                body = _safety_report_markdown(request.data)
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

        metrics = getattr(req.app.state, "metrics", None)
        lock = getattr(req.app.state, "metrics_lock", None)
        if metrics and lock:
            with lock:
                metrics["report_requests_total"] = metrics.get("report_requests_total", 0) + 1

        return ReportResponse(
            report_id=report_id,
            report_type=request.report_type,
            format=request.format,
            generated_at=now,
            title=title,
            content=content,
            metadata={
                "agent": "clinical-trial-intelligence-agent",
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
            {"id": "markdown", "name": "Markdown", "extension": ".md", "mime": "text/markdown", "description": "Human-readable clinical trial report"},
            {"id": "json", "name": "JSON", "extension": ".json", "mime": "application/json", "description": "Structured data export"},
            {"id": "pdf", "name": "PDF", "extension": ".pdf", "mime": "application/pdf", "description": "Printable clinical trial report"},
            {"id": "fhir", "name": "FHIR R4", "extension": ".json", "mime": "application/fhir+json", "description": "HL7 FHIR R4 ResearchStudy resource"},
        ],
        "report_types": [
            "protocol_summary",
            "matching_report",
            "safety_report",
            "regulatory_summary",
            "competitive_landscape",
            "eligibility_analysis",
            "diversity_report",
            "dct_assessment",
        ],
    }
