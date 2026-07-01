"""Single-cell report generation and export routes.

Provides endpoints for generating structured single-cell analysis reports
in multiple formats: Markdown, JSON, PDF, and FHIR R4 DiagnosticReport.
Supports cell type annotation reports, TME profiles, drug response
summaries, spatial analysis reports, and workflow results.

Author: Adam Jones
Date: March 2026
"""

import json
import uuid
from datetime import datetime, timezone
from typing import Optional

from fastapi import APIRouter, HTTPException, Request
from loguru import logger
from pydantic import BaseModel, Field

router = APIRouter(prefix="/v1/reports", tags=["reports"])


# =====================================================================
# Schemas
# =====================================================================

class ReportRequest(BaseModel):
    """Request to generate a single-cell analysis report."""
    report_type: str = Field(
        ...,
        description=(
            "Type: cell_type_annotation | tme_profile | drug_response | "
            "subclonal_analysis | spatial_niche | trajectory | "
            "ligand_receptor | biomarker | cart_validation | "
            "treatment_monitoring | general"
        ),
    )
    format: str = Field("markdown", pattern="^(markdown|json|pdf|fhir)$")
    patient_id: Optional[str] = None
    encounter_id: Optional[str] = None
    title: Optional[str] = None
    data: dict = Field(default={}, description="Report payload (analysis results, etc.)")
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

def _generate_markdown_header(title: str, patient_id: Optional[str] = None, encounter_id: Optional[str] = None) -> str:
    """Standard markdown report header."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        f"# {title}",
        "",
        f"**Generated:** {now}",
        "**Agent:** Single-Cell Intelligence Agent v1.0.0",
    ]
    if patient_id:
        lines.append(f"**Patient ID:** {patient_id}")
    if encounter_id:
        lines.append(f"**Encounter ID:** {encounter_id}")
    lines.extend(["", "---", ""])
    return "\n".join(lines)


def _cell_type_annotation_markdown(data: dict) -> str:
    """Format cell type annotation results as markdown."""
    lines = [
        "## Cell Type Annotation Results",
        "",
    ]
    cell_types = data.get("cell_types", [])
    if cell_types:
        lines.extend([
            "| Cell Type | Compartment | Proportion | Confidence |",
            "|-----------|-------------|------------|------------|",
        ])
        for ct in cell_types:
            if isinstance(ct, dict):
                lines.append(
                    f"| {ct.get('cell_type', 'N/A')} | "
                    f"{ct.get('compartment', 'N/A')} | "
                    f"{ct.get('proportion', 'N/A')} | "
                    f"{ct.get('confidence', 'N/A')} |"
                )
        lines.append("")

    strategy = data.get("strategy_used", "N/A")
    lines.append(f"**Strategy:** {strategy}")
    lines.append("")

    recs = data.get("recommendations", [])
    if recs:
        lines.append("## Recommendations")
        lines.append("")
        for rec in recs:
            lines.append(f"- {rec}")
        lines.append("")

    return "\n".join(lines)


def _tme_profile_markdown(data: dict) -> str:
    """Format TME profiling results as markdown."""
    lines = [
        "## Tumor Microenvironment Profile",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
        f"| TME Classification | **{data.get('tme_classification', 'N/A')}** |",
        f"| Immune Score | **{data.get('immune_score', 'N/A')}** |",
        f"| Stromal Score | **{data.get('stromal_score', 'N/A')}** |",
        f"| Immune Phenotype | **{data.get('immune_phenotype', 'N/A')}** |",
        "",
    ]

    therapy = data.get("therapy_prediction", {})
    if therapy:
        lines.append("### Therapy Predictions")
        lines.append("")
        for k, v in therapy.items():
            lines.append(f"- **{k}:** {v}")
        lines.append("")

    recs = data.get("recommendations", [])
    if recs:
        lines.append("## Recommendations")
        lines.append("")
        for rec in recs:
            lines.append(f"- {rec}")
        lines.append("")

    return "\n".join(lines)


def _generate_fhir_diagnostic_report(data: dict, title: str, patient_id: Optional[str]) -> dict:
    """Generate a FHIR R4 DiagnosticReport resource."""
    now = datetime.now(timezone.utc).isoformat()
    resource = {
        "resourceType": "DiagnosticReport",
        "id": str(uuid.uuid4()),
        "status": "final",
        "category": [{"coding": [{"system": "http://loinc.org", "code": "LP7796-8", "display": "Genomics"}]}],
        "code": {"text": title},
        "effectiveDateTime": now,
        "issued": now,
        "conclusion": data.get("interpretation", data.get("summary", "")),
        "meta": {
            "lastUpdated": now,
            "source": "single-cell-intelligence-agent",
        },
    }
    if patient_id:
        resource["subject"] = {"reference": f"Patient/{patient_id}"}
    return resource


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/generate", response_model=ReportResponse)
async def generate_report(request: ReportRequest, req: Request):
    """Generate a formatted single-cell analysis report."""
    report_id = str(uuid.uuid4())[:12]
    now = datetime.now(timezone.utc).isoformat()
    title = request.title or f"Single-Cell {request.report_type.replace('_', ' ').title()}"

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
                "encounter_id": request.encounter_id,
                "data": request.data,
            }, indent=2)

        elif request.format == "pdf":
            # Generate real PDF using reportlab with NVIDIA dark theme
            try:
                from reportlab.lib.pagesizes import letter
                from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
                from reportlab.lib.units import inch
                from reportlab.lib.colors import HexColor
                from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
                import io

                buffer = io.BytesIO()
                doc_pdf = SimpleDocTemplate(buffer, pagesize=letter,
                                           topMargin=0.75*inch, bottomMargin=0.75*inch,
                                           leftMargin=0.75*inch, rightMargin=0.75*inch)
                styles = getSampleStyleSheet()
                navy = HexColor("#1B2333")
                teal = HexColor("#1AAFCC")
                HexColor("#76B900")

                title_style = ParagraphStyle("SCTitle", parent=styles["Title"],
                                             textColor=navy, fontSize=18, spaceAfter=12)
                heading_style = ParagraphStyle("SCHeading", parent=styles["Heading2"],
                                               textColor=teal, fontSize=13, spaceAfter=8)
                body_style = ParagraphStyle("SCBody", parent=styles["Normal"],
                                            fontSize=10, leading=14, spaceAfter=6)
                footer_style = ParagraphStyle("SCFooter", parent=styles["Normal"],
                                              fontSize=8, textColor=HexColor("#999999"))

                elements = []
                elements.append(Paragraph(title, title_style))
                elements.append(Paragraph("Single-Cell Intelligence Agent — HCLS AI Factory", heading_style))
                elements.append(Spacer(1, 12))

                # Metadata section
                if request.patient_id:
                    elements.append(Paragraph(f"<b>Patient ID:</b> {request.patient_id}", body_style))
                if request.encounter_id:
                    elements.append(Paragraph(f"<b>Encounter ID:</b> {request.encounter_id}", body_style))
                elements.append(Paragraph(f"<b>Report Type:</b> {request.report_type}", body_style))
                elements.append(Paragraph(f"<b>Generated:</b> {now}", body_style))
                elements.append(Spacer(1, 12))

                # Data sections
                for key, value in request.data.items():
                    elements.append(Paragraph(key.replace("_", " ").title(), heading_style))
                    if isinstance(value, list):
                        for item in value:
                            if isinstance(item, dict):
                                for k, v in item.items():
                                    elements.append(Paragraph(f"<b>{k}:</b> {v}", body_style))
                                elements.append(Spacer(1, 6))
                            else:
                                elements.append(Paragraph(f"&bull; {item}", body_style))
                    elif isinstance(value, dict):
                        for k, v in value.items():
                            elements.append(Paragraph(f"<b>{k}:</b> {v}", body_style))
                    else:
                        elements.append(Paragraph(str(value), body_style))
                    elements.append(Spacer(1, 8))

                # Footer
                elements.append(Spacer(1, 24))
                elements.append(Paragraph(
                    "This report was generated by the Single-Cell Intelligence Agent, "
                    "part of the HCLS AI Factory platform. For research and clinical decision "
                    "support only — not a standalone diagnostic.", footer_style))

                doc_pdf.build(elements)
                content = buffer.getvalue().decode("latin-1")
            except ImportError:
                # Fallback to formatted text if reportlab not available
                lines = [f"{'='*60}", f"  {title}", f"{'='*60}", ""]
                if request.patient_id:
                    lines.append(f"Patient ID: {request.patient_id}")
                if request.encounter_id:
                    lines.append(f"Encounter ID: {request.encounter_id}")
                lines.append(f"Report Type: {request.report_type}")
                lines.append(f"Generated: {now}")
                lines.append("")
                for key, value in request.data.items():
                    lines.append(f"--- {key.replace('_', ' ').title()} ---")
                    if isinstance(value, list):
                        for item in value:
                            lines.append(f"  * {item}")
                    elif isinstance(value, dict):
                        for k, v in value.items():
                            lines.append(f"  {k}: {v}")
                    else:
                        lines.append(f"  {value}")
                    lines.append("")
                lines.append("Generated by Single-Cell Intelligence Agent — HCLS AI Factory")
                content = "\n".join(lines)

        else:  # markdown
            header = _generate_markdown_header(title, request.patient_id, request.encounter_id)
            if request.report_type == "cell_type_annotation":
                body = _cell_type_annotation_markdown(request.data)
            elif request.report_type == "tme_profile":
                body = _tme_profile_markdown(request.data)
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
                "agent": "single-cell-intelligence-agent",
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
            {"id": "markdown", "name": "Markdown", "extension": ".md", "mime": "text/markdown", "description": "Human-readable single-cell analysis report"},
            {"id": "json", "name": "JSON", "extension": ".json", "mime": "application/json", "description": "Structured data export"},
            {"id": "pdf", "name": "PDF", "extension": ".pdf", "mime": "application/pdf", "description": "Printable single-cell report"},
            {"id": "fhir", "name": "FHIR R4", "extension": ".json", "mime": "application/fhir+json", "description": "HL7 FHIR R4 DiagnosticReport resource"},
        ],
        "report_types": [
            "cell_type_annotation",
            "tme_profile",
            "drug_response",
            "subclonal_analysis",
            "spatial_niche",
            "trajectory",
            "ligand_receptor",
            "biomarker",
            "cart_validation",
            "treatment_monitoring",
            "general",
        ],
    }
