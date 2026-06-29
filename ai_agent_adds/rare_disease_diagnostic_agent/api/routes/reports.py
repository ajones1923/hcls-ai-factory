"""Rare disease diagnostic report generation and export routes.

Provides endpoints for generating structured rare disease diagnostic reports
in multiple formats: Markdown, JSON, PDF, FHIR R4, and GA4GH Phenopacket.
Supports differential diagnosis reports, variant interpretation reports,
phenotype analysis, and therapeutic option summaries.

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
    """Request to generate a rare disease diagnostic report."""
    report_type: str = Field(
        ...,
        description=(
            "Type: differential_diagnosis | variant_interpretation | "
            "phenotype_analysis | therapeutic_summary | trial_eligibility | "
            "case_summary | gene_panel | natural_history"
        ),
    )
    format: str = Field("markdown", pattern="^(markdown|json|pdf|fhir|phenopacket)$")
    patient_id: Optional[str] = None
    disease_id: Optional[str] = None
    title: Optional[str] = None
    data: dict = Field(default={}, description="Report payload (diagnosis results, variant data, etc.)")
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

def _generate_markdown_header(title: str, disease_id: Optional[str] = None, patient_id: Optional[str] = None) -> str:
    """Standard markdown report header."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        f"# {title}",
        "",
        f"**Generated:** {now}",
        "**Agent:** Rare Disease Diagnostic Agent v1.0.0",
    ]
    if disease_id:
        lines.append(f"**Disease ID:** {disease_id}")
    if patient_id:
        lines.append(f"**Patient ID:** {patient_id}")
    lines.extend(["", "---", ""])
    return "\n".join(lines)


def _differential_diagnosis_markdown(data: dict) -> str:
    """Format differential diagnosis results as markdown."""
    lines = []

    phenotype_summary = data.get("phenotype_summary", "")
    if phenotype_summary:
        lines.extend([
            "## Phenotype Summary",
            "",
            phenotype_summary,
            "",
        ])

    differential = data.get("differential", [])
    if differential:
        lines.extend([
            "## Differential Diagnosis",
            "",
            "| Rank | Disease | ID | Confidence | Phenotype Overlap | Gene Match |",
            "|------|---------|-----|------------|-------------------|------------|",
        ])
        for i, dx in enumerate(differential, 1):
            if isinstance(dx, dict):
                lines.append(
                    f"| {i} | {dx.get('disease_name', 'N/A')} | "
                    f"{dx.get('disease_id', 'N/A')} | "
                    f"**{dx.get('confidence', 0):.1%}** | "
                    f"{dx.get('phenotype_overlap', 0):.1%} | "
                    f"{'Yes' if dx.get('gene_match') else 'No'} |"
                )
        lines.append("")

    recommendations = data.get("recommendations", [])
    if recommendations:
        lines.append("## Recommendations")
        lines.append("")
        for rec in recommendations:
            lines.append(f"- {rec}")
        lines.append("")

    return "\n".join(lines)


def _variant_interpretation_markdown(data: dict) -> str:
    """Format variant interpretation results as markdown."""
    lines = [
        "## Variant Classification",
        "",
        "| Field | Value |",
        "|-------|-------|",
        f"| Gene | **{data.get('gene', 'N/A')}** |",
        f"| Variant | **{data.get('variant', 'N/A')}** |",
        f"| Classification | **{data.get('classification', 'VUS')}** |",
        f"| Pathogenicity Score | **{data.get('pathogenicity_score', 'N/A')}** |",
        f"| ClinVar | {data.get('clinvar_significance', 'Not found')} |",
        f"| AlphaMissense | {data.get('alphamissense_score', 'N/A')} |",
        "",
    ]

    criteria = data.get("acmg_criteria", [])
    if criteria:
        lines.extend([
            "## ACMG/AMP Criteria Applied",
            "",
            "| Code | Strength | Met | Evidence |",
            "|------|----------|-----|----------|",
        ])
        for c in criteria:
            if isinstance(c, dict):
                met = "Yes" if c.get("met") else "No"
                lines.append(
                    f"| {c.get('code', '')} | {c.get('strength', '')} | "
                    f"{met} | {c.get('evidence', '')[:80]} |"
                )
        lines.append("")

    return "\n".join(lines)


def _generate_fhir_diagnostic_report(data: dict, title: str, patient_id: Optional[str]) -> dict:
    """Generate a FHIR R4 DiagnosticReport resource."""
    now = datetime.now(timezone.utc).isoformat()
    resource = {
        "resourceType": "DiagnosticReport",
        "id": str(uuid.uuid4()),
        "status": "final",
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "81247-9",
                "display": "Master HL7 genetic variant reporting panel",
            }],
            "text": title,
        },
        "issued": now,
        "conclusion": data.get("phenotype_summary", data.get("classification", "")),
        "meta": {
            "lastUpdated": now,
            "source": "rare-disease-diagnostic-agent",
        },
    }
    if patient_id:
        resource["subject"] = {"reference": f"Patient/{patient_id}"}
    return resource


def _generate_phenopacket(data: dict, title: str, patient_id: Optional[str]) -> dict:
    """Generate a GA4GH Phenopacket v2 resource."""
    now = datetime.now(timezone.utc).isoformat()
    phenopacket = {
        "id": str(uuid.uuid4()),
        "subject": {
            "id": patient_id or "unknown",
        },
        "phenotypicFeatures": [],
        "interpretations": [],
        "metaData": {
            "created": now,
            "createdBy": "rare-disease-diagnostic-agent",
            "phenopacketSchemaVersion": "2.0",
            "resources": [
                {
                    "id": "hp",
                    "name": "Human Phenotype Ontology",
                    "url": "http://purl.obolibrary.org/obo/hp.owl",
                    "version": "2024-04-26",
                    "namespacePrefix": "HP",
                    "iriPrefix": "http://purl.obolibrary.org/obo/HP_",
                },
            ],
        },
    }

    # Add phenotypes from data
    phenotypes = data.get("phenotypes", [])
    for p in phenotypes:
        if isinstance(p, dict):
            phenopacket["phenotypicFeatures"].append({
                "type": {
                    "id": p.get("id", ""),
                    "label": p.get("label", ""),
                },
            })
        elif isinstance(p, str):
            phenopacket["phenotypicFeatures"].append({
                "type": {"id": p, "label": p},
            })

    return phenopacket


# =====================================================================
# Endpoints
# =====================================================================

@router.post("/generate", response_model=ReportResponse)
async def generate_report(request: ReportRequest, req: Request):
    """Generate a formatted rare disease diagnostic report."""
    report_id = str(uuid.uuid4())[:12]
    now = datetime.now(timezone.utc).isoformat()
    title = request.title or f"Rare Disease {request.report_type.replace('_', ' ').title()}"

    try:
        if request.format == "fhir":
            fhir_resource = _generate_fhir_diagnostic_report(
                request.data, title, request.patient_id,
            )
            content = json.dumps(fhir_resource, indent=2)

        elif request.format == "phenopacket":
            phenopacket = _generate_phenopacket(
                request.data, title, request.patient_id,
            )
            content = json.dumps(phenopacket, indent=2)

        elif request.format == "json":
            content = json.dumps({
                "report_id": report_id,
                "title": title,
                "type": request.report_type,
                "generated": now,
                "disease_id": request.disease_id,
                "patient_id": request.patient_id,
                "data": request.data,
            }, indent=2)

        elif request.format == "pdf":
            # Generate real PDF using reportlab
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

                title_style = ParagraphStyle("RDTitle", parent=styles["Title"],
                                             textColor=navy, fontSize=18, spaceAfter=12)
                heading_style = ParagraphStyle("RDHeading", parent=styles["Heading2"],
                                               textColor=teal, fontSize=13, spaceAfter=8)
                body_style = ParagraphStyle("RDBody", parent=styles["Normal"],
                                            fontSize=10, leading=14, spaceAfter=6)
                footer_style = ParagraphStyle("RDFooter", parent=styles["Normal"],
                                              fontSize=8, textColor=HexColor("#999999"))

                elements = []
                elements.append(Paragraph(title, title_style))
                elements.append(Paragraph("Rare Disease Diagnostic Agent — HCLS AI Factory", heading_style))
                elements.append(Spacer(1, 12))

                if request.patient_id:
                    elements.append(Paragraph(f"<b>Patient ID:</b> {request.patient_id}", body_style))
                if request.disease_id:
                    elements.append(Paragraph(f"<b>Disease ID:</b> {request.disease_id}", body_style))
                elements.append(Paragraph(f"<b>Report Type:</b> {request.report_type}", body_style))
                elements.append(Paragraph(f"<b>Generated:</b> {now}", body_style))
                elements.append(Spacer(1, 12))

                for key, value in request.data.items():
                    elements.append(Paragraph(key.replace("_", " ").title(), heading_style))
                    if isinstance(value, list):
                        for item in value:
                            if isinstance(item, dict):
                                for k, v in item.items():
                                    elements.append(Paragraph(f"<b>{k}:</b> {v}", body_style))
                                elements.append(Spacer(1, 6))
                            else:
                                elements.append(Paragraph(f"• {item}", body_style))
                    elif isinstance(value, dict):
                        for k, v in value.items():
                            elements.append(Paragraph(f"<b>{k}:</b> {v}", body_style))
                    else:
                        elements.append(Paragraph(str(value), body_style))
                    elements.append(Spacer(1, 8))

                elements.append(Spacer(1, 24))
                elements.append(Paragraph(
                    "This report was generated by the Rare Disease Diagnostic Agent, "
                    "part of the HCLS AI Factory platform. For research and clinical decision "
                    "support only — not a standalone diagnostic.", footer_style))

                doc_pdf.build(elements)
                content = buffer.getvalue().decode("latin-1")
            except ImportError:
                # Fallback to formatted text if reportlab not available
                lines = [f"{'='*60}", f"  {title}", f"{'='*60}", ""]
                if request.patient_id:
                    lines.append(f"Patient ID: {request.patient_id}")
                if request.disease_id:
                    lines.append(f"Disease ID: {request.disease_id}")
                lines.append(f"Report Type: {request.report_type}")
                lines.append(f"Generated: {now}")
                lines.append("")
                for key, value in request.data.items():
                    lines.append(f"--- {key.replace('_', ' ').title()} ---")
                    if isinstance(value, list):
                        for item in value:
                            lines.append(f"  • {item}")
                    elif isinstance(value, dict):
                        for k, v in value.items():
                            lines.append(f"  {k}: {v}")
                    else:
                        lines.append(f"  {value}")
                    lines.append("")
                lines.append("Generated by Rare Disease Diagnostic Agent — HCLS AI Factory")
                content = "\n".join(lines)

        else:  # markdown
            header = _generate_markdown_header(title, request.disease_id, request.patient_id)
            if request.report_type == "differential_diagnosis":
                body = _differential_diagnosis_markdown(request.data)
            elif request.report_type == "variant_interpretation":
                body = _variant_interpretation_markdown(request.data)
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
                "agent": "rare-disease-diagnostic-agent",
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
            {"id": "markdown", "name": "Markdown", "extension": ".md", "mime": "text/markdown", "description": "Human-readable rare disease diagnostic report"},
            {"id": "json", "name": "JSON", "extension": ".json", "mime": "application/json", "description": "Structured data export"},
            {"id": "pdf", "name": "PDF", "extension": ".pdf", "mime": "application/pdf", "description": "Printable diagnostic report"},
            {"id": "fhir", "name": "FHIR R4", "extension": ".json", "mime": "application/fhir+json", "description": "HL7 FHIR R4 DiagnosticReport resource"},
            {"id": "phenopacket", "name": "GA4GH Phenopacket", "extension": ".json", "mime": "application/json", "description": "GA4GH Phenopacket v2 for rare disease data exchange"},
        ],
        "report_types": [
            "differential_diagnosis",
            "variant_interpretation",
            "phenotype_analysis",
            "therapeutic_summary",
            "trial_eligibility",
            "case_summary",
            "gene_panel",
            "natural_history",
        ],
    }
