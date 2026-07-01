"""Neurology report generation and export routes.

Provides endpoints for generating structured neurology reports in
multiple formats: Markdown, JSON, PDF, and FHIR R4 DiagnosticReport.
Supports stroke triage reports, cognitive assessments, epilepsy
classifications, tumor grading summaries, and clinical scale reports.

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
    """Request to generate a neurology report."""
    report_type: str = Field(
        ...,
        description=(
            "Type: stroke_triage | cognitive_assessment | epilepsy_classification | "
            "tumor_grading | ms_assessment | parkinsons_assessment | headache_report | "
            "neuromuscular_evaluation | scale_summary | general"
        ),
    )
    format: str = Field("markdown", pattern="^(markdown|json|pdf|fhir)$")
    patient_id: Optional[str] = None
    encounter_id: Optional[str] = None
    title: Optional[str] = None
    data: dict = Field(default={}, description="Report payload (scale results, triage data, etc.)")
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
        "**Agent:** Neurology Intelligence Agent v1.0.0",
    ]
    if patient_id:
        lines.append(f"**Patient ID:** {patient_id}")
    if encounter_id:
        lines.append(f"**Encounter ID:** {encounter_id}")
    lines.extend(["", "---", ""])
    return "\n".join(lines)


def _stroke_triage_markdown(data: dict) -> str:
    """Format stroke triage results as markdown."""
    lines = [
        "## Stroke Triage Assessment",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
        f"| NIHSS Score | **{data.get('nihss_score', 'N/A')}** |",
        f"| NIHSS Severity | **{data.get('nihss_severity', 'N/A')}** |",
        f"| tPA Eligible | **{data.get('tpa_eligible', 'N/A')}** |",
        f"| Thrombectomy Eligible | **{data.get('thrombectomy_eligible', 'N/A')}** |",
        f"| Urgency | **{data.get('urgency_level', 'N/A')}** |",
        "",
    ]

    recs = data.get("recommendations", [])
    if recs:
        lines.append("## Recommendations")
        lines.append("")
        for rec in recs:
            lines.append(f"- {rec}")
        lines.append("")

    guidelines = data.get("guidelines_cited", [])
    if guidelines:
        lines.append("## Guidelines Cited")
        lines.append("")
        for gl in guidelines:
            lines.append(f"- {gl}")
        lines.append("")

    return "\n".join(lines)


def _scale_summary_markdown(data: dict) -> str:
    """Format clinical scale results as markdown."""
    lines = [
        "## Clinical Scale Assessment",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
        f"| Scale | **{data.get('scale_name', 'N/A')}** |",
        f"| Score | **{data.get('total_score', 'N/A')}** / {data.get('max_score', 'N/A')} |",
        f"| Severity | **{data.get('severity_category', 'N/A')}** |",
        "",
        f"**Interpretation:** {data.get('interpretation', 'N/A')}",
        "",
    ]

    recs = data.get("recommendations", [])
    if recs:
        lines.append("## Recommendations")
        lines.append("")
        for rec in recs:
            lines.append(f"- {rec}")
        lines.append("")

    return "\n".join(lines)


def _cognitive_assessment_markdown(data: dict) -> str:
    """Format cognitive/dementia assessment as markdown."""
    lines = [
        "## Cognitive Assessment",
        "",
    ]
    if data.get("moca_score") is not None:
        lines.extend([
            f"**MoCA Score:** {data.get('moca_score')} / 30",
            f"**Interpretation:** {data.get('moca_interpretation', 'N/A')}",
            "",
        ])

    diff = data.get("differential_diagnosis", [])
    if diff:
        lines.extend([
            "### Differential Diagnosis",
            "",
            "| Diagnosis | Likelihood | Rationale |",
            "|-----------|------------|-----------|",
        ])
        for d in diff:
            if isinstance(d, dict):
                lines.append(f"| {d.get('diagnosis', 'N/A')} | {d.get('likelihood', 'N/A')} | {d.get('rationale', '')} |")
        lines.append("")

    workup = data.get("recommended_workup", [])
    if workup:
        lines.append("### Recommended Workup")
        lines.append("")
        for w in workup:
            lines.append(f"- {w}")
        lines.append("")

    return "\n".join(lines)


def _generate_fhir_diagnostic_report(data: dict, title: str, patient_id: Optional[str]) -> dict:
    """Generate a FHIR R4 DiagnosticReport resource."""
    now = datetime.now(timezone.utc).isoformat()
    resource = {
        "resourceType": "DiagnosticReport",
        "id": str(uuid.uuid4()),
        "status": "final",
        "category": [{"coding": [{"system": "http://loinc.org", "code": "LP7839-6", "display": "Neurology"}]}],
        "code": {"text": title},
        "effectiveDateTime": now,
        "issued": now,
        "conclusion": data.get("interpretation", data.get("summary", "")),
        "meta": {
            "lastUpdated": now,
            "source": "neurology-intelligence-agent",
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
    """Generate a formatted neurology report."""
    report_id = str(uuid.uuid4())[:12]
    now = datetime.now(timezone.utc).isoformat()
    title = request.title or f"Neurology {request.report_type.replace('_', ' ').title()}"

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
            content = f"[PDF generation placeholder] Title: {title} | Data keys: {list(request.data.keys())}"

        else:  # markdown
            header = _generate_markdown_header(title, request.patient_id, request.encounter_id)
            if request.report_type == "stroke_triage":
                body = _stroke_triage_markdown(request.data)
            elif request.report_type in ("cognitive_assessment", "dementia_evaluation"):
                body = _cognitive_assessment_markdown(request.data)
            elif request.report_type == "scale_summary":
                body = _scale_summary_markdown(request.data)
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
                "agent": "neurology-intelligence-agent",
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
            {"id": "markdown", "name": "Markdown", "extension": ".md", "mime": "text/markdown", "description": "Human-readable neurology report"},
            {"id": "json", "name": "JSON", "extension": ".json", "mime": "application/json", "description": "Structured data export"},
            {"id": "pdf", "name": "PDF", "extension": ".pdf", "mime": "application/pdf", "description": "Printable neurology report"},
            {"id": "fhir", "name": "FHIR R4", "extension": ".json", "mime": "application/fhir+json", "description": "HL7 FHIR R4 DiagnosticReport resource"},
        ],
        "report_types": [
            "stroke_triage",
            "cognitive_assessment",
            "epilepsy_classification",
            "tumor_grading",
            "ms_assessment",
            "parkinsons_assessment",
            "headache_report",
            "neuromuscular_evaluation",
            "scale_summary",
            "general",
        ],
    }
