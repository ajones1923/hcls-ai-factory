"""Report generation routes.

Provides endpoints to generate unified clinical reports for a given
patient, in multiple formats (PDF, Markdown, JSON).

Author: Adam Jones
Date: February 2026
"""

from __future__ import annotations

import io
import logging
import time

from fastapi import APIRouter, HTTPException

logger = logging.getLogger(__name__)
from fastapi.responses import JSONResponse, PlainTextResponse, StreamingResponse
from pydantic import BaseModel

from src.metrics import record_report_generated, record_pipeline_stage

router = APIRouter(prefix="/api", tags=["reports"])


# ── Response schema ──────────────────────────────────────────────────

class ReportMeta(BaseModel):
    """Metadata about a generated report."""

    patient_id: str
    format: str
    generated_at: str
    sections: int = 0
    processing_time_ms: float = 0.0


# ── Helpers ──────────────────────────────────────────────────────────

def _build_report_data(patient_id: str) -> dict:
    """Assemble a report payload for the given patient.

    In a full deployment this would query patient-specific variant data,
    RAG evidence, and drug-discovery results.  For now it returns a
    structured placeholder that downstream formatters can consume.
    """
    return {
        "patient_id": patient_id,
        "title": f"CAR-T Intelligence Report -- Patient {patient_id}",
        "sections": [
            {
                "heading": "Variant Summary",
                "body": "Variant analysis results will be populated from the genomics pipeline.",
            },
            {
                "heading": "Literature Evidence",
                "body": "Cross-collection RAG evidence will be synthesised here.",
            },
            {
                "heading": "Clinical Trial Matches",
                "body": "Matched clinical trials from ClinicalTrials.gov data.",
            },
            {
                "heading": "Recommended Constructs",
                "body": "CAR construct recommendations based on target antigen profile.",
            },
            {
                "heading": "Safety Profile",
                "body": "Toxicity and safety considerations from real-world evidence.",
            },
        ],
    }


def _render_markdown(data: dict) -> str:
    """Render the report dict as Markdown text."""
    lines = [f"# {data['title']}", ""]
    for section in data.get("sections", []):
        lines.append(f"## {section['heading']}")
        lines.append("")
        lines.append(section["body"])
        lines.append("")
    return "\n".join(lines)


def _render_pdf_bytes(data: dict) -> bytes:
    """Render the report dict as a simple PDF.

    Uses reportlab if available; falls back to a plain-text stub wrapped
    in a PDF container.
    """
    try:
        from reportlab.lib.pagesizes import letter
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
        from reportlab.lib.styles import getSampleStyleSheet

        buf = io.BytesIO()
        doc = SimpleDocTemplate(buf, pagesize=letter)
        styles = getSampleStyleSheet()
        story = [Paragraph(data["title"], styles["Title"]), Spacer(1, 20)]
        for section in data.get("sections", []):
            story.append(Paragraph(section["heading"], styles["Heading2"]))
            story.append(Spacer(1, 8))
            story.append(Paragraph(section["body"], styles["BodyText"]))
            story.append(Spacer(1, 14))
        doc.build(story)
        return buf.getvalue()
    except ImportError:
        # Fallback: return Markdown content as bytes
        return _render_markdown(data).encode("utf-8")


# ── Endpoints ────────────────────────────────────────────────────────

@router.get("/reports/{patient_id}")
async def generate_report(patient_id: str):
    """Generate a unified clinical report (defaults to JSON)."""
    t0 = time.perf_counter()

    try:
        data = _build_report_data(patient_id)
        elapsed_ms = (time.perf_counter() - t0) * 1000
        record_report_generated("json")
        record_pipeline_stage("report_generation", elapsed_ms / 1000)

        return JSONResponse(
            content={
                "meta": {
                    "patient_id": patient_id,
                    "format": "json",
                    "sections": len(data.get("sections", [])),
                    "processing_time_ms": round(elapsed_ms, 1),
                },
                "report": data,
            }
        )
    except Exception as exc:
        logger.error(f"Report generation failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")


@router.get("/reports/{patient_id}/{fmt}")
async def generate_report_format(patient_id: str, fmt: str):
    """Generate a report in a specific format (pdf, markdown, json)."""
    fmt = fmt.lower()
    if fmt not in ("pdf", "markdown", "json"):
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported format '{fmt}'. Use pdf, markdown, or json.",
        )

    t0 = time.perf_counter()

    try:
        data = _build_report_data(patient_id)
        elapsed_ms = (time.perf_counter() - t0) * 1000

        record_report_generated(fmt)
        record_pipeline_stage("report_generation", elapsed_ms / 1000)

        if fmt == "json":
            return JSONResponse(content={"meta": {"patient_id": patient_id, "format": fmt}, "report": data})

        if fmt == "markdown":
            md_text = _render_markdown(data)
            return PlainTextResponse(content=md_text, media_type="text/markdown")

        if fmt == "pdf":
            pdf_bytes = _render_pdf_bytes(data)
            return StreamingResponse(
                io.BytesIO(pdf_bytes),
                media_type="application/pdf",
                headers={
                    "Content-Disposition": f'attachment; filename="report_{patient_id}.pdf"'
                },
            )

    except HTTPException:
        raise
    except Exception as exc:
        logger.error(f"Report generation failed: {exc}")
        raise HTTPException(status_code=500, detail="Internal processing error")
