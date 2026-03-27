"""Export Pharmacogenomics Intelligence Agent query results to Markdown, JSON, PDF, and FHIR R4.

Provides five public functions:
  - export_markdown() -- human-readable report with evidence tables, alerts, and drug interaction matrix
  - export_json()     -- machine-readable structured data using Pydantic serialization
  - export_pdf()      -- styled PDF report via reportlab Platypus (includes PGx Passport format)
  - export_fhir_r4()  -- FHIR R4 DiagnosticReport Bundle with PGx-specific Observations (LOINC 69548-6)

PGx-specific additions vs CAR-T:
  - Alerts section in all exports (CRITICAL/WARNING/INFO with gene, drug, phenotype, recommendation)
  - Drug interaction matrix in markdown/PDF
  - FHIR Genomics profile with PGx-specific LOINC codes

Author: Adam Jones
Date: March 2026
"""

import io
import json
import re
import uuid
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from .models import ComparativeResult, CrossCollectionResult, PGxAlert, SearchHit


VERSION = "1.0.0"


# ===================================================================
# PUBLIC API
# ===================================================================


def generate_filename(extension: str) -> str:
    """Generate a timestamped filename for export.

    Args:
        extension: File extension without dot (e.g. "md", "json")

    Returns:
        Filename like pgx_query_20260312T143025Z.md
    """
    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"pgx_query_{ts}.{extension}"


def export_markdown(
    query: str,
    response_text: str,
    evidence: Optional[CrossCollectionResult] = None,
    comp_result: Optional[ComparativeResult] = None,
    filters_applied: Optional[dict] = None,
    alerts: Optional[List[PGxAlert]] = None,
) -> str:
    """Export a query result as a Markdown report.

    Args:
        query: The user's original question
        response_text: The LLM-generated response
        evidence: CrossCollectionResult (for standard queries)
        comp_result: ComparativeResult (for comparative queries)
        filters_applied: Dict of sidebar filters that were active
        alerts: List of PGxAlert clinical decision support alerts

    Returns:
        Complete Markdown report as a string
    """
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    filters_str = _format_filters(filters_applied)

    lines = [
        "# Pharmacogenomics Intelligence Report",
        "",
        f"**Query:** {query}",
        f"**Generated:** {timestamp}",
        f"**Filters:** {filters_str}",
        "",
    ]

    # Alerts section
    if alerts:
        lines.extend([
            "---",
            "",
            "## Clinical Alerts",
            "",
        ])
        for alert in alerts:
            level = alert.alert_level.value.upper()
            icon = {"CRITICAL": "!!!", "WARNING": "!!", "INFO": "i"}.get(level, "i")
            lines.append(
                f"- **[{level}]** **{alert.gene}** / **{alert.drug}** "
                f"({alert.phenotype}): {alert.recommendation}"
            )
        lines.append("")

    lines.extend([
        "---",
        "",
        "## Response",
        "",
        response_text,
        "",
        "---",
        "",
    ])

    # Evidence section
    if comp_result and comp_result.total_hits > 0:
        lines.append("## Evidence Sources (Comparative)")
        lines.append("")
        lines.append(f"### {comp_result.entity_a}")
        lines.append("")
        lines.extend(_format_evidence_section(comp_result.evidence_a))
        lines.append("")
        lines.append(f"### {comp_result.entity_b}")
        lines.append("")
        lines.extend(_format_evidence_section(comp_result.evidence_b))

        if comp_result.comparison_context:
            lines.append("")
            lines.append("## Knowledge Graph Context")
            lines.append("")
            lines.append(comp_result.comparison_context)

        # Search metrics
        lines.extend(["", "---", "", "## Search Metrics", ""])
        lines.append("| Metric | Value |")
        lines.append("|--------|-------|")
        lines.append(f"| Total Results | {comp_result.total_hits} |")
        lines.append(f"| {comp_result.entity_a} Results | {comp_result.evidence_a.hit_count} |")
        lines.append(f"| {comp_result.entity_b} Results | {comp_result.evidence_b.hit_count} |")
        lines.append(f"| Search Time | {comp_result.total_search_time_ms:.0f}ms |")

    elif evidence and evidence.hit_count > 0:
        lines.append("## Evidence Sources")
        lines.append("")
        lines.extend(_format_evidence_section(evidence))

        if evidence.knowledge_context:
            lines.append("")
            lines.append("## Knowledge Graph Context")
            lines.append("")
            lines.append(evidence.knowledge_context)

        # Search metrics
        lines.extend(["", "---", "", "## Search Metrics", ""])
        lines.append("| Metric | Value |")
        lines.append("|--------|-------|")
        lines.append(f"| Total Results | {evidence.hit_count} |")
        lines.append(f"| Collections Searched | {evidence.total_collections_searched} |")
        lines.append(f"| Search Time | {evidence.search_time_ms:.0f}ms |")

    # Footer
    lines.extend(["", "---", ""])
    lines.append(f"*Generated by HCLS AI Factory -- Pharmacogenomics Intelligence Agent v{VERSION}*")
    lines.append("")

    return "\n".join(lines)


def export_json(
    query: str,
    response_text: str,
    evidence: Optional[CrossCollectionResult] = None,
    comp_result: Optional[ComparativeResult] = None,
    alerts: Optional[List[PGxAlert]] = None,
) -> str:
    """Export a query result as structured JSON.

    Uses Pydantic .model_dump() for proper serialization of evidence models.

    Args:
        query: The user's original question
        response_text: The LLM-generated response
        evidence: CrossCollectionResult (for standard queries)
        comp_result: ComparativeResult (for comparative queries)
        alerts: List of PGxAlert clinical decision support alerts

    Returns:
        Pretty-printed JSON string
    """
    is_comparative = comp_result is not None and comp_result.total_hits > 0
    timestamp = datetime.now(timezone.utc).isoformat()

    data: Dict[str, Any] = {
        "report_type": "pgx_intelligence_query",
        "version": VERSION,
        "generated_at": timestamp,
        "query": query,
        "response": response_text,
        "is_comparative": is_comparative,
    }

    # Alerts
    if alerts:
        data["alerts"] = [
            {
                "alert_level": a.alert_level.value,
                "gene": a.gene,
                "drug": a.drug,
                "phenotype": a.phenotype,
                "recommendation": a.recommendation,
                "evidence_pmids": a.evidence_pmids,
            }
            for a in alerts
        ]

    if is_comparative:
        data["comparative"] = {
            "entity_a": comp_result.entity_a,
            "entity_b": comp_result.entity_b,
            "evidence_a": comp_result.evidence_a.model_dump(),
            "evidence_b": comp_result.evidence_b.model_dump(),
            "comparison_context": comp_result.comparison_context,
        }
        data["search_metrics"] = {
            "total_results": comp_result.total_hits,
            "entity_a_results": comp_result.evidence_a.hit_count,
            "entity_b_results": comp_result.evidence_b.hit_count,
            "search_time_ms": round(comp_result.total_search_time_ms, 1),
        }
    elif evidence:
        data["evidence"] = evidence.model_dump()
        data["search_metrics"] = {
            "total_results": evidence.hit_count,
            "collections_searched": evidence.total_collections_searched,
            "search_time_ms": round(evidence.search_time_ms, 1),
        }

    return json.dumps(data, indent=2, default=str)


# ===================================================================
# PRIVATE HELPERS
# ===================================================================


def _format_filters(filters_applied: Optional[dict]) -> str:
    """Format the sidebar filters for display."""
    if not filters_applied:
        return "None"
    parts = []
    for key, value in filters_applied.items():
        if value and value not in ("All Genes", "All Drugs"):
            parts.append(f"{key}: {value}")
    return ", ".join(parts) if parts else "None"


def _format_citation_link(collection: str, record_id: str) -> str:
    """Format a clickable citation link."""
    if collection in ("Evidence", "ClinicalEvidence", "pgx_clinical_evidence") and record_id.isdigit():
        return f"[PMID {record_id}](https://pubmed.ncbi.nlm.nih.gov/{record_id}/)"
    if collection in ("ClinicalTrial", "pgx_clinical_trials") and record_id.upper().startswith("NCT"):
        return f"[{record_id}](https://clinicaltrials.gov/study/{record_id})"
    if collection in ("FDALabel", "pgx_fda_labels"):
        return f"[FDA: {record_id}](https://www.accessdata.fda.gov/scripts/cder/daf/)"
    return record_id


def _format_evidence_section(evidence: CrossCollectionResult) -> list:
    """Format all evidence for a CrossCollectionResult, grouped by collection."""
    lines: list = []
    by_coll = evidence.hits_by_collection()
    for coll_name, hits in by_coll.items():
        lines.extend(_format_evidence_table(hits, coll_name))
        lines.append("")
    return lines


def _format_evidence_table(hits: list, collection_name: str) -> list:
    """Format a Markdown table for hits from a single collection.

    Uses collection-specific columns to surface the most relevant metadata.
    """
    lines = [f"### {collection_name} ({len(hits)} results)", ""]

    if collection_name in ("pgx_gene_reference", "GeneReference"):
        lines.append("| # | ID | Score | Gene | Allele | Function | Activity Score |")
        lines.append("|---|-----|-------|------|--------|----------|----------------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            gene = m.get("gene", "")
            allele = m.get("star_allele", "")
            func = m.get("function_status", "")[:30]
            ascore = m.get("activity_score", "")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {gene} | {allele} | {func} | {ascore} |")

    elif collection_name in ("pgx_drug_guidelines", "DrugGuideline"):
        lines.append("| # | ID | Score | Gene | Drug | Phenotype | Action | Alert |")
        lines.append("|---|-----|-------|------|------|-----------|--------|-------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            gene = m.get("gene", "")
            drug = m.get("drug", "")[:25]
            pheno = m.get("phenotype", "")[:20]
            action = m.get("clinical_action", "")
            alert = m.get("alert_level", "")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {gene} | {drug} | {pheno} | {action} | {alert} |")

    elif collection_name in ("pgx_drug_interactions", "DrugInteraction"):
        lines.append("| # | ID | Score | Drug | Gene | Type | Effect |")
        lines.append("|---|-----|-------|------|------|------|--------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            drug = m.get("drug", "")[:25]
            gene = m.get("gene", "")
            itype = m.get("interaction_type", "")
            effect = m.get("effect_description", "")[:40]
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {drug} | {gene} | {itype} | {effect} |")

    elif collection_name in ("pgx_hla_hypersensitivity", "HLA", "HLAHypersensitivity"):
        lines.append("| # | ID | Score | HLA Allele | Drug | Reaction | Severity | Screening |")
        lines.append("|---|-----|-------|------------|------|----------|----------|-----------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            hla = m.get("hla_allele", "")
            drug = m.get("drug", "")[:20]
            reaction = m.get("reaction_type", "")[:20]
            severity = m.get("severity", "")
            screening = "Required" if m.get("screening_mandatory") else "Recommended"
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {hla} | {drug} | {reaction} | {severity} | {screening} |")

    elif collection_name in ("pgx_phenoconversion", "Phenoconversion"):
        lines.append("| # | ID | Score | Enzyme | Precipitant | Strength | Effect |")
        lines.append("|---|-----|-------|--------|-------------|----------|--------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            enzyme = m.get("affected_enzyme", "")
            precip = m.get("precipitant_drug", "")[:25]
            strength = m.get("interaction_type", "")
            effect = m.get("effect_on_phenotype", "")[:35]
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {enzyme} | {precip} | {strength} | {effect} |")

    elif collection_name in ("pgx_dosing_algorithms", "DosingAlgorithm"):
        lines.append("| # | ID | Score | Drug | Genes | Algorithm | Validation |")
        lines.append("|---|-----|-------|------|-------|-----------|------------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            drug = m.get("drug", "")[:20]
            genes = m.get("genes_involved", "")[:15]
            algo = m.get("algorithm_name", "")[:25]
            valid = m.get("accuracy_metrics", "")[:20]
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {drug} | {genes} | {algo} | {valid} |")

    elif collection_name in ("pgx_clinical_evidence", "Evidence", "ClinicalEvidence"):
        lines.append("| # | ID | Score | Source | Title | Gene-Drug | Year |")
        lines.append("|---|-----|-------|--------|-------|-----------|------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            link = _format_citation_link(collection_name, hit.id)
            title = m.get("title", "")[:50]
            gene_drug = f"{m.get('gene', '')}-{m.get('drug', '')}" if m.get("gene") else ""
            year = m.get("year", "")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {link} | {title} | {gene_drug} | {year} |")

    elif collection_name in ("pgx_population_data", "PopulationData"):
        lines.append("| # | ID | Score | Gene | Allele | Population | Frequency | N |")
        lines.append("|---|-----|-------|------|--------|------------|-----------|---|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            gene = m.get("gene", "")
            allele = m.get("star_allele", "")
            pop = m.get("population", "")[:15]
            freq = m.get("allele_frequency", "")
            n = m.get("sample_size", "")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {gene} | {allele} | {pop} | {freq} | {n} |")

    elif collection_name in ("pgx_fda_labels", "FDALabel"):
        lines.append("| # | ID | Score | Drug | Gene | Section | Type | Requirement |")
        lines.append("|---|-----|-------|------|------|---------|------|-------------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            drug = m.get("drug", "")[:20]
            gene = m.get("gene", "")
            section = m.get("labeling_section", "")[:20]
            ltype = m.get("label_type", "")[:20]
            req = m.get("requirement_level", "")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {drug} | {gene} | {section} | {ltype} | {req} |")

    elif collection_name in ("pgx_drug_alternatives", "DrugAlternative"):
        lines.append("| # | ID | Score | Primary Drug | Gene | Phenotype | Alternative | Evidence |")
        lines.append("|---|-----|-------|-------------|------|-----------|-------------|----------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            primary = m.get("primary_drug", "")[:20]
            gene = m.get("gene", "")
            pheno = m.get("phenotype", "")[:15]
            alt = m.get("alternative_drug", "")[:20]
            evid = m.get("evidence_level", "")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {primary} | {gene} | {pheno} | {alt} | {evid} |")

    elif collection_name in ("pgx_clinical_trials", "Trial", "PGxClinicalTrial"):
        lines.append("| # | NCT ID | Score | Source | Phase | Status | Gene-Drug |")
        lines.append("|---|--------|-------|--------|-------|--------|-----------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            link = _format_citation_link(collection_name, hit.id)
            phase = m.get("phase", "")
            status = m.get("status", "")
            gene_drug = f"{m.get('gene', '')}-{m.get('drug', '')}" if m.get("gene") else ""
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {link} | {phase} | {status} | {gene_drug} |")

    elif collection_name in ("genomic_evidence", "Genomic"):
        lines.append("| # | ID | Score | Gene | Consequence | Impact | Clinical Significance |")
        lines.append("|---|-----|-------|------|-------------|--------|----------------------|")
        for i, hit in enumerate(hits[:10], 1):
            m = hit.metadata
            gene = m.get("gene", "")
            consequence = m.get("consequence", "")[:25]
            impact = m.get("impact", "")
            clin_sig = m.get("clinical_significance", "")[:25]
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {gene} | {consequence} | {impact} | {clin_sig} |")

    else:
        # Generic fallback
        lines.append("| # | ID | Score | Text |")
        lines.append("|---|-----|-------|------|")
        for i, hit in enumerate(hits[:10], 1):
            text = hit.text[:80].replace("|", "\\|")
            lines.append(f"| {i} | {hit.id} | {hit.score:.3f} | {text} |")

    return lines


# ===================================================================
# PDF EXPORT
# ===================================================================

try:
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.lib.units import inch
    from reportlab.platypus import (
        BaseDocTemplate,
        CondPageBreak,
        Frame,
        HRFlowable,
        KeepTogether,
        NextPageTemplate,
        PageBreak,
        PageTemplate,
        Paragraph,
        Spacer,
        Table,
        TableStyle,
    )
    _HAS_REPORTLAB = True
except ImportError:
    _HAS_REPORTLAB = False

# -- PDF color palette (only available if reportlab installed) --
if _HAS_REPORTLAB:
    _NVIDIA_GREEN = colors.HexColor("#76B900")
    _NVIDIA_GREEN_DARK = colors.HexColor("#5A8F00")
    _DARK_BG = colors.HexColor("#1B1B2F")
    _DARK_BG2 = colors.HexColor("#162447")
    _TABLE_ALT = colors.HexColor("#F7F9FC")
    _TABLE_BORDER = colors.HexColor("#E2E8F0")
    _LIGHT_GRAY = colors.HexColor("#64748B")
    _MEDIUM_GRAY = colors.HexColor("#94A3B8")
    _TEXT_PRIMARY = colors.HexColor("#1E293B")
    _TEXT_SECONDARY = colors.HexColor("#475569")
    _INFO_BG = colors.HexColor("#F0FDF4")
    _INFO_BORDER = colors.HexColor("#BBF7D0")
    _SECTION_BG = colors.HexColor("#F8FAFC")
    _WHITE = colors.white

# Alert level colors
_ALERT_COLORS = {
    "critical": colors.HexColor("#DC2626"),  # red
    "warning": colors.HexColor("#D97706"),   # amber
    "info": colors.HexColor("#2563EB"),      # blue
}

# Collection accent colors for evidence table headers
_COLLECTION_COLORS = {
    "pgx_gene_reference":       colors.HexColor("#059669"),  # emerald
    "pgx_drug_guidelines":      colors.HexColor("#2563EB"),  # blue
    "pgx_drug_interactions":    colors.HexColor("#DC2626"),  # red
    "pgx_hla_hypersensitivity": colors.HexColor("#E11D48"),  # rose
    "pgx_phenoconversion":      colors.HexColor("#7C3AED"),  # violet
    "pgx_dosing_algorithms":    colors.HexColor("#0891B2"),  # cyan
    "pgx_clinical_evidence":    colors.HexColor("#4F46E5"),  # indigo
    "pgx_population_data":      colors.HexColor("#CA8A04"),  # yellow
    "pgx_clinical_trials":      colors.HexColor("#9333EA"),  # purple
    "pgx_fda_labels":           colors.HexColor("#0D9488"),  # teal
    "pgx_drug_alternatives":    colors.HexColor("#D97706"),  # amber
    "pgx_patient_profiles":     colors.HexColor("#059669"),  # emerald
    "pgx_implementation":       colors.HexColor("#64748B"),  # slate
    "pgx_education":            colors.HexColor("#2563EB"),  # blue
    "genomic_evidence":         colors.HexColor("#9333EA"),  # purple
}

_PAGE_W, _PAGE_H = letter
_MARGIN_L = 0.65 * inch
_MARGIN_R = 0.65 * inch
_MARGIN_T = 0.85 * inch
_MARGIN_B = 0.75 * inch
_CONTENT_W = _PAGE_W - _MARGIN_L - _MARGIN_R


# -- Page decoration callbacks --

def _first_page(canvas, doc):
    """Draw branded header banner on the first page."""
    canvas.saveState()
    _h = 62
    canvas.setFillColor(_DARK_BG)
    canvas.rect(0, _PAGE_H - _h, _PAGE_W, _h, fill=1, stroke=0)
    canvas.setFillColor(_NVIDIA_GREEN)
    canvas.rect(0, _PAGE_H - _h - 3, _PAGE_W, 3, fill=1, stroke=0)
    canvas.setFillColor(_WHITE)
    canvas.setFont("Helvetica-Bold", 18)
    canvas.drawString(_MARGIN_L, _PAGE_H - 42, "Pharmacogenomics Intelligence Report")
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(_NVIDIA_GREEN)
    canvas.drawRightString(_PAGE_W - _MARGIN_R, _PAGE_H - 30, f"v{VERSION}")
    canvas.setFillColor(_MEDIUM_GRAY)
    canvas.drawRightString(_PAGE_W - _MARGIN_R, _PAGE_H - 42, "HCLS AI Factory")
    _draw_footer(canvas, doc)
    canvas.restoreState()


def _later_pages(canvas, doc):
    """Draw running header + footer on continuation pages."""
    canvas.saveState()
    _h = 28
    canvas.setFillColor(_DARK_BG)
    canvas.rect(0, _PAGE_H - _h, _PAGE_W, _h, fill=1, stroke=0)
    canvas.setFillColor(_NVIDIA_GREEN)
    canvas.rect(0, _PAGE_H - _h - 2, _PAGE_W, 2, fill=1, stroke=0)
    canvas.setFillColor(_WHITE)
    canvas.setFont("Helvetica-Bold", 9)
    canvas.drawString(_MARGIN_L, _PAGE_H - 19, "Pharmacogenomics Intelligence Report")
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(_MEDIUM_GRAY)
    canvas.drawRightString(_PAGE_W - _MARGIN_R, _PAGE_H - 19, "HCLS AI Factory")
    _draw_footer(canvas, doc)
    canvas.restoreState()


def _draw_footer(canvas, doc):
    """Branded footer with green bar, version, and page number."""
    y = _MARGIN_B - 20
    canvas.setStrokeColor(_NVIDIA_GREEN)
    canvas.setLineWidth(1.5)
    canvas.line(_MARGIN_L, y + 14, _PAGE_W - _MARGIN_R, y + 14)
    canvas.setFont("Helvetica", 7)
    canvas.setFillColor(_LIGHT_GRAY)
    canvas.drawString(
        _MARGIN_L, y,
        f"HCLS AI Factory  |  Pharmacogenomics Intelligence Agent v{VERSION}",
    )
    canvas.drawRightString(_PAGE_W - _MARGIN_R, y, f"Page {doc.page}")


# -- Style builder --

def _build_pdf_styles() -> dict:
    """Build custom ParagraphStyles for the PDF report."""
    base = getSampleStyleSheet()
    return {
        "h2": ParagraphStyle(
            "PDFH2", parent=base["Heading2"],
            fontSize=14, leading=18,
            textColor=_DARK_BG, fontName="Helvetica-Bold",
            spaceBefore=16, spaceAfter=8,
            borderPadding=(0, 0, 2, 0),
            keepWithNext=1,
        ),
        "h2_green": ParagraphStyle(
            "PDFH2G", parent=base["Heading2"],
            fontSize=13, leading=17,
            textColor=_NVIDIA_GREEN_DARK, fontName="Helvetica-Bold",
            spaceBefore=14, spaceAfter=6,
            keepWithNext=1,
        ),
        "h3": ParagraphStyle(
            "PDFH3", parent=base["Heading3"],
            fontSize=11, leading=14,
            textColor=_TEXT_PRIMARY, fontName="Helvetica-Bold",
            spaceBefore=10, spaceAfter=4,
            keepWithNext=1,
        ),
        "body": ParagraphStyle(
            "PDFBody", parent=base["BodyText"],
            fontSize=10, leading=15, spaceAfter=7,
            textColor=_TEXT_PRIMARY,
        ),
        "body_sm": ParagraphStyle(
            "PDFBodySm", parent=base["BodyText"],
            fontSize=9, leading=13, spaceAfter=5,
            textColor=_TEXT_SECONDARY,
        ),
        "meta_label": ParagraphStyle(
            "PDFMetaL", parent=base["BodyText"],
            fontSize=8, leading=11,
            textColor=_LIGHT_GRAY, fontName="Helvetica-Bold",
            spaceAfter=1,
        ),
        "meta_value": ParagraphStyle(
            "PDFMetaV", parent=base["BodyText"],
            fontSize=10, leading=14,
            textColor=_TEXT_PRIMARY,
            spaceAfter=6,
        ),
        "footer": ParagraphStyle(
            "PDFFooter", parent=base["BodyText"],
            fontSize=8, textColor=_LIGHT_GRAY, alignment=TA_CENTER,
            spaceBefore=12,
        ),
        "summary_num": ParagraphStyle(
            "SummaryNum", parent=base["BodyText"],
            fontSize=18, leading=22, fontName="Helvetica-Bold",
            textColor=_NVIDIA_GREEN_DARK, alignment=TA_CENTER,
        ),
        "summary_label": ParagraphStyle(
            "SummaryLabel", parent=base["BodyText"],
            fontSize=7, leading=10,
            textColor=_LIGHT_GRAY, alignment=TA_CENTER,
        ),
        "alert_text": ParagraphStyle(
            "AlertText", parent=base["BodyText"],
            fontSize=9, leading=13,
            textColor=_TEXT_PRIMARY,
            spaceAfter=4,
        ),
    }


# -- Utility helpers --

def _pdf_escape(text: str) -> str:
    """Escape text for reportlab Paragraph XML."""
    return (text
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;"))


def _trunc(text, max_len: int = 50) -> str:
    """Truncate and escape text for PDF table cells."""
    s = str(text)
    if len(s) > max_len:
        s = s[:max_len - 1] + "\u2026"
    return _pdf_escape(s)


# -- Alerts card for PDF --

def _build_alerts_card(alerts: List[PGxAlert], styles: dict) -> list:
    """Build a styled alerts section for the PDF."""
    if not alerts:
        return []

    flowables = [
        Paragraph("Clinical Alerts", styles["h2"]),
    ]

    for alert in alerts:
        level = alert.alert_level.value
        accent = _ALERT_COLORS.get(level, _ALERT_COLORS["info"])
        accent_hex = (f"#{int(accent.red*255):02x}{int(accent.green*255):02x}"
                      f"{int(accent.blue*255):02x}") if hasattr(accent, 'red') else "#2563EB"

        alert_text = (
            f'<font color="{accent_hex}"><b>[{level.upper()}]</b></font> '
            f'<b>{_pdf_escape(alert.gene)}</b> / <b>{_pdf_escape(alert.drug)}</b> '
            f'({_pdf_escape(alert.phenotype)}): '
            f'{_pdf_escape(alert.recommendation)}'
        )
        flowables.append(Paragraph(alert_text, styles["alert_text"]))

    flowables.append(Spacer(1, 10))
    return flowables


# -- Query info card --

def _build_query_card(query: str, timestamp: str, filters_str: str,
                      styles: dict) -> list:
    """Build a styled info card for query metadata."""
    card_data = [
        [Paragraph('<font color="#64748B"><b>QUERY</b></font>', styles["meta_label"]),
         Paragraph('<font color="#64748B"><b>GENERATED</b></font>', styles["meta_label"]),
         Paragraph('<font color="#64748B"><b>FILTERS</b></font>', styles["meta_label"])],
        [Paragraph(_pdf_escape(query), styles["meta_value"]),
         Paragraph(timestamp, styles["meta_value"]),
         Paragraph(_pdf_escape(filters_str), styles["meta_value"])],
    ]
    card = Table(card_data, colWidths=[_CONTENT_W * 0.50, _CONTENT_W * 0.25,
                                       _CONTENT_W * 0.25])
    card.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), _INFO_BG),
        ("BOX", (0, 0), (-1, -1), 0.75, _INFO_BORDER),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 10),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LINEBELOW", (0, 0), (-1, 0), 0, _INFO_BG),
    ]))
    return [card, Spacer(1, 14)]


# -- Summary card --

def _build_summary_card(evidence=None, comp_result=None, alerts=None) -> list:
    """Build a compact summary card showing key search metrics."""
    styles = _build_pdf_styles()

    if comp_result and comp_result.total_hits > 0:
        total = comp_result.total_hits
        colls = "Comparative"
        time_ms = comp_result.total_search_time_ms
    elif evidence and evidence.hit_count > 0:
        total = evidence.hit_count
        colls = str(evidence.total_collections_searched)
        time_ms = evidence.search_time_ms
    else:
        return []

    alert_count = len(alerts) if alerts else 0

    cells = [
        [Paragraph(f"<b>{total}</b>", styles["summary_num"]),
         Paragraph(f"<b>{colls}</b>", styles["summary_num"]),
         Paragraph(f"<b>{time_ms:.0f}ms</b>", styles["summary_num"]),
         Paragraph(f"<b>{alert_count}</b>", styles["summary_num"])],
        [Paragraph("Total Results", styles["summary_label"]),
         Paragraph("Collections", styles["summary_label"]),
         Paragraph("Search Time", styles["summary_label"]),
         Paragraph("Alerts", styles["summary_label"])],
    ]
    w = _CONTENT_W / 4
    tbl = Table(cells, colWidths=[w, w, w, w])
    tbl.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), _SECTION_BG),
        ("BOX", (0, 0), (-1, -1), 0.75, _TABLE_BORDER),
        ("LINEAFTER", (0, 0), (0, -1), 0.5, _TABLE_BORDER),
        ("LINEAFTER", (1, 0), (1, -1), 0.5, _TABLE_BORDER),
        ("LINEAFTER", (2, 0), (2, -1), 0.5, _TABLE_BORDER),
        ("TOPPADDING", (0, 0), (-1, 0), 10),
        ("BOTTOMPADDING", (0, -1), (-1, -1), 10),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
    ]))
    return [tbl, Spacer(1, 14)]


# -- Evidence PDF table --

def _build_pdf_evidence_table(hits: list, collection_name: str) -> list:
    """Build a styled evidence table for hits from a single collection."""
    styles = _build_pdf_styles()
    accent = _COLLECTION_COLORS.get(collection_name, _NVIDIA_GREEN)
    _accent_hex = (f"#{int(accent.red*255):02x}{int(accent.green*255):02x}"
                   f"{int(accent.blue*255):02x}") if hasattr(accent, 'red') else "#76B900"

    heading_text = (f'<font color="{_accent_hex}">'
                    f'\u275A</font>  '
                    f'{_pdf_escape(collection_name)} '
                    f'<font color="#94A3B8">({len(hits)} results)</font>')
    flowables = [Paragraph(heading_text, styles["h3"])]

    # Build generic rows (collection-specific formatting mirrors markdown)
    header = ["#", "ID", "Score", "Text"]
    rows = [header]
    for i, hit in enumerate(hits[:10], 1):
        rows.append([
            str(i), _trunc(hit.id, 24), f"{hit.score:.3f}",
            _trunc(hit.text, 70),
        ])

    body_style = ParagraphStyle("Cell", fontSize=8, leading=10,
                                textColor=_TEXT_PRIMARY)
    header_style = ParagraphStyle("CellH", fontSize=8, leading=10,
                                  textColor=_WHITE, fontName="Helvetica-Bold")
    para_rows = []
    for r_idx, row in enumerate(rows):
        style = header_style if r_idx == 0 else body_style
        para_rows.append([Paragraph(str(c), style) for c in row])

    table = Table(para_rows, repeatRows=1)
    cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), accent),
        ("TEXTCOLOR", (0, 0), (-1, 0), _WHITE),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("GRID", (0, 0), (-1, -1), 0.3, _TABLE_BORDER),
        ("BOX", (0, 0), (-1, -1), 0.75, accent),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
    ]
    for r_idx in range(1, len(para_rows)):
        if r_idx % 2 == 0:
            cmds.append(("BACKGROUND", (0, r_idx), (-1, r_idx), _TABLE_ALT))

    table.setStyle(TableStyle(cmds))
    flowables.append(table)
    flowables.append(Spacer(1, 12))

    if len(rows) <= 7:
        return [KeepTogether(flowables)]
    return flowables


# -- Markdown to flowables converter --

def _md_to_flowables(text: str, styles: dict) -> list:
    """Convert markdown-ish response text to reportlab flowables."""
    flowables = []
    normalized = text.strip()
    normalized = re.sub(r'(?<!\n)\n(#{1,3} )', r'\n\n\1', normalized)
    normalized = re.sub(r'(^#{1,3} [^\n]+)\n(?!\n)', r'\1\n\n', normalized,
                        flags=re.MULTILINE)
    normalized = re.sub(r'(?<!\n)\n(> )', r'\n\n\1', normalized)
    normalized = re.sub(r'(?<!\n)\n(-{3,})(?=\n|$)', r'\n\n\1', normalized)
    normalized = re.sub(r'(?<!\n)\n(\|)', r'\n\n\1', normalized)

    blocks = re.split(r'\n{2,}', normalized.strip())

    for block in blocks:
        block = block.strip()
        if not block:
            continue

        if re.match(r'^-{3,}$', block):
            flowables.append(Spacer(1, 4))
            flowables.append(HRFlowable(width="100%", thickness=0.5,
                                        color=_TABLE_BORDER))
            flowables.append(Spacer(1, 4))
            continue

        if block.startswith("### "):
            heading_line = block.split("\n")[0][4:]
            heading_line = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>',
                                  _pdf_escape(heading_line))
            flowables.append(Paragraph(heading_line, styles["h3"]))
            rest = "\n".join(block.split("\n")[1:]).strip()
            if rest:
                flowables.extend(_md_to_flowables(rest, styles))
            continue

        if block.startswith("## "):
            heading_line = block.split("\n")[0][3:]
            flowables.append(Paragraph(_pdf_escape(heading_line), styles["h2"]))
            rest = "\n".join(block.split("\n")[1:]).strip()
            if rest:
                flowables.extend(_md_to_flowables(rest, styles))
            continue

        if block.startswith("# "):
            heading_line = block.split("\n")[0][2:]
            flowables.append(Paragraph(f'<b>{_pdf_escape(heading_line)}</b>',
                                       styles["h2"]))
            rest = "\n".join(block.split("\n")[1:]).strip()
            if rest:
                flowables.extend(_md_to_flowables(rest, styles))
            continue

        if block.startswith("> "):
            quote_style = ParagraphStyle(
                "Quote", parent=styles["body"],
                leftIndent=18, borderLeftWidth=3,
                borderLeftColor=_NVIDIA_GREEN,
                borderPadding=(6, 0, 0, 8),
                textColor=_TEXT_SECONDARY,
                fontName="Helvetica-Oblique",
                fontSize=9, leading=13,
            )
            clean = re.sub(r'^> ?', '', block, flags=re.MULTILINE)
            clean = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', _pdf_escape(clean))
            clean = clean.replace("\n", "<br/>")
            flowables.append(Paragraph(clean, quote_style))
            continue

        escaped = _pdf_escape(block)
        escaped = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', escaped)
        escaped = re.sub(
            r'\[(.+?)\]\((.+?)\)',
            r'<a href="\2" color="#2563EB"><u>\1</u></a>',
            escaped,
        )
        lines = escaped.split("\n")
        processed = []
        for line in lines:
            if line.startswith("- "):
                processed.append(
                    f'<font color="#76B900">\u2022</font>  {line[2:]}')
            elif re.match(r'^\d+\.\s', line):
                num_match = re.match(r'^(\d+\.)\s(.*)', line)
                if num_match:
                    processed.append(
                        f'<font color="#76B900"><b>{num_match.group(1)}</b>'
                        f'</font>  {num_match.group(2)}')
                else:
                    processed.append(line)
            else:
                processed.append(line)
        escaped = "<br/>".join(processed)
        flowables.append(Paragraph(escaped, styles["body"]))

    return flowables


# -- Section divider --

def _section_divider() -> list:
    """A styled section divider with spacing."""
    return [
        Spacer(1, 8),
        HRFlowable(width="100%", thickness=0.75, color=_TABLE_BORDER,
                    spaceAfter=4, spaceBefore=4),
        Spacer(1, 4),
    ]


# ===================================================================
# MAIN PDF EXPORT
# ===================================================================

def export_pdf(
    query: str,
    response_text: str,
    evidence: Optional[CrossCollectionResult] = None,
    comp_result: Optional[ComparativeResult] = None,
    alerts: Optional[List[PGxAlert]] = None,
) -> bytes:
    """Export a query result as a professionally styled PDF report.

    Uses reportlab Platypus with branded NVIDIA/HCLS themed styling.
    Includes PGx-specific alerts section.

    Args:
        query: The user's original question
        response_text: The LLM-generated response
        evidence: CrossCollectionResult (for standard queries)
        comp_result: ComparativeResult (for comparative queries)
        alerts: List of PGxAlert clinical decision support alerts

    Returns:
        PDF file content as bytes

    Raises:
        ImportError: If reportlab is not installed.
    """
    if not _HAS_REPORTLAB:
        raise ImportError("reportlab is required for PDF export. Install with: pip install reportlab")

    buffer = io.BytesIO()

    frame_first = Frame(
        _MARGIN_L, _MARGIN_B,
        _CONTENT_W, _PAGE_H - _MARGIN_T - _MARGIN_B - 20,
        id="first",
    )
    frame_later = Frame(
        _MARGIN_L, _MARGIN_B,
        _CONTENT_W, _PAGE_H - _MARGIN_T - _MARGIN_B + 10,
        id="later",
    )
    doc = BaseDocTemplate(
        buffer, pagesize=letter,
        leftMargin=_MARGIN_L, rightMargin=_MARGIN_R,
        topMargin=_MARGIN_T, bottomMargin=_MARGIN_B,
        title="Pharmacogenomics Intelligence Report",
        author="HCLS AI Factory",
    )
    doc.addPageTemplates([
        PageTemplate(id="First", frames=[frame_first], onPage=_first_page),
        PageTemplate(id="Later", frames=[frame_later], onPage=_later_pages),
    ])

    styles = _build_pdf_styles()
    elements = []

    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    filters_str = "None"

    # -- Query info card --
    elements.append(Spacer(1, 6))
    elements.extend(_build_query_card(query, timestamp, filters_str, styles))

    # -- Evidence summary card --
    elements.extend(_build_summary_card(evidence, comp_result, alerts))

    # -- Alerts --
    if alerts:
        elements.extend(_build_alerts_card(alerts, styles))

    # Switch to continuation page template
    elements.append(NextPageTemplate("Later"))

    # -- Response --
    elements.append(Paragraph("Response", styles["h2"]))
    elements.extend(_md_to_flowables(response_text, styles))
    elements.extend(_section_divider())

    # -- Evidence tables --
    if comp_result and comp_result.total_hits > 0:
        elements.append(CondPageBreak(2.5 * inch))
        elements.append(Paragraph("Evidence Sources (Comparative)", styles["h2"]))

        elements.append(Paragraph(
            f'<font color="#2563EB">\u25B6</font>  {_pdf_escape(comp_result.entity_a)}',
            styles["h3"],
        ))
        for coll_name, hits in comp_result.evidence_a.hits_by_collection().items():
            elements.extend(_build_pdf_evidence_table(hits, coll_name))

        elements.append(Paragraph(
            f'<font color="#7C3AED">\u25B6</font>  {_pdf_escape(comp_result.entity_b)}',
            styles["h3"],
        ))
        for coll_name, hits in comp_result.evidence_b.hits_by_collection().items():
            elements.extend(_build_pdf_evidence_table(hits, coll_name))

    elif evidence and evidence.hit_count > 0:
        elements.append(CondPageBreak(2.5 * inch))
        section_heading = Paragraph("Evidence Sources", styles["h2"])
        colls = list(evidence.hits_by_collection().items())
        first_flowables = _build_pdf_evidence_table(colls[0][1], colls[0][0])
        elements.append(KeepTogether([section_heading] + first_flowables))
        for coll_name, hits in colls[1:]:
            elements.extend(_build_pdf_evidence_table(hits, coll_name))

    elements.append(Spacer(1, 20))
    doc.build(elements)
    return buffer.getvalue()


# ===================================================================
# FHIR R4 EXPORT
# ===================================================================

# LOINC codes for PGx-specific observations
_FHIR_LOINC_CODES: Dict[str, tuple] = {
    "pgx_report": ("51969-4", "Genetic analysis report"),
    "variant_assessment": ("69548-6", "Genetic variant assessment"),
    "metabolizer_status": ("53040-2", "Genetic variant assessment of predicted drug metabolism phenotype"),
    "hla_typing": ("13303-3", "HLA typing"),
    "drug_metabolism": ("55233-1", "Genetic analysis - Drug metabolism"),
}


def export_fhir_r4(
    query: str,
    response_text: str,
    evidence: Optional[CrossCollectionResult] = None,
    alerts: Optional[List[PGxAlert]] = None,
    patient_id: str = "anonymous",
) -> str:
    """Export a PGx query result as a FHIR R4 Bundle (JSON string).

    Creates a FHIR Bundle of type 'collection' containing:
      - Patient resource stub
      - DiagnosticReport with PGx analysis summary
      - Observation resources for each PGx alert (LOINC 69548-6: Genetic variant assessment)
      - Observation resources for key evidence findings

    Args:
        query: The user's original question
        response_text: The LLM-generated response
        evidence: CrossCollectionResult from search
        alerts: List of PGxAlert clinical decision support alerts
        patient_id: Patient identifier (default 'anonymous')

    Returns:
        FHIR R4 Bundle as a JSON string
    """
    now_iso = datetime.now(timezone.utc).isoformat()
    entries: List[Dict[str, Any]] = []

    # --- Patient resource ---
    patient_fullurl = f"urn:uuid:{uuid.uuid4()}"
    entries.append({
        "fullUrl": patient_fullurl,
        "resource": {
            "resourceType": "Patient",
            "id": patient_id,
            "identifier": [{
                "system": "urn:hcls-ai-factory:pgx-intelligence",
                "value": patient_id,
            }],
            "active": True,
        },
    })

    # --- PGx Alert Observations (LOINC 69548-6) ---
    observation_refs: List[str] = []

    if alerts:
        for alert in alerts:
            obs_fullurl = f"urn:uuid:{uuid.uuid4()}"
            observation_refs.append(obs_fullurl)

            # Map alert level to FHIR interpretation
            interpretation_code = {
                "critical": ("AA", "Critical abnormal"),
                "warning": ("A", "Abnormal"),
                "info": ("N", "Normal"),
            }.get(alert.alert_level.value, ("N", "Normal"))

            observation: Dict[str, Any] = {
                "fullUrl": obs_fullurl,
                "resource": {
                    "resourceType": "Observation",
                    "id": str(uuid.uuid4()),
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
                            "code": _FHIR_LOINC_CODES["variant_assessment"][0],
                            "display": _FHIR_LOINC_CODES["variant_assessment"][1],
                        }],
                        "text": f"PGx Alert: {alert.gene} / {alert.drug}",
                    },
                    "subject": {"reference": patient_fullurl},
                    "effectiveDateTime": now_iso,
                    "interpretation": [{
                        "coding": [{
                            "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
                            "code": interpretation_code[0],
                            "display": interpretation_code[1],
                        }],
                    }],
                    "valueString": alert.recommendation[:500],
                    "component": [
                        {
                            "code": {"text": "Gene"},
                            "valueString": alert.gene,
                        },
                        {
                            "code": {"text": "Drug"},
                            "valueString": alert.drug,
                        },
                        {
                            "code": {"text": "Phenotype"},
                            "valueString": alert.phenotype,
                        },
                        {
                            "code": {"text": "Alert Level"},
                            "valueString": alert.alert_level.value,
                        },
                    ],
                },
            }
            entries.append(observation)

    # --- Evidence Observations ---
    if evidence:
        for hit in evidence.hits[:20]:
            obs_fullurl = f"urn:uuid:{uuid.uuid4()}"
            observation_refs.append(obs_fullurl)

            observation = {
                "fullUrl": obs_fullurl,
                "resource": {
                    "resourceType": "Observation",
                    "id": str(uuid.uuid4()),
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
                            "code": _FHIR_LOINC_CODES["pgx_report"][0],
                            "display": _FHIR_LOINC_CODES["pgx_report"][1],
                        }],
                        "text": f"PGx Evidence: {hit.collection}",
                    },
                    "subject": {"reference": patient_fullurl},
                    "effectiveDateTime": now_iso,
                    "valueString": hit.text[:500],
                    "component": [],
                },
            }

            # Add PGx-specific metadata as components
            if hit.metadata.get("gene"):
                observation["resource"]["component"].append({
                    "code": {"text": "Gene"},
                    "valueString": hit.metadata["gene"],
                })
            if hit.metadata.get("drug"):
                observation["resource"]["component"].append({
                    "code": {"text": "Drug"},
                    "valueString": hit.metadata["drug"],
                })
            if hit.metadata.get("phenotype"):
                observation["resource"]["component"].append({
                    "code": {"text": "Phenotype"},
                    "valueString": hit.metadata["phenotype"],
                })
            if hit.metadata.get("hla_allele"):
                observation["resource"]["component"].append({
                    "code": {"text": "HLA Allele"},
                    "valueString": hit.metadata["hla_allele"],
                })
            if hit.metadata.get("clinical_action"):
                observation["resource"]["component"].append({
                    "code": {"text": "Clinical Action"},
                    "valueString": hit.metadata["clinical_action"],
                })

            # Relevance score
            observation["resource"]["component"].append({
                "code": {"text": "Relevance Score"},
                "valueQuantity": {
                    "value": hit.score,
                    "unit": "{score}",
                    "system": "http://unitsofmeasure.org",
                    "code": "{score}",
                },
            })

            entries.append(observation)

    # --- DiagnosticReport resource ---
    report_fullurl = f"urn:uuid:{uuid.uuid4()}"
    diagnostic_report: Dict[str, Any] = {
        "resourceType": "DiagnosticReport",
        "id": str(uuid.uuid4()),
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
                "code": _FHIR_LOINC_CODES["variant_assessment"][0],
                "display": _FHIR_LOINC_CODES["variant_assessment"][1],
            }],
            "text": "Pharmacogenomics Intelligence Agent Report",
        },
        "subject": {"reference": patient_fullurl},
        "effectiveDateTime": now_iso,
        "issued": now_iso,
        "performer": [{
            "display": "HCLS AI Factory - Pharmacogenomics Intelligence Agent",
        }],
        "result": [{"reference": ref} for ref in observation_refs],
        "conclusion": response_text[:2000],
    }

    # Add extension for query context
    diagnostic_report["extension"] = [{
        "url": "urn:hcls-ai-factory:pgx:query",
        "valueString": query,
    }]

    # Add alert summary as presentedForm
    if alerts:
        alert_summary = "; ".join(
            f"[{a.alert_level.value.upper()}] {a.gene}/{a.drug}: {a.recommendation[:100]}"
            for a in alerts
        )
        diagnostic_report["presentedForm"] = [{
            "contentType": "text/plain",
            "data": None,  # Would be base64 in production
            "title": "PGx Clinical Alerts Summary",
        }]

    entries.append({
        "fullUrl": report_fullurl,
        "resource": diagnostic_report,
    })

    # --- Build Bundle ---
    bundle: Dict[str, Any] = {
        "resourceType": "Bundle",
        "id": str(uuid.uuid4()),
        "type": "collection",
        "timestamp": now_iso,
        "meta": {
            "profile": [
                "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report",
            ],
        },
        "entry": entries,
    }

    return json.dumps(bundle, indent=2)
